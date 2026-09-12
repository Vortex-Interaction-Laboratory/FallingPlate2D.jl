# WaterLily simulation in the plate frame, coupled to the rigid-body dynamics.

struct FallingPlateSim{S<:Section}
    sim::Simulation
    state::RigidState
    dyn::PlateDynamics
    section::S
    Ufar::Ref{SVector{2,Float64}}    # far-field velocity, grid frame (= −v_body)
    afluid::Ref{SVector{2,Float64}}  # body force on the fluid, grid frame (= −a_body)
    x₀::SVector{2,Float64}           # plate centre in the grid
    L::Float64
    U::Float64
    ν::Float64
    g::Float64
    Δt_max::Float64
    free::NTuple{3,Bool}             # (X, Y, θ) degrees of freedom released
    uexit::Vector{Float64}           # convective-exit x-velocity profile (all y indices)
    convective_exit::Bool
    Da::Float64                      # Darcy number K/c²; 0 = impermeable
    gust::Gust
end

"""
    falling_plate(; χ, ρ, Re, L, dom, xbody, shape, φ₀, T, Δt_max)
    falling_plate(; χ, m, Ga, Da, e, I_O, L, gust, ...)

Build a plate of chord `L` cells and chord-to-thickness ratio `χ`, released
from rest at chord angle `φ₀` from the horizontal, in a domain of `dom .* L`
cells with the plate centre `xbody` chords from the inflow, the lateral
direction periodic and a convective exit downstream. `Da > 0` makes the plate
permeable (Darcy–Brinkman, see `darcy.jl`).

Two conventions are available for the weight.

  - **Buoyant plate** (default): density ratio `ρ`, `Re = c·U/ν` with
    `U = sqrt(2(ρ − 1) f g) = 1`, weight `(m − A)g` — the convention of
    Andersen et al. (2005).
  - **Permeable plate**, used when `m` is given: the convention of Pandi
    et al. (JFM 2026) and Bose & Viola (2023). Lengths are scaled on the plate
    length, velocities on the gravitational velocity `u_g ≡ (m̂ĝ/ρ̂l̂²)^(1/2)`, so
    that the nondimensional gravity is `g = 1/m`; buoyancy is null because the
    solid volume of a filamentous body is negligible; and `Ga ≡ u_g·c/ν` (pass
    `Ga`, or equivalently `Re`, both being `1/ν`). The centre of gravity may
    be offset from the plate centre by `e` plate lengths along the body normal,
    with `I_O` the moment of inertia about the plate centre. `gust` applies the
    gust of `gust.jl` as a transient change of the gravitational force.

The convective exit is implemented here through the boundary function rather
than with WaterLily's `exitBC=true`: WaterLily recomputes and flux-corrects its
exit only in the predictor, so when the body force changes within a step (as
it does with a coupled body acceleration) the corrector and projection see an
incompatible flux, and the error feeds back through the plate force and grows.
`update_exit!` recomputes the profile once per step from the start-of-step
velocity, then every `BC!` call in the step sees the same, flux-consistent
values. `exit=:waterlily` uses WaterLily's `exitBC=true`, `exit=:dirichlet` a
plain Dirichlet exit (both for diagnostics).
"""
function falling_plate(; χ=14.0, ρ=2.7, Re=2090.0, L=80, dom=(8, 8), xbody=3,
                        shape=:ellipse, φ₀=0.1, T=Float64, Δt_max=0.5, free=(true, true, true),
                        perdir=(2,), exit=:own, Da=0.0,
                        m=nothing, Ga=nothing, e=0.0, I_O=nothing, gust=NOGUST)
    sec = shape == :ellipse ? Ellipse(L, χ) : RoundedPlate(L, χ)
    U = 1.0
    Ga === nothing || (Re = Ga)
    ν = U*L/Re
    if m === nothing                                   # buoyant plate (Andersen convention)
        g = U^2/(2(ρ - 1)*thickness(sec))
        dyn = PlateDynamics(sec, ρ, g)
    else                                               # permeable plate (Pandi convention)
        g = U^2/(m*L)                                  # g = 1/m, in grid units
        mg, eg = m*L^2, e*L                            # mass per unit span and offset, grid units
        IO = I_O === nothing ? m*(1 + 1/χ^2)/12 + m*e^2 : I_O
        IG = (IO - m*e^2)*L^4                          # inertia about the centre of gravity
        IG > 0 || throw(ArgumentError("I_O = \$IO is below m·e² = \$(m*e^2)"))
        dyn = PlateDynamics(sec, 1.0, g; e=eg, A=0.0, Is=IG, m=mg)
    end
    θ₀ = π/2 + φ₀
    state = RigidState(θ₀)
    Ufar = Ref(SA[0.0, 0.0]); afluid = Ref(SA[0.0, 0.0])
    x₀ = SA{T}[xbody*L, dom[2]*L/2]
    body = AutoBody(sdf(sec), RigidMap(x₀, T(θ₀)))
    N1, N2 = dom[1]*L + 2, dom[2]*L + 2           # array sizes including ghosts
    uexit = zeros(N2)
    xexit = N1 - 2.5                                # loc(1, I)[1] = N1 - 2 on the exit face
    exit in (:own, :waterlily, :dirichlet) || throw(ArgumentError("exit must be :own, :waterlily or :dirichlet"))
    uBC = exit == :own ? ((i, x, t) -> (i == 1 && x[1] > xexit) ? uexit[round(Int, x[2] + 1.5)] : Ufar[][i]) :
                         ((i, x, t) -> Ufar[][i])
    gfun = (i, x, t) -> afluid[][i]
    sim = Simulation((dom[1]*L, dom[2]*L), uBC, L; U, ν, g=gfun, body, T,
                     perdir, exitBC=(exit == :waterlily))
    FallingPlateSim(sim, state, dyn, sec, Ufar, afluid, x₀, Float64(L), U, ν, g, Δt_max, free,
                    uexit, exit == :own, Float64(Da), gust)
end

"""
    update_exit!(fp)

First-order upwind convection of the start-of-step exit profile at the mean
far-field speed, then a uniform shift so that the exit flux equals the inflow
flux at the end of the step. Call after the far field has been updated and
before `sim_step!`.
"""
function update_exit!(fp::FallingPlateSim)
    fp.convective_exit || return
    u = fp.sim.flow.u; N1, N2 = size(fp.sim.flow.p)
    Δt = fp.sim.flow.Δt[end]; Ux = fp.Ufar[][1]
    Uc = min(max(Ux, 0.0)*Δt, 1.0)                  # upwind CFL number, clipped
    ue = fp.uexit
    @inbounds for j in 2:N2-1
        ue[j] = u[N1, j, 1] - Uc*(u[N1, j, 1] - u[N1-1, j, 1])
    end
    m = sum(@view ue[2:N2-1])/(N2 - 2)
    @inbounds for j in 2:N2-1
        ue[j] += Ux - m
    end
    ue[1] = ue[N2-1]; ue[N2] = ue[2]
    nothing
end

"""
    rigid_velocity_fill!(flow, body, x₀, v, ω)

Set the imposed body velocity `flow.V` on every face lying fully inside the
body (`μ₀ == 0`) to the rigid-body value `v + ω × (x − x₀)`.

WaterLily fills `V` only within `2 + ϵ` cells of the surface and leaves it zero
deeper inside. For a moving body that makes the interior velocity field
discontinuous, with a divergence where the pressure equation has no
coefficients; the unremovable residual produces a spurious interior pressure
of O(10⁴) that leaks into the surface band and corrupts the force and moment
(wrong sign, wrong magnitude). Filling the interior with the exact rigid
velocity makes it divergence-free and restores a smooth interior pressure.
Call after `measure!` and before the flow step.
"""
function rigid_velocity_fill!(flow, body, x₀::SVector{2}, v::SVector{2}, ω::Real)
    μ₀, V = flow.μ₀, flow.V
    T = eltype(V); N = size(flow.p); t = sum(flow.Δt)
    @inbounds for i in 1:2, I in CartesianIndices(N)
        if μ₀[I, i] == 0
            x = WaterLily.loc(i, I, T)
            WaterLily.sdf(body, x, t) < 0 || continue     # skip the zero-mobility boundary faces
            r = x - x₀
            V[I, i] = i == 1 ? v[1] - ω*r[2] : v[2] + ω*r[1]
        end
    end
    nothing
end

"""
    prepare_flow!(fp)

Measure the body, fill the interior rigid velocity, apply the Darcy mobility,
and refresh the Poisson operator: everything `sim_step!(remeasure=true)` would
do, plus the two corrections. Call after the body map has been updated.
"""
function prepare_flow!(fp::FallingPlateSim)
    sim = fp.sim
    measure!(sim.flow, sim.body; t=sum(sim.flow.Δt), ϵ=sim.ϵ)
    rigid_velocity_fill!(sim.flow, sim.body, fp.x₀, zero(fp.x₀), fp.state.ω)   # body centre fixed in the frame
    apply_darcy!(sim.flow, fp.Da, fp.L)
    WaterLily.update!(sim.pois)
    nothing
end

"Force on the body (grid frame) and moment about its centre from the current flow."
function body_loads(fp::FallingPlateSim)
    F = -SVector{2,Float64}(WaterLily.total_force(fp.sim))
    M = -WaterLily.total_moment(fp.x₀, fp.sim)[1]
    F, M
end

"""
    step!(fp::FallingPlateSim) -> (F_grid, M, Δt)

One coupled step: read the loads, advance the body, update the far-field
velocity, the frame body force and the body map, then advance the flow.
"""
function step!(fp::FallingPlateSim)
    sim, s = fp.sim, fp.state
    Δt = min(sim.flow.Δt[end], fp.Δt_max); sim.flow.Δt[end] = Δt
    ag = gust_acceleration(fp.gust, sim_time(sim))/fp.L*fp.gust.dir   # grid units
    F, M = body_loads(fp)
    step!(s, fp.dyn, to_lab(F), M, Δt; free=fp.free, gust=ag)
    fp.Ufar[] = -to_grid(centre_velocity(s, fp.dyn))
    fp.afluid[] = -to_grid(centre_acceleration(s, fp.dyn))
    update_exit!(fp)
    sim.body = setmap(sim.body; θ=s.θ, ω=s.ω)
    prepare_flow!(fp)
    sim_step!(sim; remeasure=false)
    F, M, Δt
end

const CSV_HEADER = "t,X,Y,phi,vX,vY,omega,FX,FY,M,dt,pois,Xg,Zearth,uG"

"""
    run!(fp; tend, out=nothing, every=0.05, report=1.0, ωmax=10, vmax=4, Δt_min=1e-3)

Advance to `tend` convective times `c/U`. Writes a CSV row every `every`
convective times to `out` (nondimensional: lengths by `c`, velocities by `U`,
forces by `U²c`, moment by `U²c²`, angular velocity by `U/c`), prints one line
per `report` convective times, and aborts if the state is not finite, if
`|ω| c/U > ωmax`, if `|v| > vmax·U` (a plate diving edge-on at several `U`
makes the realised `Re` meaningless and the CFL step tiny), or if
`Δt < Δt_min` cells. Returns `:done` or `:aborted`.
"""
function run!(fp::FallingPlateSim; tend, out=nothing, every=0.05, report=1.0,
              ωmax=10.0, vmax=4.0, Δt_min=1e-3)
    sim, s, L, U = fp.sim, fp.state, fp.L, fp.U
    io = out === nothing ? nothing : open(out, "w")
    io === nothing || println(io, CSV_HEADER)
    next_out, next_rep, t0, n, status = 0.0, report, time(), 0, :done
    while sim_time(sim) < tend
        F, M, Δt = step!(fp); n += 1
        t = sim_time(sim)
        ω⁺ = s.ω*L/U
        if io !== nothing && t ≥ next_out
            Fl = to_lab(F)
            dG = gust_displacement(fp.gust, t)
            @printf(io, "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.4f,%d,%.5f,%.5f,%.5f\n",
                t, s.X[1]/L, s.X[2]/L, chord_angle(s.θ), s.v[1]/U, s.v[2]/U, ω⁺,
                Fl[1]/(U^2*L), Fl[2]/(U^2*L), M/(U^2*L^2), Δt, last(sim.pois.n),
                s.X[1]/L + dG*fp.gust.dir[1], s.X[2]/L + dG*fp.gust.dir[2],
                gust_velocity(fp.gust, t))
            next_out += every
        end
        if t ≥ next_rep
            @printf("t=%6.2f  X/c=%7.2f Y/c=%8.2f  φ=%7.2f  vX=%6.3f vY=%6.3f  ωc/U=%6.3f  Δt=%.3f  %5.1f s\n",
                t, s.X[1]/L, s.X[2]/L, chord_angle(s.θ), s.v[1]/U, s.v[2]/U, ω⁺, Δt, time() - t0)
            flush(stdout); next_rep += report
        end
        bad = !all(isfinite, s.v) || !isfinite(s.ω) || abs(ω⁺) > ωmax || Δt < Δt_min ||
              sqrt(sum(abs2, s.v)) > vmax*U
        if bad
            u = sim.flow.u; N = size(sim.flow.p)
            mag = [sqrt(sum(abs2, @view u[I, :])) for I in CartesianIndices(N)]
            I = argmax(mag)
            @warn "aborting" t v=s.v ω=ω⁺ Δt umax=mag[I] at=Tuple(I) grid=N body=Tuple(round.(Int, fp.x₀))
            status = :aborted; break
        end
    end
    io === nothing || close(io)
    @printf("%s: %d steps, %.1f s wall, %.1f ms/step\n", status, n, time() - t0, 1e3*(time() - t0)/max(n, 1))
    status
end
