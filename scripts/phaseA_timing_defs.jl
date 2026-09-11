# Phase A: verify the WaterLily API and measure wall time on this laptop.
#
# Freely rotating ellipse (translation locked), Float64, AutoBody SDF + RigidMap,
# implicit added mass, run for `tend` convective times on two grids.
# Run:  julia -t 4 --project=. scripts/phaseA_timing.jl
using WaterLily, StaticArrays, Printf

# Approximate signed distance to an ellipse with semi-axes a, b (first-order
# normalised implicit function; adequate for BDIM's 2-cell kernel).
function ellipse_sdf(a, b)
    (x, t) -> begin
        f = sqrt((x[1]/a)^2 + (x[2]/b)^2) - 1
        g = sqrt((x[1]/a^2)^2 + (x[2]/b^2)^2)
        f / g
    end
end

function make_sim(L; χ=4, Re=200, dom=(8, 16), T=Float64, U=1.0, θ₀=0.1)
    a, b = T(L/2), T(L/(2χ))
    x₀ = SA{T}[dom[1]*L/2, dom[2]*L/4]        # upstream quarter of the domain
    body = AutoBody(ellipse_sdf(a, b), RigidMap(x₀, T(θ₀)))
    sim = Simulation((dom[1]*L, dom[2]*L), (zero(T), T(U)), L; ν=U*L/Re, body, T)
    return sim, x₀, a, b
end

# Advance to `tend` convective times, integrating the rotation ODE with the
# added moment of inertia treated implicitly (as in TwoD_FreeRotatingEllipse.jl).
function run!(sim, x₀, a, b; tend, ρ_s=2.0)
    T = eltype(sim.flow.p)
    Is = ρ_s*π*a^3*b/4                  # solid moment of inertia (ρ_f = 1)
    Ia = π/8*(a^2 - b^2)^2              # added moment of inertia (Lamb)
    θ, ω, α = sim.body.map.θ, zero(T), zero(T)
    steps = 0; next_report = 1.0; t0 = time()
    while sim_time(sim) < tend
        Δt = sim.flow.Δt[end]
        M = -WaterLily.pressure_moment(x₀, sim.flow, sim.body)[1]
        α = (M + Ia*α)/(Is + Ia)
        ω += Δt*α; θ += Δt*ω
        sim.body = setmap(sim.body; θ=T(θ), ω=T(ω))
        sim_step!(sim; remeasure=true)
        steps += 1
        if sim_time(sim) >= next_report
            @printf("    t=%5.1f  steps=%6d  Δt=%.3f  pois_iter=%d  θ=%7.3f  ω=%7.3f  wall=%6.1f s\n",
                sim_time(sim), steps, Δt, last(sim.pois.n), θ, ω*sim.L/sim.U, time()-t0); flush(stdout)
            next_report += 1.0
        end
    end
    F = WaterLily.total_force(sim)
    return (; steps, θ, ω=ω*sim.L/sim.U, F=F ./ (sim.L*sim.U^2))
end

