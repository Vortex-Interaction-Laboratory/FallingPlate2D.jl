using FallingPlate2D, WaterLily, StaticArrays, Test
FP = FallingPlate2D

@testset "sections" begin
    e = Ellipse(32, 4)
    @test area(e) ≈ π*16*4
    @test inertia(e) ≈ π/4*16*4*(16^2 + 4^2)
    @test added_mass(e) == (π*4^2, π*16^2, π/8*(16^2 - 4^2)^2)
    Ia, Ir = andersen_Istar(Ellipse(112, 14), 2.7)
    @test Ir ≈ 0.164 atol = 2e-3          # Andersen et al. (2005) table 1, β = 1/14
    p = RoundedPlate(32, 4)
    @test area(p) ≈ 4*16*4 - (4 - π)*16
    @test FP.sdf(e)(SA[0.0, 0.0], 0.0) < 0 && FP.sdf(e)(SA[20.0, 0.0], 0.0) > 0
    @test FP.sdf(p)(SA[0.0, 0.0], 0.0) ≈ -4 && FP.sdf(p)(SA[0.0, 8.0], 0.0) ≈ 4
end

@testset "rigid step from rest" begin
    d = PlateDynamics(Ellipse(32, 4), 2.7, 0.01)
    s = RigidState(π/2)                      # horizontal plate: normal along Y
    FP.step!(s, d, SA[0.0, 0.0], 0.0, 0.5)
    @test s.a[1] ≈ 0 atol = 1e-12
    @test s.a[2] ≈ -(d.m - d.A)*0.01/(d.m + d.m22)
    @test s.α == 0
end

@testset "frames" begin
    @test to_lab(to_grid(SA[1.0, 2.0])) == SA[1.0, 2.0]
    @test to_grid(SA[0.0, -1.0]) == SA[-1.0, 0.0]  # falling (−Y) is −x in the grid: relative flow +x
    @test abs(FP.chord_dir(π/2)[1]) ≈ 1 && abs(FP.chord_dir(π/2)[2]) < 1e-12   # θ = π/2 is a horizontal chord
    @test FP.chord_dir(1.0)' * FP.normal_dir(1.0) ≈ 0 atol = 1e-12
end

@testset "drag opposes relative motion" begin
    # held plate in a stream: force along the stream
    fp = falling_plate(; χ=4, ρ=2.7, Re=200, L=16, dom=(8, 8), φ₀=0.0, free=(false, false, false))
    fp.Ufar[] = SA[1.0, 0.0]
    for _ in 1:40; sim_step!(fp.sim; remeasure=false); end
    F, _ = FP.body_loads(fp)
    @test F[1] > 0
    # plate moving through still fluid with the interior velocity filled: force against the motion
    L = 16; sec = Ellipse(L, 4); V = SA[0.5, 0.0]; x₀ = SA[3L*1.0, 4L*1.0]
    sim = Simulation((8L, 8L), (0.0, 0.0), L; U=1.0, ν=L/200, body=AutoBody(FP.sdf(sec), RigidMap(x₀, π/2; V)), T=Float64)
    for _ in 1:60
        sim.flow.Δt[end] = min(sim.flow.Δt[end], 0.5); x₀ += sim.flow.Δt[end]*V
        sim.body = setmap(sim.body; x₀, V); measure!(sim); rigid_velocity_fill!(sim.flow, sim.body, x₀, V, 0.0)
        sim_step!(sim; remeasure=false)
    end
    @test -WaterLily.total_force(sim)[1] < 0
    @test maximum(abs, sim.flow.p) < 10          # no spurious interior pressure
end

@testset "coupled fall stays finite" begin
    fp = falling_plate(; χ=4, ρ=2.7, Re=200, L=16, dom=(8, 8), φ₀=0.1)
    @test run!(fp; tend=3.0, report=100.0) == :done
    @test fp.state.v[2] < -0.3                   # falling
    @test all(isfinite, fp.sim.flow.p)
end

@testset "Darcy: drag falls and seepage rises with Da" begin
    CD = Float64[]; useep = Float64[]
    for Da in (0.0, 1e-3, 1e-2)
        fp = falling_plate(; χ=4, ρ=2.7, Re=200, L=16, dom=(8, 8), φ₀=0.0, free=(false, false, false), Da)
        fp.Ufar[] = SA[1.0, 0.0]; fp.uexit .= 1.0
        for _ in 1:60
            fp.sim.flow.Δt[end] = min(fp.sim.flow.Δt[end], 0.5)
            prepare_flow!(fp); sim_step!(fp.sim; remeasure=false)
        end
        F, _ = FP.body_loads(fp); I = CartesianIndex(round.(Int, fp.x₀ .+ 1.5)...)
        push!(CD, 2F[1]/fp.L); push!(useep, fp.sim.flow.u[I, 1])
        @test all(isfinite, fp.sim.flow.p)
    end
    @test CD[1] > CD[2] > CD[3] > 0
    @test useep[1] ≈ 0 atol = 1e-6
    @test 0 < useep[2] < useep[3] < 1.1
end

@testset "gust profiles" begin
    for iso in (true, false)
        g = Gust(; dir=[0, 1], GR=1.0, t₀=15.0, T_G=2.0, isolated=iso)
        @test gust_velocity(g, 14.0) == 0 && gust_acceleration(g, 14.0) == 0
        @test gust_velocity(g, 16.0) ≈ (iso ? 1.0 : 0.5)        # isolated peaks at mid-period
        @test gust_velocity(g, 18.0) ≈ (iso ? 0.0 : 1.0)        # isolated returns the fluid to rest
        dt = 1e-4; ts = 14:dt:19
        @test sum(gust_acceleration(g, t) for t in ts)*dt ≈ (iso ? 0.0 : 1.0) atol = 1e-6
        @test sum(gust_velocity(g, t) for t in ts)*dt ≈ gust_displacement(g, 19.0) atol = 1e-4
    end
    @test gust_velocity(Gust(; GR=1.0, T_G=2.0, isolated=true), 17.0) ≈ 0 atol = 1e-12
    @test peak_acceleration(Gust(; GR=1.0, T_G=2.0, isolated=true)) ≈ π/2
end

@testset "gust does not displace a body that ignores the fluid" begin
    # With no fluid force the isolated gust must leave the earth-frame altitude
    # unchanged: the body keeps its velocity and the air returns to rest.
    g = Gust(; dir=[0, 1], GR=1.0, t₀=2.0, T_G=2.0, isolated=true)
    d = PlateDynamics(Ellipse(32, 10), 1.0, 0.04; A=0.0)
    Z = map((NOGUST, g)) do gu
        s = RigidState(π/2); Δt = 1e-3
        for n in 1:20_000
            t = n*Δt
            FP.step!(s, d, SA[0.0, 0.0], 0.0, Δt; gust=gust_acceleration(gu, t)*gu.dir)
        end
        s.X[2] + gust_displacement(gu, 20.0)
    end
    @test Z[2] ≈ Z[1] rtol = 1e-3
end

@testset "bottom-heavy plate is pendulum-stable" begin
    # The fluid force acts at the plate centre, the weight at the offset centre
    # of gravity: a tilted plate must feel a restoring moment.
    d = PlateDynamics(Ellipse(32, 10), 1.0, 0.04; e=1.16*32, A=0.0)
    for φ in (0.2, -0.2)
        s = RigidState(π/2 + φ)
        FP.step!(s, d, SA[0.0, 1.0], 0.0, 1e-6)     # upward drag, no moment about the centre
        @test sign(s.α) == -sign(φ)                  # rotates the chord back towards horizontal
    end
end

@testset "gravity modulation reaches the requested extremes" begin
    g, T_G, t₀, A = 1.25, 2.0, 5.0, 0.5
    gm = gravity_modulation(; amplitude=A, T_G, t₀, g)
    ratio(t) = 1 + gust_acceleration(gm, t)/g
    @test ratio(t₀ - 0.1) ≈ 1                       # quiescent before onset
    @test ratio(t₀ + T_G/4) ≈ 1 - A   atol=1e-12    # down to 0.5 g
    @test ratio(t₀ + T_G/2) ≈ 1       atol=1e-12    # back to g
    @test ratio(t₀ + 3T_G/4) ≈ 1 + A  atol=1e-12    # up to 1.5 g
    @test ratio(t₀ + T_G + 0.1) ≈ 1                 # and back, no net impulse
    @test gust_velocity(gm, t₀ + T_G) ≈ 0 atol=1e-12
end

@testset "checkpoint restart reproduces the uninterrupted run" begin
    mk() = falling_plate(; χ=10.0, m=0.8, e=1.16, I_O=1.15, Da=1e-4, Ga=150.0,
                         L=16, dom=(5, 5), xbody=2.5, shape=:plate)
    f = tempname() * ".jld2"
    a = mk(); run!(a; tend=0.3, report=1e9); save_checkpoint(a, f)
    run!(a; tend=0.6, report=1e9)
    b = mk(); load_checkpoint!(b, f); run!(b; tend=0.6, report=1e9)
    @test a.state.X ≈ b.state.X
    @test a.state.v ≈ b.state.v
    @test a.state.θ ≈ b.state.θ
    @test a.state.ω ≈ b.state.ω
    @test sum(a.sim.flow.Δt) ≈ sum(b.sim.flow.Δt)
    rm(f, force=true)
end
