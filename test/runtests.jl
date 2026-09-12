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
