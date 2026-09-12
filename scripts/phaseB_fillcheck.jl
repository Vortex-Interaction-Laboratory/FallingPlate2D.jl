using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
L = 32; sec = Ellipse(L, 4); ν = L/500; θ = π/2; x0 = SA[3L*1.0, 4L*1.0]
println("--- (B') plate moving +0.5 through still fluid, interior velocity filled ---")
V = SA[0.5, 0.0]; x₀ = x0
sim = Simulation((8L, 8L), (0.0, 0.0), L; U=1.0, ν, body=AutoBody(FP.sdf(sec), RigidMap(x0, θ; V)), T=Float64)
for n in 1:150
    Δt = min(sim.flow.Δt[end], 0.5); sim.flow.Δt[end] = Δt
    global x₀ = x₀ + Δt*V; sim.body = setmap(sim.body; x₀, V)
    measure!(sim); rigid_velocity_fill!(sim.flow, sim.body, x₀, V, 0.0); sim_step!(sim; remeasure=false)
end
j = round(Int, x₀[2] + 1.5); i0 = round(Int, x₀[1] + 1.5)
for di in -8:4:8; @printf("   %+3d  p=%9.4f  μ₀=%.2f\n", di, sim.flow.p[i0+di, j], sim.flow.μ₀[i0+di, j, 1]); end
F = -SVector{2}(WaterLily.total_force(sim))
@printf("   F_on_body=(%8.3f,%8.3f)  → C_D = %.2f  (expect Fx<0; fixed-plate C_D at U=1 was 1.67)\n", F..., -2F[1]/(0.25L))

println("\n--- constant rotation in still fluid, interior filled: expect M·ω < 0 ---")
for ωstar in (0.5, 2.0)
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.0, free=(false,false,false))
    s = fp.sim; ω = ωstar/fp.L; θr = fp.state.θ; tp = 0.0
    while sim_time(s) < 3.0
        Δt = min(s.flow.Δt[end], 0.5); s.flow.Δt[end] = Δt; θr += Δt*ω
        s.body = setmap(s.body; θ=θr, ω); measure!(s); rigid_velocity_fill!(s.flow, s.body, fp.x₀, SA[0.0,0.0], ω); sim_step!(s; remeasure=false)
        if sim_time(s) - tp ≥ 0.75
            Mp = -WaterLily.pressure_moment(fp.x₀, s)[1]; Mv = -WaterLily.viscous_moment(fp.x₀, s)[1]
            @printf("  ωc/U=%.1f t=%.2f  M_p=%9.4f  M_v=%9.4f  M/(U²c²)=%9.5f\n", ωstar, sim_time(s), Mp, Mv, (Mp+Mv)/fp.L^2); tp = sim_time(s)
        end
    end
end

println("\n--- 3-DOF coupled fall (T3), χ=4 ρ=2.7 Re=500, 12 c/U ---")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1)
run!(fp; tend=12.0, out="out/checkT3.csv", every=0.05, report=1.0)
summarise("out/checkT3.csv"; window=6.0, ν=fp.ν, L=fp.L)
