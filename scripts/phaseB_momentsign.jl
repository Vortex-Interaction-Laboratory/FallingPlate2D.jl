# Sign and magnitude of the fluid moment on a plate rotating at constant ω in still fluid.
using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
for ωstar in (0.5, 2.0)
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.0, free=(false,false,false))
    sim = fp.sim; ω = ωstar*fp.U/fp.L; θ = fp.state.θ
    @printf("\nω c/U = %.1f (tip speed %.2f U). Expect M·ω < 0 after the start-up transient.\n", ωstar, ω*fp.L/2)
    t_prev = 0.0
    while sim_time(sim) < 3.0
        Δt = min(sim.flow.Δt[end], 0.5); sim.flow.Δt[end] = Δt
        θ += Δt*ω
        sim.body = setmap(sim.body; θ, ω)
        sim_step!(sim; remeasure=true)
        if sim_time(sim) - t_prev ≥ 0.5
            Mp = -WaterLily.pressure_moment(fp.x₀, sim)[1]; Mv = -WaterLily.viscous_moment(fp.x₀, sim)[1]
            @printf("  t=%.2f  M_p=%9.4f  M_v=%9.4f  M/(U²c²)=%9.5f  Δt=%.3f\n", sim_time(sim), Mp, Mv, (Mp+Mv)/fp.L^2, Δt)
            t_prev = sim_time(sim)
        end
    end
end
