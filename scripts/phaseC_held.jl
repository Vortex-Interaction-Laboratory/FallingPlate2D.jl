# Held broadside plate in a unit stream for a range of Darcy numbers:
# drag should fall and the seepage velocity through the plate rise with Da.
using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
@printf("%-8s %8s %10s %10s %12s\n", "Da", "C_D", "u_seep/U", "screen/Δ", "pois iters")
for Da in (0.0, 1e-5, 1e-4, 1e-3, 1e-2)
    fp = falling_plate(; χ=4, ρ=2.7, Re=200, L=32, dom=(8,8), φ₀=0.0, free=(false,false,false), Da)
    sim = fp.sim; fp.Ufar[] = SA[1.0, 0.0]; fp.uexit .= 1.0
    while sim_time(sim) < 4.0
        sim.flow.Δt[end] = min(sim.flow.Δt[end], 0.5)
        prepare_flow!(fp); sim_step!(sim; remeasure=false)
    end
    F, _ = FP.body_loads(fp)
    I = CartesianIndex(round.(Int, fp.x₀ .+ 1.5)...)
    @printf("%-8.0e %8.3f %10.4f %10.2f %12d\n", Da, 2F[1]/fp.L, sim.flow.u[I, 1], screening_length(Da, fp.L), last(sim.pois.n))
end
