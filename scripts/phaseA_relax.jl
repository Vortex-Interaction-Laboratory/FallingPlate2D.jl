# Does over-estimating the added inertia in the implicit term (Ia_eff = k·Ia) stabilise
# the χ=10, Re=1000, ρ_s=2 case?  At convergence the k·Ia·α terms cancel, so k only
# affects stability and the transient lag, not the solution.
using WaterLily, StaticArrays, Printf
include(joinpath(@__DIR__, "phaseA_timing_defs.jl"))
function relax(L; χ, Re, ρ_s, k, tend, dom=(8,16))
    sim, x₀, a, b = make_sim(L; χ, Re, dom)
    Is = ρ_s*π*a^3*b/4; Ia = π/8*(a^2-b^2)^2; Ie = k*Ia
    θ, ω, α = sim.body.map.θ, 0.0, 0.0
    ωmax = 0.0; n = 0; blew = false
    while sim_time(sim) < tend
        Δt = sim.flow.Δt[end]
        M = -WaterLily.pressure_moment(x₀, sim.flow, sim.body)[1]
        α = (M + Ie*α)/(Is + Ie); ω += Δt*α; θ += Δt*ω
        sim.body = setmap(sim.body; θ, ω)
        sim_step!(sim; remeasure=true); n += 1
        ωmax = max(ωmax, abs(ω)*sim.L/sim.U)
        if abs(ω)*sim.L/sim.U > 20 || Δt < 1e-3; blew = true; break; end
    end
    @printf("L=%3d χ=%2d Re=%4d ρ_s=%4.1f Is/Ia=%5.2f k=%3.1f  t_end=%5.2f steps=%5d  |ω*|max=%7.2f  θ_end=%7.2f  %s\n",
        L, χ, Re, ρ_s, Is/Ia, k, sim_time(sim), n, ωmax, θ, blew ? "UNSTABLE" : "stable"); flush(stdout)
end
relax(16; χ=4, Re=100, ρ_s=2.0, k=1.0, tend=0.3)          # compile
relax(80; χ=10, Re=1000, ρ_s=2.0, k=2.0, tend=2.0)
relax(80; χ=10, Re=1000, ρ_s=2.0, k=4.0, tend=2.0)
