# Which parameter drives the rotational coupling instability seen at L=80, χ=10, Re=1000, ρ_s=2?
using WaterLily, StaticArrays, Printf
include(joinpath(@__DIR__, "phaseA_timing_defs.jl"))
function stab(L; χ, Re, ρ_s, tend, dom=(8,16))
    sim, x₀, a, b = make_sim(L; χ, Re, dom)
    Is = ρ_s*π*a^3*b/4; Ia = π/8*(a^2-b^2)^2
    θ, ω, α = sim.body.map.θ, 0.0, 0.0
    ωmax = 0.0; n = 0; blew = false
    while sim_time(sim) < tend
        Δt = sim.flow.Δt[end]
        M = -WaterLily.pressure_moment(x₀, sim.flow, sim.body)[1]
        α = (M + Ia*α)/(Is + Ia); ω += Δt*α; θ += Δt*ω
        sim.body = setmap(sim.body; θ, ω)
        sim_step!(sim; remeasure=true); n += 1
        ωmax = max(ωmax, abs(ω)*sim.L/sim.U)
        if abs(ω)*sim.L/sim.U > 20 || Δt < 1e-3; blew = true; break; end
    end
    @printf("L=%3d χ=%2d Re=%4d ρ_s=%4.1f Is/Ia=%5.2f  t_end=%5.2f steps=%5d  |ω*|max=%7.2f  θ_end=%7.2f  %s\n",
        L, χ, Re, ρ_s, Is/Ia, sim_time(sim), n, ωmax, θ, blew ? "UNSTABLE" : "stable"); flush(stdout)
end
stab(16; χ=4, Re=100, ρ_s=2.0, tend=0.3)                 # compile
println("--- L=80 variants (2 c/v each unless unstable) ---")
stab(80; χ=10, Re=1000, ρ_s=20.0, tend=2.0)   # heavy plate: is it a mass-ratio (added-mass mismatch) instability?
stab(80; χ=10, Re=100,  ρ_s=2.0,  tend=2.0)   # low Re: is it a Reynolds-number effect?
stab(80; χ=4,  Re=1000, ρ_s=2.0,  tend=2.0)   # thick body: is it thinness / cells across thickness?
