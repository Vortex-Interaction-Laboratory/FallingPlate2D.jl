# Per-step diagnosis of the check-mode grid (L=80, χ=10, 640×1280).
using WaterLily, StaticArrays, Printf
include(joinpath(@__DIR__, "phaseA_timing_defs.jl"))
function diag(L; χ, Re, nsteps, dom=(8,16), every=1)
    sim, x₀, a, b = make_sim(L; χ, Re, dom)
    T = Float64
    Is = 2π*a^3*b/4; Ia = π/8*(a^2-b^2)^2
    θ, ω, α = sim.body.map.θ, 0.0, 0.0
    @printf("L=%d χ=%d grid=%s  Re=%d\n", L, χ, size(sim.flow.p), Re)
    for n in 1:nsteps
        t0 = time()
        tm = @elapsed M = -WaterLily.pressure_moment(x₀, sim.flow, sim.body)[1]
        α = (M + Ia*α)/(Is + Ia); ω += sim.flow.Δt[end]*α; θ += sim.flow.Δt[end]*ω
        sim.body = setmap(sim.body; θ=T(θ), ω=T(ω))
        tmeas = @elapsed measure!(sim)
        tstep = @elapsed sim_step!(sim; remeasure=false)
        umax = maximum(abs, sim.flow.u)
        (n % every == 0 || n ≤ 3) && @printf("n=%4d  t=%6.3f  Δt=%.4f  pois=%3d  |u|max=%8.3f  M=%9.3f  θ=%8.3f  ω*=%8.4f  wall/step %.3f s\n",
            n, sim_time(sim), sim.flow.Δt[end], last(sim.pois.n), umax, M, θ, ω*sim.L/sim.U, time()-t0); flush(stdout)
    end
end
diag(16; χ=4, Re=100, nsteps=3)          # compile
println("--- L=80 ---"); diag(80; χ=10, Re=1000, nsteps=400, every=20)
