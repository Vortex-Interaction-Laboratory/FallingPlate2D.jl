# Phase A timing driver; definitions in phaseA_timing_defs.jl
using WaterLily, StaticArrays, Printf
include(joinpath(@__DIR__, "phaseA_timing_defs.jl"))

function bench(L; χ, Re, tend, dom=(8, 16))
    sim, x₀, a, b = make_sim(L; χ, Re, dom)
    t = @elapsed r = run!(sim, x₀, a, b; tend)
    N = prod(size(sim.flow.p))
    @printf("L=%3d χ=%2d grid=%4d×%4d Re=%4d  tend=%5.1f  steps=%5d  wall=%7.1f s  %6.1f ms/step  %5.1f ns/cell/step  θ=%6.3f ω=%6.3f  F=(%6.3f,%6.3f)  → 100 c/v ≈ %5.1f min\n",
        L, χ, size(sim.flow.p)..., Re, tend, r.steps, t, 1e3*t/r.steps, 1e9*t/r.steps/N, r.θ, r.ω, r.F..., 100/tend*t/60)
    flush(stdout)
end

println("Julia ", VERSION, ", WaterLily ", pkgversion(WaterLily), ", threads = ", Threads.nthreads(), ", T = Float64"); flush(stdout)
# warm-up (compilation) on a tiny case
bench(16; χ=4, Re=100, tend=0.5)
println("--- timings ---")
bench(32; χ=4,  Re=200,  tend=5.0)    # toy   mode: Δ = c/32,  256×512
bench(80; χ=10, Re=1000, tend=1.0)    # check mode: Δ = c/80,  640×1280 (per-step cost only)
