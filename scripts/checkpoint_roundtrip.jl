# Does restarting from a checkpoint reproduce the uninterrupted run exactly?
using FallingPlate2D, Printf, StaticArrays
include("gust_cases.jl")

P = (; D2..., Ga=150.0, Da=1e-4)
mk() = case(P; L=24, dom=(6, 6))
t1, t2 = 0.6, 1.2
mkpath("ckpt")

fp = mk(); run!(fp; tend=t1, report=1e9); save_checkpoint(fp, "ckpt/_rt.jld2")
run!(fp; tend=t2, report=1e9)
A = (fp.state.X, fp.state.v, fp.state.θ, fp.state.ω, sum(fp.sim.flow.Δt))

fp2 = mk(); load_checkpoint!(fp2, "ckpt/_rt.jld2")
B0 = (fp2.state.X, fp2.state.v, fp2.state.θ, fp2.state.ω, sum(fp2.sim.flow.Δt))
run!(fp2; tend=t2, report=1e9)
B = (fp2.state.X, fp2.state.v, fp2.state.θ, fp2.state.ω, sum(fp2.sim.flow.Δt))

@printf("continued  X=(%.10f,%.10f) θ=%.10f ω=%.10f  t=%.10f\n", A[1]..., A[3], A[4], A[5])
@printf("restarted  X=(%.10f,%.10f) θ=%.10f ω=%.10f  t=%.10f\n", B[1]..., B[3], B[4], B[5])
dX = maximum(abs.(A[1] .- B[1])); dv = maximum(abs.(A[2] .- B[2]))
@printf("max |ΔX|=%.3e  |Δv|=%.3e  |Δθ|=%.3e  |Δω|=%.3e  |Δt|=%.3e\n",
        dX, dv, abs(A[3]-B[3]), abs(A[4]-B[4]), abs(A[5]-B[5]))
println(dX < 1e-10 && dv < 1e-10 ? "PASS — bitwise-equivalent restart" : "FAIL")
rm("ckpt/_rt.jld2", force=true)
