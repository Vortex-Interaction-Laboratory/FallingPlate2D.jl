# Free-falling permeable plate at Ga = 150, Da = 1e-4.
#
#   lap1  release from rest, run to terminal velocity, save the state
#   end1  continue lap1 for 10 convective times at constant gravity
#   end2  continue lap1 for 10 convective times with the gravitational force
#         modulated over one period of two convective times: down to 0.5 g,
#         back to g, up to 1.5 g, back to g.
#
# Geometry and inertia are the established D2 case (the 2D equivalent of the
# 3D permeable disk: χ=10, m=0.80, e=1.16, I_O=1.15); only Ga and Da are set
# by the request.
#
#   julia -t N --project=. scripts/gravity_lap.jl {lap1|end1|end2} [L]
using FallingPlate2D, Printf
include("gust_cases.jl")

mode = length(ARGS) ≥ 1 ? ARGS[1] : "lap1"
L    = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 64
P    = (; D2..., Ga=150.0, Da=1e-4)
tlap, tadd = 20.0, 10.0          # to terminal velocity, then the comparison window
g    = 1/P.m                     # nondimensional gravity
mkpath("out"); mkpath("ckpt")

if mode == "lap1"
    fp = case(P; L)
    header(fp, "lap1  (release from rest → terminal velocity)")
    run!(fp; tend=tlap, out="out/lap1.csv", every=0.02, report=1.0)
    save_checkpoint(fp, "ckpt/lap1.jld2")
    println("saved ckpt/lap1.jld2 at t = ", round(tlap, digits=3))

elseif mode in ("end1", "end2")
    gust = mode == "end1" ? NOGUST :
           gravity_modulation(; amplitude=0.5, T_G=2.0, t₀=tlap, g)
    fp = case(P; L, gust)
    header(fp, mode == "end1" ? "end1  (constant gravity)" : "end2  (modulated gravity)")
    load_checkpoint!(fp, "ckpt/lap1.jld2")
    @printf("restarted at t=%.3f  Y/c=%.4f  vY=%.4f\n",
            sum(fp.sim.flow.Δt[1:end-1])*fp.U/fp.L, fp.state.X[2]/fp.L, fp.state.v[2]/fp.U)
    run!(fp; tend=tlap+tadd, out="out/$mode.csv", every=0.02, report=1.0)
    save_checkpoint(fp, "ckpt/$mode.jld2")
    println("saved ckpt/$mode.jld2")
else
    error("mode must be lap1, end1 or end2")
end
