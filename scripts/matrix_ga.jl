# Galilei-number sweep at the largest amplitude, Da = 1e-5, rotation locked.
# Same branch time and window as the Ga = 150 campaign: lap to t = 30, branch
# t = 30 -> 40, so the results are directly comparable rather than converged.
#
#   matrix_ga.jl lap <Ga>
#   matrix_ga.jl run <Ga> const
#   matrix_ga.jl run <Ga> v down|up        vertical   fluctuation of g_y
#   matrix_ga.jl run <Ga> h pos|neg        horizontal fluctuation of g_x
using FallingPlate2D, Printf, StaticArrays
include("gust_cases.jl")
const TLAP, TADD, LGRID, DA, A = 30.0, 10.0, 64, 1e-5, 2.5
build(Ga; gust=NOGUST) = falling_plate(; χ=D2.χ, m=D2.m, e=D2.e, I_O=D2.I_O, Da=DA, Ga,
                                        L=LGRID, dom=(16,16), xbody=8, φ₀=0.0, shape=:plate,
                                        gust, Δt_max=0.5, free=(true,true,false))
mode = ARGS[1]; Ga = parse(Float64, ARGS[2]); G = @sprintf("Ga%d", round(Int, Ga))
mkpath("out"); mkpath("ckpt")

if mode == "lap"
    fp = build(Ga); header(fp, "lap $G Da=1e-5 → t=$TLAP")
    run!(fp; tend=TLAP, out="out/galap_$G.csv", every=0.05, report=5.0) == :done || error("aborted")
    save_checkpoint(fp, "ckpt/galap_$G.jld2"); println("saved ckpt/galap_$G.jld2")
else
    g = 1/D2.m
    if ARGS[3] == "const"
        gust, name = NOGUST, "g_$(G)_const"
    else
        dir, ord = ARGS[3], ARGS[4]
        e, amp = dir == "v" ? (SA[0.0,1.0], ord == "down" ?  A : -A) :
                              (SA[1.0,0.0], ord == "pos"  ?  A : -A)
        gust = gravity_modulation(; amplitude=amp, T_G=2.0, t₀=TLAP, g, dir=e)
        name = "g_$(G)_$(dir)_$(ord)"
    end
    fp = build(Ga; gust); header(fp, name)
    load_checkpoint!(fp, "ckpt/galap_$G.jld2")
    st = run!(fp; tend=TLAP+TADD, out="out/$name.csv", every=0.02, report=2.0, vmax=8.0, ωmax=20.0)
    save_checkpoint(fp, "ckpt/$name.jld2"); println("$name: $st")
end
