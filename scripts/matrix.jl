# Da × gravity-modulation-amplitude matrix at Ga = 150.
#
#   matrix.jl lap <Da>                     lap to the terminal limit cycle, save checkpoint
#   matrix.jl run <Da> const               continue 10 c/u_g at constant gravity
#   matrix.jl run <Da> <A> down|up         continue 10 c/u_g with the gravity modulated
#                                          over one period T_G = 2: "down" falls to
#                                          (1-A)g first then rises to (1+A)g, "up" the
#                                          reverse. A is the semi-amplitude in units of g.
using FallingPlate2D, Printf
include("gust_cases.jl")

const TLAP, TADD, LGRID = 50.0, 10.0, 64
tag(Da) = replace(@sprintf("%.0e", Da), "e-0" => "e-")
params(Da) = (; D2..., Ga=150.0, Da)

mode = ARGS[1]; Da = parse(Float64, ARGS[2]); T = tag(Da)
mkpath("out"); mkpath("ckpt")

if mode == "lap"
    fp = case(params(Da); L=LGRID)
    header(fp, "lap Da=$T → t=$TLAP")
    st = run!(fp; tend=TLAP, out="out/mlap_$T.csv", every=0.05, report=5.0)
    st == :done || error("lap Da=$T aborted")
    save_checkpoint(fp, "ckpt/mlap_$T.jld2")
    println("saved ckpt/mlap_$T.jld2")

elseif mode == "run"
    g = 1/params(Da).m
    if ARGS[3] == "const"
        gust, name = NOGUST, "m_$(T)_const"
    else
        A = parse(Float64, ARGS[3]); ord = ARGS[4]
        ord in ("down", "up") || error("order must be down or up")
        # amplitude>0 lowers gravity first; <0 raises it first
        gust = gravity_modulation(; amplitude=(ord == "down" ? A : -A), T_G=2.0, t₀=TLAP, g)
        name = @sprintf("m_%s_A%.1f_%s", T, A, ord)
    end
    fp = case(params(Da); L=LGRID, gust)
    header(fp, name)
    load_checkpoint!(fp, "ckpt/mlap_$T.jld2")
    st = run!(fp; tend=TLAP+TADD, out="out/$name.csv", every=0.02, report=2.0,
              vmax=8.0, ωmax=20.0)          # relaxed: A≥1.5 reverses gravity
    println("$name: $st")
else
    error("mode must be lap or run")
end
