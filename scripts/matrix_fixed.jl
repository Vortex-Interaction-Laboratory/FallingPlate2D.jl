# As matrix.jl, but with the angular position constrained: the plate is held
# horizontal (φ = 0) and cannot rotate, free = (X, Y, θ) = (true, true, false).
#
# This removes the pendulum/flutter mode, so any lasting displacement after the
# modulation must come from the translational mechanism — wake and added-mass
# memory — rather than from the phase of a rotational oscillation.
#
#   matrix_fixed.jl lap <Da>
#   matrix_fixed.jl run <Da> const | <A> down|up
using FallingPlate2D, Printf
include("gust_cases.jl")

const TLAP, TADD, LGRID = 30.0, 10.0, 64
tag(Da) = replace(@sprintf("%.0e", Da), "e-0" => "e-")
build(Da; gust=NOGUST) = falling_plate(; χ=D2.χ, m=D2.m, e=D2.e, I_O=D2.I_O, Da, Ga=150.0,
                                        L=LGRID, dom=(16, 16), xbody=8, φ₀=0.0, shape=:plate,
                                        gust, Δt_max=0.5, free=(true, true, false))
mode = ARGS[1]; Da = parse(Float64, ARGS[2]); T = tag(Da)
mkpath("out"); mkpath("ckpt")

if mode == "lap"
    fp = build(Da)
    header(fp, "FIXED-θ lap Da=$T → t=$TLAP")
    run!(fp; tend=TLAP, out="out/flap_$T.csv", every=0.05, report=5.0) == :done || error("lap aborted")
    save_checkpoint(fp, "ckpt/flap_$T.jld2")
    println("saved ckpt/flap_$T.jld2")
else
    g = 1/D2.m
    if ARGS[3] == "const"
        gust, name = NOGUST, "f_$(T)_const"
    else
        A = parse(Float64, ARGS[3]); ord = ARGS[4]
        gust = gravity_modulation(; amplitude=(ord == "down" ? A : -A), T_G=2.0, t₀=TLAP, g)
        name = @sprintf("f_%s_A%.1f_%s", T, A, ord)
    end
    fp = build(Da; gust)
    header(fp, name)
    load_checkpoint!(fp, "ckpt/flap_$T.jld2")
    st = run!(fp; tend=TLAP+TADD, out="out/$name.csv", every=0.02, report=2.0, vmax=8.0, ωmax=20.0)
    println("$name: $st")
end
