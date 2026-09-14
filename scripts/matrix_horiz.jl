# As matrix_fixed.jl (rotation locked, φ = 0), but the gravity oscillation is
# HORIZONTAL: the vertical component of gravity stays constant at g and a
# transient sideways force of peak magnitude A·g is added over one period.
#
# The base state is independent of the perturbation, so the laps and the
# constant-gravity references of the fixed-angle campaign are reused.
#
#   matrix_horiz.jl <Da> <A> pos|neg      pos = pushed along +X first
using FallingPlate2D, Printf, StaticArrays
include("gust_cases.jl")
const TLAP, TADD, LGRID = 30.0, 10.0, 64
tag(Da) = replace(@sprintf("%.0e", Da), "e-0" => "e-")

Da = parse(Float64, ARGS[1]); A = parse(Float64, ARGS[2]); ord = ARGS[3]
ord in ("pos", "neg") || error("order must be pos or neg")
T = tag(Da); g = 1/D2.m
gust = gravity_modulation(; amplitude=(ord == "pos" ? A : -A), T_G=2.0, t₀=TLAP, g,
                            dir=SA[1.0, 0.0])
fp = falling_plate(; χ=D2.χ, m=D2.m, e=D2.e, I_O=D2.I_O, Da, Ga=150.0, L=LGRID,
                     dom=(16, 16), xbody=8, φ₀=0.0, shape=:plate, gust, Δt_max=0.5,
                     free=(true, true, false))
name = @sprintf("h_%s_A%.1f_%s", T, A, ord)
header(fp, name)
load_checkpoint!(fp, "ckpt/flap_$T.jld2")         # reuse the fixed-angle terminal state
st = run!(fp; tend=TLAP+TADD, out="out/$name.csv", every=0.02, report=2.0, vmax=8.0, ωmax=20.0)
save_checkpoint(fp, "ckpt/$name.jld2")
println("$name: $st")
