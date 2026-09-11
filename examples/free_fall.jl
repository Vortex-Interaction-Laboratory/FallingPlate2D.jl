# Free fall of a plate section; defaults reproduce the aluminium plate of Andersen,
# Pesavento & Wang (2005), table 1, β = 1/14: χ = 14, ρ* = 2.7, gravitational
# Reynolds number c·sqrt(2(ρ*−1) f g)/ν = 2090 (their Re on the mean descent
# speed is 1147, i.e. <vY> ≈ 0.55 U).
#
#   julia -t 4 --project=. examples/free_fall.jl [χ ρ Re L tend name]
using FallingPlate2D, Printf
χ    = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 14.0
ρ    = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 2.7
Re   = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 2090.0
L    = length(ARGS) ≥ 4 ? parse(Int, ARGS[4])     : 80
tend = length(ARGS) ≥ 5 ? parse(Float64, ARGS[5]) : 60.0
name = length(ARGS) ≥ 6 ? ARGS[6] : @sprintf("chi%g_rho%g_Re%g_L%d", χ, ρ, Re, L)

fp = falling_plate(; χ, ρ, Re, L, dom=(8, 8), xbody=3, φ₀=0.1)
Ia, Ir = andersen_Istar(fp.section, ρ)
@printf("FallingPlate2D free fall: χ=%g ρ*=%g Re_g=%g  L=%d cells (%.1f across the thickness)  grid %s\n",
        χ, ρ, Re, L, L/χ, size(fp.sim.flow.p))
@printf("I* = %.3f (actual ellipse), %.3f (Andersen's rectangular formula);  g = %.4g, ν = %.4g (grid units)\n", Ia, Ir, fp.g, fp.ν)
mkpath("out")
out = "out/$name.csv"
status = run!(fp; tend, out, every=0.05, report=1.0)
println("status = $status, output = $out")
summarise(out; window=min(30.0, tend/2), ν=fp.ν, L=fp.L)
