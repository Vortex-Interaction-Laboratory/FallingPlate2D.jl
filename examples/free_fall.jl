# Free fall of a plate section. Defaults reproduce the configuration of the
# direct numerical simulations of Andersen, Pesavento & Wang (2005): an
# ellipse with β = 1/8 (χ = 8) whose density ratio ρ* = 2.52 gives the
# dimensionless moment of inertia I* = 0.16 of their β = 1/14 aluminium
# plate (table 1), with a gravitational Reynolds number
# c·sqrt(2(ρ*−1) f g)/ν = 2090 so that Re on the mean descent speed is
# ca. 1147 if <vY> ≈ 0.55 U as in the experiment.
# (χ = 14, ρ* = 2.7 gives I* = 0.097 for an ellipse and an edge-on dive to
# 6 U; see BRIEF section 8.)
#
#   julia -t 4 --project=. examples/free_fall.jl [χ ρ Re L tend name]
using FallingPlate2D, Printf
χ    = length(ARGS) ≥ 1 ? parse(Float64, ARGS[1]) : 8.0
ρ    = length(ARGS) ≥ 2 ? parse(Float64, ARGS[2]) : 2.52
Re   = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 2090.0
L    = length(ARGS) ≥ 4 ? parse(Int, ARGS[4])     : 64
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
