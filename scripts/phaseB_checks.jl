# Phase B sanity checks (cheap, L = 32).
using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D

println("--- T1: drag sign on a held plate, broadside to a unit stream ---")
fp = falling_plate(; χ=4, ρ=2.7, Re=200, L=32, dom=(8,8), φ₀=0.0)
fp.Ufar[] = SA[1.0, 0.0]                        # stream along +x (grid)
while sim_time(fp.sim) < 3.0
    sim_step!(fp.sim; remeasure=false)
end
F, M = FP.body_loads(fp)
@printf("F_grid = (%.3f, %.3f)  →  Fx/(U²c) = %.3f  (must be > 0: drag along the stream; C_D = %.2f)\n",
        F..., F[1]/fp.L, 2F[1]/fp.L)

println("\n--- T2: numerical added mass from the first steps of a straight fall ---")
fp = falling_plate(; χ=4, ρ=10.0, Re=200, L=32, dom=(8,8), φ₀=0.0)
d = fp.dyn
a_lamb = (d.m - d.A)*fp.g/(d.m + d.m22)
@printf("Lamb: m=%.1f A=%.1f m22=%.1f  → a_Y expected = %.5f\n", d.m, d.A, d.m22, -a_lamb)
for n in 1:12
    F, M, Δt = FP.step!(fp)
    aY = fp.state.a[2]
    m22num = (d.m - d.A)*fp.g/abs(aY) - d.m
    @printf("n=%2d  Δt=%.3f  a_Y=%.5f  vY=%.4f  inferred m22_num/m22 = %.3f\n", n, Δt, aY, fp.state.v[2], m22num/d.m22)
end

println("\n--- T3: short coupled free fall, χ=4, ρ=2.7, Re=500, 6 c/U ---")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1)
Ia, Ir = andersen_Istar(fp.section, 2.7)
@printf("I* (actual ellipse) = %.3f, I* (rectangular formula) = %.3f\n", Ia, Ir)
status = run!(fp; tend=6.0, out="out/checkT3.csv", every=0.05, report=1.0)
summarise("out/checkT3.csv"; window=3.0, ν=fp.ν, L=fp.L)
