using FallingPlate2D, StaticArrays, Printf
println("=== 1-DOF (Y) fall, WaterLily exitBC=true, flow now +x ===")
# temporarily build with WaterLily's own exit to see whether it is stable when used the right way round
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, free=(false,true,false), exit=:waterlily)
run!(fp; tend=8.0, report=2.0)
println("\n=== 3-DOF fall χ=4 ρ=2.7 Re=500, own exit, 12 c/U ===")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1)
run!(fp; tend=12.0, out="out/checkT3.csv", every=0.05, report=1.0)
summarise("out/checkT3.csv"; window=6.0, ν=fp.ν, L=fp.L)
