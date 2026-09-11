using FallingPlate2D
println("=== 1-DOF (Y) coupled fall, own convective exit ===")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, free=(false,true,false))
run!(fp; tend=8.0, out="out/isolate.csv", every=0.05, report=1.0)
println("\n=== 3-DOF coupled fall (T3 repeat), own convective exit ===")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1)
run!(fp; tend=10.0, out="out/checkT3.csv", every=0.05, report=1.0)
summarise("out/checkT3.csv"; window=4.0, ν=fp.ν, L=fp.L)
