using FallingPlate2D
println("=== 1-DOF (Y) coupled fall, exitBC=false (Dirichlet exit) ===")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, free=(false,true,false), exitBC=false)
run!(fp; tend=8.0, out="out/isolate.csv", every=0.05, report=1.0)
