using FallingPlate2D
for Da in (0.0, 1e-3)
    println("\n=== 3-DOF fall χ=4 ρ=2.7 Re=500, Da=$Da, 8 c/U ===")
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, Da)
    run!(fp; tend=8.0, out="out/phaseC_Da$(Da).csv", every=0.05, report=2.0)
    summarise("out/phaseC_Da$(Da).csv"; window=4.0, ν=fp.ν, L=fp.L)
end
