# Toy Darcy sweep in free fall: χ=4, ρ*=2.7, Re_g=500, L=32, 20 c/U each.
using FallingPlate2D, Printf
rows = []
for Da in (0.0, 1e-4, 1e-3, 1e-2)
    println("\n=== Da = $Da ===")
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, Da)
    st = run!(fp; tend=20.0, out=@sprintf("out/sweep_Da%.0e.csv", Da), every=0.05, report=5.0)
    r = summarise(@sprintf("out/sweep_Da%.0e.csv", Da); window=10.0, ν=fp.ν, L=fp.L)
    push!(rows, (Da, st, r))
end
println("\n=== summary: χ=4 ρ*=2.7 Re_g=500 L=32, means over 10–20 c/U ===")
@printf("%-8s %-8s %8s %8s %8s %8s\n", "Da", "status", "<vY>/U", "<vX>/U", "<ω>c/U", "Re_vY")
for (Da, st, r) in rows
    @printf("%-8.0e %-8s %8.3f %8.3f %8.3f %8.0f\n", Da, st, r.vY, r.vX, r.ω, r.Re)
end
