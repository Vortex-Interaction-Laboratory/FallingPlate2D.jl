using FallingPlate2D, Printf
cases = [("V1: Y only, periodic y",            (; free=(false,true,false))),
         ("V2: Y only, Dirichlet y (perdir=())",(; free=(false,true,false), perdir=())),
         ("V4: Y only, periodic y, ρ=10",       (; free=(false,true,false), ρ=10.0))]
for (name, kw) in cases
    println("\n=== $name ===")
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, kw...)
    run!(fp; tend=8.0, out="out/isolate.csv", every=0.05, report=1.0)
end
