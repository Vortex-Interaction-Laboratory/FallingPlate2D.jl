using FallingPlate2D, Printf
for (name, kw) in [("A: Δt_max=0.2", (; Δt_max=0.2)),
                   ("B: rotation locked", (; free=(true,true,false))),
                   ("C: translation locked, stream U=0 (rotation only, at rest)", (; free=(false,false,true)))]
    println("\n=== variant $name ===")
    fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, kw...)
    run!(fp; tend=6.0, out="out/isolate.csv", every=0.05, report=1.0)
end
