# Trajectories of the Darcy sweep at the Andersen inertia, one panel per Da.
#   julia --project=examples examples/plot_sweep.jl
using Plots, FallingPlate2D
FP = FallingPlate2D
default(fontfamily="Computer Modern", linewidth=1.0, framestyle=:box, grid=false, dpi=150, legend=false)
tags = ["0", "1e-6", "1e-5", "3e-5", "1e-4"]
ps = []
for tag in tags
    r = FP.readrun("out/sweep_Re1900_Da$tag.csv")
    p = plot(r.X, r.Y; color=:black, xlabel="X/c", ylabel=tag == "0" ? "Y/c" : "", aspect_ratio=:equal,
             title="Da = $tag", titlefontsize=10, xlims=(-25, 10), ylims=(-45, 2))
    for k in 1:10:length(r.t)
        x, y, φ = r.X[k], r.Y[k], r.phi[k]
        plot!(p, [x - 0.5cos(φ), x + 0.5cos(φ)], [y - 0.5sin(φ), y + 0.5sin(φ)]; color=:red)
    end
    push!(ps, p)
end
plot(ps...; layout=(1, 5), size=(1800, 620), bottom_margin=5Plots.mm, left_margin=3Plots.mm)
savefig("out/sweep_Re1900.png"); println("saved out/sweep_Re1900.png")
