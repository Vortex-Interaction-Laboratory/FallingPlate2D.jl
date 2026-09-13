# Approach to the terminal state: the lap run (t = 0 → 50) at Ga=150, Da=1e-4.
#   julia --project=examples examples/plot_lap.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
a = FP.readrun("out/lap1.csv"); b = FP.readrun("out/lap1_ext.csv")
sel = b.t .> 20.0
t  = vcat(a.t, b.t[sel]);  X = vcat(a.X, b.X[sel]);  Y = vcat(a.Y, b.Y[sel])
vY = vcat(a.vY, b.vY[sel]); φ = vcat(a.phi, b.phi[sel])
Tp = 5.8                                   # measured pendulum period

# cycle means and amplitudes
edges = collect(8.0:Tp:t[end]); cm = Float64[]; ca = Float64[]; ct = Float64[]
for lo in edges[1:end-1]
    s = (t .≥ lo) .& (t .< lo+Tp); any(s) || continue
    push!(ct, lo+Tp/2); push!(cm, sum(vY[s])/count(s))
    push!(ca, (maximum(φ[s]) - minimum(φ[s]))/2)
end

default(fontfamily="Computer Modern", linewidth=1.4, framestyle=:box, grid=false, dpi=200)
pa = plot(X, Y; color=:black, label=false, xlabel="\$X/c\$", ylabel="\$Y/c\$",
          aspect_ratio=:equal, title="trajectory", titlefontsize=10)
for k in 1:120:length(t)
    plot!(pa, [X[k]-0.5cos(φ[k]), X[k]+0.5cos(φ[k])], [Y[k]-0.5sin(φ[k]), Y[k]+0.5sin(φ[k])];
          color=:firebrick, linewidth=1, label=false)
end
pb = plot(t, vY; color=:navy, label="\$v_Y\$", xlabel="\$t\\,u_g/c\$", ylabel="\$v_Y/u_g\$", legend=:topright)
plot!(pb, ct, cm; color=:black, marker=:circle, markersize=3, linewidth=1.8, label="cycle mean")
hline!(pb, [cm[end]]; color=:gray, linestyle=:dash, linewidth=1, label=@sprintf("terminal %.3f", cm[end]))
pc = plot(t, φ; color=:seagreen, label="\$\\varphi\$", xlabel="\$t\\,u_g/c\$", ylabel="\$\\varphi\$ [rad]", legend=:topleft)
plot!(pc, ct, ca; color=:black, marker=:circle, markersize=3, linewidth=1.8, label="amplitude")
plot!(pc, ct, -ca; color=:black, linewidth=1.8, label=false)
plot(pa, pb, pc; layout=@layout([a{0.22w} b c]), size=(1500, 450),
     left_margin=7Plots.mm, bottom_margin=7Plots.mm)
savefig("out/lap1_approach.png"); println("saved out/lap1_approach.png")
@printf("terminal cycle mean vY = %.4f u_g;  φ amplitude = %.4f rad;  period ≈ %.1f c/u_g\n", cm[end], ca[end], Tp)
@printf("last two cycle means: %.5f, %.5f (change %+.2e)\n", cm[end-1], cm[end], cm[end]-cm[end-1])
