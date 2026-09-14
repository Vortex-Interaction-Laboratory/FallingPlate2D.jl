# Rotation locked (φ = 0) against the free-rotating case, Ga = 150.
#   julia --project=examples examples/plot_fixed.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
DAS = ("1e-5", "1e-4"); AS = (0.5, 1.5, 2.5)
COL = Dict(0.5=>:navy, 1.5=>:firebrick, 2.5=>:darkorange)
ip(t, y, q) = (i = searchsortedfirst(t, q); i ≤ 1 ? y[1] : i > length(t) ? y[end] :
               y[i-1] + (y[i]-y[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
default(fontfamily="Computer Modern", linewidth=1.4, framestyle=:box, grid=false, dpi=200,
        legendfontsize=6)
panels = Any[]; res = Dict()
for D in DAS
    c = FP.readrun("out/f_$(D)_const.csv")
    p = plot(xlabel="\$t\\,u_g/c\$", ylabel=(D == DAS[1] ? "\$\\Delta Y/c\$" : ""),
             title="rotation locked, \$Da=10^{$(D[end-1:end])}\$", titlefontsize=9, legend=:topleft)
    vspan!(p, [30, 32]; color=:gray, alpha=0.15, label=""); hline!(p, [0]; color=:gray, linewidth=0.8, label="")
    for A in AS, o in ("down", "up")
        r = FP.readrun(@sprintf("out/f_%s_A%.1f_%s.csv", D, A, o))
        ts = range(30, 40; length=400)
        d = [ip(r.t, r.Y, q) - ip(c.t, c.Y, q) for q in ts]
        plot!(p, ts, d; color=COL[A], linestyle=(o == "down" ? :solid : :dash),
              label=(o == "down" ? "\$A=$A\$" : ""))
        res[(D, A, o)] = d[end]
    end
    push!(panels, p)
end
# residual vs A: locked (both orderings positive) against free (antisymmetric)
freeres = Dict(("1e-5",0.5,"down")=>-0.0743, ("1e-5",1.5,"down")=>0.5773, ("1e-5",2.5,"down")=>0.6171,
               ("1e-5",0.5,"up")=>0.0656,   ("1e-5",1.5,"up")=>0.5664,  ("1e-5",2.5,"up")=>0.9312,
               ("1e-4",0.5,"down")=>0.0414, ("1e-4",1.5,"down")=>0.0972,("1e-4",2.5,"down")=>0.1198,
               ("1e-4",0.5,"up")=>-0.0439,  ("1e-4",1.5,"up")=>-0.0737, ("1e-4",2.5,"up")=>-0.1428)
pr = plot(xlabel="\$A\$", ylabel="residual \$\\Delta Y/c\$", title="lasting displacement",
          titlefontsize=9, legend=:topleft)
hline!(pr, [0]; color=:gray, linewidth=0.8, label="")
for (D, m) in zip(DAS, (:circle, :square))
    plot!(pr, collect(AS), [res[(D,A,"down")] for A in AS]; color=:black, marker=m, label="locked $D, down")
    plot!(pr, collect(AS), [res[(D,A,"up")] for A in AS]; color=:black, marker=m, linestyle=:dash, label="locked $D, up")
    plot!(pr, collect(AS), [freeres[(D,A,"down")] for A in AS]; color=:firebrick, marker=m, label="free $D, down")
    plot!(pr, collect(AS), [freeres[(D,A,"up")] for A in AS]; color=:firebrick, marker=m, linestyle=:dash, label="free $D, up")
end
push!(panels, pr)
pb = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$v_Y/u_g\$", title="base states", titlefontsize=9, legend=:right)
for (D, col) in zip(DAS, (:seagreen, :navy))
    c = FP.readrun("out/f_$(D)_const.csv"); plot!(pb, c.t, c.vY; color=col, label="locked \$10^{$(D[end-1:end])}\$")
end
push!(panels, pb)
plot(panels...; layout=(1,4), size=(1800,400), left_margin=8Plots.mm, bottom_margin=8Plots.mm)
savefig("out/fixed.png"); println("saved out/fixed.png")
