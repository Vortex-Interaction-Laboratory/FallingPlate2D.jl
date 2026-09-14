# Horizontal gravity perturbation against the vertical one, rotation locked, Ga = 150.
#   julia --project=examples examples/plot_horiz.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
DAS = ("1e-5", "1e-4"); AS = (0.5, 1.5, 2.5)
COL = Dict(0.5=>:navy, 1.5=>:firebrick, 2.5=>:darkorange)
ip(t,y,q) = (i = searchsortedfirst(t,q); i ≤ 1 ? y[1] : i > length(t) ? y[end] :
             y[i-1] + (y[i]-y[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
default(fontfamily="Computer Modern", linewidth=1.4, framestyle=:box, grid=false, dpi=200,
        legendfontsize=6)
P = Any[]
for D in DAS
    c = FP.readrun("out/f_$(D)_const.csv")
    p = plot(xlabel="\$t\\,u_g/c\$", ylabel=(D==DAS[1] ? "\$\\Delta X/c\$" : ""),
             title="horizontal, \$Da=10^{$(D[end-1:end])}\$", titlefontsize=9, legend=:topleft)
    vspan!(p,[30,32]; color=:gray, alpha=0.15, label=""); hline!(p,[0]; color=:gray, linewidth=0.8, label="")
    for A in AS, o in ("pos","neg")
        r = FP.readrun(@sprintf("out/h_%s_A%.1f_%s.csv", D, A, o))
        ts = range(30,40; length=400)
        plot!(p, ts, [ip(r.t,r.X,q)-ip(c.t,c.X,q) for q in ts]; color=COL[A],
              linestyle=(o=="pos" ? :solid : :dash), label=(o=="pos" ? "\$A=$A\$" : ""))
    end
    push!(P, p)
end
# vertical, same axes scale, for direct amplitude comparison
for D in DAS
    c = FP.readrun("out/f_$(D)_const.csv")
    p = plot(xlabel="\$t\\,u_g/c\$", ylabel=(D==DAS[1] ? "\$\\Delta Y/c\$" : ""),
             title="vertical, \$Da=10^{$(D[end-1:end])}\$", titlefontsize=9, legend=:topleft)
    vspan!(p,[30,32]; color=:gray, alpha=0.15, label=""); hline!(p,[0]; color=:gray, linewidth=0.8, label="")
    for A in AS, o in ("down","up")
        r = FP.readrun(@sprintf("out/f_%s_A%.1f_%s.csv", D, A, o))
        ts = range(30,40; length=400)
        plot!(p, ts, [ip(r.t,r.Y,q)-ip(c.t,c.Y,q) for q in ts]; color=COL[A],
              linestyle=(o=="down" ? :solid : :dash), label=(o=="down" ? "\$A=$A\$" : ""))
    end
    ylims!(p, ylims(P[1]))                       # same scale as the horizontal panels
    push!(P, p)
end
plot(P...; layout=(1,4), size=(1800,400), left_margin=8Plots.mm, bottom_margin=8Plots.mm)
savefig("out/horiz.png"); println("saved out/horiz.png")
