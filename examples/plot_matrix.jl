# Da x semi-amplitude x ordering matrix at Ga = 150.
#   julia --project=examples examples/plot_matrix.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
DAS = ("1e-5", "1e-4", "1e-3"); AS = (0.5, 1.5, 2.5)
COL = Dict(0.5=>:navy, 1.5=>:firebrick, 2.5=>:darkorange)
ip(t, y, q) = (i = searchsortedfirst(t, q); i ≤ 1 ? y[1] : i > length(t) ? y[end] :
               y[i-1] + (y[i]-y[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
default(fontfamily="Computer Modern", linewidth=1.4, framestyle=:box, grid=false, dpi=200,
        legendfontsize=6)

panels = Any[]
resid = Dict{Tuple{String,Float64},Tuple{Float64,Float64}}()
for D in DAS
    c = FP.readrun("out/m_$(D)_const.csv")
    p = plot(xlabel="\$t\\,u_g/c\$", ylabel=(D == DAS[1] ? "\$\\Delta Y/c\$" : ""),
             title="\$Da=10^{$(D[end-1:end])}\$", titlefontsize=9, legend=:topright)
    vspan!(p, [50, 52]; color=:gray, alpha=0.15, label="")
    hline!(p, [0]; color=:gray, linewidth=0.8, label="")
    for A in AS, o in ("down", "up")
        r = FP.readrun(@sprintf("out/m_%s_A%.1f_%s.csv", D, A, o))
        ts = range(50, 60; length=400)
        d = [ip(r.t, r.Y, q) - ip(c.t, c.Y, q) for q in ts]
        plot!(p, ts, d; color=COL[A], linestyle=(o == "down" ? :solid : :dash),
              label=(o == "down" ? "\$A=$A\$" : ""))
        resid[(D, A)] = (o == "down" ? (d[end], get(resid, (D,A), (0.0,0.0))[2]) :
                                       (get(resid, (D,A), (0.0,0.0))[1], d[end]))
    end
    push!(panels, p)
end

# base states
pb = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$v_Y/u_g\$", title="base states (constant \$g\$)",
          titlefontsize=9, legend=:right)
for (D, col) in zip(DAS, (:seagreen, :navy, :firebrick))
    c = FP.readrun("out/m_$(D)_const.csv")
    plot!(pb, c.t, c.vY; color=col, label="\$10^{$(D[end-1:end])}\$")
end
push!(panels, pb)

plot(panels...; layout=(1, 4), size=(1800, 400), left_margin=8Plots.mm, bottom_margin=8Plots.mm)
savefig("out/matrix.png"); println("saved out/matrix.png")

@printf("\n%-8s %5s %12s %12s %12s\n", "Da", "A", "peak down", "residual dn", "residual up")
for D in DAS, A in AS
    c = FP.readrun("out/m_$(D)_const.csv")
    r = FP.readrun(@sprintf("out/m_%s_A%.1f_down.csv", D, A))
    ts = range(50, 60; length=400); d = [ip(r.t,r.Y,q) - ip(c.t,c.Y,q) for q in ts]
    pk = d[argmax(abs.(d))]
    @printf("%-8s %5.1f %12.4f %12.5f %12.5f\n", D, A, pk, resid[(D,A)][1], resid[(D,A)][2])
end
