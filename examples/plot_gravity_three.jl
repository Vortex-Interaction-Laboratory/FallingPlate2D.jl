# Three permeable plates branched from the same state at t = 50 (Ga=150, Da=1e-4):
#   end1  constant gravity
#   end2  gravity down first:  g → 0.5g → g → 1.5g → g
#   end3  gravity up first:    g → 1.5g → g → 0.5g → g
#   julia --project=examples examples/plot_gravity_three.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D

_a = FP.readrun("out/lap1.csv"); _b = FP.readrun("out/lap1_ext.csv"); _s = _b.t .> _a.t[end]
lap = (; t=vcat(_a.t,_b.t[_s]), X=vcat(_a.X,_b.X[_s]), Y=vcat(_a.Y,_b.Y[_s]),
         phi=vcat(_a.phi,_b.phi[_s]), vY=vcat(_a.vY,_b.vY[_s]), omega=vcat(_a.omega,_b.omega[_s]))
c1, c2, c3 = FP.readrun("out/end1.csv"), FP.readrun("out/end2.csv"), FP.readrun("out/end3.csv")
m, T_G, t₀ = 0.80, 2.0, 50.0; g = 1/m
gd = gravity_modulation(; amplitude= 0.5, T_G, t₀, g)
gu = gravity_modulation(; amplitude=-0.5, T_G, t₀, g)
ratio(gm, t) = 1 + gust_acceleration(gm, t)/g

interp(t, y, q) = (i = searchsortedfirst(t, q); i ≤ 1 ? y[1] : i > length(t) ? y[end] :
                   y[i-1] + (y[i]-y[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
ts = range(t₀, min(c1.t[end], c2.t[end], c3.t[end]); length=600)
ΔY2 = [interp(c2.t,c2.Y,q) - interp(c1.t,c1.Y,q) for q in ts]
ΔY3 = [interp(c3.t,c3.Y,q) - interp(c1.t,c1.Y,q) for q in ts]

default(fontfamily="Computer Modern", linewidth=1.6, framestyle=:box, grid=false, dpi=200)
CC, C2, C3, SH = :navy, :firebrick, :darkorange, :gray65
band = (t₀, t₀ + T_G)

pa = plot(xlabel="\$X/c\$", ylabel="\$Y/c\$", legend=:bottomleft, aspect_ratio=:equal,
          legendfontsize=7)
plot!(pa, lap.X, lap.Y; color=SH, label="shared, \$t<50\$")
plot!(pa, c1.X, c1.Y; color=CC, label="constant")
plot!(pa, c2.X, c2.Y; color=C2, linestyle=:dash, label="down first")
plot!(pa, c3.X, c3.Y; color=C3, linestyle=:dashdot, label="up first")

pb = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$g_\\mathrm{eff}/g\$", legend=:topright, legendfontsize=7)
vspan!(pb, [band...]; color=:gray, alpha=0.15, label="")
plot!(pb, ts, fill(1.0, length(ts)); color=CC, label="constant")
plot!(pb, ts, ratio.(Ref(gd), ts); color=C2, linestyle=:dash, label="down first")
plot!(pb, ts, ratio.(Ref(gu), ts); color=C3, linestyle=:dashdot, label="up first")

pc = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$v_Y/u_g\$", legend=:bottomright, legendfontsize=7)
vspan!(pc, [band...]; color=:gray, alpha=0.15, label="")
plot!(pc, c1.t, c1.vY; color=CC, label="constant")
plot!(pc, c2.t, c2.vY; color=C2, linestyle=:dash, label="down first")
plot!(pc, c3.t, c3.vY; color=C3, linestyle=:dashdot, label="up first")

pd = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$\\Delta Y/c\$", legend=:topright, legendfontsize=7)
vspan!(pd, [band...]; color=:gray, alpha=0.15, label="")
plot!(pd, ts, ΔY2; color=C2, linestyle=:dash, label="down first")
plot!(pd, ts, ΔY3; color=C3, linestyle=:dashdot, label="up first")
plot!(pd, ts, ΔY2 .+ ΔY3; color=:black, linewidth=1.2, label="sum")
hline!(pd, [0]; color=:gray, linewidth=0.8, label="")

plot(pa, pb, pc, pd; layout=@layout([a{0.19w} b c d]), size=(1750, 430),
     left_margin=7Plots.mm, bottom_margin=8Plots.mm)
savefig("out/gravity_three.png"); println("saved out/gravity_three.png")

@printf("\nbranch point t=%.1f: Y/c=%.4f vY=%.4f φ=%+.4f ω=%+.4f\n",
        lap.t[end], lap.Y[end], lap.vY[end], lap.phi[end], lap.omega[end])
@printf("%-14s %10s %10s %10s\n", "", "Y(60)/c", "ΔY/c", "vY(60)")
for (n, r, d) in (("constant", c1, 0.0), ("down first", c2, ΔY2[end]), ("up first", c3, ΔY3[end]))
    @printf("%-14s %10.4f %+10.5f %10.5f\n", n, interp(r.t,r.Y,ts[end]), d, interp(r.t,r.vY,ts[end]))
end
@printf("\npeak ΔY: down first %+0.5f at t=%.2f;  up first %+0.5f at t=%.2f\n",
        maximum(abs, ΔY2)*sign(ΔY2[argmax(abs.(ΔY2))]), ts[argmax(abs.(ΔY2))],
        maximum(abs, ΔY3)*sign(ΔY3[argmax(abs.(ΔY3))]), ts[argmax(abs.(ΔY3))])
@printf("asymmetry at t=60: ΔY(down) + ΔY(up) = %+.5f c  (zero if the response were linear)\n",
        ΔY2[end] + ΔY3[end])
