#   julia --project=examples examples/plot_gust.jl [L]
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
L = length(ARGS) ≥ 1 ? ARGS[1] : "64"
load(n) = FP.readrun("out/$n.csv")
ref, up, dn = load("D2_nogust_L$L"), load("D2_up_L$L"), load("D2_down_L$L")
function at(r, c, t)
    i = searchsortedfirst(r.t, t); i ≤ 1 && return getproperty(r, c)[1]
    i > length(r.t) && return getproperty(r, c)[end]
    y = getproperty(r, c); y[i-1] + (y[i]-y[i-1])*(t-r.t[i-1])/(r.t[i]-r.t[i-1])
end
ts = 10:0.05:minimum(r.t[end] for r in (ref, up, dn))
zu = [at(up, :Zearth, t) - at(ref, :Zearth, t) for t in ts]
zd = [at(dn, :Zearth, t) - at(ref, :Zearth, t) for t in ts]
default(fontfamily="Computer Modern", linewidth=1.6, framestyle=:box, grid=false, dpi=160)
band = (15, 17)
p1 = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$Z_G\$", legend=:topleft)
vspan!(p1, [band...]; color=:gray, alpha=0.18, label="")
plot!(p1, ts, zu; color=:firebrick, label="up-gust, \$GR=+1\$")
plot!(p1, ts, zd; color=:navy, label="down-gust, \$GR=-1\$")
plot!(p1, ts, zu .+ zd; color=:black, linestyle=:dash, label="sum")
hline!(p1, [0]; color=:gray, linewidth=0.8, label="")
p2 = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$u_Z/u_g\$", legend=:bottomright)
vspan!(p2, [band...]; color=:gray, alpha=0.18, label="")
for (r, c, l) in ((ref, :gray40, "quiescent"), (up, :firebrick, "up-gust"), (dn, :navy, "down-gust"))
    sel = r.t .≥ 10
    plot!(p2, r.t[sel], r.vY[sel] .+ r.uG[sel]; color=c, label=l)   # earth-frame vertical velocity
end
plot(p1, p2; layout=(1, 2), size=(1150, 430), left_margin=6Plots.mm, bottom_margin=6Plots.mm)
savefig("out/gust_L$L.png"); println("saved out/gust_L$L.png")
