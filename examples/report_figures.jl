# Figures for the narrative report. Vector PDF, VOILAb conventions.
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D
AS = (0.5, 1.5, 2.5); COL = Dict(0.5=>:black, 1.5=>:firebrick, 2.5=>:royalblue)
LS  = Dict(0.5=>:solid, 1.5=>:dash, 2.5=>:dashdot)
ip(t,y,q) = (i = searchsortedfirst(t,q); i ≤ 1 ? y[1] : i > length(t) ? y[end] :
             y[i-1] + (y[i]-y[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
default(fontfamily="Computer Modern", linewidth=1.6, framestyle=:box, grid=false,
        background_color=:white, legendfontsize=8, guidefontsize=11, tickfontsize=9)

# ---- figure 1: free to rotate, vertical perturbation, Da = 1e-4 -------------
c = FP.readrun("out/m_1e-4_const.csv"); ts = range(50, 60; length=400)
p1 = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$\\Delta Y/c\$", legend=:topright,
          xticks=50:2:60, size=(560,400))
vspan!(p1, [50,52]; color=:gray, alpha=0.12, label="")
hline!(p1, [0]; color=:gray, linewidth=0.7, label="")
for A in AS, (o,lab) in (("down","\$g_y\$ reduced first"), ("up","\$g_y\$ raised first"))
    r = FP.readrun(@sprintf("out/m_1e-4_A%.1f_%s.csv", A, o))
    plot!(p1, ts, [ip(r.t,r.Y,q)-ip(c.t,c.Y,q) for q in ts]; color=COL[A], linestyle=LS[A],
          linealpha=(o=="down" ? 1.0 : 0.45),
          label=(A==0.5 ? lab : (o=="down" ? "\$A=$A\$" : "")))
end
savefig(p1, "report/fig1.pdf")

# ---- figure 2: angular position fixed, vertical perturbation ---------------
function panel(pre, D, t0, ylab, os, labs)
    c = FP.readrun("out/f_$(D)_const.csv"); ts = range(t0, t0+10; length=400)
    p = plot(xlabel="\$t\\,u_g/c\$", ylabel=ylab, legend=:topleft, xticks=t0:2:t0+10)
    vspan!(p, [t0,t0+2]; color=:gray, alpha=0.12, label="")
    hline!(p, [0]; color=:gray, linewidth=0.7, label="")
    for A in AS, (o,lab) in zip(os, labs)
        r = FP.readrun(@sprintf("out/%s_%s_A%.1f_%s.csv", pre, D, A, o))
        plot!(p, ts, [ip(r.t,r.Y,q)-ip(c.t,c.Y,q) for q in ts]; color=COL[A], linestyle=LS[A],
              linealpha=(o==os[1] ? 1.0 : 0.45),
              label=(A==0.5 ? lab : (o==os[1] ? "\$A=$A\$" : "")))
    end
    p
end
va = panel("f","1e-5",30.0,"\$\\Delta Y/c\$",("down","up"),("\$g_y\$ reduced first","\$g_y\$ raised first"))
vb = panel("f","1e-4",30.0,"",("down","up"),("\$g_y\$ reduced first","\$g_y\$ raised first"))
ylims!(va,(-0.55,0.70)); ylims!(vb,(-0.55,0.70))
savefig(plot(va,vb; layout=(1,2), size=(1000,400), left_margin=6Plots.mm, bottom_margin=6Plots.mm),
        "report/fig2.pdf")

ha = panel("h","1e-5",30.0,"\$\\Delta Y/c\$",("pos","neg"),("\$g_x>0\$ first","\$g_x<0\$ first"))
hb = panel("h","1e-4",30.0,"",("pos","neg"),("\$g_x>0\$ first","\$g_x<0\$ first"))
ylims!(ha,(-0.3,6.2)); ylims!(hb,(-0.3,6.2))
savefig(plot(ha,hb; layout=(1,2), size=(1000,400), left_margin=6Plots.mm, bottom_margin=6Plots.mm),
        "report/fig3.pdf")

# ---- figure 4: magnitude comparison ---------------------------------------
p4 = plot(xlabel="\$A\$", ylabel="\$\\Delta Y/c\$", yscale=:log10, legend=:bottomright,
          xticks=(collect(AS), ["0.5","1.5","2.5"]), yticks=[1e-3,1e-2,1e-1,1e0,1e1], size=(560,420))
for (D,mk) in (("1e-5",:circle), ("1e-4",:square))
    cf = FP.readrun("out/f_$(D)_const.csv")
    hv = Float64[]; vv = Float64[]
    for A in AS
        h = [ (r = FP.readrun(@sprintf("out/h_%s_A%.1f_%s.csv", D, A, o));
               ip(r.t,r.Y,40.0)-ip(cf.t,cf.Y,40.0)) for o in ("pos","neg") ]
        v = [ (r = FP.readrun(@sprintf("out/f_%s_A%.1f_%s.csv", D, A, o));
               ip(r.t,r.Y,40.0)-ip(cf.t,cf.Y,40.0)) for o in ("down","up") ]
        push!(hv, sum(h)/2); push!(vv, sum(v)/2)
    end
    plot!(p4, collect(AS), hv; color=:royalblue, marker=mk, markersize=5, linestyle=:solid,
          label="\$g_x\$, \$Da=10^{$(D[end-1:end])}\$")
    plot!(p4, collect(AS), vv; color=:black, marker=mk, markersize=5, linestyle=:dash,
          label="\$g_y\$, \$Da=10^{$(D[end-1:end])}\$")
end
savefig(p4, "report/fig4.pdf")
println("figures written to report/")
