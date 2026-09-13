# Trajectories of two permeable plates branched from the same state at t = 20:
# one at constant gravity, one whose gravitational force is modulated over a
# single period of two convective times (0.5 g → g → 1.5 g → g).
#   julia --project=examples examples/plot_gravity_lap.jl
using Plots, FallingPlate2D, Printf
FP = FallingPlate2D

_a  = FP.readrun("out/lap1.csv")            # t = 0 → 20
_b  = FP.readrun("out/lap1_ext.csv")        # t = 20 → 50 (continued from the checkpoint)
_s  = _b.t .> _a.t[end]
lap = (; t=vcat(_a.t,_b.t[_s]), X=vcat(_a.X,_b.X[_s]), Y=vcat(_a.Y,_b.Y[_s]),
         phi=vcat(_a.phi,_b.phi[_s]), vX=vcat(_a.vX,_b.vX[_s]), vY=vcat(_a.vY,_b.vY[_s]),
         omega=vcat(_a.omega,_b.omega[_s]))
c1  = FP.readrun("out/end1.csv")      # constant gravity
c2  = FP.readrun("out/end2.csv")      # modulated gravity
m, T_G, t₀ = 0.80, 2.0, 50.0
g   = 1/m
gm  = gravity_modulation(; amplitude=0.5, T_G, t₀, g)
gratio(t) = 1 + gust_acceleration(gm, t)/g

# join the shared history to each branch so the trajectories are complete
cat2(a, b, s) = (vcat(getproperty(a, s), getproperty(b, s)))
T1, X1, Y1, V1 = cat2(lap,c1,:t), cat2(lap,c1,:X), cat2(lap,c1,:Y), cat2(lap,c1,:vY)
T2, X2, Y2, V2 = cat2(lap,c2,:t), cat2(lap,c2,:X), cat2(lap,c2,:Y), cat2(lap,c2,:vY)

interp(t, x, q) = (i = searchsortedfirst(t, q); i ≤ 1 ? x[1] : i > length(t) ? x[end] :
                   x[i-1] + (x[i]-x[i-1])*(q-t[i-1])/(t[i]-t[i-1]))
ts = range(t₀, min(c1.t[end], c2.t[end]); length=400)
ΔY = [interp(T2,Y2,q) - interp(T1,Y1,q) for q in ts]
ΔX = [interp(T2,X2,q) - interp(T1,X1,q) for q in ts]

default(fontfamily="Computer Modern", linewidth=1.6, framestyle=:box, grid=false, dpi=200)
const CONST, MOD, SH = :navy, :firebrick, :gray65
band = (t₀, t₀ + T_G)

# (a) trajectory
pa = plot(xlabel="\$X/c\$", ylabel="\$Y/c\$", legend=:bottomleft, aspect_ratio=:equal)
plot!(pa, lap.X, lap.Y; color=SH, label="shared, \$t<50\$")
plot!(pa, c1.X, c1.Y; color=CONST, label="constant \$g\$")
plot!(pa, c2.X, c2.Y; color=MOD, linestyle=:dash, label="modulated \$g\$")
for (r, col) in ((c1, CONST), (c2, MOD)), k in 1:50:length(r.t)
    x, y, φ = r.X[k], r.Y[k], r.phi[k]
    plot!(pa, [x-0.5cos(φ), x+0.5cos(φ)], [y-0.5sin(φ), y+0.5sin(φ)]; color=col, linewidth=1, label=false)
end

# (b) vertical velocity, with the modulation window shaded
pb = plot(xlabel="\$t\\,u_g/c\$", ylabel="\$v_Y/u_g\$", legend=:bottomleft)
vspan!(pb, [band...]; color=:gray, alpha=0.15, label="")
plot!(pb, T1, V1; color=CONST, label="constant \$g\$")
plot!(pb, T2, V2; color=MOD, linestyle=:dash, label="modulated \$g\$")
pbt = twinx(pb)
plot!(pbt, ts, gratio.(ts); color=:black, linestyle=:dot, linewidth=1.2,
      ylabel="\$g_\\mathrm{eff}/g\$", label=false, legend=false, grid=false)

# (c) separation of the two plates
pc = plot(xlabel="\$t\\,u_g/c\$", ylabel="separation \$/c\$", legend=:bottomleft)
vspan!(pc, [band...]; color=:gray, alpha=0.15, label="")
plot!(pc, ts, ΔY; color=:black, label="\$\\Delta Y\$")
plot!(pc, ts, ΔX; color=SH, linestyle=:dash, label="\$\\Delta X\$")
hline!(pc, [0]; color=:gray, linewidth=0.8, label="")

plot(pa, pb, pc; layout=@layout([a{0.28w} b c]), size=(1500, 470),
     left_margin=7Plots.mm, bottom_margin=7Plots.mm, right_margin=7Plots.mm)
savefig("out/gravity_lap.png"); println("saved out/gravity_lap.png")

φ₀ = lap.phi[end]; ω₀ = lap.omega[end]
@printf("\nbranch point t=%.1f:  vY=%.4f  Y/c=%.3f  φ=%+.4f rad  ω=%+.4f (pendulum phase at onset)\n",
        lap.t[end], lap.vY[end], lap.Y[end], φ₀, ω₀)
@printf("at t=%.0f:  Y const=%.4f  Y mod=%.4f  ΔY=%+.4f c\n", ts[end],
        interp(T1,Y1,ts[end]), interp(T2,Y2,ts[end]), ΔY[end])
@printf("final vY:  const=%.4f  mod=%.4f  (difference %+.2e)\n",
        interp(T1,V1,ts[end]), interp(T2,V2,ts[end]), interp(T2,V2,ts[end])-interp(T1,V1,ts[end]))
