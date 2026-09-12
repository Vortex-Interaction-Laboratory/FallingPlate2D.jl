# Altitude gain Z_G of the vertical-gust study.
#   Z_G(t) ≡ Z_earth,gust(t) − Z_earth,quiescent(t)
# exactly as in Pandi et al.: the altitude difference between a plate that
# experiences the gust and one that falls in quiescent fluid. It is zero for a
# body that does not interact with the fluid, and zero before the gust onset.
using FallingPlate2D, Printf
FP = FallingPlate2D

"Linear interpolation of column `c` of run `r` onto time `t`."
function at(r, c, t)
    i = searchsortedfirst(r.t, t)
    i ≤ 1 && return getproperty(r, c)[1]
    i > length(r.t) && return getproperty(r, c)[end]
    t0, t1 = r.t[i-1], r.t[i]; y = getproperty(r, c)
    y[i-1] + (y[i] - y[i-1])*(t - t0)/(t1 - t0)
end

"Z_G(t) for a gust run against its quiescent reference."
ZG(gust, ref, t) = at(gust, :Zearth, t) - at(ref, :Zearth, t)

load(name) = FP.readrun("out/$name.csv")

function report(L; tend=30.0, t₀=15.0)
    ref, up, dn = load("D2_nogust_L$L"), load("D2_up_L$L"), load("D2_down_L$L")
    @printf("\n=== 2D permeable plate, χ=10 m=0.80 e=1.16 I_O=1.15 Da=2.4e-4 Ga=150, L=%d ===\n", L)
    u_t = -at(ref, :vY, t₀)
    @printf("terminal velocity  u_t/u_g = %.4f   (3D disk: 1.50)   Re = u_t·c/ν = %.0f  (3D disk: 225)\n",
            u_t, u_t*150)
    @printf("attitude at gust onset: φ = %+.4f rad, ω c/u_g = %+.4f  (plate is broadside and steady)\n",
            at(ref, :phi, t₀), at(ref, :omega, t₀))
    @printf("bookkeeping check, Z_G just before the gust: up %+.2e, down %+.2e (must be 0)\n",
            ZG(up, ref, t₀ - 0.5), ZG(dn, ref, t₀ - 0.5))
    println("\n   t      Z_G(up)   Z_G(down)     sum       antisym")
    for t in (16.0, 17.0, 18.0, 20.0, 22.0, 25.0, 27.5, tend)
        t > min(up.t[end], dn.t[end], ref.t[end]) && continue
        zu, zd = ZG(up, ref, t), ZG(dn, ref, t)
        @printf("%6.1f  %+9.4f  %+9.4f  %+9.4f  %+9.4f\n", t, zu, zd, zu + zd, (zu - zd)/2)
    end
    zu, zd = ZG(up, ref, tend), ZG(dn, ref, tend)
    @printf("\nAt t = %.0f:  Z_G(up) = %+.4f,  Z_G(down) = %+.4f,  SUM = %+.4f  (%s altitude gain)\n",
            tend, zu, zd, zu + zd, zu + zd > 0 ? "NET POSITIVE" : "NET NEGATIVE")
    (; L, u_t, zu, zd, sum=zu + zd)
end

function validation(; tend=35.0)
    ref, gst = load("V_nogust_L64"), load("V_gust_L64")
    println("\n=== Validation: Bose & Viola (2023), Ga=50, m=0.5, Da=1e-4, χ=10 ===")
    @printf("transverse sustained gust, GR=1, t_G=1, t₀=15\n")
    @printf("terminal velocity u_t/u_g = %.4f   (paper: u_t < 1, the plate rises in a frame falling at u_g)\n",
            -at(ref, :vY, 15.0))
    for t in (20.0, 25.0, 30.0, 35.0)
        t > min(ref.t[end], gst.t[end]) && continue
        @printf("  Z_G(t=%.0f) = %+.4f\n", t, ZG(gst, ref, t))
    end
    z = ZG(gst, ref, min(tend, ref.t[end], gst.t[end]))
    @printf("Z_G at t=35: %.3f   Bose & Viola (2023): 1.4   → %+.0f%%\n", z, 100(z - 1.4)/1.4)
    z
end

isfile("out/V_gust_L64.csv") && validation()
res = [report(L) for L in (64, 96) if isfile("out/D2_down_L$L.csv")]
if length(res) == 2
    a, b = res
    @printf("\n=== Grid verification (L=%d → %d) ===\n", a.L, b.L)
    @printf("%-12s %10s %10s %10s\n", "", "L=$(a.L)", "L=$(b.L)", "change")
    for (k, lab) in ((:u_t, "u_t/u_g"), (:zu, "Z_G(up)"), (:zd, "Z_G(down)"), (:sum, "sum"))
        x, y = getproperty(a, k), getproperty(b, k)
        @printf("%-12s %10.4f %10.4f %+10.4f\n", lab, x, y, y - x)
    end
end
