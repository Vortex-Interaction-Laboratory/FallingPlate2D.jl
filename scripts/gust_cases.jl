# Shared definitions for the vertical-gust study.
using FallingPlate2D, Printf
FP = FallingPlate2D

# 2D equivalent of the 3D permeable disk (see REPORT): χ=10, m=0.80, e=1.16,
# I_O=1.15, Da=2.4e-4, Ga=150.
const D2 = (; χ=10.0, m=0.80, e=1.16, I_O=1.15, Da=2.4e-4, Ga=150.0)
# Validation case of Bose & Viola (2023) fig. 4a: Ga=50, m=0.5, Da=1e-4, χ=10.
const V  = (; χ=10.0, m=0.50, e=0.0,  I_O=nothing, Da=1e-4,   Ga=50.0)

function case(p; L=64, dom=(16, 16), φ₀=0.1, gust=NOGUST, shape=:plate, Δt_max=0.5)
    kw = p.I_O === nothing ? (;) : (; I_O=p.I_O)
    falling_plate(; χ=p.χ, m=p.m, e=p.e, Da=p.Da, Ga=p.Ga, L, dom, xbody=dom[1]/2,
                  φ₀, shape, gust, Δt_max, kw...)
end

function header(fp, name)
    @printf("%s: χ=%g m=%g e=%g Da=%g Ga=%g | g=%.4f (=1/m) ν=%.4g grid %s, %.1f cells/thickness\n",
        name, fp.section.a*2/(fp.section.b*2), fp.dyn.m/fp.L^2, fp.dyn.e/fp.L, fp.Da, fp.U*fp.L/fp.ν,
        fp.g*fp.L, fp.ν/fp.L, size(fp.sim.flow.p), 2fp.section.b)
    fp.gust.GR == 0 || @printf("   gust: dir=%s GR=%g t₀=%g T_G=%g %s, peak u̇_G=%.3f (= %.2f g)\n",
        fp.gust.dir, fp.gust.GR, fp.gust.t₀, fp.gust.T_G,
        fp.gust.isolated ? "isolated" : "sustained", peak_acceleration(fp.gust),
        peak_acceleration(fp.gust)/(fp.g*fp.L))
    flush(stdout)
end
