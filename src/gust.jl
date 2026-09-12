# Gust model of Pandi et al. (JFM 2026) and Bose & Viola (2023).
#
# The flow is solved in a frame that moves with the gust, so that the fluid is
# at rest at infinity and the uniform force associated with the gust, like
# gravity, is absorbed into the pressure. For a body with negligible solid
# volume (a permeable, filamentous plate) that force exerts no buoyancy, and
# the gust is felt by the body only as the transient force −m u̇_G added to the
# weight: the gust is a transient change of the gravitational force.
#
# Gust velocity, following eq. (5) of Bose & Viola (2023): a cosine ramp from
# zero to `GR·u_g` over the gust period `T_G`, beginning at `t₀`. A `sustained`
# gust holds `GR` thereafter (their case). An `isolated` gust applies the same
# cosine ramp up and back down, so the fluid returns to rest: the net impulse
# of the gust force is then zero, and the altitude gain `Z_G` stays finite for
# a vertical gust.

"""
    Gust(; dir, GR, t₀, T_G, isolated)

`dir` is the gust direction in the lab frame (`[1,0]` transverse, `[0,1]`
upward), `GR` the gust ratio (gust speed over gravitational velocity, negative
for a gust along `−dir`), `t₀` the onset time and `T_G` the gust period, both
in convective times `c/u_g`.
"""
struct Gust
    dir::SVector{2,Float64}
    GR::Float64
    t₀::Float64
    T_G::Float64
    isolated::Bool
end
function Gust(; dir=SA[0.0, 1.0], GR=0.0, t₀=15.0, T_G=2.0, isolated=true)
    d = SVector{2,Float64}(dir)
    Gust(d/sqrt(sum(abs2, d)), GR, t₀, T_G, isolated)
end
const NOGUST = Gust(; GR=0.0)

"Gust speed `u_G(t)/u_g`."
function gust_velocity(g::Gust, t)
    τ = (t - g.t₀)/g.T_G
    τ ≤ 0 && return 0.0
    if g.isolated
        τ ≥ 1 && return 0.0
        return g.GR/2*(1 - cos(2π*τ))
    else
        τ ≥ 1 && return g.GR
        return g.GR/2*(1 - cos(π*τ))
    end
end

"Gust acceleration `u̇_G(t)`, the transient change of gravity felt by the body."
function gust_acceleration(g::Gust, t)
    τ = (t - g.t₀)/g.T_G
    (τ ≤ 0 || τ ≥ 1) && return 0.0
    g.isolated ? g.GR*π/g.T_G*sin(2π*τ) : g.GR*π/(2g.T_G)*sin(π*τ)
end

"Displacement of the fluid, `∫u_G dt`, which carries the body with it."
function gust_displacement(g::Gust, t)
    τ = (t - g.t₀)/g.T_G
    τ ≤ 0 && return 0.0
    if g.isolated
        τ ≥ 1 && return g.GR*g.T_G/2
        return g.GR*g.T_G/2*(τ - sin(2π*τ)/(2π))
    else
        τ ≥ 1 && return g.GR*g.T_G*(τ - 0.5)
        return g.GR*g.T_G/2*(τ - sin(π*τ)/π)
    end
end

"Peak gust acceleration."
peak_acceleration(g::Gust) = g.isolated ? abs(g.GR)*π/g.T_G : abs(g.GR)*π/(2g.T_G)
