# Three-degree-of-freedom rigid-body dynamics of the plate with the added
# mass treated implicitly.
#
# Lab frame: X horizontal, Y up. The flow is solved in a frame translating
# with the plate centre, so the far fluid accelerates with the frame and there
# is no uniform pressure gradient at infinity. The fluid force computed by
# WaterLily after a step taken with body acceleration `a_old` is therefore
#     F = F_vort − Mₐ·a_old ,
# and the lab-frame equation of motion  (m + Mₐ)·a = F_vort + (m − A)·g
# becomes the relaxation
#     (m + Mₐ)·a_new = F + Mₐ·a_old + (m − A)·g ,
# whose fixed point is exact and whose error amplification is
# |Mₐ,num − Mₐ|/(m + Mₐ). Rotation is treated the same way with Iₛ, Iₐ.
#
# The centre of gravity G may be offset from the plate centre O by a distance
# `e` along the body normal, so that the body is bottom-heavy (a pappus with a
# seed below it, say). Translation is then integrated for G and rotation about
# G, which is exact and needs no coupling terms:
#     m·a_G = F + m·g_eff ,   I_G·α = M_G = M_O − s × F ,   s = r_G − r_O .
# The flow frame follows O, whose velocity and acceleration follow from those
# of G by rigid-body kinematics.

"Lab-frame kinematic state of the plate (cells, cells per time unit, radians)."
mutable struct RigidState
    X::SVector{2,Float64}   # position of the centre
    v::SVector{2,Float64}   # velocity
    a::SVector{2,Float64}   # acceleration (last accepted)
    θ::Float64              # body angle, anticlockwise
    ω::Float64
    α::Float64
end
RigidState(θ₀) = RigidState(zero(SVector{2,Float64}), zero(SVector{2,Float64}),
                            zero(SVector{2,Float64}), θ₀, 0.0, 0.0)

"Inertial properties (fluid density 1) and gravity vector."
struct PlateDynamics
    m::Float64      # mass
    A::Float64      # displaced fluid mass (0 when buoyancy is neglected)
    Is::Float64     # moment of inertia about the centre of gravity
    m11::Float64    # added mass along the chord
    m22::Float64    # added mass along the normal
    Ia::Float64     # added moment of inertia
    g::SVector{2,Float64}
    e::Float64      # centre of gravity offset from the plate centre, along the normal
end
function PlateDynamics(s::Section, ρ, g; e=0.0, A=area(s), Is=ρ*inertia(s), m=ρ*area(s))
    m11, m22, Ia = added_mass(s)
    PlateDynamics(m, A, Is, m11, m22, Ia, SA[0.0, -g], e)
end

# grid (x, y) = (Y, −X): the plate falls towards −x, the relative far-field
# flow is +x and leaves through the exit at x = xmax; z is unchanged.
to_grid(v::SVector{2}) = SA[v[2], -v[1]]
to_lab(w::SVector{2}) = SA[-w[2], w[1]]

"Chord and normal unit vectors in the lab frame for body angle θ (body x-axis is (cos θ, sin θ) in the grid)."
chord_dir(θ) = SA[-sin(θ), cos(θ)]
normal_dir(θ) = SA[-cos(θ), -sin(θ)]
chord_angle(θ) = θ - π/2      # φ, angle of the chord line from the horizontal, anticlockwise

"In-plane cross product ω ẑ × v."
perp(v::SVector{2}) = SA[-v[2], v[1]]

"Vector from the plate centre to the centre of gravity, lab frame."
cog_offset(s::RigidState, d::PlateDynamics) = d.e*normal_dir(s.θ)

"Velocity and acceleration of the plate centre, given those of the centre of gravity."
centre_velocity(s::RigidState, d::PlateDynamics) = s.v - s.ω*perp(cog_offset(s, d))
function centre_acceleration(s::RigidState, d::PlateDynamics)
    r = cog_offset(s, d)
    s.a - s.α*perp(r) + s.ω^2*r
end

"Added-mass tensor in the lab frame."
function added_mass_tensor(d::PlateDynamics, θ)
    ec, en = chord_dir(θ), normal_dir(θ)
    d.m11*ec*ec' + d.m22*en*en'
end

"""
    step!(s::RigidState, d::PlateDynamics, F, M, Δt; gust=zero(SVector{2,Float64}))

Advance the body by `Δt` given the lab-frame fluid force `F` on the body and
the moment `M` about the **plate centre** (anticlockwise positive), both as
computed by the flow solver after the previous step. `gust` is the gust
acceleration `u̇_G`, which the body perceives as a transient change of the
gravitational acceleration, `g_eff = g − u̇_G`. Semi-implicit Euler in time.
"""
function step!(s::RigidState, d::PlateDynamics, F::SVector{2}, M::Real, Δt;
               free=(true, true, true), gust=zero(SVector{2,Float64}))
    Ma = added_mass_tensor(d, s.θ)
    mask = SA[free[1] ? 1.0 : 0.0, free[2] ? 1.0 : 0.0]
    r = cog_offset(s, d)
    MG = M - (r[1]*F[2] - r[2]*F[1])          # moment about the centre of gravity
    geff = d.g - gust
    s.a = mask .* ((d.m*I + Ma) \ (F + Ma*s.a + (d.m - d.A)*geff))
    s.α = free[3] ? (MG + d.Ia*s.α)/(d.Is + d.Ia) : 0.0
    s.v += Δt*s.a;  s.X += Δt*s.v
    s.ω += Δt*s.α;  s.θ += Δt*s.ω
    s
end
