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
    m::Float64      # mass = ρ·A
    A::Float64      # section area (displaced fluid mass)
    Is::Float64     # moment of inertia about the centre
    m11::Float64    # added mass along the chord
    m22::Float64    # added mass along the normal
    Ia::Float64     # added moment of inertia
    g::SVector{2,Float64}
end
function PlateDynamics(s::Section, ρ, g)
    m11, m22, Ia = added_mass(s)
    PlateDynamics(ρ*area(s), area(s), ρ*inertia(s), m11, m22, Ia, SA[0.0, -g])
end

# grid (x, y) = (−Y, X)
to_grid(v::SVector{2}) = SA[-v[2], v[1]]
to_lab(w::SVector{2}) = SA[w[2], -w[1]]

"Chord and normal unit vectors in the lab frame for body angle θ."
chord_dir(θ) = SA[sin(θ), -cos(θ)]
normal_dir(θ) = SA[cos(θ), sin(θ)]
chord_angle(θ) = θ - π/2      # φ, from the horizontal

"Added-mass tensor in the lab frame."
function added_mass_tensor(d::PlateDynamics, θ)
    ec, en = chord_dir(θ), normal_dir(θ)
    d.m11*ec*ec' + d.m22*en*en'
end

"""
    step!(s::RigidState, d::PlateDynamics, F, M, Δt)

Advance the body by `Δt` given the lab-frame fluid force `F` on the body and
the moment `M` about its centre (anticlockwise positive), both as computed by
the flow solver after the previous step. Semi-implicit Euler in time.
"""
function step!(s::RigidState, d::PlateDynamics, F::SVector{2}, M::Real, Δt;
               free=(true, true, true))
    Ma = added_mass_tensor(d, s.θ)
    mask = SA[free[1] ? 1.0 : 0.0, free[2] ? 1.0 : 0.0]
    s.a = mask .* ((d.m*I + Ma) \ (F + Ma*s.a + (d.m - d.A)*d.g))
    s.α = free[3] ? (M + d.Ia*s.α)/(d.Is + d.Ia) : 0.0
    s.v += Δt*s.a;  s.X += Δt*s.v
    s.ω += Δt*s.α;  s.θ += Δt*s.ω
    s
end
