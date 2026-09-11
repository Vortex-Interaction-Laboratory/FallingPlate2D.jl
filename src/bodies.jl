# Plate sections: signed-distance functions and inertial properties.
# Body coordinates: x₁ along the chord, x₂ along the thickness, origin at the
# centre. Lengths in cells; fluid density 1. `a` is the half-chord, `b` the
# half-thickness, so the chord-to-thickness ratio is χ = a/b.

abstract type Section end

"Elliptical section with semi-axes `a` (half-chord) and `b` (half-thickness)."
struct Ellipse <: Section
    a::Float64
    b::Float64
end
Ellipse(L::Real, χ::Real) = Ellipse(L/2, L/(2χ))

"Rectangle `2a × 2b` with semicircular ends of radius `b` (a capsule)."
struct RoundedPlate <: Section
    a::Float64
    b::Float64
end
RoundedPlate(L::Real, χ::Real) = RoundedPlate(L/2, L/(2χ))

chord(s::Section) = 2s.a
thickness(s::Section) = 2s.b

"""
    sdf(s::Section) -> (x, t) -> distance

Signed distance in body coordinates, negative inside. The ellipse uses the
first-order normalised implicit function, which is exact to within a fraction
of a cell for χ ≤ 20 and adequate for the 2-cell BDIM kernel.
"""
sdf(e::Ellipse) = (x, t) -> begin
    f = sqrt((x[1]/e.a)^2 + (x[2]/e.b)^2) - 1
    g = sqrt((x[1]/e.a^2)^2 + (x[2]/e.b^2)^2)
    f/g
end
sdf(p::RoundedPlate) = (x, t) -> begin
    h = p.a - p.b                       # half-length of the straight part
    sqrt(max(abs(x[1]) - h, 0)^2 + x[2]^2) - p.b
end

"Section area (mass per unit density)."
area(e::Ellipse) = π*e.a*e.b
area(p::RoundedPlate) = 4p.a*p.b - (4 - π)*p.b^2

"Second moment of area about the centre (moment of inertia per unit density)."
inertia(e::Ellipse) = π/4*e.a*e.b*(e.a^2 + e.b^2)
function inertia(p::RoundedPlate)
    h, b = p.a - p.b, p.b
    rect = 4h*b*(h^2 + b^2)/3
    semi = π*b^4/4 + π*b^2*h^2/2 + 4h*b^3/3   # each semicircular end, about the centre
    rect + 2semi
end

"""
    added_mass(s) -> (m₁₁, m₂₂, Iₐ)

Potential-flow added mass along the chord, along the normal, and added moment
of inertia about the centre, per unit fluid density (Lamb, ellipse). The
rounded plate uses the ellipse values with the same `a, b` as an estimate; the
implicit coupling only needs an estimate.
"""
added_mass(s::Section) = (π*s.b^2, π*s.a^2, π/8*(s.a^2 - s.b^2)^2)

"""
    andersen_Istar(s, ρ) -> (actual, rectangular)

Dimensionless moment of inertia `I* = 32 I/(π ρ_f c⁴)` of Andersen, Pesavento
& Wang (2005): first with the actual section inertia, then with the
rectangular-plate formula `(8/3π) ρ β (1 + β²)`, `β = 1/χ`, which is how the
values in their table 1 were computed.
"""
function andersen_Istar(s::Section, ρ)
    c = chord(s); β = s.b/s.a
    32ρ*inertia(s)/(π*c^4), 8/(3π)*ρ*β*(1 + β^2)
end
