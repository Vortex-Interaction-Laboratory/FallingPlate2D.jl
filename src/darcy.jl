# Permeable (Darcy–Brinkman) plate through the BDIM mobility.
#
# In a Brinkman medium the fluid inside the solid obeys
#     ∂u/∂t + … = −∇p + ν∇²u − σ (u − V),   σ = ν/K,
# with K the permeability and V the solid velocity. Treated implicitly over a
# step, the update  u(1 + Δt σ) = u* − Δt∇p + Δt σ V  reads
#     u = μ (u* − Δt∇p) + (1 − μ) V,   μ = 1/(1 + Δt σ),
# which is exactly WaterLily's BDIM blend with mobility μ₀ = μ: the same
# coefficient also enters the pressure Poisson operator, so pressure acts on
# the pore fluid with the reduced mobility. A permeable body is therefore
# obtained by raising the interior mobility from 0 to μ_p = 1/(1 + Δt ν/K),
# blended across the BDIM kernel:  μ₀ ← μ₀ + (1 − μ₀) μ_p.  Da → 0 recovers
# the impermeable BDIM body exactly; Da → ∞ removes the body.
# The second-order BDIM correction μ₁ models a sharp no-slip interface and is
# zeroed for a permeable body. Darcy number Da = K/c².

"""
    apply_darcy!(flow, Da, L)

Raise the BDIM mobility inside the body to the implicit Brinkman value for
Darcy number `Da` (permeability `K = Da·L²`) and zero `μ₁`. No-op for `Da ≤ 0`.
Call after `measure!` and the interior velocity fill, before `update!(pois)`.
"""
function apply_darcy!(flow, Da::Real, L::Real)
    Da > 0 || return nothing
    σ = flow.ν/(Da*L^2)
    μp = 1/(1 + flow.Δt[end]*σ)
    μ₀ = flow.μ₀
    @inbounds for Ii in eachindex(μ₀)
        μ₀[Ii] += (1 - μ₀[Ii])*μp
    end
    # restore zero mobility on the domain boundary faces (WaterLily sets these
    # after measure!; raising them lets pressure flux leak through the walls)
    T = eltype(μ₀); D = ndims(flow.p)
    WaterLily.BC!(μ₀, zeros(SVector{D,T}), false, flow.perdir)
    flow.μ₁ .= 0
    nothing
end

"Brinkman screening length `c·√Da` in cells, which the grid must resolve."
screening_length(Da, L) = L*sqrt(Da)
