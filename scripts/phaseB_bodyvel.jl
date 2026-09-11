using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
# (a) imposed body velocity at the plate tips for prescribed rotation ω > 0 (anticlockwise)
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.0, free=(false,false,false))
sim = fp.sim; ω = 0.02; θ = fp.state.θ            # θ₀ = π/2: chord along grid y
sim.body = setmap(sim.body; θ, ω); measure!(sim)
x₀ = fp.x₀
for r in (SA[0.0, 12.0], SA[0.0, -12.0], SA[6.0, 0.0])        # tip +y, tip -y, and a point off the face
    I = CartesianIndex(round.(Int, x₀ + r .+ 1.5)...)          # cell whose centre is at x₀ + r
    Vnum = SA[sim.flow.V[I, 1], sim.flow.V[I, 2]]
    d = FP.sdf(fp.section)(fp.sim.body.map(x₀ + r, 0.0), 0.0)
    @printf("r=%s  V_num=(%8.4f,%8.4f)  ω×r=(%8.4f,%8.4f)  sdf=%6.2f\n", r, Vnum..., (ω*SA[-r[2], r[1]])..., d)
end
# (b) plate translating at V = (0.5, 0) through still fluid: force must oppose the motion
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.0, free=(false,false,false))
sim = fp.sim; V = SA[0.5, 0.0]; x₀ = fp.x₀
for n in 1:200
    Δt = min(sim.flow.Δt[end], 0.5); sim.flow.Δt[end] = Δt
    global x₀ = x₀ + Δt*V
    sim.body = setmap(sim.body; x₀, V)
    sim_step!(sim; remeasure=true)
    n % 50 == 0 && (F = -SVector{2}(WaterLily.total_force(sim)); @printf("n=%3d  x₀=%s  F_on_body=(%8.3f,%8.3f)  (expect Fx < 0)\n", n, round.(x₀, digits=1), F...))
end
