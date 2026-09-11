# Force on a plate translating at V=(0.5,0) through still fluid: RigidMap vs classic moving map.
using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
L = 32; sec = Ellipse(L, 4); ν = L/500; V = SA[0.5, 0.0]; x0 = SA[3L*1.0, 4L*1.0]; θ = π/2
function run(body, name; nsteps=150)
    sim = Simulation((8L, 8L), (0.0, 0.0), L; U=1.0, ν, body, T=Float64)
    x₀ = x0
    for n in 1:nsteps
        Δt = min(sim.flow.Δt[end], 0.5); sim.flow.Δt[end] = Δt
        x₀ = x₀ + Δt*V
        body isa AutoBody && body.map isa RigidMap && (sim.body = setmap(sim.body; x₀, V))
        sim_step!(sim; remeasure=true)
        if n % 50 == 0
            F = -SVector{2}(WaterLily.total_force(sim))
            # imposed velocity just ahead of the front face (chord along y, so the front face is at +x)
            I = CartesianIndex(round.(Int, x₀ + SA[5.0, 0.0] .+ 1.5)...)
            @printf("%-9s n=%3d  x₀=%6.1f  V_num ahead=(%7.3f,%7.3f)  F_on_body=(%8.3f,%8.3f)  expect Fx<0\n",
                    name, n, x₀[1], sim.flow.V[I,1], sim.flow.V[I,2], F...)
        end
    end
end
run(AutoBody(FP.sdf(sec), RigidMap(x0, θ; V)), "RigidMap")
R = SA[cos(θ) sin(θ); -sin(θ) cos(θ)]
run(AutoBody(FP.sdf(sec), (x, t) -> R*(x - x0 - V*t)), "classic")
