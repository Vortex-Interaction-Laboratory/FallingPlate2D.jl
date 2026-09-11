using FallingPlate2D, WaterLily, StaticArrays, Printf
FP = FallingPlate2D
L = 32; sec = Ellipse(L, 4); ν = L/500; θ = π/2; x0 = SA[3L*1.0, 4L*1.0]
function profile(sim, xc, label)
    j = round(Int, xc[2] + 1.5); i0 = round(Int, xc[1] + 1.5)
    @printf("%s: p along the motion axis (offset from centre, p, μ₀):\n", label)
    for di in -10:2:10
        i = i0 + di
        @printf("   %+3d  p=%9.4f  μ₀=%.2f\n", di, sim.flow.p[i, j], sim.flow.μ₀[i, j, 1])
    end
    F = -SVector{2}(WaterLily.pressure_force(sim)); Fv = -SVector{2}(WaterLily.viscous_force(sim))
    @printf("   F_p=(%8.3f,%8.3f)  F_v=(%8.3f,%8.3f)\n", F..., Fv...)
end
# (A) fixed plate, stream +0.5 along x (plate broadside: chord along y)
simA = Simulation((8L, 8L), (0.5, 0.0), L; U=1.0, ν, body=AutoBody(FP.sdf(sec), RigidMap(x0, θ)), T=Float64)
for n in 1:150; simA.flow.Δt[end] = min(simA.flow.Δt[end], 0.5); sim_step!(simA; remeasure=false); end
profile(simA, x0, "(A) fixed plate in stream +0.5: expect p>0 upstream (−x side), F_x>0")
# (B) plate moving at +0.5 through still fluid
V = SA[0.5, 0.0]; x₀ = x0
simB = Simulation((8L, 8L), (0.0, 0.0), L; U=1.0, ν, body=AutoBody(FP.sdf(sec), RigidMap(x0, θ; V)), T=Float64)
for n in 1:150
    Δt = min(simB.flow.Δt[end], 0.5); simB.flow.Δt[end] = Δt
    global x₀ = x₀ + Δt*V; simB.body = setmap(simB.body; x₀, V); sim_step!(simB; remeasure=true)
end
profile(simB, x₀, "(B) plate moving +0.5 in still fluid: expect p>0 ahead (+x side), F_x<0")
