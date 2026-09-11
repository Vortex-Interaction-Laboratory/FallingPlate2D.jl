# Does the convective exit blow up under an accelerating far field, with the plate held?
using FallingPlate2D, WaterLily, StaticArrays, Printf
const Tramp = 4.0   # ramp duration, c/U
ramp(t) = 0.6*min(t/Tramp, 1.0)          # far-field speed (t in c/U)
function held(fp, name; tend=9.0)
    sim = fp.sim; t0 = time(); n = 0; ok = true
    while sim_time(sim) < tend
        t = sim_time(sim)
        fp.Ufar[] = SA[ramp(t), 0.0]
        fp.afluid[] = SA[t < Tramp ? 0.6/(Tramp*fp.L) : 0.0, 0.0]   # dU/dt in grid units
        sim.flow.Δt[end] = min(sim.flow.Δt[end], 0.5)
        sim_step!(sim; remeasure=false); n += 1
        umax = maximum(abs, sim.flow.u)
        if umax > 5 || sim.flow.Δt[end] < 1e-3
            I = argmax([sqrt(sum(abs2, @view sim.flow.u[I,:])) for I in CartesianIndices(size(sim.flow.p))])
            @printf("%s: BLEW UP at t=%.2f  umax=%.1f at %s\n", name, t, umax, Tuple(I)); ok = false; break
        end
    end
    ok && @printf("%s: stable to t=%.1f, umax=%.2f, %d steps, %.0f s\n", name, tend, maximum(abs, sim.flow.u), n, time()-t0)
end
println("(i)  Ref-based far field + explicit g, exitBC=true")
held(falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, free=(false,false,false)), "(i)")
println("(iii) Ref-based far field + explicit g, exitBC=false (Dirichlet exit)")
fp = falling_plate(; χ=4, ρ=2.7, Re=500, L=32, dom=(8,8), φ₀=0.1, free=(false,false,false))
L = fp.L; body = fp.sim.body
sim3 = Simulation((8L, 8L), (i,x,t)->fp.Ufar[][i], L; U=1.0, ν=fp.ν, g=(i,x,t)->fp.afluid[][i], body, T=Float64, perdir=(2,), exitBC=false)
held(FallingPlate2D.FallingPlateSim(sim3, fp.state, fp.dyn, fp.section, fp.Ufar, fp.afluid, fp.x₀, L, 1.0, fp.ν, fp.g, 0.5, fp.free, ), "(iii)")
println("(ii) native uBC(i,x,t) ramp, WaterLily computes g by autodiff, exitBC=true")
Lr = 32.0
uBC(i,x,t) = i==1 ? 0.6*min(t/(Tramp*Lr), 1.0) : 0.0
sim2 = Simulation((8Lr, 8Lr), uBC, Lr; U=1.0, ν=fp.ν, body, T=Float64, perdir=(2,), exitBC=true)
t0 = time(); n = 0
while sim_time(sim2) < 9.0
    sim2.flow.Δt[end] = min(sim2.flow.Δt[end], 0.5); sim_step!(sim2; remeasure=false); global n += 1
    if maximum(abs, sim2.flow.u) > 5 || sim2.flow.Δt[end] < 1e-3
        @printf("(ii): BLEW UP at t=%.2f umax=%.1f\n", sim_time(sim2), maximum(abs, sim2.flow.u)); break
    end
end
sim_time(sim2) ≥ 9.0 && @printf("(ii): stable to t=9, umax=%.2f, %d steps, %.0f s\n", maximum(abs, sim2.flow.u), n, time()-t0)
