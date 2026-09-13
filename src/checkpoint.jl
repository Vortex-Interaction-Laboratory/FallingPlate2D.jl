# Restart files, so that one run can be branched at a saved instant.
#
# Stored: the velocity and pressure fields, the accumulated time-step history
# (so `sim_time` is preserved exactly), the rigid-body state, and the frame
# bookkeeping (far-field velocity, frame body force, convective-exit profile).
# Not stored: the BDIM fields `μ₀`, `μ₁` and `V`, which are functions of the
# body position and are rebuilt by `prepare_flow!` on load.

"""
    save_checkpoint(fp::FallingPlateSim, fname) -> fname

Write a restart file for `fp`.
"""
function save_checkpoint(fp::FallingPlateSim, fname)
    d = dirname(fname); isempty(d) || mkpath(d)
    s = fp.state
    jldsave(fname;
        p = Array(fp.sim.flow.p), u = Array(fp.sim.flow.u), Δt = collect(fp.sim.flow.Δt),
        X = Vector(s.X), v = Vector(s.v), a = Vector(s.a), θ = s.θ, ω = s.ω, α = s.α,
        uexit = copy(fp.uexit), Ufar = Vector(fp.Ufar[]), afluid = Vector(fp.afluid[]),
        grid = collect(size(fp.sim.flow.p)), L = fp.L, Da = fp.Da, ν = fp.ν, g = fp.g)
    fname
end

"""
    load_checkpoint!(fp::FallingPlateSim, fname) -> fp

Restore a checkpoint into `fp`, which must have been built with the same grid
and the same physical parameters. The **gust need not match**: branching one
saved state under different gravity histories is the intended use.
"""
function load_checkpoint!(fp::FallingPlateSim, fname)
    jldopen(fname) do o
        size(fp.sim.flow.p) == Tuple(o["grid"]) || throw(ArgumentError(
            "grid $(size(fp.sim.flow.p)) does not match the checkpoint $(Tuple(o["grid"]))"))
        for (nm, have, want) in (("L", fp.L, o["L"]), ("Da", fp.Da, o["Da"]), ("ν", fp.ν, o["ν"]))
            isapprox(have, want) || @warn "checkpoint $nm differs" have want
        end
        fp.sim.flow.p .= o["p"]
        fp.sim.flow.u .= o["u"]
        empty!(fp.sim.flow.Δt); append!(fp.sim.flow.Δt, o["Δt"])
        s = fp.state
        s.X = SVector{2,Float64}(o["X"]); s.v = SVector{2,Float64}(o["v"])
        s.a = SVector{2,Float64}(o["a"]); s.θ = o["θ"]; s.ω = o["ω"]; s.α = o["α"]
        fp.uexit .= o["uexit"]
        fp.Ufar[] = SVector{2,Float64}(o["Ufar"])
        fp.afluid[] = SVector{2,Float64}(o["afluid"])
    end
    fp.sim.body = setmap(fp.sim.body; θ=fp.state.θ, ω=fp.state.ω)
    prepare_flow!(fp)
    fp
end
