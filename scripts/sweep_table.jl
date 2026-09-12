# Summary table of the Darcy sweep at the Andersen inertia (Re_g = 1900, L = 64).
using FallingPlate2D, Printf
FP = FallingPlate2D
cases = [("0", 0.0), ("1e-6", 1e-6), ("1e-5", 1e-5), ("3e-5", 3e-5), ("1e-4", 1e-4)]
ν = 64/1900; L = 64.0
@printf("%-6s %8s %8s %8s %8s %8s %8s %8s  %s\n", "Da", "Re·Da", "t_end", "<vY>/U", "<vX>/U", "<ω>c/U", "revs", "Re_vY", "regime")
for (tag, Da) in cases
    f = "out/sweep_Re1900_Da$tag.csv"; isfile(f) || continue
    r = FP.readrun(f); t1 = r.t[end]; w = min(30.0, 0.6t1); sel = r.t .≥ t1 - w
    vX, vY, ω = FP.mean(r.vX[sel]), FP.mean(r.vY[sel]), FP.mean(r.omega[sel])
    φu = FP.unwrap(r.phi); revs = (φu[end] - φu[findfirst(sel)])/(2π)
    σφ = sqrt(FP.mean((r.phi[sel] .- FP.mean(r.phi[sel])).^2))
    regime = abs(revs) > 0.5 ? "tumbling" : σφ > 0.2 ? "fluttering" : "steady broadside"
    @printf("%-6s %8.4f %8.1f %8.3f %8.3f %8.3f %8.2f %8.0f  %s\n", tag, 1900Da, t1, vY, vX, ω, revs, abs(vY)*L/ν, regime)
end
