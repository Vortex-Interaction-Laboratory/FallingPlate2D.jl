# Trajectory and time series from a run! CSV.
#   julia --project=examples examples/plot_run.jl out/name.csv
using Plots, FallingPlate2D
file = length(ARGS) ≥ 1 ? ARGS[1] : "out/andersen_dns_b8_L64.csv"
r = FallingPlate2D.readrun(file)
default(fontfamily="Computer Modern", linewidth=1.2, legend=:topright, framestyle=:box, grid=false, dpi=150)
# (a) trajectory with the chord drawn every 0.5 c/U
p1 = plot(r.X, r.Y; color=:black, label=false, xlabel="X/c", ylabel="Y/c", aspect_ratio=:equal)
for k in 1:10:length(r.t)
    x, y, φ = r.X[k], r.Y[k], r.phi[k]
    plot!(p1, [x - 0.5cos(φ), x + 0.5cos(φ)], [y - 0.5sin(φ), y + 0.5sin(φ)]; color=:red, label=false)
end
# (b) time series
φu = FallingPlate2D.unwrap(r.phi)
p2 = plot(r.t, φu/(2π); label="φ/2π", xlabel="t U/c", ylabel="")
plot!(p2, r.t, r.vX; label="vX/U")
plot!(p2, r.t, r.vY; label="vY/U")
plot!(p2, r.t, r.omega; label="ω c/U")
plot(p1, p2; layout=(1, 2), size=(1200, 520), left_margin=5Plots.mm, bottom_margin=5Plots.mm)
out = replace(file, ".csv" => ".png"); savefig(out); println("saved ", out)
