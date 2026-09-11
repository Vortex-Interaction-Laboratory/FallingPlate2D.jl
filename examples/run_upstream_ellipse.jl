# Phase A check: run the maintainers' TwoD_FreeRotatingEllipse.jl UNMODIFIED
# in the `examples/` environment (WaterLily + ParametricBodies + Plots).
# Run:  julia -t 4 --project=examples examples/run_upstream_ellipse.jl
using Downloads
url = "https://raw.githubusercontent.com/WaterLily-jl/WaterLily-Examples/master/examples/TwoD_FreeRotatingEllipse.jl"
path = Downloads.download(url, joinpath(mktempdir(), "TwoD_FreeRotatingEllipse.jl"))
println("Running upstream example from ", url)
t = @elapsed include(path)     # runs main() (tend=50) and main(Re=250,thk=0.12,pivot=0.25,tend=15)
println("Upstream example completed in ", round(t/60, digits=1), " min (includes compilation and gif rendering)")
