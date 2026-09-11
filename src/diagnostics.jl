# Post-processing of the CSV written by `run!`.

"Read the CSV written by `run!` into a NamedTuple of vectors."
function readrun(file)
    lines = readlines(file)
    names = Symbol.(split(lines[1], ','))
    data = [parse.(Float64, split(l, ',')) for l in lines[2:end]]
    NamedTuple{Tuple(names)}(Tuple([d[j] for d in data] for j in eachindex(names)))
end

"Unwrap an angle series so it is continuous (no jumps by 2π)."
function unwrap(φ)
    out = copy(φ)
    for i in 2:length(out)
        d = out[i] - out[i-1]
        out[i] -= 2π*round(d/(2π))
    end
    out
end

"""
    summarise(file; window=30.0, ν, L)

Means over the last `window` convective times of the horizontal and vertical
velocity and angular velocity, the realised Reynolds number on the mean
descent speed, and the ratios compared by Andersen et al. (2005). Prints ten
lines and returns a NamedTuple.
"""
function summarise(file; window=30.0, ν=nothing, L=nothing)
    r = readrun(file)
    t1 = r.t[end]; sel = r.t .≥ t1 - window
    vX, vY, ω = mean(r.vX[sel]), mean(r.vY[sel]), mean(r.omega[sel])
    φu = unwrap(r.phi); revs = (φu[end] - φu[findfirst(sel)])/(2π)
    regime = abs(revs) > 0.5 ? "tumbling" : "fluttering/steady"
    Re = (ν === nothing || L === nothing) ? NaN : abs(vY)*L/ν
    @printf("window       : t = %.1f – %.1f c/U (%d samples)\n", t1 - window, t1, count(sel))
    @printf("<vX>/U       : %8.4f\n", vX)
    @printf("<vY>/U       : %8.4f\n", vY)
    @printf("<ω>c/U       : %8.4f\n", ω)
    @printf("revolutions  : %8.2f  → %s\n", revs, regime)
    @printf("<vX>/|<vY>|  : %8.3f   (Andersen β=1/14: 2.20)\n", vX/abs(vY))
    @printf("<ω>c/|<vY>|  : %8.3f   (Andersen β=1/14: 0.857)\n", ω/abs(vY))
    @printf("Re on <vY>   : %8.1f   (Andersen β=1/14: 1147)\n", Re)
    @printf("descent      : %.1f chords over %.1f c/U\n", -r.Y[end], t1)
    @printf("Δt (min/mean): %.3f / %.3f cells\n", minimum(r.dt), mean(r.dt))
    (; vX, vY, ω, revs, Re)
end

mean(x) = sum(x)/length(x)
