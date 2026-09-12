# Vertical-gust study: all cases, run sequentially in one process.
#   julia -t 4 --project=. scripts/gust_runs.jl [L]
include(joinpath(@__DIR__, "gust_cases.jl"))
mkpath("out")
L = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 64
const DOM = (12, 12)

function go(name, p; gust=NOGUST, tend=30.0, kw...)
    file = "out/$name.csv"
    if isfile(file) && countlines(file) > 100
        println("skip $name (exists)"); flush(stdout); return
    end
    fp = case(p; L, dom=DOM, gust, kw...)
    header(fp, name)
    st = run!(fp; tend, out=file, every=0.02, report=5.0)
    @printf("%s: %s\n\n", name, st); flush(stdout)
end

# Production: the 2D equivalent of the permeable disk. Isolated vertical gust,
# GR = ±1, T_G = 2 (two of the paper's cosine ramps), onset at t₀ = 15.
up   = Gust(; dir=[0, 1], GR= 1.0, t₀=15.0, T_G=2.0, isolated=true)
down = Gust(; dir=[0, 1], GR=-1.0, t₀=15.0, T_G=2.0, isolated=true)
go("D2_nogust_L$L", D2)
go("D2_up_L$L",     D2; gust=up)
go("D2_down_L$L",   D2; gust=down)

# Validation: Bose & Viola (2023) fig. 4a, Ga=50, m=0.5, Da=1e-4, transverse
# sustained gust GR=1, t_G=1, t₀=15 → Z_G ≈ 1.4 at t=35.
if L == 64
    tg = Gust(; dir=[1, 0], GR=1.0, t₀=15.0, T_G=1.0, isolated=false)
    go("V_nogust_L$L", V; tend=35.0)
    go("V_gust_L$L",   V; gust=tg, tend=35.0)
end
