using JLD2, Printf, Statistics
# Reader for the dt battery: failure rates at dt=0.05 against the dt=0.1 baseline,
# and the three-step Trotter series of the plateau (caches under data/local/controls/).
# Run:  julia --project=. scripts/analysis/dtreport.jl
ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const C = joinpath(ROOT, "data", "local", "controls")
loose(x) = 0.05 < x < 0.20
base = Dict((0.0,18.0)=>(1,20), (0.0,22.0)=>(2,20), (0.0,24.0)=>(4,20),
            (0.1,11.0)=>(2,10), (0.1,12.0)=>(0,10), (0.3,8.0)=>(4,10), (0.3,9.0)=>(6,10))
println("=== failure rate: dt=0.05 vs dt=0.1 baseline ===")
@printf("%-6s %-5s %-14s %-14s %s\n","p","T","dt=0.1","dt=0.05","med(good) dt=0.05")
for f in sort(filter(x->startswith(x,"seedens_dt005_"), readdir(C)))
    m = match(r"seedens_dt005_p([0-9.]+)_T([0-9.]+)\.jld2", f); m===nothing && continue
    P, T = parse(Float64,m[1]), parse(Float64,m[2])
    r = load(joinpath(C,f),"res"); ps=[v.plateau for v in values(r)]
    g = filter(loose, ps)
    b = get(base,(P,T),(-1,0))
    @printf("%-6.1f %-5.0f %-14s %-14s %s\n", P, T,
            b[1]>=0 ? "$(b[1])/$(b[2])" : "--",
            "$(count(!loose,ps))/$(length(ps))",
            isempty(g) ? "--" : @sprintf("%.4f", median(g)))
end
if isfile(joinpath(C,"trotter_fine.jld2"))
    tf = load(joinpath(C,"trotter_fine.jld2"),"res")
    println("\n=== plateau vs dt (Richardson check; dt=0.1 / 0.05 from earlier ladders) ===")
    for k in sort(collect(keys(tf)))
        @printf("p=%.1f T=%-4.0f  dt=0.025: plateau=%.4f  (%s@%d, %.0fs)\n",
                k[1], k[2], tf[k].plateau, tf[k].reason, tf[k].niters, tf[k].elapsed)
    end
end
