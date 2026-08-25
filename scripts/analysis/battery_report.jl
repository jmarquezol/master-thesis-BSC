using JLD2, Printf, Statistics

# Prints the three tables the failure battery produces: the seed-ensemble spread per rung
# (seedens_* caches), the exact diagnostics on the corrected column (dense_bulk5 cache), and
# the cutoff scan (cutrerun_* caches), all under data/local/controls/.
# Run:  julia --project=. scripts/analysis/battery_report.jl  [detail]
ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CACHE = joinpath(ROOT, "data", "local", "controls")

mad(v) = 1.4826 * median(abs.(v .- median(v)))
caches(prefix) = sort(filter(f -> startswith(f, prefix), readdir(CACHE)))

# A seed fails when its plateau cannot be a Renyi-2 temporal plateau at all. Two bands around
# pi*c/16 = 0.0982 make the sensitivity of the rate visible, the way fit windows are reported
# elsewhere. Ensemble dispersion is not used: it under-counts when failures set the median and
# over-counts when the survivors agree to four decimals.
loose(x) = 0.05 < x < 0.20
tight(x) = 0.07 < x < 0.15

function seedens_table()
    println("=== seed ensembles (single-vector route) ===")
    @printf("%-26s %-4s %-11s %-11s %-11s %-10s %s\n",
            "cache", "n", "fail loose", "fail tight", "med(good)", "spread", "rigidity")
    for f in caches("seedens_")
        res = load(joinpath(CACHE, f), "res")
        ps = [v.plateau for v in values(res)]
        rg = [v.rigidity for v in values(res)]
        good = filter(loose, ps)
        n = length(ps)
        @printf("%-26s %-4d %-11s %-11s %-10s %-10.2e %.2e\n", f, n,
                "$(count(!loose, ps))/$n", "$(count(!tight, ps))/$n",
                isempty(good) ? "   --  " : @sprintf("%.4f", median(good)),
                mad(ps), median(rg))
    end
end

function seedens_detail()
    println("\n=== per-seed plateaus ===")
    for f in caches("seedens_")
        res = load(joinpath(CACHE, f), "res")
        ks = sort(collect(keys(res)))
        println(f)
        println("  plateau: ", join([@sprintf("%.4f", res[k].plateau) for k in ks], " "))
        println("  reason : ", join([string(res[k].reason, "@", res[k].niters) for k in ks], " "))
    end
end

function dense_table()
    file = joinpath(CACHE, "dense_bulk5.jld2")
    isfile(file) || return
    res = load(file, "res")
    println("\n=== exact dense diagnostics, bulk5 column ===")
    @printf("%-5s %-5s %-6s %-8s %-10s %-10s %-10s %-10s %s\n",
            "p", "T", "dim", "|mu0|", "|mu1/mu0|", "gap01", "rigidity", "cond(V)", "supp(trunc)")
    for k in sort(collect(keys(res)), by = x -> (x[1], x[2]))
        v = res[k]
        @printf("%-5.1f %-5.1f %-6d %-8.4f %-10.5f %-10.3e %-10.3e %-10.3e %.2e\n",
                v.p, v.T, v.dim, abs(v.mu0), v.ratio, v.gap01, v.r0, v.condV,
                v.pert["trunc"].suppression)
    end
end

function cutoff_table()
    rows = [load(joinpath(CACHE, f), "res") for f in caches("cutrerun_")]
    isempty(rows) && return
    println("\n=== cutoff scan (block route, seeded) ===")
    @printf("%-10s %-5s %-9s %-5s %-7s %-11s %-8s %-9s %s\n",
            "p", "T", "cutoff", "seed", "niters", "reason", "time_s", "g", "plateau")
    for v in sort(rows, by = x -> (x.p, x.T, x.cutoff, x.seed))
        label = string(v.p, v.bc == "Up" ? "(fixed)" : "")
        @printf("%-10s %-5.1f %-9.0e %-5d %-7d %-11s %-8.0f %-9.5f %.4f\n",
                label, v.T, v.cutoff, v.seed, v.niters, v.reason, v.elapsed, v.g, v.plateau)
    end
end

seedens_table()
dense_table()
cutoff_table()
"detail" in ARGS && seedens_detail()
