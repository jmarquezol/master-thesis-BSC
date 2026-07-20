# cluster/merge_sweeps.jl — assemble the per-label sweep_*.jld2 files (v3, one per arm, written
# race-free) into the single results/data/cluster/warm_sweep.jld2 that the notebooks (NB7/NB9
# cluster-ready cells) read. Run after the arms land:  julia --project=. cluster/merge_sweeps.jl
#
# Idempotent and non-destructive: it reads every sweep_*.jld2, unions their (label,T) entries, and
# writes warm_sweep.jld2. A pre-existing warm_sweep.jld2 is folded in too (older entries kept unless
# a sweep_* file has a newer one for the same key — sweep_* wins).

using JLD2, Printf

const CLUSTER_DIR = joinpath(@__DIR__, "..", "results", "data", "cluster")
const OUT         = joinpath(CLUSTER_DIR, "warm_sweep.jld2")

merged = Dict{Tuple{String,Float64},Any}()

# Build warm_sweep.jld2 PURELY from the per-label sweep_*.jld2 files (the source of truth). We do
# NOT seed from any existing warm_sweep.jld2 — that would drag in stale/clobbered v2 data. The
# per-label files persist, so rebuilding from them is idempotent.
sweep_files = sort(filter(f -> startswith(basename(f), "sweep_") && endswith(f, ".jld2"),
                          readdir(CLUSTER_DIR; join=true)))
isempty(sweep_files) && error("no sweep_*.jld2 files in $CLUSTER_DIR — nothing to merge")
for f in sweep_files
    d = load(f, "done")
    for (k, v) in d
        merged[k] = v
    end
    @printf("  + %-28s (%d entries)\n", basename(f), length(d))
end

jldsave(OUT; done=merged)

# summary by label
labels = sort(unique(first.(keys(merged))))
println("\nmerged → $(OUT)")
for lab in labels
    ts   = sort([k[2] for k in keys(merged) if k[1] == lab])
    okts = [t for t in ts if haskey(merged[(lab, t)], :theta)]
    errs = [t for t in ts if haskey(merged[(lab, t)], :error)]
    @printf("  %-14s : %2d ok  %s%s\n", lab, length(okts), "T=$(Int.(okts))",
            isempty(errs) ? "" : "  ERR@$(Int.(errs))")
end
