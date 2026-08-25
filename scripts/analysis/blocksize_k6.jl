# Block-size control: the same gap ladder as gap_ladder.jl, run with a block of 6 states instead
# of 4. The gap is a physical quantity, so it must not depend on how many states the block carries.
# Where the two agree, the k=4 ladder is resolving the tower properly; where they part company, the
# smaller block is cutting through it, which is what the wiggles in the k=4 ladder are.
#
# Writes: data/local/nb3_gap_k6.jld2, keyed by (:k6, p, T)
# Run:    julia --project=. scripts/analysis/blocksize_k6.jl
#
# This recomputes physics and takes hours. The finished cache ships with the repository. Every
# point is saved as it finishes, so an interrupted run resumes.
#
# The shipped cache also holds two older entry types, (:exact, p) and (:seed, n), from controls
# that have since been redone properly: the dense spectra are superseded by the corrected-column
# ones in notebook 06, and the seeds by the seed ensembles there. This script leaves them alone.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra
BLAS.set_num_threads(2)

# k=6 costs far more than k=4, so each coupling is followed only as far as it stayed affordable
const LADDERS = [(0.0, collect(2.0:1.0:7.0)),
                 (0.1, collect(2.0:1.0:6.0)),
                 (0.3, collect(2.0:1.0:4.0)),
                 (0.5, collect(2.0:1.0:4.0))]
const CACHE = joinpath(ROOT, "data", "local", "nb3_gap_k6.jld2")

done = isfile(CACHE) ? load(CACHE, "done") : Dict{Any,Any}()

for (p, Ts) in LADDERS
    prevL = prevR = nothing
    for T in Ts
        key = (:k6, p, T)
        if haskey(done, key)
            @printf("p=%.1f T=%.0f cached\n", p, T)
            prevL = prevR = nothing        # a cached point stores no vectors to warm-start from
            continue
        end

        mpo, scaffold = build_alcaraz_tmpo(T; p=p, lambda=1.0, dt=0.1, nbeta=4, MPO_alg="VD2")
        sites = siteinds(scaffold)
        seedL = prevL === nothing ? nothing : MPS[pad_tmps(w, sites) for w in prevL]
        seedR = prevR === nothing ? nothing : MPS[pad_tmps(w, sites) for w in prevR]

        theta, L, R, info = block_transfer_eigs(mpo, scaffold;
            k=6, maxdim=64, maxdims=collect(2:2:64),
            cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
            itermax=1000, eps_conv=1e-6, trunc_mode=:rtm,
            n_track=2, stuck_after=120, seedL=seedL, seedR=seedR)

        done[key] = (theta=collect(theta), reason=string(info[:reason]), niters=info[:niters])
        jldsave(CACHE; done=done)
        @printf("p=%.1f T=%.0f  |mu0|=%.6f  %s in %d iters\n",
                p, T, maximum(abs.(theta)), info[:reason], info[:niters])
        prevL, prevR = L, R
        GC.gc()
    end
end

# ── the comparison this control exists for ──
# lambda0 and its +/- partner are the two largest moduli; the gap is the next one down.
gap_of(theta) = begin
    order = sortperm(abs.(theta), rev=true)
    rest = order[3:end]
    isempty(rest) ? 0.0 : maximum(abs.(theta[rest])) / abs(theta[order[1]])
end

k4 = load(joinpath(ROOT, "data", "local", "nb55_pgap.jld2"), "D")
println()
@printf("%-5s %-4s %-10s %-10s %-9s\n", "p", "T", "gap k=4", "gap k=6", "diff")
for (p, Ts) in LADDERS, T in Ts
    haskey(done, (:k6, p, T)) || continue
    g4, g6 = k4[p][T].gap, gap_of(done[(:k6, p, T)].theta)
    @printf("%-5.1f %-4.0f %-10.4f %-10.4f %-9.1e\n", p, T, g4, g6, abs(g4 - g6))
end
println("\nwrote ", relpath(CACHE, ROOT))
