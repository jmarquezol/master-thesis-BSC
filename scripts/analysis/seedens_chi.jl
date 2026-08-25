ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra, Statistics

# Does a larger bond-dimension cap recover the failed entropy rungs? Same seeds as the chi=64
# ensembles (identical Random.seed! scheme), cap doubled to 128.
# usage: julia seedens_chi.jl P NSEEDS T1 [T2 ...]
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))
const P      = parse(Float64, ARGS[1])
const NSEEDS = parse(Int, ARGS[2])
const TS     = parse.(Float64, ARGS[3:end])
const CACHE  = joinpath(ROOT, "data", "local", "controls")
plateau(prof) = mean(imag.(prof)[max(1, end ÷ 2 - 1):end ÷ 2 + 2])

for T in TS
    out = joinpath(CACHE, "seedens_chi128_p$(P)_T$(T).jld2")
    res = isfile(out) ? load(out, "res") : Dict{Int,Any}()
    println("\n=== chi=128  p=$P T=$T  ($(length(res)) of $NSEEDS cached)")
    flush(stdout)
    for s in 1:NSEEDS
        haskey(res, s) && continue
        Random.seed!(round(Int, 1_000_000P) + 1000 * round(Int, 10T) + s)   # same seed as chi=64
        t0 = time()
        r = run_pm_diagnosed(T; p=P, maxdim=128, cutoff=1e-12, nbeta=4,
                             itermax=4000, stuck_after=300, column=:bulk5)
        s2 = collect(ITransverse.gen_renyi2(r.L, r.R))
        res[s] = (; seed=s, chi=128, plateau=plateau(s2[3:end-2]),
                    rigidity=1 / (norm(r.L) * norm(r.R)), lambda0=r.lambda0,
                    reason=r.reason, niters=r.niters, elapsed=time() - t0, s2=s2)
        jldsave(out; res)
        @printf("  seed %2d  plateau=%8.4f  |mu0|=%.6f  (%s @%d)  %.0f s\n",
                s, res[s].plateau, abs(r.lambda0), r.reason, r.niters, res[s].elapsed)
        flush(stdout)
    end
end
println("\nSEEDENS-CHI-DONE")
