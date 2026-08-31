ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra
# Convergence history of the single-vector iteration for fig:conv: the per-iteration change ds
# under RTM (which levels off above the tolerance and ends on the no-improvement rule) and RDM
# (which reaches it). Same rung, same seed, production settings; only the truncation differs.
# The seed matches run 1 of the seedens ensemble at this rung.
# Writes: data/local/conv_history_sv.jld2
# Run:    julia --project=. scripts/analysis/conv_history_sv.jl
BLAS.set_num_threads(2)
const P, T, SEEDN = 0.0, 12.0, 1
const OUT = joinpath(ROOT, "data", "local", "conv_history_sv.jld2")
res = isfile(OUT) ? load(OUT, "res") : Dict{Symbol,Any}()
for (mode, alg) in ((:rtm, "RTM"), (:rdm, "densitymatrix"))
    haskey(res, mode) && continue
    Random.seed!(round(Int, 1_000_000P) + 1000 * round(Int, 10T) + SEEDN)
    el = @elapsed r = run_pm_diagnosed(T; p=P, maxdim=64, cutoff=1e-12, nbeta=4,
                                       itermax=4000, stuck_after=300, column=:bulk5, alg=alg)
    res[mode] = (ds=collect(r.ds_hist), chi=collect(r.chi_hist), reason=String(r.reason),
                 niters=r.niters, lambda0=r.lambda0, elapsed=el)
    jldsave(OUT; res=res)
    @printf("%s: %s after %d iters, chi=%d, lambda0=%.10f%+.10fim, %.0f s\n",
            mode, r.reason, r.niters, maximum(r.chi_hist), real(r.lambda0), imag(r.lambda0), el)
end
if haskey(res, :rtm) && haskey(res, :rdm)
    d = abs(res[:rtm].lambda0 - res[:rdm].lambda0)
    @printf("|lambda0(RTM) - lambda0(RDM)| = %.2e  (rel %.2e)\n", d, d / abs(res[:rtm].lambda0))
end
println("CONV-SV-DONE")
