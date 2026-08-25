ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using Printf, JLD2, LinearAlgebra, Statistics, Random

# A third time step for the plateau, so that with the dt = 0.1 and dt = 0.05 ladders there are
# three points to extrapolate in dt^2 rather than two. One seeded run per rung at dt = 0.025, with
# nbeta raised to 16 so that beta0 stays at 0.2. Read by dtreport.jl.
#
# Only the first rung (p = 0, T = 12) is in the shipped cache: at this step the chain is four times
# longer than at dt = 0.1 and the remaining rungs were not worth their cost. Running this script
# will start computing them.
#
# Usage:  julia --project=. scripts/analysis/trotter_fine.jl
# Writes: data/local/controls/trotter_fine.jld2
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

plateau(prof) = mean(imag.(prof)[max(1, end ÷ 2 - 1):end ÷ 2 + 2])
out = joinpath(ROOT, "data", "local", "controls", "trotter_fine.jld2")
res = isfile(out) ? load(out, "res") : Dict{Any,Any}()

for (P, T) in ((0.0, 12.0), (0.0, 16.0), (0.1, 6.0), (0.1, 8.0), (0.1, 10.0))
    haskey(res, (P, T)) && continue
    Random.seed!(700_000_000 + round(Int, 1_000_000P) + round(Int, 10T))
    t0 = time()
    r = run_pm_diagnosed(T; p=P, dt=0.025, nbeta=16, maxdim=64, cutoff=1e-12,
                         itermax=6000, stuck_after=300, column=:bulk5)
    s2 = collect(ITransverse.gen_renyi2(r.L, r.R))
    res[(P, T)] = (; dt=0.025, nbeta=16, plateau=plateau(s2[9:end-8]), lambda0=r.lambda0,
                     rigidity=1 / (norm(r.L) * norm(r.R)), reason=r.reason, niters=r.niters,
                     elapsed=time() - t0, s2=s2)
    jldsave(out; res)
    @printf("p=%.1f T=%-4.0f dt=0.025  plateau=%.4f  |mu0|=%.6f  (%s@%d)  %.0f s\n",
            P, T, res[(P, T)].plateau, abs(r.lambda0), r.reason, r.niters, res[(P, T)].elapsed)
    flush(stdout)
end
println("TROTTER-FINE-DONE")
