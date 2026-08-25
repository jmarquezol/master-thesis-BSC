ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Random, LinearAlgebra, Statistics, LsqFit

# Seed ensembles for the single-vector (entropy) route: the same rung repeated from independent
# reproducible seeds, so the median of the physical runs is an ensemble reading rather than one
# draw from the attractor lottery. One cache per rung, so lanes never share a file.
#
# With --dt005 the same rungs run at half the time step (nbeta doubled to hold beta0 fixed), which
# answers whether the failure rate moves with the Trotter step. The seeds are offset so the two
# ensembles draw independently.
#
# usage: julia --project=. scripts/analysis/seedens.jl [--dt005] P NSEEDS T1 [T2 ...]
# writes: data/local/controls/seedens[_dt005]_p<P>_T<T>.jld2
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

const FINE   = !isempty(ARGS) && ARGS[1] == "--dt005"
const REST    = FINE ? ARGS[2:end] : ARGS
const P      = parse(Float64, REST[1])
const NSEEDS = parse(Int, REST[2])
const TS     = parse.(Float64, REST[3:end])
const CACHE  = joinpath(ROOT, "data", "local", "controls")

# halving dt doubles nbeta so that beta0 = nbeta*dt/2 stays at 0.2, and doubles the number of
# cooling bonds trimmed from each end of the profile
const DT     = FINE ? 0.05 : 0.10
const NBETA  = FINE ? 8 : 4
const TRIM   = NBETA ÷ 2
const TAG    = FINE ? "seedens_dt005" : "seedens"
const OFFSET = FINE ? 500_000_000 : 0

W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]
function chord_c(prof, T)
    re = real.(prof); n = length(re)
    ts = range(T / (n + 1), T - T / (n + 1), length=n); bulk = (n ÷ 4):(3n ÷ 4)
    return 8 * curve_fit(lin, W.(collect(ts[bulk]), T), re[bulk], [0.06, 0.5]).param[1]
end
plateau(prof) = mean(imag.(prof)[max(1, end ÷ 2 - 1):end ÷ 2 + 2])

for T in TS
    out = joinpath(CACHE, "$(TAG)_p$(P)_T$(T).jld2")
    res = isfile(out) ? load(out, "res") : Dict{Int,Any}()
    println("\n=== dt=$DT p=$P T=$T  ($(length(res)) of $NSEEDS cached)")
    flush(stdout)
    for s in 1:NSEEDS
        haskey(res, s) && continue
        Random.seed!(OFFSET + round(Int, 1_000_000P) + 1000 * round(Int, 10T) + s)
        t0 = time()
        r = run_pm_diagnosed(T; p=P, maxdim=64, cutoff=1e-12, dt=DT, nbeta=NBETA,
                             itermax=4000, stuck_after=300, column=:bulk5)
        s2 = collect(ITransverse.gen_renyi2(r.L, r.R))
        prof = s2[(TRIM + 1):(end - TRIM)]
        res[s] = (; seed=s, dt=DT, nbeta=NBETA, c_Re=chord_c(prof, T), plateau=plateau(prof),
                    rigidity=1 / (norm(r.L) * norm(r.R)), lambda0=r.lambda0,
                    reason=r.reason, niters=r.niters, elapsed=time() - t0, s2=s2)
        jldsave(out; res)
        @printf("  seed %2d  plateau=%7.4f  |mu0|=%.6f  rig=%.2e  (%s @%d)  %.0f s\n",
                s, res[s].plateau, abs(r.lambda0), res[s].rigidity, r.reason, r.niters, res[s].elapsed)
        flush(stdout)
    end
end
println("\nSEEDENS-DONE")
