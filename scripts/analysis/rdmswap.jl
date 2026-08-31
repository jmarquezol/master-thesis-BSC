ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using Printf, JLD2, LinearAlgebra, Statistics, Random, LsqFit

# Scheme-swap test: rerun rungs with density-matrix (RDM) truncation instead of RTM, same
# reproducible seeds as the seedens RTM ensemble. If the same rungs fail either way, the
# truncation scheme is not the mechanism; if failures move, it is implicated.
# usage: julia rdmswap.jl P NSEEDS T1 [T2 ...]
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

const P      = parse(Float64, ARGS[1])
const NSEEDS = parse(Int, ARGS[2])
const TS     = parse.(Float64, ARGS[3:end])

W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]
function chord_c(prof, T)
    re = real.(prof); n = length(re)
    ts = range(T / (n + 1), T - T / (n + 1), length=n); bulk = (n ÷ 4):(3n ÷ 4)
    return 8 * curve_fit(lin, W.(collect(ts[bulk]), T), re[bulk], [0.06, 0.5]).param[1]
end
plateau(prof) = mean(imag.(prof)[max(1, end ÷ 2 - 1):end ÷ 2 + 2])

out = joinpath(ROOT, "data", "local", "controls", "rdmswap_p$(P).jld2")
res = isfile(out) ? load(out, "res") : Dict{Any,Any}()
for T in TS, s in 1:NSEEDS
    haskey(res, (T, s)) && continue
    Random.seed!(round(Int, 1_000_000P) + 1000 * round(Int, 10T) + s)   # matches seedens
    t0 = time()
    r = run_pm_diagnosed(T; p=P, maxdim=64, cutoff=1e-12, nbeta=4,
                         itermax=4000, stuck_after=300, column=:bulk5, alg="densitymatrix")
    s2 = collect(ITransverse.gen_renyi2(r.L, r.R))
    prof = s2[3:end-2]
    res[(T, s)] = (; seed=s, c_Re=chord_c(prof, T), plateau=plateau(prof),
                     rigidity=1 / (norm(r.L) * norm(r.R)), lambda0=r.lambda0,
                     reason=r.reason, niters=r.niters, elapsed=time() - t0, s2=s2)
    jldsave(out; res)
    @printf("rdm p=%.1f T=%-5.1f seed %d  c_Re=%8.4f  plateau=%7.4f  |mu0|=%.6f  (%s @%d)  %.0f s\n",
            P, T, s, res[(T, s)].c_Re, res[(T, s)].plateau, abs(r.lambda0),
            r.reason, r.niters, res[(T, s)].elapsed)
    flush(stdout)
end
println("RDMSWAP-DONE")
