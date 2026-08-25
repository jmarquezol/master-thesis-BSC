# Single-vector power-method ladders: the local entropy arms used by the appendices.
#
# Each arm is one ladder in T at fixed coupling, run with the plain power method on a single
# vector rather than the k-state block. It returns the leading eigenvalue and the Renyi-2 profile
# at every rung, which is what the entropy route to c is read from. These are local runs, so they
# support appendix material only; the numbers in the main text come from the cluster arms.
#
# The archive held one small script per arm, all with the same body. They are collected here as a
# table of arms, so adding one is a single line.
#
#   p00          the integrable point, the long control ladder (c = 1/2 is known there)
#   p00_ext      its continuation to higher T
#   p00_dt005    the same points at half the time step, the Trotter control
#   p00_chi128   selected rungs at double the bond dimension, the truncation control
#   p01          the p = 0.1 ladder
#   p01_ext      its continuation
#   p01_chi128   its rungs at double the bond dimension
#
# Usage:  julia --project=. scripts/analysis/svpm_ladder.jl <arm>
# Writes: data/local/svpm_<arm>.jld2
#
# This recomputes physics and a long arm takes days. The finished caches ship with the repository.
# Every rung is saved as it finishes and a rung that throws is recorded in `failed`, so a restart
# moves past it instead of looping on it. To push an arm further, raise the end of its T range.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LinearAlgebra
BLAS.set_num_threads(parse(Int, get(ENV, "LANE_BLAS", "2")))

# arm => (p, chi, dt, nbeta, T ladder)
const ARMS = Dict(
    "p00"        => (0.0,  64, 0.10, 4, collect(2.0:1.0:24.0)),
    "p00_ext"    => (0.0,  64, 0.10, 4, collect(25.0:1.0:29.0)),
    "p00_dt005"  => (0.0,  64, 0.05, 8, collect(2.0:2.0:24.0)),
    "p00_chi128" => (0.0, 128, 0.10, 4, collect(8.0:2.0:20.0)),
    "p01"        => (0.1,  64, 0.10, 4, collect(2.0:1.0:12.0)),
    "p01_ext"    => (0.1,  64, 0.10, 4, collect(13.0:1.0:16.0)),
    "p01_chi128" => (0.1, 128, 0.10, 4, [10.0, 11.0, 12.0, 13.0, 14.0]),
)

isempty(ARGS) && error("pick an arm: " * join(sort(collect(keys(ARMS))), ", "))
const ARM = ARGS[1]
haskey(ARMS, ARM) || error("unknown arm $ARM; pick one of " * join(sort(collect(keys(ARMS))), ", "))
const P, CHI, DT, NBETA, TS = ARMS[ARM]
const CACHE = joinpath(ROOT, "data", "local", "svpm_$(ARM).jld2")

stored = isfile(CACHE) ? load(CACHE) : Dict{String,Any}()
res    = get(stored, "res", Dict{Float64,Any}())
failed = get(stored, "failed", Float64[])

@printf("arm %s: p=%.1f chi=%d dt=%.3g nbeta=%d, T=%g..%g\n", ARM, P, CHI, DT, NBETA, first(TS), last(TS))
for T in TS
    (haskey(res, T) || T in failed) && (@printf("T=%-5.1f cached\n", T); continue)
    t0 = time()
    try
        r = run_pm_diagnosed(T; p=P, maxdim=CHI, cutoff=1e-12, dt=DT, nbeta=NBETA,
                             itermax=4000, stuck_after=300, column=:bulk5)
        s2 = ITransverse.gen_renyi2(r.L, r.R)
        # 1/(|L||R|) measures how far the two fixed points have drifted apart; it is recorded as a
        # health signal, and a fall in it is where a local ladder is cut (see entropy_c.jl).
        rigidity = 1.0 / (norm(r.L) * norm(r.R))
        res[T] = (; lambda0=r.lambda0, s2=collect(s2), rigidity, reason=r.reason,
                  niters=r.niters, chi=CHI, dt=DT, nbeta=NBETA, elapsed=time() - t0)
        jldsave(CACHE; res, failed)
        @printf("T=%-5.1f |mu0|=%.6f arg=%.4f  rigidity=%.3g  (%s @%d)  %.0f s\n",
                T, abs(r.lambda0), angle(r.lambda0), rigidity, r.reason, r.niters, time() - t0)
    catch err
        push!(failed, T)                  # record it so a restart moves on instead of looping
        jldsave(CACHE; res, failed)
        @printf("RUNG-FAIL T=%.1f  %s\n", T, sprint(showerror, err))
    end
    flush(stdout)
end
println("arm ", ARM, " done")
