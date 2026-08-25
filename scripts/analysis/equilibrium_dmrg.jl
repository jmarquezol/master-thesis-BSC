# Equilibrium ground-state calculations of Section 4.1: the three DMRG caches behind the
# central-charge figures.
#
#   nb4_cvsp.jld2       c(p) from a mid-chain fit at N = 300, on a coarse grid of couplings
#   nb4_fss.jld2        mid-chain S1 and S2 at N = 40..200, for the finite-size scaling of fig:cft_L
#   nb4_chord_N400.jld2 full entanglement profile at N = 400, for fig:cft_chord
#
# Run:  julia --project=. scripts/analysis/equilibrium_dmrg.jl
#
# This recomputes physics and is slow (DMRG up to N = 400 with maxdim 1000). Each cache is written
# as soon as its part finishes and existing caches are skipped, so an interrupted run resumes. The
# finished caches ship with the repository; nothing here needs rerunning to reproduce the thesis.

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, LsqFit

const LAMBDA = 1.0
const LOCAL = joinpath(ROOT, "data", "local")

# one ground state, with the sweep schedule used throughout
function ground_state(N, p; maxdim=[10, 20, 100, 200, 400], nsweeps=15)
    sites = siteinds("S=1/2", N)
    _, psi = dmrg(alcaraz_H(sites, LAMBDA, p), randomMPS(sites, 10);
                  nsweeps=nsweeps, maxdim=maxdim, cutoff=[1e-10],
                  noise=[1e-4, 1e-5, 1e-6, 0.0], outputlevel=0)
    return psi
end

# von Neumann and Renyi-2 entropies from the Schmidt spectrum at the middle bond
function mid_chain_entropies(psi, N)
    mid = div(N, 2)
    orthogonalize!(psi, mid)
    _, Sv, _ = svd(psi[mid], (linkind(psi, mid - 1), siteind(psi, mid)))
    probs = [Sv[n, n]^2 for n in 1:dim(Sv, 1) if Sv[n, n]^2 > 1e-12]
    return -sum(pr * log(pr) for pr in probs), -log(sum(pr^2 for pr in probs))
end

# ── 1. c(p) over a coarse grid of couplings ──────────────────────────────────
lin(x, q) = @. q[1] * x + q[2]
cvsp_cache = joinpath(LOCAL, "nb4_cvsp.jld2")
if !isfile(cvsp_cache)
    N = 300
    ps, cs = Float64[], Float64[]
    for p in 0.0:0.2:2.0
        S = compute_vn_entropy(ground_state(N, p))
        # fit the middle half, on same-parity cuts so the even-odd oscillation drops out
        quarter = div(N, 4)
        lo = iseven(quarter) ? quarter : quarter + 1
        cuts = lo:2:div(3N, 4)
        xv = [log((N / pi) * sin(l * pi / N)) / 6 for l in cuts]
        fit = curve_fit(lin, xv, S[cuts], [1.0, 1.0])
        push!(ps, p); push!(cs, fit.param[1])
        @printf("c(p) sweep: p=%.1f  c=%.3f\n", p, fit.param[1])
    end
    jldsave(cvsp_cache; p=ps, c=cs)
    println("wrote nb4_cvsp.jld2")
end

# ── 2. finite-size scaling of the mid-chain entropies ────────────────────────
fss_cache = joinpath(LOCAL, "nb4_fss.jld2")
if !isfile(fss_cache)
    N_values = [40, 60, 80, 100, 120, 140, 160, 180, 200]
    fss_data = Dict{Float64,Any}()
    for p in [0.1, 0.5]
        S1, S2 = Float64[], Float64[]
        for N in N_values
            s1, s2 = mid_chain_entropies(ground_state(N, p), N)
            push!(S1, s1); push!(S2, s2)
            @printf("fss p=%.1f  N=%3d  S1=%.4f  S2=%.4f\n", p, N, s1, s2)
        end
        fss_data[p] = (vn=S1, r2=S2)
    end
    jldsave(fss_cache; data=fss_data)
    println("wrote nb4_fss.jld2")
end

# ── 3. full entanglement profile at N = 400 ──────────────────────────────────
chord_cache = joinpath(LOCAL, "nb4_chord_N400.jld2")
N_chord = 400
chord_data = isfile(chord_cache) ? load(chord_cache, "data") : Dict{Float64,Vector{Float64}}()
for p in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0, 1.5]
    haskey(chord_data, p) && continue
    @printf("chord profile at p=%.1f (N=%d) ...\n", p, N_chord)
    chord_data[p] = compute_vn_entropy(
        ground_state(N_chord, p; maxdim=[20, 100, 200, 400, 800, 1000], nsweeps=50))
    jldsave(chord_cache; data=chord_data)      # checkpoint after every coupling
end

println("equilibrium caches complete")
