# Extending the Eq. (3) fitting window past the reach of the block iteration.
#
# The fit consumes only the leading eigenvalue, and the single-vector iteration returns exactly
# that at a fraction of the cost. So the ladder can be continued with local single-vector rungs.
#
# The difficulty is that at those times independent seeds land on genuinely different fixed points,
# so an average over the ensemble is meaningless. We pick one seed per rung by continuity of the
# phase advance, which is the same criterion Section sec:reach applies across evolution times,
# here applied across seeds. The criterion never refers to the CFT prediction.
#
# A rung is only kept if the fit residual does not degrade, which is what rejects p=0.5 at T=7.
#
# Reads:  data/cluster/sweep_*.jld2          the block-iteration arms
#         data/local/controls/seedens_p*.jld2  the single-vector seed ensembles
# Run:    julia --project=. scripts/analysis/extend_window.jl

ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

const CL  = joinpath(ROOT, "data", "cluster")
const LOC = joinpath(ROOT, "data", "local", "controls")
const V   = Dict(0.0 => 2.000, 0.1 => 2.670, 0.3 => 3.967, 0.5 => 5.212)  # v_inf from equilibrium

# the block-iteration windows of Table tab:cp, where the phase increment leaves its constant rate
const BLOCK_WINDOW = Dict(0.0 => 20.0, 0.1 => 9.5, 0.3 => 7.0, 0.5 => 6.0)
const MAX_WINDOW   = Dict(0.0 => 22.0, 0.1 => 17.0, 0.3 => 10.0, 0.5 => 6.0)

# ── loading ──────────────────────────────────────────────────────────────────
arm(file, label) = Dict(k[2] => v for (k, v) in load(joinpath(CL, file), "done")
                        if k[1] == label && !haskey(v, :error))

function cluster_arm(p)
    p == 0.0 && return merge(arm("sweep_rtm_eigs_p0.0.jld2", "rtm_eigs_p0.0"),
                             arm("sweep_rtm_eigs_p0.0_fine.jld2", "rtm_eigs_p0.0_fine"),
                             arm("sweep_rtm_eigs_p0.0_fineb.jld2", "rtm_eigs_p0.0_fineb"))
    p == 0.1 && return merge(arm("sweep_rtm_eigs_p0.1_bulk.jld2", "rtm_eigs_p0.1_bulk"),
                             arm("sweep_rtm_eigs_p0.1_fine_bulk.jld2", "rtm_eigs_p0.1_fine_bulk"),
                             arm("sweep_tower_p0.1_bulk.jld2", "tower_p0.1_bulk"),
                             arm("sweep_rtm_eigs_p0.1_fineb_bulk.jld2", "rtm_eigs_p0.1_fineb_bulk"))
    p == 0.3 && return merge(arm("sweep_ent_p0.3_bulk.jld2", "ent_p0.3_bulk"),
                             arm("sweep_tower_p0.3_bulk.jld2", "tower_p0.3_bulk"),
                             arm("sweep_rtm_eigs_p0.3_fine_bulk.jld2", "rtm_eigs_p0.3_fine_bulk"),
                             arm("sweep_rtm_eigs_p0.3_fineb_bulk.jld2", "rtm_eigs_p0.3_fineb_bulk"))
    return merge(arm("sweep_ent_p0.5_bulk.jld2", "ent_p0.5_bulk"),
                 arm("sweep_tower_p0.5_bulk.jld2", "tower_p0.5_bulk"),
                 arm("sweep_rtm_eigs_p0.5_fine_bulk.jld2", "rtm_eigs_p0.5_fine_bulk"),
                 arm("sweep_rtm_eigs_p0.5_fineb_bulk.jld2", "rtm_eigs_p0.5_fineb_bulk"))
end

# every seed ensemble available for this coupling, as T => vector of mu0
function seed_ensembles(p)
    out = Dict{Float64,Vector{ComplexF64}}()
    for f in readdir(LOC)
        m = match(Regex("^seedens_p$(p)_T([0-9.]+)\\.jld2\$"), f)
        m === nothing && continue
        res = load(joinpath(LOC, f), "res")
        out[parse(Float64, m.captures[1])] = [r.lambda0 for r in values(res)]
    end
    return out
end

# ── the phase ladder ─────────────────────────────────────────────────────────
phase_of(mu) = angle(-mu)

# unwrap so the ladder is continuous rather than folded into (-pi, pi]
function unwrap(phases)
    out = copy(phases)
    for i in 2:length(out)
        d = out[i] - out[i-1]
        out[i] = out[i-1] + d - 2pi * round(d / (2pi))
    end
    return out
end

# put a bare phase on the same branch as a reference value
same_branch(phi, reference) = phi - 2pi * round((phi - reference) / (2pi))

# Straight line through the last few accepted points, extrapolated to T. Two points is the
# default because it uses only the local rate of the phase advance, which is the quantity the
# criterion is about; longer baselines average over the curvature the fit is trying to measure.
function predict_phase(Ts, phis, T; nback=2)
    n = min(nback, length(Ts))
    x = Ts[end-n+1:end]
    y = phis[end-n+1:end]
    n == 1 && return y[1]
    slope = (n * sum(x .* y) - sum(x) * sum(y)) / (n * sum(abs2, x) - sum(x)^2)
    intercept = mean(y) - slope * mean(x)
    return intercept + slope * T
end

# ── the Eq. (3) fit ──────────────────────────────────────────────────────────
eq3(x, q) = @. q[1] - pi / x + q[2] / x^2

function fit_c(p, Ts, phis)
    y = phis ./ Ts
    f = curve_fit(eq3, Ts, y, [0.1, -0.01])
    c = 24 * V[p] * abs(f.param[2]) / pi
    rms = sqrt(mean(abs2, y .- eq3(Ts, f.param)))
    return c, rms
end

# ── per coupling: build the block ladder, then extend it seed by seed ────────
# A rung is only usable if there are enough seeds to choose between. With one run there is no
# selection at all, only acceptance of whatever it returned, so such rungs are not admitted.
const MIN_SEEDS = 3

function extend(p; rms_tolerance=3.0, nback=2, verbose=true)
    a = cluster_arm(p)
    block_Ts = sort([T for T in keys(a) if T <= BLOCK_WINDOW[p]])
    Ts = copy(block_Ts)
    phis = unwrap([phase_of(a[T].theta_phys) for T in block_Ts])
    mus = ComplexF64[a[T].theta_phys for T in block_Ts]

    c_block, rms_block = fit_c(p, Ts, phis)
    verbose && @printf("\np = %.1f\n", p)
    verbose && @printf("  block iteration   T <= %-5g %2d pts   c = %.4f   rms = %.1e\n",
            BLOCK_WINDOW[p], length(Ts), c_block, rms_block)

    ens = seed_ensembles(p)
    # T = 24 at p = 0 is excluded from the thesis: its ensemble loses 4 of 20 runs and its median
    # plateau turns back upwards, so the p = 0 sequences end at T = 22 (see MAX_WINDOW).
    candidates = sort([T for T in keys(ens) if T > BLOCK_WINDOW[p] && T <= MAX_WINDOW[p]])
    isempty(candidates) && return (; p, Ts, phis, mus, c=c_block, rms=rms_block, window=BLOCK_WINDOW[p])

    accepted_window = BLOCK_WINDOW[p]
    c_ext, rms_ext = c_block, rms_block
    for T in candidates
        predicted = predict_phase(Ts, phis, T; nback=nback)
        # the seed whose phase best continues the ladder, never the one closest to a prediction
        seed_phases = [same_branch(phase_of(mu), predicted) for mu in ens[T]]
        misses = abs.(seed_phases .- predicted)
        best = argmin(misses)

        if length(ens[T]) < MIN_SEEDS
            verbose && @printf("  T = %-5g %2d seeds   too few to select from, stopping here\n",
                               T, length(ens[T]))
            break
        end
        trial_Ts = [Ts; T]
        trial_phis = [phis; seed_phases[best]]
        c_try, rms_try = fit_c(p, trial_Ts, trial_phis)

        keep = rms_try <= rms_tolerance * rms_ext
        # A rung is only usable if one seed clearly continues the ladder and the others do not.
        # Where that separation collapses to order one, no seed is on the branch.
        separation = median(misses) / misses[best]
        verbose && @printf("  T = %-5g %2d seeds   best %.1e  typical %.1e  (x%-5.0f)   c = %.4f   rms = %.1e   %s\n",
                T, length(ens[T]), misses[best], median(misses), separation, c_try, rms_try,
                keep ? "keep" : "REJECT")
        keep || break                       # the ladder stops at the first rung that fails

        Ts, phis = trial_Ts, trial_phis
        push!(mus, ens[T][best])
        c_ext, rms_ext, accepted_window = c_try, rms_try, T
    end

    verbose && @printf("  extended          T <= %-5g %2d pts   c = %.4f   rms = %.1e\n",
            accepted_window, length(Ts), c_ext, rms_ext)
    return (; p, Ts, phis, mus, c=c_ext, rms=rms_ext, window=accepted_window)
end

results = [extend(p) for p in (0.0, 0.1, 0.3, 0.5)]

# Write the selected ladder so the figure and the audit consume the same rungs rather than
# repeating the selection. Re-running this script after the seed ensembles grow refreshes both.
const OUT = joinpath(ROOT, "data", "local", "extended_rungs.jld2")
selected = Dict(r.p => (Ts=r.Ts, phis=r.phis, mus=r.mus, window=r.window, c=r.c, rms=r.rms)
                for r in results)
jldsave(OUT; selected=selected)
println("\nwrote ", relpath(OUT, ROOT))

println("\n", "="^74)
@printf("%-5s %-14s %-14s %-9s %-9s\n", "p", "block window", "extended", "c block", "c ext")
for r in results
    a = cluster_arm(r.p)
    bTs = sort([T for T in keys(a) if T <= BLOCK_WINDOW[r.p]])
    cb, _ = fit_c(r.p, bTs, unwrap([phase_of(a[T].theta_phys) for T in bTs]))
    @printf("%-5.1f T<=%-11g T<=%-11g %-9.4f %-9.4f\n",
            r.p, BLOCK_WINDOW[r.p], r.window, cb, r.c)
end

# How much does the answer depend on the extrapolation baseline? This is the selector's own
# systematic, and it belongs with the number.
println("\nsensitivity of the extended value to the extrapolation baseline")
@printf("%-5s %-9s %-9s %-9s %-9s\n", "p", "2 points", "3", "4", "5")
for p in (0.0, 0.1, 0.3, 0.5)
    cs = [extend(p; nback=n, verbose=false).c for n in 2:5]
    @printf("%-5.1f %-9.4f %-9.4f %-9.4f %-9.4f   spread %.4f\n",
            p, cs..., maximum(cs) - minimum(cs))
end
