# cluster/wall_scan_cluster.jl — headless batch driver for the BSC master sweep.
#
# Usage:  julia --project=. cluster/wall_scan_cluster.jl <mode> [p] [Tmax]
#   <mode> ∈ {"rtm", "rdm", "cutoff", "psweep"}
#   "rtm"/"rdm"/"cutoff" take no extra args — they reproduce the original p=0.1 triple exactly
#   (T up to 14, or 12 for rdm). "psweep" takes two required extra args, <p> <Tmax>, and always
#   uses the RTM truncation route (the only one worth the cost of a p-sweep — see the RTM-vs-RDM
#   cost comparison in NB9): julia wall_scan_cluster.jl psweep 0.3 20
#
# Runs ONE full, WARM-STARTED T ladder for the requested (mode, p) configuration and writes
# crash-safe, per-T results to results/data/cluster/warm_sweep.jld2 (keyed by (label, T), so every
# job — the original p=0.1 triple plus any number of psweep jobs — accumulates into the same file
# without clobbering each other).
#
# v2 (July 2026) changes from the first cluster pass, driven by what the cold-started array sweep
# (cluster/{rtm,rdm,cutoff}_array/, one independent SLURM task per T, no warm start at all) showed:
# the dome/wall inflated at T≈6 instead of T≈10, and the naive rank-based selector picked the WRONG
# eigenvalue branch at T=12/13 (|θ_phys| jumped to 1.78/1.98, then snapped back to 1.55 at T=14).
# Two fixes, both load-bearing:
#   1. Selection now uses `pick_phys_continuity` (src/transverse_tools.jl) — searches the WHOLE
#      block for the closest-by-value member, not just the top-2 by modulus.
#   2. The solver runs through `block_transfer_eigs_adaptive` (src), which escalates k=4→6 whenever
#      the block has no tower member besides λ0 (NB5's finding: at larger p a k=4 block can be λ0
#      plus three -λ0-type partners, with λ1 entirely absent).
# And because cold RDM T=9 alone took 20.6h, a sequential warm ladder WILL outlive a single SLURM
# walltime cap — so this version also checkpoints the converged blocks to disk after every T
# (`cluster/checkpoint_<label>.jld2`, gitignored) and resumes truly warm (not just continuity-
# anchored) if the process is killed and the same script is resubmitted.

ENV["GKSwstype"] = "100"   # headless GR backend (src/thesislib.jl unconditionally `using Plots`)

include(joinpath(@__DIR__, "..", "src", "thesislib.jl"))

using LinearAlgebra
try
    @eval using MKL          # Intel MKL BLAS — much faster on the cluster's Xeon nodes
catch err
    @warn "MKL unavailable, falling back to the default BLAS (fine for local dry-runs)" err
end
using Printf

BLAS.set_num_threads(16)   # BLAS-bound workload; SLURM also sets OPENBLAS_NUM_THREADS

# Model / sweep constants — MUST match NB7's master sweep and the desktop χ-scan.
const P_NNN  = 0.1
const LAMBDA = 1.0
const DT     = 0.1
const NBETA  = 4

const CLUSTER_DIR   = joinpath(@__DIR__, "..", "results", "data", "cluster")
const CLUSTER_CACHE = joinpath(CLUSTER_DIR, "warm_sweep.jld2")

# ── nbeta trimming (verbatim from NB7/NB8): the first/last nbeta/2 bonds of a gen_renyi2 profile
#    are imaginary-time cooling, not physical real-time cuts.
function trim_dome(profile, nbeta)
    half = nbeta ÷ 2
    return collect(profile[(half + 1):(end - half)])
end

# ── phase rigidity of a bi-normalized pair (block_transfer_eigs rescales each pair so that
#    ⟨L_j|R_j⟩ = 1, so the rigidity r_j = |⟨L|R⟩|/(‖L‖‖R‖) is just 1/(‖L‖‖R‖)).
function phase_rigidity(Lj::MPS, Rj::MPS)
    return 1.0 / (norm(Lj) * norm(Rj))
end

# ── general χ/ε scan driver, warm-started and checkpointed. `pick_phys_continuity`,
#    `classify_tower`, `tower_gap`, and `block_transfer_eigs_adaptive` all come from
#    src/transverse_tools.jl (included above via thesislib.jl) — the promoted, authoritative
#    versions, not local copies.
function run_wall_scan(; chi::Int, label::String,
        Ts=collect(2.0:1.0:14.0),
        p_nnn::Float64=P_NNN,
        cutoff=1e-14, cutoffs=[fill(1e-12, 40); 1e-14],   # v3: tighter (was 1e-8/1e-10) — stabilises
                                                          # branch selection at large T (fixed the
                                                          # v2 T=10 wrong-branch jump).
        trunc_mode=:rtm, basis=:eig,
        itermax=8000, stuck_after=400,
        k=4, k_retry=4,   # v3: NO k=6 escalation. The escalation was a workaround for classifying the
                          # (physical) π-displaced odd partners as "no tower member"; x1 needs only
                          # λ0 + the smallest partner (both in k=4), and all tower/x1 selection is
                          # post-processing. Keeping k=4 also lightens memory (OOM). k_retry=k=4 ⇒
                          # the adaptive wrapper never bumps k.
        cachefile=joinpath(CLUSTER_DIR, "sweep_$(label).jld2"),   # v3: PER-LABEL cache. Every arm
                          # writes its OWN file — this ends the concurrent-write race on the single
                          # warm_sweep.jld2 that clobbered most of the v2 compute. Merge with
                          # cluster/merge_sweeps.jl for the notebooks.
        checkpointfile=joinpath(@__DIR__, "checkpoint_$(label).jld2"))

    mkpath(dirname(cachefile))
    done = isfile(cachefile) ? load(cachefile, "done") : Dict{Tuple{String,Float64},Any}()

    previous_L = nothing
    previous_R = nothing
    previous_phys = nothing

    for T in Ts
        already_done = haskey(done, (label, T)) && !haskey(done[(label, T)], :error)
        if already_done
            previous_phys = done[(label, T)].theta_phys
            previous_L = nothing
            previous_R = nothing
            continue
        end

        # About to do real work with no warm blocks in memory yet: this is either the very first
        # T of a fresh process, or the first T after skipping past everything already cached on a
        # resubmission. Either way, try to recover the last checkpoint written for this label so
        # the warm start survives a SLURM walltime kill, not just the continuity anchor.
        if previous_L === nothing && isfile(checkpointfile)
            ckpt = load(checkpointfile, "checkpoint")
            if ckpt.label == label && haskey(done, (label, ckpt.T)) &&
                    !haskey(done[(label, ckpt.T)], :error)
                previous_L = ckpt.L
                previous_R = ckpt.R
                previous_phys = done[(label, ckpt.T)].theta_phys
                @info "[$label] warm-resumed from checkpoint at T=$(ckpt.T)"
            end
        end

        try
            elapsed = @elapsed begin
            mpo, scaffold = build_alcaraz_tmpo(T; p=p_nnn, lambda=LAMBDA, dt=DT, nbeta=NBETA, MPO_alg="VD2")
            site_list = siteinds(scaffold)

            if previous_L === nothing
                seedL = nothing
                seedR = nothing
            else
                seedL = MPS[]
                seedR = MPS[]
                for converged_state in previous_L
                    push!(seedL, pad_tmps(converged_state, site_list))
                end
                for converged_state in previous_R
                    push!(seedR, pad_tmps(converged_state, site_list))
                end
            end

            theta, L, R, info = block_transfer_eigs_adaptive(mpo, scaffold;
                k=k, k_retry=k_retry, anchor=previous_phys,
                maxdim=chi, maxdims=collect(2:2:chi),
                cutoff=cutoff, cutoffs=cutoffs,
                itermax=itermax, eps_conv=1e-6, trunc_mode=trunc_mode, basis=basis,
                n_track=2, stuck_after=stuck_after,
                seedL=seedL, seedR=seedR)
            k_actual = length(theta)

            i0 = pick_phys_continuity(theta, previous_phys)

            # v3 validity guard: reject a collapsed/overflowed solve BEFORE caching it as a success.
            # At strong frustration the top of the spectrum becomes (near-)exactly degenerate, and
            # the block solver can return a spurious null Ritz vector (|θ0|≈0) or overflow (|θ0|→∞);
            # v2 cached these as "success", so a resubmit SKIPPED them and warm-start-cascaded the
            # collapse. Throwing here routes it to the catch: cached as :error (hence retried, not
            # skipped) and the warm chain is broken (previous_L/R reset), stopping the cascade.
            lam0_mag = abs(theta[i0])
            if !isfinite(lam0_mag) || lam0_mag < 1e-6 || lam0_mag > 50.0
                error("collapsed/invalid |θ0|=$lam0_mag (near-exact top degeneracy or overflow)")
            end

            dphi, cls = classify_tower(theta; i0=i0)
            gap = tower_gap(theta; i0=i0)
            s2_base = trim_dome(ITransverse.gen_renyi2(L[i0], R[i0]), NBETA)

            rigidity = Float64[]
            for j in 1:k_actual
                push!(rigidity, phase_rigidity(L[j], R[j]))
            end
            end # @elapsed

            done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
                i0=i0, theta_phys=theta[i0],
                dphi=dphi, cls=string.(cls), tower_gap=gap,
                k_used=info[:k_used], escalated=info[:escalated],
                s2_base=s2_base, peak=maximum(real.(s2_base)), rigidity=rigidity,
                reason=string(info[:reason]), niters=info[:niters], elapsed=elapsed)

            previous_phys = theta[i0]
            previous_L = L
            previous_R = R

            # Overwrite the (single, per-label) checkpoint with this T's converged blocks. Only
            # the most recent rung is ever needed to resume, so disk usage stays bounded.
            jldsave(checkpointfile; checkpoint=(label=label, T=T, L=L, R=R))

            rigidity_strings = String[]
            for r in rigidity
                push!(rigidity_strings, @sprintf("%.2g", r))
            end
            @info @sprintf("[%s] T=%.1f  %s@%d  k=%d%s  |θ0|=%.4f  gap=%.3f  peak=%.4f  r=[%s]  %.0fs",
                label, T, info[:reason], info[:niters], info[:k_used],
                info[:escalated] ? "(esc)" : "", abs(theta[i0]), gap,
                maximum(real.(s2_base)), join(rigidity_strings, ","), elapsed)
        catch err
            @warn "[$label] T=$T failed: $err"
            done[(label, T)] = (error=string(err),)
            previous_L = nothing
            previous_R = nothing
        end

        jldsave(cachefile; done=done)
        GC.gc()
    end

    n_ok = 0
    for T in Ts
        if haskey(done, (label, T)) && !haskey(done[(label, T)], :error)
            n_ok += 1
        end
    end
    if n_ok == length(Ts) && isfile(checkpointfile)
        rm(checkpointfile)
        @info "[$label] ladder complete — checkpoint removed"
    end
    println("[$label] cache: $cachefile  ($n_ok/$(length(Ts)) points done)")
    return done
end

# ── entry point: dispatch on the command-line mode ──────────────────────────────────────────────
mode = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia wall_scan_cluster.jl <rdm|psweep> [p] [Tmax]")

# v3: the `cutoff` mode is GONE — its tight cutoff (1e-12/1e-14) is now the default for every arm
# (see run_wall_scan). The p=0.1 arm is just `psweep 0.1 14` (no more standalone `rtm` label).
const RDM_LADDER = collect(2.0:1.0:12.0)   # RDM is the separate, run-LATER diagnostic (NB9: no
                                            # physical gain over RTM for 4-11x cost); own cache via Fix A.

if mode == "rdm"
    run_wall_scan(chi=64, label="rdm64", trunc_mode=:rdm, Ts=RDM_LADDER, p_nnn=P_NNN)
elseif mode == "psweep"
    # A p-sweep arm through the RTM route. Usage: julia wall_scan_cluster.jl psweep <p> <Tmax>
    # Inherits the tight default cutoff and k=4 (no escalation). Writes its own sweep_rtm_p<p>.jld2.
    length(ARGS) >= 3 || error("psweep needs two extra args: julia wall_scan_cluster.jl psweep <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    run_wall_scan(chi=64, label="rtm_p$(p_val)", Ts=collect(2.0:1.0:Tmax), trunc_mode=:rtm, p_nnn=p_val)
else
    error("unknown mode \"$mode\" — expected one of: rdm, psweep")
end
