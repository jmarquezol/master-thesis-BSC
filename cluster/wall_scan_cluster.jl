# Batch driver for the wall-scan sweeps. One process walks a ladder of evolution times T for a given
# (mode, p), warm-starting each T from the last and checkpointing as it goes, so a walltime-killed
# job resumes where it stopped when resubmitted.
#
#   julia --project=. cluster/wall_scan_cluster.jl <mode> [p] [Tmax]
#
# modes:
#   rtm / rdm / cutoff        the original p=0.1 truncation comparison (no extra args)
#   psweep    <p> <Tmax>      full run (with entropy), RTM, at coupling p
#   eigsweep  <p> <Tmax>      eigenvalues only (no entropy) — goes past the wall
#   betascan  <p> <Tmax>      full run, repeated over nbeta=2..16 — the β0 regulator scan
#
# Branch selection and the k=4→6 escalation live in src/transverse_tools.jl (pick_phys_continuity,
# block_transfer_eigs_adaptive); both matter once the spectrum gets near-degenerate at larger p.

ENV["GKSwstype"] = "100"   # headless GR backend (src/thesislib.jl unconditionally `using Plots`)

include(joinpath(@__DIR__, "..", "src", "thesislib.jl"))

using LinearAlgebra, Printf
BLAS.set_num_threads(Sys.CPU_THREADS)   # BLAS-bound workload; SLURM also sets OPENBLAS_NUM_THREADS

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
        nbeta::Int=NBETA,
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        trunc_mode=:rtm, basis=:eig,
        eigvals_only::Bool=false,
        itermax=8000, stuck_after=400,
        k=4, k_retry=6,
        cachefile=CLUSTER_CACHE,
        checkpointfile=joinpath(@__DIR__, "checkpoint_$(label).jld2"))

    # Eigenvalue-only mode (basis=:schur, no eigenvector post-processing): the leading spectrum
    # (dual-unitarity circle, Eq.(3) c, Eq.(4) x1, tower gaps) is a Rayleigh-quotient-like quantity
    # that survives the entanglement-barrier wall, whereas the entropy needs the eigenVECTORS and
    # does not. Skipping the vector work (gen_renyi2 + phase rigidity, both of which require the
    # bi-normalized pairs that :schur/eigvals_only do NOT return) makes each T cheaper and well-
    # conditioned, so the spectral ladder reaches larger T than the full-eigenvector runs do.
    if eigvals_only
        basis = :schur
    end

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
            mpo, scaffold = build_alcaraz_tmpo(T; p=p_nnn, lambda=LAMBDA, dt=DT, nbeta=nbeta, MPO_alg="VD2")
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
                eigvals_only=eigvals_only,
                n_track=2, stuck_after=stuck_after,
                seedL=seedL, seedR=seedR)
            k_actual = length(theta)

            i0 = pick_phys_continuity(theta, previous_phys)
            dphi, cls = classify_tower(theta; i0=i0)
            gap = tower_gap(theta; i0=i0)

            # The entropy (gen_renyi2) and phase rigidity both need the bi-normalized eigenVECTOR
            # pairs, which eigvals_only=:schur does NOT return — so skip them entirely in that mode
            # and store empty placeholders. Everything above (theta, tower classification, gap) needs
            # only the eigenvalues and is valid.
            if eigvals_only
                s2_base = Float64[]
                rigidity = Float64[]
            else
                s2_base = trim_dome(ITransverse.gen_renyi2(L[i0], R[i0]), nbeta)
                rigidity = Float64[]
                for j in 1:k_actual
                    push!(rigidity, phase_rigidity(L[j], R[j]))
                end
            end
            end # @elapsed

            peak = isempty(s2_base) ? NaN : maximum(real.(s2_base))   # no entropy in eigvals-only mode
            done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
                i0=i0, theta_phys=theta[i0],
                dphi=dphi, cls=string.(cls), tower_gap=gap,
                k_used=info[:k_used], escalated=info[:escalated],
                s2_base=s2_base, peak=peak, rigidity=rigidity,
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
                peak, join(rigidity_strings, ","), elapsed)
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
mode = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia wall_scan_cluster.jl <rtm|rdm|cutoff|psweep|eigsweep|betascan> [p] [Tmax]")

const FULL_LADDER    = collect(2.0:1.0:14.0)
const RTM_FULL_LADDER = collect(2.0:1.0:20.0)  # rtm alone now matches the psweep arms' T=20 reach
const RDM_LADDER     = collect(2.0:1.0:12.0)   # cold T=9 alone took 20.6h; two points past the warm
                                                # wall suffice — extend Ts + resubmit if ever needed.

if mode == "rtm"
    run_wall_scan(chi=64, label="rtm64_full", Ts=RTM_FULL_LADDER, p_nnn=P_NNN)
elseif mode == "rdm"
    run_wall_scan(chi=64, label="rdm_p0.1", trunc_mode=:rdm, Ts=RDM_LADDER, p_nnn=P_NNN,
        cachefile=joinpath(CLUSTER_DIR, "sweep_rdm_p0.1.jld2"))
elseif mode == "cutoff"
    run_wall_scan(chi=64, label="cut_tight", cutoffs=[fill(1e-10, 40); 1e-12], Ts=FULL_LADDER, p_nnn=P_NNN)
elseif mode == "psweep"
    # A p-sweep job, always through the RTM route (NB9's cost comparison: RDM buys nothing
    # physical for ~4-11x the cost, so it is not worth extending to the p-sweep at all).
    # Usage: julia wall_scan_cluster.jl psweep <p> <Tmax>
    length(ARGS) >= 3 || error("psweep needs two extra args: julia wall_scan_cluster.jl psweep <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    run_wall_scan(chi=64, label="rtm_p$(p_val)", Ts=collect(2.0:1.0:Tmax), trunc_mode=:rtm, p_nnn=p_val)
elseif mode == "eigsweep"
    # EIGENVALUE-ONLY p-sweep arm: identical configuration to `psweep` (RTM route, χ=64, same strict
    # cutoff schedule, warm-started + checkpointed) EXCEPT it runs the block solver in Schur /
    # eigvals-only mode — it computes only the leading spectrum and skips the eigenVECTOR work
    # (entropy + phase rigidity). This is the right tool for the spectral route (dual unitarity, the
    # Eq.(3) central charge, the Eq.(4) boundary exponent, and the tower gaps), all of which are
    # Rayleigh-quotient-like and survive the entanglement-barrier wall — so this arm reaches larger T
    # than the full-eigenvector `psweep` runs, which stall once the eigenvector conditioning collapses.
    # Written to its OWN per-p cache (sweep_rtm_eigs_p<p>.jld2) so it never clobbers the full runs.
    # Usage: julia wall_scan_cluster.jl eigsweep <p> <Tmax>
    length(ARGS) >= 3 || error("eigsweep needs two extra args: julia wall_scan_cluster.jl eigsweep <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    run_wall_scan(chi=64, label="rtm_eigs_p$(p_val)", Ts=collect(2.0:1.0:Tmax),
        trunc_mode=:rtm, p_nnn=p_val, eigvals_only=true,
        cachefile=joinpath(CLUSTER_DIR, "sweep_rtm_eigs_p$(p_val).jld2"))
elseif mode == "betascan"
    # Regulator scan: same full RTM run, repeated for a few values of the imaginary-time cooling
    # nbeta (β0 = nbeta·dt/2 = 0.1, 0.2, 0.3, 0.4 at dt=0.1). Lets us see how the CFT read depends on
    # β0 — too small dirties the boundary, too large inflates the finite-time correction (ε2=2β0/T).
    # Modest T, so it finishes. Usage: julia wall_scan_cluster.jl betascan <p> <Tmax>
    length(ARGS) >= 3 || error("betascan needs two extra args: julia wall_scan_cluster.jl betascan <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    betacache = joinpath(CLUSTER_DIR, "sweep_beta_p$(p_val).jld2")
    for nb in (2, 4, 6, 8, 10, 12, 14, 16)   # β0 = 0.1 … 0.8; cached rungs are skipped on rerun
        run_wall_scan(chi=64, label="beta_p$(p_val)_nb$(nb)", Ts=collect(2.0:1.0:Tmax),
            trunc_mode=:rtm, p_nnn=p_val, nbeta=nb, cachefile=betacache)
    end
else
    error("unknown mode \"$mode\" — expected one of: rtm, rdm, cutoff, psweep, eigsweep, betascan")
end
