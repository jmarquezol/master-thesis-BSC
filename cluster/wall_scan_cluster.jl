# Batch driver for the wall-scan sweeps. One process walks a ladder of evolution times T for a given
# (mode, p), warm-starting each T from the last and checkpointing as it goes, so a walltime-killed
# job resumes where it stopped when resubmitted.
#
#   julia --project=. cluster/wall_scan_cluster.jl <mode> [p] [Tmax]
#
# modes:
#   preflight                 build one small tMPO and exit — cheap check that the env still works
#   rtm / rdm / cutoff        the original p=0.1 truncation comparison (no extra args)
#   psweep    <p> <Tmax> [dT]  full run (with entropy), RTM, at coupling p; dT defaults to 1
#   eigsweep  <p> <Tmax> [dT]  eigenvalues only (no entropy) — goes past the wall; dT defaults to 1
#   betascan  <p> <Tmax>      full run, repeated over nbeta=2..16 — the β0 regulator scan
#   betawall  <p> <nbeta> <Tmax>   full run at ONE β0 on a long ladder — does β0 move the wall?
#
# Branch selection and the k escalation live in src/transverse_tools.jl; both matter once the
# spectrum gets near-degenerate at larger p.

ENV["GKSwstype"] = "100"   # headless GR backend (src/thesislib.jl unconditionally `using Plots`)

include(joinpath(@__DIR__, "..", "src", "thesislib.jl"))

using LinearAlgebra, Printf
BLAS.set_num_threads(Sys.CPU_THREADS)   # BLAS-bound workload; SLURM also sets OPENBLAS_NUM_THREADS

# Model / sweep constants — these match NB7's master sweep and the desktop χ-scan.
const P_NNN  = 0.1
const LAMBDA = 1.0
const DT     = 0.1
const NBETA  = 4

const CLUSTER_DIR   = joinpath(@__DIR__, "..", "results", "data", "cluster")
const CLUSTER_CACHE = joinpath(CLUSTER_DIR, "warm_sweep.jld2")

# ── the first/last nbeta/2 bonds of a gen_renyi2 profile are imaginary-time cooling, not
#    physical real-time cuts, so trim them before reading a dome.
function trim_dome(profile, nbeta)
    half = nbeta ÷ 2
    return collect(profile[(half + 1):(end - half)])
end

# ── phase rigidity of a bi-normalized pair (block_transfer_eigs rescales each pair so that
#    ⟨L_j|R_j⟩ = 1, so the rigidity r_j = |⟨L|R⟩|/(‖L‖‖R‖) is just 1/(‖L‖‖R‖)).
function phase_rigidity(Lj::MPS, Rj::MPS)
    return 1.0 / (norm(Lj) * norm(Rj))
end

# ── general χ/ε scan driver, warm-started and checkpointed. The selector and escalation helpers
#    come from src/transverse_tools.jl.
function run_wall_scan(; chi::Int, label::String,
        Ts=collect(2.0:1.0:14.0),
        p_nnn::Float64=P_NNN,
        nbeta::Int=NBETA,
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        trunc_mode=:rtm, basis=:eig,
        eigvals_only::Bool=false,
        itermax=8000, stuck_after=150,
        k=4, k_retry=6,
        cachefile=CLUSTER_CACHE,
        checkpointfile=joinpath(@__DIR__, "checkpoints", "checkpoint_$(label).jld2"))

    # Spectrum-only mode: de-mix on the Schur basis and skip the eigenvector work. θ is a Rayleigh
    # quotient so it survives the wall even where the vectors do not.
    if eigvals_only
        basis = :schur
    end

    mkpath(dirname(cachefile))
    mkpath(dirname(checkpointfile))
    done = isfile(cachefile) ? load(cachefile, "done") : Dict{Tuple{String,Float64},Any}()

    previous_L = nothing
    previous_R = nothing
    previous_phys = nothing
    consecutive_failures = 0

    for T in Ts
        already_done = haskey(done, (label, T)) && !haskey(done[(label, T)], :error)
        if already_done
            previous_phys = done[(label, T)].theta_phys
            previous_L = nothing
            previous_R = nothing
            continue
        end

        # No warm blocks in memory: either a fresh process, or the first T after skipping the
        # cached ones on a resubmission. Recover the last checkpoint so the warm start survives.
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

            i0, recovered = pick_phys_robust(theta, previous_phys)
            dphi, cls = classify_tower(theta; i0=i0)
            gap = tower_gap(theta; i0=i0)

            # Entropy and rigidity both need the bi-normalized pairs, which :schur does not return.
            if eigvals_only
                s2_base = ComplexF64[]
                s2_all = Vector{ComplexF64}[]
                rigidity = Float64[]
            else
                # one dome per member: i0 is only reliable across the whole ladder, not rung by rung
                s2_all = Vector{ComplexF64}[]
                rigidity = Float64[]
                for j in 1:k_actual
                    push!(s2_all, trim_dome(ITransverse.gen_renyi2(L[j], R[j]), nbeta))
                    push!(rigidity, phase_rigidity(L[j], R[j]))
                end
                s2_base = s2_all[i0]
            end
            end # @elapsed

            peak = isempty(s2_base) ? NaN : maximum(real.(s2_base))   # no entropy in eigvals-only mode
            done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
                i0=i0, theta_phys=theta[i0],
                dphi=dphi, cls=string.(cls), tower_gap=gap,
                k_used=info[:k_used], escalated=info[:escalated],
                s2_base=s2_base, s2_all=s2_all, peak=peak, rigidity=rigidity,
                reason=string(info[:reason]), niters=info[:niters], elapsed=elapsed)

            recovered && (previous_phys = theta[i0])
            previous_L = L
            previous_R = R

            # Only the most recent rung is needed to resume, so one file per label is enough.
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
            consecutive_failures += 1
        else
            consecutive_failures = 0
        end

        jldsave(cachefile; done=done)
        GC.gc()

        # A broken environment fails on every rung; stop rather than error the whole ladder.
        if consecutive_failures >= 2
            @error "[$label] aborting at T=$T after two consecutive failures — check the environment"
            break
        end
    end

    n_ok = 0
    for T in Ts
        if haskey(done, (label, T)) && !haskey(done[(label, T)], :error)
            n_ok += 1
        end
    end
    # Keep the checkpoint even when the ladder finishes, so resubmitting with a larger Tmax later
    # resumes warm instead of cold-starting the first new T (a cold restart is what broke the dome
    # at T≈6 in the first array sweep). One file per label, overwritten each rung. A stale one is
    # harmless: the resume above also checks the label and that the rung is in this cache.
    if n_ok == length(Ts) && isfile(checkpointfile)
        ckpt_T = load(checkpointfile, "checkpoint").T
        @info "[$label] ladder complete — checkpoint kept at T=$(ckpt_T) for a later extension"
    end
    println("[$label] cache: $cachefile  ($n_ok/$(length(Ts)) points done)")
    return done
end

# ── entry point: dispatch on the command-line mode ──────────────────────────────────────────────
mode = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia wall_scan_cluster.jl <preflight|rtm|rdm|cutoff|psweep|eigsweep|betascan|betawall> [p] [nbeta] [Tmax]")

const FULL_LADDER    = collect(2.0:1.0:14.0)
const RTM_FULL_LADDER = collect(2.0:1.0:20.0)  # rtm alone now matches the psweep arms' T=20 reach
const RDM_LADDER     = collect(2.0:1.0:12.0)   # cold T=9 alone took 20.6h; two points past the warm
                                                # wall suffice — extend Ts + resubmit if ever needed.

if mode == "preflight"
    # Cheap environment check: same tMPO call the ladder makes on every rung.
    mpo, scaffold = build_alcaraz_tmpo(2.0; p=0.1, lambda=LAMBDA, dt=DT, nbeta=NBETA, MPO_alg="VD2")
    println("preflight OK — tMPO built, $(length(scaffold)) sites, maxlinkdim $(maxlinkdim(mpo))")
elseif mode == "rtm"
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
    length(ARGS) >= 3 || error("psweep needs two extra args: julia wall_scan_cluster.jl psweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    dT    = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    run_wall_scan(chi=64, label="rtm_p$(p_val)", Ts=collect(2.0:dT:Tmax), trunc_mode=:rtm, p_nnn=p_val,
        cachefile=joinpath(CLUSTER_DIR, "sweep_rtm_p$(p_val).jld2"))
elseif mode == "entsweep"
    # Entropy arm from T=2, one dome per block member. Own cache: the old psweep runs stored only
    # the selected dome, so their branch choice cannot be revisited.
    # Usage: julia wall_scan_cluster.jl entsweep <p> <Tmax> [dT]
    length(ARGS) >= 3 || error("entsweep needs two extra args: julia wall_scan_cluster.jl entsweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    dT    = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    run_wall_scan(chi=64, label="ent_p$(p_val)", Ts=collect(2.0:dT:Tmax), trunc_mode=:rtm, p_nnn=p_val,
        cachefile=joinpath(CLUSTER_DIR, "sweep_ent_p$(p_val).jld2"))
elseif mode == "towerscan"
    # Deep k=8 block at small T: the tower figure needs more members than the k=4 arms carry.
    # Usage: julia wall_scan_cluster.jl towerscan <p> <Tmax>
    length(ARGS) >= 3 || error("towerscan needs two extra args: julia wall_scan_cluster.jl towerscan <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    run_wall_scan(chi=64, label="tower_p$(p_val)", Ts=collect(2.0:1.0:Tmax), k=8, k_retry=10,
        trunc_mode=:rtm, p_nnn=p_val, eigvals_only=true,
        cachefile=joinpath(CLUSTER_DIR, "sweep_tower_p$(p_val).jld2"))
elseif mode == "eigsweep"
    # Eigenvalues only: same configuration as `psweep` but in Schur/eigvals-only mode, skipping the
    # eigenvector work (entropy, rigidity). The spectrum is Rayleigh-quotient-like and survives the
    # wall, so this arm reaches larger T than the full runs — it is the right tool for dual
    # unitarity, the Eq.(3) central charge, Eq.(4), and the tower gaps. Own cache per p.
    # Usage: julia wall_scan_cluster.jl eigsweep <p> <Tmax>
    length(ARGS) >= 3 || error("eigsweep needs two extra args: julia wall_scan_cluster.jl eigsweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    # dT < 1 fills half-integer rungs into the SAME cache. The branch tracker follows the physical
    # eigenvalue by predicting its phase, and the phase advance per rung grows with the velocity, so
    # at larger p a unit ladder advances by more than pi and the branch can no longer be identified.
    dT    = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0
    run_wall_scan(chi=64, label="rtm_eigs_p$(p_val)", Ts=collect(2.0:dT:Tmax),
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
elseif mode == "betawall"
    # Does the regulator move the wall? The β0 scan showed the modulus gaps grow linearly with β0
    # (Eq. 14 of the 2026 paper, confirmed at R²≈0.99), and it is those gaps closing that ends the
    # eigenvector route — so a larger β0 might buy reach. One β0 per job, long ladder (the opposite
    # shape to betascan). Expect a modest shift at best: the enhancement carries a 1/T².
    #
    # Full eigenvector run on purpose — the wall is the entropy dome inflating, so eigvals_only
    # would answer a different question. Own cache, and the ladder starts at T=2 rather than
    # resuming betascan at T=7, because that ladder is finished and its first new rung would start
    # cold — which is what corrupted the dome in the first array sweep.
    # Usage: julia wall_scan_cluster.jl betawall <p> <nbeta> <Tmax>
    length(ARGS) >= 4 || error("betawall needs three extra args: julia wall_scan_cluster.jl betawall <p> <nbeta> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    nb    = parse(Int, ARGS[3])
    Tmax  = parse(Float64, ARGS[4])
    run_wall_scan(chi=64, label="betawall_p$(p_val)_nb$(nb)", Ts=collect(2.0:1.0:Tmax),
        trunc_mode=:rtm, p_nnn=p_val, nbeta=nb,
        cachefile=joinpath(CLUSTER_DIR, "sweep_betawall_p$(p_val).jld2"))
else
    error("unknown mode \"$mode\" — expected one of: preflight, rtm, rdm, cutoff, psweep, eigsweep, betascan, betawall")
end
