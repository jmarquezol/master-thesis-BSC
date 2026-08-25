# Walks a ladder of evolution times T for one (mode, p), saving after every rung and resuming from
# its checkpoint when resubmitted.
#
#   julia --project=.. wall_scan_cluster.jl <mode> <p> <Tmax> [dT]
#
#   eigsweep    eigenvalues only
#   entsweep    entropy, one dome per block member
#   towerscan   k=8 block for the boundary dimensions
#   preflight   build one tMPO and exit
#
# Also implemented, not in current use: psweep, ksector, betascan, betawall, rtm, rdm, cutoff.
# Env: WALL_COLUMN (legacy3|bulk5), WALL_DT, WALL_BLAS_THREADS, WALL_RETRIES, WALL_MU_JUMP.

ENV["GKSwstype"] = "100"   # headless GR backend (src/thesislib.jl unconditionally `using Plots`)

include(joinpath(@__DIR__, "..", "src", "thesislib.jl"))

using LinearAlgebra, Printf
# cap the threads so parallel jobs do not oversubscribe the cores
BLAS.set_num_threads(parse(Int, get(ENV, "WALL_BLAS_THREADS", string(Sys.CPU_THREADS))))


const P_NNN  = 0.1
const LAMBDA = 1.0

# NBETA follows DT so that beta0 = nbeta*dt/2 stays 0.2
const DT    = parse(Float64, get(ENV, "WALL_DT", "0.1"))

# only multiples of DT are realisable: T=2.25 at DT=0.1 gets 22 sites and runs at 2.2
check_dT(dT) = isapprox(dT / DT, round(dT / DT); atol=1e-9) ? dT :
    error("dT=$dT is not a multiple of the Trotter step DT=$DT")
const NBETA = round(Int, 0.4 / DT)

# retries for a rung whose leading modulus jumps away from the previous one; 0 disables
const RETRIES = parse(Int, get(ENV, "WALL_RETRIES", "2"))
const MU_JUMP = parse(Float64, get(ENV, "WALL_MU_JUMP", "0.005"))

const COLUMN = Symbol(get(ENV, "WALL_COLUMN", "legacy3"))
COLUMN in (:legacy3, :bulk5) || error("WALL_COLUMN must be legacy3 or bulk5")

const CLUSTER_DIR   = joinpath(@__DIR__, "..", "data", "cluster")
const CLUSTER_CACHE = joinpath(CLUSTER_DIR, "warm_sweep.jld2")

# the first/last nbeta/2 bonds are imaginary-time cooling, not real-time cuts
function trim_dome(profile, nbeta)
    half = nbeta ÷ 2
    return collect(profile[(half + 1):(end - half)])
end

# pairs come back bi-normalized, so the rigidity reduces to 1/(‖L‖‖R‖)
function phase_rigidity(Lj::MPS, Rj::MPS)
    return 1.0 / (norm(Lj) * norm(Rj))
end

# ladder driver; selector and k escalation live in src/transverse_tools.jl
function run_wall_scan(; chi::Int, label::String,
        Ts=collect(2.0:1.0:14.0),
        p_nnn::Float64=P_NNN,
        nbeta::Int=NBETA,
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        trunc_mode=:rtm, basis=:eig,
        eigvals_only::Bool=false,
        itermax=8000, stuck_after=150, eps_conv=1e-6, itermin_floor=20,
        k=4, k_retry=6,
        ksector::Union{Nothing,Tuple{Vector{Float64},Int}}=nothing,
        cachefile=CLUSTER_CACHE,
        checkpointfile=nothing)

    if COLUMN === :bulk5
        label = label * "_bulk"
        cachefile = replace(cachefile, r"\.jld2$" => "_bulk.jld2")
    end
    if DT != 0.1
        # chains of different length must never share a warm start
        label = label * "_dt$(DT)"
        cachefile = replace(cachefile, r"\.jld2$" => "_dt$(DT).jld2")
    end
    checkpointfile === nothing &&
        (checkpointfile = joinpath(@__DIR__, "checkpoints", "checkpoint_$(label).jld2"))

    # Schur basis, no eigenvector work: θ survives where the vectors do not
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

        # nothing warm in memory: recover the checkpoint so a resubmission stays warm
        if previous_L === nothing && isfile(checkpointfile)
            ckpt = load(checkpointfile, "checkpoint")
            # seed forwards only: pad_tmps cannot shrink a vector
            if ckpt.label == label && ckpt.T < T && haskey(done, (label, ckpt.T)) &&
                    !haskey(done[(label, ckpt.T)], :error)
                previous_L = ckpt.L
                previous_R = ckpt.R
                previous_phys = done[(label, ckpt.T)].theta_phys
                @info "[$label] warm-resumed from checkpoint at T=$(ckpt.T)"
            end
        end

        try
            elapsed = @elapsed begin
            mpo, scaffold = build_alcaraz_tmpo(T; p=p_nnn, lambda=LAMBDA, dt=DT, nbeta=nbeta, MPO_alg="VD2", column=COLUMN)
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

            # a wrong fixed point returns a finite theta that no exception catches; a jump in
            # |theta_0| from the previous rung means this rung, so redo it from a fresh seed
            local theta, L, R, info
            for attempt in 1:(1 + RETRIES)
                if ksector === nothing
                    theta, L, R, info = block_transfer_eigs_adaptive(mpo, scaffold;
                        k=k, k_retry=k_retry, anchor=previous_phys,
                        maxdim=chi, maxdims=collect(2:2:chi),
                        cutoff=cutoff, cutoffs=cutoffs,
                        itermax=itermax, eps_conv=eps_conv, itermin=itermin_floor, trunc_mode=trunc_mode, basis=basis,
                        eigvals_only=eigvals_only,
                        n_track=2, stuck_after=stuck_after,
                        seedL=seedL, seedR=seedR)
                else
                    # truncation leaks sectors, so project every iteration, not just the seeds
                    Rd, sector_sign = ksector
                    theta, L, R, info = block_transfer_eigs(mpo, scaffold;
                        k=k, maxdim=chi, maxdims=collect(2:2:chi),
                        cutoff=cutoff, cutoffs=cutoffs,
                        itermax=itermax, eps_conv=eps_conv, itermin=itermin_floor, trunc_mode=trunc_mode, basis=basis,
                        eigvals_only=eigvals_only,
                        n_track=2, stuck_after=stuck_after,
                        seedL=seedL, seedR=seedR,
                        project=psi -> project_ksector(psi, Rd, sector_sign))
                    info = merge(info, Dict(:k_used => k, :escalated => false))
                end
                jumped = previous_phys !== nothing &&
                    abs(abs(theta[pick_phys_robust(theta, previous_phys)[1]]) - abs(previous_phys)) >
                        MU_JUMP * abs(previous_phys)
                jumped || break
                if attempt > RETRIES
                    @warn @sprintf("[%s] T=%.2f still jumping after %d attempts, keeping it flagged",
                                   label, T, attempt)
                    break
                end
                @warn @sprintf("[%s] T=%.2f attempt %d: |theta0| jumped from %.6f, retrying with a fresh seed",
                               label, T, attempt, abs(previous_phys))
                seedL = nothing          # an independent draw, not a warm start from the same basin
                seedR = nothing
            end
            k_actual = length(theta)

            i0, recovered = pick_phys_robust(theta, previous_phys)
            dphi, cls = classify_tower(theta; i0=i0)
            gap = tower_gap(theta; i0=i0)

            # entropy and rigidity need the bi-normalized pairs, which :schur does not return
            if eigvals_only
                s2_base = ComplexF64[]
                s2_all = Vector{ComplexF64}[]
                rigidity = Float64[]
            else
                # one dome per member: i0 is reliable across the ladder, not rung by rung
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
            kcharge = ksector === nothing ? Float64[] :
                [real(overlap_noconj(L[j], apply_ksign(R[j], ksector[1])) /
                      overlap_noconj(L[j], R[j])) for j in 1:k_actual]
            done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
                i0=i0, theta_phys=theta[i0],
                dphi=dphi, cls=string.(cls), tower_gap=gap,
                k_used=info[:k_used], escalated=info[:escalated],
                s2_base=s2_base, s2_all=s2_all, peak=peak, rigidity=rigidity, kcharge=kcharge,
                reason=string(info[:reason]), niters=info[:niters], elapsed=elapsed)

            recovered && (previous_phys = theta[i0])
            previous_L = L
            previous_R = R

            # only the most recent rung is needed to resume
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

        # a broken environment fails on every rung; stop rather than grind
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
    # keep the checkpoint after the ladder finishes so a larger Tmax resumes warm rather than cold
    if n_ok == length(Ts) && isfile(checkpointfile)
        ckpt_T = load(checkpointfile, "checkpoint").T
        @info "[$label] ladder complete — checkpoint kept at T=$(ckpt_T) for a later extension"
    end
    println("[$label] cache: $cachefile  ($n_ok/$(length(Ts)) points done)")
    return done
end

# entry point
mode = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia wall_scan_cluster.jl <preflight|rtm|rdm|cutoff|psweep|eigsweep|betascan|betawall> [p] [nbeta] [Tmax]")

const FULL_LADDER    = collect(2.0:1.0:14.0)
const RTM_FULL_LADDER = collect(2.0:1.0:20.0)  # rtm alone now matches the psweep arms' T=20 reach
const RDM_LADDER     = collect(2.0:1.0:12.0)   # cold T=9 alone took 20.6h; two points past the warm
                                                # wall suffice — extend Ts + resubmit if ever needed.

if mode == "preflight"
    # same tMPO call the ladder makes on every rung
    mpo, scaffold = build_alcaraz_tmpo(2.0; p=0.1, lambda=LAMBDA, dt=DT, nbeta=NBETA, MPO_alg="VD2", column=COLUMN)
    println("preflight OK — column=$(COLUMN), tMPO built, $(length(scaffold)) sites, " *
            "temporal site dim $(dim(siteind(mpo, 2))), maxlinkdim $(maxlinkdim(mpo))")
elseif mode == "rtm"
    run_wall_scan(chi=64, label="rtm64_full", Ts=RTM_FULL_LADDER, p_nnn=P_NNN)
elseif mode == "rdm"
    run_wall_scan(chi=64, label="rdm_p0.1", trunc_mode=:rdm, Ts=RDM_LADDER, p_nnn=P_NNN,
        cachefile=joinpath(CLUSTER_DIR, "sweep_rdm_p0.1.jld2"))
elseif mode == "cutoff"
    run_wall_scan(chi=64, label="cut_tight", cutoffs=[fill(1e-10, 40); 1e-12], Ts=FULL_LADDER, p_nnn=P_NNN)
elseif mode == "psweep"
    # usage: psweep <p> <Tmax> [dT]
    length(ARGS) >= 3 || error("psweep needs two extra args: julia wall_scan_cluster.jl psweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    dT    = check_dT(length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0)
    run_wall_scan(chi=64, label="rtm_p$(p_val)", Ts=collect(2.0:dT:Tmax), trunc_mode=:rtm, p_nnn=p_val,
        cachefile=joinpath(CLUSTER_DIR, "sweep_rtm_p$(p_val).jld2"))
elseif mode == "fork"
    # fork <src_label> <dst_label>: seed a second interleaved chain from an existing arm's
    # frontier, so it warm-starts from the same checkpoint without sharing cache or checkpoint.
    # Idempotent: does nothing if the destination cache already exists.
    length(ARGS) >= 3 || error("fork needs: julia wall_scan_cluster.jl fork <src_label> <dst_label>")
    sfx(l) = begin
        l2 = COLUMN === :bulk5 ? l * "_bulk" : l
        DT != 0.1 ? l2 * "_dt$(DT)" : l2
    end
    slab, dlab = sfx(ARGS[2]), sfx(ARGS[3])
    scache = joinpath(CLUSTER_DIR, "sweep_$(slab).jld2")
    dcache = joinpath(CLUSTER_DIR, "sweep_$(dlab).jld2")
    sckpt  = joinpath(@__DIR__, "checkpoints", "checkpoint_$(slab).jld2")
    dckpt  = joinpath(@__DIR__, "checkpoints", "checkpoint_$(dlab).jld2")
    if isfile(dcache)
        @info "fork: $dcache already exists, nothing to do"
    else
        isfile(scache) || error("fork: source cache $scache not found")
        isfile(sckpt) || error("fork: source checkpoint $sckpt not found")
        src_done = load(scache, "done")
        ck = load(sckpt, "checkpoint")
        ck.label == slab || error("fork: checkpoint label $(ck.label) does not match $slab")
        haskey(src_done, (slab, ck.T)) || error("fork: checkpoint rung T=$(ck.T) missing from source cache")
        rec = src_done[(slab, ck.T)]
        haskey(rec, :error) && error("fork: source frontier rung T=$(ck.T) is an error record")
        dst_done = Dict{Tuple{String,Float64},Any}((dlab, ck.T) => (; rec..., label=dlab))
        jldsave(dcache; done=dst_done)
        jldsave(dckpt; checkpoint=(label=dlab, T=ck.T, L=ck.L, R=ck.R))
        @info "fork: seeded $dlab from $slab at frontier T=$(ck.T)"
    end
elseif mode == "entsweep"
    # usage: entsweep <p> <Tmax> [dT] — stores every dome, so the branch choice stays revisable
    length(ARGS) >= 3 || error("entsweep needs two extra args: julia wall_scan_cluster.jl entsweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    dT    = check_dT(length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0)
    # optional interleaving: [T0] starts the ladder above the frontier, [label] names the chain
    T0    = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 2.0
    lbl   = length(ARGS) >= 6 ? ARGS[6] : "ent_p$(p_val)"
    run_wall_scan(chi=64, label=lbl, Ts=collect(T0:dT:Tmax), trunc_mode=:rtm, p_nnn=p_val,
        cachefile=joinpath(CLUSTER_DIR, "sweep_$(lbl).jld2"))
elseif mode == "towerscan"
    # usage: towerscan <p> <Tmax> [dT] — k=8, for the boundary dimensions
    length(ARGS) >= 3 || error("towerscan needs two extra args: julia wall_scan_cluster.jl towerscan <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])

    dT    = check_dT(length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0)
    # chi=48 suffices here; itermin=80 stops the loose tolerance accepting a rung while the
    # truncation schedule is still ramping
    run_wall_scan(chi=48, label="tower_p$(p_val)", Ts=collect(2.0:dT:Tmax), k=8, k_retry=10,
        trunc_mode=:rtm, p_nnn=p_val, eigvals_only=true, eps_conv=1e-5, itermin_floor=80,
        cachefile=joinpath(CLUSTER_DIR, "sweep_tower_p$(p_val).jld2"))
elseif mode == "eigsweep"
    # usage: eigsweep <p> <Tmax> [dT] — no eigenvector work, so it reaches larger T
    length(ARGS) >= 3 || error("eigsweep needs two extra args: julia wall_scan_cluster.jl eigsweep <p> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    # a unit ladder advances the phase by more than pi at large p, where the branch can no longer
    # be identified; a finer dT keeps the advance resolvable
    dT    = check_dT(length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 1.0)
    # its own label, so it never shares a cache or checkpoint with the unit ladder
    T0  = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 2.0
    lbl = length(ARGS) >= 6 ? ARGS[6] :
          (dT == 1.0 ? "rtm_eigs_p$(p_val)" : "rtm_eigs_p$(p_val)_fine")
    run_wall_scan(chi=64, label=lbl, Ts=collect(T0:dT:Tmax),
        trunc_mode=:rtm, p_nnn=p_val, eigvals_only=true,
        cachefile=joinpath(CLUSTER_DIR, "sweep_$(lbl).jld2"))
elseif mode == "ksector"
    # usage: ksector <p> <plus|minus> <Tmax> [dT] — one K sector, so the branch is the sector label
    length(ARGS) >= 4 || error("ksector needs <p> <plus|minus> <Tmax> [dT]")
    p_val = parse(Float64, ARGS[2])
    sector_sign = ARGS[3] in ("plus", "+1", "+") ? 1 :
                  ARGS[3] in ("minus", "-1", "-") ? -1 : error("sector must be plus or minus")
    Tmax = parse(Float64, ARGS[4])
    dT = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1.0
    Rd = ksector_signs(p_val; dt=DT, nbeta=NBETA, column=COLUMN)
    lbl = "ksec_p$(p_val)_" * (sector_sign == 1 ? "plus" : "minus")
    run_wall_scan(chi=64, label=lbl, Ts=collect(2.0:dT:Tmax), p_nnn=p_val,
        eigvals_only=true, k=2, ksector=(Rd, sector_sign),
        cachefile=joinpath(CLUSTER_DIR, "sweep_$(lbl).jld2"))

elseif mode == "betascan"
    # usage: betascan <p> <Tmax> — the same run repeated over nbeta, to see how the read depends
    # on the regulator: too small dirties the boundary, too large inflates the finite-time term
    length(ARGS) >= 3 || error("betascan needs two extra args: julia wall_scan_cluster.jl betascan <p> <Tmax>")
    p_val = parse(Float64, ARGS[2])
    Tmax  = parse(Float64, ARGS[3])
    betacache = joinpath(CLUSTER_DIR, "sweep_beta_p$(p_val).jld2")
    for nb in (2, 4, 6, 8, 10, 12, 14, 16)   # β0 = 0.1 … 0.8; cached rungs are skipped on rerun
        run_wall_scan(chi=64, label="beta_p$(p_val)_nb$(nb)", Ts=collect(2.0:1.0:Tmax),
            trunc_mode=:rtm, p_nnn=p_val, nbeta=nb, cachefile=betacache)
    end
elseif mode == "betawall"
    # usage: betawall <p> <nbeta> <Tmax> — one regulator value on a long ladder. Full eigenvector
    # run on purpose: eigvals_only would answer a different question.
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
