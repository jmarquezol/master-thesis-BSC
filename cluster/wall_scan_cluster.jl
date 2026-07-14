# cluster/wall_scan_cluster.jl — headless batch driver for the BSC master sweep.
#
# Usage:  julia --project=. cluster/wall_scan_cluster.jl <mode>
#   <mode> ∈ {"rtm", "rdm", "cutoff"}
#
# Runs ONE full T=2..14 ladder of `run_wall_scan` for the requested truncation configuration and
# writes crash-safe, per-T results to results/data/nb13_wallscan_cluster.jld2 (keyed by
# (label, T), so all three modes accumulate into the same file without clobbering each other or
# the pre-existing local nb13_wallscan.jld2 from the desktop χ-scan).
#
# This script is a self-contained copy of the logic in NBs/8_eigvec_robustness.ipynb's
# NB13-HELPERS / NB13-WALLSCAN cells (only the four functions run_wall_scan actually needs — the
# notebook's other E3/E4/E6/E7 helpers are not needed here). Kept as a copy rather than promoted
# into src/ to avoid touching the already-validated notebook; promote if this needs a third home.

# Headless plotting backend: src/thesislib.jl unconditionally `using Plots`, and the GR backend
# tries to open a display unless told not to. Must be set BEFORE the include.
ENV["GKSwstype"] = "100"

include(joinpath(@__DIR__, "..", "src", "thesislib.jl"))

using LinearAlgebra, Printf
# Our workload is BLAS-bound (ITensors tensor contractions), not Threads.@threads-parallel.
# The SLURM scripts also set OPENBLAS_NUM_THREADS; this call makes the script correct even when
# run outside SLURM (e.g. the local smoke test in the README).
BLAS.set_num_threads(Sys.CPU_THREADS)

# Model / sweep constants — MUST match NB7's master sweep and the desktop χ-scan so every series
# in the eventual comparison shares the same physical setup.
const P_NNN  = 0.1     # the frustrated Alcaraz point (the wall at T≈10 for the χ=64 baseline)
const LAMBDA = 1.0     # quench TO criticality
const DT     = 0.1
const NBETA  = 4       # β0 = 0.2 conformal cooling (CLAUDE.md §13)

const CLUSTER_CACHE = joinpath(@__DIR__, "..", "results", "data", "nb13_wallscan_cluster.jld2")

# ── continuity selector (verbatim from NB7/NB8): which of the two leading Ritz values is "the
#    physical λ0"? Follow the smooth curve: keep whichever is closest to the previous T's value.
function pick_phys(theta_values, previous_physical_value)
    if previous_physical_value === nothing
        return (1, 2)
    end

    distance_to_first  = abs(theta_values[1] - previous_physical_value)
    distance_to_second = abs(theta_values[2] - previous_physical_value)

    if distance_to_first <= distance_to_second
        return (1, 2)
    else
        return (2, 1)
    end
end

# ── nbeta trimming (verbatim from NB7/NB8): the first/last nbeta/2 bonds of a gen_renyi2 profile
#    are imaginary-time cooling, not physical real-time cuts.
function trim_dome(profile, nbeta)
    half = nbeta ÷ 2
    return collect(profile[(half + 1):(end - half)])
end

# ── phase rigidity of a bi-normalized pair. block_transfer_eigs rescales each pair so that
#    ⟨L_j|R_j⟩ = 1, so the rigidity r_j = |⟨L|R⟩|/(‖L‖‖R‖) is just 1/(‖L‖‖R‖).
#    r_j → 0 signals approach to an exceptional point (eigenvector coalescence).
function phase_rigidity(Lj::MPS, Rj::MPS)
    return 1.0 / (norm(Lj) * norm(Rj))
end

# ── general χ/ε scan driver (verbatim from NB8's NB13-WALLSCAN cell). Cache keyed by (label, T);
#    each label gets its own independent warm-start chain over its own Ts.
function run_wall_scan(; chi::Int, label::String,
        Ts=collect(2.0:1.0:14.0),
        cutoff=1e-12, cutoffs=[fill(1e-8, 40); 1e-10],
        trunc_mode=:rtm, basis=:eig,
        itermax=8000, stuck_after=400,
        cachefile=CLUSTER_CACHE)

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

        try
            elapsed = @elapsed begin
            mpo, scaffold = build_alcaraz_tmpo(T; p=P_NNN, lambda=LAMBDA, dt=DT, nbeta=NBETA, MPO_alg="VD2")
            site_list = siteinds(scaffold)

            # warm start from this label's own previous T (independent chain per label)
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

            theta, L, R, info = block_transfer_eigs(mpo, scaffold;
                k=4, maxdim=chi, maxdims=collect(2:2:chi),
                cutoff=cutoff, cutoffs=cutoffs,
                itermax=itermax, eps_conv=1e-6, trunc_mode=trunc_mode, basis=basis,
                n_track=2, stuck_after=stuck_after,
                seedL=seedL, seedR=seedR)
            k = length(theta)

            i0, ip = pick_phys(theta, previous_phys)
            s2_base = trim_dome(ITransverse.gen_renyi2(L[i0], R[i0]), NBETA)

            rigidity = Float64[]
            for j in 1:k
                push!(rigidity, phase_rigidity(L[j], R[j]))
            end
            end # @elapsed

            done[(label, T)] = (label=label, T=T, chi=chi, theta=collect(theta),
                i0=i0, ip=ip, theta_phys=theta[i0], theta_partner=theta[ip],
                s2_base=s2_base, peak=maximum(real.(s2_base)), rigidity=rigidity,
                reason=string(info[:reason]), niters=info[:niters], elapsed=elapsed)

            previous_phys = theta[i0]
            previous_L = L
            previous_R = R

            rigidity_strings = String[]
            for r in rigidity
                push!(rigidity_strings, @sprintf("%.2g", r))
            end
            @info @sprintf("[%s] T=%.1f  %s@%d  |θ0|=%.4f  peak=%.4f  r=[%s]  %.0fs",
                label, T, info[:reason], info[:niters], abs(theta[i0]),
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
    println("[$label] cache: $cachefile  ($n_ok/$(length(Ts)) points done)")
    return done
end

# ── entry point: dispatch on the command-line mode ──────────────────────────────────────────────
mode = length(ARGS) >= 1 ? ARGS[1] : error("usage: julia wall_scan_cluster.jl <rtm|rdm|cutoff>")

const FULL_LADDER = collect(2.0:1.0:14.0)

if mode == "rtm"
    run_wall_scan(chi=64, label="rtm64_full", Ts=FULL_LADDER)
elseif mode == "rdm"
    run_wall_scan(chi=64, label="rdm64", trunc_mode=:rdm, Ts=FULL_LADDER)
elseif mode == "cutoff"
    run_wall_scan(chi=64, label="cut_tight", cutoffs=[fill(1e-10, 40); 1e-12], Ts=FULL_LADDER)
else
    error("unknown mode \"$mode\" — expected one of: rtm, rdm, cutoff")
end
