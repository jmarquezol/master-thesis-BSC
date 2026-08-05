# ──────────────────────────────────────────────────────────────────────────────
# Transverse-contraction toolkit
# ──────────────────────────────────────────────────────────────────────────────


# `overlap_noconj` is the bilinear overlap ⟨L|R⟩ we are interested in
overlap_lr(L::MPS, R::MPS) = overlap_noconj(L, R)

"""
    build_tmpo(mp, scheme, target_T; dt, nbeta, init_state) → (mpo, scaffold)

Model-agnostic tMPO builder. Returns the rotated forward tMPO (the spatial transfer
matrix) and a structured seed tMPS on the same time-site indices (the scaffold). The
scaffold gets overwritten with random complex tensors before any power method.
`init_state` is the single-site boundary state name (e.g. "X+"). Works for any
`ModelParams`/`ExpHRecipe` pair that has an `ITransverse.expH` method.
"""
function build_tmpo(mp::ModelParams, scheme::ExpHRecipe, target_T::Float64;
        dt::Float64=0.1, nbeta::Int=0, init_state::String="X+")
    Nsteps      = round(Int, target_T / dt) + nbeta          # time-steps + nbeta cooling sites
    s           = mp.phys_site
    init        = complex(state(s, init_state))              # single-site LEFT/RIGHT boundary
    tp          = tMPOParams(mp=mp, dt=dt, nbeta=nbeta, scheme=scheme, dbeta=-im*dt, bl=init)
    b           = FwtMPOBlocks(tp)                           # rotated bulk/boundary tensors
    # the spatial MPO's VIRTUAL bond => the temporal PHYSICAL dimension (read it dynamically)
    spatial_bond_dim = dim(inds(b.Wc, "Site,time")[1])
    time_sites  = addtags(siteinds(spatial_bond_dim, Nsteps; conserve_qns=false), "time")
    mpo         = fw_tMPO(b, time_sites, tr=init)            # the transfer matrix (an MPO)
    scaffold    = fw_tMPS(b, time_sites; tr=init, LR=:right) # structured seed tMPS (index skeleton + dims)
    return mpo, scaffold
end

"""
    build_alcaraz_tmpo(target_T; p, lambda, dt, nbeta, MPO_alg) → (mpo, scaffold)

Alcaraz-specific thin wrapper around `build_tmpo` (boundary |X+⟩). Kept for the
existing Alcaraz notebooks/sweeps.
"""
function build_alcaraz_tmpo(target_T::Float64;
        p::Float64=0.1, lambda::Float64=1.0, dt::Float64=0.1,
        nbeta::Int=0, MPO_alg::String="VD2")
    recipe = Dict("WI"=>AlcarazWI(), "WII"=>AlcarazWII(), "VD2"=>AlcarazVD2())[MPO_alg]
    return build_tmpo(AlcarazParams(lambda=lambda, p=p), recipe, target_T; dt=dt, nbeta=nbeta)
end

# Weighted sum Σ coeffs[i]*vecs[i]. A plain directsum of k vectors peaks at k times the input bond,
# which is the memory and time peak of a block iteration. `accdim` caps the running sum; the last
# addition stays exact so the caller still chooses the final compression.
function lincomb_mps(coeffs::AbstractVector, vecs::AbstractVector{MPS};
                     cutoff::Float64=1e-12, maxdim::Int=256, do_truncate::Bool=true,
                     accdim::Int=maxdim)
    acc = coeffs[1] * vecs[1]
    last = lastindex(vecs)
    for i in 2:last
        acc = +(acc, coeffs[i] * vecs[i]; alg="directsum")
        i < last && truncate!(acc; cutoff=cutoff, maxdim=accdim)
    end
    do_truncate && truncate!(acc; cutoff=cutoff, maxdim=maxdim)
    return acc
end

"""
    pad_tmps(src::MPS, target_sites::Vector{<:Index}; tailχ=4) → MPS

Remap a converged tMPS `source` (on a shorter time-site set) onto `target_sites` (the
`siteinds` of a LONGER-T scaffold) so it can warm-start the longer-T (block) power method via
`seedL`/`seedR`. The leading `length(source)` sites reuse `source`'s tensors (re-indexed to the new
site indices, bond structure intact); the extra trailing time sites are filled with a random
complex tail of bond dim `tailχ`, fanned out from the last shared site. The PM only has to relax
the freshly-added tail, saving the bulk of the iteration count along a T-ladder.
"""
function pad_tmps(src::MPS, target_sites::Vector{<:Index}; tailχ::Int=4)
    Ns = length(src)
    Nt = length(target_sites)

    Nt >= Ns || error("pad_tmps: target ($Nt sites) shorter than source ($Ns sites)")
    ssrc = siteinds(src)            
    out  = Vector{ITensor}(undef, Nt)              # pre-allocate the Nt tensors of the new MPS
    for i in 1:Ns
        # Reuse each converged tensor verbatim, only RELABELLING its physical leg (`old => new`)
        # The learned bond structure is kept intact
        out[i] = replaceind(src[i], ssrc[i] => target_sites[i])
    end
    if Nt > Ns
        jl       = Index(tailχ, "Link,l=$Ns")       # a fresh bond Index of small dimension tailχ
        out[Ns] *= randomITensor(ComplexF64, jl)    # attach it to the last shared site (`*` contracts)
        prev     = jl
        for i in (Ns+1):Nt
            rl = i < Nt ? Index(tailχ, "Link,l=$i") : nothing   # last site has no right bond
            # each new tensor carries: left bond `prev`, its physical site, and (unless last) a right bond `rl`
            out[i] = isnothing(rl) ? randomITensor(ComplexF64, prev, target_sites[i]) :
                                     randomITensor(ComplexF64, prev, target_sites[i], rl)
            prev = rl
        end
    end
    return normalize(MPS(out))              
end

"""
    block_transfer_eigs(mpo, scaffold; k, maxdim, cutoff, itermax, eps_conv,
                        n_track, cond_thresh, maxdims) → (theta, L, R, info)

Block power method for the leading `k` eigenvalues of the non-Hermitian transfer operator `mpo`,
with separate left/right bases and a non-conjugating (overlap_noconj) oblique Rayleigh-Ritz.

Options worth knowing:
  • `maxdims`/`cutoffs` — optional per-iteration schedules (cheap early, tight later); `nothing`
    uses the fixed `maxdim`/`cutoff` throughout.
  • `trunc_mode` — `:rtm` truncates each (L,R) pair jointly on |R⟩⟨L|: fewer states, but its
    non-Hermitian SVD suffers at the gap closing. `:rdm` (alias `:naive`) truncates L and R
    separately on their own density matrices: loses the L–R coupling, stays well conditioned.
  • `basis` — `:eig` de-mixes onto the Ritz eigenvectors (cond ~ 1/gap); `:schur` uses an
    orthonormalised basis instead, which survives a near-degenerate cluster.
  • `seedL`/`seedR` — warm starts on `siteinds(scaffold)`, padded with random vectors if short.
    See `pad_tmps` for reusing a converged pair across T.
  • `eigvals_only` — spectrum only; forces `:schur` and skips the bi-normalisation.

Stops at `eps_conv` (`reason="converged"`), or after `stuck_after` iterations without improvement
(`reason="stuck"`). info keys: :niters, :reason, :condS, :condS_hist, :dtheta, :theta, :theta_eigen.
"""
function block_transfer_eigs(mpo::MPO, scaffold::MPS;
        k::Int=4, maxdim::Int=256, cutoff::Float64=1e-12,
        itermax::Int=300, eps_conv::Float64=1e-8, n_track::Int=2,
        cond_thresh::Float64=1e10,
        maxdims::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        cutoffs::Union{Nothing,AbstractVector{<:Real}}=nothing,
        trunc_mode::Symbol=:rtm, itermin::Int=20, stuck_after::Int=100,
        accdim::Int=0,
        seedL::Union{Nothing,AbstractVector{MPS}}=nothing,
        seedR::Union{Nothing,AbstractVector{MPS}}=nothing,
        basis::Symbol=:eig,
        eigvals_only::Bool=false)

    # Spectrum-only mode. The eigenvalues survive the wall even when the eigenvectors do not: θ is a
    # Rayleigh quotient, so its error is O(δ²) in the vector error δ. Callers who only want θ (the λ0
    # circle, gap sweeps, the Eq.(3)/(4) reads) can therefore run much coarser. This flag is not a
    # cheaper algorithm — it just de-mixes on the orthonormal Schur basis instead of the 1/gap-
    # conditioned eigenvector one, and skips the bi-normalisation that only the entropy needs. Pass a
    # smaller maxdim/cutoff yourself to actually get the speed-up.
    if eigvals_only
        basis = :schur
    end

    # cap on the running sum inside lincomb_mps; 0 means "same as maxdim"
    accdim = accdim > 0 ? accdim : maxdim

    # SETUP
    sit  = siteinds(scaffold)

    # mpo true TRANSPOSE: 
    mpoT = swapprime(mpo, 0, 1) # swaps the bra/ket legs w/o complex-conjugating, so ⟨L|mpoT|R⟩ = ⟨mpo·L | R⟩

    # helper functions giving the bond-cap / cutoff to use at iteration `it`
    md_at(it)  = maxdims === nothing ? maxdim  : Int(maxdims[min(it, length(maxdims))])
    cut_at(it) = cutoffs === nothing ? cutoff : Float64(cutoffs[min(it, length(cutoffs))])

    # A fresh random complex MPS on the right sites
    rand_mps() = normalize(complex.(randomMPS(sit, linkdims=2k)))

    # Build the initial block of k vectors: 
    # use warm-start seeds `s` if given (padding with random ones if fewer than k were supplied), otherwise all-random
    seed_block(s) = s === nothing ? MPS[rand_mps() for _ in 1:k] :
        MPS[i <= length(s) ? normalize(complex.(s[i])) : rand_mps() for i in 1:k]   # if i <= length(s), there is a seed available
    R = seed_block(seedR)         # k right vectors |R_1..R_k⟩
    L = seed_block(seedL)         # k left  vectors ⟨L_1..L_k|

    # Bookkeeping for the iteration. NaN+0im = "not computed yet" (complex Not-a-Number).
    theta       = fill(NaN + 0im, k)   # current eigenvalue estimates (the "Ritz values")
    theta_prev  = fill(NaN + 0im, k)   # previous iteration's, to measure convergence Δθ
    dtheta_hist = Float64[]            # history of Δθ
    condS_hist  = Float64[]            # history of the overlap-matrix condition number
    condS_last  = NaN
    reason      = "maxiter"            # why we stopped (usually overwritten)
    niters      = 0
    best_dtheta = Inf                  # best Δθ seen so far (for the "stuck" early-stop)
    iters_noimp = 0                    # consecutive iterations with no improvement

    # The power method on k vectors at once: apply the transfer matrix to the block, then project
    # onto that k-dimensional subspace and solve a small k×k eigenproblem (Rayleigh-Ritz) to read
    # the eigenvalues and re-mix the block.
    for it in 1:itermax
        niters = it
        md  = md_at(it)                                 # per-iteration bond-dim cap (ramp or fixed)
        cut = cut_at(it)                                # per-iteration cutoff (schedule or fixed)

        # Apply the transfer matrix (and its transpose) to every right (left) vector
        # `applyn` = apply an MPO to an MPS (then truncate)
        AR  = MPS[applyn(mpo,  R[j]; cutoff=cut, maxdim=md) for j in 1:k]   # |AR_j⟩ = mpo |R_j⟩
        ATL = MPS[applyn(mpoT, L[j]; cutoff=cut, maxdim=md) for j in 1:k]   # ⟨ATL_j| = ⟨L_j| mpo

        # Build two small k×k matrices ("pencil") that represent the operator inside our subspace:
        #   S = overlaps ⟨L_i|R_j⟩ (the subspace "metric")
        #   M = ⟨L_i| mpo |R_j⟩ (the operator)
        S = Matrix{ComplexF64}(undef, k, k)
        M = Matrix{ComplexF64}(undef, k, k)
        for i in 1:k, j in 1:k
            S[i, j] = overlap_noconj(L[i], R[j])
            M[i, j] = overlap_noconj(L[i], AR[j])
        end
        condS_last = cond(S)                            # condition number: how close S is to singular
        push!(condS_hist, condS_last)                   # (IF large => two vectors nearly parallel => trouble)

        # M v = θ S v, mapped to the ordinary problem (S⁻¹M) v = θ v. pinv, not inv: S goes
        # near-singular as the gap closes.
        pS = pinv(S; rtol=1e-12)
        W  = pS * M
        Fr = eigen(W)                                   # Fr.values = θ's, Fr.vectors = right coeffs
        permr = sortperm(abs.(Fr.values); rev=true)     # sort by |θ| descending (largest first)
        theta = Fr.values[permr]                        # the eigenvalue estimates this iteration
        VR    = Fr.vectors[:, permr]                    # matching right mixing-coefficients (columns)

        # Left coefficients from the same decomposition: uᵀM = θ uᵀS means uᵀS is a row of VR⁻¹,
        # so u_j = pinv(S)ᵀ(VR⁻¹)ᵀ[:,j] — exactly paired with theta[j] and bi-orthogonal by
        # construction. A second eigen() with nearest-value matching mispairs inside a cluster.
        VL = transpose(pS) * transpose(pinv(VR; rtol=1e-12))

        if basis === :schur
            # cond(VR) ~ 1/gap near a degenerate cluster, so de-mixing onto the eigenvector basis
            # amplifies noise. QR-orthonormalise the |θ|-sorted coefficients instead: unitary, and
            # the leading columns still span the same subspaces.
            VR = Matrix(qr(VR).Q)
            VL = Matrix(qr(VL).Q)
        elseif basis !== :eig
            error("block_transfer_eigs: unknown basis=$(basis) (use :eig or :schur)")
        end

        # "De-mix": turn the abstract eigen-coefficients (VR, VL columns) back into actual MPS, by
        # taking those linear combinations of the applied vectors AR / ATL. This replaces the old
        # block with k cleanly-separated (approximate) eigenvectors, ready for the next iteration.
        if trunc_mode === :rtm
            # Truncate each (L_j,R_j) pair jointly on the transition matrix |R_j⟩⟨L_j| — no
            # conjugation. Fewer, cleaner states, but the non-Hermitian SVD is ill-conditioned
            # right at the gap closing. Directsum without truncating: truncate_sweep
            # orthogonalises its own inputs.
            Rnew = MPS[lincomb_mps(VR[:, j], AR;  do_truncate=false, accdim=accdim) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; do_truncate=false, accdim=accdim) for j in 1:k]
            for j in 1:k
                res = truncate_sweep(Lnew[j], Rnew[j]; cutoff=cut, maxdim=md)   # joint pair truncation
                Lnew[j], Rnew[j] = res.L, res.R
            end
        elseif trunc_mode === :rdm || trunc_mode === :naive
            # Truncate every L_j and R_j independently on its own density matrix |v⟩⟨v*|. Throws
            # away the L–R coupling, but it is Hermitian and positive, so it stays well conditioned
            # where the RTM SVD scatters. (:naive is the old alias for this route.)
            Rnew = MPS[lincomb_mps(VR[:, j], AR;  cutoff=cut, maxdim=md, accdim=accdim) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; cutoff=cut, maxdim=md, accdim=accdim) for j in 1:k]
        else
            error("block_transfer_eigs: unknown trunc_mode=$(trunc_mode) (use :rtm, :rdm, or :naive)")
        end
        # Re-normalise each new vector (or replace it with a fresh random vector if norm -> 0, Inf)
        for j in 1:k
            nr = norm(Rnew[j]); Rnew[j] = (isfinite(nr) && nr > 1e-300) ? normalize(Rnew[j]) : rand_mps()
            nl = norm(Lnew[j]); Lnew[j] = (isfinite(nl) && nl > 1e-300) ? normalize(Lnew[j]) : rand_mps()
        end
        R, L = Rnew, Lnew                               # adopt the refreshed block for next iteration

        # Convergence / stopping checks (only tracking the leading n_track eigenvalues)
        ntr = min(n_track, k)
        if it > 1 && all(isfinite, theta_prev[1:ntr])
            # Match each previous θ to its nearest current one before differencing: near-degenerate
            # values swap |θ|-sort order between iterations, and an index-wise difference would read
            # that swap as a jump and call it stuck too early.
            dtheta = 0.0
            usedc  = falses(k)
            for j in 1:ntr
                best, bestd = 0, Inf
                for m in 1:k
                    usedc[m] && continue
                    d = abs(theta[m] - theta_prev[j])
                    if isfinite(d) && d < bestd
                        bestd, best = d, m
                    end
                end
                best != 0 && (usedc[best] = true)
                dtheta = max(dtheta, bestd)
            end
            push!(dtheta_hist, dtheta)
            if dtheta < eps_conv                      
                reason = "converged"
                break                                 
            end
            if it >= itermin                    
                if dtheta < best_dtheta
                    best_dtheta = dtheta
                    iters_noimp = 0    # new best, so we reset the "patience" counter
                else
                    iters_noimp += 1                    # no improvement => tick the counter
                end
                if iters_noimp > stuck_after            # plateaued too long => give up
                    reason = "stuck"
                    break
                end
            end
        end
        theta_prev = copy(theta)                        # remember this step's θ for the next comparison

        # If S is badly conditioned, two block directions have nearly merged.
        # => replace last eigenvector (index k) with a fresh random one that is orthogonal to all the others (Gram-Schmidt orthogonalization)
        if condS_last > cond_thresh && k >= 2
            jb = k
            r = rand_mps(); l = rand_mps()
            for a in 1:(k-1)    # take random vector r (l) and subtract away any component it shared with R[1], ..., R[k-1] (L[1], ..., L[k-1])
                # bi-orthogonal projection: the (L_a,R_a) pairs are not bi-normalized during the
                # iteration, so the projection coefficient must be divided by ⟨L_a|R_a⟩
                den = overlap_noconj(L[a], R[a])
                abs(den) < 1e-14 && continue
                r = lincomb_mps([1.0, -overlap_noconj(L[a], r)/den], MPS[r, R[a]]; cutoff=cutoff, maxdim=md)
                l = lincomb_mps([1.0, -overlap_noconj(R[a], l)/den], MPS[l, L[a]]; cutoff=cutoff, maxdim=md)
            end
            R[jb] = normalize(r); L[jb] = normalize(l)
            reason = (reason == "converged") ? reason : "refreshed"
        end
    end


    theta_eigen = copy(theta)         # keep the raw eigenvalues

    # Bi-orthonormalise each pair to ⟨L_j|R_j⟩ = 1. Only the vectors need this (it is what makes
    # the phase rigidity read r_j = 1/‖L_j‖‖R_j‖), so the spectrum-only path skips it.
    if !eigvals_only
        for j in 1:k
            ov = overlap_noconj(L[j], R[j])
            if abs(ov) > 1e-10
                L[j] = (1 / sqrt(ov)) * L[j]
                R[j] = (1 / sqrt(ov)) * R[j]
            end
        end
    end

    info = Dict(:niters => niters, :reason => reason,
                :condS => condS_last, :condS_hist => condS_hist,
                :dtheta => dtheta_hist, :theta => theta,
                :theta_eigen => theta_eigen)
    return theta, L, R, info        # eigenvalues, left block, right block, diagnostics
end

# single-vector LR power method + leading-eigenvalue diagnostic
"""
    run_pm_diagnosed(target_T; p, lambda, dt, maxdim, cutoff, eps_converged,
                     nbeta, MPO_alg, alg, itermax, stuck_after) → NamedTuple

Single-vector `powermethod_lr` wrapper. Returns the bi-normalized (L,R), the
leading Rayleigh-quotient eigenvalue λ₀, the tMPO (reuse it for block_transfer_eigs),
and convergence diagnostics. 
"""
function run_pm_diagnosed(target_T::Float64;
        p::Float64=0.1, lambda::Float64=1.0, dt::Float64=0.1,
        maxdim::Int=256, cutoff::Float64=1e-14, eps_converged::Float64=1e-6,
        nbeta::Int=0, MPO_alg::String="VD2", alg::String="RTM",
        itermax::Int=5000, stuck_after::Int=200)

    mpo, scaffold = build_alcaraz_tmpo(target_T; p=p, lambda=lambda, dt=dt,
                                        nbeta=nbeta, MPO_alg=MPO_alg)
    seed_mps = deepcopy(scaffold)
    for i in eachindex(seed_mps)    # random seed to avoid subdominant-sector trap
        seed_mps[i] = randomITensor(ComplexF64, inds(seed_mps[i]))
    end
    normalize!(seed_mps)                    

    # Power Method params
    pm_params = PMParams(;
        truncp        = (; cutoff=cutoff, maxdim=maxdim, alg=alg),
        opt_method    = :nosym,        
        cutoffs       = [cutoff],
        maxdims       = 2:2:maxdim,
        itermax       = itermax,
        eps_converged = eps_converged,
        normalization = "overlap",
        stuck_after   = stuck_after,
        compute_fidelity = false)

    # Run PM
    psi_L, psi_R, pm_info = ITransverse.powermethod_lr(seed_mps, mpo, mpo, pm_params)

    # Convergence diagnostics
    ds_hist  = pm_info[:ds]                   # ds = per-step change in singular values
    chi_hist = pm_info[:chi]                  # bond dimension used each step
    niters   = length(ds_hist)
    final_ds = isempty(ds_hist) ? NaN : last(ds_hist)
    stuck    = isempty(ds_hist) || final_ds > eps_converged
    reason   = (!stuck) ? "converged" : (niters >= itermax) ? "maxiter" : "stuck"

    # Leading eigenvalue: λ₀ = ⟨L|mpo|R⟩ / ⟨L|R⟩
    lr_overlap_raw = overlap_lr(psi_L, psi_R)
    lambda0 = expval_LR(psi_L, mpo, psi_R) / lr_overlap_raw

    # Bi-normalise so ⟨L|R⟩=1
    c = sqrt(lr_overlap_raw)
    psi_L = (1/c) * psi_L
    psi_R = (1/c) * psi_R

    return (L=psi_L, R=psi_R, mpo=mpo, scaffold=scaffold, lambda0=lambda0,
            niters=niters, stuck=stuck, reason=reason, final_ds=final_ds,
            ds_hist=ds_hist, chi_hist=chi_hist)
end

# generalized temporal entropies from a converged power method:
"""
    compute_entropies(mp::ModelParams, target_T; scheme, dt, cutoff, maxdim, alg,
                      eps_converged, nbeta, use_block_pm, k_block)
        → NamedTuple(bonds, re, im, L, R, mpo)

Builds the forward tMPO for model `mp` with exponentiation `scheme` (e.g. AlcarazVD2()) and
returns the Rényi-2 temporal-entropy profile S₂(t) = -log Tr(T_t²) per internal bond (Re and Im
separately). Initial state |X+⟩ (free BC).

The boundary tMPS (L,R) come from one of two power methods:
  • `use_block_pm=false` (DEFAULT): single-vector `powermethod_lr` from a random seed. Cheap, but
    STOPS CONVERGING once the transfer-matrix gap closes (the entanglement barrier), so the
    profile is only trustworthy for short/intermediate T.
  • `use_block_pm=true`: the oblique block (subspace) method `block_transfer_eigs` with `k_block`
    Ritz vectors, taking the leading (already bi-orthonormal) pair. Stays well-behaved through the
    degeneracy and recovers the conformal dome deeper into the barrier (at higher cost).
"""
function compute_entropies(mp::ModelParams, target_T::Float64;
        scheme::ExpHRecipe, dt::Float64=0.1, cutoff::Float64=1e-12, maxdim::Int=64,
        alg::String="RTM", eps_converged::Float64=1e-6, nbeta::Int=4,
        use_block_pm::Bool=false, k_block::Int=2,
        maxdims::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        cutoffs::Union{Nothing,AbstractVector{<:Real}}=nothing,
        trunc_mode::Symbol=:rtm, init_state::String="X+", basis::Symbol=:eig,
        itermax::Int=8000, stuck_after::Int=2000, seed::Union{Nothing,MPS}=nothing)

    Ntime_steps = round(Int, target_T / dt)
    Nsteps      = Ntime_steps + nbeta
    s           = mp.phys_site
    init        = complex(state(s, init_state))

    tp = tMPOParams(mp=mp, dt=dt, nbeta=nbeta, scheme=scheme, dbeta=-im*dt, bl=init)
    b  = FwtMPOBlocks(tp)
    spatial_bond_dim = dim(inds(b.Wc, "Site,time")[1])
    time_sites = addtags(siteinds(spatial_bond_dim, Nsteps; conserve_qns=false), "time")

    mpo       = fw_tMPO(b, time_sites, tr=init)
    start_mps = fw_tMPS(b, time_sites; tr=init, LR=:right)

    if use_block_pm
        # (A) robust block PM through gap closing
        _, L_vecs, R_vecs, info = block_transfer_eigs(mpo, start_mps;
            k=k_block, maxdim=maxdim, cutoff=cutoff, itermax=itermax, eps_conv=eps_converged,
            maxdims=maxdims, cutoffs=cutoffs, trunc_mode=trunc_mode, basis=basis)
        # Warm only if method is stuck or reached max iterations
        info[:reason] in ("maxiter", "stuck") && @warn "block PM did not strictly converge at T=$target_T (reason=$(info[:reason]))"
        psi_L, psi_R = L_vecs[1], R_vecs[1]            # take the leading (dominant) pair
    else
        # (B) cheap single-vector method
        if seed === nothing
            for i in eachindex(start_mps)              # random seed (avoids the subdominant-sector trap)
                start_mps[i] = randomITensor(ComplexF64, inds(start_mps[i]))
            end
        else
            start_mps = pad_tmps(seed, siteinds(start_mps))   # warm-start: prev-T fixed point onto these time-sites
        end
        normalize!(start_mps)

        pm_params = PMParams(;
            truncp = (; cutoff=cutoff, maxdim=maxdim, alg=alg),
            opt_method = :nosym,
            cutoffs = cutoffs === nothing ? [cutoff] : cutoffs,
            maxdims = maxdims === nothing ? (2:2:maxdim) : maxdims,
            itermax = itermax,
            eps_converged = eps_converged,
            normalization = "overlap",
            stuck_after = stuck_after,
            compute_fidelity = false)

        psi_L, psi_R, _ = ITransverse.powermethod_lr(start_mps, mpo, mpo, pm_params)

        nrm   = overlap_noconj(psi_L, psi_R)
        psi_L = (1/sqrt(nrm)) * psi_L
        psi_R = (1/sqrt(nrm)) * psi_R
    end

    # Rényi-2 temporal entropy
    s2 = ITransverse.gen_renyi2(psi_L, psi_R)

    return (; bonds = 1:length(s2), re = real.(s2), im = imag.(s2),
            L = psi_L, R = psi_R, mpo = mpo)
end

"""
    plot_entropy_profiles(mp, target_times; scheme, dt, ...) → Plots.Plot

Plots Re(S₂) and Im(S₂) temporal-entropy profiles for a list of target times.
"""
function plot_entropy_profiles(mp::ModelParams, target_times::Vector{Float64};
        scheme::ExpHRecipe, dt::Float64=0.1, cutoff::Float64=1e-12, maxdim::Int=64,
        alg::String="RTM", eps_converged::Float64=1e-6, nbeta::Int=4,
        use_block_pm::Bool=false, k_block::Int=2,
        maxdims::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        cutoffs::Union{Nothing,AbstractVector{<:Real}}=nothing, trunc_mode::Symbol=:rtm)

    # Two empty plot for real and imaginary parts
    plt_real = plot(title="Re(S₂)", xlabel="temporal cut t/T", ylabel="Re(S₂)",
                    legend=:outerright, grid=true, framestyle=:box)
    plt_imag = plot(title="Im(S₂)", xlabel="temporal cut t/T", ylabel="Im(S₂)",
                    legend=:outerright, grid=true, framestyle=:box)
    n = length(target_times)
    cr = cgrad(:viridis, n, categorical=true)          # n distinct colours along a gradient (real)
    ci = cgrad(:plasma,  n, categorical=true)          # a second palette for the imaginary panel

    @showprogress "entropy profiles ($(typeof(scheme)))..." for (i, T) in enumerate(target_times)
        res = compute_entropies(mp, T; scheme=scheme, dt=dt, cutoff=cutoff, maxdim=maxdim,
                                alg=alg, eps_converged=eps_converged, nbeta=nbeta,
                                use_block_pm=use_block_pm, k_block=k_block,
                                maxdims=maxdims, cutoffs=cutoffs, trunc_mode=trunc_mode)
        x = range(0.0, 1.0, length=length(res.re))     # rescale the bond index to t/T ∈ [0,1]
        lab = "T = $(round(T, digits=1))"
        plot!(plt_real, x, res.re, label=lab, lw=2, color=cr[i])
        plot!(plt_imag, x, res.im, label=lab, lw=2, color=ci[i])
    end
    return plot(plt_real, plt_imag, layout=(1, 2), size=(1200, 450), margin=5Plots.mm)
end

# TDVP Schrödinger Loschmidt amplitude L(T)=⟨ψ0|U(T)|ψ0⟩ (crash-safe)
"""
    tdvp_loschmidt_amplitude(N, target_times; p, lambda, dt, cutoff, maxdim, cachefile)
      → Dict{Float64, NamedTuple}

Evolves |X+⟩^N with TDVP on the Alcaraz Hamiltonian and records the complex Loschmidt
amplitude at each target time.
"""
function tdvp_loschmidt_amplitude(N::Int, target_times::Vector{Float64};
        p::Float64=0.1, lambda::Float64=1.0, dt::Float64=0.05,
        cutoff::Float64=1e-12, maxdim::Int=256,
        cachefile::Union{String,Nothing}=nothing)

    # Cache file path (build a default name from p and N if none was given)
    cf   = isnothing(cachefile) ? "results/data/tdvp_loschmidt_p$(p)_N$(N).jld2" : cachefile
    done = isfile(cf) ? load(cf, "done") : Dict{Float64,Any}()   # resume from disk, or start fresh

    # Set up real-space problem
    sites = siteinds("S=1/2", N) 
    psi0  = complex(MPS(sites, "X+"))          
    os    = alcaraz_opsum(N, lambda, p)           
    H     = MPO(os, sites)                          

    sorted_Ts  = sort(target_times)
    missing_Ts = [T for T in sorted_Ts if !haskey(done, T)]   # which targets still need computing
    if isempty(missing_Ts)
        @info "All target T values already cached."
        return done
    end

    # Main Loop
    psi_t = deepcopy(psi0)
    current_t = 0.0
    for T in sorted_Ts
        steps = round(Int, (T - current_t) / dt)
        for _ in 1:steps
            # One TDVP time-step
            psi_t = tdvp(H, -im * dt, psi_t; cutoff=cutoff, maxdim=maxdim, nsite=2)
            normalize!(psi_t)
        end
        current_t = T
        haskey(done, T) && (@info "T=$T (cached, evolved through)"; continue)
        # Loschmidt amplitude G = ⟨ψ0|ψ(T)⟩
        G    = inner(psi0, psi_t); absG = abs(G)
        # Store amplitude, its modulus, the "rate" -log|G|/N (intensive), and the bond dim reached
        done[T] = (G=G, absG=absG, rate=-log(max(absG, 1e-50))/N, maxchi=maxlinkdim(psi_t))
        jldsave(cf; done)
        @info "T=$T (NEW)  |G|=$(round(absG,digits=5))  χ=$(maxlinkdim(psi_t))"
        GC.gc()
    end
    return done
end

# generic crash-safe sweep: f(T) per T, checkpoint after each
"""
    crashsafe_sweep(f, Ts; cachefile) → done::Dict

Calls `done[T] = f(T)` for each T (sorted), saving a JLD2 checkpoint after every T
and skipping already-cached T. Use a results/data/ path for `cachefile`.
"""
function crashsafe_sweep(f::Function, Ts; cachefile::String)
    done = isfile(cachefile) ? load(cachefile, "done") : Dict{Float64, Any}()
    for T in sort(collect(Ts))
        haskey(done, T) && continue                    # already computed => skip (the "resume")
        try
            done[T] = f(T)                             # the actual computation
        catch err
            # If f(T) throws, DON'T abort the whole sweep -> record the error and keep going
            @warn "T=$T failed: $err"
            done[T] = (error=string(err),)
        end
        jldsave(cachefile; done)
        GC.gc()
    end
    return done
end

"""
    thesis_plot_theme!()

Set the global `Plots.default` values used for MANUSCRIPT-QUALITY figures (July 2026 thesis pass):
larger fonts, no in-figure titles (the caption does that job in the thesis), thicker lines, a boxed
frame, and a consistent canvas/dpi. Owning notebook figure cells call this once before building a
figure destined for `thesis/imgs/`; notebook-only diagnostic figures can skip it. Idempotent.
"""
function thesis_plot_theme!()
    Plots.default(
        # NOTE: no fontfamily override — GR's Computer Modern breaks the Unicode glyphs (λ₀, Δφ, π)
        # used throughout the axis labels; the default sans font renders them all correctly.
        guidefontsize  = 14,   # axis labels
        tickfontsize   = 12,
        legendfontsize = 11,
        linewidth      = 2.5,
        markersize     = 6,
        framestyle     = :box,
        grid           = true,
        size           = (800, 480),
        dpi            = 200,
        margin         = 5Plots.mm,
        bottom_margin  = 10Plots.mm,   # 14pt axis labels need the room; 5mm clips them
        left_margin    = 10Plots.mm,
    )
    return nothing
end

# One colour and marker per frustration value, shared by every thesis figure.
const P_COLOR  = Dict(0.0 => :dodgerblue, 0.1 => :crimson, 0.3 => :seagreen, 0.5 => :darkorange)
const P_MARKER = Dict(0.0 => :circle, 0.1 => :square, 0.3 => :diamond, 0.5 => :utriangle)

# Save a row of subplots as one figure under results/imgs/
function plot_panels(panels...; filename::String, title::String="",
                     fig_size::Tuple{Int,Int}=(500*length(panels), 480))
    mkpath("results/imgs")
    plt = plot(panels...; layout=(1, length(panels)), size=fig_size,   # one row, N columns
               plot_title=title, margin=5Plots.mm)
    savefig(plt, joinpath("results", "imgs", filename))
    return plt
end

# ── towers, partners, and picking the physical λ0 ───────────────────────────────────────────────
#
# The spectrum comes in ± pairs: λ0 and a partner at nearly the same modulus, about π away in phase.
# So we classify by phase rather than by modulus rank — the tower is everything within π/2 of λ0
# (its own descendants), the partners are the rest.
#
# Two things we got wrong first and fixed (NB5, NB9):
#   - Pick λ0 by largest modulus. Tracking it by nearest complex value looks natural but cascades:
#     the phase winds as ~a·v·T, so the rule latches onto a stationary partner and every later rung
#     re-anchors on that mistake. Continuity is only used to break a sub-percent modulus tie.
#   - The leading block can be λ0 plus nothing but partners, λ1 missing entirely (p≳0.3). Hence the
#     k→6→8 escalation in block_transfer_eigs_adaptive.

# Phase convention used throughout this project: Im(λ) = arg(-τ), i.e. phase measured from -θ.
phase_of(z) = angle(-z)

# Phase difference φ(a) - φ(b), wrapped into (-π, π].
phase_difference(a, b) = mod(phase_of(a) - phase_of(b) + π, 2π) - π

"""
    classify_tower(theta; i0=argmax(abs.(theta)))

Classify every member of a transfer-matrix spectrum `theta` relative to the physical eigenvalue at
index `i0`: `:tower` if within π/2 in phase (λ0's own descendants), `:partner` otherwise (the
π-shifted -λ0-type copies). Returns `(dphi, cls)`, both vectors indexed like `theta`.
"""
function classify_tower(theta; i0::Int=argmax(abs.(theta)))
    dphi = [phase_difference(theta[j], theta[i0]) for j in eachindex(theta)]
    cls  = [abs(d) < pi / 2 ? :tower : :partner for d in dphi]
    return dphi, cls
end

"""
    tower_gap(theta; i0=argmax(abs.(theta)))

The physical gap |λ1|/|λ0|, where λ1 is the largest-modulus TOWER member (excluding i0 itself) —
never a -λ0 partner. Returns `NaN` if the block contains no tower member besides i0 (λ1 is missing
from this block; the caller needs a bigger k — see `block_transfer_eigs_adaptive`).
"""
function tower_gap(theta; i0::Int=argmax(abs.(theta)))
    _, cls = classify_tower(theta; i0=i0)
    tower_moduli = [abs(theta[j]) for j in eachindex(theta) if j != i0 && cls[j] === :tower]
    isempty(tower_moduli) && return NaN
    return maximum(tower_moduli) / abs(theta[i0])
end

"""
    pick_phys_continuity(theta, previous_phys)

Select the physical eigenvalue λ0 = the DOMINANT (largest-modulus) member of `theta`. The physical
branch is the leading eigenvalue of the transfer matrix; its ±π partner sits at *slightly smaller*
modulus for the resolvable frustration range (p≲1), so modulus-dominance is the robust discriminant.
`previous_phys` is used ONLY to break a genuine modulus TIE (top two within 1%, the exact ±π
near-degeneracy), where the physical branch is the closer of the pair in the complex plane.

This deliberately does NOT anchor on the previous complex value in general: because Im(λ0)≈a·v·T the
phase winds ~a·v per T-step, so a nearest-complex-value rule drifts onto a stationary partner and
that slip CASCADES down the ladder (see section note (A) above; the p≥0.5 cluster failure). Pass
`previous_phys=nothing` on the first rung. Returns the index `i0`.
"""
function pick_phys_continuity(theta, previous_phys)
    mags  = abs.(theta)
    order = sortperm(mags, rev=true)
    i1    = order[1]                                   # dominant eigenvalue = physical λ0
    (previous_phys === nothing || length(order) < 2) && return i1
    i2 = order[2]
    if abs(mags[i1] - mags[i2]) / mags[i1] < 0.01     # genuine ±π near-degeneracy → continuity tiebreak
        return abs(theta[i1] - previous_phys) <= abs(theta[i2] - previous_phys) ? i1 : i2
    end
    return i1
end

"""
    phys_lambda0_suspect(theta, i0, previous_phys; tol=0.30) -> Bool

`true` if the selected |λ0| deviates by more than fractional `tol` from the previous rung's — a
non-convergence / wrong-branch red flag (e.g. the p=0.8, T=4 blow-up to |θ0|≈8.5, or a collapse
toward 0). Callers can use it to re-run that point with more iterations / a larger k rather than
trust it. Returns `false` on the first rung (`previous_phys === nothing`).
"""
function phys_lambda0_suspect(theta, i0, previous_phys; tol::Float64=0.30)
    previous_phys === nothing && return false
    return abs(abs(theta[i0]) - abs(previous_phys)) / abs(previous_phys) > tol
end

"""
    pick_phys_robust(theta, reference; tol=0.05) -> (i0, recovered)

`pick_phys_continuity` restricted to candidates within `tol` of `reference`, the last accepted |λ0|.
A non-converged rung can carry a spurious Ritz value above the true λ0, and plain modulus-dominance
then returns that artefact; the real λ0 is usually still in the block. `tol=0.05` separates the
physical drift (<0.5% per rung) from the contamination (>13%).

`recovered=false` means nothing in the block matches `reference` — drop the rung and keep chaining
`reference` from the last accepted one.
"""
function pick_phys_robust(theta, reference; tol::Float64=0.05)
    reference === nothing && return (pick_phys_continuity(theta, reference), true)
    keep = [i for i in eachindex(theta) if abs(abs(theta[i]) - abs(reference)) / abs(reference) <= tol]
    isempty(keep) && return (pick_phys_continuity(theta, reference), false)
    return (keep[pick_phys_continuity(theta[keep], reference)], true)
end

"""
    block_transfer_eigs_adaptive(mpo, scaffold; k=4, k_retry=8, anchor=nothing,
                                  seedL=nothing, seedR=nothing, kwargs...)

k-adaptive wrapper around `block_transfer_eigs`. Runs at `k` first (optionally warm-started between
T-rungs via `seedL`/`seedR`); while the leading block has no tower member besides the physical λ0
(`tower_gap` returns `NaN` — λ1 absent, only ±π partners present), it **escalates k in steps of 2**
(k→6→8…) up to `k_retry`, warm-seeding each bump with the just-converged block (cheap: the extra
slots are padded with random vectors). A single k=4→6 bump is not always enough at larger frustration
(p≈0.8 needs k=8), which is why this loops rather than doing one retry. Escalation is decided fresh
every call, so a better-conditioned rung drops back to the cheap `k` on its own. Returns
`(theta, L, R, info)` like `block_transfer_eigs`, with `info[:k_used]` and `info[:escalated]` added.
"""
function block_transfer_eigs_adaptive(mpo::MPO, scaffold::MPS;
        k::Int=4, k_retry::Int=8, anchor=nothing,
        seedL::Union{Nothing,AbstractVector{MPS}}=nothing,
        seedR::Union{Nothing,AbstractVector{MPS}}=nothing,
        kwargs...)
    theta, L, R, info = block_transfer_eigs(mpo, scaffold; k=k, seedL=seedL, seedR=seedR, kwargs...)
    k_used = k

    while k_used < k_retry && isnan(tower_gap(theta; i0=pick_phys_continuity(theta, anchor)))
        k_new = min(k_used + 2, k_retry)
        @info "block_transfer_eigs_adaptive: no tower member at k=$k_used, escalating to k=$k_new (warm-seeded)"
        theta, L, R, info = block_transfer_eigs(mpo, scaffold; k=k_new, seedL=L, seedR=R, kwargs...)
        k_used = k_new
    end

    info = merge(info, Dict(:k_used => k_used, :escalated => k_used > k))
    return theta, L, R, info
end
