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

# Weighted sum of several MPS:  result = Σ_i coeffs[i] * vecs[i]
#   adding two MPS exactly GROWS the bond dimension (it stacks them), so we
#   add them with the exact "directsum" algorithm and only THEN compress once with an SVD-truncation.
#   `do_truncate=false` skips that final compression and returns the RAW (uncompressed) directsum
function lincomb_mps(coeffs::AbstractVector, vecs::AbstractVector{MPS};
                     cutoff::Float64=1e-12, maxdim::Int=256, do_truncate::Bool=true)
    acc = coeffs[1] * vecs[1]
    for i in 2:lastindex(vecs)
        acc = +(acc, coeffs[i] * vecs[i]; alg="directsum")
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

Block (subspace) power method for the leading `k` eigenvalues of the non-Hermitian
transfer operator `mpo`, with separate left/right bases and a fully non-conjugating
(overlap_noconj) oblique Rayleigh-Ritz.

`maxdims` is an OPTIONAL bond-dim ramp: when given a vector, iteration `it` truncates to
`maxdims[min(it,end)]` (cheap early iterations, then grow to the cap); when `nothing`
(default) the fixed `maxdim` is used at every iteration. 

`cutoffs` is the analogous OPTIONAL per-iteration cutoff schedule (looser early, 
tighter late); `nothing` (default) ⇒ fixed `cutoff`.

`trunc_mode` selects how each de-mixed Ritz pair is truncated — the block analogues of the two
truncations `powermethod_lr` offers via `truncp.alg`:
  • `:rtm` (DEFAULT): RTM truncation. Combine exactly (directsum), then truncate each matched (L,R)
    pair JOINTLY on its bilinear transition matrix |R⟩⟨L| via `truncate_sweep` — the same
    non-conjugating RTM route as `truncp.alg="RTM"`. Optimal for the ⟨L|R⟩ overlap and keeps fewer
    states, but its SVD of a non-Hermitian object is ill-conditioned exactly at the gap closing.
  • `:rdm` (= `:naive`, historical alias): RDM truncation. Truncate each L and R INDEPENDENTLY, each
    on its own Hermitian reduced density matrix |v⟩⟨v*| (the conjugating SVD inside `truncate!`) —
    the block analogue of `truncp.alg="densitymatrix"`. Discards the L–R coupling but stays
    well-conditioned (a positive Hermitian spectrum) through the near-degeneracy.

`itermin`/`stuck_after` give an early-stop: after `itermin` iters, break with `reason="stuck"`
once the tracked Δθ fails to improve for `stuck_after` consecutive iters

`eps_conv` gives the strict `reason="converged"` break.

`seedL`/`seedR` are OPTIONAL warm starts: vectors of MPS on `siteinds(scaffold)` used as the
initial left/right bases (any shorter than `k` is padded with random vectors). Defaults `nothing`
generates fully random seeds. A converged pair reused here (see `pad_tmps` for cross-T) converges fast.

`basis` selects the coefficient basis the block is de-mixed onto each iteration:
  • `:eig` (DEFAULT): the Ritz EIGENVECTOR basis — each (L_j,R_j) is an approximate eigenpair.
    Ill-conditioned when eigenvalues cluster (cond of the eigenvector matrix ~ 1/gap).
  • `:schur`: an ordered SCHUR-like basis (QR-orthonormalized eigen coefficients, leading-|θ|
    first). The coefficient rotations are unitary, so the block stays well-conditioned through a
    near-degenerate cluster; the leading pair still spans the dominant eigenvector, and θ is read
    from the projected pencil as usual.

info keys: :niters, :reason, :condS (final), :condS_hist, :dtheta, :theta, :theta_eigen.
"""
function block_transfer_eigs(mpo::MPO, scaffold::MPS;
        k::Int=4, maxdim::Int=256, cutoff::Float64=1e-12,
        itermax::Int=300, eps_conv::Float64=1e-8, n_track::Int=2,
        cond_thresh::Float64=1e10,
        maxdims::Union{Nothing,AbstractVector{<:Integer}}=nothing,
        cutoffs::Union{Nothing,AbstractVector{<:Real}}=nothing,
        trunc_mode::Symbol=:rtm, itermin::Int=20, stuck_after::Int=100,
        seedL::Union{Nothing,AbstractVector{MPS}}=nothing,
        seedR::Union{Nothing,AbstractVector{MPS}}=nothing,
        basis::Symbol=:eig,
        eigvals_only::Bool=false)

    # EIGENVALUES-ONLY FAST PATH (opt-in; default false keeps the full behaviour, byte-for-byte).
    #
    # Why this is safe and useful: the transfer matrix here is strongly non-normal, so near the
    # entanglement-barrier wall the individual eigenVECTORS become ill-defined (their phase rigidity
    # ⟨L|R⟩/‖L‖‖R‖ → 0). The eigenVALUES do NOT: an eigenvalue read off a Rayleigh quotient
    # θ = ⟨L|E|R⟩/⟨L|R⟩ has error O(δ²) in the eigenvector error δ (it is quadratically insensitive),
    # which is exactly why θ stays accurate to ~1e-8 while the vectors scramble. So a caller who wants
    # ONLY the spectrum — the λ0 circle / dual-unitarity plots, the gap-closing sweeps, the Eq.(3)/(4)
    # eigenvalue c-extractions, none of which need entropy-quality vectors — can afford a coarser
    # truncation and still get trustworthy θ's.
    #
    # This flag does NOT invent a cheaper algorithm (the cost is the 2k MPO applications + the block
    # de-mixing, both unavoidable for the eigenvalues too). It just bundles the two settings that make
    # a coarse run well-conditioned, and nothing else:
    #   • it de-mixes on the ORTHONORMAL (Schur) block basis rather than the eigenvector basis (whose
    #     conditioning ~1/gap blows up as the cluster closes), so the iteration never amplifies noise;
    #   • it skips the final ⟨L_j|R_j⟩ = 1 bi-normalisation below, which only prepares the VECTORS for
    #     entropy use — work this path does not need.
    # The realised speed-up is the CALLER's to take, by ALSO passing a smaller `maxdim`/`cutoff` (or a
    # `maxdims` ramp) than an entropy run would need; the flag simply makes doing so safe and
    # self-documenting. The returned `theta` is computed by the identical Rayleigh–Ritz step as always.
    if eigvals_only
        basis = :schur
    end

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

    # MAIN LOOP
    # The idea (like the ordinary power method, but for k vectors at once): repeatedly apply the
    # transfer matrix to a block of k vectors; the block rotates toward the k LARGEST eigenvectors.
    # Each step we project the operator onto our small k-dim subspace and solve a tiny k×k eigenproblem
    # (the "Rayleigh-Ritz" step) to read off eigenvalue estimates and re-mix the block cleanly.
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

        # Solve the generalised eigenproblem: M v = θ S v 
        # we map it to the ordinary problem (S⁻¹M) v = θ v
        # and we use the PSEUDO-inverse pinv(S) which behaves well even when S is near-singular (gap closing) 
        pS = pinv(S; rtol=1e-12)
        W  = pS * M
        Fr = eigen(W)                                   # Fr.values = θ's, Fr.vectors = right coeffs
        permr = sortperm(abs.(Fr.values); rev=true)     # sort by |θ| descending (largest first)
        theta = Fr.values[permr]                        # the eigenvalue estimates this iteration
        VR    = Fr.vectors[:, permr]                    # matching right mixing-coefficients (columns)

        # LEFT coefficients from the SAME decomposition (no second eigen, no matching heuristic):
        # a left pencil eigenvector obeys uᵀM = θ uᵀS ⇔ (uᵀS) W = θ (uᵀS), i.e. uᵀS is a left
        # eigenvector of W = the matching ROW of VR⁻¹. Hence u_j = pinv(S)ᵀ (VR⁻¹)ᵀ[:,j] — paired
        # with theta[j] EXACTLY and bi-orthogonal (uᵢᵀ S vⱼ = δᵢⱼ) by construction, where the old
        # nearest-value matching of a second independent eigen() could mispair the left/right
        # vectors inside a near-degenerate cluster (the tricritical failure mode).
        VL = transpose(pS) * transpose(pinv(VR; rtol=1e-12))

        if basis === :schur
            # Near a degenerate cluster the eigenvector basis itself is ill-conditioned
            # (cond(VR) ~ 1/gap), so de-mixing onto it amplifies noise. Instead rotate the block
            # onto ordered-Schur-like bases: QR-orthonormalize the |θ|-sorted eigen coefficients.
            # The rotations are unitary (well-conditioned always); the leading j columns span the
            # same leading-j eigenvector subspaces, and column 1 is still the dominant eigenvector.
            VR = Matrix(qr(VR).Q)
            VL = Matrix(qr(VL).Q)
        elseif basis !== :eig
            error("block_transfer_eigs: unknown basis=$(basis) (use :eig or :schur)")
        end

        # "De-mix": turn the abstract eigen-coefficients (VR, VL columns) back into actual MPS, by
        # taking those linear combinations of the applied vectors AR / ATL. This replaces the old
        # block with k cleanly-separated (approximate) eigenvectors, ready for the next iteration.
        if trunc_mode === :rtm
            # :rtm — RTM truncation: the block analogue of powermethod_lr's truncp.alg="RTM"
            # (tlrcontract(::Algorithm"RTM")). Truncate each matched (L_j,R_j) pair JOINTLY on its
            # bilinear transition matrix |R_j⟩⟨L_j| (NO conjugation) via `truncate_sweep`. Keeps only
            # the states that matter for the physical ⟨L|R⟩ structure (fewer states, cleaner), but the
            # non-Hermitian SVD is ill-conditioned right at the gap closing.
            # De-mix with the EXACT directsum (do_truncate=false): truncate_sweep orthogonalizes its
            # inputs itself, so a pre-truncate! here would be redundant.
            Rnew = MPS[lincomb_mps(VR[:, j], AR;  do_truncate=false) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; do_truncate=false) for j in 1:k]
            for j in 1:k
                res = truncate_sweep(Lnew[j], Rnew[j]; cutoff=cut, maxdim=md)   # joint pair truncation
                Lnew[j], Rnew[j] = res.L, res.R
            end
        elseif trunc_mode === :rdm || trunc_mode === :naive
            # :rdm — RDM truncation: the block analogue of powermethod_lr's truncp.alg="densitymatrix"
            # (the generic tlrcontract fallback → tcontract(::Algorithm"densitymatrix")). Truncate every
            # L_j and R_j INDEPENDENTLY, each on its own Hermitian reduced density matrix |v⟩⟨v*| (the
            # conjugating SVD inside `truncate!`). Discards the L–R coupling, but it is a positive
            # Hermitian eigenproblem, so it stays well-conditioned through the near-degeneracy where the
            # RTM SVD scatters. (:naive is the historical alias for this same per-vector route.)
            Rnew = MPS[lincomb_mps(VR[:, j], AR;  cutoff=cut, maxdim=md) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; cutoff=cut, maxdim=md) for j in 1:k]
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
            # How much did the tracked eigenvalues move since last step? Match each tracked
            # PREVIOUS θ to its nearest CURRENT θ (continuity) before differencing: members of a
            # ±pair or a near-degenerate cluster swap |θ|-sort order between iterations, and a
            # raw index-wise difference would register that swap as a spurious Δθ jump
            # (fake non-convergence → premature "stuck").
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
                # bi-orthogonal projection: the (L_a,R_a) pairs are NOT bi-normalized during the
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

    # Bi-orthonormalise each pair so that ⟨L_j|R_j⟩ = 1 exactly.
    # Skipped on the eigenvalues-only path: this step only conditions the VECTORS for entropy use
    # (e.g. so the phase rigidity reads r_j = 1/‖L_j‖‖R_j‖); the eigenvalues in `theta` are final
    # already, so an eigenvalues-only caller pays nothing for a normalisation it will not read.
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
    )
    return nothing
end

# Save a row of subplots as one figure under results/imgs/
function plot_panels(panels...; filename::String, title::String="",
                     fig_size::Tuple{Int,Int}=(500*length(panels), 480))
    mkpath("results/imgs")
    plt = plot(panels...; layout=(1, length(panels)), size=fig_size,   # one row, N columns
               plot_title=title, margin=5Plots.mm)
    savefig(plt, joinpath("results", "imgs", filename))
    return plt
end

# ────────────────────────────────────────────────────────────────────────────────────────────────
# PHASE-BASED SPECTRUM CLASSIFICATION AND CONTINUITY (added July 2026, promoted from NB3/NB5/NB8)
#
# The transfer-matrix spectrum comes in ± pairs: a physical eigenvalue λ0 and a "-λ0 partner" at
# nearly the same modulus but phase shifted by ≈π. Selecting "the physical branch" or "the gap" by
# MODULUS RANK ALONE (comparing only the two largest-|θ| values) silently breaks once frustration
# grows: extra partners crowd into the leading k-block and a rank-based pick can grab a partner
# instead of the genuine λ1 (NB5: at p=0.3, T=4 the whole k=4 block is λ0 plus three partners, with
# λ1 entirely absent), or lose continuity across a cold-started T-ladder (NB9). The fix is to
# classify by PHASE relative to λ0: the "tower" is every eigenvalue within π/2 of λ0's phase (its
# own CFT descendants), the "partners" are everything else (the π-shifted copies). This is exactly
# NB3/NB5's finding and NB8's `classify_block`, promoted here as the one authoritative
# implementation — new drivers should use this, not re-derive a modulus-rank selector.
#
# TWO REFINEMENTS from the July-2026 cluster p-sweep (NB9), both baked into the functions below:
#
#  (A) SELECT λ0 BY MODULUS-DOMINANCE, NOT COMPLEX CONTINUITY.  An earlier version of
#      `pick_phys_continuity` searched the WHOLE block for the member closest in the COMPLEX PLANE to
#      the previous rung's λ0. That drifts and then CASCADES: Im(λ0)≈a·v·T, so the physical λ0's phase
#      winds by ~a·v every T-step (faster at larger v, i.e. larger p). A nearest-complex-value rule
#      then latches onto a stationary near-modulus partner, and because the next rung re-anchors on
#      that slip, one wrong pick corrupts the ENTIRE ladder below it (observed p≥0.5: |λ0| collapsing
#      to 0.16 at T=4, then partners for all T>4). Modulus is the STABLE quantity — the physical λ0 is
#      the dominant (largest-|θ|) eigenvalue and its ±π partner is subdominant for the resolvable
#      range (p≲1) — so we pick by |θ| and use continuity ONLY to break a genuine <1% modulus tie
#      (the exact ±π near-degeneracy). `phys_lambda0_suspect` flags a non-convergence blow-up/collapse.
#  (B) ESCALATE k UNTIL A TOWER MEMBER IS CAPTURED.  At larger p the leading block can be λ0 + ONLY
#      π-partners (tower_gap=NaN, λ1 absent). One bump k=4→6 is not always enough (p≈0.8), so
#      `block_transfer_eigs_adaptive` now steps k→6→8 (warm-seeded) until λ1 appears or k_retry is hit.
# ────────────────────────────────────────────────────────────────────────────────────────────────

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
