# --------------------------------------------------------------------------
# Transverse-contraction toolkit
# --------------------------------------------------------------------------
# `overlap_noconj` is the bilinear overlap ⟨L|R⟩ we are interested in
overlap_lr(L::MPS, R::MPS) = overlap_noconj(L, R)

"""
    bulk_fwtmpoblocks(tp; nsites=5) → FwtMPOBlocks

Build U on 5 sites and take the bulk tensor to construct the transfer matrix column.
"""
function bulk_fwtmpoblocks(tp::tMPOParams; nsites::Int=5)
    b3 = FwtMPOBlocks(tp)
    ss = [addtags(sim(tp.mp.phys_site), "Site") for _ in 1:nsites]
    U = ITransverse.expH(ss, tp.mp, tp.scheme; dt=tp.dt)
    Uim = ITransverse.expH(ss, tp.mp, tp.scheme; dt=tp.dbeta)
    mid = (nsites + 1) ÷ 2
    L, Lim = linkinds(U), linkinds(Uim)
    dim(L[mid-1]) == dim(L[mid]) ||
        error("middle tensor is still edge-affected (bond $(dim(L[mid-1]))≠$(dim(L[mid]))); increase nsites")
    icP = ss[mid]
    time_P = sim(L[mid-1], tags="Site,time")
    time_vL = sim(icP, tags="Link,time")
    time_vR = sim(icP', tags="Link,time")
    Wc = replaceinds(U[mid], (L[mid-1], L[mid], icP, icP'), (time_P', time_P, time_vL, time_vR))
    Wcim = ITensors.permute(replaceinds(Uim[mid], (Lim[mid-1], Lim[mid], icP, icP'),
            (time_P', time_P, time_vL, time_vR)), inds(Wc)...)
    return FwtMPOBlocks(b3; Wc=Wc, Wc_im=Wcim, iL=time_vL, iR=time_vR, iP=time_P, iPs=time_P')
end

"""
    build_tmpo(mp, scheme, target_T; dt, nbeta, init_state, column) → (mpo, scaffold)

Model-agnostic tMPO builder: the rotated forward tMPO plus a seed tMPS on the same
time sites (overwritten with random tensors before any power method). `column=:legacy3`
is ITransverse's 3-site extraction (exact only for NN models); `:bulk5` is the corrected
bulk tensor for NNN models.
"""
function build_tmpo(mp::ModelParams, scheme::ExpHRecipe, target_T::Float64;
    dt::Float64=0.1, nbeta::Int=0, init_state::String="X+",
    init_state_top::String=init_state, column::Symbol=:legacy3)
    Nsteps = round(Int, target_T / dt) + nbeta          # time-steps + nbeta cooling sites
    s = mp.phys_site
    init = complex(state(s, init_state))              # bottom (t=0) temporal boundary
    init_top = complex(state(s, init_state_top))          # top boundary; defaults to the same
    tp = tMPOParams(mp=mp, dt=dt, nbeta=nbeta, scheme=scheme, dbeta=-im * dt, bl=init)
    column in (:legacy3, :bulk5) || error("unknown column mode $column")
    b = column === :bulk5 ? bulk_fwtmpoblocks(tp) : FwtMPOBlocks(tp)
    # the spatial MPO's VIRTUAL bond => the temporal PHYSICAL dimension (read it dynamically)
    spatial_bond_dim = dim(inds(b.Wc, "Site,time")[1])
    time_sites = addtags(siteinds(spatial_bond_dim, Nsteps; conserve_qns=false), "time")
    mpo = fw_tMPO(b, time_sites, bl=init, tr=init_top)   # the transfer matrix (an MPO)
    # legacy: structured fw_tMPS seed; bulk5: random seed (fw_tMPS needs the legacy edge tensors)
    scaffold = column === :bulk5 ?
               normalize(complex.(randomMPS(time_sites; linkdims=4))) :
               fw_tMPS(b, time_sites; tr=init, LR=:right)
    return mpo, scaffold
end

"""
    build_alcaraz_tmpo(target_T; p, lambda, dt, nbeta, MPO_alg) → (mpo, scaffold)

Alcaraz-specific thin wrapper around `build_tmpo` (boundary |X+⟩). Kept for the
existing Alcaraz notebooks/sweeps.
"""
function build_alcaraz_tmpo(target_T::Float64;
    p::Float64=0.1, lambda::Float64=1.0, dt::Float64=0.1,
    nbeta::Int=0, MPO_alg::String="VD2", column::Symbol=:legacy3,
    init_state::String="X+", init_state_top::String=init_state)
    recipe = Dict("WI" => AlcarazWI(), "WII" => AlcarazWII(), "VD2" => AlcarazVD2())[MPO_alg]
    return build_tmpo(AlcarazParams(lambda=lambda, p=p), recipe, target_T;
        dt=dt, nbeta=nbeta, column=column, init_state=init_state, init_state_top=init_state_top)
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
    out = Vector{ITensor}(undef, Nt)              # pre-allocate the Nt tensors of the new MPS
    for i in 1:Ns
        # Reuse each converged tensor verbatim, only RELABELLING its physical leg (`old => new`)
        # The learned bond structure is kept intact
        out[i] = replaceind(src[i], ssrc[i] => target_sites[i])
    end
    if Nt > Ns
        jl = Index(tailχ, "Link,l=$Ns")       # a fresh bond Index of small dimension tailχ
        out[Ns] *= randomITensor(ComplexF64, jl)    # attach it to the last shared site (`*` contracts)
        prev = jl
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
    eigvals_only::Bool=false,
    project::Union{Nothing,Function}=nothing)

    # Spectrum-only mode: de-mix on the orthonormal Schur basis instead
    if eigvals_only
        basis = :schur
    end

    # cap on the running sum inside lincomb_mps; 0 means "same as maxdim"
    accdim = accdim > 0 ? accdim : maxdim

    # SETUP
    sit = siteinds(scaffold)

    # mpo true TRANSPOSE: 
    mpoT = swapprime(mpo, 0, 1) # swaps the bra/ket legs w/o complex-conjugating, so ⟨L|mpoT|R⟩ = ⟨mpo·L | R⟩

    # helper functions giving the bond-cap / cutoff to use at iteration `it`
    md_at(it) = maxdims === nothing ? maxdim : Int(maxdims[min(it, length(maxdims))])
    cut_at(it) = cutoffs === nothing ? cutoff : Float64(cutoffs[min(it, length(cutoffs))])

    # A fresh random complex MPS on the right sites
    rand_mps() = normalize(complex.(randomMPS(sit, linkdims=2k)))

    # Build the initial block of k vectors: 
    # use warm-start seeds `s` if given (padding with random ones if fewer than k were supplied), otherwise all-random
    seed_block(s) = s === nothing ? MPS[rand_mps() for _ in 1:k] :
                    MPS[i <= length(s) ? normalize(complex.(s[i])) : rand_mps() for i in 1:k]   # if i <= length(s), there is a seed available
    R = seed_block(seedR)         # k right vectors |R_1..R_k⟩
    L = seed_block(seedL)         # k left  vectors ⟨L_1..L_k|

    # Bookkeeping for the iteration
    theta = fill(NaN + 0im, k)   # current eigenvalue estimates (the "Ritz values")
    theta_prev = fill(NaN + 0im, k)   # previous iteration's, to measure convergence Δθ
    dtheta_hist = Float64[]            # history of Δθ
    condS_hist = Float64[]            # history of the overlap-matrix condition number
    condS_last = NaN
    reason = "maxiter"            # why we stopped (usually overwritten)
    niters = 0
    best_dtheta = Inf                  # best Δθ seen so far (for the "stuck" early-stop)
    iters_noimp = 0                    # consecutive iterations with no improvement

    # The power method on k vectors at once: apply the transfer matrix to the block, then project
    # onto that k-dimensional subspace and solve a small k×k eigenproblem (Rayleigh-Ritz) to read
    # the eigenvalues and re-mix the block.
    for it in 1:itermax
        niters = it
        md = md_at(it)                                 # per-iteration bond-dim cap (ramp or fixed)
        cut = cut_at(it)                                # per-iteration cutoff (schedule or fixed)

        # Apply the transfer matrix (and its transpose) to every right (left) vector
        # `applyn` = apply an MPO to an MPS (then truncate)
        AR = MPS[applyn(mpo, R[j]; cutoff=cut, maxdim=md) for j in 1:k]   # |AR_j⟩ = mpo |R_j⟩
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
        W = pS * M
        Fr = eigen(W)                                   # Fr.values = θ's, Fr.vectors = right coeffs
        permr = sortperm(abs.(Fr.values); rev=true)     # sort by |θ| descending (largest first)
        theta = Fr.values[permr]                        # the eigenvalue estimates this iteration
        VR = Fr.vectors[:, permr]                    # matching right mixing-coefficients (columns)

        # Left coefficients from the same decomposition
        VL = transpose(pS) * transpose(pinv(VR; rtol=1e-12))

        if basis === :schur
            # QR de-mixing to avoid noise amplification when eigenvectors are nearly parallel
            VR = Matrix(qr(VR).Q)
            VL = Matrix(qr(VL).Q)
        elseif basis !== :eig
            error("block_transfer_eigs: unknown basis=$(basis) (use :eig or :schur)")
        end

        # "De-mix": turn the abstract eigen-coefficients (VR, VL columns) back into actual MPS, by
        # taking those linear combinations of the applied vectors AR / ATL
        if trunc_mode === :rtm
            # Truncate each (L_j,R_j) pair jointly on the transition matrix |R_j⟩⟨L_j|, no conjugation, cheaper
            Rnew = MPS[lincomb_mps(VR[:, j], AR; do_truncate=false, accdim=accdim) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; do_truncate=false, accdim=accdim) for j in 1:k]
            for j in 1:k
                res = truncate_sweep(Lnew[j], Rnew[j]; cutoff=cut, maxdim=md)   # joint pair truncation
                Lnew[j], Rnew[j] = res.L, res.R
            end
        elseif trunc_mode === :rdm || trunc_mode === :naive
            # Truncate every L_j and R_j independently on its own density matrix |v⟩⟨v*|. Throws
            # away the L–R coupling, but it is Hermitian and positive, so it stays well conditioned
            # where the RTM SVD scatters
            Rnew = MPS[lincomb_mps(VR[:, j], AR; cutoff=cut, maxdim=md, accdim=accdim) for j in 1:k]
            Lnew = MPS[lincomb_mps(VL[:, j], ATL; cutoff=cut, maxdim=md, accdim=accdim) for j in 1:k]
        else
            error("block_transfer_eigs: unknown trunc_mode=$(trunc_mode) (use :rtm, :rdm, or :naive)")
        end
        # Re-normalise each new vector (or replace it with a fresh random vector if norm -> 0, Inf)
        for j in 1:k
            nr = norm(Rnew[j])
            Rnew[j] = (isfinite(nr) && nr > 1e-300) ? normalize(Rnew[j]) : rand_mps()
            nl = norm(Lnew[j])
            Lnew[j] = (isfinite(nl) && nl > 1e-300) ? normalize(Lnew[j]) : rand_mps()
        end
        # symmetry projection: truncation does not preserve a global sector, so leaked components
        # of the other sector grow back under iteration unless removed every step
        if project !== nothing
            for j in 1:k
                Rnew[j] = normalize(project(Rnew[j]))
                Lnew[j] = normalize(project(Lnew[j]))
            end
        end
        R, L = Rnew, Lnew                               # adopt the refreshed block for next iteration

        # Convergence / stopping checks (only tracking the leading n_track eigenvalues)
        ntr = min(n_track, k)
        if it > 1 && all(isfinite, theta_prev[1:ntr])
            # Match each previous θ to its nearest current one before differencing
            dtheta = 0.0
            usedc = falses(k)
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
            r = rand_mps()
            l = rand_mps()
            for a in 1:(k-1)    # take random vector r (l) and subtract away any component it shared with R[1], ..., R[k-1] (L[1], ..., L[k-1])
                # bi-orthogonal projection: the (L_a,R_a) pairs are not bi-normalized during the
                # iteration, so the projection coefficient must be divided by ⟨L_a|R_a⟩
                den = overlap_noconj(L[a], R[a])
                abs(den) < 1e-14 && continue
                r = lincomb_mps([1.0, -overlap_noconj(L[a], r) / den], MPS[r, R[a]]; cutoff=cutoff, maxdim=md)
                l = lincomb_mps([1.0, -overlap_noconj(R[a], l) / den], MPS[l, L[a]]; cutoff=cutoff, maxdim=md)
            end
            R[jb] = normalize(r)
            L[jb] = normalize(l)
            reason = (reason == "converged") ? reason : "refreshed"
        end
    end


    theta_eigen = copy(theta)         # keep the raw eigenvalues

    # Bi-orthonormalise each pair to ⟨L_j|R_j⟩ = 1. Only the vectors need this, so the
    # spectrum-only path skips it. After normalising, `rigidity` = 1/‖L_j‖‖R_j‖ measures how far
    # the pair is from biorthogonal: it is an internal diagnostic
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
    itermax::Int=5000, stuck_after::Int=200, column::Symbol=:legacy3,
    seed::Union{Nothing,MPS}=nothing)

    mpo, scaffold = build_alcaraz_tmpo(target_T; p=p, lambda=lambda, dt=dt,
        nbeta=nbeta, MPO_alg=MPO_alg, column=column)
    # warm start from a converged vector of the previous rung when given
    if seed === nothing
        seed_mps = deepcopy(scaffold)
        for i in eachindex(seed_mps)
            seed_mps[i] = randomITensor(ComplexF64, inds(seed_mps[i]))
        end
    else
        seed_mps = pad_tmps(seed, siteinds(scaffold))
    end
    normalize!(seed_mps)

    # Power Method params
    pm_params = PMParams(;
        truncp=(; cutoff=cutoff, maxdim=maxdim, alg=alg),
        opt_method=:nosym,
        cutoffs=[cutoff],
        maxdims=2:2:maxdim,
        itermax=itermax,
        eps_converged=eps_converged,
        normalization="overlap",
        stuck_after=stuck_after,
        compute_fidelity=false)

    # Run PM
    psi_L, psi_R, pm_info = ITransverse.powermethod_lr(seed_mps, mpo, mpo, pm_params)

    # Convergence diagnostics
    ds_hist = pm_info[:ds]                   # ds = per-step change in singular values
    chi_hist = pm_info[:chi]                  # bond dimension used each step
    niters = length(ds_hist)
    final_ds = isempty(ds_hist) ? NaN : last(ds_hist)
    stuck = isempty(ds_hist) || final_ds > eps_converged
    reason = (!stuck) ? "converged" : (niters >= itermax) ? "maxiter" : "stuck"

    # Leading eigenvalue: λ₀ = ⟨L|mpo|R⟩ / ⟨L|R⟩
    lr_overlap_raw = overlap_lr(psi_L, psi_R)
    lambda0 = expval_LR(psi_L, mpo, psi_R) / lr_overlap_raw

    # Bi-normalise so ⟨L|R⟩=1
    c = sqrt(lr_overlap_raw)
    psi_L = (1 / c) * psi_L
    psi_R = (1 / c) * psi_R

    return (L=psi_L, R=psi_R, mpo=mpo, scaffold=scaffold, lambda0=lambda0,
        niters=niters, stuck=stuck, reason=reason, final_ds=final_ds,
        ds_hist=ds_hist, chi_hist=chi_hist)
end

# imaginary plateau of a Renyi-2 profile, cooling bonds trimmed
function plateau_im(s2, nbeta::Int=4)
    prof = s2[nbeta÷2+1:end-nbeta÷2]
    return mean(imag.(prof)[max(1, end ÷ 2 - 1):end÷2+2])
end

"""
    ensemble_profile(p, T; nbeta=4) → NamedTuple or nothing

Representative Renyi-2 profile for one rung of the single-vector seed ensembles. 
A seed counts as physical when its imaginary plateau falls inside a certain band, and 
among those we return the profile whose plateau is closest to their median.

Returns `nothing` when the rung has no cache or when no seed is physical.
"""
function ensemble_profile(p::Real, T::Real; nbeta::Int=4,
    root=normpath(joinpath(@__DIR__, "..")))
    file = joinpath(root, "data", "local", "controls",
        "seedens_p$(float(p))_T$(float(T)).jld2")
    isfile(file) || return nothing
    res = load(file, "res")
    good = [r for r in values(res) if 0.05 < r.plateau < 0.20]
    isempty(good) && return nothing

    plateaus = [r.plateau for r in good]
    pick = good[argmin(abs.(plateaus .- median(plateaus)))]
    trim = nbeta ÷ 2
    return (s2=collect(pick.s2)[(trim+1):(end-trim)], plateau=pick.plateau,
        nseeds=length(res), ngood=length(good))
end

# two runs agree if both the plateau and the leading modulus match
runs_agree(a, b; tol, mu_tol) =
    abs(a.plateau - b.plateau) <= tol * abs(a.plateau) &&
    abs(abs(a.lambda0) - abs(b.lambda0)) <= mu_tol * abs(a.lambda0)

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
    cf = isnothing(cachefile) ?
         normpath(joinpath(@__DIR__, "..", "data", "local", "tdvp_loschmidt_p$(p)_N$(N).jld2")) : cachefile
    done = isfile(cf) ? load(cf, "done") : Dict{Float64,Any}()   # resume from disk, or start fresh

    # Set up real-space problem
    sites = siteinds("S=1/2", N)
    psi0 = complex(MPS(sites, "X+"))
    os = alcaraz_opsum(N, lambda, p)
    H = MPO(os, sites)

    sorted_Ts = sort(target_times)
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
        G = inner(psi0, psi_t)
        absG = abs(G)
        # Store amplitude, its modulus, the "rate" -log|G|/N (intensive), and the bond dim reached
        done[T] = (G=G, absG=absG, rate=-log(max(absG, 1e-50)) / N, maxchi=maxlinkdim(psi_t))
        jldsave(cf; done)
        @info "T=$T (NEW)  |G|=$(round(absG,digits=5))  χ=$(maxlinkdim(psi_t))"
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
        guidefontsize=14,   # axis labels
        tickfontsize=12,
        legendfontsize=11,
        linewidth=2.5,
        markersize=6,
        framestyle=:box,
        grid=true,
        foreground_color_legend=nothing,   # no box around the legend
        background_color_legend=nothing,
        size=(800, 480),
        dpi=200,
        margin=5Plots.mm,
        bottom_margin=6Plots.mm,    # just clear of the axis label; 10mm left a blank band
        # between the label and the caption in every figure
        left_margin=10Plots.mm,
    )
    return nothing
end

"""
    thesis_size(frac; aspect, panels)

Canvas size for a figure included at `frac` of the text width.
"""
function thesis_size(frac::Real; aspect::Real=0.62, textwidth::Real=468)
    target = 0.607
    w = round(Int, frac * textwidth / target)
    return (w, round(Int, w * aspect))
end

# One colour and marker per coupling, shared by every thesis figure.
const P_COLOR = Dict(0.0 => :dodgerblue, 0.1 => :crimson, 0.3 => :seagreen,
    0.5 => :darkorange, 1.0 => :purple)
const P_MARKER = Dict(0.0 => :circle, 0.1 => :square, 0.3 => :diamond,
    0.5 => :utriangle, 1.0 => :star5)

"""
    save_thesis_figure(plt, name)

Save a thesis figure under `figures/` in the three formats we keep: PNG for quick viewing,
PDF for `\\includegraphics` in LaTeX, and SVG for editing in Inkscape.
"""
function save_thesis_figure(plt, name::AbstractString)
    dir = normpath(joinpath(@__DIR__, "..", "figures"))
    mkpath(dir)
    for ext in ("png", "pdf", "svg")
        savefig(plt, joinpath(dir, "$name.$ext"))
    end
    println("wrote figures/$name.{png,pdf,svg}")
    return plt
end

# ---- identifying the physical branch and the tower -------------------------
# The block iteration returns k unlabelled Ritz values. Picking the largest modulus is not
# reliable: the operator is non-normal, so an unconverged subspace can return values outside the
# spectrum. We anchor on continuity in T instead. pick_phys_robust keeps the candidates whose
# modulus lies within `tol` of the previously accepted λ0 and selects among those; if none
# qualifies, the caller discards the time point. The remaining values are then ordered by phase
# distance from λ0, which is the order the conformal tower predicts.
#
# classify_tower makes that split at a phase distance of π/2. Values beyond it are not tower
# members and do not enter the gaps. If the block contains no tower member at all, tower_gap
# returns NaN and block_transfer_eigs_adaptive raises k.

# Phase convention used throughout: Im(λ) = arg(-θ), so phases are measured from -θ.
phase_of(z) = angle(-z)

# Phase difference φ(a) - φ(b), wrapped into (-π, π].
phase_difference(a, b) = mod(phase_of(a) - phase_of(b) + π, 2π) - π

"""
    classify_tower(theta; i0=argmax(abs.(theta)))

Classify every member of a transfer-matrix spectrum `theta` relative to the physical eigenvalue at
index `i0`: `:tower` if within π/2 in phase, so λ0 and its descendants, `:partner` for anything
further away, which does not belong to the tower and is excluded from the gaps.
Returns `(dphi, cls)`, both vectors indexed like `theta`.
"""
function classify_tower(theta; i0::Int=argmax(abs.(theta)))
    dphi = [phase_difference(theta[j], theta[i0]) for j in eachindex(theta)]
    cls = [abs(d) < pi / 2 ? :tower : :partner for d in dphi]
    return dphi, cls
end

"""
    tower_gap(theta; i0=argmax(abs.(theta)))

The physical gap |λ1|/|λ0|, where λ1 is the largest-modulus TOWER member (excluding i0 itself).
Returns `NaN` when the block holds no tower member besides i0, meaning λ1 has not been captured
and the caller needs a larger k.
"""
function tower_gap(theta; i0::Int=argmax(abs.(theta)))
    _, cls = classify_tower(theta; i0=i0)
    tower_moduli = [abs(theta[j]) for j in eachindex(theta) if j != i0 && cls[j] === :tower]
    isempty(tower_moduli) && return NaN
    return maximum(tower_moduli) / abs(theta[i0])
end

"""
    tower_dims(theta, T, v; i0=argmax(abs.(theta)))

Convert the phase gaps of a spectrum into boundary dimensions x_i - x_0 = v·T·|Δφ_i|/π,
sorted ascending, with the i0 entry (zero gap) dropped.
"""
function tower_dims(theta, T::Real, v::Real; i0::Int=argmax(abs.(theta)))
    gaps = [abs(phase_difference(theta[j], theta[i0])) for j in eachindex(theta) if j != i0]
    return sort(v .* T .* gaps ./ pi)
end

"""
    pick_phys_continuity(theta, previous_phys)

Low-level helper: returns the index of the largest-modulus member of `theta`, using
`previous_phys` only to break a tie when the top two are within 1% of each other, in which case it
takes the closer of the two in the complex plane.

This is not the selection rule used in production. On its own, largest modulus is unreliable for a
non-normal operator, because an unconverged subspace can return values above the true spectrum. The
rule the results use is `pick_phys_robust`, which restricts the candidates by continuity in T before
calling this. Pass `previous_phys=nothing` on the first rung. Returns the index `i0`.
"""
function pick_phys_continuity(theta, previous_phys)
    mags = abs.(theta)
    order = sortperm(mags, rev=true)
    i1 = order[1]                                   # dominant eigenvalue = physical λ0
    (previous_phys === nothing || length(order) < 2) && return i1
    i2 = order[2]
    if abs(mags[i1] - mags[i2]) / mags[i1] < 0.01     # near-degenerate moduli: break the tie by continuity
        return abs(theta[i1] - previous_phys) <= abs(theta[i2] - previous_phys) ? i1 : i2
    end
    return i1
end

"""
    pick_phys_robust(theta, reference; tol=0.05) -> (i0, recovered)

The selection rule used for every result: keep only the candidates whose modulus lies within `tol`
of `reference`, the last accepted |λ0|, and choose among those. A rung that has not converged can
carry a spurious Ritz value above the true λ0, which plain modulus-dominance would return; the real
λ0 is usually still in the block.

`recovered=false` means nothing in the block matches `reference` => drop the rung and keep chaining
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
(`tower_gap` returns `NaN`, so λ1 was not captured), it escalates k in steps of 2
(k→6→8…) up to `k_retry`, warm-seeding each bump with the just-converged block. 
Returns `(theta, L, R, info)` like `block_transfer_eigs`, with `info[:k_used]` and `info[:escalated]` added.
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
