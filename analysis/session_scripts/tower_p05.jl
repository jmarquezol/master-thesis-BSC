include("../../src/thesislib.jl")
using ITensors, ITensorMPS, ITransverse, LinearAlgebra, Printf

# Production blocks come from a 3-site build; the corrected ones take the middle tensor of a
# 5-site build, which is a genuine bulk site once the NNN term is present.
function blocks_pair(p; dt=0.1, nbeta=4, nsites=5)
    mp     = AlcarazParams(lambda=1.0, p=p)
    recipe = AlcarazVD2()
    init   = complex(state(mp.phys_site, "X+"))
    tp     = tMPOParams(mp=mp, dt=dt, nbeta=nbeta, scheme=recipe, dbeta=-im*dt, bl=init)
    b_prod = FwtMPOBlocks(tp)

    ss  = [addtags(sim(mp.phys_site), "Site") for _ in 1:nsites]
    U   = ITransverse.expH(ss, mp, recipe; dt=tp.dt)
    Uim = ITransverse.expH(ss, mp, recipe; dt=tp.dbeta)
    mid = (nsites + 1) ÷ 2
    L, Lim = linkinds(U), linkinds(Uim)
    icP = ss[mid]
    time_P  = sim(L[mid-1], tags="Site,time")
    time_vL = sim(icP,  tags="Link,time")
    time_vR = sim(icP', tags="Link,time")
    Wc   = replaceinds(U[mid],   (L[mid-1], L[mid], icP, icP'),     (time_P', time_P, time_vL, time_vR))
    Wcim = permute(replaceinds(Uim[mid], (Lim[mid-1], Lim[mid], icP, icP'),
                               (time_P', time_P, time_vL, time_vR)), inds(Wc)...)
    b_true = FwtMPOBlocks(b_prod; Wc=Wc, Wc_im=Wcim, iL=time_vL, iR=time_vR, iP=time_P, iPs=time_P')
    return b_prod, b_true, tp, init
end

function spectrum(b, tp, init, T; k=8, maxdim=40)
    Nsteps = round(Int, T / tp.dt) + tp.nbeta
    ts  = addtags(siteinds(dim(b.iP), Nsteps; conserve_qns=false), "time")
    mpo = fw_tMPO(b, ts; tr=init)
    theta, _, _, info = block_transfer_eigs(mpo, randomMPS(ts; linkdims=2);
                                            k=k, maxdim=maxdim, cutoff=1e-12,
                                            eigvals_only=true, itermax=200)
    ord = sortperm(abs.(theta), rev=true)
    return theta[ord], info[:reason]
end

function report(label, theta, reason)
    dphi, cls = classify_tower(theta)
    npart = count(c -> c === :partner, cls)
    @printf("   %-11s reason=%-9s partners=%d/%d  tower_gap=%.4f\n",
            label, reason, npart, length(cls), tower_gap(theta))
    for j in eachindex(theta)
        @printf("      |mu|=%.5f  dphi/pi=%+.4f  %s\n",
                abs(theta[j]), dphi[j] / pi, cls[j])
    end
end

for p in (0.5,)
    for T in (2.0, 3.0)
        b_prod, b_true, tp, init = blocks_pair(p)
        @printf("p=%.1f T=%.1f   phys dim %d (production) vs %d (bulk fix)\n",
                p, T, dim(b_prod.iP), dim(b_true.iP))
        th_p, r_p = spectrum(b_prod, tp, init, T)
        report("production", th_p, r_p)
        th_t, r_t = spectrum(b_true, tp, init, T)
        report("bulk fix", th_t, r_t)
        flush(stdout)
    end
end
