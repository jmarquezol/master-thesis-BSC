include("../../src/thesislib.jl")
using ITensors, ITensorMPS, ITransverse, LinearAlgebra, Printf

# Production blocks come from a 3-site build; this rebuilds the repeated column from the middle
# tensor of a 5-site build, which is a genuine bulk site for a NNN model.
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
    Wcim = replaceinds(Uim[mid], (Lim[mid-1], Lim[mid], icP, icP'), (time_P', time_P, time_vL, time_vR))
    Wcim = permute(Wcim, inds(Wc)...)

    b_true = FwtMPOBlocks(b_prod; Wc=Wc, Wc_im=Wcim, iL=time_vL, iR=time_vR, iP=time_P, iPs=time_P')
    return b_prod, b_true, tp, init
end

function leading(mpo, scaffold; k=4, maxdim=64)
    theta, _, _, info = block_transfer_eigs(mpo, scaffold; k=k, maxdim=maxdim, cutoff=1e-12,
                                            eigvals_only=true, itermax=200)
    ord = sortperm(abs.(theta), rev=true)
    return theta[ord], info[:reason]
end

for p in (0.1, 0.5)
    for T in (1.0, 2.0, 3.0)
        b_prod, b_true, tp, init = blocks_pair(p)
        Nsteps = round(Int, T / tp.dt) + tp.nbeta

        ts_p = addtags(siteinds(dim(b_prod.iP), Nsteps; conserve_qns=false), "time")
        mpo_p = fw_tMPO(b_prod, ts_p; tr=init)

        ts_t = addtags(siteinds(dim(b_true.iP), Nsteps; conserve_qns=false), "time")
        mpo_t = fw_tMPO(b_true, ts_t; tr=init)

        th_p, r_p = leading(mpo_p, randomMPS(ts_p; linkdims=2))
        th_t, r_t = leading(mpo_t, randomMPS(ts_t; linkdims=2))

        @printf("p=%.1f T=%.1f  phys dim %d vs %d\n", p, T, dim(b_prod.iP), dim(b_true.iP))
        @printf("   production |mu0|=%.6f  arg=%.6f   (%s)\n", abs(th_p[1]), angle(th_p[1]), r_p)
        @printf("   bulk fix   |mu0|=%.6f  arg=%.6f   (%s)\n", abs(th_t[1]), angle(th_t[1]), r_t)
        @printf("   relative difference in |mu0| = %.3e,  in arg = %.3e\n",
                abs(abs(th_p[1]) - abs(th_t[1])) / abs(th_t[1]),
                abs(angle(th_p[1]) - angle(th_t[1])))
        flush(stdout)
    end
end
