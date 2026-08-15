include("../../src/thesislib.jl")
using ITensors, ITensorMPS, ITransverse, LinearAlgebra, Printf, JLD2, LsqFit, Statistics

# Same construction as the tower comparison: the production column comes from a 3-site build, the
# corrected one from the middle tensor of a 5-site build, which is a genuine bulk site for NNN.
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

# Renyi-2 profile from the dominant left/right pair, the same quantity the thesis plots.
function entropy_profile(b, tp, init, T; k=4, maxdim=64)
    Nsteps = round(Int, T / tp.dt) + tp.nbeta
    ts  = addtags(siteinds(dim(b.iP), Nsteps; conserve_qns=false), "time")
    mpo = fw_tMPO(b, ts; tr=init)
    _, Lv, Rv, info = block_transfer_eigs(mpo, randomMPS(ts; linkdims=2);
                                          k=k, maxdim=maxdim, cutoff=1e-12,
                                          itermax=300, eps_conv=1e-6, trunc_mode=:rtm)
    return ITransverse.gen_renyi2(Lv[1], Rv[1]), info[:reason]
end

W(t, T) = log((2T / pi) * sin(pi * t / T))
lin(x, q) = q[1] .* x .+ q[2]

# c from the real part: slope against the chord variable on the middle half of the cuts, times 8.
function chord_c(s2, T)
    re = real.(s2); n = length(re)
    ts = range(T/(n+1), T - T/(n+1), length=n)
    bulk = (n÷4):(3n÷4)
    fit = curve_fit(lin, W.(ts[bulk], T), re[bulk], [0.06, 0.5])
    return 8 * fit.param[1]
end

# the imaginary part is flat in the middle; take its central value
plateau(s2) = mean(imag.(s2)[max(1, length(s2)÷2 - 1):(length(s2)÷2 + 2)])

out = "../session_caches/entropy_bulkfix.jld2"
results = isfile(out) ? load(out, "results") : Dict{Tuple{Float64,Symbol,Float64},Any}()

println("Renyi-2 temporal entropies, production column vs corrected bulk column")
println("c_Re from the chord slope, Im plateau against the CFT target pi*c/16 = 0.0982 at c=1/2\n")
@printf("%-5s %-11s %-5s %-9s %-9s %-9s\n", "p", "column", "T", "c_Re", "Im plateau", "reason")
flush(stdout)

# p=0.3 and p=0.5 are the decisive rows: production never forms a conformal dome there
# (dome c at T=2: 1.53 at p=0.3; none at p=0.5). If the corrected column reads ~1/2, the
# p>=0.3 failure story is an artefact of the column; if it reproduces the failure, it stands.
schedule = ((0.1, 2.0:1.0:8.0), (0.3, 2.0:1.0:4.0), (0.5, 2.0:1.0:3.0), (0.0, 2.0:1.0:8.0))
for (p, Ts) in schedule
    b_prod, b_true, tp, init = blocks_pair(p)
    for T in Ts
        for (name, b) in ((:production, b_prod), (:corrected, b_true))
            haskey(results, (p, name, T)) && continue
            s2, reason = entropy_profile(b, tp, init, T)
            results[(p, name, T)] = (; s2, reason)
            @printf("%-5.1f %-11s %-5.1f %-9.4f %-9.4f %-9s\n",
                    p, name, T, chord_c(s2, T), plateau(s2), reason)
            flush(stdout)
            jldsave(out; results)
        end
    end
end

# Does a + b*T^(-1/2) describe each part? Fit both and report the extrapolated value.
@. sqrt_model(x, q) = q[1] + q[2] * x
println("\nfinite-time fit  y = a + b*T^(-1/2)  over the rungs above\n")
@printf("%-5s %-11s %-12s %-9s %-9s %-9s\n", "p", "column", "quantity", "a", "b", "max resid")
for p in (0.1, 0.0), name in (:production, :corrected)
    Ts = sort([T for (pp, nn, T) in keys(results) if pp == p && nn == name])
    isempty(Ts) && continue
    x = [T^(-0.5) for T in Ts]
    for (label, f) in (("c_Re", s2 -> chord_c(s2, 0.0)), ("Im plateau", plateau))
        y = label == "c_Re" ? [chord_c(results[(p, name, T)].s2, T) for T in Ts] :
                              [plateau(results[(p, name, T)].s2) for T in Ts]
        fit = curve_fit(sqrt_model, x, y, [0.5, 0.1])
        resid = maximum(abs.(y .- sqrt_model(x, fit.param)))
        @printf("%-5.1f %-11s %-12s %-9.4f %-9.4f %-9.4f\n",
                p, name, label, fit.param[1], fit.param[2], resid)
    end
end
