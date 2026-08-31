# The uncertainties of Table tab:cp, both routes.
#
# Statistical: each rung's set of physical runs is redrawn with repetition and the ladder refit;
# quoted as the standard deviation over the replicas. For the entropy route the redraw moves the
# median plateau of each rung; for the spectral route it changes which admissible seed continues
# the phase ladder (candidates within twice the best continuity miss, after the 0.5% modulus cut).
# Systematic: the standard deviation of the refit over windows shortened by one or two rungs at
# either end. The two are combined in quadrature.
#
# Reads:  data/cluster/sweep_*.jld2, data/local/controls/seedens_p*.jld2
# Run:    julia --project=. scripts/analysis/table1_errors.jl
ROOT = normpath(joinpath(@__DIR__, "..", ".."))
using JLD2, Printf, Statistics, LsqFit, Random
CL=joinpath(ROOT,"data","cluster"); C=joinpath(ROOT,"data","local","controls"); loose(x)=0.05<x<0.20
arm(f,l)=Dict(k[2]=>v for (k,v) in load(joinpath(CL,f),"done") if k[1]==l && !haskey(v,:error))
ph(z)=angle(-z); dphw(x,y)=mod(x-y+pi,2pi)-pi
@. eq3(x,q)=q[1]-pi/x+q[2]/x^2
@. bc(x,q)=q[1]+q[2]/sqrt(x)
V=Dict(0.0=>2.000,0.1=>2.670,0.3=>3.967,0.5=>5.212)
Random.seed!(20260825)

# ---------------- entropy route ----------------
# The lower end of a window is set by requiring the strip to be wide enough in conformal units for
# the asymptotic form to apply: every rung entering a fit has v*T >~ 8. Only p = 0 is affected, where
# v = 2 makes T = 2 and 3 strips of 4 and 6 lattice units; the p = 0 plateau also peaks at T = 4, so
# the monotone decay the fit describes only starts there.
function erungs(P,cut,lo=0.0)
    pts=Tuple{Float64,Vector{Float64}}[]
    for f in readdir(C)
        startswith(f,"seedens_p"*P*"_T") || continue
        T=parse(Float64,replace(replace(f,"seedens_p"*P*"_T"=>""),".jld2"=>""))
        lo<=T<=cut || continue
        g=filter(loose,[v.plateau for v in values(load(joinpath(C,f),"res"))])
        length(g)>=3 && push!(pts,(T,g))
    end
    sort(pts,by=x->x[1])
end
ec(ps)=16*curve_fit(bc,[p[1] for p in ps],[median(p[2]) for p in ps],[0.1,0.05]).param[1]/pi
println("ENTROPY ROUTE   (statistical = seed resampling, systematic = window variants)")
@printf("%-5s %-9s %-11s %-11s %s\n","p","c","stat sd","syst sd","combined")
for (P,lo,cut) in (("0.0",4.,22.),("0.1",4.,13.),("0.3",2.,9.),("0.5",2.,5.))
    ps=erungs(P,cut,lo); base=ec(ps)
    stat=std([ec([(T,[g[rand(1:length(g))] for _ in 1:length(g)]) for (T,g) in ps]) for _ in 1:400])
    variants=[ec(ps[1:end-1]), ec(ps[2:end]), ec(ps[2:end-1]), ec(ps[1:end-2])]
    syst=std(vcat(base,variants))
    @printf("%-5s %-9.4f %-11.4f %-11.4f %.3f\n",P,base,stat,syst,sqrt(stat^2+syst^2))
end
# ---------------- spectral route ----------------
cfg=Dict(
 0.0=>([("sweep_rtm_eigs_p0.0.jld2","rtm_eigs_p0.0"),("sweep_rtm_eigs_p0.0_fine.jld2","rtm_eigs_p0.0_fine"),("sweep_rtm_eigs_p0.0_fineb.jld2","rtm_eigs_p0.0_fineb")],20.,[22.]),
 0.1=>([("sweep_rtm_eigs_p0.1_bulk.jld2","rtm_eigs_p0.1_bulk"),("sweep_rtm_eigs_p0.1_fine_bulk.jld2","rtm_eigs_p0.1_fine_bulk")],9.5,[11.,12.,13.,14.,15.,16.,17.]),
 0.3=>([("sweep_ent_p0.3_bulk.jld2","ent_p0.3_bulk"),("sweep_rtm_eigs_p0.3_fine_bulk.jld2","rtm_eigs_p0.3_fine_bulk")],7.,[8.,9.,10.]),
 0.5=>([("sweep_ent_p0.5_bulk.jld2","ent_p0.5_bulk"),("sweep_rtm_eigs_p0.5_fine_bulk.jld2","rtm_eigs_p0.5_fine_bulk")],6.,Float64[]))
function sbuild(p; rnd=false)
    arms,cut,ex=cfg[p]; cl=merge([arm(f,l) for (f,l) in arms]...)
    Ts=Float64[]; mus=ComplexF64[]
    for T in sort(collect(keys(cl))); T<=cut && (push!(Ts,T); push!(mus,cl[T].theta_phys)); end
    rate=(ph(mus[end])-ph(mus[end-1]))/(Ts[end]-Ts[end-1])
    for T in ex
        fp=joinpath(C,"seedens_p$(p)_T$(T).jld2"); isfile(fp) || continue
        vs=collect(values(load(fp,"res"))); isempty(vs) && continue
        dT=T-Ts[end]; pred=ph(mus[end])+rate*dT; pm=abs(mus[end])
        ok=[v for v in vs if abs(abs(v.lambda0)-pm)/pm<0.005]; pool=isempty(ok) ? vs : ok
        errs=[abs(dphw(ph(v.lambda0),pred)) for v in pool]; b=minimum(errs)
        cand=findall(e->e<=2*b,errs); i = rnd ? rand(cand) : cand[1]
        push!(Ts,T); push!(mus,pool[i].lambda0); rate=(ph(mus[end])-ph(mus[end-1]))/dT
    end
    phi=[ph(m) for m in mus]; for i in 2:length(phi); phi[i]=phi[i-1]+dphw(phi[i],phi[i-1]); end
    Ts,phi./Ts
end
sc(p,Ts,y,idx)=24*V[p]*abs(curve_fit(eq3,Ts[idx],y[idx],[0.1,-0.01]).param[2])/pi
println("\nSPECTRAL ROUTE")
@printf("%-5s %-9s %-11s %-11s %s\n","p","c","stat sd","syst sd","combined")
for p in (0.0,0.1,0.3,0.5)
    Ts,y=sbuild(p); n=length(Ts); base=sc(p,Ts,y,1:n)
    stat = isempty(cfg[p][3]) ? 0.0 : std([(T2=sbuild(p;rnd=true); sc(p,T2[1],T2[2],1:length(T2[1]))) for _ in 1:150])
    variants=[sc(p,Ts,y,1:n-1), sc(p,Ts,y,2:n), sc(p,Ts,y,3:n), sc(p,Ts,y,2:n-1)]
    syst=std(vcat(base,variants))
    @printf("%-5.1f %-9.4f %-11.4f %-11.4f %.3f\n",p,base,stat,syst,sqrt(stat^2+syst^2))
end
