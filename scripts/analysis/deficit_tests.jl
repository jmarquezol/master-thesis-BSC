ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit
# Candidate origins for the deficit in the fitted 1/T^2 coefficient C of Eq.(3), at p=0 where
# c=1/2 and v=2 are known exactly, so C_pred = (pi/v)(x0 - c/24) = -0.032725 for x0=0.
const CL = joinpath(ROOT, "data", "cluster")
arm(f,l)=Dict(k[2]=>v for (k,v) in load(joinpath(CL,f),"done") if k[1]==l && !haskey(v,:error))
ph(z)=angle(-z); dphw(a,b)=mod(a-b+pi,2pi)-pi
function unwrapped(mus)
    p=[ph(m) for m in mus]; for i in 2:length(p); p[i]=p[i-1]+dphw(p[i],p[i-1]); end; p
end
@. eq3_2(x,q) = q[1] - pi/x + q[2]/x^2              # the production form
@. eq3_3(x,q) = q[1] - pi/x + q[2]/x^2 + q[3]/x^4   # one order further (D/T^3 in Im lambda0)
const CPRED = -pi/2*(0.5/24)

spec = merge(arm("sweep_rtm_eigs_p0.0.jld2","rtm_eigs_p0.0"), arm("sweep_rtm_eigs_p0.0_fine.jld2","rtm_eigs_p0.0_fine"))
Ts = sort(collect(keys(spec))); phi = unwrapped([spec[T].theta_phys for T in Ts]); y = phi ./ Ts

println("predicted C = ", round(CPRED, digits=5), "   (c = 1/2, v = 2, x0 = 0)")
println("\n=== test 1: does adding the next order in 1/T move C towards it? ===")
for Tmin in (2.0, 3.0, 4.0)
    k = findall(T -> T >= Tmin, Ts)
    f2 = curve_fit(eq3_2, Ts[k], y[k], [0.1,-0.01]); f3 = curve_fit(eq3_3, Ts[k], y[k], [0.1,-0.01,0.0])
    @printf("  T>=%.0f (%2d pts)  two-term C=%.5f (%.0f%% of pred)   three-term C=%.5f (%.0f%%)\n",
            Tmin, length(k), f2.param[2], 100 * f2.param[2] / CPRED, f3.param[2], 100 * f3.param[2] / CPRED)
end

println("\n=== test 2: does the regulator beta0 bias C? ===")
b = load(joinpath(CL,"sweep_beta_p0.0.jld2"),"done")
nb_of(k)=parse(Int, match(r"nb(\d+)",k[1]).captures[1])
for nb in sort(unique(nb_of.(keys(b))))
    ks = sort([k for k in keys(b) if nb_of(k)==nb], by=k->k[2])
    Tb = [k[2] for k in ks]
    length(Tb) < 5 && continue
    yb = unwrapped([b[k].theta_phys for k in ks]) ./ Tb
    f = curve_fit(eq3_2, Tb, yb, [0.1,-0.01])
    @printf("  nbeta=%-3d beta0=%.2f  T=%g..%g (%d pts)  C=%.5f  (%.0f%% of pred)\n",
            nb, nb*0.05, first(Tb), last(Tb), length(Tb), f.param[2], 100 * f.param[2] / CPRED)
end

println("\n=== test 3: is there a non-universal constant hiding in B? ===")
@. eq3_3free(x,q) = q[1] + q[2]/x + q[3]/x^2 + q[4]/x^4
for Tmin in (2.0, 3.0, 4.0)
    k = findall(T -> T >= Tmin, Ts)
    f2f = curve_fit(eq3_2, Ts[k], y[k], [0.1,-0.01])
    ff  = curve_fit((x,q)->q[1] .+ q[2]./x .+ q[3]./x.^2, Ts[k], y[k], [0.1,-pi,-0.01])
    f3f = curve_fit(eq3_3free, Ts[k], y[k], [0.1,-pi,-0.01,0.0])
    @printf("  T>=%.0f  two-term B free: B/pi=%+.4f c=%.3f | three-term B free: B/pi=%+.4f c=%.3f | three-term B pinned: c=%.3f\n",
            Tmin, ff.param[2]/pi, 24*2*abs(ff.param[3])/pi,
            f3f.param[2]/pi, 24*2*abs(f3f.param[3])/pi,
            24*2*abs(curve_fit(eq3_3, Ts[k], y[k], [0.1,-0.01,0.0]).param[2])/pi)
end
