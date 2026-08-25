ROOT = normpath(joinpath(@__DIR__, "..", ".."))
isdefined(Main, :thesis_plot_theme!) || include(joinpath(ROOT, "src", "thesislib.jl"))
using JLD2, Printf, Statistics, LsqFit

# The one coherent entropy-route analysis behind app:re and the app:errors AICc paragraph.
# Per rung: uncorrected chord fit, corrected free-c fit, and the pinned c=1/2 consistency fit
# with the Bou-Comas (iw)^(-1/2) correction and the exact beta0 regulator terms of their Eq. 51.
# Across T: the plateau extrapolation under four correction models compared by AICc.
const CL = joinpath(ROOT, "data", "cluster")
const SC = joinpath(ROOT, "data", "local")
const BETA0 = 0.2

arm(f, l) = Dict(k[2] => v for (k, v) in load(joinpath(CL, f), "done")
                 if k[1] == l && !haskey(v, :error))

# local svpm caches store untrimmed profiles; trim nbeta/2 bonds per end
function local_ladder(path; nbeta=4)
    isfile(path) || return Dict{Float64,Any}()
    d = load(path, "res")
    k = nbeta ÷ 2
    return Dict(T => (s2_base=r.s2[(k + 1):(end - k)], rigidity=r.rigidity, islocal=true) for (T, r) in d)
end
merge_ladders(ds...) = merge(ds...)

wchord(t, T) = (2T / pi) * sin(pi * t / T)
plat4(s2) = mean(imag.(s2)[max(1, length(s2) ÷ 2 - 1):(length(s2) ÷ 2 + 2)])

# cut coordinates of a trimmed profile: t_k = k*dt exactly, dt = T/(n+1)
function cuts(n, T)
    dt = T / (n + 1)
    abs(dt - 0.1) < 1e-9 || abs(dt - 0.05) < 1e-9 || error("unexpected dt=$dt at T=$T, n=$n")
    return collect(dt:dt:(n * dt))
end
window(n, frac) = (lo = max(1, round(Int, n * (1 - frac) / 2) + 1); lo:(n + 1 - lo))

# imaginary part of Eq. 51 at c=1/2 with the exact regulator terms (eta1 = beta0/T, eta2 = 2beta0/T)
function im_cft(t, T)
    eta1, eta2 = BETA0 / T, 2BETA0 / T
    return (1 / 16) * (pi / 2 - eta2 - (eta1 - eta2 * t / T) * cot(pi * t / T))
end

lsq(X, y) = X \ y                       # linear least squares, columns of X = regressors

function rung_fits(s2, T)
    n  = length(s2)
    ts = cuts(n, T)
    Wv = log.(wchord.(ts, T))
    cv = wchord.(ts, T) .^ (-0.5)
    re = real.(s2); im_ = imag.(s2)

    m50 = window(n, 0.5)
    q = lsq([ones(length(m50)) Wv[m50]], re[m50])
    c_unc = 8q[2]

    c90 = window(n, 0.9)
    X3 = [ones(length(c90)) Wv[c90] cv[c90]]
    q3 = lsq(X3, re[c90])
    resid3 = X3 * q3 .- re[c90]
    # standard error of the slope from the linear model covariance
    covm = inv(X3' * X3) * sum(abs2, resid3) / (length(c90) - 3)
    c_corr, se_corr = 8q3[2], 8 * sqrt(covm[2, 2])

    # pinned c = 1/2: Re residual against s0 + a_re w^{-1/2}; Im residual against a_im w^{-1/2}
    yre = re[c90] .- Wv[c90] ./ 16
    q2  = lsq([ones(length(c90)) cv[c90]], yre)
    rms_pin  = sqrt(mean(abs2, [ones(length(c90)) cv[c90]] * q2 .- yre))
    rms_free = sqrt(mean(abs2, resid3))
    yim  = im_[c90] .- im_cft.(ts[c90], T)
    a_im = lsq(reshape(cv[c90], :, 1), yim)[1]
    return (c_unc=c_unc, c_corr=c_corr, se_corr=se_corr,
            s0=q2[1], a_re=q2[2], a_im=a_im, rms_ratio=rms_pin / rms_free)
end

robust_keep(ys) = (med = median(ys); mad = 1.4826 * median(abs.(ys .- med));
                   mad == 0 ? trues(length(ys)) : abs.(ys .- med) .<= 4mad)

# ---- correction models for the plateau extrapolation, compared by AICc ---------------------
@. m_const(x, q) = q[1] + 0 * x
@. m_sqrt(x, q)  = q[1] + q[2] / sqrt(x)
@. m_lin(x, q)   = q[1] + q[2] / x
@. m_two(x, q)   = q[1] + q[2] / x + q[3] / x^2
const MODELS = [("const", m_const, [0.1], 1), ("a+b/sqrt(T)", m_sqrt, [0.1, 0.05], 2),
                ("a+b/T", m_lin, [0.1, 0.05], 2), ("a+b/T+c/T^2", m_two, [0.1, 0.05, 0.0], 3)]

function extrapolate(Ts, ys, tag)
    rows = []
    for (name, m, q0, k) in MODELS
        length(Ts) <= k && continue
        f = curve_fit(m, Ts, ys, q0)
        push!(rows, (name, f.param[1], sqrt(mean(abs2, m(Ts, f.param) .- ys)), k, f.param))
    end
    @printf("  %s over T=%g..%g (%d rungs)\n", tag, first(Ts), last(Ts), length(Ts))
    # residuals are only indicative across models with different parameter counts; the quoted
    # correction is the predicted T^(-1/2), the others test whether the answer depends on it
    @printf("    %-13s %-6s %-9s %-9s %s\n", "model", "npar", "a", "rms", "c=16a/pi")
    for r in rows
        @printf("    %-13s %-6d %-9.4f %-9.2e %.3f\n", r[1], r[4], r[2], r[3], 16r[2] / pi)
    end
    sq = findfirst(r -> r[1] == "a+b/sqrt(T)", rows)
    sq === nothing && return
    b = rows[sq][5][2]
    @printf("    sqrtT model amplitude b=%.4f -> implied a_im = b*sqrt(2/pi) = %.4f\n", b, b * sqrt(2 / pi))
end

# ---- the ladders ---------------------------------------------------------------------------
lad = Dict{String,Any}()
lad["cluster p=0 (merged, round1 + unused bulk)"] =
    merge(arm("sweep_ent_p0.0_bulk.jld2", "ent_p0.0_bulk"), arm("sweep_ent_p0.0.jld2", "ent_p0.0"))
lad["cluster p=0.1"] = arm("sweep_ent_p0.1_bulk.jld2", "ent_p0.1_bulk")
lad["cluster p=0.3"] = arm("sweep_ent_p0.3_bulk.jld2", "ent_p0.3_bulk")
lad["cluster p=0.5"] = arm("sweep_ent_p0.5_bulk.jld2", "ent_p0.5_bulk")
lad["local p=0 chi64"] = merge_ladders(local_ladder(joinpath(SC, "svpm_p00.jld2")),
                                       local_ladder(joinpath(SC, "svpm_p00_ext.jld2")))
lad["local p=0 dt005"] = local_ladder(joinpath(SC, "svpm_p00_dt005.jld2"); nbeta=8)
lad["local p=0.1 chi64"] = merge_ladders(local_ladder(joinpath(SC, "svpm_p01.jld2")),
                                         local_ladder(joinpath(SC, "svpm_p01_ext.jld2")))

const ORDER = ["cluster p=0 (merged, round1 + unused bulk)", "cluster p=0.1", "cluster p=0.3",
               "cluster p=0.5", "local p=0 chi64", "local p=0 dt005", "local p=0.1 chi64"]

for name in ORDER
    d = lad[name]
    isempty(d) && continue
    Ts = sort(collect(keys(d)))
    pls = [plat4(d[T].s2_base) for T in Ts]
    if haskey(first(values(d)), :islocal)
        # local ladders break intermittently, not at a wall: cut at an anomalous rigidity fall
        # (bc_fits criterion), then drop per-rung outliers below
        rigs = [d[T].rigidity for T in Ts]
        ratios = [rigs[i + 1] / rigs[i] for i in 1:(length(rigs) - 1)]
        typical = median(ratios[1:min(4, end)])
        cut = length(Ts)
        for i in eachindex(ratios)
            ratios[i] < typical / 5 && (cut = i; break)
        end
        lastT = Ts[cut]
    else
        # cluster arms end where the plateau leaves the family of its predecessors
        lastT = Ts[end]
        for i in 3:length(Ts)
            if abs(pls[i] - median(pls[1:i-1])) > 0.05
                lastT = Ts[i-1]; break
            end
        end
    end
    keepT = [T for T in Ts if T <= lastT]
    if haskey(first(values(d)), :islocal)
        keepB = robust_keep([plat4(d[T].s2_base) for T in keepT])
        droppedB = keepT[.!keepB]
        keepT = keepT[keepB]
        isempty(droppedB) || println("  [$name] broken rungs dropped by 4xMAD: T = ", join(droppedB, ", "))
    end
    fits = Dict(T => rung_fits(d[T].s2_base, T) for T in keepT)
    println("\n=== $name  (clean to T=$lastT, $(length(keepT)) of $(length(Ts)) rungs) ===")
    @printf("  %-5s %-8s %-8s %-16s %-8s %-8s %-8s %s\n",
            "T", "plat4", "c_unc", "c_corr(se)", "a_re", "a_im", "rmsP/F", "")
    for T in keepT
        f = fits[T]
        @printf("  %-5g %-8.4f %-8.3f %-7.3f(%.3f)%3s %-8.4f %-8.4f %-8.2f\n",
                T, plat4(d[T].s2_base), f.c_unc, f.c_corr, f.se_corr, "", f.a_re, f.a_im, f.rms_ratio)
    end
    # outlier-filter the two T-series independently before extrapolating
    keepP = robust_keep([plat4(d[T].s2_base) for T in keepT])
    TsP = keepT[keepP]
    dropped = setdiff(keepT, TsP)
    isempty(dropped) || println("  plateau outliers dropped: T = ", join(dropped, ", "))
    length(TsP) >= 4 && extrapolate(TsP, [plat4(d[T].s2_base) for T in TsP], "Im plateau")
end
println("\nENTROPY-C-DONE")
