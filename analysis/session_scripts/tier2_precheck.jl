include("../../src/thesislib.jl")
using SparseArrays, LinearAlgebra, Printf, JLD2, LsqFit

# ── check 1: Trotter drift at p = 0.5 / 1.0 / 1.5, dt = 0.1 vs 0.05 (App C.4 diagnostic) ──
println("Trotter drift: |psi| after 8 VD2 layers on |X+>, N=100\n")
@printf("  %-5s %-8s %-12s %-12s\n", "p", "dt", "after 1", "after 8")
for p in (0.5, 1.0, 1.5), dt in (0.1, 0.05)
    sites = siteinds("S=1/2", 100)
    U = expH_alcaraz(sites, 1.0, p; dt=dt, mpo_alg="VD2")
    psi = MPS(sites, "X+")
    n1 = NaN
    for layer in 1:8
        psi = apply(U, psi; cutoff=1e-12, maxdim=128)
        layer == 1 && (n1 = norm(psi))
    end
    @printf("  %-5.1f %-8.2f %-12.6f %-12.6f\n", p, dt, n1, norm(psi))
    flush(stdout)
end

# ── check 2: v_inf at p = 1.0 and 1.5 from the N = 10..18 ring ladder (momentum-only route) ──
function bits_of(s, N); digits(s - 1, base=2, pad=N); end
index_of(b) = 1 + sum(b[k] * 2^(k-1) for k in eachindex(b))

function ring_H(N, lambda, p)
    rows = Int[]; cols = Int[]; vals = Float64[]
    for s in 1:2^N
        b = bits_of(s, N); spin(i) = 1 - 2b[i]
        d = 0.0
        for i in 1:N; d -= spin(i) * spin(mod1(i+1, N)); end
        for i in 1:N; d -= p * spin(i) * spin(mod1(i+2, N)); end
        push!(rows, s); push!(cols, s); push!(vals, d)
        for i in 1:N
            f = copy(b); f[i] = 1 - f[i]
            push!(rows, index_of(f)); push!(cols, s); push!(vals, -lambda)
        end
        for i in 1:N
            f = copy(b); f[i] = 1 - f[i]; j = mod1(i+1, N); f[j] = 1 - f[j]
            push!(rows, index_of(f)); push!(cols, s); push!(vals, -p * lambda)
        end
    end
    sparse(rows, cols, vals, 2^N, 2^N)
end

function shift_op(N)
    rows = Int[]; cols = Int[]
    for s in 1:2^N
        b = bits_of(s, N)
        push!(rows, index_of([b[mod1(i-1, N)] for i in 1:N])); push!(cols, s)
    end
    sparse(rows, cols, ones(2^N), 2^N, 2^N)
end

function ring_velocity(N, p; howmany=14)
    H = ring_H(N, 1.0, p)
    vals, vecs, _ = ITensorMPS.eigsolve(H, randn(2^N), howmany, :SR; ishermitian=true)
    E0 = vals[1]; Tr = shift_op(N)
    ex = vals[2:howmany] .- E0
    V  = vecs[2:howmany]
    # degenerate levels come back as arbitrary mixtures; diagonalize T inside each cluster
    moms = zeros(Int, length(ex))
    i = 1
    while i <= length(ex)
        j = i
        while j < length(ex) && abs(ex[j+1] - ex[i]) < 1e-8; j += 1; end
        C = hcat(V[i:j]...)
        _, W = eigen(C' * (Tr * C))
        B = C * W
        for (kk, col) in enumerate(eachcol(B))
            moms[i + kk - 1] = round(Int, angle(col' * (Tr * col)) * N / (2pi))
        end
        i = j + 1
    end
    z = minimum(ex[moms .== 0]); nz = findall(!=(0), moms)
    isempty(nz) && error("no momentum-carrying state in the lowest $howmany levels at N=$N, p=$p")
    i = nz[argmin(abs.(moms[nz]))]
    return N * (ex[i] - z) / (2pi * abs(moms[i]))
end

println("\nv_inf extrapolation at the tier-2 couplings\n")
cache = "../../results/data/nb4_velocity_sizes.jld2"
vN = load(cache, "v")
@. lin(x, q) = q[1] + q[2] * x
for p in (1.0, 1.5)
    for N in (10, 12, 14, 16, 18)
        haskey(vN, (p, N)) && continue
        vN[(p, N)] = ring_velocity(N, p)
        jldsave(cache; v=vN)
        @printf("  ring p=%.1f N=%d  v=%.4f\n", p, N, vN[(p, N)]); flush(stdout)
    end
    x = [1 / N^2 for N in (10, 12, 14, 16, 18)]
    y = [vN[(p, N)] for N in (10, 12, 14, 16, 18)]
    f = curve_fit(lin, x, y, [2.0, -3.0])
    r = maximum(abs.(y .- lin(x, f.param)))
    @printf("  p=%.1f  v_inf = %.5f  C = %.2f  max resid %.2e\n", p, f.param[1], f.param[2], r)
end

# ── check 3: equilibrium chord fit at p = 1.0 and 1.5, N = 400 (extends the NB4 cache) ──
println("\nequilibrium chord fits, N=400\n")
using ITensorMPS: dmrg, randomMPS
chord_cache = "../../results/data/nb4_chord_N400.jld2"
chord = load(chord_cache, "data")
function vn_profile(psi)
    N = length(psi); S = Float64[]
    for b in 1:N-1
        orthogonalize!(psi, b)
        ii = b == 1 ? (siteind(psi, b),) : (linkind(psi, b-1), siteind(psi, b))
        _, sv, _ = svd(psi[b], ii)
        s = 0.0
        for n in 1:dim(sv, 1)
            pr = sv[n, n]^2
            pr > 1e-12 && (s -= pr * log(pr))
        end
        push!(S, s)
    end
    S
end
@. chordmodel(x, q) = (q[1]/6)*x + q[2]
for p in (1.0, 1.5)
    if !haskey(chord, p)
        sites = siteinds("S=1/2", 400)
        _, psi = dmrg(MPO(alcaraz_opsum(400, 1.0, p), sites), randomMPS(sites, 10);
                      nsweeps=50, maxdim=[20,100,200,400,800,1000], cutoff=[1e-12],
                      noise=[1e-4,1e-5,1e-6,0.0], outputlevel=0)
        chord[p] = vn_profile(psi)
        jldsave(chord_cache; data=chord)
    end
    S = chord[p]; bonds = 1:399
    xv = log.((800/pi) .* sin.(pi .* bonds ./ 400))
    bulk = findall(b -> 100 <= b <= 300, bonds)
    f = curve_fit(chordmodel, xv[bulk], S[bulk], [0.5, 1.0])
    @printf("  p=%.1f  c = %.4f\n", p, f.param[1]); flush(stdout)
end
println("PRECHECK-DONE")
