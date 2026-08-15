include("../../src/thesislib.jl")
using Printf, LinearAlgebra, JLD2

# K-projected ladders on the corrected column at p=0.5: where does each sector's leader sit at
# production scale? The k=8 mixed blocks found no odd-sector state in the top 8; forcing the
# solver into one sector at a time cannot miss it.
p = 0.5
signs = ksector_signs(p; column=:bulk5)
out = Dict{Tuple{Int,Float64},Any}()
cache = "../session_caches/ksec_bulk.jld2"

for T in (2.0, 3.0), sgn in (+1, -1)
    mpo, scaffold = build_alcaraz_tmpo(T; p=p, nbeta=4, column=:bulk5)
    sit = siteinds(scaffold)
    theta, _, _, info = block_transfer_eigs(mpo, scaffold; k=2, maxdim=48, cutoff=1e-12,
        eigvals_only=true, itermax=400,
        project=psi -> project_ksector(psi, signs, sgn))
    ord = sortperm(abs.(theta), rev=true)
    out[(sgn, T)] = (; theta=theta[ord], reason=info[:reason])
    @printf("p=%.1f T=%.1f sector %+d : |mu|=%.5f, %.5f  arg=%.4f, %.4f  (%s)\n",
            p, T, sgn, abs.(theta[ord])..., angle.(theta[ord])..., info[:reason])
    flush(stdout)
    jldsave(cache; out)
end

println("\nsector gap (odd leader / even leader) and relative phase:")
for T in (2.0, 3.0)
    a = out[(+1, T)].theta[1]; b = out[(-1, T)].theta[1]
    lam0, sub = abs(a) >= abs(b) ? (a, b) : (b, a)
    @printf("T=%.1f  ratio=%.4f  dphi/pi=%+.4f\n", T, abs(sub/lam0), angle(sub/lam0)/pi)
end
println("KSEC-BULK-DONE")
