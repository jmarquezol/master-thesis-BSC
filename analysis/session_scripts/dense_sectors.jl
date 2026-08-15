include("../../src/thesislib.jl")
using Printf, LinearAlgebra

# Dense sector census: 3 temporal sites (nbeta=0), every eigenvalue classified by exact K charge.
# Answers whether the corrected column has an odd-sector family near the top of the spectrum.
function census(p, column)
    signs = ksector_signs(p; column=column)
    mpo, _ = build_alcaraz_tmpo(0.3; p=p, nbeta=0, column=column)
    c = mpo[1] * mpo[2] * mpo[3]
    phys = sort([s for s in inds(c) if plev(s) == 0], by=dim)
    cc = combiner(phys...); rc = combiner(prime.(phys)...)
    M = Matrix(rc * c * cc, combinedind(rc), combinedind(cc))
    K = kron(reverse([Diagonal(real.(signs)) for _ in 1:3])...)   # site order of the combiner
    F = eigen(M)
    charges = [real(v' * K * v) / real(v' * v) for v in eachcol(F.vectors)]
    return F.values, charges
end

for p in (0.5, 0.1), column in (:bulk5, :legacy3)
    vals, ch = census(p, column)
    mixed = count(c -> abs(c) < 0.99, ch)
    ord = sortperm(abs.(vals), rev=true)
    lam0 = vals[ord[1]]
    @printf("\np=%.1f %-8s dim=%d  mixed eigenvectors=%d\n", p, column, length(vals), mixed)
    @printf("   top of each sector (|mu|, dphi/pi, charge):\n")
    shown = Dict(1 => 0, -1 => 0)
    for i in ord
        s = ch[i] > 0 ? 1 : -1
        shown[s] >= 5 && continue
        shown[s] += 1
        dphi = angle(vals[i] / lam0) / pi
        @printf("   %+d  |mu|=%.5f  dphi/pi=%+.4f\n", s, abs(vals[i]), dphi)
    end
    flush(stdout)
end
println("CENSUS-DONE")
