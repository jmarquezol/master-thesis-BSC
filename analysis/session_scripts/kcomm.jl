include("../../src/thesislib.jl")
using Printf, LinearAlgebra

# is K an exact symmetry of the corrected column, or did the census mix degenerate clusters?
for p in (0.5, 0.1), column in (:bulk5, :legacy3)
    signs = ksector_signs(p; column=column)
    mpo, _ = build_alcaraz_tmpo(0.3; p=p, nbeta=0, column=column)
    c = mpo[1] * mpo[2] * mpo[3]
    phys = sort([s for s in inds(c) if plev(s) == 0], by=dim)
    cc = combiner(phys...); rc = combiner(prime.(phys)...)
    M = Matrix(rc * c * cc, combinedind(rc), combinedind(cc))
    K = kron(reverse([Diagonal(real.(signs)) for _ in 1:3])...)
    @printf("p=%.1f %-8s |EK-KE|/|E| = %.2e   (d=%d, signs %s)\n",
            p, column, norm(M*K - K*M) / norm(M), length(signs),
            join(Int.(round.(real.(signs))), ""))
    flush(stdout)
end
println("KCOMM-DONE")
