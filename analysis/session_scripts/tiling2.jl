include("../../src/thesislib.jl")
using ITensors, ITensorMPS, LinearAlgebra, Printf
ITensors.disable_warn_order()

function dense_mpo_matrix(tensors, sites)
    c = tensors[1]
    for i in 2:length(tensors); c = c * tensors[i]; end
    cc = combiner(sites...); rc = combiner(prime.(sites)...)
    d = rc * c * cc
    return Matrix(d, combinedind(rc), combinedind(cc))
end

# Control for the tiling test: repeat the middle tensor of a 5-site build (a genuine bulk tensor,
# both its links interior) instead of the 3-site one. If this reproduces the honest 6-site
# operator, the tiling procedure is sound and the 3-site tensor is the thing at fault.
N = 6
for p in (0.1, 0.5), alg in ("WII", "VD2")
    sites = [Index(2, "S=1/2,Site,n=$i") for i in 1:N]
    Ufull = expH_alcaraz(sites, 1.0, p; dt=0.1, mpo_alg=alg)
    A = dense_mpo_matrix([Ufull[i] for i in 1:N], sites)

    s5 = sites[1:5]
    U5 = expH_alcaraz(s5, 1.0, p; dt=0.1, mpo_alg=alg)
    L = linkinds(U5)                      # bonds 1..4 of the 5-site chain
    lmid, rmid = L[2], L[3]               # the interior links of U5[3]
    dim(lmid) == dim(rmid) || @printf("p=%.1f %s: interior links differ (%d vs %d)\n", p, alg, dim(lmid), dim(rmid))

    extra = Index(dim(rmid), "splice")
    T = [U5[1], U5[2], U5[3],
         replaceinds(U5[3], (lmid, rmid, s5[3], s5[3]'), (rmid, extra, sites[4], sites[4]')),
         replaceinds(U5[4], (L[3], s5[4], s5[4]'), (extra, sites[5], sites[5]')),
         replaceinds(U5[5], (L[4], s5[5], s5[5]'), (L[4], sites[6], sites[6]'))]
    B = dense_mpo_matrix(T, sites)
    @printf("p=%.1f %-4s : bulk-tensor tiling  |A-B|/|A| = %.3e   (interior link dim %d)\n",
            p, alg, norm(A - B) / norm(A), dim(lmid))
end
