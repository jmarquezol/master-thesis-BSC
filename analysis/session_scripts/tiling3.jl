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

# Same test as before, now including p=0, where the model is nearest-neighbour only.
N = 6
for p in (0.0, 0.05, 0.1, 0.3, 0.5), alg in ("VD2",)
    sites = [Index(2, "S=1/2,Site,n=$i") for i in 1:N]
    Ufull = expH_alcaraz(sites, 1.0, p; dt=0.1, mpo_alg=alg)
    s3 = sites[1:3]
    U3 = expH_alcaraz(s3, 1.0, p; dt=0.1, mpo_alg=alg)
    l1, l2 = linkinds(U3)
    links = [Index(dim(l1), "tile,l=$i") for i in 1:(N-1)]
    T = Vector{ITensor}(undef, N)
    T[1] = replaceinds(U3[1], (l1,), (links[1],))
    for i in 2:(N-1)
        T[i] = replaceinds(U3[2], (l1, l2, s3[2], s3[2]'), (links[i-1], links[i], sites[i], sites[i]'))
    end
    T[N] = replaceinds(U3[3], (l2, s3[3], s3[3]'), (links[N-1], sites[N], sites[N]'))
    A = dense_mpo_matrix([Ufull[i] for i in 1:N], sites)
    B = dense_mpo_matrix(T, sites)
    @printf("p=%.2f %-4s : temporal phys dim = %2d (bulk MPO bond %2d)   |A-B|/|A| = %.3e\n",
            p, alg, dim(l1), maxlinkdim(Ufull), norm(A - B) / norm(A))
end
