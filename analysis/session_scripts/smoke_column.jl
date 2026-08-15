include("../../src/thesislib.jl")
using Printf, LinearAlgebra

# builds both modes and checks dims + p=0 tensor-level null
for p in (0.0, 0.1)
    m3, s3 = build_alcaraz_tmpo(1.0; p=p, nbeta=4, column=:legacy3)
    m5, s5 = build_alcaraz_tmpo(1.0; p=p, nbeta=4, column=:bulk5)
    d3 = dim(siteind(m3, 2)); d5 = dim(siteind(m5, 2))
    @printf("p=%.1f  site dim legacy=%d bulk5=%d  sites %d/%d  scaffold ok %s\n",
            p, d3, d5, length(m3), length(m5), length(s5) == length(m5))
    if p == 0.0
        A = Array(m3[3], inds(m3[3])...); B = Array(m5[3], inds(m5[3])...)
        @printf("       p=0 null on the bulk tensor: |legacy - bulk5| = %.2e\n", norm(A - B))
    end
end
println("SMOKE-OK")
