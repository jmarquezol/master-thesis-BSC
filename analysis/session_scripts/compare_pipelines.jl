using JLD2, Printf

# pipeline comparison: single-vector lambda0 against the eigenvalue-only block ladder at p=0.1
sv = load("../session_caches/svpm_p01.jld2", "res")
bk = load("../../results/data/cluster/sweep_rtm_eigs_p0.1_bulk.jld2", "done")
v = 2.67016
println("single-vector vs eigenvalue-only block, corrected column, p=0.1\n")
@printf("%-4s %-12s %-12s %-11s %-8s\n", "T", "|mu0| sv", "|mu0| block", "rel diff", "x1(block)")
for T in sort(collect(keys(sv)))
    key = ("rtm_eigs_p0.1_bulk", T)
    haskey(bk, key) || continue
    e = bk[key]
    m_sv = abs(sv[T].lambda0); m_bk = abs(e.theta_phys)
    i0 = e.i0
    gaps = [abs(angle(t / e.theta[i0])) for (i, t) in enumerate(e.theta) if i != i0]
    x1 = v * T * minimum(gaps) / pi
    @printf("%-4.0f %-12.6f %-12.6f %-11.2e %-8.3f\n", T, m_sv, m_bk, abs(m_sv - m_bk) / m_bk, x1)
end
println("PIPELINE-COMPARE-DONE")
