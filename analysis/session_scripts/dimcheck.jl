include("../../src/thesislib.jl")
sites = siteinds("S=1/2", 10)
for p in (0.0, 0.1)
    H = MPO(alcaraz_opsum(10, 1.0, p), sites)
    println("p=", p, "  D_W(H MPO) = ", maxlinkdim(H))
    for alg in ("WI", "WII", "VD2")
        U = expH_alcaraz(sites, 1.0, p; dt=0.1, mpo_alg=alg)
        println("   ", alg, ": maxlinkdim(U) = ", maxlinkdim(U), "   linkdims = ", linkdims(U))
    end
end
for p in (0.0, 0.1)
    mpo, scaffold = build_alcaraz_tmpo(2.0; p=p, dt=0.1, nbeta=4, MPO_alg="VD2")
    println("tMPO p=", p, "  temporal phys dims = ", unique(dim.(siteinds(mpo; plev=0))), "  n temporal sites = ", length(mpo))
end
