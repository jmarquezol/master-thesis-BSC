# Universal properties from quantum many-body dynamics

Master's thesis — **Joaquín G. Márquez Olguín**, supervised by **Stefano Carignano** (Barcelona Supercomputing Center), 2026.

## What this is about

At a critical point, the real-time dynamics of a simple initial state carries the universal data of the underlying conformal field theory: the central charge and the boundary operator content. The route is the Loschmidt echo. After a Wick rotation the echo becomes the partition function of a CFT on a strip, and transverse contraction reads it numerically: the space–time tensor network is contracted along the spatial direction, which reduces the whole evolution to the leading eigenvalues and eigenvectors of a single transfer matrix. Because the quench is critical, the temporal entanglement grows only logarithmically, so the contraction stays efficient at times where conventional evolution is limited by the entanglement barrier.

This programme was established for integrable chains (Carignano & Tagliacozzo; Bou-Comas et al.). The thesis asks whether it survives the loss of integrability, using a self-dual ANNNI-type chain: the transverse-field Ising model with a next-nearest-neighbour coupling of strength `p`,

```
H = -Σ_i [ σᶻ_i σᶻ_{i+1} + λ σˣ_i + p (σᶻ_i σᶻ_{i+2} + λ σˣ_i σˣ_{i+1}) ]
```

which is genuinely interacting for any `p > 0` but remains critical and in the Ising universality class up to `p ≲ 1.5`. The equilibrium checks confirm this directly (the ground-state central charge stays at the Ising value across the whole range), and the dynamical measurements at the integrable point reproduce every conformal prediction. The frustrated couplings are being recomputed on a corrected construction of the transfer-matrix column; those results, from the cluster sweeps in `cluster/`, complete the story.

## What is in the repository

- `thesis/` — the LaTeX manuscript. The source of truth is edited in Overleaf and synced here as `BSC-TFMvN.zip`.
- `NBs/` — the work as a notebook series, in order: model and benchmarks, the time-evolution MPO, temporal entropies, the equilibrium DMRG checks and the sound velocity, the spectrum and central-charge measurements, and the validation of the corrected transfer-matrix column (`13_bulk_column`).
- `src/` — the Julia library everything uses (`include("src/thesislib.jl")`): model Hamiltonians, the temporal-MPO construction including the corrected bulk-column extraction, the block power method, and the entropy routines.
- `cluster/` — the SLURM jobs for MareNostrum, with their own README.
- `results/` — cached data and figures; every figure regenerates from the notebook that owns it.
- `other_models/` — exploratory notebooks (XXZ, tricritical), off the main line.
- `ITensorExpMPOv2.jl/` — a fork of [tipfom/ITensorExpMPO.jl](https://github.com/tipfom/ITensorExpMPO.jl); all upstream work is @tipfom's, and this thesis adds the second-order VD2 kernel (Van Damme et al.) so the NNN model evolves at genuine second order.

## Running it

```julia
julia --project=.
include("src/thesislib.jl")

mpo, scaffold = build_alcaraz_tmpo(4.0; p=0.1, nbeta=4, column=:bulk5)   # corrected transfer matrix
theta, L, R, info = block_transfer_eigs(mpo, scaffold; k=4)              # leading eigenpairs
mp = AlcarazParams(lambda=1.0, p=0.1)
compute_entropies(mp, 4.0; scheme=AlcarazVD2(), nbeta=4, column=:bulk5)  # Rényi-2 temporal entropy
```

`column=:bulk5` selects the corrected bulk-column extraction; the long sweeps cache to `results/data/` as they go, so an interrupted run resumes where it stopped.

## Credits and references

- [ITransverse.jl](https://github.com/starsfordummies/ITransverse.jl) (Stefano Carignano, BSC) — the transverse-contraction library this work builds on.
- [ITensors.jl](https://github.com/ITensor/ITensors.jl) — the tensor-network foundation.
- Carignano & Tagliacozzo, arXiv:2405.14706 — the framework and the integrable benchmark.
- Bou-Comas et al., arXiv:2607.08649 — conformal data from Loschmidt echoes, finite-time corrections.
- Van Damme, Haegeman, McCulloch & Vanderstraeten, SciPost Phys. 17, 135 (2024) — the VD2 construction.
- Alcaraz et al. — the ANNNI-type model and its Ising-class finite-size scaling.
