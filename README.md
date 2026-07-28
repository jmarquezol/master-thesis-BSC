# Universal properties from quantum many-body dynamics

Master's thesis — **Joaquín G. Márquez Olguín**, supervised by **Stefano Carignano** (Barcelona Supercomputing Center).

## What's this about?

The short version: you can read the *universal* fingerprints of a critical quantum chain — things like its central charge and its spectrum of scaling dimensions — straight out of its **real-time dynamics**, without ever computing a ground state.

The idea is an old duality dressed up for tensor networks. Take a simple product state, quench it to a critical point, and look at the Loschmidt echo (the amplitude to return to where you started). Rotate real time into imaginary time and that echo becomes the partition function of a 2D conformal field theory on a strip — so all the universal CFT data is sitting right there inside the dynamics. **Transverse contraction** is what gets at it numerically: we treat time as if it were a spatial direction, build the 2D space-time network, and contract it *sideways*, along the original space axis. The whole time evolution then reduces to the eigenvalues of a single transfer matrix, and because the quench is to a critical point the temporal entanglement grows only logarithmically — so this stays cheap exactly where ordinary time evolution would choke on the entanglement barrier.

All of this has been worked out for **integrable** chains (the Ising chain, in the papers by Carignano & Tagliacozzo and by Bou-Comas et al.). The question this thesis asks is whether any of it survives once you break integrability. We test it on an **ANNNI-type chain** — the transverse-field Ising model with an extra next-nearest-neighbour coupling of strength `p`:

```
H = -Σ_i [ σᶻ_i σᶻ_{i+1} + p σᶻ_i σᶻ_{i+2} + λ σˣ_i + p λ σˣ_i σˣ_{i+1} ]
```

It's self-dual, genuinely interacting for any `p > 0`, but it stays critical and in the Ising universality class. The answer turns out to be yes: the temporal central charge stays Ising (`c(p=0.1) ≈ 0.47 ± 0.05`), the boundary operator content comes out right — including a nice surprise, an odd-parity tower that the interaction pushes out to momentum `π` — and we map out where and why the method eventually runs into a wall.

## What's in here

- `thesis/` — the LaTeX manuscript.
- `NBs/` — the story told in notebooks, roughly in order: building and benchmarking the time-evolution MPO, the temporal entropies and the block power method that computes them, the equilibrium DMRG check that the model really is Ising-class, and then the main result (central charge and operator content as a function of `p`).
- `src/` — the Julia library everything sits on (`include("src/thesislib.jl")`): the model Hamiltonians, the temporal-MPO construction, the block power method, and the entropy routines.
- `cluster/` — SLURM scripts for the heavier sweeps that run on the BSC cluster.
- `other_models/` — exploratory notebooks on the XXZ chain and the tricritical point, kept off to the side from the main line.
- `results/` — cached data and figures, each regenerable by the notebook that produced it.

## Running it

Activate the project environment and load the library:

```julia
julia --project=.
include("src/thesislib.jl")
```

and the main drivers are one call away:

```julia
mp = AlcarazParams(lambda=1.0, p=0.1)
compute_entropies(mp, 4.0; scheme=AlcarazVD2(), nbeta=4)      # Rényi-2 temporal entropy
mpo, scaffold = build_alcaraz_tmpo(4.0; p=0.1, nbeta=4)       # the temporal MPO
theta, L, R, info = block_transfer_eigs(mpo, scaffold; k=4)   # leading transfer eigenvalues
```

The long runs cache to `results/data/` as they go, so an interrupted sweep just picks up where it left off.

## Credits

- **[ITransverse.jl](https://github.com/starsfordummies/ITransverse.jl)** (Stefano Carignano, BSC) — the transverse-contraction library; every transverse algorithm used here comes from it.
- **ITensorExpMPO** (`ITensorExpMPOv2.jl/`) — a fork of [tipfom/ITensorExpMPO.jl](https://github.com/tipfom/ITensorExpMPO.jl). All of it is [@tipfom](https://github.com/tipfom)'s work except the second-order VD2 kernel added for this thesis (the Van Damme *et al.* construction, so the next-nearest-neighbour model evolves at genuine 2nd order).
- Built on top of **[ITensors.jl](https://github.com/ITensor/ITensors.jl)**.

## References

- Carignano & Tagliacozzo, *Loschmidt echo, emerging dual unitarity and scaling of generalized temporal entropies after quenches to the critical point*, arXiv:2405.14706.
- Bou-Comas, Carignano, Cerezo-Roquebrún, Lopez & Tagliacozzo, *Extracting conformal data from Loschmidt echoes after critical quenches*, arXiv:2607.08649.
- Van Damme, Haegeman, McCulloch & Vanderstraeten, *Efficient higher-order matrix product operators for time evolution*, SciPost Phys. **17**, 135 (2024) — the VD2 construction.
- Alcaraz — the original ANNNI-type model and its Ising-class finite-size scaling.
