# Universal properties from quantum many-body dynamics

Master's thesis — **Joaquín G. Márquez Olguín**, supervised by **Stefano Carignano** (Barcelona Supercomputing Center), 2026.

## What this is about

At a critical point, the real-time dynamics of a simple initial state carries the universal data of the underlying conformal field theory: the central charge and the boundary operator content. The route is the Loschmidt echo. After a Wick rotation the echo becomes the partition function of a CFT on a strip, and transverse contraction reads it numerically: the space–time tensor network is contracted along the spatial direction, which reduces the whole evolution to the leading eigenvalues and eigenvectors of a single transfer matrix. Because the quench is critical, the temporal entanglement grows only logarithmically, so the contraction stays efficient at times where conventional evolution is limited by the entanglement barrier.

This programme was established for integrable chains (Carignano & Tagliacozzo; Bou-Comas et al.). The thesis asks whether it survives the loss of integrability, using a self-dual ANNNI-type chain: the transverse-field Ising model with a next-nearest-neighbour coupling of strength `p`,

```
H = -Σ_i [ σᶻ_i σᶻ_{i+1} + λ σˣ_i + p (σᶻ_i σᶻ_{i+2} + λ σˣ_i σˣ_{i+1}) ]
```

which is genuinely interacting for any `p > 0` but remains critical and in the Ising universality class up to `p ≲ 1.5`. The answer is yes: the equilibrium checks place the model in the Ising class across the whole range, and the dynamical measurements — on a corrected construction of the transfer-matrix column (`column=:bulk5`), from the cluster sweeps in `cluster/` — return the central charge through two independent routes and seven members of the boundary operator spectrum at every coupling studied.

## What is in the repository

- `thesis/` — the LaTeX manuscript (compiled with tectonic; figures included as PDF).
- `notebooks/` — the guided tour, six notebooks in reading order: the model and its equilibrium properties, the method and its validation, the temporal entropies and the wall, the spectral route to the central charge, the boundary operator spectrum, and the numerical controls. Every cell loads shipped caches and calls library or script code — nothing heavy runs in a notebook.
- `src/` — the Julia library everything uses (`include("src/thesislib.jl")`): model Hamiltonians, the temporal-MPO construction including the corrected bulk-column extraction, the block power method, and the entropy routines.
- `scripts/figures/` — one script per thesis figure; `make_all.jl` regenerates every plot as PNG + PDF + SVG under `figures/` and syncs the PDFs into `thesis/imgs/`.
- `scripts/analysis/` — the numbers behind the thesis tables, split in two kinds. **Readers** take the shipped caches and print a result: `cluster_audit.jl` recomputes every value in the results table and shows where each one comes from, and `entropy_c.jl`, `eq3_windows.jl`, `chi_check.jl`, `deficit_tests.jl`, `mixedbc_analysis.jl`, `battery_report.jl` and `dtreport.jl` do the same for the individual controls. **Producers** are the runs that made those caches — the equilibrium ED and DMRG, the entropy and boundary ladders, the exponential-MPO benchmark, and the robustness controls (seed ensembles, cutoff scan, Trotter step, bond dimension, block size, warm starts). Every cache under `data/local/` has one, named in its header.
- `data/cluster/` — the production caches from the MareNostrum sweeps; `data/local/` — the local caches (equilibrium DMRG/ED, validation runs, controls).
- `figures/` — the generated figures (SVG is the Inkscape-editable version).
- `cluster/` — the SLURM jobs, with their own README.
- `presentation/` — the beamer draft for the defense.
- `ITensorExpMPOv2.jl/` — a fork of [tipfom/ITensorExpMPO.jl](https://github.com/tipfom/ITensorExpMPO.jl); all upstream work is @tipfom's, and this thesis adds the second-order VD2 kernel (Van Damme et al.) so the NNN model evolves at genuine second order.

## Reproducing the results

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # Manifest pins ITransverse to the tested commit

julia --project=. scripts/analysis/cluster_audit.jl   # every number in the results table, with its provenance
julia --project=. scripts/analysis/extend_window.jl   # the fitting windows, and the rungs selected for them
julia --project=. scripts/figures/make_all.jl         # every thesis plot, PNG + PDF + SVG
```

Then read the notebooks in order (Julia 1.12 kernel); they execute in minutes because everything they show is cached.

Regenerating a cache is a different matter — those are the long runs. Each producer says in its header what it writes and roughly what it costs, checkpoints as it goes, and skips whatever is already on disk, so an interrupted run resumes and a finished one costs nothing:

```bash
julia --project=. scripts/analysis/equilibrium_velocity.jl      # sound velocity, exact diagonalisation
julia --project=. scripts/analysis/bulk_column_mu0.jl           # the transfer-matrix column test
julia --project=. scripts/analysis/svpm_ladder.jl p00           # one entropy arm (see the script for the list)
julia --project=. scripts/analysis/mixedbc_ladder.jl upup       # one boundary pair
julia --project=. scripts/analysis/cutrerun.jl                  # the cutoff control, whole grid
```

The production sweeps behind the main text are not local runs at all: they are the SLURM jobs in `cluster/`, which write to `data/cluster/`.

One convention matters everywhere: for an NNN model the transfer-matrix column must be built from a five-site patch (`build_alcaraz_tmpo(...; column=:bulk5)`). The legacy three-site extraction silently drops a memory channel — notebook 2 demonstrates the difference.

## Credits and references

- [ITransverse.jl](https://github.com/starsfordummies/ITransverse.jl) (Stefano Carignano, BSC) — the transverse-contraction library this work builds on.
- [ITensors.jl](https://github.com/ITensor/ITensors.jl) — the tensor-network foundation.
- Carignano & Tagliacozzo, arXiv:2405.14706 — the framework and the integrable benchmark.
- Bou-Comas et al., arXiv:2607.08649 — conformal data from Loschmidt echoes, finite-time corrections.
- Van Damme, Haegeman, McCulloch & Vanderstraeten, SciPost Phys. 17, 135 (2024) — the VD2 construction.
- Alcaraz et al. — the ANNNI-type model and its Ising-class finite-size scaling.
