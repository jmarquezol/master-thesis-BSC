# Universal properties from quantum many-body dynamics

Master's thesis by Joaquín G. Márquez Olguín, supervised by Stefano Carignano at the Barcelona
Supercomputing Center, 2026.

## What the work is about

At a critical point, the real-time dynamics of a simple initial state carries the universal data of
the underlying conformal field theory: the central charge and the boundary operator content. We
reach it through the Loschmidt echo. After a Wick rotation the echo becomes the partition function
of a CFT on a strip, and transverse contraction evaluates it numerically: the space-time tensor
network is contracted along the spatial direction, so the evolution reduces to the leading
eigenvalues and eigenvectors of a single transfer matrix. Since the quench is critical, the temporal
entanglement grows only logarithmically, so the contraction remains efficient at times where
conventional evolution is already limited by the entanglement barrier.

This was established for integrable chains by Carignano and Tagliacozzo and by Bou-Comas et al. The
thesis asks whether it survives the loss of integrability. The model is a self-dual ANNNI-type
chain, the transverse-field Ising model with a next-nearest-neighbour coupling of strength `p`,

```
H = -Σ_i [ σᶻ_i σᶻ_{i+1} + λ σˣ_i + p (σᶻ_i σᶻ_{i+2} + λ σˣ_i σˣ_{i+1}) ]
```

which is interacting for any `p > 0` and stays critical and in the Ising universality class up to
`p ≈ 1.5`. The answer is that the signatures survive. The equilibrium checks place the model in the
Ising class over the whole range, and the dynamical measurements return the central charge through
two independent routes, together with seven members of the boundary operator spectrum, at every
coupling studied.

## How the repository is organised

| folder | what it holds |
|---|---|
| `thesis/` | the LaTeX manuscript and its figures, compiled with tectonic |
| `notebooks/` | six notebooks, in reading order, that walk through the whole project |
| `src/` | the Julia library everything else uses |
| `scripts/` | one script per figure, and the analysis that produced the numbers |
| `data/` | the cached results, so nothing long has to be rerun |
| `cluster/` | the submission scripts for the production runs |
| `defense/` | the slides for the defence |
| `figures/` | the generated plots, in PNG, PDF and SVG |

The notebooks are the place to start. They follow the structure of the thesis: the model and its
equilibrium properties, the method and its validation, the temporal entropies, the spectral route to
the central charge, the boundary operator spectrum, and the numerical controls. Each one computes at
least one result from scratch, at a small enough size to run in seconds, and then uses the cached
production data for the full sweeps, so the reader can see how a number is made without waiting for
it.

`src/thesislib.jl` loads the library: the model Hamiltonians, the temporal-MPO construction with the
corrected bulk-column extraction, the block power method, and the entropy routines.

`scripts/` has its own README explaining what every script does and whether it reads a cache or
produces one. `cluster/` likewise, for the runs that need more than a workstation.

`data/local/` holds everything computed on a workstation, including the seed ensembles under
`data/local/controls/`. `data/cluster/` holds the production sweeps.

## Reproducing the results

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # the Manifest pins ITransverse to the tested commit

julia --project=. scripts/analysis/cluster_audit.jl   # every number in the results table, with its provenance
julia --project=. scripts/analysis/extend_window.jl   # the fitting windows and the rungs chosen for them
julia --project=. scripts/figures/make_all.jl         # every thesis plot, as PNG, PDF and SVG
```

Then read the notebooks in order, with a Julia 1.12 kernel. They take a few minutes, since the long
calculations are cached.

Rebuilding a cache is a different matter, because those are the long runs. Each producer script says
in its header what it writes and roughly what it costs, checkpoints as it goes, and skips whatever
is already on disk, so an interrupted run resumes and a finished one costs nothing:

```bash
julia --project=. scripts/analysis/equilibrium_velocity.jl   # sound velocity, by exact diagonalisation
julia --project=. scripts/analysis/bulk_column_mu0.jl        # the transfer-matrix column test
julia --project=. scripts/analysis/svpm_ladder.jl p00        # one entropy arm; run without arguments for the list
julia --project=. scripts/analysis/mixedbc_ladder.jl upup    # one boundary pair
julia --project=. scripts/analysis/cutrerun.jl               # the cutoff control, over the whole grid
```

The production sweeps behind the main text were not run locally. They were submitted on MareNostrum
from `cluster/` and write into `data/cluster/`.

One convention applies everywhere. For a model with next-nearest-neighbour terms the
transfer-matrix column must be built from a five-site patch, `build_alcaraz_tmpo(...;
column=:bulk5)`. The three-site extraction is exact only for nearest-neighbour models and otherwise
drops a memory channel without any error being raised. Notebook 2 shows the difference.

## Credits and references

- [ITransverse.jl](https://github.com/starsfordummies/ITransverse.jl), by Stefano Carignano (BSC),
  is the transverse-contraction library this work builds on.
- [ITensors.jl](https://github.com/ITensor/ITensors.jl) is the tensor-network foundation.
- `ITensorExpMPOv2.jl/` is a fork of
  [tipfom/ITensorExpMPO.jl](https://github.com/tipfom/ITensorExpMPO.jl). All the upstream work is
  @tipfom's; this thesis adds the second-order VD2 kernel so that the NNN model evolves at genuine
  second order.
- Carignano and Tagliacozzo, arXiv:2405.14706, for the framework and the integrable benchmark.
- Bou-Comas et al., arXiv:2607.08649, for conformal data from Loschmidt echoes and the finite-time
  corrections.
- Van Damme, Haegeman, McCulloch and Vanderstraeten, SciPost Phys. 17, 135 (2024), for the VD2
  construction.
- Alcaraz et al., for the ANNNI-type model and its Ising-class finite-size scaling.
