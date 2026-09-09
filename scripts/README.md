# Scripts

Two folders. `figures/` has one script per thesis plot. `analysis/` has everything that produces or
reads a cached result.

Every script runs from the repository root:

```
julia --project=. scripts/analysis/cluster_audit.jl
```

## How to read the analysis folder

Scripts come in two kinds, and the difference matters if you are deciding what to run.

**Readers** take the caches that ship with the repository and print a result. They are fast, they
change nothing, and they are the ones to run first.

**Producers** are the calculations that made those caches. They are slow. Each one checkpoints as it
goes and skips whatever is already on disk, so an interrupted run resumes and a finished one costs
nothing. You only need them if you want to rebuild a cache from scratch.

Some scripts do both: they compute their cache the first time and print the report every time.

## The results table

These are the ones behind Table 1 and the central-charge numbers.

| script | kind | what it does |
|---|---|---|
| `cluster_audit.jl` | reader | recomputes every number in the results table from the cluster caches alone, and prints where each came from |
| `table1_errors.jl` | reader | the uncertainties: seed resampling for the entropy route, window variants for both |
| `extend_window.jl` | reader | chooses the Eq. (3) fitting windows, and extends them past the reach of the block iteration |
| `eq3_windows.jl` | reader | window sensitivity of the spectral fits, and the branch constant free against pinned |
| `entropy_c.jl` | reader | the entropy route: per-rung chord fits, the regulator terms, and the model comparison behind the plateau extrapolation |

## Equilibrium inputs

| script | kind | what it does |
|---|---|---|
| `equilibrium_velocity.jl` | producer | sound velocity from exact diagonalisation on rings of 10 to 18 sites, extrapolated in 1/N² |
| `equilibrium_dmrg.jl` | producer | the three DMRG caches behind the equilibrium central-charge figures |

## The evolution operator and the transfer-matrix column

| script | kind | what it does |
|---|---|---|
| `mpo_order_check.jl` | producer | norm drift and convergence order of W¹, W² and VD2 — the two MPO comparison tables |
| `mpo_dt_tail.jl` | producer | the same order measurement at smaller time steps, where the exponent settles |
| `rate_benchmark.jl` | producer | the Loschmidt rate from repeated MPO application against a TDVP reference |
| `bulk_column_mu0.jl` | producer | the three-site against five-site column, both measured against exact Krylov evolution |
| `blockpm_validation.jl` | producer | the block method against dense exact diagonalisation at short times |

## Ladders: the runs themselves

| script | kind | what it does |
|---|---|---|
| `svpm_ladder.jl` | producer | the local single-vector entropy arms used by the appendices |
| `seedens.jl` | producer | seed ensembles: the same rung repeated from independent seeds. These are what the quoted entropy uncertainty is computed from |
| `mixedbc_ladder.jl` | producer | the boundary-pair ladders at p=0, for the four pairings of free and fixed |
| `mixedbc_analysis.jl` | reader | the four pairs: their towers, and the absolute tower position from the leading phase |
| `mixedbc_consensus.jl` | producer | agreement between two independent fixed-boundary ladders, gap by gap rather than on the modulus alone |
| `ising_x1.jl` | producer | the boundary exponent at the exactly solvable point, the benchmark of notebook 2 |
| `gap_ladder.jl` | producer | the coarse map of where the window ends, at four couplings |

## Controls

Every knob that could have set the answer, turned and measured. Notebook 6 runs all of these.

| script | kind | what it varies |
|---|---|---|
| `battery_report.jl` | reader | prints the seed spreads, the cutoff scan and the dense diagnostics together |
| `chi_check.jl` | reader | bond dimension 64 against 128, on both estimators |
| `seedens_chi.jl` | producer | whether a doubled cap recovers rungs that failed at 64 |
| `dtreport.jl` | reader | Trotter step 0.05 against the 0.1 baseline |
| `trotter_fine.jl` | producer | a third time step, so the plateau can be extrapolated in δt² from three points |
| `cutrerun.jl` | producer | the singular-value cutoff, on the block path |
| `blocksize_k6.jl` | producer | six states in the block instead of four |
| `warmcold.jl` | both | warm start against cold start, same seed |
| `mode_agreement.jl` | both | the two eigensolver bases against each other |
| `costcheck.jl` | both | whether compressing the partial direct sums saves anything |
| `rdmswap.jl` | producer | density-matrix truncation instead of transition-matrix, same seeds |
| `conv_history.jl` | producer | how the block iteration converges, or stops without converging |
| `conv_history_sv.jl` | producer | the same for the single-vector iteration |
| `spectral_seeds.jl` | producer | failure statistics of the block route near the window edge |
| `dense_bulk5.jl` | producer | small transfer matrices built densely, for the conditioning argument |
| `deficit_tests.jl` | reader | candidate origins of the deficit in the fitted Eq. (3) coefficient |
| `analysis_pack.jl` | reader | the remaining catalogue verdicts, over existing caches |

## Figures

`figures/` holds one script per thesis plot. Each defines `make_<name>()` returning the plot, so a
notebook can `include` the file and call the function; run standalone, it writes PNG, PDF and SVG
into `figures/` and syncs the PDF into `thesis/imgs/`.

```
julia --project=. scripts/figures/make_all.jl        # every plot at once
```

`fig_cft_L.jl`, `fig_chord_equilibrium.jl`, `fig_velocity_vs_p.jl` and `fig_velocity_extrap.jl` are
the equilibrium figures; `fig_domes.jl`, `fig_domes_hi.jl` and `fig_re_c.jl` the temporal entropies;
`fig_spectral.jl`, `fig_beta0.jl`, `fig_tower.jl` and `fig_bc_pairs.jl` the spectral route and the
boundary spectrum; `fig_conv.jl` the convergence control.

## Where the data lives

`data/local/` holds everything computed on a workstation, with the seed ensembles under
`data/local/controls/`. `data/cluster/` holds the production sweeps; those are submitted from
`cluster/`, which has its own README.
