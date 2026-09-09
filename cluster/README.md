# Cluster sweeps

The production runs behind the results chapter. Each job walks a ladder of evolution times `T` for
one coupling `p` and saves after every rung, so if it hits the walltime you resubmit the same script
and it carries on from where it stopped.

Three kinds of run, all driven by `wall_scan_cluster.jl`:

| mode | method | what it gives |
|---|---|---|
| `eigsweep` | block iteration, eigenvalues only | the central charge, from the phase of the leading eigenvalue |
| `towerscan` | block iteration, `k=8` | the boundary operator dimensions |
| `entsweep` | single-vector power method | the temporal entropies, and the leading eigenvalue where the `eigsweep` arm does not reach |

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have it, `git checkout main && git pull`. Do not update the packages: `Manifest.toml`
pins ITransverse to the commit everything was run with.

## Submitting

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"
for f in block/*.slurm single_vector/*.slurm; do sbatch $f; done
```

The scripts are split by method: `block/` holds the `eigsweep` and `towerscan` arms,
`single_vector/` the `entsweep` ones. Submit from `cluster/` as above, or from inside either
subdirectory — each script steps back up to `cluster/` before doing anything, so both work.

Every arm is split into an `_a` and a `_b` job that advance the same ladder half a rung apart. The
`_a` job continues the existing chain. The `_b` job first calls `fork`, which seeds a second chain
from the same checkpoint, and then advances on the offset grid; that call does nothing if the chain
already exists, so resubmitting is safe. The two grids together give the `dT=0.5` ladder, and the
analysis merges arms by label, so nothing downstream needs to know about the split.

If the queue is tight, the spectral arms are the ones to protect: `eigs_p0.0`, then `eigs_p0.1`,
then the rest.

## What must not be deleted

- `data/cluster/` — the results, and what each job reads to skip rungs already done
- `cluster/checkpoints/` — what lets a resubmission resume warm

A resubmitted job uses both, so losing either means recomputing every rung from `T=2`. Neither is in
the repository: the caches are shipped under `data/cluster/`, the checkpoints are local state and are
gitignored.

## Environment

`WALL_COLUMN` selects the transfer-matrix column: `bulk5` builds the evolution operator on five sites
and takes the true middle tensor, which is what an NNN model needs. It is set in every submission script except
the two `eigsweep` ones at `p=0`, where the chain is nearest-neighbour only and the three-site
extraction is already exact. The driver appends `_bulk` to the cache label when `bulk5` is on, so the
two conventions can never overwrite each other.

`dT` must be a multiple of the Trotter step, 0.5 at `dt=0.1` and 0.25 at `dt=0.05`. The driver
refuses anything else.

`WALL_RETRIES` (default 2) recomputes a rung from a fresh seed if its leading eigenvalue jumps away
from the previous one. `WALL_RETRIES=0` turns that off.

`WALL_BLAS_THREADS` caps the thread count. Time one rung at 4, 10 and 40 before filling the queue: if
it stops scaling early, several smaller jobs beat one large one.

## Layout

```
cluster/
  wall_scan_cluster.jl        the driver, all modes
  block/                      eigsweep (8 files) and towerscan (4)
  single_vector/              entsweep (8 files)
  checkpoints/                resume state, not in the repository
  logs/                       job output, not in the repository
```

The two tower scripts at `p=0` and `p=0.1` are reconstructions: those arms ran in an earlier
round and their submission files were not kept. Each says so in its header, and the arguments
reproduce the shipped caches.

## Where the output lands

- results: `data/cluster/` — the folder to copy back
- checkpoints: `cluster/checkpoints/`
- job logs: `cluster/logs/`
