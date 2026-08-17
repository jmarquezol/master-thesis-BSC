# Cluster sweeps

Each job walks a ladder of evolution times T for one coupling p, saving after every rung. If it hits
the walltime, resubmit the same script and it carries on from where it stopped.

Three routes per coupling:

- `eigsweep` — central charge from the eigenvalue phase
- `entsweep` — temporal entropy
- `towerscan` — boundary dimensions

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have it: `git checkout main && git pull`. Do not update the packages —
`Manifest.toml` pins ITransverse to a specific commit.

## Submitting

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"

for p in 0.0 0.1 0.3 0.5 1.0; do sbatch submit_eigs_p${p}_bulk.slurm;  done
for p in 0.0 0.1 0.3 0.5 1.0; do sbatch submit_ent_p${p}_bulk.slurm;   done
for p in 0.0 0.1 0.3 0.5 1.0; do sbatch submit_tower_p${p}_bulk.slurm; done
```

The spectral arms matter most if the queue has to be cut short.

## Do not delete

- `results/data/cluster/` — the results, and what each job reads to skip work already done
- `cluster/checkpoints/` — what makes a resubmission resume warm

A resubmitted job reuses both, so deleting either means recomputing every rung from T=2.

## Notes

`dT` must be a multiple of the Trotter step: 0.5 at dt=0.1, 0.25 at dt=0.05. The driver refuses
anything else.

`WALL_RETRIES` (default 2) recomputes a rung from a fresh seed if its leading eigenvalue jumps away
from the previous one. `WALL_RETRIES=0` turns this off.

`WALL_BLAS_THREADS` caps the threads. Worth timing one rung at 4, 10 and 40 before filling the
queue: if it stops scaling early, several smaller jobs beat one large one.

Scripts under `archive/` are not to be run.

## Where things land

- results: `results/data/cluster/` — the folder to send back
- checkpoints: `cluster/checkpoints/`
- logs: `cluster/logs/`
