# Wall-scan sweeps on the cluster

These jobs run the transfer-matrix sweeps behind the temporal-CFT analysis. Each one walks a
ladder of evolution times `T` and saves after every rung, so if a job hits the walltime you just
resubmit the same script and it resumes where it stopped.

`wall_scan_cluster.jl` is the driver; the `submit_*.slurm` files are thin wrappers. Everything in
this folder is a job to run now. Finished ones are in `archive/done/`.

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have the repository: `git checkout main && git pull`. Do not update the packages —
`Manifest.toml` pins ITransverse to a specific commit and floating it has broken runs before.

No checkpoint files are needed this round: every job starts from `T = 2` and writes its own
checkpoint as it goes.

## Submitting

From inside `cluster/`:

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"

sbatch submit_rtm_eigs_p0.0_fine.slurm   # 1. the most important one
sbatch submit_ksec_p0.3_plus.slurm       # 2.
sbatch submit_ksec_p0.3_minus.slurm      # 3.
sbatch submit_ksec_p0.5_plus.slurm       # 4.
sbatch submit_ksec_p0.5_minus.slurm      # 5.
sbatch submit_tower_p0.1.slurm           # 6.
sbatch submit_tower_p0.3.slurm           # 7.
sbatch submit_tower_p0.5.slurm           # 8.
```

Every script runs the preflight itself, so a broken environment fails in a minute instead of
erroring every rung. All jobs write their own cache, so they can all run at the same time.

## What the jobs are

**`rtm_eigs_p0.0_fine`** — half-integer eigenvalue rungs at the integrable point, to `T = 18`.
The p=0.1 fine and unit ladders disagree by twenty per cent on the central charge, and this is the
control at the coupling where the answer is known. The most valuable job in the queue.

**`ksec_*`** — eigenvalue ladders confined to one symmetry sector each, to `T = 18`. The transfer
matrix has an exact Z2 symmetry whose two sectors hold the two eigenvalue families, and the
branch-tracking failures at `p ≥ 0.3` were the two families crossing in modulus inside one run.
Confining the run to one sector fixes the branch before the run starts. Each rung prints a sector
charge that should read ±1.000. A rung reporting `stuck` past `T ≈ 4` is normal — the eigenvalue
is still accurate there — and a resubmission resumes from the checkpoint.

**`tower_*`** — k=8 blocks at small `T`, for the tower figure. Short jobs.

## Where the results go

Results land in `results/data/cluster/sweep_<label>.jld2` — not in `data/`, which has confused
people before. Push them whenever; it is safe at any point:

```
git add results/data/cluster/*.jld2 && git commit -m "cluster sweep progress" && git push
```

## Watching

```
squeue -u $USER
tail -f logs/ksec_p0.3_plus-<jobid>.out
```

Each finished rung logs one line with the time, the leading eigenvalue and the timing.

## Already run

The entropy arms (`sweep_ent_p*`), the mixed eigenvalue arms (`sweep_rtm_eigs_p*` and their
half-integer `_fine` versions at p = 0.1–0.5), and the beta scans are done; their data is in
`results/data/cluster/` and their submit scripts in `archive/done/`. The mixed eigenvalue arms at
`p ≥ 0.3` are superseded by the `ksec` ladders, which resolve the same spectrum separated by
family.
