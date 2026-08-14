# Wall-scan sweeps on the cluster

Each job walks a ladder of evolution times T and saves after every rung. If it hits the walltime,
resubmit the same script: it resumes from its checkpoint. To extend a finished ladder, raise the
Tmax at the bottom of its script and resubmit.

All jobs use the corrected transfer-matrix column (August 14 fix; the old column was missing a
memory channel at p != 0). Their results go to separate `*_bulk.jld2` caches. The superseded
scripts are in `archive/legacy_column/` — do not run them, and cancel any still in the queue.
The p = 0 job is unaffected by the fix.

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have the repository: `git checkout main && git pull`. Do not update the packages —
`Manifest.toml` pins ITransverse to a specific commit.

## Submitting

From inside `cluster/`, in this order:

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"

sbatch submit_rtm_eigs_p0.0_fine.slurm    # p=0 control, highest priority

sbatch submit_ent_p0.1_bulk.slurm         # entropy arms
sbatch submit_ent_p0.3_bulk.slurm
sbatch submit_ent_p0.5_bulk.slurm
sbatch submit_ent_p1.0_bulk.slurm
sbatch submit_ent_p1.5_bulk.slurm

sbatch submit_tower_p0.1_bulk.slurm       # k=8 blocks for the boundary-tower figure
sbatch submit_tower_p0.3_bulk.slurm
sbatch submit_tower_p0.5_bulk.slurm
sbatch submit_tower_p1.0_bulk.slurm
sbatch submit_tower_p1.5_bulk.slurm

sbatch submit_ent_p1.0_bulk_dt0.1.slurm   # Trotter controls: same arms at the standard dt=0.1
sbatch submit_ent_p1.5_bulk_dt0.1.slurm
sbatch submit_tower_p1.0_bulk_dt0.1.slurm
sbatch submit_tower_p1.5_bulk_dt0.1.slurm
```

The p = 1.0 and 1.5 arms use a finer Trotter step (set inside each script; the driver adjusts the
cooling to keep beta0 = 0.2). Their `_dt0.1` twins run the standard step on purpose, to measure
how much the step matters at these couplings. All eight are exploratory: nothing depends on them,
so if the queue is tight they go last.

## Where things land

- results: `results/data/cluster/` — this is the folder to send back
- checkpoints: `cluster/checkpoints/` — never delete these
- logs: `cluster/logs/`
