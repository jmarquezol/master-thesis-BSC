# Wall-scan sweeps on the cluster

Each job walks a ladder of evolution times T and saves after every rung. If a job hits the
walltime, resubmit the same script: it resumes where it stopped.

`wall_scan_cluster.jl` is the driver; the `submit_*.slurm` files are thin wrappers. Everything in
this folder is a job to run now. Finished or superseded ones are in `archive/`.

## What changed (August 14)

We found that the transfer-matrix column used so far was not the true bulk tensor of our model
once the NNN interaction is on: it was missing a memory channel. The corrected construction is
selected with `WALL_COLUMN=bulk5`, which the `*_bulk.slurm` scripts set themselves. Their results
go to separate `*_bulk.jld2` caches, so nothing mixes with earlier data.

Because of this, the old p != 0 scripts are superseded and live in `archive/legacy_column/`.
Please do not run them; if any are already in the queue, cancel them. The p = 0 job is unaffected
(both constructions are identical there).

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If you already have the repository: `git checkout main && git pull`. Do not update the packages —
`Manifest.toml` pins ITransverse to a specific commit and floating it has broken runs before.

## Submitting

From inside `cluster/`:

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK ... site dim 7"

sbatch submit_rtm_eigs_p0.0_fine.slurm    # highest priority, unchanged from before

sbatch submit_ent_p0.1_bulk.slurm         # entropy arms, corrected column
sbatch submit_ent_p0.3_bulk.slurm
sbatch submit_ent_p0.5_bulk.slurm

sbatch submit_ksec_p0.1_plus_bulk.slurm   # sector eigenvalue ladders, corrected column
sbatch submit_ksec_p0.1_minus_bulk.slurm
sbatch submit_ksec_p0.3_plus_bulk.slurm
sbatch submit_ksec_p0.3_minus_bulk.slurm
sbatch submit_ksec_p0.5_plus_bulk.slurm
sbatch submit_ksec_p0.5_minus_bulk.slurm

sbatch submit_tower_p0.1_bulk.slurm       # deep k=8 blocks for the tower figure
sbatch submit_tower_p0.5_bulk.slurm
```

The corrected column is more expensive (site dimension 13 instead of 7, roughly 2-4x per rung),
so these first ladders stop at T=12 (towers at T=8). To extend one later, raise the Tmax at the
bottom of its script and resubmit — it continues from its checkpoint.

## Where things land

- results: `results/data/cluster/sweep_<label>_bulk.jld2`  (this is the folder to send back)
- checkpoints: `cluster/checkpoints/` — never delete these; a cold restart corrupts the entropy
  ladders
- logs: `cluster/logs/`

If a job reports "nothing saved in data/", the results are in `results/data/cluster/`, not in a
top-level `data/` folder.
