# Wall-scan sweeps on the cluster

These jobs run the transfer-matrix sweeps behind the temporal-CFT analysis. Each one walks a ladder
of evolution times `T`, warm-starting every rung from the previous one and checkpointing as it goes.
If a job hits the walltime, resubmit the same script — it resumes.

`wall_scan_cluster.jl` is the driver; the `submit_*.slurm` files are thin wrappers. Everything in
this folder still needs computation; finished jobs are in `done/`.

## Setup

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
unzip checkpoints_for_bsc.zip -d cluster/
```

Do not update the packages — `Manifest.toml` pins ITransverse to a specific commit and floating it
has broken runs before. If you already have the repository, `git checkout main && git pull`.

Then submit from inside `cluster/`:

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"
sbatch submit_rtm_eigs_p0.1.slurm                    # and the rest
```

Every submit script runs the preflight first, so a broken environment fails in a minute instead of
erroring every rung. If it fails, run `Pkg.instantiate()` again and check that `Manifest.toml` still
pins ITransverse at `f10aee05`.

## Checkpoints

`checkpoints_for_bsc.zip` holds three checkpoints, one per job that resumes:

| checkpoint | job it resumes | at |
|---|---|---|
| `checkpoint_rtm_eigs_p0.1.jld2` | `submit_rtm_eigs_p0.1` (and `_fine`) | `T = 17` |
| `checkpoint_rtm_eigs_p0.3.jld2` | `submit_rtm_eigs_p0.3` (and `_fine`) | `T = 13` |
| `checkpoint_rtm_eigs_p0.5.jld2` | `submit_rtm_eigs_p0.5` (and `_fine`) | `T = 15` |

Unpack so they land in `cluster/checkpoints/`, and check before submitting:

```
ls cluster/checkpoints/
```

Never delete them. The other jobs — the entropy sweeps and the tower scan — start at `T = 2` and
write their own.

## The jobs

**1. Eigenvalue sweeps.** The main arms, and the slowest. `p = 0` is already complete to `T = 20`;
these three resume from their checkpoints and target the same.

```
sbatch submit_rtm_eigs_p0.1.slurm      # from T=17
sbatch submit_rtm_eigs_p0.3.slurm      # from T=13
sbatch submit_rtm_eigs_p0.5.slurm      # from T=15
```

**2. Half-integer ladders.** The physical eigenvalue is followed by predicting its phase and taking
the nearest member of the block. A larger sound velocity packs the members closer, and the
prediction error grows with the step, so at `p ≥ 0.3` the branch is lost after a few rungs. Halving
the step fixes both. New rungs merge into the same cache, so the analysis sees one denser ladder.
Worth doing at `p = 0.1` too, since that is where the results are quoted.

```
sbatch submit_rtm_eigs_p0.1_fine.slurm
sbatch submit_rtm_eigs_p0.3_fine.slurm
sbatch submit_rtm_eigs_p0.5_fine.slurm
```

**3. Entropy sweeps.** New runs, cold from `T = 2`. They compute the temporal entropy as well as the
spectrum, so they cost more per rung — `p = 0` has no wall and its later rungs take hours each. At
`p ≥ 0.3` the entropy stops meaning anything a few rungs in, so those arms mainly record the wall.

They replace an earlier set that stored only one entropy profile per rung, for whichever eigenvalue
the driver picked there. That pick can only be judged with the whole ladder in hand, so it could not
be corrected afterwards. These store one profile per block member and leave the choice to the
analysis.

```
sbatch submit_ent_p0.0.slurm
sbatch submit_ent_p0.1.slurm
sbatch submit_ent_p0.3.slurm
sbatch submit_ent_p0.5.slurm
```

**4. Deep block.** The sweeps above keep four eigenvalues per rung — enough to follow the physical
branch, not enough to show the tower around it. This keeps eight, at small `T`, so it is short.

```
sbatch submit_tower_p0.1.slurm
```

Everything writes to `results/data/cluster/` under its own filename, so nothing overwrites anything.

## Watching and pushing

```
squeue -u $USER
tail -f logs/eigs_p0.1-<jobid>.out
```

Each finished `T` logs a line like
`[rtm_eigs_p0.1] T=9.0 converged@61 k=4 |θ0|=1.5443 gap=0.731 ...`.

Push results whenever — it is safe at any point:

```
git add results/data/cluster/*.jld2 && git commit -m "cluster sweep progress" && git push
```

## What's in `done/`

Kept for reference. `submit_rtm.slurm`, `submit_cutoff.slurm` and `submit_rdm.slurm` are the
`p = 0.1` truncation comparison, which answered its question. `submit_rtm_p*.slurm` are the earlier
entropy arms, replaced by the `submit_ent_*` jobs. `submit_rtm_p1.0/1.5.slurm` fall outside the
`0 ≤ p ≤ 0.5` range of the thesis. `submit_rtm_eigs_p0.0.slurm` is the finished `p = 0` ladder, and
`submit_beta_*.slurm` the regulator scans. `submit_betawall_*.slurm` asked whether a larger β₀ pushes
the wall out; the regulator scan already answers it, so they were not run — the rigidity collapses at
the same rate for every β₀ between 0.1 and 0.8. Superseded data is in
`results/data/cluster/archive/`.

## Output format

`JLD2.load(path, "done")` gives a `Dict` keyed by `(label, T)`. Each entry is a named tuple with
`theta` (the spectrum), `theta_phys`/`i0` (the physical eigenvalue and its index), `tower_gap`,
`k_used`, and — for the entropy runs — `rigidity`, `s2_all` (one Rényi-2 profile per block member)
and `s2_base`, the one at `i0`.
