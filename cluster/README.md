# Wall-scan sweeps on the cluster

This folder runs the transfer-matrix sweeps that feed the temporal-CFT analysis (NB7, NB9). Each
job is a single Julia process that walks a ladder of evolution times `T`, warm-starting each `T`
from the previous one and checkpointing as it goes, so a job that hits the walltime just needs to be
resubmitted — it resumes where it stopped.

The driver is `wall_scan_cluster.jl`; the `submit_*.slurm` files are thin wrappers around it.

## What to run now

The August 2026 round finished four eigenvalue ladders and two entropy ones. What is left is
finishing the ladders that ran out of walltime, plus a new kind of job that did not exist before.

**1. Resume the eigenvalues-only sweeps.** `p = 0` is complete (`T = 20`). The other three stopped
short: `p = 0.1` at `T = 17`, `p = 0.3` at `T = 13`, `p = 0.5` at `T = 15`. Resubmitting the same
script resumes from the cached rungs. These are the slowest jobs we have — a single rung at
`p = 0.3, T = 14` ran over seven hours on a 16-core desktop without finishing — so expect a couple
of resubmissions each.

```
sbatch submit_rtm_eigs_p0.1.slurm
sbatch submit_rtm_eigs_p0.3.slurm
sbatch submit_rtm_eigs_p0.5.slurm
```

**2. Half-integer ladders at `p = 0.3` and `0.5`.** New. We follow the physical eigenvalue by
predicting the phase it should have at the next rung. That advance grows with the sound velocity,
and at `p ≥ 0.3` a ladder spaced `ΔT = 1` advances by more than `π` between rungs — ambiguous
modulo `2π`, so the branch is lost after a few rungs and the extraction stops working. Halving the
spacing resolves the winding. The driver now takes an optional third argument for the step, and the
new rungs merge into the same cache, so the analysis just sees one denser ladder.

```
sbatch submit_rtm_eigs_p0.3_fine.slurm
sbatch submit_rtm_eigs_p0.5_fine.slurm
```

**3. The entropy arms.** `psweep 0.1` stopped at `T = 8` and `psweep 0.5` at `T = 6`. Both compute
the temporal entropy as well as the spectrum, so they are heavier per rung, and past the wall
(`T ≈ 9` at `p = 0.1`, `T ≈ 2` at `p = 0.5`) the entropy stops meaning anything. Worth having as a
record of where the wall sits, but lower priority than 1 and 2. These are the jobs that genuinely
need their checkpoint — a cold restart distorts the dome.

```
sbatch submit_rtm_p0.1.slurm
sbatch submit_rtm_p0.5.slurm
```

**Do not resubmit the RDM job.** It reached `T = 12` and answered its question: the dome breaks at
`T ≈ 10` under RDM truncation, the same place as RTM, and the four rigidities are all ≈ 5·10⁻⁶
there, so no eigenvector in the block survives. Extending it would buy more rungs on the wrong side
of the wall. See NB9 §3.

The β₀ jobs (`submit_beta_*.slurm`, `submit_betawall_*.slurm`) are done and answered — leave them.

Everything writes to `results/data/cluster/` under its own filename, so nothing overwrites anything.

## Running it

```
git checkout main && git pull
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # from the repo root, once
# then sbatch the jobs above from inside cluster/
```

Every submit script now runs `wall_scan_cluster.jl preflight` first. That builds one small tMPO and
exits, which is the same call the long run makes on every rung. It exists because of how the first
round ended: four jobs hit the 48 h walltime, and on requeue two of them could no longer dispatch
`build_alcaraz_tmpo` — a `MethodError` on the VD2 kernel from a stale precompile. The ladder then
marked every remaining `T` as errored in a few seconds. The preflight turns that into a one-minute
failure, and the driver now also stops after two consecutive failures instead of burning the rest of
the ladder. Run `Pkg.instantiate()` again if the preflight fails, and check that `Manifest.toml`
still pins ITransverse at `f10aee05`.

## Checkpoints

Checkpoints live in `cluster/checkpoints/checkpoint_<label>.jld2` and are gitignored, so they never
travel with a `git push` — they have to be copied across by hand. They are what makes a resubmission
resume *warm*; without one the first new rung starts cold, and a cold restart is what corrupted the
entropy dome in the first array sweep. Never delete them.

All the checkpoints we have are shipped as a zip alongside the repo. Unpack it so the files land in
`cluster/checkpoints/`, then check they are there before submitting:

```
unzip checkpoints.zip -d cluster/
ls cluster/checkpoints/          # expect checkpoint_rtm_eigs_p0.{1,3,5}.jld2 and checkpoint_rtm_p0.{3,5}.jld2
```

The eigenvalue jobs run correctly without them — they compute no entropy, so a cold first rung only
costs iterations. The entropy jobs in step 3 do need theirs.

Check progress with `squeue -u $USER` and `tail -f logs/eigs_p0.1-<jobid>.out`. Each finished `T`
logs a line like `[rtm_eigs_p0.1] T=9.0 converged@61 k=4 |θ0|=1.5443 gap=0.731 ...`.

If a job stops at the walltime, just `sbatch` the same script again — it resumes from the last
checkpoint. Expect the `T = 20` arms to need a couple of resubmissions.

Push results back whenever, it's safe at any point:
```
git add results/data/cluster/*.jld2 && git commit -m "cluster sweep progress" && git push
```

## What's already done

The 2026-07 full-eigenvector RTM runs (`p = 0 … 1.5`) stop at different `T` per `p` because the
eigenvector route hits the wall at different times — that is the physics, not a failed job. The ones
that have since been superseded or fall outside the `0 ≤ p ≤ 0.5` thesis range moved to
`results/data/cluster/archive/`; the duplicate copy under `cluster/data/` is gone, and the raw
per-`T` array-job outputs are in `cluster/archive/` (NB9 §1 still merges them into `cold_sweep.jld2`).

The August round delivered `sweep_rtm_eigs_p{0.0,0.1,0.3,0.5}.jld2`, `sweep_rdm_p0.1.jld2`, and the
`p = 0.3` entropy arm. That last one arrived as `warm_sweep.jld2`, because `psweep` used to write
every `p` into the driver's shared default cache; it now writes `sweep_rtm_p<p>.jld2` like
`eigsweep` does, and the file was renamed to match.

## The other jobs

`submit_rtm.slurm` and `submit_cutoff.slurm` are the original `p=0.1` truncation comparison (bond
dimension / cutoff); `submit_rtm_p*.slurm` are the full-eigenvector p-sweep. Those are all run and
archived — leave them unless you want to regenerate something.

## Output format

`JLD2.load(path, "done")` gives a `Dict` keyed by `(label, T)`. Each entry is a named tuple with
`theta` (the spectrum), `theta_phys`/`i0` (the physical eigenvalue and its index), `tower_gap`,
`k_used`, and — for the full runs only — `s2_base` (the Rényi-2 dome) and `rigidity`.
