# Wall-scan sweeps on the cluster

This folder runs the transfer-matrix sweeps that feed the temporal-CFT analysis (NB7, NB9). Each
job is a single Julia process that walks a ladder of evolution times `T`, warm-starting each `T`
from the previous one and checkpointing as it goes, so a job that hits the walltime just needs to be
resubmitted — it resumes where it stopped.

The driver is `wall_scan_cluster.jl`; the `submit_*.slurm` files are thin wrappers around it.

## What to run now

The first round finished in August 2026. What is left is finishing the ladders that ran out of
walltime, plus one arm that was never submitted.

**1. Resume the eigenvalues-only sweeps.** `p = 0` is complete (`T = 20`). The other three stopped
short: `p = 0.1` at `T = 17`, `p = 0.3` at `T = 13`, `p = 0.5` at `T = 15`. Resubmitting the same
script resumes from the cached rungs.

```
sbatch submit_rtm_eigs_p0.1.slurm
sbatch submit_rtm_eigs_p0.3.slurm
sbatch submit_rtm_eigs_p0.5.slurm
```

**2. The two entropy arms.** `psweep 0.3` reached `T = 14` and resumes; `psweep 0.5` has never run
and is the only source for the `p = 0.5` entropy dome, since the eigsweeps carry no entropy and the
2026-07 `sweep_rtm_p0.5` run is unusable.

```
sbatch submit_rtm_p0.3.slurm
sbatch submit_rtm_p0.5.slurm
```

**Do not resubmit the RDM job.** It reached `T = 12` and answered its question: the dome breaks at
`T ≈ 10` under RDM truncation, the same place as RTM, and the four rigidities are all ≈ 5·10⁻⁶
there, so no eigenvector in the block survives. Extending it would buy more rungs on the wrong side
of the wall. See NB9 §3.

**3. β₀ regulator scan, `p = 0` and `0.1`.** Each job reruns the same full sweep for a few amounts of
imaginary-time cooling (`nbeta = 2,4,…,16`, i.e. `β₀ = 0.1 … 0.8`), up to `T = 10` — to check how
much the extracted central charge and boundary exponent depend on that regulator. Short jobs.

```
sbatch submit_beta_p0.0.slurm
sbatch submit_beta_p0.1.slurm
```

**4. Does β₀ move the wall?** The β₀ scan above showed the modulus gaps between the leading
eigenvalues grow *linearly* with β₀ — and it's those gaps closing that ends the eigenvector route
(the "wall"). So a bigger β₀ should, in principle, buy more reach. These two jobs test it directly:
the same long ladder (`T` up to 14) run at `nbeta = 12` (β₀ = 0.6) and at our usual `nbeta = 4`
(β₀ = 0.2) as the control, so the only difference between them is the regulator. Compare where the
entropy dome breaks in each. Worth knowing: the gap enhancement falls off as β₀/T², so the honest
expectation is a modest shift rather than an escape — a null result is still a useful answer.

```
sbatch submit_betawall_nb12.slurm
sbatch submit_betawall_nb4.slurm
```

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

Present: `rtm_eigs_p0.1` (`T = 17`), `rtm_eigs_p0.3` (`T = 13`), `betawall_p0.1_nb4`.
Missing: `rtm_eigs_p0.5` (`T = 15`), `rtm_p0.3` (`T = 14`), `rdm_p0.1` (`T = 12`).

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
