# Wall-scan sweeps on the cluster

This folder runs the transfer-matrix sweeps that feed the temporal-CFT analysis (NB7, NB9). Each
job is a single Julia process that walks a ladder of evolution times `T`, warm-starting each `T`
from the previous one and checkpointing as it goes, so a job that hits the walltime just needs to be
resubmitted — it resumes where it stopped.

The driver is `wall_scan_cluster.jl`; the `submit_*.slurm` files are thin wrappers around it.

## What to run now

Two things, in parallel:

**1. Eigenvalues-only sweeps, `p = 0, 0.1, 0.3, 0.5`, up to `T = 20`.** These skip the eigenvectors
(so no entropy), which is what lets them run far past the wall — we only need the spectrum here (dual
unitarity, the Eq. 3 central charge, the boundary exponent, the tower gaps).

```
sbatch submit_rtm_eigs_p0.0.slurm
sbatch submit_rtm_eigs_p0.1.slurm
sbatch submit_rtm_eigs_p0.3.slurm
sbatch submit_rtm_eigs_p0.5.slurm
```

**2. One RDM check.** A single full run (with entropy) at `p = 0.1`, up to `T = 12` (about where the
RTM `p = 0.1` run reached), using RDM truncation instead of RTM. The question is just whether the wall
moves when we switch truncation at a frustrated point — compare against `sweep_rtm_p0.1.jld2`. RDM is a
few times slower than RTM.

```
sbatch submit_rdm.slurm
```

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

Check progress with `squeue -u $USER` and `tail -f logs/eigs_p0.1-<jobid>.out`. Each finished `T`
logs a line like `[rtm_eigs_p0.1] T=9.0 converged@61 k=4 |θ0|=1.5443 gap=0.731 ...`.

If a job stops at the walltime, just `sbatch` the same script again — it resumes from the last
checkpoint. Expect the `T = 20` arms to need a couple of resubmissions.

Push results back whenever, it's safe at any point:
```
git add results/data/cluster/*.jld2 && git commit -m "cluster sweep progress" && git push
```

## What's already done

The full-eigenvector RTM runs (`p = 0 … 1.5`) finished earlier and are saved in `cluster/data/` and
`results/data/cluster/sweep_rtm_p*.jld2`. They stop at different `T` per `p` because the eigenvector
route hits the wall at different times (`p=0.1` at `T≈8`, `p=0.3` at `T≈7`, …) — that's the physics,
not a failed job. The eigenvalues-only reruns above are the follow-up that pushes past those walls.

## The other jobs

`submit_rtm.slurm` and `submit_cutoff.slurm` are the original `p=0.1` truncation comparison (bond
dimension / cutoff); `submit_rtm_p*.slurm` are the full-eigenvector p-sweep. Those are all run and
archived — leave them unless you want to regenerate something.

## Output format

`JLD2.load(path, "done")` gives a `Dict` keyed by `(label, T)`. Each entry is a named tuple with
`theta` (the spectrum), `theta_phys`/`i0` (the physical eigenvalue and its index), `tower_gap`,
`k_used`, and — for the full runs only — `s2_base` (the Rényi-2 dome) and `rigidity`.
