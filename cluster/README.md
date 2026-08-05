# Wall-scan sweeps on the cluster

This folder runs the transfer-matrix sweeps that feed the temporal-CFT analysis (NB7, NB9). Each
job is a single Julia process that walks a ladder of evolution times `T`, warm-starting each `T`
from the previous one and checkpointing as it goes, so a job that hits the walltime just needs to be
resubmitted — it resumes where it stopped.

The driver is `wall_scan_cluster.jl`; the `submit_*.slurm` files are thin wrappers around it. Every
script sitting in this folder still needs computation. The ones that have already run are in
`done/`.

## What to run now

**1. Eigenvalue-only sweeps.** `p = 0` is complete up to `T = 20`. The others stopped short:
`p = 0.1` at `T = 17`, `p = 0.3` at `T = 13`, `p = 0.5` at `T = 15`. Resubmitting resumes from the
cached rungs. These are the slowest jobs, so expect a few resubmissions each.

```
sbatch submit_rtm_eigs_p0.1.slurm
sbatch submit_rtm_eigs_p0.3.slurm
sbatch submit_rtm_eigs_p0.5.slurm
```

**2. Half-integer ladders.** The physical eigenvalue is followed by predicting the phase it should
have at the next rung and taking the nearest member of the block. Two things set whether that works.
The block members are separated in phase by roughly `π x /(v T)`, so a larger sound velocity crowds
them together — they sit about 2.5 times closer at `p = 0.5` than at `p = 0`. And the prediction is
an extrapolation, so its error grows with the step. At `p ≥ 0.3` the two meet and the branch is lost
after a few rungs. Halving the step halves the extrapolation error and separates them again. The new
rungs merge into the same cache, so the analysis sees one denser ladder.

Worth doing at `p = 0.1` as well: the branch is currently trackable only to `T = 11` there, and that
is the coupling all the quantitative results are quoted at.

```
sbatch submit_rtm_eigs_p0.1_fine.slurm
sbatch submit_rtm_eigs_p0.3_fine.slurm
sbatch submit_rtm_eigs_p0.5_fine.slurm
```

**3. Entropy sweeps.** `p = 0` stopped at `T = 14`, `p = 0.1` at `T = 8`, `p = 0.5` at `T = 7`. These
compute the temporal entropy as well as the spectrum, so they cost more per rung. The `p = 0` one is
the only arm without a wall, and its later rungs run for several hours each. At the other two the
entropy stops being meaningful past the wall, so those are worth having mainly as a record of where
it sits. All three need their checkpoint: a cold restart distorts the dome.

```
sbatch submit_rtm_p0.0.slurm
sbatch submit_rtm_p0.1.slurm
sbatch submit_rtm_p0.5.slurm
```

**4. Does β₀ move the wall?** The regulator scan showed the modulus gaps between the leading
eigenvalues grow linearly with β₀, and it is those gaps closing that ends the eigenvector route.
A larger β₀ should therefore buy some reach. These two run the same ladder at `nbeta = 12`
(β₀ = 0.6) and at the usual `nbeta = 4` (β₀ = 0.2) as the control, so the only difference is the
regulator; compare where the entropy dome breaks in each. The enhancement falls off as β₀/T², so
expect a modest shift rather than an escape. Only two rungs have ever run, so this one is open.

```
sbatch submit_betawall_nb12.slurm
sbatch submit_betawall_nb4.slurm
```

Everything writes to `results/data/cluster/` under its own filename, so nothing overwrites anything.

## Running it

First time, from anywhere:

```
git clone https://github.com/jmarquezol/master-thesis-BSC.git
cd master-thesis-BSC
module load julia/1.12.0
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Do not update the packages: `Manifest.toml` pins ITransverse to a specific commit, and floating it
has broken runs before. If you already have the repository, `git checkout main && git pull` is
enough. Then unpack the checkpoints as described below, and submit from inside `cluster/`:

```
cd cluster
julia --project=.. wall_scan_cluster.jl preflight    # must print "preflight OK"
sbatch submit_rtm_eigs_p0.1.slurm                    # and the rest, in the order above
```

Every submit script runs `wall_scan_cluster.jl preflight` first. It builds one small tMPO and exits,
the same call each rung makes, so a broken environment fails in a minute instead of marking every
`T` as errored. The driver also stops after two consecutive failures rather than burning the rest of
the ladder. If the preflight fails, run `Pkg.instantiate()` again and check that `Manifest.toml`
still pins ITransverse at `f10aee05`.

## Checkpoints

Checkpoints live in `cluster/checkpoints/checkpoint_<label>.jld2`. They are gitignored, so they do
not travel with the repository and come as a separate zip. Unpack it so the files land in that
directory, and check they are there before submitting:

```
unzip checkpoints_for_bsc.zip -d cluster/
ls cluster/checkpoints/
```

They are what makes a resubmission resume warm. Never delete them. The eigenvalue jobs still run
correctly without one, since they compute no entropy and only the first rung starts cold; the
entropy jobs do need theirs.

Check progress with `squeue -u $USER` and `tail -f logs/eigs_p0.1-<jobid>.out`. Each finished `T`
logs a line like `[rtm_eigs_p0.1] T=9.0 converged@61 k=4 |θ0|=1.5443 gap=0.731 ...`.

If a job stops at the walltime, just `sbatch` the same script again — it resumes from the last
checkpoint. Expect the `T = 20` arms to need a couple of resubmissions.

Push results back whenever, it's safe at any point:
```
git add results/data/cluster/*.jld2 && git commit -m "cluster sweep progress" && git push
```

## What's already done

The full-eigenvector RTM runs stop at a different `T` for each `p`, because the eigenvector route
hits the wall at different times. That is the physics, not a failed job. Runs that were superseded,
or that fall outside the `0 ≤ p ≤ 0.5` range of the thesis, are in
`results/data/cluster/archive/`; the raw per-`T` array-job outputs are in `cluster/archive/`.

## The jobs in `done/`

Kept for reference, and to regenerate something if needed. `submit_rtm.slurm`, `submit_cutoff.slurm`
and `submit_rdm.slurm` are the `p = 0.1` truncation comparison; the RDM one answered its question,
that the dome breaks at the same `T` as under RTM truncation with every rigidity in the block at
≈ 5·10⁻⁶, so extending it would only add rungs past the wall (NB9 §3). `submit_rtm_p0.3.slurm` is the `p = 0.3` entropy arm and
`submit_rtm_p1.0/1.5.slurm` are outside the `0 ≤ p ≤ 0.5` range of the thesis;
`submit_rtm_eigs_p0.0.slurm` is the complete `p = 0` eigenvalue ladder, and `submit_beta_*.slurm`
the regulator scans.

## Output format

`JLD2.load(path, "done")` gives a `Dict` keyed by `(label, T)`. Each entry is a named tuple with
`theta` (the spectrum), `theta_phys`/`i0` (the physical eigenvalue and its index), `tower_gap`,
`k_used`, and — for the full runs only — `s2_base` (the Rényi-2 dome) and `rigidity`.
