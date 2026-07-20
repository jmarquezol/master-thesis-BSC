# Master wall-scan sweep — BSC cluster rerun (v2, warm-started)

## What's in this master job

**8 independent SLURM jobs, all running in parallel**, each a self-contained warm-started T ladder:

| job | truncation | p (NNN coupling) | T range |
|---|---|---|---|
| `submit_rtm.slurm` | RTM | 0.1 | 2..**20** |
| `submit_rdm.slurm` | RDM | 0.1 | 2..12 |
| `submit_cutoff.slurm` | RTM, tighter cutoff | 0.1 | 2..14 |
| `submit_rtm_p0.0.slurm` | RTM | 0.0 | 2..20 |
| `submit_rtm_p0.3.slurm` | RTM | 0.3 | 2..20 |
| `submit_rtm_p0.5.slurm` | RTM | 0.5 | 2..20 |
| `submit_rtm_p1.0.slurm` | RTM | 1.0 | 2..20 |
| `submit_rtm_p1.5.slurm` | RTM | 1.5 | 2..20 |

The first three are the original ε-scan (bond dimension/cutoff/truncation-algorithm independence,
see NB9) — `submit_rtm.slurm`'s ladder was extended from T=14 to T=20 so it now doubles as the
p=0.1 reference for the p-sweep arms below it. The five `submit_rtm_p*.slurm` jobs are new: a
p-sweep at fixed χ=64, always through the RTM route (NB9 found RDM buys no physical improvement
for 4-11x the cost, so it's not worth the p-sweep's extra jobs), reaching T=20 to see well past the
wall for every p. All eight write into the same `results/data/cluster/warm_sweep.jld2`, keyed by
`(label, T)` (labels: `rtm64_full`, `rdm64`, `cut_tight`, `rtm_p0.0`, `rtm_p0.3`, `rtm_p0.5`,
`rtm_p1.0`, `rtm_p1.5`), so nothing clobbers anything else — submit all eight at once.

The driver (`wall_scan_cluster.jl`) takes the p-sweep's `p` and `Tmax` as command-line args now
(`julia wall_scan_cluster.jl psweep <p> <Tmax>`); the three original modes (`rtm`/`rdm`/`cutoff`)
are unchanged and need no extra args. **The checkpoint/resume mechanism was verified locally this
session** with an explicit kill-mid-ladder-then-resubmit test: the resumed process logged
`warm-resumed from checkpoint at T=2.0` and continued correctly from there — this matters a lot
here, since T=20 is far enough that several of these jobs (especially the higher-p arms, which may
also trigger the k=4→6 escalation) will very likely need more than one 48h walltime window.

## Why a rerun (background — the original 3-job ε-scan)

Your first pass (`cluster/{rtm,rdm,cutoff}_array/`, one independent SLURM task per T, no warm
start between T's) got real, useful data — but analyzing it turned up two things worth fixing
before we treat the dome/wall-location numbers as final:

1. **No warm start ⇒ the dome inflated at T≈6 instead of T≈10.** Phase rigidity (the eigenvector-
   conditioning number) was identical to our warm-started local baseline at every T — so the
   *physics* conclusion (RDM doesn't help; the wall is truncation-independent) survives — but the
   entropy/dome broke ~4 steps earlier than it should, because a cold random restart at each T has
   nothing biasing it onto the same physical branch as the near-degenerate cluster tightens.
2. **The rank-based selector picked the wrong eigenvalue** at T=12/13 in the `rtm` array (`|θ_phys|`
   jumped to 1.78/1.98, then snapped back to 1.55 at T=14) — a selection failure, not physics.

Both are fixed in `cluster/wall_scan_cluster.jl` now: it warm-starts every T from the previous one
(seeding **and** continuity-anchoring), it checkpoints the converged blocks to disk after every T
so a walltime-killed job resumes truly warm on resubmission (not just cold-with-an-anchor), and the
selector (`pick_phys_continuity`, from `src/transverse_tools.jl`) searches the *whole* block instead
of only the top-2 by modulus. It also auto-escalates k=4→6 whenever a block has no tower member
besides λ0 (the `p≥0.3` failure mode found separately in NB5) — not needed at p=0.1, but free
insurance. **Nothing from your first pass is discarded** — it's valuable data for exactly the
comparison in point 1 above, and stays where it is.

## Steps

1. **Push anything not yet pushed** from your current cold-array run first — nothing already
   computed should be lost:
   ```
   git add cluster/*_array/worker_results_T*.jld2
   git commit -m "cold array sweep: whatever finished"
   git push
   ```

2. **Cancel the remaining cold jobs** (workers and their pending `_collect` dependency jobs):
   ```
   squeue -u $USER | grep wallscan
   scancel <job-ids>
   ```

3. **Pull the updated repo** (new `cluster/wall_scan_cluster.jl`, ready-to-run `submit_*.slurm`
   with your account/qos/module already filled in, plus the `src/transverse_tools.jl` additions the
   driver depends on) and re-instantiate — the Manifest changed (your MKL addition is already
   merged in, this just re-syncs the lock to it):
   ```
   git checkout main && git pull
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

4. **Submit all eight** (parameters are already filled in from your working `*_array/submit.sh` —
   `--account=bsc21`, `--qos=gp_bsccase`, `--cpus-per-task=16`, `julia/1.12.0` — only touch them if
   anything about your allocation has changed):
   ```
   sbatch cluster/submit_rtm.slurm
   sbatch cluster/submit_rdm.slurm
   sbatch cluster/submit_cutoff.slurm
   sbatch cluster/submit_rtm_p0.0.slurm
   sbatch cluster/submit_rtm_p0.3.slurm
   sbatch cluster/submit_rtm_p0.5.slurm
   sbatch cluster/submit_rtm_p1.0.slurm
   sbatch cluster/submit_rtm_p1.5.slurm
   ```
   These are single sequential jobs (not job arrays) — one per (truncation, p) configuration, each
   running its own warm-started T ladder internally. `rdm` is capped at T=12 (cold T=9 alone took
   ~20.6h under `:rdm`); `cutoff` stays at T=14; `rtm` and all five `rtm_p*` arms go to T=20.

5. **Monitor:**
   ```
   squeue -u $USER
   tail -f cluster/logs/wallscan_rtm-<jobid>.out           # etc. for rdm / cutoff
   tail -f cluster/logs/wallscan_rtm_p0.3-<jobid>.out       # etc. for the other p-sweep arms
   ```
   Each completed T logs a line like
   `[rtm64_full] T=9.0  stuck@1069  k=4  |θ0|=1.5443  gap=0.731  peak=0.3331  r=[...]  4533s`
   The p-sweep arms' label includes p, e.g. `[rtm_p0.3] T=9.0  ...`.

6. **If a job hits its walltime cap (2 days)**, it stops mid-ladder — no data is lost (every
   completed T is cached, and the last converged blocks are checkpointed to
   `cluster/checkpoint_<label>.jld2`, gitignored). Just resubmit the **same** `sbatch` command: it
   picks up from the last checkpoint, truly warm, not a cold restart. **Expect this to happen at
   least once** for the T=20 arms (`rtm` and all five `rtm_p*`) — reaching T=20 is a lot further
   than the original T=14, and this exact resume path was verified locally before this rerun (a
   kill-mid-ladder test produced the log line `warm-resumed from checkpoint at T=2.0` and then
   continued correctly).

7. **Push results back periodically** (safe at any point — the cache is always internally
   consistent):
   ```
   git add results/data/cluster/warm_sweep.jld2
   git commit -m "warm cluster sweep: progress update"
   git push
   ```

## What the output looks like

`results/data/cluster/warm_sweep.jld2` holds a `Dict` keyed by `(label, T)` — e.g.
`("rdm64", 9.0)` — where each value is a named tuple with `theta`, `theta_phys`, `i0`, `dphi`/`cls`
(the phase classification per spectrum member), `tower_gap`, `k_used`/`escalated` (did this T need
the k=4→6 retry?), `s2_base` (the trimmed Rényi-2 dome), `peak`, `rigidity`, `reason`
(`"converged"` or `"stuck"`), `niters`, and `elapsed`. Loadable with `JLD2.load(path, "done")`.

---

## v3 (2026-07-20) — corrected campaign after the v2 post-mortem

The v2 run mostly failed for infrastructure reasons, not physics (logs archived in
`cluster/archive_v2/copy_logs_v2/`, the clobbered cache in `cluster/archive_v2/warm_sweep_v2.jld2`).
Fixes in this version:

- **Per-label caches (the critical fix).** Every arm now writes its OWN
  `results/data/cluster/sweep_<label>.jld2`. v2 had all 8 jobs writing the SAME `warm_sweep.jld2`
  concurrently → a race clobbered most of the 2-day compute (only 16 of ~60 points survived).
  Assemble the notebooks' `warm_sweep.jld2` afterward with `julia --project=. cluster/merge_sweeps.jl`.
- **`--mem=120G`** in every script (v2 had none → 4 OOM kills). Raise toward `--mem=0` (whole node) if
  OOM recurs.
- **Tight cutoff `[1e-12…; 1e-14]`** is now the default (was 1e-8/1e-10). It fixed the v2 wrong-branch
  jump at T=10 AND converges far better (e.g. p=0.1 T=3: 32 iters vs 461 stuck). The separate
  `cutoff` mode is gone.
- **k=4, no k=6 escalation.** The escalation was a workaround for misclassifying the physical
  π-displaced partners; x1/tower selection is post-processing. Also lighter on memory.
- **Validity guard.** A solve with |θ0|<1e-6 / non-finite / >50 is cached as `:error` (retried, not
  skipped) and breaks the warm chain — stops the v2 large-p null-collapse cascade.
- **p arms:** {0.0, 0.1, 0.3, 0.5} core + {0.8, 1.0, 1.5} reach, all `psweep <p> 14`. (p=1.0/1.5 no
  longer null-collapse with the tight cutoff, but their |θ0| is large/growing — treat as exploratory.)
- **RDM:** run SEPARATELY and AFTER the RTM arms land clean (`submit_rdm.slurm`, own cache, T=12).
  WATCH its first ~10 min: v2's rdm produced zero output for 2 days — if silent again, kill & report.

**Submit order:** the 7 `submit_rtm_p*.slurm` (in parallel) → when clean, `merge_sweeps.jl` → then
`submit_rdm.slurm`. Each job is independent, crash-safe (resubmit resumes from its checkpoint), and
memory-capped.
