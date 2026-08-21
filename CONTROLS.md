# Controls and tests — the convergence investigation (2026-08-18 to 08-20)

Context file for Claude sessions. Records every numerical parameter and mechanism tested in the
failure-mechanism investigation, what was measured, and what entered the thesis. Everything here
is reproducible: every run stores its `Random.seed!` value, every claim traces to a cache in
`analysis/session_caches/` (or `results/data/cluster/`) and a script in
`analysis/session_scripts/`. Three readers print the tables:

```bash
cd analysis/session_scripts
julia --project=../.. battery_report.jl        # seed ensembles (dual band), dense, cutoff scan
julia --project=../.. dtreport.jl              # dt=0.05 ensembles vs dt=0.1 baseline
julia --project=../.. cluster_audit.jl         # every thesis number from cluster caches alone
```

Thesis anchors that carry these results: `app:failures` (tab:failrate, tab:rdm),
`app:conditioning` (tab:dense_cond), `app:errors` (controls paragraphs), `tab:blockpm_validation`
(corrected-column row), plus edited passages in `numerics.tex` (sec:reach), `P3.tex` (scheme
promise), `results.tex` (window + velocity paragraphs) and the conclusions.

---

## The failure mode, established

The single-vector (entropy) route fails by a DISCRETE, SEED-SELECTED event: at fixed (p, T) the
operator is identical for every seed, yet some seeds converge to unphysical fixed points while
the rest agree with neighbouring times. The failure probability rises with T at every coupling
and its onset moves earlier as p grows. Measured rates (loose band (0.05, 0.20) around
pi*c/16 = 0.0982; tight band (0.07, 0.15) shifts only marginal rungs):

| p | onset | rates (loose) |
|---|---|---|
| 0   | T=18 | 0/20 for T=8..17; 1/20 (18); 2/20 (22); 4/20 (24) |
| 0.1 | T~11 | 2/10 (11); 0/10 (12, but 2/10 tight); 0/10 (13) |
| 0.3 | T=7  | 0/10 (4..6); 2/10 (7); 4/10 (8); 6/10 (9, 10) |
| 0.5 | T=4  | 0/10 (3); 3/10 (4); 7/10 (5); 10/10 (6) |

Survivors stay correct through mid-ladder, but at the far end BOTH modes coexist: at p=0.3 T>=9
the surviving runs themselves spread by x2.5, and at p=0.1 T=12 (dt=0.05) survivors split into
two clusters, one below pi*c/16 (unphysical). Near the onset, only repetition distinguishes a
degraded survivor from a healthy rung. Never use ensemble dispersion (4x robust deviation) as the
failure metric — it undercounts when failures set the median and overcounts when survivors agree
to 4 decimals.

Two production events explained by the ensembles: the notorious p=0 T=17 plateau 0.724 was NOT
reproduced in 20 seeds (rarer than 1-in-20 luck); T=18, previously "good", produced -1.75 on one
of 20 seeds. There is no boundary between good and bad rungs, only a rate.

## Parameter-by-parameter verdicts

### Trotter step dt — accuracy control, NOT a failure cause
- Failure rate does not move: 19/100 fails at dt=0.1 vs 16/100 at dt=0.05 over matched rungs,
  direction flipping rung to rung. The wall does not move (failures persist at dt=0.05 at p=0.1
  T=12, p=0.3 T=8-9, p=0 T=24). Caches `seedens_dt005_*` (nbeta=8, trim 4/end, seeds offset 5e8).
- Bias is real and O(dt^2): dense operator scan (`dense_dtscan.jl`) shows |mu0| shifts shrinking
  ~x6 per halving; production plateau shifts <1% inside the window (ensemble medians p=0 T=18:
  0.1144 vs 0.1147).
- RETRACTED: the old "18% Trotter shift at the boundary rung" (p=0.1 T=11) was a single-run
  artefact — at the boundary, single-run dt comparisons measure the attractor lottery. app:errors
  sentence updated accordingly (approved 2026-08-20).
- Rigidity at matched T is dt-independent (p=0 T=18: 9.57e-6 vs 9.67e-6) or slightly worse at
  finer dt (p=0.1). Exact-operator rigidity converges to a dt-independent value; cond(V) grows
  with finer dt. Finer dt does not buy conditioning or reach.
- Cost of halving dt at fixed T: x1.2-1.5 (chi is set by physical T, not chain length; iteration
  counts unchanged). dt=0.025 (656 sites) produced a NaN blow-up — very long chains make the
  iteration MORE fragile.

### Truncation cutoff (block path) — affects the convergence path, not the answer
36-point seeded scan (`cutrerun_*`): the modulus gap g is identical to FIVE decimals across
cutoffs {1e-8, 1e-10, 1e-12} and seeds; plateau moves <1.5%. 1e-8 mostly terminates `stuck`,
tighter cutoffs mostly `converged` — but seed-dependent IN BOTH DIRECTIONS (1e-10 stuck on seed 2
at p=0.3 T=3; 1e-8 converged on seed 2 at p=0.5 T=3). Fixed-boundary control (p=0, |Up>):
every cutoff converges in 17-25 iterations, seeds irrelevant. No "tighter cutoff fixes it" claim
survives. In thesis: app:errors fourth-control paragraph (honest version, approved).

### Bond dimension chi — not binding
chi = 32/64/128 give identical niters, reason, g, plateau at fixed seed (`chi_test_*`): the
temporal states reach bond dimensions 6-16 and never touch the cap. Earlier controls: plateau
chi-independent to ~2%, |mu0| to 1e-4, but the REAL-part chord slope moves by up to 0.12 —
one reason the plateau, not the slope, is the quoted estimator.

### Truncation scheme RTM vs RDM — same physics, different verdict, x12-33 cost
8 matched rungs (`rdmswap_*`, tab:rdm): RDM reproduces the RTM plateau and |mu0| everywhere,
including T=24 at p=0 and T=11 at p=0.1 (12.1 h for that single seed). RDM reports `converged`
at every rung; RTM reports `stuck` at essentially every T>~10 — for correct and failed runs
alike. The RTM stagnation is its truncation noise floor, not the transfer matrix. Consequence:
the `stuck/converged` flag under RTM carries NO information about correctness (the only
RTM-converged ensemble in ~350 runs is p=0.5 T=3, 10/10). Fulfils P3's scheme-agreement promise.

### Phase rigidity — correlate, not control
Falls geometrically with T, faster with p, but does not determine failure: p=0.1 T=13 has the
LOWEST rigidity in the battery (7.9e-10) and 0/10 failures, while p=0.5 T=6 fails 10/10 at
rigidity three orders HIGHER (1.9e-6). Within one coupling rigidity is monotone while failures
are probabilistic. dt-independent at matched T.

### Eigenvalues vs eigenvectors — the robust/fragile split
|mu0| reproduces across seeds at every usable time (~1e-4..1e-6 relative); eigenvector
functionals (entropies) do not. Supports the eigenvalue-only mode and the repetition procedure.
Consensus caveat measured: two seeds CAN agree on a wrong value (p=0.1 T=13 accepted plateau
0.1843 -> c=0.939; p=0.5 T=7 accepted with spread 0.0000 at a wrong value) — reproducibility
alone is insufficient, smoothness in T is also required (as sec:reach states).

### Conditioning — set by the COUPLING, not by the spectral gap (refuted mechanism)
Exact diagonalisation of the corrected bulk5 column (`dense_bulk5.jl`, tab:dense_cond): at
matched dimension ~2190, cond(V) = 3.0e2 (p=0) vs 2.2e11 / 2.1e10 / 9.6e9 (p=0.1/0.3/0.5) while
the leading gap is comparable or LARGER at p!=0 and rigidity is healthy (0.12-0.36) everywhere.
This refutes the old conclusions claim that emergent dual unitarity's closing gap makes
de-mixing ill conditioned (removed 2026-08-19; also removed from app:blockpm:schur). The
perturbation-theory asymmetry (gap denominator in mixing, none in the eigenvalue shift) stands
as motivation for the eigenvalue-only mode. Bound looseness: first-order |mu0| shifts sit 1-4
orders BELOW ||dE||/r0 for truncation-like and random perturbations — eigenvalues are more
stable than the bound guarantees. Block PM validated against the exact bulk5 operator to 3.6e-7
(`validate_bulk5.jl`, row added to tab:blockpm_validation).

### vT scaling of finite-time corrections — holds for dimensions, fails for the plateau
x1 deviation at matched T falls with p (0.023 / 0.0052 / 0.0028 at p=0.1/0.3/0.5, T=2) ✓.
Entropy-plateau offset GROWS with p at fixed vT (0.0270 / 0.0320 / 0.0351 / 0.0372 at vT~16).
results.tex velocity paragraph scoped to the spectral observables accordingly.

## Open questions (marked open in app:failures and the conclusions — do not overstate)
1. WHY the coupling controls the eigenvector-basis conditioning (no mechanism).
2. Why eigenvalue stability so far exceeds the 1/rigidity bound (structured suppression).
3. No single-run criterion separates degraded survivors from healthy rungs near the onset.
4. Whether ANY truncation setting moves the failure onset (cutoff scan was small-T block-path;
   RDM failure statistics are n=1 per rung).
5. Block-route (spectral) failure statistics never ensembled; spectral windows rest on the
   phase-increment criterion.

## Cluster round 2 (2026-08-19, Pau) — health checks all green
15 arms, 0 error records, 0 escalations, |mu0| jumps <=0.18%, tower T=2 guard passes at all
five couplings (2e-6..4e-5 vs ent-arm). Arms are WALLTIME-limited, not physics-limited at p=0
(ent reached T=10.5 of 24) — resubmission safe and encouraged, priority ent/eigs/tower at p=0,
then p=0.1; p>=0.3 already reached their physics windows. Audit c = 0.451/0.466/0.494/0.527.
Figures regenerated from round 2 (`fig_tower.jl`, `fig_spectral.jl`, `fig_domes.jl` — merge
lists mirror cluster_audit.jl; p=0 tower panel now the k=8 cluster arm).

## Verification work also done this session
- Appendix C (JW / four-fermion) fully verified: every equation exact as a matrix identity,
  spin vs fermion spectra identical at 4 parameter sets (`jw_check.jl`); wick1954 added for the
  rotation (wick1950 = the theorem); opening sentence scoped; degree-preservation clause added.

## Gotchas for anyone rerunning
- Trim entropy profiles by nbeta/2 bonds per end: [3:end-2] at nbeta=4, [5:end-4] at nbeta=8,
  [9:end-8] at nbeta=16. Cluster caches are pre-trimmed; local svpm/seedens caches are NOT.
- Seeding schemes: seedens = 1e6*p + 1000*round(10T) + s; dt=0.05 ensembles offset +5e8;
  block cutrerun offset +100*seed. Production lanes before 2026-08-18 are NOT reproducible
  (global RNG, unseeded).
- `maxdims`/`cutoffs` ramps are INERT on the single-vector path (only the block path honours
  them); svpm production ran at fixed cutoff=1e-12.
- `pm_itercheck!` only counts stuck iterations once ds < 0.05 (hidden eps_max gate).
- Failure bands: report loose AND tight; quote counts ("two of ten seeds"), never bare
  percentages at n=10.
