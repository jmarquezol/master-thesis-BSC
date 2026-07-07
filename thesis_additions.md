# Thesis additions — synthesis of the July 2026 campaign
### Meeting brief + everything new + where it goes in `thesisdraft.md`

*(Companion to `barrier_section.md` (results chapter draft) and `blockpm_methods.md` (methods
§5.4 draft), both already written in the thesis's style and ready to paste. This document is the
map: what we did, what we found, what we refuted, what it means, and exactly where each piece
belongs in the manuscript. Every number below comes from an executed notebook output; figures are
referenced by filename in `results/imgs/`.)*

---

## 1. Meeting brief — six things to tell Stefano

1. **The thesis has a stronger frame now.** Instead of one result ("temporal Ising universality
   survives NNN frustration"), the narrative is **the entanglement barrier across four models** —
   Ising, ANNNI/Alcaraz, XXZ, tricritical O'Brien–Fendley — with a *three-dial law* for where the
   method's reach ends: (i) central charge, (ii) symmetry content of the quench, (iii) MPO
   symmetry. Each model isolates a different dial. Draft results chapter: `barrier_section.md`.

2. **The Alcaraz headline stands and is now better defended**: temporal central charge
   **c(p=0.1) = 0.47 ± 0.05** (clean window T=4..9, Rényi-2 chord slope, calibrated against p=0),
   boundary exponent x₁ ≈ 0.50 for both p=0 and p=0.1 (free-BC Ising). The solver behind it is now
   validated against exact diagonalization (point 5).

3. **XXZ is fully mapped, and the result is not what we expected in June.** The wall at T≈4 is
   NOT frustration and NOT parity oscillations: it is the **exact Z₂ degeneracy of the Néel
   quench** (two Néel states → two exactly degenerate transfer fixed points, unmasked at T≈4 as a
   4-fold band). The robust extraction gives **c_eff ≈ 0.75 at Δ=1** (stable T=3..8; marginal-log
   suspect) and **≈0.9–1.0 at Δ=0.5**, from the *calibrated* Im S₂ → πc/12 estimator. The famous
   parity oscillations are a **null result** in the temporal entropy (amplitude ~10⁻³).

4. **The tricritical point is resolved honestly**: equilibrium λ_c located between 0.42 and 0.43
   by an entropy-free finite-size gap (the S(N/2) estimator *overshoots* c — known
   corrections-to-scaling at TCI points, not a bug; a velocity-free energy check is consistent
   with c=7/10). Dynamically there is essentially **no clean window**: c=7/10 merges the leading
   transfer spectrum into a band almost immediately. The earlier "solver failure" was disproved —
   the eigenvalues were always right; the wobble was the physical-λ₀ *selector* hopping inside a
   genuinely degenerate band.

5. **The machinery is now validated against ground truth.** The block power method reproduces
   exact dense diagonalization of small transfer matrices on all four models (errors
   1.4×10⁻⁸ / 4.6×10⁻¹³ / 4.2×10⁻¹⁰ / 2.3×10⁻¹⁰ for Alcaraz-VD2 / Alcaraz-WII / Ising-Murg /
   tricritical). First check in the project that does not route through another power method.

6. **New this week — the asymmetry experiment (NB12), resolved (with a self-correction worth
   telling).** We built the first genuinely symmetric XXZ propagator (Murg-type, exact bond-2
   commuting Pauli layers; the package's SymSVD is demonstrably not symmetric). Stress-testing the
   plan surfaced a subtlety: the XXZ propagator is an exactly symmetric *operator*
   (‖U−Uᵀ‖~10⁻¹⁸, inevitable since H is real-symmetric) but is **not manifestly symmetric as an
   MPO tensor** — its YY layer carries σʸ (transpose-antisymmetric), and σˣ,σʸ,σᶻ cannot all be
   symmetric in any single-site basis. So the symmetric-Takagi solver is **structurally
   unavailable to XXZ** (a result in itself: it's a privilege of Ising-like models), and our
   direct Takagi runs were invalid. **The conclusion — the wall does not move — is instead
   carried by two gauge-free routes**: (a) the *construction-independence test* — the
   independent Murg construction run through the ordinary two-sided solver reproduces VD2 to
   3–4 digits, wall and all (two completely different constructions, one solver, same physics);
   and (b) the spectrum bridge (same Z₂ band in the symmetric tMPO). Plus a benign control: the
   pure TFIM (Alcaraz p=0) through our *same* generic pipeline keeps a flat |λ₀|≈1.49 with **no
   band** out to T=12 — the XXZ band at T≈4 is a property of the quench, not a broken pipeline. **So the wall is physics,
   not our adaptation of the package** — quench symmetry (dial ii) dominates MPO symmetry (dial
   iii). What *does* give Ising its T=14 is genuinely its symmetric MPO + Takagi + small dₜ +
   single-sector quench + low c — a privileged combination XXZ is denied on every count.

**One-sentence summary for the meeting:** *we turned a single-model result into a validated,
falsification-driven map of where and why transverse contraction reaches its limit, with the
limit's position explained by physics (central charge + quench symmetry), not by implementation.*

---

## 2. The reframe and the three-dial law

**Why the reframe.** Every model we touched eventually failed in the same characteristic way
(entropy dome inflates; leading transfer eigenvalues collapse onto a band of equal moduli). Chasing
each failure produced, cumulatively, something more valuable than any single c-value: a *unified,
mechanistic account* of the method's horizon. That account is the thesis's methodological
contribution; the individual c-extractions are its physics payloads.

**The mechanism (one paragraph).** A quench *to* a critical point drives the transfer matrix
toward a rescaled unitary — emergent dual unitarity. A unitary has no modulus hierarchy, so the gap
ratio |μ₁/μ₀| → 1 is the *defining approach* to the conformal regime, not a numerical accident.
Non-Hermitian linear algebra then bites: eigen*values* stay well-conditioned as the gap closes, but
eigen*vectors* — what the generalized entropies are built from — degrade with condition number
~1/gap. The better the conformal physics, the less well-defined the object we want to measure.
**The barrier is not an obstacle in front of the physics; it is the physics, seen from the
numerical side.**

**The reach table** (the thesis's central table):

| model | MPO symmetry | c | quench symmetry | degeneracy type | clean reach |
|---|---|---|---|---|---|
| Ising | symmetric (Takagi) | 1/2 | symmetric (\|X+⟩) | none until barrier | T ≈ 14 |
| Alcaraz p=0.1 | asymmetric | 1/2 | symmetric (\|X+⟩) | dynamical (barrier) | T ≈ 10 |
| XXZ Néel | **both tested** | 1 | **breaking** (Néel) | **exact, from T=0** | T ≈ 4 (either MPO) |
| tricritical | asymmetric | 7/10 | symmetric (\|X+⟩) | band (charge-driven) | ≲ 2 |

**The three dials.** (i) Larger c ⇒ faster gap closing ⇒ earlier band (tricritical is the extreme).
(ii) A symmetry-*breaking* initial state contributes an **exact** degeneracy present from T=0, only
masked by a decaying cat-state splitting (XXZ) — worse than any dynamically-developing
near-degeneracy. (iii) A symmetric MPO unlocks the better-conditioned Takagi route (Ising's whole
advantage over Alcaraz at equal c) — **except** when dial (ii) is exact: NB12 shows symmetrizing
XXZ's MPO leaves the reach unchanged. **The punchiest cross-model sentence: an exact symmetry
degeneracy of the quench is a worse enemy than frustration** (XXZ walls at 4, frustrated Alcaraz
at 10).

---

## 3. Results, model by model

### 3.1 Ising — the control (NB6)
- Reproduces Carignano–Tagliacozzo Eq. (6) with the symmetric Murg + powermethod_sym + Takagi n→1
  machinery: Re S on the c=1/2, s₀≈0.3 chord (peak 0.37→0.44→0.47 at T=4,8,12), Im S ≈ 0.13–0.14 ≈
  π/24. Clean out to **T=14**.
- λ₀(T) traces a circle of near-constant radius in the complex plane — emergent dual unitarity as
  a spectral statement. Figures: `cft_ising_validation.png`, `ising_lambda0_circle.png`.

### 3.2 Alcaraz / ANNNI-type — the headline (NB4, NB5, NB7)
- Equilibrium: DMRG gives c ≈ 1/2 across p (three independent reads). Figures: `p_dependence.png`,
  `cft_L.png`.
- **Headline: c(p=0.1) = 0.47 ± 0.05** from the clean-window (T=4..9) Rényi-2 chord slope,
  per-T values 0.43, 0.46, 0.56, 0.48, 0.46, 0.42, calibrated against p=0. Temporal Ising
  universality **survives NNN frustration**. Figures: `temporal_entropy_profiles.png`,
  `temporal_chord_fit.png`.
- Boundary exponent (Eq. 4): x₁(p=0)=0.502, x₁(p=0.1)=0.497 — free-BC Ising, robust.
- The wall at T≈10: eigenvalues stay reproducible (~10⁻⁵) but the eigenvector is ill-conditioned;
  every repair attempted (subspace projectors, continuity projection, seeds) fails for the
  structural reason of §2. Frustration closes the gap faster than pure Ising (partner-filtered gap
  0.35→0.98 over T=1..6). Figures: `gap_closing_wall.png`, `alcaraz_gap_dualunitarity.png`,
  `block_pm_ising_vs_alcaraz.png`.

### 3.3 XXZ Néel — the symmetry wall (NB8, NB9, NB12)
- **Why Néel is forced** (NB8 §0, NB9 §0 — pedagogical sections written for direct thesis reuse):
  XXZ has no transverse field, so every flip-symmetric product state is nearly an eigenstate
  (measured: |X+⟩ gives Re S₂ ≈ 0, χ=4). The only entangling product quench breaks the Z₂
  spin-flip symmetry — and its two Néel copies hand the transfer matrix an exactly degenerate
  leading pair from T=0.
- Equilibrium: DMRG c ≈ 1 across Δ ∈ [−1,1]. Figure: `xxz_c_equilibrium.png`.
- **The wall anatomy** (NB9 §4): at T≤3 the dominant eigenvalue is isolated (T=1 spectrum
  [0.98, 0.35, 0.35, 0.20] — note the degenerate *subleading* pair, the sector fingerprint);
  between T=3 and 4 a **4-fold band** (2 sectors × ± partners) locks in within 3%, exactly where
  the Re dome jumps discontinuously (Δ=1 peak 0.25→0.81). Figure: `xxz_vs_alcaraz_gap.png`.
- **The dome inflation is deterministic, not random mixing**: 3-seed test identical to 4 digits
  (0.2500×3 at T=4; 0.8891×3 at T=6). The truncated (nonlinear) iteration has a unique attractor
  that is no longer a pure eigenvector.
- **c extraction**: Re chord unusable past the jump (fits c≈6–10); the robust estimator is
  **Im S₂ → πc/12** — *calibrated*, not derived (see §4.3): **c_eff ≈ 0.75 ± 0.05 at Δ=1**
  (T=3..8, remarkably stable; marginal-operator logs at the Heisenberg point are the suspected
  bias) and **≈ 0.9–1.0 at Δ=0.5** in its clean windows. Figure: `xxz_entropy_profiles.png`.
- **Oscillations: null result** — staggered amplitude ~10⁻³ at both Δ, three orders below the
  dome. The supervisor's conjecture was tested and cleanly bounded. Figure: `xxz_oscillations.png`.
- **The intra-sector barrier is genuinely slow** (k=6 resolves past the symmetry cluster):
  |θ₅/θ₁| closes 0.14→0.75 over T=1..5 vs Alcaraz's 0.35→0.97 — "NN + unfrustrated ⇒ slow
  barrier" is TRUE in the corrected observable; the reach is killed by the symmetry cluster
  sitting on top of a slow barrier. Figure: `xxz_intrasector_gap.png`.
- **The asymmetry experiment (NB12) — resolved, see §4.4/§5**: symmetric Murg MPO built and
  verified; the wall is unchanged at T≈4, established gauge-free by (a) the
  **construction-independence test** — the independent Murg construction (= symmetric order-2 to 4 digits) through the two-sided
  solver reproduces VD2 to 3–4 digits including the T=4→5 inflation, and (b) the spectrum bridge
  (same Z₂ band, T=4: [0.864,0.829,0.829,0.826] vs [0.893,0.885,0.867,0.867]). The *direct*
  symmetric-Takagi test was found **gauge-invalid** (the σʸ obstruction: XXZ has no manifestly
  symmetric single-site-Murg MPO) — a self-correction, not a retraction; the conclusion (dial ii
  dominates dial iii) stands on the gauge-free routes.

### 3.4 Tricritical O'Brien–Fendley — the charge wall (NB10, NB11)
- **Equilibrium success**: the entropy-free finite-size gap brackets **λ_c between 0.42 and 0.43**
  (N·gap plateaus vs grows); gaplessness ends at a *point* ⇒ genuine TCI, no floating phase.
  S(N/2) *overshoots* c (0.79→1.19, never 7/10) — the known corrections-to-scaling pathology of
  tricritical points (calibrated: same method gives 0.50 on the Ising line); the velocity-free
  ground-state-energy estimator (v cancels in the ratio, only x₁=3/40 enters) is consistent with
  c = 7/10. Figures: `tricritical_c_equilibrium.png`, `tricritical_gap.png`.
- **Dynamics: no clean window, and that is the finding.** The leading spectrum is a band almost
  immediately (5 eigenvalues within ~15% by T=3 at k=6; top four within 4% by T=4.5). At the wall
  point T=4.5 the old solver stuck at 527 iterations and the fixed, ground-truth-validated solver
  also failed to converge (killed after 7.8 h) — **not negotiable by better internals**.
- The earlier "non-converged oscillating gap" claim was **disproved point by point**: new and old
  eigenvalues agree to 3–4 digits wherever both converge; the |λ₀| wobble was the pick_phys
  *selector* reporting different band members (at T=3 it picked the second-largest θ). Inside a
  band, "the physical λ₀" is not a well-posed quantity — report |θ₁| and the band spread.
- The cross-model ordering survives with the selector-free reading: **Ising < Alcaraz <
  tricritical in time-to-close** — higher c / frustration closes the barrier faster.

---

## 4. Methods contributions (and where each goes in the thesis)

### 4.1 The block power method → new §5.4 (paste `blockpm_methods.md`)
Oblique (Petrov–Galerkin) Rayleigh–Ritz on k left + k right temporal MPS; **non-conjugating**
pencils S=⟨L|R⟩, M=⟨L|E|R⟩ (bilinear, per §5.1's biorthogonality); pinv(S)·M instead of
eigen(M,S) near singular S; exact left-right pairing from ONE decomposition
(u_j = pinv(S)ᵀ(V⁻¹)ᵀe_j — guarantees pairing and biorthogonality by construction); RTM vs RDM
truncation of the de-mixed pairs; warm-starting along the T-ladder; the ±λ₀ partner structure and
the partner-filtered physical gap. Already written in draft style with local equation numbers
(B1)–(B9).

### 4.2 The validation methodology → new Appendix B (novel emphasis for a master's thesis)
- **Exact dense ground truth**: contract small tMPOs (N_s ≤ ~12 sites depending on temporal
  dimension) to dense matrices, diagonalize exactly, compare — the only check that does not route
  through another power method. Accuracies achieved: 10⁻⁸–10⁻¹³ across all four models (NB11 V2).
- **Planted-spectrum synthetic tests**: the same Rayleigh–Ritz iteration on random non-normal
  matrices with known eigenvalues (± pairs, tunable 4-fold clusters) — isolates the algorithm from
  truncation. Honest outcome: the suspected pairing defect performs identically to the fix in
  exact arithmetic — the fixes are strictly-better implementation, not the explanation of any
  observed failure.
- **Seed tests as degeneracy diagnostics**: N independent random seeds; seed-independence ⇒
  unique attractor; seed-dependence would ⇒ unresolved degenerate manifold. Used twice with
  decisive results (asymmetric XXZ: deterministic inflation; symmetric XXZ: unique attractor).
- **Pre-registered decision branches**: for each experiment we wrote down, before running, what
  each possible outcome would mean (NB12 §4b/§5 are explicit examples). Cheap discipline, kept us
  honest twice when the data contradicted the favored hypothesis.

### 4.3 The estimator toolkit → additions to §5.2 + Appendix B
- **Partner-filtered gap**: the transfer spectrum comes in near ±pairs; the naive |θ₂|/|θ₁|
  measures the pair *splitting*, not the spectral gap. Filter the −λ₀ partner (nearest to −θ_{i₀}).
- **Selector-free band reading**: inside a near-degenerate band, *any* "physical λ₀" selection
  hops between members and manufactures spurious oscillations — report |θ₁| (moduli-sorted) and
  the band spread instead. (This single lesson dissolved the entire "tricritical solver failure".)
- **Cluster-filtered intra-sector gap**: with an exact symmetry degeneracy, k must exceed the
  cluster size; |θ₅|/|θ₁| at k=6 is the true barrier observable for XXZ.
- **Im S₂ calibration**: C–T Eq. (6)'s iπc/12 offset is derived only for n→1. The naive Rényi
  continuation predicts πc/16 for S₂; measurement on the p=0 Ising line (c=1/2 known three
  independent ways) gives 0.127–0.131 ≈ π/24 = πc/12 — the offset is **empirically
  n-independent**, so c = 12·Im S₂/π, with a ~10% systematic (the offset drifts +12% at p=0.1).
  No XXZ data entered the calibration; it was applied blind.

### 4.4 The symmetric Murg-XXZ propagator + the operator-vs-tensor symmetry distinction → new §4.3
The rotated Néel-frame two-site term decomposes as SxSx − SySy − ΔSzSz: three **mutually
commuting layers**, each an exact bond-2 Murg cos/sin factorization (the Ising construction
generalized to any Pauli). Two kernels: order=2 palindrome e^{ZZ/2}e^{YY/2}e^{XX}e^{YY/2}e^{ZZ/2}
(d_t=32) and order=1 single sandwich e^{ZZ}e^{YY}e^{XX} (d_t=8, 16× cheaper). Validated: order-1
echo error vs TDVP 1.1×10⁻⁵ at dt=0.05 (better than the palindrome's own 2.6×10⁻⁵, both at VD2's
level); entropy cross-check at (Δ=0.5,T=4) matches order-2 to 4 digits (χ=8, peak 0.8632, both).

**The subtle and important part (worth a boxed remark in the thesis): "symmetric operator" ≠
"manifestly-symmetric MPO."** Every propagator we use is an exactly transpose-symmetric *operator*
(‖U−Uᵀ‖~10⁻¹⁸, since our H are real-symmetric — so even "asymmetric" VD2 is a symmetric operator).
But the symmetric-solver machinery (`powermethod_sym` + Takagi) needs the MPO *tensor* to be
symmetric under a *raw* index swap, with no gauge search. Ising's Murg MPO is manifestly symmetric
(its coupling σˣ is transpose-symmetric); XXZ's is **not** (its YY layer stores σʸ, which is
transpose-antisymmetric), and — the structural punchline — **no single-site basis makes σˣ,σʸ,σᶻ
simultaneously symmetric** (three mutually anticommuting operators; at most two symmetrizable). So
the symmetric-Takagi route is *structurally unavailable* to XXZ, and to any model requiring three
anticommuting couplings. This is why our direct Takagi runs (NB12 §4) were invalid, and why the
dial-(iii) conclusion had to be — and was — established by the gauge-free construction-independence
test (§3.3). Also: ITransverse's SymSVD builder fails even the operator-symmetric-construction bar
(normdiff 0.07–0.45), which is what motivated building the Murg version from scratch.

### 4.5 Performance engineering → short paragraphs in §5.4 / Appendix B
- **RTM joint truncation of the de-mixed (L,R) pairs is worth 97×** (27 201 s → 280 s at the
  Alcaraz T=6 reference; converges at χ=9 where per-vector RDM needs χ=44). Confirmed from the
  other side in July: RDM at tricritical T=3 ran >3 h vs RTM's 775 s.
- **Warm-starting the T-ladder** (`pad_tmps`: re-index converged tensors onto the longer time-site
  set, random tail) — the largest iteration-count saver in every sweep.
- **A trap worth documenting**: `powermethod_sym` applies its bond-dimension schedule *per
  iteration*, so a 2:2:64 ramp truncates a warm-started χ~15 vector to χ=2 at step 1 — warm rungs
  need a fixed cap. (Cost us one confused afternoon; now in the notebook and the methods text.)
- Kernel choice is a physics decision: VD2 (2nd order) for NNN models; WII is genuinely 1st order
  there but effectively 2nd order — and ~5× cheaper — for strictly-NN models like XXZ (used as an
  independent cross-check kernel, reproducing the VD2 story point for point).

---

## 5. What we believed, and what the data did to it

*(The campaign's epistemic spine. Each row: hypothesis → the test that decided it → what replaced
it. Recommend keeping a condensed version in the thesis — it is the honest shape of the work.)*

| # | We believed… | Decisive test | What's true instead |
|---|---|---|---|
| 1 | The \|λ₀\|≈2.2 spike at T≈6 is a **DQPT** (June) | block PM k=4 | Single-vector PM non-convergence artifact; \|λ₀\| flat ≈0.895; no DQPT for quenches *to* criticality |
| 2 | **Z₂ sector projection** will fix the near-degeneracy (June) | parity analysis | \|X+⟩ boundary is purely even; the odd sector is structurally dead; degeneracy lives *within* the even sector |
| 3 | XXZ's Re-chord failure is **parity oscillations** | staggered-component measurement | Oscillations are null (~10⁻³); the failure is dome inflation at the closing gap |
| 4 | The inflated XXZ dome is a **random Z₂ mixture** | 3-seed test | Seed-independent to 4 digits — a *deterministic* attractor of the truncated iteration; same wall as Alcaraz |
| 5 | The tricritical block PM is **buggy** (oscillating gap) | NB11: dense ground truth + old-vs-new comparison | Eigenvalues were always right (3–4 digit agreement); the wobble was the λ₀-*selector* hopping inside a real band |
| 6 | The **greedy left/right pairing** caused the cluster noise | planted-spectrum dense test | In exact arithmetic old and new pairing are identical; fixes kept as strictly-better, but the wall is physics + truncation |
| 7 | The direct **symmetric-Takagi run** settles dial (iii) — it converges, is seed-independent, and its warm/cold results agree, so "the symmetric route also walls" | operator-vs-tensor-symmetry check (§3b), prompted by the supervisor's stress-test | **Those runs were GAUGE-INVALID.** XXZ's tMPO is a symmetric *operator* (‖U−Uᵀ‖~10⁻¹⁸) but not a manifestly-symmetric *MPO* — its σʸ layer is transpose-antisymmetric and σˣ,σʸ,σᶻ cannot all be symmetrized in any single-site basis. So `powermethod_sym` is structurally unavailable to XXZ; the runs are discarded. A genuine self-correction. |
| 8 | (the real question) **MPO symmetry might extend XXZ's reach** (dial iii) | **construction-independence test** (symmetric-construction propagator through the *two-sided* solver) + spectrum bridge + a TFIM control | Symmetric construction reproduces VD2 to 3–4 digits *including* the T=4→5 inflation; the same Z₂ band sits in the symmetric tMPO; the pure TFIM reaches T≈12 through the *same generic pipeline*. **The wall is physics, not our adaptation of the package** — dial (ii) dominates dial (iii). |

Also refuted along the way: the naive πc/16 guess for the Rényi-2 imaginary offset (measurement:
πc/12, n-independent — row zero of the calibration story), and the June claim "XXZ's gap closes
slower than Alcaraz's" as read from the naive ratio (true only before the symmetry cluster forms;
the *intra-sector* barrier is what's genuinely slower).

**Why row 7 matters most for the meeting.** It is the cleanest illustration of the campaign's
method: a claim that *looked* solid (converged, reproducible, self-consistent) was overturned by
asking "symmetric in *what* sense?" — and the overturning did not weaken the conclusion, it
*strengthened* it (the construction-independence test is a better argument than the gauge-broken
Takagi run ever was). This is worth telling Stefano verbatim.

---

## 6. Exact mapping onto `thesisdraft.md`

Current draft: §1 Intro · §2 ANNNI model · §3 Transverse contraction (3.1–3.4) · §4 MPO of
e^{−iHt} (4.1–4.2) · §5 Non-Hermitian transfer matrix & DPTs (5.1–5.3) · §6 Conclusions ·
Appendix A (TN fundamentals).

| Addition | Goes where | Source material |
|---|---|---|
| The critical XXZ chain & the Néel quench (phase diagram, K(Δ), v(Δ), symmetries, initial-state-as-boundary-state, sublattice rotation, the two-sector structure) | **new §2.2** (ANNNI becomes §2.1) | NB8 §0 + NB9 §0 (written as prose already) |
| The tricritical O'Brien–Fendley chain (model, phase diagram, λ_c program, why S(N/2) overshoots) | **new §2.3** | NB10 intro + equilibrium sections |
| The block power method (full derivation) | **new §5.4** | `blockpm_methods.md` — paste as-is, fix cross-refs (B1)–(B9) |
| Selector-free band reading + partner filter | fold into **§5.2** (eigenvalues vs singular values) | §4.3 above / NB10–NB11 verdicts |
| Symmetric Murg construction (Ising + the XXZ generalization, order-1 insight) | **new §4.3** | §4.4 above + NB12 §1–1b |
| The results chapter: "The entanglement barrier across four models" | **new §6** (Conclusions becomes §7) | `barrier_section.md` — paste as-is; its bracketed cross-refs [tricritical], [gap-closing], [limits] need renumbering |
| Conclusions rewritten around the three-dial map + the honesty narrative | **§7 (was §6)** | §§1–2, 5 of this document |
| Validation protocols & calibrations (dense ground truth, planted spectra, seed tests, Im-S₂ calibration, kernel cross-checks, performance) | **new Appendix B** | §§4.2–4.3, 4.5 above |
| Notebook ↔ section ↔ figure ↔ cache pointer table | **new Appendix C** (half a page) | §8 below |

Suggested writing order (least → most dependent): §2.2, §2.3 (self-contained) → §4.3 → §5.4 →
Appendix B → the results chapter → Conclusions.

---

## 7. Discussion agenda for the supervisor

1. **Framing**: is "the entanglement barrier across four models / three-dial law" the right
   headline for the thesis? (And: is the cross-model map + validation methodology worth a
   methods-oriented paper after the thesis?)
2. **The one untested escape**: *sector-resolved contraction* for exact quench degeneracies —
   classify the degenerate pair by the Z₂ flip operator and contract within one symmetry sector.
   Unlike the June Z₂-projection idea (dead for Alcaraz because its boundary is single-sector),
   XXZ's degeneracy IS inter-sector, so projection could genuinely work there. Worth a chapter
   note as future work, or an actual attempt before September?
3. **The σʸ obstruction and symmetric-MPO tooling** (the meaty methods discussion). We found that
   XXZ has **no single-site-Murg gauge** in which its tMPO is manifestly transpose-symmetric — σʸ
   is antisymmetric, and σˣ,σʸ,σᶻ can't all be symmetrized at once — so `powermethod_sym`/Takagi is
   *structurally unavailable* to XXZ (and to any model needing three anticommuting couplings). This
   is a clean statement about *which models the symmetric machinery is a privilege for*. Two
   follow-on questions: (a) **the one remaining loose end** — the XXZ tMPO IS a symmetric operator,
   so a symmetric MPO representation exists in *some* (non-single-site) gauge; is it worth finding
   it to run a genuinely gauge-correct Takagi solver, or is the construction-independence test
   enough to close the dial-(iii) question for the thesis? (My reading: enough — the band is in the
   symmetric tMPO regardless, so Takagi would face the same ill-conditioned eigenvector; but you
   may want the formal check.) (b) Package notes for Stefano: ITransverse's SymSVD builder is not
   actually left-right symmetric (normdiff 0.07–0.45) — report upstream? And the checker only tests
   the *raw*-swap gauge, silently invalidating `powermethod_sym` on operators that are symmetric in
   a non-trivial gauge — worth a gauge-aware check and/or a Takagi-gauge compressor for the d_t=32
   palindrome.
4. **Remaining sweeps vs deadline**: Δ-sweep (c_eff(Δ) across the Luttinger line), Alcaraz p-sweep
   (locate p* if it exists), dt-convergence bounds. Our reading: none changes the thesis's
   conclusions; all are polish. Agree/deprioritize?
5. **The Δ=1 discrepancy**: c_eff ≈ 0.75 vs 1 at the Heisenberg point — we attribute it to
   marginal-operator logs (same pathology as equilibrium fits there) plus the ~10% calibration
   systematic. Is a dedicated check (e.g., Δ=0 free-fermion benchmark, exactly solvable) worth it?
6. **How much of the refutation narrative (§5) belongs in the thesis** vs stays in the notebooks?
   Our proposal: a condensed table in the Conclusions + full detail in Appendix B.

---

## 8. Quick-reference appendix

**Central-charge results (all with their estimator and window):**

| model | c (temporal) | estimator / window | c (equilibrium, DMRG) |
|---|---|---|---|
| Ising p=0 | 1/2 reproduced (Re chord + Im ≈ π/24) | n→1 Takagi, T≤14 | 0.50 |
| Alcaraz p=0.1 | **0.47 ± 0.05** | Rényi-2 chord, T=4..9 | ≈0.5 for p≤2 |
| XXZ Δ=0.5 | ≈0.9–1.0 | Im S₂ (calibrated), clean windows | ≈1 |
| XXZ Δ=1.0 | 0.75 ± 0.05 (log-suspect) | Im S₂ (calibrated), T=3..8 | ≈1.05 (overshoot, marginal logs) |
| tricritical λ≈0.42 | — (no window) | — | consistent w/ 7/10 (energy method); S(N/2) overshoots |

**Validation accuracies (NB11/NB12):** block PM vs exact dense: 1.4e-8 (Alcaraz VD2), 4.6e-13
(Alcaraz WII), 4.2e-10 (Ising Murg), 2.3e-10 (tricritical). Regression vs master sweep: 2e-5 (T=3),
2e-3 = pair-splitting scale (T=6). Murg-XXZ echo vs TDVP: 2.6e-5 (order 2), 1.1e-5 (order 1);
VD2 reference 4.0e-5; WII 2.6e-3. Seed spreads: 0.0000 (symmetric, both T); 4-digit identity
(asymmetric).

**Notebook inventory (1–12):**

| NB | one-line role | key figures |
|---|---|---|
| 1 | model intro; naive TEBD vs TDVP; barrier first seen | — |
| 2 | exp-MPO kernels WI/WII/VD2, Schrödinger benchmark | — |
| 3 | temporal Rényi-2 intro; single-vector → block PM; dt appendix | — |
| 4 | equilibrium c(p) ≈ 1/2 (DMRG, 3 reads) | p_dependence, cft_L |
| 5 | spectrum & limits: block PM, dual unitarity, THE WALL | block_pm_ising_vs_alcaraz, alcaraz_gap_dualunitarity, gap_closing_wall |
| 6 | Ising control: Eq.(6) + λ₀ circle, T=14 | cft_ising_validation, ising_lambda0_circle |
| 7 | **Alcaraz headline** c=0.47±0.05 | temporal_entropy_profiles, temporal_chord_fit |
| 8 | XXZ model map (§0) + Néel quench validation | xxz_c_equilibrium |
| 9 | XXZ temporal result: Z₂ story, Im-c, null oscillations, intra-sector gap | xxz_entropy_profiles, xxz_oscillations, xxz_vs_alcaraz_gap, xxz_intrasector_gap |
| 10 | tricritical: λ_c located; charge-driven band (corrected reading) | tricritical_c_equilibrium, tricritical_gap |
| 11 | block-PM validation (V1–V4) + performance §P | — |
| 12 | the asymmetry experiment — resolved: dial ii > dial iii | — |

**The two paste-ready drafts**: `barrier_section.md` (results chapter, ~150 lines) and
`blockpm_methods.md` (methods §5.4, ~210 lines) — both in the draft's voice, equations numbered
locally, placeholders marked with brackets where thesis cross-references must be filled in.
