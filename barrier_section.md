# The entanglement barrier across four models

*(drop-in thesis section, written July 2026 in the style of thesisdraft.md; companion to the
methods section in blockpm_methods.md. Sources: NBs 5–7 (Alcaraz/Ising), 8–9 (XXZ), 10
(tricritical), 11 (method validation), 12 (the symmetric-XXZ experiment, marked [P2] where its
outcome is pending).)*

## 1. One phenomenon, four disguises

Every model studied in this thesis was quenched to a critical point, and in every one of them the
transverse contraction eventually failed in the same characteristic way: the temporal-entropy dome,
after growing conformally for a while, suddenly inflated beyond any CFT chord, while the transfer
matrix's leading eigenvalues collapsed onto a band of nearly equal moduli. We call the underlying
phenomenon the **entanglement barrier**, and this section argues three things: that it is *one*
phenomenon appearing in four disguises; that its position in each model is fixed by physics —
central charge, symmetry content of the quench, and the symmetry of the MPO — and not by any
algorithmic deficiency; and that within its limits, universal information remains extractable.

The mechanism, in the language of Section 5 [non-Hermitian transfer matrix]: after a quench *to* a
critical point, conformal invariance drives the transfer matrix toward a rescaled unitary — emergent
dual unitarity. A unitary has no modulus hierarchy: all eigenvalues share $|\mu|$. The gap ratio
$|\mu_1/\mu_0|\to1$ is therefore not a numerical accident but the *defining approach* to the
conformal regime. The catch is a basic fact of non-Hermitian linear algebra: while eigen*values*
remain well-conditioned as the gap closes, the individual eigen*vectors* — precisely what the
generalized entropies are built from — degrade with condition number $\sim1/\mathrm{gap}$. The
better the conformal physics, the less well-defined the object we want to measure. The barrier is
thus a physics/numerics double image: the *same* gap closing that confirms dual unitarity destroys
the eigenvector.

A methodological point established in this work (notebook 11): this attribution is not a guess. The
block power method was validated against **exact dense diagonalization** of small transfer matrices
on all four models (agreement $10^{-8}$–$10^{-13}$), and the improved solver reproduces every
converged point of the earlier campaigns. When the entropy breaks, it is not the solver.

## 2. The four manifestations

**Ising (symmetric MPO, $c=1/2$) — the control, reach $T\approx14$.** The integrable Ising chain
enjoys every advantage at once: the Murg construction gives a left–right *symmetric* tMPO, enabling
the Autonne–Takagi diagonalization of the RTM — far better conditioned near degeneracy than any
two-sided method; the $|X^+\rangle$ quench is $\mathbb Z_2$-symmetric, so the boundary fixed point is
unique; and $c=1/2$ closes the gap slowly. Result: clean reproduction of the CFT predictions
(Re $S$ on the $c=1/2$ chord, Im $S\to\pi c/24\cdot2=\pi/24$, the $\lambda_0$ circle) out to $T=14$.
This is the benchmark every other model degrades from.

**Alcaraz / ANNNI-type (asymmetric MPO, $c=1/2$, frustrated) — the dynamical wall, reach
$T\approx10$.** Losing MPO symmetry forces the two-sided power method; the NNN frustration
accelerates the gap closing relative to Ising. The near-degeneracy here is *dynamical*: it develops
gradually as the barrier forms (partner-filtered gap $0.35\to0.98$ over $T=1..6$), leaving a wide
clean window ($T=4..9$) from which the headline of this thesis was extracted —
$c(p{=}0.1)=0.47\pm0.05$, temporal Ising universality surviving NNN frustration. Past $T\approx10$
the dome inflates irrecoverably; every attempted repair (subspace projectors, continuity tracking,
seed strategies) fails for the reason of Section 1, and the eigenvalue route becomes ambiguous once
three moduli agree to $\sim10^{-2}$.

**XXZ Néel quench (asymmetric MPO, $c=1$, exact $\mathbb Z_2$ degeneracy) — the symmetry wall, reach
$T\approx4$.** The XXZ chain has no transverse field, so every symmetric product state is nearly an
eigenstate and generates no temporal entanglement: one is *forced* to quench from the
symmetry-breaking Néel state, whose two translation-related copies hand the transfer matrix an
**exactly degenerate** leading pair from $T=0$. At small $T$ a finite cat-state splitting keeps one
combination safely dominant (the spectrum at $T=1$ is $[0.98, 0.35, 0.35, 0.20]$ — note the
degenerate *subleading* pair, the fingerprint of the two sectors); as $T$ grows the splitting decays,
and between $T=3$ and $4$ a **4-fold band** (two sectors × their $\pm$ partners) locks in within 3%.
Exactly there the Re dome jumps discontinuously (peak $0.25\to0.81$ at $\Delta=1$) and chord fits
return nonsense. A seed test shows the inflated fixed point is *deterministic* (identical to 4 digits
across independent random seeds) — the truncated, nonlinear iteration has a unique attractor that is
no longer a pure eigenvector; no seeding or sector-picking strategy can recover it. Two features
survive the wall: the imaginary offset Im $S_2\to\pi c/12$ (calibrated on the $p=0$ Ising line, where
it is empirically $n$-independent) gives $c_{\rm eff}\approx0.75$ at $\Delta=1$ (stable over $T=3..8$;
marginal-operator logs at the Heisenberg point are the suspected bias) and $\approx0.9$–$1.0$ at
$\Delta=0.5$; and the **intra-sector barrier**, resolved with a $k=6$ block as $|\theta_5/\theta_1|$,
closes *far slower* than Alcaraz's ($0.14\to0.75$ vs $0.35\to0.97$ over $T=1..5$) — confirming the
"nearest-neighbour, unfrustrated ⇒ slow barrier" expectation. The reach is destroyed not by the
barrier but by the symmetry cluster sitting on top of it: **for the transverse method, an exact
symmetry degeneracy of the quench is a worse enemy than frustration.**

**Tricritical / O'Brien–Fendley (asymmetric MPO, $c=7/10$) — the charge wall, essentially no
window.** At the tricritical point the larger central charge closes the gap fastest of all: exact
diagonalization already shows the subleading band clustering at $T=0.2$; by $T=3$ *five* Ritz values
agree within 15%; by $T=4.5$ the top four agree within 4% and the old iteration finally stuck — and
the *fixed*, ground-truth-validated solver, cold-started at that same point, ran 7.8 hours without
converging before being terminated: the wall is not negotiable by better internals. The earlier
suspicion that the noisy tricritical gap data
reflected a solver bug was disproved point by point: the fixed solver reproduces the old eigenvalues
to 3–4 digits wherever both converge, and the apparent $|\lambda_0(T)|$ oscillation is a *selector*
artifact — the continuity rule for "the physical $\lambda_0$" hopping between band members that
genuinely are equivalent. Inside a band, "the dominant eigenvector" is not merely ill-conditioned;
the question itself loses meaning.

## 3. The unified picture

| model | MPO symmetry | $c$ | quench symmetry | degeneracy type | clean reach |
|---|---|---|---|---|---|
| Ising | symmetric (Takagi) | 1/2 | symmetric ($|X^+\rangle$) | none until barrier | $T\approx14$ |
| Alcaraz $p=0.1$ | asymmetric | 1/2 | symmetric ($|X^+\rangle$) | dynamical (barrier) | $T\approx10$ |
| XXZ Néel | asymmetric | 1 | **breaking** (Néel) | **exact, from $T=0$** | $T\approx4$ |
| tricritical | asymmetric | 7/10 | symmetric ($|X^+\rangle$) | band (charge-driven) | $\lesssim2$ |

Three independent dials set the wall: **(i) central charge** — larger $c$ means more temporal
entanglement per unit $T$, a faster-closing gap, an earlier band (tricritical worst); **(ii) quench
symmetry content** — a symmetry-breaking initial state contributes an exact degeneracy that no
window precedes, only a decaying cat splitting masks (XXZ); **(iii) MPO symmetry** — a symmetric
construction unlocks Takagi diagonalization, whose conditioning near degeneracy is qualitatively
better (Ising's entire advantage over Alcaraz at equal $c$).

## 4. Is there a way out?

Tested and failed (each for the structural reason of Section 1, not for lack of tuning):
subspace-projector entropies (add classical mixing entropy), continuity projection across $T$
(drifts), sector selection by seed (the attractor is unique), better power-method internals (the
solver was never wrong — validated against exact ground truth), larger blocks (the band grows
faster than $k$).

Partial escapes, used in this thesis: extract $c$ **before** the wall (Alcaraz's clean window); use
the **imaginary offset** of $S_2$, which is first-order insensitive to the eigenvector breakdown
(XXZ); use **eigenvalue-route** observables — the physical gap, $|\lambda_0|$ flatness, the
$\lambda_0$ circle, the boundary exponent $x_1$ — which remain well-conditioned arbitrarily deep
(all models).

The decisive escape experiment (notebook 12): a **symmetric MPO for XXZ**. The rotated two-site
term decomposes as $S^xS^x - S^yS^y - \Delta S^zS^z$, three mutually commuting layers each
admitting an *exact* bond-2 symmetric (Murg-type) factorization, so a reflection-symmetric
second-order propagator $e^{ZZ/2}e^{YY/2}e^{XX}e^{YY/2}e^{ZZ/2}$ exists — unlike the package's
SymSVD attempt, which is demonstrably not symmetric. **Built and verified**: the propagator
reproduces the exact one-step evolution of a small chain, and its Néel-quench echo matches the
independent TDVP benchmark to $2.6\times10^{-5}$ — the same accuracy as the production VD2 kernel.
An independent WII cross-check of the *asymmetric* pipeline (cheaper and, for this
strictly-nearest-neighbour model, still second order) reproduces the notebook-9 story point for
point, including the eventual wall — a second, independent confirmation that dials (i)+(ii)
alone, without touching MPO symmetry, already explain the asymmetric result.

**Isolating dial (iii) itself gave an inconclusive but suggestive result.** Feeding the verified
symmetric propagator into `powermethod_sym` + the Autonne–Takagi $n\to1$ entropy, cold-started
independently at each $T$ (no warm-starting across the ladder — an implementation gap relative to
the asymmetric sweeps, left as follow-up), the extracted $c$ **does not stabilize**: $1.50, 0.25,
0.96, 1.21, 0.48$ over $T=2..6$, with no discernible trend, even though the bond dimension stays
modest throughout ($\chi=5\to15$, comparable growth to the asymmetric case, no blowup). Two
readings are possible. Either the noise is purely a seeding artifact and a warm-started rerun
would recover a clean signal — in which case the question of whether Takagi conditioning moves the
wall remains genuinely open; or, more suggestively, the *erratic, cold-start-dependent* result is
itself informative: unlike Ising, whose symmetric fixed point is unique and reproducible from any
seed, XXZ's symmetric route lands on a *different* answer at every independent draw — exactly the
signature expected if the exact $\mathbb Z_2$ Néel degeneracy persists **regardless of MPO
symmetry**, i.e. dial (ii) dominates dial (iii) here. This tentative reading is consistent with
the rest of this section's argument (the degeneracy is a property of the quench, not the
contraction scheme) but was not confirmed with a controlled warm-started sweep before this
thesis's numerical campaign closed; it is the natural next experiment.

Finally, the honest framing. Emergent dual unitarity — the physics this method was built to see —
*is* the closing of the gap. A contraction scheme whose accuracy requires an open gap therefore
carries its own horizon: the barrier is not an obstacle in front of the physics; it is the physics,
seen from the numerical side. Mapping where the horizon sits for each model, and why, is this
thesis's methodological contribution.
