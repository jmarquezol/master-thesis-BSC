# Appendix rewrite — context handoff

Written 2026-08-18 for a fresh Claude session whose task is to continue rewriting the appendices
(`thesis/appendix.tex`, sections A–G) in the style the author has established across the main
text. Read this file fully before touching anything. The author (Joaquín) reviews and approves
every change; nothing lands without his sign-off.

---

## 1. The task and the workflow

Rework `appendix.tex` section by section to the target style, the way §6 and §7 were reworked:

1. The author marks passages with inline bracket comments `[...]` in the .tex — these are the
   work queue. If a section has no comments, ask before restyling it wholesale.
2. Analyse the comments, then propose the changes INLINE IN CHAT (quote current text, show
   proposed text, one block per change). Wait for his approval.
3. Only after approval, implement. He sometimes answers with adjustments — implement those.
4. After every batch of edits: run the integrity checks and compile (recipes in §6). Report
   plainly what passed.

Two hard process rules learned this session:

- **The tree drifts.** The author edits the .tex files continuously (LaTeX Workshop, real time).
  Always re-read the current file immediately before proposing or applying an edit. Never edit
  from a stale read. His rewording of your earlier text is final — do not "fix it back".
- **Edit mechanics.** Match on long unique strings and assert exactly one match before writing.
  Never run a regex spanning `\caption{` to `\label{...}` with dot-matches-all (it matches from
  the first caption in the file and has destroyed a section once). Never blanket-replace short
  numeric strings ("0.3" also lives inside `aspect=0.36`).
- **His saves silently revert your edits.** He keeps `appendix.tex` open in VS Code; when he saves
  from a buffer loaded before your write, your changes vanish with no error on either side. This
  happened five times in the 2026-08-19/20 session, each time reverting a different subset. A clean
  compile is **not** evidence your edit is on disk. After every batch verify with `grep -F` on a
  distinctive string from each change plus an `md5`, and when he reports one thing missing, audit
  *all* recent changes — the lost subset is never only what he noticed.
- **Deleting a sentence strands its bracket comment.** Twice a comment was left behind attached to
  the previous sentence when the sentence it annotated was removed. When you delete text in a
  commented region, check for a trailing `[`.
- **Sweep for comments inline, not just at line start.** `grep -nE '^\s*\['` misses comments in the
  middle of a paragraph. Use a regex over the whole line and filter LaTeX false positives
  (`\includegraphics[...]`, `\sinh[...]`, `\mathrm{round}[...]`, `\documentclass[...]`).

Claude never commits or pushes. The author owns the repo history.

## 2. Style — where it is defined and what it means in practice

Authoritative documents, in order: `GEMINI.md` (full 8-section guide, repo root),
`.agents/rules/thesis_style.md` (condensed same rules), `CLAUDE.md` §2.1 (extra prose rules).
Gold-standard reviewed text to imitate: `thesis/main.tex` (abstract, intro, conclusions),
`thesis/cft.tex`, `thesis/P3.tex`, `thesis/P2.tex`, `thesis/P4.tex`, `thesis/numerics.tex`,
and `thesis/results.tex` §7 (reworked 2026-08-18 and then hand-edited by the author — his
current wording there is the freshest sample of his voice).

The rules that actually bite, distilled from his review feedback:

- British English. "We" for aims, choices, methods, results. Plain, precise, confident.
- Paragraphs: one purpose each; topic sentence; Problem → Mechanism → Consequence → Relevance.
- No bold or `\emph` in prose. (Figure/table captions DO keep the established bold lead-in:
  `\caption{\textbf{Short title.} Sentence...}` — that is the house convention, keep it.)
- Very few em-dashes; when used, `---` with no surrounding spaces. En-dash for ranges: `2--20`.
- No meta-commentary ("it is worth noting", "we stress", "Crucially/Notably/Importantly"),
  no cleft constructions ("What fails is..."), no rhetorical section titles.
- Banned AI-isms (full list in GEMINI.md §7): "the bridge is", "bring the theory to bear",
  "not a free lunch", "the signal is spoiled", "delve", "testament to", etc. Also banned by
  review experience: "rather than assumed in order to...", stacked emphasis words
  ("also", "exactly", "precisely"), dense multi-idea sentences.
- Citations: attach `~\cite{key}` to the claim as it is made ("as done for the integrable
  chain~\cite{carignano_tagliacozzo2025}", "was demonstrated in Ref.~\cite{...}").
  NEVER make the reference the subject: "The comparison is the one performed in
  Ref.~\cite{...}" was explicitly rejected.
- Repetition is the most common complaint. If a mechanism is explained in the main text, the
  appendix refers back ("for the reason given in Section~\ref{...}") — it does not re-explain.
  Conversely the conclusions must be self-contained: main ideas stated, not deferred to
  appendix pointers.
- Every claim carries a figure, a table, or a printed number; otherwise it is cut or hedged as
  open. Honesty over confidence: what is measured vs inferred is separated explicitly; open
  problems are stated as outlook, never papered over with weak arguments.
- No "pending the cluster" hedges in final text. A number that must come from the cluster is
  the macro `\OPEN` (defined in main.tex) with the surrounding prose written as final.
- LaTeX conventions: `\ac{...}`/`\acp{...}` for acronyms (list in main.tex ~line 77),
  `Eq.~\eqref{...}`, `Eqs.~\eqref{...} and~\eqref{...}`, `Section~\ref{...}`,
  `Figure~\ref{...}`, `Table~\ref{...}`, `Appendix~\ref{...}`, `Ref.~\cite{...}`.
- Words: "frustration" is banned in the thesis (author dislikes it) — say "the \ac{NNN}
  interaction/coupling". The π-displaced tower is never mentioned anywhere — do not reopen.
- Notebook markdown/code style rules (if you touch NBs/): CLAUDE.md §2.2–2.3.
- **Depth is free, density is not.** Appendices have no length limit and every asserted choice must
  say what was tried and what the alternatives gave, with numbers, so a viva jury has no unanswered
  question. But he cuts hard when the same idea is said twice or when a paragraph is a list of
  numbers: prefer one number that carries the argument over four that decorate it, and state the
  conclusion the numbers are for. Both complaints arrived in the same session.
- **Do not overclaim to fill a gap.** Where an identification is not verified, say so
  ("We do not identify the operator responsible for it"). He chose that over naming it. Equally, do
  not justify a fit parameter with physics the thesis never introduces — an appeal to the Ising bulk
  operator content was cut for exactly this reason.
- `Ref.~\cite{...}` is acceptable only as the object of a preposition ("demonstrated in
  Ref.~\cite{...}"). "the Casimir formula of Ref.~\cite{...}" was rejected twice. Safest form is to
  attach the citation to the claim with no "Ref." at all.

## 3. Provenance rule for numbers

Main-text numbers come from cluster runs only; missing ones are `\OPEN`. Appendices MAY use
local test runs. Local caches live in `analysis/session_caches/`; cluster caches in
`results/data/cluster/` (only `_bulk` suffixed files are on the corrected transfer-matrix
column — legacy files are not to be quoted). Never write a conclusion before the code has run
and you have seen its printed output; quote the numbers a script actually printed.

## 4. Appendix map and current state (2026-08-20)

The order was changed on 2026-08-20 so the appendices follow the order the main text first
cites them. Lettering is now:

| | section | main labels | first cited from | state |
|---|---|---|---|---|
| A | Conformal-field-theory derivations | `app:cft`, `app:cft:maps`, `app:cft:spectrum` | §2 `cft.tex` | **reworked** |
| B | Truncation schemes | `app:tensor_networks`, `app:sub:rdm`, `app:sub:rtm_truncation` | §3 `P3.tex` | **reworked** |
| C | Why the NNN coupling breaks integrability | `app:pairing`, `sec:annni:jw` | §4 `P2.tex` | **reworked** |
| D | Equilibrium measurements | `app:equilibrium`, `app:fss`, `app:velocity` | §4 `P2.tex` | **reworked** |
| E | Construction of the time-evolution MPO | `app:mpo` (+ code, outline, vd2, comparison) | §5 `P4.tex` | not yet reworked |
| F | Implementation of the block power method | `app:blockpm` (+ pairing, schur, validation, cost, library, selector, conditioning, failures) | §6 `numerics.tex` | not yet reworked |
| G | Fitting conventions and supporting profiles | `app:extraction`, `app:x1`, `app:highp`, `app:errors` | §7 `results.tex` | not yet reworked |

**Deleted**: the old "Universal versus non-universal constants" appendix. Both of its results
were unused — the entropy calibration `c(p) ≈ ½ slope(p)/slope(0)` (the thesis now takes `c`
from the imaginary plateau) and the gap-ratio trick (every `x_i` comes from Eq. (12) with the
measured `v_∞`). Its one surviving idea, that the non-universal normalisation cancels in phase
differences, is now in `app:cft:spectrum` as the `fβ + f_s` terms and in a self-contained
sentence closing §2 of `cft.tex`. Do not reinstate it.

**Equation numbering**: `\numberwithin{equation}{section}` sits immediately after `\appendix`,
so appendix equations are A.1, B.3, … while main-text equations keep plain numbers.

## 4b. What was done in the 2026-08-19/20 session

Appendices A–D were reworked. Highlights that must not be undone:

**A — CFT derivations.** Fixed a factor-of-two error: the twist operator was labelled a scaling
dimension `Δ_n = (c/12)(n−1/n)` but used with exponent `ℓ^{−4Δ_n}`, which contradicted the next
equation. It is now `x_n` with exponent `ℓ^{−2x_n}`, consistent with the appendix's own two-point
function. The entropy coefficient is now derived end to end (replica → twist operator → bulk
`2x_n` → edge-anchored `x_n` by images → chord `W` → continuation `ℓ→ivT`, which produces the
`iπc/24(1+1/n)` term of `eq:cft:sgen` and hence the `πc/16` plateau). The strip quantisation is
derived (Schwarzian → `⟨T⟩ = −π²c/24β²` → Casimir → `Ĥ = (π/β)(L₀ − c/24)`), and `E = e^{−aĤ}`
is verified rather than asserted (`μ_i = exp[−πg(x_i − c/24)]` reproduces `eq:cft:transfer`).
Chiral sector, primary operator and spin `s = h − h̄` are defined for a non-expert reader.
`cft.tex` Eq. (11) gained the missing `B/T` term, which its own "in practice we fit" sentence
already implied.

**B — Truncation schemes.** Restructured: intro (Schmidt + Eckart–Young–Mirsky, deliberately
uncited), B.1 the RDM scheme, B.2 the RTM scheme and its limits. The standalone gauge/canonical
-forms subsection was removed; the gauge argument is now derived at tensor level
(`A^{[n]} → A^{[n]}X`, conjugate carries `X†` → congruence; independent `⟨L|` carries `X⁻¹` →
similarity), opening with the invariance of `⟨L|R⟩`. **A wrong claim was fixed**: the appendix
said the two schemes have "comparable computational cost", contradicting P3's "substantially
cheaper". It now says the per-truncation cost is the same at a given bond dimension and the
saving comes from the RTM needing a smaller `χ`, with the RTM rank bounded by the operator
entanglement (`carignano_marimon_tagliacozzo2024`, verified).

**C — integrability.** Retitled "Why the NNN coupling breaks integrability", subsection removed,
duplication removed, and the Jordan–Wigner algebra now shown step by step in five equations. The
key sentence is that NN and NNN are the *same* calculation, differing only in that the string
reaches two sites further and leaves two factors instead of one.

**D — Equilibrium measurements.** The mid-chain fits now include the finite-size correction and
the **corrected values are the quoted ones**: `c = 0.500` (von Neumann) and `0.508 / 0.505`
(Rényi-2), against `0.511/0.514` and `0.584/0.588` uncorrected. The velocity subsection was
rewritten to describe the routine that actually produces `tab:velocity` (`ground_and_gap`), which
uses **momentum only** — parity is a verified check, not a selection; the `1/N²` form is justified
against `1/N` and `1/N³` at `p=0`; the XY-chain validation is the main argument, with its
Hamiltonian, dispersion and `tab:xy`; and the Casimir check shows its formula.

**Bibliography**: four verified additions — `cardy1984`, `blote_cardy_nightingale1986`,
`cardy2004bcft`, `cardy_calabrese2010`. Each was checked by reading the paper, not the abstract.
Two traps found this way, both of which had been wrong in the thesis or nearly written into it:
`cardy1984` contains **no** strip Hamiltonian (that is `cardy2004bcft` Eq. 27), and `calabrese2010`
is about **Luttinger liquids** with exponent `K/2α` — it does not support the `ℓ^{−x/n}` correction
used in `app:fss`, which is `cardy_calabrese2010`. `calabrese2010` survives only in P2, where it
supports a different claim. If that citation reappears in `app:fss`, it is a regression.

## 5. COORDINATION — freeze lifted (2026-08-20)

The failure-mechanism investigation is complete and all its thesis edits have landed:
`app:failures` (seed-ensemble, RDM and dense-conditioning tables), the `app:conditioning`
rewrite, the validation-table caveat and corrected-column row, the compressed production
options, the cutoff-control paragraph in `app:errors`, and the corresponding main-text
changes in `numerics.tex`, `P3.tex`, `results.tex` and the conclusions. Every appendix
section is now fair game for the style rework. Do not weaken the measured statements in
`app:failures`/`app:conditioning`/`app:errors` when restyling: the numbers and their
hedges were set deliberately (what is demonstrated vs open).

## 5b. Figures, and the scripts behind them

You can **read PNGs directly** with the Read tool — use it to look at a figure before and after
changing it. PDFs cannot be rendered (no poppler), but PNGs display. Three figure defects in this
session were only caught by looking: a clipped legend, colliding tick labels, and a missing panel
label the caption already referred to.

House conventions, verified across every `analysis/session_scripts/fig_*.jl`:

- `msw=0` on every marker (no stroke). `fig_cft_L.jl` was the only script missing it and its
  markers looked different from the rest of the thesis.
- `framestyle=:box` — already set globally by `thesis_plot_theme!`, so setting it again is a no-op.
- Canvas from `thesis_size(frac; aspect)`; **shrink the canvas, never raise the font**.
- Multi-panel figures put the legend in **its own subplot column** (`framestyle=:none`,
  `legend=:left`, entries drawn from `[NaN]` points, layout `@layout([a b c{0.16w}])`). An
  `:outerright` legend steals width from its own panel and makes panels unequal.
- Panel labels as `title="(a)", titlelocation=:left, titlefontsize=11`.
- Every thesis figure must regenerate from its owning notebook cell. When you edit a figure,
  update **both** the session script and the notebook cell, and keep them identical.

Scripts written this session, all runnable with `julia --project=../..` from
`analysis/session_scripts/`:

| script | what it does |
|---|---|
| `fss_corrections.jl` | the mid-chain correction fits behind D's `c = 0.500 / 0.508 / 0.505` |
| `fss_plotcheck.jl` | how much the correction moves the plotted curves (answer: <10⁻³) |
| `fig_cft_L.jl` | rebuilds `cft_L.png` (Fig 17), two panels + legend column |
| `fig_velocity_extrap.jl` | rebuilds `nb4_velocity_extrapolation.png` (Fig 18), two panels |
| `xy_check.jl` | recomputes the XY validation in `tab:xy` from `nb7_xy_control.jld2` |

**Reading papers.** Sci-Hub is not used. When a paper is paywalled, ask him — he drops the PDF on
his Desktop. There is no `pdftotext`; extract with Python by zlib-decompressing the PDF streams and
pulling the strings out of the `Tj`/`TJ` operators. That worked for every paper this session.

## 6. Build and checks

Compile (no system TeX; tectonic only):
```bash
cd thesis
PATH="$HOME/.venvs/thesis/bin:$PATH" tectonic -X compile main.tex --outfmt pdf -Z shell-escape --keep-intermediates
```
`--keep-intermediates` matters: without it main.aux on disk is stale and page numbers read
from it are wrong. minted needs `-Z shell-escape` and pygmentize (venv). LaTeX Workshop in the
IDE compiles with tectonic on save — the author watches the PDF live while you edit.

Integrity after every batch (from `thesis/`):
```bash
for r in $(grep -ho 'ref{[^}]*}' *.tex | sed 's/ref{//;s/}//' | sort -u); do
  grep -qr "label{$r}" *.tex || echo "BROKEN REF: $r"; done
for k in $(grep -ho 'cite{[^}]*}' *.tex | sed 's/cite{//;s/}//' | tr ',' '\n' | tr -d ' ' | sort -u); do
  grep -q "{$k," biblio.bib || echo "UNDEF CITE: $k"; done
grep -ho 'label{[^}]*}' *.tex | sort | uniq -d
```
Plus begin/end balance and even `$` counts per changed file. If you delete a label or
equation, grep the whole tree for references to it first (an `eq:res:c_entropy_p0` deletion
once orphaned an appendix reference).

## 7. Physics conventions you will meet in the appendices

- Model: self-dual ANNNI-type chain, λ=1 critical; field on σ^x, order parameter on σ^z;
  initial state |X+>^N (free boundary, x1=1/2); |Up> is the fixed boundary (x1=2).
- δt=0.1, N_β=4 (β₀=0.2), VD2 MPO in production; temporal chain length T/δt+N_β.
- Chord variable W defined once at `eq:chord` in cft.tex — use W, never redefine.
- Temporal cuts are boundary cuts: chord slope c/8 at n=2, not c/3. Im plateau target πc/16.
- λ_i=log(−μ_i): constant −π branch offset, pinned in Eq.(3) fits.
- Entropy profiles: trim N_β/2 bonds per end before fitting (cluster caches are pre-trimmed
  by `trim_dome`; local svpm caches are NOT — trim `[3:end-2]` at nbeta=4).
- Couplings in §7: p = 0, 0.1, 0.3, 0.5 only. Equilibrium/velocity material runs to p=1.5 —
  do not restrict it.
- No error bars anywhere (follows Bou-Comas practice); windows and controls in `app:errors`
  play that role.

## 7b. Open items carried into the next session

Small, all verified as still present on 2026-08-20:

- **Duplicated paragraph.** The sweep-locality paragraph (local environment at the current stage,
  translation invariance regenerating discarded components) appears in both Appendix B and
  Appendix F. A bracket comment in F flags it. He said leave both and decide when F is revised.
- **`app:highp` caption is stale**: it ends "The $p=1.0$ profiles use the halved Trotter step"
  while `fig:domes_hi` shows $p=0$, $0.3$ and $0.5$.
- **Two stale cross-references.** Appendix G opens with "the two boundary conditions of
  Section~\ref{sec:cft:temporal}"; the boundary conditions are introduced in §2.1,
  `sec:cft:boundary`. And `defense_prep.ipynb` still says "Appendix B reaches the same formula",
  which is Appendix A now.
- **Optional, needs his approval** (he has seen the proposal and not decided): a clause in P2
  linking the offset in the fixed-`N` profile fit to the correction measured in `app:fss`, worded
  so it is clear the 0.511 → 0.500 shift belongs to the *mid-chain* fit and does not correct the
  main text's 0.502.

## 8. Quick orientation to the rest of the repo

`CLAUDE.md` (repo root, gitignored) is the full project state file — long, but §0 (working
state), §2 (writing rules) and §6 (practical notes) are the parts that matter for this task.
`cluster/README.md` describes the cluster queue. `src/transverse_tools.jl` is the library the
appendix's algorithm text describes (block_transfer_eigs, run_pm_diagnosed,
run_pm_consensus). The thesis is edited both here and in Overleaf: if the author says "I am
editing in Overleaf", give replacement text inline in chat instead of touching files, and if
he uploads `BSC-TFMvN.zip`, that zip is canonical — sync before editing (procedure in
CLAUDE.md §6).

---

## 9. Starting prompt for the fresh session

Paste this as the first message:

> Read `APPENDIX_REWRITE.md` in the repo root first, then `GEMINI.md`,
> `.agents/rules/thesis_style.md` and section 2.1 of `CLAUDE.md`. Then read `thesis/cft.tex`,
> `thesis/P3.tex` and `thesis/results.tex`, and skim Appendices A to D of `thesis/appendix.tex`,
> which were reworked to the target style and are the closest sample of where the rest should land.
>
> Your task is to continue the appendix rework in `thesis/appendix.tex`. Appendices A to D are
> done. Still to do, in order: **E** (Construction of the time-evolution MPO), **F** (Implementation
> of the block power method) and **G** (Fitting conventions and supporting profiles). Section 4b of
> the handoff lists what was changed in A to D and what must not be undone.
>
> Work the way section 1 of that file describes. Analyse the whole section first, then propose the
> changes to me inline in chat, quoting the current text and your proposed replacement, one block
> per change, and wait for my approval before editing any file. I mark passages I want changed with
> inline bracket comments `[...]`; if a section has none, tell me and propose what you would change
> rather than restyling it wholesale.
>
> Three things I care about, learned the hard way in the last session. Check that the text describes
> the code that actually produced the numbers, not a neighbouring routine — that error was found
> twice. Verify a citation by reading the paper before you attach it to a claim; if it is paywalled,
> ask me and I will put the PDF on my Desktop. And after every batch, verify your edits are still on
> disk with `grep -F` and an md5, not just a clean compile, because my editor sometimes reverts them.
>
> Start with Appendix E and tell me what you find before proposing anything.
