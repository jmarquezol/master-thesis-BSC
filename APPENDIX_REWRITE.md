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

## 3. Provenance rule for numbers

Main-text numbers come from cluster runs only; missing ones are `\OPEN`. Appendices MAY use
local test runs. Local caches live in `analysis/session_caches/`; cluster caches in
`results/data/cluster/` (only `_bulk` suffixed files are on the corrected transfer-matrix
column — legacy files are not to be quoted). Never write a conclusion before the code has run
and you have seen its printed output; quote the numbers a script actually printed.

## 4. Appendix map and current state (2026-08-18, ~583 lines)

| section | label | lines (approx) | state |
|---|---|---|---|
| Conformal maps + entropy coefficients | `app:cft` (`app:cft:maps`, `app:cft:spectrum`) | 6–79 | not recently reviewed |
| Universal-data extraction | `app:universal` / `sec:numerics:universal` | 81–112 | not recently reviewed |
| TN gauge/truncation | `app:tensor_networks` (canonical_forms, truncation, rtm_truncation) | 114–211 | not recently reviewed |
| Pairing / JW | `app:pairing`, `sec:annni:jw` | 216–257 | not recently reviewed |
| Extraction + equilibrium | `app:extraction` (`app:fss`, `app:velocity`, `app:highp`, `app:errors`) | 259–~370 | partly reworked, see below |
| MPO construction | `app:mpo` (code, outline, vd2, comparison) | ~375–455 | not recently reviewed |
| Block PM details | `app:blockpm` (pairing, schur, validation, cost, library, selector, conditioning) | ~460–end | see freeze list |

Recent changes already made (2026-08-18), do not undo:
- `app:errors`: the Re-vs-Im estimator paragraph and `tab:reim` were DELETED (user decision).
  A "fourth control" cutoff paragraph was REMOVED pending the author's review of fresh test
  output — its text is archived at the investigation session's scratchpad. DO NOT re-add it.
- `app:extraction` gained an opening paragraph on the free/fixed boundary exponents
  (x1=0.498 / 1.996) moved out of §7.2, above `fig:x1_p0`.
- The branch-constant sentence in `app:errors` says "at every coupling".

## 5. COORDINATION — sections frozen for the investigation session

A parallel session is running a failure-mechanism investigation (seed ensembles, RDM swap,
cutoff rerun, dense corrected-column test) and will, after the author reviews its dossier,
edit these places itself:

- `app:conditioning` (end of appendix.tex) — will gain new numerical material.
- `app:errors` — the cutoff paragraph question.
- `app:blockpm:validation` — evolution-time caveat for the exact-diag table.
- `app:selector` — possibly.
- Also (main text, not yours): conclusions in `main.tex`, `numerics.tex` seed paragraph,
  `P3.tex` RTM/RDM promise, `results.tex` velocity claims, `app:highp` stale p=1.0 caption.

**Do not rewrite the content of these while that work is open** — style-only fixes there
should wait or be cleared with the author. Everything else in the appendix is fair game.
If the author says the investigation is finished, this freeze is lifted.

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
> `.agents/rules/thesis_style.md` and section 2.1 of `CLAUDE.md`. Then read
> `thesis/cft.tex`, `thesis/P3.tex` and `thesis/results.tex` so you have the voice of the
> thesis in front of you before writing anything.
>
> Your task is to rewrite the appendices in `thesis/appendix.tex` in that same voice. I have
> marked the passages I want changed with inline bracket comments `[...]` in the file — those
> are the work queue. Do not touch the sections listed as frozen in section 5 of
> `APPENDIX_REWRITE.md`.
>
> Work the way section 1 of that file describes: analyse all the bracket comments in the
> section we are on, then propose the changes to me inline in chat, quoting the current text
> and your proposed replacement, one block per change. Wait for my approval before editing any
> file. Once I approve, implement, then run the integrity checks and compile with tectonic and
> tell me what passed.
>
> Start with `\section{...}` [name the appendix section you want done first], and tell me if
> you cannot find bracket comments in it.
