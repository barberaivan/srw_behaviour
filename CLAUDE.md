# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this project is

Analysis code for a study of **southern right whale (*Eubalaena australis*)
behaviour as a function of kelp gull (*Larus dominicanus*) attacks** at Península
Valdés, Argentina. The work supports a manuscript in preparation. Two analyses:

- `behaviour/` — behavioural-state (multinomial hidden-state) Stan models for
  mothers and calves, plus behaviour imputation and predictions.
- `attack/` — kelp-gull attack **occurrence** and **intensity** models
  (mothers and calves). See `attack/Hoja de ruta` for the model roadmap (Spanish).
  **Not in active development** — see below.
- `plots/` — figures combining mothers and calves results.

### The attack analysis is not in active development

The current plan is **not to include the attack side in the paper**, so `attack/`
is dormant: keep it in the repo for reference, but do not extend, refactor or
re-run it, and do not spend effort tidying it, unless explicitly asked. Active
work happens in `behaviour/` and `plots/`.

## Task workflow — read this first

Most tasks are written in the **manuscript itself**, a Word file kept in the
Insync folder (outside this repo):

```
/home/ivan/Insync/Whales/Behaviour/Paper/srw_behaviour_paper.docx
```

**This is the working manuscript.** The Google Docs draft it came from has been
abandoned — do not read or write it, and ignore the `.gddoc` pointer sitting next
to the file above. If a task refers to "the doc", "el doc", "the draft" or "el
paper", it means this `.docx`.

Instructions addressed to you are embedded inline in it, in the form:

```
[Claude, task-description]
```

**Always in square brackets, always naming Claude first** — the convention the
manuscript itself states, adding that plain square brackets are not used for
anything else, so a `[…]` is unambiguous. The instructions are usually Spanish
even though the surrounding text is English.

When asked to work from the draft, read the file, find the `[Claude, …]` markers,
and act on them.

### Reading and writing the .docx

A `.docx` is a zip. Read it by extracting the text from `word/document.xml`
(paragraphs are `<w:p>`, text runs are `<w:t>`); `word/media/` holds its images.
There is no need for the Google Drive connector any more.

**Do not rewrite the `.docx` in place unless explicitly asked.** Editing the XML
by hand risks damaging formatting, images, comments and tracked changes in a file
the co-authors are working on. By default, deliver the result of a `[Claude, …]`
task in the terminal (or as a repo file), quoting the marker it answers so it is
obvious where the text belongs, and let the user paste it in — using Word's
tracked changes or a colour if the edit needs to be visible. Never claim an edit
was made in the manuscript when the text was only printed here.

### Language: read Spanish, think and write English

**The paper is written in English.** The instructions in the draft, and parts of
its text, are often in Spanish (e.g. `[Claude, hacé tal cosa]`) — so are the
notes and roadmaps in this repo, and *la tesis* below.

Always **think in English**. When something is read in Spanish, translate it to
English first, and act on the English. Manuscript text produced for the paper is
always English, whatever the language of the instruction that asked for it. Reply
to the user in the language they wrote in, but never let Spanish source material
carry through into the manuscript.

### Reference for the writing — "la tesis"

When writing or revising manuscript text, the best reference for tone, structure
and background is **María Piotto's undergraduate thesis** (*tesina de grado*,
UNC, 2021), *"Efectos de los ataques de gaviotas cocineras (Larus dominicanus) en
el comportamiento, la mortalidad y éxito reproductivo…"*, directed by Carina F.
Marón and co-directed by Mariano Sironi:

```
/home/ivan/Insync/Whales/Behaviour/Paper/Manuscrito_MPiotto_ACTUALIZADO_tesis-grado.docx
```

In the manuscript the user refers to it as **"la tesis"** (or similar) —
that is this file. It covers the same study system, and **what we are doing now
is almost, but not exactly, the same**: follow it for framing, terminology and
references, but always check the current analysis before carrying over a
methodological claim, since the models here have moved on.

It lives in the Insync folder, outside the repo (~12 MB, Spanish, `.docx`). Read
it by extracting the text from `word/document.xml` inside the `.docx` (it is a
zip) rather than trying to open it directly. Reference it whenever it helps.

The same folder, `~/Insync/Whales/Behaviour/Paper/`, holds secondary material
worth consulting when useful: several older idea/results drafts (`behaviour
ideas y resultados 2022-*.doc`, `borrador ideas 2025-02-14.doc`, `effects and
time 2022-05-06.doc`) and `behaviour reclassification.xlsx`. These are older than
the current analysis — treat them as history, not as specification. The
`srw_behaviour_paper.gddoc` in that folder is the dead Google Docs pointer;
ignore it.

## Repository layout & the data/model store (important)

This repo is **code only**. Raw data and heavy generated objects live **outside
git** in an Insync-synced store, linked in as two **gitignored symlinks**:

- `data/`   → raw, **private** whale/gull records (the authors do not share these
  publicly). Never commit anything under `data/` or anything derived that still
  contains raw observations.
- `models/` → heavy generated objects (~7 GB): posterior sample arrays, prediction
  samples, z-samples, fitted probabilities, attack `p_sim`, and any
  raw-observation-bearing object (Stan data lists, data subsets, imputed-behaviour
  files).

First-time setup: `./setup.sh /path/to/srw_behaviour-store` (creates the symlinks;
saves the path to gitignored `.local-paths`). See `README.md` for details.

### Where outputs go
- **Heavy or raw-data-bearing** objects → `models/<analysis>/…` (in the store).
- **Small model outputs** (parameter summaries, prediction tables) → committed in
  each analysis' `files/` subfolder (`behaviour/files/…`) or next to the attack
  scripts. Figures (`.png`) are committed.
- Decide by **privacy first, then size (>~5 MB)**. When in doubt, treat it as
  private and put it in `models/`.

## Conventions

- **Working directory is the repo root.** Open `srw_behaviour.Rproj` (or
  `setwd()` to the repo root). All paths are relative to root: read raw data via
  `read.csv("data/…")`, heavy objects via `readRDS("models/…")`. Do **not**
  hardcode absolute `/home/...` or old `Insync/Whales/...` paths.
- **Serialized R objects use `saveRDS()`/`readRDS()` with a `.rds` extension.**
  (Historically some were saved with a `.R` extension — those have been renamed;
  do not reintroduce `.R` for data objects. `.R` is for scripts only.)
- Stan models (`.stan`) need **cmdstanr**/**rstan**; fitting is computationally
  heavy — don't re-run models unless asked. Prefer loading existing `models/…`
  objects.
- Match the surrounding code's style (2-space indent, existing naming). Notes and
  roadmaps in the repo are in Spanish; code/comments may be mixed.

## Housekeeping

- Open threads for the manuscript are tracked in `ROADMAP.md` at the repo root.
  The paper is iterative — co-authors ask for things as it goes — so items there
  get done **as needed**, not as a checklist to burn down. It also records
  decisions already settled, so the reasoning is not rediscovered later.
- Attack-specific open decisions live in `attack/README.md` (that analysis is
  dormant — see above).
- Commit/push only when asked. The remote is SSH
  (`git@github.com:barberaivan/srw_behaviour.git`); pushing does not prompt.
