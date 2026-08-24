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

Most tasks are written in the **Google Docs paper draft** we use as the working
manuscript:

<https://docs.google.com/document/d/1FaTxfsi65lnKZ0DcWaMkfEj42fSK54IWUQDJYMcHneQ/edit>

Instructions addressed to you are embedded inline in that document in the form:

```
[Claude, task-description]
```

**Always in square brackets, always naming Claude first** — this is the
convention the document itself states.

When asked to work from the draft, read the doc (via the Google Drive
integration), find the `[Claude, …]` markers, and act on them.

### Language: read Spanish, think and write English

**The paper is written in English.** The instructions in the draft, and parts of
its text, are often in Spanish (e.g. `[Claude, hacé tal cosa]`) — so are the
notes and roadmaps in this repo, and *la tesis* below.

Always **think in English**. When something is read in Spanish, translate it to
English first, and act on the English. Manuscript text produced for the paper is
always English, whatever the language of the instruction that asked for it. Reply
to the user in the language they wrote in, but never let Spanish source material
carry through into the manuscript.

### The doc is read-only for Claude

The Google Drive connector can **read** the document (`read_file_content`, which
also returns comment threads) but **cannot write to it**. Its `update_file` tool
only changes metadata — title and parent folder — and there is no access to the
Google Docs API, so Claude cannot insert text, let alone apply colours,
highlights or Suggesting-mode edits inside the doc.

Therefore: never claim an edit was made in the doc. Deliver the result of a
`[Claude, …]` task in the terminal (or as a repo file), formatted so it can be
pasted into the draft, and quote the marker being answered so it is obvious which
task the text belongs to. The person then pastes it in — and, if edits need to be
visible, does so with Google Docs' own **Suggesting** mode or a colour, which
tracks authorship far better than anything Claude could mark by hand.

### Reference for the writing — "la tesis"

When writing or revising manuscript text, the best reference for tone, structure
and background is **María Piotto's undergraduate thesis** (*tesina de grado*,
UNC, 2021), *"Efectos de los ataques de gaviotas cocineras (Larus dominicanus) en
el comportamiento, la mortalidad y éxito reproductivo…"*, directed by Carina F.
Marón and co-directed by Mariano Sironi:

```
/home/ivan/Insync/Whales/Behaviour/Paper/Manuscrito_MPiotto_ACTUALIZADO_tesis-grado.docx
```

In the Google Docs draft the user refers to it as **"la tesis"** (or similar) —
that is this file. It covers the same study system, and **what we are doing now
is almost, but not exactly, the same**: follow it for framing, terminology and
references, but always check the current analysis before carrying over a
methodological claim, since the models here have moved on.

It lives in the Insync folder, outside the repo (~12 MB, Spanish, `.docx`). Read
it by extracting the text from `word/document.xml` inside the `.docx` (it is a
zip) rather than trying to open it directly. Reference it whenever it helps.

The same folder, `~/Insync/Whales/Behaviour/Paper/`, holds secondary material
worth consulting when useful: `srw_behaviour_paper.gddoc` (the local pointer to
the Google Docs draft above), several older idea/results drafts (`behaviour ideas
y resultados 2022-*.doc`, `borrador ideas 2025-02-14.doc`, `effects and time
2022-05-06.doc`) and `behaviour reclassification.xlsx`. These are older than the
current analysis — treat them as history, not as specification.

**Access requirement:** reading the doc needs the Google Drive connector to be
authorized as a Google account that the doc is shared with. The doc is owned by a
different Google account, so the user must (a) share it with the Google account
the connector signs in as — `ivanbarbera93@gmail.com` — and (b) authorize the
connector as that same account (`/mcp`). There is no "Claude" Google identity; the
file is read as whichever Google account approves the OAuth flow. If a Drive call
returns a token/authorization error, ask the user to re-authorize with `/mcp`
(restarting the session may be needed for a refreshed token to take effect).

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
