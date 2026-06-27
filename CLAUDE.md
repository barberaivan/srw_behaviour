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
- `plots/` — figures combining mothers and calves results.

## Task workflow — read this first

Most tasks are written in the **Google Docs paper draft** we use as the working
manuscript:

<https://docs.google.com/document/d/1FaTxfsi65lnKZ0DcWaMkfEj42fSK54IWUQDJYMcHneQ/edit>

Instructions addressed to you are embedded inline in that document in the form:

```
[Claude, task-description]
```

When asked to work from the draft, read the doc (via the Google Drive
integration), find the `[Claude, …]` markers, and act on them.

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

- Pending tidy-up tasks (obsolete model variants, name mismatches) are tracked in
  `CLEANUP.md`.
- Commit/push only when asked. The remote is SSH
  (`git@github.com:barberaivan/srw_behaviour.git`); pushing does not prompt.
