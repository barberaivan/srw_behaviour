# srw_behaviour

Analysis code for a study of **southern right whale (*Eubalaena australis*)
behaviour as a function of kelp gull (*Larus dominicanus*) attacks**, at Península
Valdés, Argentina. The repository accompanies the manuscript in preparation.

It contains the R and Stan scripts for the two main analyses:

- `behaviour/` — behavioural-state models for mothers and calves (multinomial
  hidden-state models in Stan), plus behaviour imputation and prediction.
- `attack/` — kelp-gull attack **occurrence** and **intensity** models
  (mothers and calves). **Dormant** — the attack side is currently not planned
  for the paper; see [`attack/README.md`](attack/README.md).
- `plots/` — figures for the manuscript, combining mothers and calves results.

Small model outputs (parameter summaries, prediction tables) are committed under
each analysis' `files/` subfolder (or alongside the attack scripts). Figures
(`.png`) are committed too. The heavy and private artifacts are **not** in git —
see below.

---

## Repository structure

```
srw_behaviour/
├── behaviour/          the behavioural-state analysis (active)
│   ├── files/          small committed outputs: parameter summaries, prediction tables
│   └── figures/        figures produced by the analysis scripts themselves
├── plots/              manuscript figures combining mothers and calves
├── attack/             gull-attack occurrence and intensity models (dormant)
├── data/    → symlink  raw, private records (outside git; see setup below)
├── models/  → symlink  heavy generated objects, ~7 GB (outside git)
├── tools/              helper scripts (currently: inserting figures into the manuscript .docx)
├── CLAUDE.md           working agreement and guidance for Claude Code
├── ROADMAP.md          open threads: what is in flight and what still needs doing
└── setup.sh            links the data/model store into the repo (run once)
```

### What each script does

**`behaviour/` — the analysis in the paper**

| Script | What it does |
|---|---|
| `behaviour_mothers_analysis.R` | The mothers pipeline end to end: reads the raw records, recodes behaviour into the 8 categories, builds the Stan data, fits the model, imputes the unobserved behaviours, and computes the predictions and attack effects. |
| `behaviour_mothers_model.stan` | The mothers model: categorical likelihood on behavioural transitions, with the latent disturbance `z` and year splines. |
| `behaviour_calves_analysis.R` | The same pipeline for calves (9 years instead of 19). |
| `behaviour_calves_model.stan` | The calves model — as the mothers one, but with three separate increase parameters and smaller spline bases. |
| `behaviour_calves_assessing prior for beta.R` | Prior check that settled `beta ~ Normal(0, 1.5)` for calves, comparing calves' and mothers' betas in shared years. |
| `behaviour_period_predictions.R` | Averages the attack-scenario predictions over the manuscript's periods (1995, 2004-2010, 2011-2013, 2014-2019, 2020-2025) and over the years shared by mothers and calves, always within each posterior sample. Feeds Figure 4. |
| `behaviour_calves_obspred.R` | The calves' behaviour simulated at the estimated `z` (their "observed attack" prediction), and the four effect types computed from it. The calves' half of the section that the analysis script left un-run; feeds Figure 5. |
| `plots_script_mother calf.R` | Exploratory mother-vs-calf figures (attack scenarios, total effects, `z` dynamics). |

**`plots/` — manuscript figures**

| Script | What it does |
|---|---|
| `plots_script.R` | The earlier, exploratory version of the manuscript figures: attack scenarios, behaviour by year (disturbed / undisturbed / observed), the effects on mothers and calves, and the attack figures. Also defines `theme_mine()`. |
| `theme_paper.R` | Sourced by the figure scripts below: the behaviour coding, `separate_behav()`, the palettes and `theme_mine()`, defined once. |
| `figure_2_etogram.R` | **Figure 2**: stacked-bar etogram by year, mothers and calves, from the model-imputed behaviour. Currently produces `figure_2_etogram_v01.png`. |
| `figure_3_behaviour_timeseries.R` | **Figure 3**: behaviour by year under the disturbed, undisturbed and observed conditions, mothers and calves. Two layouts: `figure_3_behaviour_timeseries_faceted_v01.png` (separate blocks, patchwork) and `..._dodged_v01.png` (same panels). |
| `figure_4_attack_scenarios_periods.R` | **Figure 4**: behaviour under persistent attacks and after their cessation, by period. Two layouts: `figure_4_attack_scenarios_patchwork_v01.png` and `..._nested_v01.png`. |
| `figure_5_effects_summary.R` | **Figure 5**: the potential, short-term, long-term and total effects by year, mothers and calves keyed by colour. Produces `figure_5_effects_summary_v01.png`. |
| `figure_6_attacks_behaviour_mortality.R` | **Figure 6**: attacks, grouped behaviour classes and calf mortality through the years. Two layouts: `figure_6_attacks_behaviour_mortality_mothers_v01.png` and `..._both_v01.png`. Its mortality panel is an empty placeholder — see the note for Meri in the script. |

**`attack/` — dormant.** `mothers/`, `calves/` and `mothers2/` hold the
occurrence and intensity models (`*_occurrence.*`, `*_count.*`) for each target,
`mothers_spline years/` an earlier spline variant, `plots_attack.R` the figures
combining them, and `occurrence models trials.*` the exploratory fits. See
[`attack/README.md`](attack/README.md) and `attack/Hoja de ruta` (Spanish).

**`tools/`**

| Script | What it does |
|---|---|
| `insert_figures_into_docx.py` | Inserts figures and teal notes into the manuscript `.docx` under the house rules in `CLAUDE.md`: only new paragraphs, no tracked changes, and a check that every original paragraph survives byte-identical and in order. Refuses to run twice over the same file. |

**Root.** `Supplementary information 1.qmd` (and its rendered
`Supplementary-information-1.pdf`) is the full model specification cited by the
manuscript's Methods.

### Figure naming

Manuscript figures are versioned in the file name, `<name>_vNN.png` with a
zero-padded two-digit number (`_v01`, `_v09`, `_v11`, …). Co-authors revise
figures a lot, so each new version handed to the manuscript gets a new number and
earlier files are not overwritten — the file name alone tells you which version
is the one embedded in the `.docx`. Alternative code for different versions may
be kept in the plotting scripts.

---

## Getting started (first-time setup)

This repository holds **code only**. Two kinds of files live **outside git**, in an
Insync-synced **store** folder, and are symlinked back into the repo:

| Symlink   | Holds                                                                                  |
|-----------|----------------------------------------------------------------------------------------|
| `data/`   | the **raw** whale/gull records — *private*; the authors do not share them publicly       |
| `models/` | heavy generated objects: posterior sample arrays, prediction samples, z-samples, fitted probabilities, attack `p_sim`, and anything that embeds the raw observations (Stan data, data subsets, imputed-behaviour files). ~7 GB. |

> **Why?** Code belongs in git/GitHub (which versions and backs it up); the private
> data must never be published there, and multi-GB posterior arrays would bloat the
> history. Keeping them apart also avoids the sync conflicts that arise when a
> cloud-sync tool and git fight over the same files.

Setup is **three steps**:

### 1. Get the code

```bash
git clone https://github.com/barberaivan/srw_behaviour.git
cd srw_behaviour
```

### 2. Get the store

The raw data is private — **request it from the authors.** You will receive a
folder named **`srw_behaviour-store`** containing a `data/` subfolder (and,
optionally, a `models/` subfolder with the precomputed heavy objects; if it is
absent, `setup.sh` creates an empty one and the scripts regenerate it). Put the
store somewhere stable and remember the path.

### 3. Link the store into the repo

From the repo root, run `setup.sh` **once**, giving it the path to the store:

```bash
./setup.sh /full/path/to/srw_behaviour-store
```

This creates two symlinks (`data/` and `models/`) pointing into the store, and
remembers the path (in a local, gitignored `.local-paths` file), so any later
re-run is just `./setup.sh` with no argument.

To confirm it worked:

```bash
ls data/      # raw csvs: srw_behaviour_data*.csv, Monitoreo*.csv …
ls models/    # heavy .rds objects (empty on a fresh store — regenerated by the scripts)
```

Both symlinks and `.local-paths` are gitignored, so nothing machine-specific is
ever committed.

> *Symlinks require Linux or macOS (or Windows with WSL / Developer Mode enabled).*

---

## Running the analyses

Open `srw_behaviour.Rproj` in RStudio — this sets the working directory to the repo
root, which all paths assume. Then run the scripts under `behaviour/` and `attack/`.
Scripts read raw data through `data/…` and read/write heavy objects through
`models/…`, e.g.:

```r
bdata <- read.csv("data/srw_behaviour_data_1995-2025.csv")
bm    <- readRDS("models/behaviour/behaviour_mothers_model_samples.rds")
```

A plain clone (with only `data/` linked) can refit everything; if the store also
ships `models/`, you can reproduce figures without the long model runs. The Stan
models require **cmdstanr**/**rstan** and are computationally heavy.

> **Heads-up for contributors:** because data and heavy objects are outside git,
> uncommitted code is backed up nowhere — `git commit && git push` often, work on
> one machine at a time, and `git pull` before you start.

---

See [`ROADMAP.md`](ROADMAP.md) for open threads on the manuscript, and the
decisions already settled. The `attack/` analysis is currently **not in active
development** — see [`attack/README.md`](attack/README.md).
