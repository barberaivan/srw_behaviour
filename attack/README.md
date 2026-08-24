# `attack/` — kelp-gull attack models

Occurrence and intensity (count) models for kelp gull attacks on southern right
whale mothers and calves at Península Valdés.

## Status: not in active development

**The current plan is not to include the attack side in the paper.** This folder
is kept for reference — the models were fitted and their results stand — but it
is dormant: do not extend, refactor or re-run it, and do not spend effort tidying
it, unless someone explicitly asks. Active work happens in `behaviour/` and
`plots/`.

The paper is iterative and co-authors may ask for attack results at any point, so
nothing here is deleted pre-emptively. The open questions below are recorded so
that whoever picks this up again — possibly months from now — does not have to
rediscover them.

## Layout

| Folder / file | What it is |
|---|---|
| `calves/` | calves attack occurrence + count models (`*_occurrence.*`, `*_count.*`) |
| `mothers/` | mothers attack models, occurrence **without** a year spline |
| `mothers_spline years/` | mothers attack models, occurrence **with** a single year spline (`model_gaviots_*_02`) |
| `mothers2/` | a deliberately simple mothers count model, tried as an alternative |
| `occurrence models trials.*` | early exploratory occurrence-model trials |
| `Hoja de ruta` | the running model roadmap and decision log (Spanish) |
| `plots_attack.R`, `attack_figure_all.png` | combined attack figures |

## Open decisions, if this is ever resumed

These were noticed while reorganising the repo (2026-06-27). None of them block
anything today; they only matter if the attack analysis returns to the paper.

### 1. Which mothers occurrence model is final — `mothers/` or `mothers_spline years/`?

Two variants of the mothers attack-occurrence model coexist: one **without** a
year spline (`mothers/`) and one **with** a single spline whose predictions are
averaged over years (`mothers_spline years/`).

The `Hoja de ruta` entry of **2022-10-19** leans towards keeping the **single
spline model**: it is cleaner to fit one model and present the prediction
averaged across years. The noted drawback is that this average is over the
*observed* years rather than unconditional on year — judged acceptable there,
since little varied between years.

**To resolve:** pick one, delete the other (code + its objects under
`models/attack/…`), and say which one won here and in `Hoja de ruta`.

### 2. `mothers2/` can probably go

A deliberately simple count model tried as an alternative; **it did not fit
better**. The `Hoja de ruta` entry of 2022-05-14 says plainly: *"En mothers2 puse
un modelo de count muy simple queriendo ver si ajustaba mejor. No funcionó. (Se
puede borrar luego)."*

**To resolve:** safe to delete — both `attack/mothers2/` and its store objects
under `models/attack/mothers2/`.

### 3. Small outputs sit next to the scripts, not in a `files/` subfolder

`behaviour/` keeps its small committed outputs (parameter summaries, prediction
tables) in `behaviour/files/`. The attack analyses instead leave theirs
(`*_summary.rds`, `*_prediction table.rds`) beside the scripts in each model
folder.

**To resolve:** if this analysis is revived, consider giving `attack/` a `files/`
subfolder to match `behaviour/`. This is cosmetic — it would mean updating the
paths in every `attack/**/*.R` script, so it is not worth doing while dormant.

### 4. A follow-ID ordering bug flagged in `Hoja de ruta`

The roadmap carries an uppercase note that the follow ID needs to be ordered
properly for the posterior check: correct in `calves_occurrence`, possibly wrong
in the others — with the caveat that the author suspected it may have no real
effect. Worth verifying before trusting any attack posterior check.
