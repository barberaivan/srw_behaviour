# Cleanup — pending tidy-up

Deferred items noticed while reorganising the repo (2026-06-27). None block the
analysis; resolve when convenient.

## Obsolete / superseded code to decide on

- **`attack/mothers2/`** — a deliberately simple count model tried as an
  alternative; it did not fit better. The `attack/Hoja de ruta` notes say
  *"Se puede borrar luego."* Safe to delete (code + its store objects under
  `models/attack/mothers2/`).
- **`attack/mothers/` vs `attack/mothers_spline years/`** — two mothers-attack
  model variants (no-spline vs single spline averaged over years). The notes
  (2022-10-19) lean toward keeping the **single spline model** and dropping the
  no-spline one. Pick the final model and remove the other.
- **`behaviour/behaviour_calves_analysis_noYear.R`** +
  **`behaviour/behaviour_calves_model_noYear.stan`** — the older calves version
  without year effects, superseded by the year-aware model. Remove once the
  final calves model is confirmed.
- **`--lowB` mother model files** (`behaviour_mothers_model_samples--lowB.rds`,
  `behaviour_mothers_model_data--lowB.rds`, in the store) — a low-basis spline
  run kept for comparison; drop once the full-basis run is final.

## Already moved out of the repo (for reference)

Into `~/Insync/Whales/Behaviour/`:
- `ISEC/` — conference (ISEC) speech, slides and figures.
- `to remove/` — old "FULL LIKELIHOOD" model variants and earlier imputations,
  superseded by the current models. Review and delete from there eventually.

## Misc

- The raw data + heavy objects now live in
  `~/Insync/Whales/Behaviour/srw_behaviour-store/{data,models}` (symlinked in via
  `setup.sh`). Heavy `.R` files that were actually `saveRDS()` output were renamed
  to `.rds` and the novice name suffixes (`_R object`, `_Robject`, `_R list`,
  `_temporal`) were stripped.
- `plots/plots_script.R` reads a few `…attack scenarios table object*.rds` names
  that differ from what `behaviour_*_analysis.R` now writes
  (`…attack scenarios table.rds`); reconcile these names when revisiting the
  combined plots.
- Consider giving `attack/` analyses a `files/` subfolder for their small
  committed outputs, to match `behaviour/files/` (currently they sit next to the
  scripts).
