# ROADMAP

Open threads for the manuscript. The paper is **iterative** — co-authors ask for
this or that as it goes — so this is not a checklist to burn down. Things get
done **as needed**; an item sitting here unresolved is not a problem, it is just
a note to whoever needs it next.

## Open

### Attack-scenario `_avg` prediction tables are missing

`plots/plots_script.R` (Figure 1) and `behaviour/plots_script_mother calf.R` read
three prediction tables:

```
behaviour/files/behaviour_mothers_predictions_attack scenarios table.rds       # exists
behaviour/files/behaviour_mothers_predictions_attack scenarios table_avg.rds   # MISSING
behaviour/files/behaviour_calves_predictions_attack scenarios table_avg.rds    # MISSING
```

The two `_avg` tables — predictions averaged over the calves years (2013,
2015–2018) — have never been generated, so **Figure 1 cannot be produced until
they are**. They are written by:

- `behaviour/behaviour_mothers_analysis.R` (~line 1647)
- `behaviour/behaviour_calves_analysis.R` (~line 1599)

in the section *"The same plot but averaging years"*. This is only
post-processing of the existing `…predictions_attack scenario samples.rds`
objects in the store — **not** a model refit — so it is cheap: load the
pred_samples object and run that section.

### Attack analysis

Not in active development (the plan is to leave it out of the paper). Its open
decisions — which mothers occurrence model is final, whether `mothers2/` can go,
output-folder layout, a flagged follow-ID bug — are documented in
[`attack/README.md`](attack/README.md). Nothing there needs doing unless a
co-author asks for attack results.

## Background — already settled

Kept so the reasoning is not rediscovered later.

- **Behaviour models always include year.** The no-year calves variant
  (`behaviour_calves_analysis_noYear.R` + `behaviour_calves_model_noYear.stan`)
  was deleted on 2026-08-24.
- **The mothers behaviour model uses `Bz = 5`** for the z-parameter year spline
  (4 non-linear basis functions + intercept, matching the tmat spline). The
  earlier low-basis run (`Bz = 3`), whose objects were suffixed `--lowB`, was
  deleted from the store on 2026-08-24. The final run carries the plain names.
- **Prediction-table names were reconciled** on 2026-08-24: the plot scripts read
  `…attack scenarios table.rds` / `…table_avg.rds`, which is what the current
  analysis scripts write. (They previously read older `…table object*.rds` names
  produced by the superseded no-year script.)
- **Data and heavy objects live outside git**, in
  `~/Insync/Whales/Behaviour/srw_behaviour-store/{data,models}`, symlinked in by
  `setup.sh`. During that move, heavy `.R` files that were really `saveRDS()`
  output were renamed to `.rds`, and the old name suffixes (`_R object`,
  `_Robject`, `_R list`, `_temporal`) were stripped.
- **Moved out of the repo** into `~/Insync/Whales/Behaviour/`: `ISEC/` (conference
  speech, slides and figures) and `to remove/` (old "FULL LIKELIHOOD" model
  variants and earlier imputations, superseded by the current models — can be
  deleted from there eventually).
