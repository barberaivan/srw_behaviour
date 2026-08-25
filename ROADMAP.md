# ROADMAP

Open threads for the manuscript. The paper is **iterative** — co-authors ask for
this or that as it goes — so this is not a checklist to burn down. Things get
done **as needed**; an item sitting here unresolved is not a problem, it is just
a note to whoever needs it next.

## Open

### Pick one version of Figures 3 and 6 — Iván and Meri

Both layouts of each are in the `.docx` so co-authors can choose, with Iván's
vote written next to them: version B for Figure 3, version A for Figure 6. Once
picked, the loser's code can be dropped from the script.

| Figure | Script | Versions |
|---|---|---|
| 3 | `plots/figure_3_behaviour_timeseries.R` | `..._faceted_v01.png` (A, separate blocks), `..._dodged_v01.png` (B, same panels) |
| 6 | `plots/figure_6_attacks_behaviour_mortality.R` | `..._mothers_v02.png` (A, mothers only), `..._both_v02.png` (B, mothers and calves) |

Figure 4 is already resolved: the manuscript carries a single image, the
patchwork version with panels A and B.

The `.docx` still embeds the Figure 6 `_v01` images while the repo has `_v02`.
They differ only in the wording of the empty panel C placeholder, so this
resolves itself the moment that panel is redrawn with the real series.

### Figure 6, panel C: which calf mortality series? — Meri

The panel is drawn empty on purpose, with the question written into the `.docx`
in teal next to the figure: should it be the dead calves counted by the Southern
Right Whale Health Monitoring Program, as in Piotto et al. (2024)? And do we have
mortality data up to 2025, or does the series stop earlier? The panel draws
itself as soon as the numbers are put in the `mortality` data frame in
`plots/figure_6_attacks_behaviour_mortality.R`.

### Attack analysis

Not in active development (the plan is to leave it out of the paper). Its open
decisions — which mothers occurrence model is final, whether `mothers2/` can go,
output-folder layout, a flagged follow-ID bug — are documented in
[`attack/README.md`](attack/README.md). Nothing there needs doing unless a
co-author asks for attack results.

## Background — already settled

Kept so the reasoning is not rediscovered later.

### The manuscript file

- **The draft is renamed as it is passed around**,
  `srw_behaviour_paper_<initials>_<date>.docx`. As of 2026-08-25 it is
  `srw_behaviour_paper_IB_2026-08-25.docx`, sent to Meri, and the name changes
  again on the way back. Resolve it by listing the folder and taking the most
  recent non-`BACKUP` file; both scripts in `tools/` do this, and `CLAUDE.md`
  says so.
- **`tools/insert_figures_into_docx.py` and
  `tools/insert_results_text_into_docx.py`** hold the XML machinery for writing
  into the `.docx` under the house rules (teal, new paragraphs only, every
  original paragraph verified byte-identical). Reuse them for the next batch
  rather than rediscovering the image XML, the EMU sizing and the verification.
  The second one also carries the Results text that was handed to the
  manuscript, so the wording is versioned with the code that produced its
  numbers.

### Figures and predictions

- **Where they come from.** All the `plots/figure_N_*.R` scripts run from the
  repo root in a few seconds; only `behaviour/behaviour_calves_obspred.R`
  (~10 min, ~4 GB) and `behaviour/behaviour_period_predictions.R` (~25 s, ~3 GB)
  are heavy. They post-process the fitted objects in `models/`; no model is
  refitted.
- **Careful with `plots/plots_script.R`.** Its section headers use their own
  numbering, which does not match the manuscript's: its "Figure 1" is the
  attack-scenario plot, its "Figure 2" the mothers' behaviour time series and its
  "Figure 3" the effects summary. Figures 3 to 6 of the paper are the
  `plots/figure_N_*.R` scripts; `plots_script.R` is the earlier, exploratory
  version of the same material.
- **Figure 4 periods keep only the years actually sampled**, so 2014-2019 is
  2015-2019, 2020-2025 is 2021, 2023 and 2025, and the calves' 2011-2013 is 2013
  alone. Multi-year periods are averaged within each posterior sample before the
  posterior is summarised. Curves are cut at 70 min although the tables run to
  150.
- **Figure 5 leaves out the calves' long-term and total effects.** They *were*
  computed (`behaviour/files/behaviour_calves_predictions_total effects.rds`
  holds all four) but use 2013 as the baseline year instead of 1995, so they are
  not the same quantity as the mothers'.
- **The calves' "observed" short-term effect** needs the behaviour predicted at
  the estimated `z`, which the analysis script never ran.
  `behaviour/behaviour_calves_obspred.R` runs it, advancing all 8,000 posterior
  samples together instead of one at a time. The draws are not identical to what
  the original triple loop would have produced, but they come from the same
  distribution, and the resulting distribution matches the observed proportions
  year by year, as the mothers' does.
- **Figure 6, panel A counts attacks on the mothers' dataset**, which covers the
  whole 1995-2025 series (attacks on calves were recorded from the start; it is
  the calves' *behaviour* that begins in 2013). Occurrence and intensity share a
  panel through a secondary axis, and lines break at the years without sampling.
- **Figure 6, panel B shows rest and travel split by exposure**, four curves.
  In-place activity is not shown, so they do not add up to one; the caption says
  so.

### Models and the data store

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
- **The attack-scenario `_avg` prediction tables now exist.** They had never
  been generated, which blocked Figure 1 of `plots/plots_script.R`. They are
  written by `behaviour/behaviour_period_predictions.R`, together with the
  period tables Figure 4 needs.
- **Data and heavy objects live outside git**, in
  `~/Insync/Whales/Behaviour/srw_behaviour-store/{data,models}`, symlinked in by
  `setup.sh`. During that move, heavy `.R` files that were really `saveRDS()`
  output were renamed to `.rds`, and the old name suffixes (`_R object`,
  `_Robject`, `_R list`, `_temporal`) were stripped.
- **Moved out of the repo** into `~/Insync/Whales/Behaviour/`: `ISEC/` (conference
  speech, slides and figures) and `to remove/` (old "FULL LIKELIHOOD" model
  variants and earlier imputations, superseded by the current models — can be
  deleted from there eventually).
