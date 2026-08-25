# ROADMAP

Open threads for the manuscript. The paper is **iterative** — co-authors ask for
this or that as it goes — so this is not a checklist to burn down. Things get
done **as needed**; an item sitting here unresolved is not a problem, it is just
a note to whoever needs it next.

## Done, waiting on Iván

### Results text + three figure captions — Claude, 2026-08-25

Acting on the `[Claude, ...]` markers still open in the manuscript's Results:

| Marker (paragraph) | What it asks |
|---|---|
| "write the results sections" | 1-2 paragraphs per figure, placed **before** each figure |
| Figure 4 caption | "complete this image caption" |
| Figure 5 caption | "explain again here the meaning of each effect, very succinct, small" |
| Figure 6 caption | "make this caption a bit more complete" |

The Figure 2 and Figure 3 captions carry explicit *do not edit* markers, so they
are left alone.

**Where the numbers come from.** Everything quoted in the text was recomputed
from the committed prediction tables, not read off the figures:
`behaviour/files/behaviour_{mothers,calves}_predictions_total effects.rds`
(the four δ effects per year), `…predictions disturbed and undisturbed by
year.rds` (the eight category probabilities per year and scenario),
`…attack scenarios table_periods.rds` (the time courses of Figure 4, which run to
150 min even though the figure is cut at 70), and the imputed-behaviour data in
`models/behaviour/behaviour_{mothers,calves}_data with imputed behaviours.rds`
(attack occurrence and intensity per year, and the grouped proportions of
Figure 6B).

**Decisions taken:**

- **The text is inserted into the `.docx` in teal**, following the house rules,
  with `tools/insert_results_text_into_docx.py` (same machinery as the figure
  insertion). A pre-edit backup is kept next to the manuscript.
- **Captions are inserted as new paragraphs, not edits.** The house rule forbids
  touching a paragraph a co-author wrote, and the captions are theirs, so each
  completed caption goes in as a new teal paragraph right below the original,
  ready for Iván to swap in. The `[Claude, ...]` markers stay where they are.
- **`plots/figure_6_attacks_behaviour_mortality.R` had an em-dash** in the panel C
  placeholder annotation, which the new no-em-dash rule forbids in figures. It
  was changed to a colon and the figures rebuilt as `_v02`. The em-dashes left in
  the plotting scripts are all in comments, which the rule does not cover.

**What went into the file.** 38 paragraphs, all teal, all new: five Results
blocks (one per figure, each *before* its figure, with an underlined subsection
heading in the Methods style) and the three captions. The install was verified
programmatically: all 208 original paragraphs survive byte-identical and in
order, the XML parses, the image count is unchanged (7), and the inserted text
carries no em-dash. The pre-edit file is
`~/Insync/Whales/Behaviour/Paper/srw_behaviour_paper_BACKUP_2026-08-25_pre-results.docx`.

**Doubts for Iván, none of them blocking:**

- **Subsection headings in the Results were my call.** The marker said "the
  results sections", so each figure's block got an underlined heading. Delete
  them if the Results should run as continuous prose.
- **The Figure 4 marker had already been half-resolved.** The manuscript now
  carries a single Figure 4 image and no version note, so Iván had already
  picked one; the caption I wrote describes what is in the file (the patchwork
  version, A and B panels).
- **The 1995 recovery result is the shakiest claim in the text.** Recovery from
  attacks is three to four times slower from 2011 onwards than it was in 1995,
  but 1995 is one year with the widest intervals in Figure 4. The text says so.
  This is the single most interesting result and also the one most likely to be
  challenged.
- **The dissociation claim is the backbone of the Results.** The long-term and
  total effects rise with year (Spearman 0.67 and 0.62) and are uncorrelated
  with the attack pressure of the year (0.01 and 0.12), while attacks on mothers
  actually *decline* over the series (-0.62). The text puts this to the
  discussion rather than explaining it. If the intended reading is different,
  the last paragraph of the Results is the one to rewrite.

**The manuscript was renamed** by Iván after this work went in, to
`srw_behaviour_paper_IB_2026-08-25.docx`, and sent to Meri. Initials plus date is
the convention when the draft is passed to a co-author, so the name changes
again on the way back: resolve it by listing the folder, never by hardcoding it
(both scripts in `tools/` now do this, and `CLAUDE.md` says so).

**Still open:** Iván and Meri still have to pick one version of Figures 3 and 6;
the Figure 6 panel C mortality series is still Meri's call. Note that the `.docx`
still embeds the Figure 6 `_v01` images while the repo now has `_v02` (they
differ only in the panel C placeholder wording), which resolves itself when the
panel is redrawn with the real series.

## Done, waiting on Iván and Meri

### Results figures (Figures 3-6) — Claude, 2026-08-24

Every `[Claude, ...]` marker in the manuscript's **Results** section has been
acted on. The figures are in `plots/`, their code is in the repo, and both the
figures and a short teal note for each were inserted into the `.docx`
(the working manuscript, at the time named `srw_behaviour_paper.docx`; the pre-edit file is
kept next to it as `srw_behaviour_paper_BACKUP_2026-08-24_pre-figures.docx`).
Where a marker asked for two layouts, both are in the file so co-authors can
choose — **Iván/Meri pick one per figure, and the loser's code can then be
dropped from the script.**

| Figure | Script | Files |
|---|---|---|
| 3 | `plots/figure_3_behaviour_timeseries.R` | `figure_3_behaviour_timeseries_faceted_v01.png`, `..._dodged_v01.png` |
| 4 | `plots/figure_4_attack_scenarios_periods.R` | `figure_4_attack_scenarios_patchwork_v01.png`, `..._nested_v01.png` |
| 5 | `plots/figure_5_effects_summary.R` | `figure_5_effects_summary_v01.png` |
| 6 | `plots/figure_6_attacks_behaviour_mortality.R` | `figure_6_attacks_behaviour_mortality_mothers_v01.png`, `..._both_v01.png` |

New non-figure code these needed: `behaviour/behaviour_period_predictions.R`
(period averages of the attack-scenario predictions) and
`behaviour/behaviour_calves_obspred.R` (the calves' behaviour simulated at the
estimated `z`, which the analysis script never ran). `plots/theme_paper.R` holds
the theme and behaviour coding the four scripts share, and
`tools/insert_figures_into_docx.py` is the script that put the figures into the
`.docx` — reuse it for the next batch rather than rewriting the XML machinery.

**Committed and pushed** on 2026-08-25 as `9f03950`, "Add the Results figures
(Figures 3-6) and the code behind them". The store objects
(`models/behaviour/behaviour_calves_model_samples_predict_observed_attack*.rds`)
live outside git, so they exist only on this machine.

**Where the figures come from, if a number has to be checked.** All four scripts
run from the repo root and take a few seconds each, except
`behaviour/behaviour_calves_obspred.R` (~10 min, ~4 GB) and
`behaviour/behaviour_period_predictions.R` (~25 s, ~3 GB). They only
post-process the fitted objects in `models/`; no model is refitted.

**Careful with `plots/plots_script.R`.** Its section headers use their own
numbering, which no longer matches the manuscript's: its "Figure 1" is the
attack-scenario plot, its "Figure 2" the mothers' behaviour time series and its
"Figure 3" the effects summary. Figures 3 to 6 of the paper are the
`plots/figure_N_*.R` scripts listed above; `plots_script.R` is the earlier,
exploratory version of the same material.

**Decisions taken along the way — Iván should check them:**

- **Figure 3.** Version A shares the y scale row by row between the mothers and
  the calves blocks, so the two halves can be read against each other; the calves
  block carries no y axis text. Version B keys the condition by colour and
  mother/calf by point shape *and* line type.
- **Figure 4.** Periods keep only the years actually sampled, so 2014-2019 is
  2015-2019 and 2020-2025 is 2021, 2023 and 2025. Calves get three periods
  (2011-2013 = 2013 alone). Curves are cut at 70 min, as in the other figures.
- **Figure 5.** The calves' long-term and total effects *were* computed
  (`behaviour/files/behaviour_calves_predictions_total effects.rds` holds all
  four), but they use 2013 as the baseline year instead of 1995, so they are not
  the same quantity as the mothers' and the figure leaves them out — as the
  marker asked.
- **Figure 5, calves' "observed" short-term effect.** It needs the behaviour
  predicted at the estimated `z`, and the calves' version of that simulation had
  never been run. `behaviour/behaviour_calves_obspred.R` runs it, advancing all
  8,000 posterior samples together instead of one at a time (11,829 iterations
  instead of 11,829 x 8,000). The draws are therefore not identical to what the
  original triple loop would have produced, but they come from the same
  distribution, and the resulting distribution matches the observed proportions
  year by year, exactly as the mothers' does.
- **Figure 6, panel A.** Attacks are counted on the *mothers'* dataset, which
  covers the whole 1995-2025 series (attacks on calves were recorded from the
  start; it is the calves' *behaviour* that begins in 2013). Occurrence and
  intensity share a panel through a secondary axis, and the lines break at the
  years without sampling instead of jumping the 1996-2003 gap.
- **Figure 6, panel B.** Following the marker, it shows rest and travel (slow +
  fast) split by exposure — four curves. **In-place activity is therefore not
  shown, and the four curves do not add up to one.** Say so in the caption, or
  ask for a fifth and sixth curve.

**Open for Meri — Figure 6, panel C (calf mortality).** The panel is drawn empty
on purpose, with the question written into the `.docx` in teal next to the
figure: which mortality series should go there — the dead calves counted by the
Southern Right Whale Health Monitoring Program, as in Piotto et al. (2024)? And
do we have mortality data up to 2025, or does the series stop earlier? The panel
draws itself as soon as the numbers are put in the `mortality` data frame in
`plots/figure_6_attacks_behaviour_mortality.R`.

**Not done, and not asked for.** The `{Primer resultado muy general: etograma.
Describir el etograma}` note that opens the Results is in curly brackets, so it
is a pending task rather than an instruction to Claude; no text was written for
it. Say the word and the paragraph gets written.

## Open

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
