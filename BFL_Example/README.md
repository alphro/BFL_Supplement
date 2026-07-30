# BFL reproducible example: PHMRC, target = AP

Runs the **BFL** package on the real **PHMRC** verbal-autopsy data, with **AP**
as the target site and the other five PHMRC sites (Mexico, Bohol, Dar, Pemba, UP)
as sources. Each source fits a local LCVA model and shares only its summary; BFL
aggregates them onto AP. Reproduces the paper's two LCVA baselines
(**local-self**, **local-avg**) and BFL variants (**base**, **domain**,
**partial**, **mix**) across all three within-target label-shift cases --- **no
shift** (20% of AP labeled), **mild**, and **severe** --- then shows the
package's S3 methods (`print`, `summary`, `plot`). `base` uses no target labels,
so it is shift-invariant and only run in the no-shift case.

## Files (self-contained)

- `BFL_AP_example.Rmd` --- the walkthrough. Knit it to reproduce the AP numbers.
- `helpers.R` --- everything else: the two PHMRC data-access helpers plus the
  LCVA glue (loads sites, fits LCVA per source, scores the target, assembles
  `local_summaries` for `run_BFL()`).
- `phmrc_clean.csv` --- the bundled PHMRC data.

## Dependencies

Only two packages, both not on CRAN:

```r
install.packages("remotes")
remotes::install_github("richardli/LCVA")   # LCVA fitting
# install BFL from its package directory:
# R CMD INSTALL /path/to/BFLpkg/BFL
```

Everything else (the PHMRC data and the data-access helpers) is bundled in this
folder, so there is nothing external to point at.

## Running

Run **from this folder** so the relative path to `phmrc_clean.csv` resolves:

```r
rmarkdown::render("BFL_AP_example.Rmd")
```

**Speed / memory note:** the LCVA fits (one per source, `Nitr = 2000`, 3 chains)
are the slow part and can take several minutes. Prediction with
`return_likelihood = TRUE` allocates an `Nitr x N x C` array, so the predict draws
are set to 2000 (drop 1000) to stay within a laptop's memory; the paper used 4000.
For a quick smoke test, lower `LCVA_TRAIN_ARGS$Nitr` / `LCVA_PRED_ARGS$Nitr` (e.g.
to 500) in `helpers.R`. If you still hit a "vector memory limit" error, lower the
predict draws further or raise the limit with `mem.maxVSize()`.

## What it does, step by step

1. Load AP (target) + the five sources; `prepare_ap(shift = ...)` builds the
   labeled/unlabeled split for the requested case (no-shift = 20% labeled random;
   mild = Dirichlet prevalence shift on a stacked frame; severe = per-cause
   Beta(0.2, 0.2) unlabeled fractions).
2. Fit each source's local LCVA model and score it on the AP rows -> `posterior_phi`.
3. Run the paper variants through the package:
   - **base** --- sources only, no target labels (no-shift only).
   - **domain** --- sources + a target self-model on the unlabeled rows; the
     labeled truth drives the CSMF correction (`Y_add`).
   - **partial** --- observed labels fed into the model directly.
   - **mix** --- labeled rows split per cause between the self-model and the model.
   The shift cases flip `label_shift = TRUE` and score the unlabeled subset via
   `eval_idx`; the input builders handle the severe self-model placement.
4. Two LCVA **baselines** (`local_self`, `local_avg`) are reported alongside, as
   in the paper (these do not use `run_BFL`).
5. `score_BFL()` reports Top-1, balanced, and CSMF accuracy for each; the three
   cases (no-shift, mild, severe) each print a comparison table.
6. `print` / `summary` / `plot` on the `BFL` fit and `BFL_score` objects.
