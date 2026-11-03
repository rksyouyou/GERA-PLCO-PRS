# Archived DCA code: manuscript version with DCA

This directory preserves the DCA implementation used in the manuscript alternative `v7_with_DCA`. DCA was removed from the current V7 manuscript on 2026-11-03. The code is unchanged and is not part of the current manuscript workflow. The archive contains source code only.

Run from the repository root (or V36 `R/current/` for the local mirror):

```r
source("fit_models.R")
source("optimism.R")
source("archive/with_DCA/decision_curve.R")
```

Sourcing defines functions only. Supply authorized local inputs and call the functions explicitly to run an analysis.

`decision_curve()` calculates net benefit and standardized net benefit for Base, PRS, Full, All and None at observed follow-up, across 1%–20% thresholds in steps of 0.1 percentage points. Selection uses `p >= threshold`. The archived manuscript retains a vertical 10% reference line. `decision_curve_optimism(data, models, indices)` subtracts mean training-minus-original raw net benefit across the same 500 participant resamples as AUC correction, then divides by original GERA AN prevalence. All/None reference curves and PLCO external-validation curves remain unchanged. Corrected GERA curves are point estimates without DCA CIs or P values.

- Zhang et al. *Ann Transl Med*. 2018;6:308. Decision curve analysis: a technical note. https://doi.org/10.21037/atm.2018.07.02
