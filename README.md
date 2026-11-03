# GERA–PLCO PRS

R source code and fixed GERA model parameters for the current three-model GERA–PLCO manuscript (2026-11-03; DCA removed). GERA is the development population and PLCO is independent external validation. The repository contains no participant records, genotypes, PRS observations, imputations, fitted-model RDS files, bootstrap draws, or study output files. Supply authorized study data locally and keep inputs, checkpoints and outputs outside this repository.

## Model definitions

All models use a complementary log-log link and a nondecreasing SCOP-spline, `s(time_years, bs = "mpi", k = 7)`.

| Code key | Manuscript label | Predictors in addition to the time smooth | Historical V36 key |
| --- | --- | --- | --- |
| `Base` | Model 1 | Baseline colonoscopy findings | `Base` |
| `PRS` | Model 2 | Baseline findings, PRS, and baseline findings × PRS | `Model3` |
| `Full` | Full model | The PRS model plus age, sex, BMI, and colorectal cancer family history | `Full` |

The previous age/sex and established-risk-factor-only models are omitted from this release. Code uses descriptive keys to avoid collisions with historical model numbering. The primary manuscript comparison is Full versus Base. That contrast adds established risk factors and PRS jointly. Full versus PRS assesses the addition of established risk factors to the PRS model. It does not isolate PRS beyond established risk factors. All three pairwise AUC contrasts are retained. Table 3 presents a lower-triangular matrix with P-versus-Model-1 (Base) and P-versus-Model-2 (PRS) columns: PRS versus Base in the Model 2 row, and Full versus Base and PRS in the Full row. Dashes denote self-comparisons or comparisons shown elsewhere. The AUC functions already return all three contrasts; this presentation change does not alter their calculation.

## Files

| File | Purpose |
| --- | --- |
| `model_parameters.R` | Fixed coefficients, factor levels, spline knots and centering constants for the three models |
| `predict.R` | Prediction using fixed GERA parameters |
| `fit_models.R` | Clinical imputation, cohort-specific PRS standardization and optional development fitting |
| `performance.R` | Fixed-model AUC, paired participant bootstrap CIs and Wald comparisons |
| `optimism.R` | The original 500-resample GERA AUC point correction, with refitting |
| `two_stage_bootstrap.R` | Independent, resumable outer refits and corrected-AUC CI/P summaries |
| `clinical_utility.R` | Categorical NRI at 10% and NNS for the highest-risk 10% and 20% |
| `calibration.R` | Grouped 5-year calibration, with optional GERA point correction |
| `point_score.R` | Integer scores, hypothetical profile grid and 3/5/7/10-year risk summaries |

Sourcing files defines functions without fitting models or running bootstrap analyses. Optional analysis functions execute only when explicitly called.

## Fixed-model prediction

Run from the repository directory:

```r
source("predict.R")
parameters <- load_parameters()
# newdata is an authorized, locally held data.frame.
risk <- predict_an(newdata, model = "Full", parameters = parameters)
p <- predict_models(newdata, parameters)
```

| Column | Required values |
| --- | --- |
| `age_at_base` | Age in years, centered at 50 internally |
| `sex` | `Female` or `Male` |
| `bmi` | `<25`, `25-29.9`, or `30+` |
| `colo_fh` | Colorectal cancer family history: `No` or `Yes` |
| `basecolo_out` | `Normal`, `Aden1-2`, `Aden3+`, or `AA` |
| `prs_z` | Technically processed PRS, standardized in the relevant final cohort |
| `time_years` | Positive follow-up interval in years |

Aden1-2 and Aden3+ denote 1–2 and at least 3 nonadvanced conventional adenomas. AA denotes advanced adenoma. Only predictors in the selected model are required. GERA received no additional A/O kit correction. PLCO platform correction precedes standardization. GERA and PLCO PRS values were standardized separately in their final populations, without within-platform standardization. The prediction function does not compute a PRS from genotypes or standardize a single patient. GERA coefficients stay fixed for PLCO.

`prepare_cohort()` expects `Y` (0/1 AN at the observed follow-up examination), `bs_time` (days), `prs` (technically processed PRS) and the baseline/established-risk-factor columns above. It imputes age, sex, BMI and family history with `missForest`, then standardizes PRS. In bootstrap analyses these imputations and PRS preprocessing remain fixed. Cohort eligibility and upstream genetic processing require restricted inputs and are outside this release. The final populations contain 5,337 participants/397 AN in GERA, 2,513/209 in PLCO main and 1,876/170 in PLCO sensitivity. GERA excludes baseline CSSPs. PLCO retains missing FSG intervals under the 180-day eligibility rule, and the sensitivity analysis removes Hyperplastic Polyp and Benign Polyp, NOS only from Normal.

## AUC point estimates and inference

AUC discriminates AN detected at the observed follow-up colonoscopy. Calibration evaluates 5-year predictions separately.

GERA corrected point estimates subtract mean optimism from the original apparent AUC. Mean optimism is the average training-sample minus original-cohort AUC across 500 paired participant bootstrap resamples. All three models use the same sampled rows, with coefficients, the time smooth and smoothing parameters refitted. Clinical imputations and PRS preprocessing remain fixed.

```r
source("fit_models.R")
source("optimism.R")
# Deliberate refitting only. gera must already be prepared.
models <- fit_models(gera)
point_result <- optimism_bootstrap(gera, models, B = 500L, seed = 123L)
point <- setNames(point_result$AUC$corrected_AUC, point_result$AUC$model)
```

GERA 95% CIs use type-7 percentile intervals from a two-stage participant bootstrap: 5,000 completed outer samples, each with 200 inner samples. Within an outer sample, refit the three models and calculate their apparent AUCs. For each paired inner sample, refit the models and evaluate AUC on the inner training sample and its parent outer sample. Subtract mean inner optimism from each outer apparent AUC. Approximate two-sided Wald tests divide the **original 500-resample corrected point difference** by the **SD of paired corrected differences across outer samples**. This SD is not divided by the square root of the outer count. All tests are unadjusted for multiple comparisons. The two-stage CI follows Noma et al.; the paired Wald comparison is an approximate extension.

`two_stage_outer()` runs one outer sample, with checkpoints after each outer model and every five completed inner samples. Use unique global IDs/streams across repeated submissions. The study distributed 2,000 outer IDs per batch across 800 independent single-core jobs, using IDs `task, task + 800, ...`. Each batch used disjoint L'Ecuyer-CMRG streams. `outer_streams()` prepares streams in memory. A caller can pass a supplied frozen stream for exact resumption. Failed fits raise errors; retain their IDs, diagnostics and checkpoints outside the repository. Do not launch a new ID as a replacement for an existing stream while resuming it.

`two_stage_summary(point, corrected_draws, global_ids)` uses up to 5,000 successful rows in ascending global-ID order and reports `INTERIM` when fewer are available. On the first `FINAL` summary, save its `selected_ids` and reuse them through `selected_ids = ...` thereafter. Keep the complete attempted-ID/failure audit. The adopted study summary selected 5,000 successfully completed outer samples; inference therefore conditions on fitting success. Neither selection nor replacement is based on an AUC or P value. Full server draws and the frozen selection manifest are not distributed; this release does not claim independent local verification of those files or exact reproduction of screenshot-reported endpoints.

PLCO uses 10,000 paired participant resamples with GERA models, predictions and preprocessing fixed. `auc_bootstrap()` supplies type-7 percentile CIs and approximate two-sided Wald tests from the SD of paired AUC differences. It can also compute apparent GERA inference, which must be kept separate from corrected inference.

## Clinical utility

```r
source("performance.R")
source("clinical_utility.R")
# d and p have the same participant row order; p has Base, PRS and Full columns.
auc <- auc_bootstrap(d$Y, p)
nri <- nri_bootstrap(d$Y, p)
nns <- nns_bootstrap(d$Y, p)
```

Categorical NRI uses the 10% threshold at each participant's observed follow-up, comparing PRS and Full against Base. Its event, non-event and total components use 10,000 fixed-model paired resamples and type-7 percentile CIs. NNS selects the highest-risk 10% and 20%, includes ties at the type-7 selection cutoff, recomputes cutoffs within resamples and uses type-1 percentile CIs. Fixed-time kernel NRI and score-to-surveillance-interval rules are not included. Seeds make new runs repeatable, but exact Monte Carlo endpoints depend on the original archived streams. The original NRI stream seeded once before GERA then PLCO; pass `seed = NULL` to preserve that stream across calls.

DCA source code is preserved in [`archive/with_DCA/`](archive/with_DCA/README.md) for the manuscript alternative with DCA. It is excluded from the current manuscript workflow.

## Calibration and scoring

`calibration_5y()` evaluates Full-model risk in five cohort-specific predicted-risk quintiles, fitting a time-only cloglog SCAM within each group. It returns apparent group estimates with approximate link-scale Wald intervals transformed to probabilities. `calibration_optimism(data, model, indices)` additionally corrects GERA group points using the same 500 resamples. Original internal quintile cutoffs stay fixed, with endpoints extended to 0 and 1. Each refitted Full model determines new group memberships. Refit group-specific time-only models separately in the training and original cohorts. Average the training-minus-original differences in estimated-minus-mean-predicted probability and subtract them from original group estimates, keeping the original predicted x coordinates. This is an adaptation of the bootstrap calibration framework, not a calibration intercept/slope estimate. Apparent CIs remain attached to apparent points only. The auxiliary fitting fallback checks the stable high-smoothing limit explicitly; unresolved fits stop.

`point_score()` applies each term separately and rounds after scaling AA to 10 points. `profile_grid()` creates 1,008 hypothetical combinations. `score_risk()` averages model predictions within integer scores at 3, 5, 7 and 10 years. Figure 1A/B uses hypothetical profiles; Figure 1C and Supplementary Table 5 use real participant score groups. No score-to-follow-up-interval assignment is implemented. Supplementary Figure 3 is calibration and Supplementary Figure 4 is absolute risk (previously Figure 5).

## Software and references

Prediction, model fitting and calibration require `scam` and `mgcv`. Optional imputation additionally requires `missForest`. Original GERA models used R 4.5.1, scam 1.2-21 and mgcv 1.9-4. Server inference used the documented R 4.4 runtime with those same scam/mgcv versions; preserve exact runtime metadata with each local run. Core resampling and summarization use base R and `parallel`.

- Noma et al. *Stat Med*. 2021;40:5691–5701. Two-stage bootstrap intervals for optimism-corrected discrimination. https://doi.org/10.1002/sim.9148
- Harrell FE Jr, Lee KL, Mark DB. *Stat Med*. 1996;15:361–387. Bootstrap validation and calibration framework. https://doi.org/10.1002/(SICI)1097-0258(19960229)15:4%3C361::AID-SIM168%3E3.0.CO;2-4
