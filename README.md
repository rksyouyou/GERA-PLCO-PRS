# GERA–PLCO PRS

R functions for the final SCAM risk models and their evaluation in the GERA–PLCO study of metachronous advanced neoplasia (AN).

The repository contains source code and fixed model parameters only. No participant records, genotypes, PRS observations, imputed data, fitted-model objects, bootstrap samples, or study result files are included.

## Files

| File | Purpose |
| --- | --- |
| `model_parameters.R` | Fixed coefficients and the minimal time-spline specification for five GERA models |
| `predict.R` | Absolute-risk predictions using the fixed GERA models |
| `fit_models.R` | Clinical imputation, cohort-level PRS standardization, and optional development fitting |
| `performance.R` | AUC estimates, paired participant bootstrap intervals, and two-sided Wald comparisons |
| `optimism.R` | GERA bootstrap optimism correction with all five models refitted |
| `clinical_utility.R` | Observed-follow-up categorical NRI and NNS |
| `calibration.R` | Five-year calibration with a time-only SCAM in each predicted-risk quintile |
| `point_score.R` | Integer scores and score-level risk summaries for the defined profile grid |

## Prediction

Install `scam` and its dependencies. Run from the repository directory:

```r
source("predict.R")
parameters <- load_parameters()
# newdata is a locally supplied data.frame with the columns defined below.
risk <- predict_an(newdata, model = "Full", parameters = parameters)
all_models <- predict_models(newdata, parameters)
```

| Column | Required values |
| --- | --- |
| `age_at_base` | Age in years; centered at 50 internally |
| `sex` | `Female` or `Male` |
| `bmi` | `<25`, `25-29.9`, or `30+` |
| `colo_fh` | Family history of colorectal cancer: `No` or `Yes` |
| `basecolo_out` | `Normal`, `Aden1-2`, `Aden3+`, or `AA` |
| `prs_z` | Calibrated PRS, standardized within the relevant final analysis population |
| `time_years` | Follow-up interval in years |

Aden1-2 and Aden3+ denote 1–2 and at least 3 nonadvanced conventional adenomas, respectively; AA denotes advanced adenoma. Prediction requires complete values for the predictors in the selected model. A PRS must be calibrated and placed on the study scale before use. The function does not calculate a PRS from genotypes or standardize a single patient. GERA and PLCO were standardized separately across their final study populations, without within-platform standardization.

The five models use a complementary log-log link and a nondecreasing SCOP-spline, `s(time_years, bs = "mpi", k = 7)`. Base includes baseline findings; Model1 adds age and sex; Model2 also adds BMI and family history. Model3 includes baseline findings, PRS, and their interaction. Full adds age, sex, BMI, and family history to Model3. GERA coefficients remain fixed for PLCO evaluation.

`model_parameters.R` contains coefficients, category levels, knots, and the centering constants needed for prediction. It contains no training rows or subject-level model components. Predictions use the transformed SCAM coefficients and the package's spline prediction method. The functions have been checked against the saved final models in the development and validation populations and at 3, 5, 7, and 10 years.

## Analysis with locally held data

The code starts from the eligible study populations and an already technically corrected/calibrated PRS. Upstream genotype processing, source-system extraction, and cohort selection require restricted inputs and are outside this code-only release. In the study, GERA excluded all baseline CSSPs; PLCO used the 180-day FSG limit with missing intervals retained. The PLCO sensitivity population removed Hyperplastic Polyp and Benign Polyp, NOS only from the Normal group.

`prepare_cohort()` takes `Y` (0/1 examination AN status), `bs_time` (days), `prs` (calibrated PRS), and the clinical/baseline columns above. It imputes age, sex, BMI, and family history with `missForest` and standardizes PRS separately in each supplied population. It does not change cohort eligibility. Use `fit_models()` only to refit the development models deliberately; fixed-model prediction does not refit anything.

```r
source("performance.R")
source("clinical_utility.R")
# d is a locally held, prepared cohort; preserve the same participant row order.
p <- predict_models(d, parameters)
auc <- auc_bootstrap(d$Y, p)
nri <- nri_bootstrap(d$Y, p)
nns <- nns_bootstrap(d$Y, p)
```

AUC and NRI use 10,000 paired participant resamples and type-7 percentile intervals. AUC comparisons use the standard deviation of paired bootstrap differences for two-sided Wald tests. NRI uses the 10% probability threshold at each participant's observed follow-up examination. Its event, non-event, and total components compare Full and Models1–3 with Base. NNS selects the highest-risk 10% or 20%, includes ties at the type-7 quantile cutoff, and uses type-1 percentile intervals. The same sampled rows are used for all models within each bootstrap replicate; models and preprocessing remain fixed.

The seed arguments make new bootstrap runs repeatable. Exact Monte Carlo endpoints can differ from the manuscript's archived draws, which are not distributed. For the original NRI stream across GERA then PLCO, set the seed once and call both cohorts with `seed = NULL`.

## Bootstrap optimism correction

The GERA analysis additionally uses 500 ordinary participant bootstrap samples shared by all five models (seed 123). Each model is refitted in each sample, including the time smooth and smoothing parameter, with clinical imputations and PRS preprocessing held fixed. Optimism is the mean difference between a refitted model's bootstrap-sample AUC and its AUC in the original cohort. The corrected AUC subtracts this mean from the original apparent AUC.

```r
source("fit_models.R")
source("optimism.R")
# gera is a locally held, prepared development cohort.
models <- fit_models(gera)
validation <- optimism_bootstrap(gera, models, B = 500L, seed = 123L)
validation$AUC
validation$differences
```

This is a separate model-refitting analysis from `auc_bootstrap()`, which holds predictions fixed. Corrected AUCs and paired differences are reported as point estimates. Their Monte Carlo standard errors quantify simulation precision; they are not statistical standard errors or confidence intervals. The original apparent-AUC inference and PLCO external-validation results are retained separately.

The fitting routine checks inner and outer convergence. If default fitting fails, it retries the same model and independently checks the single smoothing parameter using the same binomial UBRE criterion. Unresolved fits stop the calculation rather than silently excluding a sample. Setting `workers` above 1 enables parallel refits on macOS/Linux; Windows uses sequential refits. Functions return results in memory and do not write participant data or output files.

## Calibration and point score

`calibration_5y()` fits time-only SCAMs within five quintiles of Full-model five-year predicted risk, with a link-scale Wald interval transformed back to probability. `point_score()` applies each term's contribution separately and rounds after scaling the AA coefficient to 10 points. `profile_grid()` creates the study's 1,008 hypothetical combinations in memory; `score_risk()` averages model predictions within each integer score. These profiles are not participants. No functions write files automatically.

## Dependencies and access

Prediction and calibration require `scam` (and its dependency `mgcv`); optional clinical imputation requires `missForest`. Model parameters were fitted with R 4.5.1 and scam 1.2-21. Individual-level study data are not distributed through this repository and remain subject to the relevant institutional access requirements. Keep all local input and output files outside the repository.
