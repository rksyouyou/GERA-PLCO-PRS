# Term-wise integer score; the advanced-adenoma contribution is 10 points.
point_score <- function(newdata, parameters = load_parameters()) {
  spec <- parameters$Full
  d <- newdata
  if (!"time_years" %in% names(d)) d$time_years <- 5
  d <- prepare_predictors(d, spec)
  X <- model.matrix(as.formula(spec$formula), d,
    contrasts.arg = lapply(spec$levels, function(x) contr.treatment(x, base = 1)))
  terms <- setdiff(colnames(X), "(Intercept)")
  b <- spec$coefficients[terms]
  rowSums(round(sweep(X[, terms, drop = FALSE], 2,
                      b * 10/spec$coefficients["basecolo_outAA"], "*")))
}

profile_grid <- function() {
  cutoffs <- c(0, .05, .15, .30, .70, .85, .95, 1)
  midpoints <- qnorm((head(cutoffs, -1) + tail(cutoffs, -1))/2)
  expand.grid(age_at_base = c(55, 65, 75), sex = c("Female", "Male"),
              bmi = c("<25", "25-29.9", "30+"), colo_fh = c("No", "Yes"),
              basecolo_out = c("Normal", "Aden1-2", "Aden3+", "AA"),
              prs_z = midpoints, stringsAsFactors = FALSE)
}

score_risk <- function(parameters = load_parameters(), years = c(3, 5, 7, 10)) {
  d <- profile_grid()
  d$risk_score <- point_score(d, parameters)
  do.call(rbind, lapply(years, function(t) {
    d$time_years <- t
    d$probability <- predict_an(d, parameters = parameters)
    a <- aggregate(probability ~ risk_score, d, mean)
    a$years <- t
    a[, c("years", "risk_score", "probability")]
  }))
}
