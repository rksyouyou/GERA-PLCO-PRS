# Development models; input must already satisfy the study eligibility criteria.
prepare_cohort <- function(data, seed = 123L) {
  vars <- c("age_at_base", "sex", "bmi", "colo_fh")
  required <- c(vars, "Y", "bs_time", "basecolo_out", "prs")
  stopifnot(all(required %in% names(data)))
  d <- data[, required, drop = FALSE]
  levels <- list(sex = c("Female", "Male"), bmi = c("<25", "25-29.9", "30+"),
                 colo_fh = c("No", "Yes"),
                 basecolo_out = c("Normal", "Aden1-2", "Aden3+", "AA"))
  for (v in names(levels)) {
    if (any(!is.na(d[[v]]) & !as.character(d[[v]]) %in% levels[[v]]))
      stop("Invalid category in ", v)
    d[[v]] <- factor(d[[v]], levels = levels[[v]])
  }
  stopifnot(!anyNA(d$Y), all(d$Y %in% 0:1), !anyNA(d$basecolo_out),
            all(is.finite(d$prs)), all(is.finite(d$bs_time)), all(d$bs_time > 0))
  set.seed(seed)
  if (anyNA(d[, vars])) {
    if (!requireNamespace("missForest", quietly = TRUE)) stop("Install missForest.")
    d[, vars] <- missForest::missForest(d[, vars], verbose = FALSE)$ximp
  }
  stopifnot(!anyNA(d), sd(d$prs) > 0)
  d$prs_z <- as.numeric(scale(d$prs))
  d$age_centered <- d$age_at_base - 50
  d$time_years <- d$bs_time / 365.25
  d
}

fit_models <- function(gera) {
  if (!requireNamespace("scam", quietly = TRUE)) stop("Install scam.")
  s <- mgcv::s
  formulas <- list(
    Base = Y ~ s(time_years, bs = "mpi", k = 7) + basecolo_out,
    Model1 = Y ~ s(time_years, bs = "mpi", k = 7) + basecolo_out + age_centered + sex,
    Model2 = Y ~ s(time_years, bs = "mpi", k = 7) + basecolo_out + age_centered + sex + bmi + colo_fh,
    Model3 = Y ~ s(time_years, bs = "mpi", k = 7) + basecolo_out * prs_z,
    Full = Y ~ s(time_years, bs = "mpi", k = 7) + basecolo_out * prs_z + age_centered + sex + bmi + colo_fh)
  fits <- lapply(formulas, function(f)
    scam::scam(f, data = gera, family = binomial("cloglog")))
  stopifnot(all(vapply(fits, function(f)
    isTRUE(f$conv) && identical(f$bfgs.info$conv, "Full convergence"), logical(1))))
  fits
}
