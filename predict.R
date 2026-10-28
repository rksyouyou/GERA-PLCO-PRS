# Fixed-model prediction using standardized PRS and complete clinical inputs.
load_parameters <- function(path = "model_parameters.R") {
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  env$scam_parameters
}

prepare_predictors <- function(newdata, parameters) {
  required <- all.vars(as.formula(parameters$formula))
  required <- unique(c(setdiff(required, "age_centered"), "time_years"))
  d <- as.data.frame(newdata)
  if ("age_centered" %in% all.vars(as.formula(parameters$formula))) {
    if (!"age_at_base" %in% names(d)) stop("age_at_base is required, in years.")
    d$age_centered <- d$age_at_base - 50
  }
  if (!all(required %in% names(d)))
    stop("Missing predictors: ", paste(setdiff(required, names(d)), collapse = ", "))
  if (!nrow(d)) stop("Provide at least one row.")
  for (v in names(parameters$levels)) {
    if (anyNA(d[[v]]) || any(!as.character(d[[v]]) %in% parameters$levels[[v]]))
      stop("Invalid or missing category in ", v)
    d[[v]] <- factor(d[[v]], levels = parameters$levels[[v]])
  }
  for (v in intersect(c("age_at_base", "prs_z", "time_years"), names(d))) {
    if (!is.numeric(d[[v]]) || any(!is.finite(d[[v]])))
      stop(v, " must contain finite numeric values.")
  }
  if (any(d$time_years <= 0)) stop("time_years must be positive.")
  d
}

predict_an <- function(newdata, model = "Full", parameters = load_parameters(),
                       type = c("response", "link")) {
  type <- match.arg(type)
  if (!model %in% names(parameters)) stop("Unknown model.")
  if (!requireNamespace("scam", quietly = TRUE)) stop("Install the scam package.")
  spec <- parameters[[model]]
  d <- prepare_predictors(newdata, spec)
  contrast <- lapply(spec$levels, function(x) contr.treatment(x, base = 1))
  X <- model.matrix(as.formula(spec$formula), d, contrasts.arg = contrast)
  sm <- structure(c(spec$time, list(by = "NA", dim = 1L)),
                  class = c("mpi.smooth", "mgcv.smooth"))
  # The mpi method requires a matrix even when predicting one observation.
  basis_data <- if (nrow(d) == 1L) d[c(1L, 1L), , drop = FALSE] else d
  time_basis <- mgcv::PredictMat(sm, basis_data)[seq_len(nrow(d)), , drop = FALSE]
  beta <- spec$coefficients_transformed
  stopifnot(identical(colnames(X), names(beta)[seq_len(ncol(X))]),
            ncol(X) + ncol(time_basis) == length(beta))
  eta <- as.numeric(cbind(X, time_basis) %*% beta)
  if (type == "link") eta else binomial("cloglog")$linkinv(eta)
}

predict_models <- function(newdata, parameters = load_parameters()) {
  values <- vapply(names(parameters), function(m) predict_an(newdata, m, parameters),
                   numeric(nrow(newdata)))
  matrix(values, nrow = nrow(newdata), dimnames = list(NULL, names(parameters)))
}
