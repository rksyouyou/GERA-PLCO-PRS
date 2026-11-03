# Net benefit for binary outcomes observed at each participant's examination.
decision_curve <- function(y, predictions, thresholds = (10:200)/1000) {
  p <- as.matrix(predictions)
  stopifnot(is.numeric(y) || is.logical(y), length(y) == nrow(p),
    !anyNA(y), all(y %in% 0:1), any(y == 1), any(y == 0),
    is.numeric(p), ncol(p) > 0L, all(is.finite(p)), all(p >= 0 & p <= 1),
    !is.null(colnames(p)), !anyDuplicated(colnames(p)),
    all(nzchar(colnames(p))), !any(colnames(p) %in% c("All", "None")),
    is.numeric(thresholds), length(thresholds) > 0L,
    all(is.finite(thresholds)), all(thresholds > 0 & thresholds < 1),
    !anyDuplicated(thresholds), !is.unsorted(thresholds))
  n <- length(y)
  events <- sum(y)
  prevalence <- mean(y)
  summarize <- function(model, score = NULL) {
    do.call(rbind, lapply(thresholds, function(pt) {
      high <- if (model == "All") rep(TRUE, n) else if (model == "None")
        rep(FALSE, n) else score >= pt
      tp <- sum(high & y == 1)
      fp <- sum(high & y == 0)
      nb <- tp/n - fp/n * pt/(1 - pt)
      data.frame(model = model, threshold = pt, N = n, AN = events,
        prevalence = prevalence, selected_N = sum(high), TP = tp, FP = fp,
        FN = events - tp, TN = n - events - fp, NB = nb, sNB = nb/prevalence)
    }))
  }
  do.call(rbind, c(lapply(colnames(p), function(m) summarize(m, p[, m])),
    list(summarize("All"), summarize("None"))))
}

# Shared participant indices should be the same ones used for the AUC correction.
# Source optimism.R first. This returns corrected point curves, without CIs/P.
decision_curve_optimism <- function(data, models, indices,
                                    thresholds = (10:200)/1000) {
  keys <- c("Base", "PRS", "Full")
  stopifnot(all(keys %in% names(models)), length(indices) >= 2L)
  models <- models[keys]
  forms <- lapply(models, formula)
  s <- mgcv::s
  for (j in seq_along(forms)) environment(forms[[j]]) <- environment()
  predict_all <- function(fits, d) sapply(fits, function(f)
    as.numeric(predict(f, d, type = "response")))
  apparent <- decision_curve(data$Y, predict_all(models, data), thresholds)
  rows <- apparent$model %in% keys
  op <- matrix(NA_real_, sum(rows), length(indices))
  for (b in seq_along(indices)) {
    ix <- indices[[b]]
    stopifnot(length(ix) == nrow(data), all(ix %in% seq_len(nrow(data))))
    db <- data[ix, , drop = FALSE]
    fits <- lapply(forms, function(f) .optimism_fit(f, db)$fit)
    po <- predict_all(fits, data)
    train <- decision_curve(db$Y, po[ix, , drop = FALSE], thresholds)
    test <- decision_curve(data$Y, po, thresholds)
    op[, b] <- train$NB[rows] - test$NB[rows]
  }
  out <- apparent[, c("model", "threshold", "prevalence")]
  out$apparent_NB <- apparent$NB
  out$mean_optimism_NB <- 0
  out$mean_optimism_NB[rows] <- rowMeans(op)
  out$corrected_NB <- out$apparent_NB - out$mean_optimism_NB
  out$corrected_sNB <- out$corrected_NB / mean(data$Y)
  out$replicates <- length(indices)
  out
}
