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
