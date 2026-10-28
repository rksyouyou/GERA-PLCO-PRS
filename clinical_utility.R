# Observed-follow-up categorical NRI and number needed to screen.
# Source performance.R before using these functions.
nri_bootstrap <- function(y, predictions, threshold = .10,
                          B = 10000L, seed = 123L) {
  p <- validate_predictions(y, predictions)
  models <- c("Full", "Model1", "Model2", "Model3")
  stopifnot(all(c("Base", models) %in% colnames(p)), B >= 2L,
            threshold > 0, threshold < 1)
  high <- p >= threshold
  movement <- sapply(models, function(m) as.integer(high[, m]) - as.integer(high[, "Base"]))
  if (!is.null(seed)) set.seed(seed)
  draws <- array(NA_real_, c(B, length(models), 3L))
  for (b in seq_len(B)) {
    i <- sample.int(length(y), length(y), replace = TRUE)
    yy <- y[i]; mm <- movement[i, , drop = FALSE]
    ev <- colMeans(mm[yy == 1, , drop = FALSE])
    ne <- -colMeans(mm[yy == 0, , drop = FALSE])
    draws[b, , ] <- cbind(ev, ne, ev + ne)
  }
  if (anyNA(draws)) stop("A bootstrap sample lacked events or non-events.")
  do.call(rbind, lapply(seq_along(models), function(j) {
    m <- movement[, j]
    pr <- c(mean(m[y == 1] == 1), mean(m[y == 1] == -1),
            mean(m[y == 0] == 1), mean(m[y == 0] == -1))
    ev <- pr[1] - pr[2]; ne <- pr[4] - pr[3]
    ci <- apply(draws[, j, ], 2, quantile, probs = c(.025, .975), type = 7)
    data.frame(comparison = paste(models[j], "vs Base"),
      P_up_AN = pr[1], P_down_AN = pr[2], P_up_noAN = pr[3], P_down_noAN = pr[4],
      event_NRI = ev, event_lower = ci[1, 1], event_upper = ci[2, 1],
      nonevent_NRI = ne, nonevent_lower = ci[1, 2], nonevent_upper = ci[2, 2],
      total_NRI = ev + ne, total_lower = ci[1, 3], total_upper = ci[2, 3])
  }))
}

nns_bootstrap <- function(y, predictions, fractions = c(.10, .20),
                          B = 10000L, seed = 123L) {
  p <- validate_predictions(y, predictions)
  stopifnot(B >= 2L, all(fractions > 0 & fractions < 1))
  nns <- function(y, score, fraction) {
    selected <- score >= quantile(score, 1 - fraction, type = 7, names = FALSE)
    events <- sum(y[selected])
    if (events) sum(selected)/events else Inf
  }
  if (!is.null(seed)) set.seed(seed)
  draws <- array(NA_real_, c(B, ncol(p), length(fractions)))
  for (b in seq_len(B)) {
    i <- sample.int(length(y), length(y), replace = TRUE)
    for (j in seq_len(ncol(p))) for (k in seq_along(fractions))
      draws[b, j, k] <- nns(y[i], p[i, j], fractions[k])
  }
  do.call(rbind, lapply(seq_len(ncol(p)), function(j)
    do.call(rbind, lapply(seq_along(fractions), function(k) {
      ci <- quantile(draws[, j, k], c(.025, .975), type = 1, names = FALSE)
      data.frame(model = colnames(p)[j], target_selection_fraction = fractions[k],
                 NNS = nns(y, p[, j], fractions[k]), lower = ci[1], upper = ci[2])
    }))))
}
