# Rank AUC and paired participant bootstrap with predictions held fixed.
validate_predictions <- function(y, predictions) {
  p <- as.matrix(predictions)
  stopifnot(length(y) == nrow(p), !anyNA(y), all(y %in% 0:1),
            length(unique(y)) == 2L, all(is.finite(p)), all(p >= 0 & p <= 1),
            !is.null(colnames(p)), !anyDuplicated(colnames(p)))
  p
}

auc_plan <- function(p) {
  o <- order(p)
  list(order = o, ends = cumsum(rle(p[o])$lengths))
}

weighted_auc <- function(y, weights, plan) {
  i <- plan$order
  a <- weights[i] * y[i]; b <- weights[i] * (1 - y[i])
  if (!sum(a) || !sum(b)) return(NA_real_)
  ca <- cumsum(a)[plan$ends]; cb <- cumsum(b)[plan$ends]
  sum(diff(c(0, ca)) * (cb - diff(c(0, cb))/2)) / (sum(a) * sum(b))
}

auc_bootstrap <- function(y, predictions, B = 10000L, seed = 123L) {
  p <- validate_predictions(y, predictions)
  stopifnot(B >= 2L, ncol(p) >= 2L)
  plans <- lapply(seq_len(ncol(p)), function(j) auc_plan(p[, j]))
  observed <- vapply(plans, function(q) weighted_auc(y, rep(1, length(y)), q), 0)
  set.seed(seed)
  draws <- matrix(NA_real_, B, ncol(p), dimnames = list(NULL, colnames(p)))
  for (b in seq_len(B)) {
    w <- tabulate(sample.int(length(y), length(y), replace = TRUE), nbins = length(y))
    draws[b, ] <- vapply(plans, function(q) weighted_auc(y, w, q), 0)
  }
  if (anyNA(draws)) stop("A bootstrap sample lacked events or non-events.")
  ci <- apply(draws, 2, quantile, probs = c(.025, .975), type = 7)
  estimates <- data.frame(model = colnames(p), AUC = observed,
                          lower = ci[1, ], upper = ci[2, ], row.names = NULL)
  pairs <- combn(seq_len(ncol(p)), 2)
  contrasts <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(k) {
    i <- pairs[1, k]; j <- pairs[2, k]
    delta <- observed[j] - observed[i]
    d <- draws[, j] - draws[, i]; se <- sd(d)
    limits <- quantile(d, c(.025, .975), type = 7, names = FALSE)
    data.frame(comparator = colnames(p)[i], model = colnames(p)[j],
      delta_AUC = delta, lower = limits[1], upper = limits[2], bootstrap_SE = se,
      P_two_sided = if (se > 0) 2 * pnorm(abs(delta/se), lower.tail = FALSE) else NA_real_)
  }))
  list(AUC = estimates, comparisons = contrasts)
}
