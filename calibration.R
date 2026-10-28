# Five-year calibration in quintiles of fixed Full-model risk.
calibration_5y <- function(data, parameters = load_parameters()) {
  stopifnot(!anyNA(data$Y), all(data$Y %in% 0:1))
  d <- data
  at5 <- d; at5$time_years <- 5
  p <- predict_an(at5, parameters = parameters)
  breaks <- unique(quantile(p, seq(0, 1, .2), type = 7))
  if (length(breaks) != 6L) stop("Five distinct risk quintiles are required.")
  group <- cut(p, breaks, include.lowest = TRUE, labels = FALSE)
  s <- mgcv::s
  do.call(rbind, lapply(1:5, function(g) {
    dd <- d[group == g, ]
    f <- scam::scam(Y ~ s(time_years, bs = "mpi", k = 7), data = dd,
                    family = binomial("cloglog"))
    converged <- isTRUE(f$conv) && identical(f$bfgs.info$conv, "Full convergence")
    if (!converged) {
      f <- scam::scam(Y ~ s(time_years, bs = "mpi", k = 7), data = dd,
        family = binomial("cloglog"), optimizer = c("optim", "newton"),
        optim.method = c("L-BFGS-B", "fd"),
        control = list(optim = list(maxit = 1000, factr = 10000)))
      converged <- isTRUE(f$conv) && f$optim.info$conv %in% c("Full convergence", "Successful completion")
    }
    if (!converged) stop("Quintile model did not converge.")
    pp <- predict(f, data.frame(time_years = c(3, 5, 7, 10)), type = "link", se.fit = TRUE)
    pp$fit <- pp$fit[2]; pp$se.fit <- pp$se.fit[2]
    ci <- f$family$linkinv(as.numeric(pp$fit) + c(-1, 1) * qnorm(.975) * as.numeric(pp$se.fit))
    data.frame(group = g, mean_predicted_5y = mean(p[group == g]),
      estimated_5y = f$family$linkinv(as.numeric(pp$fit)), lower = ci[1], upper = ci[2])
  }))
}
