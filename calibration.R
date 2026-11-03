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

# Source optimism.R first. Original approximate CIs are not corrected CIs.
calibration_correction_summary <- function(apparent, errors) {
  errors <- as.matrix(errors)
  stopifnot(nrow(apparent) == 5L, nrow(errors) == 5L,
    ncol(errors) >= 2L, all(is.finite(errors)))
  out <- apparent
  out$mean_optimism_error <- rowMeans(errors)
  out$corrected_5y <- out$estimated_5y - out$mean_optimism_error
  out$optimism_MCSE <- apply(errors, 1L, sd)/sqrt(ncol(errors))
  out$replicates <- ncol(errors)
  out
}

calibration_optimism <- function(data, model, indices, apparent = NULL) {
  stopifnot(length(indices) >= 2L)
  d <- data
  at5 <- d; at5$time_years <- 5
  p0 <- as.numeric(predict(model, at5, type = "response"))
  cuts <- c(0, unname(quantile(p0, c(.2, .4, .6, .8), type = 7)), 1)
  stopifnot(!anyDuplicated(cuts))
  evaluate <- function(dd, pp) {
    groups <- cut(pp, cuts, include.lowest = TRUE, labels = FALSE)
    stopifnot(!anyNA(groups), all(tabulate(groups, 5L) > 0L))
    do.call(rbind, lapply(1:5, function(g) {
      rows <- which(groups == g)
      dg <- dd[rows, , drop = FALSE]
      stopifnot(sum(dg$Y) > 0, sum(1 - dg$Y) > 0,
        min(dg$time_years) < 5, max(dg$time_years) > 5)
      a <- .calibration_aux_fit(dg)
      # SCAM needs multiple rows when constructing the spline prediction basis.
      pp5 <- predict(a$fit, data.frame(time_years = c(3, 5, 7, 10)),
                     type = "link", se.fit = TRUE)
      eta <- as.numeric(pp5$fit)[2]
      se <- as.numeric(pp5$se.fit)[2]
      risk <- a$fit$family$linkinv(eta)
      ci <- a$fit$family$linkinv(eta + c(-1, 1) * qnorm(.975) * se)
      data.frame(group = g, mean_predicted_5y = mean(pp[rows]),
        estimated_5y = risk, lower = ci[1], upper = ci[2],
        error = risk - mean(pp[rows]), optimizer = a$optimizer, warnings = a$warnings)
    }))
  }
  original <- evaluate(d, p0)
  if (is.null(apparent)) apparent <- original[, c("group", "mean_predicted_5y",
    "estimated_5y", "lower", "upper")]
  stopifnot(nrow(apparent) == 5L,
    max(abs(apparent$mean_predicted_5y - original$mean_predicted_5y)) < 1e-10,
    max(abs(apparent$estimated_5y - original$estimated_5y)) < 1e-10)
  form <- formula(model)
  s <- mgcv::s
  environment(form) <- environment()
  op <- matrix(NA_real_, 5L, length(indices))
  diagnostics <- vector("list", length(indices))
  for (b in seq_along(indices)) {
    ix <- indices[[b]]
    stopifnot(length(ix) == nrow(d), all(ix %in% seq_len(nrow(d))))
    db <- d[ix, , drop = FALSE]
    a <- .optimism_fit(form, db)
    p5 <- as.numeric(predict(a$fit, at5, type = "response"))
    training <- evaluate(db, p5[ix])
    test <- evaluate(d, p5)
    op[, b] <- training$error - test$error
    diagnostics[[b]] <- list(full_optimizer = a$optimizer,
      full_warnings = a$warnings, training = training, original_test = test)
  }
  list(summary = calibration_correction_summary(apparent, op),
    cuts = cuts, optimism = op, diagnostics = diagnostics)
}

# Checked auxiliary fitting, including the explicitly verified smoothing limit.
.calibration_aux_fit <- function(dd) {
 s <- mgcv::s
 auxform <- Y ~ s(time_years, bs="mpi", k=7)
 warnings <- character()
 capture <- function(...)withCallingHandlers(tryCatch(scam(...),error=function(e)e),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 f <- capture(auxform,data=dd,family=binomial('cloglog'));optimizer <- 'bfgs/newton'
 if(!.optimism_fit_ok(f)) {
  f <- capture(auxform,data=dd,family=binomial('cloglog'),optimizer=c('optim','newton'),optim.method=c('L-BFGS-B','fd'),control=list(optim=list(maxit=1000,factr=10000)))
  optimizer <- 'optim L-BFGS-B/newton'
 }
 if(!.optimism_fit_ok(f)) {
  a <- tryCatch(.optimism_fit(auxform,dd),error=function(e)e)
  if(inherits(a,'error')) {
   if(!identical(conditionMessage(a),'Scalar optimum reached upper search boundary; requires inspection'))stop(a)
   # Explicitly verify the stable high-smoothing limit; reject numerical nonconvergence.
   logs <- c(10,12.5,15,17.5,20)
   limit_fits <- lapply(logs,function(z)capture(auxform,data=dd,family=binomial('cloglog'),sp=exp(z)))
   stopifnot(all(vapply(limit_fits,function(q)inherits(q,'scam')&&isTRUE(q$conv)&&all(is.finite(coef(q))),logical(1))))
   risks <- sapply(limit_fits,function(q)as.numeric(predict(q,data.frame(time_years=c(3,5,7,10)),type='response')))
   aics <- vapply(limit_fits,AIC,numeric(1))
   stopifnot(max(apply(risks,1,function(x)diff(range(x))))<1e-5,diff(range(aics))<1e-3)
   f <- limit_fits[[3]];attr(f,'scalar_smoothing_verified') <- TRUE
   optimizer <- 'verified stable high-smoothing limit; log(sp)=15'
   warnings <- c(warnings,'Scalar upper boundary inspected; stable predictions and AIC verified at log(sp)=10 through 20')
  } else {f <- a$fit;optimizer <- paste('additional fallback:',a$optimizer);warnings <- c(warnings,a$warnings)}
 }
 stopifnot(.optimism_fit_ok(f))
 grid <- sort(unique(c(seq(min(dd$time_years),max(dd$time_years),length.out=101),5)))
 curve <- as.numeric(predict(f,data.frame(time_years=grid),type='response'))
 stopifnot(all(is.finite(curve)),min(diff(curve))>=-1e-10)
 list(fit=f,optimizer=optimizer,warnings=paste(unique(warnings),collapse=' | '))
}
