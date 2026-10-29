# Bootstrap optimism correction for the five prespecified SCAM models.
# Clinical imputations and PRS preprocessing are held fixed during resampling.
suppressPackageStartupMessages(library(scam))
.optimism_auc <- function(y, p) {
    stopifnot(all(is.finite(p)), all(y %in% c(0, 1)), sum(y) > 0, sum(1 - y) > 0)
    n1 <- sum(y)
    n0 <- sum(1 - y)
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}

.optimism_outer_ok <- function(f) {
    if (isTRUE(attr(f, "scalar_smoothing_verified"))) 
        return(TRUE)
    if (!is.null(f$bfgs.info)) 
        return(identical(f$bfgs.info$conv, "Full convergence"))
    if (!is.null(f$optim.info)) 
        return(identical(f$optim.info$conv, "Successful completion"))
    FALSE
}

.optimism_fit_ok <- function(f) inherits(f, "scam") && isTRUE(f$conv) && .optimism_outer_ok(f) && all(is.finite(coef(f)))

.optimism_fit <- function(formula, data) {
    warnings <- character()
    capture_fit <- function(...) withCallingHandlers(tryCatch(scam(...), error = function(e) e), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    elapsed <- proc.time()[3]
    f <- capture_fit(formula, data = data, family = binomial("cloglog"))
    default_status <- if (inherits(f, "error")) 
        conditionMessage(f)
    else paste(f$conv, f$bfgs.info$conv)
    optimizer <- "bfgs/newton"
    if (!.optimism_fit_ok(f)) {
        sp_default <- if (inherits(f, "scam")) 
            f$sp
        else numeric()
        f <- capture_fit(formula, data = data, family = binomial("cloglog"), optimizer = c("optim", "newton"), 
            optim.method = c("Nelder-Mead", "fd"), control = scam.control(optim = list(maxit = 1000, reltol = 1e-10)))
        optimizer <- "optim Nelder-Mead/newton"
        sp_retry <- if (inherits(f, "scam")) 
            f$sp
        else numeric()
        scalar_fit <- function(logsp) {
            q <- capture_fit(formula, data = data, family = binomial("cloglog"), sp = exp(logsp))
            if (!inherits(q, "scam") || !isTRUE(q$conv) || !all(is.finite(coef(q)))) 
                return(NULL)
            q
        }
        criterion <- function(logsp) {
            q <- scalar_fit(logsp)
            if (is.null(q)) 
                return(Inf)
            AIC(q)/nrow(data) - 1
        }
        positive_sp <- c(sp_default, sp_retry)
        positive_sp <- positive_sp[is.finite(positive_sp) & positive_sp > 0]
        grid <- sort(unique(c(seq(-30, 20, by = 2), log(positive_sp))))
        values <- vapply(grid, criterion, numeric(1))
        best <- which.min(values)
        if (!is.finite(values[best])) 
            stop("No converged scalar-smoothing candidate")
        candidates <- list(scalar_fit(grid[best]))
        if (best > 1 && best < length(grid)) {
            opt <- optimize(criterion, grid[c(best - 1, best + 1)], tol = 1e-06)
            candidates[[length(candidates) + 1L]] <- scalar_fit(opt$minimum)
        }
        else if (best == length(grid)) {
            stop("Scalar optimum reached upper search boundary; requires inspection")
        }
        zero <- scalar_fit(-Inf)
        if (!is.null(zero)) 
            candidates[[length(candidates) + 1L]] <- zero
        candidates <- Filter(Negate(is.null), candidates)
        f <- candidates[[which.min(vapply(candidates, AIC, numeric(1)))]]
        attr(f, "scalar_smoothing_verified") <- TRUE
        optimizer <- "scalar UBRE grid/Brent/newton"
    }
    if (!.optimism_fit_ok(f)) 
        stop("Unresolved fit: ", default_status, "; retry: ", if (inherits(f, "error")) 
            conditionMessage(f)
        else paste(f$conv, f$optim.info$conv))
    list(fit = f, optimizer = optimizer, default_status = default_status, warnings = paste(unique(warnings), collapse = " | "), 
        seconds = unname(proc.time()[3] - elapsed))
}



optimism_summary <- function(apparent, draws, model_names = names(apparent)) {
    required <- c("replicate", "model", "training_AUC", "original_test_AUC", "optimism")
    stopifnot(all(required %in% names(draws)), !anyDuplicated(draws[, c("replicate", "model")]),
              all(is.finite(apparent)), setequal(unique(draws$model), model_names))
    B <- max(draws$replicate)
    stopifnot(B >= 2L, nrow(draws) == B * length(model_names),
              all(is.finite(draws$optimism)),
              max(abs(draws$optimism - draws$training_AUC + draws$original_test_AUC)) < 1e-12)
    op <- vapply(model_names, function(m) {
        q <- draws[draws$model == m, ]
        q <- q[order(q$replicate), ]
        stopifnot(identical(as.integer(q$replicate), seq_len(B)))
        q$optimism
    }, numeric(B))
    auc <- data.frame(model = model_names, apparent_AUC = unname(apparent[model_names]),
        mean_optimism = colMeans(op), corrected_AUC = unname(apparent[model_names]) - colMeans(op),
        optimism_MCSE = apply(op, 2, sd)/sqrt(B), replicates = B, row.names = NULL)
    pairs <- combn(model_names, 2)
    differences <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
        a <- pairs[1, j]; b <- pairs[2, j]; x <- op[, b] - op[, a]
        delta <- apparent[[b]] - apparent[[a]]
        data.frame(comparator = a, model = b, apparent_delta_AUC = delta,
            mean_delta_optimism = mean(x), corrected_delta_AUC = delta - mean(x),
            delta_optimism_MCSE = sd(x)/sqrt(B), replicates = B)
    }))
    list(AUC = auc, differences = differences)
}

optimism_bootstrap <- function(data, models, B = 500L, seed = 123L, workers = 1L) {
    model_names <- c("Base", "Model1", "Model2", "Model3", "Full")
    stopifnot(is.list(models), all(model_names %in% names(models)),
              length(B) == 1L, B >= 2L, B == as.integer(B),
              length(workers) == 1L, workers >= 1L, workers == as.integer(workers))
    models <- models[model_names]
    stopifnot(all(vapply(models, function(f) .optimism_fit_ok(f) &&
        f$family$link == "cloglog" && length(f$smooth) == 1L &&
        inherits(f$smooth[[1]], "mpi.smooth") && f$smooth[[1]]$bs.dim == 7L, logical(1))))
    forms <- lapply(models, formula)
    for (j in seq_along(forms)) environment(forms[[j]]) <- environment()
    s <- mgcv::s
    variables <- unique(unlist(lapply(forms, all.vars)))
    stopifnot(all(variables %in% names(data)))
    d <- data[, variables, drop = FALSE]
    stopifnot(nrow(d) > 1L, !anyNA(d), all(d$Y %in% 0:1), length(unique(d$Y)) == 2L)
    apparent <- vapply(models, function(f) .optimism_auc(d$Y,
        as.numeric(predict(f, d, type = "response"))), 0)
    old_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
        do.call(RNGkind, as.list(old_kind))
        if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    indices <- replicate(B, sample.int(nrow(d), nrow(d), replace = TRUE), simplify = FALSE)
    one <- function(b) {
        ix <- indices[[b]]; db <- d[ix, , drop = FALSE]
        do.call(rbind, lapply(model_names, function(m) {
            a <- .optimism_fit(forms[[m]], db); f <- a$fit
            pb <- as.numeric(predict(f, db, type = "response"))
            po <- as.numeric(predict(f, d, type = "response"))
            stopifnot(all(pb >= 0 & pb <= 1), all(po >= 0 & po <= 1),
                      max(abs(pb - po[ix])) < 1e-10)
            training <- .optimism_auc(db$Y, pb); original <- .optimism_auc(d$Y, po)
            data.frame(replicate = b, model = m, n = nrow(db), events = sum(db$Y),
                unique_participants = length(unique(ix)), training_AUC = training,
                original_test_AUC = original, optimism = training - original,
                inner_converged = isTRUE(f$conv), outer_converged = .optimism_outer_ok(f),
                optimizer = a$optimizer, smoothing_parameter = unname(f$sp),
                total_edf = sum(f$edf), warnings = a$warnings)
        }))
    }
    if (workers > 1L && .Platform$OS.type != "windows") {
        ans <- parallel::mclapply(seq_len(B), one, mc.cores = workers,
                                 mc.preschedule = FALSE, mc.set.seed = FALSE)
    } else ans <- lapply(seq_len(B), one)
    stopifnot(all(vapply(ans, is.data.frame, logical(1))))
    draws <- do.call(rbind, ans)
    stopifnot(all(draws$inner_converged & draws$outer_converged))
    summary <- optimism_summary(apparent, draws)
    list(config = list(B = as.integer(B), seed = seed, models = model_names,
        formulas = lapply(forms, function(f) paste(deparse(f), collapse = " ")),
        R = R.version.string, scam = as.character(packageVersion("scam")),
        preprocessing = "Clinical imputations and PRS preprocessing held fixed"),
        apparent = apparent, draws = draws, AUC = summary$AUC,
        differences = summary$differences)
}
