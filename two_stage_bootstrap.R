# Source fit_models.R and optimism.R before deliberately running outer refits.
# No analysis runs when this file is sourced. Keep checkpoints outside the repo.
outer_streams <- function(n, seed) {
  stopifnot(length(n) == 1L, n >= 1L, n == as.integer(n))
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv)
    else if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(seed)
  streams <- vector("list", n)
  streams[[1]] <- .Random.seed
  if (n > 1L) for (i in 2:n)
    streams[[i]] <- parallel::nextRNGStream(streams[[i - 1L]])
  streams
}

.bootstrap_hash <- function(x) {
  path <- tempfile()
  on.exit(unlink(path))
  saveRDS(x, path, version = 2)
  unname(tools::md5sum(path))
}

.bootstrap_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  on.exit(unlink(tmp))
  saveRDS(x, tmp, version = 2)
  if (!file.rename(tmp, path)) stop("Checkpoint rename failed.")
}

# One global outer ID per stream. Disjoint IDs can run as independent jobs.
# Fits within both levels share participant indices across all three models.
two_stage_outer <- function(data, stream, global_id, B_inner = 200L,
                            checkpoint = NULL) {
  forms <- model_formulas()
  keys <- names(forms)
  d <- data[, unique(unlist(lapply(forms, all.vars))), drop = FALSE]
  stopifnot(!anyNA(d), length(unique(d$Y)) == 2L,
            B_inner >= 2L, B_inner == as.integer(B_inner),
            length(stream) == 7L, global_id >= 1L)
  cfg <- list(protocol = "three-model-two-stage-v1", B_inner = as.integer(B_inner),
    global_id = global_id, stream = stream, data_hash = .bootstrap_hash(d),
    fit_hash = .bootstrap_hash(list(body(.optimism_fit), body(.optimism_fit_ok),
      body(.optimism_outer_ok), body(.optimism_auc),
      lapply(forms, function(f) paste(deparse(f), collapse = " ")))),
    R = R.version.string, scam = as.character(packageVersion("scam")),
    mgcv = as.character(packageVersion("mgcv")))
  if (!is.null(checkpoint)) {
    dir.create(dirname(checkpoint), recursive = TRUE, showWarnings = FALSE)
    lock <- paste0(checkpoint, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) stop("Checkpoint is locked.")
    on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  }
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  assign(".Random.seed", stream, .GlobalEnv)
  ix <- sample.int(nrow(d), nrow(d), replace = TRUE)
  od <- d[ix, , drop = FALSE]
  z <- if (!is.null(checkpoint) && file.exists(checkpoint)) readRDS(checkpoint) else
    list(config = cfg, global_id = global_id, outer_indices = ix, completed = 0L,
      apparent = setNames(rep(NA_real_, 3L), keys),
      training = matrix(NA_real_, B_inner, 3L, dimnames = list(NULL, keys)),
      test = matrix(NA_real_, B_inner, 3L, dimnames = list(NULL, keys)),
      diagnostics = list())
  stopifnot(identical(z$config, cfg), identical(z$outer_indices, ix))
  save_progress <- function() if (!is.null(checkpoint)) .bootstrap_atomic(z, checkpoint)
  diagnostic <- function(a, k, m) list(inner = k, model = m,
    optimizer = a$optimizer, warnings = a$warnings, sp = a$fit$sp,
    converged = .optimism_fit_ok(a$fit))
  for (m in keys) if (is.na(z$apparent[[m]])) {
    a <- .optimism_fit(forms[[m]], od)
    z$apparent[[m]] <- .optimism_auc(od$Y, predict(a$fit, od, type = "response"))
    z$diagnostics[[length(z$diagnostics) + 1L]] <- diagnostic(a, 0L, m)
    save_progress()
  }
  inner_seed <- stream
  for (k in seq_len(B_inner)) {
    inner_seed <- parallel::nextRNGSubStream(inner_seed)
    if (k <= z$completed) next
    assign(".Random.seed", inner_seed, .GlobalEnv)
    jj <- sample.int(nrow(od), nrow(od), replace = TRUE)
    id <- od[jj, , drop = FALSE]
    # Commit an inner iteration only after all three model fits succeed.
    train <- test <- setNames(numeric(3L), keys)
    diagnostics <- list()
    for (m in keys) {
      a <- .optimism_fit(forms[[m]], id)
      p <- as.numeric(predict(a$fit, od, type = "response"))
      train[[m]] <- .optimism_auc(id$Y, p[jj])
      test[[m]] <- .optimism_auc(od$Y, p)
      diagnostics[[m]] <- diagnostic(a, k, m)
    }
    z$training[k, ] <- train
    z$test[k, ] <- test
    z$diagnostics <- c(z$diagnostics, unname(diagnostics))
    z$completed <- as.integer(k)
    if (k %% 5L == 0L || k == B_inner) save_progress()
  }
  z$corrected <- z$apparent - colMeans(z$training - z$test)
  z$ok <- all(is.finite(z$corrected)) && all(vapply(z$diagnostics,
    function(q) isTRUE(q$converged), logical(1)))
  stopifnot(z$ok)
  save_progress()
  z
}

# point must be the ORIGINAL 500-resample corrected estimates, not outer means.
# Only completed successful rows may enter corrected_draws. Selection is explicit.
two_stage_summary <- function(point, corrected_draws, global_ids,
                               target = 5000L, selected_ids = NULL) {
  keys <- c("Base", "PRS", "Full")
  x <- as.matrix(corrected_draws)
  stopifnot(all(keys %in% colnames(x)), all(keys %in% names(point)),
    nrow(x) == length(global_ids), !anyDuplicated(global_ids),
    all(is.finite(x)), all(is.finite(point[keys])),
    target >= 2L, target == as.integer(target))
  if (is.null(selected_ids)) selected_ids <- head(sort(global_ids), target)
  stopifnot(length(selected_ids) >= 2L, length(selected_ids) <= target,
    !anyDuplicated(selected_ids), all(selected_ids %in% global_ids))
  x <- x[match(selected_ids, global_ids), keys, drop = FALSE]
  ci <- apply(x, 2, quantile, c(.025, .975), type = 7)
  pairs <- combn(keys, 2L)
  comparisons <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    a <- pairs[1, j]; b <- pairs[2, j]
    delta <- point[[b]] - point[[a]]
    differences <- x[, b] - x[, a]
    se <- sd(differences)
    limits <- quantile(differences, c(.025, .975), type = 7, names = FALSE)
    data.frame(comparator = a, model = b, delta_AUC = delta,
      lower = limits[1], upper = limits[2], bootstrap_SE = se,
      P_two_sided = if (se > 0) 2 * pnorm(-abs(delta/se)) else NA_real_)
  }))
  list(status = if (nrow(x) == target) "FINAL" else "INTERIM",
    outer_n = nrow(x), selected_ids = selected_ids,
    AUC = data.frame(model = keys, AUC = unname(point[keys]),
      lower = ci[1, ], upper = ci[2, ], row.names = NULL),
    comparisons = comparisons)
}
