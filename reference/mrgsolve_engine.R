# ============================================================================
# mrgsolve engine: the same interface as mrg_solve() in shared/ode_engine.R.
#
# Apps source this file when they run locally. If mrgsolve is installed and
# the model compiles, simulations go through mrgsolve; otherwise the app
# keeps the base-R engine. The browser build never includes this file.
# tests/ use it as the reference.
# ============================================================================

mrgsolve_ready <- function() requireNamespace("mrgsolve", quietly = TRUE)

.mrg_cache <- new.env()

mrgsolve_model <- function(path) {
  key <- normalizePath(path, mustWork = TRUE)
  if (is.null(.mrg_cache[[key]])) {
    .mrg_cache[[key]] <- mrgsolve::mread(tools::file_path_sans_ext(basename(path)),
                                         project = dirname(path), quiet = TRUE)
  }
  .mrg_cache[[key]]
}

#' Same contract as mrg_solve(); `model` is the object from mrg_read(), whose
#' path points at the .cpp file. Tolerances passed by the caller (tuned for
#' the base-R solver) are ignored: mrgsolve is fast enough to always run
#' tight, which also keeps it a clean reference.
mrg_solve_mrgsolve <- function(model, P = data.frame(row.names = 1), events = NULL, times,
                               init = NULL, ..., rtol_ref = 1e-8, atol_ref = 1e-12) {
  rtol <- rtol_ref
  atol <- atol_ref
  mod <- mrgsolve_model(model$path)
  if (!nrow(P)) P <- data.frame(row.names = 1)
  m <- nrow(P)
  idata <- P
  idata$ID <- seq_len(m)
  ev <- mrg_expand_events(events, m)
  data <- NULL
  if (nrow(ev)) {
    data <- do.call(rbind, lapply(seq_len(m), function(i) {
      e <- ev[is.na(ev$ID) | ev$ID == i, , drop = FALSE]
      if (!nrow(e)) return(NULL)
      data.frame(ID = i, time = e$time, cmt = e$cmt, amt = e$amt, rate = e$rate, evid = 1)
    }))
    data <- data[order(data$ID, data$time), ]
  }
  if (!is.null(init)) mod <- mrgsolve::init(mod, as.list(init))
  mod <- mrgsolve::update(mod, rtol = rtol, atol = atol, maxsteps = 1e7)
  times <- sort(unique(times))
  out <- if (is.null(data)) {
    mrgsolve::mrgsim_df(mod, idata = idata, tgrid = times, end = -1, add = times)
  } else {
    mrgsolve::mrgsim_df(mod, data = data, idata = idata, add = times, end = -1, obsonly = TRUE)
  }
  out <- out[!duplicated(out[, c("ID", "time")], fromLast = TRUE), ]
  out <- out[out$time %in% times, ]
  nt <- length(times)
  states <- array(NA_real_, c(nt, length(model$cmt), m), dimnames = list(NULL, model$cmt, NULL))
  for (k in seq_along(model$cmt)) states[, k, ] <- matrix(out[[model$cmt[k]]], nrow = nt)
  res <- list(time = times, states = states)
  for (nm in model$capture) res[[nm]] <- matrix(out[[nm]], nrow = nt)
  res
}
