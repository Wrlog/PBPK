# ============================================================================
# mrgsolve engine: runs models/pbpk.cpp through the same interface as
# pbpk_simulate() in R/pbpk_engine.R.
#
# The app sources this file when it runs locally. If mrgsolve is installed
# and the model compiles, simulations go through mrgsolve; otherwise the
# app keeps the matrix-exponential engine. The browser build never
# includes this directory, since a browser can't compile C++.
#
# tests/test_engine_vs_mrgsolve.R uses it as the reference.
# ============================================================================

mrgsolve_ready <- function() {
  requireNamespace("mrgsolve", quietly = TRUE)
}

.pbpk_mod <- NULL

pbpk_mrgsolve_model <- function(model_dir = "models") {
  if (is.null(.pbpk_mod)) {
    .pbpk_mod <<- mrgsolve::mread("pbpk", project = model_dir, quiet = TRUE)
  }
  .pbpk_mod
}

#' Parameter vector for models/pbpk.cpp from a pbpk_parameters() list
pbpk_to_mrgsolve_param <- function(p) {
  tissues <- c("lung", "brain", "heart", "kidney", "muscle", "skin", "adipose",
               "bone", "spleen", "gut", "liver", "rest")
  flows <- c("brain", "heart", "kidney", "muscle", "skin", "adipose", "bone",
             "spleen", "gut", "rest")
  c(
    stats::setNames(p$V[c(tissues, "art", "ven")], paste0("V_", c(tissues, "art", "ven"))),
    stats::setNames(p$Q[flows], paste0("Q_", flows)),
    Q_ha = p$Q[["hepatic_artery"]],
    stats::setNames(p$Kp[tissues], paste0("KP_", tissues)),
    BP = p$BP, FU = p$fu, CLINT = p$CLint, CLR = p$CLr,
    KA = p$ka, FA = p$fa, FG = p$fg
  )
}

#' Same contract as pbpk_simulate(), plasma only
pbpk_simulate_mrgsolve <- function(plist, regimen, times, model_dir = "models",
                                   rtol = 1e-8, atol = 1e-12) {
  mod <- pbpk_mrgsolve_model(model_dir)
  n <- length(plist)
  amt <- rep_len(regimen$amt, n)
  t_end <- max(times)
  reg <- expand_regimen(regimen, t_end)

  idata <- as.data.frame(do.call(rbind, lapply(plist, pbpk_to_mrgsolve_param)))
  idata$ID <- seq_len(n)

  cmt <- if (identical(regimen$route, "oral")) "LUMEN" else "VEN"
  ev_df <- do.call(rbind, lapply(seq_len(n), function(s) {
    data.frame(
      ID = s, time = reg$times, amt = amt[s], cmt = cmt, evid = 1,
      rate = if (reg$inf) amt[s] / reg$inf_dur else 0
    )
  }))

  out <- mrgsolve::mrgsim_df(
    mrgsolve::update(mod, rtol = rtol, atol = atol, maxsteps = 1e6),
    data = ev_df, idata = idata, add = sort(unique(times)), end = -1,
    obsonly = TRUE, carry_out = character(0), recover = character(0)
  )
  plasma <- matrix(out$CP, nrow = length(times), ncol = n)
  list(time = times, plasma = plasma, states = NULL)
}
