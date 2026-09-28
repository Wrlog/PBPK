# ============================================================================
# Whole-body PBPK engine: exact solution by matrix exponential.
#
# The model (models/pbpk.cpp) is perfusion-limited with linear hepatic and
# renal clearance and first-order oral absorption. That makes it a linear,
# time-invariant ODE system
#
#     dx/dt = A x + b u(t)
#
# where u is the infusion rate, which is piecewise constant. Over any step of
# length h with constant u the solution is exact:
#
#     [x; u](t + h) = expm(M h) [x; u](t),   M = [A b; 0 0]
#
# so the engine never integrates anything numerically. It needs no compiled
# code, which is what lets the app run in a browser under webR. The same
# model is run through mrgsolve by reference/mrgsolve_engine.R, and
# tests/test_engine_vs_mrgsolve.R fails the build if the two disagree.
#
# Units: amount mg, volume L, flow and clearance L/h, time h, so
# concentrations are mg/L.
# ============================================================================

STATES <- c("LUMEN", "LUNG", "BRAIN", "HEART", "KIDNEY", "MUSCLE", "SKIN",
            "ADIPOSE", "BONE", "SPLEEN", "GUT", "LIVER", "REST", "ART", "VEN")
N_STATE <- length(STATES)
U_IDX <- N_STATE + 1          # augmented state carrying the infusion rate

# Tissues whose venous outflow goes straight to the vena cava.
TO_VENOUS <- c("brain", "heart", "kidney", "muscle", "skin", "adipose", "bone", "rest")

#' Everything the ODEs need for one subject
#'
#' @param drug a DRUGS entry (or an edited copy): logP, pKa, type, fu, BP,
#'   clint, pathway, renal_factor, ka, fa, fg
#' @param clint_mult,co_mult,gfr_mult individual multipliers (variability)
pbpk_parameters <- function(drug, wt, age_y, kp = NULL,
                            clint_mult = 1, co_mult = 1, gfr_mult = 1) {
  phys <- physiology(wt, age_y)
  if (co_mult != 1) {
    phys$Q <- phys$Q * co_mult
    phys$CO <- phys$CO * co_mult
  }
  if (is.null(kp)) {
    kp <- kp_rodgers_rowland(drug$logP, drug$pKa, drug$type, drug$fu, drug$BP)
  }

  # Intrinsic clearance scales with liver mass and with enzyme maturation.
  liver_ratio <- phys$V[["liver"]] / (VOL_FRAC[["liver"]] * REF_WT)
  clint <- drug$clint * liver_ratio * enzyme_maturation(age_y, drug$pathway) * clint_mult

  list(
    wt = wt, age = age_y,
    V = phys$V, Q = phys$Q, CO = phys$CO, Kp = kp,
    BP = drug$BP, fu = drug$fu,
    CLint = clint,
    CLr = phys$GFR * gfr_mult * drug$fu * drug$renal_factor,
    GFR = phys$GFR * gfr_mult,
    ka = drug$ka, fa = drug$fa, fg = drug$fg
  )
}

#' Rate matrix A (15 x 15, amounts) for one subject
pbpk_matrix <- function(p) {
  A <- matrix(0, N_STATE, N_STATE, dimnames = list(STATES, STATES))
  V <- p$V; Q <- p$Q; Kp <- p$Kp; BP <- p$BP
  S <- function(t) toupper(t)

  # First-order rate at which amount leaves tissue t in venous blood:
  # Q * C_out, with C_out = C_tissue * BP / Kp.
  kout <- function(t) Q[[t]] * BP / (Kp[[t]] * V[[t]])

  # Arterial supply of tissue t.
  supply <- function(t, q) {
    A[S(t), "ART"] <<- A[S(t), "ART"] + q / V[["art"]]
  }

  for (t in TO_VENOUS) {
    supply(t, Q[[t]])
    A[S(t), S(t)] <- A[S(t), S(t)] - kout(t)
    A["VEN", S(t)] <- A["VEN", S(t)] + kout(t)
  }

  # Splanchnic tissues drain to the liver through the portal vein.
  for (t in c("spleen", "gut")) {
    supply(t, Q[[t]])
    A[S(t), S(t)] <- A[S(t), S(t)] - kout(t)
    A["LIVER", S(t)] <- A["LIVER", S(t)] + kout(t)
  }

  # Liver: hepatic artery in, total liver flow out, and metabolism of the
  # unbound drug (C_liver * fu / Kp_liver) at the intrinsic clearance.
  supply("liver", Q[["hepatic_artery"]])
  A["LIVER", "LIVER"] <- A["LIVER", "LIVER"] - kout("liver") -
    p$CLint * p$fu / (Kp[["liver"]] * V[["liver"]])
  A["VEN", "LIVER"] <- A["VEN", "LIVER"] + kout("liver")

  # Kidney: filtration of the plasma leaving the kidney (C_kidney / Kp).
  A["KIDNEY", "KIDNEY"] <- A["KIDNEY", "KIDNEY"] - p$CLr / (Kp[["kidney"]] * V[["kidney"]])

  # Lung takes the whole cardiac output from the venous pool.
  A["LUNG", "VEN"] <- p$CO / V[["ven"]]
  A["LUNG", "LUNG"] <- -kout("lung")
  A["ART", "LUNG"] <- kout("lung")
  A["ART", "ART"] <- -p$CO / V[["art"]]
  A["VEN", "VEN"] <- -p$CO / V[["ven"]]

  # Oral absorption: the lumen empties at ka; a fraction fa * fg of what
  # leaves reaches the gut wall intact, the rest is lost.
  A["LUMEN", "LUMEN"] <- -p$ka
  A["GUT", "LUMEN"] <- p$ka * p$fa * p$fg

  A
}

#' Matrix exponential by scaling and squaring with a [6/6] Pade approximant
expm_pade <- function(M) {
  n <- nrow(M)
  nrm <- max(rowSums(abs(M)))
  s <- if (nrm > 0.5) ceiling(log2(nrm / 0.5)) else 0
  M <- M / 2^s

  q <- 6
  cf <- numeric(q)
  cf[1] <- 0.5
  for (k in 2:q) cf[k] <- cf[k - 1] * (q - k + 1) / (k * (2 * q - k + 1))

  I <- diag(n)
  X <- M
  N <- I + cf[1] * M
  D <- I - cf[1] * M
  for (k in 2:q) {
    X <- M %*% X
    N <- N + cf[k] * X
    D <- D + (-1)^k * cf[k] * X
  }
  E <- solve(D, N)
  for (i in seq_len(s)) E <- E %*% E
  E
}

#' Dose times, amounts per subject and the infusion schedule of a regimen
#'
#' @param regimen list(route = "oral" | "iv_bolus" | "iv_inf", amt (mg,
#'   scalar or one per subject), interval (h), n_doses, inf_dur (h))
expand_regimen <- function(regimen, t_end) {
  times <- (seq_len(regimen$n_doses) - 1) * regimen$interval
  times <- times[times <= t_end + 1e-9]
  inf <- identical(regimen$route, "iv_inf") && regimen$inf_dur > 0
  list(times = times, inf = inf, inf_dur = if (inf) regimen$inf_dur else 0)
}

#' Simulate one or more subjects
#'
#' @param plist list of parameter sets from pbpk_parameters()
#' @param times output times (h)
#' @param keep_states also return every compartment's concentration
#' @return list(time, plasma = matrix [time x subject], states = array
#'   [time x state x subject] when keep_states)
pbpk_simulate <- function(plist, regimen, times, keep_states = FALSE) {
  n <- length(plist)
  times <- sort(unique(times))
  t_end <- max(times)
  reg <- expand_regimen(regimen, t_end)
  amt <- rep_len(regimen$amt, n)

  # Every point where something changes, plus every output time. The state
  # is propagated exactly between consecutive points.
  grid <- sort(unique(round(c(0, times, reg$times,
                              if (reg$inf) reg$times + reg$inf_dur), 10)))
  grid <- grid[grid <= t_end + 1e-9]
  h <- diff(grid)
  h_key <- round(h, 9)
  h_unique <- unique(h_key)

  m <- U_IDX
  # One propagator per subject per distinct step length, stored as columns
  # of (m*m) x n so a step for every subject is a handful of vector ops.
  E <- lapply(h_unique, function(hh) matrix(0, m * m, n))
  v_ven <- numeric(n); bp <- numeric(n); vols <- matrix(0, N_STATE, n)
  for (s in seq_len(n)) {
    p <- plist[[s]]
    A <- pbpk_matrix(p)
    M <- matrix(0, m, m)
    M[1:N_STATE, 1:N_STATE] <- A
    M[which(STATES == "VEN"), U_IDX] <- 1
    for (k in seq_along(h_unique)) E[[k]][, s] <- as.vector(expm_pade(M * h_unique[k]))
    v_ven[s] <- p$V[["ven"]]
    bp[s] <- p$BP
    vols[, s] <- c(1, p$V[c("lung", "brain", "heart", "kidney", "muscle", "skin",
                            "adipose", "bone", "spleen", "gut", "liver", "rest",
                            "art", "ven")])
  }

  dose_cmt <- if (identical(regimen$route, "oral")) 1 else which(STATES == "VEN")
  out_idx <- match(round(times, 10), grid)
  plasma <- matrix(NA_real_, length(times), n)
  states <- if (keep_states) array(NA_real_, c(length(times), N_STATE, n)) else NULL

  X <- matrix(0, m, n)
  record <- function(j) {
    o <- which(out_idx == j)
    if (length(o)) {
      plasma[o, ] <<- X[which(STATES == "VEN"), ] / (v_ven * bp)
      if (keep_states) states[o, , ] <<- X[1:N_STATE, , drop = FALSE] / vols
    }
  }

  for (j in seq_along(grid)) {
    tj <- grid[j]
    # Bolus doses (oral or IV) land at their dose time before recording.
    if (!reg$inf) {
      hit <- sum(abs(reg$times - tj) < 1e-9)
      if (hit) X[dose_cmt, ] <- X[dose_cmt, ] + hit * amt
    }
    record(j)
    if (j == length(grid)) break

    rate <- if (reg$inf) {
      sum(reg$times <= tj + 1e-9 & tj < reg$times + reg$inf_dur - 1e-9)
    } else 0
    X[U_IDX, ] <- rate * amt / max(reg$inf_dur, 1e-12)

    Ek <- E[[match(h_key[j], h_unique)]]
    Xn <- matrix(0, m, n)
    for (k in seq_len(m)) {
      rows <- ((k - 1) * m + 1):(k * m)
      Xn <- Xn + Ek[rows, , drop = FALSE] * rep(X[k, ], each = m)
    }
    X <- Xn
  }

  if (keep_states) dimnames(states) <- list(NULL, STATES, NULL)
  list(time = times, plasma = plasma, states = states)
}

#' Exact secondary parameters from the rate matrix
#'
#' Clearance and volume come from the moments of the unit IV-bolus response
#' (AUC = c' (-A)^-1 b, AUMC = c' A^-2 b), the terminal half-life from the
#' slowest eigenvalue of the disposition system, and oral bioavailability
#' from the ratio of oral to IV AUC. None of them depend on a time grid.
pbpk_metrics <- function(p) {
  A <- pbpk_matrix(p)
  ven <- which(STATES == "VEN")
  c_pl <- numeric(N_STATE); c_pl[ven] <- 1 / (p$V[["ven"]] * p$BP)
  b_iv <- numeric(N_STATE); b_iv[ven] <- 1
  b_po <- numeric(N_STATE); b_po[1] <- 1

  y_iv <- solve(-A, b_iv)
  auc_iv <- sum(c_pl * y_iv)
  aumc_iv <- sum(c_pl * solve(-A, y_iv))
  auc_po <- sum(c_pl * solve(-A, b_po))

  cl <- 1 / auc_iv
  disp <- A[-1, -1]
  lambda_z <- min(abs(Re(eigen(disp, only.values = TRUE)$values)))

  fub <- p$fu / p$BP
  q_h <- p$Q[["liver"]]
  e_h <- fub * p$CLint / (q_h + fub * p$CLint)

  list(
    CL = cl,
    Vss = cl * aumc_iv / auc_iv,
    t_half = log(2) / lambda_z,
    F_oral = auc_po / auc_iv,
    E_H = e_h,
    CLint = p$CLint,
    CLr = p$CLr,
    GFR = p$GFR
  )
}
