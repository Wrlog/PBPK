# ============================================================================
# Virtual populations, exposure metrics and pediatric scaling.
# ============================================================================

cv_to_sd <- function(cv_percent) sqrt(log(1 + (cv_percent / 100)^2))

#' A population of one age and weight with log-normal variability
#'
#' Variability sits on the three things that differ most between people of
#' the same size: hepatic intrinsic clearance, cardiac output and GFR.
#' Partition coefficients are drug properties and are shared.
make_population <- function(drug, wt, age_y, n, cv_clint = 0, cv_co = 0,
                            cv_gfr = 0, seed = 1) {
  set.seed(seed)
  kp <- kp_rodgers_rowland(drug$logP, drug$pKa, drug$type, drug$fu, drug$BP)
  draw <- function(cv) exp(stats::rnorm(n, 0, cv_to_sd(cv)))
  m_clint <- draw(cv_clint)
  m_co <- draw(cv_co)
  m_gfr <- draw(cv_gfr)
  lapply(seq_len(n), function(i) {
    pbpk_parameters(drug, wt, age_y, kp = kp,
                    clint_mult = m_clint[i], co_mult = m_co[i], gfr_mult = m_gfr[i])
  })
}

#' Median and prediction-interval bands across subjects at each time
summarise_profiles <- function(time, conc) {
  q <- apply(conc, 1, stats::quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
             na.rm = TRUE, names = FALSE)
  if (is.null(dim(q))) q <- matrix(q, nrow = 5)
  data.frame(time = time, q05 = q[1, ], q25 = q[2, ], med = q[3, ],
             q75 = q[4, ], q95 = q[5, ])
}

trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

#' AUC by the linear-up / log-down trapezoidal rule
#'
#' Falling segments are integrated as exponentials, which is exact for a
#' mono-exponential decline and keeps the steep drop after an IV bolus from
#' being overestimated on a coarse output grid.
auc_lin_log <- function(x, y) {
  y1 <- utils::head(y, -1); y2 <- utils::tail(y, -1); dx <- diff(x)
  down <- y2 < y1 & y2 > 0 & y1 > 0
  seg <- dx * (y1 + y2) / 2
  seg[down] <- dx[down] * (y1[down] - y2[down]) / log(y1[down] / y2[down])
  sum(seg)
}

#' Exposure over the final dosing interval, per subject
#'
#' The window runs from the last dose to one interval later (or the end of
#' the simulation, whichever comes first).
exposure_metrics <- function(time, conc, regimen) {
  t_end <- max(time)
  doses <- expand_regimen(regimen, t_end)$times
  start <- max(doses)
  end <- min(start + regimen$interval, t_end)
  if (regimen$n_doses == 1) end <- t_end
  w <- time >= start - 1e-9 & time <= end + 1e-9
  tw <- time[w]
  out <- t(apply(conc[w, , drop = FALSE], 2, function(cw) {
    c(cmax = max(cw), tmax = tw[which.max(cw)] - start,
      cmin = cw[length(cw)], auc = auc_lin_log(tw, cw))
  }))
  list(window = c(start, end), per_subject = as.data.frame(out))
}

#' Typical subjects across childhood, for the Age scaling tab
#'
#' Weight follows the growth-chart median for each age. Everything reported
#' is exact (from pbpk_metrics), not read off a simulated curve.
#'
#' @param ref_dose_mgkg adult reference dose in mg/kg
#' @param oral whether AUC should include oral bioavailability
age_scaling_table <- function(drug, ages, ref_dose_mgkg, oral, adult_age = 30, adult_wt = 70) {
  kp <- kp_rodgers_rowland(drug$logP, drug$pKa, drug$type, drug$fu, drug$BP)
  one <- function(age, wt) {
    m <- pbpk_metrics(pbpk_parameters(drug, wt, age, kp = kp))
    f <- if (oral) m$F_oral else 1
    data.frame(age = age, wt = wt, CL = m$CL, CL_kg = m$CL / wt,
               Vss_kg = m$Vss / wt, t_half = m$t_half, F = f,
               auc_per_mgkg = f * wt / m$CL,
               enzyme = enzyme_maturation(age, drug$pathway),
               gfr_mat = gfr_maturation(age))
  }
  adult <- one(adult_age, adult_wt)
  tab <- do.call(rbind, lapply(ages, function(a) one(a, typical_weight(a))))
  tab$auc <- tab$auc_per_mgkg * ref_dose_mgkg
  # Linear kinetics: the dose that matches the adult AUC scales with the
  # ratio of AUC per mg/kg.
  tab$dose_match <- ref_dose_mgkg * adult$auc_per_mgkg / tab$auc_per_mgkg
  list(table = tab, adult = adult)
}

AGE_LABELS <- c("1 week" = 7 / 365, "1 month" = 1 / 12, "3 months" = 0.25,
                "1 year" = 1, "2 years" = 2, "6 years" = 6, "12 years" = 12,
                "Adult" = 30)
