# The DDI and antibody models must reproduce what their papers report.
# Uses the browser engine only.
#
# Run from the repository root: Rscript tests/test_published_behaviour.R

suppressMessages(library(shiny))
for (f in c("ode_engine.R", "app_helpers.R", "rifampicin.R")) source(file.path("shared", f))

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}

# --- Rifampicin -> midazolam (Asaumi 2018, Figure 5) -------------------------
m <- mrg_read("models/ddi_rifampicin_midazolam.cpp")
mdz_auc <- function(rif_mg, route) {
  t_v <- if (rif_mg > 0) 7 * 24 else 0
  victim <- if (route == "po") data.frame(time = t_v + MDZ_TLAG, cmt = "MDZ_GUT", amt = 3, rate = 0)
            else data.frame(time = t_v, cmt = "MDZ_CENT", amt = 1 / (MDZ_VC_PER_KG * 70), rate = 0)
  tt <- ddi_times(t_v, 24)
  r <- mrg_solve(m, data.frame(BW = 70), rbind(rif_events(rif_mg, 7), victim), tt, rtol = 1e-5, atol = 1e-10)
  w <- tt >= t_v
  auc_trap(tt[w], r$MDZ[w, 1])
}
po0 <- mdz_auc(0, "po"); iv0 <- mdz_auc(0, "iv")
r600 <- mdz_auc(600, "po") / po0
r75 <- mdz_auc(75, "po") / po0
riv <- mdz_auc(600, "iv") / iv0
check(r600 > 0.04 && r600 < 0.15, sprintf("rifampicin 600 mg cuts oral midazolam AUC by ~90%% (ratio %.3f)", r600))
check(r75 > r600 && r75 < 0.6, sprintf("the induction is dose-dependent (75 mg: %.2f)", r75))
check(riv > 0.35 && riv < 0.6, sprintf("IV midazolam AUC ratio ~0.45 (hepatic induction only; %.2f)", riv))

# --- Rifampicin -> glibenclamide (Asaumi 2018, Figure 6) ------------------------
m <- mrg_read("models/ddi_rifampicin_glibenclamide.cpp")
glb_auc <- function(scenario) {
  t_v <- switch(scenario, ctrl = 0, iv = 0, po = 8 * 24)
  ev <- rbind(switch(scenario, ctrl = NULL, iv = rif_events(600, 1, "iv", 0), po = rif_events(600, 7)),
              data.frame(time = t_v + GLB_TLAG, cmt = "GLB_GUT", amt = 1.25, rate = 0))
  tt <- ddi_times(t_v, 48)
  r <- mrg_solve(m, data.frame(BW = 70), ev, tt, rtol = 1e-5, atol = 1e-10)
  w <- tt >= t_v
  auc_trap(tt[w], r$GLB[w, 1])
}
ctrl <- glb_auc("ctrl")
iv <- glb_auc("iv") / ctrl
po <- glb_auc("po") / ctrl
check(abs(iv - 2.08) < 0.25, sprintf("single IV rifampicin: glibenclamide AUC ratio ~2.08 as published (%.2f)", iv))
check(abs(po - 0.50) < 0.1, sprintf("after 7 days oral rifampicin: AUC ratio ~0.50 as published (%.2f)", po))

# --- Antibody mPBPK (Cao 2013/2014) ------------------------------------------------
m <- mrg_read("models/mab_mpbpk_tmdd.cpp")
nmol <- function(mg) mg * 1e6 / 148000
ev <- rbind(data.frame(time = 0, cmt = "CENT", amt = nmol(280), rate = nmol(280) / 1.5, ii = 0, addl = 0),
            data.frame(time = 168, cmt = "CENT", amt = nmol(140), rate = nmol(140) / 0.5, ii = 168, addl = 10))
tt <- seq(0, 168 * 12, by = 12)
r <- mrg_solve(m, data.frame(row.names = 1), ev, tt, rtol = 1e-5, atol = 1e-9, nonneg = TRUE)
trough <- r$CP_UGML[which(tt == 168 * 11), 1]
check(trough > 35 && trough < 90, sprintf("trastuzumab 4/2 mg/kg weekly: week-12 trough 35-90 ug/mL (%.0f)", trough))
auc_per_dose <- function(mgkg) {
  x <- mrg_solve(m, data.frame(row.names = 1), data.frame(time = 0, cmt = "CENT", amt = nmol(mgkg * 70)),
                 seq(0, 168 * 12, by = 12), rtol = 1e-5, atol = 1e-10, nonneg = TRUE)
  auc_trap(x$time, x$CP_UGML[, 1]) / mgkg
}
check(auc_per_dose(0.3) < 0.5 * auc_per_dose(10), "target-mediated clearance: low doses are cleared much faster")
lin <- mrg_solve(m, data.frame(TMDD = 0), data.frame(time = 0, cmt = "CENT", amt = nmol(700)), tt, rtol = 1e-5, atol = 1e-10)
check(max(lin$ISF_LEAKY_UGML) > max(lin$ISF_TIGHT_UGML), "leaky tissues reach higher interstitial concentrations than tight ones")

# --- CAR-T PBPK-PD, mouse (Singh 2020) ---------------------------------------------
m <- mrg_read("models/cart_pbpk_pd_singh2020.cpp")
cart <- function(cells, days = 28, P = data.frame(VTUMOR0 = 0.05)) {
  ev <- if (cells > 0) data.frame(time = 0, cmt = "C_Blood", amt = cells / 0.944) else NULL
  tt <- seq(0, 24 * days, by = 6)
  list(t = tt / 24, r = mrg_solve(m, P, ev, tt, rtol = 1e-6, atol = 1e-6))
}
ctrl <- cart(0); hi <- cart(1e7); lo <- cart(1e5)
n <- length(ctrl$t)
check(ctrl$r$TumorVolume[n, 1] / ctrl$r$TumorVolume[1, 1] > 5,
      sprintf("the untreated xenograft grows from 50 to %.0f mm3 in 28 days", ctrl$r$TumorVolume[n, 1]))
check(hi$r$TumorVolume[n, 1] < 0.05 * ctrl$r$TumorVolume[n, 1],
      sprintf("10 million CAR-T cells clear the tumour (%.1f mm3 vs %.0f untreated)", hi$r$TumorVolume[n, 1], ctrl$r$TumorVolume[n, 1]))
check(lo$r$TumorVolume[n, 1] > 0.5 * ctrl$r$TumorVolume[n, 1], "0.1 million cells do not control it: the response is dose-dependent")
i <- which.max(hi$r$CARTtumor[, 1])
check(hi$t[i] > 3 && max(hi$r$CARTtumor[, 1]) > 10 * hi$r$CARTtumor[which.min(abs(hi$t - 1)), 1],
      sprintf("CAR-T cells expand in the tumour, peaking on day %.0f", hi$t[i]))
lung <- max(hi$r$states[, "C_E_Lung", 1]); brain <- max(hi$r$states[, "C_E_Brain", 1])
check(lung > 100 * brain, "cells accumulate in lung far more than in brain, as the fitted transmigration rates imply")
# affinity saturates: 0.1 nM is no better than 1 nM, but 1000 nM is worse
kd_run <- function(kd) {
  P <- data.frame(KOFF = 1.08e-12 * kd * 6.022e11, VTUMOR0 = 0.05)
  cart(1e7, P = P)$r$TumorVolume[n, 1]
}
v <- vapply(c(0.1, 1, 1000), kd_run, 0)
check(abs(v[1] - v[2]) < 0.01 * ctrl$r$TumorVolume[n, 1] && v[3] > 10 * v[2],
      "affinity saturates: 0.1 and 1 nM give the same result, 1000 nM much less killing")

cat("PASS\n")
