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

# --- Bispecific antibodies, generalized minimal PBPK (Spinosa 2026) --------------------
source(file.path("apps", "bispecific", "R", "cases.R"))
m <- mrg_read("models/bispecific_mpbpk_spinosa2026.cpp")
bs <- function(P, ev, tt, case = "tce") mrg_solve(m, P, ev, tt, rtol = 1e-4, atol = case_tol(case)$atol, nonneg = TRUE)
# Case 1: anti-IL-13/IL-17 (BITS7201A), 750 mg IV (Fig. 3d-e: ~80% and ~90% on day 28)
P <- case_setup(m, "soluble")
r <- bs(P, dose_events(P, 750, "mg", "IV", 1, 28), c(0.01, 28), "soluble")
il17 <- 100 - r$TN_sR2_cen[2, 1]; il13 <- 100 - r$TN_sR1_cen[2, 1]
check(il17 > 70 && il17 < 90 && il13 > 85 && il13 < 99,
      sprintf("BITS7201A 750 mg IV neutralises IL-17AA by %.0f%% and IL-13 by %.0f%% on day 28 (paper ~80%% and ~90%%)", il17, il13))
# Case 2: mosunetuzumab in monkeys; T cells clear much of the drug even at 1 mg/kg (Fig. 5e-f)
P <- case_setup(m, "tce"); tt <- seq(0.1, 21, by = 0.1)
auc <- function(P) { r <- bs(P, dose_events(P, 1, "mg/kg", "IV", 1, 7), tt); auc_trap(tt, r$D1_cen_ugml[, 1]) }
a0 <- auc(P); P3 <- P; P3$rpc_mR2_cen <- 0; P20 <- P; P20$rpc_mR1_cen <- 0
a3 <- auc(P3) / a0; a20 <- auc(P20) / a0
check(a3 > 1.4 && a20 < 1.25 && a3 > a20,
      sprintf("mosunetuzumab 1 mg/kg: exposure x%.1f without CD3 binding, only x%.2f without CD20", a3, a20))
# Case 3: CD3 affinity decides how much drug T cells clear (Fig. 6b-c)
share <- function(kd) {
  P <- case_setup(m, "affinity", cd3_kd = kd)
  s <- bs(P, dose_events(P, 0.1, "mg/kg", "IV", 1, 7), c(0.01, 21))$states[2, , 1]
  s[c("CL_D1_NS_nmol", "CL_D1_TMDD_R2_cen_nmol")] / sum(s[c("CL_D1_NS_nmol", "CL_D1_TMDD_R1_cen_nmol", "CL_D1_TMDD_R2_cen_nmol")])
}
s40 <- share(40); s400 <- share(400)
check(s40[2] > 0.5 && s400[1] > 0.5,
      sprintf("0.1 mg/kg: %.0f%% cleared through CD3 at 40 nM; %.0f%% nonspecific at 400 nM", 100 * s40[2], 100 * s400[1]))
# Case 4: tumour-targeted bispecific; avidity sets tumour but not blood occupancy (Fig. 7)
ro <- function(chi, fm, dose) {
  P <- case_setup(m, "cis", chi_e = chi, format = fm)
  r <- bs(P, dose_events(P, dose, "mg/kg", "IV", 1, 7), c(0.01, 21))
  c(blood = r$RO_mR1_cen[2, 1], tumour = r$RO_mR1_eff[2, 1])
}
hi <- ro(1000, "bispecific", 0.1); lo <- ro(1, "bispecific", 0.1); bv <- ro(1000, "bivalent", 0.1)
check(hi["tumour"] > 50 && hi["blood"] < 40 && lo["tumour"] < lo["blood"] + 5 && abs(hi["blood"] - lo["blood"]) < 5,
      sprintf("0.1 mg/kg: avidity 1000 gives %.0f%% in tumour vs %.0f%% in blood; without avidity %.0f%% vs %.0f%%",
              hi["tumour"], hi["blood"], lo["tumour"], lo["blood"]))
check(bv["blood"] > bv["tumour"], sprintf("a bivalent 10 nM antibody engages blood (%.0f%%) more than tumour (%.0f%%)", bv["blood"], bv["tumour"]))

cat("PASS\n")
