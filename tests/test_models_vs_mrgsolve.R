# The DDI and antibody models, solved by the browser engine
# (shared/ode_engine.R, at the apps' tolerances) and by mrgsolve compiling
# the same model file, must agree. (The small-molecule model has its own
# exact engine and its own test: tests/test_small_molecule_vs_mrgsolve.R.)
#
# Run from the repository root: Rscript tests/test_models_vs_mrgsolve.R

suppressMessages(library(shiny))
for (f in c("ode_engine.R", "app_helpers.R", "rifampicin.R")) source(file.path("shared", f))
source(file.path("reference", "mrgsolve_engine.R"))
if (!mrgsolve_ready()) stop("mrgsolve is not installed")

rel <- function(a, b) max(abs(a - b)) / max(abs(b))
worst <- c()
report <- function(name, d) { cat(sprintf("  %-42s %.2e\n", name, d)); worst[name] <<- d }

# Rifampicin 600 mg x 7 days, then midazolam 3 mg PO and 1 mg IV
m <- mrg_read("models/ddi_rifampicin_midazolam.cpp")
t_v <- 7 * 24
ev <- rbind(rif_events(600, 7), data.frame(time = t_v + MDZ_TLAG, cmt = "MDZ_GUT", amt = 3, rate = 0))
tt <- ddi_times(t_v, 24)
a <- mrg_solve(m, data.frame(BW = 70), ev, tt, rtol = 1e-5, atol = 1e-10)
b <- mrg_solve_mrgsolve(m, data.frame(BW = 70), ev, tt)
for (v in c("MDZ", "RIF", "CYP3A_LIVER", "CYP3A_GUT")) report(paste("midazolam DDI", v), rel(a[[v]], b[[v]]))

# Rifampicin IV + glibenclamide, and 7 days oral rifampicin then glibenclamide on day 9
m <- mrg_read("models/ddi_rifampicin_glibenclamide.cpp")
for (sc in c("iv", "po")) {
  t_v <- if (sc == "iv") 0 else 8 * 24
  ev <- rbind(if (sc == "iv") rif_events(600, 1, "iv", 0) else rif_events(600, 7),
              data.frame(time = t_v + GLB_TLAG, cmt = "GLB_GUT", amt = 1.25, rate = 0))
  tt <- ddi_times(t_v, 48)
  a <- mrg_solve(m, data.frame(BW = 70), ev, tt, rtol = 1e-5, atol = 1e-10)
  b <- mrg_solve_mrgsolve(m, data.frame(BW = 70), ev, tt)
  for (v in c("GLB", "RIF", "CYP2C9_LIVER", "OATP_ACTIVITY")) report(sprintf("glibenclamide DDI (%s) %s", sc, v), rel(a[[v]], b[[v]]))
}

# Antibody: trastuzumab 4 -> 2 mg/kg weekly (TMDD), and a linear antibody
m <- mrg_read("models/mab_mpbpk_tmdd.cpp")
nmol <- function(mg) mg * 1e6 / 148000
ev <- rbind(data.frame(time = 0, cmt = "CENT", amt = nmol(280), rate = nmol(280) / 1.5, ii = 0, addl = 0),
            data.frame(time = 168, cmt = "CENT", amt = nmol(140), rate = nmol(140) / 0.5, ii = 168, addl = 10))
tt <- seq(0, 168 * 14, by = 12) + 0.01
for (tm in 1:0) {
  a <- mrg_solve(m, data.frame(TMDD = tm), ev, tt, rtol = 1e-5, atol = 1e-9, nonneg = TRUE)
  b <- mrg_solve_mrgsolve(m, data.frame(TMDD = tm), ev, tt)
  report(sprintf("antibody plasma (TMDD=%d)", tm), rel(a$CP_UGML, b$CP_UGML))
  report(sprintf("antibody leaky ISF (TMDD=%d)", tm), rel(a$ISF_LEAKY_UGML, b$ISF_LEAKY_UGML))
}

# CAR-T PBPK-PD (mouse): four dose levels against a BCMA+ xenograft
m <- mrg_read("models/cart_pbpk_pd_singh2020.cpp")
cells <- c(1e5, 1e6, 5e6, 1e7)
ev <- do.call(rbind, lapply(seq_along(cells), function(i)
  data.frame(time = 0, cmt = "C_Blood", amt = cells[i] / 0.944, ID = i)))
P <- data.frame(VTUMOR0 = rep(0.05, length(cells)))
tt <- seq(2, 672, by = 6)   # mrgsolve reports pre-dose at time 0
a <- mrg_solve(m, P, ev, tt, rtol = 1e-6, atol = 1e-6)
b <- mrg_solve_mrgsolve(m, P, ev, tt)
for (v in c("TumorVolume", "CARTblood", "CARTtumor", "CplxPT")) report(paste("cart-pbpk", v), rel(a[[v]], b[[v]]))

# Bispecific mPBPK (Spinosa 2026): soluble targets SC, a T-cell engager weekly, cis-binding
source(file.path("apps", "bispecific", "R", "cases.R"))
m <- mrg_read("models/bispecific_mpbpk_spinosa2026.cpp")
runs <- list(
  soluble = list(case_setup(m, "soluble"), function(P) dose_events(P, 300, "mg", "SC", 1, 28), c("total_D1_cen_ugml", "TN_sR1_cen", "TN_sR2_lea")),
  tce = list(case_setup(m, "tce"), function(P) dose_events(P, 0.1, "mg/kg", "IV", 3, 7), c("D1_cen_ugml", "RO_mR1_cen", "RO_mR2_cen")),
  cis = list(case_setup(m, "cis", chi_e = 1000), function(P) dose_events(P, 0.3, "mg/kg", "IV", 1, 7), c("D1_eff_ugml", "RO_mR1_cen", "RO_mR1_eff"))
)
tt <- seq(0.5, 28, by = 0.5) + 0.01   # off the dose times: mrgsolve reports pre-dose there
for (nm in names(runs)) {
  P <- runs[[nm]][[1]]; ev <- runs[[nm]][[2]](P)
  a <- mrg_solve(m, P, ev, tt, rtol = 1e-4, atol = case_tol(nm)$atol, nonneg = TRUE)
  b <- mrg_solve_mrgsolve(m, P, ev, tt)
  for (v in runs[[nm]][[3]]) report(paste("bispecific", nm, v), rel(a[[v]], b[[v]]))
}

cat(sprintf("comparisons: %d, worst: %.2e\n", length(worst), max(worst)))
if (!all(is.finite(worst)) || max(worst) > 2e-3) stop("browser engine disagrees with mrgsolve")
cat("PASS\n")
