# ============================================================================
# The case studies of Spinosa et al. 2026, as parameter sets for the model.
#
# data/spinosa2026_casestudies.csv is the paper's "Model parameters.xlsx"
# (Data S1), one sheet per case study, in the gQSPSim naming of the authors'
# project: a _D1 suffix marks the drug's parameters. case_params() maps those
# names onto the model file's.
# ============================================================================

.case_env <- new.env()

# read on first use: Shiny sources R/ before app.R, when the shared helpers
# are not loaded yet
case_table <- function() {
  if (is.null(.case_env$table)) {
    f <- file.path("data", "spinosa2026_casestudies.csv")
    if (!file.exists(f)) f <- file.path("apps", "bispecific", "data", "spinosa2026_casestudies.csv")   # tests, from the root
    .case_env$table <- utils::read.csv(f, stringsAsFactors = FALSE)
  }
  .case_env$table
}

case_params <- function(model, sheet) {
  tab <- case_table()
  x <- tab[tab$case == sheet, ]
  nm <- sub("_D1$", "", x$parameter)
  nm[nm == "F"] <- "F_sc"
  keep <- nm %in% names(model$param) & !grepl("_D2$", x$parameter)
  as.data.frame(stats::setNames(as.list(x$value[keep]), nm[keep]))
}

# Effective avidity of the second arm (Data S1): chi_e = chi / (Vchi * cells/mL).
# The apps set chi_e directly, so the avidity parameter is chi_e * Vchi * cells.
avidity_param <- function(model, chi_e, cells_per_ml) chi_e * model$param[["Vchi"]] * cells_per_ml

CASES <- list(
  soluble = list(
    label = "Soluble targets: anti-IL-13/IL-17 (BITS7201A)",
    sheet = 1, kind = "soluble", species = "Healthy volunteers, 70 kg",
    R1 = "IL-13", R2 = "IL-17AA",
    note = "Case study 1. Both targets are soluble cytokines made and cleared in blood, leaky and tight tissue; the bispecific holds them in slowly cleared complexes, so total cytokine rises while free cytokine is neutralised. Calibrated to the single-ascending-dose trial (Staton et al. 2019)."
  ),
  tce = list(
    label = "T-cell engager: mosunetuzumab (CD20 x CD3)",
    sheet = 2, kind = "trans", species = "Cynomolgus monkey, 3.12 kg",
    R1 = "CD20 (B cells)", R2 = "CD3 (T cells)",
    note = "Case study 2. The drug bridges CD20 on B cells and CD3 on T cells (trans-binding); both targets sit in blood. CD3 and CD20 capacities were estimated from anti-gD/CD3 and obinutuzumab PK, then mosunetuzumab's nonspecific clearance was fitted."
  ),
  affinity = list(
    label = "T-cell engager: CD3 affinity",
    sheet = 3, kind = "trans", species = "Cynomolgus monkey, 3.12 kg",
    R1 = "CD20 (B cells)", R2 = "CD3 (T cells)",
    note = "Case study 3. The same system with the fitted CD3 (19 nM) and CD20 (4.4 nM) capacities, comparing a moderate (40 nM) with a low (400 nM) CD3 affinity: how much of the drug is cleared through T cells."
  ),
  cis = list(
    label = "Tumour-targeted bispecific (cis) vs bivalent antibody",
    sheet = 4, kind = "cis", species = "Patients, 70 kg",
    R1 = "therapeutic receptor", R2 = "targeting receptor",
    note = "Case study 4. Cells carry both receptors everywhere; the targeting receptor is 100-fold higher on tumour cells (the efficacy space). One arm binds the targeting receptor tightly (0.05 nM), the other the therapeutic receptor weakly (100 nM), so avidity steers binding to cells that carry both."
  )
)

#' Parameters for a case with the app's choices applied
case_setup <- function(model, case, cd3_kd = 40, chi_e = 1000, format = "bispecific") {
  cs <- CASES[[case]]
  P <- case_params(model, cs$sheet)
  if (case == "affinity") {
    if (is.na(cd3_kd)) P$rpc_mR2_cen <- 0 else P$Kd_R2 <- cd3_kd
  }
  if (case == "cis") {
    cells <- P$mR1_cellspermL_cen
    P$FLAG_Avidity <- 1
    P$avid_R1 <- avidity_param(model, chi_e, cells)
    P$avid_R2 <- avidity_param(model, 1, cells)          # no avidity gain for the targeting arm
    if (format == "bivalent") {
      # Table S4: a bivalent, monospecific antibody against the therapeutic
      # receptor only (10 nM), same clearance and the same avidity.
      P$FLAG_CIS_BINDING <- 0
      P$FLAG_BIVALENT_BINDING <- 1
      P$Kd_R1 <- 10
      P$kon_R2 <- 0
    }
  }
  P
}

#' Solver tolerances: the soluble cytokines sit near 1e-5 nM, so their case
#' needs a much smaller absolute tolerance than the receptor cases (nM)
case_tol <- function(case) list(rtol = 1e-4, atol = if (case == "soluble") 1e-12 else 1e-8)

#' Doses in the model's units: IV bolus in nM (D1_cen) or SC in mg (D1_ext_mg)
dose_events <- function(P, amount, unit, route, n, every_days, ID = NULL) {
  mg <- if (unit == "mg/kg") amount * P$BW else amount
  ev <- if (route == "SC") {
    data.frame(time = 0, cmt = "D1_ext_mg", amt = mg, ii = every_days, addl = n - 1)
  } else {
    data.frame(time = 0, cmt = "D1_cen", amt = mg / P$BW / P$V_cen / P$MWab * 1e9, ii = every_days, addl = n - 1)
  }
  if (!is.null(ID)) ev$ID <- ID
  ev
}
