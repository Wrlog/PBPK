# ============================================================================
# Dosing helpers for the rifampicin DDI models (Asaumi et al. 2018).
#
# Oral doses go into the gut lumen after the lag time estimated in the
# paper; IV doses are infused into the blood concentration state, so the
# amount is divided by the blood volume. Times are in hours.
# ============================================================================

RIF_TLAG <- 0.255        # h, rifampicin oral lag time (Table 1)
MDZ_TLAG <- 0.031        # h, midazolam oral lag time (Table 1)
GLB_TLAG <- 0.773        # h, glibenclamide oral lag time (Table 1, beta 0.2)
VB_PER_KG <- 0.0743      # L/kg, blood volume (Table S1)
MDZ_VC_PER_KG <- 0.571   # L/kg, midazolam central volume (Table 1)

#' Rifampicin dose records
#'
#' @param route "po" (once daily for `days`) or "iv" (single infusion)
#' @param start time of the first dose (h)
rif_events <- function(dose_mg, days = 1, route = "po", start = 0, bw = 70, inf_h = 0.5) {
  if (dose_mg <= 0 || days <= 0) return(NULL)
  if (route == "iv") {
    amt <- dose_mg / (VB_PER_KG * bw)
    return(data.frame(time = start, cmt = "RIF_CENT", amt = amt, rate = amt / inf_h))
  }
  data.frame(time = start + (seq_len(days) - 1) * 24 + RIF_TLAG, cmt = "RIF_LUMEN", amt = dose_mg, rate = 0)
}

auc_trap <- function(t, y) sum(diff(t) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

#' Dense output around a victim dose, sparse elsewhere
ddi_times <- function(t_victim, horizon = 24, before = 0) {
  sort(unique(c(seq(0, max(t_victim, 0), by = 2),
                t_victim + c(seq(0, 2, by = 0.05), seq(2.25, horizon, by = 0.25)))))
}
