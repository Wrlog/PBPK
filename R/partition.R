# ============================================================================
# Tissue:plasma partition coefficients by the Rodgers & Rowland method.
#
#   Rodgers T, Leahy D, Rowland M. J Pharm Sci 2005;94:1259  (bases, pKa >= 7)
#   Rodgers T, Rowland M.        J Pharm Sci 2006;95:1238  (acids, neutrals,
#                                                            weak bases)
#
# A tissue is modelled as extracellular water, intracellular water, neutral
# lipids and neutral phospholipids, plus binding to acidic phospholipids
# (moderate-to-strong bases) or to albumin / lipoproteins (everything else).
# The ionisation terms use the pH of plasma, intracellular water and red
# cells.
#
# Tissue composition is the rat data set tabulated by Rodgers & Rowland, as
# distributed with Utsey et al. (Drug Metab Dispos 2020;48:903) in
# github.com/metrumresearchgroup/PBPK_PC. Using rat composition for human
# tissues is the usual convention for this method.
#
# One deliberate difference from that repository's script: there the
# albumin/lipoprotein binding term is multiplied by the intracellular
# ionisation factor X, which removes it entirely for neutrals. The published
# 2006 equation has no such factor, so it is omitted here.
# ============================================================================

RR_TISSUE <- data.frame(
  tissue = c("adipose", "bone", "brain", "gut", "heart", "kidney", "liver",
             "lung", "muscle", "skin", "spleen"),
  f_nl = c(0.853, 0.017, 0.039, 0.038, 0.014, 0.012, 0.014, 0.022, 0.010, 0.060, 0.0077),
  f_np = c(0.0016, 0.0017, 0.0015, 0.0125, 0.0111, 0.0242, 0.0240, 0.0128, 0.0072, 0.0044, 0.0113),
  f_ew = c(0.135, 0.100, 0.162, 0.282, 0.320, 0.273, 0.161, 0.336, 0.118, 0.382, 0.207),
  f_iw = c(0.017, 0.346, 0.620, 0.475, 0.456, 0.483, 0.573, 0.446, 0.630, 0.291, 0.579),
  ap   = c(0.40, 0.67, 0.40, 2.41, 2.25, 5.03, 4.56, 3.91, 1.53, 1.32, 3.18),
  ar   = c(0.049, 0.100, 0.048, 0.158, 0.157, 0.130, 0.086, 0.212, 0.064, 0.277, 0.097),
  lr   = c(0.068, 0.050, 0.041, 0.141, 0.160, 0.137, 0.161, 0.168, 0.059, 0.096, 0.207),
  stringsAsFactors = FALSE
)

RR_RBC    <- list(f_nl = 0.0017, f_np = 0.0029, f_iw = 0.603, ap = 0.5)
RR_PLASMA <- list(f_nl = 0.0023, f_np = 0.0013)

PH_PLASMA <- 7.4
PH_IW     <- 7.0
PH_RBC    <- 7.22

#' Rodgers & Rowland Kp values
#'
#' @param logP octanol:water partition coefficient (log10)
#' @param pKa  ignored for neutrals
#' @param type "neutral", "acid" or "base" (monoprotic)
#' @param fu   fraction unbound in plasma
#' @param BP   blood:plasma concentration ratio
#' @return named vector of tissue:plasma Kp for every entry of TISSUES.
#'   "rest" is the mean of the non-adipose tissues.
kp_rodgers_rowland <- function(logP, pKa = NA, type = c("neutral", "acid", "base"),
                               fu, BP) {
  type <- match.arg(type)

  P    <- 10^logP
  P_vo <- 10^(1.115 * logP - 1.35)   # vegetable oil:water, used for adipose

  # X: ionised/unionised ratio in intracellular water; Y: in plasma;
  # Z: in red cells. All zero for a neutral compound.
  if (type == "neutral") {
    X <- 0; Y <- 0; Z <- 0
  } else if (type == "acid") {
    X <- 10^(PH_IW - pKa); Y <- 10^(PH_PLASMA - pKa); Z <- 10^(PH_RBC - pKa)
  } else {
    X <- 10^(pKa - PH_IW); Y <- 10^(pKa - PH_PLASMA); Z <- 10^(pKa - PH_RBC)
  }

  strong_base <- type == "base" && pKa >= 7

  t <- RR_TISSUE
  p_t <- ifelse(t$tissue == "adipose", P_vo, P)
  lipid <- (p_t * t$f_nl + (0.3 * p_t + 0.7) * t$f_np) / (1 + Y)
  water <- t$f_ew + (1 + X) / (1 + Y) * t$f_iw

  if (strong_base) {
    # Affinity for acidic phospholipids, back-calculated from the blood:plasma
    # ratio through the red-cell partition coefficient.
    kpu_bc <- (HCT - 1 + BP) / (HCT * fu)
    ka_ap <- (kpu_bc - (1 + Z) / (1 + Y) * RR_RBC$f_iw -
                (P * RR_RBC$f_nl + (0.3 * P + 0.7) * RR_RBC$f_np) / (1 + Y)) *
      (1 + Y) / (RR_RBC$ap * Z)
    binding <- max(ka_ap, 0) * t$ap * X / (1 + Y)
  } else {
    # Plasma protein binding not explained by plasma lipids, carried into the
    # tissue in proportion to its albumin (acids, weak bases) or lipoprotein
    # (neutrals) content relative to plasma.
    ka_pr <- 1 / fu - 1 - (P * RR_PLASMA$f_nl + (0.3 * P + 0.7) * RR_PLASMA$f_np) / (1 + Y)
    ratio <- if (type == "neutral") t$lr else t$ar
    binding <- max(ka_pr, 0) * ratio
  }

  kp <- (water + lipid + binding) * fu
  names(kp) <- t$tissue
  kp[["rest"]] <- mean(kp[names(kp) != "adipose"])
  kp[TISSUES]
}
