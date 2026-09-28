# ============================================================================
# Physiology: organ volumes and blood flows, and how they scale with age.
#
# The adult reference is a 70 kg subject. Organ volumes are fixed fractions
# of body weight and blood flows are fixed fractions of cardiac output, which
# scales allometrically on weight. Values are rounded from the ICRP adult
# reference male (ICRP Publication 89, 2002) and Brown et al. (Toxicol Ind
# Health 1997;13:407), with a density of 1 kg/L throughout.
#
# Children get two maturation functions on top of the size scaling:
#   - glomerular filtration: Rhodin et al., Pediatr Nephrol 2009;24:67
#   - hepatic enzyme activity: a Hill function of postnatal age per pathway
#
# This is first-order pediatric scaling. Keeping organ fractions fixed
# ignores, for example, the much larger brain fraction in an infant; a full
# pediatric physiology (PK-Sim, Simcyp) varies every organ with age.
# ============================================================================

REF_WT <- 70

# Fraction of body weight. "rest" is whatever is left over.
VOL_FRAC <- c(
  lung    = 0.0076,
  brain   = 0.0200,
  heart   = 0.0047,
  kidney  = 0.0044,
  muscle  = 0.4000,
  skin    = 0.0370,
  adipose = 0.2140,
  bone    = 0.0856,
  spleen  = 0.0026,
  gut     = 0.0171,
  liver   = 0.0257,
  art     = 0.0257,
  ven     = 0.0514
)

# Fraction of cardiac output. Spleen and gut drain into the portal vein, so
# total liver blood flow is hepatic artery + spleen + gut. "rest" takes the
# remainder so that venous return equals cardiac output.
FLOW_FRAC <- c(
  brain          = 0.120,
  heart          = 0.040,
  kidney         = 0.190,
  muscle         = 0.170,
  skin           = 0.050,
  adipose        = 0.050,
  bone           = 0.050,
  spleen         = 0.030,
  gut            = 0.140,
  hepatic_artery = 0.065
)

CO_REF  <- 390                   # cardiac output at 70 kg (L/h)
GFR_REF <- 121.2 * 60 / 1000     # adult GFR, 121.2 mL/min per 70 kg (L/h)
HCT     <- 0.45                  # haematocrit, used by the partition model

# The tissues that get their own compartment, in the order used everywhere.
TISSUES <- c("lung", "brain", "heart", "kidney", "muscle", "skin", "adipose",
             "bone", "spleen", "gut", "liver", "rest")

#' Fraction of adult GFR reached at a given postnatal age
#'
#' Sigmoid in postmenstrual age (weeks) with TM50 = 47.7 weeks and a Hill
#' coefficient of 3.4, assuming a term birth at 40 weeks.
gfr_maturation <- function(age_y) {
  pma_wk <- age_y * 52.18 + 40
  pma_wk^3.4 / (47.7^3.4 + pma_wk^3.4)
}

# Enzyme ontogeny as a fraction of adult activity per gram of liver:
#   f(age) = f_birth + (1 - f_birth) * age^hill / (tm50^hill + age^hill)
# with age in postnatal years. The shapes follow the published pattern for
# each enzyme (CYP1A2 slowest, CYP2D6 fastest), but the numbers are rounded
# illustrations rather than fitted values. Check them against a current
# ontogeny source before quantitative pediatric use.
ONTOGENY <- list(
  "Mature (no ontogeny)" = c(f_birth = 1.00, tm50 = 1.00, hill = 1.0),
  "CYP3A4"               = c(f_birth = 0.10, tm50 = 0.50, hill = 1.0),
  "CYP1A2"               = c(f_birth = 0.05, tm50 = 0.60, hill = 1.5),
  "CYP2D6"               = c(f_birth = 0.20, tm50 = 0.10, hill = 1.0),
  "Generic hepatic"      = c(f_birth = 0.30, tm50 = 0.50, hill = 1.0)
)

enzyme_maturation <- function(age_y, pathway) {
  o <- ONTOGENY[[pathway]]
  if (is.null(o)) stop("unknown ontogeny pathway: ", pathway)
  a <- pmax(age_y, 0)
  o[["f_birth"]] + (1 - o[["f_birth"]]) * a^o[["hill"]] / (o[["tm50"]]^o[["hill"]] + a^o[["hill"]])
}

# Rough median body weight by age, used only to prefill the weight input
# when the age changes. Interpolated from growth-chart medians (sexes pooled).
WT_FOR_AGE <- data.frame(
  age = c(0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 25),
  wt  = c(3.4, 6.0, 7.8, 9.6, 12.2, 16.3, 20.5, 25.6, 32, 40, 50, 60, 67, 70)
)

typical_weight <- function(age_y) {
  stats::approx(WT_FOR_AGE$age, WT_FOR_AGE$wt, xout = age_y, rule = 2)$y
}

#' Organ volumes (L), blood flows (L/h) and GFR (L/h) for one subject
physiology <- function(wt, age_y = 30) {
  v <- VOL_FRAC * wt
  v[["rest"]] <- wt * (1 - sum(VOL_FRAC))

  co <- CO_REF * (wt / REF_WT)^0.75
  q <- FLOW_FRAC * co
  q[["rest"]] <- co * (1 - sum(FLOW_FRAC))
  q[["liver"]] <- q[["hepatic_artery"]] + q[["spleen"]] + q[["gut"]]
  q[["lung"]] <- co

  list(
    wt = wt,
    age = age_y,
    V = v,
    Q = q,
    CO = co,
    GFR = GFR_REF * (wt / REF_WT)^0.75 * gfr_maturation(age_y)
  )
}
