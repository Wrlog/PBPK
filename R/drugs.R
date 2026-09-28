# ============================================================================
# Example compounds.
#
# One of each class the partition model treats differently: a weak base, a
# neutral, a strong base and an acid. logP, pKa, fu and B:P are the inputs
# used by Utsey et al. (Drug Metab Dispos 2020;48:903) except where noted.
# CLint is the unbound hepatic intrinsic clearance, chosen so the adult
# model lands near the reported plasma clearance; it is a calibration, not
# an in vitro value. ka, Fa and Fg are rounded literature-style values.
#
# These are illustrations for exploring the model, not validated compound
# files.
# ============================================================================

DRUGS <- list(
  midazolam = list(
    label = "Midazolam", class = "Weak base",
    type = "base", logP = 3.13, pKa = 6.0, fu = 0.032, BP = 0.664,
    clint = 1400, pathway = "CYP3A4", renal_factor = 0,
    ka = 3.0, fa = 0.88, fg = 0.59,
    route = "oral", dose = 7.5, dose_basis = "flat", interval = 24, n_doses = 1,
    duration = 24,
    note = "CYP3A4 probe with high hepatic and gut extraction; adult plasma CL about 25 L/h."
  ),
  caffeine = list(
    label = "Caffeine", class = "Neutral",
    # Caffeine is essentially unionised at physiological pH, so it is treated
    # as a neutral here rather than as the base used in Utsey et al.
    type = "neutral", logP = -0.07, pKa = NA, fu = 0.681, BP = 0.98,
    clint = 9, pathway = "CYP1A2", renal_factor = 0.02,
    ka = 2.5, fa = 1, fg = 1,
    route = "oral", dose = 150, dose_basis = "flat", interval = 24, n_doses = 1,
    duration = 48,
    note = "Low-extraction CYP1A2 substrate; the slow CYP1A2 maturation makes neonatal clearance an order of magnitude lower per kg."
  ),
  metoprolol = list(
    label = "Metoprolol", class = "Strong base",
    type = "base", logP = 2.15, pKa = 9.7, fu = 0.879, BP = 1.52,
    clint = 120, pathway = "CYP2D6", renal_factor = 0.1,
    ka = 1.0, fa = 1, fg = 1,
    route = "oral", dose = 100, dose_basis = "flat", interval = 12, n_doses = 6,
    duration = 72,
    note = "CYP2D6 substrate with extensive first-pass extraction; binding to acidic phospholipids gives it a large volume."
  ),
  thiopental = list(
    label = "Thiopental", class = "Acid",
    type = "acid", logP = 2.9, pKa = 7.5, fu = 0.13, BP = 1.0,
    clint = 130, pathway = "Generic hepatic", renal_factor = 0,
    ka = 1.0, fa = 1, fg = 1,
    route = "iv_bolus", dose = 4, dose_basis = "mgkg", interval = 24, n_doses = 1,
    duration = 24,
    note = "The classic redistribution example: brain concentrations peak within minutes and fall as the drug moves into muscle and fat."
  )
)

drug_choices <- stats::setNames(names(DRUGS),
                                vapply(DRUGS, function(d) sprintf("%s (%s)", d$label, tolower(d$class)), ""))
