# Properties the model must have whatever the compound. Base R only.
#
# Run from the repository root: Rscript tests/test_model.R

for (f in c("physiology", "partition", "drugs", "pbpk_engine", "simulate")) {
  source(file.path("R", paste0(f, ".R")))
}

check <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  cat("  ok:", msg, "\n")
}

# --- Physiology ------------------------------------------------------------

phys <- physiology(70, 30)
check(abs(sum(phys$V) - 70) < 1e-9, "organ volumes add up to body weight")
venous_return <- sum(phys$Q[c(TO_VENOUS, "liver")])
check(abs(venous_return - phys$CO) < 1e-9, "venous return equals cardiac output")
check(abs(gfr_maturation(1000) - 1) < 1e-3, "GFR maturation reaches adult values")
check(all(vapply(names(ONTOGENY), function(k) abs(enzyme_maturation(200, k) - 1) < 0.01, TRUE)),
      "every enzyme ontogeny curve reaches adult activity")

# --- Partition coefficients ---------------------------------------------------

for (d in names(DRUGS)) {
  x <- DRUGS[[d]]
  kp <- kp_rodgers_rowland(x$logP, x$pKa, x$type, x$fu, x$BP)
  check(all(is.finite(kp)) && all(kp > 0), sprintf("%s: every Kp is positive", d))
}
kp_lo <- kp_rodgers_rowland(1, NA, "neutral", 0.5, 1)
kp_hi <- kp_rodgers_rowland(4, NA, "neutral", 0.5, 1)
check(kp_hi[["adipose"]] > 10 * kp_lo[["adipose"]], "a more lipophilic neutral partitions further into fat")

# --- Mass balance ------------------------------------------------------------------

# With every elimination route switched off, the dose has nowhere to go.
closed <- DRUGS$thiopental
closed$clint <- 0; closed$renal_factor <- 0; closed$fa <- 1; closed$fg <- 1
p <- pbpk_parameters(closed, 70, 30)
r <- pbpk_simulate(list(p), list(route = "oral", amt = 100, interval = 12, n_doses = 2,
                                 inf_dur = 0), c(0, 5, 30), keep_states = TRUE)
vols <- c(1, p$V[c("lung", "brain", "heart", "kidney", "muscle", "skin", "adipose",
                   "bone", "spleen", "gut", "liver", "rest", "art", "ven")])
total <- r$states[3, , 1] %*% vols
check(abs(total - 200) < 1e-6, "no elimination: every mg dosed is still in the body")

# --- Linearity and bioavailability ------------------------------------------------

p <- pbpk_parameters(DRUGS$midazolam, 70, 30)
tt <- seq(0, 24, by = 0.5)
one <- pbpk_simulate(list(p), list(route = "iv_bolus", amt = 1, interval = 24, n_doses = 1, inf_dur = 0), tt)$plasma
two <- pbpk_simulate(list(p), list(route = "iv_bolus", amt = 2, interval = 24, n_doses = 1, inf_dur = 0), tt)$plasma
check(max(abs(two - 2 * one)) < 1e-12, "concentrations are proportional to dose")

m <- pbpk_metrics(p)
# The IV bolus falls tenfold within minutes, so the trapezoid needs a fine
# grid early on to integrate it accurately.
tt <- c(seq(0, 2, by = 1e-4), seq(2.05, 400, by = 0.05))
iv <- pbpk_simulate(list(p), list(route = "iv_bolus", amt = 1, interval = 1e3, n_doses = 1, inf_dur = 0), tt)$plasma[, 1]
po <- pbpk_simulate(list(p), list(route = "oral", amt = 1, interval = 1e3, n_doses = 1, inf_dur = 0), tt)$plasma[, 1]
check(abs(trapz(tt, po) / trapz(tt, iv) - m$F_oral) < 1e-3, "simulated oral/IV AUC ratio matches F from the rate matrix")
check(abs(1 / trapz(tt, iv) - m$CL) / m$CL < 1e-3, "simulated clearance matches CL from the rate matrix")

fh <- 1 - m$E_H
check(abs(m$F_oral - p$fa * p$fg * fh) < 1e-6, "F = Fa x Fg x Fh with well-stirred hepatic extraction")

# Well-stirred liver: with no renal clearance, plasma CL is Q fub CLint / (Q + fub CLint) x BP.
fub <- p$fu / p$BP
cl_ws <- p$Q[["liver"]] * fub * p$CLint / (p$Q[["liver"]] + fub * p$CLint) * p$BP
check(abs(m$CL - cl_ws) / cl_ws < 1e-9, "hepatic clearance follows the well-stirred model")

# --- Pediatrics -------------------------------------------------------------------

caf <- age_scaling_table(DRUGS$caffeine, c(7 / 365, 1, 12), 1, oral = FALSE)$table
check(all(diff(caf$t_half) < 0), "caffeine half-life shortens as CYP1A2 matures")

cat("PASS\n")
