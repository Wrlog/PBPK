# Check the matrix-exponential engine against mrgsolve running
# models/pbpk.cpp. If they agree, the browser build computes what the
# mrgsolve model computes.
#
# Run from the repository root: Rscript tests/test_engine_vs_mrgsolve.R

for (f in c("physiology", "partition", "drugs", "pbpk_engine", "simulate")) {
  source(file.path("R", paste0(f, ".R")))
}
source(file.path("reference", "mrgsolve_engine.R"))
if (!mrgsolve_ready()) stop("mrgsolve is not installed")

# Off-grid output times, so no comparison lands exactly on a dose time
# (where one engine reports the pre-dose value and the other the post-dose).
times <- seq(0, 36, by = 0.25) + 0.0137

set.seed(7)
worst <- 0
cases <- 0
for (d in names(DRUGS)) {
  for (route in c("oral", "iv_bolus", "iv_inf")) {
    # An adult, an infant and a random subject with variability.
    plist <- list(
      pbpk_parameters(DRUGS[[d]], 70, 30),
      pbpk_parameters(DRUGS[[d]], 7.8, 0.5),
      pbpk_parameters(DRUGS[[d]], runif(1, 15, 90), runif(1, 2, 60),
                      clint_mult = exp(rnorm(1, 0, 0.4)),
                      co_mult = exp(rnorm(1, 0, 0.1)),
                      gfr_mult = exp(rnorm(1, 0, 0.2)))
    )
    reg <- list(route = route, amt = c(100, 11, 60), interval = 8,
                n_doses = 4, inf_dur = 1.5)
    got <- pbpk_simulate(plist, reg, times)$plasma
    want <- pbpk_simulate_mrgsolve(plist, reg, times)$plasma
    rel <- max(abs(got - want)) / max(abs(want))
    cat(sprintf("  %-11s %-9s max relative difference %.2e\n", d, route, rel))
    worst <- max(worst, rel)
    cases <- cases + 1
  }
}

cat(sprintf("cases compared: %d\n", cases))
cat(sprintf("worst relative difference vs mrgsolve: %.3e\n", worst))
if (!is.finite(worst) || worst > 1e-6) {
  stop(sprintf("matrix-exponential engine disagrees with mrgsolve (%.3e)", worst))
}
cat("PASS\n")
