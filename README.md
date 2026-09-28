# PBPK Simulator

A whole-body physiologically based pharmacokinetic (PBPK) model in R, with an
interactive Shiny dashboard. You pick a compound or enter your own
physicochemical and clearance data, set the subject's age and weight and the
dosing regimen, and it predicts plasma and tissue concentrations, exposure, and
how clearance changes from birth to adulthood.

The model is written for [mrgsolve](https://mrgsolve.org) (`models/pbpk.cpp`).
The browser version solves the same equations exactly with a matrix
exponential, and the test suite checks the two against each other before every
deploy.

This is for research and teaching only. The compound files are illustrations
built from published inputs and calibrated clearances, the pediatric ontogeny
curves are simplified, and nothing here is validated for clinical use.

**Browser version:** <https://wrlog.github.io/PBPK/> (runs entirely in the
browser through WebAssembly, so there's nothing to install).

## Model

Twelve perfusion-limited tissues (lung, brain, heart, kidney, muscle, skin,
adipose, bone, spleen, gut wall, liver and the rest of the body), plus arterial
and venous blood and a gut lumen for oral doses. Spleen and gut drain into the
liver through the portal vein, and the lung takes the whole cardiac output.

- **Distribution.** Tissue:plasma partition coefficients come from the
  mechanistic Rodgers & Rowland equations (J Pharm Sci 2005;94:1259 and
  2006;95:1238). Moderate-to-strong bases bind acidic phospholipids, with the
  affinity back-calculated from the blood:plasma ratio. Acids and weak bases
  bind albumin and neutrals bind lipoproteins. Tissue composition is the
  Rodgers & Rowland data set distributed with Utsey et al. (Drug Metab Dispos
  2020;48:903).
- **Clearance.** Unbound hepatic intrinsic clearance acts on the unbound
  concentration in the liver, and renal clearance is fu × GFR × a renal
  factor. With perfusion-limited kinetics this gives exactly the well-stirred
  liver model, which the tests confirm.
- **Absorption.** First-order from the lumen, with the fraction absorbed (Fa)
  and the gut-wall availability (Fg), so oral F = Fa × Fg × Fh.
- **Physiology.** Organ volumes are fixed fractions of body weight and flows
  are fixed fractions of cardiac output, rounded from the ICRP adult reference
  (ICRP 89) and Brown et al. (1997). Cardiac output scales with weight^0.75.
- **Children.** GFR matures with postmenstrual age (Rhodin et al., Pediatr
  Nephrol 2009;24:67). Hepatic intrinsic clearance scales with liver mass and
  a Hill-type ontogeny curve for the clearing enzyme (CYP3A4, CYP1A2, CYP2D6
  or a generic hepatic profile). The ontogeny parameters are rounded
  illustrations of the published shapes, not fitted values.
- **Variability.** Log-normal on intrinsic clearance, cardiac output and GFR.

This is first-order pediatric scaling. Organ fractions don't change with age
(an infant's brain is a much larger share of body weight than an adult's), so
it isn't a replacement for a full pediatric physiology like PK-Sim or Simcyp.

## Example compounds

| Compound | Class | What it shows |
|---|---|---|
| Midazolam | Weak base | High hepatic and gut extraction through CYP3A4; oral F about 0.3 |
| Caffeine | Neutral | Low extraction through CYP1A2; neonatal half-life more than 10× the adult one |
| Metoprolol | Strong base | Acidic phospholipid binding gives a large volume; CYP2D6 first pass |
| Thiopental | Acid | Redistribution: brain peaks within minutes, then drug moves to muscle and fat |

logP, pKa, fu and B:P are the inputs used by Utsey et al. 2020, except that
caffeine is treated as a neutral. CLint is calibrated so each adult lands near
its reported plasma clearance. With these inputs the adult model gives
midazolam CL 25.6 L/h, Vss 1.45 L/kg and F 0.30; caffeine t½ 4.0 h and Vss
0.47 L/kg; thiopental Vss 2.0 L/kg.

## Dashboard

| Tab | What's on it |
|---|---|
| Plasma exposure | Population median and 50%/90% prediction intervals; Cmax, AUC, clearance and half-life; exposure over the final dosing interval |
| Tissues | Every tissue's concentration against plasma, the Kp values, and each tissue's share of the volume of distribution |
| Age scaling | Clearance per kg from birth to 18 years, profiles by age after the same mg/kg dose, and the mg/kg dose that matches adult AUC at each age |
| Drug & model setup | Physicochemistry, binding, clearance, absorption, ontogeny pathway and variability, plus the Kp table and the maturation curves |

Clearance, Vss, terminal half-life and oral F are exact. They come from the
moments and eigenvalues of the model's rate matrix, not from simulated curves.

## Two engines

`models/pbpk.cpp` holds the structural model for mrgsolve. Organ volumes, flows
and Kp values are computed in R (`R/physiology.R`, `R/partition.R`) and passed
in as parameters, so the same code feeds both engines.

mrgsolve compiles C++, which a browser can't do. The model is linear and
time-invariant, though, so between dosing events the exact solution is a
matrix exponential. `R/pbpk_engine.R` uses that approach, in base R with no
compiled code, and it's what runs under webR.

Run locally with mrgsolve installed, the app simulates through mrgsolve. The
About tab shows which engine is in use.

`tests/test_engine_vs_mrgsolve.R` runs all four compounds by three routes
through both engines, for an adult, an infant and a random subject on
different doses, and fails if they differ by more than 1e-6. The current
worst case is about 2e-8.

## Running it locally

```r
install.packages(c("shiny", "shinydashboard", "DT", "ggplot2", "scales"))
install.packages("mrgsolve")   # optional; needs a C++ toolchain (Rtools on Windows)
shiny::runApp()
```

Tests, from the repository root:

```sh
Rscript tests/test_model.R                # base R only
Rscript tests/test_engine_vs_mrgsolve.R   # needs mrgsolve
```

## Using it from the console

```r
for (f in c("physiology", "partition", "drugs", "pbpk_engine", "simulate"))
  source(file.path("R", paste0(f, ".R")))

# A 2-year-old, 12 kg, given oral midazolam 0.5 mg/kg
p <- pbpk_parameters(DRUGS$midazolam, wt = 12, age_y = 2)
pbpk_metrics(p)                     # CL, Vss, t_half, F_oral, E_H, ...

sim <- pbpk_simulate(list(p),
  regimen = list(route = "oral", amt = 6, interval = 24, n_doses = 1, inf_dur = 0),
  times = seq(0, 24, by = 0.1), keep_states = TRUE)

# The same thing through mrgsolve
source("reference/mrgsolve_engine.R")
pbpk_simulate_mrgsolve(list(p), list(route = "oral", amt = 6, interval = 24,
                                     n_doses = 1, inf_dur = 0), seq(0, 24, by = 0.1))
```

For your own compound, copy one of the `DRUGS` entries in `R/drugs.R` and
change its fields.

## Files

- `app.R`: the Shiny dashboard
- `models/pbpk.cpp`: the mrgsolve model
- `R/physiology.R`: organ volumes, flows, GFR and enzyme maturation
- `R/partition.R`: Rodgers & Rowland partition coefficients
- `R/drugs.R`: example compounds
- `R/pbpk_engine.R`: rate matrix, matrix-exponential solver and exact secondary parameters
- `R/simulate.R`: virtual populations, exposure metrics and age scaling
- `R/theme.R`: styling
- `reference/mrgsolve_engine.R`: runs `models/pbpk.cpp` through the same interface
- `tests/`: model properties and the engine-vs-mrgsolve check
- `.github/workflows/shinylive.yml`: tests, WebAssembly export and GitHub Pages deploy

## License

MIT
