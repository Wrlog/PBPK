$PROB
// Whole-body, perfusion-limited PBPK model.
//
// Twelve tissue compartments plus arterial and venous blood and a gut lumen
// for oral dosing. Spleen and gut drain into the liver through the portal
// vein; every other tissue drains into the venous pool, which passes
// through the lung to the arterial pool.
//
// Elimination:
//   - liver: unbound intrinsic clearance acting on C_liver * fu / Kp_liver
//   - kidney: clearance CLR acting on the plasma leaving the kidney
//   - gut wall: only a fraction FA * FG of absorbed drug reaches the gut
//
// This file is the structural model only. Organ volumes, flows and Kp
// values are computed in R (R/physiology.R, R/partition.R) and passed in
// as parameters, so the same physiology and partition code feeds both this
// model and the matrix-exponential engine in R/pbpk_engine.R.
//
// Defaults are the 70 kg adult with midazolam.
//
// Units: mg, L, L/h, h. Concentrations in mg/L.

$PARAM @annotated
V_lung    : 0.532  : Lung volume (L)
V_brain   : 1.400  : Brain volume (L)
V_heart   : 0.329  : Heart volume (L)
V_kidney  : 0.308  : Kidney volume (L)
V_muscle  : 28.00  : Muscle volume (L)
V_skin    : 2.590  : Skin volume (L)
V_adipose : 14.98  : Adipose volume (L)
V_bone    : 5.992  : Bone volume (L)
V_spleen  : 0.182  : Spleen volume (L)
V_gut     : 1.197  : Gut volume (L)
V_liver   : 1.799  : Liver volume (L)
V_rest    : 7.294  : Rest-of-body volume (L)
V_art     : 1.799  : Arterial blood volume (L)
V_ven     : 3.598  : Venous blood volume (L)

Q_brain   : 46.80  : Brain blood flow (L/h)
Q_heart   : 15.60  : Heart blood flow (L/h)
Q_kidney  : 74.10  : Kidney blood flow (L/h)
Q_muscle  : 66.30  : Muscle blood flow (L/h)
Q_skin    : 19.50  : Skin blood flow (L/h)
Q_adipose : 19.50  : Adipose blood flow (L/h)
Q_bone    : 19.50  : Bone blood flow (L/h)
Q_spleen  : 11.70  : Spleen blood flow (L/h)
Q_gut     : 54.60  : Gut blood flow (L/h)
Q_ha      : 25.35  : Hepatic artery blood flow (L/h)
Q_rest    : 37.05  : Rest-of-body blood flow (L/h)

KP_lung    : 1.28 : Lung:plasma partition coefficient
KP_brain   : 1.71 : Brain:plasma partition coefficient
KP_heart   : 0.88 : Heart:plasma partition coefficient
KP_kidney  : 0.94 : Kidney:plasma partition coefficient
KP_muscle  : 0.58 : Muscle:plasma partition coefficient
KP_skin    : 2.81 : Skin:plasma partition coefficient
KP_adipose : 3.67 : Adipose:plasma partition coefficient
KP_bone    : 0.83 : Bone:plasma partition coefficient
KP_spleen  : 0.57 : Spleen:plasma partition coefficient
KP_gut     : 1.89 : Gut:plasma partition coefficient
KP_liver   : 0.98 : Liver:plasma partition coefficient
KP_rest    : 1.25 : Rest-of-body:plasma partition coefficient

BP    : 0.664 : Blood:plasma concentration ratio
FU    : 0.032 : Fraction unbound in plasma
CLINT : 1400  : Unbound hepatic intrinsic clearance (L/h)
CLR   : 0     : Renal clearance of plasma leaving the kidney (L/h)
KA    : 3.0   : First-order absorption rate constant (1/h)
FA    : 0.88  : Fraction absorbed from the lumen
FG    : 0.59  : Fraction escaping gut-wall metabolism

$CMT @annotated
LUMEN   : Gut lumen, oral doses (mg)
LUNG    : Lung (mg)
BRAIN   : Brain (mg)
HEART   : Heart (mg)
KIDNEY  : Kidney (mg)
MUSCLE  : Muscle (mg)
SKIN    : Skin (mg)
ADIPOSE : Adipose (mg)
BONE    : Bone (mg)
SPLEEN  : Spleen (mg)
GUT     : Gut wall (mg)
LIVER   : Liver (mg)
REST    : Rest of body (mg)
ART     : Arterial blood, IV doses do not go here (mg)
VEN     : Venous blood, IV doses go here (mg)

$MAIN
double Q_liver = Q_ha + Q_spleen + Q_gut;
double CO = Q_brain + Q_heart + Q_kidney + Q_muscle + Q_skin + Q_adipose +
            Q_bone + Q_rest + Q_liver;

$ODE
// Concentration in venous blood leaving each tissue: C_tissue * BP / Kp.
double Cart = ART / V_art;
double Cven = VEN / V_ven;

double Oout_lung    = LUNG    / V_lung    * BP / KP_lung;
double Oout_brain   = BRAIN   / V_brain   * BP / KP_brain;
double Oout_heart   = HEART   / V_heart   * BP / KP_heart;
double Oout_kidney  = KIDNEY  / V_kidney  * BP / KP_kidney;
double Oout_muscle  = MUSCLE  / V_muscle  * BP / KP_muscle;
double Oout_skin    = SKIN    / V_skin    * BP / KP_skin;
double Oout_adipose = ADIPOSE / V_adipose * BP / KP_adipose;
double Oout_bone    = BONE    / V_bone    * BP / KP_bone;
double Oout_spleen  = SPLEEN  / V_spleen  * BP / KP_spleen;
double Oout_gut     = GUT     / V_gut     * BP / KP_gut;
double Oout_liver   = LIVER   / V_liver   * BP / KP_liver;
double Oout_rest    = REST    / V_rest    * BP / KP_rest;

double Cu_liver     = LIVER  / V_liver  * FU / KP_liver;
double Cp_kidney    = KIDNEY / V_kidney / KP_kidney;

dxdt_LUMEN   = -KA * LUMEN;
dxdt_LUNG    = CO * (Cven - Oout_lung);
dxdt_BRAIN   = Q_brain   * (Cart - Oout_brain);
dxdt_HEART   = Q_heart   * (Cart - Oout_heart);
dxdt_KIDNEY  = Q_kidney  * (Cart - Oout_kidney) - CLR * Cp_kidney;
dxdt_MUSCLE  = Q_muscle  * (Cart - Oout_muscle);
dxdt_SKIN    = Q_skin    * (Cart - Oout_skin);
dxdt_ADIPOSE = Q_adipose * (Cart - Oout_adipose);
dxdt_BONE    = Q_bone    * (Cart - Oout_bone);
dxdt_SPLEEN  = Q_spleen  * (Cart - Oout_spleen);
dxdt_GUT     = Q_gut     * (Cart - Oout_gut) + KA * FA * FG * LUMEN;
dxdt_LIVER   = Q_ha * Cart + Q_spleen * Oout_spleen + Q_gut * Oout_gut -
               Q_liver * Oout_liver - CLINT * Cu_liver;
dxdt_REST    = Q_rest    * (Cart - Oout_rest);
dxdt_ART     = CO * (Oout_lung - Cart);
dxdt_VEN     = Q_brain * Oout_brain + Q_heart * Oout_heart +
               Q_kidney * Oout_kidney + Q_muscle * Oout_muscle +
               Q_skin * Oout_skin + Q_adipose * Oout_adipose +
               Q_bone * Oout_bone + Q_rest * Oout_rest +
               Q_liver * Oout_liver - CO * Cven;

$TABLE
double CP = VEN / V_ven / BP;

$CAPTURE @annotated
CP : Venous plasma concentration (mg/L)
