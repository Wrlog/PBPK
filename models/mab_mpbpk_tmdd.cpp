$PROB
Second-generation minimal PBPK model for monoclonal antibodies, with
target-mediated drug disposition in the interstitial fluid

// Cao Y, Balthasar JP, Jusko WJ. Second-generation minimal physiologically-
//   based pharmacokinetic model for monoclonal antibodies.
//   J Pharmacokinet Pharmacodyn 2013;40:597-607.
// Cao Y, Jusko WJ. Incorporating target-mediated drug disposition in a
//   minimal physiologically-based pharmacokinetic model for monoclonal
//   antibodies. J Pharmacokinet Pharmacodyn 2014;41:375-387 ("Model B",
//   pTMDD: quasi-steady-state binding to targets in tight and leaky ISF).
//
// Plasma exchanges with the interstitial fluid of two lumped tissue groups
// by convection through vascular pores: "tight" (continuous endothelium:
// muscle, skin, adipose, brain) with reflection coefficient SIGMA1, and
// "leaky" (fenestrated/discontinuous: liver, kidney, heart, ...) with
// SIGMA2. ISF drains through lymph (reflection SIGMA_L = 0.2) back to
// plasma. Nonspecific clearance is from plasma (CLp). With TMDD on, the
// target is synthesised (ksyn) and degraded (kdeg) in both ISF spaces,
// binds the antibody (quasi-steady-state constant Kss) and the complex is
// internalised (kint); only free antibody drains to lymph.
//
// Physiology for a 70 kg adult (Cao 2013/2014): plasma 2.6 L, ISF 15.6 L,
// lymph 5.2 L, lymph flow 2.9 L/day, available ISF fraction Kp 0.8 (IgG1).
// Volumes and flows scale linearly with body weight.
//
// Defaults: trastuzumab in patients (Cao & Jusko 2014, Table 2, pTMDD).
// Units: time in hours, antibody in nmol and nM (dose mg x 1e6 / MW),
// target in nM.

$PARAM @annotated
BW       : 70      : Body weight (kg)
MW       : 148000  : Molecular weight (g/mol)
Vp_70    : 2.6     : Plasma volume at 70 kg (L)
ISF_70   : 15.6    : Interstitial fluid volume at 70 kg (L)
Vlymph_70 : 5.2    : Lymph volume at 70 kg (L)
L_70     : 0.120833 : Total lymph flow at 70 kg (L/h; 2.9 L/day)
Kp       : 0.8     : Fraction of ISF available to the antibody
SIGMA1   : 0.95    : Vascular reflection coefficient, tight tissues
SIGMA2   : 0.512   : Vascular reflection coefficient, leaky tissues
SIGMAL   : 0.2     : Lymphatic reflection coefficient
CLp_kg   : 0.092   : Nonspecific plasma clearance (mL/h/kg)
TMDD     : 1       : 1 = target binding on, 0 = linear antibody
KSS      : 0.99    : Quasi-steady-state constant (nM)
KSYN     : 0.376   : Target synthesis rate (nM/h)
KDEG     : 0.0117  : Free target degradation rate (1/h)
KINT     : 0.0117  : Complex internalisation rate (1/h)

$CMT @annotated
CENT      : Antibody in plasma (nmol)
TIGHT     : Total antibody in tight ISF, free + bound (nmol)
LEAKY     : Total antibody in leaky ISF, free + bound (nmol)
LYMPH     : Antibody in lymph (nmol)
RT_TIGHT  : Total target in tight ISF (nM)
RT_LEAKY  : Total target in leaky ISF (nM)

$MAIN
double Vp = Vp_70 * BW / 70;
double ISF = ISF_70 * BW / 70;
double Vlymph = Vlymph_70 * BW / 70;
double L = L_70 * BW / 70;
double L1 = 0.33 * L;
double L2 = 0.67 * L;
double Vtight = 0.65 * ISF * Kp;
double Vleaky = 0.35 * ISF * Kp;
double CLp = CLp_kg * BW / 1000;
double R0 = TMDD > 0.5 ? KSYN / KDEG : 0.0;
RT_TIGHT_0 = R0;
RT_LEAKY_0 = R0;

$ODE
double Cp = CENT / Vp;
double Clymph = LYMPH / Vlymph;

// free antibody in each ISF space from total antibody and total target
double CtotT = TIGHT / Vtight;
double bT = CtotT - RT_TIGHT - KSS;
double dT = sqrt(bT * bT + 4.0 * KSS * CtotT);
double CfreeT = bT >= 0 ? 0.5 * (bT + dT) : 2.0 * KSS * CtotT / (dT - bT);
double CtotL = LEAKY / Vleaky;
double bL = CtotL - RT_LEAKY - KSS;
double dL = sqrt(bL * bL + 4.0 * KSS * CtotL);
double CfreeL = bL >= 0 ? 0.5 * (bL + dL) : 2.0 * KSS * CtotL / (dL - bL);
double ART = CtotT - CfreeT;
double ARL = CtotL - CfreeL;

dxdt_CENT = -(1 - SIGMA1) * L1 * Cp - (1 - SIGMA2) * L2 * Cp + L * Clymph - CLp * Cp;
dxdt_TIGHT = (1 - SIGMA1) * L1 * Cp - (1 - SIGMAL) * L1 * CfreeT - KINT * ART * Vtight;
dxdt_LEAKY = (1 - SIGMA2) * L2 * Cp - (1 - SIGMAL) * L2 * CfreeL - KINT * ARL * Vleaky;
dxdt_LYMPH = (1 - SIGMAL) * L1 * CfreeT + (1 - SIGMAL) * L2 * CfreeL - L * Clymph;
dxdt_RT_TIGHT = TMDD * (KSYN - KDEG * RT_TIGHT - (KINT - KDEG) * ART);
dxdt_RT_LEAKY = TMDD * (KSYN - KDEG * RT_LEAKY - (KINT - KDEG) * ARL);

$TABLE
double CP_UGML = CENT / Vp * MW / 1e6;                         // plasma (ug/mL)
double RO_LEAKY = RT_LEAKY > 0 ? ARL / RT_LEAKY : 0.0;         // occupancy, leaky ISF
double RO_TIGHT = RT_TIGHT > 0 ? ART / RT_TIGHT : 0.0;         // occupancy, tight ISF
double ISF_LEAKY_UGML = CfreeL * MW / 1e6;                     // free, leaky ISF (ug/mL)
double ISF_TIGHT_UGML = CfreeT * MW / 1e6;                     // free, tight ISF (ug/mL)

$CAPTURE CP_UGML RO_LEAKY RO_TIGHT ISF_LEAKY_UGML ISF_TIGHT_UGML
