$PROB
Rifampicin -> midazolam DDI: CYP3A induction in liver and gut

// Asaumi R, Toshimoto K, Tobe Y, et al. Comprehensive PBPK model of
// rifampicin for quantitative prediction of complex drug-drug
// interactions: CYP3A/2C9 induction and OATP inhibition effects.
// CPT Pharmacometrics Syst Pharmacol 2018;7:186-196 (open access).
//
// Transcribed from the model code in the article's Supplementary Text
// ("Model code of midazolam in Napp in the case of representative DDI with
// rifampicin"). Parameter values are from Supplementary Tables S1-S2 and the
// optimised values of the article's Table 1.
//
// Rifampicin: blood, 5 hepatic extracellular + 5 hepatocyte compartments
// (dispersion-type liver) with saturable OATP uptake and passive diffusion,
// muscle, skin, adipose, serosa, and a segregated-flow gut (lumen,
// enterocytes, mucosal blood). Its unbound hepatocyte and enterocyte
// concentrations induce UGT (its own metabolism, auto-induction) and CYP3A
// through an enzyme-turnover model (Emax, EC50,u, kdeg).
// Midazolam: blood, 5 liver compartments, muscle, skin, adipose, portal
// vein and gut lumen (Qgut model). CYP3A induction raises its hepatic and
// intestinal intrinsic clearance in proportion to fm.
//
// Concentrations are total (bound + unbound) blood concentrations in mg/L
// (= ug/mL); gut lumen states are amounts (mg). Time in hours.
// IV doses go into the blood concentration state, so divide mg by the
// blood / central volume (Vb, Vcentral) before dosing.

$PARAM @annotated
BW          : 70      : Body weight (kg)

// physiology (Table S1), per kg
Qh_kg       : 1.24    : Hepatic blood flow (L/h/kg)
Qm_kg       : 0.642   : Muscle blood flow (L/h/kg)
Qs_kg       : 0.257   : Skin blood flow (L/h/kg)
Qa_kg       : 0.223   : Adipose blood flow (L/h/kg)
Qser_kg     : 0.274   : Serosal blood flow (L/h/kg)
Qvilli_kg   : 0.257   : Villous (mucosal) blood flow (L/h/kg)
Qpv_kg      : 0.531   : Portal vein blood flow (L/h/kg)
Vb_kg       : 0.0743  : Blood volume (L/kg)
Vh_kg       : 0.0174  : Hepatocyte volume (L/kg)
Vi_kg       : 0.0067  : Hepatic extracellular volume (L/kg)
Vm_kg       : 0.429   : Muscle volume (L/kg)
Vs_kg       : 0.111   : Skin volume (L/kg)
Va_kg       : 0.143   : Adipose volume (L/kg)
Vser_kg     : 0.00893 : Serosa volume (L/kg)
Vmucb_kg    : 0.00099 : Mucosal blood volume (L/kg)
Vent_kg     : 0.00739 : Enterocyte volume (L/kg)
Vpv_kg      : 0.001   : Portal vein volume (L/kg)

// enzyme turnover (Table S1)
kdeg3Ah     : 0.0158  : CYP3A degradation, hepatocytes (1/h)
kdeg3Ag     : 0.0288  : CYP3A degradation, enterocytes (1/h)
kdegUGTh    : 0.0158  : UGT degradation, hepatocytes (1/h)
kdegUGTg    : 0.0288  : UGT degradation, enterocytes (1/h)

// rifampicin (Table S2; optimised values from Table 1)
KaRif       : 37.6    : Absorption rate (1/h) [Table 1]
FaRif       : 1       : Fraction absorbed
FgRif       : 0.943   : Intestinal availability
fbRif       : 0.0778  : Unbound fraction in blood
fhRif       : 0.0814  : Unbound fraction in hepatocytes
fgRif       : 0.115   : Unbound fraction in enterocytes (fE)
KpScaleRif  : 6.65    : Scaling factor on in-silico Kp [Table 1]
KpmRif      : 0.0947  : Kp muscle
KpsRif      : 0.326   : Kp skin
KpaRif      : 0.0629  : Kp adipose
KpserRif    : 0.200   : Kp serosa
fbCLintall_kg : 0.251 : fB x overall hepatic intrinsic clearance (L/h/kg) [Table 1]
RdifRif     : 0.129   : PSdif,inf / PSact,inf
betaRif     : 0.5     : CLint,all / (PSact,inf + PSdif,inf) (0.2, 0.5 or 0.8 in the paper)
gammaRif    : 0.778   : PSdif,inf / PSdif,eff
KmuOATP     : 0.146   : Unbound Km for OATP uptake (ug/mL) [Table 1]
PSdifent_kg : 0.161   : Enterocyte basolateral diffusion clearance (L/h/kg) [Table 1]
CLr_kg      : 0.011   : Renal clearance (L/h/kg)
fmUGTh      : 0.759   : Fraction metabolised by UGT, liver
fmUGTg      : 0.759   : Fraction metabolised by UGT, gut
EC50u       : 0.0526  : Unbound EC50 for induction (ug/mL) [Table 1]
EmaxUGT     : 1.34    : Emax, UGT auto-induction [Table 1]
Emax3A      : 4.57    : Emax, CYP3A induction [Table 1]

// midazolam (Table S2; optimised values from Table 1)
KaMid       : 1.29    : Absorption rate (1/h) [Table 1]
FaMid       : 1       : Fraction absorbed
Vc_kg_mid   : 0.571   : Central volume (L/kg) [Table 1]
KpSFMid     : 0.201   : Scaling factor on in-silico Kp [Table 1]
KphMid      : 6.96    : Kp liver
KpmMid      : 4.09    : Kp muscle
KpsMid      : 20.4    : Kp skin
KpaMid      : 34.4    : Kp adipose
fbCLinth_kg : 0.469   : fB x hepatic intrinsic clearance (L/h/kg) [Table 1]
fgCLintg_kg : 0.107   : fE x enterocyte intrinsic clearance (L/h/kg)
CLperm_kg   : 0.151   : Gut permeability clearance, Qgut model (L/h/kg)
fm3Ah       : 0.93    : Fraction metabolised by CYP3A, liver
fm3Ag       : 1       : Fraction metabolised by CYP3A, gut
CLrMid_kg   : 0       : Renal clearance (L/h/kg)

$CMT @annotated
MDZ_CENT  : Midazolam blood (mg/L)
MDZ_LIV1  : Midazolam liver 1 (mg/L)
MDZ_LIV2  : Midazolam liver 2 (mg/L)
MDZ_LIV3  : Midazolam liver 3 (mg/L)
MDZ_LIV4  : Midazolam liver 4 (mg/L)
MDZ_LIV5  : Midazolam liver 5 (mg/L)
MDZ_MUS   : Midazolam muscle (mg/L)
MDZ_SKIN  : Midazolam skin (mg/L)
MDZ_ADI   : Midazolam adipose (mg/L)
MDZ_PV    : Midazolam portal vein (mg/L)
MDZ_GUT   : Midazolam gut lumen (mg)
RIF_CENT  : Rifampicin blood (mg/L)
RIF_HE1   : Rifampicin hepatic extracellular 1 (mg/L)
RIF_HC1   : Rifampicin hepatocytes 1 (mg/L)
RIF_HE2   : Rifampicin hepatic extracellular 2 (mg/L)
RIF_HC2   : Rifampicin hepatocytes 2 (mg/L)
RIF_HE3   : Rifampicin hepatic extracellular 3 (mg/L)
RIF_HC3   : Rifampicin hepatocytes 3 (mg/L)
RIF_HE4   : Rifampicin hepatic extracellular 4 (mg/L)
RIF_HC4   : Rifampicin hepatocytes 4 (mg/L)
RIF_HE5   : Rifampicin hepatic extracellular 5 (mg/L)
RIF_HC5   : Rifampicin hepatocytes 5 (mg/L)
RIF_SER   : Rifampicin serosa (mg/L)
RIF_MUS   : Rifampicin muscle (mg/L)
RIF_SKIN  : Rifampicin skin (mg/L)
RIF_ADI   : Rifampicin adipose (mg/L)
RIF_LUMEN : Rifampicin gut lumen (mg)
RIF_ENT   : Rifampicin enterocytes (mg/L)
RIF_MUCB  : Rifampicin mucosal blood (mg/L)
UGT_H1    : UGT activity ratio, hepatocytes 1
UGT_H2    : UGT activity ratio, hepatocytes 2
UGT_H3    : UGT activity ratio, hepatocytes 3
UGT_H4    : UGT activity ratio, hepatocytes 4
UGT_H5    : UGT activity ratio, hepatocytes 5
UGT_G     : UGT activity ratio, enterocytes
CYP3A_H1  : CYP3A activity ratio, hepatocytes 1
CYP3A_H2  : CYP3A activity ratio, hepatocytes 2
CYP3A_H3  : CYP3A activity ratio, hepatocytes 3
CYP3A_H4  : CYP3A activity ratio, hepatocytes 4
CYP3A_H5  : CYP3A activity ratio, hepatocytes 5
CYP3A_G   : CYP3A activity ratio, enterocytes

$MAIN
double Qh = Qh_kg * BW;
double Qm = Qm_kg * BW;
double Qs = Qs_kg * BW;
double Qa = Qa_kg * BW;
double Qser = Qser_kg * BW;
double Qvilli = Qvilli_kg * BW;
double Qhart = Qh - Qser - Qvilli;
double Qpv = Qpv_kg * BW;
double VbRif = Vb_kg * BW;
double VbMid = Vc_kg_mid * BW;
double Vh = Vh_kg * BW;
double Vi = Vi_kg * BW;
double Vm = Vm_kg * BW;
double Vs = Vs_kg * BW;
double Va = Va_kg * BW;
double Vser = Vser_kg * BW;
double Vmucb = Vmucb_kg * BW;
double Vent = Vent_kg * BW;
double Vpv = Vpv_kg * BW;

// rifampicin hepatic clearances from the hybrid parameters
double CLintallRif = fbCLintall_kg * BW / fbRif;
double VmaxOATP = 1.0 / (1 + RdifRif) * CLintallRif / betaRif * KmuOATP;
double PSdifinf = RdifRif / (1 + RdifRif) * CLintallRif / betaRif;
double PSdifeff = RdifRif / (1 + RdifRif) * CLintallRif / betaRif / gammaRif;
double CLintRif = RdifRif / (1 + RdifRif) * CLintallRif / (1 - betaRif) / gammaRif;
double PSdifent = PSdifent_kg * BW;
double QgutRif = fgRif * PSdifent * Qvilli / (Qvilli + fbRif * PSdifent);
double CLmetgRif = (QgutRif * (1.0 / FgRif - 1) - (1 - FaRif) * fgRif * PSdifent * 20) / fgRif;
double CLrRif = CLr_kg * BW;

// midazolam
double fbCLinthMid = fbCLinth_kg * BW;
double fgCLintgMid = fgCLintg_kg * BW;
double CLpermMid = CLperm_kg * BW;
double QgutMid = Qvilli * CLpermMid / (Qvilli + CLpermMid);
double CLrMid = CLrMid_kg * BW;
double KPLIV = KpSFMid * KphMid;

UGT_H1_0 = 1; UGT_H2_0 = 1; UGT_H3_0 = 1; UGT_H4_0 = 1; UGT_H5_0 = 1; UGT_G_0 = 1;
CYP3A_H1_0 = 1; CYP3A_H2_0 = 1; CYP3A_H3_0 = 1; CYP3A_H4_0 = 1; CYP3A_H5_0 = 1; CYP3A_G_0 = 1;

$ODE
// ---- midazolam ----
dxdt_MDZ_CENT = 1.0 / VbMid * (Qh * MDZ_LIV5 / KPLIV - (Qh - Qpv) * MDZ_CENT +
    Qm * (MDZ_MUS / (KpSFMid * KpmMid) - MDZ_CENT) + Qs * (MDZ_SKIN / (KpSFMid * KpsMid) - MDZ_CENT) +
    Qa * (MDZ_ADI / (KpSFMid * KpaMid) - MDZ_CENT) - Qpv * MDZ_CENT - CLrMid * MDZ_CENT);
dxdt_MDZ_LIV1 = 5.0 / (Vh + Vi) * ((Qh - Qpv) * MDZ_CENT + Qpv * MDZ_PV - Qh * MDZ_LIV1 / KPLIV -
    fbCLinthMid * (1 + fm3Ah * (CYP3A_H1 - 1)) / 5 * MDZ_LIV1 / KPLIV);
dxdt_MDZ_LIV2 = 5.0 / (Vh + Vi) * (Qh * (MDZ_LIV1 - MDZ_LIV2) - fbCLinthMid * (1 + fm3Ah * (CYP3A_H2 - 1)) / 5 * MDZ_LIV2) / KPLIV;
dxdt_MDZ_LIV3 = 5.0 / (Vh + Vi) * (Qh * (MDZ_LIV2 - MDZ_LIV3) - fbCLinthMid * (1 + fm3Ah * (CYP3A_H3 - 1)) / 5 * MDZ_LIV3) / KPLIV;
dxdt_MDZ_LIV4 = 5.0 / (Vh + Vi) * (Qh * (MDZ_LIV3 - MDZ_LIV4) - fbCLinthMid * (1 + fm3Ah * (CYP3A_H4 - 1)) / 5 * MDZ_LIV4) / KPLIV;
dxdt_MDZ_LIV5 = 5.0 / (Vh + Vi) * (Qh * (MDZ_LIV4 - MDZ_LIV5) - fbCLinthMid * (1 + fm3Ah * (CYP3A_H5 - 1)) / 5 * MDZ_LIV5) / KPLIV;
dxdt_MDZ_MUS = 1.0 / Vm * Qm * (MDZ_CENT - MDZ_MUS / (KpSFMid * KpmMid));
dxdt_MDZ_SKIN = 1.0 / Vs * Qs * (MDZ_CENT - MDZ_SKIN / (KpSFMid * KpsMid));
dxdt_MDZ_ADI = 1.0 / Va * Qa * (MDZ_CENT - MDZ_ADI / (KpSFMid * KpaMid));
dxdt_MDZ_PV = 1.0 / Vpv * (Qpv * (MDZ_CENT - MDZ_PV) +
    KaMid * QgutMid / (QgutMid + fgCLintgMid * (1 + fm3Ag * (CYP3A_G - 1))) * MDZ_GUT);
dxdt_MDZ_GUT = -KaMid / FaMid * MDZ_GUT;

// ---- rifampicin ----
double UP1 = fbRif * (VmaxOATP / (KmuOATP + fbRif * RIF_HE1) + PSdifinf) * RIF_HE1;
double UP2 = fbRif * (VmaxOATP / (KmuOATP + fbRif * RIF_HE2) + PSdifinf) * RIF_HE2;
double UP3 = fbRif * (VmaxOATP / (KmuOATP + fbRif * RIF_HE3) + PSdifinf) * RIF_HE3;
double UP4 = fbRif * (VmaxOATP / (KmuOATP + fbRif * RIF_HE4) + PSdifinf) * RIF_HE4;
double UP5 = fbRif * (VmaxOATP / (KmuOATP + fbRif * RIF_HE5) + PSdifinf) * RIF_HE5;

dxdt_RIF_CENT = 1.0 / VbRif * (Qh * RIF_HE5 + Qm * (RIF_MUS / (KpScaleRif * KpmRif) - RIF_CENT) +
    Qs * (RIF_SKIN / (KpScaleRif * KpsRif) - RIF_CENT) + Qa * (RIF_ADI / (KpScaleRif * KpaRif) - RIF_CENT) -
    (Qhart + Qser + Qvilli) * RIF_CENT - CLrRif * RIF_CENT);
dxdt_RIF_HE1 = 5.0 / Vi * (Qhart * RIF_CENT + Qvilli * RIF_MUCB + Qser * RIF_SER / (KpScaleRif * KpserRif) -
    Qh * RIF_HE1 + (fhRif * PSdifeff * RIF_HC1 - UP1) / 5);
dxdt_RIF_HC1 = 5.0 / Vh * (UP1 - fhRif * PSdifeff * RIF_HC1 - fhRif * CLintRif * (1 + fmUGTh * (UGT_H1 - 1)) * RIF_HC1) / 5;
dxdt_RIF_HE2 = 5.0 / Vi * (Qh * (RIF_HE1 - RIF_HE2) + (fhRif * PSdifeff * RIF_HC2 - UP2) / 5);
dxdt_RIF_HC2 = 5.0 / Vh * (UP2 - fhRif * PSdifeff * RIF_HC2 - fhRif * CLintRif * (1 + fmUGTh * (UGT_H2 - 1)) * RIF_HC2) / 5;
dxdt_RIF_HE3 = 5.0 / Vi * (Qh * (RIF_HE2 - RIF_HE3) + (fhRif * PSdifeff * RIF_HC3 - UP3) / 5);
dxdt_RIF_HC3 = 5.0 / Vh * (UP3 - fhRif * PSdifeff * RIF_HC3 - fhRif * CLintRif * (1 + fmUGTh * (UGT_H3 - 1)) * RIF_HC3) / 5;
dxdt_RIF_HE4 = 5.0 / Vi * (Qh * (RIF_HE3 - RIF_HE4) + (fhRif * PSdifeff * RIF_HC4 - UP4) / 5);
dxdt_RIF_HC4 = 5.0 / Vh * (UP4 - fhRif * PSdifeff * RIF_HC4 - fhRif * CLintRif * (1 + fmUGTh * (UGT_H4 - 1)) * RIF_HC4) / 5;
dxdt_RIF_HE5 = 5.0 / Vi * (Qh * (RIF_HE4 - RIF_HE5) + (fhRif * PSdifeff * RIF_HC5 - UP5) / 5);
dxdt_RIF_HC5 = 5.0 / Vh * (UP5 - fhRif * PSdifeff * RIF_HC5 - fhRif * CLintRif * (1 + fmUGTh * (UGT_H5 - 1)) * RIF_HC5) / 5;
dxdt_RIF_SER = 1.0 / Vser * Qser * (RIF_CENT - RIF_SER / (KpScaleRif * KpserRif));
dxdt_RIF_MUS = 1.0 / Vm * Qm * (RIF_CENT - RIF_MUS / (KpScaleRif * KpmRif));
dxdt_RIF_SKIN = 1.0 / Vs * Qs * (RIF_CENT - RIF_SKIN / (KpScaleRif * KpsRif));
dxdt_RIF_ADI = 1.0 / Va * Qa * (RIF_CENT - RIF_ADI / (KpScaleRif * KpaRif));
dxdt_RIF_LUMEN = -KaRif / FaRif * RIF_LUMEN + fgRif * PSdifent * 20 * RIF_ENT;
dxdt_RIF_ENT = 1.0 / Vent * (KaRif * RIF_LUMEN + fbRif * PSdifent * RIF_MUCB -
    fgRif * (PSdifent * 21 + CLmetgRif * (1 + fmUGTg * (UGT_G - 1))) * RIF_ENT);
dxdt_RIF_MUCB = 1.0 / Vmucb * (Qvilli * (RIF_CENT - RIF_MUCB) + fgRif * PSdifent * RIF_ENT - fbRif * PSdifent * RIF_MUCB);

// ---- enzyme induction by unbound rifampicin ----
double IH1 = fhRif * RIF_HC1 / (fhRif * RIF_HC1 + EC50u);
double IH2 = fhRif * RIF_HC2 / (fhRif * RIF_HC2 + EC50u);
double IH3 = fhRif * RIF_HC3 / (fhRif * RIF_HC3 + EC50u);
double IH4 = fhRif * RIF_HC4 / (fhRif * RIF_HC4 + EC50u);
double IH5 = fhRif * RIF_HC5 / (fhRif * RIF_HC5 + EC50u);
double IG  = fgRif * RIF_ENT / (fgRif * RIF_ENT + EC50u);
dxdt_UGT_H1 = kdegUGTh * (1 + EmaxUGT * IH1 - UGT_H1);
dxdt_UGT_H2 = kdegUGTh * (1 + EmaxUGT * IH2 - UGT_H2);
dxdt_UGT_H3 = kdegUGTh * (1 + EmaxUGT * IH3 - UGT_H3);
dxdt_UGT_H4 = kdegUGTh * (1 + EmaxUGT * IH4 - UGT_H4);
dxdt_UGT_H5 = kdegUGTh * (1 + EmaxUGT * IH5 - UGT_H5);
dxdt_UGT_G  = kdegUGTg * (1 + EmaxUGT * IG - UGT_G);
dxdt_CYP3A_H1 = kdeg3Ah * (1 + Emax3A * IH1 - CYP3A_H1);
dxdt_CYP3A_H2 = kdeg3Ah * (1 + Emax3A * IH2 - CYP3A_H2);
dxdt_CYP3A_H3 = kdeg3Ah * (1 + Emax3A * IH3 - CYP3A_H3);
dxdt_CYP3A_H4 = kdeg3Ah * (1 + Emax3A * IH4 - CYP3A_H4);
dxdt_CYP3A_H5 = kdeg3Ah * (1 + Emax3A * IH5 - CYP3A_H5);
dxdt_CYP3A_G  = kdeg3Ag * (1 + Emax3A * IG - CYP3A_G);

$TABLE
double MDZ = MDZ_CENT * 1000;       // midazolam blood (ng/mL)
double RIF = RIF_CENT;              // rifampicin blood (ug/mL)
double CYP3A_LIVER = (CYP3A_H1 + CYP3A_H2 + CYP3A_H3 + CYP3A_H4 + CYP3A_H5) / 5;
double CYP3A_GUT = CYP3A_G;

$CAPTURE MDZ RIF CYP3A_LIVER CYP3A_GUT
