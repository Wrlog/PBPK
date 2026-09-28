"""Build models/ddi_rifampicin_glibenclamide.cpp.

The rifampicin part (physiology, liver, gut, UGT and CYP3A induction) is
taken from models/ddi_rifampicin_midazolam.cpp so the two DDI models share
one rifampicin; CYP2C9 induction and the glibenclamide model are added,
transcribed from the glibenclamide Napp code in the Supplementary Text of
Asaumi et al. 2018 (CPT Pharmacometrics Syst Pharmacol 7:186).

Usage: python tools/build_glibenclamide_ddi.py
"""

SRC = "models/ddi_rifampicin_midazolam.cpp"
DST = "models/ddi_rifampicin_glibenclamide.cpp"

src = open(SRC).read()
head, rest = src.split("$PARAM @annotated", 1)
params, rest = rest.split("$CMT @annotated", 1)
cmts, rest = rest.split("$MAIN", 1)
main, rest = rest.split("$ODE", 1)
ode, rest = rest.split("$TABLE", 1)

params = params.split("// midazolam (Table S2; optimised values from Table 1)")[0]
params = params.replace(
    "kdegUGTg    : 0.0288  : UGT degradation, enterocytes (1/h)\n",
    "kdegUGTg    : 0.0288  : UGT degradation, enterocytes (1/h)\n"
    "kdeg2C9h    : 0.00666 : CYP2C9 degradation, hepatocytes (1/h)\n")
params = params.replace(
    "Emax3A      : 4.57    : Emax, CYP3A induction [Table 1]\n",
    "Emax3A      : 4.57    : Emax, CYP3A induction [Table 1]\n"
    "Emax2C9     : 2.41    : Emax, CYP2C9 induction [Table 1]\n"
    "KiuOATP     : 0.186   : Unbound Ki of rifampicin for OATP uptake (ug/mL; 0.226 uM)\n")
params += """// glibenclamide (Table S2; optimised values from Table 1, beta = 0.2)
KaGlb       : 0.445   : Absorption rate (1/h) [Table 1]
FaGlb       : 1       : Fraction absorbed
fbGlb       : 0.000774 : Unbound fraction in blood
fhGlb       : 0.0221  : Unbound fraction in hepatocytes
KpSFGlb     : 0.568   : Scaling factor on in-silico Kp [Table 1, beta 0.2]
KpmGlb      : 0.104   : Kp muscle
KpsGlb      : 0.447   : Kp skin
KpaGlb      : 0.0795  : Kp adipose
fbCLintallGlb_kg : 0.123 : fB x overall hepatic intrinsic clearance (L/h/kg) [Table 1]
RdifGlb     : 0.246   : PSdif,inf / PSact,inf
betaGlb     : 0.2     : CLint,all / (PSact,inf + PSdif,inf) (the paper's preferred value)
gammaGlb    : 0.24    : PSdif,inf / PSdif,eff
fmh2C9Glb   : 0.85    : Fraction metabolised by CYP2C9, liver (with beta 0.2)
fmh3AGlb    : 0.15    : Fraction metabolised by CYP3A, liver (with beta 0.2)
fmg3AGlb    : 1       : Fraction metabolised by CYP3A, gut
fgCLintgGlb_kg : 0.0131 : fE x enterocyte intrinsic clearance (L/h/kg)
CLpermGlb_kg : 0.103  : Gut permeability clearance, Qgut model (L/h/kg)
CLrGlb_kg   : 0       : Renal clearance (L/h/kg)

"""

glb_cmts = [
    "GLB_CENT  : Glibenclamide blood (mg/L)",
]
for i in range(1, 6):
    glb_cmts.append("GLB_HE%d   : Glibenclamide hepatic extracellular %d (mg/L)" % (i, i))
    glb_cmts.append("GLB_HC%d   : Glibenclamide hepatocytes %d (mg/L)" % (i, i))
glb_cmts += [
    "GLB_MUS   : Glibenclamide muscle (mg/L)",
    "GLB_SKIN  : Glibenclamide skin (mg/L)",
    "GLB_ADI   : Glibenclamide adipose (mg/L)",
    "GLB_PV    : Glibenclamide portal vein (mg/L)",
    "GLB_GUT   : Glibenclamide gut lumen (mg)",
]
rif_cmts = [l for l in cmts.strip("\n").splitlines() if l.strip() and not l.startswith("MDZ_")]
cyp2c9 = ["CYP2C9_H%d  : CYP2C9 activity ratio, hepatocytes %d" % (i, i) for i in range(1, 6)]
cmts = "\n".join(glb_cmts + rif_cmts + cyp2c9) + "\n\n"

main = main.split("// midazolam")[0].replace("double VbMid = Vc_kg_mid * BW;\n", "") + """double VbGlb = Vb_kg * BW;
double CLintallGlb = fbCLintallGlb_kg * BW / fbGlb;
double PSactinfGlb = 1.0 / (1 + RdifGlb) * CLintallGlb / betaGlb;
double PSdifinfGlb = RdifGlb / (1 + RdifGlb) * CLintallGlb / betaGlb;
double PSdifeffGlb = RdifGlb / (1 + RdifGlb) * CLintallGlb / betaGlb / gammaGlb;
double CLinthGlb = RdifGlb / (1 + RdifGlb) * CLintallGlb / (1 - betaGlb) / gammaGlb;
double CLpermGlb = CLpermGlb_kg * BW;
double QgutGlb = Qvilli * CLpermGlb / (Qvilli + CLpermGlb);
double fgCLintgGlb = fgCLintgGlb_kg * BW;
double CLrGlb = CLrGlb_kg * BW;

UGT_H1_0 = 1; UGT_H2_0 = 1; UGT_H3_0 = 1; UGT_H4_0 = 1; UGT_H5_0 = 1; UGT_G_0 = 1;
CYP3A_H1_0 = 1; CYP3A_H2_0 = 1; CYP3A_H3_0 = 1; CYP3A_H4_0 = 1; CYP3A_H5_0 = 1; CYP3A_G_0 = 1;
CYP2C9_H1_0 = 1; CYP2C9_H2_0 = 1; CYP2C9_H3_0 = 1; CYP2C9_H4_0 = 1; CYP2C9_H5_0 = 1;

"""

rif_ode = "// ---- rifampicin ----" + ode.split("// ---- rifampicin ----")[1]
rif_ode = rif_ode.replace(
    "dxdt_CYP3A_G  = kdeg3Ag * (1 + Emax3A * IG - CYP3A_G);",
    "dxdt_CYP3A_G  = kdeg3Ag * (1 + Emax3A * IG - CYP3A_G);\n" +
    "\n".join("dxdt_CYP2C9_H%d = kdeg2C9h * (1 + Emax2C9 * IH%d - CYP2C9_H%d);" % (i, i, i) for i in range(1, 6)))

glb = ["// ---- glibenclamide: OATP uptake competitively inhibited by unbound rifampicin ----"]
for i in range(1, 6):
    glb.append("double UPG%d = fbGlb * (PSactinfGlb / (1 + fbRif * RIF_HE%d / KiuOATP) + PSdifinfGlb) * GLB_HE%d;" % (i, i, i))
glb.append("dxdt_GLB_CENT = 1.0 / VbGlb * (Qh * GLB_HE5 - (Qh - Qpv) * GLB_CENT - Qpv * GLB_CENT +\n"
           "    Qm * (GLB_MUS / (KpSFGlb * KpmGlb) - GLB_CENT) + Qs * (GLB_SKIN / (KpSFGlb * KpsGlb) - GLB_CENT) +\n"
           "    Qa * (GLB_ADI / (KpSFGlb * KpaGlb) - GLB_CENT) - CLrGlb * GLB_CENT);")
glb.append("dxdt_GLB_HE1 = 5.0 / Vi * ((Qh - Qpv) * GLB_CENT + Qpv * GLB_PV - Qh * GLB_HE1 + (fhGlb * PSdifeffGlb * GLB_HC1 - UPG1) / 5);")
for i in range(1, 6):
    if i > 1:
        glb.append("dxdt_GLB_HE%d = 5.0 / Vi * (Qh * (GLB_HE%d - GLB_HE%d) + (fhGlb * PSdifeffGlb * GLB_HC%d - UPG%d) / 5);"
                   % (i, i - 1, i, i, i))
    glb.append("dxdt_GLB_HC%d = 5.0 / Vh * (UPG%d - fhGlb * PSdifeffGlb * GLB_HC%d - fhGlb * CLinthGlb * "
               "(1 + fmh3AGlb * (CYP3A_H%d - 1) + fmh2C9Glb * (CYP2C9_H%d - 1)) * GLB_HC%d) / 5;" % (i, i, i, i, i, i))
glb += [
    "dxdt_GLB_MUS = 1.0 / Vm * Qm * (GLB_CENT - GLB_MUS / (KpSFGlb * KpmGlb));",
    "dxdt_GLB_SKIN = 1.0 / Vs * Qs * (GLB_CENT - GLB_SKIN / (KpSFGlb * KpsGlb));",
    "dxdt_GLB_ADI = 1.0 / Va * Qa * (GLB_CENT - GLB_ADI / (KpSFGlb * KpaGlb));",
    "dxdt_GLB_PV = 1.0 / Vpv * (Qpv * (GLB_CENT - GLB_PV) +\n"
    "    KaGlb * QgutGlb / (QgutGlb + fgCLintgGlb * (1 + fmg3AGlb * (CYP3A_G - 1))) * GLB_GUT);",
    "dxdt_GLB_GUT = -KaGlb / FaGlb * GLB_GUT;",
    "",
]

hdr = """$PROB
Rifampicin -> glibenclamide: OATP inhibition plus CYP2C9/3A induction

// Asaumi R, Toshimoto K, Tobe Y, et al. Comprehensive PBPK model of
// rifampicin for quantitative prediction of complex drug-drug
// interactions: CYP3A/2C9 induction and OATP inhibition effects.
// CPT Pharmacometrics Syst Pharmacol 2018;7:186-196 (open access).
//
// Generated by tools/build_glibenclamide_ddi.py. Glibenclamide is
// transcribed from the glibenclamide model code in the article's
// Supplementary Text, with values from Supplementary Tables S1-S2, Table 1
// (beta = 0.2 column) and the OATP inhibition constant quoted in the
// Results (Ki,u 0.226 uM, from Yoshikado et al. 2016).
//
// Glibenclamide enters hepatocytes by OATP (saturable in rifampicin's
// case; for glibenclamide linear, inhibited by unbound rifampicin in the
// same liver zone) and by passive diffusion, and is metabolised by CYP2C9
// and CYP3A. Rifampicin does two opposite things: acutely it blocks uptake
// (exposure up); over days it induces CYP2C9/3A (exposure down). The
// rifampicin part is the same as models/ddi_rifampicin_midazolam.cpp, with
// CYP2C9 induction added.
//
// Concentrations are total blood concentrations in mg/L; gut lumen states
// are amounts (mg). Time in hours. IV doses go into the blood
// concentration state (divide mg by the blood volume).

"""

tbl = """$TABLE
double GLB = GLB_CENT * 1000;       // glibenclamide blood (ng/mL)
double RIF = RIF_CENT;              // rifampicin blood (ug/mL)
double CYP3A_LIVER = (CYP3A_H1 + CYP3A_H2 + CYP3A_H3 + CYP3A_H4 + CYP3A_H5) / 5;
double CYP2C9_LIVER = (CYP2C9_H1 + CYP2C9_H2 + CYP2C9_H3 + CYP2C9_H4 + CYP2C9_H5) / 5;
double OATP_ACTIVITY = 1.0 / (1 + fbRif * RIF_HE1 / KiuOATP);

$CAPTURE GLB RIF CYP3A_LIVER CYP2C9_LIVER OATP_ACTIVITY
"""

out = (hdr + "$PARAM @annotated" + params + "$CMT @annotated\n" + cmts + "$MAIN" + main +
       "$ODE\n" + "\n".join(glb) + "\n" + rif_ode.rstrip() + "\n\n" + tbl)
open(DST, "w").write(out)
print("wrote", DST)
