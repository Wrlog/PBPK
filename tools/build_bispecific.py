"""Build models/bispecific_mpbpk_spinosa2026.cpp from the SimBiology model.

Usage: python tools/build_bispecific.py
Input: tools/bispecific_spinosa2026_source.json, written by
tools/extract_simbiology.py from the project in the Supporting Information
(Data S1, PSP-2025-0197-s04.sbproj) of Spinosa et al., "A Generalized Minimal
PBPK-PD Model of Bispecific Antibodies", CPT Pharmacometrics Syst Pharmacol
2026;15:e70167 (open access).

The project has one SimBiology compartment of volume 1 with unit conversion
off, so every species' derivative is the plain sum of the rates of the
reactions it takes part in; the volume conversions between the minimal-PBPK
spaces (central, tight, leaky, lymph, efficacy, safety) are already written
into the rates. Rules carry over as they are: initial assignments go to
$MAIN, repeated assignments to $ODE (and $TABLE for outputs), rate rules
become their own compartments. The two MATLAB functions the rates call
(transport_between_compartments.m, binding_on_off.m, also in Data S1) are
expanded inline.

Default parameter values are the project's own with its active "nominal"
variants applied; the case studies of the paper (Model parameters.xlsx, in
tools/bispecific_spinosa2026_casestudies.csv) are applied by the app.

Left out: the dosing events. The project doses mg/kg or mg into bookkeeping
species that an event converts into a bolus of D1_cen; here the app gives
that bolus directly (nM in D1_cen, or mg in D1_ext_mg for subcutaneous).
Renamed: F (subcutaneous bioavailability) to F_sc.
"""
import ast
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = json.load(open(ROOT / 'tools' / 'bispecific_spinosa2026_source.json'))
OUT = ROOT / 'models' / 'bispecific_mpbpk_spinosa2026.cpp'

RENAME = {'F': 'F_sc'}
# species that only exist for the dosing events
EVENT_ONLY = {'D1_cen_mpk', 'D1_cen_mg'}
EVENT_PARAMS = {'dose_on_off', 'dose_time'}
OUTPUTS = ['D1_cen_ugml', 'total_D1_cen_ugml', 'D1_lea_ugml', 'D1_tig_ugml', 'D1_eff_ugml', 'D1_saf_ugml',
           'RO_mR1_cen', 'RO_mR2_cen', 'RO_mR1_lea', 'RO_mR2_lea', 'RO_mR1_tig', 'RO_mR2_tig',
           'RO_mR1_eff', 'RO_mR2_eff', 'TN_sR1_cen', 'TN_sR2_cen', 'TN_sR1_lea', 'TN_sR2_lea',
           'TN_sR1_tig', 'TN_sR2_tig', 'total_sR1_cen', 'total_sR2_cen', 'trimer_per_mR1_cell',
           'trimer_per_mR2_cell', 'CL_D1_tot_nmol', 'CL_D1_tot_cen_nmol']


def rn(x):
    return RENAME.get(x, x)


def side(s):
    out = {}
    for term in s.split('+'):
        term = term.strip()
        if not term or term == 'null':
            continue
        m = re.match(r'(\d+)\s+(\w+)$', term)
        k, name = (int(m.group(1)), m.group(2)) if m else (1, term)
        out[name] = out.get(name, 0) + k
    return out


def to_c(expr):
    """SimBiology expression -> C, with the two helper functions inlined."""
    expr = re.sub(r'\[([A-Za-z_]\w*)\]', r'\1', expr)     # SimBiology [quoted] names
    tree = ast.parse(expr.replace('^', '**'), mode='eval').body

    def emit(n):
        if isinstance(n, ast.BinOp):
            l, r = emit(n.left), emit(n.right)
            if isinstance(n.op, ast.Pow):
                return 'pow(%s, %s)' % (l, r)
            op = {ast.Add: '+', ast.Sub: '-', ast.Mult: '*', ast.Div: '/'}[type(n.op)]
            return '(%s %s %s)' % (l, op, r)
        if isinstance(n, ast.UnaryOp):
            return '(-%s)' % emit(n.operand) if isinstance(n.op, ast.USub) else emit(n.operand)
        if isinstance(n, ast.Constant):
            return repr(float(n.value))
        if isinstance(n, ast.Name):
            return rn(n.id)
        if isinstance(n, ast.Call):
            a = [emit(x) for x in n.args]
            if n.func.id == 'transport_between_compartments':
                # rate = C * L * (1 - sig) ./ V   (transport_between_compartments.m)
                C, sig, L, V = a
                return '(%s * %s * (1.0 - %s) / %s)' % (C, L, sig, V)
            if n.func.id == 'binding_on_off':
                # kon_obs = chi*kon; koff = kon*KD; rate = v*(kon_obs*D*R - koff*Complex)  (binding_on_off.m)
                D, R, Cx, kon, KD, v = a[:6]
                chi = a[6] if len(a) > 6 else '1.0'
                return '(%s * (%s * %s * %s * %s - %s * %s * %s))' % (v, chi, kon, D, R, kon, KD, Cx)
            f = {'max': 'fmax', 'min': 'fmin'}.get(n.func.id, n.func.id)
            return '%s(%s)' % (f, ', '.join(a))
        raise ValueError(ast.dump(n))
    return emit(tree)


def names(expr):
    return set(re.findall(r'(?<![0-9.])[A-Za-z_][A-Za-z0-9_]*', expr))


def parameters():
    p = {x['name']: x['value'] for x in SRC['parameters']}
    for v in SRC['variants']:
        if v['active']:
            p.update(v['values'])
    return p


def main():
    p = parameters()
    species = [s['name'] for s in SRC['species']]
    rules = {'initialAssignment': {}, 'repeatedAssignment': {}, 'rate': {}}
    for r in SRC['rules']:
        assert r['active']
        lhs, rhs = [x.strip() for x in r['rule'].split('=', 1)]
        rules[r['type']][lhs] = rhs
    ia, ra, rr = rules['initialAssignment'], rules['repeatedAssignment'], rules['rate']

    rx = []
    for r in SRC['reactions']:
        assert r['active']
        lhs, rhs = re.split(r'<->|->', r['reaction'])
        net = side(rhs)
        for k, v in side(lhs).items():
            net[k] = net.get(k, 0) - v
        net = {k: v for k, v in net.items() if v}
        assert not set(net) & set(ra), r
        rx.append((r['name'], net, r['rate']))

    states = [s for s in species if s not in ra and s not in EVENT_ONLY]
    derived = set(ia) - set(states)                 # parameters set by initial assignment
    for s in states:
        assert any(s in net for _, net, _ in rx) or s in rr, 'orphan species %s' % s

    # repeated assignments in dependency order
    order, done = [], set()
    while len(order) < len(ra):
        for n in sorted(set(ra) - done):
            if not (names(ra[n]) & (set(ra) - done - {n})):
                order.append(n); done.add(n)
    # initial assignments to parameters, in dependency order
    ia_order, done = [], set()
    while len(ia_order) < len(derived):
        for n in sorted(derived - done):
            if not (names(ia[n]) & (derived - done - {n})):
                ia_order.append(n); done.add(n)

    used = set()
    for _, _, rate in rx:
        used |= names(rate)
    for d in (ia, ra, rr):
        for v in d.values():
            used |= names(v)
    funcs = {'transport_between_compartments', 'binding_on_off', 'log', 'exp', 'time'}
    pars = sorted((used - set(species) - set(ra) - derived - funcs - EVENT_PARAMS) & set(p), key=str.lower)
    missing = used - set(p) - set(species) - set(ra) - derived - funcs
    assert not missing, missing

    L = ['$PROB\nGeneralized minimal PBPK model for bispecific antibodies: soluble, trans- and cis-binding targets\n']
    L.append('''// Spinosa P, Joslyn L, Ramanujan S, Gadkar K, Hosseini I. A Generalized
// Minimal PBPK-PD Model of Bispecific Antibodies: Case Studies and
// Applications in Drug Development. CPT Pharmacometrics Syst Pharmacol
// 2026;15:e70167 (open access).
//
// GENERATED by tools/build_bispecific.py from the authors' SimBiology project
// (Data S1, PSP-2025-0197-s04.sbproj). Do not edit by hand.
//
// A minimal PBPK model (Cao & Jusko 2014) with blood (cen), leaky (lea) and
// tight (tig) tissue, lymph (lym) and two optional tissues for efficacy (eff,
// e.g. a tumour) and safety (saf). In every space except lymph the drug D1
// binds two targets, R1 and R2, which may be soluble (sR) or on cells (mR).
// A second arm binds either the other target on another cell (trans, a T-cell
// engager), the other target on the same cell (cis, with an avidity factor),
// or the same target again (bivalent antibody). Targets are made at a constant
// rate and degraded first-order; complexes internalise at the target
// degradation rate. Drug-soluble target complexes clear like free drug and,
// optionally, traffic with it.
//
// Units: time day; drug and targets nM in each space (the SimBiology
// compartment has volume 1); volumes mL/kg; flows mL/day/kg. Doses: IV bolus
// of dose(mg/kg) / V_cen / MWab * 1e9 nM into D1_cen; subcutaneous doses in
// mg into D1_ext_mg (first-order absorption kabs, bioavailability F_sc).
// Defaults: the project's nominal variants (a linear IgG, no targets).
''')
    L.append('$PARAM')
    for k in pars:
        L.append('%s = %r' % (rn(k), float(p[k])))
    L.append('\n$CMT ' + ' '.join(states))
    L.append('\n$MAIN')
    for n in ia_order:
        L.append('double %s = %s;' % (n, to_c(ia[n])))
    for s in states:
        if s in ia:
            L.append('%s_0 = %s;' % (s, to_c(ia[s])))
    L.append('\n$ODE')
    for n in order:
        L.append('double %s = %s;' % (n, to_c(ra[n])))
    # each reaction rate once, then every species sums the rates it takes part in
    for i, (name, _, rate) in enumerate(rx):
        L.append('double v%d = %s;   // %s' % (i + 1, to_c(rate), name))
    for s in states:
        terms = []
        for i, (_, net, _) in enumerate(rx):
            if s in net:
                c = net[s]
                terms.append(('+ ' if c > 0 else '- ') + ('%d * ' % abs(c) if abs(c) != 1 else '') + 'v%d' % (i + 1))
        if s in rr:
            terms.append('+ ' + to_c(rr[s]))
        L.append('dxdt_%s = %s;' % (s, ' '.join(terms).lstrip('+ ') if terms else '0'))
    L.append('\n$TABLE')
    for n in order:
        L.append('%s = %s;' % (n, to_c(ra[n])))
    L.append('\n$CAPTURE ' + ' '.join(OUTPUTS))
    with open(OUT, 'w', newline='\n') as f:
        f.write('\n'.join(L) + '\n')
    print('%d states, %d reactions, %d rate rules, %d assignments, %d derived, %d parameters' % (
        len(states), len(rx), len(rr), len(order), len(ia_order), len(pars)))


if __name__ == '__main__':
    main()
