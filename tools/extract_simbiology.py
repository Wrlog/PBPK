"""Extract a SimBiology model from a .sbproj file, without MATLAB.

Usage: python tools/extract_simbiology.py <project.sbproj> <out.json>

A .sbproj file is a zip archive whose simbiodata.mat holds the model as a
MATLAB class object (SimBiology.Model). Recent MATLAB versions save class
objects in the "MCOS" layout: the variable is only a reference, and the
objects themselves live in a hidden byte array at the end of the MAT file.
That array is itself a small MAT stream holding one cell array:
  - cell 0: metadata - the property and class names, which class each object
    belongs to, and for each object its property list as (name, kind, value)
    triples. kind 0: value is a name index, kind 1: value indexes the cell
    array (offset by 2), kind 2: value is stored inline (logicals).
  - the other cells: the property values.
Object references inside property values are uint32 arrays starting with
0xDD000000, followed by the dimensions, the object numbers and the class.

This reads that layout and writes the species, parameters, rules, reactions
(with their Active flags), events, variants and doses as plain JSON. Used to
rebuild the generalized bispecific model of Spinosa et al. (CPT Pharmacometrics
Syst Pharmacol 2026;15:e70167) from the project in its Supporting Information.
"""
import json
import struct
import sys
import zipfile
import zlib

FMT = {1: 'b', 2: 'B', 3: 'h', 4: 'H', 5: 'i', 6: 'I', 7: 'f', 9: 'd', 12: 'q', 13: 'Q', 16: 'B', 17: 'H', 18: 'I'}
SZ = {'b': 1, 'B': 1, 'h': 2, 'H': 2, 'i': 4, 'I': 4, 'f': 4, 'd': 8, 'q': 8, 'Q': 8}
REF = 0xDD000000


def decompress(mat, start=128):
    out, p = b'', start
    while p < len(mat):
        t, n = struct.unpack_from('<II', mat, p)
        d = mat[p + 8:p + 8 + n]
        out += zlib.decompress(d) if t == 15 else struct.pack('<II', t, n) + d
        p += 8 + n
    return out


def tag(b, p):
    t, n = struct.unpack_from('<II', b, p)
    if t >> 16:                                   # small data element
        return t & 0xffff, t >> 16, p + 4, p + 8
    return t, n, p + 8, p + 8 + ((n + 7) // 8) * 8


def vals(b, t, n, dp):
    f = FMT.get(t, 'B')
    return list(struct.unpack_from('<%d%s' % (n // SZ[f], f), b, dp))


def matrix(b, p, end):
    """One miMATRIX element -> a Python value (cells become lists)."""
    if end - p < 8:
        return None
    t, n, dp, p = tag(b, p); flags = vals(b, t, n, dp)[0]; cls = flags & 0xff
    t, n, dp, p = tag(b, p)                                     # dimensions
    t, n, dp, p = tag(b, p)                                     # array name
    if cls == 1:                                                # cell
        items = []
        while p < end:
            t, n, dp, nxt = tag(b, p)
            items.append(matrix(b, dp, dp + n) if t == 14 else None)
            p = nxt
        return items
    if cls in (2, 3):                                           # struct / object
        if cls == 3:
            t, n, dp, p = tag(b, p)
        t, n, dp, p = tag(b, p); width = vals(b, t, n, dp)[0]
        t, n, dp, p = tag(b, p); raw = bytes(vals(b, t, n, dp))
        fields = [raw[i:i + width].split(b'\0')[0].decode() for i in range(0, len(raw), width)]
        out, i = {}, 0
        while p < end:
            t, n, dp, nxt = tag(b, p)
            if t == 14 and fields:
                out[fields[i % len(fields)]] = matrix(b, dp, dp + n)
            p = nxt; i += 1
        return out
    if cls == 17:                                               # opaque (MCOS wrapper)
        t, n, dp, p = tag(b, p)                                 # class name
        while p < end:
            t, n, dp, nxt = tag(b, p)
            if t == 14:
                return matrix(b, dp, dp + n)
            p = nxt
        return None
    if p >= end:
        return []
    t, n, dp, p = tag(b, p)
    v = vals(b, t, n, dp)
    if cls == 4:                                                # char
        return ''.join(chr(c) for c in v)
    return v if len(v) != 1 else v[0]


def objects(project):
    """All objects of the saved model: {number: (class, {property: value})}."""
    b = decompress(zipfile.ZipFile(project).read('simbiodata.mat'))
    p, elems = 0, []
    while p < len(b):
        t, n = struct.unpack_from('<II', b, p)
        elems.append((p + 8, n)); p += 8 + n
    dp, n = elems[-1]                                           # the hidden workspace
    for _ in range(3):
        t_, n_, d_, dp = tag(b, dp)
    t_, n_, d_, dp = tag(b, dp)
    inner = bytes(vals(b, t_, n_, d_))[8:]
    t, n = struct.unpack_from('<II', inner, 0)
    wrapper = matrix(inner, 8, 8 + n)
    cells = wrapper['MCOS']
    md = bytes(cells[0])
    ver, nnames = struct.unpack_from('<II', md, 0)
    off = struct.unpack_from('<8I', md, 8)
    p, names = 40, []
    while len(names) < nnames:
        e = md.index(b'\0', p); names.append(md[p:e].decode('latin1')); p = e + 1
    name = lambda i: names[i - 1] if i else ''
    classes = [name(struct.unpack_from('<4I', md, q)[1]) for q in range(off[0], off[1], 16)]
    objs = [struct.unpack_from('<6I', md, q) for q in range(off[2], off[3], 24)]
    p, blocks = off[3], []
    while p < off[4]:
        k = struct.unpack_from('<I', md, p)[0]; p += 4
        blocks.append([struct.unpack_from('<3I', md, p + 12 * j) for j in range(k)]); p += 12 * k
        p = (p + 7) // 8 * 8
    out = {}
    for pos, (c, _, _, _, blk, _) in enumerate(objs):
        if pos == 0:
            continue
        props = {}
        for a, kind, v in blocks[blk]:
            props[name(a)] = name(v) if kind == 0 else (cells[v + 2] if kind == 1 else v)
        out[pos] = (classes[c], props)
    return out


def refs(v):
    if isinstance(v, list) and len(v) >= 3 and v[0] == REF:
        nd = v[1]; count = 1
        for d in v[2:2 + nd]:
            count *= d
        return v[2 + nd:2 + nd + count]
    return []


def strip(x):
    return x.replace('minPBPK.', '') if isinstance(x, str) else x


def main(project, out):
    O = objects(project)
    model = next(p for c, p in O.values() if c == 'Model')
    comp = [O[i][1] for i in refs(model['Compartments'])]
    cfg = [O[i][1] for i in refs(model['ConfigSets'])]
    opts = [O[i][1] for c in cfg for i in refs(c.get('CompileOptions'))]
    doc = {
        'model': model['Name'],
        'compartments': [{'name': c['Name'], 'value': c['Value']} for c in comp],
        'compile_options': {k: opts[0].get(k) for k in ('DefaultSpeciesDimension', 'DimensionalAnalysis', 'UnitConversion')},
        'species': [], 'parameters': [], 'rules': [], 'reactions': [], 'events': [], 'variants': [], 'doses': [],
    }
    for c in comp:
        for i in refs(c['Species']):
            s = O[i][1]
            doc['species'].append({'name': s['Name'], 'value': s['Value'], 'units': s['Units'],
                                   'boundary': bool(s.get('BoundaryCondition')), 'constant': bool(s.get('Constant'))})
    for i in refs(model['Parameters']):
        s = O[i][1]
        doc['parameters'].append({'name': s['Name'], 'value': s['Value'], 'units': s['Units'],
                                  'constant': bool(s.get('Constant', 1))})
    for i in refs(model['Rules']):
        s = O[i][1]
        doc['rules'].append({'type': s['RuleType'], 'rule': strip(s['Rule']), 'active': bool(s.get('Active', 1))})
    for i in refs(model['Reactions']):
        s = O[i][1]
        law = refs(s.get('KineticLaw'))
        local = [O[j][1]['Name'] for k in law for j in refs(O[k][1].get('Parameters'))]
        assert not local, 'reaction-scoped parameters are not handled: %s' % local
        doc['reactions'].append({'name': s['Name'], 'reaction': s['Reaction'], 'rate': strip(s['ReactionRate']),
                                 'active': bool(s.get('Active', 1))})
    for i in refs(model['Events']):
        s = O[i][1]
        fcns = s['EventFcns']
        fcns = [fcns[k] for k in sorted(fcns, key=int)] if isinstance(fcns, dict) else [fcns]
        while len(fcns) == 1 and isinstance(fcns[0], list):
            fcns = fcns[0]
        doc['events'].append({'trigger': strip(s['Trigger']), 'functions': [strip(f) for f in fcns],
                              'active': bool(s.get('Active', 1))})
    for i in refs(model['Variants']):
        s = O[i][1]
        content = s['Content']
        items = [content[k] for k in sorted(content, key=lambda x: int(x))] if isinstance(content, dict) else content
        values = {x[1]: x[3] for x in items}
        doc['variants'].append({'name': s['Name'], 'active': bool(s.get('Active')), 'values': values})
    for i in refs(model['Doses']):
        s = O[i][1]
        doc['doses'].append({k: strip(s.get(k)) for k in ('Name', 'TargetName', 'Amount', 'Rate', 'Time',
                                                          'StartTime', 'Interval', 'RepeatCount') if k in s})
    json.dump(doc, open(out, 'w'), indent=1, default=float)
    print('%s: %d species, %d parameters, %d rules, %d reactions, %d events, %d variants, %d doses' % (
        doc['model'], len(doc['species']), len(doc['parameters']), len(doc['rules']), len(doc['reactions']),
        len(doc['events']), len(doc['variants']), len(doc['doses'])))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
