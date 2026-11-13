#!/usr/bin/env python3
"""Generate ../aquarium233.js from Bert Chan's Lenia "aquarium" collection (233.json).

Usage: python3 tools/gen_aquarium.py [path/to/Lenia/Python/found/233.json]

Takes the entries from 5N7KKM to HAESRE (file lines 2-76: 5N7KKM to GDNQYX, then HGMLH8, PJ1EF7,
KZ2BLD, XEH4YR, HAESRE), skips "G6G6CR metamorpha beaked" and "VT049W gyrotachus"
(its own kernels: a different species), renames VT049W's first "para" to "paralimus", KZ2BLD's
unnamed second entry to "virgatus" and its "gyrans tachus" to "gyrans", HGMLH8's "parapedes" to
"tardus" and "para" to "tachus", numbers a repeated name
within a code ("auctus", "auctus 2"), and writes one pattern per entry: id "CODE-name", display
name "CODE name". The self-replicating forms ("R" in the file, and XEH4YR auctus) get " 💥" after
their display name; their id keeps "-replicating" (e.g. VT049W-gyrans-replicating), so URLs and
the aliases in core.js stay valid.
Kernels are stored once per rule code; an entry whose kernels differ from its code's first
entry would carry its own (none does now).

tools/extra_patterns.json adds patterns found in Cognoship itself (code, name, R, cells), each
placed after the last entry of its code and using that code's rule.
"""
import json, re, sys, os

SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/github/Lenia/Python/found/233.json')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'aquarium233.js')
# Left out: one form, and VT049W gyrotachus, which is a different species: its kernel 0 has
# m = 0.302 instead of VT049W's 0.272.
SKIP = {('G6G6CR', 'metamorpha beaked'), ('VT049W', 'gyrotachus')}
SKIP_IDS = {'XEH4YR-auctus-2', 'HAESRE-zigzag-2'}   # repeated names, left out (after numbering)
Z18A9R_NAMES = ['gyrans', 'paragyrans', 'antigyrans']   # unnamed in the file (lines 32-34)
KEYS = ['b', 'm', 's', 'h', 'r', 'c']

lines = open(SRC).read().split('\n')
entries = [json.loads(l.rstrip().rstrip(',')) for l in lines[1:76] if l.strip()]
# Renamed after the steps above (so VT049W's first "para" still becomes "paralimus")
POST_RENAME = {('G6G6CR', 'fortified'): 'munitus', ('VT049W', 'para'): 'paragyrans',
               ('KH97WU', 'antipara'): 'antigyrans', ('KH97WU', 'dual'): 'impetus+gyrans',
               ('HAESRE', 'anti'): 'antigyrans', ('HAESRE', 'para'): 'paragyrans',
               ('VT049W', 'rectus+gyrans'): 'pararectus'}
BOOM = {('XEH4YR', 'auctus')}   # self-replicating, though not marked R in the file
RENAME = {('KZ2BLD', ''): 'virgatus', ('KZ2BLD', 'gyrans tachus'): 'gyrans',
          ('HGMLH8', 'parapedes'): 'tardus', ('HGMLH8', 'para'): 'tachus'}

# Object detection per rule (Cognoship creatures.js): merge saddle, chosen by hand in tune mode
# (?tune=1). 1.0 = distinct potential peaks never merge; lower = bumps of one creature merge.
DETECT = {'5N7KKM': {'saddle': 0.95}, 'VT049W': {'saddle': 0.95}, 'G6G6CR': {'saddle': 0.9}, 'Z18A9R': {'saddle': 0.9},
          'KZ2BLD': {'saddle': 0.9}, 'XEH4YR': {'saddle': 0.9}}
DETECT_DEFAULT = {'saddle': 1.0}

def clean(code, name, k):
    name = name.lstrip('*').strip()
    if name.startswith('M '): name = name[2:]          # "**M limus" (5N7KKM)
    name = name.replace('OMIre', '').strip()
    if code == 'Z18A9R': name = Z18A9R_NAMES[k]
    name = re.sub(r'\s+R$', ' (replicating)', name)   # "vagus R" -> "vagus (replicating)"
    return RENAME.get((code, name), name)

# Changes to a code's rule (tools/rule_tweaks.json, shared with tools/from_screenshot.py): one
# kernel parameter set to a new value, e.g. Z18A9R kernel 3's m, so transgyrans is stable
TWEAKS = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rule_tweaks.json')))

rules, pats, seen, zk = {}, [], {}, 0
for e in entries:
    code = e['code'].split('<')[0]
    k = zk if code == 'Z18A9R' else 0
    if code == 'Z18A9R': zk += 1
    name = clean(code, e['name'], k)
    if (code, name) in SKIP: continue
    if (code, name) == ('VT049W', 'para') and 'VT049W-paralimus' not in seen:
        name = 'paralimus'                              # the first of VT049W's two "para" entries
    name = POST_RENAME.get((code, name), name)
    p0 = e['params'][0]
    ks = [{key: p[key] for key in KEYS} for p in e['params']]
    for t in TWEAKS:
        if t['code'] == code:
            assert ks[t['kernel']][t['key']] == t['was'], (code, t)
            ks[t['kernel']][t['key']] = t['value']
    assert all(p['kn'] == 1 and p['gn'] == 1 for p in e['params']), e['name']
    base = re.sub(r'[ ()]+', '-', name).strip('-')
    pid = f'{code}-{base}'
    if pid in seen:                                     # a repeated name (different cells): number it
        n = 2
        while f'{pid}-{n}' in seen: n += 1
        pid, name = f'{pid}-{n}', f'{name} {n}'
    seen[pid] = 1
    if pid in SKIP_IDS: continue
    if name.endswith(' (replicating)') or (code, name) in BOOM:   # display: an emoji in place
        name = name.replace(' (replicating)', '') + ' 💥'
    pat = {'id': pid, 'name': f'{code} {name}', 'rule': code, 'R': p0['R'], 'cells': e['cells']}
    if code not in rules: rules[code] = {'T': p0['T'], 'kernels': ks}
    elif ks != rules[code]['kernels']: pat['kernels'] = ks
    pats.append(pat)

EXTRA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'extra_patterns.json')
for x in json.load(open(EXTRA)) if os.path.exists(EXTRA) else []:
    assert x['code'] in rules, x['code']
    pid = f"{x['code']}-{re.sub(r'[ ()]+', '-', x['name']).strip('-')}"
    assert pid not in seen, pid
    seen[pid] = 1
    last = max(i for i, p in enumerate(pats) if p['rule'] == x['code'])
    pats.insert(last + 1, {'id': pid, 'name': f"{x['code']} {x['name']}", 'rule': x['code'], 'R': x['R'], 'cells': x['cells']})

def kjs(k): return '{' + ', '.join(f'{key}:{json.dumps(k[key])}' for key in KEYS) + '}'
out = ['"use strict";',
       '// GENERATED by tools/gen_aquarium.py from Lenia/Python/found/233.json -- do not edit by hand.',
       '// Rules (T, detect = object-detection settings, kernels) per code, and one pattern per entry (R = native radius, cells = RLE per channel).',
       'const AQUARIUM233 = {', '  rules: {']
for code, r in rules.items():
    det = ', '.join(f'{k}: {v}' for k, v in DETECT.get(code, DETECT_DEFAULT).items())
    out.append(f'    {json.dumps(code)}: {{ T: {r["T"]}, detect: {{ {det} }}, kernels: [')
    out += [f'      {kjs(k)},' for k in r['kernels']]
    out.append('    ] },')
out += ['  },', '  patterns: [']
for p in pats:
    extra = ''
    if 'kernels' in p: extra = ', kernels: [\n' + ''.join(f'      {kjs(k)},\n' for k in p['kernels']) + '    ]'
    out.append(f'    {{ id: {json.dumps(p["id"])}, name: {json.dumps(p["name"])}, rule: {json.dumps(p["rule"])}, R: {p["R"]}{extra},')
    out.append(f'      cells: {json.dumps(p["cells"])} }},')
out += ['  ],', '};', "if (typeof module !== 'undefined') module.exports = AQUARIUM233;", '']
open(OUT, 'w').write('\n'.join(out))
for p in pats: print(p['id'], '|', p['name'], '| R', p['R'], '| own kernels' if 'kernels' in p else '')
print(len(pats), 'patterns,', len(rules), 'rules ->', os.path.normpath(OUT))
