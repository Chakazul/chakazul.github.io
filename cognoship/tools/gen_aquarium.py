#!/usr/bin/env python3
"""Generate ../aquarium233.js from Bert Chan's Lenia "aquarium" collection (233.json).

Usage: python3 tools/gen_aquarium.py [path/to/Lenia/Python/found/233.json]

Takes the entries from 5N7KKM to HAESRE (file lines 2-76: 5N7KKM to GDNQYX, then HGMLH8, PJ1EF7,
KZ2BLD, XEH4YR, HAESRE), skips "G6G6CR metamorpha beaked" and "VT049W gyrotachus", gives HAESRE zigzag its second entry's cells (CELLS_FROM),
(its own kernels: a different species), renames VT049W's first "para" to "rectogyrans" (was "paralimus"), KZ2BLD's
unnamed second entry to "virgatus" and its "gyrans tachus" to "limus", HGMLH8's "parapedes" to
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
# A repeated entry whose cells replace an earlier one's (same params; the earlier one keeps its id,
# name and place): HAESRE zigzag's first entry starts mid-turn, hard to clone; its second (file
# line 76) is on a straight stretch (2026-11-19)
CELLS_FROM = {'HAESRE-zigzag-2': 'HAESRE-zigzag'}
Z18A9R_NAMES = ['gyrans', 'paragyrans', 'antigyrans']   # unnamed in the file (lines 32-34)
KEYS = ['b', 'm', 's', 'h', 'r', 'c']

lines = open(SRC).read().split('\n')
entries = [json.loads(l.rstrip().rstrip(',')) for l in lines[1:76] if l.strip()]
# Renamed after the steps above (VT049W's first "para" is named apart: vt_para)
POST_RENAME = {('G6G6CR', 'fortified'): 'munitus', ('VT049W', 'para'): 'paragyrans',
               ('KH97WU', 'antipara'): 'antigyrans', ('KH97WU', 'dual'): 'impetus+gyrans',
               ('HAESRE', 'anti'): 'antigyrans', ('HAESRE', 'para'): 'paragyrans',
               ('VT049W', 'rectus+gyrans'): 'sequerectus'}   # (was rectogyrans until 2026-11-18)
# A code's morphs in this order in the lab's picker (the rest after them, as they come)
MORPH_ORDER = {'VT049W': ['rectus', 'gyrans', 'tachus', 'lithos', 'paragyrans', 'pararectus', 'rectogyrans', 'sequerectus']}
BOOM = {('XEH4YR', 'auctus')}   # self-replicating, though not marked R in the file
RENAME = {('KZ2BLD', ''): 'virgatus', ('KZ2BLD', 'gyrans tachus'): 'limus',
          ('HGMLH8', 'parapedes'): 'paratardus', ('HGMLH8', 'para'): 'paratachus'}

# Object detection per rule (Cognoship creatures.js), chosen by eye in the lab's Tuning mode:
# saddle, the merge of one creature's bumps into a part (1.0 = distinct potential peaks never
# merge); wholeSaddle, the second merge of parts into wholes (twins, trios; 0.6 keeps the twins of
# the measurement in output/whole_saddle whole, GDNQYX paratardus included; a loose twin, gyrans
# orbiting a lithos, is whole only while its parts touch, at any value).
DETECT = {'VT049W': {'saddle': 1.0},   # (was 0.95: 1.0 keeps pararectus's two rectus apart, two parts of a whole)
          'KZ2BLD': {'saddle': 0.95}, 'G6G6CR': {'saddle': 0.9},
          'KH97WU': {'wholeSaddle': 0.4},
          'Z18A9R': {'wholeSaddle': 0.2}}   # (transgyrans: its gyrans part between crossings; whole 54 % of the time at 0.6, 99 % at 0.3 as
                                            # measured, ~2/3 in the lab after the kernel 4 tweak; always at 0.2. Touching
                                            # Z18A9R creatures already join at 0.3, 90 % of contacts; 99.9 % at 0.2)
DETECT_DEFAULT = {'saddle': 1.0, 'wholeSaddle': 0.6}

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
vt_para = False   # VT049W's first "para" named yet
for e in entries:
    code = e['code'].split('<')[0]
    k = zk if code == 'Z18A9R' else 0
    if code == 'Z18A9R': zk += 1
    name = clean(code, e['name'], k)
    if (code, name) in SKIP: continue
    if (code, name) == ('VT049W', 'para') and not vt_para:
        name, vt_para = 'rectogyrans', True             # the first of VT049W's two "para" entries (paralimus until 2026-11-18)
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
    if pid in CELLS_FROM:
        tgt = next(p for p in pats if p['id'] == CELLS_FROM[pid])
        assert tgt['R'] == p0['R'] and ks == (tgt.get('kernels') or rules[code]['kernels']), pid
        tgt['cells'] = e['cells']
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

for code, order in MORPH_ORDER.items():   # (stable: the others keep their places after them)
    idx = [i for i, p in enumerate(pats) if p['rule'] == code]
    base = lambda p: p['name'].split(' ', 1)[1].replace(' 💥', '')
    rank = lambda p: order.index(base(p)) if base(p) in order else len(order)
    for i, p in zip(idx, sorted((pats[i] for i in idx), key=rank)): pats[i] = p

def kjs(k): return '{' + ', '.join(f'{key}:{json.dumps(k[key])}' for key in KEYS) + '}'
out = ['"use strict";',
       '// GENERATED by tools/gen_aquarium.py from Lenia/Python/found/233.json -- do not edit by hand.',
       '// Rules (T, detect = object-detection settings, kernels) per code, and one pattern per entry (R = native radius, cells = RLE per channel).',
       'const AQUARIUM233 = {', '  rules: {']
# The morphs recorded as wholes (twins, trios) in a species' long history, assets/stats/<code>.json
# ("level": "whole"), for the game's species book (their thumbnails show more of the picture);
# rerun after recording new stats
STATS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'assets', 'stats')
def wholes(code):
    p = os.path.join(STATS, f'{code}.json')
    if not os.path.exists(p): return []
    return sorted({sg['morph'] for sg in json.load(open(p))['segs'] if sg.get('level') == 'whole'})

for code, r in rules.items():
    det = ', '.join(f'{k}: {v}' for k, v in {**DETECT_DEFAULT, **DETECT.get(code, {})}.items())
    out.append(f'    {json.dumps(code)}: {{ T: {r["T"]}, detect: {{ {det} }}, wholes: {json.dumps(wholes(code))}, kernels: [')
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
