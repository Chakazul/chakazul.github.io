#!/usr/bin/env python3
"""Add a pattern captured in Cognoship (',' the ship, '.' the whole board, in the lab or the game:
core.js capture) to tools/extra_patterns.json, and regenerate aquarium233.js.

A capture holds the exact cell values, so unlike tools/from_screenshot.py there is nothing to read
back or settle: it is the pattern as it was. Usage (from Cognoship/):
  python3 tools/add_pattern.py "~/Downloads/Z18A9R unnamed 2026-11-18-10-20-30.json" --name trigyrans [--check]
  --name   the morph's name (default: the one in the file, the detected morph or "unnamed")
  --code   the species code (default: the file's)
  --check  first run it --steps steps (500) under the code's rule on Cognoship's 384x256 torus (as
           from_screenshot.py's verify: needs numpy, scipy, e.g. the complexity venv) and report
           the number of creatures and the mass drift; basic, on purpose: the lab is where it is
           tried out (see from_screenshot.py's docstring)
  --no-gen don't run tools/gen_aquarium.py after
An entry of the same code and name is replaced."""
import argparse, json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def check(entry, steps, src):
    import numpy as np
    sys.path.insert(0, HERE)
    import from_screenshot as fs
    sys.path.insert(0, os.path.join(HERE, '..', 'output', 'transgyrans', 'tune'))
    rows = []
    for c in entry['cells']:   # decode the RLE (as aquarium233.js does)
        rows.append(rle(c))
    h, w = max(len(r) for r in rows), max(max((len(x) for x in r), default=0) for r in rows)
    A0 = np.zeros((3, h, w))
    for k, r in enumerate(rows):
        for y, row in enumerate(r): A0[k, y, :len(row)] = row
    if h > fs.H or w > fs.W: sys.exit(f'too large for the board: {w}x{h}')
    ps, T = fs.load_rule(src, entry['code'])
    fs.CFG['R'] = entry['R']
    Ks, dt = fs.kernels(ps, entry['R']), 1 / T
    A = fs.place(A0)
    m0, hist = A.sum(), []
    for t in range(steps + 1):
        if t % 100 == 0: hist.append((t, A.sum() / m0, fs.parts(A, ps, Ks)))
        if t < steps: A = fs.step(A, ps, Ks, dt)
    n = sorted({h[2] for h in hist})
    drift = max(abs(h[1] - 1) for h in hist)
    print(f'check: {steps} steps, creatures {n}, mass drift {drift:.2%} '
          + ('(holds)' if len(n) == 1 and 0 not in n and drift <= 0.015 else '(changes: look at it in the lab)'))


def rle(s):
    """Lenia RLE -> rows of 0..1 values"""
    rows, row, n, i = [], [], '', 0
    while i < len(s):
        c = s[i]
        if c.isdigit(): n += c; i += 1; continue
        k = int(n) if n else 1; n = ''
        if c in '$!':
            rows.append(row); row = []
            if c == '$': rows += [[] for _ in range(k - 1)]
            i += 1; continue
        if c == '.': v = 0; i += 1
        elif 'p' <= c <= 'y': v = 25 + (ord(c) - ord('p')) * 24 + ord(s[i + 1]) - ord('A'); i += 2
        else: v = ord(c) - ord('A') + 1; i += 1
        row += [v / 255] * k
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('json', help="a capture (',' or '.' in Cognoship)")
    ap.add_argument('--name'); ap.add_argument('--code')
    ap.add_argument('--check', action='store_true'); ap.add_argument('--steps', type=int, default=500)
    ap.add_argument('--no-gen', action='store_true')
    ap.add_argument('--src', default=os.path.expanduser('~/github/Lenia/Python/found/233.json'))
    a = ap.parse_args()
    e = json.load(open(os.path.expanduser(a.json)))
    entry = {'code': a.code or e['code'], 'name': a.name or e['name'], 'R': e['R'], 'source': e.get('source', ''), 'cells': e['cells']}
    if entry['name'] in ('unnamed', 'board') and not a.name:
        print(f'note: the name is "{entry["name"]}"; give one with --name')
    if a.check: check(entry, a.steps, a.src)
    path = os.path.join(HERE, 'extra_patterns.json')
    extra = json.load(open(path)) if os.path.exists(path) else []
    extra = [x for x in extra if (x['code'], x['name']) != (entry['code'], entry['name'])] + [entry]
    json.dump(extra, open(path, 'w'), indent=1)
    print(f'added {entry["code"]} {entry["name"]} to {path}')
    if not a.no_gen:
        subprocess.run([sys.executable, os.path.join(HERE, 'gen_aquarium.py')], check=True)


if __name__ == '__main__':
    main()
