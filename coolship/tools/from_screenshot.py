#!/usr/bin/env python3
"""Rebuild a Coolship pattern from a screenshot of the world view, and optionally add it to
tools/extra_patterns.json (then rerun tools/gen_aquarium.py).

Needs numpy, scipy and Pillow, e.g. the complexity venv:
  ~/github/complexity/.venv/bin/python tools/from_screenshot.py "../assets/screenshots/Z18A9R trio.png" \\
      --code Z18A9R --name trio [--add]

Steps (the ones used for Z18A9R trio):
  1. read the cell grid off the screenshot: the pixel period (screen px per cell) from the
     spectrum of the image's edges, its phase from where the edges fall, then one median value
     per cell, RGB = the three channels; dim grey cells (the zoom frame, drawn on top) and
     values below --floor are dropped, and the result is cropped;
  2. a screenshot is only close to a real state: its colours look brighter than the cell
     values (value ** 1/gamma, gamma ~ 1.4-2.2, from the display / screenshot colour handling),
     and small errors push the creature off its attractor. So try it over gammas (--gammas) and
     scales (--scales), and let each settle under the code's rule (LeniaF.py's, at --R, on
     Coolship's 384x256 torus) for --settle steps (300), then run --verify more steps (500).
     The check is short on purpose: a form may be unstable or turn into another in the long
     run and still be worth a pattern;
  3. a pair passes when settling keeps 85-115 % of the seed's mass, and over the verify run the
     number of creatures (found as Coolship does: basins of the potential K*A) never changes
     and the mass stays within 1.5 %. Among the passes, the one most like the screenshot wins:
     the correlation of per-channel radial profiles around the CoM (rotation-free), with the
     settled state shown as the screenshot would show it;
  4. its settled state, recentred and cropped, is the pattern, at R = --R.
Writes into --out (default: /tmp/from_screenshot_<name>): seed.png (the grid as read),
scales.png (each pair at the end of the verify run; rows = gammas, columns = scales),
verify.png (the chosen pattern at the start, middle and end of the verify run) and pattern.json.
"""
import argparse, itertools, json, os, re, sys
from fractions import Fraction
from multiprocessing import Pool

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
H, W, B = 256, 384, 4   # Coolship's board and block size


# ---- rule (a port of LeniaF.py for old-format 233.json entries) ------------------------------
def load_rule(src, code):
    lines = open(src).read().split('\n')
    for line in lines:
        line = line.strip().rstrip(',')
        if line.startswith('{') and json.loads(line)['code'].split('<')[0] == code:
            e = json.loads(line)
            break
    else:
        sys.exit(f'no entry with code {code} in {src}')
    ps = [dict(p) for p in e['params']]
    T, dn = ps[0]['T'], [0.0, 0.0, 0.0]
    for p in ps: dn[p['c'][1]] += p['h']
    for p in ps: p['h'] = round(p['h'] / dn[p['c'][1]], 3)   # free_h
    if max(p['h'] for p in ps) <= 0.5:
        for p in ps: p['h'] *= 2
        T *= 2
    return ps, T


def kernels(ps, R):
    y, x = np.mgrid[0:H, 0:W]
    D = np.hypot((x - W // 2) / R, (y - H // 2) / R)
    out = []
    for p in ps:
        b = np.array([float(Fraction(f)) for f in p['b'].split(',')])
        Br = len(b) * D / p['r']
        r = np.minimum(Br % 1, 1)
        k = (D < p['r']) * (r > 0) * (r < 1) * (4 * r * (1 - r)) ** 4 * b[np.minimum(Br.astype(int), len(b) - 1)]
        out.append(np.fft.fft2(np.fft.ifftshift(k / k.sum())))
    return out


def step(A, ps, Ks, dt):
    F = [np.fft.fft2(a) for a in A]
    D = np.zeros_like(A)
    for p, K in zip(ps, Ks):
        U = np.real(np.fft.ifft2(K * F[p['c'][0]]))
        D[p['c'][1]] += dt * p['h'] * (np.maximum(0, 1 - (U - p['m']) ** 2 / (9 * p['s'] ** 2)) ** 4 * 2 - 1)
    return np.clip(A + D, 0, 1)


# ---- helpers -----------------------------------------------------------------------------------
YY, XX = np.mgrid[0:H, 0:W]


def com(m):
    ay = np.angle((m * np.exp(2j * np.pi * YY / H)).sum())
    ax = np.angle((m * np.exp(2j * np.pi * XX / W)).sum())
    return np.array([ay * H, ax * W]) / (2 * np.pi)


def recentre(A):
    c = com(A.sum(0)).round().astype(int)
    return np.roll(A, (H // 2 - c[0], W // 2 - c[1]), (1, 2))


def parts(A, ps, Ks):
    """Creatures as Coolship's creatures.js finds them: per 4x4 block the potential K*A (mean
    over kernels) at the block centre; basins of its peaks above 0.05 (each block points to its
    highest 8-neighbour, on the torus), merged across saddles >= 0.8 x the lower peak, those
    smaller than a disk of radius R/4 dropped."""
    F = [np.fft.fft2(a) for a in A]
    P = np.mean([np.real(np.fft.ifft2(K * F[p['c'][0]])) for p, K in zip(ps, Ks)], 0)[B // 2::B, B // 2::B]
    bh, bw = P.shape
    idx = np.arange(bh * bw).reshape(bh, bw)
    bu, bi = P.copy(), idx.copy()
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            u, i = np.roll(P, (dy, dx), (0, 1)), np.roll(idx, (dy, dx), (0, 1))
            hi = u > bu
            bu, bi = np.where(hi, u, bu), np.where(hi, i, bi)
    ptr = bi.ravel()
    while not np.array_equal(nxt := ptr[ptr], ptr): ptr = nxt
    lab = np.where(P > 0.05, ptr.reshape(bh, bw), -1)
    pairs = {}
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        o, uo = np.roll(lab, (-dy, -dx), (0, 1)), np.roll(P, (-dy, -dx), (0, 1))
        m = (lab >= 0) & (o >= 0) & (lab != o)
        for x, y, h in zip(np.minimum(lab[m], o[m]), np.maximum(lab[m], o[m]), np.minimum(P[m], uo[m])):
            pairs[(x, y)] = max(pairs.get((x, y), -1), h)
    parent = {}
    def find(i):
        while parent.get(i, i) != i: i = parent[i]
        return i
    Pf = P.ravel()
    for (x, y), h in sorted(pairs.items(), key=lambda kv: -kv[1]):
        rx, ry = find(x), find(y)
        if rx != ry and h >= 0.8 * min(Pf[rx], Pf[ry]):
            if Pf[rx] >= Pf[ry]: parent[ry] = rx
            else: parent[rx] = ry
    roots = [find(i) for i in lab.ravel() if i >= 0]
    _, counts = np.unique(roots, return_counts=True)
    return int((counts * B * B >= np.pi * (CFG['R'] / 4) ** 2).sum())


def place(S):
    h, w = S.shape[1:]
    A = np.zeros((3, H, W))
    A[:, H // 2 - h // 2:H // 2 - h // 2 + h, W // 2 - w // 2:W // 2 - w // 2 + w] = S
    return A


def crop(A):
    A = np.where(A < 1 / 255, 0, A)
    ys, xs = np.nonzero(A.sum(0) > 0)
    return A[:, ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def val2ch(v):
    if v == 0: return '.'
    if v < 25: return chr(ord('A') + v - 1)
    return chr(ord('p') + (v - 25) // 24) + chr(ord('A') + (v - 25) % 24)


def cells2rle(a):
    rows = [''.join((str(n) if n > 1 else '') + val2ch(v)
                    for v, g in itertools.groupby(r) for n in [len(list(g))])
            for r in np.rint(a * 255).astype(int).tolist()]
    return '$'.join(rows) + '!'


def save_png(path, A, z=3):
    Image.fromarray((np.kron(np.clip(A, 0, 1).transpose(1, 2, 0), np.ones((z, z, 1))) * 255).astype(np.uint8)).save(path)


def window(A, r=64):
    return recentre(A)[:, H // 2 - r:H // 2 + r, W // 2 - r:W // 2 + r]


# ---- 1. the cell grid --------------------------------------------------------------------------
def read_grid(path, px=None, floor=0.03, keep_grey=False):
    im = np.asarray(Image.open(path).convert('RGB')).astype(float) / 255
    lum = im.sum(2)
    gx, gy = np.abs(np.diff(lum, axis=1)).sum(0), np.abs(np.diff(lum, axis=0)).sum(1)
    if px is None:   # strongest edge period between 2 and 16 px, over both axes
        cands = []
        for g in (gx, gy):
            f = np.abs(np.fft.rfft(g - g.mean()))
            k = np.arange(len(f))
            ok = (k > 0) & (len(g) / np.maximum(k, 1) >= 2) & (len(g) / np.maximum(k, 1) <= 16)
            cands.append(len(g) / k[ok][np.argmax(f[ok])])
        px = float(np.mean(cands))

    def score(g, ph, per):
        return g[np.round(np.arange(ph, len(g), per)).astype(int).clip(0, len(g) - 1)].sum()
    _, per, ph = max((score(gx, ph, per) + score(gy, ph, per), per, ph)
                     for per in np.arange(px * 0.97, px * 1.03, 0.01) for ph in np.arange(0, per, 0.1))
    ny, nx = int((im.shape[0] - ph) // per), int((im.shape[1] - ph) // per)
    A = np.zeros((3, ny, nx))
    inset = max(0.6, per * 0.15)   # skip the cell's blurred edges
    for j in range(ny):
        for i in range(nx):
            x0, y0 = ph + i * per + 1, ph + j * per + 1
            xs = slice(int(np.ceil(x0 + inset)), int(np.floor(x0 + per - 1 - inset)) + 1)
            ys = slice(int(np.ceil(y0 + inset)), int(np.floor(y0 + per - 1 - inset)) + 1)
            A[:, j, i] = np.median(im[ys, xs].reshape(-1, 3), 0)
    A[A < floor] = 0
    if not keep_grey:   # the zoom frame: dim grey, all channels equal
        A[:, (A.max(0) - A.min(0) < 0.06) & (A.max(0) < 0.35)] = 0
    return crop(A), per, ph


# ---- 2-3. settle each (gamma, scale) -----------------------------------------------------------
CFG = {}


def radial(A, z=1.0, rmax=None):
    """Per channel, the mean value at each distance (in cells / z) from the total CoM: a
    rotation-free fingerprint of the shape. A: (3, h, w), not on the torus."""
    m = A.sum(0)
    yy, xx = np.mgrid[0:A.shape[1], 0:A.shape[2]]
    cy, cx = (m * yy).sum() / m.sum(), (m * xx).sum() / m.sum()
    r = np.hypot(yy - cy, xx - cx) / z
    rmax = rmax or int(r[m > 0].max()) + 2
    k = np.minimum(r.astype(int), rmax)
    prof = np.stack([np.bincount(k.ravel(), c.ravel(), rmax + 1) / np.maximum(np.bincount(k.ravel(), None, rmax + 1), 1)
                     for c in A])
    return prof[:, :rmax]


def similarity(raw, final, g, z):
    """Correlation of the radial profiles of the screenshot (as displayed) and the settled state
    shown the way the screenshot shows it (values ** (1/g)), sizes matched by z."""
    a = radial(raw)
    b = radial(crop(window(final, 100)) ** (1 / g), z, a.shape[1])
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


def run_job(job):
    g, z = job
    ps, T, R, raw, settle, verify = (CFG[k] for k in ('ps', 'T', 'R', 'raw', 'settle', 'verify'))
    dt, Ks = 1 / T, kernels(ps, R)
    S = np.stack([ndimage.zoom(c, z, order=1) for c in raw ** g]).clip(0, 1)
    if max(S.shape[1:]) > min(H, W) // 2:
        return g, z, None, None, [], 0.0
    A = place(S)
    start = A.sum()
    for _ in range(settle):
        A = step(A, ps, Ks, dt)
    settled, hist = A.copy(), []
    for t in range(verify + 1):
        if t % 100 == 0:
            hist.append((t, float(A.sum()) / start, parts(A, ps, Ks)))
        if t < verify:
            A = step(A, ps, Ks, dt)
    sim = similarity(raw, A, g, z) if A.sum() > 1 else 0.0
    return g, z, settled, A, hist, sim


def verdict(hist):
    """(passed, reason, mass kept through settling)."""
    if not hist or hist[0][1] < 0.01:
        return False, 'died', 0
    n = {h[2] for h in hist}
    mass = np.array([h[1] for h in hist])
    kept, drift = float(mass[0]), float(np.abs(mass / mass[0] - 1).max())
    info = f'creatures {sorted(n)}, kept {kept:.0%} of the seed mass, then drift {drift:.2%}'
    ok = len(n) == 1 and 0 not in n and drift <= DRIFT and abs(kept - 1) <= KEEP
    return ok, info, kept


KEEP = 0.15    # a pass keeps 85-115 % of the seed's mass through settling,
DRIFT = 0.015  # and then stays within 1.5 %


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('png')
    ap.add_argument('--code', required=True, help='rule code in 233.json, e.g. Z18A9R')
    ap.add_argument('--name', required=True, help='form name, e.g. trio')
    ap.add_argument('--R', type=int, default=18, help='the R the screenshot was taken at (default 18)')
    ap.add_argument('--px', type=float, help='screen px per cell (default: detected)')
    ap.add_argument('--gammas', default='1,1.4,1.8,2.2',
                    help='screenshot value v -> cell value v ** gamma; screenshots look brighter than the cells')
    ap.add_argument('--scales', default='0.86,0.9,0.95,1.0,1.05', help='seed scales to try')
    ap.add_argument('--settle', type=int, default=300, help='steps to settle off the screenshot (default 300)')
    ap.add_argument('--verify', type=int, default=500,
                    help='steps that must hold after settling (default 500); short on purpose, since a form '
                         'may be unstable or turn into another later and still be worth a pattern')
    ap.add_argument('--floor', type=float, default=0.03, help='screenshot values below this are 0')
    ap.add_argument('--keep-grey', action='store_true', help="don't drop dim grey cells")
    ap.add_argument('--src', default=os.path.expanduser('~/github/Lenia/Python/found/233.json'))
    ap.add_argument('--out', help='folder for images and pattern.json')
    ap.add_argument('--add', action='store_true', help='append to tools/extra_patterns.json')
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    a = ap.parse_args()

    out = a.out or os.path.join('/tmp', 'from_screenshot_' + re.sub(r'\W+', '_', a.name))
    os.makedirs(out, exist_ok=True)
    raw, per, ph = read_grid(a.png, a.px, a.floor, a.keep_grey)
    save_png(os.path.join(out, 'seed.png'), raw, 5)
    print(f'grid: {per:.2f} px per cell, phase {ph:.1f}; seed {raw.shape[2]}x{raw.shape[1]} cells')

    ps, T = load_rule(a.src, a.code)
    CFG.update(ps=ps, T=T, R=a.R, raw=raw, settle=a.settle, verify=a.verify)
    gammas = [float(x) for x in a.gammas.split(',')]
    scales = [float(x) for x in a.scales.split(',')]
    jobs = [(g, z) for g in gammas for z in scales]
    print(f'settling {len(jobs)} (gamma, scale) pairs: {a.settle} + {a.verify} steps each, at R = {a.R}...', flush=True)
    with Pool(min(a.jobs, len(jobs))) as pool:
        res = pool.map(run_job, jobs)

    best, rows = None, []
    for g in gammas:
        row = []
        for rg, z, settled, final, hist, sim in res:
            if rg != g: continue
            ok, why, kept = verdict(hist)
            print(f'  gamma {g:.1f} scale {z:.2f}: {"PASS" if ok else "fail"}  shape {sim:+.2f}  {why}')
            row.append(window(final) if final is not None else np.zeros((3, 128, 128)))
            if ok and (best is None or sim > best[5]):
                best = (g, z, settled, final, hist, sim)
        rows.append(np.concatenate(row, 2))
    save_png(os.path.join(out, 'scales.png'), np.concatenate(rows, 1), 1)
    if best is None:
        sys.exit(f'nothing passed; see {out}/scales.png (rows: gammas, columns: scales; '
                 'try other --gammas, --scales, --px or --R)')

    g, z, settled, final, hist, sim = best
    if sim < 0.8:
        print(f'WARNING: the best pass looks unlike the screenshot (shape {sim:+.2f}); check verify.png')
    mass = np.array([h[1] for h in hist]); drift = float(np.abs(mass / mass[0] - 1).max())
    cells = crop(recentre(settled))
    mid = settled.copy()
    dt, Ks = 1 / T, kernels(ps, a.R)
    for _ in range(a.verify // 2): mid = step(mid, ps, Ks, dt)
    save_png(os.path.join(out, 'verify.png'), np.concatenate([window(settled), window(mid), window(final)], 2), 2)
    entry = {'code': a.code, 'name': a.name, 'R': a.R,
             'source': f'Reverse-engineered from {os.path.basename(a.png)} (a Coolship screenshot at R = {a.R}) '
                       f'by tools/from_screenshot.py: {per:.2f} px per cell, values ** {g}, scaled by {z}, settled '
                       f'for {a.settle} steps under the LeniaF rule; {hist[-1][2]} creature(s), mass within '
                       f'{drift:.2%} for {a.verify} more steps, shape similarity {sim:.2f}.',
             'cells': [cells2rle(c) for c in cells]}
    json.dump(entry, open(os.path.join(out, 'pattern.json'), 'w'), indent=1)
    print(f'chose gamma {g}, scale {z}: {cells.shape[2]}x{cells.shape[1]} cells, mass {cells.sum():.0f}, shape {sim:.2f}')
    print(f'images: {out}/seed.png, {out}/scales.png, {out}/verify.png')
    if a.add:
        path = os.path.join(HERE, 'extra_patterns.json')
        extra = json.load(open(path)) if os.path.exists(path) else []
        extra = [x for x in extra if (x['code'], x['name']) != (a.code, a.name)] + [entry]
        json.dump(extra, open(path, 'w'), indent=1)
        print(f'added to {path}; now run: python3 tools/gen_aquarium.py')


if __name__ == '__main__':
    main()
