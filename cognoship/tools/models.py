#!/usr/bin/env python3
"""Maneuver experiments: which action, in which condition, gives which outcome.

Runs Cognoship's Lenia rule headless on the GPU (torch), thousands of boards at once, and
writes a morph model: for each morph, how to steer it (turn clockwise / anticlockwise) and how to turn
it into another morph, and for pairs of morphs, how a clone placed beside a ship makes a twin.
For game AI (bots steering towards the player, autopilot, transforming at will) and for reading.

Needs torch (CUDA) with numpy, scipy, matplotlib, Pillow; e.g. the evocosmos-drama venv:
  PY=~/github/evocosmos-drama/.venv/bin/python
  $PY tools/models.py steer --code 5N7KKM --morph limus            # G mask, named spots, all phases
  $PY tools/models.py clone --code GDNQYX --pair gyrans,gyrans     # clone sweep of a pair
  $PY tools/models.py book  --code 5N7KKM                          # models/<code>.md + .json
  $PY tools/models.py show  --code 5N7KKM --morph limus --spot 'right edge' --hold 30  # GIF

What is the same as the game (so an action found here works there):
  - the rule: aquarium233.js (tweaks applied), LeniaF's free_h, polynomial kernels and growth,
    R = 18 (patterns resampled from their native R as rule.js formCells), dt = 1 / T;
  - the mask: sim.glsl's cognitive mask, U = K*(A(1-B)) / K*(1-B), a disc of radius 4 cells
    (core.js: round(2.5 R / 12)) over one channel, following the creature's centre of mass;
    held in the heading-up frame (Upward on): its offset turns with core.js's slowHeading;
  - the kinematics (speed, turn rate, heading): core.js updateKinematics, ported;
  - the morph names: detect.js, ported, against assets/stats/<code>.json (parts and wholes).
What differs: a 128 x 128 torus (192 x 192 for clones) with only the experiment on it, FFT in
float32 (the game's default is FFT too); creatures are found by thresholded connected regions
(not core.js's potential basins), tracked by nearest centre.

Frames and signs (board: x = column right, y = row down, as seen on screen):
  - a morph's frame is its heading (the velocity's direction, core.js slowHeading) if it moves,
    else its body axis (the B channel's CoM to the G channel's, core.js 'axis');
  - an offset (fwd, side) in R units: fwd along the frame, side to its right as seen on screen;
  - a turn is in degrees, + = clockwise as seen (to the right), - = to the left.

Outputs: raw results in output/models/<code>/ (gitignored); the book in models/.
"""
import argparse, itertools, json, math, os, re, subprocess, sys, time
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'output', 'models')
BOOK = os.path.join(ROOT, 'models')

R_GAME = 18
MASK_R = 4                       # core.js: round(2.5 * R / 12) at R = 18

# core.js kinematics constants (in units of R and unit time)
SPEED_MIN_R = 0.005
VEL_TAU, ARROW_TAU, OMEGA_TAU = 15.0, 2.5, 5.0
JUMP_MAX, LAG_MAX = math.pi / 4, math.pi / 3

CH = 'RGB'

# ---- the rule and the patterns (aquarium233.js) -------------------------------------------------
_AQ = None


def aquarium():
    global _AQ
    if _AQ is None:
        js = "const A=require(process.argv[1]);process.stdout.write(JSON.stringify(A))"
        cache = os.path.join(OUT, 'aquarium233.json')   # (a machine without node, e.g. a GPU server: written by the sync)
        try:
            _AQ = json.loads(subprocess.run(['node', '-e', js, os.path.join(ROOT, 'aquarium233.js')],
                                            capture_output=True, text=True, check=True).stdout)
        except FileNotFoundError:
            _AQ = json.load(open(cache))
    return _AQ


def frac(s):
    a, *b = s.split('/')
    return float(a) / float(b[0]) if b else float(a)


def rule_of(code):
    r = aquarium()['rules'][code]
    ks = [dict(k) for k in r['kernels']]
    T = r['T']
    dn = [0.0, 0.0, 0.0]
    for k in ks: dn[k['c'][1]] += k['h']
    for k in ks: k['h'] = math.floor(k['h'] / dn[k['c'][1]] * 1000 + 0.5) / 1000   # rule.js leniaF
    if max(k['h'] for k in ks) <= 0.5:
        for k in ks: k['h'] *= 2
        T *= 2
    return {'code': code, 'kernels': ks, 'T': T, 'dt': 1 / T, 'wholes': r.get('wholes', [])}


def ch2val(c):
    if c in '.b': return 0
    if c == 'o': return 255
    if len(c) == 1: return ord(c) - 65 + 1
    return (ord(c[0]) - 112) * 24 + (ord(c[1]) - 65 + 25)


def rle2cells(st):
    rows, row, last, count = [], [], '', ''
    st = st.rstrip('!') + '$'
    for ch in st:
        if ch.isdigit(): count += ch
        elif ch in 'pqrstuvwxy@': last = ch
        else:
            n = int(count) if count else 1
            if ch == '$':
                rows.append(row); rows.extend([] for _ in range(n - 1)); row = []
            else:
                row.extend([ch2val(last + ch) / 255] * n)
            last, count = '', ''
    w = max(len(r) for r in rows)
    return np.array([r + [0] * (w - len(r)) for r in rows], float)


def resample(a, z):
    h, w = a.shape
    H2, W2 = max(1, math.floor(h * z + 0.5)), max(1, math.floor(w * z + 0.5))
    ys = np.minimum((np.arange(H2) / z).astype(int), h - 1)
    xs = np.minimum((np.arange(W2) / z).astype(int), w - 1)
    return a[ys][:, xs]


def pattern(code, morph, R=R_GAME):
    """The morph's cells (3, h, w) at radius R, as rule.js formCells."""
    pats = [p for p in aquarium()['patterns'] if p['rule'] == code]
    # by id (CODE-morph), else by its name (VT049W-gyrans-replicating is 'VT049W gyrans 💥')
    named = lambda p: p['name'].replace(code + ' ', '').replace(' \U0001F4A5', '').strip() == morph
    for p in [p for p in pats if p['id'].split('-', 1)[1] == morph] or [p for p in pats if named(p)]:
            ch = [resample(rle2cells(s), R / p['R']) for s in p['cells']]
            h, w = max(c.shape[0] for c in ch), max(c.shape[1] for c in ch)
            out = np.zeros((3, h, w))
            for i, c in enumerate(ch): out[i, :c.shape[0], :c.shape[1]] = c
            return out
    sys.exit(f'no pattern {code}-{morph}; have: ' + ', '.join(p['id'] for p in aquarium()['patterns'] if p['rule'] == code))


def morphs_of(code):
    return [p['id'].split('-', 1)[1] for p in aquarium()['patterns'] if p['rule'] == code]


def transform(cells, angle, flip=False):
    """cells turned by angle (rad, clockwise as seen) after an optional left-right flip, about the
    CoM, bilinear (as clone.js build); returns (3, S, S) with the CoM at the centre."""
    m = cells.sum(0)
    h, w = m.shape
    yy, xx = np.mgrid[0:h, 0:w]
    cy, cx = (m * yy).sum() / m.sum(), (m * xx).sum() / m.sum()
    S = int(math.ceil(math.hypot(h, w))) + 4
    c, s = math.cos(angle), math.sin(angle)
    M = np.array([[c, -s], [s, c]]) @ (np.diag([-1, 1]) if flip else np.eye(2))
    inv = np.linalg.inv(M)
    oy, ox = np.mgrid[0:S, 0:S] - (S - 1) / 2
    sx = inv[0, 0] * ox + inv[0, 1] * oy + cx
    sy = inv[1, 0] * ox + inv[1, 1] * oy + cy
    x0, y0 = np.floor(sx).astype(int), np.floor(sy).astype(int)
    fx, fy = sx - x0, sy - y0

    def at(k, y, x):
        ok = (y >= 0) & (y < h) & (x >= 0) & (x < w)
        return np.where(ok, cells[k][np.clip(y, 0, h - 1), np.clip(x, 0, w - 1)], 0)
    return np.array([(at(k, y0, x0) * (1 - fx) + at(k, y0, x0 + 1) * fx) * (1 - fy)
                     + (at(k, y0 + 1, x0) * (1 - fx) + at(k, y0 + 1, x0 + 1) * fx) * fy for k in range(3)])


def axis_angle(cells):
    """Board angle of the B -> G body axis of a small form (no wrap)."""
    h, w = cells.shape[1:]
    yy, xx = np.mgrid[0:h, 0:w]
    cg = [(cells[k] * yy).sum() / cells[k].sum() for k in (1, 2)], [(cells[k] * xx).sum() / cells[k].sum() for k in (1, 2)]
    return math.atan2(cg[0][0] - cg[0][1], cg[1][0] - cg[1][1])


# ---- the simulation (batched, torch) ----------------------------------------------------------
class Sim:
    def __init__(self, rule, H, W, R=R_GAME, device='cuda'):
        import torch
        self.t = torch
        self.H, self.W, self.R, self.dev = H, W, R, device
        ks = rule['kernels']
        y, x = np.mgrid[0:H, 0:W]
        D = np.hypot((x - W // 2) / R, (y - H // 2) / R)
        Kf = []
        for k in ks:   # rule.js shell, kn = 1 (polynomial quad4)
            b = np.array([frac(f) for f in k['b'].split(',')])
            Br = len(b) * D / k['r']
            f = np.minimum(Br % 1, 1)
            kk = (D < k['r']) * (f > 0) * (f < 1) * (4 * f * (1 - f)) ** 4 * b[np.minimum(Br.astype(int), len(b) - 1)]
            Kf.append(np.fft.rfft2(np.fft.ifftshift(kk / kk.sum())))
        self.K = torch.tensor(np.array(Kf), dtype=torch.complex64, device=device)
        self.src = torch.tensor([k['c'][0] for k in ks], device=device)
        self.dst = torch.tensor([k['c'][1] for k in ks], device=device)
        f32 = dict(dtype=torch.float32, device=device)
        self.m = torch.tensor([k['m'] for k in ks], **f32)[None, :, None, None]
        self.s = torch.tensor([k['s'] for k in ks], **f32)[None, :, None, None]
        self.h = torch.tensor([k['h'] for k in ks], **f32)[None, :, None, None]
        self.dt = rule['dt']
        self.Y = torch.arange(H, **f32)[:, None]
        self.X = torch.arange(W, **f32)[None, :]
        self.cy = torch.cos(2 * math.pi * torch.arange(H, **f32) / H)
        self.sy = torch.sin(2 * math.pi * torch.arange(H, **f32) / H)
        self.cx = torch.cos(2 * math.pi * torch.arange(W, **f32) / W)
        self.sx = torch.sin(2 * math.pi * torch.arange(W, **f32) / W)

    def step(self, A, vis=None):
        t = self.t
        F = t.fft.rfft2(A if vis is None else A * vis)
        U = t.fft.irfft2(F[:, self.src] * self.K, s=(self.H, self.W))
        if vis is not None:
            N = t.fft.irfft2(t.fft.rfft2(vis)[:, self.src] * self.K, s=(self.H, self.W))
            U = t.where(N > 1e-6, U / N.clamp_min(1e-6), t.zeros_like(U))
        q = (1 - (U - self.m) ** 2 / (9 * self.s ** 2)).clamp_min(0)
        G = self.h * ((q * q) ** 2 * 2 - 1)
        Dl = t.zeros_like(A).index_add_(1, self.dst, G)
        return (A + self.dt * Dl).clamp(0, 1)

    def disc(self, cx, cy, r):
        """(B, H, W) discs on the torus, centres cx, cy (B,) tensors, radius r (B,) or float"""
        t = self.t
        dx = (self.X[None] - cx[:, None, None] + self.W / 2) % self.W - self.W / 2
        dy = (self.Y[None] - cy[:, None, None] + self.H / 2) % self.H - self.H / 2
        r = r if isinstance(r, float) else r[:, None, None]
        return (dx * dx + dy * dy <= r * r).float()

    def measure(self, A):
        """per board: [com_x, com_y, mR, mG, mB, gx, gy, bx, by] (CoMs on the torus, cells)"""
        t = self.t
        ms = A.sum((2, 3))
        ry, rx = A.sum(3), A.sum(2)     # (B, 3, H), (B, 3, W)
        tot_y, tot_x = ry.sum(1), rx.sum(1)

        def ang(v, c, s, n):
            return (t.atan2(v @ s, v @ c) % (2 * math.pi)) * n / (2 * math.pi)
        out = t.stack([ang(tot_x, self.cx, self.sx, self.W), ang(tot_y, self.cy, self.sy, self.H),
                       ms[:, 0], ms[:, 1], ms[:, 2],
                       ang(rx[:, 1], self.cx, self.sx, self.W), ang(ry[:, 1], self.cy, self.sy, self.H),
                       ang(rx[:, 2], self.cx, self.sx, self.W), ang(ry[:, 2], self.cy, self.sy, self.H)], 1)
        return out.cpu().numpy().astype(float)


def jdefault(o):
    if isinstance(o, np.generic): return o.item()
    if isinstance(o, np.ndarray): return o.tolist()
    raise TypeError(type(o))



def mode(xs):
    """the most common of xs; a tie goes to the one seen first (not set order, which changes from run
    to run with Python's string hashing: models rebuilt the same way each time)"""
    return Counter(xs).most_common(1)[0][0]

def wrapd(d, n): return (d + n / 2) % n - n / 2
def wrapa(a): return (a + math.pi) % (2 * math.pi) - math.pi


class Kin:
    """core.js updateKinematics, vectorised over n creatures. v in cells per unit time."""
    def __init__(self, n, R, T):
        self.R, self.T = R, T
        z = lambda: np.zeros(n)
        self.age, self.omega, self.omegaAge, self.speed, self.velAng = z(), z(), z(), z(), z()
        self.fast, self.vel = np.zeros((n, 2)), np.zeros((n, 2))
        self.fastAng, self.slow = np.full(n, np.nan), np.full(n, np.nan)

    def copy_rows(self, idx):
        k = Kin(len(idx), self.R, self.T)
        for a in ('age', 'omega', 'omegaAge', 'speed', 'velAng', 'fast', 'vel', 'fastAng', 'slow'):
            setattr(k, a, getattr(self, a)[idx].copy())
        return k

    def update(self, v):
        T, smin = self.T, SPEED_MIN_R * self.R
        self.age += T
        ease = lambda tau: np.maximum(1 - math.exp(-T / tau), T / self.age)
        a = ease(ARROW_TAU)
        self.fast += a[:, None] * (v - self.fast)
        fs = np.hypot(self.fast[:, 0], self.fast[:, 1])
        fAng = np.where(fs > smin / 2, np.arctan2(self.fast[:, 1], self.fast[:, 0]), np.nan)
        both = ~np.isnan(fAng) & ~np.isnan(self.fastAng)
        self.omegaAge = np.where(both, self.omegaAge + T, self.omegaAge)
        b = np.maximum(1 - math.exp(-T / OMEGA_TAU), T / np.maximum(self.omegaAge, 1e-9))
        d = np.where(both, (fAng - self.fastAng + math.pi) % (2 * math.pi) - math.pi, 0)
        self.omega = np.where(both, self.omega + b * (d / T - self.omega), self.omega * math.exp(-T / OMEGA_TAU))
        self.fastAng = fAng
        wT = self.omega * T
        mx, my = (1 - (1 - a) * np.cos(wT)) / a, (1 - a) * np.sin(wT) / a
        tt = np.clip(np.arctan2(my, mx), -LAG_MAX, LAG_MAX)
        g = np.minimum(2, np.hypot(mx, my))
        c, s = np.cos(tt), np.sin(tt)
        self.vel = np.stack([(self.fast[:, 0] * c - self.fast[:, 1] * s) * g, (self.fast[:, 0] * s + self.fast[:, 1] * c) * g], 1)
        self.speed = np.hypot(self.vel[:, 0], self.vel[:, 1])
        self.velAng = np.arctan2(self.vel[:, 1], self.vel[:, 0])
        vs = np.hypot(v[:, 0], v[:, 1])
        ang = np.where(vs > smin / 2, np.arctan2(v[:, 1], v[:, 0]), np.nan)
        new = np.isnan(self.slow)
        pred = (self.slow + self.omega * T + math.pi) % (2 * math.pi) - math.pi
        r = np.clip((ang - pred + math.pi) % (2 * math.pi) - math.pi, -JUMP_MAX, JUMP_MAX)
        pulled = (pred + ease(VEL_TAU) * r + math.pi) % (2 * math.pi) - math.pi
        self.slow = np.where(new, ang, np.where(np.isnan(ang), pred, pulled))


# ---- morph detection (detect.js, ported) ------------------------------------------------------
MASS0, SPEED0, OMEGA0 = 0.01, 0.002, 0.01
WIDTHS = np.array([0.03, 0.03, 0.03, 0.3, 0.3])
FUZZ, NEFF, P_MIN, FIT_MIN, CLAMP, REF_MAX, WARM, MIN_N, WIN = 0.75, 2, 0.9, -6, -20, 250, 10, 25, 50


def feat(x):
    x = np.asarray(x, float)
    return np.log(x + np.array([MASS0, MASS0, MASS0, SPEED0, OMEGA0]))


def detect_models(code):
    """{'part': model, 'whole': model} from assets/stats/<code>.json"""
    path = os.path.join(ROOT, 'assets', 'stats', f'{code}.json')
    if not os.path.exists(path): return {}
    obj = json.load(open(path))
    idx = [obj['keys'].index(k) for k in ('mR', 'mG', 'mB', 'speed', 'omega')]
    rows, segs = [], []   # flatten, as StatsChart.log.load + push
    for sg in obj['segs']:
        lvl = 'whole' if sg.get('level') == 'whole' else 'part'
        for n, row in enumerate(sg['data']):
            gap = (n == 0 and sg.get('gap', False)) or not rows
            fz = bool(sg.get('fuzzy'))
            last = segs[-1] if segs else None
            if gap or last['morph'] != sg['morph'] or last['fuzzy'] != fz or last['level'] != lvl:
                segs.append({'start': len(rows), 'morph': sg['morph'], 'gap': gap, 'fuzzy': fz, 'level': lvl})
            rows.append([np.nan if row[i] is None else row[i] for i in idx])
    rows = np.array(rows, float)
    models = {}
    for level in ('part', 'whole'):
        by, fuzzy = {}, {}
        for i, sg in enumerate(segs):
            if sg['level'] != level: continue
            end = segs[i + 1]['start'] if i + 1 < len(segs) else len(rows)
            g = i
            while g > 0 and not segs[g]['gap']: g -= 1
            p0 = segs[g]['start'] + WARM
            key = ('~' if sg['fuzzy'] else '') + sg['morph']
            for j in range(max(sg['start'], p0), end):
                x = rows[j]
                if np.all(x >= 0): by.setdefault(key, []).append(x)
            if sg['fuzzy']: fuzzy[key] = True
        ms = []
        for key, pts in by.items():
            if len(pts) < MIN_N: continue
            pts = np.array(pts)
            stride = max(1, len(pts) / REF_MAX)
            ii, q = [], 0.0
            while q < len(pts): ii.append(int(q)); q += stride
            ref = feat(pts[ii])
            widths = WIDTHS
            if key in fuzzy:
                fs = feat(pts)
                widths = np.maximum(WIDTHS, FUZZ * fs.std(0))
            ms.append({'name': key.lstrip('~'), 'fuzzy': key in fuzzy, 'ref': ref, 'inv': 1 / widths,
                       'logNorm': np.log(widths).sum()})
        if ms: models[level] = ms
    return models


def classify(model, xs):
    """detect.js classify on samples xs (n, 5) (oldest first, already past the warm-up):
    (name or None, best, p, fit)"""
    zs = feat(xs[-WIN:])
    if len(zs) < MIN_N: return None, None, 0.0, -99
    L = []
    for mp in model:
        d = (((zs[:, None, :] - mp['ref'][None]) * mp['inv']) ** 2).sum(2)
        sc = np.maximum(CLAMP, np.log(np.exp(-d / 2).mean(1) + 1e-300))
        L.append(sc.mean())
    L = np.array(L)
    N = L - np.array([mp['logNorm'] for mp in model])
    b = int(N.argmax())
    e = np.exp(NEFF * (N - N[b]))
    p = e[b] / e.sum()
    known = p >= P_MIN and L[b] >= FIT_MIN
    return (model[b]['name'] if known else None), model[b]['name'], float(p), float(L[b])


# ---- bodies: regions on the torus -------------------------------------------------------------
def label_torus(m, thr, grow):
    """connected regions of m > thr (dilated by grow cells) on the torus: labels (H, W), n"""
    from scipy import ndimage
    mask = m > thr
    if grow > 0:
        mask = ndimage.binary_dilation(mask, iterations=grow) if not mask.all() else mask
        # wrap the dilation: also dilate a rolled copy (cheap way to cross edges)
        sh = (m.shape[0] // 2, m.shape[1] // 2)
        mask |= np.roll(ndimage.binary_dilation(np.roll(m > thr, sh, (0, 1)), iterations=grow), (-sh[0], -sh[1]), (0, 1))
    lab, n = ndimage.label(mask)
    parent = list(range(n + 1))

    def find(i):
        while parent[i] != i: parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for a, b in itertools.chain(zip(lab[0], lab[-1]), zip(lab[:, 0], lab[:, -1])):
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb: parent[ra] = rb
    # (a vector remap: a board shattered into thousands of bits took hours with list.index)
    rt = np.array([find(i) for i in range(1, n + 1)], int)
    roots, inv = np.unique(rt, return_inverse=True)
    remap = np.zeros(n + 1, int)
    remap[1:] = inv + 1
    return remap[lab], len(roots)


def bodies_of(A, R, debris):
    """bodies of a board (3, H, W) numpy: list of (cx, cy, mass/R^2), heaviest first, debris dropped"""
    m = A.sum(0)
    lab, n = label_torus(m, 0.02, 2)
    H, W = m.shape
    yy, xx = np.mgrid[0:H, 0:W]
    # (each body's mass and torus centre in one pass over the board: bincount by label; a board
    # shattered into thousands of bits was a pass per bit, hours)
    L, w = lab.ravel(), m.ravel()
    mass = np.bincount(L, weights=w, minlength=n + 1) / R ** 2
    ey = np.bincount(L, weights=w * np.cos(2 * np.pi * yy.ravel() / H), minlength=n + 1) + 1j * np.bincount(L, weights=w * np.sin(2 * np.pi * yy.ravel() / H), minlength=n + 1)
    ex = np.bincount(L, weights=w * np.cos(2 * np.pi * xx.ravel() / W), minlength=n + 1) + 1j * np.bincount(L, weights=w * np.sin(2 * np.pi * xx.ravel() / W), minlength=n + 1)
    out = [(np.angle(ex[i]) * W / (2 * np.pi) % W, np.angle(ey[i]) * H / (2 * np.pi) % H, mass[i])
           for i in range(1, n + 1) if mass[i] >= debris]
    return sorted(out, key=lambda b: -b[2])


# ---- one batch of runs ------------------------------------------------------------------------
SERIES = ('x', 'y', 'mass', 'mR', 'mG', 'mB', 'head', 'speed', 'omega', 'velAng')
KMAX = 4          # bodies tracked per board in the outcome window
OUT_WIN = WARM + WIN + 10


def run_batch(sim, A0, kin0, com0, plan, steps, R, models, debris, log_every=0):
    """Run boards A0 (B, 3, H, W) torch, with kinematics kin0 (Kin) and unwrapped CoMs com0 (B, 2).
    plan: dict of per-board arrays: action ('mask' / 'erase' / 'draw' / 'none'), ch, r, dur,
    off (B, 2) in cells (fwd, side), frame ('vel' / 'axis'), axis0 (B,) the start axis angle.
    Returns per-board time series and the outcome classification."""
    t = sim.t
    B = A0.shape[0]
    A = A0.clone()
    kin = kin0
    com = com0.copy()
    meas = sim.measure(A)
    T_unit = sim.dt
    ser = {k: np.zeros((B, steps), np.float32) for k in SERIES}
    ch = t.tensor(plan['ch'], device=sim.dev)
    onehot1 = t.nn.functional.one_hot(ch, 3).float()[:, :, None, None]
    onehot2 = t.nn.functional.one_hot(t.tensor(plan.get('ch2', plan['ch']), device=sim.dev), 3).float()[:, :, None, None]
    rr = t.tensor(plan['r'], dtype=t.float32, device=sim.dev)
    act = np.array(plan['action'])
    dur = np.array(plan['dur'])
    off = np.array(plan['off'], float)
    mass0 = meas[:, 2:5].sum(1)
    # an optional second hold (channel ch2, else the same), from step at2 for dur2 steps at off2: a combo
    at2 = np.array(plan.get('at2', [10 ** 9] * B)); dur2 = np.array(plan.get('dur2', [0] * B))
    off2 = np.array(plan.get('off2', [(0.0, 0.0)] * B), float)
    off1 = off
    frame2_axis = np.array(plan.get('frame2', ['vel'] * B)) == 'axis'
    for s in range(steps):
        second = (s >= at2) & (s < at2 + dur2)
        holding = ((s < dur) | second) & (act != 'none')
        off = np.where(second[:, None], off2, off1)
        onehot = t.where(t.tensor(second, device=sim.dev)[:, None, None, None], onehot2, onehot1)
        # (the second hold in a still morph's frame: its body axis now, B's CoM to G's)
        use_axis = np.where(second, frame2_axis, np.array(plan['frame']) == 'axis')
        ax_now = np.arctan2(wrapd(meas[:, 6] - meas[:, 8], sim.H), wrapd(meas[:, 5] - meas[:, 7], sim.W))
        # (held fixed from the second hold's start, as the sweeps hold a still morph's spot: its axis then)
        ax2 = np.where(s == at2, ax_now, ax2) if s else np.where(s == at2, ax_now, np.zeros(B))
        plan_axis0 = np.where(second, ax2, plan['axis0'])
        frame_ang = np.where(use_axis, plan_axis0, np.where(np.isnan(kin.slow), plan_axis0, kin.slow))
        if holding.any():
            c, sn = np.cos(frame_ang), np.sin(frame_ang)
            # (fwd, side): fwd along the frame, side to its right as seen (y down): right = (-fy, fx)
            dx = off[:, 0] * c - off[:, 1] * sn
            dy = off[:, 0] * sn + off[:, 1] * c
            cx = t.tensor(meas[:, 0] + dx, dtype=t.float32, device=sim.dev)
            cy = t.tensor(meas[:, 1] + dy, dtype=t.float32, device=sim.dev)
            hold_t = t.tensor(holding, device=sim.dev)[:, None, None, None]
            disc = sim.disc(cx, cy, rr)[:, None] * onehot * hold_t   # (B, 3, H, W)
            is_mask = t.tensor(act == 'mask', device=sim.dev)[:, None, None, None]
            vis = 1 - disc * is_mask
            A = sim.step(A, vis if (holding & (act == 'mask')).any() else None)
            er = disc * t.tensor(act == 'erase', device=sim.dev)[:, None, None, None]
            dr = disc * t.tensor(act == 'draw', device=sim.dev)[:, None, None, None]
            A = t.maximum(A * (1 - er), dr)
        else:
            A = sim.step(A)
        m2 = sim.measure(A)
        d = np.stack([wrapd(m2[:, 0] - meas[:, 0], sim.W), wrapd(m2[:, 1] - meas[:, 1], sim.H)], 1)
        com += d
        meas = m2
        kin.update(d / T_unit)
        ser['x'][:, s], ser['y'][:, s] = com[:, 0], com[:, 1]
        ser['mass'][:, s] = meas[:, 2:5].sum(1) / R ** 2
        ser['mR'][:, s], ser['mG'][:, s], ser['mB'][:, s] = (meas[:, 2:5] / R ** 2).T
        ser['head'][:, s] = np.where(np.isnan(kin.slow), 0, kin.slow)
        ser['speed'][:, s] = kin.speed / R
        ser['omega'][:, s] = kin.omega
        ser['velAng'][:, s] = kin.velAng
        if log_every and s % log_every == 0:
            print(f'    step {s}/{steps}', flush=True)
    # outcome window: rerun the last OUT_WIN steps? No: track bodies from a snapshot OUT_WIN
    # steps before the end. Simpler and exact: keep stepping OUT_WIN more steps from the end.
    out = outcome(sim, A, R, models, debris, mass0 / R ** 2)
    return ser, out, A


def outcome(sim, A, R, models, debris, mass0):
    """Track each board's bodies for OUT_WIN steps from A, and name them (detect.js)."""
    t = sim.t
    B = A.shape[0]
    An = A.cpu().numpy()
    bods = [bodies_of(An[i], R, debris) for i in range(B)]
    nb0 = np.array([len(b) for b in bods])
    K = KMAX
    cen = np.zeros((B, K, 2)); valid = np.zeros((B, K), bool)
    for i, bl in enumerate(bods):
        for k, (x, y, _) in enumerate(bl[:K]): cen[i, k] = (x, y); valid[i, k] = True
    unw = cen.copy()                                         # unwrapped body centres
    kins = Kin(B * K, R, sim.dt)
    wkin = Kin(B, R, sim.dt)
    wcom = None
    feats = np.full((B, K, OUT_WIN, 5), np.nan); wfeats = np.full((B, OUT_WIN, 5), np.nan)
    womega = np.zeros((B, OUT_WIN)); bomega = np.zeros((B, K, OUT_WIN))
    lim2 = (2 * R) ** 2
    for s in range(OUT_WIN):
        A = sim.step(A)
        cx = t.tensor(cen[..., 0], dtype=t.float32, device=sim.dev)[:, :, None, None]
        cy = t.tensor(cen[..., 1], dtype=t.float32, device=sim.dev)[:, :, None, None]
        dx = (sim.X[None, None] - cx + sim.W / 2) % sim.W - sim.W / 2
        dy = (sim.Y[None, None] - cy + sim.H / 2) % sim.H - sim.H / 2
        d2 = dx * dx + dy * dy
        d2 = t.where(t.tensor(valid, device=sim.dev)[:, :, None, None], d2, t.full_like(d2, 1e12))
        near = d2.argmin(1, keepdim=True)
        own = (t.zeros_like(d2).scatter_(1, near, 1.0)) * (d2 <= lim2).float()    # (B, K, H, W)
        m = A.sum(1, keepdim=True)                                                 # (B, 1, H, W)
        w = own * m
        mass = w.sum((2, 3)).clamp_min(1e-9)
        mx = (w * dx).sum((2, 3)) / mass
        my = (w * dy).sum((2, 3)) / mass
        chm = t.einsum('bkhw,bchw->bkc', own, A)
        mx, my, chm = mx.cpu().numpy(), my.cpu().numpy(), chm.cpu().numpy()
        cen[..., 0] = (cen[..., 0] + mx) % sim.W
        cen[..., 1] = (cen[..., 1] + my) % sim.H
        unw[..., 0] += mx; unw[..., 1] += my
        v = np.stack([mx, my], -1).reshape(B * K, 2) / sim.dt
        kins.update(v)
        sp = kins.speed.reshape(B, K) / R; om = kins.omega.reshape(B, K)
        feats[:, :, s] = np.concatenate([chm / R ** 2, sp[..., None], np.abs(om)[..., None]], -1)
        bomega[:, :, s] = om
        # the whole: every valid body together
        wm = (chm.sum(2) * valid)
        tot = wm.sum(1).clip(1e-9)
        wc = (unw * wm[..., None]).sum(1) / tot[:, None]
        if wcom is not None:
            wkin.update((wc - wcom) / sim.dt)
        wcom = wc
        wfeats[:, s] = np.concatenate([(chm * valid[..., None]).sum(1) / R ** 2, (wkin.speed / R)[:, None], np.abs(wkin.omega)[:, None]], -1)
        womega[:, s] = wkin.omega
    An = A.cpu().numpy()
    results = []
    for i in range(B):
        nb1 = len(bodies_of(An[i], R, debris))
        tot = float(np.nansum(feats[i, :, -1, :3][valid[i]]))
        r = {'bodies': int(nb0[i]), 'bodies_end': nb1, 'mass': round(tot, 3),
             'mass_rel': round(tot / max(mass0[i], 1e-9), 3)}
        if nb0[i] == 0 or tot < 0.25 * mass0[i]:
            r['label'] = 'dissolved'
        elif nb0[i] > K or tot > 4 * mass0[i]:
            r['label'] = 'multiplied' if nb0[i] > K else 'exploded'
        else:
            names, parts = [], []
            for k in range(min(K, nb0[i])):
                xs = feats[i, k, WARM:]
                pn = classify(models['part'], xs) if 'part' in models else (None, None, 0, -99)
                om = float(bomega[i, k, WARM:].mean())
                parts.append({'name': pn[0], 'best': pn[1], 'p': round(pn[2], 3), 'mass': round(float(xs[:, :3].sum(1).mean()), 3),
                              'speed': round(float(xs[:, 3].mean()), 4), 'omega': round(om, 4)})
            wx = wfeats[i, WARM + 1:]
            wn = classify(models['whole'], wx) if 'whole' in models else (None, None, 0, -99)
            whole = {'name': wn[0], 'best': wn[1], 'p': round(wn[2], 3), 'mass': round(float(wx[:, :3].sum(1).mean()), 3),
                     'speed': round(float(wx[:, 3].mean()), 4), 'omega': round(float(womega[i, WARM + 1:].mean()), 4)}
            r['parts'], r['whole'] = parts, whole
            if nb0[i] == 1:
                p = parts[0]
                # one region: a single morph, or a twin fused into one region
                if p['name'] and not (whole['name'] and whole['p'] > p['p'] and whole['name'] != p['name']):
                    r['label'] = p['name']
                elif whole['name']:
                    r['label'] = whole['name']
                else:
                    r['label'] = '???'
            else:
                spread = max(math.hypot(wrapd(a[0] - b[0], sim.W), wrapd(a[1] - b[1], sim.H))
                             for a, b in itertools.combinations([cen[i, k] for k in range(min(K, nb0[i]))], 2))
                r['spread'] = round(spread / R, 2)
                if whole['name'] and spread < 3.5 * R:
                    r['label'] = whole['name']
                else:
                    r['label'] = ' + '.join(sorted((p['name'] or '???') for p in parts))
            if nb1 != nb0[i] and r['label'] not in ('dissolved',):
                r['unsettled'] = True
        results.append(r)
    return results


# ---- start states -----------------------------------------------------------------------------
def place(sim, cells, quarter):
    f = np.rot90(cells, quarter % 4, (1, 2)).copy()
    b = np.zeros((3, sim.H, sim.W))
    h, w = f.shape[1:]
    y0, x0 = sim.H // 2 - h // 2, sim.W // 2 - w // 2
    b[:, y0:y0 + h, x0:x0 + w] = f
    return b


def cycle_period(mass, vx, vy):
    """The morph's own cycle in steps: the first strong autocorrelation peak of its mass pulse,
    else of its velocity (a circler's loop); 48 if neither shows one."""
    for sig in (mass, None):
        if sig is None:
            z = (vx - vx.mean()) + 1j * (vy - vy.mean())
            ac = np.real(np.array([np.vdot(z[:-l], z[l:]) for l in range(1, len(z) // 2)]))
            ac /= max(np.vdot(z, z).real, 1e-12)
        else:
            x = sig - sig.mean()
            if x.std() < 1e-6 * max(abs(sig.mean()), 1e-9): continue
            ac = np.array([np.dot(x[:-l], x[l:]) for l in range(1, len(x) // 2)]) / np.dot(x, x)
        for l in range(5, len(ac) - 1):
            if ac[l] > ac[l - 1] and ac[l] >= ac[l + 1] and ac[l] > 0.3:
                return l + 1
    return 48


PRE = WIN + 10     # stats kept from before the action, so the name box has a history


_MIRROR = {}


def start_set(sim, code, morph, starts, settle, R, mirror=None):
    """Settled copies of the morph. starts: [(quarter turns, extra steps)]: each copy is the
    pattern turned by quarter x 90 degrees (exact) and run settle + extra steps.
    Returns a dict: A (n, 3, H, W) torch, kin (Kin), com (n, 2) unwrapped, axis (n,) angles,
    pre (n, PRE, 5) the stats [mR, mG, mB, speed, |omega|] of the PRE steps before,
    mass (n, settle+max extra) total mass series, phase (n,) steps since the last mass peak.
    The book's convention: clockwise = +, and every morph starts turning +: a pattern that settles
    turning anticlockwise is mirrored (mirror=None: decided once per morph by a probe run)."""
    t = sim.t
    if mirror is None:
        key = (code, morph, R)
        if key not in _MIRROR:
            _MIRROR[key] = start_set(sim, code, morph, [(0, 0)], settle, R, mirror=False)['omega_end'][0] < -OMEGA_SIDE
        mirror = _MIRROR[key]
    cells = pattern(code, morph, R)
    if mirror: cells = cells[:, :, ::-1].copy()
    rots = sorted({q % 4 for q, _ in starts})
    A = t.tensor(np.array([place(sim, cells, q) for q in rots]), dtype=t.float32, device=sim.dev)
    nr = len(rots)
    total = settle + max(e for _, e in starts)
    meas = sim.measure(A)
    com = meas[:, :2].copy()
    kin = Kin(nr, R, sim.dt)
    hist = np.zeros((nr, total, 5)); mass = np.zeros((nr, total)); vel = np.zeros((nr, total, 2)); om = np.zeros((nr, total))
    snaps = {}
    for s in range(total):
        A = sim.step(A)
        m2 = sim.measure(A)
        d = np.stack([wrapd(m2[:, 0] - meas[:, 0], sim.W), wrapd(m2[:, 1] - meas[:, 1], sim.H)], 1)
        com += d; meas = m2
        kin.update(d / sim.dt)
        hist[:, s, :3] = m2[:, 2:5] / R ** 2
        hist[:, s, 3] = kin.speed / R; hist[:, s, 4] = np.abs(kin.omega); om[:, s] = kin.omega
        mass[:, s] = m2[:, 2:5].sum(1); vel[:, s] = d
        for i, (q, e) in enumerate(starts):
            if s + 1 == settle + e:
                r = rots.index(q % 4)
                snaps[i] = (A[r].clone(), kin.copy_rows([r]), com[r].copy(), m2[r].copy(), r, s + 1)
    out = {'A': t.stack([snaps[i][0] for i in range(len(starts))])}
    K = Kin(len(starts), R, sim.dt)
    for a in ('age', 'omega', 'omegaAge', 'speed', 'velAng', 'fast', 'vel', 'fastAng', 'slow'):
        setattr(K, a, np.concatenate([getattr(snaps[i][1], a) for i in range(len(starts))]))
    out['kin'] = K
    out['com'] = np.array([snaps[i][2] for i in range(len(starts))])
    M = np.array([snaps[i][3] for i in range(len(starts))])
    out['axis'] = np.arctan2(wrapd(M[:, 6] - M[:, 8], sim.H), wrapd(M[:, 5] - M[:, 7], sim.W))
    out['pre'] = np.array([hist[snaps[i][4], snaps[i][5] - PRE:snaps[i][5]] for i in range(len(starts))])
    out['period'] = cycle_period(mass[0, settle // 2:], vel[0, settle // 2:, 0], vel[0, settle // 2:, 1])
    # phase: steps since the cycle's main mass peak (the largest in the period before settling ends)
    P = out['period']
    ph = []
    for i in range(len(starts)):
        r, n = snaps[i][4], snaps[i][5]
        if P <= 0: ph.append(0); continue   # (no cycle found)
        w0 = max(0, settle - P)   # (a cycle longer than the settling: the peak of what there is)
        p0 = w0 + int(np.argmax(mass[r, w0:settle]))
        ph.append((n - 1 - p0) % P)
    out['phase'] = np.array(ph)
    out['omega_end'] = om[:, max(0, settle - 100):settle].mean(1)
    out['mirrored'] = bool(mirror)
    return out


def made_start_set(sim, code, morph, made, starts, settle, R, models):
    """Starts of a morph made by an action from another (made: 'parent:ch:spot:hold:wait', e.g.
    HGMLH8 'paratardus:R:front-left edge:90:40': a short-lived morph whose captured pattern wanders
    off, made fresh as the game makes it): the parent's starts, the mask held at the spot (the
    parent model's offset) in its frame, then wait steps; as start_set's dict, only the starts
    named morph at the end ('kept': their indices)."""
    t = sim.t
    parent, ch, spot, hold, wait = made.split(':')
    hold, wait = int(hold), int(wait)
    book = json.load(open(os.path.join(BOOK, f'{code}.json')))
    pm = book['morphs'][parent]
    off = next(e['off'] for e in pm['actions'] if e['spot'] == spot)
    S = start_set(sim, code, parent, starts, settle, R)
    A, K, com = S['A'], S['kin'], S['com'].copy()
    n = A.shape[0]
    meas = sim.measure(A)
    onehot = t.nn.functional.one_hot(t.full((n,), CH.index(ch), device=sim.dev), 3).float()[:, :, None, None]
    total = hold + wait
    hist = np.concatenate([S['pre'], np.zeros((n, total, 5))], 1)
    om = np.zeros((n, total))
    for s in range(total):
        if s < hold:
            fa = np.where(np.isnan(K.slow), S['axis'], K.slow) if pm['frame'] == 'vel' else S['axis']
            c, sn = np.cos(fa), np.sin(fa)
            dx = off[0] * R * c - off[1] * R * sn; dy = off[0] * R * sn + off[1] * R * c
            cx = t.tensor(meas[:, 0] + dx, dtype=t.float32, device=sim.dev); cy = t.tensor(meas[:, 1] + dy, dtype=t.float32, device=sim.dev)
            disc = sim.disc(cx, cy, t.full((n,), float(MASK_R), device=sim.dev))[:, None] * onehot
            A = sim.step(A, 1 - disc)
        else:
            A = sim.step(A)
        m2 = sim.measure(A)
        d = np.stack([wrapd(m2[:, 0] - meas[:, 0], sim.W), wrapd(m2[:, 1] - meas[:, 1], sim.H)], 1)
        com += d; meas = m2
        K.update(d / sim.dt)
        j = PRE + s
        hist[:, j, :3] = m2[:, 2:5] / R ** 2; hist[:, j, 3] = K.speed / R; hist[:, j, 4] = np.abs(K.omega); om[:, s] = K.omega
    ser = {k: hist[:, :, i] for i, k in enumerate(('mR', 'mG', 'mB', 'speed', 'omega'))}
    named = named_series(ser, models)
    kept = [i for i in range(n) if named[i, -1] == morph]
    print(f'  made from {parent} ({ch} {spot} {hold}, then {wait} steps): {len(kept)}/{n} named {morph} '
          f'({dict(Counter(x or "???" for x in named[:, -1]))})', flush=True)
    out = {'A': A[kept], 'com': com[kept], 'pre': hist[kept, -PRE:], 'kept': kept,
           'axis': np.arctan2(wrapd(meas[kept, 6] - meas[kept, 8], sim.H), wrapd(meas[kept, 5] - meas[kept, 7], sim.W)),
           'period': S['period'], 'phase': S['phase'][kept], 'omega_end': om[kept, -min(100, total):].mean(1), 'mirrored': S['mirrored']}
    Kk = Kin(len(kept), R, sim.dt)
    for at in ('age', 'omega', 'omegaAge', 'speed', 'velAng', 'fast', 'vel', 'fastAng', 'slow'):
        setattr(Kk, at, getattr(K, at)[kept].copy())
    out['kin'] = Kk
    return out


def start_states(sim, code, morph, reps, settle, R, models=None, debris=None):
    """(for the clone sweep) reps copies: turned k x 90 degrees, settled settle + 6k steps"""
    s = start_set(sim, code, morph, [(k % 4, 6 * k) for k in range(reps)], settle, R)
    return s['A'], s['kin'], s['com'], s['axis']


def stack_kin(kins, idx, R, dt):
    K = Kin(len(idx), R, dt)
    for a in ('age', 'omega', 'omegaAge', 'speed', 'velAng', 'fast', 'vel', 'fastAng', 'slow'):
        setattr(K, a, getattr(kins, a)[idx].copy())
    return K


def named_series(series, models, pre=None, device='cuda'):
    """What detect.js would name at each step, from the whole board's stats series (one body):
    (N, T) array of names ('' = ???), windows of the WIN newest samples (none before WIN)."""
    import torch
    X = np.stack([series['mR'], series['mG'], series['mB'], series['speed'], np.abs(series['omega'])], -1)
    P = 0
    if pre is not None:
        P = pre.shape[1]
        X = np.concatenate([pre, X], 1)
    N, T, _ = X.shape
    Z = torch.tensor(feat(np.maximum(X, 0)), dtype=torch.float32, device=device).reshape(N * T, 5)
    out = np.full((N, T), '', object)
    best_known = np.full((N, T), -np.inf)
    for level, model in models.items():
        Ls = []
        for mp in model:
            ref = torch.tensor(mp['ref'], dtype=torch.float32, device=device)
            inv = torch.tensor(mp['inv'], dtype=torch.float32, device=device)
            sc = torch.empty(N * T, device=device)
            for c0 in range(0, N * T, 200000):
                d = (((Z[c0:c0 + 200000, None, :] - ref[None]) * inv) ** 2).sum(2)
                sc[c0:c0 + 200000] = torch.log(torch.exp(-d / 2).mean(1) + 1e-30).clamp_min(CLAMP)
            sc = sc.reshape(N, T)
            cs = torch.cumsum(sc, 1)
            L = torch.full_like(sc, -99.0)
            L[:, WIN - 1:] = (cs[:, WIN - 1:] - torch.cat([torch.zeros(N, 1, device=device), cs[:, :-WIN]], 1)) / WIN
            Ls.append(L)
        L = torch.stack(Ls, -1)                                    # (N, T, M)
        Nn = L - torch.tensor([mp['logNorm'] for mp in model], device=device)
        b = Nn.argmax(-1)
        e = torch.exp(NEFF * (Nn - Nn.max(-1, keepdim=True).values))
        p = (e.max(-1).values / e.sum(-1)).cpu().numpy()
        Lb = L.gather(-1, b[..., None])[..., 0].cpu().numpy()
        b = b.cpu().numpy()
        known = (p >= P_MIN) & (Lb >= FIT_MIN)
        names = np.array([mp['name'] for mp in model], object)[b]
        take = known & (Lb > best_known)
        out[take] = names[take]
        best_known[take] = Lb[take]
    out[:, :WIN - 1] = ''
    return out[:, P:]


def recovery(named, label, start):
    """steps from start until the series names `label` and keeps naming it to the end; None if never"""
    if not label or label in ('???', 'dissolved', 'exploded', 'multiplied') or ' + ' in label: return None
    bad = np.nonzero(named[start:] != label)[0]
    if len(bad) == 0: return 0
    r = int(bad[-1]) + 1
    return r if start + r < len(named) else None


def circmean(a, axis=-1):
    return np.arctan2(np.sin(a).mean(axis), np.cos(a).mean(axis))


# ---- the steer experiment (masks on one morph) -----------------------------------------------
DIRS = ['front', 'front-right', 'right', 'tail-right', 'tail', 'tail-left', 'left', 'front-left']
DISTS = {'body': 0.45, 'edge': 0.9, 'outside': 1.3}


def morph_extent(code, morph, R):
    """the pattern's width at R: the larger side of its cells above 0.05 (any channel)"""
    c = pattern(code, morph, R).sum(0)
    ys, xs = np.nonzero(c > 0.05)
    return int(max(ys.max() - ys.min(), xs.max() - xs.min()) + 1)


def named_spots(R, dists=('body', 'edge'), scale=1.0):
    """the nameable spots: centre, and 8 directions at each distance named (body 0.45 R, edge 0.9 R,
    outside 1.3 R; x scale for a larger morph, e.g. a twin); (name, (fwd, side) cells)"""
    out = [('centre', (0.0, 0.0))]
    for dn, dd in ((d, DISTS[d] * scale) for d in dists):
        for k, nm in enumerate(DIRS):
            a = math.radians(45 * k)          # 0 = ahead, 90 = right
            out.append((f'{nm} {dn}', (round(dd * R * math.cos(a), 2), round(dd * R * math.sin(a), 2))))
    return out


def sequence(named, min_run=5):
    """the name box over time, as a list of names ('???' for none), runs shorter than min_run dropped"""
    runs = []
    for v in named:
        v = v or '???'
        if runs and runs[-1][0] == v: runs[-1][1] += 1
        else: runs.append([v, 1])
    runs = [r for r in runs if r[1] >= min_run] or runs[-1:]
    out = []
    for v, _ in runs:
        if not out or out[-1] != v: out.append(v)
    return out


def cmd_steer(a):
    import torch
    R = a.R
    rule = rule_of(a.code)
    models = detect_models(a.code)
    sims = {}
    base = a
    for morph in base.morph.split(','):
        t0 = time.time()
        a = argparse.Namespace(**vars(base))
        # --like: the morph's sweep again as its model's (steer-<morph>.json): the same starts and
        # settings, other channels (e.g. RB), saved under a tag; the book merges it into the model
        like = None
        if a.like:
            ld = json.load(open(os.path.join(OUT, a.code, f'steer-{morph}.json')))
            like = ld['meta']
            a.settle, a.observe, a.durs = like['settle'], like['observe'], like['durs']
            rb = next((r for r in ld['runs'] if r['spot'] == 'right body'), None)   # (the oldest sweeps did not record their scale: from a spot's place)
            a.spot_scale = like.get('spot_scale') or (round(math.hypot(*rb['off']) / (DISTS['body'] * R), 2) if rb else 1.0)
            a.size = like.get('size', 128)
            a.dists = ','.join(sorted({n.split(' ')[-1] for n in like['spots'] if n != 'centre'}, key=['body', 'edge', 'outside'].index))
            a.period, a.phase_signal, a.age, a.states = like.get('swept_period', 0) or 0, like.get('phase_signal', 'omega'), like.get('age_sweep', False), 0
            a.tag = a.tag or 'more'
        # auto (0): the spots' scale by the morph's width against a 5N7KKM single's (29 cells at R 18),
        # so 'edge' is at its edge; the board by its width (room for the spots outside and a turn)
        ext = morph_extent(a.code, morph, R)
        spot_scale = a.spot_scale if a.spot_scale > 0 else float(np.clip(round(ext / (29 * R / 18), 1), 0.7, 2.5))
        size = a.size if a.size > 0 else (128 if ext <= 40 else 160 if ext <= 52 else 192)
        if size not in sims: sims[size] = Sim(rule, size, size, R)
        sim = sims[size]
        a.spot_scale_now, a.size_now = spot_scale, size
        print(f'{a.code} {morph}: {ext} cells wide, spots x {spot_scale}, board {size}', flush=True)
        made = a.made or (like or {}).get('made')
        probe = start_set(sim, a.code, made.split(':')[0] if made else morph, [(0, 0)], a.settle + 240, R)   # (made: the parent's cycle)
        period = a.period or probe['period']   # (--period: an oscillator's swing, e.g. HAESRE zigzag's 115, longer than its mass pulse)
        offs = sorted({round(k * period / a.phases) for k in range(a.phases)})    # distinct phases
        # other headings: settled longer by whole cycles (a turning morph has turned meanwhile), and
        # turned by a quarter (exact, so only a still morph's frame differs from the grid's)
        starts = [(q, (0 if a.age or a.no_shift else q * 13 * period) + o) for q in range(a.angles) for o in offs]   # (--age, --no-shift: only turned, not settled longer: a short-lived or unstable morph)   # (--age: a short-lived morph, at each age, only turned)
        states = None
        if like:   # (its model's sweep's starts, cycle and states, exactly)
            starts, period = [tuple(x) for x in like['starts']], like['period']
            if like.get('states'): states = {k: v for k, v in like['states'].items() if k != 'start_state'}
        if a.states:   # a chaotic morph: starts drawn from each of its states in a long run (see kin_states)
            states, starts = kin_states(sim, a.code, morph, a.settle, a.probe, R, a.states, a.phases * a.angles)
        if like and like.get('made'): a.made = like['made']
        if a.made:   # (made fresh from its parent by an action: the parent's starts, those that became it)
            S = made_start_set(sim, a.code, morph, a.made, starts, a.settle, R, models)
            starts = [starts[i] for i in S['kept']]
            if not like: period = S['period']
        else:
            S = start_set(sim, a.code, morph, starts, a.settle, R)
        A0, K0, C0, AX0 = S['A'], S['kin'], S['com'], S['axis']
        if states:   # each start's state as its own kinematics show it (a chaotic run need not repeat exactly)
            states['start_state'] = [state_of(states, K0.speed[i] / R, K0.omega[i]) for i in range(len(starts))]
            print(f'  starts per state: {np.bincount(states["start_state"], minlength=a.states).tolist()}', flush=True)
        moving = bool(np.median(K0.speed / R) > SPEED_MIN_R * 1.5)
        frame = 'vel' if moving else 'axis'
        spots = named_spots(R, a.dists.split(','), spot_scale) if a.spots == 'named' else [
            (f'{f},{s}', (f, s)) for f in range(-40, 41, a.spacing) for s in range(-40, 41, a.spacing) if math.hypot(f, s) <= a.extent * R]
        if a.extra_spots:   # (named spots of one's own, (fwd, side) in R, e.g. XEH4YR rectus's sweet spot for cor)
            extra = [(nm.strip(), tuple(float(v) * R for v in xy.split(','))) for nm, xy in (x.split(':') for x in a.extra_spots.split(';'))]
            spots = extra if a.only_extra else spots + extra
        runs = []
        for st in range(len(starts)):
            runs.append(dict(start=st, spot='-', action='none', ch=0, r=float(MASK_R), dur=0, off=(0.0, 0.0)))
            for ch, dur, (nm, off) in itertools.product(a.channels, a.durs, spots):
                runs.append(dict(start=st, spot=nm, action='mask', ch=CH.index(ch), r=float(MASK_R), dur=dur, off=off))
        runs.sort(key=lambda r: (r['dur'], r['start']))
        steps = max(a.durs) + a.observe
        print(f'{a.code} {morph}: frame {frame}, cycle {period} steps, {len(starts)} starts, {len(spots)} spots, '
              f'{len(runs)} trials of {steps} + {OUT_WIN} steps', flush=True)
        res, series = [], {k: [] for k in SERIES}
        for b0 in range(0, len(runs), a.batch):
            chunk = runs[b0:b0 + a.batch]
            idx = [r['start'] for r in chunk]
            plan = {k: [r[k] for r in chunk] for k in ('action', 'ch', 'r', 'dur', 'off')}
            plan['frame'] = [frame] * len(chunk)
            plan['axis0'] = AX0[idx]
            ser, out, _ = run_batch(sim, A0[idx], stack_kin(K0, idx, R, sim.dt), C0[idx], plan, steps, R, models, a.debris)
            for r, o in zip(chunk, out): res.append({**r, **o})
            for k in series: series[k].append(ser[k])
            el = time.time() - t0
            print(f'  {len(res)}/{len(runs)} trials, {el:.0f} s, ~{el / len(res) * (len(runs) - len(res)):.0f} s left', flush=True)
        series = {k: np.concatenate(v) for k, v in series.items()}
        # the name box over time (with the stats from before the action as its history)
        named = named_series(series, models, pre=S['pre'][[r['start'] for r in res]])
        ctrl = {r['start']: i for i, r in enumerate(res) if r['action'] == 'none'}
        for i, r in enumerate(res):
            r['seq'] = sequence(named[i])
            r['named_end'] = named[i, -1] or '???'
            r['recover'] = recovery(named[i], r['named_end'], r['dur'])
            r['phase'] = int(S['phase'][r['start']])
            om = series['omega'][i, -100:]
            r['omega_end'] = round(float(om.mean()), 4)
        for i, r in enumerate(res):
            c = ctrl[r['start']]
            h0 = float(series['head'][c, 0]) if frame == 'vel' else float(AX0[r['start']])
            when = r['dur'] + max(r['recover'] if r['recover'] is not None else a.short, 30) + 5
            for tag, tt in (('rec', when), ('long', steps - 1)):
                tt = min(tt, steps - 1)
                w = slice(max(0, tt - 10), tt + 1)
                if frame == 'vel':
                    dh = circmean(series['velAng'][i, w]) - circmean(series['velAng'][c, w])
                    r[f'turn_{tag}'] = round(math.degrees(wrapa(dh)), 1)
                    r[f'speed_{tag}'] = round(float(series['speed'][i, w].mean() / max(series['speed'][c, w].mean(), 1e-6)), 3)
                dx = series['x'][i, w].mean() - series['x'][c, w].mean()
                dy = series['y'][i, w].mean() - series['y'][c, w].mean()
                cs, sn = math.cos(h0), math.sin(h0)
                r[f'shift_{tag}'] = [round(float(dx * cs + dy * sn) / R, 2), round(float(-dx * sn + dy * cs) / R, 2)]
        d = os.path.join(OUT, a.code)
        os.makedirs(d, exist_ok=True)
        meta = {'code': a.code, 'morph': morph, 'kind': 'steer', 'R': R, 'dt': sim.dt, 'size': size, 'frame': frame, 'spot_scale': spot_scale,
                'period': period, 'starts': starts, 'phases': a.phases, 'angles': a.angles, 'settle': a.settle,
                'observe': a.observe, 'steps': steps, 'durs': a.durs, 'spots': [n for n, _ in spots],
                'channels': a.channels, 'control': [res[ctrl[k]] for k in sorted(ctrl)],
                'speed': float(np.median(K0.speed / R)), 'mirrored': S['mirrored'],
                'omega': float(np.median([res[ctrl[k]]['omega_end'] for k in ctrl])), 'seconds': round(time.time() - t0),
                # each start's motion as the game's kinematics see it (an oscillator's phase from them:
                # the arrow's direction off the heading-up heading, and the turn rate)
                'start_kin': [[float(K0.velAng[i]), float(K0.slow[i]), float(K0.omega[i]), float(K0.speed[i] / R)] for i in range(len(starts))],
                'phase_signal': a.phase_signal, 'swept_period': a.period, 'age_sweep': bool(a.age), 'like': bool(like), 'made': a.made or None, **({'states': states} if states else {})}
        nm = morph + (f'-{a.tag}' if a.tag else '')   # (--tag: an experiment's own files, not the model's)
        json.dump({'meta': meta, 'runs': res}, open(os.path.join(d, f'steer-{nm}.json'), 'w'), default=jdefault)
        np.savez_compressed(os.path.join(d, f'steer-{nm}-series.npz'), named=named.astype(str),
                            **{k: v for k, v in series.items()})
        np.save(os.path.join(d, f'steer-{nm}-start.npy'), A0.cpu().numpy().astype(np.float16))
        print(f'  saved {d}/steer-{nm}.json ({time.time() - t0:.0f} s)', flush=True)
        torch.cuda.empty_cache()


# ---- the clone experiment (a clone placed beside a ship) --------------------------------------
def cmd_clone(a):
    import torch
    t = torch
    R = a.R
    rule = rule_of(a.code)
    sim = Sim(rule, a.size, a.size, R)
    models = detect_models(a.code)
    for pair in a.pair.split(';'):
        ship, cl = pair.split(',')
        t0 = time.time()
        A0, K0, C0, AX0 = start_states(sim, a.code, ship, a.reps, a.settle, R, models, a.debris)
        moving = bool(np.median(K0.speed / R) > SPEED_MIN_R * 1.5)
        frame = 'vel' if moving else 'axis'
        F0 = np.where(moving & ~np.isnan(K0.slow), K0.slow, AX0) if moving else AX0
        cells = pattern(a.code, cl, R)
        runs = []
        for rep in range(a.reps):
            runs.append(dict(rep=rep, dist=0.0, bearing=0, turn=0, flip=False, none=True))
            for d, bear, turn, flip in itertools.product(a.dists, range(0, 360, 360 // a.bearings), range(0, 360, 45), (False, True)):
                runs.append(dict(rep=rep, dist=d, bearing=bear, turn=turn, flip=flip, none=False))
        # Each clone's own axis angle once flipped, for the turn convention:
        # the clone's body axis points at frame + turn (degrees, + clockwise as seen).
        ax_n, ax_f = axis_angle(cells), axis_angle(transform(cells, 0, True))
        steps = a.observe
        print(f'{a.code} {ship} + clone {cl}: frame {frame}, {len(runs)} runs of {steps} + {OUT_WIN} steps', flush=True)
        res = []
        for b0 in range(0, len(runs), a.batch):
            chunk = runs[b0:b0 + a.batch]
            boards = A0[[r['rep'] for r in chunk]].clone()
            com = C0[[r['rep'] for r in chunk]].copy()
            mnow = sim.measure(boards)
            for j, r in enumerate(chunk):
                if r['none']: continue
                fr = F0[r['rep']]
                ang = fr + math.radians(r['turn']) - (ax_f if r['flip'] else ax_n)
                f = transform(cells, ang, r['flip'])
                bang = fr + math.radians(r['bearing'])     # bearing: 0 = ahead, 90 = right
                px = mnow[j, 0] + r['dist'] * R * math.cos(bang)
                py = mnow[j, 1] + r['dist'] * R * math.sin(bang)
                S = f.shape[1]
                ys = (np.arange(S) + int(round(py)) - S // 2) % sim.H
                xs = (np.arange(S) + int(round(px)) - S // 2) % sim.W
                ft = t.tensor(f, dtype=t.float32, device=sim.dev)
                sub = boards[j][:, ys][:, :, xs]
                boards[j][:, ys[:, None], xs[None, :]] = t.maximum(sub, ft)
            plan = {'action': ['none'] * len(chunk), 'ch': [0] * len(chunk), 'r': [0.0] * len(chunk), 'dur': [0] * len(chunk),
                    'off': [(0, 0)] * len(chunk), 'frame': [frame] * len(chunk), 'axis0': AX0[[r['rep'] for r in chunk]]}
            kin = stack_kin(K0, [r['rep'] for r in chunk], R, sim.dt)
            ser, out, _ = run_batch(sim, boards, kin, com, plan, steps, R, models, a.debris)
            for r, o in zip(chunk, out): res.append({**r, **o})
            el = time.time() - t0
            print(f'  {len(res)}/{len(runs)} runs, {el:.0f} s, ~{el / len(res) * (len(runs) - len(res)):.0f} s left', flush=True)
        d = os.path.join(OUT, a.code)
        os.makedirs(d, exist_ok=True)
        meta = {'code': a.code, 'ship': ship, 'clone': cl, 'kind': 'clone', 'R': R, 'size': a.size, 'frame': frame,
                'reps': a.reps, 'settle': a.settle, 'observe': a.observe, 'dists': a.dists, 'bearings': a.bearings,
                'seconds': round(time.time() - t0)}
        json.dump({'meta': meta, 'runs': res}, open(os.path.join(d, f'clone-{ship}+{cl}.json'), 'w'), default=jdefault)
        np.save(os.path.join(d, f'clone-{ship}+{cl}-start.npy'), A0.cpu().numpy().astype(np.float16))
        print(f'  saved {d}/clone-{ship}+{cl}.json ({time.time() - t0:.0f} s)', flush=True)
        torch.cuda.empty_cache()


# ---- the morph model ---------------------------------------------------------------------------
def steps_s(n): return f'{n} steps ({n / 30:.2g} s)'


def rec_s(e):
    """steps after the hold until the morph is named again (the ship is '???', not steerable, till then)"""
    return f"{e['recover']} (≤{e['recover_max']})" if 'recover' in e else 'not always'


def where(off, R):
    """(fwd, side) cells -> words"""
    f, s = off[0] / R, off[1] / R
    if abs(f) < 0.12 and abs(s) < 0.12: return 'centre'
    parts = []
    if abs(f) >= 0.12: parts.append(('front' if f > 0 else 'tail'))
    if abs(s) >= 0.12: parts.append(('right' if s > 0 else 'left'))
    return '-'.join(parts)


def group_runs(runs, keys):
    g = {}
    for r in runs: g.setdefault(tuple(r[k] if not isinstance(r[k], list) else tuple(r[k]) for k in keys), []).append(r)
    return g


DEAD = ('dissolved', 'exploded', 'multiplied')
OMEGA_SIDE = 0.01     # rad per unit time: below, a morph has no turning direction


def side(om): return 0 if abs(om) < OMEGA_SIDE else (1 if om > 0 else -1)


def effect(r, base, chir):
    """what one trial did, in the book's words"""
    if r.get('replicated'): return f"replicate: {r['replicated']}"
    if r['label'] in DEAD: return 'dissolve'
    if ' + ' in r['label']:   # several bodies at the end (the name box follows only one of them)
        parts = r['label'].split(' + ')
        # copies of one morph, with the mass grown (not a twin coming apart): it replicated
        if '???' not in parts and len(set(parts)) == 1 and r.get('mass_rel', 1) >= 1.5: return f'replicate: {parts[0]}'
        return 'split: ' + r['label']
    end = r['named_end'] if r['named_end'] != '???' else r['label']
    if end in ('???',) or ' + ' in end: return 'unsettled' if end == '???' else 'split: ' + end
    sd = side(r['omega_end'])
    if end == base:
        return 'reverse' if chir and sd == -chir else 'keep'
    if chir and sd: return f'into {end}' + (' (same turn)' if sd == chir else ' (reversed)')
    return f'into {end}'


def phase_words(ok, bad, period):
    """when it works: the phases (steps since the cycle's mass peak) where it always does"""
    if not bad: return 'any phase'
    if not ok: return 'no phase'
    always = sorted(set(ok) - set(bad))
    some = sorted(set(ok) & set(bad))
    w = f"phase {','.join(map(str, always))} of {period}" if always else ''
    if some: w += (' (sometimes ' if w else 'sometimes at ') + ','.join(map(str, some)) + (')' if w else '')
    return w


def load_ends(path):
    """a sweep's series, only the last 10 steps of what end_feat reads (each array decompressed once)"""
    z = np.load(path)
    out = {k: np.asarray(z[k])[:, -10:] for k in ('mR', 'mG', 'mB', 'speed', 'omega')}
    m = np.asarray(z['mass'])
    out['peak'] = m.max(1) / np.maximum(m[:, 0], 1e-9)   # (the most mass at any time, of the start's: a replication swells first)
    return out


def end_feat(z, i):
    """a trial's end in the detector's features (log mR, mG, mB, speed, |omega|), its last 10 steps"""
    return feat([float(z[k][i, -10:].mean()) for k in ('mR', 'mG', 'mB')] + [float(z['speed'][i, -10:].mean()), float(np.abs(z['omega'][i, -10:]).mean())])


def control_refs(d):
    """each swept morph's own settled end (its untouched copies): name -> features"""
    refs = {}
    for f in os.listdir(d):
        mt = re.match(r'steer-([a-z+]+)-series\.npz$', f)
        if not mt: continue
        z = load_ends(os.path.join(d, f))
        runs = json.load(open(os.path.join(d, f'steer-{mt.group(1)}.json')))['runs']
        ii = [i for i, r in enumerate(runs) if r['action'] == 'none']
        if ii and len(runs) == len(z['mB']): refs[mt.group(1)] = np.mean([end_feat(z, i) for i in ii], 0)
    return refs


PHASE_BINS = 12   # an oscillator's swing cut into this many phases
PHASE_GAIN_MIN, TURN_SPREAD_MIN = 0.05, 10.0   # the phase is worth telling if it predicts the outcomes this much better than chance, or spreads the turns this much (deg)


def swing_phases(start_kin, starts, period, dt, signal='omega'):
    """an oscillator's phase at each start (0 .. 2 pi) as the game can tell it: the phase of its turn
    rate (core kin.omega, a sine over the swing): omega ~ sin(phase). The starts are spread over the
    swing (offsets o), so a sine fitted to their turn rates against o gives it. models.js swingBin
    fits the same sine to the ship's recent turn rates."""
    w0 = 2 * math.pi / (period * dt)
    t = np.array([o % period for _, o in starts]) * dt
    om = np.array([k[2] if signal == 'omega' else k[3] for k in start_kin])   # (the turn rate, or the speed)
    X = np.stack([np.sin(w0 * t), np.cos(w0 * t), np.ones_like(t)], 1)
    a, b, _ = np.linalg.lstsq(X, om, rcond=None)[0]    # om ~ a sin(w0 t) + b cos(w0 t) = A sin(w0 t + phi)
    return (w0 * t + math.atan2(b, a)) % (2 * math.pi)


REP_PEAK = 2.0   # a trial whose mass ever reached this many times its start replicated (PJ1EF7 rectus: 2.8-7.4; into tachus 1.2), however it ended


REF_TOL = np.array([0.15, 0.15, 0.15, 1.5, 1.5])   # an unnamed end within this of a morph's own settled end (log units) is that morph


def steer_model(data, a=None, refs=None, z=None):
    meta, runs = data['meta'], data['runs']
    # starts that were not the morph to begin with (an unstable one turned into another while settling,
    # e.g. PJ1EF7 lithos into tachus): their trials left out (their untouched copy's first name)
    bad = {r['start'] for r in runs if r['action'] == 'none' and r['seq'] and r['seq'][0] != meta['morph']}
    if bad and len(bad) < len({r['start'] for r in runs}):
        keep = [i for i, r in enumerate(runs) if r['start'] not in bad]
        print(f"  {meta['morph']}: {len(bad)} starts were not {meta['morph']} at first: left out", flush=True)
        runs = [runs[i] for i in keep]
        if z is not None: z = {k: v[keep] for k, v in z.items()}
        data = {**data, 'runs': runs}
    # an unnamed end (???) that is a swept morph's own settled state (e.g. KH97WU dot, named in the
    # lab but not by its recorded stats here): that morph
    if refs and z is not None and len(z['mB']) == len(runs):
        for i, r in enumerate(runs):
            if r['named_end'] == '???' and r['label'] == '???':
                x = end_feat(z, i)
                hit = [n for n, ref in refs.items() if np.all(np.abs(x - ref) <= REF_TOL)]
                if hit: r['named_end'] = r['label'] = min(hit, key=lambda n: np.abs(x - refs[n]).sum())
    # a replication: the mass swelled to REP_PEAK times (its copies may merge, fade or collide after);
    # the copies: the morph most named at the ends of the morph's replications
    if z is not None and 'peak' in z and len(z['peak']) == len(runs):
        # (not one that grew into a heavier morph and stays it, e.g. VT049W rectus -> paragyrans: it
        # ends as one named body holding the swollen mass)
        grown = lambda i, r: ' + ' not in r['label'] and r['label'] not in DEAD and r['label'] != '???' and r.get('mass_rel', 1) >= 0.8 * z['peak'][i]
        reps = [r for i, r in enumerate(runs) if r['action'] != 'none' and z['peak'][i] >= REP_PEAK and not grown(i, r)]
        names = [p for r in reps for p in r['label'].split(' + ') if p not in ('???',) and p not in DEAD]
        copy = mode(names) if names else '???'
        for r in reps: r['replicated'] = copy
    R, frame, morph, period = meta['R'], meta['frame'], meta['morph'], meta.get('period', 0)
    ctrl = [r for r in runs if r['action'] == 'none'] or meta['control']   # (as named above)
    base = mode(c['named_end'] for c in ctrl)
    chir = side(np.median([c['omega_end'] for c in ctrl]))
    groups = group_runs([r for r in runs if r['action'] != 'none'], ('ch', 'spot', 'dur'))
    # an oscillator swept over its swing (--period, start_kin): each start's phase bin, as the game will tell it
    sbin, nbins = None, PHASE_BINS
    if meta.get('states'):   # a chaotic morph: its states
        sbin, nbins = meta['states']['start_state'], len(meta['states']['centroids'])
    elif meta.get('age_sweep'):   # a short-lived morph: its age (steps since it formed, settle + offset), in bins over the swept span
        nbins = meta['phases']
        sbin = [min(nbins - 1, round((o % period) / period * nbins)) for _, o in meta['starts']]   # (the offsets are at k / nbins of it)
    elif meta.get('start_kin') and period and meta.get('swept_period', period):   # (swept over a cycle: --period; HAESRE zigzag's sweep predates the key)
        ph = swing_phases(meta['start_kin'], meta['starts'], period, meta['dt'], meta.get('phase_signal', 'omega'))
        sbin = [int(p / (2 * math.pi) * PHASE_BINS) % PHASE_BINS for p in ph]
    rows = []
    for (ch, spot, dur), rs in groups.items():
        effs = [effect(r, base, chir) for r in rs]
        main = mode(effs)
        mine = [r for r, e in zip(rs, effs) if e == main]
        seqs = [' > '.join(r['seq']) for r in mine]
        rec = [r['recover'] for r in mine]
        e = {'ch': CH[ch], 'spot': spot, 'hold': dur, 'off': [round(rs[0]['off'][0] / R, 2), round(rs[0]['off'][1] / R, 2)],
             'effect': main, 'p': round(len(mine) / len(rs), 2), 'n': len(rs),
             'effects': {x: effs.count(x) for x in sorted(set(effs))},
             'seq': mode(seqs), 'always_named': all(r['seq'] == [base] for r in mine),
             'unnamed_max': (max(v for v in rec if v is not None) if any(v is not None for v in rec) else None),
             'phases': phase_words([r['phase'] for r in mine], [r['phase'] for r in rs if r not in mine], period)}
        # the runner-up outcome (in at least 2 trials, a tenth of them): [effect, p]
        others = sorted(((effs.count(x), x) for x in set(effs) if x != main), reverse=True)
        if others and others[0][0] >= max(2, 0.1 * len(rs)): e['alt'] = [others[0][1], round(others[0][0] / len(rs), 2)]
        if main.startswith('replicate: '):   # how many copies (bodies at the end)
            nb = [len(r['label'].split(' + ')) if r['label'] not in DEAD else 0 for r in mine]
            e['copies'] = [min(nb), max(nb)]
            # what it becomes first, the short-lived morph that replicates (VT049W: gyrans), if named on the way
            cp = main[11:]
            via = [next((n for n in r['seq'][1:] if n not in (base, '???', cp)), None) for r in mine]
            via = [v for v in via if v]
            if len(via) >= len(mine) / 2: e['via'] = mode(via)
        if frame == 'vel' and main in ('keep', 'reverse'):
            tv = np.radians([r['turn_rec'] for r in mine])
            e['turn'] = round(math.degrees(circmean(tv)), 1)
            e['turn_sd'] = round(float(np.degrees(np.sqrt(-2 * np.log(max(1e-9, abs(np.exp(1j * tv).mean())))))), 1)
        if sbin is not None:   # per phase bin: the main outcome there, its share, its turn
            bp = []
            for b in range(nbins):
                xs = [(r, x) for r, x in zip(rs, effs) if sbin[r['start']] == b]
                if not xs: bp.append(None); continue
                ex = [x for _, x in xs]
                mb = mode(ex)
                o = [mb, round(ex.count(mb) / len(ex), 2), len(ex)]
                oth = sorted(((ex.count(y), y) for y in set(ex) if y != mb), reverse=True)
                e.setdefault('alt_by_phase', [None] * nbins)[b] = [oth[0][1], round(oth[0][0] / len(ex), 2)] if oth and oth[0][0] >= 2 else None
                if frame == 'vel' and mb in ('keep', 'reverse'):
                    tv = np.radians([r['turn_rec'] for r, x in xs if x == mb])
                    o.append(round(math.degrees(circmean(tv)), 1))
                    o.append(round(float(np.degrees(np.sqrt(-2 * np.log(max(1e-9, abs(np.exp(1j * tv).mean())))))), 1))   # (its spread there)
                bp.append(o)
            e['by_phase'] = bp
        e['_effs'] = [(r['phase'] if sbin is None else sbin[r['start']], x) for r, x in zip(rs, effs)]
        passing = {n for r in mine for n in r['seq'] if n not in (base, '???') and n != main.replace('into ', '').split(' (')[0]}
        if passing: e['passing'] = sorted(passing)
        # for the map
        e['outcome'] = base if main == 'keep' else f'{base}, reversed' if main == 'reverse' else main.replace('into ', '').replace('dissolve', 'dissolved')
        e['action'], e['r'] = 'mask', MASK_R
        e['where'] = spot
        rows.append(e)
    good = lambda e: e['p'] >= 0.75
    menu = {
        'turn': sorted([e for e in rows if e['effect'] == 'keep' and good(e) and 'turn' in e], key=lambda e: e['turn']),
        'reverse': sorted([e for e in rows if e['effect'] == 'reverse' and good(e)], key=lambda e: (-e['p'], e['hold'])),
        'into': {}, 'dissolve': sorted([e for e in rows if e['effect'] == 'dissolve' and good(e)], key=lambda e: (-e['p'], e['hold'])),
        'passing': {}, 'unreliable': [e for e in rows if not good(e)],
        'split': sorted([e for e in rows if e['effect'].startswith(('split: ', 'replicate: ')) and good(e)], key=lambda e: (e['effect'], -e['p'], e['hold']))}
    for e in rows:
        if e['effect'].startswith('into ') and good(e): menu['into'].setdefault(e['effect'][5:], []).append(e)
        for n in e.get('passing', []): menu['passing'].setdefault(n, []).append(e)
    for v in menu['into'].values(): v.sort(key=lambda e: (-e['p'], not e['always_named'], e['hold']))
    # does the phase matter? where an action has no main outcome, does knowing the phase predict it
    # better than the same number of random groups would (the purity of the outcomes per group)?
    rng = np.random.default_rng(0)
    gain = []
    for e in menu['unreliable']:
        ph, ef = zip(*e['_effs'])
        def purity(keys):
            by = {}
            for k, x in zip(keys, ef): by.setdefault(k, []).append(x)
            return sum(max(v.count(y) for y in set(v)) for v in by.values()) / len(ef)
        gain.append(purity(ph) - np.mean([purity(rng.permutation(ph)) for _ in range(20)]))
    phase_gain = round(float(np.mean(gain)), 3) if gain else 0.0
    # swept by phase or state: is it worth telling the phase? if it changes the outcomes (beyond
    # chance) or the turns (a typical action's turn spread over the phases)
    if sbin is not None:
        spreads = [max(t) - min(t) for e in rows if e.get('by_phase') for t in [[x[3] for x in e['by_phase'] if x and len(x) > 3]] if len(t) >= 3]
        within = [x[4] for e in rows if e.get('by_phase') for x in e['by_phase'] if x and len(x) > 4]
        turn_spread = round(float(np.median(spreads)), 1) if spreads else 0.0
        turn_within = round(float(np.median(within)), 1) if within else 0.0
    for e in rows: e.pop('_effs')
    # (the turns: spread over the phases this much, and well beyond their spread within one: a chaotic
    # morph's turns, HGMLH8 metamorpha's, vary as much within a state)
    if sbin is not None: swing_extra = {'phase_useful': bool(phase_gain >= PHASE_GAIN_MIN or (turn_spread >= TURN_SPREAD_MIN and turn_spread >= 2 * turn_within)),
                                        'phase_turn_spread': turn_spread, 'phase_turn_within': turn_within}
    swing = ({'phase_bins': nbins, 'phase_kind': 'age', 'age_span': period, 'age_settle': meta['settle'], 'start_bins': sbin} if meta.get('age_sweep') else
             {'phase_bins': nbins, 'phase_kind': 'state', 'states': {**{k: meta['states'][k] for k in ('centroids', 'scale', 'share')}, 'names': state_names(meta['states']['centroids'])}, 'start_bins': sbin}
             if meta.get('states') else
             {'phase_bins': nbins, 'phase_kind': 'swing', 'phase_signal': meta.get('phase_signal', 'omega'), 'swing_period': period, 'start_bins': sbin}) if sbin is not None else {}
    if sbin is not None:
        if a is not None and (a.code, morph) in PHASE_USEFUL: swing_extra['phase_useful'] = True
        swing.update(swing_extra)
    return {**swing, 'morph': morph, 'base': base, 'frame': frame, 'chirality': chir, 'period': period,
            'phase_gain': phase_gain, 'phase_matters': phase_gain >= 0.1, 'mirrored': meta.get('mirrored', False),
            'omega': round(meta['omega'], 4),
            'control': [c['named_end'] for c in ctrl], 'menu': menu, 'atlas': rows,
            'trials': len(runs), 'starts': len(meta.get('starts', []))}


def name_s(e):
    if e['always_named']: return 'always named'
    return f"{e['seq']}" + (f" (??? ≤ {e['unnamed_max']})" if e['unnamed_max'] is not None else '')


def eff_s(ef, rt):
    """an effect in signs: limus + (turning clockwise), gyrans − (anticlockwise), ..."""
    c = rt['chirality']
    sg = lambda k: '+' if k > 0 else '−'
    if ef == 'keep': return f"{rt['base']} {sg(c)}" if c else rt['base']
    if ef == 'reverse': return f"{rt['base']} {sg(-c)}"
    if ef.endswith('(same turn)'): return ef[5:-12] + ' ' + sg(c)
    if ef.endswith('(reversed)'): return ef[5:-11] + ' ' + sg(-c)
    return ef.replace('into ', '').replace('dissolve', 'dissolved')


def deg_s(x):
    v = int(round(x))
    return f'{v:+d}°' if v else '0°'


def act_s(e): return f"{e['ch']} {e['spot']}, {e['hold'] / 30:g} s"


def steer_md(rt):
    m = rt['menu']
    ph = rt['phase_matters']
    PH, PS = (' phase |', '---|') if ph else ('', '')
    pc = lambda e: f" {e['phases']} |" if ph else ''
    turning = {1: 'turns clockwise (+)', -1: 'turns anticlockwise (−)', 0: 'no turning direction'}[rt['chirality']]
    L = [f"## {rt['morph']}", '',
         f"Untouched: named {', '.join(sorted(set(rt['control'])))}; turn rate {rt['omega']:+.3f} rad/s ({turning})"
         f"{'; the pattern mirrored so it turns +' if rt['mirrored'] else ''}; frame: {'heading' if rt['frame'] == 'vel' else 'body axis (still)'}. "
         f"{rt['trials']} trials ({rt['starts']} starts: phases × angles). "
         + (f"Cycle {rt['period']} steps, and the phase matters (outcomes per phase {rt['phase_gain']:+.2f} purer than chance): see the phase column."
            if ph else f"Cycle {rt['period']} steps, but the phase does not matter (outcomes per phase {rt['phase_gain']:+.2f} purer than chance), so it is left out.")
         + (f" Its body axis spins {rt['spin']:+.3f} rad per unit time." if rt.get('spin') else '')
         + (f" Oscillating (heading wobble {rt.get('wobble_deg')}°, every {rt.get('wobble_period')} steps)." if rt.get('oscillating') else '')
         + (' Drifting (a slow spin).' if rt.get('drift') else '')
         + (f" Swept by state ({', '.join(f'{n} {sh}' for n, sh in zip(rt['states']['names'], rt['states']['share']))} steps of a long run): each action's outcome per state in the JSON (by_phase)." if rt.get('phase_kind') == 'state'
            else f" Swept over its {rt['swing_period']}-step {'speed pulse' if rt.get('phase_signal') == 'speed' else 'swing'}: each action's outcome per phase ({rt['phase_bins']} bins, from the {'speed' if rt.get('phase_signal') == 'speed' else 'turn rate'}'s cycle) in the JSON (by_phase)." if rt.get('phase_kind') == 'swing' else '')
         + (f" Chaotic: {rt['chaotic']}." if rt.get('chaotic') else '') + (f" Unstable: {rt['unstable']}." if rt.get('unstable') else '')
         + (f" {rt['note'][0].upper() + rt['note'][1:]}." if rt.get('note') else '') + (' Not for maneuvers in the game.' if rt.get('no_maneuver') else ''), '']
    if m['turn']:
        L += [f"**Turn** (stays {rt['base']}{', still turning +' if rt['chirality'] else ''}), from anticlockwise to clockwise:", '',
              '| turn | action | name box | p |' + PH, '|---|---|---|---|' + PS]
        for e in m['turn']:
            L.append(f"| {deg_s(e['turn'])} ± {abs(e['turn_sd']):.0f} | {act_s(e)} | {name_s(e)} | {e['p']:.2f} |" + pc(e))
        L.append('')
    if m['reverse']:
        L += ['**Change direction** (stays ' + rt['base'] + ', now turning −: its mirror image):', '',
              '| turn | action | name box | p |' + PH, '|---|---|---|---|' + PS]
        for e in sorted(m['reverse'], key=lambda e: e.get('turn', 0)):
            L.append(f"| {deg_s(e.get('turn', 0))} | {act_s(e)} | {name_s(e)} | {e['p']:.2f} |" + pc(e))
        L.append('')
    for to, es in sorted(m['into'].items()):
        to = to.replace('(same turn)', '(turning +)').replace('(reversed)', '(turning −)')
        L += [f"**Into {to}:**", '', '| action | name box | p |' + PH, '|---|---|---|' + PS]
        for e in es[:6]:
            L.append(f"| {act_s(e)} | {name_s(e)} | {e['p']:.2f} |" + pc(e))
        if len(es) > 6: L.append(f'| *{len(es)} actions in all* | | |' + (' |' if ph else ''))
        L.append('')
    if m.get('split'):
        L += ['**Split or replicate** (several ships at the end; replicate: copies of one morph, the mass grown, e.g. through an unnamed ??? stage such as PJ1EF7 vagus):', '',
              '| outcome | action | name box | p |' + PH, '|---|---|---|---|' + PS]
        for e in m['split'][:12]:
            o = e['effect'].replace('split: ', 'split into ').replace('replicate: ', 'replicates into ') + (f" ×{e['copies'][0]}–{e['copies'][1]}" if 'copies' in e else '')
            L.append(f"| {o} | {act_s(e)} | {name_s(e)} | {e['p']:.2f} |" + pc(e))
        L.append('')
    if m['dissolve']:
        L += ['**Dissolve:** ' + '; '.join(f"{act_s(e)} (p {e['p']:.2f})" for e in m['dissolve'][:6]), '']
    if m['passing']:
        L += ['**Seen in passing** (named for a while, then gone): ' + '; '.join(
            f"{n}: {', '.join(act_s(e) for e in es[:3])}" for n, es in sorted(m['passing'].items())), '']
    if m['unreliable']:
        why = 'it varies with the phase' if ph else 'it varies from trial to trial, not with the phase'
        L += [f'**Unreliable** (no outcome in 3 of 4 trials; {why}; + turning clockwise, − anticlockwise):', '', '| action | outcomes |' + (' when the main one happens |' if ph else ''),
              '|---|---|' + PS]
        for e in sorted(m['unreliable'], key=lambda e: (e['spot'], e['hold'])):
            L.append(f"| {act_s(e)} | {', '.join(f'{eff_s(k, rt)} {v}' for k, v in e['effects'].items())} |" + (f" {eff_s(e['effect'], rt)}: {e['phases']} |" if ph else ''))
        L.append('')
    L += [f"Map: `models/{{code}}-{rt['morph']}.png`", '']
    return L


def clone_model(data):
    meta, runs = data['meta'], data['runs']
    ctrl = [r['label'] for r in runs if r['none']]
    groups = group_runs([r for r in runs if not r['none']], ('dist', 'bearing', 'turn', 'flip'))
    outs, atlas = {}, []
    for k, rs in groups.items():
        labels = [r['label'] for r in rs]
        top = mode(labels)
        e = {'dist': k[0], 'bearing': k[1], 'turn': k[2], 'flip': k[3], 'outcome': top, 'p': round(labels.count(top) / len(rs), 2),
             'n': len(rs), 'outcomes': {l: labels.count(l) for l in sorted(set(labels))}}
        atlas.append(e)
        if e['p'] >= 0.5 and top not in ('dissolved',):
            outs.setdefault(top, []).append(e)
    for es in outs.values(): es.sort(key=lambda e: (-e['p'], e['dist']))
    return {'ship': meta['ship'], 'clone': meta['clone'], 'frame': meta['frame'], 'control': ctrl,
            'outcomes': outs, 'count': {o: len(es) for o, es in outs.items()}, 'atlas': atlas}


def fmt_off(e): return f"({e['off'][0]:+.2f}, {e['off'][1]:+.2f}) {e['where']}"


OSC_WOBBLE, OSC_PERIOD, OSC_SPEED, SPIN_MIN = 15.0, 300, 0.04, 0.05   # oscillating: a real mover (R per unit time), its heading wobbling (deg) with a cycle of at most OSC_PERIOD steps (not a slow drift); spinning (rad per unit time; a still morph's axis drifts less)


# What watching them in the lab shows that the measurements don't (the author's notes): drifting
# (a slow spin / turn), oscillating (regular: a swing with a steady cycle), chaotic (phases that
# come and go at unpredictable times), unstable (turns into another in time), and no_maneuver (the
# game doesn't steer by its model)
G6_NOTE = 'amoeba-like: red dots round its blue body push or pull it, themselves appearing, vanishing and moving at random; maneuvers not well defined'
NOTES = {
    ('VT049W', 'sequerectus'): {'drift': True},
    ('KH97WU', 'limus'): {'drift': True},
    ('Z18A9R', 'antigyrans'): {'note': 'its measured spin rate is unreliable (B and G centres nearly together)'},
    # a still morph pushed out of place, in absolute directions (any side; not by its axis): a G mask on
    # one side, from the sweeps' lithos trials (the mask's board direction against where it went).
    # convert: it becomes the mover `into`, heading `dev` degrees off straight away from the mask (0:
    # away; 80: sideways, to either side); shove: it stays itself and slides `move` R away. dist: the
    # spot distance (body, edge), hold: steps, p: how often.
    ('VT049W', 'lithos'): {'pushes': [{'kind': 'convert', 'into': 'tachus', 'dist': 'edge', 'hold': 60, 'p': 0.95, 'dev': 0}]},
    ('GDNQYX', 'lithos'): {'pushes': [{'kind': 'convert', 'into': 'tardus', 'dist': 'body', 'hold': 60, 'p': 0.84, 'dev': 0}]},
    ('Z18A9R', 'lithos'): {'pushes': [{'kind': 'convert', 'into': 'tardus', 'dist': 'edge', 'hold': 90, 'p': 0.94, 'dev': 80},
                                      {'kind': 'shove', 'dist': 'body', 'hold': 90, 'p': 0.96, 'move': 0.6}]},
    ('PJ1EF7', 'lithos'): {'unstable': 'often turns into tachus after a long stay, at an unpredictable time'},
    ('KH97WU', 'impetus'): {'chaotic': 'mostly in a still phase; now and then, at unpredictable times, a moving phase, then still again'},
    ('HGMLH8', 'metamorpha'): {'chaotic': 'switches unpredictably among a tachus phase (as half a paratachus), a tardus phase (as half a paratardus) and a dashing phase', 'oscillating': False},
    ('GDNQYX', 'rectus'): {'oscillating': True}, ('KZ2BLD', 'virgatus'): {'oscillating': True}, ('HGMLH8', 'paratardus'): {'oscillating': True},
    ('HAESRE', 'zigzag'): {'oscillating': True}, ('Z18A9R', 'transgyrans'): {'oscillating': True},
    **{('G6G6CR', n): {'chaotic': G6_NOTE, 'no_maneuver': True} for n in ('rectus', 'limus', 'tardus', 'lithos', 'munitus')},
    ('G6G6CR', 'tardus'): {'chaotic': G6_NOTE + '; moves very slowly (slower than limus)', 'no_maneuver': True},
}


# a species' replications all go through this short-lived morph first (named only now and then in
# the sweeps: it mostly shows as ???): VT049W, one gyrans from a single, two from a twin
REPLICATE_VIA = {'VT049W': 'gyrans'}


# still, round morphs (no axis worth steering by): their actions are told in absolute directions,
# a push from one side; the book's push table (push_table) says what each distance and hold does
# round (near radial symmetry: a push from any side does the same, the other way) or dihedral (two
# axes of symmetry, longer one way: a push depends on its angle to the long axis, along / diagonal /
# across, the same at the 4 mirror places; either end of the long axis taken as up)
# a swept oscillator whose phase is told although its gain overall is small: a few phases do what
# no other does (HAESRE zigzag: replicates at phases 7-8 of 12, G front edge)
PHASE_USEFUL = {('HAESRE', 'zigzag')}
LITHOS_TYPE = {('5N7KKM', 'lithos'): 'round', ('VT049W', 'lithos'): 'round', ('Z18A9R', 'lithos'): 'round',
               ('HAESRE', 'lithos'): 'round', ('KH97WU', 'dot'): 'round',
               ('GDNQYX', 'lithos'): 'dihedral', ('PJ1EF7', 'lithos'): 'dihedral', ('KH97WU', 'impetus'): 'dihedral',
               ('5N7KKM', 'paralithos'): 'dihedral'}


def long_axis(m):
    """the direction a 2D field is longest along (its mass's principal axis), rad, mod pi; on the torus
    (centred on its centre of mass first). models.js / core longAxis does the same."""
    H, W = m.shape
    cx, cy = torus_com(m)
    mm = np.roll(np.roll(m, H // 2 - int(round(cy)), 0), W // 2 - int(round(cx)), 1)
    ys, xs = np.mgrid[0:H, 0:W]
    w = mm.sum()
    mx, my = (mm * xs).sum() / w, (mm * ys).sum() / w
    ixx, iyy, ixy = (mm * (xs - mx) ** 2).sum() / w, (mm * (ys - my) ** 2).sum() / w, (mm * (xs - mx) * (ys - my)).sum() / w
    return 0.5 * math.atan2(2 * ixy, ixx - iyy)


SIDE_NAMES = ('along', 'diagonal', 'across')   # (a push's angle to the long axis: 0-22.5, 22.5-67.5, 67.5-90 deg)


def torus_com(m):
    """the centre of mass (x, y) of a 2D field on the torus"""
    H, W = m.shape
    cy = np.angle((m.sum(1) * np.exp(2j * np.pi * np.arange(H) / H)).sum()) * H / (2 * np.pi)
    cx = np.angle((m.sum(0) * np.exp(2j * np.pi * np.arange(W) / W)).sum()) * W / (2 * np.pi)
    return cx % W, cy % H


def push_table(d, code, morph, sym='round'):
    """A still, round morph pushed from one side, whatever the side (its sweep's trials, the mask's
    board direction from the start boards' B -> G axis): per spot distance and hold, the main outcome
    (keep, into X, dissolve), its share p, how far it went (R, net over the observation), and its
    exit's angle off straight away from the mask (|deg|, median), with the share within 30 degrees."""
    meta, runs = d['meta'], d['runs']
    sp = os.path.join(OUT, code, f'steer-{morph}-start.npy'); zp = os.path.join(OUT, code, f'steer-{morph}-series.npz')
    if not (os.path.exists(sp) and os.path.exists(zp)): return None
    A0 = np.load(sp).astype(np.float32); H, W = A0.shape[2:]
    ax, la = [], []
    for b in A0:
        gx, gy = torus_com(b[1]); bx, by = torus_com(b[2])
        ax.append(math.atan2(wrapd(gy - by, H), wrapd(gx - bx, W)))
        la.append(long_axis(b.sum(0)))
    z = np.load(zp); x, y = np.asarray(z['x']), np.asarray(z['y']); head = np.asarray(z['head'])
    R = meta['R']
    groups = {}
    bad = {r['start'] for r in runs if r['action'] == 'none' and r['seq'] and r['seq'][0] != morph}   # (not the morph at first: left out)
    for i, r in enumerate(runs):
        if r['action'] == 'none' or abs(r['off'][0]) + abs(r['off'][1]) < 1e-6 or r['start'] in bad: continue
        end = r['named_end'] if r['named_end'] != '???' else r['label']
        eff = 'dissolve' if end in DEAD else 'keep' if end == morph else f'into {end}'
        dx, dy = x[i, -1] - x[i, 0], y[i, -1] - y[i, 0]
        dist = math.hypot(dx, dy) / R
        # (the spot's frame as run_batch held it, mid-hold: the heading if known by then (head is 0 while
        # it is not: a still start), else the start's B -> G axis, as run_batch falls back to)
        hv = float(head[i, max(0, r['dur'] // 2)])
        f0 = hv if meta['frame'] == 'vel' and hv != 0.0 else ax[r['start']]
        md = f0 + math.atan2(r['off'][1], r['off'][0])
        dev = abs(math.degrees(wrapa(math.atan2(dy, dx) - md - math.pi))) if dist > 0.3 else None
        side = None
        if sym == 'dihedral':   # (the mask's angle to the long axis, folded by the symmetry: 0 along .. 90 across)
            ph = abs(math.degrees(wrapa(2 * (md - la[r['start']]))) / 2)
            side = SIDE_NAMES[0 if ph < 22.5 else 1 if ph < 67.5 else 2]
        groups.setdefault((r['spot'].split(' ')[-1], r['dur'], side), []).append((eff, dist, dev, math.hypot(*r['off']) / R))
    out = []
    for (dist_name, hold, side), v in sorted(groups.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2] or '')):
        effs = [e for e, *_ in v]
        main = mode(effs)
        mine = [t for t in v if t[0] == main]
        devs = [t[2] for t in mine if t[2] is not None]
        out.append({**({'side': side} if side else {}), 'n': len(v), 'dist': dist_name, 'r': round(float(np.median([t[3] for t in v])), 3), 'hold': hold, 'effect': main,
                    'p': round(len(mine) / len(v), 2), 'move': round(float(np.median([t[1] for t in mine])), 1),
                    'dev': round(float(np.median(devs))) if devs else None, 'away': round(float(np.mean(np.array(devs) < 30)), 2) if devs else None})
    return out


def allowed_morphs(code):
    """the species' morphs in sandbox.txt's [add_allow_species] (the ones the game detects), or None"""
    for line in open(os.path.join(ROOT, 'sandbox.txt')):
        m = re.match(r'\[add_allow_species:[^:]*:' + code + r':([^\]]*)\]', line.strip())
        if m: return [x.strip() for x in m.group(1).split(',')]
    return None


def rename_series(runs, z, models):
    """named_series of a sweep's saved series, the stats from before the action (not saved) stood in
    for by its start's control run (the action 'none'): its first WIN steps, the same morph
    untouched (a window's score is a mean over its samples, so their order does not matter)"""
    ctrl = {r['start']: i for i, r in enumerate(runs) if r['action'] == 'none'}
    X = np.stack([z['mR'], z['mG'], z['mB'], z['speed'], np.abs(z['omega'])], -1)
    pre = X[[ctrl[r['start']] for r in runs], :WIN]
    return named_series(z, models, pre=pre)


def cmd_rename(a):
    """the sweeps' name box again with the current detection (e.g. after WIN changed), from their
    saved series: each run's seq, named_end and recover (its turns and shifts, measured once named
    again, are kept). Then book."""
    d = os.path.join(OUT, a.code)
    models = detect_models(a.code)
    for f in sorted(os.listdir(d)):
        mt = re.match(r'steer-([a-z+]+(?:-[a-z0-9]+)?)\.json$', f)
        sp = os.path.join(d, f'steer-{mt.group(1)}-series.npz') if mt else None
        if not mt or not os.path.exists(sp): continue
        data = json.load(open(os.path.join(d, f)))
        z = dict(np.load(sp, allow_pickle=True))
        named = rename_series(data['runs'], z, models)
        ch = 0
        for i, r in enumerate(data['runs']):
            old = (r['seq'], r['named_end'], r['recover'])
            r['seq'] = sequence(named[i])
            r['named_end'] = named[i, -1] or '???'
            r['recover'] = recovery(named[i], r['named_end'], r['dur'])
            ch += old != (r['seq'], r['named_end'], r['recover'])
        ctrl = {r['start']: r for r in data['runs'] if r['action'] == 'none'}
        data['meta']['control'] = [ctrl[k] for k in sorted(ctrl)]
        data['meta']['win'] = WIN
        json.dump(data, open(os.path.join(d, f), 'w'), default=jdefault)
        z['named'] = named.astype(str)
        np.savez_compressed(sp, **z)
        print(f'{a.code} {mt.group(1)}: {len(data["runs"])} runs, {ch} changed', flush=True)


def cmd_book(a):
    d = os.path.join(OUT, a.code)
    allowed = allowed_morphs(a.code)
    props = json.load(open(os.path.join(d, 'props.json'))) if os.path.exists(os.path.join(d, 'props.json')) else {}
    refs = control_refs(d)
    os.makedirs(BOOK, exist_ok=True)
    book = {'code': a.code, 'R': R_GAME, 'mask_r': MASK_R, 'conventions': {
        'off': '(fwd, side) in R from the CoM: fwd along the frame, side to its right as seen (board y down)',
        'frame': "'vel' = heading (core.js slowHeading, Upward on), 'axis' = B->G body axis",
        'turn': 'degrees, + clockwise as seen (right), measured once the name box names it again (at least 30 steps after the hold), against an untouched copy from the same start',
        'spot': 'centre, or a direction (front, front-right, ... in the frame) on the body (0.45 R) or at its edge (0.9 R)',
        'phase': 'steps since the main peak of the total mass in its cycle (0 .. period - 1), when the action starts',
        'hold': 'steps (game speed 1 = 30 steps a second)', 'p': 'share of the start phases and angles with this outcome'},
        'morphs': {}, 'clones': [],
        # the twins and triples (detected as wholes: several parts together)
        'wholes': sorted(m['name'] for m in detect_models(a.code).get('whole', [])),
        # the morph replications go through (the lab's graph puts it in the middle)
        **({'hub': REPLICATE_VIA[a.code]} if a.code in REPLICATE_VIA else {})}
    L = [f'# {a.code} morph model', '',
         'Generated by `tools/models.py book` from headless experiments (see the tool\'s docstring). '
         f'R = {R_GAME}; mask disc radius {MASK_R} cells (the game\'s). ',
         'Offsets `(fwd, side)` are in R from the centre of mass: + fwd = ahead, + side = to the right as seen. '
         'Turns in degrees, + = clockwise (right); every morph is started turning + (mirrored if its pattern turns −). Hold in steps (30 steps = 1 s at speed 1). '
         '`p` = share of the trials (start phases and angles) that gave this outcome. '
         'Actions: mask the channel at a spot for a hold; spots are the centre, or a direction in the ship\'s frame '
         '(front, front-right, ...) on the body (0.45 R from the centre) or at its edge (0.9 R). '
         '*Name box*: what the game\'s detection shows, start to end (??? = not named: not steerable by this book). '
         'Turns are measured once it is named again, against an untouched copy from the same start. '
         'A mirror-image ship (turning the other way) takes the mirrored action: left and right swapped, turn negated.', '']
    files = sorted(f for f in os.listdir(d) if f.endswith('.json'))
    for f in files:
        mt = re.match(r'steer-([a-z+]+)\.json$', f)
        if not mt: continue
        if allowed and mt.group(1) not in allowed: print(f'  {mt.group(1)}: not in [add_allow_species] (kept; the lab dims it)')
        zp = os.path.join(d, f'steer-{mt.group(1)}-series.npz')
        data = json.load(open(os.path.join(d, f)))
        main = data
        zz = load_ends(zp) if os.path.exists(zp) else None
        # sweeps of the same starts with other mask channels (steer --like): their trials merged in
        for g in sorted(os.listdir(d)):
            mg = re.match(rf'steer-{re.escape(mt.group(1))}-([a-z0-9]+)\.json$', g)
            if not mg: continue
            d2 = json.load(open(os.path.join(d, g)))
            if not d2['meta'].get('like'): continue
            keep = [i for i, r in enumerate(d2['runs']) if r['action'] != 'none']
            z2p = os.path.join(d, f"steer-{mt.group(1)}-{mg.group(1)}-series.npz")
            data = {**data, 'runs': data['runs'] + [d2['runs'][i] for i in keep]}
            if zz is not None and os.path.exists(z2p):
                z2 = load_ends(z2p)
                zz = {k: np.concatenate([zz[k], z2[k][keep]]) for k in zz}
            print(f"  {mt.group(1)}: + {d2['meta']['channels']} masks ({len(keep)} trials, {g})", flush=True)
        rt = steer_model(data, a, refs, zz)
        if (a.code, rt['morph']) in LITHOS_TYPE:   # (absolute directions: what a push from a side does)
            rt['lithos_type'] = LITHOS_TYPE[(a.code, rt['morph'])]   # ('round' / 'dihedral')
            rt['push_table'] = push_table(main, a.code, rt['morph'], rt['lithos_type'])   # (the main sweep's: G masks)
        rbp = os.path.join(d, f"robust-{rt['morph']}.json")   # (robust: each reliable action's share a cell or a few steps off)
        if os.path.exists(rbp):
            rob = {(x['ch'], x['spot'], x['hold']): x for x in json.load(open(rbp))['actions']}
            for e in rt['atlas']:
                x = rob.get((e['ch'], e['spot'], e['hold']))
                if x:
                    e['p_robust'] = x['p_robust']
                    if 'p_fresh' in x: e['p_fresh'] = x['p_fresh']
            print(f"  {rt['morph']}: + robustness of {len(rob)} actions", flush=True)
        pp = props.get(rt['morph'])
        if pp:   # its motion over a long run (props): a still spinner's spin, an oscillator's wobble
            rt['speed'] = pp['speed']   # (R per unit time, over a long run)
            rt['spin'] = pp['spin'] if abs(pp['spin']) >= SPIN_MIN else 0.0   # (of the start the sweep used: mirrored if it was)
            travels = pp['speed'] >= OSC_SPEED   # (not a slow drift, whose heading is mostly noise, e.g. PJ1EF7 lithos; a zigzag nets little)
            rt['oscillating'] = travels and pp['wobble_deg'] >= OSC_WOBBLE and 0 < pp['wobble_period'] <= OSC_PERIOD
            rt['wobble_deg'], rt['wobble_period'], rt['mass_pulse'] = pp['wobble_deg'], pp['wobble_period'], pp['mass_pulse']
        rt['whole'] = rt['morph'] in book['wholes']
        rt.update(json.loads(json.dumps(NOTES.get((a.code, rt['morph']), {}))))
        for pu in rt.get('pushes', []):   # (its distance in R, as the spot's)
            sp = next((e for e in rt['atlas'] if e['spot'] == f"right {pu['dist']}"), None)
            if sp: pu['r'] = round(math.hypot(*sp['off']), 3)
        if a.code in REPLICATE_VIA:
            for e in rt['atlas']:
                if e['effect'].startswith('replicate: '): e['via'] = REPLICATE_VIA[a.code]
        book['morphs'][rt['morph']] = rt
        L += [x.replace('{code}', a.code) for x in steer_md(rt)]
        if not a.no_maps: steer_map(json.load(open(os.path.join(d, f))), rt, a.code)
    for f in files:
        if not f.startswith('clone-'): continue
        rt = clone_model(json.load(open(os.path.join(d, f))))
        book['clones'].append(rt)
        L += [f"## clone: {rt['ship']} + {rt['clone']}", '',
              'A copy of the clone\'s pattern pasted at `dist` R from the ship\'s centre, at `bearing` (0 = ahead, 90 = right), '
              'its body axis turned `turn`° from the ship\'s frame (+ clockwise), mirrored if `flip`.', '',
              '| outcome | positions | best: dist, bearing, turn, flip | p |', '|---|---|---|---|']
        for o, es in sorted(rt['outcomes'].items(), key=lambda kv: -len(kv[1])):
            bests = '; '.join(f"{e['dist']} R, {e['bearing']}°, {e['turn']}°{', flip' if e['flip'] else ''}" for e in es[:3])
            L.append(f"| {o} | {len(es)} | {bests} | {es[0]['p']:.2f} |")
        L.append('')
    # composite maneuvers tried (composite runs): sequences of maneuvers, their results per value of
    # {w} (files composite-<name>.json; strategy-<name>.json from before the name changed)
    book['composites'] = []
    for f in sorted(os.listdir(d)):
        if f.startswith(('composite-', 'strategy-')) and f.endswith('.json'):
            st = json.load(open(os.path.join(d, f)))
            book['composites'].append({k: st[k] for k in ('name', 'morph', 'steps', 'vals', 'note')} | {'results': [
                {k: r[k] for k in ('val', 'p', 'time', 'tries', 'turn', 'turn_sd', 'travel', 'fails')} for r in st['results']]})
            L += [f"## composite: {st['name']} ({st['morph']})", '', f"{st['note']}. Steps: `{st['steps']}`.", '',
                  '| {w} | p | steps | tries | net turn | travel vs untouched | fails |', '|---|---|---|---|---|---|---|']
            for r in st['results']:
                L.append(f"| {r['val'] or ''} | {r['p']:.2f} | {r['time'] or ''} | {r['tries'] or ''} | {'' if r['turn'] is None else f"{r['turn']:+.0f}° ± {r['turn_sd']:.0f}"} | {'' if r['travel'] is None else f"{r['travel']:+.1f} R"} | {', '.join(f'{k} {v}' for k, v in r['fails'].items())} |")
            L.append('')
    slim = json.loads(json.dumps(book))
    if not a.atlas:
        # per morph, every action tried (the lab's Model mode draws them): its spot, hold, main
        # outcome and share, the turn, and the name box
        keep = ('ch', 'spot', 'hold', 'off', 'effect', 'p', 'turn', 'always_named', 'unnamed_max', 'seq', 'copies', 'by_phase', 'via', 'alt', 'alt_by_phase', 'p_robust', 'p_fresh')
        for m in slim['morphs'].values(): m['actions'] = [{k: e[k] for k in keep if k in e} for e in m.pop('atlas')]   # (the action level)
        for c in slim['clones']: c.pop('atlas')
    json.dump(slim, open(os.path.join(BOOK, f'{a.code}.json'), 'w'), indent=1)
    open(os.path.join(BOOK, f'{a.code}.md'), 'w').write('\n'.join(L) + '\n')
    print(f'wrote models/{a.code}.md and models/{a.code}.json')


def steer_map(data, rt, code):
    """One picture per morph, a panel per (channel, hold); the creature in its frame behind
    (heading up: yellow arrow; its turning direction, + clockwise: the blue curved arrow).
    Each spot tried: an arrow if it stays the morph, pointing to its new heading; a square if it
    becomes something else, lettered with the morph, grey if it dissolves. Blue: its turning
    direction kept (or none to keep), orange: it turns the other way after. Faded: the outcome
    holds in under 3 of 4 trials."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    meta = data['meta']
    R = meta['R']
    start = np.load(os.path.join(OUT, code, f"steer-{meta['morph']}-start.npy")).astype(float)[0]
    ser = np.load(os.path.join(OUT, code, f"steer-{meta['morph']}-series.npz"))
    i0 = next(i for i, r in enumerate(data['runs']) if r['action'] == 'none' and r['start'] == 0)
    H, W = start.shape[1:]
    yy, xx = np.mgrid[0:H, 0:W]
    cm = lambda m: (np.angle((m * np.exp(2j * np.pi * yy / H)).sum()) * H / (2 * np.pi), np.angle((m * np.exp(2j * np.pi * xx / W)).sum()) * W / (2 * np.pi))
    if meta['frame'] == 'vel': ang = float(ser['head'][i0, 0])
    else:
        g, b = cm(start[1]), cm(start[2])
        ang = math.atan2(wrapd(g[0] - b[0], H), wrapd(g[1] - b[1], W))
    cy, cx = cm(start.sum(0))
    atlas = rt['atlas']
    reach = max(math.hypot(*e['off']) for e in atlas) * R
    ext = reach + MASK_R + 8
    g = np.linspace(-ext, ext, 240)
    S, F = np.meshgrid(g, g[::-1])   # cols: side right; rows: fwd up
    bx = cx + F * math.cos(ang) - S * math.sin(ang)
    by = cy + F * math.sin(ang) + S * math.cos(ang)
    img = np.stack([start[k][np.round(by).astype(int) % H, np.round(bx).astype(int) % W] for k in range(3)], -1)
    img = np.clip(img, 0, 1) * 0.45
    SAME, FLIP, GONE, SPLIT, TWIN = '#2a8cff', '#ff8c00', '#808080', '#ff3df2', '#2ecc40'
    wholes = {m['name'] for m in detect_models(code).get('whole', [])}
    code_morph = rt['morph']
    base, chir = rt['base'], rt['chirality'] or 1     # starts turn + (mirrored if not)

    def square(e):
        """(colour, letter) of an outcome that is not the morph"""
        ef = e['effect']
        if ef == 'dissolve': return GONE, ''
        if ef.startswith('split: '):     # splits into two (or more) ships: magenta, their initials (LL, GL)
            return SPLIT, ''.join(n[0].upper() if n[0].isalpha() else '?' for n in ef[7:].split(' + '))
        if ef.startswith('replicate: '):   # replicates into several copies: magenta, its initial and + (T+)
            return SPLIT, ef[11].upper() + '+'
        name = ef.replace('into ', '').split(' (')[0]
        if name in wholes and code_morph not in wholes: return TWIN, name   # a single into a twin / triple: green
        return (FLIP if ef.endswith('(reversed)') else SAME), name
    holds = sorted({e['hold'] for e in atlas})
    chans = sorted({e['ch'] for e in atlas})
    P = 5.2
    fig, axs = plt.subplots(len(chans), len(holds), figsize=(P * len(holds), P * len(chans) + 2.2), squeeze=False, facecolor='#111')
    E = ext / R
    names = {}
    for i, ch in enumerate(chans):
        for j, du in enumerate(holds):
            ax = axs[i, j]
            ax.set_facecolor('#000')
            ax.imshow(img, extent=(-E, E, -E, E))
            for e in atlas:
                if e['ch'] != ch or e['hold'] != du: continue
                f, s = e['off']
                alpha = 1.0 if e['p'] >= 0.75 else 0.35
                if e['effect'] in ('keep', 'reverse'):
                    a = math.radians(e.get('turn', 0.0))    # the new heading, from up, + clockwise as seen
                    col = SAME if e['effect'] == 'keep' else FLIP
                    L = 0.13
                    dx, dy = L * math.sin(a), L * math.cos(a)
                    ax.add_patch(FancyArrowPatch((s - dx, f - dy), (s + dx, f + dy), arrowstyle='simple,head_width=1.6,head_length=1.2,tail_width=0.7',
                                                 mutation_scale=11, fc=col, ec='#000', lw=0.8, alpha=alpha))
                else:
                    col, name = square(e)
                    ax.plot(s, f, 's', ms=16, color=col, alpha=alpha, mec='#000', mew=0.8)
                    if name:
                        mark = name if col == SPLIT else name[0].upper()
                        if col != SPLIT: names.setdefault(mark, set()).add(name)
                        ax.text(s, f, mark, ha='center', va='center', fontsize=10 if len(mark) == 1 else 7.5, fontweight='bold', color='#000', alpha=alpha)
            ax.set_xlim(-E, E); ax.set_ylim(-E, E)
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f'{ch} mask, hold {du} steps ({du / 30:g} s)', color='#ddd', fontsize=12)
            # heading (yellow) and turning direction (+ clockwise: blue)
            ax.add_patch(FancyArrowPatch((0, E * 0.72), (0, E * 0.95), arrowstyle='-|>', mutation_scale=18, color='#ffe600', lw=2))
            if rt['chirality']:
                c = rt['chirality']
                ax.add_patch(FancyArrowPatch((-0.16 * c, E * 0.80), (0.16 * c, E * 0.80), connectionstyle=f'arc3,rad={-0.6 * c}',
                                             arrowstyle='-|>', mutation_scale=14, color='#ffe600', lw=2))
    # legend
    arrow = lambda c: plt.Line2D([], [], marker='^', ls='', ms=12, mec='#000', color=c)
    sq = lambda c: plt.Line2D([], [], marker='s', ls='', ms=12, mec='#000', color=c)
    hs = [arrow(SAME), arrow(FLIP), sq(SAME), sq(FLIP), sq(TWIN), sq(SPLIT), sq(GONE)]
    lab = [f'stays {base}, same turning direction; arrow: its new heading', f'stays {base}, turns the other way (its mirror image)',
           'becomes the lettered morph, same turning direction', 'becomes the lettered morph, turning the other way', 'becomes a twin / triple (from a single)', 'splits into two ships (LL: limus + limus), or replicates (T+: copies of tachus)', 'dissolved']
    key = '   '.join(f"{k}: {', '.join(sorted(v))}" for k, v in sorted(names.items()))
    fig.legend(hs, lab, loc='lower center', ncol=3, facecolor='#222', labelcolor='#ddd', fontsize=11, frameon=False)
    fig.suptitle(f"{code} {base}{' (mirrored pattern)' if rt.get('mirrored') else ''}: where the mask is held → what it does\n"
                 f"blue = turning direction kept (or none), orange = turned the other way · faded: under 3 of 4 trials · yellow: heading and turning"
                 + (f" · letters: {key}" if key else ''), color='#eee', fontsize=13)
    fig.tight_layout(rect=(0, 0.09, 1, 0.93))
    fig.savefig(os.path.join(BOOK, f"{code}-{rt['morph']}.png"), dpi=110, facecolor='#111')
    plt.close(fig)


# ---- show one action as a GIF ------------------------------------------------------------------
def outcome_hits(effect, morph, label):
    """a board's name is an action's outcome: keep / reverse (the morph), into X, dissolve, split: X
    (its parts), replicate (several bodies, or multiplied / exploded); turning not judged"""
    if effect in ('keep', 'reverse'): return label == morph
    if effect == 'dissolve': return label == 'dissolved'
    if effect.startswith('into '): return label == effect.split(' (')[0][5:]
    if effect.startswith('split: '): return label == effect[7:]
    if effect.startswith('replicate'): return ' + ' in label or label in ('multiplied', 'exploded')
    return False


def cmd_robust(a):
    """Robustness of a morph's reliable actions (p >= --pmin, not keep): each tried again with its mask
    a cell off (ahead, behind, either side) and its hold --dhold steps shorter or longer, from the
    sweep's starts (--starts of them): p_robust, the share of those near misses with its outcome (the
    lab and the game place a mask a cell or so off and hold it a step or two off; a near-kill
    knife-edge, e.g. KH97WU gyrans B centre 15 -> limus, shows as low). Saved as
    output/models/<code>/robust-<morph>.json; book adds p_robust to the actions."""
    book = json.load(open(os.path.join(BOOK, f'{a.code}.json')))
    data = json.load(open(os.path.join(OUT, a.code, f'steer-{a.morph}.json')))
    meta, m = data['meta'], book['morphs'][a.morph]
    R = meta['R']
    acts = [e for e in m['actions'] if e['p'] >= a.pmin and e['effect'] != 'keep']
    seen, todo = set(), []
    for e in sorted(acts, key=lambda e: -e['p']):
        k = (e['ch'], e['spot'], e['hold'])
        if k not in seen: seen.add(k); todo.append(e)
    print(f"{a.code} {a.morph}: {len(todo)} actions at p >= {a.pmin}", flush=True)
    sim = Sim(rule_of(a.code), meta.get('size', 128), meta.get('size', 128), R)
    models = detect_models(a.code)
    st = meta['starts'][::max(1, len(meta['starts']) // a.starts)][:a.starts]
    S = start_set(sim, a.code, a.morph, st, meta['settle'], R)
    ns = len(st)
    jit = [(dx, dy, dh) for dx, dy in ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)) for dh in (-a.dhold, 0, a.dhold) if (dx, dy, dh) != (0, 0, 0)]
    res = []
    for e in todo:
        boards = [(j, s_) for j in jit for s_ in range(ns)]
        idx = [s_ for _, s_ in boards]
        plan = {'action': ['mask'] * len(idx), 'ch': [CH.index(e['ch'])] * len(idx), 'r': [float(MASK_R)] * len(idx),
                'dur': [max(1, e['hold'] + j[2]) for j, _ in boards],
                'off': [(e['off'][0] * R + j[0], e['off'][1] * R + j[1]) for j, _ in boards],
                'frame': [m['frame']] * len(idx), 'axis0': np.array([S['axis'][s_] for s_ in idx]) if 'axis' in S else np.zeros(len(idx))}
        steps = e['hold'] + a.dhold + meta['observe']
        _, out, _ = run_batch(sim, S['A'][idx].clone(), stack_kin(S['kin'], idx, R, sim.dt), S['com'][idx].copy(), plan, steps, R, models, a.debris)
        ok = [outcome_hits(e['effect'], a.morph, o['label']) for o in out]
        by = {f"{'%+d' % j[0] if j[0] else '0'},{'%+d' % j[1] if j[1] else '0'} hold {e['hold'] + j[2]}": round(float(np.mean(ok[k * ns:(k + 1) * ns])), 2) for k, j in enumerate(jit)}
        r = {'ch': e['ch'], 'spot': e['spot'], 'hold': e['hold'], 'effect': e['effect'], 'p': e['p'], 'p_robust': round(float(np.mean(ok)), 2),
             'by': by, 'fails': dict(Counter(o['label'] for o, k in zip(out, ok) if not k).most_common(4))}
        res.append(r)
        print(f"  {e['ch']} {e['spot']} {e['hold']}: {e['effect']} p {e['p']:.2f} -> robust {r['p_robust']:.2f}  fails {r['fails']}", flush=True)
    json.dump({'code': a.code, 'morph': a.morph, 'pmin': a.pmin, 'dhold': a.dhold, 'starts': ns, 'actions': res},
              open(os.path.join(OUT, a.code, f'robust-{a.morph}.json'), 'w'), default=jdefault)


def cmd_verify(a):
    """Re-test the actions the lab offers (effect not keep / dissolve, p >= 0.75; a replication from
    0.5) of each of a species' morphs, from fresh starts (the sweep's start ages + --shift steps:
    states the sweep never saw) and as the lab places them: exactly, the heading off by +-dang
    degrees, the mask a cell off (ahead, behind, either side), the hold --dhold steps off. An
    outcome counts as the board's name at the end (a replication: two or more bodies then, not
    only a swelling). Saved as output/models/<code>/robust-<morph>.json (p_robust: the share over
    all of these; p_fresh: placed exactly); book adds them, models.js takes the lower p."""
    book = json.load(open(os.path.join(BOOK, f'{a.code}.json')))
    models = detect_models(a.code)
    morphs = a.morph.split(',') if a.morph else sorted(book['morphs'])
    for morph in morphs:
        fp = os.path.join(OUT, a.code, f'steer-{morph}.json')
        if not os.path.exists(fp): continue
        meta = json.load(open(fp))['meta']
        m, R = book['morphs'][morph], meta['R']
        seen, todo = set(), []
        for e in sorted(m['actions'], key=lambda e: -e['p']):
            if e['effect'] in ('keep', 'dissolve'): continue
            if e['p'] < (0.5 if e['effect'].startswith('replicate') else a.pmin): continue
            k = (e['ch'], e['spot'], e['hold'])
            if k not in seen: seen.add(k); todo.append(e)
        if not todo: continue
        t0 = time.time()
        sim = Sim(rule_of(a.code), meta.get('size', 128), meta.get('size', 128), R)
        st = meta['starts'][::max(1, len(meta['starts']) // a.starts)][:a.starts]
        fixed_age = meta.get('age_sweep') or meta.get('made')
        st = [(q, e if fixed_age else e + a.shift) for q, e in st]
        S = (made_start_set(sim, a.code, morph, meta['made'], st, meta['settle'], R, models) if meta.get('made')
             else start_set(sim, a.code, morph, st, meta['settle'], R, mirror=meta.get('mirrored')))
        ns = S['A'].shape[0]
        dr = math.radians(a.dang)
        # (rotation of the offset, cells ahead, cells right, hold change)
        jit = [(0, 0, 0, 0), (dr, 0, 0, 0), (-dr, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0),
               (0, 0, 0, a.dhold), (0, 0, 0, -a.dhold)]
        res = []
        for hold in sorted({e['hold'] for e in todo}):
            es = [e for e in todo if e['hold'] == hold]
            boards = [(ei, j, s_) for ei in range(len(es)) for j in range(len(jit)) for s_ in range(ns)]
            ok = []
            for b0 in range(0, len(boards), a.batch):
                ch_ = boards[b0:b0 + a.batch]
                idx = [s_ for _, _, s_ in ch_]
                offs = []
                for ei, j, _ in ch_:
                    f, sd = es[ei]['off'][0] * R, es[ei]['off'][1] * R
                    th, df, ds, _ = jit[j]
                    offs.append((f * math.cos(th) - sd * math.sin(th) + df, f * math.sin(th) + sd * math.cos(th) + ds))
                plan = {'action': ['mask'] * len(idx), 'ch': [CH.index(es[ei]['ch']) for ei, _, _ in ch_], 'r': [float(MASK_R)] * len(idx),
                        'dur': [max(1, hold + jit[j][3]) for _, j, _ in ch_], 'off': offs,
                        'frame': [m['frame']] * len(idx), 'axis0': np.array([S['axis'][s_] for s_ in idx])}
                _, out, _ = run_batch(sim, S['A'][idx].clone(), stack_kin(S['kin'], idx, R, sim.dt), S['com'][idx].copy(), plan,
                                      hold + a.dhold + meta['observe'], R, models, a.debris)
                ok += [(outcome_hits(es[ei]['effect'], morph, o['label']), o['label']) for (ei, _, _), o in zip(ch_, out)]
            for ei, e in enumerate(es):
                mine = ok[ei * len(jit) * ns:(ei + 1) * len(jit) * ns]
                hits = [h for h, _ in mine]
                exact = hits[:ns]
                r = {'ch': e['ch'], 'spot': e['spot'], 'hold': e['hold'], 'effect': e['effect'], 'p': e['p'],
                     'p_fresh': round(float(np.mean(exact)), 2), 'p_robust': round(float(np.mean(hits)), 2),
                     'by': {f"{'rot %+d' % round(math.degrees(j[0])) if j[0] else ''}{' fwd %+d' % j[1] if j[1] else ''}{' side %+d' % j[2] if j[2] else ''}{' hold %+d' % j[3] if j[3] else ''}".strip() or 'exact':
                            round(float(np.mean(hits[k * ns:(k + 1) * ns])), 2) for k, j in enumerate(jit)},
                     'fails': dict(Counter(lb for h, lb in mine if not h).most_common(4))}
                res.append(r)
                flag = '  <<' if r['p_robust'] < 0.5 else ''
                print(f"  {morph} {e['ch']} {e['spot']} {e['hold']}: {e['effect']} p {e['p']:.2f} -> fresh {r['p_fresh']:.2f} robust {r['p_robust']:.2f}  fails {r['fails']}{flag}", flush=True)
        json.dump({'code': a.code, 'morph': morph, 'pmin': a.pmin, 'dhold': a.dhold, 'dang': a.dang, 'shift': a.shift, 'starts': ns,
                   'kind': 'verify', 'actions': res}, open(os.path.join(OUT, a.code, f'robust-{morph}.json'), 'w'), default=jdefault)
        print(f"== {a.code} {morph}: {len(res)} actions, {ns} starts, {time.time() - t0:.0f} s", flush=True)


def cmd_show(a):
    """Replay one action on a fresh settled copy (start --start of the sweep's starts) and save a
    GIF following the ship, with the mask disc and the name box as the game would show it."""
    from PIL import Image, ImageDraw
    import torch
    data = json.load(open(os.path.join(OUT, a.code, f'steer-{a.morph}.json')))
    meta = data['meta']
    R = meta['R']
    spots = dict(named_spots(R, tuple(DISTS), meta.get('spot_scale', 1.0)))
    off = spots[a.spot]
    rule = rule_of(a.code)
    size = a.size or meta.get('size', 128)   # (the sweep's board, unless --size)
    sim = Sim(rule, size, size, R)
    models = detect_models(a.code)
    q, extra = meta['starts'][a.start]
    S = start_set(sim, a.code, a.morph, [(q, extra)], meta['settle'], R)
    A, kin = S['A'], S['kin']
    meas = sim.measure(A)
    ch = CH.index(a.ch)
    oh = torch.zeros(1, 3, 1, 1, device=sim.dev); oh[0, ch] = 1
    total = a.hold + a.after
    imgs, feats = [], []
    for s in range(total):
        fa = S['axis'][0] if meta['frame'] == 'axis' or np.isnan(kin.slow[0]) else kin.slow[0]
        cx = meas[0, 0] + off[0] * math.cos(fa) - off[1] * math.sin(fa)
        cy = meas[0, 1] + off[0] * math.sin(fa) + off[1] * math.cos(fa)
        if s < a.hold:
            disc = sim.disc(torch.tensor([cx], device=sim.dev), torch.tensor([cy], device=sim.dev), float(MASK_R))[:, None]
            A = sim.step(A, 1 - disc * oh)
        else:
            A = sim.step(A)
        m2 = sim.measure(A)
        d = np.stack([wrapd(m2[:, 0] - meas[:, 0], sim.W), wrapd(m2[:, 1] - meas[:, 1], sim.H)], 1)
        meas = m2
        kin.update(d / sim.dt)
        feats.append([*(m2[0, 2:5] / R ** 2), kin.speed[0] / R, kin.omega[0]])
        if (s + 1) % 100 == 0:   # (the bodies on the board now: an explosion into copies shows here)
            bs = bodies_of(A[0].cpu().numpy(), R, a.debris)
            print(f'  step {s + 1}: {len(bs)} bodies, masses {[round(b[2], 1) for b in bs]}', flush=True)
        if s % a.every == 0:
            img = (np.clip(A[0].cpu().numpy(), 0, 1).transpose(1, 2, 0) * 255).astype(np.uint8)
            sh = (sim.H // 2 - int(meas[0, 1]), sim.W // 2 - int(meas[0, 0]))
            imgs.append((s, np.roll(img, sh, (0, 1)), (wrapd(cx - meas[0, 0], sim.W) + sim.W // 2, wrapd(cy - meas[0, 1], sim.H) + sim.H // 2)))
    F = np.array(feats)
    ser = {'mR': F[None, :, 0], 'mG': F[None, :, 1], 'mB': F[None, :, 2], 'speed': F[None, :, 3], 'omega': F[None, :, 4]}
    named = named_series(ser, models, pre=S['pre'])[0]
    z = 3
    frames = []
    for s, img, (ox, oy) in imgs:
        im = Image.fromarray(img).resize((sim.W * z, sim.H * z), Image.NEAREST)
        dr = ImageDraw.Draw(im)
        if s < a.hold:
            rr = MASK_R * z
            dr.ellipse((ox * z - rr, oy * z - rr, ox * z + rr, oy * z + rr), outline=(255, 255, 255))
        nm = named[s] or '???'
        dr.text((4, 4), f"{a.morph}: {a.ch} mask {a.spot}, {a.hold / 30:g} s   step {s}", fill=(255, 255, 0))
        dr.text((4, 18), f"name box: {nm}", fill=(255, 255, 255) if nm != '???' else (255, 80, 80))
        frames.append(im)
    tag = f"{a.ch}-{a.spot.replace(' ', '_')}-{a.hold}-s{a.start}"
    path = os.path.join(OUT, a.code, f'show-{a.morph}-{tag}.gif')
    frames[0].save(path, save_all=True, append_images=frames[1:], duration=1000 // 30 * a.every, loop=0)
    print(path, ' > '.join(sequence(named)))


def kin_run(sim, code, morph, settle, steps, R):
    """a settled copy run on for steps, as the game's kinematics see it: [speed / R, omega] a step"""
    S = start_set(sim, code, morph, [(0, 0)], settle, R)
    A, kin = S['A'], S['kin']
    meas = sim.measure(A)
    out = []
    for _ in range(steps):
        A = sim.step(A)
        m2 = sim.measure(A)
        d = np.stack([wrapd(m2[:, 0] - meas[:, 0], sim.W), wrapd(m2[:, 1] - meas[:, 1], sim.H)], 1)
        meas = m2
        kin.update(d / sim.dt)
        out.append([kin.speed[0] / R, kin.omega[0]])
    return np.array(out)


STATE_SMOOTH = 15   # steps: the kinematics smoothed so before telling the state (detection noise)


def kin_states(sim, code, morph, settle, steps, R, K, n):
    """a chaotic morph's states: k-means (K) of its long run's [speed / R, |omega|] (smoothed, each
    scaled by its spread), named by speed and turning; and n starts drawn from them, as many from
    each state as it allows, at times where it stays in that state a while (spread over the run)"""
    F = kin_run(sim, code, morph, settle, steps, R)
    ker = np.ones(STATE_SMOOTH) / STATE_SMOOTH
    X = np.stack([np.convolve(F[:, 0], ker, 'same'), np.convolve(np.abs(F[:, 1]), ker, 'same')], 1)
    scale = X.std(0) + 1e-9
    Z = X / scale
    rng = np.random.default_rng(0)
    C = Z[rng.choice(len(Z), K, replace=False)]
    for _ in range(50):
        lab = np.argmin(((Z[:, None] - C[None]) ** 2).sum(2), 1)
        C = np.array([Z[lab == k].mean(0) if (lab == k).any() else C[k] for k in range(K)])
    order = np.argsort(C[:, 0])   # slowest first
    C = C[order]; lab = np.argsort(order)[lab]
    cen = C * scale
    names = state_names(cen)
    per = max(1, n // K)
    starts = []
    for k in range(K):
        ok = [t for t in range(10, len(lab) - 10) if (lab[t - 10:t + 11] == k).all()]
        if not ok: continue
        pick = [ok[int(i)] for i in np.linspace(0, len(ok) - 1, min(per, len(ok)))]
        starts += [(0, t + 1) for t in pick]   # (the run's step t is t + 1 steps after settling)
    st = {'centroids': cen.tolist(), 'scale': scale.tolist(), 'names': names, 'share': np.bincount(lab, minlength=K).tolist(),
          'smooth': STATE_SMOOTH}
    print(f"  states: {[(nm, [round(x, 3) for x in c], sh) for nm, c, sh in zip(names, cen, st['share'])]}", flush=True)
    return st, starts


def state_names(cen):
    """names for states ([speed / R, |omega|] centroids, slowest first): still; then slow, fast by rank; turning if it turns"""
    out = []
    moving = [i for i, (sp, _) in enumerate(cen) if sp >= SPEED_MIN_R * 1.5]
    for i, (sp, om) in enumerate(cen):
        if i not in moving: out.append('still'); continue
        r = moving.index(i)
        out.append(('turning' if om > 0.06 else 'slow' if r == 0 and len(moving) > 1 else 'fast' if r == len(moving) - 1 else 'moving'))
    return out


def state_of(st, speed, omega):
    """the nearest state to [speed / R, |omega|] (models.js stateBin does the same)"""
    z = np.array([speed, abs(omega)]) / np.array(st['scale'])
    c = np.array(st['centroids']) / np.array(st['scale'])
    return int(np.argmin(((c - z) ** 2).sum(1)))


def acf_period(x, lo=8):
    """the first strong autocorrelation peak of x (steps), and its height; (0, 0) if none"""
    x = x - x.mean()
    if x.std() < 1e-9: return 0, 0.0
    ac = np.array([np.dot(x[:-l], x[l:]) for l in range(1, len(x) // 3)]) / np.dot(x, x)
    for l in range(lo, len(ac) - 1):
        if ac[l] > ac[l - 1] and ac[l] >= ac[l + 1] and ac[l] > 0.2: return l + 1, float(ac[l])
    return 0, 0.0


def cmd_cycles(a):
    """Each morph's cycles as the game's kinematics see them over a long run: the period and
    strength of its turn rate's and its speed's autocorrelation, and their spread."""
    rule = rule_of(a.code)
    for morph in a.morph.split(','):
        ext = morph_extent(a.code, morph, R_GAME)
        size = 128 if ext <= 40 else 160 if ext <= 52 else 192
        F = kin_run(Sim(rule, size, size, R_GAME), a.code, morph, a.settle, a.steps, R_GAME)
        po, ho = acf_period(F[:, 1]); ps, hs = acf_period(F[:, 0])
        print(f"{a.code} {morph}: turn rate {F[:, 1].mean():+.3f} +- {F[:, 1].std():.3f} (period {po}, acf {ho:.2f}); "
              f"speed {F[:, 0].mean():.3f} +- {F[:, 0].std():.3f} (period {ps}, acf {hs:.2f})", flush=True)
        np.save(os.path.join(OUT, a.code, f'kin-{morph}.npy'), F)


def best_into(book, src, dst):
    """the book's best action turning src into dst (any turning): (spot, hold, off fwd/side in R, ch, p), or None"""
    es = [e for e in book['morphs'][src]['actions'] if e['effect'].split(' (')[0] == f'into {dst}' and e['p'] >= 0.5]
    if not es: return None
    e = sorted(es, key=lambda e: (-round(e['p'] * 20), e['hold']))[0]
    return e


def cmd_combo(a):
    """A composite maneuver: the ship's morph A into B (its best action), a wait, then B back into A
    (B's best action, as B's frame), over start phases x angles and several waits; the net turn
    against an untouched copy, what it ends as, and how far it went. E.g. GDNQYX rectus > gyrans >
    rectus (gyrans spins meanwhile: the wait sets the angle)."""
    book = json.load(open(os.path.join(BOOK, f'{a.code}.json')))
    A_, B_ = a.morph, a.via
    if a.first_effect:   # (the first action by its outcome, e.g. VT049W rectus 'replicate: tachus', which goes through gyrans)
        es = [e for e in book['morphs'][A_]['actions'] if e['effect'].startswith(a.first_effect)]
        e1 = sorted(es, key=lambda e: (-round(e['p'] * 20), e['hold']))[0] if es else None
    else: e1 = best_into(book, A_, B_)
    # the second: B's best action back into A, or (--second all) every spot of B's model at each of --holds2
    if a.second == 'all':
        seen_spots = {}
        for e in book['morphs'][B_]['actions']: seen_spots.setdefault(e['spot'], e['off'])
        seconds = [{'spot': sp, 'off': off, 'hold': h, 'ch': c, 'p': None} for c in a.chs2 for sp, off in seen_spots.items() for h in map(int, a.holds2.split(','))]
    else: seconds = [best_into(book, B_, A_)]
    if not e1 or not seconds[0]: sys.exit(f'no action {A_} -> {B_} or {B_} -> {A_} in the model')
    e2 = seconds[0]
    R = R_GAME
    print(f'{a.code} {A_} > {B_} > {A_}: {e1["ch"]} {e1["spot"]} {e1["hold"]} (p {e1["p"]}), then {e2["ch"]} {e2["spot"]} {e2["hold"]} (p {e2["p"]})', flush=True)
    ext = morph_extent(a.code, A_, R)
    size = a.size or (128 if ext <= 40 else 160 if ext <= 52 else 192)
    sim = Sim(rule_of(a.code), size, size, R)
    models = detect_models(a.code)
    probe = start_set(sim, a.code, A_, [(0, 0)], a.settle + 240, R)
    period = probe['period']
    starts = [(q, q * 13 * period + round(k * period / a.phases)) for q in range(a.angles) for k in range(a.phases)]
    S = start_set(sim, a.code, A_, starts, a.settle, R)
    waits = [int(x) for x in a.waits.split(',')]
    steps = e1['hold'] + max(waits) + max(e['hold'] for e in seconds) + a.observe
    rows = [(st, w, k) for st in range(len(starts)) for w in [None] + waits for k in ([0] if w is None else range(len(seconds)))]   # (None: untouched)
    ch = CH.index(e1['ch'])
    res = []
    for b0 in range(0, len(rows), a.batch):
        chunk = rows[b0:b0 + a.batch]
        idx = [st for st, _, _ in chunk]
        plan = {'action': ['none' if w is None else 'mask' for _, w, _ in chunk], 'ch': [ch] * len(chunk), 'r': [float(MASK_R)] * len(chunk),
                'dur': [e1['hold']] * len(chunk), 'off': [(e1['off'][0] * R, e1['off'][1] * R)] * len(chunk),
                'frame': [book['morphs'][A_]['frame']] * len(chunk), 'axis0': S['axis'][idx],
                'ch2': [CH.index(seconds[k]['ch']) for _, _, k in chunk],
                'at2': [10 ** 9 if w is None else e1['hold'] + w for _, w, _ in chunk], 'dur2': [seconds[k]['hold'] for _, _, k in chunk],
                'off2': [(seconds[k]['off'][0] * R, seconds[k]['off'][1] * R) for _, _, k in chunk], 'frame2': [book['morphs'][B_]['frame']] * len(chunk)}
        ser, out, _ = run_batch(sim, S['A'][idx], stack_kin(S['kin'], idx, R, sim.dt), S['com'][idx], plan, steps, R, models, a.debris)
        for (st, w, k), o, i in zip(chunk, out, range(len(chunk))):
            res.append({'start': st, 'wait': w, 'second': k, 'label': o['label'], 'velAng': float(circmean(ser['velAng'][i, -10:])),
                        'x': float(ser['x'][i, -1] - ser['x'][i, 0]), 'y': float(ser['y'][i, -1] - ser['y'][i, 0])})
    ctrl = {r['start']: r for r in res if r['wait'] is None}
    if a.second == 'all':   # per wait: the second actions that end best, and the ends overall
        summary = []
        for w in waits:
            rs = [r for r in res if r['wait'] == w]
            print(f'  wait {w}: ends {Counter(r["label"] for r in rs).most_common(6)}', flush=True)
            for k, e in enumerate(seconds):
                ks = [r for r in rs if r['second'] == k]
                c = Counter(r['label'] for r in ks)
                summary.append({'wait': w, 'ch': e['ch'], 'spot': e['spot'], 'hold': e['hold'], 'n': len(ks), 'ends': dict(c)})
            for tgt in sorted({r['label'] for r in rs if ' + ' not in r['label'] and r['label'] not in ('???',) + DEAD}):
                best = max((sm for sm in summary if sm['wait'] == w), key=lambda sm: sm['ends'].get(tgt, 0) / max(1, sm['n']))
                pv = best['ends'].get(tgt, 0) / max(1, best['n'])
                if pv >= 0.5: print(f'    -> {tgt}: best {best["ch"]} {best["spot"]} {best["hold"]} (p {pv:.2f})', flush=True)
        json.dump({'code': a.code, 'morph': A_, 'via': B_, 'first': e1, 'seconds': seconds, 'summary': summary, 'runs': res},
                  open(os.path.join(OUT, a.code, f'combo-{A_}-{B_}-all.json'), 'w'), default=jdefault)
        return
    print(f'  wait (steps) | ends as {A_} | net turn: median, spread (deg) | travel vs untouched (R)')
    summary = []
    for w in waits:
        rs = [r for r in res if r['wait'] == w]
        ok = [r for r in rs if r['label'] == A_]
        tv = [math.degrees(wrapa(r['velAng'] - ctrl[r['start']]['velAng'])) for r in ok]
        med = math.degrees(circmean(np.radians(tv))) if tv else float('nan')
        sd = float(np.degrees(np.sqrt(-2 * np.log(max(1e-9, abs(np.exp(1j * np.radians(tv)).mean())))))) if tv else float('nan')
        trav = np.median([(math.hypot(r['x'], r['y']) - math.hypot(ctrl[r['start']]['x'], ctrl[r['start']]['y'])) / R for r in rs])
        ends = Counter(r['label'] for r in rs).most_common(3)
        print(f'  {w:4d} | {len(ok)}/{len(rs)} {ends} | {med:+6.0f} ± {sd:3.0f} | {trav:+.1f}', flush=True)
        summary.append({'wait': w, 'ok': len(ok), 'n': len(rs), 'turn': med, 'turn_sd': sd, 'ends': ends, 'travel': float(trav)})
    json.dump({'code': a.code, 'morph': A_, 'via': B_, 'first': e1, 'second': e2, 'summary': summary, 'runs': res},
              open(os.path.join(OUT, a.code, f'combo-{A_}-{B_}.json'), 'w'), default=jdefault)


def best_action(book, morph, target):
    """the model's best action (pooled) from morph for outcome target ('into X', any turning;
    'keep'; 'dissolve'), or None: highest p (in 0.05), then the shorter hold"""
    m = book['morphs'].get(morph)
    if not m: return None
    hit = lambda ef: ef.split(' (')[0] == target or ef == target or ef.startswith(target)
    es = [e for e in m['actions'] if hit(e['effect']) and e['p'] >= 0.1]
    es += [{**e, 'effect': e['alt'][0], 'p': e['alt'][1]} for e in m['actions'] if e.get('alt') and hit(e['alt'][0]) and e['alt'][1] >= 0.1]   # (its runner-up: e.g. XEH4YR rectus -> tachus)
    return sorted(es, key=lambda e: (-round(e['p'] * 20), e['hold']))[0] if es else None


def effect_hits(ef, target):
    """an action's outcome is the one wanted: 'into X' (any turning), 'keep', 'dissolve', or by its start ('replicate')"""
    return ef.split(' (')[0] == target or ef == target or ef.startswith(target)


def board_phase(sig, period, dt):
    """an oscillator's phase now (rad) from its last swing of signal (turn rates or speeds, a step
    each): the sine fitted, as swing_phases / models.js swingBin (omega ~ sin(phase))"""
    y = np.asarray(sig[-period:], float)
    tt = np.arange(len(y)) * dt
    w0 = 2 * math.pi / (period * dt)
    X = np.stack([np.sin(w0 * tt), np.cos(w0 * tt), np.ones_like(tt)], 1)
    aa, bb, _ = np.linalg.lstsq(X, y, rcond=None)[0]
    return (w0 * tt[-1] + math.atan2(bb, aa)) % (2 * math.pi)


def phase_action(m, target, psi, dt, cost=0.1):
    """an oscillator's action for an outcome, by phase (its model's by_phase, alt_by_phase): the
    (action, phase bin, steps to wait) best for reliability less `cost` a game second of waiting
    (the reaction time, as behaviors.js reach), the wait to the middle of that phase's bin"""
    nb, P = m['phase_bins'], m['swing_period']
    w0 = 2 * math.pi / (P * dt)
    best, bs = None, -1e9
    for e in m['actions']:
        for b in range(nb):
            for x in ((e.get('by_phase') or [None] * nb)[b], (e.get('alt_by_phase') or [None] * nb)[b]):
                if not x or not effect_hits(x[0], target) or x[1] < 0.3: continue
                wait = int(round(((b + 0.5) * 2 * math.pi / nb - psi) % (2 * math.pi) / (w0 * dt)))
                sc = x[1] - cost * wait / 30 - 0.01 * e['hold'] / 30
                if sc > bs: bs, best = sc, ({**e, 'effect': x[0], 'p': x[1]}, b, wait)
    return best


def age_action(m, target, age):
    """a short-lived morph swept by age (VT049W gyrans): its action for an outcome at an age not yet
    passed (its model's by_phase: bin k the sweep's starts at age_settle + k / bins of age_span), the
    most reliable (in 0.05), then the shortest wait, then the shortest hold: (action, bin, steps to
    wait); None once older than the span (or none hits)"""
    nb, span, st = m['phase_bins'], m['age_span'], m['age_settle']
    now = (age - st) / span * nb
    if now > nb - 0.5: return None
    best, bk = None, None
    for e in m['actions']:
        for b in range(max(0, round(now)), nb):
            for x in ((e.get('by_phase') or [None] * nb)[b], (e.get('alt_by_phase') or [None] * nb)[b]):
                if not x or len(x) < 3 or x[2] < 3 or not effect_hits(x[0], target) or x[1] < 0.1: continue
                wait = max(0, int(round(st + b * span / nb - age)))
                k = (-round(x[1] * 20), wait, e['hold'])
                if bk is None or k < bk: bk, best = k, ({**e, 'effect': x[0], 'p': x[1]}, b, wait)
    return best


def parse_steps(text):
    """'into gyrans; wait 60; into rectus x3' -> [('into', 'gyrans', tries), ('wait', 60, 0), ...]"""
    out = []
    for part in [p.strip() for p in text.split(';') if p.strip()]:
        w = part.split()
        if w[0] == 'wait': out.append(('wait', int(w[1]), 0, None))
        elif w[0] == 'into':   # into X [xN] [from Y]: Y's model when it is not named (a short-lived stage, e.g. VT049W gyrans)
            frm = w[w.index('from') + 1] if 'from' in w else None
            out.append(('into', w[1], int(w[2][1:]) if len(w) > 2 and w[2].startswith('x') else 1, frm))
        elif w[0] == 'do':   # do <effect> [from Y]: as into's from (e.g. 'do split: tachus + tachus from gyrans')
            frm = w[w.index('from') + 1] if 'from' in w else None
            out.append(('do', ' '.join(w[1:w.index('from')] if frm else w[1:]), 1, frm))
        elif w[0] == 'hold':   # hold <spot words | fwd,side in R> <steps> [ch]: a mask held there, no model
            chs = w[-1] if w[-1] in CH else 'G'
            body = w[1:] if w[-1] not in CH else w[1:-1]   # (the spot's words, then the steps)
            out.append(('hold', (' '.join(body[:-1]), int(body[-1]), chs), 1, None))   # (the action for an outcome, e.g. 'do replicate: tachus', whatever it becomes)
        else: sys.exit(f'unknown step {part!r} (into X [xN] | do <effect> | wait N)')
    return out


def isolate(sim, A, com, idx, R, debris):
    """--isolate: boards idx cleared of all but their heaviest body (copies, split-off parts, debris),
    so the next action acts in a clear world, aimed at that body (a mask's place is from the board's
    CoM); com (unwrapped) moved with the board's CoM"""
    An = A.cpu().numpy()
    m0 = sim.measure(A)
    keep = np.ones(An.shape, np.float32)
    for i in idx:
        lab, nb = label_torus(An[i].sum(0), 0.02, 2)
        if nb <= 1: continue
        mass = np.bincount(lab.ravel(), weights=An[i].sum(0).ravel(), minlength=nb + 1)
        k = 1 + int(np.argmax(mass[1:]))
        keep[i] = (lab == k)[None].astype(np.float32)
    A = A * sim.t.tensor(keep, device=A.device)
    m1 = sim.measure(A)
    com = com + np.stack([wrapd(m1[:, 0] - m0[:, 0], sim.W), wrapd(m1[:, 1] - m0[:, 1], sim.H)], 1)
    return A, com


def candidates(m, target, age=None):
    """a morph's actions for an outcome, best first (as best_action: p in 0.05, then the hold; an
    age-swept morph's for ages to come, then the shortest wait): [(action, wait)], each spot and hold once"""
    out = []
    if m.get('phase_useful') and m.get('phase_kind') == 'age' and age is not None:
        nb, span, st = m['phase_bins'], m['age_span'], m['age_settle']
        now = (age - st) / span * nb
        if now > nb - 0.5: return []
        for e in m['actions']:
            for b in range(max(0, round(now)), nb):
                for x in ((e.get('by_phase') or [None] * nb)[b], (e.get('alt_by_phase') or [None] * nb)[b]):
                    if x and len(x) >= 3 and x[2] >= 3 and effect_hits(x[0], target) and x[1] >= 0.1:
                        w = max(0, int(round(st + b * span / nb - age)))
                        out.append(((-round(x[1] * 20), w, e['hold']), {**e, 'effect': x[0], 'p': x[1]}, w))
    else:
        for e in m['actions']:
            for ef, p in ((e['effect'], e['p']), tuple(e['alt']) if e.get('alt') else (None, 0)):
                if ef and effect_hits(ef, target) and p >= 0.1: out.append(((-round(p * 20), 0, e['hold']), {**e, 'effect': ef, 'p': p}, 0))
    seen, res = set(), []
    for _, e, w in sorted(out, key=lambda x: x[0]):
        k = (e['ch'], e['spot'], e['hold'])
        if k not in seen: seen.add(k); res.append((e, w))
    return res


def verify_acts(a, sim, A, kin, com, acts, info, kind, arg, last, R, models):
    """--verify K: each acting board's K best candidates (its chosen action first) simulated from the
    board as it is now (a copy each, as the step would run them: hold, then settle_step); the first
    whose outcome is the step's done instead (the most reliable that works here); none: the chosen"""
    rows, cand = [], []
    for i, (e0, off0, fr, ax, w0) in acts.items():
        cur, m, mir, ax, target, age = info[i]
        cs = [(e0, w0)] + [c for c in candidates(m, target, age) if (c[0]['ch'], c[0]['spot'], c[0]['hold']) != (e0['ch'], e0['spot'], e0['hold'])]
        for e, w in cs[:a.verify]:
            rows.append(i); cand.append((i, (e, (e['off'][0] * R, e['off'][1] * R * mir), m['frame'], ax, w)))
    nr = len(rows)
    plan = {'action': ['mask'] * nr, 'ch': [CH.index(c[1][0]['ch']) for c in cand], 'r': [float(MASK_R)] * nr,
            'dur': [0] * nr, 'off': [(0.0, 0.0)] * nr, 'frame': [c[1][2] for c in cand], 'axis0': np.array([c[1][3] for c in cand]),
            'at2': [c[1][4] for c in cand], 'dur2': [c[1][0]['hold'] for c in cand], 'off2': [c[1][1] for c in cand], 'frame2': [c[1][2] for c in cand]}
    dur = max(c[1][4] + c[1][0]['hold'] for c in cand) + a.settle_step
    _, out, _ = run_batch(sim, A[rows].clone(), kin.copy_rows(rows), com[rows].copy(), plan, dur, R, models, a.debris)
    targets = (a.target.split(',') if a.target and last else None)
    def ok(lb, cur):
        if targets: return any(lb == t or (t == 'many' and (' + ' in lb or lb in ('multiplied', 'exploded'))) for t in targets)
        return lb == arg if kind == 'into' else outcome_hits(arg, cur, lb)
    new = dict(acts)
    for i in acts:
        for (j, c), o in zip(cand, out):
            if j == i and ok(o['label'], info[i][0]): new[i] = c; break
    a.verify_stats = getattr(a, 'verify_stats', Counter())
    a.verify_stats.update('changed' if new[i][0] is not acts[i][0] else 'kept' for i in acts)
    return new


def cmd_composite(a):
    """A composite maneuver: a sequence of maneuvers, closed loop as a bot would do it: each 'into X' step
    takes the model's best action from whatever morph each board is now (mirrored for one turning
    the other way), holds it, lets it settle (--settle-step), and names what it became; failed
    tries are tried again (into X xN); 'wait N' does nothing for N steps. Over the morph's starts
    (phases x angles), with an untouched copy of each: success (each board's end named the last
    step's morph), the steps it took, the net turn and the travel against the untouched copy, and
    where the failures ended. --vals: values for a {w} in the steps (e.g. the wait between, which
    sets a turn's angle), each run. Saved as output/models/<code>/composite-<name>.json (once strategy-<name>.json: the
    command was called strategy, still an alias)."""
    book = json.load(open(os.path.join(BOOK, f'{a.code}.json')))
    R = R_GAME
    ext = morph_extent(a.code, a.morph, R)
    size = a.size or (128 if ext <= 40 else 160 if ext <= 52 else 192)
    sim = Sim(rule_of(a.code), size, size, R)
    models = detect_models(a.code)
    m_age = book['morphs'].get(a.morph, {})
    if m_age.get('phase_kind') == 'age':   # (an age-swept start: as its sweep, young: settled age_settle, the starts over its span)
        a.settle, a.no_shift = m_age['age_settle'], True
        print(f"  {a.morph}: by age: settled {a.settle} steps, starts over its {m_age['age_span']}-step span", flush=True)
    probe = start_set(sim, a.code, a.morph, [(0, 0)], a.settle + 240, R)
    period = probe['period']
    starts = [(q, (0 if a.no_shift else q * 13 * period) + round(k * period / a.phases)) for q in range(a.angles) for k in range(a.phases)]
    S = start_set(sim, a.code, a.morph, starts, a.settle, R)
    n = len(starts)
    results = []
    for val in ((a.vals.split('|') if '|' in a.vals else a.vals.split(',')) if a.vals else [None]):   # ('|' between values that hold commas: fwd,side)
        steps = parse_steps(a.steps.replace('{w}', val) if val is not None else a.steps)
        idx = list(range(n)) * 2   # (n boards doing it, n untouched)
        A = S['A'][idx].clone(); kin = stack_kin(S['kin'], idx, R, sim.dt); com = S['com'][idx].copy()
        label = [a.morph] * (2 * n); t_done = [None] * n; tries = [0] * n; path = [[a.morph] for _ in range(n)]
        x0 = com.copy(); t = 0; ok_all = [True] * n
        first_head = None
        # each board's birth as an age-swept morph (VT049W gyrans): a start's from its settling; one
        # reached on the way (e.g. rectus replicating through gyrans) from the start of the action
        # that made it (t_act: each board's last action's start)
        birth = [-(a.settle + starts[i % n][1]) if m_age.get('phase_kind') == 'age' else None for i in range(2 * n)]
        t_act = [0] * (2 * n)
        hist = {'omega': [], 'speed': [], 'mass': []}   # (each board's turn rate, speed, mass, a step each: an oscillator's phase; a replication's swell)
        def rec(ser):
            for k in hist: hist[k].append(ser[k])
        sig = lambda k: np.concatenate(hist[k], 1)
        m0 = book['morphs'].get(a.morph, {})
        if m0.get('phase_useful') and m0.get('phase_kind') == 'swing':   # (a swing first: each board's phase known before it acts)
            none = {'action': ['none'] * (2 * n), 'ch': [0] * (2 * n), 'r': [float(MASK_R)] * (2 * n), 'dur': [0] * (2 * n),
                    'off': [(0.0, 0.0)] * (2 * n), 'frame': ['vel'] * (2 * n), 'axis0': np.zeros(2 * n)}
            ser, out, A = run_batch(sim, A, kin, com, none, m0['swing_period'], R, models, a.debris)
            com = np.stack([ser['x'][:, -1], ser['y'][:, -1]], 1); t += m0['swing_period']; rec(ser)
        for kind, arg, ntry, frm in steps:
            if kind == 'wait' and arg <= 0: continue
            if kind == 'hold':   # (a mask held at a spot of the board's morph now, in its frame, for N steps)
                spot, nst, chn = arg
                if re.match(r'^-?[0-9.]+,-?[0-9.]+$', spot): off = tuple(float(v) * R for v in spot.split(','))
                else: off = next((tuple(v * R for v in e['off']) for e in book['morphs'][a.morph]['actions'] if e['spot'] == spot), None)   # (as the model's spot: its scale)
                if off is None: sys.exit(f'unknown spot {spot!r}')
                fr = book['morphs'][a.morph]['frame']
                plan = {'action': ['mask'] * n + ['none'] * n, 'ch': [CH.index(chn)] * (2 * n), 'r': [float(MASK_R)] * (2 * n), 'dur': [nst] * (2 * n),
                        'off': [off] * (2 * n), 'frame': [fr] * (2 * n), 'axis0': np.array([math.atan2(wrapd(mm[6] - mm[8], sim.H), wrapd(mm[5] - mm[7], sim.W)) for mm in sim.measure(A)])}
                ser, out, A = run_batch(sim, A, kin, com, plan, nst + a.settle_step, R, models, a.debris)
                com = np.stack([ser['x'][:, -1], ser['y'][:, -1]], 1); t += nst + a.settle_step; rec(ser)
                for i in range(2 * n): label[i] = out[i]['label']
                for i in range(n): tries[i] += 1; path[i].append(label[i])
                continue
            if kind == 'wait':
                plan = {'action': ['none'] * (2 * n), 'ch': [0] * (2 * n), 'r': [float(MASK_R)] * (2 * n), 'dur': [0] * (2 * n),
                        'off': [(0.0, 0.0)] * (2 * n), 'frame': ['vel'] * (2 * n), 'axis0': np.zeros(2 * n)}
                ser, out, A = run_batch(sim, A, kin, com, plan, arg, R, models, a.debris)
                com = np.stack([ser['x'][:, -1], ser['y'][:, -1]], 1); t += arg; rec(ser)
                if first_head is None: first_head = ser['velAng'][:, 0]
                continue
            for _ in range(max(1, ntry)):
                todo = [i for i in range(n) if ok_all[i] and (kind == 'do' or label[i] != arg)]
                if not todo: break
                if a.isolate: A, com = isolate(sim, A, com, todo, R, a.debris)
                acts, info = {}, {}
                meas = sim.measure(A)
                for i in todo:
                    cur = label[i] if label[i] in book['morphs'] or not frm else frm
                    label_i = cur
                    target = arg if kind == 'do' else f'into {arg}'
                    m = book['morphs'].get(label_i, {})
                    e, wait = best_action(book, cur, target), 0
                    # an oscillator whose phase matters: its phase now (its history), the action and the
                    # phase that are best less the wait for it; the mask held when that phase comes
                    if m.get('phase_useful') and m.get('phase_kind') == 'swing' and hist['omega'] and sig('omega').shape[1] >= m['swing_period']:
                        psi = board_phase(sig('omega' if m.get('phase_signal', 'omega') == 'omega' else 'speed')[i], m['swing_period'], sim.dt)
                        pa = phase_action(m, target, psi, sim.dt)
                        if pa: e, _, wait = pa
                    if m.get('phase_useful') and m.get('phase_kind') == 'age':   # (its age now: the action for an age to come)
                        if birth[i] is None: birth[i] = t_act[i]
                        pa = age_action(m, target, t - birth[i])
                        e, wait = (pa[0], pa[2]) if pa else (None, 0)
                    if not e: continue   # (no action from what it is now: not acted on; the step's end judges it)
                    mir = -1 if m['chirality'] and np.sign(kin.omega[i]) == -m['chirality'] else 1
                    ax = math.atan2(wrapd(meas[i, 6] - meas[i, 8], sim.H), wrapd(meas[i, 5] - meas[i, 7], sim.W))
                    acts[i] = (e, (e['off'][0] * R, e['off'][1] * R * mir), m['frame'], ax, wait)
                    info[i] = (cur, m, mir, ax, target, t - birth[i] if birth[i] is not None else None)
                if not acts: break
                if a.verify > 1:   # (check first: the model's K best candidates simulated from each board as it is now; the first that works done)
                    acts = verify_acts(a, sim, A, kin, com, acts, info, kind, arg, (kind, arg) == tuple(steps[-1][:2]), R, models)
                for i in acts: t_act[i] = t
                dur = max(x[4] + x[0]['hold'] for x in acts.values()) + a.settle_step
                # (each board's mask held from its wait on: run_batch's second hold)
                plan = {'action': ['mask' if i in acts else 'none' for i in range(2 * n)],
                        'ch': [CH.index(acts[i][0]['ch']) if i in acts else 0 for i in range(2 * n)], 'r': [float(MASK_R)] * (2 * n),
                        'dur': [0] * (2 * n), 'off': [(0.0, 0.0)] * (2 * n),
                        'frame': [acts[i][2] if i in acts else 'vel' for i in range(2 * n)],
                        'axis0': np.array([acts[i][3] if i in acts else 0.0 for i in range(2 * n)]),
                        'at2': [acts[i][4] if i in acts else 10 ** 9 for i in range(2 * n)],
                        'dur2': [acts[i][0]['hold'] if i in acts else 0 for i in range(2 * n)],
                        'off2': [acts[i][1] if i in acts else (0.0, 0.0) for i in range(2 * n)],
                        'frame2': [acts[i][2] if i in acts else 'vel' for i in range(2 * n)]}
                ser, out, A = run_batch(sim, A, kin, com, plan, dur, R, models, a.debris)
                com = np.stack([ser['x'][:, -1], ser['y'][:, -1]], 1); t += dur; rec(ser)
                if first_head is None: first_head = ser['velAng'][:, 0]
                for i in range(2 * n): label[i] = out[i]['label']
                for i in acts:
                    tries[i] += 1; path[i].append(label[i])
                    if label[i] == arg and t_done[i] is None and (kind, arg) == tuple(steps[-1][:2]): t_done[i] = t
            for i in range(n):
                if ok_all[i] and kind == 'into' and label[i] != arg: ok_all[i] = False
        # the last step's morph at the end, the net turn and travel against the untouched copy
        fin = run_batch(sim, A, kin, com, {'action': ['none'] * (2 * n), 'ch': [0] * (2 * n), 'r': [float(MASK_R)] * (2 * n), 'dur': [0] * (2 * n),
                                           'off': [(0.0, 0.0)] * (2 * n), 'frame': ['vel'] * (2 * n), 'axis0': np.zeros(2 * n)}, 30, R, models, a.debris)[0]
        va = np.array([circmean(fin['velAng'][i, -10:]) for i in range(2 * n)])
        trav = np.hypot(fin['x'][:, -1] - x0[:, 0], fin['y'][:, -1] - x0[:, 1]) / R
        targets = (a.target or (steps[-1][1] if steps[-1][0] == 'into' else a.morph)).split(',')
        # ('many': replicated, the copies left at the end, or the mass swollen to REP_PEAK times its start on
        # the way, as the models judge a replication: the copies may merge again)
        mass = sig('mass') if hist['mass'] else None
        swell = (mass.max(1) / np.maximum(mass[:, 0], 1e-9)) if mass is not None else np.zeros(2 * n)
        hit = lambda i, lb: any(lb == t or (t == 'many' and (lb in ('exploded', 'multiplied') or ' + ' in lb or swell[i] >= REP_PEAK)) for t in targets)
        ok = [i for i in range(n) if ok_all[i] and hit(i, label[i])]
        turns = [math.degrees(wrapa(va[i] - va[i + n])) for i in ok]
        res = {'val': val, 'steps': a.steps if val is None else a.steps.replace('{w}', val), 'n': n, 'ok': len(ok), 'p': round(len(ok) / n, 2),
               'time': int(np.median([t] * len(ok))) if ok else None, 'tries': round(float(np.mean([tries[i] for i in ok])), 1) if ok else None,
               'turn': round(math.degrees(circmean(np.radians(turns))), 1) if turns else None,
               'turn_sd': round(float(np.degrees(np.sqrt(-2 * np.log(max(1e-9, abs(np.exp(1j * np.radians(turns)).mean())))))), 1) if turns else None,
               'travel': round(float(np.median(trav[ok] - trav[[i + n for i in ok]])), 1) if ok else None,
               'fails': dict(Counter(label[i] for i in range(n) if i not in ok).most_common(4)),
               'paths': dict(Counter(' > '.join(path[i]) for i in range(n)).most_common(4))}
        if a.verify > 1: res['verify'] = dict(getattr(a, 'verify_stats', {})); a.verify_stats = Counter(); print(f"  (verify {a.verify}: tries {res['verify']})", flush=True)
        results.append(res)
        print(f"  {res['steps']}: p {res['p']} ({res['ok']}/{n}), {res['time']} steps, tries {res['tries']}, turn {res['turn']} ± {res['turn_sd']}, "
              f"travel {res['travel']:+} R vs untouched, fails {res['fails']}" if res['ok'] else f"  {res['steps']}: p 0, fails {res['fails']}", flush=True)
    name = a.name or re.sub(r'[^a-z0-9]+', '-', f'{a.morph} {a.steps}'.lower()).strip('-')
    slug = re.sub(r'[^a-z0-9]+', '-', name.lower()).strip('-')
    old = os.path.join(OUT, a.code, f'strategy-{slug}.json')
    if os.path.exists(old): os.remove(old)   # (the same composite run again: its new file replaces the old-named one)
    json.dump({'code': a.code, 'morph': a.morph, 'name': name, 'steps': a.steps, 'vals': a.vals, 'note': a.note, 'results': results},
              open(os.path.join(OUT, a.code, f"composite-{slug}.json"), 'w'), default=jdefault)   # (a name in words, e.g. 'limus turns sharply': its file a slug)


def cmd_props(a):
    """Each morph's motion over a long run, settled: its speed, the spin of its body axis (a still
    spinner: antigyrans), the wobble of its heading round its mean turn (an oscillator: zigzag),
    its period and the pulse of its mass. Saved as output/models/<code>/props.json (book reads it)."""
    rule = rule_of(a.code)
    names = a.morph.split(',') if a.morph else [m.group(1) for f in sorted(os.listdir(os.path.join(OUT, a.code))) for m in [re.match(r'steer-([a-z+]+)\.json$', f)] if m]
    out = {}
    for morph in names:
        ext = morph_extent(a.code, morph, R_GAME)
        size = 128 if ext <= 40 else 160 if ext <= 52 else 192
        sim = Sim(rule, size, size, R_GAME)
        S = start_set(sim, a.code, morph, [(0, 0)], a.settle, R_GAME)
        A = S['A']
        meas = sim.measure(A)
        com, ax, mass = [np.zeros(2)], [], []
        for _ in range(a.steps):
            A = sim.step(A)
            m2 = sim.measure(A)
            com.append(com[-1] + [wrapd(m2[0, 0] - meas[0, 0], sim.W), wrapd(m2[0, 1] - meas[0, 1], sim.H)])
            ax.append(math.atan2(wrapd(m2[0, 6] - m2[0, 8], sim.H), wrapd(m2[0, 5] - m2[0, 7], sim.W)))
            mass.append(m2[0, 2:5].sum()); meas = m2
        com, mass = np.array(com), np.array(mass)
        T = a.steps * sim.dt
        L = 10   # heading: the displacement over L steps
        d = com[L:] - com[:-L]
        speed = float(np.hypot(*d.T).mean() / (L * sim.dt) / R_GAME)
        moving = speed > SPEED_MIN_R * 1.5
        hd = np.unwrap(np.arctan2(d[:, 1], d[:, 0]))
        tt = np.arange(len(hd)) * sim.dt
        k = np.polyfit(tt, hd, 1) if moving else [0.0, 0.0]
        res = hd - np.polyval(k, tt) if moving else np.zeros(1)
        wobble = float(np.degrees(res.std()))
        per = 0
        if moving and wobble > 1:
            f = np.abs(np.fft.rfft(res - res.mean())); f[0] = 0
            i = int(np.argmax(f)); per = int(round(len(res) / i)) if i else 0
        axu = np.unwrap(np.array(ax))
        spin = float((axu[-1] - axu[0]) / T)
        net = float(np.hypot(*(com[-1] - com[0])) / T / R_GAME)   # (travel: a jitter to and fro has none)
        o = {'speed': round(speed, 4), 'net_speed': round(net, 4), 'moving': bool(moving), 'turn': round(float(k[0]), 4), 'spin': round(spin, 4),
             'wobble_deg': round(wobble, 1), 'wobble_period': per, 'mass_pulse': round(float(mass.std() / mass.mean()), 3),
             'mirrored': S['mirrored']}
        out[morph] = o
        print(a.code, morph, o, flush=True)
    path = os.path.join(OUT, a.code, 'props.json')
    old = json.load(open(path)) if os.path.exists(path) else {}
    old.update(out)
    json.dump(old, open(path, 'w'), indent=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--code', required=True)
    common.add_argument('--R', type=int, default=R_GAME)
    common.add_argument('--settle', type=int, default=300)
    common.add_argument('--debris', type=float, default=0.3, help='bodies lighter than this (mass / R^2) are ignored')
    s = sub.add_parser('steer', parents=[common])
    s.add_argument('--morph', required=True, help='comma-separated')
    s.add_argument('--channels', default='G', help='mask channels, e.g. G or B or RGB')
    s.add_argument('--durs', type=lambda v: [int(x) for x in v.split(',')], default=[15, 30, 60], help='holds, steps (30 = 1 s)')
    s.add_argument('--spots', default='named', help="'named' (centre + 8 directions x --dists) or 'grid'")
    s.add_argument('--dists', default='body,edge', help='named spot distances: body (0.45 R), edge (0.9 R), outside (1.3 R)')
    s.add_argument('--spot-scale', type=float, default=1.0, help='the spot distances x this, for a larger morph (a twin: 1.6); 0: by its width')
    s.add_argument('--spacing', type=int, default=4, help='grid spacing, cells')
    s.add_argument('--extent', type=float, default=1.6, help='grid within this many R of the CoM')
    s.add_argument('--phases', type=int, default=12, help='start phases over the morph\'s cycle')
    s.add_argument('--period', type=int, default=0, help='the cycle (steps) to spread the start phases over (0: the detected one)')
    s.add_argument('--phase-signal', default='omega', choices=['omega', 'speed'], help="an oscillator's phase from its turn rate's cycle or its speed's")
    s.add_argument('--like', action='store_true', help="as the morph's model sweep (its starts and settings exactly), other --channels, saved under --tag (default 'more'): merged into the model by book")
    s.add_argument('--no-shift', action='store_true', help='the turned copies not settled whole cycles longer (an unstable morph, e.g. PJ1EF7 lithos, decays meanwhile)')
    s.add_argument('--extra-spots', default='', help="spots of one's own, 'name:fwd,side;...' in R (e.g. 'centre-back 0.2:-0.2,0')")
    s.add_argument('--only-extra', action='store_true', help='only the --extra-spots (e.g. a --like sweep of a few more spots, merged by book)')
    s.add_argument('--made', default='', help="starts made from another morph by an action: 'parent:ch:spot:hold:wait' (a short-lived morph whose captured pattern wanders off; the parent's model gives the spot)")
    s.add_argument('--tag', default='', help="an experiment's name: saved as steer-<morph>-<tag> (the book leaves it out)")
    s.add_argument('--age', action='store_true', help='a short-lived morph (VT049W gyrans): starts over its first --period steps of life (its age), not over cycles')
    s.add_argument('--states', type=int, default=0, help='a chaotic morph: starts from this many states of a long run (k-means of speed and turn rate)')
    s.add_argument('--probe', type=int, default=6000, help='--states: the long run, steps')
    s.add_argument('--angles', type=int, default=4, help='start angles (quarter turns)')
    s.add_argument('--observe', type=int, default=400)
    s.add_argument('--short', type=int, default=60, help='turn measured this many steps after the hold')
    s.add_argument('--size', type=int, default=128)
    s.add_argument('--batch', type=int, default=384)
    c = sub.add_parser('clone', parents=[common])
    c.add_argument('--pair', required=True, help='ship,clone[;ship,clone...]')
    c.add_argument('--dists', type=lambda v: [float(x) for x in v.split(',')], default=[1.0, 1.4, 1.8, 2.2, 2.6])
    c.add_argument('--bearings', type=int, default=16)
    c.add_argument('--reps', type=int, default=2)
    c.add_argument('--observe', type=int, default=1200)
    c.add_argument('--size', type=int, default=192)
    c.add_argument('--batch', type=int, default=192)
    b = sub.add_parser('book', parents=[common])
    b.add_argument('--no-maps', action='store_true')
    b.add_argument('--atlas', action='store_true', help='keep every action\'s outcome in the JSON')
    sg = sub.add_parser('composite', aliases=['strategy'], parents=[common])
    sg.add_argument('--morph', required=True, help='the morph it starts as')
    sg.add_argument('--steps', required=True, help="'into gyrans; wait {w}; into rectus x3'")
    sg.add_argument('--vals', default='', help='values for {w} in the steps, comma-separated: a run each')
    sg.add_argument('--name', default=''); sg.add_argument('--note', default=''); sg.add_argument('--target', default='', help="success: ends as one of these, comma-separated; 'many': replicated (default: the last into step)")
    sg.add_argument('--settle-step', type=int, default=60, help='steps after a hold before naming what it became')
    sg.add_argument('--phases', type=int, default=6); sg.add_argument('--angles', type=int, default=4)
    sg.add_argument('--no-shift', action='store_true'); sg.add_argument('--size', type=int, default=0)
    sg.add_argument('--verify', type=int, default=0, help="check first: each try, the model's K best actions simulated from the board as it is now, the first that works done (K 1: off)")
    sg.add_argument('--isolate', action='store_true', help='before each action, clear each board of all but its heaviest body (copies, split-off parts): a clear world, the mask aimed at it')
    co = sub.add_parser('combo', parents=[common])
    co.add_argument('--morph', required=True); co.add_argument('--via', required=True)
    co.add_argument('--waits', default='0,15,30,45,60,90,120,150', help='steps between the two holds')
    co.add_argument('--first-effect', default='', help="the first action by its outcome (default: into --via), e.g. 'replicate: tachus'")
    co.add_argument('--second', default='best', choices=['best', 'all'], help="the second action: the model's best back, or every spot of --via's at each of --holds2")
    co.add_argument('--holds2', default='15,30')
    co.add_argument('--chs2', default='G', help='the second hold\'s channels (--second all), e.g. RGB')
    co.add_argument('--phases', type=int, default=4); co.add_argument('--angles', type=int, default=4)
    co.add_argument('--observe', type=int, default=300); co.add_argument('--size', type=int, default=0); co.add_argument('--batch', type=int, default=384)
    cy = sub.add_parser('cycles', parents=[common])
    cy.add_argument('--morph', required=True)
    cy.add_argument('--steps', type=int, default=4000)
    pr = sub.add_parser('props', parents=[common])
    pr.add_argument('--morph', default='', help="comma-separated (default: the morphs swept)")
    pr.add_argument('--steps', type=int, default=1500)
    sub.add_parser('rename', parents=[common])
    rb = sub.add_parser('robust', parents=[common])
    rb.add_argument('--morph', required=True); rb.add_argument('--pmin', type=float, default=0.75)
    rb.add_argument('--dhold', type=int, default=3, help='the hold this many steps shorter and longer')
    rb.add_argument('--starts', type=int, default=8, help="how many of the sweep's starts")
    vf = sub.add_parser('verify', parents=[common])
    vf.add_argument('--morph', default='', help='comma-separated (default: all with a sweep)')
    vf.add_argument('--pmin', type=float, default=0.75)
    vf.add_argument('--dhold', type=int, default=3); vf.add_argument('--dang', type=float, default=15)
    vf.add_argument('--shift', type=int, default=37, help="steps added to the sweep's start ages (fresh states)")
    vf.add_argument('--starts', type=int, default=12); vf.add_argument('--batch', type=int, default=384)
    w = sub.add_parser('show', parents=[common])
    w.add_argument('--morph', required=True)
    w.add_argument('--spot', required=True, help="e.g. 'right edge', 'tail-left body', 'centre'")
    w.add_argument('--hold', type=int, default=30)
    w.add_argument('--ch', default='G')
    w.add_argument('--start', type=int, default=0, help="which of the sweep's starts (phase x angle)")
    w.add_argument('--after', type=int, default=300)
    w.add_argument('--every', type=int, default=2)
    w.add_argument('--size', type=int, default=0, help="the board (0: the sweep's)")
    a = ap.parse_args()
    {'steer': cmd_steer, 'clone': cmd_clone, 'book': cmd_book, 'show': cmd_show, 'props': cmd_props, 'cycles': cmd_cycles, 'combo': cmd_combo, 'composite': cmd_composite, 'strategy': cmd_composite, 'robust': cmd_robust, 'verify': cmd_verify, 'rename': cmd_rename}[a.cmd](a)


if __name__ == '__main__':
    main()
