"use strict";
// ============================================================================================
//  Finding and following creatures (Phase 1b, PLAN.md 4.5). Pure functions, no WebGL.
//
//  Input is blocks.glsl's coarse map: per BLOCK x BLOCK patch, (mass, sum m*dx, sum m*dy, P),
//  P = the potential K * A at the patch centre (mean over kernels), how the rule sees the world.
//  Segmentation follows ~/github/complexity (measures/individuality.py: basins, merge_shallow,
//  detect), on the block grid and the torus:
//    1. every block points to its highest 8-neighbour (itself if none is higher); following the
//       pointers ends at a peak of P, and each peak's catchment basin is one candidate;
//       blocks with P <= threshold are background;
//    2. neighbouring basins merge when the saddle between them (the highest point of their
//       common border) is at least `saddle` x the lower peak: bumps on one creature, not two;
//    3. basins smaller than minArea cells are dropped.
//  Unlike connected components of the mass, this separates touching creatures (a colony, or a
//  replicator just after it divides) whose potentials still have separate peaks.
//  A creature's CoM is exact: rebuilt from the per-block mass moments, with the basin unwrapped
//  across the board edges.
// ============================================================================================

const Creatures = (function () {
  const NB8 = [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, -1], [1, -1], [-1, 1]];
  const DEFAULTS = { threshold: 0.05, saddle: 0.8, minArea: 0 };

  // Basin labels per block (the peak's block index), -1 for background (step 1).
  function basins(P, BW, BH, thr) {
    const n = BW * BH, ptr = new Int32Array(n);
    for (let b = 0; b < n; b++) {
      const x = b % BW, y = (b - x) / BW;
      let best = b, bu = P[b];
      for (const [dx, dy] of NB8) {
        const nb = ((y + dy + BH) % BH) * BW + (x + dx + BW) % BW;
        if (P[nb] > bu) { bu = P[nb]; best = nb; }
      }
      ptr[b] = best;
    }
    for (let changed = true; changed;) {   // pointer jumping to the peaks
      changed = false;
      for (let b = 0; b < n; b++) {
        const q = ptr[ptr[b]];
        if (q !== ptr[b]) { ptr[b] = q; changed = true; }
      }
    }
    const lab = new Int32Array(n);
    for (let b = 0; b < n; b++) lab[b] = P[b] > thr ? ptr[b] : -1;
    return lab;
  }

  // Basins merged across saddles >= saddle x the lower peak, highest saddles first; a merged
  // basin keeps the higher peak (step 2). Relabels lab in place.
  function mergeShallow(P, lab, BW, BH, saddle) {
    const n = BW * BH, pairs = new Map();
    for (let b = 0; b < n; b++) {
      const la = lab[b];
      if (la < 0) continue;
      const x = b % BW, y = (b - x) / BW;
      for (const [dx, dy] of [[1, 0], [0, 1], [1, 1], [-1, 1]]) {
        const nb = ((y + dy + BH) % BH) * BW + (x + dx + BW) % BW, lb = lab[nb];
        if (lb < 0 || lb === la) continue;
        const key = la < lb ? la * n + lb : lb * n + la, h = Math.min(P[b], P[nb]);
        if (!(pairs.get(key) >= h)) pairs.set(key, h);
      }
    }
    if (!pairs.size) return;
    const parent = new Map();
    const find = i => {
      let r = i;
      while (parent.has(r) && parent.get(r) !== r) r = parent.get(r);
      while (i !== r) { const nx = parent.get(i); parent.set(i, r); i = nx; }
      return r;
    };
    for (const [key, h] of [...pairs].sort((u, v) => v[1] - u[1])) {
      const rx = find(Math.floor(key / n)), ry = find(key % n);   // labels are peak blocks
      if (rx === ry || h < saddle * Math.min(P[rx], P[ry])) continue;
      if (P[rx] >= P[ry]) parent.set(ry, rx); else parent.set(rx, ry);
    }
    for (let b = 0; b < n; b++) if (lab[b] >= 0) lab[b] = find(lab[b]);
  }

  // blk: Float32Array of BW*BH*4. opts: { threshold, saddle, minArea (cells) }.
  // Returns { comps, labels, total }:
  //   comps[i] = { mass, col, row, radius }  (radius = radius of gyration, cells)
  //   labels[b] = index into comps, or -1 (background, or a basin smaller than minArea)
  function find(blk, BW, BH, B, W, H, opts) {
    const { threshold, saddle, minArea } = { ...DEFAULTS, ...opts };
    const n = BW * BH;
    const P = new Float32Array(n);
    let total = 0;
    for (let b = 0; b < n; b++) { P[b] = blk[b * 4 + 3]; total += blk[b * 4]; }
    const lab = basins(P, BW, BH, threshold);
    mergeShallow(P, lab, BW, BH, saddle);

    const labels = new Int32Array(n).fill(-1);
    const seen = new Uint8Array(n);
    const ux = new Int32Array(n), uy = new Int32Array(n);
    const comps = [], stack = [], members = [];
    for (let b0 = 0; b0 < n; b0++) {
      if (seen[b0] || lab[b0] < 0) continue;
      const id = lab[b0];
      seen[b0] = 1;
      ux[b0] = b0 % BW; uy[b0] = Math.floor(b0 / BW);
      stack.length = 0; members.length = 0;
      stack.push(b0);
      let m = 0, sx = 0, sy = 0;
      while (stack.length) {   // the basin's blocks, unwrapped (8-connected, as the basins are)
        const b = stack.pop();
        members.push(b);
        const bm = blk[b * 4];
        m += bm;
        sx += bm * ux[b] * B + blk[b * 4 + 1];
        sy += bm * uy[b] * B + blk[b * 4 + 2];
        const x = b % BW, y = (b - x) / BW;
        for (const [dx, dy] of NB8) {
          const nx = (x + dx + BW) % BW, ny = (y + dy + BH) % BH, nb = ny * BW + nx;
          if (seen[nb] || lab[nb] !== id) continue;
          seen[nb] = 1;
          ux[nb] = ux[b] + dx; uy[nb] = uy[b] + dy;
          stack.push(nb);
        }
      }
      if (members.length * B * B < minArea || m <= 1e-6) continue;
      const cx = sx / m, cy = sy / m;
      let r2 = 0;
      for (const b of members) {
        const dx = ux[b] * B + B / 2 - cx, dy = uy[b] * B + B / 2 - cy;
        r2 += blk[b * 4] * (dx * dx + dy * dy);
      }
      const idx = comps.length;
      for (const b of members) labels[b] = idx;
      comps.push({ mass: m, col: ((cx % W) + W) % W, row: ((cy % H) + H) % H, radius: Math.sqrt(r2 / m) });
    }
    return { comps, labels, total };
  }

  const wrapD = (d, n) => d - n * Math.round(d / n);

  function dist(c, col, row, W, H) {
    return Math.hypot(wrapD(c.col - col, W), wrapD(c.row - row, H));
  }

  function labelAt(found, col, row, BW, BH, B) {
    const bx = ((Math.floor(col / B) % BW) + BW) % BW, by = ((Math.floor(row / B) % BH) + BH) % BH;
    return found.labels[by * BW + bx];
  }

  // The component at (col, row), or else the nearest one within `gate` cells; -1 if none.
  // Used for clicks and for tracking. On a split the predicted point usually lies inside one
  // half (it keeps the focus); if it falls in the gap, the nearer half wins.
  function pick(found, col, row, gate, BW, BH, B, W, H) {
    const at = labelAt(found, col, row, BW, BH, B);
    if (at >= 0) return at;
    let best = -1, bd = gate;
    found.comps.forEach((c, i) => {
      const d = dist(c, col, row, W, H);
      if (d < bd) { bd = d; best = i; }
    });
    return best;
  }

  return { find, pick, dist };
})();

if (typeof module !== 'undefined') module.exports = Creatures;
