"use strict";
// ============================================================================================
//  Rule maths for the GPU engine, kept free of WebGL so it can be checked in node.
//
//  Pure functions, no WebGL: RLE decoding of Lenia cells, nearest-neighbour resampling, and
//  the packed kernel tables that sim.glsl consumes. Packing: the (up to) 16 kernels are four
//  vec4 groups. Per tap, sim.glsl does ONE fetch of the board (rgb = the three channels) and
//  forms each group's source values with a 3x4 selection matrix (SRC[g] * rgb), so a kernel's
//  source channel costs no extra fetch. Growth goes back to channels through 4x3 matrices
//  (DST[g] * (h * G)), then is divided by the per-channel sum of h (LeniaNDKC's normalisation).
// ============================================================================================

const Rule = (function () {
  const GROUPS = 4;               // 4 x vec4 = 16 kernel slots
  const SLOTS = GROUPS * 4;

  // ---- RLE (LeniaNDKC Board.rle2cells, 2D only) ------------------------------------------
  function ch2val(c) {
    if (c === '.' || c === 'b') return 0;
    if (c === 'o') return 255;
    if (c.length === 1) return c.charCodeAt(0) - 65 + 1;
    return (c.charCodeAt(0) - 112) * 24 + (c.charCodeAt(1) - 65 + 25);
  }

  function rle2cells(st) {
    const rows = [];
    let row = [], last = '', count = '';
    st = st.replace(/!$/, '') + '$';
    for (const ch of st) {
      if (ch >= '0' && ch <= '9') count += ch;
      else if ('pqrstuvwxy@'.includes(ch)) last = ch;
      else {
        const n = count ? parseInt(count, 10) : 1;
        if (ch === '$') {
          rows.push(row);
          for (let i = 1; i < n; i++) rows.push([]);
          row = [];
        } else {
          const v = ch2val(last + ch) / 255;
          for (let i = 0; i < n; i++) row.push(v);
        }
        last = ''; count = '';
      }
    }
    const w = Math.max(...rows.map(r => r.length));
    return rows.map(r => r.concat(new Array(w - r.length).fill(0)));
  }

  // Nearest-neighbour zoom, same index rule as the numpy port (floor(i / z)).
  function resample(A, z) {
    const h = A.length, w = A[0].length;
    const H2 = Math.max(1, Math.round(h * z)), W2 = Math.max(1, Math.round(w * z));
    const out = [];
    for (let y = 0; y < H2; y++) {
      const sy = Math.min(Math.floor(y / z), h - 1), r = [];
      for (let x = 0; x < W2; x++) r.push(A[sy][Math.min(Math.floor(x / z), w - 1)]);
      out.push(r);
    }
    return out;
  }

  // Decode a form's three channel RLEs, resample to radius R, and pad all to one shape.
  // A form rotated by `angle` radians about its centre, nearest-neighbour (each output cell
  // takes the source cell its centre maps back onto), in a box just large enough.
  function rotate(form, angle) {
    const c = Math.cos(angle), s = Math.sin(angle);
    const w = Math.ceil(Math.abs(form.w * c) + Math.abs(form.h * s));
    const h = Math.ceil(Math.abs(form.w * s) + Math.abs(form.h * c));
    const cells = form.cells.map(src => {
      const out = [];
      for (let y = 0; y < h; y++) {
        const r = new Array(w).fill(0);
        for (let x = 0; x < w; x++) {
          const dx = x + 0.5 - w / 2, dy = y + 0.5 - h / 2;
          const sx = Math.floor(c * dx + s * dy + form.w / 2), sy = Math.floor(-s * dx + c * dy + form.h / 2);
          if (sx >= 0 && sy >= 0 && sx < form.w && sy < form.h) r[x] = src[sy][sx];
        }
        out.push(r);
      }
      return out;
    });
    return { w, h, cells };
  }

  function formCells(rles, nativeR, R) {
    const z = R / nativeR;
    const ch = rles.map(s => resample(rle2cells(s), z));
    const h = Math.max(...ch.map(c => c.length));
    const w = Math.max(...ch.map(c => c[0].length));
    return { h, w, cells: ch.map(c => {
      const out = [];
      for (let y = 0; y < h; y++) {
        const r = [];
        for (let x = 0; x < w; x++) r.push(y < c.length && x < c[0].length ? c[y][x] : 0);
        out.push(r);
      }
      return out;
    }) };
  }

  // ---- kernels --------------------------------------------------------------------------
  function frac(s) {
    const [a, b] = s.split('/');
    return b === undefined ? parseFloat(a) : parseFloat(a) / parseFloat(b);
  }

  // Kernel shell, LeniaF.kernel_shell. d = distance / R. Core kn (LeniaF numbering):
  // 1 = polynomial quad4 (4f(1-f))^4, 2 = exponential bump exp(4 - 1/(f(1-f))).
  function shell(d, k, kn = 1) {
    const b = k.b.split(',').map(frac), B = b.length;
    if (d >= k.r) return 0;
    const Br = B * d / k.r;
    const f = Math.min(Br % 1, 1);
    if (f <= 0 || f >= 1) return 0;
    const core = kn === 1 ? Math.pow(4 * f * (1 - f), 4) : Math.exp(4 - 1 / (f * (1 - f)));
    return core * b[Math.min(Math.floor(Br), B - 1)];
  }

  // LeniaF's free_h, applied when it loads an old-format entry (no 'model'): each h becomes
  // round(h / sum of h over kernels with the same destination, 3); if all are <= 1/2, h and T
  // are both doubled. The update is then A += (1/T) * sum h G, with no further normalisation.
  function leniaF(kernels, T) {
    const Dn = [0, 0, 0];
    for (const k of kernels) Dn[k.c[1]] += k.h;
    let ks = kernels.map(k => ({ ...k, h: Math.round(k.h / Dn[k.c[1]] * 1000) / 1000 }));
    if (Math.max(...ks.map(k => k.h)) <= 0.5) { ks = ks.map(k => ({ ...k, h: k.h * 2 })); T *= 2; }
    return { kernels: ks, T };
  }

  // Tap table: every integer offset where any kernel is non-zero. weights[t*SLOTS + k] is
  // kernel k's normalised weight at tap t. Slots >= kernels.length stay zero (h = 0 too).
  function buildTaps(kernels, R, kn = 1) {
    const span = Math.ceil(R);
    const offs = [], raw = [];
    const sums = new Float64Array(SLOTS);
    for (let dy = -span; dy <= span; dy++) {
      for (let dx = -span; dx <= span; dx++) {
        const d = Math.sqrt(dx * dx + dy * dy) / R;
        const w = new Float64Array(SLOTS);
        let any = false;
        kernels.forEach((k, i) => { w[i] = shell(d, k, kn); if (w[i] > 0) any = true; });
        if (!any) continue;
        offs.push([dx, dy]); raw.push(w);
        for (let i = 0; i < SLOTS; i++) sums[i] += w[i];
      }
    }
    const weights = new Float32Array(offs.length * SLOTS);
    raw.forEach((w, t) => {
      for (let i = 0; i < SLOTS; i++) weights[t * SLOTS + i] = sums[i] > 0 ? w[i] / sums[i] : 0;
    });
    return { count: offs.length, offs, weights };
  }

  // Selection matrices and growth parameters, packed per vec4 group.
  //   src[g]: 12 floats, column-major mat3x4 (3 columns = channels, 4 rows = slots)
  //   dst[g]: 12 floats, column-major mat4x3 (4 columns = slots, 3 rows = channels)
  function buildParams(kernels) {
    const src = [], dst = [], m = [], s = [], h = [];
    const hsum = [0, 0, 0];
    for (let g = 0; g < GROUPS; g++) {
      const S = new Float32Array(12), D = new Float32Array(12);
      const mg = [0, 0, 0, 0], sg = [1, 1, 1, 1], hg = [0, 0, 0, 0];
      for (let j = 0; j < 4; j++) {
        const k = kernels[g * 4 + j];
        if (!k) continue;
        S[k.c[0] * 4 + j] = 1;           // column c0, row j
        D[j * 3 + k.c[1]] = 1;           // column j, row c1
        mg[j] = k.m; sg[j] = k.s; hg[j] = k.h;
        hsum[k.c[1]] += k.h;
      }
      src.push(S); dst.push(D); m.push(mg); s.push(sg); h.push(hg);
    }
    return { src, dst, m, s, h, hsum };
  }

  return { GROUPS, SLOTS, rle2cells, resample, formCells, rotate, shell, leniaF, buildTaps, buildParams };
})();

if (typeof module !== 'undefined') module.exports = Rule;
