"use strict";
// ============================================================================================
//  Morph detection (Detect button): which morph of the species a creature is, from its short
//  history compared with the long history's pieces of each morph.
//
//  Each morph's long-history samples are a cloud of points in a feature space: log(R, G, B mass + MASS0),
//  log(speed + SPEED0) and log(turn rate + OMEGA0), each divided by a fixed kernel width (WIDTHS).
//  Logs make a mass difference relative, and the offsets make speed and turn rate most sensitive
//  near zero, so a still creature (lithos) and a slowly drifting one (tardus) are far apart even
//  where their masses are alike. Widths are fixed, not from the data's spread: a still morph's
//  cloud is nearly a point, and what has to be covered is how much two copies differ (a copy added
//  at another angle has masses ~1-2% off), not how much one copy varies. A creature's recent
//  samples are scored against each cloud with a Gaussian kernel density (a sample's score = log of
//  the mean kernel over the cloud's points), and the mean score per morph becomes a probability
//  with a softmax. Only where the points lie counts, not their order in time, so it does not
//  matter where in its cycle a creature is.
//
//  Fuzzy morphs (pieces marked "fuzzy" in the long-history file, e.g. VT049W cytos: the
//  nearly-still, slowly shape-shifting offspring of the replicating gyrans) vary too much from one
//  creature to the next for the fixed widths. Their kernel width per feature is taken from their
//  own spread instead (FUZZ x its standard deviation, at least the fixed width), so their cloud is
//  broad. To compare clouds of different widths, the softmax uses each morph's normalised density
//  (score minus the sum of the log widths): a broad cloud pays for its breadth, so where a precise
//  morph fits, it wins, and a fuzzy one is named only where no precise one is near. Its name is
//  shown with a "~" (e.g. "~cytos"). Tested on VT049W's six cytos pieces, each left out of the
//  reference in turn: every 50-step window of it named ~cytos (FUZZ 0.5: 67%); rectus, lithos and
//  tachus still named right, and each of them, left out of the reference, stays ??? rather than
//  ~cytos (at FUZZ 1.25 cytos swallows lithos and tachus).
//
//  Unknown (null) when: too few samples since the last gap; the best morph's probability is under
//  P_MIN; or the samples fit no cloud (mean score under FIT_MIN, e.g. a morph never recorded).
//  Tested offline on the saved histories (each morph's first half against its second): 91% of
//  50-step windows named, 99.3% of those right; the misses are the shapeshifter (HGMLH8
//  metamorpha), which really passes through its neighbours' shapes. A morph left out of the
//  reference is wrongly named 10% of the time. Lithos and tardus (still vs slowly drifting, masses
//  alike) stay apart with every mass 5% off; with the earlier scaling (each stat by its spread
//  over the species), a lithos 2% low in red mass was taken for tardus or unknown.
// ============================================================================================

const MorphDetect = (function () {
  const F = ['mR', 'mG', 'mB', 'speed', 'omega'];   // total mass left out: it is their sum
  const MASS0 = 0.01;      // g: so an empty channel (a two-channel morph, e.g. XEH4YR cor) has a log
  const SPEED0 = 0.002;    // m/s: a still creature's measured speed is ~0.0001-0.0003
  const OMEGA0 = 0.01;     // rad/s
  const WIDTHS = [0.03, 0.03, 0.03, 0.3, 0.3];   // kernel width per feature (masses: 3%)
  const FUZZ = 0.75;        // a fuzzy morph's kernel width: this x its own spread, per feature
  // The feature vector of a sample x (values in F order), not yet divided by the widths.
  const feat = x => [Math.log(x[0] + MASS0), Math.log(x[1] + MASS0), Math.log(x[2] + MASS0),
                     Math.log(x[3] + SPEED0), Math.log(x[4] + OMEGA0)];
  const usable = x => x.every(v => v >= 0);   // NaN (no focus) fails too
  const NEFF = 2;          // softmax gain: the window's samples are far from independent
  const P_MIN = 0.9;       // a morph is named only at this probability or more
  const FIT_MIN = -6;      // mean score below this: fits no morph (a lone ref point scores ~-5.5)
  const CLAMP = -20;       // a sample's score is floored here, so one stray sample cannot veto
  const REF_MAX = 250;     // reference points per morph (strided beyond that)
  const WARM = 10;         // samples skipped after a gap: the speed and turn filters settle
  const MIN_N = 25;        // samples needed (after WARM) to name a morph
  const WIN = 50;          // newest samples used

  // HSV hue in degrees of (r, g, b) masses (as statschart.js's derived hue); NaN if grey
  function hueOf(r, g, b) {
    const mx = Math.max(r, g, b), d = mx - Math.min(r, g, b);
    if (!(d > 1e-9)) return NaN;
    const h = mx === r ? (g - b) / d : mx === g ? 2 + (b - r) / d : 4 + (r - g) / d;
    return (h * 60 + 360) % 360;
  }

  // The model of a long history (StatsChart.log): per morph its reference points, kernel widths
  // (as 1 / width), whether it is fuzzy, and its mean hue. null if it has fewer than two morphs
  // (nothing to choose between).
  function build(long) {
    // keyed by name, "~" before a fuzzy one: a fuzzy piece named like a precise morph stays apart
    const by = {}, fuzzy = {};
    long.segs.forEach((sg, i) => {
      const end = i + 1 < long.segs.length ? long.segs[i + 1].start : long.len;
      const key = (sg.fuzzy ? '~' : '') + sg.morph;
      // samples since the piece's gap: find where it started, to skip the warm-up
      let g = i; while (g > 0 && !long.segs[g].gap) g--;
      const p0 = long.segs[g].start + WARM;
      for (let j = Math.max(sg.start, p0); j < end; j++) {
        const x = F.map(k => long.get(k, j));
        if (usable(x)) (by[key] ||= []).push(x);
      }
      if (sg.fuzzy) fuzzy[key] = true;
    });
    const names = Object.keys(by).filter(m => by[m].length >= MIN_N);
    if (names.length < 2) return null;
    const morphs = names.map(key => {
      const name = key.replace(/^~/, ''), pts = by[key], stride = Math.max(1, pts.length / REF_MAX), ref = [];
      for (let s = 0; s < pts.length; s += stride) ref.push(feat(pts[Math.floor(s)]));
      let widths = WIDTHS;
      if (fuzzy[key]) {   // from the morph's own spread, per feature
        const fs = pts.map(feat);
        widths = WIDTHS.map((w, q) => {
          const mu = fs.reduce((a, f) => a + f[q], 0) / fs.length;
          return Math.max(w, FUZZ * Math.sqrt(fs.reduce((a, f) => a + (f[q] - mu) ** 2, 0) / fs.length));
        });
      }
      const m = q => pts.reduce((a, x) => a + x[q], 0) / pts.length;
      return { name, fuzzy: !!fuzzy[key], ref: Float64Array.from(ref.flat()), n: ref.length,
               inv: widths.map(w => 1 / w), logNorm: widths.reduce((a, w) => a + Math.log(w), 0),
               hue: hueOf(m(0), m(1), m(2)) };
    });
    return { morphs, len: long.len };
  }

  // Score of one sample's features z against a morph's cloud (distances in its kernel widths).
  function score(mp, z) {
    const r = mp.ref, w = mp.inv, D = F.length;
    let s = 0;
    for (let i = 0; i < r.length; i += D) {
      let d = 0;
      for (let q = 0; q < D; q++) { const e = (z[q] - r[i + q]) * w[q]; d += e * e; }
      s += Math.exp(-d / 2);
    }
    return Math.max(CLAMP, Math.log(s / mp.n + 1e-300));
  }

  // Classify a history (ring: get(k, j), brk(j), len): { name, fuzzy, label, p, hue, probs } with
  // name null if unknown (label: the name, "~" before a fuzzy one), or null if there are too few
  // samples yet.
  function classify(model, hist) {
    if (!model) return null;
    const zs = [];
    let since = 0;   // valid samples back to the last gap
    for (let j = hist.len - 1; j >= 0; j--) {
      const x = F.map(k => hist.get(k, j));
      if (usable(x)) { since++; if (zs.length < WIN) zs.push(feat(x)); }
      if (hist.brk(j)) break;
    }
    // the warm-up: drop the oldest samples if they are within WARM of the gap
    const keep = Math.min(zs.length, since - WARM);
    if (keep < MIN_N) return null;
    zs.length = keep;
    const L = model.morphs.map(mp => zs.reduce((a, z) => a + score(mp, z), 0) / keep);
    const N = L.map((v, i) => v - model.morphs[i].logNorm);   // normalised: comparable across widths
    const b = N.indexOf(Math.max(...N));
    const e = N.map(v => Math.exp(NEFF * (v - N[b]))), Z = e.reduce((a, v) => a + v, 0);
    const probs = e.map(v => v / Z), p = probs[b], mp = model.morphs[b];
    const known = p >= P_MIN && L[b] >= FIT_MIN;
    return { name: known ? mp.name : null, fuzzy: mp.fuzzy, label: known ? (mp.fuzzy ? '~' : '') + mp.name : null,
             p, fit: L[b], hue: mp.hue, best: mp.name,
             probs: Object.fromEntries(model.morphs.map((m, i) => [(m.fuzzy ? '~' : '') + m.name, probs[i]])) };
  }

  return { build, classify, hueOf };
})();
if (typeof module !== 'undefined') module.exports = MorphDetect;
