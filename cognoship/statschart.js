"use strict";
// ============================================================================================
//  Stats of the focus (stats mode, key A): the two histories and their phase-space chart.
//
//  Time (t) is relative: the short history's from its oldest sample, each long-history piece's
//  from its first (also when loaded from a file).
//  Short history: a ring of the last SHORT_N samples (one per sim step), always recorded; a step
//  with no focus is an empty sample (NaN, not drawn), so the line runs out after the focus dies. Drawn as a line fading from the brightest colour (newest) to a darker one,
//  with a dot at the newest point.
//  Long history: unlimited, filled by Record (button, key S), 300 steps a press. Each sample also
//  carries the morph selected at the time. Drawn in a dark colour, each morph labelled once with
//  its name, placed where it overlaps no other label.
//  Both: a change of focus, or (long only) each Record, starts a new piece
//  with a gap in the line.
//
//  Each axis shows one stat. "Channel masses" is three series (R, G, B) and "total + channel
//  masses" four, so a chart can hold several trajectories, one per colour. Axes with several
//  series are paired by colour (R with R, ...); a single-series axis is shared by all of them.
//  3D is an orthographic view of the unit cube, rotated by dragging.
// ============================================================================================

const StatsChart = (function () {
  const KEYS = ['t', 'mass', 'mR', 'mG', 'mB', 'speed', 'omega'];
  const SHORT_N = 50;
  // name: on the axes; menu: in the axis pickers (adds the detail). Units, a spaceship's:
  // 1 m = R cells, 1 s = one unit time (T steps), 1 g = mass filling 1 m² at cell value 1.
  const STATS = {
    t:     { name: 'time',           unit: 's',     series: [['t', 'w']] },
    mass:  { name: 'total mass',     unit: 'g',     series: [['mass', 'w']] },
    chan:  { name: 'channel masses', unit: 'g',     series: [['mR', 'r'], ['mG', 'g'], ['mB', 'b']], menu: 'channel masses (R G B)' },
    both:  { name: 'total + channel masses', unit: 'g', series: [['mass', 'w'], ['mR', 'r'], ['mG', 'g'], ['mB', 'b']] },
    hue:   { name: 'hue',            unit: '°',     series: [['hue', 'w']] },
    speed: { name: 'speed',          unit: 'm/s',   series: [['speed', 'w']] },
    omega: { name: 'turn rate',      unit: 'rad/s', series: [['omega', 'w']] },
  };
  const OPTIONS = Object.keys(STATS);
  const C_OPTIONS = ['none', ...OPTIONS];   // the colour axis: also none (each trajectory its own colour)
  // per colour: [newest, oldest of the short history, long history]
  const COLORS = {
    w: [[255, 255, 255], [110, 110, 110], [70, 70, 70]],
    r: [[255, 70, 70], [140, 38, 38], [96, 30, 30]],
    g: [[80, 255, 110], [34, 140, 52], [26, 92, 38]],
    b: [[90, 150, 255], [40, 72, 150], [30, 48, 104]],
  };
  // Derived stats: computed from a history's samples when drawn, not recorded.
  // hue = the HSV hue of (R, G, B) masses, in degrees 0–360 (0 red, 120 green, 240 blue); empty
  // when the three are equal (grey). It wraps: just past red it jumps from ~360 to ~0.
  const DERIVED = {
    hue: (hist, j) => {
      const r = hist.get('mR', j), g = hist.get('mG', j), b = hist.get('mB', j);
      const mx = Math.max(r, g, b), d = mx - Math.min(r, g, b);
      if (!(d > 1e-9)) return NaN;
      const h = mx === r ? (g - b) / d : mx === g ? 2 + (b - r) / d : 4 + (r - g) / d;
      return (h * 60 + 360) % 360;
    },
  };
  const value = (hist, k, j) => DERIVED[k] ? DERIVED[k](hist, j) : hist.get(k, j);
  // ---- outliers ---------------------------------------------------------------------------
  // A sample is an outlier when it is far from its time neighbours (up to OUT_W on each side,
  // not across a gap) compared with how much they vary among themselves: |v − median| >
  // OUT_K × 1.4826 × MAD of the neighbours (≈ OUT_K standard deviations). A spike of one or two
  // steps (a detection glitch) is caught; a lasting jump is not, since after it the neighbours
  // on the far side agree with it. Outliers are left out of the axis ranges and of the lines (the
  // line joins their neighbours). The newest sample has neighbours on one side only, so a real
  // jump there counts as an outlier for a step or two, until the samples after it arrive.
  const OUT_W = 5, OUT_K = 8;   // tested: catches 1–2 step spikes, not a lasting jump; ~0.1% of noise
  const nbBuf = new Float64Array(2 * OUT_W), dvBuf = new Float64Array(2 * OUT_W);
  const medBuf = new Float64Array(2 * OUT_W);
  function median(a, n) {   // insertion sort: n <= 2 * OUT_W
    for (let i = 0; i < n; i++) {
      const v = a[i]; let j = i - 1;
      while (j >= 0 && medBuf[j] > v) { medBuf[j + 1] = medBuf[j]; j--; }
      medBuf[j + 1] = v;
    }
    return n % 2 ? medBuf[(n - 1) / 2] : (medBuf[n / 2 - 1] + medBuf[n / 2]) / 2;
  }
  function isOutlier(hist, k, j) {
    const v = value(hist, k, j);
    if (!isFinite(v)) return false;
    let n = 0;
    for (let d = 1; d <= OUT_W && j - d >= 0; d++) {   // back, stopping at a gap
      if (hist.brk(j - d + 1)) break;
      const u = value(hist, k, j - d); if (isFinite(u)) nbBuf[n++] = u;
    }
    for (let d = 1; d <= OUT_W && j + d < hist.len; d++) {   // forward
      if (hist.brk(j + d)) break;
      const u = value(hist, k, j + d); if (isFinite(u)) nbBuf[n++] = u;
    }
    if (n < 3) return false;
    const m = median(nbBuf, n);
    for (let q = 0; q < n; q++) dvBuf[q] = Math.abs(nbBuf[q] - m);
    return Math.abs(v - m) > OUT_K * 1.4826 * median(dvBuf, n) + 1e-9 * (Math.abs(m) + 1);
  }
  // Outlier flags of stat k over a whole history (Uint8Array, by j). For the long history they are
  // kept (hist.outl, cleared with it) and only the last OUT_W samples redone as it grows.
  function outliers(hist, k) {
    if (!hist.outl) {
      const f = new Uint8Array(hist.len);
      for (let j = 0; j < hist.len; j++) f[j] = isOutlier(hist, k, j);
      return f;
    }
    let c = hist.outl[k] ||= { n: 0, f: new Uint8Array(1024), fin: 0, lo: Infinity, hi: -Infinity };
    if (c.f.length < hist.len) { const g = new Uint8Array(2 * hist.len); g.set(c.f); c.f = g; }
    for (let j = Math.max(0, c.n - OUT_W); j < hist.len; j++) c.f[j] = isOutlier(hist, k, j);
    c.n = hist.len;
    return c.f;
  }
  // [min, max] of stat k over a history, outliers (flags f) left out. For the long history, the
  // samples whose flags are final (all but the last OUT_W) are folded into a kept min/max once.
  function rangeOf(hist, k, f) {
    let a = Infinity, b = -Infinity, j0 = 0;
    const c = hist.outl && hist.outl[k];
    if (c) {
      for (const end = Math.max(0, hist.len - OUT_W); c.fin < end; c.fin++) {
        if (f[c.fin]) continue;
        const v = value(hist, k, c.fin);
        if (v < c.lo) c.lo = v; if (v > c.hi) c.hi = v;
      }
      a = c.lo; b = c.hi; j0 = c.fin;
    }
    for (let j = j0; j < hist.len; j++) {
      if (f[j]) continue;
      const v = value(hist, k, j);
      if (v < a) a = v; if (v > b) b = v;
    }
    return [a, b];
  }
  // Axis ranges as shown: a range grows at once (nothing drawn falls outside) and shrinks by
  // RANGE_EASE of the way per frame (no jumps when an extreme point ages out).
  const RANGE_EASE = 0.08;
  let shown = [];   // per axis: { stat, dim, lo, hi }
  const OTHER_TAIL = 10;     // every-object mode: samples drawn of each object besides the focus
  const HIST_MAX = 8000;   // long-history points drawn per trajectory (strided beyond that)
  const rgb = c => `rgb(${c[0]},${c[1]},${c[2]})`;
  // CSS colour from hue (degrees; NaN: grey), saturation and value, 0..1
  function hsv(hDeg, sat, val) {
    if (!Number.isFinite(hDeg)) { const g = Math.round(val * 160); return `rgb(${g},${g},${g})`; }
    const f = n => { const k = (n + hDeg / 60) % 6; return Math.round(255 * val * (1 - sat * Math.max(0, Math.min(k, 4 - k, 1)))); };
    return `rgb(${f(5)},${f(3)},${f(1)})`;
  }
  const mix = (a, b, f) => a.map((v, i) => Math.round(v + (b[i] - v) * f));
  const fmt = v => v === 0 ? '0' : Math.abs(v) >= 1000 || Math.abs(v) < 0.01 ? v.toExponential(1) : v.toPrecision(3);

  // ---- histories ------------------------------------------------------------------------
  // Short: a ring of the last cap samples. get(k, j), j = 0 the oldest; brk(j): gap before j.
  function ring(cap = SHORT_N) {
    const data = Object.fromEntries(KEYS.map(k => [k, new Float32Array(cap)])), brks = new Uint8Array(cap);
    const at = j => (h.head - h.len + j + cap) % cap;
    // t is read relative to the oldest sample kept (tBase), so the window's time runs from 0
    let tBase = 0;
    const h = {
      cap, len: 0, head: 0,
      clear() { h.len = 0; h.head = 0; tBase = 0; },
      push(s, gap) {
        for (const k of KEYS) data[k][h.head] = s[k];
        brks[h.head] = gap ? 1 : 0;
        h.head = (h.head + 1) % cap; h.len = Math.min(h.len + 1, cap);
        tBase = 0;
        for (let j = 0; j < h.len; j++) { const v = data.t[at(j)]; if (Number.isFinite(v)) { tBase = v; break; } }
      },
      get: (k, j) => k === 't' ? data.t[at(j)] - tBase : data[k][at(j)],
      raw: (k, j) => data[k][at(j)],   // t not made relative (absolute sim time)
      brk: j => brks[at(j)] === 1,
    };
    return h;
  }

  // Long: grows without limit. segs: pieces of one focus and one morph, { start, morph, gap, fuzzy }
  // (gap: a new focus or a new recording, so the line is not joined to the previous piece; fuzzy:
  // set by hand in the file, for a loose morph matched by its own spread, see detect.js).
  function log() {
    let cap = 4096, data, brks, t0 = 0;
    const h = {
      len: 0, segs: [],
      clear() {
        cap = 4096; h.len = 0; h.segs = []; h.outl = {};
        data = Object.fromEntries(KEYS.map(k => [k, new Float32Array(cap)])); brks = new Uint8Array(cap);
      },
      // level: 'whole' for a twin or trio recorded as one (core.js wholes), else the default (a part)
      push(s, gap, morph, fuzzy = false, level = 'part') {
        if (h.len === cap) {   // grow x2
          cap *= 2;
          for (const k of KEYS) { const d = new Float32Array(cap); d.set(data[k]); data[k] = d; }
          const b = new Uint8Array(cap); b.set(brks); brks = b;
        }
        const j = h.len++;
        gap = gap || j === 0;
        if (gap) t0 = s.t;   // each piece's time starts at 0
        for (const k of KEYS) data[k][j] = k === 't' ? s.t - t0 : s[k];
        brks[j] = gap ? 1 : 0;
        const last = h.segs[h.segs.length - 1];
        if (gap || last.morph !== morph || last.fuzzy !== fuzzy || last.level !== level) h.segs.push({ start: j, morph, gap, fuzzy, level });
      },
      get: (k, j) => data[k][j],
      brk: j => brks[j] === 1,
      pieces: () => h.segs.filter(sg => sg.gap).length,
      // Remove the last piece (from the last gap to the end).
      dropLastPiece() {
        let i = h.segs.length - 1;
        if (i < 0) return;
        while (i > 0 && !h.segs[i].gap) i--;   // segs[0] always starts with a gap
        h.len = h.segs[i].start;
        h.segs.length = i;
        h.outl = {};   // outlier flags and kept ranges: redone on the next draw
      },
      // As a plain object for a file: each segment with its morph, gap and samples (rounded to 4
      // significant digits, in KEYS order). load() takes it back, appending.
      toJSON() {
        const r = v => Number.isFinite(v) ? +v.toPrecision(4) : null;
        return {
          keys: KEYS,
          segs: h.segs.map((sg, i) => {
            const end = i + 1 < h.segs.length ? h.segs[i + 1].start : h.len, rows = [];
            for (let j = sg.start; j < end; j++) rows.push(KEYS.map(k => r(data[k][j])));
            return { morph: sg.morph, gap: sg.gap, ...(sg.fuzzy ? { fuzzy: true } : {}), ...(sg.level === 'whole' ? { level: 'whole' } : {}), data: rows };
          }),
        };
      },
      load(obj) {
        const idx = KEYS.map(k => obj.keys.indexOf(k));
        for (const sg of obj.segs) sg.data.forEach((row, n) => {
          const s = {};
          KEYS.forEach((k, q) => { const v = idx[q] >= 0 ? row[idx[q]] : null; s[k] = v === null ? NaN : v; });
          h.push(s, n === 0 && sg.gap, sg.morph, !!sg.fuzzy, sg.level === 'whole' ? 'whole' : 'part');
        });
      },
    };
    h.clear();
    return h;
  }

  // ---- chart ----------------------------------------------------------------------------
  // One trajectory per colour: keys[i] = the stats field on axis i.
  function trajectories(axes) {
    const multi = axes.filter(a => STATS[a].series.length > 1);
    let colors = multi.length ? STATS[multi[0]].series.map(s => s[1]) : ['w'];
    for (const a of multi) colors = colors.filter(c => STATS[a].series.some(s => s[1] === c));
    return colors.map(c => ({
      color: c,
      keys: axes.map(a => { const s = STATS[a].series; return (s.length > 1 ? s.find(x => x[1] === c) : s[0])[0]; }),
    }));
  }

  // view: { dim: 2 | 3, axes: [x, y, z], c, lock, yaw, pitch } (lock: keep every axis's range). others: short histories of the other
  // objects (every-object mode), each drawn as a short tail (OTHER_TAIL samples; all of them when time is on
  // an axis), under the focus's line.
  // hover: one of them (under the mouse in the world view), drawn with its line as well.
  function draw(canvas, short, long, view, others = [], hover = null) {
    drawChart(canvas, short, long, view, others, hover);
    // the ship view's border (core.js drawOverlay), since the chart takes its place: the panel is
    // PANEL_W canvas px wide there, so its sizes are scaled to this canvas
    const w = canvas.clientWidth, h = canvas.clientHeight, k = w / PANEL_W;
    if (!w || !h) return;
    const ctx = canvas.getContext('2d');
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(160,156,172,.6)'; ctx.lineWidth = 2 * k;
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(3 * k, k, w - 4 * k, h - 2 * k, 10 * k);
    else ctx.rect(3 * k, k, w - 4 * k, h - 2 * k);
    ctx.stroke();
  }
  const PANEL_W = 384;   // core.js PW

  // Clickable axis labels of the last draw, in CSS px: { axis, x0, y0, x1, y1 } (see hitAxis).
  let hits = [];
  // A label's box from its text, anchor, alignment and baseline (rotated: the 2D y title).
  function addHit(ctx, axis, text, x, y, align, base, rotated = false) {
    const w = ctx.measureText(text).width, hh = 14;
    let x0 = align === 'left' ? x : align === 'right' ? x - w : x - w / 2;
    let y0 = base === 'top' ? y : y - hh / 2;
    if (rotated) hits.push({ axis, x0: x - hh / 2, y0: y - w / 2, x1: x + hh / 2, y1: y + w / 2 });
    else hits.push({ axis, x0: x0 - 3, y0: y0 - 2, x1: x0 + w + 3, y1: y0 + hh + 2 });
  }
  // The axis (0, 1, 2) whose label is at (x, y), CSS px in the canvas, or -1.
  const hitAxis = (x, y) => (hits.find(h => x >= h.x0 && x <= h.x1 && y >= h.y0 && y <= h.y1) || { axis: -1 }).axis;

  function drawChart(canvas, short, long, view, others, hover) {
    hits = [];
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth, h = canvas.clientHeight;
    if (!w || !h) return;
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = '#000'; ctx.fillRect(0, 0, w, h);
    ctx.font = '11px ui-monospace, Menlo, Consolas, monospace';
    ctx.lineJoin = 'round'; ctx.lineCap = 'round';

    const axes = view.axes.slice(0, view.dim);
    // colour axis (view.c): a stat that sets each sample's hue; paired by colour like the others,
    // kept apart as t.ckey
    const cStat = view.c && view.c !== 'none' ? view.c : null;
    const trajs = trajectories(cStat ? [...axes, cStat] : axes)
      .map(t => cStat ? { color: t.color, keys: t.keys.slice(0, -1), ckey: t.keys[t.keys.length - 1] } : t);

    // outlier flags, per history and stat, this frame
    const flagCache = new Map();
    const flags = (hist, k) => {
      let m = flagCache.get(hist);
      if (!m) flagCache.set(hist, m = {});
      return m[k] ||= outliers(hist, k);
    };
    const out = (hist, t, j) => t.keys.some(k => flags(hist, k)[j]);

    // ranges: per axis, over every series on it (R, G, B share a scale) and the histories drawn,
    // outliers left out; then eased (shown). i = 3: the colour axis (no padding).
    const axisRange = (i, keys, stat, pad = 0.05) => {
      let a = Infinity, b = -Infinity;
      for (const k of keys) {
        for (const hist of [short, long, hover]) {
          if (!hist || !hist.len) continue;
          const [p, q] = rangeOf(hist, k, flags(hist, k));
          a = Math.min(a, p); b = Math.max(b, q);
        }
        for (const o of others) {   // only what is drawn of them: their tails (all of it along t)
          if (!o.len) continue;
          const f = flags(o, k);
          for (let j = axes.includes('t') ? 0 : Math.max(0, o.len - OTHER_TAIL); j < o.len; j++) {
            if (f[j]) continue;
            const v = value(o, k, j);
            if (v < a) a = v; if (v > b) b = v;
          }
        }
      }
      if (!(a <= b)) { a = 0; b = 1; }
      if (b - a < 1e-9) { const p = Math.max(Math.abs(a) * 0.05, 1e-6); a -= p; b += p; }
      const p = (b - a) * pad;
      a -= p; b += p;
      const sh = shown[i];
      if (view.lock && !view.relock && sh && sh.stat === stat && sh.dim === view.dim) return [sh.lo, sh.hi];   // locked
      if (!sh || sh.stat !== stat || sh.dim !== view.dim) shown[i] = { stat, dim: view.dim, lo: a, hi: b };
      else {
        sh.lo = a < sh.lo ? a : sh.lo + (a - sh.lo) * RANGE_EASE;
        sh.hi = b > sh.hi ? b : sh.hi + (b - sh.hi) * RANGE_EASE;
      }
      return [shown[i].lo, shown[i].hi];
    };
    const lo = [], hi = [];
    for (let i = 0; i < axes.length; i++) {
      const [a, b] = axisRange(i, new Set(trajs.map(t => t.keys[i])), axes[i]);
      lo.push(a); hi.push(b);
    }
    let cLo = 0, cHi = 1;
    if (cStat && cStat !== 'hue') [cLo, cHi] = axisRange(3, new Set(trajs.map(t => t.ckey)), cStat, 0);
    // the hue of sample j: the hue stat itself, or the colour stat scaled max -> 0° (red),
    // min -> 240° (blue); NaN (grey) where the stat is undefined
    const hueAt = (hist, t, j) => {
      const v = value(hist, t.ckey, j);
      if (!Number.isFinite(v)) return NaN;
      if (cStat === 'hue') return v;
      return 240 * (1 - Math.min(1, Math.max(0, (v - cLo) / (cHi - cLo))));
    };
    const unit = (v, i) => (v - lo[i]) / (hi[i] - lo[i]);   // 0..1

    // projection: data (one value per axis) -> canvas px
    let project, clipped = false;
    if (view.dim === 2) {
      const L = 62, R = 14, T = 30, B = 38, pw = w - L - R, ph = h - T - B;
      project = v => [L + unit(v[0], 0) * pw, T + (1 - unit(v[1], 1)) * ph];
      drawAxes2D(ctx, axes, lo, hi, L, T, pw, ph);
      drawColorLegend(ctx, view.c, cLo, cHi, 8, 8);
      ctx.save(); ctx.beginPath(); ctx.rect(L, T, pw, ph); ctx.clip(); clipped = true;   // locked axes: points may fall outside
    } else {
      const s = Math.min(w, h) * 0.3, cx = w / 2, cy = h / 2 + 6;
      const cY = Math.cos(view.yaw), sY = Math.sin(view.yaw), cP = Math.cos(view.pitch), sP = Math.sin(view.pitch);
      const proj3 = (x, y, z) => {   // x, y, z in -1..1; z is up
        const X = x * cY - y * sY, D = x * sY + y * cY;
        return [cx + s * X, cy - s * (z * cP + D * sP)];
      };
      project = v => proj3(2 * unit(v[0], 0) - 1, 2 * unit(v[1], 1) - 1, 2 * unit(v[2], 2) - 1);
      drawAxes3D(ctx, axes, lo, hi, proj3);
      drawColorLegend(ctx, view.c, cLo, cHi, 8, 8 + 15 * 3);
    }
    const pt = (hist, t, j) => project(t.keys.map(k => value(hist, k, j)));

    if (!short.len && !long.len && !others.some(o => o.len)) {
      if (clipped) ctx.restore();
      return;   // nothing to draw yet: just the axes
    }

    // long history: one path per trajectory, a piece per segment (joined to the previous piece
    // unless it starts with a gap), strided to at most HIST_MAX points
    const labels = new Map();   // morph label -> candidate points along its line (all its pieces)
    if (long.len) {
      const stride = Math.max(1, Math.ceil(long.len / HIST_MAX));
      ctx.lineWidth = 1;
      for (const t of trajs) {
        // one path per colour: the trajectory's own, or with a colour axis hue bins of 15°, each
        // step coloured by its end sample
        const paths = new Map();
        const colAt = j => cStat ? hsv(Math.round(hueAt(long, t, j) / 15) * 15, 0.75, 0.5) : rgb(COLORS[t.color][2]);
        const line = (col, a, b) => {
          let p = paths.get(col);
          if (!p) paths.set(col, p = new Path2D());
          p.moveTo(a[0], a[1]); p.lineTo(b[0], b[1]);
        };
        long.segs.forEach((sg, i) => {
          const end = i + 1 < long.segs.length ? long.segs[i + 1].start - 1 : long.len - 1;
          const a = sg.gap ? sg.start : sg.start - 1;
          let prev = null;   // an empty point (a derived stat undefined there) lifts the pen
          const key = (sg.fuzzy ? '~' : '') + sg.morph;
          if (!labels.has(key)) labels.set(key, []);
          const cand = labels.get(key);
          for (let j = a; j <= end; j++) {
            if (j !== a && j !== end && (j - a) % stride) continue;
            if (out(long, t, j)) continue;   // an outlier: the line joins its neighbours
            const q = pt(long, t, j);
            if (!isFinite(q[0]) || !isFinite(q[1])) { prev = null; continue; }
            line(colAt(j), prev || [q[0] - 0.01, q[1]], q);   // the first point: a dot
            if ((j - a) % (10 * stride) === 0) cand.push(q);
            prev = q;
          }
        });
        for (const [col, p] of paths) { ctx.strokeStyle = col; ctx.stroke(p); }
      }
    }

    // short histories: each step its own colour, darkest (oldest) to brightest (newest); the
    // other objects first, the focus on top with a larger dot
    ctx.lineWidth = 1.5;
    // tail: how many of the newest samples to draw as a line (the rest is not drawn)
    const drawShort = (hist, dot, tail = Infinity) => {
      const n = hist.len;
      const isOut = (t, j) => out(hist, t, j);
      for (const t of trajs) {
        const [bright, dark] = COLORS[t.color];
        let prev = null;   // the last point drawn (outliers are skipped, the line joins their neighbours)
        for (let j = Math.max(0, n - tail); j < n; j++) {
          if (j > 0 && hist.brk(j)) prev = null;
          if (isOut(t, j)) continue;
          const b = pt(hist, t, j);
          if (!isFinite(b[0]) || !isFinite(b[1])) { prev = null; continue; }   // no object then
          if (prev) {
            const f = j / (n - 1);
            ctx.strokeStyle = cStat ? hsv(hueAt(hist, t, j), 0.85, 0.4 + 0.6 * f) : rgb(mix(dark, bright, f));
            ctx.beginPath(); ctx.moveTo(prev[0], prev[1]); ctx.lineTo(b[0], b[1]); ctx.stroke();
          }
          prev = b;
        }
      }
      if (n) for (const t of trajs) {
        if (isOut(t, n - 1)) continue;   // the newest is an outlier: no dot this step
        const [x, y] = pt(hist, t, n - 1);
        if (!isFinite(x) || !isFinite(y)) continue;
        ctx.fillStyle = cStat ? hsv(hueAt(hist, t, n - 1), 0.85, 1) : rgb(COLORS[t.color][0]);
        ctx.beginPath(); ctx.arc(x, y, dot, 0, 2 * Math.PI); ctx.fill();
      }
    };
    // the other objects: a short tail (their last OTHER_TAIL samples), or the whole line when time is on an
    // axis (a line along t is what shows their history there)
    const otherTail = axes.includes('t') ? Infinity : OTHER_TAIL;
    for (const o of others) drawShort(o, 3.5, otherTail);
    if (hover) drawShort(hover, 3.5);   // the object under the mouse in the world view: its line too
    drawShort(short, 3.5);

    if (clipped) ctx.restore();
    // morph labels of the long history: each name once, at a point on its line where its box
    // overlaps no label already placed (points tried from the middle of its run outwards, at four
    // corners around each point); if none is free, at the middle point anyway
    ctx.font = '10px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillStyle = '#9a9a9a';
    const placed = [], LH = 12;
    const free = r => r.x >= 0 && r.y >= 0 && r.x + r.w <= w && r.y + r.h <= h &&
      !placed.some(o => r.x < o.x + o.w && o.x < r.x + r.w && r.y < o.y + o.h && o.y < r.y + r.h);
    for (const [text, cand] of labels) {
      if (!cand.length) continue;
      const tw = ctx.measureText(text).width + 2, mid = cand.length >> 1;
      const order = [mid];
      for (let d = 1; order.length < Math.min(cand.length, 60); d++) {
        if (mid + d < cand.length) order.push(mid + d);
        if (mid - d >= 0) order.push(mid - d);
      }
      const boxes = q => [[3, -LH - 2], [-tw - 3, -LH - 2], [3, 2], [-tw - 3, 2]]
        .map(([dx, dy]) => ({ x: q[0] + dx, y: q[1] + dy, w: tw, h: LH }));
      let box = null;
      for (const i of order) { box = boxes(cand[i]).find(free); if (box) break; }
      box ||= boxes(cand[mid])[0];
      placed.push(box);
      ctx.fillText(text, box.x + 1, box.y);
    }
  }

  const axisTitle = a => STATS[a].unit ? `${STATS[a].name} (${STATS[a].unit})` : STATS[a].name;

  // The colour axis's legend at (x, y): "c: <stat>" (clickable, axis 3), and a colour bar with
  // its ends: blue (min) to red (max), or the hue circle for hue.
  function drawColorLegend(ctx, c, cLo, cHi, x, y) {
    ctx.font = '11px ui-monospace, Menlo, Consolas, monospace';
    ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillStyle = '#999';
    const label = c && c !== 'none' ? `c: ${axisTitle(c)}` : 'c: none';
    ctx.fillText(label, x, y);
    addHit(ctx, 3, label, x, y, 'left', 'top');
    if (!c || c === 'none') return;
    // laid out left to right: label, min, bar, max (each measured, so nothing overlaps)
    const [a, b] = c === 'hue' ? ['0', '360'] : [fmt(cLo), fmt(cHi)];
    const ax = x + ctx.measureText(label).width + 10;
    const bx = ax + ctx.measureText(a).width + 4, bw = 60, by = y + 4, bh = 4;
    const g = ctx.createLinearGradient(bx, 0, bx + bw, 0);
    for (let k = 0; k <= 6; k++) g.addColorStop(k / 6, hsv(c === 'hue' ? 60 * k : 240 * (1 - k / 6), 0.85, 1));
    ctx.fillStyle = g; ctx.fillRect(bx, by, bw, bh);
    ctx.fillStyle = '#999'; ctx.textAlign = 'left';
    ctx.fillText(a, ax, y);
    ctx.fillText(b, bx + bw + 4, y);
  }

  function drawAxes2D(ctx, axes, lo, hi, L, T, pw, ph) {
    ctx.strokeStyle = '#222'; ctx.lineWidth = 1;
    for (const f of [0, 0.5, 1]) {   // grid at min, middle, max
      ctx.beginPath(); ctx.moveTo(L + f * pw, T); ctx.lineTo(L + f * pw, T + ph); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(L, T + f * ph); ctx.lineTo(L + pw, T + f * ph); ctx.stroke();
    }
    ctx.strokeStyle = '#555';
    ctx.beginPath(); ctx.moveTo(L, T); ctx.lineTo(L, T + ph); ctx.lineTo(L + pw, T + ph); ctx.stroke();
    ctx.fillStyle = '#999';
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const f of [0, 0.5, 1]) ctx.fillText(fmt(lo[0] + f * (hi[0] - lo[0])), L + f * pw, T + ph + 4);
    ctx.fillText(axisTitle(axes[0]), L + pw / 2, T + ph + 20);
    addHit(ctx, 0, axisTitle(axes[0]), L + pw / 2, T + ph + 20, 'center', 'top');
    addHit(ctx, 1, axisTitle(axes[1]), 12, T + ph / 2, 'center', 'middle', true);
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const f of [0, 0.5, 1]) ctx.fillText(fmt(lo[1] + f * (hi[1] - lo[1])), L - 5, T + (1 - f) * ph);
    ctx.save(); ctx.translate(12, T + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center'; ctx.fillText(axisTitle(axes[1]), 0, 0); ctx.restore();
  }

  function drawAxes3D(ctx, axes, lo, hi, proj3) {
    const C = [-1, 1];
    ctx.lineWidth = 1; ctx.strokeStyle = '#262626';
    for (const a of C) for (const b of C) {   // the cube's 12 edges
      for (const e of [[[-1, a, b], [1, a, b]], [[a, -1, b], [a, 1, b]], [[a, b, -1], [a, b, 1]]]) {
        const p = proj3(...e[0]), q = proj3(...e[1]);
        ctx.beginPath(); ctx.moveTo(p[0], p[1]); ctx.lineTo(q[0], q[1]); ctx.stroke();
      }
    }
    // the three axes from the min corner, labelled at their far ends with the stat
    const o = proj3(-1, -1, -1), ends = [[1, -1, -1], [-1, 1, -1], [-1, -1, 1]], names = ['x', 'y', 'z'];
    ctx.strokeStyle = '#666'; ctx.fillStyle = '#aaa'; ctx.textBaseline = 'middle';
    ends.forEach((e, i) => {
      const p = proj3(...e);
      ctx.beginPath(); ctx.moveTo(o[0], o[1]); ctx.lineTo(p[0], p[1]); ctx.stroke();
      ctx.textAlign = p[0] >= o[0] ? 'left' : 'right';
      const dx = p[0] >= o[0] ? 5 : -5;
      const label = `${names[i]}: ${STATS[axes[i]].name}`;
      ctx.fillStyle = '#ccc'; ctx.fillText(label, p[0] + dx, p[1]);
      addHit(ctx, i, label, p[0] + dx, p[1], ctx.textAlign, 'middle');
    });
    // legend: each axis's stat and range
    ctx.textAlign = 'left'; ctx.textBaseline = 'top'; ctx.fillStyle = '#999';
    axes.forEach((a, i) => {
      const row = `${names[i]}: ${axisTitle(a)}  ${fmt(lo[i])} … ${fmt(hi[i])}`;
      ctx.fillText(row, 8, 8 + 15 * i);
      addHit(ctx, i, row, 8, 8 + 15 * i, 'left', 'top');
    });
    ctx.fillStyle = '#666';
    ctx.fillText('drag to rotate · click a label to change its stat', 8, 8 + 15 * 4);
  }

  // Forget the axis ranges: the next draw sets them afresh from the data.
  const resetRanges = () => { shown = []; };

  // The long-history file as text: header keys aligned, one sample per line, each segment's
  // columns aligned on the decimal point.
  function statsFileText(o) {
    const q = JSON.stringify;
    const head = ['format', 'species', 'R', 'dt', 'keys'], kw = Math.max(...head.map(k => q(k).length), q('segs').length);
    const key = (k, w) => q(k).padEnd(w) + ': ';
    const L = ['{'];
    for (const k of head) if (k in o) L.push('    ' + key(k, kw) + (Array.isArray(o[k]) ? '[' + o[k].map(v => q(v)).join(', ') + ']' : q(o[k])) + ',');
    L.push('    ' + key('segs', kw) + '[');
    o.segs.forEach((sg, si) => {
      const sk = Object.keys(sg).filter(k => k !== 'data'), sw = Math.max(...[...sk, 'data'].map(k => q(k).length));
      L.push('        {');
      for (const k of sk) L.push('            ' + key(k, sw) + q(sg[k]) + ',');
      // each value split at its decimal point (or exponent / end): integer part right-aligned,
      // the rest left-aligned, per column
      const cells = sg.data.map(r => r.map(v => { const t = q(v); const m = t.search(/[.e]/); return m < 0 ? [t, ''] : [t.slice(0, m), t.slice(m)]; }));
      const n = cells.length ? cells[0].length : 0;
      const iw = [], fw = [];
      for (let c = 0; c < n; c++) { iw[c] = Math.max(...cells.map(r => r[c][0].length)); fw[c] = Math.max(...cells.map(r => r[c][1].length)); }
      L.push('            ' + key('data', sw) + '[');
      cells.forEach((r, j) => L.push('                [' + r.map(([a, b], c) => a.padStart(iw[c]) + b.padEnd(fw[c])).join(', ') + ']' + (j < cells.length - 1 ? ',' : '')));
      L.push('            ]');
      L.push('        }' + (si < o.segs.length - 1 ? ',' : ''));
    });
    L.push('    ]', '}');
    return L.join('\n') + '\n';
  }

  return { draw, hitAxis, resetRanges, ring, log, statsFileText, OPTIONS, C_OPTIONS, STATS, SHORT_N };
})();
