"use strict";
// ============================================================================================
//  The dots channel of the game ([dot:<n>|off] in script.txt), ported from Carlman: a separate
//  Lenia field of small stationary solitons (Carlman's dot: rule0, mu 0.38, sigma 0.07, at half
//  the kernel radius) that breathe under their own rule but never move. Any creature eats a dot
//  by covering it (its mass erases the dot's, see shaders/dotsim.glsl); a dot bitten deep enough
//  dissolves, a graze heals. The game keeps <n> dots on the board: an eaten one is replaced at
//  once, at a random spot away from the creatures and the other dots. The player's ship eating
//  one scores 10 (the score, shown at the top, starts from 0 when the dots are turned on).
// ============================================================================================

const Dots = (function () {
  const C = Core;
  // Carlman's dot soliton (models/solitons_direction.json rule0_mu0.3800_s0.0700_R18), 22 x 22 at R 18.
  const SOLITON = [
    [0,0,0,0,0,0,0,0.0015,0.0208,0.0451,0.0599,0.0574,0.0391,0.0142,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0.0068,0.0697,0.1637,0.2481,0.6059,0.6721,0.6601,0.5712,0.2241,0.1349,0.0475,0.0009,0,0,0,0,0],
    [0,0,0,0,0.0384,0.1698,0.64,0.9472,1,1,1,1,1,1,0.8626,0.5417,0.1276,0.0184,0,0,0,0],
    [0,0,0,0.0551,0.2293,0.8843,1,1,1,1,1,1,1,1,1,1,0.7632,0.178,0.0268,0,0,0],
    [0,0,0.0384,0.2293,0.9695,1,1,1,1,1,1,1,1,1,1,1,1,0.8295,0.1731,0.0143,0,0],
    [0,0.0068,0.1698,0.8843,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7354,0.1142,0,0],
    [0,0.0697,0.64,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.3662,0.0333,0],
    [0.0015,0.1637,0.9472,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7858,0.1061,0],
    [0.0208,0.2481,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.1836,0.0035],
    [0.0451,0.6059,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2408,0.0167],
    [0.0599,0.6721,1,1,1,1,1,1,1,1,0.9983,0.9991,1,1,1,1,1,1,1,1,0.5115,0.0258],
    [0.0574,0.6601,1,1,1,1,1,1,1,1,0.9991,0.9999,1,1,1,1,1,1,1,1,0.501,0.0245],
    [0.0391,0.5712,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2268,0.013],
    [0.0142,0.2241,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.9638,0.1633,0.0014],
    [0,0.1349,0.8626,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.7105,0.0841,0],
    [0,0.0475,0.5417,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.2158,0.0201,0],
    [0,0.0009,0.1276,0.7632,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.6274,0.0814,0,0],
    [0,0,0.0184,0.178,0.8295,1,1,1,1,1,1,1,1,1,1,1,1,0.7086,0.1279,0.0038,0,0],
    [0,0,0,0.0268,0.1731,0.7354,1,1,1,1,1,1,1,1,1,1,0.6274,0.1279,0.0104,0,0,0],
    [0,0,0,0,0.0143,0.1142,0.3662,0.7858,1,1,1,1,1,0.9638,0.7105,0.2158,0.0814,0.0038,0,0,0,0],
    [0,0,0,0,0,0,0.0333,0.1061,0.1836,0.2408,0.5115,0.501,0.2268,0.1633,0.0841,0.0201,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0.0035,0.0167,0.0258,0.0245,0.013,0.0014,0,0,0,0,0,0,0,0]
  ];
  const SOLITON_R = 18;
  const RULE = { mu: 0.38, sigma: 0.07, dt: 0.1 };   // Carlman's dt, per step
  const GONE = 0.4;      // a dot is eaten when its window holds less than this of its stamped mass
  const POINTS = 10;     // per dot the player's ship eats
  // Colour: a pale blue at first, its hue turning slowly towards cyan and on round the wheel.
  const HUE0 = 220, HUE_DEG_PER_S = -2, SAT = 0.6, LIGHT = 0.82;   // a full turn in 3 minutes
  // The field drifts as the star field (to the left by default; [stars] sets its angle and speed),
  // a bit slower than its near layer at speed 1 (draw.glsl starLayer 3.5 cells per second), with no
  // sway, by game time (core driftSec: the simulation's rate while playing, so the dots keep pace
  // with the ship at any speed and on any device; while paused, the speed-1 rate in real time),
  // while the stars go by real time. Sites are kept in the field's own frame; a site at (col, row) is at board
  // (col - offset[0], row - offset[1]).
  const DRIFT = 6;     // cells per game second (30 steps: 1 cell per 5 steps)
  let offset = [0, 0];

  let ready = false, dotR = 9, form = null, rad = 9;
  let on = false, target = 0, sites = [], score = 0, t0 = 0;
  const hooks = { score() {}, color(hue, rgb) {} };   // color: each frame while on, the dots' hue (degrees) and colour

  // The soliton resampled to `scale` of its size, bilinearly (Carlman's resizeSoliton).
  function resize(src, scale) {
    const h = src.length, w = src[0].length;
    const nh = Math.max(1, Math.round(h * scale)), nw = Math.max(1, Math.round(w * scale));
    const flat = new Float32Array(nh * nw);
    let mass = 0;
    for (let y = 0; y < nh; y++) {
      const fy = (y + 0.5) * h / nh - 0.5;
      const y0 = Math.max(0, Math.min(h - 1, Math.floor(fy))), y1 = Math.min(h - 1, y0 + 1), ty = Math.min(1, Math.max(0, fy - y0));
      for (let x = 0; x < nw; x++) {
        const fx = (x + 0.5) * w / nw - 0.5;
        const x0 = Math.max(0, Math.min(w - 1, Math.floor(fx))), x1 = Math.min(w - 1, x0 + 1), tx = Math.min(1, Math.max(0, fx - x0));
        const a = src[y0][x0] + (src[y0][x1] - src[y0][x0]) * tx, b = src[y1][x0] + (src[y1][x1] - src[y1][x0]) * tx;
        flat[y * nw + x] = a + (b - a) * ty;
        mass += flat[y * nw + x];
      }
    }
    return { flat, h: nh, w: nw, mass };
  }

  function setup() {
    if (ready) return;
    ready = true;
    dotR = Math.max(3, Math.round(C.R / 2));   // half the creatures' radius, as Carlman's dots
    form = resize(SOLITON, dotR / SOLITON_R);
    rad = Math.max(dotR, Math.ceil(Math.max(form.h, form.w) / 2) + 1);   // the window holds the whole stamp
    SimGL.setDotRule({ R: dotR, ...RULE });
  }

  // A random spot for a new dot: inside the field's frame (the stamp does not wrap), at least 3 dot radii
  // from the other dots and clear of every creature. Null if none was found this time.
  function findSpot() {
    const m = rad + 1, W = C.BW, H = C.BH, comps = C.found ? C.found.comps : [];
    for (let k = 0; k < 200; k++) {
      const col = m + Math.floor(Math.random() * (W - 2 * m)), row = m + Math.floor(Math.random() * (H - 2 * m));
      if (sites.some(s => Math.hypot(C.wrapD(s.col - col, W), C.wrapD(s.row - row, H)) < 3 * dotR)) continue;
      if (comps.some(c => Math.hypot(C.wrapD(c.col - (col - offset[0]), W), C.wrapD(c.row - (row - offset[1]), H)) < (c.radius || C.R) + C.R + rad)) continue;
      return { col, row };
    }
    return null;
  }
  function spawn() {
    const p = findSpot();
    if (!p) return false;
    SimGL.putDots(form.flat, form.h, form.w, p.row - Math.floor(form.h / 2), p.col - Math.floor(form.w / 2));
    sites.push({ col: p.col, row: p.row, mass: form.mass });
    return true;
  }
  // Wipe what is left of a dot (an eaten one's sliver, which might otherwise regrow).
  function wipe(s) {
    const e = 2 * rad + 1;
    SimGL.putDots(null, e, e, s.row - rad, s.col - rad);
  }

  // [dot:<n>]: keep n dots (the score starts from 0 if the dots were off); [dot:off]: none, no score.
  function set(n) {
    if (n === 'off') {
      on = false; target = 0; sites = [];
      SimGL.setDots(false);
      hooks.score();
      return;
    }
    setup();
    if (!on) { on = true; score = 0; sites = []; t0 = performance.now(); SimGL.setDots(true); }
    target = Math.min(n, SimGL.MAX_SITES);
    while (sites.length > target) wipe(sites.pop());
    hooks.score();
  }

  // Each frame: score and replace the eaten dots, top up to the target, turn the hue.
  function frame() {
    if (!on) return;
    const dt = C.driftSec, d = C.starDrift();
    offset = [(offset[0] - d.dir[0] * DRIFT * d.speed * dt) % C.BW, (offset[1] - d.dir[1] * DRIFT * d.speed * dt) % C.BH];
    SimGL.setDotOffset(...offset);
    const left = SimGL.readDotSites(sites.map(s => [s.col, s.row]), rad);
    const found = C.found;
    sites = sites.filter((s, i) => {
      if (left[i] >= GONE * s.mass) return true;
      // eaten: by the player's ship if it is the creature there (or the nearest one close by)
      const at = found ? Creatures.pick(found, s.col - Math.round(offset[0]), s.row - Math.round(offset[1]), C.R + rad, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
      if (at >= 0 && at === C.focusIdx) { score += POINTS; hooks.score(); }
      wipe(s);
      return false;
    });
    while (sites.length < target && spawn());
    const hue = (((HUE0 + HUE_DEG_PER_S * (performance.now() - t0) / 1000) % 360) + 360) % 360;
    const rgb = hsl(hue, SAT, LIGHT);
    SimGL.setDotColor(rgb);
    hooks.color(hue, rgb);
  }

  // HSL (h in degrees, s and l in 0..1) to RGB in 0..1.
  function hsl(h, s, l) {
    const a = s * Math.min(l, 1 - l);
    const f = n => { const k = (n + h / 30) % 12; return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)); };
    return [f(0), f(8), f(4)];
  }

  // Empty the field and put the dots in afresh (the score kept): the game's Reset.
  function clear() { if (!on) return; sites = []; SimGL.setDots(true); }

  // As if the player's ship ate a dot (?dev=1's P key): the score goes up by POINTS
  function eat() { score += POINTS; hooks.score(); }
  return { set, clear, frame, eat, hooks, hsl, get on() { return on; }, get score() { return score; }, get count() { return sites.length; } };
})();
