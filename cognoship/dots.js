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
  // Patterns that [dot:<name>:<how>] puts in the dots channel instead of the round dots, one at a
  // time and not replaced when eaten (a beacon). The channel runs one rule at a time: the dots' or
  // the pattern's. From Lenia/Python/animals.json (Lenia.py's numbering: kn 1 = exp core, gn 1 =
  // Gaussian growth, as dotsim.glsl's).
  //   Heptadentium scintillans (2D7s): checked in numpy with the dots' maths (exp core, rings
  //   1,5/6, Gaussian growth, dt 0.1): at R 18 it stays put (0.4 cells in 3,000 steps) and its mass
  //   holds within 1 %; at R 12 too; at R 9 (the round dots' R) it dies.
  const PATTERNS = {
    heptadentium: { name: 'Heptadentium scintillans', R: 18, betas: [1, 5 / 6], core: 'exp', mu: 0.22, sigma: 0.025, dt: 0.1,
      rle: '22.NpFK$11.wQwSuLrA6.pLtLwX2yOyGpB$9.tE3yOyIvOrN5.rPvMxT3yOyLuV$8.tH2yOtAqUtIwIwEtIpN2.MtJwGwXxV2yOyLuSpW$7.qCwNwWrA3.pNtPvTyFxHxN3yOyJ3yOyNuMpC$7.qWvHvArI5.pXtOxM9yOvXpC$7.pKtOvTvPvHtGrC4.pAtL8yOyMrK$7.AqS3yOyHwWtV4.XwG3yOyIvSsUyCyOwN$8.6yOyCsS3.BuP3yOxDrN2.wQyOrO$6.AvU7yOvJ4.tN2yOyLuM3.pEyOyH3.E$5.qQvGxRvLtVvJyI3yOvI4.rA2yOxDrB3.NxAyOtLSqHrQrAA$2.pVqDsLwWyOsBODSrSxFyOyN6.qTxRuD4.qBwV2yOyJuVwFwNrT$.uGtWuTxMyOuE5.2I13.XtUxP4yOyH2yO$pQxGxKyH2yON19.OtAxWyM3yOyJvR2yOpX$vD4yOxQqM19.OxD4yOyCtVpTwEyOxA$xU4yOxCvDtTsJqL17.wR3yOxStD2.sHyOxI$xJyOxLvFxRyCyK2yOyHxC17.xFxNwIsX2.DsVyNvB$uPyOvTuTxT6yOwH16.qIsAqV3.qTvIwHrU$rGwQvTvVyN7yO21.pNvTtFrXD$MtBuOvK4yOyEyHyArX20.pAxAuPW$.pAqRrR3yOuLrUqBM20.pHxDyOA$4.wCyOwHpH19.pPqKqPrQwEyOtG$4.pU2yO19.sRwUwNwAxJ2yO$4.MyFyO18.sJxJ6yOP$4.pRxTyO3.VrTrR12.qTuUyE5yOtCqHE$4.rVyDyO2.qXvCyIyOyApW11.rTvGxQyKyMyIxLvNuFqL$3.qFuLyOyApAqXvE4yOsH6.rIpK3.HrVuLvQvIvAtMvFwPsQ$3.tExAyOuDqLtFxJ3yOwVqT3.pGwNyOyAsUJ3.pBsRuMtUtLsWwWyOsJ$3.vMyJyOsXpGuFyB3yOvIU3.rQxR2yOvWqX3.qAxE2yOyA3yOK$4.2yOtQpIsFyA3yOuIA3.sDxS2yOyIuUqWpBqU7yOtP$4.vWyOxCrRsNvGuU2yOvWB3.rXxL8yOyKxJxMyG2yO$5.xKyOvOuEuGsEpF2yOqJ3.rUxE6yOyDtDqIqRsRvAwM$6.tJvIuCsDpJ2.yJyOqB2.tBxP4yOxVpO5.qV$7.PpKP4.vIxUsWtAyL3yOwIuBqOB$15.uFyF4yOxWvVuCqTB$15.rMwFyM4yOwLuAqH$16.sWwV4yOwFrXJ$18.sBvOvMtCpR!' },
  };
  const PATTERN_R_MIN = 12;   // below this, Heptadentium dies (numpy)
  // A pattern by its name: the key, or the full name or its first word, case and spaces ignored
  const norm = x => String(x).toLowerCase().replace(/\s+/g, '');
  function patternKey(name) {
    const n = norm(name);
    return Object.keys(PATTERNS).find(k => k === n || norm(PATTERNS[k].name) === n || norm(PATTERNS[k].name.split(' ')[0]) === n) ?? null;
  }
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

  let kind = null, dotR = 9, form = null, rad = 9;   // kind: 'dots' (round dots) or a PATTERNS key
  let on = false, target = 0, sites = [], score = 0, t0 = 0;
  // The scrolling world (the game's [clear_world:scroll]; null: none): { toBoard(wx, wy, margin):
  // the board cell [col, row] of a world point, or null if it is farther than half the board less
  // margin from the ship (on either axis); toWorld(col, row); extended(reach): a world point out of
  // the world view; ok(col, row, reach): a board spot away from the seam }. A site then also has
  // its world point (wx, wy), and is parked (wiped off the board, kept by its world point) while
  // too far from the ship; it comes back when the ship nears it.
  let world = null;
  let beaconsDrift = true;   // the beacons drift with the field as the round dots do (the game: not in challenge mode)
  // cells inside the board's edge less the dot's window (about 21 cells for Heptadentium) at which a
  // site comes back / goes: parked before the seam band (6 cells, cleared each frame) reaches it, as it
  // would wear it away ("eaten"); back as near the band as that allows, mostly under the game's black
  // stripe (2 R: a beacon, 42 cells across, can't hide wholly between it and the band)
  const PARK_IN = 10, PARK_OUT = 6;
  function setWorld(w) { world = w; }
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

  // The channel's rule and stamp: k = 'dots' (the round dots) or a PATTERNS key. A change of rule
  // empties the field (what is there would not live under the other rule).
  function setup(k) {
    if (kind === k) return;
    const changed = kind !== null;
    kind = k;
    if (k === 'dots') {
      dotR = Math.max(3, Math.round(C.R / 2));   // half the creatures' radius, as Carlman's dots
      form = resize(SOLITON, dotR / SOLITON_R);
      SimGL.setDotRule({ R: dotR, ...RULE });
    } else {
      const P = PATTERNS[k];
      if (!P.cells) P.cells = Rule.rle2cells(P.rle);
      dotR = Math.max(PATTERN_R_MIN, Math.round(P.R * C.R / 18));   // its own R at the game's R 18
      form = resize(P.cells, dotR / P.R);
      SimGL.setDotRule({ R: dotR, mu: P.mu, sigma: P.sigma, dt: P.dt, betas: P.betas, core: P.core });
    }
    rad = Math.max(dotR, Math.ceil(Math.max(form.h, form.w) / 2) + 1);   // the window holds the whole stamp
    if (changed && on) { sites = []; SimGL.setDots(true); }
  }

  // A spot for a new dot: inside the field's frame (the stamp does not wrap), at least 3 dot radii
  // from the other dots and clear of every creature. how: random (the default), far (of 30 such
  // spots, the one farthest from the player's ship), center (the board's centre, as it is now:
  // the field drifts). Null if none was found this time.
  function findSpot(how = 'random') {
    const m = rad + 1, W = C.BW, H = C.BH, comps = C.found ? C.found.comps : [];
    const wrapF = (v, n) => ((v % n) + n) % n;
    if (how === 'center') {
      const c = world ? world.center() : [W / 2, H / 2];   // (the scrolling world: the view's middle)
      const col = Math.round(wrapF(c[0] + offset[0], W)), row = Math.round(wrapF(c[1] + offset[1], H));
      return world ? { col, row } : { col: Math.min(W - m, Math.max(m, col)), row: Math.min(H - m, Math.max(m, row)) };   // (stamps wrap now; the torus's frame kept as it was)
    }
    const ship = how === 'far' ? C.ship() : null;
    let best = null, bestD = -1, found = 0;
    for (let k = 0; k < 400 && found < (ship ? 30 : 1); k++) {
      const col = m + Math.floor(Math.random() * (W - 2 * m)), row = m + Math.floor(Math.random() * (H - 2 * m));
      if (world && !world.ok(col - offset[0], row - offset[1], rad + PARK_IN)) continue;   // away from the scrolling world's seam
      if (sites.some(s => Math.hypot(C.wrapD(s.col - col, W), C.wrapD(s.row - row, H)) < 3 * dotR)) continue;
      if (comps.some(c => Math.hypot(C.wrapD(c.col - (col - offset[0]), W), C.wrapD(c.row - (row - offset[1]), H)) < (c.radius || C.R) + C.R + rad)) continue;
      found++;
      const d = ship ? Math.hypot(C.wrapD(ship.col - (col - offset[0]), W), C.wrapD(ship.row - (row - offset[1]), H)) : 0;
      if (d > bestD) { best = { col, row }; bestD = d; }
    }
    return best;
  }
  // (points: what the ship eating it scores; default POINTS)
  // how extended: out of the world view (the scrolling world: a parked site, put on the board as
  // the ship nears it; else as far)
  function spawn(how, points = POINTS, dist) {
    if (how === 'extended') {
      if (world && kind !== 'dots') { const [wx, wy] = world.extended(rad, dist); sites.push({ wx, wy, parked: true, mass: form.mass, points }); return true; }
      how = 'far';
    }
    let p = null;
    const at = how && how.place ? how.place(form.w, form.h) : Array.isArray(how) ? how : null;   // a placer (the scrolling world's edge), or a board spot
    if (at) p = { col: Math.round((((at[0] + offset[0]) % C.BW) + C.BW) % C.BW), row: Math.round((((at[1] + offset[1]) % C.BH) + C.BH) % C.BH) };
    if (!p) p = findSpot(how);
    if (!p) return false;
    stamp(p.col, p.row);
    const s = { col: p.col, row: p.row, mass: form.mass, points };
    if (world && kind !== 'dots') [s.wx, s.wy] = world.toWorld(p.col - offset[0], p.row - offset[1]);   // (a pattern, a beacon: kept by its world point; round dots are not: the seam takes them, others are put in)
    sites.push(s);
    return true;
  }
  // Write an h x w block (flat, or zeros if null) into the field at (row0, col0), wrapping round the
  // torus (putDots takes only blocks inside the board)
  function putWrapped(flat, h, w, row0, col0) {
    const W = C.BW, H = C.BH;
    row0 = ((row0 % H) + H) % H; col0 = ((col0 % W) + W) % W;
    for (const [y0, y1] of [[0, Math.min(h, H - row0)], [Math.min(h, H - row0), h]])
      for (const [x0, x1] of [[0, Math.min(w, W - col0)], [Math.min(w, W - col0), w]]) {
        if (y1 <= y0 || x1 <= x0) continue;
        let part = null;
        if (flat) { part = new Float32Array((y1 - y0) * (x1 - x0)); for (let y = y0; y < y1; y++) for (let x = x0; x < x1; x++) part[(y - y0) * (x1 - x0) + x - x0] = flat[y * w + x]; }
        SimGL.putDots(part, y1 - y0, x1 - x0, (row0 + y0) % H, (col0 + x0) % W);
      }
  }
  const stamp = (col, row) => putWrapped(form.flat, form.h, form.w, row - Math.floor(form.h / 2), col - Math.floor(form.w / 2));
  // Wipe what is left of a dot (an eaten one's sliver, which might otherwise regrow).
  function wipe(s) {
    const e = 2 * rad + 1;
    putWrapped(null, e, e, s.row - rad, s.col - rad);
  }
  // The scrolling world, each frame: a site too far from the ship is parked, a parked one near
  // enough put back on the board (its world point kept as it drifts with the field)
  function scrollSites() {
    if (!world) return;
    const W = C.BW, H = C.BH;
    for (const s of sites) {
      if (s.wx === undefined) {   // a round dot touching an edge: taken away, its replacement queued for the opposite edge
        const q = world.edgeOut(s.col - offset[0], s.row - offset[1], rad);
        if (q) { s.gone = true; queue.push(q); }
        continue;
      }
      if (s.parked) {
        const b = world.toBoard(s.wx, s.wy, rad + PARK_IN);
        if (!b) continue;
        s.col = Math.round((((b[0] + offset[0]) % W) + W) % W); s.row = Math.round((((b[1] + offset[1]) % H) + H) % H);
        stamp(s.col, s.row); s.parked = false;
      } else {
        [s.wx, s.wy] = world.toWorld(s.col - offset[0], s.row - offset[1]);
        if (!world.toBoard(s.wx, s.wy, rad + PARK_OUT)) { wipe(s); s.parked = true; }
      }
    }
  }

  // [dot:<n>]: keep n dots (the score starts from 0 if the dots were off); [dot:off]: none, no score.
  // how: where the first n are put (findSpot's; then the replacements too, but in a scrolling world).
  // The scrolling world: a dot touching an edge is taken away (no score) and its replacement put on
  // the opposite edge (world.edgeOut, queue), other replacements (an eaten dot's) on the edge the ship
  // heads to (world.edgeSpot)
  let refill = 'random', puts = 0, queue = [];   // puts: dots put in since [dot:<n>] (the first n by how)
  function set(n, how = 'random') {
    if (n === 'off') {
      on = false; target = 0; sites = [];
      SimGL.setDots(false);
      hooks.score();
      return;
    }
    turnOn();
    setup('dots');
    target = Math.min(n, SimGL.MAX_SITES);
    refill = how; puts = sites.length; queue = [];
    while (sites.length > target) wipe(sites.pop());
    hooks.score();
  }

  function turnOn() {
    if (on) return;
    on = true; score = 0; sites = []; t0 = performance.now(); SimGL.setDots(true);
  }
  // [dot:<pattern>[:<n>][:<how>]]: put one <pattern> (a PATTERNS name, e.g. Heptadentium) in the dots
  // channel at a spot found by `how` (random, far, center: findSpot), switching the channel to its
  // rule (the round dots go). Eaten by the ship, it scores `points` (default as a dot); it is not
  // replaced. False if the name is unknown or no spot was found.
  // dist: (extended, a scrolling world) its distance from the ship, in R
  function put(name, how = 'random', points = POINTS, dist) {
    const k = patternKey(name);
    if (!k) return false;
    turnOn();
    setup(k);
    target = 0;   // nothing topped up
    return spawn(how, points, dist);
  }

  // Each frame: score and replace the eaten dots, top up to the target, turn the hue.
  function frame() {
    if (!on) return;
    const dt = kind !== 'dots' && !beaconsDrift ? 0 : C.driftSec, d = C.starDrift();   // (beacons: stationary unless beaconsDrift)
    const cells = d.cells ?? DRIFT * d.speed;   // cells per game second ([stars:follow]: the ship's speed times its factor)
    offset = [(offset[0] - d.dir[0] * cells * dt) % C.BW, (offset[1] - d.dir[1] * cells * dt) % C.BH];
    SimGL.setDotOffset(...offset);
    scrollSites();
    limitMass();
    const live = sites.filter(s => !s.parked);
    const left = SimGL.readDotSites(live.map(s => [s.col, s.row]), rad);
    const found = C.found;
    // the ship's parts (a twin or trio ship: all of them)
    const shipParts = C.focusIdx >= 0 ? (C.objects().find(o => o.parts.includes(C.focusIdx))?.parts ?? [C.focusIdx]) : [];
    const lost = [];   // beacons gone, but not to the ship: put back nearby
    sites = sites.filter(s => {
      if (s.parked) return true;
      if (s.gone) { wipe(s); return false; }   // (not eaten: no score)
      if (left[live.indexOf(s)] >= GONE * s.mass) return true;
      // eaten: by the player's ship if it is the creature there (or the nearest one close by)
      const at = found ? Creatures.pick(found, s.col - Math.round(offset[0]), s.row - Math.round(offset[1]), C.R + rad, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
      const byShip = at >= 0 && shipParts.includes(at);
      if (byShip) { score += s.points ?? POINTS; hooks.score(); }
      wipe(s);
      if (!byShip && kind !== 'dots') lost.push(s);   // a beacon eaten by another creature (or worn away): not counted
      return false;
    });
    for (const s of lost) respawnNear(s);
    while (sites.length < target) {   // the scrolling world: the opposite edge for one taken away at an edge; once n are in, the edge ahead
      const at = world && puts >= target ? (queue.length ? world.queuedSpot(queue.shift(), rad) : world.edgeSpot(rad)) : null;
      if (!spawn(at || refill)) break;
      puts++;
    }
    const hue = (((HUE0 + HUE_DEG_PER_S * (performance.now() - t0) / 1000) % 360) + 360) % 360;
    const rgb = hsl(hue, SAT, LIGHT);
    SimGL.setDotColor(rgb);
    hooks.color(hue, rgb);
  }

  // A beacon lost (eaten by another creature, or worn away) put back near where it was: its centre
  // 2 to 4 R from there, at a random angle, clear of every creature (else as near as could be found)
  function respawnNear(s) {
    const comps = C.found ? C.found.comps : [], W = C.BW, H = C.BH;
    let best = null, bestGap = -Infinity;
    for (let k = 0; k < 60; k++) {
      const a = Math.random() * 2 * Math.PI, d = (2 + 2 * Math.random()) * C.R;
      const col = s.col + d * Math.cos(a), row = s.row + d * Math.sin(a);   // (the field's frame)
      if (world && !world.ok(col - offset[0], row - offset[1], rad + PARK_IN)) continue;   // (the scrolling world: away from its edges)
      const gap = Math.min(Infinity, ...comps.map(c => Math.hypot(C.wrapD(c.col - (col - offset[0]), W), C.wrapD(c.row - (row - offset[1]), H)) - (c.radius || C.R) - rad));
      if (gap > bestGap) { best = [col, row]; bestGap = gap; }
      if (gap > C.R) break;   // clear: take it
    }
    const [col, row] = best || [s.col, s.row];
    const t = { col: Math.round((((col % W) + W) % W)), row: Math.round((((row % H) + H) % H)), mass: form.mass, points: s.points };
    if (world) [t.wx, t.wy] = world.toWorld(t.col - offset[0], t.row - offset[1]);
    stamp(t.col, t.row);
    sites.push(t);
  }
  // Beacons left (pattern sites, on the board or not): [need:clear_beacon] waits for none
  const beaconCount = () => kind && kind !== 'dots' ? sites.length : 0;

  // The field's mass limit (as the bullets' [bullet_limit], with no tag of its own): round dots, n + 2
  // dots' worth ([dot:<n>]); beacons, 3 beacons' worth (each on the board). Over it (the field
  // grew, "exploded"), it is emptied and put back as it should be: round dots afresh (frame refills
  // them; no score), beacons stamped afresh at their sites (none lost, nor counted as eaten).
  function limitMass() {
    if (!form) return;
    const live = sites.filter(s => !s.parked);
    const limit = (kind === 'dots' ? target + 2 : 3 * Math.max(1, live.length)) * form.mass;
    const blk = SimGL.readDotBlocks();
    let total = 0;
    for (let i = 0; i < blk.length; i++) total += blk[i];
    if (total <= limit) return;
    console.info(`Cognoship: the dots channel grew to ${(total / form.mass).toFixed(1)} ${kind === 'dots' ? 'dots' : 'beacons'}' worth: put back as it was`);
    SimGL.setDots(true);   // (emptied)
    if (kind === 'dots') sites = sites.filter(s => s.parked);
    else for (const s of live) stamp(s.col, s.row);
  }

  // HSL (h in degrees, s and l in 0..1) to RGB in 0..1.
  function hsl(h, s, l) {
    const a = s * Math.min(l, 1 - l);
    const f = n => { const k = (n + h / 30) % 12; return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)); };
    return [f(0), f(8), f(4)];
  }

  // Empty the field and put the dots in afresh (the score kept): the game's Reset.
  function clear() { if (!on) return; sites = []; SimGL.setDots(true); }
  // The parked sites' world points (the scrolling world: the game points to them)
  const parked = () => sites.filter(s => s.parked).map(s => [s.wx, s.wy]);
  // The beacons on the board (pattern sites, not parked): [col, row, radius] (board cells), for the
  // game to put creatures clear of them
  const beaconCells = () => kind && kind !== 'dots' ? sites.filter(s => !s.parked).map(s => [s.col - offset[0], s.row - offset[1], rad]) : [];
  // A beacon's radius (its window), cells
  const beaconRadius = () => rad;
  // Every pattern site's world point (the scrolling world's beacons, on the board or not)
  const beacons = () => sites.filter(s => s.wx !== undefined).map(s => [s.wx, s.wy]);

  // As if the player's ship ate a dot (?dev=1's P key): the score goes up by POINTS
  function eat() { score += POINTS; hooks.score(); }
  return { set, put, patternKey, clear, frame, eat, hooks, hsl, setWorld, parked, beacons, beaconCells, beaconRadius, beaconCount,
           set beaconsDrift(v) { beaconsDrift = !!v; }, get on() { return on; }, get score() { return score; }, get count() { return sites.length; } };
})();
