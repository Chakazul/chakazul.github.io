"use strict";
// ============================================================================================
//  The bullets channel of the game ([gun] and [bullet] in script.txt): a separate Lenia field of
//  patterns (bullet_patterns.js: Scutium valvatus by default; most are gliders) under the classic
//  Lenia rule, on the GPU
//  (glsim.js setBullets, shaders/bulletsim.glsl). Where a bullet overlaps a creature, it takes
//  mass out of the creature's channels (or masks them from the creature's own kernels), and the
//  creature's mass wears the bullet down, both gradually, unlike the dots, which are eaten at
//  once: so a bullet sinks some way into a body before it is gone.
//
//  The player fires (A, fire()): a bullet is put just in front of the ship, flying the way the
//  ship moves; a ship standing still (no arrow) fires at a random angle. A ship faster than the bullets runs
//  into its own: part of the game (a fast morph is a poor gunship). At most one shot every
//  FIRE_STEPS (0.7 game seconds, in steps); and if the field holds more than maxBullets (8, [bullet_limit]) bullets' worth of mass (bullets that
//  ran into each other and grew), it is emptied.
//
//  Tuning (URL): ?bhit= mass a bullet takes out of a creature per step per unit of bullet mass
//  (default 0.2), ?bwear= bullet mass worn per bullet step per unit of creature mass (0.03),
//  ?bsub= bullet steps per sim step (5: Scutium then flies at ~2.1 cells per step), ?bcool= the least
//  time between shots, in seconds (0.7).
// ============================================================================================

const Bullets = (function () {
  const C = Core;
  // The patterns, by name (bullet_patterns.js), in the order [bullet:all] lists them.
  const KINDS = BULLET_PATTERNS, IDS = Object.keys(KINDS);
  // A name in the script matches whatever its case and spacing ("Scutium fluens", "scutium fluens",
  // "scutiumfluens"); the ids used before the names still work.
  const norm = s => String(s).toLowerCase().replace(/[\s_-]+/g, '');
  const OLD = { scutium2: 'Scutium', scutium3: 'Scutium fluens', scutium4: 'Scutium vagus', scutium5: 'Scutium serratus',
                'scutium mk ii': 'Scutium fluens', 'scutium mk iii': 'Scutium vagus', 'scutium mk iv': 'Scutium serratus',
                orbium1: 'Orbium', paraptera1: 'Paraptera' };
  const BY_NAME = new Map(IDS.map(n => [norm(n), n]).concat(Object.entries(OLD).map(([o, n]) => [norm(o), n])));
  const kindOf = s => BY_NAME.get(norm(s));
  const betasOf = b => String(b).split(',').map(x => { const [n, d] = x.split('/'); return d ? n / d : +n; });
  const num = (k, d) => { const v = parseFloat(C.q.get(k)); return Number.isFinite(v) && v >= 0 ? v : d; };
  const IDLE = 0.01;   // the field counts as empty below this many bullets' worth of mass
  const HIT = num('bhit', 0.2), WEAR = num('bwear', 0.03), SUB = Math.max(1, Math.round(num('bsub', 3)));
  const FIRE_STEPS = num('bcool', 0.7) * C.STEPS_PER_S;   // the least time between shots, in steps: ?bcool= game seconds
  let maxBullets = 8;      // more mass than this many bullets in the field: it is emptied ([bullet_limit:<x>])
  const GAP = 8;            // cells between the ship's edge and a new bullet's back edge (as most patterns had when placed by their box)
  const EDGE = 0.1;         // a stamp cell at least this is part of the bullet (its back edge, fire)
  // colour by the channels they act on: white for all three, else that channel's colour (as
  // core.js CH_COLORS)
  const COLORS = { rgb: [1, 1, 1], r: [1, 0.231, 0.231], g: [0.231, 1, 0.353], b: [0.231, 0.482, 1] };

  let on = false, kind = null, lastShot = -Infinity, mode = 2, chOpt = 'color_now', chNow = null, wearAll = false;
  let list = [];   // the patterns the player can fire ([bullet:...]), none at first; kind: the one in use
  const hooks = { change() {} };   // the gun turned on or off, the list changed, or another pattern picked
  let src = null, srcMass = 0, bR = 13, scale = 1, size = 0;

  // [gun[:<option>...]], options in any order (separated by ':' or ','): on | off (the gun), erase |
  // mask (what a bullet does to a creature), rgb | r | g | b | color_now (which channels;
  // color_now: the player's channel pick, the R G B patches, followed live by the whole field),
  // wear_color | wear_all (the creature channels that wear a bullet down: the same ones it acts
  // on, or all three). Defaults: mask, color_now, wear_color. It fires the pattern in use, from
  // the list [bullet:...] keeps.
  function setGun(options = []) {
    const o = { off: false, mode: 'mask', ch: 'color_now', wear: 'wear_color' };
    for (const t0 of options.flatMap(a => a.split(',')).map(s => s.trim()).filter(Boolean)) {
      const t = t0.toLowerCase();
      if (t === 'off') o.off = true;
      else if (t === 'on') o.off = false;
      else if (t === 'erase' || t === 'mask') o.mode = t;
      else if (['rgb', 'r', 'g', 'b', 'color_now'].includes(t)) o.ch = t;
      else if (t === 'wear_color' || t === 'wear_all') o.wear = t;
      else console.warn(`script.txt: unknown option [gun:...${t0}...], ignored`);
    }
    if (o.off) { on = false; setKind(null); hooks.change(); return; }
    mode = o.mode === 'mask' ? 2 : 1; chOpt = o.ch; wearAll = o.wear === 'wear_all';
    on = true;
    setKind(list.includes(kind) ? kind : list[0] ?? null, true);
    hooks.change();
  }
  // [bullet:<names>[:<equip>]], <names> in order, separated by ',': a pattern name adds it to the
  // list; random adds one not in it yet, at random; all adds every pattern; clear empties the
  // list. The gun switches to <equip> (added if not listed), else to the last pattern added; if
  // none was, the one in use stays if it is still listed, else it is the list's first.
  function setList(names = '', equip) {
    let last = null;
    const toks = names.split(',').map(s => s.trim()).filter(Boolean);
    const eq = equip?.trim() ? kindOf(equip.trim()) : null;   // the one to equip
    if (equip?.trim()) toks.push(equip.trim());   // added to the list too (warned if unknown)
    for (const t0 of toks) {
      const t = t0.toLowerCase(), add = n => { if (n && !list.includes(n)) { list.push(n); last = n; } };
      if (t === 'clear') list = [];
      else if (t === 'all') IDS.forEach(add);
      else if (t === 'random') { const rest = IDS.filter(n => !list.includes(n)); add(rest[Math.floor(Math.random() * rest.length)]); }
      else if (kindOf(t0)) add(kindOf(t0));
      else console.warn(`script.txt: unknown bullet [bullet:...${t0}...], ignored`);
    }
    list.sort((a, b) => IDS.indexOf(a) - IDS.indexOf(b));   // always in the patterns' own order (the README's)
    if (eq) last = eq;   // equipped even if it was listed already
    if (on) setKind(last ?? (list.includes(kind) ? kind : list[0] ?? null));
    else if (last) kind = last;   // gun off: it starts with it when turned on
    hooks.change();
  }
  // Add pattern `pick` (if not listed yet), else one not listed yet at random, quietly: the one in
  // use stays (none in use: the new one). Returns its name, or null if none could be added (the
  // game's bullet book, creating one: [new_bullet:<n>])
  function addRandom(pick) {
    const rest = IDS.filter(n => !list.includes(n));
    if (!rest.length || (pick !== undefined && !rest.includes(pick))) return null;
    const name = pick ?? rest[Math.floor(Math.random() * rest.length)];
    list.push(name);
    list.sort((a, b) => IDS.indexOf(a) - IDS.indexOf(b));
    if (!kind) { if (on) setKind(list[0]); else kind = list[0]; }
    hooks.change();
    return name;
  }
  // Use pattern `id` (null: none, the field off): its rule, and an empty field (the old bullets
  // would run under the new rule). force: re-apply the gun's options even if it is the same one.
  function setKind(id, force = false) {
    if (!on || !id) { kind = on ? null : kind; SimGL.setBullets(null); return; }
    if (id !== kind || force) {
      const fresh = id !== kind;
      kind = id;
      const K = KINDS[kind];
      const rows = Rule.rle2cells(K.rle), w = Math.max(...rows.map(r => r.length));
      src = rows.map(r => Array.from({ length: w }, (_, x) => r[x] || 0));
      srcMass = src.flat().reduce((a, b) => a + b, 0);
      bR = Math.max(4, Math.round(13 * C.R / 18));   // every bullet at R 13 (with the default creature R 18), resampled to it
      scale = bR / K.R;
      size = Math.ceil(Math.hypot(src.length, w) * scale) + 2;   // the stamp at any angle
      if (fresh) SimGL.setBullets(null);   // off then on (applyChannel): the field emptied
      SimGL.setBulletRule({ R: bR, betas: betasOf(K.b), mu: K.mu, sigma: K.sigma, dt: 1 / K.T });
    }
    chNow = null;
    applyChannel();
  }
  // The bullet book: switch to pattern `name`, if the gun has it (the field is emptied).
  function select(name) {
    if (!on || !list.includes(name) || name === kind) return false;
    setKind(name);
    hooks.change();
    return true;
  }
  // D / F, or the ◀ ▶ buttons: the previous / next pattern in the list (the field is emptied).
  function cycle(step) {
    if (!on || list.length < 2) return false;
    const i = list.indexOf(kind);
    setKind(list[(i + step + list.length) % list.length]);
    hooks.change();
    return true;
  }

  // The glider turned by `ang` (radians, board frame: col right, row down) and resampled to the
  // run radius, bilinearly, on a size x size grid centred on the source's centre.
  function stamp(ang) {
    const h = src.length, w = src[0].length, out = new Float32Array(size * size);
    const c = Math.cos(ang), s = Math.sin(ang), m = (size - 1) / 2;
    const at = (y, x) => (y >= 0 && y < h && x >= 0 && x < w ? src[y][x] : 0);
    for (let y = 0; y < size; y++)
      for (let x = 0; x < size; x++) {
        const dx = (x - m) / scale, dy = (y - m) / scale;
        const sx = c * dx + s * dy + (w - 1) / 2, sy = -s * dx + c * dy + (h - 1) / 2;   // turned back
        const x0 = Math.floor(sx), y0 = Math.floor(sy), fx = sx - x0, fy = sy - y0;
        const v = (at(y0, x0) * (1 - fx) + at(y0, x0 + 1) * fx) * (1 - fy) + (at(y0 + 1, x0) * (1 - fx) + at(y0 + 1, x0 + 1) * fx) * fy;
        out[y * size + x] = Math.min(1, Math.max(0, v));
      }
    return out;
  }

  // Write a size x size stamp with its top-left at (row0, col0), wrapping round the torus: in up
  // to four pieces, since putBullets takes only blocks inside the board.
  function putWrapped(flat, row0, col0) {
    const W = C.BW, H = C.BH;
    row0 = ((row0 % H) + H) % H; col0 = ((col0 % W) + W) % W;
    const hs = [[0, Math.min(size, H - row0)], [Math.min(size, H - row0), size]];
    const ws = [[0, Math.min(size, W - col0)], [Math.min(size, W - col0), size]];
    for (const [y0, y1] of hs) for (const [x0, x1] of ws) {
      if (y1 <= y0 || x1 <= x0) continue;
      const h = y1 - y0, w = x1 - x0, part = new Float32Array(h * w);
      for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) part[y * w + x] = flat[(y0 + y) * size + x0 + x];
      SimGL.putBullets(part, h, w, (row0 + y0) % H, (col0 + x0) % W);
    }
  }

  // Could fire now: on, and the ship moving (its arrow shows), FIRE_STEPS after the last shot (steps: no cooldown passes while paused).
  function canFire() {
    return on && !!kind && !!C.focus && C.simSteps - lastShot >= FIRE_STEPS;
  }
  // A or the Fire button: fire a bullet from the ship if canFire(). True if a bullet was put in.
  // From the front, the way the ship moves (its arrow); while it is (near) stationary, the arrow
  // hidden (core SPEED_MIN), at a random angle.
  function fire() {
    if (!canFire()) return false;
    const ship = C.focus, k = C.kin;
    lastShot = C.simSteps;
    const moving = k && k.speed >= C.SPEED_MIN;
    const ang0 = moving ? k.velAng : Math.random() * 2 * Math.PI;
    const K = KINDS[kind], dir0 = K.dir ? Math.atan2(K.dir[1], K.dir[0]) : 0;   // dir null: it stays put
    const ang = ang0 + (K.back ? Math.PI : 0);   // a bomb (back): dropped behind the ship, set off away from it
    // Placed by the bullet's own mass, not its stamp's box (a pattern stored off centre, or with a
    // wide empty margin, would come out too far, or to one side): its centre of mass on the line
    // ahead of the ship (behind, for a bomb), its edge nearest the ship (the cells >= EDGE farthest
    // back along `ang`) GAP cells from the ship's edge. The same distance for every pattern.
    const flat = stamp(ang - dir0), ux = Math.cos(ang), uy = Math.sin(ang);
    let m = 0, mx = 0, my = 0;
    for (let y = 0; y < size; y++) for (let x = 0; x < size; x++) { const v = flat[y * size + x]; m += v; mx += v * x; my += v * y; }
    const cx = m > 0 ? mx / m : (size - 1) / 2, cy = m > 0 ? my / m : (size - 1) / 2;
    let back = 0;   // the back edge, along the flight, from the centre of mass (<= 0)
    for (let y = 0; y < size; y++) for (let x = 0; x < size; x++)
      if (flat[y * size + x] >= EDGE) back = Math.min(back, (x - cx) * ux + (y - cy) * uy);
    const d = (ship.radius || C.R) + GAP - back;   // the centre of mass, ahead of the ship
    const col = ship.col + d * ux, row = ship.row + d * uy;
    putWrapped(flat, Math.round(row - cy), Math.round(col - cx));
    return true;
  }

  // Each frame: follow the channel pick (color_now); empty the field if it holds too much mass.
  // The channels the bullets act on, and their colour; for color_now, the player's pick now. A change
  // keeps the field (bullets in flight switch too).
  function applyChannel() {
    const c = chOpt === 'color_now' ? ['r', 'g', 'b'][C.channel] || 'rgb' : chOpt;
    if (c === chNow) return;
    chNow = c;
    const ch = c === 'rgb' ? [1, 1, 1] : ['r', 'g', 'b'].map(x => x === c ? 1 : 0);
    SimGL.setBullets({ mode, ch, wearCh: wearAll ? [1, 1, 1] : ch, hit: HIT, wear: WEAR, sub: SUB });
    SimGL.setBulletColor(COLORS[c]);
  }

  // Once the field has (almost) emptied, the last of it is wiped and the field goes idle (glsim:
  // no bullet steps, no bullet path in the creatures' step, no readback) until the next shot.
  function frame() {
    if (!on || !kind) return;
    applyChannel();
    if (SimGL.bulletsIdle) return;
    const blk = SimGL.readBulletBlocks();
    let total = 0;
    for (let i = 0; i < blk.length; i++) total += blk[i];
    const one = srcMass * scale * scale;   // one bullet's mass
    if (total > maxBullets * one || total < IDLE * one) SimGL.clearBullets();
  }

  // [bullet_limit:<x>]: the field's mass limit, in bullets' worth (default 8)
  function setLimit(x = 8) { maxBullets = x; }
  // Empty the bullets channel (the game's Reset, [clear_world]).
  function clear() { SimGL.clearBullets(); }   // (also with the gun off: nothing left to come back)

  return { setGun, setList, addRandom, kindOf, setLimit, clear, fire, canFire, cycle, select, frame, hooks, KINDS, get on() { return on; }, get kind() { return kind; },
           get list() { return list.slice(); },
           get color() { return COLORS[chNow] || COLORS.rgb; } };   // the bullets' colour now, rgb 0..1
})();
