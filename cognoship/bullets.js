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
//
//  World draw ([world_draw], [add_bullet:nebula|antimatter]): two more types after the patterns,
//  Nebula and Antimatter (FIELDS), which fire nothing: in use, the player paints the R G B fields
//  instead (field.js), which mask (Nebula) or erase (Antimatter) the creatures' channels; each is
//  in the book once added (fields), never in the list of patterns (nor the inventory).
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
  const num = (k, d) => { const v = parseFloat(C.q.get(k)); return Number.isFinite(v) && v >= 0 ? v : d; };
  const IDLE = 0.01;   // the field counts as empty below this many bullets' worth of mass
  const HIT = num('bhit', 0.2), WEAR = num('bwear', 0.03), SUB = Math.max(1, Math.round(num('bsub', 3)));
  const FIRE_STEPS = num('bcool', 0.7) * C.STEPS_PER_S;   // the least time between shots, in steps: ?bcool= game seconds
  let maxBullets = 8;      // more mass than this many bullets in the field: it is emptied ([bullet_limit:<x>])
  // colour by the channels they act on: white for all three, else that channel's colour (as
  // core.js CH_COLORS)
  const COLORS = { rgb: [1, 1, 1], r: [1, 0.231, 0.231], g: [0.231, 1, 0.353], b: [0.231, 0.482, 1] };

  // the world-draw types: how their field acts (glsim.js setField mode: 2 mask, 1 erase)
  const FIELDS = { 'Antimatter': { mode: 1 }, 'Nebula': { mode: 2 } }, FIELD_NAMES = Object.keys(FIELDS);
  const fields = new Set();   // the world-draw types in the book (none at first)
  const fieldKinds = () => FIELD_NAMES.filter(n => fields.has(n));
  const fieldOf = s => FIELD_NAMES.find(n => norm(n) === norm(s));   // a world-draw type's name, whatever its case
  const choices = () => list.concat(fieldKinds());   // what the book, select and cycle offer
  let on = false, kind = null, lastShot = -Infinity, mode = 2, chOpt = 'color_now', chNow = null, wearAll = false;
  let list = [];   // the patterns the player can fire ([bullet:...]), none at first; kind: the one in use
  const hooks = { change() {} };   // the gun turned on or off, the list changed, or another pattern picked
  let st = null;   // the pattern in use, as a stamp (BulletStamp)

  // [gun[:<option>...]], options in any order (separated by ':' or ','): on | off (the gun), erase |
  // mask (what a bullet does to a creature), rgb | r | g | b | color_now (which channels;
  // color_now: the player's channel pick, the R G B patches, followed live by the whole field),
  // wear_color | wear_all (the creature channels that wear a bullet down: the same ones it acts
  // on, or all three). Defaults: mask, color_now, wear_color. It fires the pattern in use, from
  // the list [bullet:...] keeps.
  // The options as given, over the defaults d ({ off, mode, ch, wear }); token: for the warnings
  function parseGun(options, d, token = 'gun') {
    const o = { ...d };
    for (const t0 of options.flatMap(a => a.split(',')).map(s => s.trim()).filter(Boolean)) {
      const t = t0.toLowerCase();
      if (t === 'off') o.off = true;
      else if (t === 'on') o.off = false;
      else if (t === 'erase' || t === 'mask') o.mode = t;
      else if (['rgb', 'r', 'g', 'b', 'color_now'].includes(t)) o.ch = t;
      else if (t === 'wear_color' || t === 'wear_all') o.wear = t;
      else console.warn(`script.txt: unknown option [${token}:...${t0}...], ignored`);
    }
    return o;
  }
  function setGun(options = []) {
    const o = parseGun(options, { off: false, mode: 'mask', ch: 'color_now', wear: 'wear_color' });
    if (o.off) { on = false; setKind(null); hooks.change(); return; }
    mode = o.mode === 'mask' ? 2 : 1; chOpt = o.ch; wearAll = o.wear === 'wear_all';
    on = true;
    setKind(choices().includes(kind) ? kind : list[list.length - 1] ?? null, true);   // (none chosen yet: the list's last)
    hooks.change();
  }
  // [bullet:<names>[:<equip>]], <names> in order, separated by ',': a pattern name adds it to the
  // list; random adds one not in it yet, at random; all adds every pattern; clear empties the
  // list (of patterns: Nebula and Antimatter stay); nebula / antimatter add that world-draw type to
  // the book. The gun switches to <equip> (added if not listed), else to the last pattern added; if
  // none was, the one in use stays if it is still listed, else it is the list's first.
  function setList(names = '', equip) {
    let last = null;
    const toks = names.split(',').map(s => s.trim()).filter(Boolean);
    const eq = equip?.trim() ? kindOf(equip.trim()) ?? fieldOf(equip.trim()) : null;   // the one to equip
    if (equip?.trim()) toks.push(equip.trim());   // added to the list too (warned if unknown)
    for (const t0 of toks) {
      const t = t0.toLowerCase(), add = n => { if (n && !list.includes(n)) { list.push(n); last = n; } };
      if (t === 'clear') list = [];
      else if (t === 'all') IDS.forEach(add);
      else if (t === 'random') { const rest = IDS.filter(n => !list.includes(n)); add(rest[Math.floor(Math.random() * rest.length)]); }
      else if (kindOf(t0)) add(kindOf(t0));
      else if (fieldOf(t0)) { fields.add(fieldOf(t0)); last = fieldOf(t0); }   // Nebula / Antimatter: into the book
      else console.warn(`script.txt: unknown bullet [bullet:...${t0}...], ignored`);
    }
    list.sort((a, b) => IDS.indexOf(a) - IDS.indexOf(b));   // always in the patterns' own order (the README's)
    if (eq) last = eq;   // equipped even if it was listed already
    if (on) setKind(last ?? (choices().includes(kind) ? kind : list[0] ?? null));
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
    if (!on || !id) { kind = on ? null : kind; SimGL.setBullets(null); SimGL.setField(null); return; }
    if (FIELDS[id]) {   // a world-draw type: no bullets, the field in their place
      kind = id; st = null; chNow = null;
      SimGL.setBullets(null);
      SimGL.setField({ mode: FIELDS[id].mode, hit: HIT });
      return;
    }
    SimGL.setField(null);
    if (!st) kind = null;   // (back from a field type: the pattern made afresh)
    if (id !== kind || force) {
      const fresh = id !== kind;
      kind = id;
      st = BulletStamp.make(KINDS[kind]);
      if (fresh) SimGL.setBullets(null);   // off then on (applyChannel): the field emptied
      SimGL.setBulletRule(st.rule);
    }
    chNow = null;
    applyChannel();
  }
  // The bullet book: switch to pattern `name`, if the gun has it (the field is emptied).
  function select(name) {
    if (!on || !choices().includes(name) || name === kind) return false;
    setKind(name);
    hooks.change();
    return true;
  }
  // D / F, or the ◀ ▶ buttons: the previous / next pattern in the list (the field is emptied).
  function cycle(step) {
    const all = choices();
    if (!on || all.length < 2) return false;
    const i = all.indexOf(kind);
    setKind(all[(i + step + all.length) % all.length]);
    hooks.change();
    return true;
  }

  // Could fire now: on, and the ship moving (its arrow shows), FIRE_STEPS after the last shot (steps: no cooldown passes while paused).
  function canFire() {
    return on && !!kind && !FIELDS[kind] && !!C.focus && C.simSteps - lastShot >= FIRE_STEPS;
  }
  // A or the Fire button: fire a bullet from the ship (a twin or trio: the whole) if canFire(). True if a bullet was put in.
  // From the front, the way the ship moves (its arrow); while it is (near) stationary, the arrow
  // hidden (core SPEED_MIN), at a random angle. ang (radians): that way instead (Autofire's aim).
  function fire(ang) {
    if (!canFire()) return false;
    const ship = C.ship(), k = ship.kin;   // a twin or trio ship: from its whole (its centre, heading and reach)
    lastShot = C.simSteps;
    const moving = k && k.speed >= C.SPEED_MIN;
    const ang0 = ang ?? (moving ? k.velAng : Math.random() * 2 * Math.PI);
    BulletStamp.fire(st, ang0, ship.col, ship.row, ship.reach || C.R, SimGL.putBullets);
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
    if (!on || !kind || FIELDS[kind]) return;
    applyChannel();
    if (SimGL.bulletsIdle) return;
    const blk = SimGL.readBulletBlocks();
    let total = 0;
    for (let i = 0; i < blk.length; i++) total += blk[i];
    const one = st.mass;   // one bullet's mass
    if (total > maxBullets * one || total < IDLE * one) SimGL.clearBullets();
  }

  // [world_draw]: names (world-draw type names) added to the book, quietly; null: both taken out
  // (one in use: back to the list's last pattern, or none). False if a name is unknown (nothing changed).
  function setWorldDraw(names) {
    if (names === null) fields.clear();
    else {
      const ks = names.map(fieldOf);
      if (ks.includes(undefined)) return false;
      ks.forEach(k => fields.add(k));
    }
    if (FIELDS[kind] && !fields.has(kind)) { if (on) setKind(list[list.length - 1] ?? null); else kind = list[list.length - 1] ?? null; }
    hooks.change();
    return true;
  }

  // [bullet_limit:<x>]: the field's mass limit, in bullets' worth (default 8)
  function setLimit(x = 8) { maxBullets = x; }
  // Empty the bullets channel (the game's Reset, [clear_world]).
  function clear() { SimGL.clearBullets(); }   // (also with the gun off: nothing left to come back)

  return { setGun, parseGun, setList, addRandom, kindOf, setLimit, clear, fire, canFire, cycle, select, frame, setWorldDraw, hooks, KINDS, HIT, WEAR, SUB,
           get maxBullets() { return maxBullets; }, get on() { return on; }, get kind() { return kind; },
           get list() { return list.slice(); },
           // the world-draw types offered now, and the one in use (null: a pattern, or none)
           get fieldKinds() { return fieldKinds(); }, get choices() { return choices(); }, get field() { return on && FIELDS[kind] ? kind : null; },
           get count() { return choices().length; },
           // (a field type: the colour of the field being painted, the player's channel pick)
           get color() { return on && FIELDS[kind] ? COLORS[['r', 'g', 'b'][C.channel]] || COLORS.rgb : COLORS[chNow] || COLORS.rgb; },
           // the gun's settings now, as the bullets field's options (the bots' bullets take them too: Hostile)
           get settings() {
             const c = chOpt === 'color_now' ? ['r', 'g', 'b'][C.channel] || 'rgb' : chOpt;
             const ch = c === 'rgb' ? [1, 1, 1] : ['r', 'g', 'b'].map(x => x === c ? 1 : 0);
             return { mode, ch, wearCh: wearAll ? [1, 1, 1] : ch, hit: HIT, wear: WEAR, sub: SUB };
           } };   // the bullets' colour now, rgb 0..1
})();

// --------------------------------------------------------------------------------------------
//  A bullet pattern as a stamp, and firing one from a creature: shared by the player's gun
//  (Bullets) and the bots' (Hostile).
// --------------------------------------------------------------------------------------------
const BulletStamp = (function () {
  const C = Core;
  const GAP = 8;     // cells between the firing creature's edge and a new bullet's back edge (as most patterns had when placed by their box)
  const EDGE = 0.1;  // a stamp cell at least this is part of the bullet (its back edge, fire)
  const betasOf = b => String(b).split(',').map(x => { const [n, d] = x.split('/'); return d ? n / d : +n; });

  // K: a pattern of bullet_patterns.js -> { K, rule (for SimGL.setBulletRule / setHostileRule),
  // size (the stamp's edge at any angle), mass (one bullet's, resampled), stamp(ang) }
  function make(K) {
    const rows = Rule.rle2cells(K.rle), w = Math.max(...rows.map(r => r.length));
    const src = rows.map(r => Array.from({ length: w }, (_, x) => r[x] || 0));
    const srcMass = src.flat().reduce((a, b) => a + b, 0);
    const bR = Math.max(4, Math.round(13 * C.R / 18));   // every bullet at R 13 (with the default creature R 18), resampled to it
    const scale = bR / K.R;
    const size = Math.ceil(Math.hypot(src.length, w) * scale) + 2;   // the stamp at any angle
    // The glider turned by `ang` (radians, board frame: col right, row down) and resampled to the
    // run radius, bilinearly, on a size x size grid centred on the source's centre.
    function stamp(ang) {
      const h = src.length, out = new Float32Array(size * size);
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
    return { K, size, mass: srcMass * scale * scale, stamp,
             rule: { R: bR, betas: betasOf(K.b), mu: K.mu, sigma: K.sigma, dt: 1 / K.T } };
  }

  // Write a size x size stamp with its top-left at (row0, col0), wrapping round the torus: in up
  // to four pieces, since put (SimGL.putBullets / putHostile) takes only blocks inside the board.
  function putWrapped(flat, size, row0, col0, put) {
    const W = C.BW, H = C.BH;
    row0 = ((row0 % H) + H) % H; col0 = ((col0 % W) + W) % W;
    const hs = [[0, Math.min(size, H - row0)], [Math.min(size, H - row0), size]];
    const ws = [[0, Math.min(size, W - col0)], [Math.min(size, W - col0), size]];
    for (const [y0, y1] of hs) for (const [x0, x1] of ws) {
      if (y1 <= y0 || x1 <= x0) continue;
      const h = y1 - y0, w = x1 - x0, part = new Float32Array(h * w);
      for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) part[y * w + x] = flat[(y0 + y) * size + x0 + x];
      put(part, h, w, (row0 + y0) % H, (col0 + x0) % W);
    }
  }

  // Fire st from a creature at (col, row) of radius `radius` (cells), flying at angle ang0
  // (radians, board frame), into the field `put` writes. Placed by the bullet's own mass, not its
  // stamp's box (a pattern stored off centre, or with a wide empty margin, would come out too
  // far, or to one side): its centre of mass on the line ahead of the creature (behind, for a
  // bomb: K.back, set off away from it), its edge nearest the creature (the cells >= EDGE
  // farthest back along the flight) GAP cells from the creature's edge. The same for every pattern.
  function fire(st, ang0, col0, row0, radius, put) {
    const K = st.K, size = st.size, dir0 = K.dir ? Math.atan2(K.dir[1], K.dir[0]) : 0;   // dir null: it stays put
    const ang = ang0 + (K.back ? Math.PI : 0);
    const flat = st.stamp(ang - dir0), ux = Math.cos(ang), uy = Math.sin(ang);
    let m = 0, mx = 0, my = 0;
    for (let y = 0; y < size; y++) for (let x = 0; x < size; x++) { const v = flat[y * size + x]; m += v; mx += v * x; my += v * y; }
    const cx = m > 0 ? mx / m : (size - 1) / 2, cy = m > 0 ? my / m : (size - 1) / 2;
    let back = 0;   // the back edge, along the flight, from the centre of mass (<= 0)
    for (let y = 0; y < size; y++) for (let x = 0; x < size; x++)
      if (flat[y * size + x] >= EDGE) back = Math.min(back, (x - cx) * ux + (y - cy) * uy);
    const d = radius + GAP - back;   // the centre of mass, ahead of the creature
    putWrapped(flat, size, Math.round(row0 + d * uy - cy), Math.round(col0 + d * ux - cx), put);
  }

  return { make, fire };
})();

// --------------------------------------------------------------------------------------------
//  The bots' bullets (challenge mode's [put_hostile]): a second bullets field (glsim.js
//  setHostile), white, acting as the player's bullets do ([gun]'s settings), the bots
//  included; its own pattern (one at a time: the last [put_hostile]'s), fired by the game
//  (game.js, the bots) through fire(). Emptied, as the player's, once it holds more than
//  maxBullets bullets' worth ([hostile_bullet_limit], default 8) or almost nothing.
// --------------------------------------------------------------------------------------------
const Hostile = (function () {
  let st = null, on = false, maxBullets = 8;   // maxBullets: [hostile_bullet_limit]
  // The bots' bullets act as the player's (the gun's settings, [gun]: mask or erase, the channels,
  // which channels wear them, followed live: the player's channel pick, color_now, too), unless
  // [hostile_gun] gives them their own (own: parseGun's { mode, ch, wear }; null: follow [gun]).
  // hold: [hostile_gun:off], the bots hold fire.
  let applied = '', own = null, hold = false;
  function settings() {
    if (!own) return Bullets.settings;
    const c = own.ch === 'color_now' ? ['r', 'g', 'b'][Core.channel] || 'rgb' : own.ch;
    const ch = c === 'rgb' ? [1, 1, 1] : ['r', 'g', 'b'].map(x => x === c ? 1 : 0);
    return { mode: own.mode === 'erase' ? 1 : 2, ch, wearCh: own.wear === 'wear_all' ? [1, 1, 1] : ch,
             hit: Bullets.HIT, wear: Bullets.WEAR, sub: Bullets.SUB };
  }
  function apply() {
    const o = settings();
    applied = JSON.stringify(o);
    SimGL.setHostile(o);
  }
  // [hostile_gun[:<options>]]: the bots' bullets' own settings, [gun]'s options (mask | erase, rgb |
  // r | g | b | color_now, wear_color | wear_all; on | off: the bots fire or hold fire), over the
  // defaults mask, color_now, wear_color; follow (or no options): as the player's [gun] again.
  function setGun(options = []) {
    const toks = options.flatMap(a => a.split(',')).map(s => s.trim().toLowerCase()).filter(Boolean);
    if (!toks.length || toks.includes('follow')) { own = null; hold = toks.includes('off'); }
    else {
      const o = Bullets.parseGun(toks.filter(t => t !== 'follow'), { off: false, mode: 'mask', ch: 'color_now', wear: 'wear_color' }, 'hostile_gun');
      hold = o.off;
      own = toks.some(t => t !== 'on' && t !== 'off') ? { mode: o.mode, ch: o.ch, wear: o.wear } : own;   // (on / off alone: the settings kept)
    }
    if (on) apply();
  }
  // Arm with pattern `name` (a bullet_patterns.js name): its rule, the field emptied. False if unknown.
  function set(name) {
    const k = Bullets.kindOf(name);
    if (!k) return false;
    if (!st || st.K !== Bullets.KINDS[k]) {
      st = BulletStamp.make(Bullets.KINDS[k]);
      SimGL.setHostile(null);
      SimGL.setHostileRule(st.rule);
    }
    apply();
    on = true;
    return true;
  }
  function off() { on = false; SimGL.setHostile(null); }
  function fire(col, row, radius, ang) { if (on && !hold) BulletStamp.fire(st, ang, col, row, radius, SimGL.putHostile); }
  function frame() {
    if (on && JSON.stringify(settings()) !== applied) apply();   // the gun changed (or the channel pick)
    if (!on || SimGL.hostileIdle) return;
    const blk = SimGL.readHostileBlocks();
    let total = 0;
    for (let i = 0; i < blk.length; i++) total += blk[i];
    if (total > maxBullets * st.mass || total < 0.01 * st.mass) SimGL.clearHostile();
  }
  function clear() { SimGL.clearHostile(); }
  function setLimit(x = 8) { maxBullets = x; }
  return { set, off, fire, frame, clear, setLimit, setGun, get on() { return on; }, get settings() { return settings(); } };
})();
