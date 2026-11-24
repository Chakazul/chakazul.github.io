"use strict";
// ============================================================================================
//  Morph models (models/<code>.json, written by tools/models.py from headless experiments
//  with the game's rule, mask and detection): which action (a mask held at a spot in the ship's
//  frame for a hold) gives which outcome, per morph. Shared by the lab's Model mode and, later,
//  the game's bots and autopilot, so both use the same condition -> action -> outcome rules.
//
//  The condition: the ship's named morph (the name box; ??? = no actions) and its turning
//  direction (the sign of its turn rate). The book's morphs all turn + (clockwise; mirrored if
//  their pattern does not): a ship turning - uses the mirror image (spots' sides swapped, turns
//  negated). Outcomes are said relative to the ship: a turn right / left as seen (clockwise /
//  anticlockwise), keeping its turning direction or reversing it; into another morph, keeping
//  or reversing; dissolve.
// ============================================================================================

const Models = (function () {
  const OMEGA_SIDE = 0.01;   // rad per unit time: a slower turn has no turning direction (as the book's)
  const P_MIN = 0.75;        // an action counts for its outcome if it gives it in 3 of 4 trials
  const P_REP = 0.5;         // a replication (rare and useful) is offered from half the trials
  const TURN_MIN = 5;        // degrees: a smaller heading change is no turn
  const TURN_BIG = 45;       // degrees: from here a turn is big
  const books = {};          // code -> book, null (none), or undefined (loading)

  // Load models/<code>.json once (then: get(code))
  function load(code, then) {
    if (code in books) { if (then && books[code] !== undefined) then(books[code]); return; }
    books[code] = undefined;
    fetch(`models/${encodeURIComponent(code)}.json`, { cache: 'no-cache' })
      .then(r => r.ok ? r.json() : null).catch(() => null)
      .then(b => { if (b) robustP(b); books[code] = b; if (then) then(b); });
  }
  const get = code => books[code];
  // An action re-tested (tools/models.py verify: from fresh starts, its mask a cell off, the heading
  // 15 degrees off, its hold a few steps off, as the lab and the game place it; a replication only
  // if the copies last): its p the lower of the two, so the menus, the graph and the bots prefer
  // robust actions; p_sweep the sweep's own
  function robustP(b) {
    for (const m of Object.values(b.morphs || {})) for (const e of m.actions || [])
      if (e.p_robust !== undefined && e.p_robust < e.p) { e.p_sweep = e.p; e.p = e.p_robust; }
  }

  // The ship's turning direction from its turn rate, kept while it turns too slowly to tell (last)
  const sign = (omega, last = 1) => Math.abs(omega) >= OMEGA_SIDE ? Math.sign(omega) : last;

  // One action seen from a ship: m the book's morph, sg the ship's turning direction (ignored for
  // a morph with none). Returns { e (the book's row), fs: [fwd, side] cells, turn (degrees as
  // seen, + clockwise; null if not measured), spot (its name as the ship sees it: left and right
  // swapped if mirrored), keeps (its turning direction kept: true; false:
  // reversed), into (the morph it becomes, or null: stays; 'dissolved'), split (the morphs of the
  // bodies it splits into, ['limus', 'limus'], '???' for one not named; or null) }
  function seen(book, m, e, sg) {
    const mir = mirrored(m, sg) ? -1 : 1, R = book.R;
    const x = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(e.effect);
    const spot = mir < 0 ? e.spot.replace(/left|right/g, w => w === 'left' ? 'right' : 'left') : e.spot;
    return { e, spot, fs: [e.off[0] * R, e.off[1] * R * mir], turn: e.turn !== undefined ? e.turn * mir : null,
             keeps: !(e.effect === 'reverse' || (x && x[2] === 'reversed')),
             into: e.effect === 'dissolve' ? 'dissolved' : x ? x[1] : ['keep', 'reverse'].includes(e.effect) || replicateOf(e.effect) ? null : e.effect,
             split: splitOf(e.effect), replicate: replicateOf(e.effect), copies: e.copies || null };
  }
  // Seen in mirror image: a mover turning the other way (sg < 0), or a still spinner (spin) spinning
  // the other way from the model's (sg < 0: the lab's spin sign, observed spin x the model's)
  function mirrored(m, sg) { return !!(m && (m.chirality || m.spin) && sg < 0); }
  // 'split: gyrans + limus' -> ['gyrans', 'limus'], else null
  const splitOf = ef => ef.startsWith('split: ') ? ef.slice(7).split(' + ') : null;
  // 'replicate: tachus' (copies of one morph, the mass grown: e.g. PJ1EF7 rectus through vagus) -> 'tachus', else null
  const replicateOf = ef => ef.startsWith('replicate: ') ? ef.slice(11) : null;
  // a twin or triple (several parts detected as one whole), by the book's wholes
  const isWhole = (book, name) => !!(book && (book.wholes || []).includes(name));
  // A split's mark: the parts' initials, e.g. LL, GL (? for one not named)
  const splitMark = parts => parts.map(n => /^[a-z]/i.test(n) ? n[0].toUpperCase() : '?').join('');

  // The outcome a seen action gives, as a key (or null: none worth offering): turn-<small|big>-
  // <right|left>[-rev], into-<morph>[-rev], split-<morph>+<morph>, dissolve
  function outcomeOf(a) {
    if (a.e.p < (a.replicate ? P_REP : P_MIN)) return null;   // (a replication, rare and useful: offered from half the trials)
    if (a.into === 'dissolved') return 'dissolve';
    if (a.replicate) return `replicate-${a.replicate}`;
    if (a.split) return a.split.includes('???') ? null : `split-${a.split.join('+')}`;
    if (a.into) return /^[a-z]+$/i.test(a.into) && a.into !== 'unsettled' ? `into-${a.into}${a.keeps ? '' : '-rev'}` : null;   // (unsettled: ??? at the end, no morph)
    if (a.turn === null || Math.abs(a.turn) < TURN_MIN) return a.keeps ? null : 'reverse';
    return `turn-${Math.abs(a.turn) >= TURN_BIG ? 'big' : 'small'}-${a.turn > 0 ? 'right' : 'left'}${a.keeps ? '' : '-rev'}`;
  }

  // Words for an outcome key (base: the ship's morph)
  function words(key, base) {
    if (key === 'dissolve') return 'Dissolve';
    if (key === 'reverse') return `Turn the other way (no heading change)`;
    if (key.startsWith('split-')) return `Split into ${key.slice(6).split('+').join(' + ')}`;
    if (key.startsWith('replicate-')) return `Replicate into several ${key.slice(10)}`;
    const rev = key.endsWith('-rev'), k = rev ? key.slice(0, -4) : key;
    const t = /^turn-(small|big)-(right|left)$/.exec(k);
    if (t) return `${t[1] === 'small' ? 'Small' : 'Big'} turn ${t[2]}${rev ? ', then turning the other way' : ''}`;
    const i = /^into-(.+)$/.exec(k);
    return `Into ${i[1]}${rev ? ', turning the other way' : ''}`;
  }
  const ORDER = k => k.startsWith('turn-') ? (k.endsWith('-rev') ? 1 : 0) : k === 'reverse' ? 1 : k.startsWith('into-') ? 2 : k.startsWith('replicate-') ? 3 : k.startsWith('split-') ? 4 : 5;
  const GROUPS = ['Turn (keeps its turning direction)', 'Turn and turn the other way after', 'Into another morph', 'Replicate (multiply)', 'Split into two', 'Dissolve'];

  // The best action among those giving one outcome: the most reliable (p, in steps of 0.05), the
  // name box kept (always named), the shortest hold, the shortest time unnamed; for a turn, then
  // the one nearest the middle of its size (small 20 degrees, big 90)
  function rank(a, key) {
    const mid = key.includes('-big-') ? 90 : 20;
    return [-Math.round(a.e.p * 20), a.e.always_named ? 0 : 1, a.e.hold, a.e.unnamed_max || 0,
            a.turn !== null ? Math.abs(Math.abs(a.turn) - mid) : 0];
  }
  const less = (x, y) => { for (let i = 0; i < x.length; i++) if (x[i] !== y[i]) return x[i] < y[i]; return false; };

  // The outcomes the ship can reach now, each with its best action: [{ key, label, group, best,
  // n (actions giving it) }], in menu order. book: the species' book; name: the ship's named morph
  // (null: ??? -> none); sg: its turning direction (sign()). A morph not in the book: none.
  function maneuvers(book, name, sg, mm = null, everyPhase = false) {   // (mm: the morph's entry as now, e.g. phased(); everyPhase: allActions)
    const m = mm || (book && name ? book.morphs[name] : null);
    if (!m) return [];
    const by = {};
    for (const e of everyPhase ? allActions(m) : m.actions) {
      const a = seen(book, m, e, sg), key = outcomeOf(a);
      if (!key) continue;
      const r = rank(a, key), o = by[key] || (by[key] = { key, n: 0, best: null, r: null });
      o.n++;
      if (!o.best || less(r, o.r)) { o.best = a; o.r = r; }
    }
    return Object.values(by).map(o => ({ key: o.key, label: words(o.key, m.base), group: GROUPS[ORDER(o.key)], best: o.best, n: o.n }))
      .sort((x, y) => ORDER(x.key) - ORDER(y.key) || x.key.localeCompare(y.key));
  }

  // The best action for one outcome now, or null if none
  const choose = (book, name, sg, key, mm = null, everyPhase = false) => (maneuvers(book, name, sg, mm, everyPhase).find(o => o.key === key) || {}).best || null;

  // An oscillator swept over its swing (the book's phase_bins, swing_period): its phase bin now,
  // as tools/models.py swing_phases: a sine fitted to its recent turn rates (hist: [[t, omega]],
  // t in time units, omega core kin.omega) over the last swing; omega ~ sin(phase). dt: time units
  // a step. null: not swept so, or not yet a swing's worth of history
  function swingBin(m, hist, dt) {
    if (m && m.phase_kind === 'state') return stateBin(m, hist, dt);
    if (m && m.phase_kind === 'age') return ageBin(m, hist, dt);
    if (!m || !m.phase_bins || !hist || hist.length < 8) return null;
    const sig = m.phase_signal === 'speed' ? 2 : 1;   // (hist: [t, omega, speed / R]: the turn rate's cycle, or the speed's)
    const P = m.swing_period * dt, w0 = 2 * Math.PI / P, tn = hist[hist.length - 1][0];
    const h = hist.filter(x => x[0] > tn - P);
    if (tn - h[0][0] < 0.8 * P) return null;
    // least squares: omega ~ a sin(w0 t) + b cos(w0 t) + c (normal equations, 3 x 3)
    const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], v = [0, 0, 0];
    for (const x of h) {
      const t = x[0], om = x[sig];
      const f = [Math.sin(w0 * t), Math.cos(w0 * t), 1];
      for (let r = 0; r < 3; r++) { v[r] += f[r] * om; for (let c = 0; c < 3; c++) M[r][c] += f[r] * f[c]; }
    }
    const det = A => A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    const D = det(M);
    if (Math.abs(D) < 1e-12) return null;
    const col = k => M.map((row, r) => row.map((x, c) => c === k ? v[r] : x));
    const a = det(col(0)) / D, b = det(col(1)) / D;
    let ph = (w0 * tn + Math.atan2(b, a)) % (2 * Math.PI);
    if (ph < 0) ph += 2 * Math.PI;
    return Math.floor(ph / (2 * Math.PI) * m.phase_bins) % m.phase_bins;
  }
  // A short-lived morph swept by age (VT049W gyrans: age_settle + offsets over age_span steps):
  // its age now, from when it was first named so (hist[0]) plus the name box's delay (NAMED_AT steps
  // of history it needs), in bins; null once older than the span (or with no history)
  const NAMED_AT = 50;
  function ageBin(m, hist, dt) {
    if (!m || !hist || !hist.length) return null;
    const age = (hist[hist.length - 1][0] - hist[0][0]) / dt + NAMED_AT;
    const x = (age - m.age_settle) / m.age_span * m.phase_bins;   // (bin k: the sweep's starts at k / bins of the span)
    return x > m.phase_bins - 0.5 ? null : Math.max(0, Math.round(x));
  }
  // A chaotic morph (the book's states: centroids of [speed / R, |omega|], their scale): the state it
  // is in now, the nearest to its recent kinematics (averaged over the last 15 steps, as the sweep's)
  function stateBin(m, hist, dt) {
    if (!m || !m.states || !hist || !hist.length) return null;
    const tn = hist[hist.length - 1][0], h = hist.filter(x => x[0] > tn - 15 * dt);
    const sp = h.reduce((a, x) => a + x[2], 0) / h.length, om = h.reduce((a, x) => a + Math.abs(x[1]), 0) / h.length;
    const { centroids: c, scale: k } = m.states;
    let best = 0, bd = Infinity;
    c.forEach((v, i) => { const d = ((v[0] - sp) / k[0]) ** 2 + ((v[1] - om) / k[1]) ** 2; if (d < bd) { bd = d; best = i; } });
    return best;
  }
  // Every action of a morph in every phase where the phase matters (each spot once per phase bin,
  // with its outcome there and phase: the bin, to wait for), else its actions as they are: what it
  // can do at all, whatever its phase now
  // An action's share in one phase bin, scaled down as its re-test scaled its whole p (for its main outcome)
  const phaseP = (e, b) => b[0] === e.effect && e.p_sweep ? b[1] * e.p / e.p_sweep : b[1];
  function allActions(m) {
    if (!m || !m.phase_useful) return m ? m.actions : [];
    const out = [];
    for (const e of m.actions) {
      if (!e.by_phase) { out.push(e); continue; }
      e.by_phase.forEach((b, i) => {
        if (!b || b[2] < 3) return;   // (too few trials there)
        const o = { ...e, effect: b[0], p: phaseP(e, b), phase: i, alt: e.alt_by_phase ? e.alt_by_phase[i] : null };
        if (b[3] !== undefined) o.turn = b[3]; else delete o.turn;
        out.push(o);
      });
    }
    return out;
  }
  // A composite step's action, as the experiments chose it (tools/models.py best_action): among the
  // morph's actions in every phase, those whose outcome, or runner-up, is the one wanted ('into X',
  // any turning; or by its start, 'replicate'), however reliable (p 0.1 on), the most reliable, then
  // the shortest hold; seen as the ship sees it (sg). null if none. A short-lived morph swept by age
  // (VT049W gyrans; now: its age bin, swingBin): only ages not yet passed, then the shortest wait
  // (as tools/models.py age_action); its age not known (or past the span): none of them
  function stepAction(book, m, target, sg, now = null) {
    if (!m) return null;
    const hits = ef => ef && (ef.split(' (')[0] === target || ef === target || ef.startsWith(target));
    const age = m.phase_kind === 'age' && m.phase_useful, known = now !== null && now !== undefined;
    const wait = e => age && e.phase !== undefined ? e.phase - now : 0;
    let best = null;
    for (const e of allActions(m)) for (const [ef, p] of [[e.effect, e.p], e.alt ? [e.alt[0], e.alt[1]] : [null, 0]]) {
      if (!hits(ef) || p < 0.1) continue;
      if (age && e.phase !== undefined && (!known || e.phase < now)) continue;   // (an age passed: it cannot come again)
      const c = { ...e, effect: ef, p }, r = Math.round(p * 20), rb = best ? Math.round(best.p * 20) : -1;
      if (!best || r > rb || (r === rb && (wait(e) < wait(best) || (wait(e) === wait(best) && e.hold < best.hold)))) best = c;
    }
    return best ? seen(book, m, best, sg) : null;
  }
  // The morph's entry with each action's outcome in phase bin `bin` (its main outcome there, share
  // and turn); the others as they are
  function phased(m, bin) {
    if (bin === null || !m) return m;
    return { ...m, phase_now: bin, actions: m.actions.map(e => {
      const b = e.by_phase && e.by_phase[bin];
      if (!b) return e;
      const o = { ...e, effect: b[0], p: phaseP(e, b), alt: e.alt_by_phase ? e.alt_by_phase[bin] : null };
      if (b[3] !== undefined) o.turn = b[3]; else delete o.turn;
      return o;
    }) };
  }

  return { OMEGA_SIDE, P_MIN, TURN_MIN, TURN_BIG, load, get, sign, seen, outcomeOf, maneuvers, choose, words, splitOf, splitMark, replicateOf, isWhole, swingBin, phased, allActions, mirrored, stepAction };
})();
