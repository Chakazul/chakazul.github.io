"use strict";
// ============================================================================================
//  Behaviours ([hostile_behavior:<morphs>:<behavior>]; the settings book's Hostile and Player
//  Behavior): how a creature steers by the morph model (models.js) and where it fires: a hostile
//  bot, going for the player's ship, or the player's ship on autopilot, going for the nearest
//  hostile. A behaviour is for one species and some of its morphs (what the model says they can
//  do); as another morph, ??? (not maneuverable) or part of a twin, a bot drifts and fires at the
//  ship as one with no behaviour, and the autopilot does nothing.
//
//  BEHAVIORS are written as rules, at the bottom of this file (see balanced). A behavior is a few
//  states; each is played as one morph (as:) and has a rule, read each time the creature decides:
//
//      attack: { as: 'limus', rule: s => s.near(4) ? 'turn' : s.approach() },
//
//  A rule gives one of:
//    'state'              go to that state (and read its rule now);
//    s.approach() ...     an action (the model's: a mask held at a spot for a while), or null if
//                         none helps now: do nothing;
//    s.once(action, 'x')  that action, then state x (x at once if there is no such action);
//    null                 nothing now.
//  With no opponent (a bot's: the player's ship dissolved; the autopilot's: no hostile left) for
//  LOST_WAIT, every behavior calms down and rests: its rest
//  state, as its morph, whatever the rules say; once one is back, the rest state's rule again.
//  Played as another morph than the state's, the creature first becomes it (the model's action
//  for that), then follows the rule. Decisions are for the moment only; the creature waits until
//  an action is over and it is settled before the next. fire: when it fires (along its motion).
//
//  What a rule can ask (s): the opponent (the ship, for a bot; the nearest hostile, for the
//  autopilot), distances in R:
//    s.near(r) / s.far(r)   the opponent nearer / farther than r R (no opponent: far)
//    s.facing(deg)          moving with the opponent within deg degrees of straight ahead
//    s.is(morph)            the creature is that morph now
//    s.turns()              turns it has spun in this state (its turn rate over the time)
//    s.time()               game seconds in this state
//    s.hunting              set by s.hunt (below), until the opponent is near
//  and do:
//    s.approach(tol)        turn towards where the opponent will be (nothing within tol degrees, 12)
//    s.flee(tol)            turn away from it
//    s.turnRound(min)       one big turn (at least min degrees, 90), the one ending nearest to away
//    s.become(morph)        turn into that morph
//    s.reach(key)           the best action for any outcome (models.js keys: 'replicate-zigzag', 'into-limus', ...)
//    s.run(name, next, opts) a composite maneuver of the model (a sequence of maneuvers, e.g. 'rectus-turn-via-gyrans'),
//                           step by step; then the state next. opts.turn: the turn wanted (picks the wait)
//    s.turnTo(off)          the turn wanted to face the opponent (degrees), for opts.turn
//    s.push(off)            a still morph that can be pushed out of place (VT049W, GDNQYX lithos): the push
//                           sending it towards the opponent (turned by off rad)
//  An oscillating, chaotic or short-lived morph whose phase matters steers by its actions as they work
//  in its phase now; become / reach take an action of another phase by waiting for it.
//    s.once(action, next)   (above)
//    s.hunt(r, next)        go to next, hunting: s.hunting is true until the opponent is nearer than r R
//
//  Under the rules: decide(ctx), aim(ctx), called by game.js with
//  ctx: { book (the species' model), name (the creature's named morph), sg (its turning direction),
//         heading (board angle of its heading, rad), omega (its turn rate, rad per time unit),
//         moveAng (of its motion, rad, or null if still), dist (cells to the opponent, or
//         Infinity), bearing (board angle to it), rel(t) (its place from the creature t time units
//         ahead, both going on as now; null: no opponent), dt (time units a step), tu (time units
//         a game second), R (cells), now (game seconds), mem (the creature's own memory) }
//  Angles are on the board: + clockwise as seen (rows grow downwards), as the model's turns.
// ============================================================================================

const Behaviors = (function () {
  const wrap = a => a - 2 * Math.PI * Math.round(a / (2 * Math.PI));
  const wrapDeg = d => d - 360 * Math.round(d / 360);
  const deg = a => a * 180 / Math.PI;

  const LOST_WAIT = 1.5;   // game s with no opponent (the player's ship dissolved) before resting
  const SETTLE = 35;
  // an action for one phase of an oscillating (or chaotic, or short-lived) morph, not its phase now:
  // decide again in a moment (game.js looks again a couple of steps later)
  const WAIT = Object.freeze({ wait: true });   // steps after a hold until its turn has taken effect (the model measures it 30+ after)

  // The heading error (degrees, + clockwise) the creature would have t time units ahead doing
  // nothing more: the heading it wants then (want(t), rad) less its heading going on turning at its rate
  const errorAt = (ctx, want, t) => deg(wrap(want(t) - (ctx.heading + ctx.omega * t)));

  // The action that brings the heading nearest to the one it wants, staying the same morph
  // (keeping its turning direction or not), judged when the action has taken effect: want(t) is
  // the wanted heading t time units ahead (the opponent moves meanwhile), the creature keeps
  // turning at its own rate, plus the model's turn (measured against an untouched copy). Only one
  // that ends nearer than doing nothing would by then; the name box kept counts for 20 degrees,
  // each second of hold for 3. null if none helps. So a small error gets a small turn, a large one
  // a large turn. minTurn (degrees): only actions turning at least this much (a big turn).
  function turnAction(ctx, want, minTurn = 0) {
    const { book, name, sg } = ctx, m = ctx.m || (book && book.morphs[name]);   // (as it works in its phase now)
    if (!m) return null;
    let best = null, bc = Infinity;
    for (const e of m.actions) {
      if (!(e.effect === 'keep' || e.effect === 'reverse') || e.p < Models.P_MIN || e.turn === undefined || Math.abs(e.turn) < minTurn) continue;
      const a = Models.seen(book, m, e, sg), t = (e.hold + SETTLE) * ctx.dt;
      const idle = errorAt(ctx, want, t), miss = Math.abs(wrapDeg(idle - a.turn));
      if (miss >= Math.abs(idle) - 3) continue;   // (no better than doing nothing, by then)
      const cost = miss + (e.always_named ? 0 : 20) + 3 * e.hold / 30;
      if (cost < bc) { bc = cost; best = a; }
    }
    return best;
  }

  // ---- what a rule sees and does (s) ----------------------------------------------------------
  // The action for an outcome (models.js key) in any phase, as reaction time matters: each action's
  // reliability less WAIT_COST a second it must wait for its phase (a fair action in this phase or the
  // next can beat a perfect one much later); it, if it can be done now, else WAIT; null if none
  const WAIT_COST = 0.1;   // p a second of waiting
  function waitSeconds(ctx, m, ph) {
    const now = ctx.m && ctx.m.phase_now;
    if (ph === undefined || ph === now) return 0;
    if (now === undefined || now === null) return null;   // (its phase not known yet)
    const bins = m.phase_bins, ahead = m.phase_kind === 'age' ? ph - now : ((ph - now) % bins + bins) % bins;
    if (ahead < 0) return null;   // (too young for it no more)
    const span = m.phase_kind === 'swing' ? m.swing_period : m.phase_kind === 'age' ? m.age_span : 30 * bins;   // (a state: about a second each)
    return ahead * span / bins / (ctx.tu / ctx.dt) ;
  }
  function reach(ctx, key) {
    const m = ctx.mAll || ctx.book.morphs[ctx.name];
    let best = null, bs = -Infinity;
    for (const e of Models.allActions(m)) {
      const a = Models.seen(ctx.book, m, e, ctx.sg);
      if (Models.outcomeOf(a) !== key) continue;
      const w = waitSeconds(ctx, m, e.phase);
      if (w === null) continue;
      const sc = e.p - WAIT_COST * w - 0.01 * e.hold / 30;   // (then the shorter hold)
      if (sc > bs) { bs = sc; best = a; }
    }
    if (!best) return null;
    return waitSeconds(ctx, m, best.e.phase) === 0 ? best : WAIT;
  }

  // A still morph pushed out of place (the model's push: VT049W lithos into tachus, GDNQYX lithos
  // into tardus): a mask held on its side opposite the way wanted (towards where the opponent will be,
  // turned by off rad), in the board's frame; null if it can't be pushed or has no opponent
  // From the model's push table (lithos-types): each reliable push (p 3 of 4) of the kind (convert: it
  // becomes a mover; shove: it stays and slides) that leaves one known way (away from the mask, or back
  // towards it; not sideways), at each place the morph allows (round: any side; dihedral: along,
  // diagonal or across its long axis, at the 4 mirror places), the one whose way out is nearest the way
  // wanted (within 60 degrees), more reliable and shorter first.
  function push(ctx, off, kind = 'convert') {
    const m = ctx.mAll || ctx.book.morphs[ctx.name], tab = m && m.push_table;
    if (!tab || !ctx.rel) return null;
    const dih = m.lithos_type === 'dihedral';
    if (dih && ctx.longAxis === null) return null;
    let best = null, bs = -Infinity;
    for (const pt of tab) {
      const into = /^into (\S+)/.exec(pt.effect);
      const ok = kind === 'convert' ? !!into : pt.effect === 'keep' && pt.move >= 0.3 && (pt.away ?? 0) >= 0.5;
      if (!ok || pt.p < Models.P_MIN || pt.dev === null || pt.dev === undefined || (pt.dev >= 45 && pt.dev < 135)) continue;
      const r = ctx.rel(pt.hold * ctx.dt), want = Math.atan2(r[1], r[0]) + off;
      const back = pt.dev >= 135;   // (it leaves towards the mask's side)
      const masks = !dih ? [want + (back ? 0 : Math.PI)]
        : (pt.side === 'along' ? [0, Math.PI] : pt.side === 'across' ? [Math.PI / 2, -Math.PI / 2] : [Math.PI / 4, -Math.PI / 4, 3 * Math.PI / 4, -3 * Math.PI / 4]).map(a => ctx.longAxis + a);
      for (const th of masks) {
        const out = th + (back ? 0 : Math.PI), miss = Math.abs(deg(wrap(out - want)));
        if (miss > 60) continue;
        const sc = pt.p - miss / 180 - 0.01 * pt.hold / 30;
        if (sc > bs) { bs = sc; const d = pt.r * ctx.R; best = { pt, th, d }; }
      }
    }
    if (!best) return null;
    const { pt, th, d } = best;
    return { e: { ch: 'G', spot: `push (${pt.dist}${pt.side ? ', ' + pt.side : ''})`, hold: pt.hold, p: pt.p, effect: pt.effect, always_named: false },
             fs: [d * Math.cos(th), d * Math.sin(th)], frame: 'board', spot: 'push', turn: null };
  }

  // ---- composite maneuvers: sequences of maneuvers (the model's, tried in experiments) -------------
  // A composite's steps ('into X [xN]', 'do <effect>', 'wait N' steps), run closed loop: each decision
  // does the step it is at (becoming X: its best action now, or WAIT for a phase), moving on when the
  // creature is X, or N steps have passed. With a {w} (e.g. the wait that sets a turn): the value tried
  // whose net turn is nearest opts.turn (degrees), else the most reliable. Done (or failed: a step's
  // morph not reached in its tries, or not the composite's morph to begin with): the state `next`.
  const compositeOf = (ctx, name) => ((ctx.book && ctx.book.composites) || []).find(s => s.name === name) || null;
  function pickVal(st, opts) {
    const rs = st.results.filter(r => r.p >= 0.5);
    if (!rs.length) return null;
    if (opts.turn !== undefined) return rs.reduce((a, r) => !a || Math.abs(wrapDeg(r.turn - opts.turn)) < Math.abs(wrapDeg(a.turn - opts.turn)) ? r : a, null);
    return rs.reduce((a, r) => !a || r.p > a.p ? r : a, null);
  }
  function run(ctx, name, next, opts = {}) {
    const { mem } = ctx, st = compositeOf(ctx, name);
    let sr = mem.comp;
    if (!sr || sr.name !== name) {
      if (!st || ctx.name !== st.morph) return next;
      const r = pickVal(st, opts);
      if (!r) return next;
      const steps = (st.steps.replace('{w}', r.val ?? '')).split(';').map(x => x.trim().split(/\s+/)).filter(w => w[0]);
      sr = mem.comp = { name, steps, i: 0, t0: ctx.now, tries: 0, expect: r };
    }
    for (let hop = 0; hop < 4; hop++) {
      const w = sr.steps[sr.i];
      if (!w) { mem.comp = null; return next; }   // (done)
      if (w[0] === 'wait') {
        if ((ctx.now - sr.t0) * ctx.tu / ctx.dt < +w[1]) return WAIT;   // (steps: game s x steps a second)
        sr.i++; sr.t0 = ctx.now; sr.tries = 0; continue;
      }
      if (w[0] === 'into') {
        if (ctx.name === w[1]) { sr.i++; sr.t0 = ctx.now; sr.tries = 0; continue; }
        const most = w[2] && w[2][0] === 'x' ? +w[2].slice(1) : 1;
        if (sr.tries >= most) { mem.comp = null; return next; }   // (failed)
        if (sr.tries && (ctx.now - sr.t0) * ctx.tu / ctx.dt < 150) return WAIT;   // (a try settling and named: the name box needs its history)
        // (the step's action as the experiments chose it, however reliable; a phase's: waited for)
        const m = ctx.mAll || ctx.book.morphs[ctx.name], sa = Models.stepAction(ctx.book, m, `into ${w[1]}`, ctx.sg, ctx.m && ctx.m.phase_now);
        const a = sa && sa.e.phase !== undefined && !(ctx.m && sa.e.phase === ctx.m.phase_now) ? WAIT : sa;
        if (a && a !== WAIT) sr.t0 = ctx.now + a.e.hold * ctx.dt / ctx.tu;
        if (a && a !== WAIT) sr.tries++;
        return a || (mem.comp = null, next);
      }
      if (w[0] === 'hold') {   // hold <spot | fwd,side in R> <steps> [ch]: a mask held there, as tried
        if (sr.tries) { sr.i++; sr.t0 = ctx.now; sr.tries = 0; continue; }
        const m = ctx.mAll || ctx.book.morphs[ctx.name], chn = ['R', 'G', 'B'].includes(w[w.length - 1]) ? w.pop() : 'G';
        const nst = +w[w.length - 1], spot = w.slice(1, -1).join(' ');
        const off = /^-?[0-9.]+,-?[0-9.]+$/.test(spot) ? spot.split(',').map(Number) : ((m.actions.find(e => e.spot === spot) || {}).off);
        if (!off) { mem.comp = null; return next; }
        sr.tries++;
        return { e: { ch: chn, spot, hold: nst, p: null, effect: 'hold', always_named: false }, fs: [off[0] * ctx.R, off[1] * ctx.R], spot, turn: null };
      }
      if (w[0] === 'do') {
        if (sr.tries) { sr.i++; sr.t0 = ctx.now; sr.tries = 0; continue; }
        const m = ctx.mAll || ctx.book.morphs[ctx.name];
        const key = (Models.maneuvers(ctx.book, ctx.name, ctx.sg, m, true).find(o => o.key.startsWith(w.slice(1).join(' ').split(':')[0])) || {}).key;
        const a = key ? reach(ctx, key) : null;
        if (a && a !== WAIT) sr.tries++;
        return a || (mem.comp = null, next);
      }
      sr.i++;
    }
    return WAIT;
  }
  const become = (ctx, morph) => ctx.name === morph ? null : reach(ctx, `into-${morph}`) || reach(ctx, `into-${morph}-rev`);

  function situation(ctx) {
    const { mem, R } = ctx, d = ctx.dist;
    // the heading wanted t ahead: to where the opponent will be, plus off (rad)
    const toward = off => t => { const r = ctx.rel(t); return Math.atan2(r[1], r[0]) + off; };
    const steer = (off, tol) => {
      if (!ctx.rel) return null;   // (no opponent)
      const want = toward(off);
      return Math.abs(errorAt(ctx, want, 15 * ctx.dt)) < tol ? null : turnAction(ctx, want);   // (the error half a second ahead)
    };
    return {
      near: r => d < r * R,
      far: r => d > r * R,
      facing: a => ctx.moveAng !== null && Number.isFinite(d) && Math.abs(deg(wrap(ctx.bearing - ctx.moveAng))) <= a,
      is: morph => ctx.name === morph,
      turns: () => (mem.spun || 0) / (2 * Math.PI),
      time: () => ctx.now - (mem.since ?? ctx.now),
      get hunting() { return !!mem.hunt; },
      approach: (tol = 12) => steer(0, tol),
      flee: (tol = 12) => steer(Math.PI, tol),
      turnRound: (min = 90) => ctx.rel ? turnAction(ctx, toward(Math.PI), min) : null,
      become: morph => become(ctx, morph),
      run: (name, next, opts) => run(ctx, name, next, opts),
      // the turn wanted to face where the opponent will be (degrees, + clockwise), for a composite's opts.turn
      turnTo: (off = 0) => ctx.rel ? wrapDeg(deg(Math.atan2(ctx.rel(30 * ctx.dt)[1], ctx.rel(30 * ctx.dt)[0]) + off - ctx.heading)) : 0,
      reach: key => reach(ctx, key),
      push: (off = 0, kind = 'convert') => push(ctx, off, kind),
      once: (action, next) => action === WAIT ? WAIT : action ? { action, next } : next,
      hunt: (r, next) => { mem.hunt = r; mem.huntTo = next; return next; },
    };
  }

  // ---- a behavior: its rules run as decide / aim -------------------------------------------------
  // def: { species: [codes], morphs: [names], start: s => 'state', states: { name: { as, rule } },
  //        fire: s => true | false }
  function behavior(def) {
    // a new state: its clock and spin from 0; hunting only into the state s.hunt named
    const enter = (mem, state, now) => { mem.mode = state; mem.since = now; mem.spun = 0; if (state !== mem.huntTo) mem.hunt = null; mem.huntTo = null; };
    return {
      species: def.species, morphs: def.morphs, states: def.states,
      decide(ctx) {
        const { mem } = ctx, s = situation(ctx);
        if (!mem.mode || !def.states[mem.mode]) enter(mem, def.start(s), ctx.now);
        if (mem.after) { const n = mem.after; mem.after = null; enter(mem, n, ctx.now); }   // (once's next, now its action is done)
        // bookkeeping: how far it has spun as its state's morph (its turn rate over the time since it last looked)
        const dt = Math.max(0, ctx.now - (mem.looked ?? ctx.now));
        mem.looked = ctx.now;
        if (ctx.name === def.states[mem.mode].as) mem.spun = (mem.spun || 0) + Math.abs(ctx.omega) * dt * ctx.tu;
        if (mem.hunt && s.near(mem.hunt)) mem.hunt = null;   // (hunting: found it)
        // no opponent (the player's ship dissolved): a moment as it was, then (after LOST_WAIT) rest (the rest state, as its
        // morph) whatever the rules say, until one is back: then rest afresh, by its rule
        if (!ctx.rel) {
          mem.lost ??= ctx.now;
          const rest = def.states.rest;
          if (rest && ctx.now - mem.lost < LOST_WAIT) return null;   // (a moment as it was)
          if (rest) {
            if (mem.mode !== 'rest' || !mem.restLost) enter(mem, 'rest', ctx.now);
            mem.after = null; mem.restLost = true;
            return rest.as && ctx.name !== rest.as ? s.become(rest.as) : null;
          }
        } else if (mem.lost != null) {
          if (mem.restLost) enter(mem, 'rest', ctx.now);
          mem.lost = null; mem.restLost = false;
        }
        for (let hop = 0; hop < 6; hop++) {   // (a few state changes at most, in one decision)
          const st = def.states[mem.mode];
          if (st.as && ctx.name !== st.as) return s.become(st.as);   // played as another morph: become it first
          const r = st.rule(s);
          if (r === WAIT) return WAIT;
          if (typeof r === 'string') {
            if (!def.states[r]) { console.warn(`behaviors.js: no state ${r}`); return null; }
            if (r === mem.mode) return null;
            enter(mem, r, ctx.now);
            continue;
          }
          if (r && r.next !== undefined && r.action) { mem.after = r.next; return r.action; }   // once(action, next): next after it
          return r || null;
        }
        return null;
      },
      // where it fires: along its motion, when fire says (null: not now, the shot kept)
      aim(ctx) { return ctx.moveAng !== null && def.fire(situation(ctx)) ? ctx.moveAng : null; },
    };
  }

  // ============================================================================================
  //  BEHAVIORS
  // ============================================================================================
  const LIST = {
    // balanced (5N7KKM limus, gyrans): attack when the opponent is within reach, turn round and
    // escape when it is too close, rest once safely away, until it comes within reach again (or
    // after two turns of resting: then hunt it from anywhere).
    balanced: behavior({
      species: ['5N7KKM'], morphs: ['limus', 'gyrans'],
      start: s => s.near(7) ? 'attack' : 'rest',
      states: {
        attack: { as: 'limus',  rule: s => s.near(4) ? 'turn'
                                         : s.far(8) && !s.hunting ? 'rest'
                                         : s.approach() },
        turn:   { as: 'limus',  rule: s => s.once(s.turnRound(90), 'escape') },
        escape: { as: 'limus',  rule: s => s.far(8) ? 'rest' : s.flee() },
        rest:   { as: 'gyrans', rule: s => s.near(7) ? 'attack'
                                         : s.turns() >= 2 ? s.hunt(7, 'attack')
                                         : null },
      },
      fire: s => s.is('limus') && s.near(8) && s.facing(15),
    }),

    // aggressive (5N7KKM limus, gyrans): no rest, little escape: attack the opponent from
    // anywhere; only when it is very close, turn round and escape, but briefly (until a little
    // away, or 3 s), then attack again (a gyrans turns into limus first).
    aggressive: behavior({
      species: ['5N7KKM'], morphs: ['limus', 'gyrans'],
      start: s => 'attack',
      states: {
        attack: { as: 'limus', rule: s => s.near(2.5) ? 'turn' : s.approach() },
        turn:   { as: 'limus', rule: s => s.once(s.turnRound(90), 'escape') },
        escape: { as: 'limus', rule: s => s.far(5) || s.time() > 3 ? 'attack' : s.flee() },
        rest:   { as: 'gyrans', rule: s => 'attack' },   // (only with no opponent, see decide)
      },
      fire: s => s.is('limus') && s.near(7) && s.facing(10),
    }),
  };

  const get = name => LIST[name] || null;
  const has = name => name in LIST;
  return { get, has, behavior, turnAction, WAIT, names: () => Object.keys(LIST) };
})();
