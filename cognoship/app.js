"use strict";
// ============================================================================================
//  Cognoship, Phases 1 and 1b: the controllable ship, a pattern picker, and several creatures
//  at once with a clickable focus (see PLAN.md sections 4 and 4.5).
//
//  CPU side: layout, input, finding and following creatures, the heading / view logic, the
//  overlay (tool patches, velocity arrow, tool outline) and the readout. The simulation and the
//  tools run on the GPU (glsim.js). Each frame the CPU reads back a coarse block-mass map and
//  labels its connected components (creatures.js); one of them is the focus, which the zoom
//  view, the arrow and the tools follow.
//
//  Coordinates. Board cells are (col, row), row 0 at the top of the screen. CoMs are in cell
//  index units (cell i's centre is at i). "View" vectors are the zoom panel's (right, down) in
//  cells; rot() maps them to board vectors, so the tool lands under the pointer in any view mode.
// ============================================================================================

(function () {
  // ---- parameters -----------------------------------------------------------------------
  const q = new URLSearchParams(location.search);
  const num = (k, d, lo, hi) => {
    const v = parseFloat(q.get(k));
    return Number.isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d;
  };
  // Kernel radius, cells. 18: below 16, gyrans R stops replicating and lithos moves (README).
  const R = Math.round(num('R', 18, 6, 30));
  const DT_OVERRIDE = q.has('dt') ? num('dt', 0.2, 0.02, 1) : null;   // else the pattern's own dt
  const FORM_ALIASES = { aquarium: 'VT049W-rectus', rectus: 'VT049W-rectus',
    replicating: 'VT049W-gyrans-replicating', gyrans: 'VT049W-gyrans-replicating', gyransR: 'VT049W-gyrans-replicating' };
  const findPattern = id => PATTERNS.find(p => p.id === id) || PATTERNS.find(p => p.id === DEFAULT_PATTERN);
  const DEBUG = q.get('debug') === '1';

  // Canvas layout, in canvas px (the canvas is 1152 x 512 and scales with CSS).
  const CW = 1152, CH = 512;
  const BW = 384, BH = 256, WS = 2;               // board cells; world view px per cell
  const PX = BW * WS;                              // panel starts here (left 2/3 is the world)
  const PW = CW - PX;                              // panel width
  // Two rows of tool patches: top = channel (R G B), bottom = action (Mask Erase Draw); zoom view between.
  const PATCH = { top: 6, bottom: 452, h: 54, w: 104, gap: 16 };
  const ZOOM = { x: PX, y: 64, edge: PW };          // zoom view square
  // Zoom slider: px per cell in the zoom view (7 slots); the default shows about 4 R across.
  const ZOOMS = [2, 3, 4, 6, 8, 12, 16];
  const nearestIdx = (list, v) => list.reduce((b, x, i) => Math.abs(x - v) < Math.abs(list[b] - v) ? i : b, 0);
  let zoomIdx = nearestIdx(ZOOMS, num('zoom', ZOOM.edge / (4 * R), 1, 32));
  let CELL_PX = ZOOMS[zoomIdx], ZOOM_CELLS = ZOOM.edge / CELL_PX;   // cells across the zoom view

  // Kinematics (speeds in cells per unit time, times in units of simulated time).
  // The Aquarium glides at ~0.072 R per unit time (t = steps / T, T = 2).
  const SPEED_MIN = 0.005 * R;   // below: near stationary, arrow hidden (the tardus forms are slow)
  const SPEED_HI = SPEED_MIN * 1.5;  // hysteresis: "moving" switches on above this
  const VEL_TAU = 15;            // heading-up smoothing, unit time (= 30 steps; pulses, detection jitter)
  const ARROW_TAU = 2.5;         // arrow smoothing, unit time (~5 steps)
  const OMEGA_TAU = 5;           // turn-rate smoothing, unit time
  const JUMP_MAX = Math.PI / 4;  // heading residuals are capped at this (detection jumps)
  const LAG_MAX = Math.PI / 3;   // cap on the turn-lag correction (60 degrees)
  const ARROW_T = 12;            // arrow = displacement over this much time, at zoom scale

  // Creatures (creatures.js): basins of the potential K * A above 0.05, merged across saddles
  // >= 0.8 x the lower peak, at least a disk of radius R/4 (the values of ~/github/complexity).
  const MIN_AREA_UNIT = Math.PI * (R / 4) ** 2;
  const DETECT = { threshold: num('thr', 0.1, 0.01, 1), saddle: num('saddle', 0.9, 0.1, 2),
                   minArea: num('area', 1, 0, 10) * MIN_AREA_UNIT };   // mutable: tuning-mode sliders
  // Page mode, cycled with A (Shift+A backwards): action (play) -> stats (the chart in place of the
  // ship view) -> action+stats (the chart left of the ship view, over part of the world view) ->
  // tuning (detection sliders). Always starts in action mode (no URL parameter).
  const MODES = ['action', 'stats', 'action+stats', 'tuning'];
  const MODE_NAMES = { action: 'Action', stats: 'Stats', 'action+stats': 'Action+Stats', tuning: 'Tuning' };
  let mode = 'action';
  const statsOn = () => mode === 'stats' || mode === 'action+stats';   // the chart is shown
  // Morph detection (Detect button, detect.js): none, the focus only (Focused: a name
  // box at the top of the ship view), or every object (Every, the default: also its circle in its morph's hue and its
  // name in the world view).
  const DETECT_MODES = ['none', 'focus', 'every'];
  const DETECT_NAMES = { none: 'None', focus: 'Focused', every: 'Every' };
  let detectMode = 'every';   // the default
  const DETECT_PER_FRAME = 4;   // every-object detection: tracks classified per frame, in turn
  const md = { model: null, len: -1, next: 0, dirty: true };   // the long history's model, rebuilt as it changes
  const TRACK_GATE = R;          // the focus may move at most this far between frames (cells)
  const CLICK_GATE = 2 * R;      // a click this close to a creature selects it

  const CH_COLORS = ['#ff3b3b', '#3bff5a', '#3b7bff'];
  const CH_FILLS = ['rgba(255,59,59,.5)', 'rgba(59,255,90,.5)', 'rgba(59,123,255,.5)'];
  const CH_NAMES = ['R', 'G', 'B'];
  // Actions: 0 mask (hide from the senses), 1 erase (set to 0), 2 draw (set to 1). Circle line:
  // dashed, solid, dotted.
  const ACTIONS = ['Mask', 'Erase', 'Draw'];
  const ACTION_DASH = [[6, 4], [], [2, 4]];

  // ---- state ----------------------------------------------------------------------------
  // Speed slider: 7 slots. Each slot runs SPEEDS[i] * STEPS_PER_SPEED steps per frame, so
  // speed 1 is one step every 4 frames (under LeniaF's rule the glider moves ~0.65 cells a step).
  const SPEEDS = [1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8], STEPS_PER_SPEED = 1 / 4;
  let running = true, speedIdx = nearestIdx(SPEEDS, num('speed', 1, 1 / 8, 8)), stepAcc = 0;
  let maskR = Math.round(num('radius', 2.5 * R / 12, 1, 8)), channel = 1, action = 0;   // selected tool: G mask (best for steering and shape-shifting)
  // Zoom-view rotation modes:
  //   board: no rotation; vel: heading = the arrow's smoothed velocity (arrow points straight up);
  //   axis: heading = the B->G body axis (stable at any speed, but not the travel direction:
  //         ~56 degrees off for the circling gyrans form, a few degrees of wobble for the Aquarium).
  // 'axis' (body axis up) is kept only for ?headup=axis; the button toggles board / heading up.
  const VIEW_LABELS = { board: 'View: board up', vel: 'View: heading up', axis: 'View: body axis up' };
  let viewMode = { '1': 'vel', vel: 'vel', axis: 'axis' }[q.get('headup')] || 'board';

  // '+' in a query string decodes to a space ("lithos+gyrans"), so map spaces back.
  const formRaw = (q.get('form') || '').replace(/ /g, '+');
  let pattern = findPattern(FORM_ALIASES[formRaw] || formRaw);
  let dt = DT_OVERRIDE ?? pattern.dt;

  let steps = 0, spawnMass = null, status = 'alive';
  let found = { comps: [], labels: new Int32Array(0), total: 0 };   // last creature scan
  let focus = null;          // { col, row, mass, radius } of the followed creature, or null
  let focusIdx = -1;         // its index in found.comps (this scan)
  let autoFocus = false;     // after a spawn: focus the creature nearest the board centre
  let kin = null;
  let masking = null;        // { id, view: [x, y] } while the pointer holds a tool
  let hover = null;          // [x, y] view offset of a mouse hovering over the ship view
  let msPerStep = 0;

  const glCanvas = document.getElementById('gl');
  const ui = document.getElementById('ui');
  const ctx = ui.getContext('2d');
  const readout = document.getElementById('readout');

  // vel/speed/velAng: the arrow's velocity (short average, turn-lag corrected), also used to
  // predict the focus position. slowHeading: the heading-up view's angle. omega: turn rate
  // (rad per unit time). age: time measured since the focus was set. See updateKinematics.
  const newKinematics = () => ({ vel: [0, 0], speed: 0, velAng: 0, fast: [0, 0], fastAng: null, omega: 0, omegaAge: 0,
                                 slowHeading: null, age: 0, moving: false, heading: null, axisHeading: null });
  function resetKinematics() { kin = newKinematics(); }

  // ---- helpers --------------------------------------------------------------------------
  const wrapD = (d, n) => d - n * Math.round(d / n);
  const deg = a => ((a * 180 / Math.PI) + 360) % 360;

  // Rotation taking view (right, down) to board (col, row). Heading-up puts heading at "up".
  function viewAngle() {
    const h = viewMode === 'vel' ? kin.heading : viewMode === 'axis' ? kin.axisHeading : null;
    return h === null ? 0 : h + Math.PI / 2;
  }
  function rot(v, phi = viewAngle()) {
    const c = Math.cos(phi), s = Math.sin(phi);
    return [v[0] * c - v[1] * s, v[0] * s + v[1] * c];
  }
  function unrot(b, phi = viewAngle()) { return rot(b, -phi); }

  // The tool circle's offset from the focus CoM (board frame), or null if no tool is held.
  function toolOffset() {
    if (!masking || !focus) return null;
    return rot(masking.view);
  }

  // ---- simulation loop ------------------------------------------------------------------
  function setPattern(id) {
    pattern = findPattern(id);
    const rule = patternRule(pattern);
    dt = DT_OVERRIDE ?? rule.dt;
    SimGL.setRule(rule.kernels, R, rule.kn);
    SimGL.setUpdate(rule.norm, dt, rule.gn);
    // merge saddle per species (aquarium233.js rules[code].detect), unless ?saddle= is given
    if (stats.species !== speciesOf(pattern)) {
      stats.recLeft = 0; syncRecord();   // a recording in progress is discarded with the long history
      stats.short.clear(); stats.long.clear(); tracks = []; stats.species = speciesOf(pattern);
      // locked axes: new ranges for the new species, locked again once its saved history has
      // loaded, or (no file) once the short history has filled (frame)
      if (chart.lock) { StatsChart.resetRanges(); chart.relock = true; }
      loadSavedStats(stats.species);
    }
    if (!q.has('saddle')) DETECT.saddle = AQUARIUM233.rules[pattern.rule]?.detect?.saddle ?? 0.9;
    const ts = $('t-saddle');
    ts.value = DETECT.saddle; $('t-saddle-v').textContent = DETECT.saddle.toFixed(2); fillRange(ts);
    syncPicker();
  }

  function spawn() {
    const f = Rule.formCells(pattern.cells, pattern.R, R);
    SimGL.seed(f, Math.floor(BH / 2 - f.h / 2), Math.floor(BW / 2 - f.w / 2));
    steps = 0; spawnMass = null; status = 'alive'; stepAcc = 0;
    focus = null; autoFocus = true;
    resetKinematics(); tracks = []; stats.short.clear();   // a new world: no short histories
    analyze(0);
  }

  // Add the pattern, rotated by a random angle, at a random spot: the emptiest of 24 random candidates (least mass within
  // the form's half-size + R of its centre). Focuses it if nothing is focused.
  function addPattern() {
    const f = Rule.rotate(Rule.formCells(pattern.cells, pattern.R, R), Math.random() * 2 * Math.PI);
    const reach = Math.max(f.w, f.h) / 2 + R, B = SimGL.BLOCK, blk = SimGL.readBlocks();
    let best = null, bestM = Infinity;
    for (let k = 0; k < 24; k++) {
      const col = Math.random() * BW, row = Math.random() * BH;
      let m = 0;
      for (let by = 0; by < SimGL.BH; by++)
        for (let bx = 0; bx < SimGL.BW; bx++) {
          const d = Math.hypot(wrapD(bx * B + B / 2 - col, BW), wrapD(by * B + B / 2 - row, BH));
          if (d < reach) m += blk[(by * SimGL.BW + bx) * 4];
        }
      if (m < bestM) { bestM = m; best = [col, row]; }
    }
    SimGL.add(f, Math.floor(best[1] - f.h / 2), Math.floor(best[0] - f.w / 2));
    if (status !== 'alive') { status = 'alive'; running = true; }
    analyze(0);
    // follow the new creature
    const i = Creatures.pick(found, best[0], best[1], 2 * R, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
    if (i >= 0) setFocus(i);
    syncButtons();
  }

  // Run n steps. The tool's centre is extrapolated along the focus velocity within the frame,
  // since the focus CoM is only measured once per frame.
  function advance(n) {
    const t0 = DEBUG ? performance.now() : 0;
    const off = toolOffset();
    for (let i = 0; i < n; i++) {
      let tool = null;
      if (off) {
        const col = focus.col + kin.vel[0] * dt * i + off[0];
        const row = focus.row + kin.vel[1] * dt * i + off[1];
        tool = { col, row, r: maskR, ch: channel, action };
      }
      SimGL.step(tool);
    }
    steps += n;
    if (DEBUG) { SimGL.finish(); msPerStep = 0.9 * msPerStep + 0.1 * (performance.now() - t0) / n; }
  }

  // Find the creatures, follow the focus, update velocity and heading. n = steps since last call.
  function analyze(n) {
    const blk = SimGL.readBlocks();
    found = Creatures.find(blk, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH, DETECT);

    // No "Exploded" check: a board full of copies just keeps running.
    if (found.comps.length === 0) { status = 'dissolved'; running = false; syncButtons(); }

    const T = n * dt;
    if (autoFocus) {
      autoFocus = false;
      setFocus(Creatures.pick(found, BW / 2, BH / 2, Infinity, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH));
      if (focus) spawnMass = focus.mass;
    } else if (focus) {
      // Track: look where the focus should be now. Inside a component -> that one (also keeps
      // the focus through a merge); else the nearest within the gate (on a split, the nearer
      // half). Nothing -> the creature died: no focus.
      const pc = focus.col + kin.vel[0] * T, pr = focus.row + kin.vel[1] * T;
      const i = Creatures.pick(found, pc, pr, TRACK_GATE, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
      focusIdx = i;
      if (i < 0) { focus = null; resetKinematics(); masking = null; stopRecord(); }
      else {
        const c = found.comps[i];
        if (T > 0) updateKinematics([wrapD(c.col - focus.col, BW) / T, wrapD(c.row - focus.row, BH) / T], T);
        focus = { col: c.col, row: c.row, mass: c.mass, radius: c.radius };
      }
    }
    if (n > 0 && tracking()) trackObjects(T);
    if (!focus) { if (n > 0) recordStats(); return; }   // no focus: the short history still drains

    // 'vel' heading: the long-averaged, turn-lag corrected heading (see updateKinematics). Near
    // stationary the last heading is held.
    if (kin.moving || kin.heading === null) kin.heading = kin.slowHeading;

    // 'axis' heading: from the B channel's CoM to the G channel's within the focus (G leads B
    // by ~5 cells in the Aquarium). Falls back to the velocity heading if the axis is lost.
    if (viewMode === 'axis') {
      const c = SimGL.readChannelComs(focus.col, focus.row, 2 * focus.radius + 4);
      const g = { row: c[4], col: c[5], m: c[6] }, b = { row: c[8], col: c[9], m: c[10] };
      const ax = [wrapD(g.col - b.col, BW), wrapD(g.row - b.row, BH)];
      const valid = c[7] > 0.5 && c[11] > 0.5 && Math.hypot(ax[0], ax[1]) > 1 &&
                    g.m > 0.05 * focus.mass && b.m > 0.05 * focus.mass;
      if (valid) kin.axisHeading = Math.atan2(ax[1], ax[0]);
      else if (kin.axisHeading === null || kin.moving) kin.axisHeading = kin.heading;
    }

    if (n > 0) recordStats();
  }

  // Saved long history of a species: assets/stats/<species>.json (from Download). Loaded into the
  // long history when the species is (none: nothing happens). Skipped if it was recorded at
  // another R or dt (the creatures behave differently there), or if recording began meanwhile.
  function loadSavedStats(sp) {
    fetch(`assets/stats/${encodeURIComponent(sp)}.json`, { cache: 'no-cache' })
      .then(r => r.ok ? r.json() : null)
      .then(obj => {
        if (!obj || stats.species !== sp || stats.long.len) return;
        if (obj.R !== R || Math.abs(obj.dt - dt) > 1e-9) {
          console.info(`Cognoship: saved stats for ${sp} skipped (recorded at R=${obj.R}, dt=${obj.dt})`);
          return;
        }
        stats.long.load(obj);
        if (chart.relock) { StatsChart.resetRanges(); chart.relock = false; }   // lock on the loaded history
      })
      .catch(() => {});
  }

  // ---- stats: time series of the focus ------------------------------------------------
  // One sample per sim step (frame() analyses after every step): t, total mass, mass per channel
  // (R, G, B), speed and turn rate, all independent of the kernel radius R, in a spaceship's units
  // (1 m = R cells, 1 s = unit time): t in s, mass in g = sum of cell values / R^2 (over the
  // focus's blocks, as found.comps.mass), speed in m/s, turn rate in rad/s.
  // Speed and turn rate are the arrow's (updateKinematics): speed = |kin.vel|, turn rate =
  // |kin.omega| in rad per unit time, how fast the arrow's direction turns, either way. Two
  // histories (statschart.js): short = the last 50 steps, always; long = unlimited, filled by
  // Record (button, key S) REC_STEPS steps at a time, each sample with the morph selected then.
  // A change of focus (spawn, reset, click, add) starts a new piece in both, and each Record in
  // the long one: a gap in the line. A change of species clears both.
  const stats = { short: StatsChart.ring(), long: StatsChart.log(), species: null,
                  recLeft: 0, shortGap: false, longGap: false,   // recLeft: steps still to record
                  recName: '', recFuzzy: false };   // the morph name being recorded; fuzzy (Shift+S, detect.js)
  const NO_SAMPLE = { t: NaN, mass: NaN, mR: NaN, mG: NaN, mB: NaN, speed: NaN, omega: NaN };
  const REC_STEPS = 300;   // one press of Record (S) records this many steps of the focus
  window.cognoshipStats = stats;   // for poking at it from the console: .short / .long .get(key, j)

  function recordStats() {
    if (!focus || focusIdx < 0) {   // no focus: an empty sample in the short history, so it ages out
      stats.short.push(NO_SAMPLE, false);
      return;
    }
    const s = sampleOf(focusIdx, kin);
    stats.short.push(s, stats.shortGap); stats.shortGap = false;
    if (stats.recLeft > 0) {
      stats.long.push(s, stats.longGap, stats.recName, stats.recFuzzy);
      stats.longGap = false;
      if (--stats.recLeft === 0) syncRecord();
    }
  }

  // Channel masses of every object of this scan: found.chan[3 * i + c], summed over its blocks.
  function channelMasses() {
    if (found.chan) return found.chan;
    const lab = found.labels, ch = SimGL.blockChannels, out = new Float64Array(3 * found.comps.length);
    for (let k = 0; k < lab.length; k++) {
      const i = lab[k];
      if (i < 0) continue;
      out[3 * i] += ch[4 * k]; out[3 * i + 1] += ch[4 * k + 1]; out[3 * i + 2] += ch[4 * k + 2];
    }
    return (found.chan = out);
  }
  // The stats sample of object i (this scan) with kinematics k, in R units.
  function sampleOf(i, k) {
    const ch = channelMasses(), m = 1 / (R * R);
    return { t: steps * dt, mass: found.comps[i].mass * m, mR: ch[3 * i] * m, mG: ch[3 * i + 1] * m, mB: ch[3 * i + 2] * m,
             speed: k.speed / R, omega: Math.abs(k.omega) };
  }

  // Every-object mode (stats panel button, ?every=1): every object besides the focus keeps a short
  // history too, but is drawn as one dot at its current stats; the focus as its fading line. On
  // a change of focus the two swap (setFocus): the new focus brings its history and filters, the
  // old one becomes a track, so its line shows again if it is focused again. Each step
  // the tracks are matched to this scan's objects: greedily, nearest predicted position first,
  // within TRACK_GATE, the focus's object excluded. A matched track updates its own kinematics
  // (the arrow's filters, as for the focus) and gets a sample; an unmatched one gets an empty
  // sample (no dot) and is dropped after TRACK_MISS_MAX steps unmatched; an unmatched object
  // starts a new track. The long history still records the focus only.
  const TRACK_MISS_MAX = 10;   // keeps a track's filters through a brief detection dropout
  let everyObject = true;   // on by default (not in the URL)
  const tracking = () => everyObject || detectMode === 'every';   // other objects keep tracks
  let worldHover = null;                        // [col, row] of a mouse over the world view, or null
  // Every-object mode: the object under the mouse in the world view (not the focus), and its
  // track: { i, tr } (tr null if it has none yet), or null.
  function hovered() {
    if (!statsOn() || !everyObject || !worldHover) return null;
    const i = Creatures.pick(found, worldHover[0], worldHover[1], R / 2, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
    if (i < 0 || i === focusIdx) return null;
    const c = found.comps[i];
    let tr = null, best = 2;   // cells: a track sits where it was last matched, on this scan's object
    for (const t of tracks) {
      const d = Math.hypot(wrapD(t.col - c.col, BW), wrapD(t.row - c.row, BH));
      if (d < best) { best = d; tr = t; }
    }
    return { i, tr };
  }
  let tracks = [];   // { col, row, kin, hist, miss }
  function trackObjects(T) {
    const comps = found.comps, used = new Uint8Array(comps.length);
    if (focusIdx >= 0) used[focusIdx] = 1;
    const pairs = [];
    tracks.forEach((tr, a) => {
      const pc = tr.col + tr.kin.vel[0] * T, pr = tr.row + tr.kin.vel[1] * T;
      comps.forEach((c, b) => {
        const d = Math.hypot(wrapD(c.col - pc, BW), wrapD(c.row - pr, BH));
        if (d < TRACK_GATE) pairs.push([d, a, b]);
      });
    });
    pairs.sort((x, y) => x[0] - y[0]);
    const done = new Uint8Array(tracks.length);
    for (const [, a, b] of pairs) {
      if (done[a] || used[b]) continue;
      done[a] = used[b] = 1;
      const tr = tracks[a], c = comps[b];
      if (T > 0) updateKinematics([wrapD(c.col - tr.col, BW) / T, wrapD(c.row - tr.row, BH) / T], T, tr.kin);
      tr.col = c.col; tr.row = c.row; tr.miss = 0;
      tr.hist.push(sampleOf(b, tr.kin), false);
    }
    tracks.forEach((tr, a) => { if (!done[a]) { tr.miss++; tr.hist.push(NO_SAMPLE, false); } });
    tracks = tracks.filter(tr => tr.miss < TRACK_MISS_MAX);
    comps.forEach((c, b) => {
      if (used[b]) return;
      const tr = { col: c.col, row: c.row, kin: newKinematics(), hist: StatsChart.ring(), miss: 0 };
      tr.hist.push(sampleOf(b, tr.kin), false);
      tracks.push(tr);
    });
  }

  // Velocity filters, all fed by v = the focus CoM's displacement / T since the last frame.
  //   Arrow: a short exponential average of v. Averaging a turning velocity makes it trail the
  //     true direction and shrink; the average's exact response to a steady turn at omega is
  //     divided out (rotation capped at LAG_MAX, gain at 2).
  //   omega: the turn rate of the arrow's average (already damped against pulses and detection
  //     jumps, and its angle changes add up exactly), averaged over OMEGA_TAU.
  //   Heading-up: v's direction averaged as an angle over VEL_TAU and predicted forward by omega
  //     each update, so a steady turn is followed with no lag. Each pull towards v is capped at
  //     JUMP_MAX, so a jump in the detected CoM barely moves the view. (Averaging v as a vector
  //     would not do: for a fast spinner it shrinks towards zero and its direction means nothing.)
  //   All averages start as running means: the first measurement counts in full, nothing ramps
  //   up from zero.
  function updateKinematics(v, T, k = kin) {   // k: the focus's, or another object's (every-object mode)
    k.age += T;
    const wrapA = a => a - 2 * Math.PI * Math.round(a / (2 * Math.PI));
    const ease = tau => Math.max(1 - Math.exp(-T / tau), T / k.age);   // EMA weight, running mean at first

    const a = ease(ARROW_TAU);
    k.fast = [k.fast[0] + a * (v[0] - k.fast[0]), k.fast[1] + a * (v[1] - k.fast[1])];

    // turn rate of the short average, while it is long enough to have a direction
    const fs = Math.hypot(k.fast[0], k.fast[1]);
    const fAng = fs > SPEED_MIN / 2 ? Math.atan2(k.fast[1], k.fast[0]) : null;
    if (fAng !== null && k.fastAng !== null) {
      k.omegaAge += T;
      const b = Math.max(1 - Math.exp(-T / OMEGA_TAU), T / k.omegaAge);
      k.omega += b * (wrapA(fAng - k.fastAng) / T - k.omega);
    } else k.omega *= Math.exp(-T / OMEGA_TAU);
    k.fastAng = fAng;

    // arrow: for v_k = e^(i omega T k) the average is v_k * a / (1 - (1 - a) e^(-i omega T));
    // divide that factor out
    const wT = k.omega * T;
    const mx = (1 - (1 - a) * Math.cos(wT)) / a, my = (1 - a) * Math.sin(wT) / a;
    const t = Math.max(-LAG_MAX, Math.min(LAG_MAX, Math.atan2(my, mx)));
    k.vel = rot(k.fast, t).map(x => x * Math.min(2, Math.hypot(mx, my)));
    k.speed = Math.hypot(k.vel[0], k.vel[1]);
    k.velAng = Math.atan2(k.vel[1], k.vel[0]);
    if (k.speed > SPEED_HI) k.moving = true;
    else if (k.speed < SPEED_MIN) k.moving = false;

    // heading-up: predict by the turn rate, then pull (capped) towards v's direction
    const ang = Math.hypot(v[0], v[1]) > SPEED_MIN / 2 ? Math.atan2(v[1], v[0]) : null;
    if (k.slowHeading === null) k.slowHeading = ang;
    else {
      k.slowHeading = wrapA(k.slowHeading + k.omega * T);
      if (ang !== null) {
        const r = Math.max(-JUMP_MAX, Math.min(JUMP_MAX, wrapA(ang - k.slowHeading)));
        k.slowHeading = wrapA(k.slowHeading + ease(VEL_TAU) * r);
      }
    }
  }

  function setFocus(i) {
    masking = null;
    md.dirty = true;   // re-detect the new focus's morph, even while paused
    if (tracking() && focus && i === focusIdx) return;   // the same creature again: keep its line
    stopRecord();   // a recording is of one creature: a new focus ends it (kept)
    stats.longGap = true;
    // every-object mode: the new focus takes over its track (short history and filters), the old
    // focus becomes a track. A track sits where it was last matched, i.e. on this scan's object.
    let taken = null;
    if (tracking() && i >= 0) {
      const c = found.comps[i];
      let best = 2;   // cells
      for (const tr of tracks) {
        const d = Math.hypot(wrapD(tr.col - c.col, BW), wrapD(tr.row - c.row, BH));
        if (d < best) { best = d; taken = tr; }
      }
    }
    if (tracking() && focus)
      tracks.push({ col: focus.col, row: focus.row, kin, hist: stats.short, miss: 0 });
    if (taken) {
      tracks = tracks.filter(tr => tr !== taken);
      stats.short = taken.hist; kin = taken.kin;
    } else {
      if (tracking() && focus) stats.short = StatsChart.ring();   // the old one went to its track
      resetKinematics();
      stats.shortGap = true;
    }
    focusIdx = i;
    if (i < 0) { focus = null; return; }
    const c = found.comps[i];
    focus = { col: c.col, row: c.row, mass: c.mass, radius: c.radius };
  }

  function frame() {
    let n = 0;
    if (running) {
      stepAcc += SPEEDS[speedIdx] * STEPS_PER_SPEED;
      n = Math.floor(stepAcc);
      stepAcc -= n;
      // analyse after every step, not once per frame, so the stats get one sample per step
      // (at most 2 steps per frame, at the top speed)
      for (let i = 0; i < n && running; i++) { advance(1); analyze(1); }
    }
    if (n === 0) analyze(0);
    if (detectMode !== 'none' && (n > 0 || md.dirty || stats.long.len !== md.len)) detectMorphs();
    const p = viewAngle(), c = Math.cos(p), s = Math.sin(p);
    SimGL.render({
      canvasW: CW, canvasH: CH, worldScale: WS,
      zoom: [ZOOM.x, ZOOM.y, ZOOM.edge, CELL_PX],
      rot: [c, s, -s, c],
      focus: focus ? [focus.col, focus.row, 1] : [0, 0, 0],
      stars: soundOn, time: (performance.now() / 1000) % 3600,   // star field follows the Sound button
    });
    drawOverlay();
    drawReadout();
    if (chart.relock && stats.short.len >= stats.short.cap) chart.relock = false;
    if (statsOn()) StatsChart.draw($('c-canvas'), stats.short, stats.long, chart, everyObject ? tracks.map(t => t.hist) : [], hovered()?.tr?.hist);
    requestAnimationFrame(frame);
  }

  // Morph detection: the result (detect.js classify) is kept on each history as .det, so it moves
  // with it when the focus changes. The focus is classified every frame, the other objects'
  // tracks DETECT_PER_FRAME per frame in turn.
  function detectMorphs() {
    md.dirty = false;
    if (stats.long.len !== md.len) { md.model = MorphDetect.build(stats.long); md.len = stats.long.len; }
    stats.short.det = focus ? MorphDetect.classify(md.model, stats.short) : null;
    if (detectMode !== 'every') return;
    for (let q = 0; q < Math.min(DETECT_PER_FRAME, tracks.length); q++) {
      const tr = tracks[md.next++ % tracks.length];
      tr.hist.det = MorphDetect.classify(md.model, tr.hist);
    }
  }
  // The detection result of object k of this scan: the focus's, or that of the track on it.
  function detOf(k, isF) {
    if (isF) return stats.short.det;
    const c = found.comps[k];
    for (const t of tracks) if (Math.hypot(wrapD(t.col - c.col, BW), wrapD(t.row - c.row, BH)) < 2) return t.hist.det;
    return null;
  }
  const hueFill = (h, a) => Number.isFinite(h) ? `hsla(${h.toFixed(0)},100%,55%,${a})` : `rgba(200,200,200,${a})`;

  // ---- overlay --------------------------------------------------------------------------
  function patchRect(i, bottom) {
    const total = 3 * PATCH.w + 2 * PATCH.gap;
    return { x: PX + (PW - total) / 2 + i * (PATCH.w + PATCH.gap), y: bottom ? PATCH.bottom : PATCH.top,
             w: PATCH.w, h: PATCH.h };
  }
  const zoomCenter = () => [ZOOM.x + ZOOM.edge / 2, ZOOM.y + ZOOM.edge / 2];

  function drawOverlay() {
    ctx.clearRect(0, 0, CW, CH);
    ctx.lineJoin = 'round';

    // panel separator
    ctx.fillStyle = '#1c1826';
    ctx.fillRect(PX, 0, 2, CH);

    // tool patches: top row the channel, bottom row the action (in the channel's colour)
    for (const bottom of [false, true]) {
      for (let i = 0; i < 3; i++) {
        const r = patchRect(i, bottom), sel = bottom ? i === action : i === channel;
        ctx.globalAlpha = sel ? 1 : 0.4;
        ctx.fillStyle = CH_COLORS[bottom ? channel : i];
        ctx.fillRect(r.x, r.y, r.w, r.h);
        ctx.globalAlpha = 1;
        if (sel) {
          ctx.strokeStyle = '#fff'; ctx.lineWidth = 4;
          ctx.strokeRect(r.x + 2, r.y + 2, r.w - 4, r.h - 4);
        }
        ctx.fillStyle = sel ? '#000' : 'rgba(0,0,0,.75)';
        ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.font = bottom ? '700 17px system-ui, sans-serif' : '700 22px system-ui, sans-serif';
        ctx.fillText(bottom ? ACTIONS[i] : CH_NAMES[i], r.x + r.w / 2, r.y + r.h / 2);
      }
    }

    const [cx, cy] = zoomCenter();

    // Tuning mode, stats mode with every object on, and Detect: Every: a circle around every
    // detected object in the world view; focus in yellow, the one under the mouse (every-object
    // mode) brighter. Detect: Every also draws each circle thick in its morph's mean hue (in place of yellow / white) and writes
    // the morph's name above it (neither if the morph is unknown).
    const detEvery = detectMode === 'every';
    if (mode === 'tuning' || (statsOn() && everyObject) || detEvery) {
      const hv = hovered();
      ctx.setLineDash([]);
      for (const [k, c] of found.comps.entries()) {
        const isF = focus && Math.hypot(wrapD(c.col - focus.col, BW), wrapD(c.row - focus.row, BH)) < 0.5;
        const isH = hv && hv.i === k;
        const x = (c.col + 0.5) * WS, y = (c.row + 0.5) * WS, r = Math.max(3, 2 * c.radius * WS);
        ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI);
        const det = detEvery ? detOf(k, isF) : null;
        ctx.strokeStyle = det && det.name ? hueFill(det.hue, 1)   // a known morph: thick, in its hue
          : isF ? 'rgba(255,212,0,.9)' : isH ? 'rgba(255,255,255,.95)' : 'rgba(255,255,255,.45)';
        ctx.lineWidth = det && det.name ? 3 : isH ? 2 : 1;
        ctx.stroke();
        if (det && det.name) {
          ctx.font = '600 11px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
          ctx.lineWidth = 3; ctx.strokeStyle = 'rgba(0,0,0,.7)'; ctx.strokeText(det.label, x, y - r - 2);
          ctx.fillStyle = '#fff'; ctx.fillText(det.label, x, y - r - 2);
        }
      }
    }

    // zoom frame around the focus in the world view (rotated with the view)
    if (focus) {
      const h = ZOOM_CELLS / 2, corners = [[-h, -h], [h, -h], [h, h], [-h, h]].map(v => {
        const b = rot(v);
        return [(focus.col + 0.5 + b[0]) * WS, (focus.row + 0.5 + b[1]) * WS];
      });
      ctx.save();
      ctx.beginPath(); ctx.rect(0, 0, PX, CH); ctx.clip();
      ctx.strokeStyle = 'rgba(255,255,255,.35)'; ctx.lineWidth = 1.5;
      ctx.beginPath(); corners.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]));
      ctx.closePath(); ctx.stroke();
      ctx.restore();
    }

    // tool outline, in both views (dashed = mask, solid = erase, dotted = draw)
    const off = toolOffset();
    if (off) {
      ctx.strokeStyle = CH_COLORS[channel]; ctx.lineWidth = 2; ctx.setLineDash(ACTION_DASH[action].map(x => x * 2 / 3));
      ctx.beginPath();
      ctx.arc((focus.col + 0.5 + off[0]) * WS, (focus.row + 0.5 + off[1]) * WS, Math.max(2, maskR * WS), 0, 2 * Math.PI);
      ctx.stroke();
    }
    // ship view: the tool circle follows the pointer while it hovers, filled while pressed
    const ring = masking ? masking.view : hover;
    if (ring && focus) {
      ctx.save();
      ctx.beginPath(); ctx.rect(ZOOM.x, ZOOM.y, ZOOM.edge, ZOOM.edge); ctx.clip();
      ctx.strokeStyle = CH_COLORS[channel]; ctx.lineWidth = 2.5; ctx.setLineDash(ACTION_DASH[action]);
      ctx.beginPath();
      ctx.arc(cx + ring[0] * CELL_PX, cy + ring[1] * CELL_PX, maskR * CELL_PX, 0, 2 * Math.PI);
      if (masking) { ctx.fillStyle = CH_FILLS[channel]; ctx.fill(); }
      ctx.stroke();
      ctx.restore();
    }

    // velocity arrow, from the zoom centre, in view space
    if (focus && kin.speed >= SPEED_MIN) {
      const v = unrot(kin.vel);
      let len = kin.speed * ARROW_T * CELL_PX;
      len = Math.min(len, ZOOM.edge / 2 - 14);
      const ux = v[0] / kin.speed, uy = v[1] / kin.speed;
      const ex = cx + ux * len, ey = cy + uy * len;
      ctx.strokeStyle = '#ffd400'; ctx.fillStyle = '#ffd400'; ctx.lineWidth = 4; ctx.lineCap = 'round';
      ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(ex - ux * 10, ey - uy * 10); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(ex, ey);
      ctx.lineTo(ex - ux * 18 - uy * 10, ey - uy * 18 + ux * 10);
      ctx.lineTo(ex - ux * 18 + uy * 10, ey - uy * 18 - ux * 10);
      ctx.closePath(); ctx.fill();
    }

    // grey round-corner border around the ship panel: the zoom view and both rows of tool patches
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(160,156,172,.6)'; ctx.lineWidth = 2;
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(PX + 3, 1, PW - 4, CH - 2, 10);
    else ctx.rect(PX + 3, 1, PW - 4, CH - 2);
    ctx.stroke();

    // small labels in the zoom view's top corners
    ctx.setLineDash([]);
    ctx.font = '12px ui-monospace, Menlo, Consolas, monospace';
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillStyle = 'rgba(200,192,214,.75)';
    ctx.fillText(VIEW_LABELS[viewMode].replace('View: ', ''), ZOOM.x + 8, ZOOM.y + 6);
    ctx.textAlign = 'right';
    ctx.fillText(`${ZOOM_CELLS}×${ZOOM_CELLS} cells`, ZOOM.x + ZOOM.edge - 8, ZOOM.y + 6);

    // Detect: the focus's morph in a box at the top of the ship view, under the channel patches
    if (detectMode !== 'none' && focus) {
      const det = stats.short.det, txt = det && det.name ? det.label : '???';
      ctx.font = '600 22px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      const w = ctx.measureText(txt).width + 24, bx = cx - w / 2, by = ZOOM.y + 8;
      ctx.fillStyle = 'rgba(40,40,46,.75)';
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(bx, by, w, 34, 7); else ctx.rect(bx, by, w, 34);
      ctx.fill();
      ctx.fillStyle = '#fff'; ctx.fillText(txt, cx, by + 17);
    }

    // no focus: say how to get one
    if (!focus && status === 'alive') {
      ctx.fillStyle = 'rgba(220,214,232,.85)'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = '15px system-ui, sans-serif';
      ctx.fillText('No creature selected.', cx, cy - 12);
      ctx.fillText('Click a creature on the left to follow it.', cx, cy + 12);
    }

    // end-of-run banner over the world
    if (status !== 'alive') {
      const msg = 'Dissolved';
      ctx.fillStyle = 'rgba(10,8,16,.85)';
      ctx.fillRect(PX / 2 - 170, CH / 2 - 42, 340, 84);
      ctx.fillStyle = '#fff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = '600 26px system-ui, sans-serif';
      ctx.fillText(msg, PX / 2, CH / 2 - 12);
      ctx.font = '14px system-ui, sans-serif'; ctx.fillStyle = '#bbb';
      ctx.fillText('Reset (Z) or Add (X) to try again', PX / 2, CH / 2 + 18);
    }
  }

  function drawReadout() {
    const n = found.comps.length;
    const st = status === 'alive' ? `<b class="ok">alive</b>` : `<b class="bad">${status}</b>`;
    // fixed-width numbers (the readout is monospace): pad with no-break spaces
    const pad = (v, n) => String(v).padStart(n, '\u00a0');
    let html = `<span>step <b>${pad(steps, 4)}</b></span><span>t <b>${pad((steps * dt).toFixed(1), 6)}</b></span>` +
      `<span>R <b>${R}</b></span><span>dt <b>${dt}</b></span>` +
      `<span>creatures <b>${pad(n, 3)}</b></span><span>${st}</span>`;
    if (focus) {
      const pct = spawnMass ? Math.round(100 * focus.mass / spawnMass) : 100;
      const head = pad(kin.speed >= SPEED_MIN ? (Math.round(deg(-kin.velAng)) % 360 + 360) % 360 : '–', 3) + '°';
      html += `<span>focus [mass <b>${pad(focus.mass.toFixed(0), 4)}</b> (${pad(pct, 3)}%) · ` +
        `speed <b>${kin.speed.toFixed(2)}</b> cells/t · heading <b>${head}</b>]</span>`;
    } else html += `<span>focus [<b>none</b>]</span>`;
    if (stats.recLeft > 0) html += `<span>rec <b class="bad">●</b> ${pad(REC_STEPS - stats.recLeft, 3)}/${REC_STEPS}` +
      ` ${stats.recFuzzy ? '~' : ''}${stats.recName}</span>`;
    html += `<span>${ACTIONS[action].toLowerCase()} [${CH_NAMES[channel]} <b>${masking ? 'on\u00a0' : 'off'}</b> r=${maskR}]</span>`;
    if (DEBUG) html += `<span>${msPerStep.toFixed(2)} ms/step · ${SimGL.tapCount} taps</span>`;
    // Tuning: the share of the board's mass inside detected objects (the rest is in background
    // blocks, under the threshold, or in basins under the min area)
    if (mode === 'tuning' && found.total > 0) {
      const inObj = found.comps.reduce((s, c) => s + c.mass, 0);
      html += `<span>in objects <b>${pad((100 * inObj / found.total).toFixed(1), 5)}%</b> of board mass</span>`;
    }
    readout.innerHTML = html;
    if (statsOn()) {
      const L = stats.long, pieces = L.segs.filter(s => s.gap).length;
      const txt = `recorded ${L.len} steps · ${pieces} piece${pieces === 1 ? '' : 's'}`;
      if ($('c-info').textContent !== txt) $('c-info').textContent = txt;
    }
  }

  // ---- input ----------------------------------------------------------------------------
  function canvasPoint(e) {
    const r = ui.getBoundingClientRect();
    return [(e.clientX - r.left) * CW / r.width, (e.clientY - r.top) * CH / r.height];
  }
  const inZoom = p => p[0] >= ZOOM.x && p[0] < ZOOM.x + ZOOM.edge && p[1] >= ZOOM.y && p[1] < ZOOM.y + ZOOM.edge;
  function viewOffset(p) {
    const [cx, cy] = zoomCenter(), h = ZOOM_CELLS / 2;
    const v = [(p[0] - cx) / CELL_PX, (p[1] - cy) / CELL_PX];
    return [Math.max(-h, Math.min(h, v[0])), Math.max(-h, Math.min(h, v[1]))];
  }

  ui.addEventListener('pointerdown', e => {
    const p = canvasPoint(e);
    if (p[0] < PX) {
      // world view: focus the creature under (or near) the click
      const i = Creatures.pick(found, p[0] / WS - 0.5, p[1] / WS - 0.5, CLICK_GATE, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
      if (i >= 0) setFocus(i);
      return;
    }
    for (const bottom of [false, true]) {
      for (let i = 0; i < 3; i++) {
        const r = patchRect(i, bottom);
        if (p[0] >= r.x && p[0] < r.x + r.w && p[1] >= r.y && p[1] < r.y + r.h) {
          if (bottom) action = i; else channel = i;
          return;
        }
      }
    }
    if (inZoom(p) && focus && !masking && e.button === 0) {
      masking = { id: e.pointerId, view: viewOffset(p) };
      ui.setPointerCapture(e.pointerId);
      e.preventDefault();
    }
  });
  // Tooltip over the tool patches (they are drawn on the canvas, so the canvas title changes).
  const CH_WORDS = ['red', 'green', 'blue'];
  function patchTip(p) {
    for (const bottom of [false, true]) for (let i = 0; i < 3; i++) {
      const r = patchRect(i, bottom);
      if (p[0] >= r.x && p[0] < r.x + r.w && p[1] >= r.y && p[1] < r.y + r.h)
        return bottom ? [
          `Mask (4): pressing on the ship view hides the ${CH_WORDS[channel]} channel there from the creature's senses`,
          `Erase (5): pressing on the ship view sets ${CH_WORDS[channel]} mass there to 0`,
          `Draw (6): pressing on the ship view sets ${CH_WORDS[channel]} mass there to 1`][i]
          : `${CH_NAMES[i]} (${i + 1}): act on the ${CH_WORDS[i]} channel`;
    }
    return '';
  }
  ui.addEventListener('pointermove', e => {
    const p = canvasPoint(e);
    const tip = patchTip(p);
    if (ui.title !== tip) ui.title = tip;
    if (masking && e.pointerId === masking.id) masking.view = viewOffset(p);
    hover = e.pointerType === 'mouse' && inZoom(p) ? viewOffset(p) : null;
    worldHover = e.pointerType === 'mouse' && p[0] < PX ? [p[0] / WS - 0.5, p[1] / WS - 0.5] : null;
  });
  ui.addEventListener('pointerleave', () => { hover = null; worldHover = null; });
  const endMask = e => { if (masking && e.pointerId === masking.id) masking = null; };
  ui.addEventListener('pointerup', endMask);
  ui.addEventListener('pointercancel', endMask);
  ui.addEventListener('lostpointercapture', endMask);
  ui.addEventListener('contextmenu', e => e.preventDefault());

  const $ = id => document.getElementById(id);
  // Briefly highlight the controls an action used (a click or its shortcut key).
  function flash(...ids) {
    for (const id of ids) {
      const el = $(id).tagName === 'INPUT' ? $(id).parentElement : $(id);   // a slider: its label
      el.classList.remove('flash'); void el.offsetWidth; el.classList.add('flash');
      clearTimeout(el._flash); el._flash = setTimeout(() => el.classList.remove('flash'), 180);
    }
  }
  function syncButtons() {
    $('play').textContent = running ? 'Pause' : 'Play';
    $('headup').textContent = viewMode === 'axis' ? 'Body axis up' : 'Heading up';
    $('headup').setAttribute('aria-pressed', viewMode !== 'board');
    const sp = SPEEDS[speedIdx];
    $('speed-v').textContent = sp < 1 ? '1/' + Math.round(1 / sp) : sp;
    $('radius-v').textContent = maskR;
    $('zoom-v').textContent = CELL_PX + '\u00d7';
    for (const el of document.querySelectorAll('input[type=range]')) fillRange(el);
  }
  // The filled part of a slider bar (CSS --pct), since the bar is drawn by hand for the theme.
  function fillRange(el) { el.style.setProperty('--pct', 100 * (el.value - el.min) / (el.max - el.min) + '%'); }
  document.addEventListener('input', e => { if (e.target.type === 'range') fillRange(e.target); });
  function togglePlay() { if (status === 'alive') wantRun = running = !running; syncButtons(); }
  // A reset keeps the play / pause state the user chose (wantRun), even if the sim stopped because
  // everything dissolved.
  let wantRun = true;
  function reset() { spawn(); running = wantRun; syncButtons(); }
  function changeSpeed(d) {
    speedIdx = Math.max(0, Math.min(SPEEDS.length - 1, speedIdx + d));
    $('speed').value = speedIdx + 1; syncButtons();
  }
  function stepOnce() { if (status !== 'alive') return; wantRun = running = false; advance(1); analyze(1); syncButtons(); }
  function toggleHeadUp() { viewMode = viewMode === 'board' ? 'vel' : 'board'; syncButtons(); }

  // Two pickers: species (a rule code, e.g. 5N7KKM) and morph (a form of it, e.g. limus). Each
  // legacy pattern has its own rule settings, so it counts as a species with one morph.
  const speciesOf = p => p.logic === 'legacy' ? p.id : p.rule;
  const SPECIES_NAMES = { '5N7KKM': 'Aquarium original', 'VT049W': 'Aquarium replicating', 'G6G6CR': 'Cilium',
    'Z18A9R': 'Aquarium spinning', 'KH97WU': 'Aquarium dancing', 'GDNQYX': 'Aquarium drifting',
    'HGMLH8': 'Aquarium shapeshifting', 'PJ1EF7': 'Aquarium multiplying', 'KZ2BLD': 'Aquarium striped',
    'XEH4YR': 'Aquarium exploding', 'HAESRE': 'Aquarium zigzagging' };
  const SPECIES_LAST = ['G6G6CR'];   // Cilium, the one non-Aquarium species, goes to the bottom
  const speciesName = p => p.logic === 'legacy' ? p.name : `${SPECIES_NAMES[p.rule] || p.rule} (${p.rule})`;
  const morphName = p => p.logic === 'legacy' ? 'rectus' : p.name.replace(p.rule + ' ', '');
  const BOOM = / \u{1F4A5}$/u;   // the self-replicating morphs' emoji: shown in the picker, not recorded
  const listed = () => {   // hidden: only via ?form=
    const l = PATTERNS.filter(p => !p.hidden || p === pattern);
    return l.filter(p => !SPECIES_LAST.includes(p.rule)).concat(l.filter(p => SPECIES_LAST.includes(p.rule)));
  };
  const morphs = () => listed().filter(p => speciesOf(p) === speciesOf(pattern));
  const spSel = $('species'), mSel = $('morph');
  function syncPicker() {
    spSel.textContent = '';
    for (const p of listed()) {
      if ([...spSel.options].some(o => o.value === speciesOf(p))) continue;
      spSel.add(new Option(speciesName(p), speciesOf(p)));
    }
    spSel.value = speciesOf(pattern);
    mSel.textContent = '';
    for (const p of morphs()) mSel.add(new Option(morphName(p), p.id));
    mSel.value = pattern.id;
  }
  // A new species changes the rule: reset the world with its first morph.
  spSel.onchange = () => {
    setPattern(listed().find(p => speciesOf(p) === spSel.value).id); reset(); spSel.blur();
    flash('species', 'reset');
  };
  // Previous / next species (Shift+C / Shift+V), wrapping round: as picking it in the list.
  function stepSpecies(d) {
    const n = spSel.options.length;
    spSel.selectedIndex = (spSel.selectedIndex + d + n) % n;
    spSel.onchange();
  }
  // A new morph of the same species only changes what Add (X) puts in; the world is kept.
  function selectMorph(id) { pattern = findPattern(id); mSel.value = pattern.id; }
  // Picking a morph in the list, even the current one, also adds it (like Add). A native select
  // fires no change for the option already chosen, so it is cleared while open and restored after.
  const clearMorph = () => { if (mSel.selectedIndex >= 0) mSel.selectedIndex = -1; };
  mSel.addEventListener('pointerdown', clearMorph);
  mSel.addEventListener('focus', clearMorph);
  const restoreMorph = () => { if (mSel.selectedIndex < 0) mSel.value = pattern.id; };
  mSel.addEventListener('blur', restoreMorph);
  mSel.addEventListener('keydown', e => { if (e.key === 'Escape') restoreMorph(); });
  mSel.onchange = () => { selectMorph(mSel.value); mSel.blur(); flash('morph', 'add'); addPattern(); };
  // Previous / next morph (C / V), within the species, wrapping round.
  function stepPattern(d) {
    const n = mSel.options.length, i = (mSel.selectedIndex + d + n) % n;
    selectMorph(mSel.options[i].value);
  }
  const prevMorph = () => { stepPattern(-1); flash('pat-prev'); };
  const nextMorph = () => { stepPattern(1); flash('pat-next'); };
  $('pat-prev').onclick = prevMorph;
  $('pat-next').onclick = nextMorph;

  // Each action with the controls it touches: Step while playing also pauses (Play/Pause).
  const doPlay = () => { flash('play'); togglePlay(); };
  const doStep = () => { flash(...(running && status === 'alive' ? ['play', 'step'] : ['step'])); stepOnce(); };
  const doReset = () => { flash('reset'); reset(); };
  const doAdd = () => { flash('add'); addPattern(); };
  const doView = () => { flash('headup'); toggleHeadUp(); };
  const doSpeed = d => { flash('speed'); changeSpeed(d); };
  $('play').onclick = doPlay;
  $('step').onclick = doStep;
  $('reset').onclick = doReset;
  $('add').onclick = doAdd;
  $('headup').onclick = doView;
  $('speed').value = speedIdx + 1;
  $('speed').oninput = e => { speedIdx = +e.target.value - 1; syncButtons(); };
  $('zoom').value = zoomIdx + 1;
  $('zoom').oninput = e => {
    zoomIdx = +e.target.value - 1;
    CELL_PX = ZOOMS[zoomIdx]; ZOOM_CELLS = ZOOM.edge / CELL_PX;
    syncButtons();
  };
  $('radius').value = maskR;
  $('radius').oninput = e => { maskR = +e.target.value; syncButtons(); };

  // Modes: tuning and stats each show their controls under the screen, above the help text. The URL follows, so a reload keeps it.
  function applyMode() {
    $('mode').textContent = 'Mode: ' + MODE_NAMES[mode];
    $('mode').setAttribute('aria-pressed', mode !== 'action');   // yellow when not Action
    $('tune').hidden = mode !== 'tuning';
    $('stats').hidden = !statsOn();
    $('c-canvas').hidden = !statsOn();
    $('c-canvas').classList.toggle('beside', mode === 'action+stats');   // left of the ship view, not over it
    $('c-move').hidden = mode !== 'action+stats';
    if (mode === 'stats') { masking = null; hover = null; }   // the chart covers the ship view
  }
  // action+stats: the arrow in the chart's top-right corner moves it to the left half of the world
  // view (then points right, to move it back)
  let chartFar = false;
  $('c-move').onclick = () => {
    chartFar = !chartFar;
    for (const id of ['c-canvas', 'c-move']) $(id).classList.toggle('far', chartFar);
    $('c-move').innerHTML = chartFar ? '&#9654;' : '&#9664;';
    $('c-move').title = chartFar ? 'Move the chart back, next to the ship view' : 'Move the chart to the left edge';
  };
  // Detect: None / Focus / Every (detectMorphs)
  $('detect').onclick = () => {
    flash('detect');
    detectMode = DETECT_MODES[(DETECT_MODES.indexOf(detectMode) + 1) % DETECT_MODES.length];
    $('detect').textContent = 'Detect: ' + DETECT_NAMES[detectMode];
    $('detect').setAttribute('aria-pressed', detectMode !== 'none');   // yellow when on
    if (!tracking()) tracks = [];
    md.dirty = true;
  };
  function nextMode(d = 1) { flash('mode'); mode = MODES[(MODES.indexOf(mode) + d + MODES.length) % MODES.length]; applyMode(); }
  $('mode').onclick = () => nextMode();
  applyMode();

  // Stats chart (statschart.js): 3D or 2D phase space, a stat per axis (stats: t mass chan both hue
  // speed omega). Drag the 3D view to rotate it.
  const chart = { dim: 3,   // 3D by default (not in the URL)
                  axes: ['mass', 'hue', 'speed'],   // defaults (not in the URL)
                  c: 'omega',   // colour axis
                  // Lock axes: keep every axis's min and max; on by default (relock: see setPattern, which
                  // also relocks on the first species's data)
                  lock: true, relock: false,
                  yaw: -2.5, pitch: 0.5 };   // z up, x to the front-left, y to the front-right
  const axisSel = ['c-x', 'c-y', 'c-z'].map($);
  function syncChart() {
    $('c-dim').textContent = chart.dim === 3 ? '3D' : '2D';
    $('c-every').textContent = everyObject ? 'Every object: ON' : 'Every object: OFF';
    $('c-every').setAttribute('aria-pressed', everyObject);
    $('c-z-l').hidden = chart.dim !== 3;
    $('c-canvas').style.touchAction = chart.dim === 3 ? 'none' : '';   // 3D: a drag rotates, not scrolls
  }
  const cSel = $('c-c');
  for (const k of StatsChart.C_OPTIONS) cSel.add(new Option(k === 'none' ? 'none' : StatsChart.STATS[k].menu || StatsChart.STATS[k].name, k));
  cSel.value = chart.c;
  cSel.onchange = () => { chart.c = cSel.value; cSel.blur(); syncChart(); };
  axisSel.forEach((sel, i) => {
    for (const k of StatsChart.OPTIONS) sel.add(new Option(StatsChart.STATS[k].menu || StatsChart.STATS[k].name, k));
    sel.value = chart.axes[i];
    sel.onchange = () => { chart.axes[i] = sel.value; sel.blur(); syncChart(); };
  });
  $('c-dim').onclick = () => { chart.dim = chart.dim === 3 ? 2 : 3; syncChart(); };
  $('c-lock').onclick = () => {
    flash('c-lock'); chart.lock = !chart.lock; chart.relock = false;
    $('c-lock').textContent = chart.lock ? 'Lock axes: ON' : 'Lock axes: OFF';
    $('c-lock').setAttribute('aria-pressed', chart.lock);
  };
  $('c-every').onclick = () => { flash('c-every'); everyObject = !everyObject; if (!tracking()) tracks = []; syncChart(); };
  // Record (S): records the next REC_STEPS steps of the focus into the long history, as a new
  // piece. Pressing it while recording stops early. It asks for the morph name to record under
  // (offered: the picker's morph). Shift+S records a fuzzy piece (detect.js), offering the last
  // fuzzy name recorded. A name need not be a morph in the picker (e.g. cytos).
  function syncRecord() {
    $('c-rec').textContent = stats.recLeft > 0 ? 'Stop' : 'Record';
    $('c-rec').setAttribute('aria-pressed', stats.recLeft > 0);
  }
  const doRecord = (fuzzy = false) => {
    flash('c-rec');
    if (stats.recLeft > 0) { stats.recLeft = 0; stats.longGap = true; }   // stop
    else {
      const last = fuzzy && [...stats.long.segs].reverse().find(sg => sg.fuzzy);
      const name = fuzzy
        ? prompt('Record a fuzzy morph (matched loosely, shown as ~name). Name:', last ? last.morph : '')
        : prompt('Record the focus as morph:', morphName(pattern).replace(BOOM, ''));
      if (!name || !name.trim()) return;   // cancelled
      stats.recName = name.trim(); stats.recFuzzy = fuzzy;
      stats.recLeft = REC_STEPS; stats.recStart = stats.long.len;
      stats.longGap = true;
      // the piece starts with the focus's short history gathered so far (since its last gap),
      // the rest of the REC_STEPS steps recorded as they come
      const h = stats.short;
      if (focus && h.len) {
        let j0 = h.len - 1;
        while (j0 > 0 && !h.brk(j0)) j0--;
        for (let j = j0; j < h.len; j++) {
          const s = Object.fromEntries(['t', 'mass', 'mR', 'mG', 'mB', 'speed', 'omega'].map(k => [k, h.raw(k, j)]));
          if (!Number.isFinite(s.mass)) continue;   // no focus at that step
          stats.long.push(s, stats.longGap, stats.recName, stats.recFuzzy);
          stats.longGap = false; stats.recLeft--;
        }
      }
    }
    syncRecord();
  };
  $('c-rec').onclick = e => doRecord(e.shiftKey);
  // A recording ends early (what it has recorded is kept) on a new focus or when the focus dies.
  function stopRecord() {
    if (stats.recLeft <= 0) return;
    stats.recLeft = 0; stats.longGap = true; syncRecord();
  }
  // Delete: removes the last piece of the long history; during a recording, the piece being
  // recorded (and stops it; nothing if it has no samples yet).
  $('c-clear').onclick = () => {
    flash('c-clear');
    if (stats.recLeft > 0) {
      stopRecord();
      if (stats.long.len > stats.recStart) stats.long.dropLastPiece();
      return;
    }
    stats.long.dropLastPiece();
  };
  // Download: the long history as <species>.json. Saved into assets/stats/ and committed, it is
  // loaded whenever that species is (loadSavedStats); new recordings add to it.
  $('c-save').onclick = () => {
    flash('c-save');
    const obj = { format: 1, species: stats.species, R, dt, ...stats.long.toJSON() };
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([StatsChart.statsFileText(obj)], { type: 'application/json' }));   // aligned, one sample a line
    a.download = `${stats.species}.json`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  syncRecord();
  // Step axis i (0 x, 1 y, 2 z, 3 colour) to the next (d = 1) or previous (d = -1) stat: a click on
  // its label, or keys D F G H (Shift: backwards).
  function cycleAxis(i, d) {
    const o = i === 3 ? StatsChart.C_OPTIONS : StatsChart.OPTIONS;
    const step = v => o[(o.indexOf(v) + d + o.length) % o.length];
    if (i === 3) { chart.c = cSel.value = step(chart.c); flash('c-c'); }
    else { chart.axes[i] = axisSel[i].value = step(chart.axes[i]); flash(axisSel[i].id); }
    syncChart();
  }
  // On the chart: drag rotates the 3D view; a click (moved < 5 px) on an axis label, or on an axis
  // row of the 3D legend, steps that axis to the next stat.
  let chartDrag = null;
  const cc = $('c-canvas');
  const chartXY = e => { const r = cc.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
  cc.addEventListener('pointerdown', e => {
    chartDrag = { id: e.pointerId, x: e.clientX, y: e.clientY, moved: 0 };
    cc.setPointerCapture(e.pointerId);
  });
  cc.addEventListener('pointermove', e => {
    if (!chartDrag || e.pointerId !== chartDrag.id) {   // hover: a hand over a label
      const c = StatsChart.hitAxis(...chartXY(e)) >= 0 ? 'pointer' : chart.dim === 3 ? 'grab' : '';
      if (cc.style.cursor !== c) cc.style.cursor = c;
      return;
    }
    const dx = e.clientX - chartDrag.x, dy = e.clientY - chartDrag.y;
    chartDrag.moved += Math.hypot(dx, dy);
    chartDrag.x = e.clientX; chartDrag.y = e.clientY;
    if (chart.dim !== 3) return;
    chart.yaw += dx * 0.01;
    chart.pitch = Math.max(-1.5, Math.min(1.5, chart.pitch + dy * 0.01));
  });
  cc.addEventListener('pointerup', e => {
    if (chartDrag && e.pointerId === chartDrag.id && chartDrag.moved < 5) {
      const i = StatsChart.hitAxis(...chartXY(e));
      if (i >= 0) cycleAxis(i, 1);
    }
    chartDrag = null;
  });
  for (const ev of ['pointercancel', 'lostpointercapture']) cc.addEventListener(ev, () => { chartDrag = null; });
  syncChart();

  // Tuning mode: detection sliders (also ?saddle= ?thr= ?area=).
  const tune = [['t-saddle', v => DETECT.saddle = v, DETECT.saddle],
                ['t-thr', v => DETECT.threshold = v, DETECT.threshold],
                ['t-area', v => DETECT.minArea = v * MIN_AREA_UNIT, DETECT.minArea / MIN_AREA_UNIT]];
  for (const [id, set, v0] of tune) {
    const el = $(id), out = $(id + '-v');
    el.value = v0; out.textContent = (+el.value).toFixed(2); fillRange(el);
    el.oninput = () => { set(+el.value); out.textContent = (+el.value).toFixed(2); analyze(0); };
  }

  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (e.code === 'Space') { e.preventDefault(); doStep(); }
    else if (k === 'Enter') { e.preventDefault(); doPlay(); }
    else if (k === 'z') doReset();
    else if (k === 'x') doAdd();
    else if (k === 'c') e.shiftKey ? stepSpecies(-1) : prevMorph();
    else if (k === 'v') e.shiftKey ? stepSpecies(1) : nextMorph();
    else if (k === 'q') doSpeed(1);
    else if (k === 'w') doSpeed(-1);
    else if (k === 'e') doView();
    else if (k === 'a') nextMode(e.shiftKey ? -1 : 1);
    else if (k === 's') doRecord(e.shiftKey);
    else if ('dfgh'.includes(k) && k.length === 1) cycleAxis('dfgh'.indexOf(k), e.shiftKey ? -1 : 1);
    else if (e.key >= '1' && e.key <= '3') channel = +e.key - 1;
    else if (e.key >= '4' && e.key <= '6') action = +e.key - 4;
  });

  // Background music, royalty free: a playlist looped with 10 s of silence between tracks. The
  // Stars button turns the music and the star field on or off together. Off by default; browsers
  // only let audio start after a user gesture, so it begins on the first click, tap or key.
  const TRACKS = ['assets/music/long-way-rexlambo.mp3',              // "Long Way" by Rexlambo
                  'assets/music/a-few-jumps-away-arthur-vyncke.mp3']; // "A Few Jumps Away" by Arthur Vyncke
  const TRACK_GAP_MS = 10000;
  const music = new Audio();
  music.volume = 0.5;
  let track = 0, gapTimer = null;
  music.src = TRACKS[track];
  music.addEventListener('ended', () => {
    track = (track + 1) % TRACKS.length;
    music.src = TRACKS[track];
    gapTimer = setTimeout(() => { gapTimer = null; if (soundOn) music.play().catch(() => {}); }, TRACK_GAP_MS);
  });
  let soundOn = false;   // off by default; the browser remembers the choice
  try { soundOn = localStorage.getItem('cognoship-sound') === 'on'; } catch (e) {}
  function syncSound() {
    $('sound').textContent = soundOn ? 'Stars: ON' : 'Stars: OFF';
    $('sound').setAttribute('aria-pressed', soundOn);
    if (!soundOn) music.pause();
    else if (!gapTimer) music.play().catch(() => {});   // in a gap between tracks, the timer starts the next
  }
  $('sound').onclick = () => {
    flash('sound'); soundOn = !soundOn; syncSound();
    try { localStorage.setItem('cognoship-sound', soundOn ? 'on' : 'off'); } catch (e) {}
  };
  const startMusic = () => { if (soundOn && music.paused && !gapTimer) music.play().catch(() => {}); };
  for (const ev of ['pointerdown', 'keydown']) document.addEventListener(ev, startMusic, { capture: true });   // retried until it plays
  syncSound();

  // theme and fullscreen
  $('theme').onclick = () => {
    const light = document.documentElement.getAttribute('data-theme') !== 'light';
    if (light) document.documentElement.setAttribute('data-theme', 'light');
    else document.documentElement.removeAttribute('data-theme');
    try { localStorage.setItem('cognoship-theme', light ? 'light' : 'dark'); } catch (e) {}
  };
  // In fullscreen the two control rows stay below the canvas; --rows is their height, so the
  // canvas takes what is left.
  function fitFs() {
    if (!document.documentElement.classList.contains('fs')) return;
    const rows = ['.controls', '.row2'].map(q => document.querySelector(q));
    const h = rows.reduce((a, el) => a + el.offsetHeight + parseFloat(getComputedStyle(el).marginTop), 0) + 8;
    document.documentElement.style.setProperty('--rows', h + 'px');
  }
  function setFs(on) { document.documentElement.classList.toggle('fs', on); fitFs(); }
  window.addEventListener('resize', fitFs);
  if (window.ResizeObserver) {   // the rows can re-wrap when the canvas width changes
    const ro = new ResizeObserver(fitFs);
    ['.controls', '.row2'].forEach(q => ro.observe(document.querySelector(q)));
  }
  $('fullscreen').onclick = async () => {
    setFs(true);
    try { await document.documentElement.requestFullscreen?.(); } catch (e) {}
    try { await screen.orientation?.lock?.('landscape'); } catch (e) {}
  };
  $('fs-exit').onclick = () => { setFs(false); if (document.fullscreenElement) document.exitFullscreen(); };
  document.addEventListener('fullscreenchange', () => { if (!document.fullscreenElement) setFs(false); });

  // ---- start ----------------------------------------------------------------------------
  resetKinematics();
  syncButtons();
  SimGL.init(glCanvas, { W: BW, H: BH, R, kernels: patternRule(pattern).kernels })
    .then(() => { setPattern(pattern.id); spawn(); requestAnimationFrame(frame); })
    .catch(err => {
      const box = $('error');
      box.hidden = false;
      box.textContent = 'Could not start: ' + err.message +
        (location.protocol === 'file:' ? ' (serve this folder over http, e.g. python3 -m http.server)' : '');
      console.error(err);
    });
})();
