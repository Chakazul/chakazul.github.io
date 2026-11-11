"use strict";
// ============================================================================================
//  Cognoship core, shared by the game (index.html, game.js) and the lab (lab.html, lab.js):
//  the controllable ship, several creatures at once with a clickable focus (see PLAN.md
//  sections 4 and 4.5), their stats histories and morph detection.
//
//  CPU side: layout, canvas input, finding and following creatures, the heading / view logic,
//  the overlay (tool patches, velocity arrow, tool outline). The simulation and the tools run on
//  the GPU (glsim.js). Each frame the CPU reads back a coarse block-mass map and labels its
//  connected components (creatures.js); one of them is the focus, which the zoom view, the
//  arrow and the tools follow.
//
//  The page (lab.js, game.js) owns its HTML controls, readout and keys. It reads and sets the
//  state through the Core object, sets Core.hooks to hear about changes and to draw after each
//  frame, Core.opt to choose what the overlay shows, then calls Core.start().
//
//  Coordinates. Board cells are (col, row), row 0 at the top of the screen. CoMs are in cell
//  index units (cell i's centre is at i). "View" vectors are the zoom panel's (right, down) in
//  cells; rot() maps them to board vectors, so the tool lands under the pointer in any view mode.
// ============================================================================================

const Core = (function () {
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
                   minArea: num('area', 1, 0, 10) * MIN_AREA_UNIT };   // mutable: the lab's tuning sliders
  // Morph detection (detect.js): none, the focus only (Focused: a name box at the top of the
  // ship view), or every object (Every, the default: also its circle in its morph's hue and its
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

  // What the page shows, asked each frame (the page replaces these):
  //   circlesAll: a circle around every object in the world view (as Detect: Every does)
  //   hoverTracks: the object under the mouse in the world view is highlighted (hovered())
  //   stars: the star field behind the board
  //   viewLabels: the ship view's corner labels (view mode, cells across)
  //   actionPatches: the bottom row of tool patches (Mask Erase Draw); without it the action
  //     stays as set (Core.action)
  //   morphItalic: detected morph names in italics (not "???")
  //   banner: the end-of-run banner's layout, { title, line, inShip }: font sizes in canvas px of
  //     "Dissolved" and of the line below it (Core.banner; 0: no line), and whether it is centred
  //     in the ship view (else the world view); null: the lab's (26, 14 in a 340 x 84 box over the world)
  const opt = { circlesAll: () => false, hoverTracks: () => false, stars: () => false,
                viewLabels: () => true, actionPatches: () => true, morphItalic: () => false,
                banner: () => null, morphLabel: det => det.label,
                // clickFocus: a click on the world view follows the creature there; noFocusText:
                // the ship view's lines when nothing is followed
                clickFocus: () => true,
                // over: show the end-of-run banner ("Dissolved"); default: nothing left alive
                over: () => status !== 'alive',
                // pauseWhenEmpty: pause when nothing is left alive (the lab; the game plays on)
                pauseWhenEmpty: () => true,
                // nameBox: show the focus's morph name box (detection runs either way)
                // autoFocus: after a spawn or a cleared world, follow the creature nearest the
                // board centre (the lab); off, only an explicit focus does (the game's put_player)
                autoFocus: () => true,
                // overFill: a colour for the ship panel's background while over() holds, or null
                overFill: () => null,
                // worldShade: 0..1, how far the world view is faded to black
                worldShade: () => 0,
                // showKeys: draw the tool patches' keys (while the page shows its key badges)
                showKeys: () => false,
                // nameDot: a colour dot (the morph's hue) before the name in the name box
                nameDot: () => false,
                // nameBoxLow: the name box at the bottom of the ship view, just above the strip of
                // controls (the game's Upward / Zoom); else at the top, under the channel patches (the lab)
                nameBoxLow: () => false,
                // zoomPress(off, e): a press on the ship view at board offset off (from the focus CoM);
                // true = taken by the page (no tool). drawGhost(ctx, at, centre, cellPx, phi, focus):
                // draw the page's own thing in the ship view (at: the hovering pointer, canvas px, or
                // null; centre: the view's centre = the focus CoM; phi: viewAngle), instead of the tool
                // circle when it returns true.
                zoomPress: () => false,
                keepHover: () => false,
                // patchErase: the selected channel's patch toggles erase (no action patches)
                patchErase: () => false,
                // frameColor: the world view's frame around the followed creature (null: as the panel border)
                frameColor: () => null,
                // emptyView: with nothing followed, the ship view shows the world's centre (unrotated)
                // instead of staying blank
                emptyView: () => false,
                // focusCircle: a yellow circle around the followed creature in the world view (as
                // Detect: Every draws it), inside its frame, also when the other circles are off
                focusCircle: () => false,
                // zoomDrag(off, e): the pointer of a press zoomPress took moves (off as zoomPress's)
                zoomDrag: () => {},
                drawGhost: () => false,
                nameBox: () => true,
                noFocusText: () => ['No creature selected.', 'Click a creature on the left to follow it.'] };
  // morphLabel(det): the name shown for a detected (known) morph, or null to show it as unknown
  // ("???", no hue): the game shows only the morphs the player has been told about.
  const shownLabel = det => det && det.name ? opt.morphLabel(det) : null;
  const ital = () => opt.morphItalic() ? 'italic ' : '';
  // The page's callbacks: change (buttons may need syncing), record (recording started or
  // stopped), pattern (setPattern ran), species (a new species: the histories were cleared),
  // statsLoaded (the species's saved long history loaded), frame (after drawing each frame).
  const hooks = { change() {}, record() {}, pattern() {}, species() {}, statsLoaded() {}, frame() {} };

  // ---- state ----------------------------------------------------------------------------
  // Speed: 7 slots. Slot i runs SPEEDS[i] * STEPS_PER_S steps per real second, by the clock, not per
  // frame, so the pace is the same at any refresh or frame rate: speed 1 is 30 steps a second (one
  // every 2 frames at 60 Hz, every 4 at 120 Hz, every frame at 30 fps; under LeniaF's rule the glider
  // moves ~0.65 cells a step). A frame runs at most MAX_STEPS_PER_FRAME steps, and counts at most
  // MAX_FRAME_S of clock time (a hiccup, a hidden tab): a device that cannot keep up plays slower,
  // smoothly, rather than piling up steps that make each frame slower still.
  const SPEEDS = [1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8];
  // A "game second": the steps speed 1 runs in a real second. The game counts its times (survive
  // time, discovery, fire cooldown, ...) in steps, so they stay matched to the creatures at any speed
  // and on a device that cannot keep up.
  const STEPS_PER_S = 30;
  const MAX_STEPS_PER_FRAME = 4, MAX_FRAME_S = 0.1;
  let simSteps = 0;     // every step run since the page loaded (steps resets with the world; this never does)
  let driftSec = 0;     // game seconds the dots drift this frame (frame(): driftTick)
  let starSec = 0;      // real seconds the star field drifts this frame
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
  // A reset keeps the play / pause state the user chose (wantRun), even if the sim stopped because
  // everything dissolved.
  let wantRun = true;

  const glCanvas = document.getElementById('gl');
  const ui = document.getElementById('ui');
  const ctx = ui.getContext('2d');

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

  // ---- species and morphs ---------------------------------------------------------------
  // Two levels: species (a rule code, e.g. 5N7KKM) and morph (a form of it, e.g. limus). Each
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
  // A new morph of the same species only changes what Add puts in; the world is kept.
  function selectMorph(id) { pattern = findPattern(id); }

  // ---- simulation loop ------------------------------------------------------------------
  let ruleLoaded = false;   // a species's rule is set (setPattern); until then nothing steps or is scanned
  function setPattern(id) {
    ruleLoaded = true;
    pattern = findPattern(id);
    const rule = patternRule(pattern);
    dt = DT_OVERRIDE ?? rule.dt;
    SimGL.setRule(rule.kernels, R, rule.kn);
    SimGL.setUpdate(rule.norm, dt, rule.gn);
    // merge saddle per species (aquarium233.js rules[code].detect), unless ?saddle= is given
    if (stats.species !== speciesOf(pattern)) {
      stats.recLeft = 0; hooks.record();   // a recording in progress is discarded with the long history
      stats.short.clear(); stats.long.clear(); tracks = []; stats.species = speciesOf(pattern);
      hooks.species();
      loadSavedStats(stats.species);
    }
    if (!q.has('saddle')) DETECT.saddle = AQUARIUM233.rules[pattern.rule]?.detect?.saddle ?? 0.9;
    hooks.pattern();
  }

  function spawn() {
    const f = Rule.formCells(pattern.cells, pattern.R, R);
    SimGL.seed(f, Math.floor(BH / 2 - f.h / 2), Math.floor(BW / 2 - f.w / 2));
    steps = 0; spawnMass = null; status = 'alive'; stepAcc = 0;
    focus = null; autoFocus = true;
    resetKinematics(); tracks = []; stats.short.clear();   // a new world: no short histories
    analyze(0);
  }

  // Put the pattern into the world (keeping what is there), by `how`:
  //   'center': unrotated, at the board centre (as a reset, but without clearing);
  //   'random': rotated by a random angle, at the emptiest of 24 random candidates (least mass
  //     within the form's half-size + R of its centre), as the lab's Add;
  //   'far': as 'random', but of the empty candidates (mass within EMPTY_TOL of the least) the one
  //     farthest from the focus (as 'random' if there is no focus).
  // focusNew: true = follow the new creature (the lab's Add); false = only if nothing is followed;
  // 'never' = never (the game's [min_morph]: a lost ship stays lost until Respawn).
  const EMPTY_TOL = 1;   // cell mass: "empty" is at most this much more than the emptiest candidate
  function putPattern(how = 'random', focusNew = true) {
    let f = Rule.formCells(pattern.cells, pattern.R, R), best;
    if (how === 'center') best = [BW / 2, BH / 2];
    else {
      f = Rule.rotate(f, Math.random() * 2 * Math.PI);
      const reach = Math.max(f.w, f.h) / 2 + R, B = SimGL.BLOCK, blk = SimGL.readBlocks();
      const cands = [];
      for (let k = 0; k < 24; k++) {
        const col = Math.random() * BW, row = Math.random() * BH;
        let m = 0;
        for (let by = 0; by < SimGL.BH; by++)
          for (let bx = 0; bx < SimGL.BW; bx++) {
            const d = Math.hypot(wrapD(bx * B + B / 2 - col, BW), wrapD(by * B + B / 2 - row, BH));
            if (d < reach) m += blk[(by * SimGL.BW + bx) * 4];
          }
        cands.push({ col, row, m });
      }
      let pick = cands.reduce((a, c) => c.m < a.m ? c : a);   // the emptiest (the first of equals)
      if (how === 'far' && focus) {
        const dist = c => Math.hypot(wrapD(c.col - focus.col, BW), wrapD(c.row - focus.row, BH));
        pick = cands.filter(c => c.m <= pick.m + EMPTY_TOL).reduce((a, c) => dist(c) > dist(a) ? c : a);
      }
      best = [pick.col, pick.row];
    }
    SimGL.add(f, Math.floor(best[1] - f.h / 2), Math.floor(best[0] - f.w / 2));
    if (status !== 'alive') { status = 'alive'; running = true; }
    analyze(0);
    // follow the new creature
    if (focusNew === true || (focusNew !== 'never' && !focus)) {
      const i = Creatures.pick(found, best[0], best[1], 2 * R, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
      if (i >= 0) setFocus(i);
    }
    hooks.change();
  }
  const addPattern = () => putPattern('random', true);   // the lab's Add

  // Empty the world: no creatures, no focus, step 0. Status stays 'alive' until a frame finds the
  // board empty (so a putPattern right after does not flash "Dissolved"); the next creature put
  // in is followed.
  function clearWorld() {
    SimGL.seed({ w: 0, h: 0, cells: [[], [], []] }, 0, 0);
    stopRecord();
    steps = 0; spawnMass = null; status = 'alive'; stepAcc = 0;
    focus = null; focusIdx = -1; autoFocus = true; masking = null;
    resetKinematics(); tracks = []; stats.short.clear();
    hooks.change();
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
    steps += n; simSteps += n;
    if (DEBUG) { SimGL.finish(); msPerStep = 0.9 * msPerStep + 0.1 * (performance.now() - t0) / n; }
  }

  // Find the creatures, follow the focus, update velocity and heading. n = steps since last call.
  function analyze(n) {
    const blk = SimGL.readBlocks();
    found = Creatures.find(blk, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH, DETECT);

    // No "Exploded" check: a board full of copies just keeps running.
    if (found.comps.length === 0) { status = 'dissolved'; if (opt.pauseWhenEmpty()) running = false; hooks.change(); }

    const T = n * dt;
    if (autoFocus && !opt.autoFocus()) autoFocus = false;
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

  // Saved long history of a species: assets/stats/<species>.json (the lab's Download). Loaded into the
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
        hooks.statsLoaded();
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
  // the lab's Record REC_STEPS steps at a time, each sample with the morph selected then.
  // A change of focus (spawn, reset, click, add) starts a new piece in both, and each Record in
  // the long one: a gap in the line. A change of species clears both.
  const stats = { short: StatsChart.ring(), long: StatsChart.log(), species: null,
                  recLeft: 0, shortGap: false, longGap: false,   // recLeft: steps still to record
                  recName: '', recFuzzy: false };   // the morph name being recorded; fuzzy (detect.js)
  const NO_SAMPLE = { t: NaN, mass: NaN, mR: NaN, mG: NaN, mB: NaN, speed: NaN, omega: NaN };
  const REC_STEPS = 300;   // one Record records this many steps of the focus
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
      if (--stats.recLeft === 0) hooks.record();
    }
  }
  // A recording ends early (what it has recorded is kept) on a new focus or when the focus dies.
  function stopRecord() {
    if (stats.recLeft <= 0) return;
    stats.recLeft = 0; stats.longGap = true; hooks.record();
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
  // The masses of all objects of this scan summed, per channel [R, G, B], in g (as the stats);
  // null before the first scan.
  function worldMass() {
    if (!found || !found.comps) return null;
    const ch = channelMasses(), m = 1 / (R * R), out = [0, 0, 0];
    for (let i = 0; i < ch.length; i++) out[i % 3] += ch[i] * m;
    return out;
  }
  // The stats sample of object i (this scan) with kinematics k, in R units.
  function sampleOf(i, k) {
    const ch = channelMasses(), m = 1 / (R * R);
    return { t: steps * dt, mass: found.comps[i].mass * m, mR: ch[3 * i] * m, mG: ch[3 * i + 1] * m, mB: ch[3 * i + 2] * m,
             speed: k.speed / R, omega: Math.abs(k.omega) };
  }

  // Every-object mode (everyObject, the lab's stats panel button): every object besides the focus keeps a short
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
  // The object under the mouse in the world view (not the focus), while opt.hoverTracks(), and its
  // track: { i, tr } (tr null if it has none yet), or null.
  function hovered() {
    if (!opt.hoverTracks() || !worldHover) return null;
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
  function setEveryObject(on) { everyObject = on; if (!tracking()) tracks = []; }

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

  // The star field's drift (setStars): angle (degrees, clockwise from left: 0 drifts left, 90 up),
  // speed and brightness (1 = the default). The drift is integrated here, so changing the speed or
  // angle moves the stars on smoothly from where they are. off wraps every 10 boards: the layers'
  // speeds (draw.glsl, 0.6 / 1.5 / 3.5) times 10 are whole numbers, so the wrap is invisible.
  // The stars go by real time (starSec), whatever the speed or pause. The dots (dots.js) go by
  // game time (driftSec): while playing, the simulation's own rate (the speed times the frame's
  // clock time, smooth, not step by step; slowed with the steps when a frame hits its cap), so they
  // keep pace with the creatures; while paused (or no species yet), at the speed-1 rate in real time,
  // so they still drift.
  const star = { angle: 0, speed: 1, bright: 1, off: [0, 0], t: 0, last: null };
  const SWAY_PERIOD = 2 * Math.PI / 0.3;   // s (draw.glsl STAR_SWAY_W)
  function driftTick() {
    const now = performance.now() / 1000, real = star.last === null ? 0 : Math.min(0.25, now - star.last);
    star.last = now;
    starSec = real;
    driftSec = real;   // while playing, frame() sets it from the steps' rate
    star.t = (star.t + starSec) % (100 * SWAY_PERIOD);
  }
  function setStars(angle = 0, speed = 1, bright = 1) { Object.assign(star, { angle, speed, bright }); }
  // the drift's unit direction (board x right, y down) and speed factor (the dots follow it too)
  function starDrift() { const a = star.angle * Math.PI / 180; return { dir: [-Math.cos(a), -Math.sin(a)], speed: star.speed }; }
  function starView() {
    const a = star.angle * Math.PI / 180, dir = [-Math.cos(a), -Math.sin(a)];
    star.off[0] = (star.off[0] + dir[0] * star.speed * starSec) % (10 * BW);
    star.off[1] = (star.off[1] + dir[1] * star.speed * starSec) % (10 * BH);
    return { stars: opt.stars(), time: star.t, starOff: star.off, starDir: dir, starSpeed: star.speed, starBright: star.bright };
  }

  function frame() {
    let n = 0;
    driftTick();
    if (!ruleLoaded) {   // started empty, no species yet (the game's script loads one)
      SimGL.render({ canvasW: CW, canvasH: CH, worldScale: WS, zoom: [ZOOM.x, ZOOM.y, ZOOM.edge, CELL_PX],
                     rot: [1, 0, 0, 1], focus: [0, 0, 0], ...starView() });
      drawOverlay(); hooks.frame(); requestAnimationFrame(frame);
      return;
    }
    if (running) {
      const sec = Math.min(MAX_FRAME_S, starSec), want = SPEEDS[speedIdx] * STEPS_PER_S * sec;
      stepAcc += want;
      const capped = stepAcc >= MAX_STEPS_PER_FRAME + 1;
      n = Math.min(MAX_STEPS_PER_FRAME, Math.floor(stepAcc));
      stepAcc = capped ? 0 : stepAcc - n;   // what the cap held back is dropped, not owed
      driftSec = (capped ? n : want) / STEPS_PER_S;   // game seconds: the dots keep pace with the steps run
      // analyse after every step, not once per frame, so the stats get one sample per step
      // (at most MAX_STEPS_PER_FRAME per frame)
      for (let i = 0; i < n && running; i++) { advance(1); analyze(1); }
    } else if (masking && focus && action !== 0) {
      // paused: erase / draw still act, at once (the mask has nothing to act on without a step)
      const off = toolOffset();
      SimGL.paint({ col: focus.col + off[0], row: focus.row + off[1], r: maskR, ch: channel, action });
    }
    if (n === 0) analyze(0);
    if (detectMode !== 'none' && (n > 0 || md.dirty || stats.long.len !== md.len)) detectMorphs();
    const p = focus ? viewAngle() : 0, c = Math.cos(p), s = Math.sin(p);
    SimGL.render({
      canvasW: CW, canvasH: CH, worldScale: WS,
      zoom: [ZOOM.x, ZOOM.y, ZOOM.edge, CELL_PX],
      rot: [c, s, -s, c],
      focus: focus ? [focus.col, focus.row, 1] : opt.emptyView() ? [BW / 2 - 0.5, BH / 2 - 0.5, 1] : [0, 0, 0],
      ...starView(),   // the star field (opt.stars: the lab's Stars button)
    });
    drawOverlay();
    hooks.frame();
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
  function setDetectMode(m) { detectMode = m; if (!tracking()) tracks = []; md.dirty = true; }
  // The detection result of object k of this scan: the focus's, or that of the track on it.
  function detOf(k, isF) {
    if (isF) return stats.short.det;
    const c = found.comps[k];
    for (const t of tracks) if (Math.hypot(wrapD(t.col - c.col, BW), wrapD(t.row - c.row, BH)) < 2) return t.hist.det;
    return null;
  }
  const PANEL_BORDER = '#605e67';   // the ship panel's border, and the world view's frame of it: opaque (its old
  // 60% grey over black), so it looks the same over anything and however often it is drawn
  function panelBorder() {
    ctx.setLineDash([]);
    ctx.strokeStyle = PANEL_BORDER; ctx.lineWidth = 2;
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(PX + 3, 1, PW - 4, CH - 2, 10);
    else ctx.rect(PX + 3, 1, PW - 4, CH - 2);
    ctx.stroke();
  }
  const hueFill = (h, a) => Number.isFinite(h) ? `hsla(${h.toFixed(0)},100%,55%,${a})` : `rgba(200,200,200,${a})`;

  // ---- overlay --------------------------------------------------------------------------
  function patchRect(i, bottom) {
    const total = 3 * PATCH.w + 2 * PATCH.gap;
    return { x: PX + (PW - total) / 2 + i * (PATCH.w + PATCH.gap), y: bottom ? PATCH.bottom : PATCH.top,
             w: PATCH.w, h: PATCH.h };
  }
  const patchRows = () => opt.actionPatches() ? [false, true] : [false];   // bottom = the action row
  const zoomCenter = () => [ZOOM.x + ZOOM.edge / 2, ZOOM.y + ZOOM.edge / 2];

  function drawOverlay() {
    ctx.clearRect(0, 0, CW, CH);
    ctx.lineJoin = 'round';

    // the ship panel's background while the run is over (opt.overFill, e.g. grey), under the rest
    const overFill = opt.over() && opt.overFill();
    if (overFill) { ctx.fillStyle = overFill; ctx.fillRect(PX, 0, PW, CH); }   // (faded with the rest: worldShade)

    // panel separator
    ctx.fillStyle = '#1c1826';
    ctx.fillRect(PX, 0, 2, CH);

    // tool patches: top row the channel, bottom row the action (in the channel's colour)
    for (const bottom of patchRows()) {
      for (let i = 0; i < 3; i++) {
        const r = patchRect(i, bottom), sel = bottom ? i === action : i === channel;
        // erase on the selected channel patch (opt.patchErase, the game's recycle): striped, "-G"
        const striped = !bottom && sel && action === 1 && opt.patchErase();
        ctx.globalAlpha = sel ? 1 : 0.4;
        ctx.fillStyle = CH_COLORS[bottom ? channel : i];
        if (striped) {
          ctx.save();
          ctx.beginPath(); ctx.rect(r.x, r.y, r.w, r.h); ctx.clip();
          ctx.fillStyle = '#111'; ctx.fillRect(r.x, r.y, r.w, r.h);
          ctx.strokeStyle = CH_COLORS[i]; ctx.lineWidth = 7;
          for (let d = -r.h; d < r.w + r.h; d += 16) { ctx.beginPath(); ctx.moveTo(r.x + d, r.y + r.h); ctx.lineTo(r.x + d + r.h, r.y); ctx.stroke(); }
          ctx.restore();
        } else ctx.fillRect(r.x, r.y, r.w, r.h);
        ctx.globalAlpha = 1;
        if (sel) {
          ctx.strokeStyle = '#fff'; ctx.lineWidth = 4;
          ctx.strokeRect(r.x + 2, r.y + 2, r.w - 4, r.h - 4);
        }
        ctx.fillStyle = sel ? '#000' : 'rgba(0,0,0,.75)';
        ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        ctx.font = bottom ? '700 17px system-ui, sans-serif' : '700 22px system-ui, sans-serif';
        const label = bottom ? ACTIONS[i] : (striped ? '\u2212' : '') + CH_NAMES[i];
        if (striped) {   // on the stripes: dark text in a light halo
          ctx.lineWidth = 5; ctx.strokeStyle = 'rgba(255,255,255,.85)'; ctx.strokeText(label, r.x + r.w / 2, r.y + r.h / 2);
        }
        ctx.fillText(label, r.x + r.w / 2, r.y + r.h / 2);
        if (!bottom && opt.showKeys()) {   // its key (1-3), as the page's key badges (style.css .kbd: 10 CSS px)
          const k = CW / (ui.getBoundingClientRect().width || CW);   // canvas px per CSS px
          ctx.font = `700 ${10 * k}px system-ui, sans-serif`;
          const t = String(i + 1), w = ctx.measureText(t).width + 8 * k, h = 14 * k;
          const bx = r.x + r.w - w - 2 * k, by = r.y + 2 * k;
          ctx.fillStyle = 'rgba(0,0,0,.78)';
          ctx.beginPath(); if (ctx.roundRect) ctx.roundRect(bx, by, w, h, 4 * k); else ctx.rect(bx, by, w, h); ctx.fill();
          ctx.fillStyle = '#fff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          ctx.fillText(t, bx + w / 2, by + h / 2);
        }
      }
    }

    const [cx, cy] = zoomCenter();

    // opt.circlesAll() (the lab's tuning mode, stats mode with every object on) and Detect: Every:
    // a circle around every detected object in the world view; focus in yellow, the one under the
    // mouse (hovered()) brighter. Detect: Every also draws each circle thick in its morph's mean
    // hue (in place of yellow / white) and writes the morph's name above it (neither if the morph is unknown).
    const detEvery = detectMode === 'every';
    if (opt.circlesAll() || detEvery) {
      const hv = hovered();
      ctx.setLineDash([]);
      for (const [k, c] of found.comps.entries()) {
        const isF = focus && Math.hypot(wrapD(c.col - focus.col, BW), wrapD(c.row - focus.row, BH)) < 0.5;
        const isH = hv && hv.i === k;
        const x = (c.col + 0.5) * WS, y = (c.row + 0.5) * WS, r = Math.max(3, 2 * c.radius * WS);
        ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI);
        const det = detEvery ? detOf(k, isF) : null, lab = shownLabel(det);
        ctx.strokeStyle = lab ? hueFill(det.hue, 1)   // a known morph: thick, in its hue
          : isF ? 'rgba(255,212,0,.9)' : isH ? 'rgba(255,255,255,.95)' : 'rgba(255,255,255,.45)';
        ctx.lineWidth = lab ? 3 : isH ? 2 : 1;
        ctx.stroke();
        if (lab) {
          ctx.font = ital() + '600 11px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
          ctx.lineWidth = 3; ctx.strokeStyle = 'rgba(0,0,0,.7)'; ctx.strokeText(lab, x, y - r - 2);
          ctx.fillStyle = '#fff'; ctx.fillText(lab, x, y - r - 2);
        }
      }
    }

    // zoom frame around the focus in the world view (rotated with the view), as the ship panel's
    // border; across the board's edge it wraps round: a copy on each side it crosses (up to four
    // near a corner)
    if (focus) {
      const h = ZOOM_CELLS / 2, corners = [[-h, -h], [h, -h], [h, h], [-h, h]].map(v => {
        const b = rot(v);
        return [(focus.col + 0.5 + b[0]) * WS, (focus.row + 0.5 + b[1]) * WS];
      });
      ctx.save();
      ctx.beginPath(); ctx.rect(0, 0, PX, CH); ctx.clip();
      ctx.setLineDash([]);
      ctx.strokeStyle = opt.frameColor() || PANEL_BORDER; ctx.lineWidth = 2;
      for (const dx of [-BW * WS, 0, BW * WS])
        for (const dy of [-BH * WS, 0, BH * WS]) {
          ctx.beginPath(); corners.forEach((p, i) => i ? ctx.lineTo(p[0] + dx, p[1] + dy) : ctx.moveTo(p[0] + dx, p[1] + dy));
          ctx.closePath(); ctx.stroke();
        }
      // opt.focusCircle: the followed creature's circle (unless the circles above drew it already)
      if (opt.focusCircle() && !(opt.circlesAll() || detEvery)) {
        const x = (focus.col + 0.5) * WS, y = (focus.row + 0.5) * WS, r = Math.max(3, 2 * focus.radius * WS);
        ctx.strokeStyle = 'rgba(255,212,0,.9)'; ctx.lineWidth = 1.5;
        for (const dx of [-BW * WS, 0, BW * WS])
          for (const dy of [-BH * WS, 0, BH * WS]) {
            ctx.beginPath(); ctx.arc(x + dx, y + dy, r, 0, 2 * Math.PI); ctx.stroke();
          }
      }
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
    let ghost = false;
    if (focus && !masking) {
      ctx.save();
      ctx.beginPath(); ctx.rect(ZOOM.x, ZOOM.y, ZOOM.edge, ZOOM.edge); ctx.clip();
      ghost = opt.drawGhost(ctx, ring ? [cx + ring[0] * CELL_PX, cy + ring[1] * CELL_PX] : null,
                            [cx, cy], CELL_PX, viewAngle(), focus);
      ctx.restore();
    }
    if (ring && focus && !ghost) {
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
    panelBorder();

    // small labels in the zoom view's top corners
    ctx.setLineDash([]);
    if (opt.viewLabels()) {
      ctx.font = '12px ui-monospace, Menlo, Consolas, monospace';
      ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      ctx.fillStyle = 'rgba(200,192,214,.75)';
      ctx.fillText(VIEW_LABELS[viewMode].replace('View: ', ''), ZOOM.x + 8, ZOOM.y + 6);
      ctx.textAlign = 'right';
      ctx.fillText(`${ZOOM_CELLS}×${ZOOM_CELLS} cells`, ZOOM.x + ZOOM.edge - 8, ZOOM.y + 6);
    }

    // Detect: the focus's morph in a box at the top of the ship view, under the channel patches (or
    // at its bottom, above the controls' strip: opt.nameBoxLow)
    if (detectMode !== 'none' && focus && opt.nameBox()) {
      const lab = shownLabel(stats.short.det), txt = lab || '???';
      ctx.font = (lab ? ital() : '') + '600 22px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      // a known morph's colour dot before its name (opt.nameDot; as the game's morph book)
      const det = stats.short.det, dot = lab && opt.nameDot(), D = 13, G = 10;
      const tw = ctx.measureText(txt).width, w = tw + 24 + (dot ? D + G : 0), bx = cx - w / 2, by = opt.nameBoxLow() ? PATCH.bottom - 34 - 10 : ZOOM.y + 16;
      ctx.fillStyle = 'rgba(40,40,46,.75)';
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(bx, by, w, 34, 7); else ctx.rect(bx, by, w, 34);
      ctx.fill();
      if (dot) {
        ctx.fillStyle = Number.isFinite(det.hue) ? `hsl(${det.hue.toFixed(0)},100%,55%)` : '#ccc';
        ctx.beginPath(); ctx.arc(bx + 12 + D / 2, by + 17, D / 2, 0, 2 * Math.PI); ctx.fill();
      }
      ctx.fillStyle = '#fff'; ctx.fillText(txt, bx + w - 12 - tw / 2, by + 17);
    }

    // no focus: say how to get one
    if (!focus && !opt.over()) {
      ctx.fillStyle = 'rgba(220,214,232,.85)'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = '15px system-ui, sans-serif';
      const lines = opt.noFocusText();
      lines.forEach((t, k) => ctx.fillText(t, cx, cy + 24 * (k - (lines.length - 1) / 2)));
    }

    // end-of-run banner over the world
    const over = opt.over(), bn = over && opt.banner();
    if (bn) {   // laid out by the page: the box fits the text
      const t = bn.title, s = bn.line || 0, msg = 'Dissolved';
      const x = bn.inShip ? PX + PW / 2 : PX / 2, maxW = bn.inShip ? PW - 24 : PX - 24;
      ctx.font = `600 ${t}px system-ui, sans-serif`;
      let w = ctx.measureText(msg).width + 1.2 * t;
      if (s) { ctx.font = `${s}px system-ui, sans-serif`; w = Math.max(w, ctx.measureText(core.banner).width + 2 * s); }
      w = Math.min(maxW, w);
      const h = 1.8 * t + (s ? 1.6 * s : 0);   // padding 0.4 t, title, [gap 0.6 s, line,] padding 0.4 t
      const top = CH / 2 - h / 2 + (bn.dy || 0);   // dy: the page may move it up or down
      ctx.fillStyle = 'rgba(10,8,16,.85)';
      ctx.beginPath();
      if (ctx.roundRect) ctx.roundRect(x - w / 2, top, w, h, 10); else ctx.rect(x - w / 2, top, w, h);
      ctx.fill();
      ctx.fillStyle = '#fff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = `600 ${t}px system-ui, sans-serif`;
      ctx.fillText(msg, x, top + 0.4 * t + t / 2);
      if (s) {
        ctx.font = `${s}px system-ui, sans-serif`; ctx.fillStyle = '#bbb';
        ctx.fillText(core.banner, x, top + 0.4 * t + t + 0.6 * s + s / 2);
      }
    } else if (over) {
      const msg = 'Dissolved';
      ctx.fillStyle = 'rgba(10,8,16,.85)';
      ctx.fillRect(PX / 2 - 170, CH / 2 - 42, 340, 84);
      ctx.fillStyle = '#fff'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.font = '600 26px system-ui, sans-serif';
      ctx.fillText(msg, PX / 2, CH / 2 - 12);
      ctx.font = '14px system-ui, sans-serif'; ctx.fillStyle = '#bbb';
      ctx.fillText(core.banner, PX / 2, CH / 2 + 18);
    }

    // the world view and the ship view's square shaded black by opt.worldShade() (0 clear ..
    // 1 hidden): everything drawn in them, "Dissolved" included (the game's scene and chapter fades);
    // the tool patches above and below the square, and the page's buttons, stay
    const shade = opt.worldShade();
    if (shade > 0) {
      ctx.fillStyle = `rgba(0,0,0,${Math.min(1, shade)})`;
      ctx.fillRect(0, 0, PX, CH);
      ctx.fillRect(ZOOM.x, ZOOM.y, ZOOM.edge, ZOOM.edge);
      panelBorder();   // its sides run through the square: drawn again, unfaded
    }
  }

  // ---- canvas input ---------------------------------------------------------------------
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
      if (!opt.clickFocus()) return;
      const i = Creatures.pick(found, p[0] / WS - 0.5, p[1] / WS - 0.5, CLICK_GATE, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH);
      if (i >= 0) setFocus(i);
      return;
    }
    for (const bottom of patchRows()) {
      for (let i = 0; i < 3; i++) {
        const r = patchRect(i, bottom);
        if (p[0] >= r.x && p[0] < r.x + r.w && p[1] >= r.y && p[1] < r.y + r.h) {
          if (bottom) action = i; else pickChannel(i);
          return;
        }
      }
    }
    if (inZoom(p) && focus && !masking && e.button === 0) {
      if (opt.zoomPress(rot(viewOffset(p)), e)) {   // taken by the page: its drag follows (opt.zoomDrag)
        pageDrag = e.pointerId; ui.setPointerCapture(e.pointerId); e.preventDefault(); return;
      }
      masking = { id: e.pointerId, view: viewOffset(p) };
      ui.setPointerCapture(e.pointerId);
      e.preventDefault();
    }
  });
  // A channel patch (or its key): that channel. With opt.patchErase (the game, which has no action
  // patches), the selected channel's patch again toggles erase (its "recycle") on and off, and
  // another channel goes back to mask.
  function pickChannel(i) {
    if (opt.patchErase()) action = i === channel && action === 0 ? 1 : 0;
    channel = i;
  }
  // Tooltip over the tool patches (they are drawn on the canvas, so the canvas title changes).
  const CH_WORDS = ['red', 'green', 'blue'];
  function patchTip(p) {
    for (const bottom of patchRows()) for (let i = 0; i < 3; i++) {
      const r = patchRect(i, bottom);
      if (p[0] >= r.x && p[0] < r.x + r.w && p[1] >= r.y && p[1] < r.y + r.h)
        return bottom ? [
          `Mask (4): pressing on the ship view hides the ${CH_WORDS[channel]} channel there from the creature's senses`,
          `Erase (5): pressing on the ship view sets ${CH_WORDS[channel]} mass there to 0`,
          `Draw (6): pressing on the ship view sets ${CH_WORDS[channel]} mass there to 1`][i]
          : opt.patchErase() ? `${CH_NAMES[i]} (${i + 1}): mask the ${CH_WORDS[i]} channel; again: erase it (recycle), striped`
          : `${CH_NAMES[i]} (${i + 1}): act on the ${CH_WORDS[i]} channel`;
    }
    return '';
  }
  ui.addEventListener('pointermove', e => {
    const p = canvasPoint(e);
    const tip = patchTip(p);
    if (ui.title !== tip) ui.title = tip;
    ui.style.cursor = tip ? 'pointer' : '';   // a hand over the patches
    if (masking && e.pointerId === masking.id) masking.view = viewOffset(p);
    if (pageDrag === e.pointerId) opt.zoomDrag(rot(viewOffset(p)), e);
    hover = e.pointerType === 'mouse' && inZoom(p) ? viewOffset(p) : null;
    worldHover = e.pointerType === 'mouse' && p[0] < PX ? [p[0] / WS - 0.5, p[1] / WS - 0.5] : null;
  });
  // (opt.keepHover(el): keep the ship-view hover point when the pointer moves onto el, e.g. a
  // button over the view)
  ui.addEventListener('pointerleave', e => { if (!opt.keepHover(e.relatedTarget)) hover = null; worldHover = null; });
  let pageDrag = null;   // the pointer of a press the page took (opt.zoomPress), while held
  const endMask = e => { if (masking && e.pointerId === masking.id) masking = null; if (pageDrag === e.pointerId) pageDrag = null; };
  ui.addEventListener('pointerup', endMask);
  ui.addEventListener('pointercancel', endMask);
  ui.addEventListener('lostpointercapture', endMask);
  ui.addEventListener('contextmenu', e => e.preventDefault());
  // The ship-view hover point from a pointer event elsewhere (a page's button over the view), so
  // what follows the pointer keeps following it there.
  function hoverAt(e) { const p = canvasPoint(e); hover = e.pointerType === 'mouse' && inZoom(p) ? viewOffset(p) : null; }
  // Drop the held tool and the hover ring (e.g. when something covers the ship view).
  function clearPointer() { masking = null; hover = null; }

  // ---- actions (the page's buttons and keys call these) ---------------------------------
  function setRunning(on) { if (status === 'alive' || !on) wantRun = running = on; hooks.change(); }
  function togglePlay() { if (status === 'alive') wantRun = running = !running; hooks.change(); }
  function reset() { spawn(); running = wantRun; hooks.change(); }
  function changeSpeed(d) { setSpeed(speedIdx + d); }
  function setSpeed(i) { speedIdx = Math.max(0, Math.min(SPEEDS.length - 1, i)); hooks.change(); }
  function stepOnce() { if (status !== 'alive') return; wantRun = running = false; advance(1); analyze(1); hooks.change(); }
  function toggleHeadUp() { viewMode = viewMode === 'board' ? 'vel' : 'board'; hooks.change(); }
  function setZoom(i) { zoomIdx = i; CELL_PX = ZOOMS[zoomIdx]; ZOOM_CELLS = ZOOM.edge / CELL_PX; hooks.change(); }
  function setMaskR(r) { maskR = r; hooks.change(); }

  // ---- start ----------------------------------------------------------------------------
  // Call once the page has set its hooks and options. Resolves when the first frame is queued.
  // empty: start with an empty world and no species (the game: its script loads one and puts the
  // creatures in); else the default (or ?form=) pattern's species, spawned at the centre (the lab).
  function start({ empty = false } = {}) {
    resetKinematics();
    hooks.change();
    return SimGL.init(glCanvas, { W: BW, H: BH, R, kernels: patternRule(pattern).kernels })
      .then(() => {
        if (empty) clearWorld(); else { setPattern(pattern.id); spawn(); }
        requestAnimationFrame(frame);
      });
  }

  const core = {
    // constants
    q, R, DEBUG, BW, BH, WS, ZOOM, SPEEDS, ZOOMS, MIN_AREA_UNIT, DETECT, DETECT_MODES, DETECT_NAMES,
    CH_NAMES, ACTIONS, SPEED_MIN, REC_STEPS, STEPS_PER_S, VIEW_LABELS, SPECIES_NAMES, BOOM,
    // live objects (mutated in place)
    stats, md, opt, hooks,
    // the end-of-run banner's second line
    banner: 'Reset or Add to try again',
    // helpers
    deg, wrapD, viewAngle, nearestIdx, findPattern, speciesOf, speciesName, morphName, listed, morphs, hovered, tracking, worldMass,
    // actions
    start, setPattern, selectMorph, spawn, addPattern, putPattern, clearWorld, analyze, setFocus, stopRecord, clearPointer, hoverAt,
    pickChannel,
    setDetectMode, setEveryObject, setRunning, togglePlay, reset, changeSpeed, setSpeed, stepOnce, toggleHeadUp,
    setZoom, setMaskR, setStars, starDrift,
  };
  // live state, read through getters (and set where a page may change it directly)
  const rw = { channel: [() => channel, v => { channel = v; }], action: [() => action, v => { action = v; }] };
  const ro = { running: () => running, status: () => status, steps: () => steps, dt: () => dt, pattern: () => pattern,
               focus: () => focus, focusIdx: () => focusIdx, kin: () => kin, spawnMass: () => spawnMass,
               found: () => found, masking: () => masking, tracks: () => tracks, speedIdx: () => speedIdx,
               maskR: () => maskR, viewMode: () => viewMode, zoomIdx: () => zoomIdx, CELL_PX: () => CELL_PX,
               ZOOM_CELLS: () => ZOOM_CELLS, detectMode: () => detectMode, everyObject: () => everyObject,
               msPerStep: () => msPerStep, simSteps: () => simSteps, driftSec: () => driftSec };
  for (const [k, [get, set]] of Object.entries(rw)) Object.defineProperty(core, k, { get, set, enumerable: true });
  for (const [k, get] of Object.entries(ro)) Object.defineProperty(core, k, { get, enumerable: true });
  return core;
})();
