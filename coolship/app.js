"use strict";
// ============================================================================================
//  Coolship, Phases 1 and 1b: the controllable ship, a pattern picker, and several creatures
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
  // Two rows of tool patches: top = channel (R G B), bottom = action (Mask Erase Insert); zoom view between.
  const PATCH = { top: 6, bottom: 452, h: 54, w: 104, gap: 16 };
  const ZOOM = { x: PX, y: 64, edge: PW };          // zoom view square
  // Zoom slider: px per cell in the zoom view (7 slots); the default shows about 4 R across.
  const ZOOMS = [2, 3, 4, 6, 8, 12, 16];
  const nearestIdx = (list, v) => list.reduce((b, x, i) => Math.abs(x - v) < Math.abs(list[b] - v) ? i : b, 0);
  let zoomIdx = nearestIdx(ZOOMS, num('zoom', ZOOM.edge / (4 * R), 1, 32));
  let CELL_PX = ZOOMS[zoomIdx], ZOOM_CELLS = ZOOM.edge / CELL_PX;   // cells across the zoom view

  // Kinematics (speeds in cells per unit time, times in units of simulated time).
  // The Aquarium glides at ~0.072 R per unit time (t = steps / T, T = 2).
  const SPEED_MIN = 0.02 * R;    // below: near stationary, arrow hidden
  const SPEED_HI = SPEED_MIN * 1.5;  // hysteresis: "moving" switches on above this
  const VEL_TAU = 7.5;           // velocity smoothing, unit time (= 15 steps; VT049W pulses)
  const ARROW_T = 12;            // arrow = displacement over this much time, at zoom scale

  // Creatures (creatures.js): basins of the potential K * A above 0.05, merged across saddles
  // >= 0.8 x the lower peak, at least a disk of radius R/4 (the values of ~/github/complexity).
  const DETECT = { threshold: 0.05, saddle: 0.8, minArea: Math.PI * (R / 4) ** 2 };
  const TRACK_GATE = R;          // the focus may move at most this far between frames (cells)
  const CLICK_GATE = 2 * R;      // a click this close to a creature selects it

  const CH_COLORS = ['#ff3b3b', '#3bff5a', '#3b7bff'];
  const CH_FILLS = ['rgba(255,59,59,.5)', 'rgba(59,255,90,.5)', 'rgba(59,123,255,.5)'];
  const CH_NAMES = ['R', 'G', 'B'];
  // Actions: 0 mask (hide from the senses), 1 erase (set to 0), 2 insert (set to 1). Circle line:
  // dashed, solid, dotted.
  const ACTIONS = ['Mask', 'Erase', 'Insert'];
  const ACTION_DASH = [[6, 4], [], [2, 4]];

  // ---- state ----------------------------------------------------------------------------
  // Speed slider: 7 slots. Each slot runs SPEEDS[i] * STEPS_PER_SPEED steps per frame, so
  // speed 1 is one step every 4 frames (under LeniaF's rule the glider moves ~0.65 cells a step).
  const SPEEDS = [1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8], STEPS_PER_SPEED = 1 / 4;
  let running = true, speedIdx = nearestIdx(SPEEDS, num('speed', 1, 1 / 8, 8)), stepAcc = 0;
  let maskR = num('radius', Math.round(2.5 * R / 12 * 2) / 2, 1, 8), channel = 1, action = 0;   // selected tool: G mask (best for steering and shape-shifting)
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
  let autoFocus = false;     // after a spawn: focus the creature nearest the board centre
  let kin = null;
  let masking = null;        // { id, view: [x, y] } while the pointer holds a tool
  let hover = null;          // [x, y] view offset of a mouse hovering over the ship view
  let msPerStep = 0;

  const glCanvas = document.getElementById('gl');
  const ui = document.getElementById('ui');
  const ctx = ui.getContext('2d');
  const readout = document.getElementById('readout');

  function resetKinematics() {
    kin = { vel: [0, 0], speed: 0, velAng: 0, moving: false, heading: null, axisHeading: null };
  }

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
    syncPicker();
  }

  function spawn() {
    const f = Rule.formCells(pattern.cells, pattern.R, R);
    SimGL.seed(f, Math.floor(BH / 2 - f.h / 2), Math.floor(BW / 2 - f.w / 2));
    steps = 0; spawnMass = null; status = 'alive'; stepAcc = 0;
    focus = null; autoFocus = true;
    resetKinematics();
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
    const hadFocus = !!focus;
    analyze(0);
    if (!hadFocus) {
      setFocus(Creatures.pick(found, best[0], best[1], 2 * R, SimGL.BW, SimGL.BH, SimGL.BLOCK, BW, BH));
    }
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
      if (i < 0) { focus = null; resetKinematics(); masking = null; }
      else {
        const c = found.comps[i];
        if (T > 0) {
          const v = [wrapD(c.col - focus.col, BW) / T, wrapD(c.row - focus.row, BH) / T];
          const a = 1 - Math.exp(-T / VEL_TAU);
          kin.vel = [kin.vel[0] + a * (v[0] - kin.vel[0]), kin.vel[1] + a * (v[1] - kin.vel[1])];
          kin.speed = Math.hypot(kin.vel[0], kin.vel[1]);
          kin.velAng = Math.atan2(kin.vel[1], kin.vel[0]);
          if (kin.speed > SPEED_HI) kin.moving = true;
          else if (kin.speed < SPEED_MIN) kin.moving = false;
        }
        focus = { col: c.col, row: c.row, mass: c.mass, radius: c.radius };
      }
    }
    if (!focus) return;

    // 'vel' heading: the same smoothed velocity the arrow shows, so the arrow always points
    // straight up. Near stationary (arrow hidden) the last heading is held.
    if (kin.moving || kin.heading === null) kin.heading = kin.velAng;

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
  }

  function setFocus(i) {
    masking = null;
    resetKinematics();
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
      if (n > 0) advance(n);
    }
    analyze(n);
    const p = viewAngle(), c = Math.cos(p), s = Math.sin(p);
    SimGL.render({
      canvasW: CW, canvasH: CH, worldScale: WS,
      zoom: [ZOOM.x, ZOOM.y, ZOOM.edge, CELL_PX],
      rot: [c, s, -s, c],
      focus: focus ? [focus.col, focus.row, 1] : [0, 0, 0],
    });
    drawOverlay();
    drawReadout();
    requestAnimationFrame(frame);
  }

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

    // tool outline, in both views (dashed = mask, solid = erase, dotted = insert)
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

    // grey round-corner border around the ship view
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(160,156,172,.6)'; ctx.lineWidth = 2;
    ctx.beginPath();
    if (ctx.roundRect) ctx.roundRect(ZOOM.x + 3, ZOOM.y + 1, ZOOM.edge - 4, ZOOM.edge - 2, 10);
    else ctx.rect(ZOOM.x + 3, ZOOM.y + 1, ZOOM.edge - 4, ZOOM.edge - 2);
    ctx.stroke();

    // small labels in the zoom view's top corners
    ctx.setLineDash([]);
    ctx.font = '12px ui-monospace, Menlo, Consolas, monospace';
    ctx.textAlign = 'left'; ctx.textBaseline = 'top';
    ctx.fillStyle = 'rgba(200,192,214,.75)';
    ctx.fillText(VIEW_LABELS[viewMode].replace('View: ', ''), ZOOM.x + 8, ZOOM.y + 6);
    ctx.textAlign = 'right';
    ctx.fillText(`${ZOOM_CELLS}×${ZOOM_CELLS} cells`, ZOOM.x + ZOOM.edge - 8, ZOOM.y + 6);

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
    let html = `<span>step <b>${steps}</b></span><span>t <b>${(steps * dt).toFixed(1)}</b></span>` +
      `<span>R <b>${R}</b></span><span>dt <b>${dt}</b></span>` +
      `<span>creatures <b>${n}</b></span><span>${st}</span>`;
    if (focus) {
      const pct = spawnMass ? Math.round(100 * focus.mass / spawnMass) : 100;
      const head = kin.speed >= SPEED_MIN ? deg(-kin.velAng).toFixed(0) + '°' : '–';
      html += `<span>focus [mass <b>${focus.mass.toFixed(0)}</b> (${pct}%) · ` +
        `speed <b>${kin.speed.toFixed(2)}</b> cells/t · heading <b>${head}</b>]</span>`;
    } else html += `<span>focus [<b>none</b>]</span>`;
    html += `<span>${ACTIONS[action].toLowerCase()} [${CH_NAMES[channel]} <b>${masking ? 'on' : 'off'}</b> r=${maskR}]</span>`;
    if (DEBUG) html += `<span>${msPerStep.toFixed(2)} ms/step · ${SimGL.tapCount} taps</span>`;
    readout.innerHTML = html;
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
          `Insert (6): pressing on the ship view sets ${CH_WORDS[channel]} mass there to 1`][i]
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
  });
  ui.addEventListener('pointerleave', () => { hover = null; });
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
  function togglePlay() { if (status === 'alive') running = !running; syncButtons(); }
  function reset() { spawn(); running = true; syncButtons(); }
  function changeSpeed(d) {
    speedIdx = Math.max(0, Math.min(SPEEDS.length - 1, speedIdx + d));
    $('speed').value = speedIdx + 1; syncButtons();
  }
  function stepOnce() { if (status !== 'alive') return; running = false; advance(1); analyze(1); syncButtons(); }
  function toggleHeadUp() { viewMode = viewMode === 'board' ? 'vel' : 'board'; syncButtons(); }

  // Two pickers: species (a rule code, e.g. 5N7KKM) and morph (a form of it, e.g. limus). Each
  // legacy pattern has its own rule settings, so it counts as a species with one morph.
  const speciesOf = p => p.logic === 'legacy' ? p.id : p.rule;
  const SPECIES_NAMES = { '5N7KKM': 'Aquarium', 'VT049W': 'Aquarium replicating', 'G6G6CR': 'Cilium',
    'Z18A9R': 'Aquarium spinning', 'KH97WU': 'Aquarium dancing', 'GDNQYX': 'Aquarium drifting' };
  const SPECIES_LAST = ['G6G6CR'];   // Cilium, the one non-Aquarium species, goes to the bottom
  const speciesName = p => p.logic === 'legacy' ? p.name : `${SPECIES_NAMES[p.rule] || p.rule} (${p.rule})`;
  const morphName = p => p.logic === 'legacy' ? 'rectus' : p.name.replace(p.rule + ' ', '');
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

  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (e.code === 'Space') { e.preventDefault(); doStep(); }
    else if (k === 'Enter') { e.preventDefault(); doPlay(); }
    else if (k === 'z') doReset();
    else if (k === 'x') doAdd();
    else if (k === 'c') prevMorph();
    else if (k === 'v') nextMorph();
    else if (k === 'q') doSpeed(1);
    else if (k === 'w') doSpeed(-1);
    else if (k === 'e') doView();
    else if (e.key >= '1' && e.key <= '3') channel = +e.key - 1;
    else if (e.key >= '4' && e.key <= '6') action = +e.key - 4;
  });

  // theme and fullscreen
  $('theme').onclick = () => {
    const light = document.documentElement.getAttribute('data-theme') !== 'light';
    if (light) document.documentElement.setAttribute('data-theme', 'light');
    else document.documentElement.removeAttribute('data-theme');
    try { localStorage.setItem('coolship-theme', light ? 'light' : 'dark'); } catch (e) {}
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
