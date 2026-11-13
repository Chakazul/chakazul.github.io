"use strict";
// ============================================================================================
//  Cognoship Lab (lab.html): the experiment page on top of the shared core (core.js). The
//  lab's controls, its readout, the modes (Action, Stats, Action+Stats, Tuning), the
//  phase-space chart (statschart.js), recording and downloading the long history, the tuning
//  sliders and the keys. The simulation, the focus, its stats and the overlay are in core.js.
// ============================================================================================

(function () {
  const C = Core;
  const { $, flash, fillRange } = Page;
  const readout = $('readout');

  // Page mode, cycled with A (Shift+A backwards): action (play) -> stats (the chart in place of the
  // ship view) -> action+stats (the chart left of the ship view, over part of the world view) ->
  // tuning (detection sliders). Always starts in action mode (no URL parameter).
  const MODES = ['action', 'stats', 'action+stats', 'tuning'];
  const MODE_NAMES = { action: 'Action', stats: 'Stats', 'action+stats': 'Action+Stats', tuning: 'Tuning' };
  let mode = 'action';
  const statsOn = () => mode === 'stats' || mode === 'action+stats';   // the chart is shown

  // What the core's overlay shows: tuning mode, and stats mode with every object on, circle
  // every object; in stats mode with every object on, the one under the mouse is highlighted.
  C.opt.circlesAll = () => mode === 'tuning' || (statsOn() && C.everyObject);
  C.opt.hoverTracks = () => statsOn() && C.everyObject;
  C.opt.stars = () => Page.soundOn;   // the star field follows the Stars button
  C.opt.wholes = () => true;   // twins and trios found and named as one (core.js), the Twin saddle slider
  C.banner = 'Reset (Z) or Spray (X) to try again';

  // Stats chart (statschart.js): 3D or 2D phase space, a stat per axis (stats: t mass chan both hue
  // speed omega). Drag the 3D view to rotate it.
  const chart = { dim: 3,   // 3D by default (not in the URL)
                  axes: ['mass', 'hue', 'speed'],   // defaults (not in the URL)
                  c: 'omega',   // colour axis
                  // Lock axes: keep every axis's min and max; on by default (relock: see the species
                  // hook, which also relocks on the first species's data)
                  lock: true, relock: false,
                  yaw: -2.5, pitch: 0.5 };   // z up, x to the front-left, y to the front-right
  const stats = C.stats;

  // ---- core hooks -------------------------------------------------------------------------
  C.hooks.change = () => syncButtons();
  C.hooks.record = () => syncRecord();
  // locked axes: new ranges for the new species, locked again once its saved history has
  // loaded, or (no file) once the short history has filled (frame)
  C.hooks.species = () => { if (chart.lock) { StatsChart.resetRanges(); chart.relock = true; } };
  C.hooks.statsLoaded = () => { if (chart.relock) { StatsChart.resetRanges(); chart.relock = false; } };   // lock on the loaded history
  C.hooks.pattern = () => {
    const ts = $('t-saddle');
    ts.value = C.DETECT.saddle; $('t-saddle-v').textContent = C.DETECT.saddle.toFixed(2); fillRange(ts);
    const tw = $('t-whole');
    tw.value = C.DETECT.wholeSaddle; $('t-whole-v').textContent = C.DETECT.wholeSaddle.toFixed(2); fillRange(tw);
    syncPicker();
  };
  C.hooks.frame = () => {
    drawReadout();
    if (chart.relock && stats.short.len >= stats.short.cap) chart.relock = false;
    const w = stats.recLeft > 0 && stats.level === 'whole' ? C.focusWhole() : null;   // recording a whole: its line
    if (statsOn()) StatsChart.draw($('c-canvas'), w ? w.hist : stats.short, stats.long, chart, C.everyObject ? C.tracks.map(t => t.hist) : [], C.hovered()?.tr?.hist);
  };

  // ---- readout ----------------------------------------------------------------------------
  function drawReadout() {
    const n = C.found.comps.length, status = C.status, focus = C.focus, kin = C.kin, steps = C.steps;
    const st = status === 'alive' ? `<b class="ok">alive</b>` : `<b class="bad">${status}</b>`;
    // fixed-width numbers (the readout is monospace): pad with no-break spaces
    const pad = (v, n) => String(v).padStart(n, ' ');
    let html = `<span>step <b>${pad(steps, 4)}</b></span><span>t <b>${pad((steps * C.dt).toFixed(1), 6)}</b></span>` +
      `<span>R <b>${C.R}</b></span><span>dt <b>${C.dt}</b></span>` +
      `<span>creatures <b>${pad(n, 3)}</b></span><span>${st}</span>`;
    if (focus) {
      const pct = C.spawnMass ? Math.round(100 * focus.mass / C.spawnMass) : 100;
      const head = pad(kin.speed >= C.SPEED_MIN ? (Math.round(C.deg(-kin.velAng)) % 360 + 360) % 360 : '–', 3) + '°';
      html += `<span>focus [mass <b>${pad(focus.mass.toFixed(0), 4)}</b> (${pad(pct, 3)}%) · ` +
        `speed <b>${kin.speed.toFixed(2)}</b> cells/t · heading <b>${head}</b>]</span>`;
    } else html += `<span>focus [<b>none</b>]</span>`;
    if (stats.recLeft > 0) html += `<span>rec <b class="bad">●</b> ${pad(C.REC_STEPS - stats.recLeft, 3)}/${C.REC_STEPS}` +
      ` ${stats.recFuzzy ? '~' : ''}${stats.recName}</span>`;
    html += `<span>${C.ACTIONS[C.action].toLowerCase()} [${C.CH_NAMES[C.channel]} <b>${C.masking ? 'on ' : 'off'}</b> r=${C.maskR}]</span>`;
    if (C.DEBUG) html += `<span>${C.msPerStep.toFixed(2)} ms/step · ${SimGL.tapCount} taps</span>`;
    // Tuning: the share of the board's mass inside detected objects (the rest is in background
    // blocks, under the threshold, or in basins under the min area)
    const found = C.found;
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

  // ---- controls ---------------------------------------------------------------------------
  function syncButtons() {
    $('play').textContent = C.running ? 'Pause' : 'Play';
    $('headup').textContent = C.viewMode === 'axis' ? 'Body axis up' : 'Heading up';
    $('headup').setAttribute('aria-pressed', C.viewMode !== 'board');
    const sp = C.SPEEDS[C.speedIdx];
    $('speed').value = C.speedIdx + 1;
    $('speed-v').textContent = sp < 1 ? '1/' + Math.round(1 / sp) : sp;
    $('radius-v').textContent = C.maskR;
    $('zoom-v').textContent = C.CELL_PX + '×';
    for (const el of document.querySelectorAll('input[type=range]')) fillRange(el);
  }

  // Two pickers: species and morph (core.js).
  const spSel = $('species'), mSel = $('morph');
  function syncPicker() {
    spSel.textContent = '';
    for (const p of C.listed()) {
      if ([...spSel.options].some(o => o.value === C.speciesOf(p))) continue;
      spSel.add(new Option(C.speciesName(p), C.speciesOf(p)));
    }
    spSel.value = C.speciesOf(C.pattern);
    mSel.textContent = '';
    for (const p of C.morphs()) mSel.add(new Option(C.morphName(p), p.id));
    mSel.value = C.pattern.id;
  }
  // A new species changes the rule: reset the world with its first morph.
  spSel.onchange = () => {
    if (clone.held) endClone();   // (another rule)
    C.setPattern(C.listed().find(p => C.speciesOf(p) === spSel.value).id); C.reset(); spSel.blur();
    flash('species', 'reset');
  };
  // Previous / next species (Shift+C / Shift+V), wrapping round: as picking it in the list.
  function stepSpecies(d) {
    const n = spSel.options.length;
    spSel.selectedIndex = (spSel.selectedIndex + d + n) % n;
    spSel.onchange();
  }
  // A new morph of the same species only changes what Add (X) puts in; the world is kept.
  function selectMorph(id) {
    C.selectMorph(id); mSel.value = C.pattern.id;
    if (clone.held) Clone.swap(clone.held, C.pattern, C.R);   // the held one changes with it
  }
  // Picking a morph in the list, even the current one, also adds it (like Add). A native select
  // fires no change for the option already chosen, so it is cleared while open and restored after.
  const clearMorph = () => { if (mSel.selectedIndex >= 0) mSel.selectedIndex = -1; };
  mSel.addEventListener('pointerdown', clearMorph);
  mSel.addEventListener('focus', clearMorph);
  const restoreMorph = () => { if (mSel.selectedIndex < 0) mSel.value = C.pattern.id; };
  mSel.addEventListener('blur', restoreMorph);
  mSel.addEventListener('keydown', e => { if (e.key === 'Escape') restoreMorph(); });
  mSel.onchange = () => {
    selectMorph(mSel.value); mSel.blur();
    if (clone.held) flash('morph', 'add');   // (holding one to add: that one is held instead)
    else { flash('morph', 'spray'); C.addPattern(); }
  };
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
  const doPlay = () => {
    if (clone.held) { clone.held.wasRunning = false; endClone(); }   // (Play while holding one: dropped)
    flash('play'); C.togglePlay();
  };
  const doStep = () => { flash(...(C.running && C.status === 'alive' ? ['play', 'step'] : ['step'])); C.stepOnce(); };
  const doReset = () => { flash('reset'); C.reset(); };
  const doSpray = () => { flash('spray'); C.addPattern(); };
  const doAdd = () => { flash('add'); startClone(); };
  const doView = () => { flash('headup'); C.toggleHeadUp(); };
  const doSpeed = d => { flash('speed'); C.changeSpeed(d); };
  $('play').onclick = doPlay;
  $('step').onclick = doStep;
  $('reset').onclick = doReset;
  $('spray').onclick = doSpray;
  $('add').onclick = doAdd;
  $('headup').onclick = doView;
  $('speed').value = C.speedIdx + 1;
  $('speed').oninput = e => C.setSpeed(+e.target.value - 1);
  $('zoom').value = C.zoomIdx + 1;
  $('zoom').oninput = e => C.setZoom(+e.target.value - 1);
  $('radius').value = C.maskR;
  $('radius').oninput = e => C.setMaskR(+e.target.value);

  // Modes: tuning and stats each show their controls under the screen, above the help text.
  function applyMode() {
    $('mode').textContent = 'Mode: ' + MODE_NAMES[mode];
    $('mode').setAttribute('aria-pressed', mode !== 'action');   // yellow when not Action
    $('tune').hidden = mode !== 'tuning';
    $('stats').hidden = !statsOn();
    $('c-canvas').hidden = !statsOn();
    $('c-canvas').classList.toggle('beside', mode === 'action+stats');   // left of the ship view, not over it
    $('c-move').hidden = mode !== 'action+stats';
    if (mode === 'stats') C.clearPointer();   // the chart covers the ship view
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
  // Detect: None / Focus / Every (core.js detectMorphs)
  $('detect').onclick = () => {
    flash('detect');
    const M = C.DETECT_MODES;
    C.setDetectMode(M[(M.indexOf(C.detectMode) + 1) % M.length]);
    $('detect').textContent = 'Detect: ' + C.DETECT_NAMES[C.detectMode];
    $('detect').setAttribute('aria-pressed', C.detectMode !== 'none');   // yellow when on
  };
  function nextMode(d = 1) { flash('mode'); mode = MODES[(MODES.indexOf(mode) + d + MODES.length) % MODES.length]; applyMode(); }
  $('mode').onclick = () => nextMode();
  applyMode();

  // ---- stats chart ------------------------------------------------------------------------
  const axisSel = ['c-x', 'c-y', 'c-z'].map($);
  function syncChart() {
    $('c-dim').textContent = chart.dim === 3 ? '3D' : '2D';
    $('c-every').textContent = C.everyObject ? 'Every object: ON' : 'Every object: OFF';
    $('c-every').setAttribute('aria-pressed', C.everyObject);
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
  $('c-every').onclick = () => { flash('c-every'); C.setEveryObject(!C.everyObject); syncChart(); };
  // Record (S), Record Whole (D), Record Fuzzy (F): record the next REC_STEPS steps into the long
  // history, as a new piece: of the focus; of the whole it is part of (a twin, a trio: core.js
  // stats.level 'whole'); or of the focus as a fuzzy morph (detect.js). Each asks for the morph name
  // (offered: the picker's morph; fuzzy: the last fuzzy name recorded; a name need not be a morph in
  // the picker, e.g. cytos). Any of them pressed while recording stops early.
  const REC_BUTTONS = { part: ['c-rec', 'Record'], whole: ['c-rec-whole', 'Record Whole'], fuzzy: ['c-rec-fuzzy', 'Record Fuzzy'] };
  let recKind = 'part';
  function syncRecord() {
    for (const [kind, [id, text]] of Object.entries(REC_BUTTONS)) {
      const on = stats.recLeft > 0 && kind === recKind;
      $(id).textContent = on ? 'Stop' : text;
      $(id).setAttribute('aria-pressed', on);
    }
  }
  const doRecord = (kind = 'part') => {
    const fuzzy = kind === 'fuzzy';
    flash(REC_BUTTONS[stats.recLeft > 0 ? recKind : kind][0]);
    if (stats.recLeft > 0) { stats.recLeft = 0; stats.longGap = true; }   // stop
    else {
      recKind = kind; stats.level = kind === 'whole' ? 'whole' : 'part';
      const last = fuzzy && [...stats.long.segs].reverse().find(sg => sg.fuzzy);
      const name = fuzzy
        ? prompt('Record a fuzzy morph (matched loosely, shown as ~name). Name:', last ? last.morph : '')
        : prompt(stats.level === 'whole' ? 'Record the whole (twin or trio) the focus is part of, as morph:' : 'Record the focus as morph:', C.morphName(C.pattern).replace(C.BOOM, ''));
      if (!name || !name.trim()) return;   // cancelled
      stats.recName = name.trim(); stats.recFuzzy = fuzzy;
      stats.recLeft = C.REC_STEPS; stats.recStart = stats.long.len;
      stats.longGap = true;
      // the piece starts with the focus's short history gathered so far (since its last gap),
      // the rest of the REC_STEPS steps recorded as they come
      const h = stats.level === 'whole' ? C.focusWhole()?.hist : stats.short;
      if (C.focus && h && h.len) {
        let j0 = h.len - 1;
        while (j0 > 0 && !h.brk(j0)) j0--;
        for (let j = j0; j < h.len; j++) {
          const s = Object.fromEntries(['t', 'mass', 'mR', 'mG', 'mB', 'speed', 'omega'].map(k => [k, h.raw(k, j)]));
          if (!Number.isFinite(s.mass)) continue;   // no focus at that step
          stats.long.push(s, stats.longGap, stats.recName, stats.recFuzzy, stats.level);
          stats.longGap = false; stats.recLeft--;
        }
      }
    }
    syncRecord();
  };
  $('c-rec').onclick = () => doRecord('part');
  $('c-rec-whole').onclick = () => doRecord('whole');
  $('c-rec-fuzzy').onclick = () => doRecord('fuzzy');
  // Delete: removes the last piece of the long history; during a recording, the piece being
  // recorded (and stops it; nothing if it has no samples yet).
  $('c-clear').onclick = () => {
    flash('c-clear');
    if (stats.recLeft > 0) {
      C.stopRecord();
      if (stats.long.len > stats.recStart) stats.long.dropLastPiece();
      return;
    }
    stats.long.dropLastPiece();
  };
  // Download: the long history as <species>.json. Saved into assets/stats/ and committed, it is
  // loaded whenever that species is (core.js loadSavedStats); new recordings add to it.
  $('c-save').onclick = () => {
    flash('c-save');
    const obj = { format: 1, species: stats.species, R: C.R, dt: C.dt, ...stats.long.toJSON() };
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([StatsChart.statsFileText(obj)], { type: 'application/json' }));   // aligned, one sample a line
    a.download = `${stats.species}.json`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };
  syncRecord();
  // Step axis i (0 x, 1 y, 2 z, 3 colour) to the next (d = 1) or previous (d = -1) stat: a click on
  // its label, or keys G H J K (Shift: backwards).
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

  // ---- tuning -----------------------------------------------------------------------------
  // Tuning mode: detection sliders (also ?saddle= ?thr= ?area=).
  const D = C.DETECT, U = C.MIN_AREA_UNIT;
  const tune = [['t-saddle', v => D.saddle = v, D.saddle],
                ['t-thr', v => D.threshold = v, D.threshold],
                ['t-area', v => D.minArea = v * U, D.minArea / U],
                ['t-whole', v => D.wholeSaddle = v, D.wholeSaddle]];
  for (const [id, set, v0] of tune) {
    const el = $(id), out = $(id + '-v');
    el.value = v0; out.textContent = (+el.value).toFixed(2); fillRange(el);
    el.oninput = () => { set(+el.value); out.textContent = (+el.value).toFixed(2); C.analyze(0); };
  }

  // ---- Add: a morph held to paste ---------------------------------------------------------
  // Add (B) holds the current morph (clone.js): it follows the pointer over the world view and the
  // ship view alike (shown in both: one board place), turned by ↻ (N) and flipped by ↔ (M) (as
  // seen in the view the pointer was last in), until a click pastes it there, or Esc / ✗ cancels.
  // A touch places it instead (moved by dragging), pasted by ✓ (or Enter). The simulation pauses
  // while it is held, then plays on as before. Picking another morph holds that one instead.
  const clone = { held: null };
  let inWorld = true;   // the pointer was last over the world view (else the ship view)
  $('ui').addEventListener('pointermove', e => {
    const r = $('ui').getBoundingClientRect();
    inWorld = (e.clientX - r.left) / r.width < C.PX / C.CW;
  });
  function startClone() {
    const was = clone.held ? clone.held.wasRunning : C.running;
    C.setRunning(false);
    clone.held = Object.assign(Clone.hold(C.pattern, C.R), { wasRunning: was });
    syncClone();
  }
  function endClone() {
    const was = clone.held && clone.held.wasRunning;
    clone.held = null; syncClone();
    if (was) C.setRunning(true);
  }
  function syncClone() {
    const H = clone.held;
    $('clone-tools').hidden = !H;
    $('clone-ok').hidden = !(H && H.at);   // ✓ once a touch has placed it
    $('add').classList.toggle('clone-on', !!H);
  }
  function flipClone() { Clone.flip(clone.held, inWorld ? 0 : C.viewAngle()); }
  // Paste it at cell at, follow what it became (a pasted creature, or the one it joined)
  function pasteClone(at) {
    Clone.paste(clone.held, at);
    C.analyze(0);
    const i = C.found ? Creatures.pick(C.found, at[0], at[1], 2 * C.R, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
    if (i >= 0) C.setFocus(i);
    endClone();
  }
  // Enter / ✓: where a touch placed it, else under the mouse (either view)
  function confirmClone() {
    const at = clone.held.at || C.pointerCell();
    if (at) { flash('clone-ok'); pasteClone(at); }
  }
  const where = () => clone.held && (clone.held.at || C.pointerCell());
  const press = (at, e) => {
    if (!clone.held) return false;
    if (e && e.pointerType !== 'mouse') { clone.held.at = at; syncClone(); return true; }
    pasteClone(at);
    return true;
  };
  const drag = (at, e) => { if (clone.held && clone.held.at && e && e.pointerType !== 'mouse') clone.held.at = at; };
  C.opt.worldPress = press;
  C.opt.worldDrag = drag;
  const fromShip = off => { const sh = C.ship(); return [sh.col + off[0], sh.row + off[1]]; };   // (off: from the ship view's centre)
  C.opt.zoomPress = (off, e) => C.focus ? press(fromShip(off), e) : false;
  C.opt.zoomDrag = (off, e) => { if (C.focus) drag(fromShip(off), e); };
  // the ghost, in both views, at one board place
  C.opt.worldGhost = () => {
    const H = clone.held, at = where();
    if (!at) return null;
    return { at, r: H.form.w * C.WS * 0.75, draw: (ctx, x, y) => Clone.draw(ctx, H, x, y, 0, C.WS) };
  };
  C.opt.drawGhost = (ctx, hoverAt, centre, cellPx, phi, focus) => {
    const H = clone.held;
    if (!H) return false;
    const at = where();
    if (!at) return true;   // nothing to show, and no tool circle either
    const b = [C.wrapD(at[0] - focus.col, C.BW), C.wrapD(at[1] - focus.row, C.BH)];
    const c = Math.cos(-phi), s = Math.sin(-phi);   // board -> view
    Clone.draw(ctx, H, centre[0] + (b[0] * c - b[1] * s) * cellPx, centre[1] + (b[0] * s + b[1] * c) * cellPx, phi, cellPx);
    return true;
  };
  // the ghost keeps following the pointer over the buttons (they sit on the world view)
  C.opt.keepHover = el => !!clone.held && !!el && $('clone-tools').contains(el);
  $('clone-tools').addEventListener('pointermove', e => { if (clone.held) C.hoverAt(e); });
  $('clone-rot').onclick = () => { if (clone.held) Clone.rotate(clone.held); };
  $('clone-flip').onclick = () => { if (clone.held) flipClone(); };
  $('clone-ok').onclick = () => { if (clone.held) confirmClone(); };
  $('clone-cancel').onclick = () => { if (clone.held) endClone(); };

  // ---- keys -------------------------------------------------------------------------------
  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (e.code === 'Space') { e.preventDefault(); doStep(); }
    else if (k === 'Enter') { e.preventDefault(); if (clone.held) confirmClone(); else doPlay(); }
    else if (k === 'Escape') { if (clone.held) { flash('clone-cancel'); endClone(); } }
    else if (k === 'z') doReset();
    else if (k === 'x') doSpray();
    else if (k === 'b') doAdd();
    else if (k === 'n') { if (clone.held) { flash('clone-rot'); Clone.rotate(clone.held); } }
    else if (k === 'm') { if (clone.held) { flash('clone-flip'); flipClone(); } }
    else if (k === 'c') e.shiftKey ? stepSpecies(-1) : prevMorph();
    else if (k === 'v') e.shiftKey ? stepSpecies(1) : nextMorph();
    else if (k === 'q') doSpeed(1);
    else if (k === 'w') doSpeed(-1);
    else if (k === 'e') doView();
    else if (k === 'a') nextMode(e.shiftKey ? -1 : 1);
    else if (k === 's') doRecord('part');
    else if (k === 'd') doRecord('whole');
    else if (k === 'f') doRecord('fuzzy');
    else if ('ghjk'.includes(k) && k.length === 1) cycleAxis('ghjk'.indexOf(k), e.shiftKey ? -1 : 1);
    else if (e.key >= '1' && e.key <= '3') C.channel = +e.key - 1;
    else if (e.key >= '4' && e.key <= '6') C.action = +e.key - 4;
  });

  // ---- start ------------------------------------------------------------------------------
  C.start().catch(Page.showError);
})();
