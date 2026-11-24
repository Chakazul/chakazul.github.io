"use strict";
// ============================================================================================
//  Cognoship Lab (lab.html): the experiment page on top of the shared core (core.js). The
//  lab's controls, its readout, the modes (Action, Stats, Action+Stats, Model, Tuning), the
//  phase-space chart (statschart.js), recording and downloading the long history, the tuning
//  sliders and the keys. The simulation, the focus, its stats and the overlay are in core.js.
// ============================================================================================

(function () {
  const C = Core;
  const { $, flash, fillRange } = Page;
  const readout = $('readout');

  // Page mode, cycled with A (Shift+A backwards): action (play) -> stats (the chart in place of the
  // ship view) -> action+stats (the chart left of the ship view, over part of the world view) ->
  // model (the morph model's actions drawn on the ship, a click does one) -> tuning (detection
  // sliders). Always starts in action mode (no URL parameter).
  const MODES = ['action', 'stats', 'action+stats', 'model', 'tuning'];
  const MODE_NAMES = { action: 'Action', stats: 'Stats', 'action+stats': 'Action+Stats', model: 'Model', tuning: 'Tuning' };
  let mode = 'action';
  let graphOn = true;   // Model mode: the morph graph shown (Graph; on at first)
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
  C.hooks.species = () => { if (chart.lock) { StatsChart.resetRanges(); chart.relock = true; } if (mode === 'model') loadBook(); };
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
    if (mode === 'model') modelInfo();
    if (chart.relock && stats.short.len >= stats.short.cap) chart.relock = false;
    const w = stats.recLeft > 0 && stats.level === 'whole' ? C.focusWhole() : null;   // recording a whole: its line
    // the other objects' lines: every part's, and every recognized twin or trio's (a whole named, as
    // its circle shows it: core wholeTracks, label set), on the board now (not the one being recorded)
    const others = C.everyObject ? C.tracks.map(t => t.hist)
      .concat(C.wholeTracks.filter(t => t.label && t.wi >= 0 && t !== w).map(t => t.hist)) : [];
    if (statsOn()) StatsChart.draw($('c-canvas'), w ? w.hist : stats.short, stats.long, chart, others, C.hovered()?.tr?.hist);
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
    stopComposite('species changed');
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
  const doReset = () => { stopComposite('reset'); flash('reset'); C.reset(); };
  // a composite running stops on a reset or another species (its ship is gone)
  function stopComposite(why) { if (mod.comp) { mod.compEnd = `${mod.comp.st.name}: stopped (${why})`; mod.comp = null; } }
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
    $('model').hidden = mode !== 'model';
    $('r-holds').hidden = mode !== 'model';
    $('r-graph').hidden = mode !== 'model' || !graphOn;
    if (mode === 'model') { loadBook(); if (C.viewMode === 'board') C.toggleHeadUp(); }   // the book's frame: heading up
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
  C.opt.zoomPress = (off, e) => {
    if (!C.focus) return false;
    if (mode === 'model' && !clone.held) { const ic = iconAt(off); if (ic) { doAction(ic); return true; } }
    return press(fromShip(off), e);
  };
  C.opt.zoomDrag = (off, e) => { if (C.focus) drag(fromShip(off), e); };
  // the ghost, in both views, at one board place
  C.opt.worldGhost = () => {
    const H = clone.held, at = where();
    if (!at) return null;
    return { at, r: H.form.w * C.WS * 0.75, draw: (ctx, x, y) => Clone.draw(ctx, H, x, y, 0, C.WS) };
  };
  C.opt.drawGhost = (ctx, hoverAt, centre, cellPx, phi, focus) => {
    const H = clone.held;
    if (!H && mode === 'model') {   // over an action: the disc it would hold there, in place of the tool circle
      const ic = hoverAt && hovered(), B = book();
      if (!ic || !B) return false;
      const v = C.unrot(ic.b);
      ctx.strokeStyle = C.CH_COLORS[C.CH_NAMES.indexOf(ic.e.ch)]; ctx.lineWidth = 2; ctx.setLineDash([6, 4]);
      ctx.beginPath(); ctx.arc(centre[0] + v[0] * cellPx, centre[1] + v[1] * cellPx, B.mask_r * cellPx, 0, 2 * Math.PI); ctx.stroke();
      ctx.setLineDash([]);
      return true;
    }
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

  // ---- model mode -------------------------------------------------------------------------
  // The morph model of the species (models/<code>.json, from tools/models.py: headless
  // experiments with the game's rule, mask and detection): for each morph, what a G mask (the
  // game's disc) held at a spot in its frame for a while did, from many starts. Drawn on the ship
  // while the name box names a morph of the book: an arrow if it stays that morph (pointing to its
  // new heading), a square lettered with the morph it becomes, grey if it dissolves; blue if its
  // turning direction stays as it is (or it has none), orange if it turns the other way after;
  // faded under 3 of 4 trials. The book's
  // morphs all turn + (mirrored if not): a ship turning − sees it mirrored (left and right swapped,
  // turns negated). A click on one does it, as the player would: the mask held there for the hold.
  // The book and the choice of an action for an outcome: models.js (shared with the game).
  const SAME = '#2a8cff', FLIP = '#ff8c00', GONE = '#808080', SPLIT = '#ff3df2', ODD = '#ff3df2', TWIN = '#2ecc40';
  const mod = { hold: 30, sign: 1, icons: [], last: '', doing: null, done: null, menu: '', outs: [] };
  const book = () => Models.get(C.speciesOf(C.pattern));
  const loadBook = () => Models.load(C.speciesOf(C.pattern));
  // The followed ship's book entry, and its turning direction (kept while it turns too slowly to tell)
  function shipModel() {
    const B = book(), det = C.focus ? C.shipDet() : null;
    const name = det && det.name ? det.name : null;   // (a fuzzy one, ~rectus, as rectus)
    const m0 = B && name ? B.morphs[name] : null;
    mod.sign = Models.sign(C.ship()?.kin.omega ?? 0, mod.sign);
    // an oscillator swept over its swing (HAESRE zigzag): its actions as they work in its phase now
    // (its recent turn rates, kept while it stays this morph: the swing's phase from them)
    const kin = C.ship()?.kin, tn = C.simSteps * C.dt;
    if (!m0 || !m0.phase_bins || !kin || mod.omName !== name) { mod.omHist = []; mod.omName = name; }
    if (m0 && m0.phase_bins && kin && (!mod.omHist.length || tn > mod.omHist[mod.omHist.length - 1][0])) {
      mod.omHist.push([tn, kin.omega, kin.speed / C.R]);
      if (m0.phase_kind !== 'age')   // (an age counts from the first)
        while (mod.omHist.length && mod.omHist[0][0] < tn - 2 * (m0.swing_period || 30) * C.dt) mod.omHist.shift();
    }
    // a still spinner (antigyrans): which way it spins now, from its body axis, against the model's
    // (mod.spinSg: -1 its mirror image); kept while it spins too slowly to tell
    if (m0 && m0.spin) {
      const ax = C.frameAngle('axis');
      if (ax !== null && mod.axPrev && tn > mod.axPrev[0]) {
        let d = ax - mod.axPrev[1]; d -= 2 * Math.PI * Math.round(d / (2 * Math.PI));
        mod.axRate = 0.8 * (mod.axRate || 0) + 0.2 * d / (tn - mod.axPrev[0]);
        if (Math.abs(mod.axRate) > 0.02) mod.spinSg = Math.sign(mod.axRate) * Math.sign(m0.spin);
      }
      if (ax !== null) mod.axPrev = [tn, ax];
    } else { mod.axPrev = null; mod.axRate = 0; mod.spinSg = 1; }
    const m = m0 && m0.phase_bins && m0.phase_useful ? Models.phased(m0, Models.swingBin(m0, mod.omHist, C.dt)) : m0;   // (only where the phase matters)
    return { B, name, m, m0 };
  }
  // The ship's sign for its model: a mover's turning direction (mod.sign), a still spinner's spin
  // against the model's (mod.spinSg), else 0; < 0: it sees the model mirrored
  const shipSg = m => !m ? 0 : m.chirality ? mod.sign : m.spin ? (mod.spinSg || 1) : 0;
  // An action's outcome: [colour, letter or null (an arrow)]: blue the turning direction kept
  // (or none to keep), orange turned the other way; green a single becoming a twin or triple
  function look(e, B, from) {
    const ef = e.effect;
    if (ef === 'keep') return [SAME, null];
    if (ef === 'reverse') return [FLIP, null];
    if (ef === 'dissolve') return [GONE, ''];
    const sp = Models.splitOf(ef);
    if (sp) return [SPLIT, Models.splitMark(sp)];   // splits into two (or more) ships: magenta, their initials (LL, GL)
    const rp = Models.replicateOf(ef);
    if (rp) return [TWIN, rp[0].toUpperCase() + '+'];   // replicates into copies: green (clones or multiplies), T+
    const m = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(ef);
    if (!m) return [ODD, '?'];
    if (Models.isWhole(B, m[1]) && !Models.isWhole(B, from)) return [TWIN, m[1][0].toUpperCase()];
    return [m[2] === 'reversed' ? FLIP : SAME, m[1][0].toUpperCase()];
  }
  function outcomeWords(e, m, sg) {
    const sn = s => s > 0 ? 'turning +' : s < 0 ? 'turning −' : '';
    const tt = e.turn !== undefined ? Math.round(e.turn * (sg < 0 ? -1 : 1)) : null;
    const ef = e.effect, t = tt === null ? '' : `, heading ${tt > 0 ? '+' : ''}${tt}°`;
    if (ef === 'keep') return `stays ${m.base}${sg ? ', ' + sn(sg) : ''}${t}`;
    if (ef === 'reverse') return `stays ${m.base}, now ${sn(-sg)}${t}`;
    if (ef === 'dissolve') return 'dissolves';
    const x = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(ef);
    if (x) return `becomes ${x[1]}${x[2] ? ', ' + sn(x[2] === 'same turn' ? sg : -sg) : ''}`;
    const sp = Models.splitOf(ef);
    if (sp) return `splits into ${sp.join(' + ')}`;
    const rp = Models.replicateOf(ef);
    if (rp) return `replicates into ${rp}${e.copies ? ` (${e.copies[0]}–${e.copies[1]} left at the end)` : ''}`;
    return ef;
  }
  const actWords = (e, spot = e.spot) => `${e.ch} ${spot}, ${e.hold / 30} s`;
  C.opt.wantAxis = () => mode === 'model';
  // Heading up turns as the ship's morph model's frame (its action spots stay put), in any mode
  C.opt.upFrame = () => {
    loadBook();
    const B = book(), det = C.focus ? C.shipDet() : null, m = B && det && det.name ? B.morphs[det.name] : null;
    return !m ? null : m.lithos_type === 'dihedral' ? 'long' : m.lithos_type ? 'board' : m.frame;   // (a lithos-type: its spots' frame)
  };
  // the long axis measured while a dihedral lithos-type is followed in Model mode (its spots' frame)
  C.opt.wantLongAxis = () => {
    if (mode !== 'model') return false;
    const B = book(), det = C.focus ? C.shipDet() : null, m = B && det && det.name ? B.morphs[det.name] : null;
    return !!(m && m.lithos_type === 'dihedral');
  };
  C.opt.drawShip = (ctx, centre, px, phi, sh) => {
    mod.icons = [];
    if (mode !== 'model' || C.auto) return;
    const { B, m } = shipModel();
    // (a still morph with no body axis, its G and B centres together, e.g. GDNQYX lithos: its spots on the board's frame)
    const a0 = C.frameAngle(m ? m.frame : 'vel') ?? (m && (m.frame === 'axis' || m.lithos_type) ? 0 : null);   // (a lithos-type: its own frames, below)
    if (!m || a0 === null) return;
    const R = B.R;
    const sg = shipSg(m), mir = Models.mirrored(m, sg) ? -1 : 1;
    const S = Math.max(7, Math.min(15, 0.1 * R * px));
    const hov = hovered();
    ctx.save();
    ctx.setLineDash([]);
    // a still morph pushed out of place (the model's push: lithos): 8 spots round it on the board's
    // frame (absolute, not its axis), each sending it the opposite way: a square with the new morph's
    // letter and a line the way it will go
    // a lithos-type morph (still and round, the model's push_table): its actions in absolute directions,
    // 8 sides x each spot distance, for the hold chosen: a push from that side. A square: it becomes
    // the lettered morph, a line the way it leaves (two when sideways, either way); an arrow: it stays
    // and slides that way; a ring: nothing much; grey: it dissolves; faded: under 3 of 4 trials
    if (m.lithos_type && m.push_table) {
      // round: the board's directions; dihedral: round its long axis (either end), 0 / 180 along it,
      // 90 / 270 across, the others diagonal (the same effect at the 4 mirror places)
      const dih = m.lithos_type === 'dihedral', base = dih ? C.frameAngle('long') ?? 0 : 0;
      const DIRS = dih ? ['end', 'end-diagonal', 'side', 'end-diagonal', 'other end', 'end-diagonal', 'side', 'end-diagonal']
                       : ['right', 'lower right', 'below', 'lower left', 'left', 'upper left', 'above', 'upper right'];
      const SIDE = k => k % 4 === 0 ? 'along' : k % 4 === 2 ? 'across' : 'diagonal';
      const rows = m.push_table.filter(t => t.hold === mod.hold);
      if (dih) {   // its long axis as measured, end to end (a check: the spots are placed on it)
        const u = C.unrot([Math.cos(base), Math.sin(base)]), L = 1.6 * R * px;
        ctx.beginPath(); ctx.moveTo(centre[0] - u[0] * L, centre[1] - u[1] * L); ctx.lineTo(centre[0] + u[0] * L, centre[1] + u[1] * L);
        ctx.setLineDash([6, 5]); ctx.lineWidth = 1.5; ctx.strokeStyle = 'rgba(255,255,255,0.6)'; ctx.globalAlpha = 1; ctx.stroke(); ctx.setLineDash([]);
      }
      for (let k = 0; k < 8; k++) for (const pt of rows.filter(t => !dih || t.side === SIDE(k))) {
        const th = base + k * Math.PI / 4, d = pt.r * R, b = [d * Math.cos(th), d * Math.sin(th)], v = C.unrot(b);
        const x = centre[0] + v[0] * px, y = centre[1] + v[1] * px, away = dih ? 'the far side' : DIRS[(k + 4) % 8];
        const into = /^into (\S+)/.exec(pt.effect), dev = pt.dev ?? 0;
        const shove = pt.effect === 'keep' && pt.move >= 0.3 && (pt.away ?? 0) >= 0.5;
        const what = into ? (dev >= 135 ? `it leaves back towards the mask (${dev}° off ${away})` : dev >= 45 ? `it leaves sideways, ${dev}° off ${away}` : `it leaves ${away}`)
          : pt.effect === 'dissolve' ? 'it dissolves' : shove ? `slides it ${pt.move} R ${away}` : 'nothing much';
        const e = { ch: 'G', spot: `${DIRS[k]} ${pt.dist}${dih ? ` (${pt.side} the long axis)` : ''} (${what})`, hold: pt.hold, p: pt.p, effect: pt.effect, off: [0, 0], always_named: false };
        mod.icons.push({ b, e, fs: b, spot: e.spot, frame: 'board' });
        const hot = hov && hov.e.spot === e.spot;
        ctx.globalAlpha = pt.p >= 0.75 ? 1 : 0.35;
        const edge = () => { ctx.lineWidth = hot ? 3 : 1.5; ctx.strokeStyle = hot ? '#fff' : '#000'; ctx.stroke(); };
        const goes = (dev >= 135 ? [180] : dev >= 45 ? [dev, -dev] : [0]).map(t => C.unrot([-Math.cos(th + t * Math.PI / 180), -Math.sin(th + t * Math.PI / 180)]));
        if (into) {
          for (const go of goes) { ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + go[0] * S * 2.2, y + go[1] * S * 2.2); ctx.strokeStyle = SAME; ctx.lineWidth = 2.5; ctx.stroke(); }
          ctx.beginPath(); ctx.rect(x - S * 0.8, y - S * 0.8, S * 1.6, S * 1.6); ctx.fillStyle = SAME; ctx.fill(); edge();
          ctx.fillStyle = '#000'; ctx.font = `bold ${Math.round(S * 1.2)}px system-ui, sans-serif`; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
          ctx.fillText(into[1][0].toUpperCase(), x, y + 1);
        } else if (pt.effect === 'dissolve') {
          ctx.beginPath(); ctx.rect(x - S * 0.8, y - S * 0.8, S * 1.6, S * 1.6); ctx.fillStyle = GONE; ctx.fill(); edge();
        } else if (shove) {
          const dd = goes[0], n = [-dd[1], dd[0]];
          const P = (f, sd) => [x + (dd[0] * f + n[0] * sd) * S, y + (dd[1] * f + n[1] * sd) * S];
          const pts = [P(1.1, 0), P(0.1, 0.75), P(0.1, 0.3), P(-1, 0.3), P(-1, -0.3), P(0.1, -0.3), P(0.1, -0.75)];
          ctx.beginPath(); pts.forEach((q, i) => i ? ctx.lineTo(...q) : ctx.moveTo(...q)); ctx.closePath(); ctx.fillStyle = SAME; ctx.fill(); edge();
        } else {
          ctx.beginPath(); ctx.arc(x, y, S * 0.45, 0, 2 * Math.PI); ctx.lineWidth = 2.5; ctx.strokeStyle = SAME; ctx.stroke();
          if (hot) { ctx.lineWidth = 1.5; ctx.strokeStyle = '#fff'; ctx.stroke(); }
        }
      }
      ctx.restore();
      return;
    }
    // one channel's actions at a time (a book of two, e.g. XEH4YR's G and B, at the same spots):
    // the selected tool channel's, else the book's first
    const chs = [...new Set(m.actions.map(e => e.ch))], ch = chs.includes(C.CH_NAMES[C.channel]) ? C.CH_NAMES[C.channel] : chs[0];
    for (const e of m.actions) {
      if (e.hold !== mod.hold || e.ch !== ch) continue;
      const fs = [e.off[0] * R, e.off[1] * R * mir];
      const b = C.frameOffset(fs, m.frame), v = C.unrot(b);
      const x = centre[0] + v[0] * px, y = centre[1] + v[1] * px;
      mod.icons.push({ b, e, fs, spot: Models.seen(B, m, e, sg).spot });
      const [col, letter] = look(e, B, m.morph);
      ctx.globalAlpha = e.p >= 0.75 ? 1 : 0.35;
      ctx.lineWidth = hov && hov.e === e ? 3 : 1.5;
      ctx.strokeStyle = hov && hov.e === e ? '#fff' : '#000';
      ctx.fillStyle = col;
      if (letter === null) {   // an arrow to the new heading
        const h = m.frame === 'vel' ? a0 + (e.turn || 0) * mir * Math.PI / 180 : a0;
        const d = C.unrot([Math.cos(h), Math.sin(h)]), n = [-d[1], d[0]];
        const P = (f, s) => [x + (d[0] * f + n[0] * s) * S, y + (d[1] * f + n[1] * s) * S];
        const pts = [P(1.1, 0), P(0.1, 0.75), P(0.1, 0.3), P(-1, 0.3), P(-1, -0.3), P(0.1, -0.3), P(0.1, -0.75)];
        ctx.beginPath(); pts.forEach((p, i) => i ? ctx.lineTo(...p) : ctx.moveTo(...p)); ctx.closePath();
        ctx.fill(); ctx.stroke();
      } else {
        ctx.beginPath(); ctx.rect(x - S * 0.8, y - S * 0.8, S * 1.6, S * 1.6); ctx.fill(); ctx.stroke();
        if (letter) {
          ctx.fillStyle = '#000'; ctx.font = `bold ${Math.round(S * (letter.length > 1 ? 0.85 : 1.2))}px system-ui, sans-serif`;
          ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.fillText(letter, x, y + 1);
        }
      }
    }
    ctx.restore();
  };
  // The icon under the mouse (board offset from the ship, cells), or near off
  function iconAt(off) {
    if (!off || !mod.icons.length) return null;
    const B = book(), lim = Math.max(3, 0.12 * (B ? B.R : 18));
    let best = null, bd = lim;
    for (const ic of mod.icons) {
      const d = Math.hypot(ic.b[0] - off[0], ic.b[1] - off[1]);
      if (d < bd) { bd = d; best = ic; }
    }
    return best;
  }
  function hovered() {
    const at = C.pointerCell(true), sh = C.ship();
    if (!at || !sh) return null;
    return iconAt([C.wrapD(at[0] - sh.col, C.BW), C.wrapD(at[1] - sh.row, C.BH)]);
  }
  // Do it: the mask held there for the hold, as the player would (core pressFor). ic: an icon or a
  // seen action (models.js), { e, fs, spot }; why: the outcome picked, if any
  function doAction(ic, why = '') {
    const { B, m } = shipModel(), sg = shipSg(m);
    if (ic.e.phase !== undefined && m.phase_now !== ic.e.phase) { mod.waitFor = { ic, why, morph: m.morph }; return; }   // (one phase's action: when it comes, modelInfo)
    mod.waitFor = null;
    mod.doing = { e: ic.e, spot: ic.spot, m, sg, why, start: C.steps };
    mod.done = null;
    // (a push: the board's frame)
    C.pressFor({ fs: ic.fs, frame: ic.frame || m.frame, steps: ic.e.hold, ch: C.CH_NAMES.indexOf(ic.e.ch), r: B.mask_r,
                 done: () => { mod.done = mod.doing; mod.doing = null; } });
  }
  // An action's details, the same wherever it is hovered (a spot on the ship, an arrow in the graph,
  // an item of the Maneuver menu), shown after the Maneuver menu: the morph, the mask (channel, spot,
  // hold), what it does (its colour), p, the name box
  const COLOUR_NAME = { [SAME]: 'blue', [FLIP]: 'orange', [GONE]: 'grey', [SPLIT]: 'magenta', [TWIN]: 'green' };
  function detailText(B, m, e, spot, sg = 1) {   // (sg < 0: a ship turning −, its turns mirrored as its spots)
    const [col] = look(e, B, m.morph), ef = e.effect, b = m.morph;
    const x = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(ef), sp = Models.splitOf(ef), rp = Models.replicateOf(ef);
    const tu = e.turn === undefined ? null : e.turn * (Models.mirrored(m, sg) ? -1 : 1);
    const tt = tu !== null && Math.abs(tu) >= 1 ? `, turn ${tu > 0 ? '+' : ''}${Math.round(tu)}°` : '';
    const whatOf = (ef, turn) => {
      const x = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(ef), sp = Models.splitOf(ef), rp = Models.replicateOf(ef);
      return ef === 'keep' ? `${b}${turn}` : ef === 'reverse' ? `${b}, change turning${turn}` : ef === 'dissolve' ? `${b} → dissolved`
        : ef === 'unsettled' ? `${b} → ???` : x ? `${b} → ${x[1]}${x[2] === 'reversed' ? ', change turning' : ''}` : sp ? `${b} → ${sp.join(' + ')}`
        : rp ? `${b} → replicates into ${rp}${e.copies && ef === e.effect ? ` (${e.copies[0]}–${e.copies[1]} left at the end)` : ''}` : `${b} → ?`;
    };
    // (the runner-up outcome, if any; the same morph, its turning not judged: said so, not as if the same)
    const same = e.alt && x && /^into \S+$/.test(e.alt[0]) && e.alt[0] === `into ${x[1]}`;
    const what = whatOf(ef, tt), alt = e.alt ? ` · else ${whatOf(e.alt[0], '')}${same ? ', its turning not known' : ''} (p ${e.alt[1].toFixed(2)})` : '';
    const chs = new Set(m.actions.map(s => s.ch));
    return `${COLOUR_NAME[col] || 'magenta'} · ${chs.size > 1 ? e.ch + ' ' : ''}${spot} · ${e.hold / 30} s${e.phase !== undefined ? ` · at phase ${e.phase + 1}/${m.phase_bins}` : ''} · p ${e.p.toFixed(2)}${e.p_sweep !== undefined ? ` (${e.p_sweep.toFixed(2)} in the sweep${e.p_fresh !== undefined ? `, ${e.p_fresh.toFixed(2)} from fresh starts` : ''}; ${e.p_robust < 0.5 ? 'fragile: ' : ''}${e.p.toFixed(2)} placed a little off)` : ''} · ${what}${alt}`;
  }
  // A morph's characteristics (hovered in the graph): twin or single, its motion (R per game
  // second; still), rotating (a mover's turning, a still one's spin, drifting), oscillating, chaotic,
  // unstable, notes; from the book (tools/models.py props and the author's notes)
  function morphText(B, n, source) {
    const m = B.morphs[n], deg = w => Math.round(Math.abs(w) * C.STEPS_PER_S * C.dt * 180 / Math.PI);
    const parts = [Models.isWhole(B, n) ? (n.startsWith('tri') ? 'triple' : 'twin') : 'single'];
    if (!m) return `${n}: ${parts[0]} · no model of its own (only reached from others)`;
    const sp = m.speed !== undefined ? m.speed * C.STEPS_PER_S * C.dt : null;
    parts.push(sp !== null ? (sp < 0.1 ? 'still' : `moves ${sp.toFixed(2)} R/s${sp < 0.45 ? ' (slow)' : ''}`) : m.frame === 'axis' ? 'still' : 'moves');
    if (m.drift) parts.push(`drifting ${m.omega > 0 ? '↻' : '↺'}${deg(m.omega) ? ` ${deg(m.omega)}°/s` : ''} (a slow spin)`);
    else if (m.chirality) parts.push(`turns ${m.omega > 0 ? '↻ clockwise' : '↺ anticlockwise'} ${deg(m.omega)}°/s`);
    else if (m.spin) parts.push(`spins ${m.spin > 0 ? '↻' : '↺'} ${deg(m.spin)}°/s`);
    if (m.oscillating) parts.push(`∿ oscillating${m.swing_period ? `, every ${m.swing_period} steps` : ''}${m.phase_useful ? ' (its actions depend on the phase)' : ''}`);
    if (m.chaotic) parts.push(`≋ chaotic: ${m.chaotic}`);
    if (m.unstable) parts.push(`! unstable: ${m.unstable}`);
    if (m.note) parts.push(m.note);
    if (m.no_maneuver) parts.push('not for maneuvers in the game');
    if (m.lithos_type) parts.push(m.lithos_type === 'dihedral' ? 'dihedral lithos-type: its actions by their angle to its long axis (along, diagonal, across; either end up)' : 'round lithos-type: its actions in absolute directions (a push from one side)');
    for (const pu of m.pushes || []) parts.push(pu.kind === 'shove'
      ? `can be shoved ${pu.move} R from any side (G mask at its ${pu.dist}, ${pu.hold / 30} s; p ${pu.p})`
      : `can be pushed into ${pu.into} from any side (G mask at its ${pu.dist}, ${pu.hold / 30} s; it leaves ${pu.dev ? `${pu.dev}° off the far side, either way` : 'the far side'}; p ${pu.p})`);
    return `${n}: ${parts.join(' · ')} · click to reset with it`;
  }
  // what each hovers now, the first shown; else the composite running (after a green ▶)
  const details = { menu: '', graph: '', ship: '', comp: '' };
  const esc = t => t.replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' })[c]);
  function showDetail(src, t) {
    details[src] = t;
    const el = $('r-outcome-act'), v = details.menu || details.graph || details.ship;
    const html = v ? esc(v) : details.comp ? esc(details.comp) : '';
    $('r-comp-run').hidden = !mod.comp;   // (a composite running: ▶ beside the Composite menu; not its end's message)
    if (el.innerHTML !== html) el.innerHTML = html;
  }
  function modelInfo() {
    const { B, name, m } = shipModel();
    let t;
    const sp = C.speciesName(C.pattern);
    compTick();
    const sr = mod.comp;
    if (sr) mod.compEnd = null;
    showDetail('comp', !sr && mod.compEnd ? mod.compEnd : sr ? `${sr.st.morph} → ${sr.st.name}${sr.r.val ? ` (w ${sr.r.val})` : ''}: ${!sr.started ? `waiting for ${sr.st.morph}` : `step ${sr.i + 1}/${sr.steps.length}, ${sr.steps[sr.i] ? sr.steps[sr.i].join(' ') : 'done'}`} · ${compText(sr)}` : '');
    if (mod.waitFor && !mod.doing) {
      const w = mod.waitFor, { m: mw } = shipModel();
      if (!mw || mw.morph !== w.morph || C.auto) mod.waitFor = null;   // (it became another morph: no more)
      else if (mw.phase_now === w.ic.e.phase) doAction(w.ic, w.why);
    }
    if (mod.comp && !mod.doing && !mod.waitFor) t = `composite ${mod.comp.st.name}${mod.comp.r.val ? ` (w ${mod.comp.r.val})` : ''}: step ${mod.comp.i + 1}/${mod.comp.steps.length}, ${mod.comp.steps[mod.comp.i] ? mod.comp.steps[mod.comp.i].join(' ') : 'done'}`;
    else if (mod.waitFor && !mod.doing) {
      const w = mod.waitFor, { m: mw } = shipModel();
      t = `${w.why ? w.why + ': ' : ''}waiting for phase ${w.ic.e.phase + 1}/${mw ? mw.phase_bins : '?'} (now ${mw && mw.phase_now !== undefined ? mw.phase_now + 1 : '?'}) to do ${actWords(w.ic.e, w.ic.spot)}`;
    } else if (mod.doing) {
      const d = mod.doing, a = C.auto;
      t = `${d.why ? d.why + ': ' : ''}doing ${actWords(d.e, d.spot)}: ${a ? d.e.hold - a.left : d.e.hold}/${d.e.hold} steps${C.running ? '' : ' (paused)'} · expected: ${outcomeWords(d.e, d.m, d.sg)} (p ${d.e.p.toFixed(2)}; name box ${d.e.always_named ? 'always named' : d.e.seq})`;
    } else if (B === undefined) t = 'loading the morph model…';
    else if (!B) t = `no morph model for ${sp} yet`;
    else if (C.detectMode === 'none') t = 'Detect is off: turn it on to see the model';
    else if (!C.focus) t = 'no ship followed';
    else if (!name) t = '???: not a known morph now, so no actions (not maneuverable)';
    else if (!m) t = `${name}: not in the morph model yet`;
    else {
      const ic = hovered();
      const sg = shipSg(m);
      const turning = m.chirality ? `turning ${sg > 0 ? '+ (clockwise)' : '− (anticlockwise: the model mirrored)'}` : 'no turning direction';
      // a model of several mask channels (XEH4YR: G and B): the one shown, and how to switch
      const chs = [...new Set(m.actions.map(e => e.ch))], chNow = chs.includes(C.CH_NAMES[C.channel]) ? C.CH_NAMES[C.channel] : chs[0];
      const chText = chs.length > 1 ? ` · models for ${chs.join(' and ')}: showing ${chNow} (switch with the ${chs.join(' / ')} buttons or 1–3)` : '';
      showDetail('ship', ic ? detailText(B, m, ic.e, ic.spot, sg) : '');
      const phase = !m.phase_bins ? '' : !m.phase_useful ? ` · swept by ${m.phase_kind === 'state' ? 'state' : 'phase'}: it hardly matters` : m.phase_kind === 'state'
        ? ` · state: ${m.phase_now !== undefined ? m.states.names[m.phase_now] : '?'} (the actions as they work in it)`
        : ` · ${m.phase_signal === 'speed' ? 'pulse' : 'swing'} phase ${m.phase_now !== undefined ? m.phase_now + 1 : '?'}/${m.phase_bins} (the actions as they work now)`;
      t = `${name}, ${turning}${phase} · ${m.actions.filter(e => e.hold === mod.hold && e.ch === chNow).length} actions held ${mod.hold / 30} s${chText} · hover one for what it does, click to do it`;
      if (!ic && mod.done) t = `${mod.done.why ? mod.done.why + ': ' : ''}did ${actWords(mod.done.e, mod.done.spot)} (expected: ${outcomeWords(mod.done.e, mod.done.m, mod.done.sg)}) · ` + t;
    }
    if (!m || C.detectMode === 'none' || !C.focus) showDetail('ship', '');
    if (t !== mod.last) { $('r-info').textContent = t; mod.last = t; }
    syncManeuvers(B, name, m);
    syncComposites(B, name);
    drawGraph();
  }
  // The Maneuver menu: what the ship can reach now (models.js maneuvers), rebuilt when that changes;
  // picking one does its best action at once, then the menu shows its prompt again
  const oSel = $('r-outcome');
  function syncManeuvers(B, name, m) {
    const m0 = B && name ? B.morphs[name] : null;   // (what it can do in any phase: an action for one phase waits for it)
    const outs = m0 && !C.auto ? Models.maneuvers(B, name, shipSg(m0), m0, true) : [];
    mod.menuCtx = { B, m, sg: shipSg(m) };
    mod.outs = outs;
    const sig = outs.map(o => o.key).join() + '|' + (C.auto ? 'busy' : m ? 'm' : 'none');
    if (sig === mod.menu || !oMenu.hidden) return;   // (not while it is open)
    mod.menu = sig;
    oSel.textContent = C.auto ? 'doing one…' : outs.length ? 'pick one… ▾' : 'none now';
    oSel.disabled = !outs.length;
    oMenu.textContent = '';
    let g = null;
    for (const o of outs) {
      if (g !== o.group) { g = o.group; const h = document.createElement('div'); h.className = 'omenu-head'; h.textContent = g; oMenu.appendChild(h); }
      const t = o.best.turn !== null && o.key.startsWith('turn-') ? ` (${o.best.turn > 0 ? '+' : ''}${Math.round(o.best.turn)}°)` : '';
      const it = document.createElement('div');
      it.className = 'omenu-item'; it.textContent = o.label + t; it.dataset.key = o.key;
      oMenu.appendChild(it);
    }
  }
  // ---- Composite menu: the model's composite maneuvers starting from the ship's morph, one item per value
  // tried (e.g. a turn's wait), with p, steps and net turn; a click runs it step by step (compTick)
  const sSel = $('r-comp'), sMenu = $('r-comp-menu');
  let compItems = [];
  // Its kinds, from its steps (val: the value tried, for {w}): composite (+, two maneuvers or more),
  // repeating (↻, a step retried: xN), a long hold (⏱, beyond the sweeps' 90 steps: no maneuver
  // does it). A single maneuver done once is none of them: the Maneuver menu has it.
  function compKinds(st, val) {
    const steps = (val ? st.steps.replace('{w}', val) : st.steps).split(';').map(w => w.trim()).filter(w => w && !w.startsWith('wait'));
    const rep = steps.some(w => /\sx([2-9]|\d\d+)\b/.test(w));
    const longHold = steps.some(w => { const m = /^hold\s.*?\s(\d+)(\s+[RGB])?$/.exec(w); return m && +m[1] > 90; });
    return (steps.length > 1 ? '+' : '') + (rep ? '↻' : '') + (longHold ? '⏱' : '');
  }
  // A composite's results worth listing: all, or for a long list of a turn's waits only the most
  // reliable small and big turn each way (big: 90 degrees or more), each with its label
  function compResults(st) {
    const rs = st.results.filter(r => r.p > 0);
    if (rs.length <= 4) return rs.map(r => ({ r, tag: r.val ? `w ${r.val}` : '' }));
    const turned = rs.filter(r => r.turn !== null && r.turn !== undefined);
    if (turned.length < 4) return rs.sort((a, b) => b.p - a.p).slice(0, 4).map(r => ({ r, tag: `w ${r.val}` }));
    const best = {};
    for (const r of turned) {
      const k = (Math.abs(r.turn) >= 90 ? 'big' : 'small') + ' turn ' + (r.turn >= 0 ? 'right' : 'left');
      if (!best[k] || r.p > best[k].p) best[k] = r;
    }
    return ['small turn left', 'big turn left', 'small turn right', 'big turn right'].filter(k => best[k])
      .map(k => ({ r: best[k], tag: `${k}, ${best[k].turn > 0 ? '+' : ''}${Math.round(best[k].turn)}°` }));
  }
  function syncComposites(B, name) {
    // (every composite of the species, but a single maneuver done once: one starting from another
    // morph resets the world with it)
    const items = B ? (B.composites || []).flatMap(st => compResults(st).map(({ r, tag }) => ({ st, r, tag, kinds: compKinds(st, r.val) })))
      .filter(x => x.kinds) : [];
    const sig = items.map(x => x.st.name + x.r.val).join();
    if (sig === mod.compSig || !sMenu.hidden) return;
    mod.compSig = sig; compItems = items;
    sSel.textContent = items.length ? 'pick one… ▾' : 'none'; sSel.disabled = !items.length;
    sMenu.textContent = '';
    items.forEach((x, k) => {
      const it = document.createElement('div'); it.className = 'omenu-item'; it.dataset.k = k;
      it.textContent = `${x.kinds} ${x.st.name}${x.tag ? ` (${x.tag})` : ''}`;   // (its kinds, its name and which result; the rest on hover)
      sMenu.appendChild(it);
    });
  }
  const compText = x => `${x.st.note ? x.st.note + ' · ' : ''}${x.r.val ? x.st.steps.replace('{w}', x.r.val) : x.st.steps} · p ${x.r.p.toFixed(2)} · ${x.r.time} steps${x.r.tries ? `, ${x.r.tries} tries` : ''}${x.r.turn !== null && x.r.turn !== undefined ? ` · net turn ${x.r.turn > 0 ? '+' : ''}${Math.round(x.r.turn)}° ± ${Math.round(x.r.turn_sd)}` : ''}${x.r.travel !== null && x.r.travel !== undefined ? ` · travel ${x.r.travel > 0 ? '+' : ''}${x.r.travel} R vs untouched` : ''}${Object.keys(x.r.fails || {}).length ? ` · fails: ${Object.entries(x.r.fails).map(([k, v]) => `${k} ${v}`).join(', ')}` : ''}`;
  sSel.onclick = e => { e.stopPropagation(); if (sMenu.hidden && compItems.length && !C.auto) sMenu.hidden = false; else { sMenu.hidden = true; showDetail('menu', ''); } };
  sMenu.addEventListener('mouseover', e => { const x = compItems[e.target.dataset && e.target.dataset.k]; showDetail('menu', x ? compText(x) : ''); });
  sMenu.addEventListener('click', e => {
    e.stopPropagation();
    const x = compItems[e.target.dataset && e.target.dataset.k];
    sMenu.hidden = true; showDetail('menu', '');
    if (!x) return;
    const steps = (x.r.val ? x.st.steps.replace('{w}', x.r.val) : x.st.steps).split(';').map(s => s.trim().split(/\s+/)).filter(w => w[0]);
    if (shipModel().name !== x.st.morph) {   // (another morph: the world reset with it, at the centre)
      const pat = C.morphs().find(p => C.morphName(p).replace(C.BOOM, '') === x.st.morph);
      if (pat) { selectMorph(pat.id); doReset(); }
    }
    mod.comp = { ...x, steps, i: 0, t0: C.steps, tries: 0, started: false };   // (after the reset, which stops one running)
  });
  document.addEventListener('click', () => { if (!sMenu.hidden) { sMenu.hidden = true; showDetail('menu', ''); } });
  // steps after a try's hold before it is judged: the morph settled and the name box named it (it
  // needs ~50 steps of history; the experiments judge ~130 steps after)
  const COMP_SETTLE = 150;
  // the running composite, each frame: its step done, or its next action (as the Maneuver menu's)
  function compTick() {
    const sr = mod.comp;
    if (!sr || mod.doing || mod.waitFor || C.auto) return;
    const { B, name, m0 } = shipModel();
    if (!sr.started) {   // (it begins once the ship is its morph: after a reset, once named)
      if (name === sr.st.morph) { sr.started = true; sr.t0 = C.steps; }
      else if (C.steps - sr.t0 > 300) { mod.comp = null; mod.compEnd = `${sr.st.name}: not started (the ship never ${sr.st.morph})`; }
      return;
    }
    for (let hop = 0; hop < 4; hop++) {
      const w = sr.steps[sr.i];
      if (!w) { mod.compEnd = `${sr.st.name}: done (${name || '???'})`; mod.comp = null; return; }
      // ??? (no name: a morph changing between two, the name box catching up): wait, up to a try's
      // allowance, then judged as failed below
      if (!B || !name) {
        if (!sr.unnamed) sr.unnamed = C.steps;
        if (C.steps - sr.unnamed < COMP_SETTLE + 150) return;
        mod.compEnd = `${sr.st.name}: stopped, the ship unnamed (???) too long`; mod.comp = null; return;
      }
      sr.unnamed = null;
      if (w[0] === 'wait') { if (C.steps - sr.t0 < +w[1]) return; sr.i++; sr.t0 = C.steps; sr.tries = 0; continue; }
      if (w[0] === 'hold') {   // hold <spot | fwd,side in R> <steps> [ch]: a mask held there, as tried
        if (sr.tries) { if (C.steps - sr.t0 < 0) return; sr.i++; sr.t0 = C.steps; sr.tries = 0; continue; }
        const ww = [...w], chn = ['R', 'G', 'B'].includes(ww[ww.length - 1]) ? ww.pop() : 'G', nst = +ww[ww.length - 1], spot = ww.slice(1, -1).join(' ');
        const off = /^-?[0-9.]+,-?[0-9.]+$/.test(spot) ? spot.split(',').map(Number) : ((m0.actions.find(e => e.spot === spot) || {}).off);
        if (!off) { mod.comp = null; return; }
        sr.tries++; sr.t0 = C.steps + nst;
        doAction({ e: { ch: chn, spot, hold: nst, p: 1, effect: 'keep', always_named: false }, fs: [off[0] * B.R, off[1] * B.R], spot }, `composite ${sr.st.name}`);
        return;
      }
      const most = w[0] === 'into' && w[2] && w[2][0] === 'x' ? +w[2].slice(1) : 1;
      if (w[0] === 'into' && name === w[1]) { sr.i++; sr.t0 = C.steps; sr.tries = 0; continue; }
      if (w[0] === 'do' && sr.tries) { sr.i++; sr.t0 = C.steps; sr.tries = 0; continue; }
      if (sr.tries && C.steps - sr.t0 < COMP_SETTLE) return;   // (a try settling and named: the name box needs its history)
      if (sr.tries >= most) { mod.compEnd = `${sr.st.name}: failed, ${name} after ${sr.tries} tries at "${w.join(' ')}"`; mod.comp = null; return; }   // (failed)
      // (the step's action as the experiments chose it: however reliable; Models.stepAction)
      const a = Models.stepAction(B, m0, w[0] === 'into' ? `into ${w[1]}` : w.slice(1, w.includes('from') ? w.indexOf('from') : w.length).join(' '), shipSg(m0), shipModel().m?.phase_now);   // (an age-swept morph: an age to come)
      if (!a) { mod.comp = null; mod.compEnd = `${sr.st.name}: stopped, no action from ${name} for "${w.join(' ')}"`; return; }
      // (one for one phase waits for it, doAction)
      sr.tries++; sr.t0 = C.steps + a.e.hold;
      mod.hold = a.e.hold; syncHold();
      doAction(a, `composite ${sr.st.name}`);
      return;
    }
  }
  // the menu: open on a click; an item hovered shows the action it would take after the menu
  // (#r-outcome-act); a click does it
  const oMenu = $('r-outcome-menu'), oAct = $('r-outcome-act');
  const outcomeAct = o => { const c = mod.menuCtx; return c && c.m ? detailText(c.B, c.m, o.best.e, o.best.spot, c.sg) : ''; };
  const closeMenu = () => { oMenu.hidden = true; showDetail('menu', ''); };
  oSel.onclick = e => { e.stopPropagation(); if (oMenu.hidden && mod.outs.length && !C.auto) oMenu.hidden = false; else closeMenu(); };
  oMenu.addEventListener('mouseover', e => {
    const o = e.target.dataset && mod.outs.find(x => x.key === e.target.dataset.key);
    showDetail('menu', o ? outcomeAct(o) : '');
  });
  oMenu.addEventListener('click', e => {
    e.stopPropagation();
    const o = e.target.dataset && mod.outs.find(x => x.key === e.target.dataset.key);
    if (!o) return;
    closeMenu();
    if (C.auto) return;
    mod.hold = o.best.e.hold; syncHold();
    doAction(o.best, o.label);
  });
  document.addEventListener('click', () => { if (!oMenu.hidden) closeMenu(); });
  for (const b of document.querySelectorAll('#r-holds [data-hold]')) b.onclick = () => { mod.hold = +b.dataset.hold; syncHold(); };
  function syncHold() { for (const b of document.querySelectorAll('#r-holds [data-hold]')) b.setAttribute('aria-pressed', +b.dataset.hold === mod.hold); }
  syncHold();

  // The morph graph (Graph): the book's morphs as nodes (dashed: one only reached, not in the book;
  // grey: dissolved), an arrow for each morph it can become (and a loop for turning the other way
  // in place) with its best action: the most reliable, then the name box kept, the shortest hold.
  // Blue keeps the turning direction, orange turns the other way; faded under 3 of 4 trials,
  // dashed and fainter under 3 of 4 but at least half. Spots as a ship turning + sees them (one
  // turning - takes the mirror image), the followed ship's own as it sees them. Its node is
  // ringed yellow; its arrows are thicker, and a click on one does its action.
  const gc = $('r-graph');
  let graphSig = '', graphHits = [], graphHover = -1;   // (graphHover: the edge under the pointer)
  let graphCh = null;   // the mask channel the graph shows (null: all; a click on G / R / B in its key)
  $('r-graph-btn').onclick = () => {
    graphOn = !graphOn; $('r-graph-btn').setAttribute('aria-pressed', graphOn);
    gc.hidden = mode !== 'model' || !graphOn; graphSig = '';
  };
  // ch (a mask channel): the same arrows (each the best of all channels), each with that channel's
  // own best action instead (chNotBest: another channel's does better) or none (chNone)
  function graphEdges(B, ch = null) {
    const edges = [];
    const better = (r, o) => { if (!o) return true; const i = r.findIndex((v, j) => v !== o.r[j]); return i >= 0 && r[i] < o.r[i]; };   // (lexicographic: the first that differs)
    for (const [name, m] of Object.entries(B.morphs)) {
      const best = {}, bestC = {};
      for (const e of Models.allActions(m)) {
        if (e.p < 0.5) continue;
        const x = /^into (\S+)(?: \((same turn|reversed)\))?$/.exec(e.effect), sp = Models.splitOf(e.effect), rp = Models.replicateOf(e.effect);
        const to = e.effect === 'dissolve' ? 'dissolved' : e.effect === 'reverse' ? name : x ? x[1] : rp ? rp : sp && !sp.includes('???') ? e.effect : null;
        if (!to || !(sp || /^[a-z]+$/i.test(to))) continue;
        const keeps = !(e.effect === 'reverse' || (x && x[2] === 'reversed'));
        const k = (rp ? 'rep:' : '') + to + (keeps ? '' : '-rev'), r = [-Math.round(e.p * 20), e.always_named ? 0 : 1, e.hold, e.unnamed_max || 0];
        const ed = { k, from: name, to, keeps, e, r, split: sp, rep: rp, twin: !!x && Models.isWhole(B, x[1]) && !Models.isWhole(B, name) };
        if (better(r, best[k])) best[k] = ed;
        if (ch && e.ch === ch && better(r, bestC[k])) bestC[k] = ed;
      }
      if (ch) for (const k of Object.keys(best)) {
        const o = best[k], c = bestC[k];
        best[k] = c ? { ...c, bestE: o.e, chNotBest: better(o.r, c) } : { ...o, bestE: o.e, chNone: true };
      }
      // a split: an arrow to each morph it splits into (magenta); a replication through a short-lived
      // morph (VT049W: gyrans): the action's arrow to it (blue from a single, magenta from a twin: one
      // gyrans or two), and that morph's own green arrow to the copies (once; not an action of its own)
      for (const ed of Object.values(best)) {
        if (ed.split) for (const n of [...new Set(ed.split)]) edges.push({ ...ed, to: n });
        else if (ed.rep && ed.e.via) {
          const two = Models.isWhole(B, name);
          edges.push({ ...ed, to: ed.e.via, rep: null, viaOf: ed.rep, split: two ? [ed.e.via, ed.e.via] : null, keeps: true });
          if (!edges.some(x => x.passive && x.from === ed.e.via && x.to === ed.rep)) edges.push({ ...ed, from: ed.e.via, owner: name, passive: true });
        } else edges.push(ed);
      }
    }
    return edges;
  }
  function drawGraph() {
    if (gc.hidden) return;
    const B = book(), { name } = shipModel();
    const sig = [C.speciesOf(C.pattern), B ? 1 : 0, name, mod.sign, mod.spinSg, gc.clientWidth, gc.clientHeight, C.auto ? 1 : 0, graphHover, graphCh].join();
    if (sig === graphSig) return;
    graphSig = sig;
    const dpr = window.devicePixelRatio || 1, W = gc.clientWidth, H = gc.clientHeight;
    gc.width = Math.round(W * dpr); gc.height = Math.round(H * dpr);
    const g = gc.getContext('2d');
    g.setTransform(dpr, 0, 0, dpr, 0, 0);
    g.clearRect(0, 0, W, H);
    graphHits = [];
    g.font = '13px system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'middle';
    if (!B) { g.fillStyle = '#bbb'; g.fillText(B === undefined ? 'loading the morph model…' : `no morph model for ${C.speciesOf(C.pattern)} yet`, W / 2, H / 2); return; }
    // the mask channels its actions were tried with (G, R, B), each in its colour; the one shown
    const chs = [...new Set(Object.values(B.morphs).flatMap(m => m.actions.map(e => e.ch)))].sort((a, b) => C.CH_NAMES.indexOf(a) - C.CH_NAMES.indexOf(b));
    const gch = chs.length > 1 && chs.includes(graphCh) ? graphCh : null;
    const edges = graphEdges(B, gch);
    const order = C.morphs().map(p => C.morphName(p).replace(/ \u{1F4A5}$/u, '')), at = n => { const i = order.indexOf(n); return i < 0 ? 999 : i; };
    const names = Object.keys(B.morphs);
    for (const e of edges) if (e.to !== 'dissolved' && !names.includes(e.to)) names.push(e.to);
    names.sort((a, b) => at(a) - at(b) || a.localeCompare(b));   // (round the ellipse as the Morph list, the others after)
    // the hub (the morph replications go through, VT049W gyrans): in the middle, off the ellipse
    const hub = B.hub && names.includes(B.hub) ? B.hub : null;
    if (hub) names.splice(names.indexOf(hub), 1);
    // nodes on a circle, dissolved in the middle
    const cx = W / 2, cy = H / 2 + 8, rx = W * 0.38, ry = H * 0.37, pos = {};   // an ellipse, as wide as the view allows
    names.forEach((n, i) => { const a = -Math.PI / 2 + 2 * Math.PI * (i + 0.5) / names.length; pos[n] = [cx + rx * Math.cos(a), cy + ry * Math.sin(a)]; });
    pos.dissolved = hub ? [cx, cy + ry * 0.6] : [cx, cy];
    if (hub) { pos[hub] = [cx, cy]; names.push(hub); }   // (drawn with the others; placed in the middle)
    const NR = 30;   // a node's radius, px
    // node shapes: a circle (a single), a rectangle as wide as its name (a twin or triple), a
    // diamond (dissolved); rim(n, u): from its centre to its border along the unit vector u
    g.font = 'italic 600 13px system-ui, sans-serif';
    const box = {};
    for (const n of [...names, 'dissolved'])
      box[n] = n === 'dissolved' ? { kind: 'diamond', d: NR * 0.95 }
             : Models.isWhole(B, n) ? { kind: 'rect', hw: Math.max(NR, g.measureText(n).width / 2 + 10), hh: NR * 0.72 }
             : { kind: 'circle', r: NR };
    const rim = (n, u) => {
      const b = box[n], ax = Math.abs(u[0]) || 1e-9, ay = Math.abs(u[1]) || 1e-9;
      return b.kind === 'rect' ? Math.min(b.hw / ax, b.hh / ay) : b.kind === 'diamond' ? b.d / (ax + ay) : b.r;
    };
    g.fillStyle = '#ddd'; g.font = '600 13px system-ui, sans-serif';
    g.fillText(`Model graph: ${C.speciesName(C.pattern)}`, W / 2, 14);   // (its key: under the Model row, lab.html)
    // its key: the channels, each clickable (more than one): all, or one channel's graph
    g.font = 'bold 12px system-ui, sans-serif';
    const opts = chs.length > 1 ? ['all', ...chs] : chs;
    const parts = ['mask channels modelled: ', ...opts.flatMap((c, i) => i ? [' · ', c] : [c])];
    const widths = parts.map(t => g.measureText(t).width), tot = widths.reduce((a, b) => a + b, 0);
    let px0 = W / 2 - tot / 2;
    g.textAlign = 'left';
    parts.forEach((t, i) => {
      const ci = C.CH_NAMES.indexOf(t), opt = i > 0 && opts.includes(t), on = opts.length > 1 && opt && (t === 'all' ? !gch : t === gch);
      g.fillStyle = ci >= 0 ? C.CH_COLORS[ci] : opt ? '#ddd' : '#aaa';
      g.font = opt ? 'bold 12px system-ui, sans-serif' : '12px system-ui, sans-serif';
      g.globalAlpha = opts.length > 1 && opt && !on ? 0.55 : 1;
      g.fillText(t, px0, 30);
      if (on) g.fillRect(px0, 38, widths[i], 2);   // (the one shown: underlined)
      g.globalAlpha = 1;
      if (opts.length > 1 && opt) graphHits.push({ rect: [px0 - 4, 20, px0 + widths[i] + 4, 42], chSel: t === 'all' ? null : t,
        detail: t === 'all' ? 'the graph by all mask channels: each arrow its best action' : `the graph by ${t} masks: each arrow's best ${t} action (dimmed: another channel's does better; faint: no ${t} action does it)` });
      px0 += widths[i];
    });
    g.textAlign = 'center';
    g.font = '11px system-ui, sans-serif'; g.fillStyle = '#999';
    // edges
    // An arrow for each edge (its whole arrow the hover and click target, details on hover); the
    // followed ship's own, clickable, thick; the one under the pointer haloed white
    const pairN = {};
    edges.forEach((ed, idx) => {
      const src = ed.owner || ed.from, own = !ed.passive && !ed.chNone && src === name, a = Models.seen(B, B.morphs[src], ed.e, own ? shipSg(B.morphs[src]) : 1);
      const label = `${ed.e.ch} ${a.spot}, ${ed.e.hold / 30} s · ${ed.e.p.toFixed(2)}`;
      const col = ed.split ? SPLIT : ed.rep || ed.twin ? TWIN : ed.keeps ? SAME : FLIP;   // (green: clones or multiplies; magenta: splits or reduces)
      const click = own && !C.auto, hot = idx === graphHover, sure = ed.e.p >= 0.75;
      const path = new Path2D(), pts = [];
      let head = null;
      const P = pos[ed.from], Q = pos[ed.to];
      if (ed.from === ed.to) {   // a loop: turns the other way in place
        const out = [P[0] - cx, P[1] - cy], L = Math.hypot(...out) || 1, u = [out[0] / L, out[1] / L];
        const c = [P[0] + u[0] * (rim(ed.from, u) + 18), P[1] + u[1] * (rim(ed.from, u) + 18)];
        path.arc(c[0], c[1], 18, 0, 2 * Math.PI);
        for (let k = 0; k <= 24; k++) pts.push([c[0] + 18 * Math.cos(k * Math.PI / 12), c[1] + 18 * Math.sin(k * Math.PI / 12)]);
      } else {
        const k = ed.from + '>' + ed.to, i = pairN[k] = (pairN[k] || 0) + 1;
        const d = [Q[0] - P[0], Q[1] - P[1]], L = Math.hypot(...d), u = [d[0] / L, d[1] / L], n = [-u[1], u[0]];
        const bend = ed.to === 'dissolved' ? 0 : 30 + 46 * (i - 1);
        const m = [(P[0] + Q[0]) / 2 + n[0] * bend, (P[1] + Q[1]) / 2 + n[1] * bend];
        const r0 = rim(ed.from, u), r1 = rim(ed.to, u) + 4;
        const s0 = [P[0] + u[0] * r0 + n[0] * bend * 0.3, P[1] + u[1] * r0 + n[1] * bend * 0.3];
        const e0 = [Q[0] - u[0] * r1 + n[0] * bend * 0.3, Q[1] - u[1] * r1 + n[1] * bend * 0.3];
        path.moveTo(...s0); path.quadraticCurveTo(m[0], m[1], ...e0);
        for (let k = 0; k <= 20; k++) {
          const t = k / 20;
          pts.push([(1 - t) ** 2 * s0[0] + 2 * t * (1 - t) * m[0] + t * t * e0[0], (1 - t) ** 2 * s0[1] + 2 * t * (1 - t) * m[1] + t * t * e0[1]]);
        }
        const t = [e0[0] - m[0], e0[1] - m[1]], tl = Math.hypot(...t) || 1, tu = [t[0] / tl, t[1] / tl], hs = click ? 1.5 : 1;   // arrowhead
        head = new Path2D();
        head.moveTo(...e0);
        head.lineTo(e0[0] - (tu[0] * 11 + tu[1] * 6) * hs, e0[1] - (tu[1] * 11 - tu[0] * 6) * hs);
        head.lineTo(e0[0] - (tu[0] * 11 - tu[1] * 6) * hs, e0[1] - (tu[1] * 11 + tu[0] * 6) * hs);
        head.closePath();
      }
      g.lineCap = 'round'; g.lineJoin = 'round';
      if (hot) {   // the halo
        g.setLineDash([]); g.globalAlpha = 0.75; g.strokeStyle = '#fff'; g.fillStyle = '#fff';
        g.lineWidth = (click ? 8 : 2.5) + 7; g.stroke(path);
        if (head) { g.lineWidth = 7; g.stroke(head); g.fill(head); }
      }
      const al = C.allowedMorphs(C.speciesOf(C.pattern)), dim = al && [ed.from, ed.to].some(n => n !== 'dissolved' && !al.includes(n)) ? 0.4 : 1;
      const chDim = ed.chNone ? 0.15 : ed.chNotBest ? 0.4 : 1;   // (one channel's graph: another channel's action does it better, or only another's does it)
      g.globalAlpha = (sure || hot ? 1 : 0.55) * (hot ? Math.max(chDim, 0.6) : dim * chDim);   // (to or from a morph the game doesn't detect: dimmed)
      const lw = click ? 8 : 2.5;
      g.setLineDash(sure ? [] : click ? [16, 12] : [6, 4]);   // (less reliable: dashed, the thick sparser)
      if (!sure) g.lineCap = 'butt';   // (round ends would close the gaps)
      g.strokeStyle = col; g.fillStyle = col; g.lineWidth = lw;
      g.stroke(path);
      g.setLineDash([]);
      if (head) g.fill(head);
      g.globalAlpha = 1; g.lineCap = 'butt';
      const key = Models.outcomeOf(a), what = ed.rep ? `replicates into ${ed.e.copies ? ed.e.copies[0] + '–' + ed.e.copies[1] : 'several'} ${ed.rep}` : ed.split ? `splits into ${ed.split.join(' + ')}` : ed.to === 'dissolved' ? 'dissolves' : `${ed.to === ed.from ? 'stays ' + ed.from : 'becomes ' + ed.to}, ${ed.keeps ? 'keeping its turning direction' : 'turning the other way'}`;
      const bestBy = ed.bestE ? `${ed.bestE.ch} masks: p ${ed.bestE.p.toFixed(2)}, hold ${ed.bestE.hold / 30} s` : '';
      const detail = ed.passive ? `${ed.from} → replicates into ${ed.to} (short-lived; reached by ${src}'s actions)`
        : ed.chNone ? `${ed.from} → ${ed.to === 'dissolved' ? 'dissolves' : ed.to}: no ${gch} mask action does it (best with ${bestBy})`
        : detailText(B, B.morphs[src], ed.e, a.spot, own ? shipSg(B.morphs[src]) : 1) + (ed.chNotBest ? ` · not the best: ${bestBy}` : '');
      graphHits.push({ pts, edge: idx, a, own, detail,
                       why: key ? Models.words(key, B.morphs[src].base) : ed.split ? `split into ${ed.split.join(' + ')} (p ${ed.e.p.toFixed(2)})` : `${ed.to === ed.from ? 'turn the other way' : 'into ' + ed.to} (p ${ed.e.p.toFixed(2)})` });
    });
    // nodes
    for (const n of [...names, 'dissolved']) {
      if (n === 'dissolved' && !edges.some(e => e.to === 'dissolved')) continue;
      const [x, y] = pos[n], inBook = !!B.morphs[n], al = C.allowedMorphs(C.speciesOf(C.pattern));
      const extra = n !== 'dissolved' && al && !al.includes(n);   // (not one the game detects: dimmed)
      g.globalAlpha = extra ? 0.4 : 1;
      const bx = box[n];
      g.beginPath();
      if (bx.kind === 'rect') g.roundRect(x - bx.hw, y - bx.hh, 2 * bx.hw, 2 * bx.hh, 8);   // a twin or triple
      else if (bx.kind === 'diamond') { g.moveTo(x, y - bx.d); g.lineTo(x + bx.d, y); g.lineTo(x, y + bx.d); g.lineTo(x - bx.d, y); g.closePath(); }
      else g.arc(x, y, bx.r, 0, 2 * Math.PI);
      g.fillStyle = n === 'dissolved' ? '#444' : '#1d1d24'; g.fill();
      // (a source, nothing becomes it: no arrow in from another morph)
      const source = n !== 'dissolved' && !edges.some(e => e.to === n && e.from !== n);
      g.setLineDash([]);
      g.lineWidth = n === name ? 3.5 : 1.5; g.strokeStyle = n === name ? '#ffd400' : '#aaa'; g.stroke();
      g.setLineDash([]);
      g.fillStyle = '#fff'; g.font = (n === 'dissolved' ? '' : 'italic ') + '600 13px system-ui, sans-serif';
      g.fillText(n, x, y);
      const pat = n === 'dissolved' ? null : C.morphs().find(p => p.id === `${C.speciesOf(C.pattern)}-${n}`)
        || C.morphs().find(p => C.morphName(p).replace(C.BOOM, '') === n);   // (by name too: VT049W-gyrans-replicating is gyrans 💥)
      const m = B.morphs[n];
      // marks under the name: rotating (↻ clockwise, ↺ anticlockwise: a mover's turning, the followed
      // ship's as it turns now; a still one's spin), oscillating (∿)
      const rot = m && m.chirality ? (n === name ? (mod.sign > 0 ? '↻ +' : '↺ −') : '↻') : m && m.spin ? (m.spin > 0 ? '↻' : '↺') : m && m.drift ? (m.omega > 0 ? '↻' : '↺') : '';
      const marks = [rot, m && m.oscillating ? '∿' : '', m && m.chaotic ? '≋' : '', m && m.unstable ? '!' : ''].filter(Boolean).join(' ');
      if (marks) { g.font = '11px system-ui, sans-serif'; g.fillStyle = '#bbb'; g.fillText(marks, x, y + 15); }
      if (n !== 'dissolved') graphHits.push({ x, y, r: bx.kind === 'rect' ? bx.hw : NR, node: pat ? pat.id : null, isNode: true, detail: morphText(B, n, source).replace(pat ? '' : ' · click to reset with it', ''),
        tip: `${n}${Models.isWhole(B, n) ? ' (a twin or triple)' : ''}${m && m.oscillating ? `, oscillating (heading wobble ${m.wobble_deg}°, every ${m.wobble_period} steps; mass pulse ${Math.round(m.mass_pulse * 100)}%)` : ''}${m && m.spin ? `, spinning ${m.spin > 0 ? 'clockwise' : 'anticlockwise'}` : ''}${m && m.drift ? ', drifting (a slow spin)' : ''}${m && m.chaotic ? `, chaotic: ${m.chaotic}` : ''}${m && m.unstable ? `, unstable: ${m.unstable}` : ''}${m && m.no_maneuver ? ' (not for maneuvers in the game)' : ''}${inBook ? '' : ', no model of its own'}: click to reset the world with it` });
      g.globalAlpha = 1;
    }
    g.font = '11px system-ui, sans-serif'; g.fillStyle = '#999';
  }
  // what is under the pointer: a morph (its circle first), else the nearest arrow within 7 px
  const segDist = (x, y, p, q) => {
    const d = [q[0] - p[0], q[1] - p[1]], L2 = d[0] * d[0] + d[1] * d[1] || 1;
    const t = Math.max(0, Math.min(1, ((x - p[0]) * d[0] + (y - p[1]) * d[1]) / L2));
    return Math.hypot(x - p[0] - t * d[0], y - p[1] - t * d[1]);
  };
  const graphHit = e => {
    const r = gc.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    const node = graphHits.find(h => (h.isNode && Math.hypot(x - h.x, y - h.y) <= h.r) || (h.rect && x >= h.rect[0] && x <= h.rect[2] && y >= h.rect[1] && y <= h.rect[3]));
    if (node) return node;
    let best = null, bd = 7;
    for (const h of graphHits) if (h.pts) for (let k = 1; k < h.pts.length; k++) {
      const dd = segDist(x, y, h.pts[k - 1], h.pts[k]);
      if (dd < bd) { bd = dd; best = h; }
    }
    return best;
  };
  gc.addEventListener('pointerleave', () => { showDetail('graph', ''); if (graphHover >= 0) { graphHover = -1; drawGraph(); } });
  gc.addEventListener('pointermove', e => {   // an arrow: highlighted, its details after the Maneuver menu; a morph: its tooltip; a hand over a clickable one
    const h = graphHit(e), hv = h && h.pts ? h.edge : -1;
    if (hv !== graphHover) { graphHover = hv; drawGraph(); }
    showDetail('graph', h ? h.detail || '' : '');   // (an arrow's action, or a morph's characteristics)
    if (gc.title) gc.title = '';
    gc.style.cursor = h && (h.node || h.rect || (h.own && !C.auto)) ? 'pointer' : '';
  });
  gc.addEventListener('click', e => {
    const h = graphHit(e);
    if (h && h.rect) { graphCh = h.chSel; graphSig = ''; drawGraph(); return; }   // a mask channel: its graph
    if (h && h.node) { selectMorph(h.node); doReset(); graphSig = ''; return; }   // a morph: reset the world with it
    if (h && h.own && !C.auto && C.focus) { mod.hold = h.a.e.hold; syncHold(); doAction(h.a, h.why); graphSig = ''; }
  });

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
