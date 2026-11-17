"use strict";
// ============================================================================================
//  World draw ([world_draw:on]; the bullet book's Nebula and Antimatter): three static fields,
//  R, G and B, that the player paints straight onto the world, in either view. On the GPU they
//  are one texture (glsim.js setField), standing in for the bullets, the bots' bullets and the
//  dots while a field type is in use: never stepped, never worn down, but each field masks
//  (Nebula) or erases (Antimatter) its own channel of the creatures, as the bullets do.
//
//  The fields live here, on the CPU (buf: W x H x 4, rgb used), and the rows a stroke touched are
//  uploaded once a frame (tick). Tools, by the R G B patches (core.js fieldTool): 0 paintbrush
//  (raises the field to the brush's coverage, anti-aliased), 1 soft brush (watercolour: a blurred
//  disc, 1 at its middle fading smoothly to 0 at its edge; also raises the field to it, so no
//  build-up), 2 eraser (lowers it to 1 - coverage). Width: 2 to 32 cells (WIDTHS); the soft brush
//  and the eraser twice that (TOOL_SCALE).
//  Undo: one level, the last stroke (or Clear); undo again redoes it.
// ============================================================================================

const Field = (function () {
  const C = Core;
  const W = C.BW, H = C.BH;
  const WIDTHS = [2, 4, 8, 16, 32];  // paintbrush widths, cells (the Width slider)
  const TOOL_SCALE = [1, 2, 2];      // per tool (paintbrush, soft brush, eraser): its width, times the slider's
  let buf = new Float32Array(W * H * 4);
  let undoBuf = null, undone = false;   // the field before the last stroke (after an undo: the stroke's)
  let widthIdx = 1;
  let stroke = null;                    // { ch, tool, at: [col, row] }
  const dirty = new Uint8Array(H);      // rows to upload
  let anyDirty = false;
  const hooks = { change() {} };        // the undo state changed

  const widthOf = tool => WIDTHS[widthIdx] * (TOOL_SCALE[tool] ?? 1);
  const radius = tool => widthOf(tool) / 2;
  // One dab of tool at (cx, cy) (cell-index units: cell i's centre at i), on channel ch, wrapping
  // round the torus.
  function dab(cx, cy, ch, tool) {
    const r = radius(tool), R = Math.ceil(r + 1), x0 = Math.round(cx), y0 = Math.round(cy);
    for (let dy = -R; dy <= R; dy++) {
      const y = (((y0 + dy) % H) + H) % H;
      let hit = false;
      for (let dx = -R; dx <= R; dx++) {
        const d = Math.hypot(x0 + dx - cx, y0 + dy - cy);
        let k;
        if (tool === 1) { const t = d / (r + 0.5); if (t >= 1) continue; k = 1 - t * t * (3 - 2 * t); }   // the blurred disc: smoothstep from 1 (middle) to 0 (edge)
        else { k = Math.min(1, r + 0.5 - d); if (k <= 0) continue; }   // coverage: a 1-cell anti-aliased edge
        const i = (y * W + (((x0 + dx) % W) + W) % W) * 4 + ch, v = buf[i];
        buf[i] = tool === 2 ? Math.min(v, 1 - k) : Math.max(v, k);
        hit = true;
      }
      if (hit) { dirty[y] = 1; anyDirty = true; }
    }
  }
  function snapshot() { undoBuf = buf.slice(); undone = false; hooks.change(); }
  function markAll() { dirty.fill(1); anyDirty = true; }

  // A stroke: begin at a press, move with the drag (dabs every third of the radius along the way; the soft brush every half cell),
  // end at the release
  function begin(at, ch, tool) {
    snapshot();
    stroke = { ch, tool, at };
    dab(at[0], at[1], ch, tool);
  }
  function move(at) {
    if (!stroke) return;
    const a = stroke.at, dx = C.wrapD(at[0] - a[0], W), dy = C.wrapD(at[1] - a[1], H);
    const gap = stroke.tool === 1 ? 0.5 : Math.max(0.5, radius(stroke.tool) / 3);   // (the soft brush closer: no ridges along it)
    const n = Math.max(1, Math.ceil(Math.hypot(dx, dy) / gap));
    for (let k = 1; k <= n; k++) {
      const p = [a[0] + dx * k / n, a[1] + dy * k / n];
      dab(p[0], p[1], stroke.ch, stroke.tool);
    }
    stroke.at = at;
  }
  function end() { stroke = null; flush(); }

  // Each frame: the rows touched since are uploaded
  function tick() { flush(); }
  function flush() {
    if (!anyDirty) return;
    anyDirty = false;
    for (let y = 0; y < H;) {
      if (!dirty[y]) { y++; continue; }
      let e = y;
      while (e < H && dirty[e]) e++;
      SimGL.putFieldRows(buf, y, e);
      y = e;
    }
    dirty.fill(0);
    let any = false;
    for (let i = 0; i < buf.length && !any; i++) any = buf[i] > 0;
    SimGL.setFieldEmpty(!any);   // (an empty field: no effect on the creatures, no cost)
  }

  // The scrolling world's seam (game.js scrollFrame, each frame; null: none): the columns and rows
  // within band of it, and all it swept over since the last frame, emptied in every channel, so what
  // was painted drifts away for good instead of coming round to the ship's front (not undoable)
  let lastSeam = null;
  function clearSeam(col, row, band) {
    if (col === null) { lastSeam = null; return; }
    const prev = lastSeam || [col, row];
    lastSeam = [col, row];
    const dc = C.wrapD(col - prev[0], W), dr = C.wrapD(row - prev[1], H);
    const cols = new Set(), rows = [];
    for (let x = Math.floor(Math.min(prev[0], prev[0] + dc) - band); x <= Math.ceil(Math.max(prev[0], prev[0] + dc) + band); x++) cols.add(((x % W) + W) % W);
    for (let y = Math.floor(Math.min(prev[1], prev[1] + dr) - band); y <= Math.ceil(Math.max(prev[1], prev[1] + dr) + band); y++) rows.push(((y % H) + H) % H);
    const wipe = (y, x) => { const i = (y * W + x) * 4; if (buf[i] || buf[i + 1] || buf[i + 2]) { buf[i] = buf[i + 1] = buf[i + 2] = 0; dirty[y] = 1; anyDirty = true; } };
    for (let y = 0; y < H; y++) for (const x of cols) wipe(y, x);
    for (const y of rows) for (let x = 0; x < W; x++) wipe(y, x);
  }

  // Clear: channel ch emptied (undoable, as a stroke)
  function clear(ch) {
    stroke = null;
    snapshot();
    for (let i = ch; i < buf.length; i += 4) buf[i] = 0;
    markAll(); flush();
  }
  // Undo the last stroke (or Clear); again: redo it. False if there is none.
  function undo() {
    if (!undoBuf) return false;
    stroke = null;
    [buf, undoBuf] = [undoBuf, buf];
    undone = !undone;
    markAll(); flush(); hooks.change();
    return true;
  }
  // Everything emptied, no undo ([clear_world], [world_draw:off])
  function reset() {
    stroke = null; undoBuf = null; undone = false;
    buf.fill(0);
    markAll(); flush(); hooks.change();
  }
  function setWidth(i) { widthIdx = Math.max(0, Math.min(WIDTHS.length - 1, i)); }

  return { begin, move, end, tick, clear, clearSeam, undo, reset, setWidth, widthOf, hooks, WIDTHS,
           get width() { return WIDTHS[widthIdx]; }, get widthIdx() { return widthIdx; },
           get canUndo() { return !!undoBuf; }, get undone() { return undone; }, get active() { return !!stroke; },
           get tool() { return stroke ? stroke.tool : null; } };   // the held stroke's tool (null: none)
})();
