"use strict";
// ============================================================================================
//  CARL maze demo, WebGL edition.
//
//  Same demo as CARL/maze_playground.html and the same trained policy; what changed is where
//  the Lenia simulation runs. The CPU version stepped the automaton in JS with an FFT-based
//  convolution, swept the whole board each step to locate the soliton, and painted the canvas
//  a pixel at a time. All of that is now GPU work in glsim.js. This file keeps the parts that
//  are genuinely sequential or genuinely CPU-bound: maze generation, episode bookkeeping, the
//  overlay drawing, the UI, and the policy network itself (onnxruntime-web, WASM).
//
//  The step order is deliberately identical to the CPU demo's agentStep(), because the policy
//  is sensitive to it: act on the current board, then step, then locate the soliton, then judge
//  the episode. What differs is only that locating the soliton and cropping its neighbourhood
//  happen on the GPU and come back together in a single readback.
// ============================================================================================

// ====================================================================================
//  CONFIG -- locked to the canonical direction run (models/meta_direction.json)
// ====================================================================================
const CFG = {
  K: 4,                       // frame stack
  R: 18,                      // kernel radius (fixed -- all 48 training solitons share it)
  netSize: 96,                 // the model was trained on 96x96 toroidal grids with no mazes --
                                // every inference call crops a 96x96 toroidal window centered on
                                // the soliton's CoM out of the (larger) maze board, so the model
                                // always sees input matching its training distribution regardless
                                // of the selected board size, and runs faster to boot.
  dt: 0.1,
  actionValue: [0.3, -0.3, 0.0],  // add / remove / no-op (output channel order)
  actionRadius: 7,
  windowSize: 4,               // steps of CoM history kept for the "current direction" arrow
  pinnedTime: 0.005,           // 50/10000 -- see meta_direction.json for why this is pinned
  costMin: 0, costMax: 5,      // true action_cost_range from meta_direction.json -- the model
                                // normalizes action_cost against this, so it must match training,
                                // independent of whatever range the UI slider exposes. The slider
                                // deliberately runs past it (to 10), which normalizes to 2.0 rather
                                // than clamping at 1.0: clamping would make the top half of the
                                // slider inert, and rescaling to the slider's own range would
                                // silently change what every existing cost value means.
  massDeathFraction: 0.3,      // soliton counts as "dead" below this fraction of its spawn mass
  massExplodeLimit: 1500,      // ...and as blown up above this absolute mass
  modelUrl: 'models/agent_direction.onnx',
  solitonsUrl: 'models/solitons_direction.json',
  thumbDir: 'assets/solitons/',        // one PNG per allowed rule, named "<rule name>.png"
  // The curated subset of update rules this demo exposes, out of the 48 the maze agent was
  // trained across. The picker renders these in *this* order, five per row, so the leading five
  // are the ones a first-time visitor sees on the top row: the solitons that survive the maze
  // most reliably.
  allowedRuleNames: [
    'rule73_mu0.2250_s0.0250_R18', 'rule74_mu0.2300_s0.0350_R18', 'rule67_mu0.2800_s0.0410_R18',
    'rule72_mu0.1650_s0.0200_R18', 'rule57_mu0.2350_s0.0360_R18',
    'rule75_mu0.2600_s0.0420_R18', 'rule55_mu0.2750_s0.0450_R18', 'rule58_mu0.2250_s0.0340_R18',
    'rule60_mu0.2850_s0.0520_R18', 'rule23_mu0.2650_s0.0390_R18', 'rule81_mu0.2700_s0.0480_R18',
    'rule68_mu0.3050_s0.0580_R18', 'rule76_mu0.2850_s0.0290_R18', 'rule33_mu0.2100_s0.0250_R18',
    'rule37_mu0.2650_s0.0330_R18', 'rule7_mu0.3200_s0.0560_R18', 'rule21_mu0.2400_s0.0240_R18',
    'rule46_mu0.2500_s0.0270_R18', 'rule63_mu0.3200_s0.0660_R18', 'rule61_mu0.2850_s0.0610_R18',
  ],
  defaultRuleName: 'rule73_mu0.2250_s0.0250_R18',  // pre-selected on load
  easyRuleCount: 5,            // leading entries flagged as easy-to-steer -- one picker row
  actionFadeSeconds: 0.85,     // how long an intervention marker takes to fade out, in sim time
};
const MA = CFG.actionValue[0]; // action magnitude (0.3)

// Board colors. The state ramp itself is plain black->white and lives in draw.glsl; these two
// are passed in as uniforms.
const WALL_RGB = [120, 110, 90], TARGET_RGB = [51, 153, 255];

// ====================================================================================
//  Simulation state
//
//  Note what is *not* here any more: the board array, the frame stack, the kernel taps and the
//  FFT scratch buffers. The board lives in GPU textures and the frame stack is a ring of them;
//  the only board-sized array left on this side is the wall mask, which is generated once per
//  maze and uploaded.
// ====================================================================================
let H = 150, W = 150, N = H * W;
let mu = 0.24, sig = 0.024, betas = [1.0, 0.5];
let wall = new Uint8Array(0);
let maze = { start: [0, 0], target: [0, 0], proximity: 12 };
let targetSigma = 3;                          // goal blob width, display only
let dir = [0, 1];                             // current target direction (dy,dx), unit length
let comHistory = [[0, 0], [0, 0], [0, 0], [0, 0]];
let lastCoM = null, lastDist = null, initialMass = 0;
// The most recent crop readback, channel-packed (frame k in channel k). Held by reference: the
// engine reuses one buffer, and it is always consumed into a tensor before the next readback
// overwrites it.
let lastCrop = null;
let steps = 0, actions = 0, goalReached = false, solitonDead = false;
let courseChanges = 0;          // target direction changes the user made this episode
let failReason = null;          // 'exploded' | 'dissolved' once solitonDead is set
let running = false, session = null, busy = false, lastMs = 0, lastAction = null, lastQ = null;
// Recent interventions, newest last: {r,c,sign,step}. Each marker fades out over
// CFG.actionFadeSeconds of *simulation* time, so the trail reads the same at any sim speed
// (and freezes mid-fade rather than vanishing when the run is paused).
let actionTrail = [];
let bank = [], currentIndex = 0;
let mazeEnabled = true;
// Who intervenes on the soliton. When true CARL's policy is never queried: the sim still
// advances, but the only actions on the board are the ones the user clicks in.
let humanActs = false;
// The user's queued action, {gx,gy,sign}, or null. At most one is ever held: a step consumes it
// exactly where agentAct() would have run, so a human turn and a CARL turn are the same turn.
let pendingAction = null;
// The 90deg rotation the soliton spawns with. Rolled once per maze, not per spawn, so Respawn
// re-runs the same starting configuration instead of quietly changing the soliton's heading.
let spawnRotation = 0;
let sps = 60, stepAcc = 0, lastT = 0, measSps = 0, rateSteps = 0, rateTime = 0;
// Higher than the CPU demo's cap: with the convolution on the GPU a frame can absorb far more
// steps before it stops keeping up, and the measured rate readout shows where the real ceiling
// lands on a given machine.
const MAX_STEPS_PER_FRAME = 16;
const MAX_TRAIL = 48;            // hard cap on the intervention trail (the fade usually ends it first)

const cv = document.getElementById('maze');
const ov = document.getElementById('overlay'), octx = ov.getContext('2d');
const $ = id => document.getElementById(id);

// ====================================================================================
//  Toroidal helper (the CoM itself is computed on the GPU -- see com.glsl)
// ====================================================================================
function toroidalDelta(ny, nx, oy, ox) {
  let dy = ny - oy, dx = nx - ox;
  if (dy > H / 2) dy -= H; else if (dy < -H / 2) dy += H;
  if (dx > W / 2) dx -= W; else if (dx < -W / 2) dx += W;
  return [dy, dx];
}

// ====================================================================================
//  Maze generation -- randomized-DFS spanning tree + loops, carved onto the pixel grid.
//  Direct port of scripts/interactive_maze_agent_demo.py :: generate_maze_dfs
// ====================================================================================
function generateMaze(bh, bw, cw, ww, edgeWall, enableWalls) {
  const cellSize = cw + ww, border = edgeWall;
  let cols = Math.floor((bw - 2 * border + ww) / cellSize), rows = Math.floor((bh - 2 * border + ww) / cellSize);
  cols = Math.max(2, cols); rows = Math.max(2, rows);
  // cw is a floor, not a target: whatever room floor()'ing cols/rows above leaves unused gets
  // spent widening the corridors, not the outer border -- otherwise a thick ww/edgeWall can
  // strand a large leftover remainder as an oversized edge margin instead of maze structure.
  const cwX = Math.max(cw, Math.floor((bw - 2 * border - (cols - 1) * ww) / cols));
  const cwY = Math.max(cw, Math.floor((bh - 2 * border - (rows - 1) * ww) / rows));
  const cellSizeX = cwX + ww, cellSizeY = cwY + ww;
  const offX = Math.max(0, border + Math.floor((bw - 2 * border - (cols * cellSizeX - ww)) / 2));
  const offY = Math.max(0, border + Math.floor((bh - 2 * border - (rows * cellSizeY - ww)) / 2));

  const visited = Array.from({ length: rows }, () => new Array(cols).fill(false));
  const edges = Array.from({ length: rows }, () => Array.from({ length: cols }, () => new Set()));
  const dirs4 = [[-1, 0], [1, 0], [0, -1], [0, 1]];
  const dstack = [[0, 0]]; visited[0][0] = true;
  while (dstack.length) {
    const [r, c] = dstack[dstack.length - 1];
    const nbrs = [];
    for (const [dr, dc] of dirs4) {
      const nr = r + dr, nc = c + dc;
      if (nr >= 0 && nr < rows && nc >= 0 && nc < cols && !visited[nr][nc]) nbrs.push([nr, nc]);
    }
    if (nbrs.length) {
      const [nr, nc] = nbrs[(Math.random() * nbrs.length) | 0];
      edges[r][c].add(nr + ',' + nc); edges[nr][nc].add(r + ',' + c);
      visited[nr][nc] = true; dstack.push([nr, nc]);
    } else dstack.pop();
  }
  const numExtra = Math.max(1, Math.floor(rows * cols / 5));
  for (let i = 0; i < numExtra; i++) {
    const r = (Math.random() * rows) | 0, c = (Math.random() * cols) | 0;
    const dd = dirs4.slice().sort(() => Math.random() - 0.5);
    for (const [dr, dc] of dd) {
      const nr = r + dr, nc = c + dc;
      if (nr >= 0 && nr < rows && nc >= 0 && nc < cols && !edges[r][c].has(nr + ',' + nc)) {
        edges[r][c].add(nr + ',' + nc); edges[nr][nc].add(r + ',' + c); break;
      }
    }
  }

  const wallArr = new Uint8Array(bh * bw).fill(enableWalls ? 1 : 0);
  if (enableWalls) {
    const carve = (y0, x0, y1, x1) => {
      y0 = Math.max(0, y0); x0 = Math.max(0, x0); y1 = Math.min(bh, y1); x1 = Math.min(bw, x1);
      for (let y = y0; y < y1; y++) { const row = y * bw; for (let x = x0; x < x1; x++) wallArr[row + x] = 0; }
    };
    for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
      const top = offY + r * cellSizeY, left = offX + c * cellSizeX;
      carve(top, left, top + cwY, left + cwX);
    }
    for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
      for (const key of edges[r][c]) {
        const [nr, nc] = key.split(',').map(Number);
        if (nr * cols + nc <= r * cols + c) continue;
        const top1 = offY + r * cellSizeY, left1 = offX + c * cellSizeX;
        const top2 = offY + nr * cellSizeY, left2 = offX + nc * cellSizeX;
        carve(Math.min(top1, top2), Math.min(left1, left2), Math.max(top1, top2) + cwY, Math.max(left1, left2) + cwX);
      }
    }
  }

  function bfsDist(sr, sc) {
    const dist = Array.from({ length: rows }, () => new Array(cols).fill(-1));
    dist[sr][sc] = 0; const q = [[sr, sc]]; let qi = 0;
    while (qi < q.length) {
      const [r, c] = q[qi++];
      for (const key of edges[r][c]) {
        const [nr, nc] = key.split(',').map(Number);
        if (dist[nr][nc] === -1) { dist[nr][nc] = dist[r][c] + 1; q.push([nr, nc]); }
      }
    }
    return dist;
  }
  const cellCenter = (r, c) => [offY + r * cellSizeY + ((cwY / 2) | 0), offX + c * cellSizeX + ((cwX / 2) | 0)];

  const d0 = bfsDist(0, 0);
  let far = [0, 0], maxD = -1;
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) if (d0[r][c] > maxD) { maxD = d0[r][c]; far = [r, c]; }
  const d1 = bfsDist(far[0], far[1]);
  let end = [0, 0], maxD2 = -1;
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) if (d1[r][c] > maxD2) { maxD2 = d1[r][c]; end = [r, c]; }

  return { wall: wallArr, start: cellCenter(far[0], far[1]), target: cellCenter(end[0], end[1]) };
}

function boardParams(size) {
  // Corridor width (cw) is floored by the widest soliton actually used here (rule74 is 40px
  // wide, models/solitons_direction.json) -- 44px leaves it a ~4px margin per side, as far as
  // this can shrink without the soliton clipping wall edges on every turn. edgeWall matches ww
  // so the outer perimeter reads as the same thickness as the inner separator walls.
  // The 150 board is the exception: a cw floor of 44 only fits 2x2 cells there, so it drops to
  // 38 to fit 3x3 instead. 100 is tighter again, and 2x2 is all it can hold; its floor is 36
  // because that is what makes the 2x2 land on the intended 9px border.
  return { cw: size <= 100 ? 36 : size <= 150 ? 38 : 44, ww: 9, edgeWall: 9 };
}

// ====================================================================================
//  Soliton placement: exact 90deg rotation (no interpolation) + place its own center of
//  mass at the maze start, then zero out anything that lands on a wall.
// ====================================================================================
function rotate90(src, h, w, k) {
  k = ((k % 4) + 4) % 4;
  if (k === 0) return { data: src.slice(), h, w };
  if (k === 2) {
    const out = new Float32Array(h * w);
    for (let i = 0; i < h * w; i++) out[h * w - 1 - i] = src[i];
    return { data: out, h, w };
  }
  const nh = w, nw = h, out = new Float32Array(nh * nw);
  for (let y = 0; y < h; y++) for (let x = 0; x < w; x++) {
    const v = src[y * w + x]; if (v === 0) continue;
    let ny, nx;
    if (k === 1) { ny = x; nx = h - 1 - y; } else { ny = w - 1 - x; nx = y; }
    out[ny * nw + nx] = v;
  }
  return { data: out, h: nh, w: nw };
}

function placeSoliton(entry) {
  mu = entry.mu; sig = entry.sigma; betas = entry.betas.slice();
  SimGL.setRule({ mu, sigma: sig, betas, R: CFG.R, dt: CFG.dt });

  const rot = rotate90(entry.flat, entry.h, entry.w, spawnRotation);
  let sy = 0, sx = 0, sm = 0;
  for (let y = 0; y < rot.h; y++) for (let x = 0; x < rot.w; x++) {
    const v = rot.data[y * rot.w + x]; if (v > 0) { sm += v; sy += v * y; sx += v * x; }
  }
  const cy = sy / sm, cx = sx / sm;

  // Built once on the CPU and uploaded; from here on the board only exists on the GPU.
  const arr = new Float32Array(N);
  const [tr, tc] = maze.start;
  const oy = Math.round(tr - cy), ox = Math.round(tc - cx);
  for (let y = 0; y < rot.h; y++) {
    const gy = ((oy + y) % H + H) % H, grow = gy * W;
    for (let x = 0; x < rot.w; x++) {
      const v = rot.data[y * rot.w + x]; if (v <= 0) continue;
      const gx = ((ox + x) % W + W) % W;
      arr[grow + gx] = v;
    }
  }
  for (let i = 0; i < N; i++) if (wall[i]) arr[i] = 0;    // applyWallCollision, at spawn
  SimGL.uploadState(arr);

  steps = 0; actions = 0; courseChanges = 0; goalReached = false; solitonDead = false; failReason = null;
  lastAction = null; lastQ = null;
  actionTrail.length = 0;
  pendingAction = null;

  // One analysis pass with no step behind it, so the spawn CoM and the policy's first input
  // window come from the same place every later step gets them from.
  SimGL.prime();
  const rb = SimGL.readback();
  lastCrop = rb.crop;
  lastCoM = rb.valid ? [rb.row, rb.col, rb.mass] : null;
  initialMass = lastCoM ? lastCoM[2] : 0;
  comHistory = Array.from({ length: CFG.windowSize }, () => lastCoM ? [lastCoM[0], lastCoM[1]] : [tr, tc]);
  if (lastCoM) {
    const ddy = maze.target[0] - lastCoM[0], ddx = maze.target[1] - lastCoM[1];
    lastDist = Math.hypot(ddy, ddx);
    const n = lastDist || 1; dir = [ddy / n, ddx / n];
  }
  render();
}

function newMaze() {
  spawnRotation = (Math.random() * 4) | 0;
  const size = +$('board').value;
  const { cw, ww, edgeWall } = boardParams(size);
  maze = generateMaze(H, W, cw, ww, edgeWall, mazeEnabled);
  maze.proximity = cw * 0.3;
  wall = maze.wall;
  SimGL.setWall(wall);
  targetSigma = Math.max(3, size * 0.02);
  placeSoliton(bank[currentIndex]);
}
function setBoard(size) {
  H = W = size; N = H * W;
  SimGL.setBoard(W, H);
  boardReady = true;
  newMaze();
}
function respawnCurrentSoliton() { placeSoliton(bank[currentIndex]); }

// ====================================================================================
//  Agent step
// ====================================================================================
// Crops a CFG.netSize toroidal window out of each stacked frame, centered on the soliton's
// current CoM. The crop itself happens in crop.glsl; this is the CPU-side copy of its origin
// math, used to map the policy's per-cell output back to board coordinates and to bound what
// the user is allowed to click in human mode.
function agentWindowOrigin(cy, cx) {
  const half = CFG.netSize >> 1;
  return [Math.round(cy) - half, Math.round(cx) - half];
}
function buildContext() {
  return Float32Array.from([
    CFG.pinnedTime, dir[0], dir[1],
    (actionCost - CFG.costMin) / (CFG.costMax - CFG.costMin),
  ]);
}

// CARL's half of a turn: pick a spot + sign from the policy. Unlike the CPU version this only
// *returns* the intervention -- applying it is a GPU pass inside SimGL.step().
async function agentAct() {
  const t0 = performance.now();
  const S = CFG.netSize, SS = S * S;
  // The crop arrives channel-packed (frame k in channel k of one RGBA texture); the model wants
  // [1,K,S,S], i.e. frame-major. This de-interleave is the whole cost of the new input path.
  const data = new Float32Array(CFG.K * SS);
  for (let k = 0; k < CFG.K; k++) {
    const base = k * SS;
    for (let i = 0; i < SS; i++) data[base + i] = lastCrop[i * 4 + k];
  }
  const feeds = {
    state:   new ort.Tensor('float32', data, [1, CFG.K, S, S]),
    context: new ort.Tensor('float32', buildContext(), [1, 4]),
  };
  const out = await session.run(feeds);
  const q = out.q.data;
  lastMs = performance.now() - t0;
  let best = 0, bv = q[0];
  for (let i = 1; i < q.length; i++) if (q[i] > bv) { bv = q[i]; best = i; }
  lastQ = q;
  const at = Math.floor(best / SS), pos = best % SS, lr = Math.floor(pos / S), lc = pos % S;
  const [oy, ox] = agentWindowOrigin(lastCoM[0], lastCoM[1]);
  const r = ((oy + lr) % H + H) % H, c = ((ox + lc) % W + W) % W;
  const sign = Math.sign(CFG.actionValue[at]);
  lastAction = { r, c, sign };
  if (sign === 0) return null;
  actions++;
  actionTrail.push({ r, c, sign, step: steps });
  if (actionTrail.length > MAX_TRAIL) actionTrail.shift();
  return { x: c, y: r, delta: sign * MA, radius: CFG.actionRadius };
}

// The user's half of a turn: hand over the one queued action, if any, and clear the queue. Runs
// at the same point in the step as agentAct(), so the action is registered on this step and the
// board is free to take the next one.
function takePendingAction() {
  const a = pendingAction;
  if (!a) return null;
  pendingAction = null;
  actions++;
  actionTrail.push({ r: ((a.gy % H) + H) % H, c: ((a.gx % W) + W) % W, sign: a.sign, step: steps });
  if (actionTrail.length > MAX_TRAIL) actionTrail.shift();
  return { x: a.gx, y: a.gy, delta: a.sign * MA, radius: CFG.actionRadius };
}

async function agentStep() {
  if (!lastCoM) return;
  let action = null;
  if (humanActs) { lastMs = 0; action = takePendingAction(); }   // no inference -- readout shows "--"
  else if (!session) return;
  else action = await agentAct();

  // Action, step, wall mask, CoM reduction and the next input crop, queued back to back with no
  // synchronization; the single readback below is the only point where the CPU waits on the GPU.
  SimGL.step(action);
  const rb = SimGL.readback();
  lastCrop = rb.crop;
  steps++;
  finishStep(rb);
}

// Bookkeeping shared by both actors: take the soliton's freshly computed center of mass and
// decide whether the episode has ended.
function finishStep(rb) {
  if (!rb.valid) {
    lastCoM = null;
    if (!goalReached && !solitonDead) { solitonDead = true; failReason = 'dissolved'; running = false; }
    return;
  }
  lastCoM = [rb.row, rb.col, rb.mass];
  comHistory.push([rb.row, rb.col]); if (comHistory.length > CFG.windowSize) comHistory.shift();
  lastDist = Math.hypot(maze.target[0] - rb.row, maze.target[1] - rb.col);
  if (!goalReached && !solitonDead) {
    if (lastDist < maze.proximity) goalReached = true;
    else if (rb.mass > CFG.massExplodeLimit) { solitonDead = true; failReason = 'exploded'; }
    else if (rb.mass < CFG.massDeathFraction * initialMass) { solitonDead = true; failReason = 'dissolved'; }
    if (goalReached || solitonDead) running = false;
  }
}

// ====================================================================================
//  Rendering -- the board itself is a shader pass (draw.glsl); everything below is the
//  overlay canvas sitting on top of it, unchanged from the CPU demo.
// ====================================================================================
function drawArrow(c, x0, y0, dx, dy, len, color) {
  const n = Math.hypot(dx, dy); if (n < 1e-6) return;
  const ux = dx / n, uy = dy / n, x1 = x0 + ux * len, y1 = y0 + uy * len;
  c.strokeStyle = color; c.fillStyle = color; c.lineWidth = Math.max(2, len * 0.07);
  c.beginPath(); c.moveTo(x0, y0); c.lineTo(x1, y1); c.stroke();
  const hs = len * 0.24, ang = Math.atan2(uy, ux);
  c.beginPath(); c.moveTo(x1, y1);
  c.lineTo(x1 - hs * Math.cos(ang - 0.4), y1 - hs * Math.sin(ang - 0.4));
  c.lineTo(x1 - hs * Math.cos(ang + 0.4), y1 - hs * Math.sin(ang + 0.4));
  c.closePath(); c.fill();
}
// The green/red discs stamped on the board each step. Each one lingers and fades out over
// CFG.actionFadeSeconds of sim time instead of blinking for a single frame, so the pattern of
// interventions stays readable at 20+ steps/sec.
function drawActionTrail(sc) {
  if (!actionTrail.length) return;
  const fadeSteps = Math.max(4, CFG.actionFadeSeconds * sps);
  const rad = CFG.actionRadius * sc, newest = actionTrail[actionTrail.length - 1].step;
  while (actionTrail.length && steps - actionTrail[0].step > fadeSteps) actionTrail.shift();
  for (const a of actionTrail) {
    const age = (steps - a.step) / fadeSteps;
    if (age > 1) continue;
    const k = 1 - age, fade = k * k;                 // ease-out: bright for a beat, then a long tail
    const px = (a.c + 0.5) * sc, py = (a.r + 0.5) * sc;
    const rgb = a.sign > 0 ? '83,185,121' : '222,73,104';
    octx.beginPath(); octx.arc(px, py, rad, 0, 7);
    octx.fillStyle = `rgba(${rgb},${(0.30 * fade).toFixed(3)})`; octx.fill();
    octx.strokeStyle = `rgba(${rgb},${(0.85 * fade + 0.06 * k).toFixed(3)})`;
    octx.lineWidth = (a.step === newest ? 3 : 2) * Math.max(1, sc * 0.5);
    octx.stroke();
    if (a.step === newest) {                     // halo so the current intervention stands out
      octx.beginPath(); octx.arc(px, py, rad + 3, 0, 7);
      octx.strokeStyle = `rgba(${rgb},.35)`; octx.lineWidth = 1.5; octx.stroke();
    }
  }
}
// The queued-but-not-yet-applied action, as a hollow dashed ring.
function drawPendingAction(sc) {
  if (!pendingAction) return;
  const { gx, gy, sign } = pendingAction, rad = CFG.actionRadius * sc;
  const rgb = sign > 0 ? '83,185,121' : '222,73,104';
  octx.beginPath(); octx.arc((gx + 0.5) * sc, (gy + 0.5) * sc, rad, 0, 7);
  octx.setLineDash([rad * 0.55, rad * 0.45]);
  octx.strokeStyle = `rgba(${rgb},.9)`;
  octx.lineWidth = 2 * Math.max(1, sc * 0.5);
  octx.stroke();
  octx.setLineDash([]);
}
let boardReady = false;   // render() can be reached from UI handlers before the first setBoard()
function render() {
  if (!boardReady) return;
  SimGL.draw({
    targetRow: maze.target[0], targetCol: maze.target[1], targetSigma,
    wallRGB: WALL_RGB, targetRGB: TARGET_RGB,
  });

  const r = cv.getBoundingClientRect();
  if (ov.width !== Math.round(r.width)) { ov.width = Math.round(r.width) || W; ov.height = Math.round(r.height) || H; }
  octx.clearRect(0, 0, ov.width, ov.height);
  const sc = ov.width / W;

  if (lastCoM) {
    const cx = lastCoM[1] * sc, cy = lastCoM[0] * sc, alen = Math.max(H, W) * 0.16 * sc;
    const cs = getComputedStyle(document.documentElement);
    // The target direction is CARL's instruction. With CARL idle nothing consumes it, so the blue
    // arrow would be a claim about an agent that isn't acting -- hide it while the user has the board.
    if (!humanActs) drawArrow(octx, cx, cy, dir[1], dir[0], alen, cs.getPropertyValue('--target').trim());
    const old = comHistory[0];
    const [vdy, vdx] = toroidalDelta(lastCoM[0], lastCoM[1], old[0], old[1]);
    drawArrow(octx, cx, cy, vdx, vdy, alen, cs.getPropertyValue('--current').trim());
  }
  drawActionTrail(sc);
  drawPendingAction(sc);

  $('r-dist').textContent = lastDist != null ? lastDist.toFixed(1) : '–';
  $('r-dist').className = (lastDist != null && lastDist < maze.proximity * 1.5) ? 'ok' : '';
  $('r-mass').textContent = lastCoM ? lastCoM[2].toFixed(0) : '0';
  $('r-step').textContent = steps;
  $('r-acts').textContent = actions;
  $('r-course').textContent = courseChanges;
  $('r-ms').textContent = lastMs ? lastMs.toFixed(1) + 'ms' : '–';
  $('r-rate').textContent = measSps ? measSps.toFixed(0) + '/s' : '–';
  updateResult();
}
// Episode outcome, shown over the board once the run stops. Rebuilt only when the outcome
// actually changes -- render() runs on every frame.
let resultKey = '';
function updateResult() {
  const el = $('result');
  const key = goalReached ? `win:${steps}:${actions}:${courseChanges}`
            : solitonDead ? `lose:${failReason}:${steps}:${actions}:${courseChanges}` : '';
  if (key === resultKey) return;
  resultKey = key;
  if (!key) { el.hidden = true; el.innerHTML = ''; return; }
  const plural = (n, w) => `<span class="s">${n} ${w}${n === 1 ? '' : 's'}</span>`;
  if (goalReached) {
    el.className = 'result win';
    el.innerHTML = `<b>Goal reached</b>${plural(steps, 'episode step')} · ${plural(actions, 'agent action')}`
                 + ` · ${plural(courseChanges, 'target direction change')}`;
  } else {
    el.className = 'result lose';
    const why = failReason === 'exploded' ? 'Fail: soliton exploded' : 'Fail: soliton dissolved';
    el.innerHTML = `<b>${why}</b>after ${plural(steps, 'episode step')} · ${plural(actions, 'agent action')}`
                 + ` · ${plural(courseChanges, 'target direction change')}`;
  }
  el.hidden = false;
}

// ====================================================================================
//  Main loop (accumulator scheduler)
// ====================================================================================
async function loop() {
  if (!running) return;
  if (!busy) {
    busy = true;
    const now = performance.now();
    if (!lastT) lastT = now;
    let elapsed = (now - lastT) / 1000; if (elapsed > 0.25) elapsed = 0.25;
    lastT = now;
    stepAcc += elapsed * sps;
    let want = Math.floor(stepAcc); stepAcc -= want;
    const n = Math.min(want, MAX_STEPS_PER_FRAME);
    try {
      for (let k = 0; k < n; k++) {
        if (goalReached || solitonDead) { running = false; break; }
        await agentStep();
      }
    } finally { busy = false; }
    rateSteps += n; rateTime += elapsed;
    if (rateTime >= 0.5) { measSps = rateSteps / rateTime; rateSteps = 0; rateTime = 0; }
    if (n > 0) render();
  }
  requestAnimationFrame(loop);
}
function setRunning(v) {
  running = v && !goalReached && !solitonDead;
  $('play').textContent = running ? 'Pause' : 'Play';
  $('play').classList.toggle('on', running);
  if (running) { lastT = 0; stepAcc = 0; rateSteps = 0; rateTime = 0; loop(); }
  else render();
}

// ====================================================================================
//  Input: click-to-steer, arrow keys, sliders/buttons/selects
// ====================================================================================
let actionCost = 2.5;
// Every user-initiated steer goes through here, so a single place counts course changes.
// Re-issuing the direction the agent is already chasing (holding an arrow key down, clicking
// straight ahead) isn't a change of course, so it isn't counted.
function steerTo(dy, dx) {
  const n = Math.hypot(dy, dx); if (n < 1e-6) return;
  const ny = dy / n, nx = dx / n;
  if (ny * dir[0] + nx * dir[1] < 0.9999) courseChanges++;
  dir = [ny, nx];
  render();
}
cv.addEventListener('click', e => {
  if (humanActs || !lastCoM) return;   // in human mode the left button acts instead; arrows still steer
  const r = cv.getBoundingClientRect();
  const px = (e.clientX - r.left) / r.width * W, py = (e.clientY - r.top) / r.height * H;
  steerTo(py - lastCoM[0], px - lastCoM[1]);
});
window.addEventListener('keydown', e => {
  const tag = document.activeElement && document.activeElement.tagName;
  if (tag === 'SELECT' || tag === 'INPUT') return;
  const map = { ArrowUp: [-1, 0], ArrowDown: [1, 0], ArrowLeft: [0, -1], ArrowRight: [0, 1] };
  const v = map[e.key];
  if (v) { steerTo(v[0], v[1]); e.preventDefault(); }
});
$('play').addEventListener('click', () => setRunning(!running));
// Step always leaves the run paused and the board one step further on, so it reads as a scrub
// rather than a nudge to a still-running sim. If the loop happens to be mid-batch when the click
// lands (busy), that in-flight step is the advance -- stepping again here would double it.
async function stepOnce() {
  setRunning(false);
  if (busy || goalReached || solitonDead) return;
  busy = true;
  try { await agentStep(); } finally { busy = false; }
  render();
}
$('step1').addEventListener('click', stepOnce);

// ------------------------------------------------------------------------------------
//  Human actor: the same intervention CARL makes -- same +-MA over CFG.actionRadius, same
//  disc, same "actions" tally -- and, like CARL, at most one per step. A click queues the
//  action rather than applying it on the spot; the next step stamps it in at exactly the
//  point in the turn CARL's own action would land.
// ------------------------------------------------------------------------------------
// CARL's reach: the policy sees a CFG.netSize window centered on the soliton and returns one q
// value per cell of exactly that window -- so a spot outside it is not a move CARL could make.
// Same bound, same origin math, for the user.
function inAgentReach(gy, gx) {
  if (!lastCoM) return false;
  const S = CFG.netSize, [oy, ox] = agentWindowOrigin(lastCoM[0], lastCoM[1]);
  return ((gy - oy) % H + H) % H < S && ((gx - ox) % W + W) % W < S;
}
let toastTimer = 0;
function showToast(msg) {
  const el = $('toast');
  el.textContent = msg;
  el.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => el.classList.remove('show'), 2000);
}
function humanAct(e, sign) {
  if (goalReached || solitonDead) return;
  if (pendingAction || busy) return;              // one action per step, no queue-jumping
  const r = cv.getBoundingClientRect();
  const gx = Math.floor((e.clientX - r.left) / r.width * W), gy = Math.floor((e.clientY - r.top) / r.height * H);
  if (!inAgentReach(gy, gx)) {
    showToast(`Out of reach — only actions close to the soliton are.`);
    return;
  }
  pendingAction = { gx, gy, sign };
  if (running) render();                        // show it queued
  else stepOnce();                              // paused: act and advance in one motion
}
cv.addEventListener('pointerdown', e => {
  if (!humanActs) return;
  if (e.button === 0) humanAct(e, 1);
  else if (e.button === 2) humanAct(e, -1);
  else if (e.button === 1) stepOnce();
  else return;
  e.preventDefault();               // no text selection, and no middle-click autoscroll
});
// Right- and middle-click are actions here, not browser gestures.
cv.addEventListener('contextmenu', e => { if (humanActs) e.preventDefault(); });
cv.addEventListener('auxclick', e => { if (humanActs) e.preventDefault(); });
$('actor').addEventListener('click', () => {
  humanActs = !humanActs;
  $('actor').textContent = humanActs ? 'You act' : 'CARL acts';
  $('actor').classList.toggle('on', humanActs);
  $('human-hint').hidden = !humanActs;
  $('legend-target').hidden = humanActs;         // no blue arrow on the board, no key for it
  cv.classList.toggle('human', humanActs);
  lastMs = 0; lastQ = null; lastAction = null; pendingAction = null;
  render();
});
$('respawn').addEventListener('click', respawnCurrentSoliton);
$('newmaze').addEventListener('click', newMaze);
$('shuffle').addEventListener('click', () => {
  currentIndex = (Math.random() * bank.length) | 0;
  updateSolitonSelection();
  respawnCurrentSoliton();
});
$('board').addEventListener('change', e => { const was = running; setRunning(false); setBoard(+e.target.value); if (was) setRunning(true); });
$('maze-toggle').addEventListener('click', () => {
  mazeEnabled = !mazeEnabled;
  $('maze-toggle').textContent = mazeEnabled ? 'On' : 'Off';
  $('maze-toggle').classList.toggle('on', mazeEnabled);
  const was = running; setRunning(false); newMaze(); if (was) setRunning(true);
});
function setActionCost(v) {
  actionCost = v;
  const el = $('v-cost');
  el.textContent = actionCost.toFixed(2);
  const ood = actionCost > CFG.costMax;
  el.classList.toggle('ood', ood);
  el.title = ood ? `Beyond the 0-${CFG.costMax} range CARL was trained on -- it is extrapolating here.` : '';
}
$('cost').addEventListener('input', e => setActionCost(+e.target.value));
$('spd').addEventListener('input', e => { sps = +e.target.value; $('v-spd').textContent = sps + '/s'; });

// ====================================================================================
//  Theme
// ====================================================================================
function applyTheme(t) {
  document.documentElement.setAttribute('data-theme', t);
  $('theme-toggle').textContent = t === 'dark' ? '☀️' : '🌙';
  $('theme-toggle').title = t === 'dark' ? 'Switch to light mode' : 'Switch to dark mode';
  try { localStorage.setItem('carl-theme', t); } catch (e) {}
}
function initTheme() {
  applyTheme(document.documentElement.getAttribute('data-theme') === 'dark' ? 'dark' : 'light');
  $('theme-toggle').addEventListener('click', () => {
    applyTheme(document.documentElement.getAttribute('data-theme') === 'dark' ? 'light' : 'dark');
    render();
  });
}

// ====================================================================================
//  Boot
// ====================================================================================
// Stops the run without going through setRunning(), which would try to render -- fail() has to
// stay callable before the board exists, e.g. when GPU init itself is what failed.
function fail(msg) {
  $('tip').innerHTML = msg;
  running = false;
  $('play').textContent = 'Play';
  $('play').classList.remove('on');
}
// No inline fallback soliton here, unlike the CPU demo: this version fetches its shaders as
// well as its model, so it cannot run from file:// under any circumstances, and the offline
// case the fallback existed for cannot arise.
async function loadSolitonBank() {
  let raw = null;
  try {
    const res = await fetch(CFG.solitonsUrl);
    if (res.ok) raw = await res.json();
  } catch (e) { /* reported below */ }
  if (!raw) { fail('Could not load the soliton bank — serve this folder over http(s), not file://.'); return false; }
  // The picker's order is CFG.allowedRuleNames', not the JSON's -- the leading entries are the
  // ones we want on the top row, so sort the fetched entries back into that order.
  raw = raw.filter(e => CFG.allowedRuleNames.includes(e.name))
           .sort((a, b) => CFG.allowedRuleNames.indexOf(a.name) - CFG.allowedRuleNames.indexOf(b.name));
  bank = raw.filter(e => Math.max(e.h, e.w) <= 80).map(e => ({
    name: e.name, mu: e.mu, sigma: e.sigma, betas: e.betas, h: e.h, w: e.w,
    flat: Float32Array.from(e.state.flat()),
  }));
  if (!bank.length) { fail('The soliton bank loaded but contained no usable entries.'); return false; }

  const grid = $('soliton-grid'); grid.innerHTML = '';
  bank.forEach((b, i) => {
    const btn = document.createElement('button');
    const easy = i < CFG.easyRuleCount;
    btn.type = 'button'; btn.className = 'soliton-thumb' + (easy ? ' easy' : ''); btn.dataset.index = i;
    btn.title = `${b.name} · β=[${b.betas.join(',')}] · μ${b.mu.toFixed(3)} · σ${b.sigma.toFixed(3)}`
              + (easy ? ' · ★ easy to control' : '');
    const img = document.createElement('img');
    img.src = CFG.thumbDir + encodeURIComponent(b.name) + '.png';
    img.alt = b.name; img.loading = 'lazy';
    // Replace just the image on error -- the star badge below is a sibling, not part of it.
    img.addEventListener('error', () => { img.replaceWith(Object.assign(document.createElement('span'), { textContent: '?' })); });
    btn.appendChild(img);
    if (easy) {
      const star = document.createElement('span');
      star.className = 'star'; star.textContent = '★'; star.setAttribute('aria-hidden', 'true');
      btn.appendChild(star);
    }
    btn.addEventListener('click', () => { currentIndex = i; updateSolitonSelection(); respawnCurrentSoliton(); });
    grid.appendChild(btn);
  });
  const def = bank.findIndex(b => b.name === CFG.defaultRuleName);
  currentIndex = def >= 0 ? def : 0;
  updateSolitonSelection();
  return true;
}
function updateSolitonSelection() {
  $('soliton-grid').querySelectorAll('.soliton-thumb').forEach(btn => {
    btn.classList.toggle('sel', +btn.dataset.index === currentIndex);
  });
}
// Forced to WASM (CPU) -- WebGPU is disabled regardless of browser support, carried over from
// the CPU demo. Note that moving inference to the WebGPU backend would not remove the readback
// in agentStep(): the simulation above lives in a WebGL2 context, and the two APIs do not share
// GPU memory, so the crop would still have to travel through the CPU to reach it.
async function loadModel() {
  ort.env.wasm.wasmPaths = 'https://cdn.jsdelivr.net/npm/onnxruntime-web@1.20.1/dist/';
  // Multi-threaded WASM needs SharedArrayBuffer, which needs the cross-origin isolation the
  // coi-serviceworker provides; on a first visit (before the service worker has taken over)
  // crossOriginIsolated is still false, so fall back to a single thread rather than let
  // onnxruntime-web throw. Overridable via ?threads=N for re-tuning on other machines; the
  // thread pool is sized once at wasm-module init, so each value needs a fresh page load.
  const threadsOverride = parseInt(new URLSearchParams(location.search).get('threads'), 10);
  ort.env.wasm.numThreads = threadsOverride > 0 ? threadsOverride
    : (window.crossOriginIsolated ? Math.min(navigator.hardwareConcurrency || 6, 6) : 1);

  // WASM pays a one-time JIT cost on its very first run. Absorb that here, before the model is
  // used for real, instead of letting it land on the user's first real simulation step.
  const warmup = async s => {
    try {
      const S = CFG.netSize, SS = S * S;
      const data = new Float32Array(CFG.K * SS);
      for (let k = 0; k < CFG.K; k++) {
        const base = k * SS;
        for (let i = 0; i < SS; i++) data[base + i] = lastCrop[i * 4 + k];
      }
      const feeds = { state: new ort.Tensor('float32', data, [1, CFG.K, S, S]), context: new ort.Tensor('float32', buildContext(), [1, 4]) };
      await s.run(feeds);
    } catch (e) { /* best-effort -- worst case the first real step pays this cost instead */ }
  };

  try {
    session = await ort.InferenceSession.create(CFG.modelUrl, { executionProviders: ['wasm'] });
    await warmup(session);
  } catch (e) {
    console.error(e);
    fail('Model failed to load — serve over http(s), not file://.');
    return;
  }
  render();
}

(async function init() {
  initTheme();
  setActionCost(+$('cost').value);
  try {
    await SimGL.init(cv, CFG.netSize);
  } catch (e) {
    console.error(e);
    fail('<b>GPU init failed:</b> ' + e.message);
    return;
  }
  if (!await loadSolitonBank()) return;
  setBoard(+$('board').value);
  setRunning(true);
  loadModel();
})();
