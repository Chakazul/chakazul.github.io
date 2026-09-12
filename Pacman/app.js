"use strict";
// ============================================================================================
//  CARL maze demo, WebGL edition.
//
//  Same demo and trained policy as CARL/maze_playground.html; what changed is where the Lenia
//  simulation runs. The CPU version stepped the automaton in JS (FFT convolution, a full-board
//  sweep for the soliton, per-pixel canvas paint) -- all of that is now GPU work in glsim.js.
//  This file keeps only what's genuinely sequential or CPU-bound: maze generation, episode
//  bookkeeping, the overlay, the UI, and the policy network (onnxruntime-web, WASM).
//
//  Step order matches the CPU demo's agentStep() exactly, since the policy is sensitive to it
//  (act, step, locate, judge) -- what differs is that locating the soliton and cropping its
//  neighbourhood happen on the GPU and come back together in one readback.
// ============================================================================================

// ====================================================================================
//  URL params and feature switches -- e.g. ?sound=0&dots=0&net=96
// ====================================================================================
// "0"/"false" turns a switch off; anything else, or absence, leaves the default.
function boolParam(name, def) {
  const v = new URLSearchParams(location.search).get(name);
  return v === null ? def : v !== '0' && v.toLowerCase() !== 'false';
}
function intParam(name, def) {
  const v = parseInt(new URLSearchParams(location.search).get(name), 10);
  return v > 0 ? v : def;
}
// A whole percentage, 0-100, as a fraction. Unlike intParam this accepts 0 as a real setting.
function pctParam(name, def) {
  const v = parseInt(new URLSearchParams(location.search).get(name), 10);
  return v >= 0 && v <= 100 ? v / 100 : def;
}
const SOUND_ENABLED = boolParam('sound', true);
// Pellet (channel-2) field. Off: placeDots() never sets a rule, so channel 2's sim pass and its
// contribution to the drawn board are both skipped -- only Pac-Man's channel moves.
const DOTS_ENABLED = boolParam('dots', true);
// Ghost (channel-3) field, same shape as the dots switch. Off: no ghost sim pass and no eat check
// -- i.e. also the switch for "nothing can kill Pac-Man but himself".
const GHOSTS_ENABLED = boolParam('ghost', true);
// Dying still plays out (jingle, pause, respawn) but never costs a life -- for testing deep levels
// without a game over cutting the run short.
const GOD_MODE = boolParam('god', false);
// WebGPU by default for the policy net (faster than WASM on the mobile devices this was tuned
// on), falling back to wasm per-op. Doesn't remove the CPU roundtrip in agentStep() -- the sim
// runs in a separate WebGL2 context sharing no memory with WebGPU -- but speeds up the net's own
// conv work. `?ep=wasm` forces CPU-only, for re-comparing on a new device.
const EXECUTION_PROVIDERS = new URLSearchParams(location.search).get('ep') === 'wasm'
  ? ['wasm'] : ['webgpu', 'wasm'];

// ====================================================================================
//  CONFIG -- locked to the canonical direction run (models/meta_direction.json)
// ====================================================================================
const CFG = {
  K: 4,                       // frame stack
  R: 18,                      // kernel radius (fixed -- all 48 training solitons share it)
  // The model was trained on 96x96 toroidal grids with no mazes; every inference crops a 96x96
  // toroidal window centred on the CoM out of the larger board, so input always matches training
  // regardless of board size. ?net=N overrides this for experiments (fully convolutional, so it
  // still runs, but off-distribution) and also sets the crop size the one readback carries.
  netSize: intParam('net', 96),
  dt: 0.1,
  // Per-channel dt multipliers -- Pac-Man's and the ghosts' own speed, independent of each other
  // and of the stationary dots channel. Each layers under glsim.js's frightened-window multiplier
  // rather than replacing it (PACMAN_FRIGHTENED_SPEED, GHOST_FRIGHTENED_SPEED); 1 leaves it as is.
  channel1Speed: 1.5,
  channel3Speed: 1.5,
  actionValue: [0.3, -0.3, 0.0],  // add / remove / no-op (output channel order)
  actionRadius: 7,
  windowSize: 4,               // steps of CoM history kept for the "current direction" arrow
  pinnedTime: 0.005,           // 50/10000 -- see meta_direction.json for why this is pinned
  // The true action_cost_range from meta_direction.json -- the model normalizes against this, so
  // it must match training regardless of the UI slider's own range. The slider deliberately runs
  // past costMax (to 10): clamping would make its top half inert, and rescaling to the slider's
  // own range would silently change what every existing cost value means.
  costMin: 0, costMax: 5,
  massDeathFraction: 0.3,      // soliton counts as "dead" below this fraction of its spawn mass
  // ...and as blown up above this absolute mass. Lowered from the CPU demo's 1500, which windowing
  // made unreachable: an undisturbed rule74 soliton peaks at 437 (146% of spawn) and one fed
  // CARL-sized mass every step plateaus at 679, so 1500 never fired. 1000 clears both peaks with
  // room to spare while staying well under a real runaway (6% of the 128-window's capacity).
  massExplodeLimit: 1000,
  modelUrl: 'models/agent_direction.onnx',
  solitonsUrl: 'models/solitons_direction.json',
  thumbDir: 'assets/solitons/',        // one PNG per allowed rule, named "<rule name>.png"
  // Curated subset of the 48 rules the maze agent was trained on, in picker order (five per row,
  // so the leading five are the ones a first-time visitor sees) -- these survive the maze best.
  allowedRuleNames: [
    'rule73_mu0.2250_s0.0250_R18', 'rule74_mu0.2300_s0.0350_R18', 'rule67_mu0.2800_s0.0410_R18',
    'rule72_mu0.1650_s0.0200_R18', 'rule57_mu0.2350_s0.0360_R18',
    'rule75_mu0.2600_s0.0420_R18', 'rule55_mu0.2750_s0.0450_R18', 'rule58_mu0.2250_s0.0340_R18',
    'rule60_mu0.2850_s0.0520_R18', 'rule23_mu0.2650_s0.0390_R18', 'rule81_mu0.2700_s0.0480_R18',
    'rule68_mu0.3050_s0.0580_R18', 'rule76_mu0.2850_s0.0290_R18', 'rule33_mu0.2100_s0.0250_R18',
    'rule37_mu0.2650_s0.0330_R18', 'rule7_mu0.3200_s0.0560_R18', 'rule21_mu0.2400_s0.0240_R18',
    'rule46_mu0.2500_s0.0270_R18', 'rule63_mu0.3200_s0.0660_R18',
  ],
  defaultRuleName: 'rule74_mu0.2300_s0.0350_R18',  // or rule73_mu0.2250_s0.0250_R18
  // Dots channel's rule: one soliton per '.' in the layout, free-running under its own growth,
  // never steered.
  channel2RuleName: 'rule0_mu0.3800_s0.0700_R18',
  // Half CFG.R, shrinking the dots to half size -- see resizeSoliton() in placeDots(), which
  // always resamples by exactly channel2R/R.
  channel2R: 9,
  // Ghost channel's rule: one soliton per 'M', respawned alongside Pac-Man. Kept at the native
  // CFG.R (unlike the dots) since a Lenia pattern only moves like itself at the radius it was
  // found at, and a ghost should be Pac-Man-sized. Free-running -- CARL never sees or steers it.
  //
  // Chosen for speed: a ghost that can't keep up isn't a threat. Measured in a corridor of this
  // maze, rule74 runs 0.479 px/step against rule73's 0.194 (rule73 wastes more motion on a wider
  // lean angle). Matching defaultRuleName is deliberate: it puts the ghost at exactly Pac-Man's
  // pace on the default pick, so it closes only by out-navigating him.
  channel3RuleName: 'rule74r_mu0.2300_s0.0350_R18',
  // Difficulty dial: how often a ghost's turn deliberately closes on Pac-Man rather than rolling
  // the dice. 0 = wanders, 1 = closes whenever a junction allows it. ?chase=N (0-100) to retune.
  chaseBias: pctParam('chase', 0.6),
  easyRuleCount: 5,            // leading entries flagged as easy-to-steer -- one picker row
  actionFadeSeconds: 0.85,     // how long an intervention marker takes to fade out, in sim time
  // Wall-clock (not sim) gap with nothing eaten before the eat_dot_0/1 "waka waka" loop stops.
  eatSoundGraceMs: 1000,
  // Board counts as cleared once the dots channel's live total mass (SimGL.readback().dotsMass)
  // falls below this -- not exactly 0, since a dot mid-erasure can leave a sliver its own growth
  // rule would sustain forever. An absolute value rather than a fraction of the spawn total: the
  // dots are themselves free-running solitons, so their total drifts with growth/decay rather
  // than only falling as they're eaten.
  dotsWinMass: 0.1,
  // A dot scores the step its own window (dotsites.glsl) drops below this fraction of what it
  // held at stamping. Wide margin either way (CPU-mirror measurement): an untouched or grazed dot
  // never falls below ~80% over its ~52-step breathing cycle, while a bitten one dissolves to 0
  // within a few steps.
  dotGoneFraction: 0.4,
  // How long a power pellet's frightened window stays open, wall-clock ms (player-perceived time,
  // so it doesn't scale with sim speed). A second pellet eaten mid-window resets rather than stacks.
  frightenedDurationMs: 15000,
  // Bonus fruit: a bonus dot appears on Pac-Man's spawn tile the moment the level's dots-eaten
  // fraction crosses each entry of bonusFruitThresholds -- twice a level at the defaults (30%,
  // 70%) -- worth BONUS_FRUIT_POINTS[level-1] instead of the usual 10. Left uneaten for
  // bonusFruitTimeout steps, it disappears instead (the next threshold, if any, still arms in its
  // own time -- see updateBonusFruit()). Configurable -- just change the numbers.
  bonusFruitThresholds: [0.3, 0.7],
  bonusFruitTimeout: 800,
};
// Pickup distance from Pac-Man's own centre of mass -- his kernel radius (CFG.R), since the fruit
// is picked up on contact rather than dissolved like an ordinary dot (see updateBonusFruit()).
const BONUS_FRUIT_RADIUS = CFG.R;
const MA = CFG.actionValue[0]; // action magnitude (0.3)
// Keyed by level (1-9, i.e. index level-1) rather than by which of a level's two fruits it is --
// see updateBonusFruit(). Array.from(), not a plain string split: every one of these is outside
// the BMP (a surrogate pair in UTF-16), so [...str]/split('') would cut them in half.
const BONUS_FRUIT_EMOJIS = Array.from('🍊🍎🍒🍓🍉🍭🍄🍩🍖');
const BONUS_FRUIT_POINTS = [100, 100, 100, 200, 500, 700, 1000, 2000, 5000];

// Board palette, 0-255, passed to draw.glsl as uniforms: background->soliton is the mass 0->1
// ramp, channels 2/3 layer over it as soliton2/soliton3, walls paint flat.
const COLORS = {
  soliton:    [255, 221, 51],   // Pac-Man yellow
  soliton2:   [255, 255, 255],  // ordinary dots, in the free-running channel-2 field
  pellet:     [255, 255, 0],    // power pellets, same field -- see placeDots()'s power mask
  soliton3:   [255, 40, 40],    // fallback for ghost mass with no owner recorded (arcade red)
  // One per ghost, in spawn order, cycled if there are more ghosts than colours. ghostblit.glsl
  // records which one owns each pixel, so it stays exact even where tiles overlap.
  ghosts: [
    [255,  40,  40],            // Blinky red
    [ 90, 160, 255],            // Inky blue
    [ 80, 230, 120],            // Clyde green
    [255, 170,  60],            // amber
    [124,  58, 237],            // purple
    [255, 105, 180],            // pink
    [188,   0, 211],            // violet
    [  0, 100,   0],            // dark green
    [130, 130, 130],            // grey
  ],
  wall:       [33, 33, 180],    // dark arcade blue
  background: [0, 0, 0],
  frightened: [0, 0, 139],      // every ghost's colour during a power pellet's frightened window
};

// ====================================================================================
//  Simulation state
//
//  What's not here any more: the board array, frame stack, kernel taps, FFT scratch -- the board
//  lives in GPU textures. The only board-sized array left here is the wall mask, generated once
//  per maze and uploaded.
// ====================================================================================
// 5 rows x 11 columns of maze cells; 550x250 is exactly 11:5 so both axes land on ~40px cells
// (see MAZE_LAYOUT below). No board-size control -- the layout fixes the cell count, and any
// other size just starves the corridors.
const BOARD_H = 250, BOARD_W = 550;
let H = BOARD_H, W = BOARD_W, N = H * W;
let mu = 0.24, sig = 0.024, betas = [1.0, 0.5];
let wall = new Uint8Array(0);
let maze = { start: [0, 0] };
let dir = [0, 1];                             // current target direction (dy,dx), unit length
let comHistory = [[0, 0], [0, 0], [0, 0], [0, 0]];
let lastCoM = null, initialMass = 0;
// Most recent crop readback, channel-packed. Held by reference (the engine reuses one buffer),
// always consumed into a tensor before the next readback overwrites it.
let lastCrop = null;
let steps = 0, actions = 0, solitonDead = false;
let courseChanges = 0;          // target direction changes the user made this episode
// Pending auto-respawn from a death, so a manual Restart/etc. during the gap can cancel it.
let deathTimer = 0;
let levelWon = false;
// Shown as 💛 in the title row. Reset only on a genuine new game (Restart, maze toggle, soliton
// picker, initial load, or game-over), not on an ordinary death-respawn (so they count down
// across deaths). A win carries the count into the next board and adds one, capped at MAX_LIVES.
const STARTING_LIVES = 3;
const MAX_LIVES = 5;
let lives = STARTING_LIVES;
// Level N spawns the maze's ghosts numbered 1..N, so the starting level is the starting ghost
// count (see placeGhosts()). Advances by one per win, resets on a genuine new game same as lives
// -- not on a death-respawn or a win -- so the game gets harder round over round, not per death.
// ?level=N to start elsewhere (also what a reset falls back to).
const STARTING_LEVEL = intParam('level', 3);
let level = STARTING_LEVEL;
// The maze layout's highest numbered ghost spawn (see MAZE_LAYOUT) -- clearing this level already
// has every ghost in play, so a further level would spawn nothing new. Winning it ends the run
// instead of quietly advancing to an identical board.
const MAX_LEVEL = 9;
let winTimer = 0;    // pending auto-restart from clearing the board, same shape as deathTimer
// Tracked per ghost (newGhost()'s `frightened`), not as one global mode: a pellet marks every
// not-already-frightened ghost, but one eaten and respawned ends its own window immediately while
// its packmates keep counting down theirs. frightenedTimer is the shared duration, restarted on
// every pellet eaten so a second pellet extends the window rather than stacking one behind it.
let frightenedTimer = 0;
// Set by frightenedTimer's callback (no fresh readback there to turn ghosts with) -- just flags
// the window ended. updateFrightened() reacts on the next step it notices, turning every *still*-
// frightened ghost back toward Pac-Man (one already normal from an earlier respawn is untouched).
let frightenedExpired = false;
// True while the sim holds for the eat_ghost.wav jingle -- set by handleGhostEaten(), cleared by
// the sound's own 'ended' event (unlike solitonDead/levelWon this has no timer of its own).
let ghostEatPause = false;
// Floating "200" (etc.) shown where the ghost died, live for exactly as long as ghostEatPause --
// set alongside it in handleGhostEaten(), cleared alongside it too. {r, c, text} or null.
let ghostEatPopup = null;
// Bonus fruit currently on the board, or null -- {r, c, emoji, points, spawnStep}, always stamped
// at maze.start. bonusFruitThresholdIdx is which of CFG.bonusFruitThresholds is still armed (0,
// then 1, then done for the level) -- it resets only where dotSites itself does, in placeDots(),
// so a death respawn (which leaves already-eaten dots eaten) doesn't re-arm a threshold the level
// already passed, and a fresh level's dots do get both thresholds back.
let bonusFruit = null;
let bonusFruitThresholdIdx = 0;
// True while the sim holds for the eat_fruit.wav jingle, same shape as ghostEatPause.
let fruitEatPause = false;
// Same idea as ghostEatPopup, for the bonus fruit -- shown where it was picked up.
let fruitEatPopup = null;
// Persists across wins and death-respawns, resets only on a genuine new game. 10/dot, 50/pellet,
// awarded once each (rb.dotSites vs CFG.dotGoneFraction) -- not off mass erased under Pac-Man's
// pixels (misses what a bitten dot sheds afterward) nor the channel's raw total (swings ~15 dots'
// worth as every dot breathes in phase).
let score = 0;
// One entry per placeDots() stamp, in SimGL.setDotSites() order: {r, c, points, spawnMass, eaten}.
// Rebuilt only with the dots field itself, so a death respawn keeps already-eaten dots counted.
let dotSites = [];
// 200/400/800/... per ghost eaten within one frightened window, arcade-style. Reset to 1 on a
// fresh pellet (even mid-window) and on every fresh spawn.
let ghostChainMultiplier = 1;
// eat_dot_0/1 "waka waka" loop: eatActive while alternating, eatToggle picks which plays next,
// eatDeadline (wall-clock) is when it's allowed to stop -- pushed out by CFG.eatSoundGraceMs per
// eaten dot. Checked only where a new play starts (playNextEatSound), not off a separate timer,
// so at most one sample plays at once.
let eatActive = false, eatToggle = 0, eatDeadline = 0;
let running = false, session = null, busy = false, lastMs = 0, lastAction = null, lastQ = null;
// True from playStartSound() until its jingle finishes (covers waiting for the unlocking gesture
// too). setRunning() refuses to unpause while set, so Play/Space can't cut the intro short.
let introPlaying = false;
// Recent interventions, newest last: {r,c,sign,step}. Fades over CFG.actionFadeSeconds of *sim*
// time, so the trail reads the same at any sim speed and freezes rather than vanishes when paused.
let actionTrail = [];
let bank = [], currentIndex = 0;
// Every usable JSON entry by name, whether or not the picker offers it -- the free-running
// channels look their rules up here rather than in `bank`.
let ruleBank = new Map();
const mazeEnabled = true;
// Chrome on/off: the control deck and CARL's overlay hide together, so the board reads as a game.
let showOverlay = true;
// 'sometimes' (default) queries CARL only for sometimesWindow steps after a spawn/steer, then
// goes idle -- much cheaper than 'always' since inference is the expensive part of a step, and
// looks the same whenever the player is actually engaged. 'human': CARL is never queried; the
// user's own clicks are the only actions.
let actorMode = 'sometimes';
let humanActs = false;                 // derived from actorMode === 'human', kept for readability
let sometimesWindow = intParam('sometimes', 100);  // ?sometimes=N -- steps CARL stays active for in 'sometimes' mode
let sometimesRemaining = 0;            // steps left in the current active window ('sometimes' mode only)
// Further thins out inference within an active window: 1 infers every step, N>1 infers only every
// Nth active step and skips the action in between -- the sim step, readback and rendering still
// run every step, just without a fresh intervention.
let inferenceStride = intParam('stride', 1);  // ?stride=N
let strideCounter = 0;                 // active steps left to skip before the next inference is due
// User's queued action, {gx,gy,sign}, or null. At most one held: a step consumes it exactly where
// agentAct() would run, so a human turn and a CARL turn are the same turn.
let pendingAction = null;
// 90deg spawn rotation, rolled once per maze (not per spawn) so Respawn repeats the same heading.
let spawnRotation = 0;
let sps = 60, stepAcc = 0, lastT = 0, measSps = 0, rateSteps = 0, rateTime = 0;
// A batch yields once it's spent this much of the frame, rather than a fixed step count. The
// fixed cap this replaces (16) assumed pure-GPU steps; with CARL acting, each also pays an
// inference (tens of ms on a phone), so 16 back to back could freeze the board for most of a
// second before the next repaint. Budgeting by wall-clock time instead runs however many steps
// fit, so a slow device just runs slower rather than in lurches. ?budget=N to tune per device.
const FRAME_BUDGET_MS = intParam('budget', 10);
const MAX_TRAIL = 48;            // hard cap on the intervention trail (the fade usually ends it first)

const cv = document.getElementById('maze');
const ov = document.getElementById('overlay'), octx = ov.getContext('2d');
const $ = id => document.getElementById(id);
const sndStart = new Audio('assets/sound/start.wav');
const sndDeath = new Audio('assets/sound/death_0.wav');
const sndIntermission = new Audio('assets/sound/intermission.wav');
const sndFright = new Audio('assets/sound/fright.wav');
const sndEatGhost = new Audio('assets/sound/eat_ghost.wav');
const sndEatFruit = new Audio('assets/sound/eat_fruit.wav');
const sndEat = [new Audio('assets/sound/eat_dot_0.wav'), new Audio('assets/sound/eat_dot_1.wav')];
// Registered once, not per-play: 'ended' only fires on a natural finish (never from
// stopEatingSound()'s pause()), so a one-shot listener re-added each play would pile up whenever
// a play gets interrupted mid-note (a death, a Restart) -- each stale listener left behind means
// an extra concurrent call the next time that element finishes naturally, playing two tracks at
// once. playNextEatSound() itself decides whether to keep going.
sndEat.forEach(snd => snd.addEventListener('ended', playNextEatSound));

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
//  Maze layout -- fixed, hand-authored (replacing the randomized-DFS generator the CARL demo
//  shipped with), but the same geometry: wide open *cells* separated by thin *walls* between
//  them. Written on a doubled grid, the standard maze-ASCII form: an RxC maze is (2R+1)x(2C+1)
//  characters, odd/odd indices are cells ((row,col) -> cell ((row-1)/2, (col-1)/2)), even/even
//  are corner posts (always wall), and everything else is a wall slot between two cells.
//
//  '+' and '-' are pure authoring aliases for '.' and ' ', kept so a hand-edited row still reads
//  as ASCII art. '^' is open floor too, except a ghost standing on it is sent north rather than
//  choosing for itself -- what gives a centre room a one-way exit.
//
//  '.'/'+' places a dot and 'O' a power pellet (both a channel-2 soliton, see placeDots()) on
//  either a cell slot (centred on the tile) or a wall slot (centred in the gap, like a pellet in
//  a corridor) -- ignored on a corner post, which is always wall.
//
//  On a cell slot: 'C' Pac-Man's spawn · '1'-'9' a numbered ghost's spawn (placeGhosts() spawns
//  ghosts numbered at or below the current `level`) · '#' solid cell · 'O' a power pellet (same
//  soliton/rule as a dot, just recoloured -- see placeDots()'s power mask) · anything else open
//  floor (a 'C' or digit is open floor too, just marking what spawns there).
//  On a wall slot: '#'/'|' wall · anything else an open passage -- including on the outer ring,
//  since the sim wraps toroidally regardless of walls, so an opening there is a real side tunnel
//  to the opposite edge, not a dead end.
//
//  Cell width (cw) is derived from the board so the maze always fills it, so keep cw above the
//  soliton's width (30-46px, models/solitons_direction.json) or turns clip it. The 5x11 maze on
//  the 550x250 board lands at cwX=40, cwY=39, wall thickness matching the old generator (see
//  MAZE_WW). Solitons wider than cw are clipped to fit -- see placeSoliton().
// ====================================================================================
const MAZE_WALL_CHARS = '#|';
// Headings in rotate90()'s k order, so (k+1)%4 is a right turn, (k+3)%4 a left one, and the
// difference between two headings is the quarter turns between them -- exact and rule-independent.
// The *absolute* heading a freshly stamped pattern travels is not (rule74's canonical heading is
// right, rule73's is left), so a ghost's heading is always read off its own motion, never assumed
// from its stamp rotation -- see ghostHeading().
const DIRS = [[0, 1], [1, 0], [0, -1], [-1, 0]];
// Cells a ghost may turn on -- the same "dot"/"open"/"ghost spawn"/"power pellet" characters
// doing double duty, since in this layout they land on exactly the 32 corner/junction cells.
const GHOST_TURN_CHARS = '+-O^123456789';
const GHOST_NORTH_CHAR = '^';   // the one that doesn't leave the choice open -- always sent north
// Wall thickness and outer border, matching the old generator (ww=9, edgeWall=9). At 5x5 cells
// that leaves cw=39, tighter than the widest curated soliton (46px) -- see placeSoliton()'s clip.
const MAZE_WW = 9, MAZE_EDGE = 9;

// 5 rows x 11 cols, mirrored left-to-right about the centre column: the right half is the left
// half reflected, except that 'C' is not duplicated -- the soliton spawns on the left only.
const MAZE_LAYOUT = [
  "####### ####### #######",
  "#O...9#+.......+#7...O#",
  "#.###.#.#######.#.###.#",
  "#.#+.+.+...3...+.+.+#.#",
  "#.#.###.##   ##.###.#.#",
  " +.4...+#1 ^ 2#+...5.+ ",
  "#.#.###.#######.###.#.#",
  "#.#+.+.+.. C ..+.+.+#.#",
  "#.###.#.#######.#.###.#",
  "#O...6#+...O...+#8...O#",
  "####### ####### #######",
];

function layoutMaze(layout, bh, bw, enableWalls) {
  const gh = layout.length, gw = Math.max(...layout.map(r => r.length));
  const rows = (gh - 1) >> 1, cols = (gw - 1) >> 1;
  if (!(gh & 1) || !(gw & 1) || rows < 1 || cols < 1) {
    console.warn(`maze layout should be (2R+1)x(2C+1) characters; got ${gh}x${gw}`);
  }
  // Anything past the end of a ragged row, or outside the grid entirely, reads as wall.
  const at = (r, c) => (r >= 0 && r < gh && c >= 0 && c < layout[r].length ? layout[r][c] : '#');
  const isWall = ch => MAZE_WALL_CHARS.includes(ch);

  const cwX = Math.max(1, Math.floor((bw - 2 * MAZE_EDGE - (cols - 1) * MAZE_WW) / cols));
  const cwY = Math.max(1, Math.floor((bh - 2 * MAZE_EDGE - (rows - 1) * MAZE_WW) / rows));
  const stepX = cwX + MAZE_WW, stepY = cwY + MAZE_WW;
  // Leftover from the floor()s above is spread as extra margin, keeping the maze centered.
  const offX = MAZE_EDGE + Math.max(0, Math.floor((bw - 2 * MAZE_EDGE - (cols * stepX - MAZE_WW)) / 2));
  const offY = MAZE_EDGE + Math.max(0, Math.floor((bh - 2 * MAZE_EDGE - (rows * stepY - MAZE_WW)) / 2));
  const cellTop = r => offY + r * stepY, cellLeft = c => offX + c * stepX;

  // Pixel centre of a slot on either parity: odd index -> that cell's middle, even index -> the
  // middle of the gap in front of cell i/2. Used only for dots -- walls and the spawn stay
  // confined to their own slot parity.
  const slotCenter = (i, cellStart, cw, count, span) => {
    if (i & 1) return cellStart((i - 1) >> 1) + (cw >> 1);
    const k = i >> 1;
    const lo = k === 0 ? 0 : cellStart(k - 1) + cw;
    const hi = k === count ? span : cellStart(k);
    return (lo + hi) >> 1;
  };
  const gridPos = (gr, gc) =>
    [slotCenter(gr, cellTop, cwY, rows, bh), slotCenter(gc, cellLeft, cwX, cols, bw)];

  // Carve-out-of-solid: start all wall, open up the cells, then open the passages between them.
  const wallArr = new Uint8Array(bh * bw).fill(enableWalls ? 1 : 0);
  const carve = (y0, x0, y1, x1) => {
    y0 = Math.max(0, y0); x0 = Math.max(0, x0); y1 = Math.min(bh, y1); x1 = Math.min(bw, x1);
    for (let y = y0; y < y1; y++) { const row = y * bw; for (let x = x0; x < x1; x++) wallArr[row + x] = 0; }
  };

  let startCell = null;
  const ghostCells = [];      // {r, c, id}, id from the layout digit -- see placeGhosts()'s level filter
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    const ch = at(2 * r + 1, 2 * c + 1);
    if (ch === 'C' && !startCell) startCell = [r, c];
    if (ch >= '1' && ch <= '9') ghostCells.push({ r, c, id: +ch });
    if (isWall(ch)) continue;                                   // solid cell: leave it filled
    if (enableWalls) carve(cellTop(r), cellLeft(c), cellTop(r) + cwY, cellLeft(c) + cwX);
  }
  // Sorted by id (not scan order): this fixes which array index -- and colour, see render() --
  // a given numbered ghost always gets, regardless of which levels include it.
  ghostCells.sort((a, b) => a.id - b.id);

  // Dots/pellets are read off the whole doubled grid, cell and wall slots alike. Two separate
  // lists rather than one tagged list, since placeDots() stamps both the same way regardless.
  const dotPos = [], powerPos = [];
  for (let gr = 0; gr < gh; gr++) for (let gc = 0; gc < layout[gr].length; gc++) {
    const ch = at(gr, gc);
    if (ch !== '.' && ch !== '+' && ch !== 'O') continue;
    if (gr % 2 === 0 && gc % 2 === 0) continue;   // corner post: always wall, can't hold a dot
    (ch === 'O' ? powerPos : dotPos).push(gridPos(gr, gc));
  }
  if (enableWalls) for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    // Only the gap itself is carved, not the union of the two cells, so a passage next to a
    // solid cell opens the connector without hollowing the cell out.
    if (c + 1 < cols && !isWall(at(2 * r + 1, 2 * c + 2)))
      carve(cellTop(r), cellLeft(c) + cwX, cellTop(r) + cwY, cellLeft(c + 1));
    if (r + 1 < rows && !isWall(at(2 * r + 2, 2 * c + 1)))
      carve(cellTop(r) + cwY, cellLeft(c), cellTop(r + 1), cellLeft(c) + cwX);
  }
  // Same idea on the outer ring: an open wall slot there carves clear to the physical board edge
  // instead of to a neighbour, becoming a real side tunnel since sim.glsl wraps toroidally anyway.
  if (enableWalls) {
    for (let r = 0; r < rows; r++) {
      if (!isWall(at(2 * r + 1, 0))) carve(cellTop(r), 0, cellTop(r) + cwY, cellLeft(0));
      if (!isWall(at(2 * r + 1, 2 * cols))) carve(cellTop(r), cellLeft(cols - 1) + cwX, cellTop(r) + cwY, bw);
    }
    for (let c = 0; c < cols; c++) {
      if (!isWall(at(0, 2 * c + 1))) carve(0, cellLeft(c), cellTop(0), cellLeft(c) + cwX);
      if (!isWall(at(2 * rows, 2 * c + 1))) carve(cellTop(rows - 1) + cwY, cellLeft(c), bh, cellLeft(c) + cwX);
    }
  }

  // Corner posts stay square from fill() above at every grid intersection. Round each one to an
  // inscribed circle, but only on the side(s) that face open floor -- a side facing a solid wall
  // bar stays square so the bar's flat end merges flush with the post instead of the circle
  // notching a chip out of its corner. A corner rounds only when BOTH adjacent sides are open.
  if (enableWalls) for (let i = 0; i <= rows; i++) for (let j = 0; j <= cols; j++) {
    const y0 = i === 0 ? 0 : cellTop(i - 1) + cwY, y1 = i === rows ? bh : cellTop(i);
    const x0 = j === 0 ? 0 : cellLeft(j - 1) + cwX, x1 = j === cols ? bw : cellLeft(j);
    const cy = (y0 + y1) / 2, cx = (x0 + x1) / 2, rad = Math.min(y1 - y0, x1 - x0) / 2;

    // Sampled one pixel outside the post's box, toroidally wrapped (the board wraps).
    const isOpenAt = (y, x) => {
      const yy = ((y % bh) + bh) % bh, xx = ((x % bw) + bw) % bw;
      return wallArr[yy * bw + xx] === 0;
    };
    const midY = Math.floor((y0 + y1 - 1) / 2), midX = Math.floor((x0 + x1 - 1) / 2);
    const openUp = isOpenAt(y0 - 1, midX), openDown = isOpenAt(y1, midX);
    const openLeft = isOpenAt(midY, x0 - 1), openRight = isOpenAt(midY, x1);

    for (let y = y0; y < y1; y++) {
      const row = y * bw, dy = y + 0.5 - cy;
      const openV = dy < 0 ? openUp : openDown;
      for (let x = x0; x < x1; x++) {
        const dx = x + 0.5 - cx;
        const openH = dx < 0 ? openLeft : openRight;
        if (openV && openH && dx * dx + dy * dy > rad * rad) wallArr[row + x] = 0;
      }
    }
  }

  if (!startCell) startCell = [0, 0];
  const cellCenter = ([r, c]) => [cellTop(r) + (cwY >> 1), cellLeft(c) + (cwX >> 1)];

  // Per-cell decision table for the ghosts: layout character, which DIRS headings can leave it,
  // and pixel centre -- built here, the one place that knows both the layout and the geometry.
  const cells = [];
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    const gr = 2 * r + 1, gc = 2 * c + 1;
    const [cy, cx] = cellCenter([r, c]);
    cells.push({ ch: at(gr, gc), open: DIRS.map(([dy, dx]) => !isWall(at(gr + dy, gc + dx))), cy, cx });
  }
  // Which cell a board pixel falls in, or -1 outside the grid. The wall gap past a cell reads as
  // that cell; callers gate on distance to the centre anyway, so the gap never satisfies them.
  const cellIndexAt = (y, x) => {
    const r = Math.floor((y - offY) / stepY), c = Math.floor((x - offX) / stepX);
    return (r >= 0 && r < rows && c >= 0 && c < cols) ? r * cols + c : -1;
  };

  return {
    wall: wallArr, start: cellCenter(startCell), dots: dotPos, power: powerPos,
    // [y, x, id] per ghost -- id is the layout digit, which placeGhosts() filters the current
    // level's spawns against.
    ghosts: ghostCells.map(g => [...cellCenter([g.r, g.c]), g.id]),
    cells, cellIndexAt,
    cellX: cwX, cellY: cwY,
  };
}

// ====================================================================================
//  Soliton placement: exact 90deg rotation, centre of mass placed at the maze start (clipped to
//  the spawn cell's box if wider than it), then anything still on a wall zeroed out.
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

// Resamples a soliton to `scale` of its native size, bilinearly. A Lenia pattern is tied to the
// kernel radius it was found at, so the dots channel's smaller radius (placeDots()) and this
// resample always move together by the same factor.
function resizeSoliton(entry, scale) {
  const h = entry.h, w = entry.w;
  const nh = Math.max(1, Math.round(h * scale)), nw = Math.max(1, Math.round(w * scale));
  const flat = new Float32Array(nh * nw);
  for (let y = 0; y < nh; y++) {
    const fy = (y + 0.5) * h / nh - 0.5;
    const y0 = Math.max(0, Math.min(h - 1, Math.floor(fy))), y1 = Math.min(h - 1, y0 + 1);
    const ty = Math.min(1, Math.max(0, fy - y0));
    for (let x = 0; x < nw; x++) {
      const fx = (x + 0.5) * w / nw - 0.5;
      const x0 = Math.max(0, Math.min(w - 1, Math.floor(fx))), x1 = Math.min(w - 1, x0 + 1);
      const tx = Math.min(1, Math.max(0, fx - x0));
      const v00 = entry.flat[y0 * w + x0], v01 = entry.flat[y0 * w + x1];
      const v10 = entry.flat[y1 * w + x0], v11 = entry.flat[y1 * w + x1];
      const v0 = v00 + (v01 - v00) * tx, v1 = v10 + (v11 - v10) * tx;
      flat[y * nw + x] = v0 + (v1 - v0) * ty;
    }
  }
  return { flat, h: nh, w: nw };
}

// Writes one soliton into `arr` centred on (tr,tc). `mask`, if given, gets a 1 at every cell the
// stamp touches -- how placeDots() tells a power pellet's footprint apart from a dot's, since
// both are the same soliton and rule.
function stampSoliton(arr, entry, tr, tc, rotation, mask) {
  const rot = rotate90(entry.flat, entry.h, entry.w, rotation);
  let sy = 0, sx = 0, sm = 0;
  for (let y = 0; y < rot.h; y++) for (let x = 0; x < rot.w; x++) {
    const v = rot.data[y * rot.w + x]; if (v > 0) { sm += v; sy += v * y; sx += v * x; }
  }
  const cy = sy / sm, cx = sx / sm;

  // Clip window sized to the spawn cell, centred on the soliton's centroid -- trims a soliton
  // wider than the cell instead of letting it overhang into a neighbouring passage (the wall
  // zeroing below only catches an overhang that lands on a wall, not on open corridor next door).
  const halfY = maze.cellY / 2, halfX = maze.cellX / 2;

  const oy = Math.round(tr - cy), ox = Math.round(tc - cx);
  for (let y = 0; y < rot.h; y++) {
    if (Math.abs(y - cy) > halfY) continue;
    const gy = ((oy + y) % H + H) % H, grow = gy * W;
    for (let x = 0; x < rot.w; x++) {
      if (Math.abs(x - cx) > halfX) continue;
      const v = rot.data[y * rot.w + x]; if (v <= 0) continue;
      const gx = ((ox + x) % W + W) % W;
      arr[grow + gx] = v;
      if (mask) mask[grow + gx] = 1;
    }
  }
}

// The dots channel: one soliton per '.'/'+' plus one per 'O' (power pellet -- same soliton and
// rule, just recoloured), all free-running under a single shared rule; CARL never touches it.
// Rebuilt whenever the agent's soliton is placed, so Restart/maze toggle restore the full set.
function placeDots() {
  dotSites = [];
  bonusFruitThresholdIdx = 0;   // a fresh dot set gets both bonus-fruit thresholds back
  if (!DOTS_ENABLED) return;
  // No dots or pellets in the layout: leave the rule unset, so SimGL.step() skips the whole
  // second convolution rather than paying for one over an empty board.
  const entry = (maze.dots.length || maze.power.length)
    ? ruleBank.get(CFG.channel2RuleName) : null;
  if (!entry) return;
  SimGL.setRule2({ mu: entry.mu, sigma: entry.sigma, betas: entry.betas, R: CFG.channel2R, dt: CFG.dt });

  const scaled = resizeSoliton(entry, CFG.channel2R / CFG.R);
  const arr = new Float32Array(N);
  // 1 wherever a power pellet landed, so draw.glsl can colour it apart from a plain dot.
  const power = new Uint8Array(N);
  for (const [r, c] of maze.dots) stampSoliton(arr, scaled, r, c, 0);
  for (const [r, c] of maze.power) stampSoliton(arr, scaled, r, c, 0, power);
  for (let i = 0; i < N; i++) if (wall[i]) { arr[i] = 0; power[i] = 0; }
  SimGL.uploadState2(arr);
  SimGL.setPowerMask(power);

  // One scoring site per stamp (see `score`). Window reaches channel2R either side: a dot's mass
  // stays within ~5px of centre and neighbours sit ~24px apart, so no window picks up another
  // dot's mass. Baseline is summed off the board just uploaded, so a wall-clipped stamp is judged
  // against what it actually started with.
  const rad = CFG.channel2R;
  const windowMass = (r, c) => {
    let m = 0;
    for (let y = r - rad; y <= r + rad; y++) {
      const row = ((y % H) + H) % H * W;
      for (let x = c - rad; x <= c + rad; x++) m += arr[row + ((x % W) + W) % W];
    }
    return m;
  };
  for (const [r, c] of maze.dots) dotSites.push({ r, c, points: 10, spawnMass: windowMass(r, c), eaten: false });
  for (const [r, c] of maze.power) dotSites.push({ r, c, points: 50, spawnMass: windowMass(r, c), eaten: false });
  SimGL.setDotSites(dotSites.map(s => [s.r, s.c]), rad);
}

// The ghost channel: one soliton per numbered spawn, free-running like the dots -- CARL neither
// sees nor steers it. It erases Pac-Man's mass on overlap (glsim.js's ghostEat), which is also
// how a ghost kills: enough erased and finishStep() reads the mass below CFG.massDeathFraction.
// Rebuilt on every spawn (death respawns included), so ghosts always restart in their house.
// Each ghost runs in a private GHOST_WINDOW-square tile rather than a shared field -- two ghosts
// sharing one would be two solitons that annihilate or blow up on contact; tiles let them pass
// through each other instead. See shaders/ghostsim.glsl.
const GHOST_WINDOW = SimGL.GHOST_WIN;   // the tile edge the engine allocates

let ghosts = [];              // one entry per 'M'; [] means there are none to judge or steer

function newGhost(r, c, cellIdx, home = [r, c]) {
  return {
    home,                        // its own numbered cell -- kept for a future level filter, unused by respawn
    y: r, x: c,                 // last known board CoM
    origin: [Math.round(r) - (GHOST_WINDOW >> 1), Math.round(c) - (GHOST_WINDOW >> 1)],
    shift: [0, 0],              // whole cells its tile slides next step, to re-centre it
    cell: cellIdx,
    // Without this the house cell (where it spawns dead-centre) would read as "just passed the
    // centre" and trigger a turn on its very first step.
    turned: true,
    prevD2: Infinity,           // squared distance to that cell's centre one step ago
    anchor: [r, c],             // CoM on entering `cell` -- the baseline its heading is read from
    // Stamped in the bank's own orientation and turned to face the house exit once it has moved
    // far enough to say which way it's going -- aiming it at stamp time would need knowing each
    // rule's canonical heading, which varies.
    aim: cellIdx >= 0 ? maze.cells[cellIdx].open.indexOf(true) : -1,
    spawnMass: 0,               // 0 disarms this ghost's death check
    // Per-ghost: a ghost eaten mid-frightened-window respawns with this false, ending its own
    // window immediately while untouched packmates keep counting down theirs. See updateFrightened().
    frightened: false,
  };
}

// Builds one ghost's tile: the soliton centred in it, masked against the maze cells it covers.
// Returns the tile's mass, which arms its death check immediately rather than waiting on a readback.
function buildGhostTile(entry, g) {
  const win = GHOST_WINDOW, half = win >> 1;
  const tile = new Float32Array(win * win);
  // stampSoliton() works in board coordinates, so stamp into a board-sized scratch and cut the
  // tile out -- one implementation of the centring/clipping, not two.
  const board = new Float32Array(N);
  stampSoliton(board, entry, g.y, g.x, 0);
  let mass = 0;
  for (let ly = 0; ly < win; ly++) {
    const by = ((g.origin[0] + ly) % H + H) % H;
    for (let lx = 0; lx < win; lx++) {
      const bx = ((g.origin[1] + lx) % W + W) % W;
      const v = wall[by * W + bx] ? 0 : board[by * W + bx];
      tile[ly * win + lx] = v;
      mass += v;
    }
  }
  return { tile, mass };
}

function pushGhostTiles() {
  SimGL.setGhostTiles(ghosts.map(g => g.origin), ghosts.map(g => g.shift),
    ghosts.map(g => g.frightened));
}

function placeGhosts() {
  // Cleared first, refilled only once ghosts are actually on the board -- spawnMass=0 arms
  // nothing, so any path leaving the channel empty must leave the death check disarmed too, or
  // "no ghost" reads as "dead ghost" and respawns every step.
  ghosts = [];
  SimGL.setGhostTiles([], []);
  if (!GHOSTS_ENABLED) return;
  // Level N spawns every numbered ghost at or below N (level 3, the start, is ghosts 1-3); a
  // level past the layout's highest number just spawns all of them, since nothing's left to exclude.
  const spawns = maze.ghosts.filter(([, , id]) => id <= level);
  const entry = spawns.length ? ruleBank.get(CFG.channel3RuleName) : null;
  if (!entry) return;      // no ghost spawn at or below this level: leave the whole channel unset
  SimGL.setRule3({ mu: entry.mu, sigma: entry.sigma, betas: entry.betas, R: CFG.R,
                   dt: CFG.dt * CFG.channel3Speed });

  for (const [r, c] of spawns) ghosts.push(newGhost(r, c, maze.cellIndexAt(r, c)));
  pushGhostTiles();        // origins must be current before any tile is written
  ghosts.forEach((g, i) => {
    const { tile, mass } = buildGhostTile(entry, g);
    g.spawnMass = mass;
    SimGL.uploadGhostTile(i, tile);
  });
}

// ------------------------------------------------------------------------------------
//  Ghost steering: 90deg rotations at junctions.
//
//  A free-running soliton already travels straight, so the only decision is what happens at a
//  corner or crossroad. Rather than steer the turn with the policy (a second inference per step),
//  the ghost's own tile is rotated a quarter turn about its CoM (rotate.glsl) -- instant, always
//  succeeds, and (unlike re-stamping the canonical pattern) doesn't read as a blink.
// ------------------------------------------------------------------------------------
// A ghost's heading, read off how it has actually moved since its anchor -- inferring it from the
// stamp rotation would be wrong for any rule not facing right in the bank. -1 until it has
// travelled far enough (6px) for the answer to be trustworthy at 0.3-0.45px/step and up to ~30deg
// off-axis.
const GHOST_HEADING_MIN_PX = 6;
function ghostHeading(g) {
  const [dy, dx] = toroidalDelta(g.y, g.x, g.anchor[0], g.anchor[1]);
  if (dy * dy + dx * dx < GHOST_HEADING_MIN_PX * GHOST_HEADING_MIN_PX) return -1;
  return Math.abs(dx) > Math.abs(dy) ? (dx > 0 ? 0 : 2) : (dy > 0 ? 1 : 3);
}

// Which of `choices` heads most directly at (sign>0) or away from (sign<0) Pac-Man from (gy,gx),
// scored against the toroidal vector to him so a side-tunnel corridor is judged on where it comes
// out. -1 if he isn't on the board. Fleeing just negates the same score, so both share one function.
function ghostChaseDir(choices, gy, gx, rb, sign = 1) {
  if (!rb.valid) return -1;
  const [dy, dx] = toroidalDelta(rb.row, rb.col, gy, gx);
  let best = -1, bestScore = -Infinity;
  for (const k of choices) {
    const score = sign * (DIRS[k][0] * dy + DIRS[k][1] * dx);
    if (score > bestScore) { bestScore = score; best = k; }
  }
  return best;
}

// A ghost dies the same two ways Pac-Man does -- dissolved, or past the explode limit -- judged
// against its own spawn mass on the same CFG thresholds. The consequence differs: no jingle, no
// held board, no episode end, no effect on the other ghosts -- it's just put back in its house.
function ghostDied(g, s) {
  if (!g.spawnMass) return false;         // not armed -- nothing to have died
  return !s.valid || s.mass > CFG.massExplodeLimit || s.mass < CFG.massDeathFraction * g.spawnMass;
}

// Puts one ghost back on the board, leaving the rest of the pack running. Rewriting its tile is
// also the whole cleanup after an explosion: a tile has hard edges, so the mess stays inside it,
// and the board-space texture is rebuilt from the tiles every step rather than accumulated.
//
// Every death sends the ghost back to symbol '1''s cell -- the one actual house in the centre
// room -- rather than its own numbered spawn, since the numbered cells are scattered around the
// maze and not all plausibly a "house". Its own numbered cell still passes through as `home`, for
// a future level filter to know where it belongs.
function respawnGhost(i, eaten) {
  const entry = ruleBank.get(CFG.channel3RuleName);
  if (!entry) return;
  const home = ghosts[i].home;
  const house = maze.ghosts.find(([, , id]) => id === 1);
  const [r, c] = house || home;
  const g = newGhost(r, c, maze.cellIndexAt(r, c), home);
  ghosts[i] = g;
  pushGhostTiles();          // its origin is back at the spawn point before the tile is written
  const { tile, mass } = buildGhostTile(entry, g);
  g.spawnMass = mass;
  SimGL.uploadGhostTile(i, tile);
}

// Flips ghost `i` by a hard 180 if it's heading the "wrong" way for mode `away` (toward Pac-Man
// when it should flee, or away when it should chase) -- left alone if already correct or if it
// hasn't moved far enough for a readable heading (ghostHeading() returns -1): never a left/right
// correction, only "no change or 180".
function turnGhost180(rb, i, away) {
  if (!rb.valid) return;
  const g = ghosts[i], s = rb.ghosts[i];
  if (!s || !s.valid) return;
  const cur = ghostHeading(g);
  if (cur < 0) return;
  const [dy, dx] = toroidalDelta(rb.row, rb.col, g.y, g.x);
  const towardPac = DIRS[cur][0] * dy + DIRS[cur][1] * dx > 0;
  if (towardPac !== away) return;      // already facing the way this mode wants
  SimGL.rotateGhost(i, Math.round(s.localCol), Math.round(s.localRow), 2);
  g.anchor = [g.y, g.x]; g.prevD2 = Infinity;
  g.turned = true;   // this cell's turn is spent -- steerGhost() shouldn't also fire one below
}

// Starts/extends the shared frightened window on a pellet eaten, ends it on the step that first
// notices frightenedExpired -- in both cases turning (turnGhost180()) only the ghosts that
// actually transition, skipping ones already frightened (window just extends) or already ended
// early via respawnGhost(). A pellet eaten now takes effect next step, the same one-step lag
// ghost death detection already has (SimGL.step() used last step's frightened flags).
function updateFrightened(rb) {
  if (rb.pelletEaten > EAT_SOUND_MIN_MASS) {
    playFrightSound();
    clearTimeout(frightenedTimer);
    frightenedTimer = setTimeout(() => { frightenedExpired = true; }, CFG.frightenedDurationMs);
    frightenedExpired = false;   // a fresh pellet supersedes any pending expiry from the old window
    for (let i = 0; i < ghosts.length; i++) {
      if (ghosts[i].frightened) continue;   // already fleeing -- window just extended
      ghosts[i].frightened = true;
      turnGhost180(rb, i, true);
    }
  }
  if (frightenedExpired) {
    frightenedExpired = false;
    for (let i = 0; i < ghosts.length; i++) {
      if (!ghosts[i].frightened) continue;  // ended early already -- nothing to do
      ghosts[i].frightened = false;
      turnGhost180(rb, i, false);
    }
  }
}

// Dots and pellets score once each, the step their own window is found empty (see `score`). The
// combo reset stays on the first bite of a pellet, same moment updateFrightened() opens the window.
function updateScore(rb) {
  if (rb.pelletEaten > EAT_SOUND_MIN_MASS) ghostChainMultiplier = 1;   // even mid-window
  // Length mismatch means a readback from before the last placeDots() -- nothing to judge.
  if (!rb.dotSites || rb.dotSites.length !== dotSites.length) return;
  for (let i = 0; i < dotSites.length; i++) {
    const s = dotSites[i];
    if (s.eaten || rb.dotSites[i] >= CFG.dotGoneFraction * s.spawnMass) continue;
    s.eaten = true;
    score += s.points;
  }
}

function steerGhosts(rb) {
  if (!GHOSTS_ENABLED || !ghosts.length) return;
  for (let i = 0; i < ghosts.length; i++) {
    const s = rb.ghosts[i];
    if (!s) continue;
    if (ghostDied(ghosts[i], s)) {
      const eaten = ghosts[i].frightened;   // dissolved while frightened -- Pac-Man ate it
      const dy = ghosts[i].y, dx = ghosts[i].x;   // captured before respawnGhost() moves it home
      respawnGhost(i, eaten);               // just this one; the rest run on
      if (eaten) {
        const pts = 200 * ghostChainMultiplier;
        score += pts;
        ghostChainMultiplier *= 2;          // next ghost in this same window is worth double
        handleGhostEaten(dy, dx, pts);
      }
      continue;
    }
    if (s.valid) steerGhost(ghosts[i], i, s, rb);
  }
  // Origins and shifts changed above; the engine needs them before the next step's sim.
  pushGhostTiles();
}

function steerGhost(g, idx, s, rb) {
  const win = GHOST_WINDOW, half = win >> 1;
  // Board position is the tile origin plus where the engine reports the soliton inside it; the
  // tile then slides by whole cells to re-centre it (effective on the next step's sim).
  g.y = ((g.origin[0] + s.localRow) % H + H) % H;
  g.x = ((g.origin[1] + s.localCol) % W + W) % W;
  const sy = Math.round(s.localRow) - half, sx = Math.round(s.localCol) - half;
  g.shift = [sy, sx];
  // Wrapped into board range, not left to accumulate: ghostblit.glsl only searches one board-width
  // for a tile's near copy, so an origin drifted further (several trips around the torus) would
  // stop matching any pixel until it wanders back into range.
  g.origin = [(((g.origin[0] + sy) % H) + H) % H, (((g.origin[1] + sx) % W) + W) % W];

  const i = maze.cellIndexAt(g.y, g.x);
  if (i !== g.cell) {
    g.cell = i; g.turned = false; g.prevD2 = Infinity;
    g.anchor = [g.y, g.x];                  // fresh baseline: this cell's run is the heading
  }

  // Pivot for any rotation below, tile-local and integer -- exact only about a whole cell.
  const pivX = Math.round(s.localCol), pivY = Math.round(s.localRow);

  // One-off turn out of the spawn house, once there's enough motion to read a heading. Everything
  // below is ordinary junction logic and doesn't apply yet.
  if (g.aim >= 0) {
    const cur = ghostHeading(g);
    if (cur < 0) return;
    if (cur !== g.aim) SimGL.rotateGhost(idx, pivX, pivY, (g.aim - cur + 4) % 4);
    g.aim = -1;
    g.anchor = [g.y, g.x];
    return;
  }

  if (i < 0 || g.turned) return;
  const cell = maze.cells[i];
  if (!GHOST_TURN_CHARS.includes(cell.ch)) return;

  // Turn at closest approach to the tile centre (distance stops falling), not on first crossing
  // into a band around it: the rotation pivots on the CoM, so turning elsewhere swings the
  // soliton into the corner. A distance threshold can't land on the centre without being small
  // enough for a faster soliton to step clean over.
  const dy = g.y - cell.cy, dx = g.x - cell.cx;
  const d2 = dy * dy + dx * dx;
  if (d2 > g.prevD2) {
    g.turned = true;
    const cur = ghostHeading(g);
    if (cur < 0) return;      // too little travel to read a heading -- leave it running straight

    // Left, straight, right -- never a reversal, so the ghost reads as patrolling. A '^' cell
    // overrides this and sends it north (a centre room's one-way exit). Only actually-open
    // candidates survive; a dead end leaves reversing as the only option.
    const wanted = cell.ch === GHOST_NORTH_CHAR ? [3] : [(cur + 3) % 4, cur, (cur + 1) % 4];
    const open = wanted.filter(k => cell.open[k]);
    let next;
    if (!open.length) next = (cur + 2) % 4;
    else {
      // CFG.chaseBias of the time it takes the opening that closes on (or, frightened, flees)
      // Pac-Man; the rest rolls. Biasing only the turns still compounds hard, since every junction
      // is another chance to correct -- the difficulty dial, not the speed.
      const chase = Math.random() < CFG.chaseBias
        ? ghostChaseDir(open, g.y, g.x, rb, g.frightened ? -1 : 1) : -1;
      next = chase >= 0 ? chase : open[(Math.random() * open.length) | 0];
    }

    if (next === cur) return;                     // straight on: there is nothing to do at all
    SimGL.rotateGhost(idx, pivX, pivY, (next - cur + 4) % 4);
    g.anchor = [g.y, g.x];                        // heading changed: the old run is no baseline
    return;
  }
  g.prevD2 = d2;
}

// resetDots is false only for the automatic respawn after a death: the dots channel keeps running
// as-is (already-eaten dots stay eaten). Every user-triggered (re)spawn restores the full set.
// `resetLives` defaults to `playIntro` -- every ordinary fresh start refills lives, except
// handleWin() clearing the board, which is a fresh spawn but not a new *game*, so it carries
// remaining lives into the next board and passes false explicitly.
function placeSoliton(entry, resetDots = true, playIntro = true, resetLives = playIntro, resetLevel = resetLives) {
  mu = entry.mu; sig = entry.sigma; betas = entry.betas.slice();
  SimGL.setRule({ mu, sigma: sig, betas, R: CFG.R, dt: CFG.dt * CFG.channel1Speed });
  // Must land before placeGhosts() below reads `level` -- respawnAfterWin() has already bumped it
  // by the time this runs.
  if (resetLevel) level = STARTING_LEVEL;

  // Built once on the CPU and uploaded; from here on the board only exists on the GPU.
  const arr = new Float32Array(N);
  const [tr, tc] = maze.start;
  stampSoliton(arr, entry, tr, tc, spawnRotation);
  for (let i = 0; i < N; i++) if (wall[i]) arr[i] = 0;    // applyWallCollision, at spawn
  SimGL.uploadState(arr);
  if (resetDots) placeDots();
  placeGhosts();     // unconditional: a ghost that just ate Pac-Man mustn't still sit on his
                     // spawn tile when he comes back

  clearTimeout(deathTimer);
  clearTimeout(winTimer);
  clearTimeout(frightenedTimer);
  stopEatingSound();
  sndEatGhost.pause(); sndEatGhost.currentTime = 0;   // in case a respawn lands mid-jingle
  sndEatFruit.pause(); sndEatFruit.currentTime = 0;
  steps = 0; actions = 0; courseChanges = 0; solitonDead = false; levelWon = false; ghostEatPause = false;
  ghostEatPopup = null;
  frightenedExpired = false;   // placeGhosts() above already gave every ghost a fresh, unfrightened newGhost()
  ghostChainMultiplier = 1;    // fresh ghosts (placeGhosts() above is unconditional): fresh combo too
  bonusFruit = null; fruitEatPause = false; fruitEatPopup = null;
  if (resetLives) lives = STARTING_LIVES;
  if (resetLives) score = 0;
  sometimesRemaining = sometimesWindow;
  strideCounter = 0;
  lastAction = null; lastQ = null;
  actionTrail.length = 0;
  pendingAction = null;

  // One analysis pass with no step behind it, so the spawn CoM and the policy's first window
  // come from the same place every later step gets them from.
  SimGL.prime();
  const rb = SimGL.readback();
  lastCrop = rb.crop;
  lastCoM = rb.valid ? [rb.row, rb.col, rb.mass] : null;
  initialMass = lastCoM ? lastCoM[2] : 0;
  comHistory = Array.from({ length: CFG.windowSize }, () => lastCoM ? [lastCoM[0], lastCoM[1]] : [tr, tc]);
  render();
  // Skipped for the automatic post-death respawn, which carries on the same run rather than
  // starting a fresh one -- no jingle, no pause.
  if (playIntro) playStartSound();
}
// Held paused through the start jingle, then runs the instant it ends. Skipped for the automatic
// death respawn (see placeSoliton()).
//
// Browsers block audio until a user gesture. Restart/respawn/etc. are already inside a click
// so they play immediately, but the very first call from page load has no gesture yet and gets
// rejected -- so instead of starting silently, show a prompt and retry once the page sees its
// actual first gesture.
function playStartSound() {
  if (!SOUND_ENABLED) { setRunning(true); return; }
  setRunning(false);
  introPlaying = true;
  const attempt = () => { sndStart.currentTime = 0; return sndStart.play(); };
  const finish = () => { introPlaying = false; setRunning(true); };
  attempt()
    .then(() => sndStart.addEventListener('ended', finish, { once: true }))
    .catch(() => {
      showStartPrompt();
      waitForGesture(() => {
        hideStartPrompt();
        attempt()
          .then(() => sndStart.addEventListener('ended', finish, { once: true }))
          .catch(finish);   // blocked even inside a gesture -- give up silently, but still start
      });
    });
}
// Keys the HTML spec excludes from a "user activation" gesture (the UI Events Modifier Keys
// table, plus Escape). Pressing only one of these keeps the game waiting rather than starting
// silently.
const NON_ACTIVATING_KEYS = new Set([
  'Alt', 'AltGraph', 'CapsLock', 'Control', 'Fn', 'FnLock', 'Hyper', 'Meta', 'NumLock', 'OS',
  'ScrollLock', 'Shift', 'Super', 'Symbol', 'SymbolLock', 'Escape',
]);
// Waits for the page's next "real" gesture, then runs `cb` once. Shared by playStartSound()'s
// autoplay-blocked fallback, showGameOverPrompt(), showWinPrompt() and showStartScreen().
function waitForGesture(cb) {
  const start = () => {
    window.removeEventListener('pointerup', onPointer);
    window.removeEventListener('keydown', onKey);
    cb();
  };
  // pointerup, not pointerdown: iOS Safari only counts a *completed* tap as the gesture that
  // unlocks audio. Listening on pointerdown would consume the listener before that, so the retry
  // fails for the same reason as the original call and the game starts silently.
  const onPointer = () => start();
  // Modifier keys (and Escape) don't count as activation -- keep listening past a bare
  // Alt/Ctrl/Shift/CapsLock instead of consuming the gesture on it.
  const onKey = e => { if (!NON_ACTIVATING_KEYS.has(e.key)) start(); };
  window.addEventListener('pointerup', onPointer, { once: true });
  window.addEventListener('keydown', onKey);
}
// Reuses the (otherwise unused) result banner for the "waiting for a gesture" notice, shown both
// on first page load and on game over.
function showStartPrompt() {
  const el = $('result');
  el.className = 'result';
  el.innerHTML = '<b>Start</b>';
  el.hidden = false;
}
function hideStartPrompt() {
  const el = $('result');
  el.hidden = true;
  el.innerHTML = '';
}
// Hold on the "Start" prompt instead of restarting under the player -- same presentation as the
// first page load. respawnCurrentSoliton() runs inside the resulting gesture, so its
// playStartSound() plays immediately rather than needing a second prompt. Shared tail of both
// showGameOverPrompt() (lives ran out) and showWinPrompt() (final level cleared) -- either way the
// next gesture leads back to a fresh game, not a continued one.
function showStartScreen() {
  setRunning(false);
  showStartPrompt();
  waitForGesture(() => { hideStartPrompt(); respawnCurrentSoliton(); });
}
// Lives ran out: hold on a "Game over" banner, then fall through to the "Start" prompt on the next
// gesture.
function showGameOverPrompt() {
  setRunning(false);
  const el = $('result');
  el.className = 'result';
  el.innerHTML = '<b>Game over</b>';
  el.hidden = false;
  waitForGesture(() => { el.hidden = true; el.innerHTML = ''; showStartScreen(); });
}
// Final level cleared: hold on a "You win!" banner instead of quietly advancing past MAX_LEVEL,
// then fall through to the "Start" prompt on the next gesture.
function showWinPrompt() {
  setRunning(false);
  const el = $('result');
  el.className = 'result';
  el.innerHTML = '<b>You win!</b>';
  el.hidden = false;
  waitForGesture(() => { el.hidden = true; el.innerHTML = ''; showStartScreen(); });
}

function newMaze() {
  spawnRotation = (Math.random() * 4) | 0;
  maze = layoutMaze(MAZE_LAYOUT, H, W, mazeEnabled);
  wall = maze.wall;
  SimGL.setWall(wall);
  placeSoliton(bank[currentIndex]);
}
function initBoard() {
  SimGL.setBoard(W, H);
  boardReady = true;
  newMaze();
}
function respawnCurrentSoliton() { placeSoliton(bank[currentIndex]); }
// Death-triggered respawn -- same spawn, but leaves the dots channel alone (see placeSoliton()).
function respawnAfterDeath() { placeSoliton(bank[currentIndex], false, false); }
// Win-triggered respawn -- full fresh spawn, but lives carry over rather than refilling (the life
// gained for the win is already applied in handleWin()), the level advances first so placeGhosts()
// spawns the next level's count, and `playIntro` is false since the intermission jingle just
// finished.
function respawnAfterWin() { level++; placeSoliton(bank[currentIndex], true, false, false); }

// ====================================================================================
//  Agent step
// ====================================================================================
// CPU-side copy of crop.glsl's origin math, used to map the policy's per-cell output back to
// board coordinates and to bound what the user may click in human mode.
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
  // The crop arrives channel-packed; the model wants frame-major [1,K,S,S]. This de-interleave is
  // the whole cost of the GPU input path.
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

// The user's half of a turn: hand over the one queued action, if any, and clear the queue.
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
  else if (actorMode === 'sometimes' && sometimesRemaining <= 0) {
    lastMs = 0; lastQ = null; lastAction = null;      // CARL idle this step -- no inference, no action
  }
  else if (!session) return;
  else if (strideCounter > 0) {
    strideCounter--;
    lastMs = 0; lastQ = null; lastAction = null;      // stride-skipped step -- no inference, no action
    if (actorMode === 'sometimes') sometimesRemaining--;
  }
  else {
    action = await agentAct();
    strideCounter = inferenceStride - 1;
    if (actorMode === 'sometimes') sometimesRemaining--;
  }

  // Action, step, reduction and the next crop, queued back to back; the readback below is the
  // only point the CPU waits on the GPU. Per-ghost frightened state rides along via
  // pushGhostTiles() at the end of steerGhosts(), not as an argument here.
  SimGL.step(action);
  const rb = SimGL.readback();
  lastCrop = rb.crop;
  steps++;
  updateFrightened(rb);
  updateScore(rb);
  steerGhosts(rb);
  finishStep(rb);
}

// Bookkeeping shared by both actors: take the soliton's freshly computed CoM and decide whether
// the episode has ended.
function finishStep(rb) {
  if (rb.eaten > EAT_SOUND_MIN_MASS) noteEating();
  if (rb.hasDots && !solitonDead && !levelWon && rb.dotsMass < CFG.dotsWinMass) {
    handleWin();
    return;
  }
  if (!rb.valid) {
    lastCoM = null;
    if (!solitonDead) handleDeath();
    return;
  }
  lastCoM = [rb.row, rb.col, rb.mass];
  comHistory.push([rb.row, rb.col]); if (comHistory.length > CFG.windowSize) comHistory.shift();
  if (!solitonDead) {
    if (rb.mass > CFG.massExplodeLimit) handleDeath(true);
    else if (rb.mass < CFG.massDeathFraction * initialMass) handleDeath();
  }
  if (!solitonDead && !levelWon) updateBonusFruit();
}

// Bonus fruit: appears on Pac-Man's own spawn tile the moment the level's dots-eaten fraction
// crosses CFG.bonusFruitThresholds[bonusFruitThresholdIdx] -- twice a level at the defaults, since
// there are two thresholds -- and disappears unresolved after CFG.bonusFruitTimeout steps. Emoji
// and points come from BONUS_FRUIT_EMOJIS/BONUS_FRUIT_POINTS at the current `level` (1-9), fixed
// at spawn so a mid-flight level change (there isn't one, but just in case) can't retarget an
// already-showing fruit. Picked up on contact -- unlike an ordinary dot it isn't part of the
// channel-2 field at all, just a CPU-tracked position judged against lastCoM, since a bonus item
// is meant to disappear the instant Pac-Man touches it rather than dissolve over several steps.
function updateBonusFruit() {
  if (bonusFruit) {
    if (steps - bonusFruit.spawnStep >= CFG.bonusFruitTimeout) bonusFruit = null;
  } else if (dotSites.length && bonusFruitThresholdIdx < CFG.bonusFruitThresholds.length) {
    const eatenFrac = dotSites.filter(s => s.eaten).length / dotSites.length;
    if (eatenFrac >= CFG.bonusFruitThresholds[bonusFruitThresholdIdx]) {
      bonusFruitThresholdIdx++;
      const [r, c] = maze.start;
      const idx = Math.min(BONUS_FRUIT_EMOJIS.length, Math.max(1, level)) - 1;
      bonusFruit = { r, c, emoji: BONUS_FRUIT_EMOJIS[idx], points: BONUS_FRUIT_POINTS[idx], spawnStep: steps };
    }
  }
  if (!bonusFruit || !lastCoM) return;
  const [dy, dx] = toroidalDelta(lastCoM[0], lastCoM[1], bonusFruit.r, bonusFruit.c);
  if (dy * dy + dx * dx > BONUS_FRUIT_RADIUS * BONUS_FRUIT_RADIUS) return;
  const { r, c, points } = bonusFruit;
  score += points;
  bonusFruit = null;
  handleBonusFruitEaten(r, c, points);
}

// Same beat as handleGhostEaten(): holds the board for the jingle's own length, no extra pause.
// (r, c) is where the fruit sat, shown as a floating score popup for as long as the board holds.
function handleBonusFruitEaten(r, c, points) {
  fruitEatPause = true;
  fruitEatPopup = { r, c, text: `${points}` };
  const finish = () => { fruitEatPause = false; fruitEatPopup = null; };
  if (!SOUND_ENABLED) { finish(); return; }
  sndEatFruit.currentTime = 0;
  sndEatFruit.play()
    .then(() => sndEatFruit.addEventListener('ended', finish, { once: true }))
    .catch(finish);
}
// Clears rounding error in the GPU reduction, not a heuristic threshold.
const EAT_SOUND_MIN_MASS = 1e-4;
// eat_dot_0/1 "waka waka" pair, alternating while dots keep getting eaten. Every eaten dot pushes
// eatDeadline out by CFG.eatSoundGraceMs; once passed, the next scheduled play sees that and stops
// rather than queuing another sample, so a sample already playing always finishes uncut and never
// overlaps a fresh start.
function noteEating() {
  if (!SOUND_ENABLED) return;
  eatDeadline = performance.now() + CFG.eatSoundGraceMs;
  if (!eatActive) { eatActive = true; playNextEatSound(); }
}
function playNextEatSound() {
  if (!eatActive) return;               // stopEatingSound() cut the chain -- nothing to continue
  if (performance.now() >= eatDeadline) { eatActive = false; return; }
  const snd = sndEat[eatToggle];
  eatToggle = 1 - eatToggle;
  snd.currentTime = 0;
  snd.play().catch(() => { eatActive = false; });   // playback blocked -- don't spin retrying forever
}
// Stops the loop outright, mid-note if need be -- used on a new life so a leftover chomp can't
// bleed into it.
function stopEatingSound() {
  eatActive = false;
  for (const snd of sndEat) { snd.pause(); snd.currentTime = 0; }
}
// Death is a beat, not a stop: while lives remain the sim just holds (loop() stops stepping while
// solitonDead) through the jingle plus a further second, then respawns on its own. Game over once
// lives run out instead -- see showGameOverPrompt().
const DEATH_PAUSE_MS = 1000;
// `exploded` (mass-runaway, see finishStep()) costs no life, since it isn't something the player
// could have steered around the way a ghost is -- same beat and jingle otherwise.
function handleDeath(exploded = false) {
  solitonDead = true;
  if (!exploded && !GOD_MODE) lives = Math.max(0, lives - 1);
  // Lives run out: hold on a "Game over" banner rather than restarting immediately -- the player
  // chooses when the next game begins.
  const next = lives > 0 ? respawnAfterDeath : showGameOverPrompt;
  if (!SOUND_ENABLED) { deathTimer = setTimeout(next, DEATH_PAUSE_MS); return; }
  sndDeath.currentTime = 0;
  // Timed off the jingle actually ending; if playback is blocked, fall back to the pause alone.
  const afterSound = () => { deathTimer = setTimeout(next, DEATH_PAUSE_MS); };
  sndDeath.play()
    .then(() => sndDeath.addEventListener('ended', afterSound, { once: true }))
    .catch(afterSound);
}

// Clearing the board is also a beat, not a stop: the sim holds through the intermission jingle
// plus a further second, then a fresh spawn -- same path as Restart, except lives carry over
// since this is a level clear, not a new game.
const WIN_PAUSE_MS = 1000;
function handleWin() {
  levelWon = true;
  // Awarded ahead of the intermission jingle so the 💛 count updates before the reward beat plays.
  lives = Math.min(MAX_LIVES, lives + 1);
  stopEatingSound();      // the last dot's chomp shouldn't bleed into the intermission jingle
  const next = level >= MAX_LEVEL ? showWinPrompt : respawnAfterWin;
  if (!SOUND_ENABLED) { winTimer = setTimeout(next, WIN_PAUSE_MS); return; }
  sndIntermission.currentTime = 0;
  const afterSound = () => { winTimer = setTimeout(next, WIN_PAUSE_MS); };
  sndIntermission.play()
    .then(() => sndIntermission.addEventListener('ended', afterSound, { once: true }))
    .catch(afterSound);
}

// Fire-and-forget, unlike the death/win/eat-ghost jingles -- the game keeps running underneath
// it, so it's just restarted from the top on every pellet eaten rather than queued or awaited.
function playFrightSound() {
  if (!SOUND_ENABLED) return;
  sndFright.currentTime = 0;
  sndFright.play().catch(() => {});
}

// The board holds for the sound's own length (same beat as handleDeath()/handleWin(), no extra
// pause after) -- the eaten ghost already went home via respawnGhost(), called just before this.
// (r, c) is where it died (captured by the caller before respawnGhost() moved it), shown as a
// floating score popup for exactly as long as the board holds.
function handleGhostEaten(r, c, points) {
  ghostEatPause = true;
  ghostEatPopup = { r, c, text: `${points}` };
  const finish = () => { ghostEatPause = false; ghostEatPopup = null; };
  if (!SOUND_ENABLED) { finish(); return; }
  sndEatGhost.currentTime = 0;
  sndEatGhost.play()
    .then(() => sndEatGhost.addEventListener('ended', finish, { once: true }))
    .catch(finish);
}

// ====================================================================================
//  Rendering -- the board is a shader pass (draw.glsl); everything below is the overlay canvas
//  on top, unchanged from the CPU demo.
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
// Green/red discs stamped on the board each step, fading out over CFG.actionFadeSeconds of sim
// time instead of blinking for a single frame, so the intervention pattern reads at 20+ steps/sec.
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
// One 👀 per ghost, centred on its CoM -- g.y/g.x are continuous (straight off the GPU tile
// reduction), so this uses the same no-offset mapping lastCoM's arrows use, not actionTrail's
// "+0.5 to recentre a rounded index" one.
function drawGhostEyes(sc) {
  if (!ghosts.length) return;
  octx.font = `${Math.max(7, (maze.cellY || 40) * sc * 0.5)}px sans-serif`;
  octx.textAlign = 'center';
  octx.textBaseline = 'bottom';
  // fillStyle isn't reset between draws -- drawActionTrail() leaves it at a low-alpha rgba(), so
  // without setting it here the eyes would inherit that leftover alpha instead of drawing opaque.
  for (const g of ghosts) {
    octx.fillStyle = g.frightened ? '#fff7' : '#ffff';
    octx.fillText('👀', g.x * sc, g.y * sc);
  }
}
// The bonus fruit: a plain dot (same look as an ordinary one) with its emoji drawn over it, both
// gone the instant updateBonusFruit() judges it eaten.
function drawBonusFruit(sc) {
  if (!bonusFruit) return;
  const px = bonusFruit.c * sc, py = bonusFruit.r * sc;
  octx.beginPath(); octx.arc(px, py, Math.max(2, sc * 3), 0, 7);
  octx.fillStyle = `rgb(${COLORS.soliton2.join(',')})`; octx.fill();
  octx.font = `${Math.max(10, (maze.cellY || 40) * sc * 0.8)}px sans-serif`;
  octx.textAlign = 'center'; octx.textBaseline = 'middle';
  octx.fillText(bonusFruit.emoji, px, py);
}
// A floating score number under wherever a ghost or the bonus fruit was just eaten, live for
// exactly as long as its own pause flag holds the board (cleared by handleGhostEaten()/handleBonusFruitEaten()
// alongside ghostEatPause/fruitEatPause, not by a timer of its own here).
function drawScorePopup(sc, p) {
  if (!p) return;
  const px = p.c * sc, py = (p.r + (maze.cellY || 40) * 0.4) * sc;   // slightly below the spot
  octx.font = `bold ${Math.max(11, (maze.cellY || 40) * sc * 0.42)}px sans-serif`;
  octx.textAlign = 'center'; octx.textBaseline = 'top';
  octx.lineWidth = Math.max(1, sc * 0.6);
  octx.strokeStyle = 'rgba(0,0,0,.7)'; octx.strokeText(p.text, px, py);
  octx.fillStyle = '#fff'; octx.fillText(p.text, px, py);
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
let boardReady = false;   // render() can be reached from UI handlers before the first initBoard()
function render() {
  if (!boardReady) return;
  // Per ghost, not all-or-nothing: an already-respawned packmate keeps its own colour while the
  // rest are still fleeing blue. Built fresh each frame since which index is which colour depends
  // on state.
  const ghostColors = ghosts.map((g, i) =>
    g.frightened ? COLORS.frightened : COLORS.ghosts[i % COLORS.ghosts.length]);
  SimGL.draw({ ...COLORS, ghosts: ghostColors });

  const r = cv.getBoundingClientRect();
  if (ov.width !== Math.round(r.width)) { ov.width = Math.round(r.width) || W; ov.height = Math.round(r.height) || H; }
  octx.clearRect(0, 0, ov.width, ov.height);
  const sc = ov.width / W;

  // Not gated on showOverlay unlike the rest below: it's part of what a ghost *is* to the player,
  // not a CARL debugging aid.
  drawGhostEyes(sc);
  drawBonusFruit(sc);
  drawScorePopup(sc, ghostEatPopup);
  drawScorePopup(sc, fruitEatPopup);

  // CARL's working-out drawn over the board -- direction arrows and intervention discs. Hidden
  // together with the control deck, leaving just the game.
  if (showOverlay) {
    if (lastCoM) {
      // One tile long: enough to read the heading, short enough not to cover the maze.
      const cx = lastCoM[1] * sc, cy = lastCoM[0] * sc, alen = (maze.cellX || 40) * sc;
      const cs = getComputedStyle(document.documentElement);
      // Hide the target arrow while CARL is idle -- nothing is acting on it, so it would claim an
      // agent is working when it isn't.
      if (!humanActs) drawArrow(octx, cx, cy, dir[1], dir[0], alen, cs.getPropertyValue('--target').trim());
      const old = comHistory[0];
      const [vdy, vdx] = toroidalDelta(lastCoM[0], lastCoM[1], old[0], old[1]);
      drawArrow(octx, cx, cy, vdx, vdy, alen, cs.getPropertyValue('--current').trim());
    }
    drawActionTrail(sc);
    drawPendingAction(sc);
  }

  $('r-mass').textContent = lastCoM ? lastCoM[2].toFixed(0) : '0';
  $('r-step').textContent = steps;
  $('r-acts').textContent = actions;
  $('r-course').textContent = courseChanges;
  $('r-ms').textContent = lastMs ? lastMs.toFixed(1) + 'ms' : '–';
  $('r-rate').textContent = measSps ? measSps.toFixed(0) + '/s' : '–';
  $('lives').textContent = '💛'.repeat(lives);
  $('score').textContent = score;
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
    // Steps the budget can't afford are dropped, not carried over (that would only lengthen the
    // next batch). Checked after a step, not before, so one step blowing the whole budget still
    // advances by one instead of stalling forever.
    let done = 0;
    try {
      for (let k = 0; k < want; k++) {
        if (solitonDead || levelWon || ghostEatPause || fruitEatPause) break;   // holding for a jingle
        await agentStep();
        done++;
        if (performance.now() - now >= FRAME_BUDGET_MS) break;
      }
    } finally { busy = false; }
    rateSteps += done; rateTime += elapsed;
    if (rateTime >= 0.5) { measSps = rateSteps / rateTime; rateSteps = 0; rateTime = 0; }
    if (done > 0) render();
  }
  requestAnimationFrame(loop);
}
function setRunning(v) {
  if (v && introPlaying) return;   // no unpausing until the start jingle has finished
  running = v && !solitonDead;
  $('play').textContent = running ? 'Pause' : 'Play';
  $('play').classList.toggle('on', running);
  if (running) { lastT = 0; stepAcc = 0; rateSteps = 0; rateTime = 0; loop(); }
  else render();
}

// ====================================================================================
//  Input: click-to-steer, arrow keys, sliders/buttons/selects
// ====================================================================================
let actionCost = 2.5;
// Every user-initiated steer goes through here. Re-issuing the direction already being chased
// (holding an arrow key, clicking straight ahead) isn't counted as a change of course.
function steerTo(dy, dx) {
  const n = Math.hypot(dy, dx); if (n < 1e-6) return;
  const ny = dy / n, nx = dx / n;
  if (ny * dir[0] + nx * dir[1] < 0.9999) courseChanges++;
  dir = [ny, nx];
  sometimesRemaining = sometimesWindow;    // a fresh steer wakes CARL up again in 'sometimes' mode
  strideCounter = 0;                       // ...and gets an immediate inference, not a stale skip
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
  if (v) { steerTo(v[0], v[1]); e.preventDefault(); return; }
  if (e.key === ' ') { setRunning(!running); e.preventDefault(); }
});
$('play').addEventListener('click', () => setRunning(!running));
// Step always leaves the run paused, one step further on. If the loop is mid-batch when the
// click lands (busy), that in-flight step is the advance -- stepping again here would double it.
async function stepOnce() {
  setRunning(false);
  if (busy || solitonDead) return;
  busy = true;
  try { await agentStep(); } finally { busy = false; }
  render();
}
$('step1').addEventListener('click', stepOnce);

// ------------------------------------------------------------------------------------
//  Human actor: the same intervention CARL makes (same +-MA over CFG.actionRadius, same disc,
//  same tally), at most one per step. A click queues the action; the next step stamps it in at
//  exactly the point CARL's own action would land.
// ------------------------------------------------------------------------------------
// CARL's reach: the policy returns one q value per cell of its CFG.netSize window, so a spot
// outside it isn't a move CARL could make. Same bound, same origin math, for the user.
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
  if (solitonDead) return;
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
const ACTOR_MODES = ['sometimes', 'always', 'human'];
const ACTOR_LABELS = { sometimes: 'CARL acts sometimes', always: 'CARL acts always', human: 'You act' };
function setActorMode(mode) {
  actorMode = mode;
  humanActs = mode === 'human';
  $('actor').textContent = ACTOR_LABELS[mode];
  $('actor').classList.toggle('on', humanActs);
  $('human-hint').hidden = !humanActs;
  $('legend-target').hidden = humanActs;         // no blue arrow on the board, no key for it
  cv.classList.toggle('human', humanActs);
  if (mode === 'sometimes') sometimesRemaining = sometimesWindow;   // fresh window on entering the mode
  strideCounter = 0;
  lastMs = 0; lastQ = null; lastAction = null; pendingAction = null;
  render();
}
$('actor').addEventListener('click', () => {
  setActorMode(ACTOR_MODES[(ACTOR_MODES.indexOf(actorMode) + 1) % ACTOR_MODES.length]);
});
$('respawn').addEventListener('click', respawnCurrentSoliton);
$('shuffle').addEventListener('click', () => {
  currentIndex = (Math.random() * bank.length) | 0;
  updateSolitonSelection();
  respawnCurrentSoliton();
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

// Chrome toggle: the deck and CARL's overlay go together. The button itself stays put, since it
// is the only way back to the controls.
function setChrome(on) {
  showOverlay = on;
  $('deck').hidden = !on;
  const b = $('ui-toggle');
  b.textContent = on ? '👁' : '🙈';
  b.title = on ? 'Hide the controls and the agent overlay' : 'Show the controls and the agent overlay';
  b.setAttribute('aria-pressed', String(!on));
  render();
}
$('ui-toggle').addEventListener('click', () => setChrome(!showOverlay));

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
// Stops the run without going through setRunning() (which would try to render) -- fail() must
// stay callable before the board exists, e.g. when GPU init itself is what failed.
function fail(msg) {
  $('tip').innerHTML = msg;
  running = false;
  $('play').textContent = 'Play';
  $('play').classList.remove('on');
}
// No inline fallback soliton, unlike the CPU demo: this version fetches shaders too, so it can't
// run from file:// under any circumstances.
async function loadSolitonBank() {
  let raw = null;
  try {
    const res = await fetch(CFG.solitonsUrl);
    if (res.ok) raw = await res.json();
  } catch (e) { /* reported below */ }
  if (!raw) { fail('Could not load the soliton bank — serve this folder over http(s), not file://.'); return false; }
  ruleBank = new Map(raw.filter(e => Math.max(e.h, e.w) <= 80).map(e => [e.name, {
    name: e.name, mu: e.mu, sigma: e.sigma, betas: e.betas, h: e.h, w: e.w,
    flat: Float32Array.from(e.state.flat()),
  }]));
  // Picker order is CFG.allowedRuleNames', not the JSON's, so the leading entries land on top row.
  bank = CFG.allowedRuleNames.map(name => ruleBank.get(name)).filter(Boolean);
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
async function loadModel() {
  ort.env.wasm.wasmPaths = 'https://cdn.jsdelivr.net/npm/onnxruntime-web@1.20.1/dist/';
  // Multi-threaded WASM needs SharedArrayBuffer, which needs cross-origin isolation -- still false
  // on a first visit before coi-serviceworker takes over, so fall back to one thread rather than
  // let onnxruntime-web throw. ?threads=N overrides (needs a fresh page load to take effect).
  const threadsOverride = parseInt(new URLSearchParams(location.search).get('threads'), 10);
  ort.env.wasm.numThreads = threadsOverride > 0 ? threadsOverride
    : (window.crossOriginIsolated ? Math.min(navigator.hardwareConcurrency || 6, 6) : 1);

  // Absorb WASM's one-time JIT cost here rather than on the user's first real sim step.
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
    session = await ort.InferenceSession.create(CFG.modelUrl, { executionProviders: EXECUTION_PROVIDERS });
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
  initBoard();     // spawns the soliton, which plays the start jingle and then starts the game
  loadModel();
})();
