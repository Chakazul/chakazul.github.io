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
//  URL params and feature switches -- e.g. ?sound=0&dots=0&net=96
// ====================================================================================
// "0"/"false" turns a switch off; anything else, or the param being absent, leaves the default.
function boolParam(name, def) {
  const v = new URLSearchParams(location.search).get(name);
  return v === null ? def : v !== '0' && v.toLowerCase() !== 'false';
}
// Positive integers only -- absent, unparseable or <= 0 all leave the default.
function intParam(name, def) {
  const v = parseInt(new URLSearchParams(location.search).get(name), 10);
  return v > 0 ? v : def;
}
// A whole percentage, 0-100, returned as a fraction. Unlike intParam this accepts 0: for the
// things it reads, "never" is a setting someone might actually want, not a missing value.
function pctParam(name, def) {
  const v = parseInt(new URLSearchParams(location.search).get(name), 10);
  return v >= 0 && v <= 100 ? v / 100 : def;
}
// Sound effects + the death beat. Off: no start jingle (the game begins the instant the soliton
// spawns instead of waiting on one), no eat_dot chomp, and a death just holds+respawns in silence
// instead of playing the death jingle first.
const SOUND_ENABLED = boolParam('sound', true);
// The pellet (channel-2) soliton field -- its rule, its stamps, its rendering. Off: placeDots()
// never calls SimGL.setRule2()/uploadState2(), so ch2Kernel stays null and sim.glsl's per-step
// channel-2 pass and its contribution to the drawn board are both skipped -- only the Pac-Man
// channel moves.
const DOTS_ENABLED = boolParam('dots', true);
// The ghost (channel-3) soliton field, same switch shape as the dots. Off: placeGhosts() never
// calls SimGL.setRule3()/uploadState3(), so ch3Kernel stays null and both halves of the channel
// are skipped -- its own per-step sim pass, and the eat check that erases Pac-Man where a ghost
// overlaps him. So this is also the switch for "nothing can kill Pac-Man but himself".
const GHOSTS_ENABLED = boolParam('ghost', true);
// Execution provider for the policy net: 'wasm' (default) or '?ep=webgpu' to try GPU compute
// instead. WebGPU can't remove the CPU roundtrip in agentStep() -- the sim runs in a separate
// WebGL2 context with no memory sharing with WebGPU, so the crop still crosses through CPU
// either way -- but it can still speed up the net's own conv work (a real 4x96x96 CNN, not a
// toy MLP) on a device with a capable, well-supported GPU. Mobile WebGPU support is newer and
// patchier than WASM SIMD+threads (older Android drivers, iOS Safari, in-app webviews), so this
// is a knob for A/B testing on real devices rather than a default change. 'webgpu' is listed
// with a 'wasm' fallback so any op the WebGPU EP doesn't support still lands on wasm.
const EXECUTION_PROVIDERS = new URLSearchParams(location.search).get('ep') === 'webgpu'
  ? ['webgpu', 'wasm'] : ['wasm'];

// ====================================================================================
//  CONFIG -- locked to the canonical direction run (models/meta_direction.json)
// ====================================================================================
const CFG = {
  K: 4,                       // frame stack
  R: 18,                      // kernel radius (fixed -- all 48 training solitons share it)
  netSize: intParam('net', 96), // the model was trained on 96x96 toroidal grids with no mazes --
                                // every inference call crops a 96x96 toroidal window centered on
                                // the soliton's CoM out of the (larger) maze board, so the model
                                // always sees input matching its training distribution regardless
                                // of the selected board size, and runs faster to boot. ?net=N
                                // overrides the window edge for experiments: the policy is fully
                                // convolutional so a different N still runs, but it is off the
                                // distribution the model was trained on, and N also sets how much
                                // the one readback per step has to carry (N*N*4 floats).
  dt: 0.1,
  // Per-channel multipliers on `dt` above -- Pac-Man's own speed and the ghosts' own speed,
  // independent of each other and of the dots channel (which has none: it never moves). Each
  // layers under glsim.js's own frightened-window multiplier rather than replacing it: Pac-Man
  // runs at dt*channel1Speed normally and dt*channel1Speed*1.25 (PACMAN_FRIGHTENED_SPEED) while
  // any ghost is frightened; a ghost runs at dt*channel3Speed normally and dt*channel3Speed*0.75
  // (GHOST_FRIGHTENED_SPEED) while it itself is. 1 leaves both exactly as they were.
  channel1Speed: 1.5,
  channel3Speed: 1.5,
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
  // ...and as blown up above this absolute mass. Lowered from the CPU demo's 1500, which windowing
  // exposed as unreachable: measured on rule74, an undisturbed soliton peaks at 437 (146% of its
  // ~298 spawn mass) and one with CARL-sized mass injected every single step plateaus at 679, so
  // 1500 never fired and a blown-up Pac-Man simply sat there. 1000 clears that forced ceiling by
  // ~1.5x and the healthy peak by ~2.3x, while staying far under what a real runaway reaches --
  // it is 6% of the 128-square window's capacity, so genuine unbounded growth crosses it early.
  // Note the window is not what caps the mass: 96 and 128 windows gave the identical 679 peak, so
  // this is the growth function's own ceiling, not a clipping artefact.
  massExplodeLimit: 1000,
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
    'rule46_mu0.2500_s0.0270_R18', 'rule63_mu0.3200_s0.0660_R18', 'rule0_mu0.3800_s0.0700_R18',
  ],
  defaultRuleName: 'rule74_mu0.2300_s0.0350_R18',  // or rule73_mu0.2250_s0.0250_R18
  // The dots channel's rule. One soliton of it is dropped on every '.' cell of the layout; they
  // then run free under their own growth parameters, with no policy steering them.
  channel2RuleName: 'rule0_mu0.3800_s0.0700_R18',
  // Half CFG.R, which shrinks the dots to half size -- see resizeSoliton() in placeDots(): the
  // pattern is always resampled by exactly channel2R/R, so this is the only number to change.
  channel2R: 9,
  // The ghost channel's rule. One soliton on every 'M' of the layout, respawned alongside Pac-Man.
  // Kept at the native CFG.R (unlike the dots) because a ghost is meant to be Pac-Man-sized and to
  // travel the corridors the way he does, and a Lenia pattern only moves like itself at the radius
  // it was found at. Free-running: CARL never sees this channel and never steers it.
  //
  // Chosen for speed -- a ghost that cannot keep up is not a threat. Measured in a corridor of
  // this maze, rule74 runs 0.479 px/step against rule73's 0.194, a wider gap than in open space
  // (0.452 vs 0.270) because rule73 wastes much more of its motion on the ~32deg lean it travels
  // at, where rule74 leans only ~11deg. It being the same rule as defaultRuleName is deliberate
  // rather than incidental: it puts the ghost at exactly Pac-Man's own pace on the default pick,
  // so it closes only when it out-navigates him, not because it simply moves faster.
  channel3RuleName: 'rule74_mu0.2300_s0.0350_R18',
  // How often a ghost's turn is a deliberate move toward Pac-Man rather than a roll of the dice --
  // the difficulty dial. At 0 it wanders and only meets him by accident; at 1 it closes on him at
  // every junction that offers the option. ?chase=N (a whole percentage, 0-100) to retune.
  chaseBias: pctParam('chase', 0.7),
  easyRuleCount: 5,            // leading entries flagged as easy-to-steer -- one picker row
  actionFadeSeconds: 0.85,     // how long an intervention marker takes to fade out, in sim time
  // How long a gap with nothing eaten is allowed before the eat_dot_0/1 "waka waka" loop stops.
  // Measured in real (wall-clock) time, not sim time, since it times a sound rather than the sim.
  eatSoundGraceMs: 1000,
  // The board counts as cleared once the dots channel's live total mass (SimGL.readback().dotsMass,
  // a real GPU reduction over channel 2 -- see runDotsAnalysis()) falls below this absolute value --
  // not exactly 0, since a dot mid-erasure can leave a sliver behind that channel 2's own growth
  // rule would otherwise sustain forever. An absolute mass, not a fraction of the mass placeDots()
  // spawned with: the dots are themselves free-running Lenia solitons, so their total drifts with
  // the rule's own growth/decay rather than only ever going down as they're eaten, and a fraction
  // of a stale spawn-time snapshot doesn't track that.
  dotsWinMass: 0.1,
  // How long a power pellet's frightened window stays open, in wall-clock ms -- same convention as
  // eatSoundGraceMs: about player-perceived time, not sim steps, so it doesn't scale with sim
  // speed. Eating another pellet while already frightened resets the clock rather than stacking.
  frightenedDurationMs: 15000,
};
const MA = CFG.actionValue[0]; // action magnitude (0.3)

// Board palette, 0-255. All of these go to draw.glsl as uniforms, so changing them here is the
// only edit needed: the board ramps from `background` at mass 0 to `soliton` at mass 1, channels
// 2 and 3 are laid over that in `soliton2`/`soliton3`, and wall cells are painted flat in `wall`.
const COLORS = {
  soliton:    [255, 221, 51],   // Pac-Man yellow
  soliton2:   [255, 255, 255],  // ordinary dots, in the free-running channel-2 field
  pellet:     [255, 255, 0],    // power pellets, same field -- see placeDots()'s power mask
  soliton3:   [255, 40, 40],    // fallback for ghost mass with no owner recorded -- arcade red,
                                // the one hue that reads as danger against the yellow and blue
  // One per ghost, in spawn order, cycled if there are more ghosts than colours. Which one applies
  // to a pixel is recorded in the composited texture by ghostblit.glsl, so it stays exact even
  // where two ghosts' tiles overlap in board space.
  ghosts: [
    [255,  40,  40],            // Blinky red
    [ 90, 160, 255],            // Inky blue
    [ 80, 230, 120],            // Clyde green
    [255, 170,  60],            // amber
    [124,  58, 237],            // purple -- same hue as the CSS --easy accent
  ],
  wall:       [33, 33, 180],    // dark arcade blue
  background: [0, 0, 0],
  frightened: [0, 0, 139],      // every ghost's colour during a power pellet's frightened window
};

// ====================================================================================
//  Simulation state
//
//  Note what is *not* here any more: the board array, the frame stack, the kernel taps and the
//  FFT scratch buffers. The board lives in GPU textures and the frame stack is a ring of them;
//  the only board-sized array left on this side is the wall mask, which is generated once per
//  maze and uploaded.
// ====================================================================================
// Wide board: 5 rows x 11 columns of maze cells. 550x250 is exactly 11:5, so the pixel aspect
// matches the tile grid and both axes land on ~40px cells (see the MAZE_LAYOUT comment below).
// There is no board-size control -- the cell count is fixed by the layout, and any other size
// just starves the corridors.
const BOARD_H = 250, BOARD_W = 550;
let H = BOARD_H, W = BOARD_W, N = H * W;
let mu = 0.24, sig = 0.024, betas = [1.0, 0.5];
let wall = new Uint8Array(0);
let maze = { start: [0, 0] };
let dir = [0, 1];                             // current target direction (dy,dx), unit length
let comHistory = [[0, 0], [0, 0], [0, 0], [0, 0]];
let lastCoM = null, initialMass = 0;
// The most recent crop readback, channel-packed (frame k in channel k). Held by reference: the
// engine reuses one buffer, and it is always consumed into a tensor before the next readback
// overwrites it.
let lastCrop = null;
let steps = 0, actions = 0, solitonDead = false;
let courseChanges = 0;          // target direction changes the user made this episode
// The pending auto-respawn from a death, so a manual Restart/New Maze/etc. during that gap can
// cancel it -- otherwise it would fire a second, unwanted respawn on top of the manual one.
let deathTimer = 0;
let levelWon = false;
// Game-over lives, shown as 💛 in the title row (see render()). Reset only on a genuine new game
// (placeSoliton()'s `resetLives` -- Restart, maze toggle, soliton picker, initial load, and the
// full restart handleDeath() falls back to once lives run out), not on the ordinary death-respawn
// in between (which is what makes them count down across deaths) or on clearing the board (which
// carries the same life count into the next board, like an arcade level clear -- see handleWin()).
const STARTING_LIVES = 3;
let lives = STARTING_LIVES;
// The level system: level N spawns the maze's ghosts numbered 1..N (see placeGhosts()), so the
// starting level is also the starting ghost count. Advances by one every win (see
// respawnAfterWin()) and resets on a genuine new game the same way lives do (placeSoliton()'s
// `resetLevel`, defaulting to `resetLives`) -- not on an ordinary death-respawn or on clearing the
// board, which is what makes the game harder round over round instead of every death. ?level=N to
// start somewhere other than 3 -- also what a reset falls back to, so it holds across a Restart.
const STARTING_LEVEL = intParam('level', 3);
let level = STARTING_LEVEL;
// The pending auto-restart from clearing the board, same shape as deathTimer.
let winTimer = 0;
// Frightened is tracked per ghost (see newGhost()'s `frightened` field), not as one global mode:
// eating a pellet marks every currently-not-already-frightened ghost, but a ghost that then gets
// eaten and respawns (respawnGhost() -> newGhost()) ends its own window immediately and comes back
// chasing, while its packmates keep counting down theirs. frightenedTimer is the single shared
// duration they all started from -- cleared and restarted on every pellet eaten while it's already
// running, so a second pellet extends the window rather than stacking a second one behind it.
let frightenedTimer = 0;
// Set by frightenedTimer's callback, which fires between agentStep() calls with no fresh `rb`
// (ghost tile-local positions) to turn ghosts with -- it only flags that the window has ended.
// updateFrightened() is what actually reacts to it, on whichever step notices it next, turning
// every *still*-frightened ghost back toward Pac-Man and clearing the flag (an already-normal one,
// from an earlier respawn, is left alone -- it ended its own window already).
let frightenedExpired = false;
// True while the sim holds for the eat_ghost.wav jingle -- set by handleGhostEaten(), cleared once
// the sound actually finishes (or fails to play at all). Checked in loop() alongside
// solitonDead/levelWon; unlike those two this one has no timer of its own, since it's the sound's
// own 'ended' event that ends the hold.
let ghostEatPause = false;
// The eat_dot_0/1 "waka waka" loop. eatActive is true while the pair is alternating; eatToggle
// picks which of the two plays next; eatDeadline (performance.now()-based) is when it's allowed
// to stop -- every eaten dot pushes it out by another full CFG.eatSoundGraceMs. The deadline is
// checked only at the one place a new play is actually started (see playNextEatSound), rather
// than off a separate timer callback that could fire mid-playback and race a restart -- that is
// what keeps this to one sample audible at a time.
let eatActive = false, eatToggle = 0, eatDeadline = 0;
let running = false, session = null, busy = false, lastMs = 0, lastAction = null, lastQ = null;
// True from the moment a spawn calls playStartSound() until its jingle actually finishes --
// covers both "waiting for the first gesture to unlock audio" and "jingle audibly playing".
// setRunning() refuses to unpause while this is set, so Play/Space/etc. can't cut the intro short.
let introPlaying = false;
// Recent interventions, newest last: {r,c,sign,step}. Each marker fades out over
// CFG.actionFadeSeconds of *simulation* time, so the trail reads the same at any sim speed
// (and freezes mid-fade rather than vanishing when the run is paused).
let actionTrail = [];
let bank = [], currentIndex = 0;
let mazeEnabled = true;
// Chrome on/off: the control deck and CARL's overlay (direction arrows, intervention discs) hide
// together, so the board can be watched as a game rather than as an instrumented demo.
let showOverlay = true;
// Who intervenes on the soliton, one of 'sometimes' (default), 'always', or 'human'. In 'human'
// mode CARL's policy is never queried: the sim still advances, but the only actions on the board
// are the ones the user clicks in. 'sometimes' is the performance mode: CARL is only queried for
// sometimesWindow steps after an episode starts or after the user steers, then goes idle (no
// inference, no action) until the next steer -- inference is the expensive part of a step, so this
// is much cheaper to run than 'always' while looking the same whenever the player is engaged.
let actorMode = 'sometimes';
let humanActs = false;                 // derived from actorMode === 'human', kept for readability
let sometimesWindow = 100;             // configurable steps CARL stays active for in 'sometimes' mode
let sometimesRemaining = 0;            // steps left in the current active window ('sometimes' mode only)
// The user's queued action, {gx,gy,sign}, or null. At most one is ever held: a step consumes it
// exactly where agentAct() would have run, so a human turn and a CARL turn are the same turn.
let pendingAction = null;
// The 90deg rotation the soliton spawns with. Rolled once per maze, not per spawn, so Respawn
// re-runs the same starting configuration instead of quietly changing the soliton's heading.
let spawnRotation = 0;
let sps = 60, stepAcc = 0, lastT = 0, measSps = 0, rateSteps = 0, rateTime = 0;
// A batch of steps yields once it has spent this much of the frame, rather than running a fixed
// number of steps. The fixed cap this replaces (16) was sized for a step that is pure GPU work,
// which is what a step costs while *you* are acting. With CARL acting every step also pays for an
// inference -- tens of milliseconds on a phone -- so 16 of them ran back to back for most of a
// second, and since render() only lands after the batch, the board visibly froze between repaints
// even though the step rate itself was tolerable. Budgeting by wall-clock time instead spends the
// same frame on however many steps actually fit: the sim runs slower on a slow device rather than
// in lurches, and the rate readout reports what was really achieved. ?budget=N to tune on a device.
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
const sndEat = [new Audio('assets/sound/eat_dot_0.wav'), new Audio('assets/sound/eat_dot_1.wav')];
// Registered once, permanently, rather than per-play: 'ended' only fires on a natural finish
// (never from stopEatingSound()'s pause()), so a one-shot listener re-added on every play would
// pile up whenever a play gets interrupted before it can fire and be removed -- e.g. a death or
// Restart landing mid-note -- and each stale listener left behind means one more concurrent call
// into playNextEatSound() the next time that same element does finish naturally, which is exactly
// what plays two tracks over each other. playNextEatSound() itself is what decides whether to
// keep going, so this listener just hands control back to it every time.
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
//  Maze layout -- a fixed, hand-authored maze, replacing the randomized-DFS generator the
//  CARL demo shipped with. The geometry it produces is the same one that generator built:
//  wide open *cells* (cw px) separated by thin *walls* (ww px), i.e. walls live between
//  cells rather than occupying cells of their own. Only the connectivity is now authored
//  instead of rolled.
//
//  That in-between geometry is why the layout is written on a doubled grid: an RxC maze is
//  (2R+1) x (2C+1) characters, the standard maze-ASCII form, which is also the shape real
//  Pac-Man level data takes -- so the 5x11 maze below is authored as 11x23 characters.
//
//        col:  0 1 2 3 4 5 6 7 8 9 10        odd index  -> cell
//    row 0     # # # # # # # # # # #         even index -> the wall slot between two cells
//    row 1     # C . . . . . . . . #         (row,col) both odd  -> cell (r,c) = ((row-1)/2, (col-1)/2)
//    row 2     # . # # # . # # # . #         both even           -> corner post, always wall
//
//  Two characters are pure authoring aliases, there so a hand-edited row stays readable as ASCII
//  art: '+' means exactly '.', and '-' means exactly ' '. Nothing reads them differently. '^' is
//  a third of the same shape -- open floor, exactly like '-' -- except to a ghost standing on it,
//  which is sent north instead of choosing for itself. That is what lets a centre room have a
//  one-way exit: without it a ghost can rattle around inside a pocket indefinitely.
//
//  '.' (or '+') places one dot, and 'O' one power pellet -- both a channel-2 soliton, see
//  placeDots() -- wherever they appear, cell slot or wall slot alike: on a cell slot it centers on
//  the tile, on a wall slot it centers on the gap between the two tiles either side (a pellet
//  sitting in a corridor, same as real Pac-Man). A '.'/'+'/'O' on a corner post (both indices even)
//  is impossible to satisfy -- corner posts are always wall -- and is ignored.
//
//  Characters:
//    on a cell slot:  'C' Pac-Man's spawn · '1'-'9' a numbered ghost's spawn (see placeGhosts()
//                     and CFG/`level` below -- the current level only spawns ghosts numbered at or
//                     below it) · '#' solid (filled) cell · 'O' a power pellet -- same channel-2
//                     soliton as a '.'/'+' dot, under the same shared rule, just recoloured (see
//                     placeDots()'s power mask): the rule fixes one equilibrium size for everything
//                     in that channel, so a bigger *stamp* would just relax back down to ordinary
//                     dot size rather than stay distinct · anything else (' ', '-', '^') is open
//                     floor, no dot. 'C' and a digit are open floor too; they only mark what spawns
//                     on the tile.
//    on a wall slot:  '#' (or '|') wall · anything else (' ', '-', '^', '.', '+', 'O') is an open
//                     passage between the two neighbouring cells -- including on the outer
//                     ring (row/col 0 and row/col 2R/2C): the sim wraps toroidally regardless
//                     of walls, so an opening there is a real Pac-Man-style side tunnel to the
//                     opposite edge, not a dead end. Corner posts (both indices even, e.g.
//                     (0,0)) are always wall, on the border same as in the interior.
//
//  Sizing: cw is derived from the board so the maze always fills it -- so the *layout* fixes
//  how many cells there are, and the board size fixes how many pixels each one gets. Keep cw
//  above the soliton's width or it clips the walls on every turn: the curated solitons run
//  30-46px across (models/solitons_direction.json, rule58 is the 46px outlier). The 5x11 maze
//  on the 550x250 board lands at cwX=40, cwY=39 -- near enough square, with wall thickness
//  matching CARL-WebGL's original generator (see MAZE_WW below). Solitons wider than cw are
//  clipped to fit -- see placeSoliton().
// ====================================================================================
const MAZE_WALL_CHARS = '#|';
// The four headings, in the order rotate90()'s k uses, so (k+1)%4 is a right turn and (k+3)%4 a
// left one, and the difference between two headings is the number of quarter turns between them.
// That *relative* part is exact and rule-independent. What is not is the absolute part: which way
// a freshly stamped pattern travels depends on how its bank entry happens to be drawn -- rule74's
// canonical heading is right, rule73's is left. So the ghost's heading is never assumed from the
// rotation it was stamped with; it is read back off its own motion (see ghostHeading()).
const DIRS = [[0, 1], [1, 0], [0, -1], [-1, 0]];
// Cells the ghost is allowed to change direction on. These are the same characters that already
// mean "dot" ('+'), "open" ('-', '^'), "ghost spawn" ('1'-'9') and "power pellet" ('O'), doing
// double duty as turn markers -- which costs nothing, because in this layout they land on exactly
// the 32 corner/junction cells and on no straight corridor cell at all.
const GHOST_TURN_CHARS = '+-O^123456789';
// ...and the one that does not leave the choice open: a ghost reaching it is sent north.
const GHOST_NORTH_CHAR = '^';
// Wall thickness and outer border, in pixels -- matches the randomized generator this replaced
// (CARL-WebGL's original ww=9, edgeWall=9). At 5x5 cells that leaves cw=39, tighter than the
// widest curated soliton (46px, models/solitons_direction.json) -- see placeSoliton()'s clip.
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
  "#O...6#+.......+#8...O#",
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
  // Whatever the floor()s above leave over is spread as extra margin, keeping the maze centered.
  const offX = MAZE_EDGE + Math.max(0, Math.floor((bw - 2 * MAZE_EDGE - (cols * stepX - MAZE_WW)) / 2));
  const offY = MAZE_EDGE + Math.max(0, Math.floor((bh - 2 * MAZE_EDGE - (rows * stepY - MAZE_WW)) / 2));
  const cellTop = r => offY + r * stepY, cellLeft = c => offX + c * stepX;

  // Pixel centre of a slot on one axis, whatever its parity: an odd index is a cell and resolves
  // to that cell's middle; an even index is the gap in front of cell i/2 and resolves to the
  // middle of that gap -- the outer margin at the two ends, the wall slot everywhere else. Used
  // only for dots -- walls and the spawn stay confined to their own slot parity, as documented
  // above.
  const slotCenter = (i, cellStart, cw, count, span) => {
    if (i & 1) return cellStart((i - 1) >> 1) + (cw >> 1);
    const k = i >> 1;
    const lo = k === 0 ? 0 : cellStart(k - 1) + cw;
    const hi = k === count ? span : cellStart(k);
    return (lo + hi) >> 1;
  };
  const gridPos = (gr, gc) =>
    [slotCenter(gr, cellTop, cwY, rows, bh), slotCenter(gc, cellLeft, cwX, cols, bw)];

  // Same carve-out-of-solid approach as the generator this replaces: start all wall, open up
  // the cells, then open the passages between them.
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
  // Sorted by id, not left in the row-major scan order above -- ghostCells.map(cellCenter) below
  // feeds maze.ghosts to placeGhosts() in this order, and that order is what fixes which array
  // index (and so which colour, see render()'s ghostColors) a given numbered ghost always gets,
  // regardless of which levels include it.
  ghostCells.sort((a, b) => a.id - b.id);

  // Dots and power pellets are read off the whole doubled grid, cell slots and wall slots alike,
  // since '.'/'+'/'O' are all valid on either (see the layout comment above). Kept as two separate
  // lists rather than one tagged list: placeDots() stamps them with the same soliton and rule, so
  // there is nothing a caller would do differently per-entry except which list it came from.
  const dotPos = [], powerPos = [];
  for (let gr = 0; gr < gh; gr++) for (let gc = 0; gc < layout[gr].length; gc++) {
    const ch = at(gr, gc);
    if (ch !== '.' && ch !== '+' && ch !== 'O') continue;
    if (gr % 2 === 0 && gc % 2 === 0) continue;   // corner post: always wall, can't hold a dot
    (ch === 'O' ? powerPos : dotPos).push(gridPos(gr, gc));
  }
  if (enableWalls) for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    // Only the gap itself is carved, not the union of the two cells -- so a passage next to a
    // solid cell opens the connector without also hollowing the cell out.
    if (c + 1 < cols && !isWall(at(2 * r + 1, 2 * c + 2)))
      carve(cellTop(r), cellLeft(c) + cwX, cellTop(r) + cwY, cellLeft(c + 1));
    if (r + 1 < rows && !isWall(at(2 * r + 2, 2 * c + 1)))
      carve(cellTop(r) + cwY, cellLeft(c), cellTop(r + 1), cellLeft(c) + cwX);
  }
  // The outer ring of wall slots (row/col 0 and row/col gh-1/gw-1) works the same way -- an
  // open one carves its cell's gap all the way out to the physical board edge, through the
  // MAZE_EDGE margin, rather than to a neighbouring cell (there isn't one). sim.glsl already
  // wraps the board toroidally regardless of walls, so a corridor opened clear to the edge on
  // both sides of the board is a real Pac-Man-style side tunnel to the opposite edge, not just
  // a dead end. Corner posts ((0,0) etc.) are never carved, same as interior ones.
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

  // Corner posts are never carved above -- they stay exactly the fixed WW-ish square the initial
  // fill() left them as, at every grid intersection (interior and border alike). Round each one
  // down to a circle inscribed in its own box only on the side(s) that face open floor; a side
  // that instead continues into a solid wall bar is left full square, so that bar's full-width
  // flat end merges flush with the post rather than butting against a circle that -- being
  // inscribed in a box exactly as wide as the bar -- would otherwise only touch that bar's edge
  // at a single tangent point, notching out the rest of the bar's corners right where they meet
  // the post (visible as a chip/pinch). A corner rounds only when BOTH of its two adjacent sides
  // are open: that is precisely "elongate the wall into the post by the circle's own radius", since
  // leaving a solid-facing quadrant untouched is the same as the bar already reaching the post's
  // centre line, where the inscribed circle is at its full-width equator.
  if (enableWalls) for (let i = 0; i <= rows; i++) for (let j = 0; j <= cols; j++) {
    const y0 = i === 0 ? 0 : cellTop(i - 1) + cwY, y1 = i === rows ? bh : cellTop(i);
    const x0 = j === 0 ? 0 : cellLeft(j - 1) + cwX, x1 = j === cols ? bw : cellLeft(j);
    const cy = (y0 + y1) / 2, cx = (x0 + x1) / 2, rad = Math.min(y1 - y0, x1 - x0) / 2;

    // Sampled one pixel outside the post's own box, toroidally wrapped (the board wraps, so the
    // post's true neighbour past a board edge is the opposite edge, not "nothing").
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

  // Per-cell decision table for the ghosts: the layout character on the cell, which of the four
  // DIRS headings can leave it, and its pixel centre. Built here because this is the one place
  // that knows both the layout characters and the pixel geometry they resolve to.
  const cells = [];
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    const gr = 2 * r + 1, gc = 2 * c + 1;
    const [cy, cx] = cellCenter([r, c]);
    cells.push({ ch: at(gr, gc), open: DIRS.map(([dy, dx]) => !isWall(at(gr + dy, gc + dx))), cy, cx });
  }
  // Which cell a board pixel falls in, or -1 outside the grid (the outer margin, where a tunnel
  // wraps). The wall gap past a cell reads as that cell rather than as a slot of its own; callers
  // gate on distance to the centre anyway, so the gap never satisfies them.
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
//  Soliton placement: exact 90deg rotation (no interpolation) + place its own center of
//  mass at the maze start, clipped to the spawn cell's cwX x cwY box if the soliton is
//  wider than that (its edges are empty or near-empty, so the clip costs it little), then
//  zero out anything that still lands on a wall.
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

// Resamples a soliton's pixel data to `scale` of its native size, bilinearly. A Lenia pattern is
// tied to the kernel radius it was found at -- shrinking just the kernel radius without also
// shrinking the pattern breaks it (it's no longer the same normalized neighbourhood the pattern
// self-organized around), so the dots channel's smaller kernel radius (see placeDots()) and this
// resample always move together, by the same factor.
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

// Writes one soliton into `arr` with its own centre of mass landing on (tr,tc). `mask`, if given,
// gets a 1 written at every cell the stamp actually touches -- how placeDots() tells a power
// pellet's footprint apart from an ordinary dot's for colouring, since both are the same soliton
// under the same rule and nothing about the mass values themselves distinguishes one from the other.
function stampSoliton(arr, entry, tr, tc, rotation, mask) {
  const rot = rotate90(entry.flat, entry.h, entry.w, rotation);
  let sy = 0, sx = 0, sm = 0;
  for (let y = 0; y < rot.h; y++) for (let x = 0; x < rot.w; x++) {
    const v = rot.data[y * rot.w + x]; if (v > 0) { sm += v; sy += v * y; sx += v * x; }
  }
  const cy = sy / sm, cx = sx / sm;

  // Clip window, centered on the soliton's own centroid, sized to the spawn cell -- anything
  // farther than half the cell width/height from (cy,cx) is dropped rather than deposited, so a
  // soliton wider than the cell is trimmed to fit instead of overhanging into a neighbouring
  // passage (the caller's wall zeroing only catches an overhang that lands on a wall, not one
  // that lands on an open corridor next door).
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

// The dots channel: one soliton on every '.'/'+' of the layout, plus one on every 'O' (a power
// pellet -- same soliton, same rule, just recoloured; see the mask below), all under a single
// shared rule. CARL never sees this channel and never acts on it -- it just runs. Rebuilt whenever
// the agent's soliton is placed, so Restart and toggling the maze restore the full set.
function placeDots() {
  if (!DOTS_ENABLED) return;
  // No dots or pellets in the layout means no second channel at all: leaving its rule unset is
  // what keeps SimGL.step() from paying for a second convolution over an empty board.
  const entry = (maze.dots.length || maze.power.length)
    ? bank.find(b => b.name === CFG.channel2RuleName) : null;
  if (!entry) return;
  SimGL.setRule2({ mu: entry.mu, sigma: entry.sigma, betas: entry.betas, R: CFG.channel2R, dt: CFG.dt });

  const scaled = resizeSoliton(entry, CFG.channel2R / CFG.R);
  const arr = new Float32Array(N);
  // 1 wherever a power pellet's stamp landed, so draw.glsl can colour it apart from a plain dot --
  // the two are otherwise the same mass under the same rule (see stampSoliton()).
  const power = new Uint8Array(N);
  for (const [r, c] of maze.dots) stampSoliton(arr, scaled, r, c, 0);
  for (const [r, c] of maze.power) stampSoliton(arr, scaled, r, c, 0, power);
  for (let i = 0; i < N; i++) if (wall[i]) { arr[i] = 0; power[i] = 0; }
  SimGL.uploadState2(arr);
  SimGL.setPowerMask(power);
}

// The ghost channel: one soliton on every 'M' of the layout, all under a single shared rule and,
// like the dots, running free -- CARL neither sees this channel nor steers it. The one thing it
// does to the game is erase Pac-Man's mass wherever it overlaps him (glsim.js's ghostEat), which
// is also how a ghost kills: enough of him erased and finishStep() reads the mass below
// CFG.massDeathFraction and calls it a death, through the ordinary death path.
//
// Rebuilt on every spawn, death respawns included, so the ghosts always start back in their house
// rather than wherever they had drifted to when Pac-Man died.
//
// Each ghost is simulated in a private GHOST_WINDOW-square tile rather than in a shared field.
// Two ghosts in one field are two Lenia solitons, and two Lenia solitons that meet annihilate or
// blow up; tiles let them pass through each other the way arcade ghosts do. A tile is a window
// onto the board -- it carries a board origin, so walls still apply -- and it slides by whole
// cells each step to keep its soliton centred, whole cells because a fractional slide would mean
// resampling the pattern away. See shaders/ghostsim.glsl.
const GHOST_WINDOW = SimGL.GHOST_WIN;   // the tile edge the engine allocates

let ghosts = [];              // one entry per 'M'; [] means there are none to judge or steer

function newGhost(r, c, cellIdx, home = [r, c]) {
  return {
    home,                        // its own numbered cell -- where a non-eaten respawn goes back to
    y: r, x: c,                 // last known board CoM
    origin: [Math.round(r) - (GHOST_WINDOW >> 1), Math.round(c) - (GHOST_WINDOW >> 1)],
    shift: [0, 0],              // whole cells its tile slides next step, to re-centre it
    cell: cellIdx,
    // The house is itself a marked cell and the ghost spawns dead on its centre, so without
    // `turned` it would count as having just passed the centre and turn on its very first step.
    turned: true,
    prevD2: Infinity,           // squared distance to that cell's centre one step ago
    anchor: [r, c],             // CoM on entering `cell` -- the baseline its heading is read from
    // Stamped in the bank's own orientation, whatever that is, and turned to face the way out of
    // the house as soon as it has moved far enough to say which way it is going -- a few px, well
    // short of the ~20 to the wall. Aiming it at stamp time would mean knowing the rule's
    // canonical heading, which varies by rule.
    aim: cellIdx >= 0 ? maze.cells[cellIdx].open.indexOf(true) : -1,
    spawnMass: 0,               // 0 disarms this ghost's death check
    // Per-ghost, not global: a ghost eaten mid-frightened-window respawns here with this false
    // (see respawnGhost(), which is just this same newGhost()), ending its own window immediately
    // while its packmates -- untouched -- keep counting down theirs. See updateFrightened().
    frightened: false,
  };
}

// Builds one ghost's tile: the soliton centred in it, masked against the maze at the board cells
// the tile currently covers. Returns the tile's mass, which is what arms its death check -- summed
// here rather than waited for from a readback, so the check is live from the step it appears.
function buildGhostTile(entry, g) {
  const win = GHOST_WINDOW, half = win >> 1;
  const tile = new Float32Array(win * win);
  // stampSoliton() works in board coordinates, so stamp into a board-sized scratch and cut the
  // tile out of it -- that keeps one implementation of the centring and clipping, not two.
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
  // Cleared first, and only refilled once ghosts are actually on the board: an entry's spawnMass
  // is what arms its death check, so every path that leaves the channel empty must also leave that
  // disarmed -- otherwise "no ghost" reads as "dead ghost" and respawns on a loop, every step.
  ghosts = [];
  SimGL.setGhostTiles([], []);
  if (!GHOSTS_ENABLED) return;
  // The level system: level N spawns every numbered ghost at or below N, so level 3 (the starting
  // level) is ghosts 1-3, level 4 adds ghost 4, and so on -- and a level past the highest number
  // the layout actually has just spawns all of them, since the filter below has nothing left to
  // exclude (see updateLevel()'s comment for "if level number is higher than the max ghost").
  const spawns = maze.ghosts.filter(([, , id]) => id <= level);
  const entry = spawns.length ? bank.find(b => b.name === CFG.channel3RuleName) : null;
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
//  A free-running soliton already travels in a straight line on its own, so the only thing that
//  needs deciding is what happens at a corner or crossroad. Rather than steer the turn with the
//  policy -- a second inference every step, which is what a step actually costs -- the ghost's own
//  tile is rotated a quarter turn about its centre of mass (rotate.glsl). The turn is then instant
//  and always succeeds, which is also how an arcade ghost turns.
//
//  Rotating rather than re-stamping is what keeps the turn from reading as a blink: the pattern
//  that comes out the far side is the one that went in, mass for mass, not a fresh copy of the
//  canonical soliton from the bank.
// ------------------------------------------------------------------------------------
// A ghost's heading, read off how it has actually moved since its anchor rather than inferred from
// the rotation it was stamped with -- which would be wrong for any rule whose bank entry does not
// happen to face right. -1 until it has travelled far enough for the answer to be trustworthy:
// these solitons run at 0.3-0.45 px/step and up to ~30deg off their nominal axis, so a couple of
// px is not a heading. Six is, and there are ~20 to cross before a heading is needed.
const GHOST_HEADING_MIN_PX = 6;
function ghostHeading(g) {
  const [dy, dx] = toroidalDelta(g.y, g.x, g.anchor[0], g.anchor[1]);
  if (dy * dy + dx * dx < GHOST_HEADING_MIN_PX * GHOST_HEADING_MIN_PX) return -1;
  return Math.abs(dx) > Math.abs(dy) ? (dx > 0 ? 0 : 2) : (dy > 0 ? 1 : 3);
}

// Which of `choices` heads most directly at (sign>0) or away from (sign<0) Pac-Man from (gy,gx).
// Scored against the toroidal vector to him, so a heading whose corridor leaves by a side tunnel
// is judged on where it actually comes out rather than on the long way round. -1 if he is not on
// the board to be chased or fled from. Fleeing reuses the exact same scoring, just negated --
// "most away" is "least toward" -- so frightened ghosts and normal ones share one function.
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

// A ghost dies the same two ways Pac-Man does -- dissolved away, or grown past the explode limit --
// judged against its own spawn mass on the same CFG thresholds. What differs is the consequence:
// no jingle, no held board, no episode end, and no effect on the other ghosts. It is simply put
// back in its own house. The explode limit stays meaningful inside a tile, which a window on
// Pac-Man's channel would not: a 96-square tile holds up to 9216 mass against the 1500 threshold,
// so a runaway still trips it long before it saturates.
function ghostDied(g, s) {
  if (!g.spawnMass) return false;         // not armed -- nothing to have died
  return !s.valid || s.mass > CFG.massExplodeLimit || s.mass < CFG.massDeathFraction * g.spawnMass;
}

// Puts one ghost back on the board, leaving the rest of the pack running. Rewriting its tile is
// also the whole of the cleanup after an explosion: a tile has hard edges, so however far the mess
// spread it is still inside that one tile, and the board-space texture everything downstream reads
// is rebuilt from the tiles every step rather than accumulated -- so there is nowhere else for
// debris to have got to, and nothing else to scrub.
//
// `eaten` (a frightened-mode kill, as opposed to an explosion) sends it to symbol '1''s cell -- the
// one actual ghost house, inside the centre room -- rather than back to its own numbered spawn: the
// numbered cells are scattered around the maze (see MAZE_LAYOUT), not all of them a "house" a ghost
// could plausibly walk out of again. Its own numbered cell is passed through as `home` regardless,
// so a later non-eaten respawn (or a future level filtering it back out and back in) still knows
// where it actually belongs.
function respawnGhost(i, eaten) {
  const entry = bank.find(b => b.name === CFG.channel3RuleName);
  if (!entry) return;
  const home = ghosts[i].home;
  const house = maze.ghosts.find(([, , id]) => id === 1);
  const [r, c] = eaten && house ? house : home;
  const g = newGhost(r, c, maze.cellIndexAt(r, c), home);
  ghosts[i] = g;
  pushGhostTiles();          // its origin is back at the spawn point before the tile is written
  const { tile, mass } = buildGhostTile(entry, g);
  g.spawnMass = mass;
  SimGL.uploadGhostTile(i, tile);
}

// Flips ghost `i` if it is currently heading the "wrong" way for the mode `away` names -- toward
// Pac-Man when it should be fleeing him, or away from him when it should be chasing again -- by a
// hard 180. Left alone if it's already heading the right way, or hasn't moved far enough since its
// last turn to have a readable heading at all (ghostHeading() returns -1): "either no change or
// turn 180", never a left/right correction. The pivot is `s.localRow/localCol`, the tile-local
// position this same readback just reported -- exactly what SimGL.rotateGhost() needs, and current
// regardless of anything steerGhost() does to `g.origin`/`g.shift` afterwards for the *next* step.
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

// Starts/extends the shared frightened window on a pellet eaten, and ends it on the step that
// first notices frightenedExpired -- in both cases turning (see turnGhost180()) only the ghosts
// that transition: a pellet eaten while some ghosts are already frightened just extends their
// clock, no re-turn, and an expiry only turns whichever ghosts are *still* frightened, skipping
// any that already ended their own window early via respawnGhost() (which resets g.frightened to
// false the moment a frightened ghost gets eaten and returns home). SimGL.step() already used the
// *previous* per-ghost frightened flags for this step's eat checks (see pushGhostTiles()), so a
// pellet eaten just now takes effect next step -- the same one-step lag ghost death detection
// already has.
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

function steerGhosts(rb) {
  if (!GHOSTS_ENABLED || !ghosts.length) return;
  for (let i = 0; i < ghosts.length; i++) {
    const s = rb.ghosts[i];
    if (!s) continue;
    if (ghostDied(ghosts[i], s)) {
      const eaten = ghosts[i].frightened;   // dissolved while frightened -- Pac-Man ate it
      respawnGhost(i, eaten);               // just this one; the rest run on
      if (eaten) handleGhostEaten();
      continue;
    }
    if (s.valid) steerGhost(ghosts[i], i, s, rb);
  }
  // Origins and shifts changed above; the engine needs them before the next step's sim.
  pushGhostTiles();
}

function steerGhost(g, idx, s, rb) {
  const win = GHOST_WINDOW, half = win >> 1;
  // The engine reports where the soliton sits inside its tile; the board position is that plus the
  // tile's origin. The tile then slides by whole cells to put it back in the middle, which is what
  // keeps it from ever reaching an edge -- the shift takes effect on the next step's sim.
  g.y = ((g.origin[0] + s.localRow) % H + H) % H;
  g.x = ((g.origin[1] + s.localCol) % W + W) % W;
  const sy = Math.round(s.localRow) - half, sx = Math.round(s.localCol) - half;
  g.shift = [sy, sx];
  // Wrapped into board range, not left to accumulate: ghostblit.glsl locates a tile by shifting
  // it at most one board-width to find the near copy, so an origin that has drifted further than
  // that (e.g. after several trips around a toroidal edge) stops matching any pixel at all -- the
  // ghost goes invisible until further wandering drifts it back into range on its own.
  g.origin = [(((g.origin[0] + sy) % H) + H) % H, (((g.origin[1] + sx) % W) + W) % W];

  const i = maze.cellIndexAt(g.y, g.x);
  if (i !== g.cell) {
    g.cell = i; g.turned = false; g.prevD2 = Infinity;
    g.anchor = [g.y, g.x];                  // fresh baseline: this cell's run is the heading
  }

  // Pivot for any rotation below, in tile-local coordinates and integer -- the rotation is only
  // exact about a whole cell.
  const pivX = Math.round(s.localCol), pivY = Math.round(s.localRow);

  // The one-off turn out of the spawn house, once there is enough motion to say which way it is
  // currently pointing. Everything below is the ordinary junction logic and does not apply yet.
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

  // Turn at the closest approach to the tile centre, not on first crossing into some band around
  // it. The distance falls while the ghost runs in and rises once it is past, so the step where it
  // stops falling is the centre crossing itself -- which is where a turn has to happen, since the
  // rotation pivots on the CoM and anywhere else swings the soliton into the corner it is turning
  // around. A distance threshold cannot do this: it fires on entry to the band, a third of a tile
  // early, and tightening it enough to land on the centre makes it small enough to step over.
  const dy = g.y - cell.cy, dx = g.x - cell.cx;
  const d2 = dy * dy + dx * dx;
  if (d2 > g.prevD2) {
    g.turned = true;
    const cur = ghostHeading(g);
    if (cur < 0) return;      // too little travel to read a heading -- leave it running straight

    // Left, straight and right -- never a reversal, so the ghost reads as patrolling rather than
    // dithering. A '^' cell overrides all of that and sends it north, which is what gives a centre
    // room a one-way exit. Whatever the candidates, only the ones that are actually open survive;
    // if none does (a dead end) reversing is all that is left.
    const wanted = cell.ch === GHOST_NORTH_CHAR ? [3] : [(cur + 3) % 4, cur, (cur + 1) % 4];
    const open = wanted.filter(k => cell.open[k]);
    let next;
    if (!open.length) next = (cur + 2) % 4;
    else {
      // CFG.chaseBias of the time it takes the opening that closes on (or, frightened, opens away
      // from) Pac-Man; the rest of the time it rolls. Biasing only the turns still compounds hard,
      // because every junction is another chance to correct -- which is why this is the difficulty
      // dial and not the speed.
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

// resetDots is false only for the automatic respawn after a death: the dots channel just keeps
// running with whatever it already had (dots already eaten stay eaten). Every user-triggered
// (re)spawn -- Restart, New Maze, the maze toggle, the soliton picker -- restores the full set.
// `resetLives` defaults to `playIntro` -- every ordinary fresh start (Restart, New Maze, the
// soliton picker, initial load, and the full restart handleDeath() falls back to once lives run
// out) refills lives back to STARTING_LIVES same as everything else about a fresh start. The one
// exception is handleWin() clearing the board: that's still a fresh spawn (dots refilled, jingle
// played) but not a new *game* -- clearing the board carries the player's remaining lives into the
// next one, same as a level clear would in the arcade original -- so it passes false explicitly.
function placeSoliton(entry, resetDots = true, playIntro = true, resetLives = playIntro, resetLevel = resetLives) {
  mu = entry.mu; sig = entry.sigma; betas = entry.betas.slice();
  SimGL.setRule({ mu, sigma: sig, betas, R: CFG.R, dt: CFG.dt * CFG.channel1Speed });
  // Resolved before placeGhosts() below, which reads `level` to decide which numbered ghosts to
  // spawn -- respawnAfterWin() has already bumped it by the time this runs, and a reset here must
  // land before that same call, not after it.
  if (resetLevel) level = STARTING_LEVEL;

  // Built once on the CPU and uploaded; from here on the board only exists on the GPU.
  const arr = new Float32Array(N);
  const [tr, tc] = maze.start;
  stampSoliton(arr, entry, tr, tc, spawnRotation);
  for (let i = 0; i < N; i++) if (wall[i]) arr[i] = 0;    // applyWallCollision, at spawn
  SimGL.uploadState(arr);
  if (resetDots) placeDots();
  placeGhosts();     // unconditional: a ghost that ate Pac-Man must not still be sitting on his
                     // spawn tile when he comes back, or the respawn is eaten on arrival

  clearTimeout(deathTimer);
  clearTimeout(winTimer);
  clearTimeout(frightenedTimer);
  stopEatingSound();
  sndEatGhost.pause(); sndEatGhost.currentTime = 0;   // in case a respawn lands mid-jingle
  steps = 0; actions = 0; courseChanges = 0; solitonDead = false; levelWon = false; ghostEatPause = false;
  frightenedExpired = false;   // placeGhosts() above already gave every ghost a fresh, unfrightened newGhost()
  if (resetLives) lives = STARTING_LIVES;
  sometimesRemaining = sometimesWindow;
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
  render();
  // Skipped for the automatic post-death respawn: that one carries on the same run rather than
  // starting a fresh one, so it gets no jingle and no pause -- the sim just keeps going the
  // instant solitonDead clears above.
  if (playIntro) playStartSound();
}
// Held paused through the start jingle -- the board is up and visible but frozen -- then the run
// begins the instant it ends. Runs after every manual (re)spawn (initial boot, Restart, maze
// toggle, soliton picker) -- not the automatic respawn after a death, see placeSoliton().
//
// Browsers refuse to play audio with sound until the page has seen a user gesture. Restart/maze
// toggle/etc. are themselves triggered from inside a click, so they're already past that gate and
// play immediately -- but the very first call, from page load, has no gesture behind it yet and
// gets rejected. Rather than starting silently in that case, show a prompt and wait for the
// page's actual first gesture (a click, a tap, any key), then retry -- that attempt is inside a
// real gesture, so it succeeds -- before starting the game.
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
      const start = () => {
        window.removeEventListener('pointerup', onPointer);
        window.removeEventListener('keydown', onKey);
        hideStartPrompt();
        attempt()
          .then(() => sndStart.addEventListener('ended', finish, { once: true }))
          .catch(finish);   // blocked even inside a gesture -- give up silently, but still start
      };
      // pointerup, not pointerdown: iOS Safari (and other strict mobile browsers) only counts a
      // *completed* tap -- touchend/pointerup/click -- as the gesture that unlocks audio, not the
      // touch-start. Listening on pointerdown consumed the one-shot listener on the down-phase,
      // attempt() failed again for the same reason as the very first (gestureless) call, and it
      // fell straight through to the catch(finish) below -- silently starting the game with no
      // jingle on a tap, while every other spawn path (Restart, etc.) plays fine because a button
      // click is a real completed gesture on any browser.
      const onPointer = () => start();
      // Modifier keys (and Escape) don't count as a real "user activation" for autoplay purposes
      // -- a bare Alt/Ctrl/Shift/CapsLock press would otherwise fall straight through to the
      // catch() above and start the game silently. Keep listening past those instead of consuming
      // the one-shot gesture on them.
      const onKey = e => { if (!NON_ACTIVATING_KEYS.has(e.key)) start(); };
      window.addEventListener('pointerup', onPointer, { once: true });
      window.addEventListener('keydown', onKey);
    });
}
// Keys the HTML spec excludes from counting as a "user activation" gesture: the UI Events
// Modifier Keys table, plus Escape (excluded separately by the activation-triggering-input-event
// definition). Pressing only one of these keeps the game waiting rather than starting silently.
const NON_ACTIVATING_KEYS = new Set([
  'Alt', 'AltGraph', 'CapsLock', 'Control', 'Fn', 'FnLock', 'Hyper', 'Meta', 'NumLock', 'OS',
  'ScrollLock', 'Shift', 'Super', 'Symbol', 'SymbolLock', 'Escape',
]);
// Reuses the (otherwise unused) result banner element for a "waiting for the first click/key"
// notice -- same spot, same look.
function showStartPrompt() {
  const el = $('result');
  el.className = 'result';
  el.innerHTML = '<b>Click, tap, or press a key to start</b>';
  el.hidden = false;
}
function hideStartPrompt() {
  const el = $('result');
  el.hidden = true;
  el.innerHTML = '';
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
// The death-triggered respawn -- same spawn, but leaves the dots channel alone (see placeSoliton).
function respawnAfterDeath() { placeSoliton(bank[currentIndex], false, false); }
// The win-triggered respawn -- a full fresh spawn like respawnCurrentSoliton(), except lives carry
// over into the next board rather than refilling (see placeSoliton()'s `resetLives`), and the
// level advances by one first, so placeGhosts() -- called from inside placeSoliton() -- spawns the
// next level's ghost count. `level` isn't clamped to the layout's highest numbered ghost here: the
// filter in placeGhosts() (`id <= level`) just has nothing left to exclude once level passes it,
// so a level with no matching digit is silently a no-op rather than needing special-casing.
function respawnAfterWin() { level++; placeSoliton(bank[currentIndex], true, true, false); }

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
  else if (actorMode === 'sometimes' && sometimesRemaining <= 0) {
    lastMs = 0; lastQ = null; lastAction = null;      // CARL idle this step -- no inference, no action
  }
  else if (!session) return;
  else {
    action = await agentAct();
    if (actorMode === 'sometimes') sometimesRemaining--;
  }

  // Action, step, wall mask, CoM reduction and the next input crop, queued back to back with no
  // synchronization; the single readback below is the only point where the CPU waits on the GPU.
  // Per-ghost frightened state rides along on ghosts[i].frightened via pushGhostTiles() (called at
  // the end of steerGhosts() below), not as an argument here -- see glsim.js's gFrightened.
  SimGL.step(action);
  const rb = SimGL.readback();
  lastCrop = rb.crop;
  steps++;
  updateFrightened(rb);
  steerGhosts(rb);
  finishStep(rb);
}

// Bookkeeping shared by both actors: take the soliton's freshly computed center of mass and
// decide whether the episode has ended.
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
}
// Ignore floating-point noise from the reduction -- SimGL.readback().eaten is an exact sum of
// erased dot mass, not a heuristic, so this only needs to clear rounding error, not tune a
// detector.
const EAT_SOUND_MIN_MASS = 1e-4;
// The eat_dot_0/1 "waka waka" pair, alternating for as long as dots keep getting eaten. Every
// eaten dot pushes eatDeadline out by another CFG.eatSoundGraceMs; once it's passed, the chain's
// next scheduled play (not a separate timer) sees that and stops instead of queuing another
// sample -- so it always finishes the sample already playing rather than cutting one off, and
// there is never a moment where a fresh start and an old tail could both be sounding at once.
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
// Stops the loop outright (mid-note if need be) -- used when a new life starts, so a chomp left
// over from the previous one can't bleed into it.
function stopEatingSound() {
  eatActive = false;
  for (const snd of sndEat) { snd.pause(); snd.currentTime = 0; }
}
// Death is a beat, not a stop: no "Fail" screen and the Play/Pause button never flips, so the
// game never visibly pauses -- the sim just holds (loop() stops stepping while solitonDead)
// through the death jingle, then a further second of silence, then either the soliton respawns on
// its own (lives left) or the whole game restarts (see the `next` pick below), same as a manual
// Restart in the latter case.
const DEATH_PAUSE_MS = 1000;
// `exploded`, when true, is a mass-runaway death (rb.mass > CFG.massExplodeLimit) rather than an
// ordinary dissolve or ghost kill -- see finishStep(). That one doesn't cost a life: it isn't
// something the player could have steered around the way running into a ghost is, so it just
// respawns the soliton in place, same beat and jingle otherwise.
function handleDeath(exploded = false) {
  solitonDead = true;
  if (!exploded) lives = Math.max(0, lives - 1);
  // Lives run out: a full restart (placeSoliton()'s playIntro path) rather than the ordinary
  // death-respawn -- refills the dots, replays the start jingle, and resets lives right back to
  // STARTING_LIVES, same as hitting Restart by hand.
  const next = lives > 0 ? respawnAfterDeath : respawnCurrentSoliton;
  if (!SOUND_ENABLED) { deathTimer = setTimeout(next, DEATH_PAUSE_MS); return; }
  sndDeath.currentTime = 0;
  // The pause is timed off the jingle actually ending, not off starting it -- if playback is
  // blocked for some reason, fall back to the pause alone rather than never respawning.
  const afterSound = () => { deathTimer = setTimeout(next, DEATH_PAUSE_MS); };
  sndDeath.play()
    .then(() => sndDeath.addEventListener('ended', afterSound, { once: true }))
    .catch(afterSound);
}

// Clearing the board is also a beat rather than a stop: the sim holds (loop() stops stepping
// while levelWon) through the intermission jingle, then a further second of silence, then a fresh
// spawn -- same path as a manual Restart (start jingle, refilled dots channel; see
// placeSoliton()), except lives carry over rather than refilling, since this is a level clear, not
// a new game -- see respawnAfterWin().
const WIN_PAUSE_MS = 1000;
function handleWin() {
  levelWon = true;
  stopEatingSound();      // the last dot's chomp shouldn't bleed into the intermission jingle
  if (!SOUND_ENABLED) { winTimer = setTimeout(respawnAfterWin, WIN_PAUSE_MS); return; }
  sndIntermission.currentTime = 0;
  const afterSound = () => { winTimer = setTimeout(respawnAfterWin, WIN_PAUSE_MS); };
  sndIntermission.play()
    .then(() => sndIntermission.addEventListener('ended', afterSound, { once: true }))
    .catch(afterSound);
}

// Fire-and-forget, unlike the death/win/eat-ghost jingles -- the game keeps running underneath it,
// so it's just restarted from the top on every pellet eaten (including one eaten mid-window, which
// only extends the frightened timer -- see updateFrightened()) rather than queued or awaited.
function playFrightSound() {
  if (!SOUND_ENABLED) return;
  sndFright.currentTime = 0;
  sndFright.play().catch(() => {});
}

// The one moment eating a ghost differs from eating a dot: the board holds (loop() stops stepping
// while ghostEatPause) for the sound's own length, same beat as handleDeath()/handleWin() but with
// no further pause tacked on afterward and no consequence beyond the hold itself -- the eaten
// ghost already went home via respawnGhost(), called just before this from steerGhosts().
function handleGhostEaten() {
  ghostEatPause = true;
  const finish = () => { ghostEatPause = false; };
  if (!SOUND_ENABLED) { finish(); return; }
  sndEatGhost.currentTime = 0;
  sndEatGhost.play()
    .then(() => sndEatGhost.addEventListener('ended', finish, { once: true }))
    .catch(finish);
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
// One 👀 per ghost, centred on its CoM -- g.y/g.x are continuous (straight off the GPU tile
// reduction, see steerGhost()), not rounded pixel indices, so this uses the same no-offset mapping
// lastCoM's own arrows use rather than actionTrail's "+0.5 to recentre a rounded index" one.
function drawGhostEyes(sc) {
  if (!ghosts.length) return;
  octx.font = `${Math.max(7, (maze.cellY || 40) * sc * 0.5)}px sans-serif`;
  octx.textAlign = 'center';
  octx.textBaseline = 'bottom';
  // octx.fillStyle isn't reset between draws -- drawActionTrail() is the only other place that
  // touches it, always to a low-alpha rgba() for its fading action markers, and never sets it back
  // afterward. Without setting it here too, fillText() would inherit that leftover alpha instead
  // of drawing opaque, which is exactly why the eyes faded whenever CARL had recently acted.
  octx.fillStyle = '#fff';
  for (const g of ghosts) octx.fillText('👀', g.x * sc, g.y * sc);
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
  // Per ghost, not all-or-nothing: only the ones currently frightened read as the frightened
  // colour, in place of their own -- an already-respawned packmate keeps chasing in its own colour
  // while the rest are still fleeing blue. Built fresh each frame rather than cycled by glsim.js
  // (as COLORS.ghosts alone would be), since which index is which colour now depends on state.
  const ghostColors = ghosts.map((g, i) =>
    g.frightened ? COLORS.frightened : COLORS.ghosts[i % COLORS.ghosts.length]);
  SimGL.draw({ ...COLORS, ghosts: ghostColors });

  const r = cv.getBoundingClientRect();
  if (ov.width !== Math.round(r.width)) { ov.width = Math.round(r.width) || W; ov.height = Math.round(r.height) || H; }
  octx.clearRect(0, 0, ov.width, ov.height);
  const sc = ov.width / W;

  // Unlike the rest of this overlay, not gated on showOverlay: it's part of what a ghost *is* to
  // the player, not a CARL debugging aid, so it stays up with the chrome hidden the same way the
  // ghosts' own board-rendered colour does.
  drawGhostEyes(sc);

  // Everything below is CARL's working-out drawn over the board -- the two direction arrows and
  // the intervention discs. Hidden together with the control deck, leaving just the game.
  if (showOverlay) {
    if (lastCoM) {
      // One tile long: enough to read the heading against the maze, short enough not to cover it.
      const cx = lastCoM[1] * sc, cy = lastCoM[0] * sc, alen = (maze.cellX || 40) * sc;
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
  }

  $('r-mass').textContent = lastCoM ? lastCoM[2].toFixed(0) : '0';
  $('r-step').textContent = steps;
  $('r-acts').textContent = actions;
  $('r-course').textContent = courseChanges;
  $('r-ms').textContent = lastMs ? lastMs.toFixed(1) + 'ms' : '–';
  $('r-rate').textContent = measSps ? measSps.toFixed(0) + '/s' : '–';
  $('lives').textContent = '💛'.repeat(lives);
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
    // Steps the budget could not afford are dropped here rather than carried over -- rolling them
    // into the next frame would only make the following batch longer, and so on downwards. The
    // budget is checked after a step, not before, so a device where one step alone blows the whole
    // budget still advances by one instead of stalling forever.
    let done = 0;
    try {
      for (let k = 0; k < want; k++) {
        if (solitonDead || levelWon || ghostEatPause) break;   // holding for a jingle
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
// Every user-initiated steer goes through here, so a single place counts course changes.
// Re-issuing the direction the agent is already chasing (holding an arrow key down, clicking
// straight ahead) isn't a change of course, so it isn't counted.
function steerTo(dy, dx) {
  const n = Math.hypot(dy, dx); if (n < 1e-6) return;
  const ny = dy / n, nx = dx / n;
  if (ny * dir[0] + nx * dir[1] < 0.9999) courseChanges++;
  dir = [ny, nx];
  sometimesRemaining = sometimesWindow;    // a fresh steer wakes CARL up again in 'sometimes' mode
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
// Step always leaves the run paused and the board one step further on, so it reads as a scrub
// rather than a nudge to a still-running sim. If the loop happens to be mid-batch when the click
// lands (busy), that in-flight step is the advance -- stepping again here would double it.
async function stepOnce() {
  setRunning(false);
  if (busy || solitonDead) return;
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
$('maze-toggle').addEventListener('click', () => {
  mazeEnabled = !mazeEnabled;
  $('maze-toggle').textContent = mazeEnabled ? 'On' : 'Off';
  $('maze-toggle').classList.toggle('on', mazeEnabled);
  newMaze();       // respawns the soliton, which plays the start jingle and then starts the game
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
$('window').addEventListener('input', e => {
  sometimesWindow = +e.target.value;
  $('v-window').textContent = sometimesWindow;
});

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
