# Lenia Pac-Man

A Pac-Man skin on `../CARL-WebGL`: same GPU engine, same trained policy, but the free-form maze
generator is replaced with a fixed, hand-authored Pac-Man-style level, and two further
free-running Lenia channels play the parts of the dots/pellets and the ghosts. Serve the folder
over http(s) and open `index.html` — it fetches its shaders and its model, so `file://` will not
work.

## What's Pac-Man about it

- **The maze.** `MAZE_LAYOUT` in `app.js` hand-authors a fixed 5×11-cell level in the standard
  maze-ASCII form (odd/odd = cell, even/even = corner post, everything else = the wall slot
  between two cells) instead of the randomized-DFS generator `CARL-WebGL` used. Openings on the
  outer ring aren't dead ends: the sim already wraps toroidally, so they become real Pac-Man-style
  side tunnels to the opposite edge. Wall corner posts are rounded into circles, but only on the
  side(s) that face open floor — a post next to a solid wall bar stays square so the bar's flat
  end meets it flush instead of notching a chip out of the corner. Two layout characters are pure
  authoring aliases, so a hand-edited row stays readable as ASCII art: `+` means exactly `.` (a
  dot) and `-` means exactly a space (open). `^` is open floor too, but a ghost standing on one is
  sent north rather than choosing for itself — that is what gives the centre room a one-way exit.
- **The dots and power pellets.** A second Lenia channel (`channel2RuleName`/`channel2R` in
  `CFG`), rendered in white, gets one small free-running soliton dropped on every `.`/`+` in the
  layout, and one on every `O` — the same soliton and rule, just recoloured (`setPowerMask()`),
  since the shared channel-2 rule fixes one equilibrium size regardless of stamp size. Neither is
  steered by the policy. `shaders/eatsum.glsl` + `shaders/eatreduce.glsl` add a block-tiled
  reduction, alongside the CoM one, that measures how much dot mass sits under Pac-Man's channel
  above `EAT_THRESHOLD` each step, split into total and power-pellet-only; `sim.glsl` erases that
  mass as part of the same step. Whether anything (or any pellet) was eaten comes back in the same
  readback as the CoM and crop, and drives the `eat_dot_0`/`eat_dot_1` "waka waka" loop
  (`noteEating()`) and the frightened window (below).
- **The ghosts.** A third Lenia channel (`channel3RuleName` in `CFG`), rendered per-ghost from
  `COLORS.ghosts`, gets one full-size free-running soliton per numbered spawn (`1`-`9`) in the
  layout, gated by the current `level` (level *N* spawns every ghost numbered at or below *N* —
  see Levels and lives, below). Like the dots it's never steered by the policy. Ordinarily it's
  the eating relationship run backwards: `sim.glsl`'s same one-fetch-and-compare erases *Pac-Man's*
  mass wherever a non-frightened ghost overlaps him, draining him below `massDeathFraction` and
  triggering the ordinary death path — no separate collision test. During a frightened window
  (below) that reverses per-ghost. Ghosts are re-placed on every spawn, including death respawns,
  so they always restart in their house. `?ghost=0` turns the whole channel off, which also
  removes the only thing that can kill Pac-Man other than his own dynamics.
- **Power pellets & frightened ghosts.** Eating an `O` starts (or, mid-window, extends) a shared
  `frightenedDurationMs` window: every not-already-frightened ghost turns to flee, slows to
  `GHOST_FRIGHTENED_SPEED` of normal, and can now be eaten by Pac-Man instead of the reverse — the
  same one-fetch erase, just pointed the other way, per ghost tile (`uFrightened[i]` in
  `ghostsim.glsl`). Any ghost that dies (dissolves or explodes, see Ghost death below) while
  frightened is credited as eaten and scores 200, doubling for every next one credited in the same
  window (200, 400, 800, ...) — there's no separate bite signal, but frightened dynamics are
  dominated by the erase-on-touch effect, so a frightened death is Pac-Man's doing in practice. It
  respawns in its house immediately, ending its own flight early while its still-frightened
  packmates keep counting down theirs (`updateFrightened()`). A death while not frightened is an
  ordinary ghost death (below): no score, same as always.
- **Ghost death.** A ghost dies the same two ways Pac-Man does — dissolved below
  `massDeathFraction` of its spawn mass, or grown past `massExplodeLimit` — judged on the mass that
  already comes back in its tile's CoM readback. The consequence differs: no jingle, no held board,
  no episode end, and no effect on the other ghosts. It is put back in its own house, alone, while
  the rest of the pack keeps running.

  Rewriting its tile is also the entire cleanup after an explosion: a tile has hard edges, so
  however far the mess spread it is still inside that one tile, and the board-space texture
  everything downstream reads is rebuilt from the tiles every step rather than accumulated, so
  there is nowhere else for debris to have got to. The explode limit also stays meaningful inside
  a tile, which is not true of every windowing scheme — 1000 is well under a 96-square tile's
  9216-cell capacity, so a runaway trips it long before it saturates. Spawn mass is summed off the
  stamp rather than waited for from a readback, so the check is armed from the step the ghost
  appears, and stays disarmed whenever the channel is empty — what stops "no ghost" reading as
  "dead ghost" and respawning on a loop.
- **Levels and lives.** Lives start at `STARTING_LIVES` (3, capped at `MAX_LIVES` 5), shown as 💛
  in the title row, and persist across death-respawns within a game — lost one per death (except
  a mass-explosion, which isn't something the player could have steered around) and reset only on
  a genuine new game (Restart, New Maze, soliton picker, or running out and choosing to play
  again). Clearing the board (all dots and pellets gone) awards a life back, advances `level` by
  one, and respawns everything — dots refilled, ghost count following the new level — after an
  intermission jingle, carrying the current life count forward rather than resetting it. `?level=N`
  picks a starting level other than the default 3 (also what a Restart falls back to). `?god=1`
  keeps the death beat (jingle, pause, respawn) but stops it from costing a life, so a game over
  never interrupts a run.
- **Pac-Man is windowed too.** His simulation is confined to a 128-square patch that follows his
  centre of mass, taking his channel from 188M tap-iterations per step to 22M. Unlike the ghosts he
  keeps an ordinary board-sized texture and only the *step* is scissored, so the policy's frame
  stack, the CoM reduction, the crop and both eat checks carry on reading a full board and needed
  no changes at all.

  128 rather than a ghost's 96, for a reason a ghost does not have: he takes interventions. CARL
  picks a cell anywhere in its 96×96 view — 48px off centre — and lays a 7px disc on it, so real
  mass can arrive 55px out; a 96 window (half-width 48) would silently clip anything past 41px.
  The other thing to check is that the crop still sees what it used to: each frame of the stack was
  windowed around the CoM it had at the time, and the crop is taken at the *current* CoM, so the
  oldest frame is read up to 49.9px from its own window centre against a half-width of 64 — 16px
  of headroom, and the soliton only spans ~23px anyway, so the outer region is empty in the
  full-board version too.

  Total per-step cost with three ghosts: 803M before any of this, 276M with the ghosts windowed,
  110M with both. The dots are now the most expensive channel on the board.
- **Several ghosts, one pass.** Each numbered spawn is simulated in a private 96×96 tile of one
  atlas texture (`shaders/ghostsim.glsl`) rather than in a shared board-sized field. Two ghosts in
  one field are two Lenia solitons, and two Lenia solitons that meet annihilate or blow up; tiles
  make them pass through each other the way arcade ghosts do. The convolution refuses to read past
  its own tile's edge, so the isolation is structural — verified by packing a neighbouring tile
  with mass right up against the shared boundary and confirming a ghost's step is bit-identical
  either way.

  A tile is a window onto the board, not a world of its own: it carries a board origin, so the wall
  lookup is a real maze lookup and ghosts still run the corridors. It slides by whole cells each
  step to keep its soliton centred — whole cells because a fractional slide means resampling, and
  resampling a soliton every step smears it away. Over 400 steps of travel the soliton stays within
  1.1px of the tile centre, against a 48px half-edge.

  Cost scales with the number of ghosts rather than the board: a tile is 9216 cells against the
  board's 137500, so four ghosts come to about a quarter of the single whole-board pass this
  replaces. `shaders/ghostblit.glsl` composites the tiles back to a board-sized texture — in one
  pass over the board rather than a draw per tile, since tiles can overlap once ghosts are close
  and separate draws would have the later one *replace* the earlier instead of adding to it. It
  records which ghost owns each pixel as it goes, which is what `draw.glsl` colours from, so the
  colours stay exact even where tiles overlap. `shaders/tilered.glsl` → `shaders/tilecom.glsl` then
  locate every ghost inside its own tile in two draws for the whole pack, landing in one texel
  each so the readback stays the single stall one ghost cost.
- **Ghost navigation.** A free-running soliton already travels in a straight line, so the only
  thing needing a decision is what happens at a corner or crossroad. Rather than steer the turn
  with the policy — a second inference on every step, and inference is what a step actually costs
  — `steerGhost()` rotates the ghost's own field a quarter turn about its centre of mass
  (`shaders/rotate.glsl`). Turn cells are marked by the layout characters `+ - O ^ 1-9`, which in
  this level land on exactly the 32 corner/junction cells and on no straight corridor. At one,
  left/straight/right are filtered to whichever are open and one is chosen — never a reversal,
  unless it is a dead end and reversing is all there is. A `^` cell overrides the choice entirely
  and sends the ghost north.
- **Chasing.** `CFG.chaseBias` (`?chase=N`, a whole percentage, default 60) is how often that
  choice is the opening that closes on Pac-Man (or, frightened, opens away from him) rather than a
  roll of the dice, scored by dot product against the toroidal vector to him so a side tunnel is
  judged on where it comes out. This is the difficulty dial, and biasing only the turns compounds,
  because every junction is another chance to correct. Share of time spent within one tile of
  Pac-Man, over 24 runs against a Pac-Man moving at the same speed (both on `rule74`), measured
  when the default was 40%:

  | `chase` | 0% | 20% | 40% | 70% | 100% |
  |---|---|---|---|---|---|
  | time within a tile | 4.2% | 5.5% | 6.8% | 15.0% | 65.6% |

  Note the cliff at the top: with the two solitons matched for speed, a ghost that always turns
  toward him never loses him again — it just shadows him, mean separation 37px. The useful range
  is roughly 0–70%.
- **Ghost speed.** `channel3RuleName` is picked for pace, since a ghost that cannot keep up is not
  a threat. Measured along a corridor of this maze: `rule74` runs 0.479 px/step against `rule73`'s
  0.194 — a wider gap than in open space (0.452 vs 0.270), because `rule73` wastes much more of
  its motion on the ~32° lean it travels at where `rule74` leans only ~11°. That it is also
  `defaultRuleName` is deliberate: it puts the ghost at exactly Pac-Man's pace on the default
  pick, so it closes by out-navigating him rather than by simply being faster.
- **Reading the heading.** The ghost's direction of travel is measured from how its CoM has
  actually moved since it entered the current cell, never inferred from the rotation it was
  stamped with. It has to be: which way a freshly stamped soliton travels depends on how its bank
  entry happens to be drawn, and that is not consistent — rule74's canonical heading is right,
  rule73's is left. Assuming one silently rotates every turn decision by a constant for any rule
  that disagrees. The same measurement is what aims the ghost out of its spawn house: it is
  stamped in whatever orientation the bank holds, then turned to face the exit once it has moved
  the few px needed to say which way it is pointing.
- **Turning on the tile centre.** The turn fires at the ghost's closest approach to the tile
  centre, found by watching its distance stop falling, rather than when that distance first drops
  under a threshold. It has to be the centre: the rotation pivots on the CoM, so turning anywhere
  else swings the soliton into the corner it is going around. A threshold cannot express this —
  it fires on *entry* to the band, a third of a tile early, and tightening it enough to land on the
  centre makes it small enough for a faster soliton to step clean over. Closest approach has no
  tolerance to tune, is speed-independent, and always fires. Measured over a 20,000-step walk with
  the soliton's real drift modelled: median 0.53px from the centre along the direction of travel
  (max 2.04px), 0 steps inside a wall, all 55 open cells visited.
- **Why rotate rather than re-stamp.** Stamping the canonical pattern at the new heading is the
  obvious way to turn a soliton, and it reads as a blink: whatever the running pattern had evolved
  into is discarded for a fresh copy from the bank. Rotating the texture keeps it. A 90° rotation
  about an *integer* pivot maps texel centres exactly onto texel centres, so the pass is a
  permutation — one `texelFetch` per cell, no filtering, no resampling blur. Verified on an evolved
  (40-step) pattern: mass in equals mass out to the bit for all three turn counts, the CoM shifts
  only ~0.5px (pivot rounding), and the result still travels at the full 0.45 px/step on the
  heading the turn count claims. Nothing else about the ghost is touched — not moved, not
  recentred, not reset — so a turn changes only which way it points.
- **Sound.** A start jingle plays on spawn (deferred, if needed, to the page's first click/tap/key,
  per browser autoplay rules); dying plays a death jingle, holds the board for `DEATH_PAUSE_MS`,
  then auto-respawns (no confirmation) while lives remain, or holds on a "Start" prompt once they
  run out — same paused, wait-for-a-gesture presentation as the very first page load. Eating a
  power pellet plays a fire-and-forget fright sound; eating a frightened ghost holds the board for
  its own jingle, same beat as a death or a win but with no extra pause after. Clearing a level
  plays an intermission jingle ahead of the next spawn. `?sound=0` disables all of it (jingles,
  chomp, fright sound, and every pause — the sim just holds silently and moves on).
- **Score**, shown centered in the title row: 10 per dot, 50 per power pellet, and 200 for a ghost
  eaten during a pellet's frightened window, doubling for every next ghost eaten in that same window
  (200, 400, 800, ...) before resetting on the next pellet. Dots and pellets are counted per dot:
  `shaders/dotsites.glsl` sums the channel-2 mass left in a small window around every dot's stamped
  position (one texel per dot, read back with everything else), and a dot scores once, the step its
  window drops below `CFG.dotGoneFraction` (40%) of what it held at stamping. That works because a
  dot is bistable — measured on a CPU mirror of `sim.glsl`, an untouched or merely grazed dot never
  falls below ~80% over its ~52-step breathing cycle, while one Pac-Man bites into dissolves all the
  way to 0 within a few steps. Two cheaper signals were tried and don't work. Summing the mass erased
  under Pac-Man (`rb.eaten`, `shaders/eatreduce.glsl`) misses the part of a bitten dot he never
  covered, which then dissolves on its own, so it caught only 25–70% of each dot and many dots never
  scored. The channel's total mass swings by ~15 dots' worth, because every dot breathes ±15% in
  phase with every other. `rb.eaten`/`rb.pelletEaten` still drive the chomp sound and the frightened
  window, which only need "was anything bitten this step". Persists across level wins and death
  respawns, resetting only on a genuine new game.
- **Bonus fruit.** Twice a level, a bonus dot appears on Pac-Man's own spawn tile the moment the
  level's dots-eaten fraction crosses each entry of `CFG.bonusFruitThresholds` (30%, then 70%),
  drawn as an ordinary dot with an emoji over it. Left uneaten for `CFG.bonusFruitTimeout` steps
  (500) it disappears instead — the level's other threshold, if not yet crossed, still arms in its
  own time. Both which emoji and how much it's worth are keyed by the current `level` (1-9):
  `🍊🍎🍒🍓🍉🍭🍄🍖💩` and 100/100/100/200/500/700/1000/2000/5000 respectively (`BONUS_FRUIT_EMOJIS` /
  `BONUS_FRUIT_POINTS` in `app.js`) — both of a level's fruits are the same, since the arcade
  original ties fruit to level rather than to which of the level's two it is. The threshold-armed
  state (`bonusFruitThresholdIdx`) resets only where `dotSites` itself does, in `placeDots()` — a
  death respawn (which leaves already-eaten dots eaten) doesn't re-arm a threshold the level
  already passed, while a fresh level's dots get both back. Unlike an ordinary dot the fruit isn't
  part of the channel-2 field at all — it's a CPU-tracked position (`bonusFruit` in `app.js`)
  judged each step against `lastCoM`, since a bonus item is meant to vanish the instant Pac-Man
  touches it rather than dissolve over several steps the way a real dot does. Eating it plays
  `eat_fruit.wav` and holds the board for its own length, same beat as an eaten ghost. Both an
  eaten ghost and an eaten fruit flash a floating score-number popup slightly below where they died
  or were picked up (`ghostEatPopup`/`fruitEatPopup`), live for exactly as long as the board holds
  for the jingle — set alongside `ghostEatPause`/`fruitEatPause` and cleared the same moment.
- **Three actor modes**, cycled by the one button: **CARL acts sometimes** (default) only queries
  the policy for a configurable window of steps after the episode starts or after you last steer,
  then goes idle — same per-step cost as **You act** the rest of the time, which matters because
  policy inference, not the GPU work, is what a step actually costs (see Differences on purpose,
  below). **CARL acts always** queries it every step, as `CARL-WebGL` always did. **You act** turns
  the policy off entirely; left/right-click add or remove mass yourself.
- **No board-size picker.** `CARL-WebGL`'s 100–250 size choices are gone — the maze layout fixes
  the cell count, and the fixed 550×250 board is sized so wall thickness and cell width match what
  the layout needs (see the sizing comment above `MAZE_LAYOUT`). Undersized cells would clip the
  soliton on every turn.

## Why (inherited from CARL-WebGL)

The CPU demo's (`../CARL/maze_playground.html`) per-step cost was dominated by its convolution. It
used an FFT to avoid an O(taps·cells) direct sum, but the board sizes it offered are not powers of
two, so every step fell onto the Bluestein path — three padded transforms per axis, in JS. On a GPU
the direct sum is the fast path instead: every cell is an independent fragment, so the ~950-tap
kernel costs nothing that matters, and the FFT is not needed at all.

Locating the soliton and painting the board were the other two whole-board JS loops. Both are now
shader passes.

## What runs where

| | CPU demo | this one |
|---|---|---|
| convolution + growth (Pac-Man, windowed; dots) | JS, FFT | `shaders/sim.glsl` |
| wall collision | JS sweep | folded into `sim.glsl` |
| intervention (add/remove mass) | JS | `shaders/action.glsl` |
| center of mass + total mass | JS sweep of the whole board | `shaders/reduce.glsl` → `shaders/com.glsl` |
| eaten dot / pellet mass | n/a (no dots channel) | `shaders/eatreduce.glsl` → `shaders/eatsum.glsl` |
| mass left per dot (scoring) | n/a (no dots channel) | `shaders/dotsites.glsl` |
| ghost center of mass (per ghost) | n/a (no ghost channel) | `shaders/tilered.glsl` → `shaders/tilecom.glsl` |
| ghost 90° turn | n/a (no ghost channel) | `shaders/rotate.glsl` |
| ghost tiles → board space | n/a (no ghost channel) | `shaders/ghostblit.glsl` |
| policy input window | JS crop of 4 stored boards | `shaders/crop.glsl` |
| board rendering (walls, all three channels) | per-pixel JS + `putImageData` | `shaders/draw.glsl` |
| **policy network** | onnxruntime-web (WASM) | onnxruntime-web (WebGPU by default, falling back to WASM per-op; `?ep=wasm` forces CPU-only) |
| maze layout, episode logic, overlay, UI, sound | JS | unchanged shape, new content (`app.js`) |

## The one synchronization point

The policy runs on the CPU, so something has to cross back from the GPU every step. Each step
queues its passes with no readback between them — action, sim, reduce, com, eatreduce+eatsum,
dotsites, crop, plus the two free-running channels' own sim passes and the ghost reduction — and the
CPU then collects everything it needs in a single `readback()`:

- the **CoM** (a 1×1 texture: row, col, mass, valid) for the episode bookkeeping and the overlay
- the **eaten-mass totals** (dot and power-pellet, a 1×1 texture), thresholded on the CPU into
  "something/a pellet was eaten"
- the **mass left per dot** (one row, one texel per dot) for the score
- the **per-ghost CoM** (one row, one texel per ghost) for the turn logic — free, in that the
  stall has already happened by the time it is read
- the **96×96×4 crop** the policy reads

The crop pass reads the CoM out of a *texture* rather than a uniform, which is what makes one stall
enough: it does not have to wait for the CPU to be told where the soliton is before it can crop
around it.

Board state is a ring of four `R32F` textures rather than the usual ping-pong pair, because the
policy is fed a 4-frame stack and all four frames must be cropped at the *same* origin. Cropping
each frame at the CoM it had when captured would show a stationary soliton and destroy the velocity
signal the policy depends on. A step overwrites the four-steps-ago frame, which is exactly the one
falling out of the stack, so four slots is the whole requirement.

The policy's own compute runs on onnxruntime-web's WebGPU execution provider by default (measured
faster than WASM on the mobile devices this was tuned on; `?ep=wasm` forces CPU-only for
re-comparing on a new device). This does not remove the synchronization point above: the sim runs
in a WebGL2 context, and WebGL2 and WebGPU share no GPU memory, so the crop still has to cross
through the CPU to reach the net either way. WebGPU only speeds up the net's own conv work.

## Fidelity to the CPU version

The step order is deliberately identical to the CPU demo's `agentStep()` — act on the current
board, then step, then locate, then judge — because the policy is sensitive to it. Beyond that:

- The kernel is built by the same arithmetic as the CPU demo's `buildKernel()`, including the
  1e-7 tap threshold and normalization by the unthresholded sum, and uploaded as a texture. So the
  weights are the ones the CPU version used, not a shader re-derivation of them.
- The shader sums `state[p+d]`; the CPU's FFT computes the `state[p-d]` form. These agree because
  the kernel is radially symmetric (verified: max asymmetry ~1e-17). Direct sum vs. FFT circular
  convolution agree to ~4e-15.
- The CoM is the same toroidal circular mean, just summed in two stages; it reproduces the direct
  global computation to ~1e-14, and the 16×16 block grid covers every cell at the board size used
  here.
- Board state is `R32F`, not the 8-bit textures the `Lenia/WebGL` demos use — the policy was
  trained on float32 states and 8-bit quantization would be a real change to its input.
- The dots and ghost channels, the power-pellet/frightened mechanic, and the level/lives system
  have no CPU-demo equivalent to match fidelity against — they're new to this fork, not a port.

## Differences on purpose

- The FPS slider goes to 480 (was 120), and a frame absorbs as many steps as fit a wall-clock
  budget (`?budget=`, 10ms) rather than a fixed count. A fixed count is wrong here because a step
  costs wildly different amounts depending on who's acting: pure GPU work while you are, plus a
  policy inference while CARL is. Sixteen of the latter ran for most of a second before the one
  repaint at the end of the batch, which read as the board freezing. The `rate` readout shows the
  rate actually achieved, which is where the real ceiling shows up.
- A third actor mode, "CARL acts sometimes" (the default), only pays for inference for a
  configurable window of steps after the episode starts or after you last steer, then goes idle —
  same cheap per-step cost as "You act" the rest of the time. This is the mode that benefits most
  from the wall-clock budget above, since its per-step cost keeps switching between the two
  extremes mid-run.
- No inline fallback soliton. The CPU demo carried one for the `file://` case; this version fetches
  shaders too, so that case cannot arise.
- Fixed hand-authored maze instead of the randomized-DFS generator, and the dots/pellets, ghosts,
  frightened window, scoring, and level/lives system — all new to this fork, not present in either
  the CPU demo or `CARL-WebGL`.

## Known limits

- Needs WebGL 2 **and** `EXT_color_buffer_float` (rendering to float textures). Both are checked
  at startup and reported in the page if missing.
- The readback is synchronous, so the sim cannot run ahead of the policy. If that becomes the
  bottleneck, the next step is PBO-based async readback (`readPixels` into a `PIXEL_PACK_BUFFER`
  plus a fence), trading one frame of latency for no stall.
- Board edge is capped at 256 by the reduction's fixed 16×16 block loop; `setBoard` throws if that
  is exceeded rather than silently missing cells.
- Up to `MAX_GHOSTS` (9) ghosts are supported, matching the layout's numbered spawns (`1`-`9`); a
  higher count would need a wider ghost atlas and the matching `#define` bumped in three shaders.
- Eaten dots and pellets don't come back until the dots channel is reset — Restart, New Maze,
  soliton change, or a level win.
