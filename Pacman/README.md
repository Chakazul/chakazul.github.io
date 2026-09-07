# Lenia Pac-Man

A Pac-Man skin on `../CARL-WebGL`: same GPU engine, same trained policy, but the free-form maze
generator is replaced with a fixed, hand-authored Pac-Man-style level, and a second, free-running
Lenia channel plays the part of the dots. Serve the folder over http(s) and open `index.html` — it
fetches its shaders and its model, so `file://` will not work.

## What's Pac-Man about it

- **The maze.** `MAZE_LAYOUT` in `app.js` hand-authors a fixed 5×11-cell level in the standard
  maze-ASCII form (odd/odd = cell, even/even = corner post, everything else = the wall slot
  between two cells) instead of the randomized-DFS generator `CARL-WebGL` used. Openings on the
  outer ring aren't dead ends: the sim already wraps toroidally, so they become real Pac-Man-style
  side tunnels to the opposite edge. Wall corner posts are rounded into circles, but only on the
  side(s) that face open floor — a post next to a solid wall bar stays square so the bar's flat
  end meets it flush instead of notching a chip out of the corner. `M`/`O` in the layout (ghost
  house, power pellets) are placed but not yet interpreted — see Known limits.
- **The dots.** A second Lenia channel (`channel2RuleName`/`channel2R` in `CFG`), rendered in
  white, gets one small free-running soliton dropped on every `.` in the layout and is never
  steered by the policy. `shaders/eatsum.glsl` + `shaders/eatreduce.glsl` add a block-tiled
  reduction, alongside the CoM one, that measures how much dot mass sits under Pac-Man's channel
  above `EAT_THRESHOLD` each step; `sim.glsl` erases that mass as part of the same step. Whether
  anything was eaten comes back in the same readback as the CoM and crop, and drives the
  `eat_dot_0`/`eat_dot_1` "waka waka" loop (`noteEating()`) — it keeps alternating for as long as
  dots keep landing within `CFG.eatSoundGraceMs` of each other, so a good run reads as one
  continuous chomp rather than discrete blips.
- **Sound and lives.** A start jingle plays on spawn (deferred, if needed, to the page's first
  click/tap/key, per browser autoplay rules); dying — mass below `massDeathFraction` of spawn mass
  or above `massExplodeLimit`, same thresholds as `CARL-WebGL` — plays a death jingle, holds the
  board for `DEATH_PAUSE_MS`, then auto-respawns the same soliton with no confirmation, so the game
  never shows a "Game Over" screen. `?sound=0` disables all of it (jingles, chomp, and the death
  pause — a death just holds silently and respawns).
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
| convolution + growth (both channels) | JS, FFT | `shaders/sim.glsl` |
| wall collision | JS sweep | folded into `sim.glsl` |
| intervention (add/remove mass) | JS | `shaders/action.glsl` |
| center of mass + total mass | JS sweep of the whole board | `shaders/reduce.glsl` → `shaders/com.glsl` |
| eaten dot mass | n/a (no dots channel) | `shaders/eatreduce.glsl` → `shaders/eatsum.glsl` |
| policy input window | JS crop of 4 stored boards | `shaders/crop.glsl` |
| board rendering (walls, both channels) | per-pixel JS + `putImageData` | `shaders/draw.glsl` |
| **policy network** | onnxruntime-web (WASM) | **unchanged** |
| maze layout, episode logic, overlay, UI, sound | JS | unchanged shape, new content (`app.js`) |

## The one synchronization point

The policy runs on the CPU, so something has to cross back from the GPU every step. Each step
queues six passes with no readback between them — action, sim, reduce, com, eatreduce+eatsum,
crop — and the CPU then collects all three results it needs in a single `readback()`:

- the **CoM** (a 1×1 texture: row, col, mass, valid) for the episode bookkeeping and the overlay
- the **eaten-mass total** (a 1×1 texture), thresholded on the CPU into "something was eaten"
- the **96×96×4 crop** the policy reads

The crop pass reads the CoM out of a *texture* rather than a uniform, which is what makes one stall
enough: it does not have to wait for the CPU to be told where the soliton is before it can crop
around it.

Board state is a ring of four `R32F` textures rather than the usual ping-pong pair, because the
policy is fed a 4-frame stack and all four frames must be cropped at the *same* origin. Cropping
each frame at the CoM it had when captured would show a stationary soliton and destroy the velocity
signal the policy depends on. A step overwrites the four-steps-ago frame, which is exactly the one
falling out of the stack, so four slots is the whole requirement.

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
- The dots channel has no CPU-demo equivalent to match fidelity against — it's new to this fork,
  not a port.

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
- Fixed hand-authored maze instead of the randomized-DFS generator, and a second, unsteered Lenia
  channel for dots — both new to this fork, not present in either the CPU demo or `CARL-WebGL`.

## Known limits

- Needs WebGL 2 **and** `EXT_color_buffer_float` (rendering to float textures). Both are checked
  at startup and reported in the page if missing.
- The readback is synchronous, so the sim cannot run ahead of the policy. If that becomes the
  bottleneck, the next step is PBO-based async readback (`readPixels` into a `PIXEL_PACK_BUFFER`
  plus a fence), trading one frame of latency for no stall.
- Board edge is capped at 256 by the reduction's fixed 16×16 block loop; `setBoard` throws if that
  is exceeded rather than silently missing cells.
- `M` (ghost house) and `O` (power pellets) in `MAZE_LAYOUT` are placed but not yet interpreted —
  `layoutMaze()` currently treats them as plain open floor. There are no ghosts and no power-pellet
  effect yet.
- No lives counter or score — death is an infinite respawn loop, and eaten dots don't come back
  (until Restart/New Maze/soliton change resets the dots channel).
