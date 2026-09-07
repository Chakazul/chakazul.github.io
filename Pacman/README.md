# CARL Maze Demo — WebGL

The same demo as `../CARL/maze_playground.html`, running the same trained policy, with the Lenia
simulation moved from JS onto the GPU. Serve the folder over http(s) and open `index.html` — it
fetches its shaders and its model, so `file://` will not work.

## Why

The CPU demo's per-step cost was dominated by its convolution. It used an FFT to avoid an
O(taps·cells) direct sum, but the board sizes offered in the UI (100/150/200/250) are not powers
of two, so every step fell onto the Bluestein path — three padded transforms per axis, in JS.
On a GPU the direct sum is the fast path instead: every cell is an independent fragment, so the
~950-tap kernel costs nothing that matters, and the FFT is not needed at all.

Locating the soliton and painting the board were the other two whole-board JS loops. Both are
now shader passes.

## What runs where

| | CPU demo | this one |
|---|---|---|
| convolution + growth | JS, FFT | `shaders/sim.glsl` |
| wall collision | JS sweep | folded into `sim.glsl` |
| intervention (add/remove mass) | JS | `shaders/action.glsl` |
| center of mass + total mass | JS sweep of the whole board | `shaders/reduce.glsl` → `shaders/com.glsl` |
| policy input window | JS crop of 4 stored boards | `shaders/crop.glsl` |
| board rendering | per-pixel JS + `putImageData` | `shaders/draw.glsl` |
| **policy network** | onnxruntime-web (WASM) | **unchanged** |
| maze generation, episode logic, overlay, UI | JS | unchanged (`app.js`) |

## The one synchronization point

The policy runs on the CPU, so something has to cross back from the GPU every step. Each step
queues five passes with no readback between them — action, sim, reduce, com, crop — and the CPU
then collects both results it needs in a single `readback()`:

- the **CoM** (a 1×1 texture: row, col, mass, valid) for the episode bookkeeping and the overlay
- the **96×96×4 crop** the policy reads

The crop pass reads the CoM out of a *texture* rather than a uniform, which is what makes one
stall enough: it does not have to wait for the CPU to be told where the soliton is before it can
crop around it.

Board state is a ring of four `R32F` textures rather than the usual ping-pong pair, because the
policy is fed a 4-frame stack and all four frames must be cropped at the *same* origin. Cropping
each frame at the CoM it had when captured would show a stationary soliton and destroy the
velocity signal the policy depends on. A step overwrites the four-steps-ago frame, which is
exactly the one falling out of the stack, so four slots is the whole requirement.

## Fidelity to the CPU version

The step order is deliberately identical to the CPU demo's `agentStep()` — act on the current
board, then step, then locate, then judge — because the policy is sensitive to it. Beyond that:

- The kernel is built by the same arithmetic as the CPU demo's `buildKernel()`, including the
  1e-7 tap threshold and normalization by the unthresholded sum, and uploaded as a texture. So
  the weights are the ones the CPU version used, not a shader re-derivation of them.
- The shader sums `state[p+d]`; the CPU's FFT computes the `state[p-d]` form. These agree because
  the kernel is radially symmetric (verified: max asymmetry ~1e-17). Direct sum vs. FFT circular
  convolution agree to ~4e-15.
- The CoM is the same toroidal circular mean, just summed in two stages; it reproduces the direct
  global computation to ~1e-14, and the 16×16 block grid covers every cell at all four board sizes.
- Board state is `R32F`, not the 8-bit textures the `Lenia/WebGL` demos use — the policy was
  trained on float32 states and 8-bit quantization would be a real change to its input.

## Differences on purpose

- The FPS slider goes to 480 (was 120), and a frame absorbs as many steps as fit a wall-clock
  budget (`?budget=`, 10ms) rather than a fixed count. A fixed count is wrong here because a step
  costs wildly different amounts in the two actor modes: pure GPU work while you are acting, plus
  a policy inference while CARL is. Sixteen of the latter ran for most of a second before the one
  repaint at the end of the batch, which read as the board freezing. The `rate` readout shows the
  rate actually achieved, which is where the real ceiling shows up.
- No inline fallback soliton. The CPU demo carried one for the `file://` case; this version
  fetches shaders too, so that case cannot arise.

## Known limits

- Needs WebGL 2 **and** `EXT_color_buffer_float` (rendering to float textures). Both are checked
  at startup and reported in the page if missing.
- The readback is synchronous, so the sim cannot run ahead of the policy. If that becomes the
  bottleneck, the next step is PBO-based async readback (`readPixels` into a `PIXEL_PACK_BUFFER`
  plus a fence), trading one frame of latency for no stall.
- Board edge is capped at 256 by the reduction's fixed 16×16 block loop; `setBoard` throws if
  that is exceeded rather than silently missing cells.
