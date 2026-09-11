# Performance notes

The game runs smoothly on a reasonably fast desktop, but stutters on weaker laptops and on
mobile — worst when CARL is acting. This is where the per-step budget actually goes, and what
can be done about it, cheapest first.

Nothing here has been implemented yet; it is a work list, not a changelog.

## Where the time goes

Two separate bottlenecks. On a weak device they compound, because a stalled readback and an
overloaded GPU each make the other's queue longer.

### 1. GPU: the convolution tap loops, and essentially nothing else

`sim.glsl` and `ghostsim.glsl` both run a `(2R+1)^2` tap loop per cell, and each iteration does
*two* `texelFetch`es — one for the kernel weight, then (if it is nonzero) one for the state. Every
other pass on the board — the reductions, the eat check, the blit, the draw — is 1–5 fetches per
cell and rounds to nothing beside this.

Counted for the fixed 550×250 board:

| pass | fragments | grid | fetches/step |
|---|---|---|---|
| Pac-Man, 128 window, R=18 | 16,384 | 37² | 38.4M |
| dots, full board, R=9 | 137,500 | 19² | **82.6M** |
| ghosts ×3, 96 tile, R=18 | 27,648 | 37² | 64.7M |
| **total, level 3 (the starting level)** | | | **185.7M** |
| **total, level 9 (all ghosts)** | | | **315.2M** |

At the default 60 steps/s that is **11.1 Gfetch/s** at level 3 and **18.9 Gfetch/s** at level 9. A
mid-range phone GPU delivers single-digit Gtexel/s in practice, so the default speed is 2–5× over
budget on mobile before the level system has even started adding ghosts — and `MAX_GHOSTS` is 9,
so a good player's reward for clearing boards is a 1.7× heavier step.

Two numbers worth carrying forward:

- **The dots are 44% of the total.** They are the largest single consumer on the board, and they
  are unsteered stationary decoration that CARL never sees.
- **Half of every tap loop is the kernel fetch.** The `if (w == 0.0) continue;` guard skips the
  *state* fetch for a zero tap, but the kernel fetch that discovered it was zero has already
  happened. At R=18, 972 of 1369 grid positions are nonzero; the other 397 cost a fetch each and
  buy nothing.

### 2. CPU: one synchronous 144 KB readback per step, five `readPixels` calls deep

`readback()` in `glsim.js` issues five separate `readPixels` calls per step — CoM (1×1), eaten mass
(1×1), ghost CoM (`MAX_GHOSTS`×1), dots CoM (1×1), and the crop. The crop is 96×96 RGBA32F =
**144 KB, read every single step, whether or not anything is going to use it**. In "You act" mode
and during the idle stretches of "CARL acts sometimes" it is read and thrown away.

On a tile-based mobile GPU a synchronous readback forces a tile resolve and a full pipeline flush.
This is very likely the largest mobile cost after the convolution itself.

## The work list

### Tier 1 — small diffs, behaviour bit-identical

1. **Read the crop only when something is about to infer.**
   `agentAct()` consumes `lastCrop` from the *previous* step's readback, so moving the crop
   `readPixels` out of `readback()` and into the top of `agentAct()` reads exactly the same pixels
   from exactly the same pass — the crop FBO still holds them, and nothing overwrites it in
   between. It removes 144 KB and a pipeline flush from every step where CARL is idle, which is
   most steps in the default actor mode and all of them in "You act". Biggest win per line changed.

2. **Skip the zero taps with per-row extents.**
   The kernel is a disc, so its nonzero taps in each row are a *contiguous* run — checked: the sum
   of per-row `[lo,hi]` extents is 973 against 972 actual nonzero taps at R=18, and 241 against 240
   at R=9, so bounding the inner loop by the extents skips exactly one zero tap more than the
   `w == 0.0` guard does and never skips a nonzero one. Pass `lo`/`hi` as a `uniform int[2R+1]` and
   loop over the run instead of the full row. Exact, not an approximation. **−18% fetches.**

3. **Bake the kernel weights into the shader instead of fetching them.**
   Generate `const float K[973] = float[](...)` into the shader source at link time and index it by
   the loop counter, which removes the kernel `texelFetch` entirely. The weights are still the ones
   `kernelData()` produces, so this stays bit-identical to the CPU demo the way the README's
   fidelity section requires. The kernel only changes when the rule changes (the soliton picker,
   the ghost/dots rule), so a shader recompile there is affordable.

   | | current | +row extents (#2) | +weights in shader (#3) |
   |---|---|---|---|
   | Pac-Man | 38.4M | 31.9M | 15.9M |
   | dots | 82.6M | 66.3M | 33.0M |
   | ghosts ×3 | 64.7M | 53.8M | 26.9M |
   | **total, level 3** | **185.7M** | **152.0M** | **75.8M** |
   | **total, level 9** | **315.2M** | **259.6M** | **129.5M** |

   **−59% at level 3.** Benchmark this one on a real phone before committing to it: some
   Mali/Adreno drivers materialize a ~1000-entry constant array badly, and if this is one of those
   cases, #2 alone still stands on its own.

4. **Collapse the four 1×1 readbacks into one.**
   A tiny gather pass that writes CoM, eaten mass, ghost CoM and dots CoM into a single
   12×1 RGBA32F texture, read with one `readPixels`. The data is trivial either way; what this
   removes is three driver round-trips per step, which are not free on mobile.

5. **Ship `ort.wasm.min.js` rather than `ort.webgpu.min.js`.**
   `index.html` loads the 333 KB WebGPU bundle and `loadModel()` then forces
   `executionProviders: ['wasm']`, so the WebGPU half is downloaded and parsed for nothing on every
   visit. (The README already explains why the WebGPU EP would not help: it does not share GPU
   memory with the WebGL2 context, so the crop would still travel through the CPU.)

### Tier 2 — a visible tradeoff, but the dots are 44% of the GPU cost

6. **Stop stepping the dots every step.**
   They are stationary — `CFG`'s own comment says the dots channel "never moves", which is why
   `sim.glsl` skips wall masking for it. Either step the channel every N steps (cost ÷N, they just
   breathe more slowly) or do not step it at all (stamp once, keep only the erase). Eating has to
   stay instant either way, so the off-steps need a cheap erase-only pass — ~3 fetches per cell
   against the 361 a full step costs, i.e. free. Suggest a `?dotsim=N` param, 1 on desktop and 3–4
   on mobile.

   A per-dot tile atlas, the way the ghosts already work, was costed and rejected: the layout has
   104 dots, so even 24² tiles come to 21.6M against the current 49.6M grid iterations — 2.3× for a
   large refactor, against ÷N for a few lines.

7. **Cap the ghost count on weak devices.**
   Level 9 is +130M fetches/step over level 3. Clamping spawns against measured `sps` keeps late
   levels playable instead of turning the player's reward into a slideshow.

8. **Auto-scale quality rather than picking one default.**
   The honest fix for "smooth here, laggy there" is a controller that watches achieved `measSps`
   against the requested `sps` and walks down dots-sim rate → target `sps` → ghost cap until they
   meet. The speed slider's tooltip currently asks the player to do this by hand ("try 15–30
   steps/sec"), which is the right diagnosis and the wrong owner.

9. **Thread count on mobile.**
   `loadModel()` sets `numThreads = min(hardwareConcurrency || 6, 6)`. Phones commonly report 8
   cores of which four are tiny, and oversubscribing them is slower than 2–4 threads. Also worth
   confirming `crossOriginIsolated` is actually true on the target phones: if `coi-serviceworker`
   has not taken over yet, the fallback is `numThreads = 1` and inference runs 4–6× slower than
   the desktop numbers would suggest.

### Tier 3 — the big one, only if Tier 1 + 2 fall short

10. **Pack 4 cells per texel (`RGBA32F` state, W/4 wide).**
    One state fetch feeds four output lanes through a 4×4 weight matrix per (dy, dx-block): per
    output texel, 1 state fetch + 4 weight fetches produce 16 MACs, against 32 fetches for the same
    16 MACs today. That is ~4.6× fewer fetches for ~1.7× more ALU, so call it 2–2.5× on a
    fetch-bound mobile GPU. The catch is that every pass reading the state texture has to change
    with it — `crop`, `reduce`, `action`, `eatreduce`, `draw`, `rotate`, `ghostblit` — so it is a
    rewrite, not an edit.

11. **PBO async readback.**
    Already noted in the README's Known limits: `readPixels` into a `PIXEL_PACK_BUFFER` plus a
    `fenceSync`, trading one step of policy latency for no stall. Worth doing *after* #1, which
    removes most of the stalls this would otherwise be fixing.

## Considered and rejected

- **`textureGather`.** Would fetch 2×2 texels per instruction, a clean 4×. Not available: WebGL 2
  is GLSL ES 3.00, and `textureGather` is ES 3.10 and up.
- **A separable (sum-of-Gaussians) approximation of the kernel.** Two 1D passes would be ~13×
  cheaper per term, but it changes the weights, and the README's fidelity section is explicit that
  the shader multiplies by the exact weights the CPU version used because the policy is sensitive
  to them. It also risks destabilizing the solitons themselves. Defensible for the free-running
  dots and ghost channels alone, where no policy is watching — but not worth the split.

## How the numbers were derived

Fragment counts are the scissor rectangles the passes actually use: 128² for Pac-Man's window
(`PAC_WIN`), 96²×`gCount` for the ghost atlas (`GHOST_WIN`), and the full 550×250 board for the
dots. Grid sizes are `(2R+1)²` for `CFG.R = 18` (Pac-Man, ghosts) and `CFG.channel2R = 9` (dots).
Nonzero and per-row-extent tap counts come from running `glsim.js`'s own `kernelData()` — same
`quad4()` arithmetic, same 1e-7 threshold — and counting. "Fetches" is
`fragments × (grid + nonzero)`: one kernel fetch per grid position, plus one state fetch per
nonzero tap.
