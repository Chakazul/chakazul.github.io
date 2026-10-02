# Performance notes

The game runs smoothly on a reasonably fast desktop, but stutters on weaker laptops and on
mobile — worst when CARL is acting. This is where the per-step budget actually goes, and what
can be done about it, cheapest first.

Started as a work list. Items marked **(done)** have since been implemented; the rest is still
a list, not a changelog. Measured device numbers are in "Device measurements" below.

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

1. **Read the crop only when something is about to infer.** **(done: `readCrop()`)**
   `agentAct()` consumes `lastCrop` from the *previous* step's readback, so moving the crop
   `readPixels` out of `readback()` and into the top of `agentAct()` reads exactly the same pixels
   from exactly the same pass — the crop FBO still holds them, and nothing overwrites it in
   between. It removes 144 KB and a pipeline flush from every step where CARL is idle, which is
   most steps in the default actor mode and all of them in "You act". Biggest win per line changed.

2. **Skip the zero taps with per-row extents.** **(done, together with #3)**
   The kernel is a disc, so its nonzero taps in each row are a *contiguous* run — checked: the sum
   of per-row `[lo,hi]` extents is 973 against 972 actual nonzero taps at R=18, and 241 against 240
   at R=9, so bounding the inner loop by the extents skips exactly one zero tap more than the
   `w == 0.0` guard does and never skips a nonzero one. Pass `lo`/`hi` as a `uniform int[2R+1]` and
   loop over the run instead of the full row. Exact, not an approximation. **−18% fetches.**

3. **Bake the kernel weights into the shader instead of fetching them.** **(done, as a uniform
   buffer rather than a const array: see the note at the end of this item)**
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

   *As implemented:* weights go in a uniform buffer object, not a generated `const` array. That
   sidesteps the const-array risk above, and UBO reads in the same order across fragments are
   the usual fast path for this, but it is still unmeasured on a phone. Each row's nonzero run is
   zero-padded to whole vec4 chunks, so the shader reads four weights at a time with constant
   component indices. The padding costs a little: up to 1068 state fetches per cell at R=18 (vs.
   973 unpadded) and 292 at R=9 (vs. 241), so the saving is −54% / −51% per cell rather than
   the table's −59%. With #6 at its phone default (`dotsim` 3), level 3 comes to ~61M fetches per
   step against 185.7M (−67%), and ~87M (−53%) on desktop, where the dots still step every
   step. Fetches are a proxy, not a timing: the device re-run is what counts. The sums are bit-identical to
   the texture path (checked in float32 over every bank rule, at wrap and tile edges), and
   `?kernel=tex` keeps the texture path for an on-device A/B.

4. **Collapse the four 1×1 readbacks into one.** **(done: `gather.glsl`, the per-dot row included)**
   A tiny gather pass that writes CoM, eaten mass, ghost CoM and dots CoM into a single
   12×1 RGBA32F texture, read with one `readPixels`. The data is trivial either way; what this
   removes is three driver round-trips per step, which are not free on mobile.

5. **Ship `ort.wasm.min.js` rather than `ort.webgpu.min.js`.**
   `index.html` loads the 333 KB WebGPU bundle and `loadModel()` then forces
   `executionProviders: ['wasm']`, so the WebGPU half is downloaded and parsed for nothing on every
   visit. (The README already explains why the WebGPU EP would not help: it does not share GPU
   memory with the WebGL2 context, so the crop would still travel through the CPU.)

### Tier 2 — a visible tradeoff, but the dots are 44% of the GPU cost

6. **Stop stepping the dots every step.** **(done: `?dotsim=N`, default 3 on phones/tablets, 1 elsewhere)**
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

9. **Thread count on mobile.** **(done: capped at 4, from the device measurements)**
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
    with it — `crop`, `reduce`, `action`, `eatreduce`, `dotsites`, `draw`, `rotate`, `ghostblit` — so it is a
    rewrite, not an edit.

11. **PBO async readback.**
    Already noted in the README's Known limits: `readPixels` into a `PIXEL_PACK_BUFFER` plus a
    `fenceSync`, trading one step of policy latency for no stall. Worth doing *after* #1, which
    removes most of the stalls this would otherwise be fixing.

## CARL inference

Everything above is about the GPU sim. Policy inference is a separate cost, and on slow devices
it is the one players actually notice: a step with CARL acting is a step without CARL plus one
`session.run()`. "CARL acts sometimes" limits how *often* that cost is paid; the items below are
about making each payment cheaper. Lowering the inference rate with `?stride=N` was tried first,
and it makes pattern maintenance visibly worse.

### What one inference costs

`models/agent_direction.onnx` is a 4-level U-Net with FiLM conditioning: the 4-float context
(time, direction, cost) goes through a small MLP that scales and shifts each conv layer's output.
It has 0.79M parameters and **0.377 GMAC per call at 96×96**, in 694 ONNX nodes, of which only
15 are convolutions. The other nodes are circular padding built from Slice/Concat (124 nodes),
14 separate FiLM MLPs (56 `Gemm`), and shape bookkeeping (`Shape`/`Gather`/`Unsqueeze`/...).

| layer group | MMAC | share |
|---|---|---|
| encoders (96² → 48² → 24²) | 90 | 24% |
| bottleneck (12²) | 32 | 8% |
| **decoders** (skips are concatenated, so each block's first conv sees 3× the channels) | **255** | **68%** |

Native onnxruntime on CPU, single thread, on the desktop the game runs smoothly on:

| variant | ms / inference |
|---|---|
| 96, basic graph optimization | 9.0 |
| 96, full optimization | 7.2 |
| 96, full optimization + static input shape | 7.2 |
| 96, dynamic int8 quantization | 6.9 |
| 64 crop | 3.3 |
| 48 crop | 2.0 |

WASM is typically 2–3× slower than native, and a phone CPU is several times slower again, so an
estimate of 50–100+ ms per inference on mobile is plausible. That is many times a whole sim step.
The cost is in arithmetic, not per-node overhead, so on CPU it only shrinks by doing fewer MACs.

### Closed-loop check of the cheap ideas

A numpy copy of the game's step order (act → Lenia step at `dt × channel1Speed` → toroidal CoM →
4-frame crop at the current CoM), driving the real ONNX policy. 160² torus, **no maze walls**, 2
solitons (`rule74`, `rule73`) × 4 directions × 500 steps per config. It is not validated against
the game, so read it as a sanity check, not a measurement:

| config | deaths / 8 |
|---|---|
| 96 crop, infer every step (current) | 0 |
| `?net=64` | 2 |
| `?net=48` | 3 |
| stride 2, no action on skipped steps (current `?stride=2`) | 0 |
| stride 4, no action on skipped steps | 0 |
| stride 2, **repeat** last action on skipped steps | 0 |
| stride 4, **repeat** last action on skipped steps | **8** (all within ~100 steps) |

Steering progress was too noisy across 8 runs to rank configs, and without walls the stride
rows can't show the maintenance loss seen in the game. One pattern was clear: at 96, every
step, the policy chose **no-op on ~90–95% of steps**. Most inference only confirms "do nothing",
which is why skipping steps is tempting, and why skipping the few that matter hurts.

### In this repo, no retraining

1. **Measure on the target devices first.** Mobile has no console, so put a diagnostics line on
   the page: which execution provider actually loaded (WebGPU, or the silent WASM fallback),
   `crossOriginIsolated`, `numThreads`, and per-step inference ms against sim+readback ms.
   Without cross-origin isolation, WASM drops to 1 thread, which alone is 4–6× slower.
   Check that `coi-serviceworker` really takes control on iOS Safari.

2. **Pipeline the sim with inference (one step of action latency).** Today a step runs one
   thing after another:

   ```
   [ infer on crop(t) ][ GPU: act + sim + reduce + crop  →  readback ]  → step t+1
   ```

   Per-step time is `infer + sim + readback`. Pipelined, as soon as the readback returns crop(t),
   start inference on it *without awaiting it*. Then immediately queue the next GPU step using
   the action from the inference that just finished, the one started on crop(t−1):

   ```
   CPU/worker:  [ infer crop(t-1) ][ infer crop(t)   ][ infer crop(t+1) ]
   GPU:              [ step t+1 + readback ][ step t+2 + readback ] ...
   ```

   Per-step time becomes `max(infer, sim + readback)`. This is not a stride: CARL still decides
   every step, but each decision is based on a state one step old. Notes:
   - It only overlaps if inference really runs off the main thread. That holds for WebGPU, or
     for WASM with `ort.env.wasm.proxy = true`, which runs it in a worker. Plain WASM computes on
     the main thread and serializes again.
   - The gain is largest when the two halves are similar in cost (up to 2×). If inference is
     80 ms against a 10 ms sim step, it saves ~11%. The other benefit is that the page keeps
     drawing frames while inference runs.
   - The one-step delay is a distribution shift: the policy was trained to act on the state it
     sees. The soliton moves under 0.5 px/step against a 7 px action disc, so it is probably
     tolerable, but it has to be checked in play. Retraining item 4 below removes the doubt.
   - Not to be confused with PBO async readback (Tier 3 #11), which hides the readback stall but
     not inference.

3. **Simplify the graph offline.** Ship a pre-optimized model made by a one-off Python script.
   The context changes only on a steer or a cost-slider move, so the 14 FiLM MLPs can run once
   per change, with their γ/β outputs fed to the conv graph as inputs. Even better, fold γ/β
   into each conv's weights and bias. Replacing the Slice/Concat circular padding with `Conv`
   zero padding removes ~124 nodes, each a full feature-map copy. That changes the network's
   math at the crop edge, where the crop is usually empty, so check action agreement against
   the original. Expected gain on CPU is 0–20%. On WebGPU every node is a separate dispatch, so
   the gain could be much larger, but that is unmeasured.

   **(Done, partly: `tools/split_model.py`.)** What was built:
   - **FiLM split out** into `agent_direction_film.onnx` (context → 28 γ/β tensors, 112 ops),
     run on WASM once per context change and cached (`filmFeeds()`). The core takes γ/β as
     inputs, which removes the MLPs, the reshapes and all the shape arithmetic.
   - **Upsamples by a constant scale of 2** instead of a size computed from the skip's shape.
     Identical at every level of this net, and it keeps the core fully convolutional.
   - Two cores: `agent_direction_core.onnx` keeps the Slice/Concat circular padding (150 ops,
     any crop size). `agent_direction_core_gather.onnx` does each circular pad as two `Gather`s
     with constant indices `[n-1, 0..n-1, 0]` (94 ops, fixed at 96×96).
   - Both reproduce the original's Q map **bit for bit**: max |Δq| = 0, same action on 400/400
     crops from closed-loop rollouts (2 solitons × 4 directions, random action costs). The
     Slice core also matches at 64, 80 and 128.
   - Native CPU, 1 thread: original 7.1 ms, Slice core 6.8 ms, **Gather core 9.9 ms**: `Gather`
     is slow on CPU (37% of that model's time). Hence the default: Gather core on WebGPU at
     96, Slice core otherwise. `?model=orig|slice|gather` overrides, and the debug Bench times
     all three on each backend.

   Rejected along the way:
   - **Zero padding.** Same action on only 22/400 crops: the net depends on the wraparound.
   - **`Pad` with `mode='wrap'`.** One op per conv instead of 6, and onnxruntime-web 1.20.1 has
     a WebGPU kernel for it. But that kernel's shader source has a stray `]`
     (`k += i32(uniforms.x_shape[…]]);`), so it would very likely fail to compile. Worth
     revisiting on a newer onnxruntime-web.
   - Folding γ/β into conv weights: the weights would become per-context inputs, re-uploaded
     every inference (0.79M floats).

4. **fp16 on WebGPU, static int8 on WASM.** Dynamic int8 bought nothing (7.2 → 6.9 ms). Static
   QDQ int8, calibrated on real crops, is where WASM speedups usually come from (~2×). Either
   precision needs an argmax-agreement check against fp32.

5. **Run the U-Net in WebGL2 shaders (the big one).** The crop is already a texture in the sim's
   context. 15 conv passes, 3 max-pools, 3 upsamples, FiLM folded into the conv weights per
   context, then an argmax pass, with only the chosen action read back. That removes the 144 KB
   crop readback, the CPU hop between the WebGL2 and WebGPU contexts, and the dependence on
   threads or WebGPU support. 0.38 GMAC of fragment-shader work competes with the sim for the
   GPU, so do Tier 1 first. It needs the same fidelity check against onnxruntime that the sim
   had against the CPU demo.

### Needs retraining or distillation (outside this repo)

1. **Distill into a smaller student.** Train it to reproduce the current model's Q maps and
   argmax on the current model's own rollouts. This is supervised learning, no RL loop. Halving
   every channel width gives ~4× fewer MACs. Adding skip connections instead of concatenating
   them shrinks the decoders, which are 68% of the cost. Best expected payoff for the effort.

2. **Coarser action map.** Output Q at 48×48 instead of 96×96. The action is a radius-7 disc,
   so 1 px placement precision buys little. This drops the last decoder level (85 MMAC, 23%) and
   shrinks the output 4×.

3. **Train for a lower decision rate, rather than skipping decisions.** The current model was
   only ever trained to act every step, so `?stride=N` just leaves it blind in between, and
   repeating its action is worse (above). Train with one decision per k Lenia steps, and add k
   to the context (as the cost already is) so a single model covers k = 1–4 and the game can pick
   k per device. Optionally allow a stronger or multi-site intervention per decision, so a
   low-rate agent keeps the same control authority.

4. **Train with one step of action latency,** so in-repo item 2 is in-distribution rather than
   an approximation.

5. **Train on a smaller window with zero padding.** The game board is not a 96² torus anyway,
   and cost scales with crop area: 64² is 2.2× cheaper, 48² is 3.7×. Running the existing model
   at those sizes (`?net=`) killed solitons above, so a retrain is the only way to get there.

6. **A tiny "act now?" gate.** With ~90%+ of decisions being no-op, a small classifier trained
   on the full policy's no-op vs act labels could decide when the full net needs to run,
   replacing the fixed `sometimes` window and `stride` with a learned trigger.

Suggested order: diagnostics on a real phone (#1), then pipelining (#2) if the measured split
says it helps. Externally: distillation plus a coarser action head, then decision rate k as a
context input. Only port the network to WebGL (#5) if a distilled model is still too slow.

## Measuring on a device

`?debug=1` loads `debug.js`, which adds a panel to the page. It is inert without the param.

- **Live line** (p50/p90 ms): `infer` is the whole policy call (de-interleave + `session.run` +
  argmax), and `run` is `session.run` alone. `readback` is the synchronous stall waiting for the
  GPU, i.e. roughly the sim's GPU cost. `submit` is the CPU cost of queueing the passes.
  `step carl` / `step idle` are whole steps with and without an inference. `frame` is the gap
  between animation frames, which shows any main-thread blocking. `rate` is achieved against
  requested steps/s.
- **Bench** pauses the game and times the policy alone. It creates a fresh session per
  execution provider (WebGPU if the browser has it, then WASM) and per policy graph (`orig`,
  `slice`, `gather`; see `?model=`), timing the first run (shader compile / JIT), session
  creation and 20 runs on the live crop. (Reports before `reportVersion` 2 timed the original
  model at 96², 64² and 48² crops instead.) It also runs a
  main-thread test: `maxMainThreadGapMs` ≈ `runMs` means inference blocks the page, so
  pipelining could not overlap it; a few ms means inference runs off-thread.
- **Sweep threads** reloads the page with `?threads=1,2,4,max` and benches each. WASM's thread
  count is fixed for the life of a page, so this can't be done in one load.
- **Copy report** gives one JSON blob: device/browser/GPU info, cross-origin isolation, the live
  stats, and every stored bench. Benches accumulate in `localStorage` across reloads until
  **Clear benches**.

Protocol per device:

1. Open `index.html?debug=1` and wait for the model to load. If **COI** shows `NO`, reload
   once (the service worker takes over on the second load). If it stays `NO`, that is a finding.
2. **Sweep threads** and wait until it says done (about a minute on a phone). No need to copy
   anything yet: the benches are stored and go into every later report.
3. Re-open the plain `?debug=1` URL, so the sweep's last `?threads=` doesn't stick. Tap **Reset
   live**, then play ~90 s in the default "CARL acts sometimes" mode, steering every few
   seconds. The live stats split steps by whether CARL inferred (`step.carl` / `step.idle`), so
   this one run gives both the full cost and the sim-only cost. "You act" adds nothing beyond
   the idle steps.
4. **Copy report** and paste it, labelled with the device name.

What the numbers decide:

- `run96` per EP across devices → whether a smaller network is mandatory. That means
  distillation, or the coarser action map.
- `blocking` → whether pipelining can work on that device, and with which EP.
- `infer` vs. `readback` on CARL steps → the pipelining payoff (`max` vs. `sum`), and whether
  the GPU tiers above matter more than inference.
- `run64` / `run48` (version-1 reports) → what a smaller-window retrain would buy.
- `run` per `model` and `ep` (version-2 reports) → which policy graph to default to per backend.
- `infer` − `run` → the JS de-interleave cost (expected to be negligible).
- `webgpuBackendUp`, `crossOriginIsolated`, `wasmThreads` → configuration problems, which
  are free to fix.

## Device measurements (2026-10-02)

Six devices, `?debug=1&sound=0`, default "CARL acts sometimes" mode at level 3 and 60 steps/s
requested, plus a thread sweep on each. Medians in ms. "GPU sim" is `gpu.readback`, the stall
waiting for a step's queued passes, which approximates the sim's GPU cost.

| device | GPU (WebGPU adapter) | WebGPU 96 | WASM 1 / 2 / 4 / 8 threads, 96 | fastest CARL | GPU sim / step | step with CARL / idle | steps/s |
|---|---|---|---|---|---|---|---|
| Windows laptop, Chrome | AMD 890M (rdna-3) | 6.7 | 23 / 12.5 / 7.0 / 7.7 | tie, ~7 | 3.6 | 10.5 / 4.3 | 59 |
| Windows laptop, Chrome | Intel Xe3-LPG | 8.1 | 21 / 12.4 / 8.8 / 6.2 | WASM 8t | 6.3 | 15.7 / 11.6 | 59 |
| MacBook, Chrome | Apple M2 (metal-3) | 8.3 | 25 / 13.2 / 15.2 / 18.2 | WebGPU | 8.6 | 17.2 / 10.4 | 60 |
| Pixel, Chrome | Mali-G715 (valhall) | 24–30 | 38 / 30 / 17 / 60 | WASM 4t | 18.5 | 44.7 / 21.3 | 22 |
| iPhone, Chrome (WebKit) | Apple GPU | ~30 | 52 / 66 / 70 / – (4 cores) | WebGPU | 20.2 | – / 20.4 ¹ | 28 |
| old MacBook, Chrome 128 ² | Intel Iris Pro (gen-7) | 47, **wrong output** | 45 / 26 / 15 / – (4 cores) | WASM 4t | 26.3 | 47.8 / 25.4 | 35 |

¹ No CARL steps were recorded in the iPhone's live run (CARL stayed idle); its CARL cost comes
from the bench. ² Old MacBook: run with `?ep=wasm` because WebGPU inference there returns wrong
actions (interventions land ~3–4 cells above Pac-Man; WASM on the same machine is correct). Kept
for completeness only, not an optimization target.

Smaller crops, at each device's best WASM thread count (64 / 48): Windows AMD 3.7 / 2.5,
Pixel 8.5 / 5.6, iPhone 24 / 15, M2 6.2 / 4.1, Intel Xe3 3.6 / 2.7. WebGPU at 64 / 48 is within
~1–3 ms of its 96 figure on every device except the old MacBook.

### Findings

1. **On phones the GPU sim is as expensive as CARL.** 18–20 ms per step on Pixel and iPhone
   before any inference, which caps them at ~50 steps/s with CARL idle. Desktops are 4–9 ms.
   The GPU work list above (Tier 1 #1–#4, Tier 2 #6) is now as important as anything on the
   CARL side, not a separate concern.

2. **WebGPU inference has a floor of ~6.5–8.5 ms, even on fast desktop GPUs, and it barely
   moves with crop size.** So its cost is per-node dispatch overhead in onnxruntime's JS, not
   arithmetic: the 694-node graph is the problem. The main-thread test agrees: WebGPU holds the
   main thread for most of each run on the Pixel (max gap 20–27 ms of ~30). The lever for WebGPU
   is the graph simplification (CARL inference #3). A smaller network or crop does not help it.

3. **WASM inference is arithmetic-bound.** 96 → 64 → 48 gives ~2× and ~3–4× on every device, so
   distillation and a smaller-window retrain pay off on WASM. WASM also blocks the main thread
   for the whole run (zero ticks in the test), so pipelining (CARL inference #2) needs inference
   moved into a worker first.

4. **The best WASM thread count is device-specific, and more threads can be slower:** 2 on M2,
   4 on Pixel, 8 on Intel Xe3, 1 on iPhone, and 8 on the Pixel is 3.5× slower than 4. The
   current default, `min(hardwareConcurrency, 6)`, oversubscribes phones with big.LITTLE cores.
   4 is the best single default: best or near-best everywhere WASM is the faster backend.

5. **No single backend is fastest everywhere.** WebGPU wins on M2 and iPhone and ties on AMD;
   WASM wins on Pixel (17 vs 25–30 ms) and Intel Xe3 (6.2 vs 8.1). WebGPU stays the right
   default: it is best or within ~2 ms on every target device except the Pixel.

6. **Pipelining pays more than first estimated.** CARL and the sim cost about the same on most
   devices (M2 8.7 vs 8.6, Pixel on WASM 17 vs 18.5), so overlapping them approaches the 2×
   ceiling on CARL steps rather than the ~11% feared when inference was assumed to dominate.

7. **The two newest desktops reach 60 steps/s, with no headroom.** A CARL step is 15.7–17.2 ms
   against a 16.7 ms budget.

8. **The WebGPU output bug is not Mac-wide.** The M2 MacBook runs WebGPU by default at full
   speed. If its steering is confirmed correct, the bug is specific to the old Intel Iris Pro
   driver.

### Revised order

1. Default WASM threads to 4. WebGPU stays the default backend. (Done.)
2. GPU sim Tier 1 (#1 crop read only when inferring, #2 row extents, #3 weights in shader, #4
   one 1×1 readback) and Tier 2 #6 (step the dots less often): targets the phones' 18–20 ms.
   (Done; not yet re-measured on the devices.)
3. Graph simplification (fold FiLM, static shapes, conv padding), mainly for WebGPU's
   per-node floor. Moved ahead of pipelining by the second round of measurements (below).
   (Done as the FiLM split + Slice/Gather cores; not yet measured on the devices.)
4. Inference in a worker, pipelined with the sim: CARL steps cost `max(infer, sim)`. Then
   distillation / a smaller window (retraining) for WASM.

### After the GPU sim fixes (same day, second round)

Kernel uniform buffer, crop read only when inferring, one gathered readback, and dots stepped
every 3rd step on phones (`dotsim` 3; desktop stays at 1). Medians in ms:

| | Windows (AMD) | iPhone | Pixel | Pixel, `?kernel=tex` |
|---|---|---|---|---|
| GPU sim per step, before → now | 3.6 → **2.3** (−37%) | 20.2 → **11.0** (−45%) | 18.5 → **9.5** (−49%) | 12.4 |
| step without CARL, before → now | 4.3 → 2.6 | 20.4 → 13.6 | 21.3 → 9.7 | 12.5 |
| CARL inference (WebGPU, incl. crop read) | 7.3 | 34.2 | 35.0 | 34.7 |
| step with CARL, before → now | 10.5 → 9.5 | – → 39.1 | 44.7 → 45.3 | 47.9 |

- Both halves earn their place. On the Pixel, the dots change alone takes the sim from 18.5 to
  12.4 ms (the `kernel=tex` run), and the uniform-buffer kernel takes it on to 9.5 ms. The UBO
  path stays the default.
- The crop read now sits inside `infer` (`infer` − `run` grew to ~2.6 ms on the Pixel and
  ~0.5–0.8 ms elsewhere). It is paid on inferring steps only, rather than every step.
- The Pixel's CARL step did not improve (44.7 → 45.3 ms): WebGPU `run` measured ~32 ms this
  round against ~29 last time, probably run-to-run noise or a warmer phone. That ate the
  sim's 9 ms saving, along with the crop read.
- `measSps` in a report is a snapshot taken at the moment of Copy, so ignore it; compare step
  times instead.

**Consequence for the order.** On phones CARL is now ~3× the sim (~33 vs ~10 ms), so pipelining
hides only ~10 ms of a 40–45 ms CARL step. Graph simplification moves ahead of it: it attacks
WebGPU's per-node overhead directly, and the iPhone has no faster backend to switch to.
Pipelining follows, and pairs best with WASM on the Pixel (17 ms at 4 threads, against ~32 ms
for WebGPU).

### Debug-tool flaws seen in these runs (fixed since)

- `Session already started` (some WebGPU bench rows): the bench can start while the game's own
  warm-up inference is still running, because `loadModel()` sets `session` before awaiting
  the warm-up. The bench should wait for in-flight runs.
- `memory access out of bounds` on the last step of a sweep (Pixel 8 threads, iPhone 4 threads):
  most likely sessions accumulating across benches. A session that errors is never released.

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

The CARL inference numbers come from loading `models/agent_direction.onnx` with the `onnx` Python
package, running shape inference at a 1×4×96×96 input, and summing `out_H × out_W × C_out × C_in ×
3 × 3` over the 15 convs. Timings are native `onnxruntime` CPU, `intra_op_num_threads = 1`, 40
runs after 5 warm-ups, on random input (latency does not depend on the values). The closed-loop
table uses an FFT Lenia built from `kernelData()`'s arithmetic, `action.glsl`'s disc rule, and
`crop.glsl`'s frame order and origin. The scripts were throwaway and are not in the repo.
