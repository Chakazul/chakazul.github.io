# Cognoship: plan

A spaceship-piloting game. The ship is a multi-channel Lenia creature (the VT049W "Aquarium"),
and the only way to steer it is **cognitive masking**: briefly hiding part of the creature from
its own convolution, as in Cool et al. 2026. Carlman steers its creature with a trained RL agent
(Cvjetko et al. 2026). Here the player is the controller and has no ML help. The skill is learning
where, when and in which colour channel to mask the creature so that it turns, changes form, or
survives.

The project has three phases. **Phase 1** builds the ship and proves that masking can steer it.
After that, the next step can be either **Phase 2** (design the game on the world view) or
**Phase 3** (the mask map and a CARL-like agent). Phase 3 is part of the game, and is also a
research agenda in its own right. Phases 2 and 3 don't depend on each other, and each can feed the
other later.

Status (revised 2026-11-13):
- Phases 1 and 1b are done and checked in the browser, on desktop and phone (2026-11-04).
- Phase 3 was chosen (2026-11-05). Its measuring half is built in the lab (`lab.html`): stats
  mode, the phase-space chart, long histories, and live morph detection (`detect.js`). Its
  press-and-outcome half (press records, mask maps, sweep, agent) is not started.
- A game was built anyway (2026-11-06 to 11-13), in a different shape from §5: a scripted story
  (`script.txt`, 7 chapters plus a final one, all 11 species) and a sandbox mode (`sandbox.txt`),
  with the morph book, dots and energy, bullets, cloning, and twins and trios followed as wholes.
  §7 and §8 mark what of §5 and §8 it covers.
- Not built: §9's performance work, and most of §8.2's bigger mechanics.

---

## 1. Requirements (paraphrased from the brief)

### Phase 1: prove that masking can steer the Aquarium

- **One page, WebGL2, multi-channel Lenia.** Model it on `../Carlman/` (GPU engine split into
  shader passes) and on `https://chakazul.github.io/Lenia/WebGL/index4.html` (shadertoy-style
  ping-pong sim plus a draw pass). Using WebGPU is allowed if its advantages outweigh its costs
  (decided in §3).
- **Canvas.** One large rectangular canvas, about the size and aspect of Carlman's, on a black
  background.
- **Left 2/3: the world.** 3-channel Lenia (R, G, B) on a torus (repeating borders). It starts
  with a single Aquarium (VT049W) in the centre. The creature must glide steadily forward and must
  not self-replicate.
- **Right 1/3: the zoom panel.**
  - The same creature in real time, heavily zoomed and visibly pixelated (nearest-neighbour).
  - The view is centred on the creature's centre of mass (CoM) each frame.
  - A **yellow arrow** shows the creature's velocity: its direction is the heading and its length
    is the speed. The arrow is hidden when the speed is below a threshold (near stationary).
- **Channel picker** at the top of the panel: three large red, green and blue patches. Clicking a
  patch selects that channel. Red is selected by default.
- **Cognitive masking.** While the mouse (or a finger) is held down on the zoomed creature, a small
  circle under the pointer becomes the occlusion mask **B** for the selected channel. That area is
  left out of the computation of the potential **U**. Releasing the button ends the masking. Under
  masking the creature may turn, change phenotype (for example to a rotating or stationary form),
  or die.
- **The mask is relative to the creature, not to the grid.** Cool's Fig. 4 fixes a small mask at
  an absolute grid position. Here the mask sits at an offset from the creature's CoM, so it moves
  with the creature. The offset is the pointer's position in the zoom panel, which is already
  CoM-centred.
- **Purpose.** Let the player experiment with when, where and in which channel to mask, to reach a
  goal (change direction, change phenotype). In effect, the creature becomes a spaceship.

### Phase 1b: choosing the pattern, and several creatures at once

- **Pattern picker.** Choose which pattern to spawn. The default is the Aquarium (the gliding
  VT049W). Add a **self-replicating Aquarium**.
- **Focus by clicking.** A self-replicating Aquarium makes many creatures, and the default one
  can too after heavy masking or erasing. Clicking a creature in the world view (left panel)
  makes it the focused ship: the zoom view, the arrow and the tools then follow it.
- **When the focused creature splits** into two, keep the focus on one of them.
- **When the focused creature dies**, nothing is focused and the ship view goes blank until
  another creature is clicked.
- **When every creature has died**, show the "Dissolved" screen.

### Phase 2: design the game (left 2/3)

Once Phase 1 shows the Aquarium can be controlled, design the game played on the world view:
which Lenia game elements to use (walls or obstacles that remove mass, walls or obstacles that
mask, other creatures, a goal, and so on) and what the objective is. §5 lists candidates. The
choice is made after Phase 1.

### Phase 3: the mask map and a CARL-like agent (part of the game, and a research agenda)

- **Define forms quantitatively.** Find the creature's forms (attractors) as clusters in a phase
  space of metrics such as mass and speed, rather than naming them by eye.
- **The mask map.** Record every mask press as (state, action) → outcome, and show what has been
  learned as maps per form and channel. This is **part of the game UI** (a panel the player can
  hide), not a debug tool. It is the player's own map of the creature, built by experimenting.
- **The analogy with CARL.** The player is learning the same kind of mapping that CARL's network
  learns: from a high-dimensional pattern to a low-dimensional intervention that reaches a goal.
  The mask map is a low-level, readable representation of that mapping, in a natural mind or an
  artificial network.
- **Later:** a learned co-pilot, and a CARL-like RL agent that steers by masking, whose behaviour
  can be compared with the players' maps. §6 has the design and the research questions.

### Constraints carried over from this repo

- Static and no build step. Shaders are fetched, so the page is served over http(s), as Carlman
  is.
- The page is also played on **mobile web**. Use `precision highp` throughout. Desktop runs cannot
  validate GPU precision (see `../Carlman/shaders/action.glsl`). Make touch input work like mouse
  input.
- No headless-browser verification. Cheap checks only: `node --check`, and throwaway node scripts
  for maths.

---

## 2. Resources

| Resource | What it gives us |
|---|---|
| `../assets/papers/Cool 2026 Lenia agnosiophobia.pdf` | The occlusion rule (below), Fig. 4 targeted-occlusion sensitivity maps (a 3×3 mask per body pixel, giving heading change, recovery time, distortion, and death or explosion), the recovery/phenotype metrics, and the dynamical-systems framing (attractor, cognitive basin, heading as a free variable). |
| `https://github.com/jessescool/lenia-umwelt` | The paper's Python/torch code: `substrate/` (Lenia update with occlusion), `metrics_and_machinery/` (Wasserstein profile distance, interventions), `experiments/sweep.py` (the Fig. 4 sweep), `environments/` (11 barrier masks: funnel, corridor, pegs, shuriken, box, capsule, chips, guidelines, membrane, noise, ring). |
| `https://jessescool.github.io/lenia-umwelt/` | The authors' JS demo (FFT, single channel): click to occlude, transient or fixed occlusions, circle/triangle/square brushes. A reference for how occlusion feels interactively. |
| `../Carlman/` (`glsim.js`, `shaders/*.glsl`, `README.md`, `PERF.md`) | A proven WebGL2 engine on this repo's targets: R32F texture ring, GPU CoM reduction (`reduce.glsl` → `com.glsl`), a crop shader that reads the CoM from a texture so there is no CPU stall, wall masking, a mobile-tested precision setup, and fullscreen/theme UI. |
| `../Lenia/WebGL/index2.html`, `shadertoy/lenia5.glsl` | Multi-channel Lenia in GLSL: per-kernel `src`/`dst` channel tables packed into `mat4`s, with VT049W as `species 0` (labelled "sometimes reproductive" there). |
| `~/github/Lenia/Python/found/233.json` (+ `LeniaNDKC.py`) | **The source used.** Chan's "aquarium" collection, which holds rule `VT049W<XFKA9P` at R = 25 with several named forms of the *same* rule: *rectus* (straight glider, our ship), *gyrans* (circling), *lithos* (stone), *tachus*, *para*, and so on. `LeniaNDKC.py` defines the exact update: per-channel growth normalised by Σh. |
| `~/github/chakazul.github.io/lenia-CE/pattern_list.js` | Another VT049W ("Amazing cell(s)", R = 12). It has slightly different `h` values and its own dt and scale quirks, so it is not used. |
| `../Cvjetko 2026`, `../Hamon 2025` (same folder) | Context only. CARL is used by Carlman. Hamon's sensorimotor agents are candidates for Phase 2's other creatures. |
| `../DESIGN.md` | Earlier game-design notes, including "Cool: herding without ML". |

### The occlusion rule (Cool et al., Eq. 2) and how it extends to 3 channels

Single channel: **U = K∗(A·(1−B)) / K∗(1−B)**. Masked cells are set to 0 in the numerator. The
denominator renormalizes by the share of kernel mass that falls on visible cells, so the
creature "trusts the visible part more" and does not just see zeros. With no mask in reach, this
is the ordinary rule.

Multi-channel: each kernel *k* reads a source channel `src[k]` and writes growth to `dst[k]`. A
mask on channel *c* (`B_c`) hides channel *c* from every kernel that reads it:

```
U_k  = Σ_taps w_k · A_src[k] · (1 − B_src[k])  /  Σ_taps w_k · (1 − B_src[k])
A_c' = clip01( A_c + dt · Σ_{k : dst[k] = c} h_k · G_k(U_k) )
```

If the denominator is 0 (the kernel is fully masked), use U = 0 or skip the growth term; choose
one and document it. The mask hides information only. It does not remove mass, so masked cells
still update from their own (unmasked) neighbourhood. That difference from Carlman's walls, which
remove mass, is what the game is built on.

---

## 3. Technology decision: WebGL2, not WebGPU

**WebGL2** wins for this project:

- Carlman's engine, its shader idioms and its mobile fixes carry straight over. A WebGPU rewrite
  reuses none of them.
- WebGL2 runs on effectively every phone browser. WebGPU on mobile is still patchy, and we
  already have an open Mac-specific WebGPU bug (Carlman's policy backend).
- The sim is small: one 3-channel board of about 384×256. A direct-convolution fragment shader
  handles it easily, so WebGPU's compute shaders, shared memory and storage buffers buy little.

WebGPU would be worth reconsidering only if we needed FFT convolution (much larger kernels or
boards) or many creatures at once. Neither is in scope.

---

## 4. Phase 1 design: the controllable ship

### 4.1 Creature and world

- Rule and cells: VT049W *rectus* from `Lenia/Python/found/233.json`, resampled from R = 25 to
  **R = 12**, with **dt = 0.2**.
- **Result of the step 0 check.** No tuning was needed: it glides for 5,000 steps, does not
  replicate, and its mass holds at about 336. The details, the shader-vs-numpy match, and a coarse
  numerical mask sweep that already shows strong steering are in README-MORE.md, "Fidelity checks".
  The original plan for the check follows.
- **Step 0 is a fidelity check.** Run it free and confirm that it glides steadily, does not
  replicate, and keeps a stable mass. VT049W is known as "sometimes reproductive". If it
  replicates, try, in this order:
  1. a different starting orientation, or an initial state taken after settling;
  2. a small σ/μ tweak on the kernel driving fission (as `lenia5.glsl` did for Z18A9R);
  3. a scale-up (R = 18 or 24, resampling the cells), which also makes the mask finer relative to
     the body. Cool upscaled creatures 4× for this reason.

  Record the variant we choose in `species.js`, with notes.
- World: 3 channels, toroidal, about **384×256** cells for the left 2/3. Exact size is chosen so
  the cell size is an integer number of screen pixels. Carlman's board is 550×250 at R = 18, and
  scaled to R = 12 this is about the same physical size relative to the creature.

### 4.2 GPU passes (`glsim.js` + `shaders/`)

| Pass | Job |
|---|---|
| `sim.glsl` | Direct convolution over a (2R+1)² window. Each tap does **one** RGBA board fetch, which serves all 15 kernels, plus one mask fetch. Kernel weights live in a small texture or UBO (15 weights per tap, precomputed and normalised on the CPU). It accumulates the numerator and, only when any mask is active, the denominator, then applies growth and clipping. The board is a ping-pong pair of `RGBA16F`/`RGBA32F` textures with `REPEAT` wrap. |
| `mask` (inside `sim.glsl`) | B_c(x) = the level mask texture (empty in Phase 1, used by Phase 2's masking walls) OR the **dynamic circle**. The circle is computed analytically from `CoM + offset`, radius and channel; the CoM is read from the CoM texture (Carlman's crop trick), so placing the mask needs no CPU round trip. Distances are torus-wrapped. |
| `reduce.glsl` → `com.glsl` | Mass and **wrap-aware** CoM. Each axis uses a mass-weighted circular mean (Σ m·cos θ, Σ m·sin θ), so the CoM stays correct when the creature crosses an edge. Output is a 1×1 texture. |
| `draw.glsl` | Left view: the board as RGB on black, optional faint outline of the mask. Right view: a crop of about 48×48 cells around the CoM, scaled to the panel with `NEAREST` (heavily pixelated). |

Readback: one tiny 1×1 CoM/mass read per rendered frame (not per step), used for the arrow and
Phase 2 logic. Phase 3's form features are added later as extra sums in the same reduction. Use
async readback (`fenceSync` plus a PBO) if it shows up in profiles.

Performance budget: 384×256 × 625 taps ≈ 61 M fetches per step. That is fine on desktop for
several steps per frame and should allow 1 to 2 steps per frame on mid-range phones. Optimisation
if needed: compute the denominator only in a band of R + r around the mask circle and level-mask
tiles.

### 4.3 Overlay and UI (`core.js`, 2D canvas over the WebGL canvas)

- **Velocity arrow.** Velocity is the CoM change per frame, torus-unwrapped, drawn from the panel
  centre and scaled by a fixed gain. It is hidden below `speedMin`; the threshold is tuned in step
  0 from the free-gliding speed.
  - **As built (rev. 2):** a short EMA (`ARROW_TAU` = 2.5 unit time, ~5 steps) that starts as a
    running mean, so the arrow has its full length from the first frame instead of ramping up.
    The EMA's response to a steady turn at rate ω, `a / (1 - (1-a) e^(-iωT))`, is divided out
    (rotation capped at 60°, gain at 2). ω is the turn rate of that short average's direction,
    averaged over `OMEGA_TAU` = 5. (Taking ω from raw step-to-step directions instead was tried:
    capping detection jumps there breaks the sum of the angle changes, and the view drifted up to
    60° on a straight glider.) On synthetic circles (periods 10–100 t) the arrow is exact once ω
    has settled. `speedMin` = 0.005 R, low enough for the *tardus* forms.
  - Rev. 1 used one EMA over 15 steps for both arrow and view. It ramped up from zero after a
    focus change (95% only after ~45 steps) and on a turning ship trailed by atan(ωτ), shrinking
    by its cosine: ~55° for a 30 t circle, near perpendicular for a 10 t one.
- **Channel patches.** Three large R/G/B squares at the top of the panel, with the selected one
  highlighted. Keys `1`, `2`, `3` select too.
- **Masking input.** `pointerdown` inside the zoom view starts masking and `pointerup`, `cancel`
  or `leave` ends it. While masking, the offset follows the pointer, converted from panel pixels
  to cells relative to the CoM. The mask circle is drawn as an outline in the selected colour. A
  slider sets the mask radius (default about 2 to 3 cells at R = 12, roughly the paper's 3×3).
- **Heading-up toggle.** The zoom view has a toggle that turns it to keep the ship pointing up.
  Off (the default) is board-relative, as in the brief. On, the crop rotates with the ship, so a
  press on the left flank stays on the left flank as the ship turns. Pointer offsets are rotated
  back into board coordinates before they reach the shader, so the mask lands under the pointer
  either way.
  - **As built (rev. 2):** the view angle is the raw velocity direction averaged *as an angle*
    over `VEL_TAU` = 15 unit time (30 steps), predicted forward by ω each update, so a steady turn
    is followed with no lag and pulses do not shake it. Each pull towards the raw direction is
    capped at 45° (`JUMP_MAX`) against detection jumps (5N7KKM paragyrans, VT049W para): with
    ±10° noise plus a 70° jump on 10% of frames, the view stays within 14° of a straight glider's
    heading. A full alpha-beta tracker (ω corrected from the same residual) was tried and dropped:
    at this smoothing its ω adapts too slowly and it slipped half a turn on a 30 t circle. A vector average was not used: for a fast spinner it shrinks towards zero
    and its direction stops meaning anything. It updates while the ship is moving (hysteresis on
    the arrow's speed: on above 1.5 × `speedMin`, off below `speedMin`) and holds the last heading
    when near stationary. A sudden turn onset (30 t circle) is trailed by up to ~45° for ~40
    steps while ω catches up. No rate limit and no freeze.
  - **What was tried and dropped:**
    - *Body axis as heading* (B-channel CoM → G-channel CoM). It exists at any speed, but it is
      not the direction of travel: it sits ~56° off for gyrans and wobbles ±5° for rectus, so the
      arrow pointed obliquely.
    - *Axis plus a learned axis-to-velocity offset.* The offset came from noisy per-frame CoM
      deltas, so while running the target kept moving.
    - *A rate limit* (90° and then 180° per second). With a moving target, the view lagged and
      only settled when paused.
    - *Freezing the view when the ship rotates in place*, with a "Follow spin" override. It
      misfired on masked turns, which slow the ship while its body turns.
  - The velocity arrow is drawn in view space; since it is the faster estimate, with heading-up
    it points about up and swings a little around it.
- **Status readout:** mass, speed, heading, and a state of *alive*, *dissolved* (mass below a
  fraction of its spawn mass) or *exploded* (mass above a limit). Phase 3 adds the form label
  (§6.1).
- Controls: Space to pause, `.` to step one frame, Respawn (Enter), a sim-speed slider, and the
  fullscreen and theme buttons as in Carlman.
- URL params (as built): `?form=`, `?R=`, `?dt=`, `?speed=`, `?radius=`, `?headup=`, `?debug=1`. `?debug=1` is for
  developer timings only. The mask map is not a debug tool.

### 4.4 Phase 1 exit criteria: is the ship controllable?

- [x] Free run: the creature glides for at least 5,000 steps with no fission and steady mass.
- [x] Manual: repeated mask presses can turn it left and right on demand.
- [x] Manual: at least one press pattern visibly changes its form (for example to rotating or
  stationary), and ideally back.

### 4.5 Phase 1b design: pattern choice and multiple creatures

- **Pattern picker.** A select box (and `?form=`) listing the patterns in `species.js`, with the
  Aquarium first.
- **The self-replicating Aquarium.** Find a seed that reliably divides, and check it in numpy
  first, as in §4.1. Candidates to try:
  - the other forms of the same rule in `found/233.json` (the *para* and *+* entries);
  - the "Amazing cell(s)" VT049W in lenia-CE;
  - the VT049W labelled "sometimes reproductive" in `Lenia/WebGL/shadertoy/lenia5.glsl`.

  Prefer a seed of the same rule, so it shares the board with the Aquarium. Record how often it
  divides and whether the copies keep dividing until the board fills: a board that fills up is an
  explosion in all but name.
- **Finding creatures.** Once per frame, read back a coarse mass map, for example 4×4-cell block
  sums (96×64 for a 384×256 board), from a reduction pass. On the CPU, label its connected
  components with torus wrap and a mass threshold. Each component above a minimum mass is a
  creature, and its CoM comes from the block sums. That is cheap (about 6k blocks) and needs no
  per-creature GPU work.
- **Focus and tracking.** The focus is one component, followed from frame to frame by nearest
  CoM within a gating distance. Clicking in the world view focuses the component under or nearest
  to the click.
  - The focused creature's CoM must reach the GPU in place of the whole-board CoM, both for the
    mask and for the zoom view. Either upload it as a uniform each frame (about one frame of lag,
    which is fine at about 0.1 cells per step), or run the existing CoM reduction over a window
    around the last focus position (Carlman's window approach), which keeps it fully on the GPU.
  - Velocity and heading reset when the focus changes.
- **Splits.** When the tracked component becomes two, keep the one nearer the previous CoM, or the
  heavier one if that is ambiguous. Velocity is kept, since the heading of either half carries on.
- **Death.** If the focused component disappears or drops below the minimum mass, set focus to
  none. The zoom view goes blank, the arrow is hidden, tools do nothing, and the panel says
  "click a creature to follow it". If no component is left on the board, show the "Dissolved"
  screen. As built, there is no "Exploded" stop: a board full of copies keeps running.

### 4.6 File layout

```
Cognoship/
  PLAN.md         this file
  index.html      the game page (game.js)
  lab.html        the lab page (lab.js): every experiment control
  style.css       styles, both pages
  page.js         page chrome, both pages: theme, fullscreen, music, helpers
  core.js         shared: canvas input, focus tracking, stats, detection, overlay
  lab.js          the lab's UI, readout, stats chart, recording, tuning
  game.js         the game's UI and HUD (game logic to come)
  forms.js        (Phase 3) phase-space features, clustering, form labels
  maskmap.js      (Phase 3) press records, map derivation, JSON import/export, map panel
  glsim.js        WebGL2 engine (passes, textures, readback)
  species.js      VT049W kernels + cells (+ any tuned variant)
  shaders/        vertex, sim, reduce, com, draw (.glsl)
  README.md       controls, URL params, fidelity notes (when built)
```

### 4.7 Stats mode and spin (stats mode built, see §7 Phase 3 and README; spin not built)

Stats mode (key `A`, `?mode=stats`) records the focus every sim step: mass, mass per channel,
the arrow's speed and its turn rate, in R units (README, "Stats mode"). The turn rate is how fast
the *direction of travel* turns. It is not the creature's spin: a creature spinning in place has
an almost stationary CoM, so its velocity, and the turn rate, carry no signal. Single-channel
Lenia has many such spinners with exact n-fold rotational symmetry (2-fold to 10-fold and more).

**Spin from complex angular moments.** About the focus CoM, with each cell at z = x + iy and mass
m, for k = 2 … 12:

  Z_k = Σ m · (z / |z|)^k = Σ m · e^{ikθ}

the k-th angular Fourier coefficient of the mass. A rotation by α multiplies Z_k by e^{ikα}, and
for an n-fold symmetric creature every Z_k is 0 except at multiples of n, so Z_n carries the
signal. Per step:

  spin = arg( Z_n(t) · conj(Z_n(t − 1)) ) / (n · dt)        (rad per unit time)

then smoothed like the turn rate (OMEGA_TAU). The product of consecutive values gives the phase
change directly, so no phase unwrapping or long-term tracking.

- **Choosing n:** the k with the largest |Z_k| / Σm (0 to 1). Record it as `fold` next to `spin`.
  When n changes between steps, skip that step's spin (the phases of different k don't compare).
- **One method for all creatures:** k = 2 is the familiar long-axis (second moment) method,
  which is why it fails for 3-fold and higher shapes (their Z_2 is 0); k = 1 is always 0 about the
  CoM. An asymmetric ship has a strong Z_2 or Z_3, so the same code gives the ships' spin too,
  without the B→G channel-CoM axis of the "body axis up" view (which needs channel CoMs at least
  a cell apart, so it doesn't apply to single-channel creatures).
- **Limits:** a turn of more than π/n per step aliases (18° per step at n = 10; spinners are
  expected to be far slower, not checked). A nearly round creature, or one with little contrast
  between its arms, has a tiny |Z_n| and a noisy phase: no spin below a cut-off, e.g.
  |Z_n|/Σm < 0.02. CoM jitter adds noise at high k; the unit weight (z/|z|)^k, not z^k, keeps the
  outer cells from dominating, and cells within about 1 cell of the CoM are left out (undefined
  angle). Pulsing changes |Z_n|, not its phase. Spiral arms (rings with different phases) still
  work: a rigid rotation shifts every ring's phase equally. What it measures is the rotation of
  the pattern, the only spin there is in Lenia (no material moves).
- **GPU cost:** one more two-stage reduction like reduce → com, summing m·cos kθ and m·sin kθ for
  k = 2 … 12 (22 numbers, per channel or summed) within the focus window around the focus CoM,
  read back in the same step as the blocks.

---

## 5. Phase 2: game design candidates (an option after Phase 1)

Elements, all based on Lenia:

- **Mass-removing walls and asteroids:** Carlman-style hard hazards that erase mass on contact.
- **Masking fields ("nebulae"):** static per-channel occlusion regions from the level mask.
  Agnosiophobia makes the ship turn away from them without player input, so they act as soft
  bumpers and guide rails. Shapes can borrow from lenia-umwelt's environments (funnel, corridor,
  pegs, ring).
- **Beacons and goals:** a target region to reach. A goal could also require arriving *as* a given
  phenotype (for example, park as a stationary form).
- **Other creatures:** single-channel solitons living in one channel (they interact through the
  cross-channel kernels), or Hamon-style sensorimotor agents as autonomous traffic or predators.
- **Mask energy:** total mask-time as fuel, which rewards precise, minimal steering.

Candidate objectives:

1. **Courier:** reach a sequence of beacons through an asteroid and nebula field, scored on time
   and fuel used.
2. **Shape-shift gates:** gates that only a given phenotype can pass (for example, a narrow slot
   needs the rotating form).
3. **Survival drift:** stay alive as long as possible as hazards close in. Nebulae help and walls
   kill.

---

## 6. Phase 3: the mask map and a CARL-like agent (an option after Phase 1)

### 6.1 Defining forms from phase space

"Gliding", "rotating" and "stationary" are only qualitative names, and the boundary between two
forms can't be drawn by hand. So forms are defined from data instead: each form is an
**attractor**, and it shows up as a cluster in a low-dimensional phase space of metrics. This is
Cool's framework made operational. Their canonical morphoset (a centre plus `d_max`) is one such
cluster region, measured with a different metric.

- **State features**, all invariant to translation and rotation. Position and heading are the free
  variables, so they must not separate forms. Candidates:
  - mass per channel (3 values);
  - CoM speed;
  - angular velocity of the body axis;
  - radius of gyration (spatial spread);
  - distances between the channel CoMs.

  The reduction pass outputs these as extra sums, so they cost nothing extra to read back.
- **Windowed over the pulse.** VT049W oscillates, so even a steady form traces a loop (a limit
  cycle) rather than sitting at a point. Each feature is taken as its mean and its amplitude over a
  sliding window of about one pulse period, which collapses each loop to a tight blob. Measure the
  period during the first check of the creature.
- **Finding clusters.** Record trajectories from free runs and from the aftermath of many presses.
  Standardise the features, then cluster them (GMM or DBSCAN, run in-browser on the recorded
  points). Hand-drawn boxes in the phase plot are a quick first pass. Each cluster gets a centre
  and a radius (for example, a Mahalanobis distance covering 99% of its points), as with `d_max`.
- **Labels:**
  - **a form:** inside a cluster's region for k consecutive windows (Cool uses k = 5);
  - **transient:** outside every region. This is the creature moving through its basin between
    attractors, the "slow, distorted recovery" that the paper links to heading changes;
  - **dissolved** or **exploded:** absolute mass thresholds;
  - **unnamed form #n:** a new, dense cluster that appears in later data. It is a discovery, and a
    candidate for a name.
- **Phase-space plot (in the map panel, §6.3):** a 2D projection (pick two axes, or PCA) with
  the clusters shaded, the live state as a dot with a fading trail, and press events marked on the
  trail. It shows a press knocking the state out of one blob and where it lands. *Built so far as
  stats mode (§4.7, milestones): 2D/3D axes picked from the stats, trails and saved recordings;
  still to come: clusters shaded, presses marked, PCA.*

### 6.2 Mask map: (state, action) → outcome

The mask map is closely analogous to CARL:

| | CARL (Carlman) | Cognoship |
|---|---|---|
| Intervention | add or remove mass | a cognitive mask |
| Input (state) | the whole pattern (a 96×96×4 crop) | the pattern's phase-space state (+ body frame) |
| Goal | survive and follow a direction | the player's experiment or objective |
| Learner | an RL network | the player's mind |

In both, a high-dimensional state has to be mapped to a low-dimensional action that reaches a
goal. The player builds that mapping in their head, as CARL's network does in its weights. The
mask map is a low-resolution, readable slice of that mapping, natural or artificial.

**Store raw records and derive the maps.** Each press stores one record:

- *state at press:* the feature vector, form label, body-frame heading, and pulse phase;
- *action:* channel, body-frame offset, radius, and duration;
- *outcome:* the state trajectory for M steps after release, the form label at the end, the
  heading change, the recovery time, and whether it survived.

Records come from two sources:

- **Passive:** every player press. The map fills in where the player keeps clicking.
- **Active sweep:** a button that fills gaps off-screen. It copies the state, holds the mask at
  each grid offset, and records the result. Python with the lenia-umwelt code is a fallback.

Form definitions will change as clusters are refined, so maps are never stored. They are
**re-derived from the raw records** each time:

- A **map** is the records filtered by (form at press time, channel), binned by body-frame offset.
  Records pressed during *transient* states stay in the data but go into no form's map.
- Each cell shows one of: heading change, the form at the end (this also shows transitions, for
  example "this spot turns a glider into a rotator"), recovery time, or survival rate. It also keeps
  a sample count and how much the outcomes vary, and is drawn faded until it has enough samples.
  If outcomes vary a lot with pulse phase, add phase as a binning dimension.
- Offsets use the body frame from §4.3, so presses made at different headings fall in the same
  cell. A form without a usable body axis falls back to board-relative and is labelled as such.
- **Export and import as JSON**, for offline analysis and for combining data across sessions and
  players.

### 6.3 In the game UI

The mask map is part of the game, not a developer tool. It can be hidden, but it is the player's
record of what they have discovered.

- **Map overlay on the zoom panel:** a heat map for the current form and selected channel, drawn
  under the creature in the body frame. It is switchable between heading change, end form,
  recovery time and survival. A **Map** button shows or hides it.
- **Map panel** (collapsible, below the canvas or as a drawer, like Carlman's Setup): the
  phase-space plot, the list of known forms with their sample counts, the active-sweep button, and
  JSON import and export.
- In fullscreen, both stay hidden unless the player turns them on.
- URL param `?map=0` starts with both hidden.

### 6.4 Research agenda

**Observation from Phase 1 play (2026-11-04).** Masking a colour seems to *push* that colour, as
if the creature were afraid of the mask. That fits the paper's "agnosiophobia" (fear of the
unknown). Erasing mass is more destructive. The two tools differ in kind, not just strength: the
mask moves the creature by what it can't sense, and the eraser by what it loses. Worth measuring
with the mask maps: per channel, compare the heading change and survival of mask presses against
erase presses at the same offsets.

The mask map turns Cognoship into an instrument as well as a game. In order:

1. **Basin structure of a multi-form creature.** Cool's Fig. 4 maps single-form creatures. With
   VT049W, the maps per form and channel show which spots are lethal, which turn the creature,
   and which move it to another form. The paper's central claim (zones that turn the creature
   border lethal zones; large turns need long, distorted recovery) can then be tested per form
   and per channel.
2. **Channel specificity.** Do the three channels act as different "senses", with distinct turn
   and transition zones? Which channel is the best rudder?
3. **Learned co-pilot.** A small forward model (k-nearest neighbours, or a tiny MLP in JS) trained
   on the records predicts a press's outcome and suggests where to press for a goal.
4. **A CARL that masks.** Train a goal-conditioned RL agent whose action is (channel, offset,
   duration) masking instead of adding or removing mass, offline in the CARL/rlca codebase, and
   ship it like Carlman's policy (onnxruntime-web). Then:
   - derive its mask map from its own rollouts, with the same record format, and compare it with
     the maps players have built (natural versus artificial mapping);
   - let it act as a rival or co-pilot in the game, or as a benchmark ("beat the agent's route").
5. **Progression.** In the game, the player's map (what they have discovered) can be a visible
   collection, or can unlock things: new forms named, new channels trusted, goals that need a
   given transition.

### 6.5 Phase 3 exit criteria

- [ ] The phase-space plot shows at least two separate form clusters, and the form label is stable
  within each one.
- [ ] The mask maps for the gliding form show consistent left and right turn zones on at least one
  channel.
- [ ] Records survive JSON export and import, and maps re-derive correctly after the form
  definitions change.

---

## 7. Milestones

Tick a box when it is done *and* checked. "(numpy)" means checked off-line only, not yet in the
browser.

**Planning**
- [x] PLAN.md: requirements, resources, technology choice, three phases.

**Phase 1: the controllable ship** (§4)
- [x] Find the creature: VT049W *rectus* from `Lenia/Python/found/233.json`.
- [x] Fidelity check (§4.1): glides steadily at R = 12, dt = 0.2, no fission (numpy).
- [x] Shader maths matches numpy, with and without a mask (node emulation of `sim.glsl`).
- [x] Numerical mask sweep: masks turn it by up to ±160°, and some spots kill it (numpy).
- [x] WebGL2 engine: 3-channel sim on a torus, per-channel CoM reduction, draw pass.
- [x] World view (left 2/3) and pixelated zoom panel centred on the CoM (right 1/3).
- [x] Yellow velocity arrow, hidden when near stationary.
- [x] R/G/B channel patches, green by default (best for steering and shape-shifting; the brief said red).
- [x] Press-and-hold masking relative to the CoM (Eq. 2, renormalised), with an outline in both
  views.
- [x] View rotation button: board up / heading up (velocity) / body axis up.
- [x] Erase-mass tool: a second row of R/G/B patches at the bottom of the panel. It removes mass
  in the circle instead of masking, for comparison (an addition to the brief).
- [x] Speed slider (7 slots 1/8–8, ¼ step per frame per unit), mask radius slider, pause/step/respawn, death
  detection, readout, theme and fullscreen.
- [x] Shaders compile and run in a desktop browser (seen in use).
- [x] Performance on a phone (`?debug=1`): about 0.2 ms/step with 1,004 taps at R = 18
  (2026-11-04), so the kernel weights stay in a texture.
- [x] Mobile shader precision checked on a device: the Aquarium behaves as on desktop (2026-11-04).
- [x] Phase 1 exit criteria (§4.4) pass (checked by hand in the browser).

**Phase 1b: pattern choice and several creatures** (§4.5)
- [x] Page text calls the glider "Aquarium (VT049W)".
- [x] Pattern picker (select box and `?form=`), defaulting to the Aquarium.
- [x] Find and check (numpy) a self-replicating Aquarium, and add it to `species.js`. It is the
  same seed under Lenia/WebGL's unnormalised "Fission cells" update at dt = 0.5. It fills the
  board within ~1,000 steps, and no slower setting exists (README, "Checks for Phase 1b").
- [x] Tame the replication rate (e.g. a σ/μ tweak, as `lenia5.glsl` did for Z18A9R), if the
  fast fill is too much to play with. Done by the switch to LeniaF (below): *VT049W gyrans
  (replicating)* holds at 3–4 creatures at R = 18, so no tweak is needed. Earlier, at dt = 0.2 (README, "Checks for Phase 1b"): the only
  replicators found (kernel 5 μ × 0.90–0.955) swell before dividing. μ × 0.955 is in the picker
  as a test pattern for comparison.
- [x] Mobile check (2026-11-04): speed is fine; `?debug=1` shows about 0.2 ms/step with 1,004
  taps at R = 18.
- [x] Coarse mass-map reduction pass and readback (`blocks.glsl`).
- [x] Connected components on the CPU, with torus wrap, giving per-creature mass and CoM (`creatures.js`).
- [x] Click a creature in the world view to focus it; zoom, arrow and tools follow the focus.
- [x] Focus tracking from frame to frame; on a split, keep one half (replayed in node).
- [x] Focused creature dies: no focus, blank ship view, tools disabled, prompt to click one.
- [x] Every creature dies: "Dissolved" screen. No "Exploded" stop: the simulation continues.
- [x] Circling form (gyrans) dropped from the picker: it is only a morph of VT049W.
- [x] Default speed 1 step per frame; Zoom slider for the ship view (2–16 px per cell).
- [x] Page text credits Cool et al. 2026, with links to the paper and the authors' demo.
- [x] README updated.
- [x] Switch to LeniaF.py's rule (polynomial kernel and growth, `free_h`, T = 2), which is what
  233.json was found with. In a numpy port it reproduces the collection: *gyrans R* is the only
  self-replicator among 5N7KKM … GDNQYX, and *lithos* is still.
- [x] Default R = 18 (below 16, *gyrans R* stops replicating and *lithos* moves). Self-replicating
  Aquarium = *gyrans R*, which holds at 3–4 creatures at R = 18. The legacy fission patterns are
  hidden.
- [x] Picker lists all 43 entries of 233.json from 5N7KKM to GDNQYX (later extended to HAESRE, see Phase 3) ("code form", id
  "code-form"), generated into `aquarium233.js` by `tools/gen_aquarium.py`; default VT049W rectus.
- [x] Creature detection as in ~/github/complexity: basins of the potential K∗A (sampled per
  block in `blocks.glsl`), merged across shallow saddles, so touching creatures separate.
- [x] Keys and buttons: Enter pause, Space step, Q/W speed, E view, Z reset, X add the pattern
  at a random spot, C/V previous/next pattern.
- [x] Z18A9R trio, rebuilt from a screenshot and settled in numpy, added after the Z18A9R forms
  (`tools/extra_patterns.json`).
- [x] `tools/from_screenshot.py`: rebuilds a pattern from a screenshot (brightness correction,
  scale search, shape check). Z18A9R lithos and tardus added with it, trio redone with it; VT049W gyrotachus dropped (a
  different species).
- [x] Checked in the browser by hand, desktop and phone, including fullscreen (2026-11-04). GPU cost
  at R = 18: about 1,000 taps per cell, 0.2 ms/step on the phone.

**Then choose one: Phase 2 or Phase 3**
- [x] Decision made: Phase 3 (2026-11-05). Phase 2 stays an option for later.

**Phase 2: the game** (§5; built as a story game instead, see the status at the top and §8)
- [x] Concept: a scripted story (`script.txt`) of chapters, one species each, with goals as
  `[need:…]` conditions; a sandbox mode (`sandbox.txt`); the morph book as the collection.
- [ ] Level mask texture in `sim.glsl` (masking walls / nebulae). Still none of §5's level
  geometry: every world is an open torus.
- [ ] Mass-removing hazards. Only the player's bullets remove mass; nothing in the world does
  (see §8.2 items 6 and 7).
- [x] Goals and flow: `[need:…]` conditions with progress, the SUCCESS box, Dissolved, Retry /
  Reset from the last checkpoint. Scoring is the energy (dots eaten), spent on bullet patterns.

**Phase 3: the mask map and agent** (§6)

Built 2026-11-05: checked in node (syntax, and the chart, tracking and file code on synthetic
data) and used in the browser during development; README "Stats mode" has the details.
- [x] Page modes, cycled by the Mode button or `A` (`?mode=`): normal, stats, tuning. The help
  text shows in every mode, under the mode's controls.
- [x] Tuning mode: detection sliders (merge saddle, threshold, min area), circles around the
  detected objects, and the share of the board's mass inside objects. Threshold lowered to 0.1
  from 0.15 (99–100% of the mass in objects, against 97–100%).
- [x] Detection and kinematics run after every sim step, not once per frame (one stats sample per
  step).
- [x] Mass per channel in the reduction (a second target of `blocks.glsl`, summed per object).
- [x] Five more species from 233.json, HGMLH8 to HAESRE (20 morphs), named "Aquarium shapeshifting / multiplying /
  striped / exploding / zigzagging (CODE)",
  before Cilium in the picker. Not yet checked (numpy or browser) as the first six were.
- [ ] Other feature sums: spread (radius of gyration), inter-channel CoM distances.
- [ ] Spin and fold count of the focus from complex angular moments, into the stats (§4.7).
- [x] Stats per step, in units independent of R (m = R cells, s = unit time, g = mass of 1 m² at
  value 1): time, total mass, R/G/B masses, speed and |turn rate| (the arrow's filters); hue
  derived from the R/G/B masses when drawn.
- [x] Phase-space plot (stats mode, `statschart.js`), in place of the ship view: 2D or 3D (drag
  to rotate), a stat per axis plus a colour axis (min blue to max red, or the hue itself), axis
  labels that cycle the stat on click, Lock axes, and outliers (one- or two-step spikes) left out
  of lines and ranges.
- [x] Short history (last 50 steps, fading) and long history (Record / `S`: 300 steps a press,
  each piece labelled with its morph; Delete deletes the last piece). Time is relative to each
  piece's start.
- [x] Every-object mode: every object tracked with its own filters and short history, drawn as a
  5-sample tail (whole line with time on an axis); hovering one in the world view shows its line;
  a change of focus swaps its history in.
- [x] Long histories saved per species: Download `<species>.json`, committed to `assets/stats/`,
  loaded with the species (skipped at another R or dt); new recordings add pieces to it.
- [ ] Post-press trajectories: mark presses (channel, tool, offset) on the recorded history, and
  record the aftermath of each press.
- [x] ~~Windowed features over the pulse period~~ Superseded: `detect.js` scores the last 50
  samples against each morph's cloud by kernel density, ignoring their order, so the pulse's loop
  doesn't matter.
- [x] Live form labels (§6.1), done differently: supervised, not clustered. Morphs are named from
  recorded long histories (`detect.js`: kernel density + softmax, fuzzy morphs, per-species merge
  saddles; offline 91% of windows named, 99.3% of them right), and twins and trios are named as
  wholes (2026-11-13). Unsupervised clustering to *discover* new forms (§6.1's "unnamed form #n")
  is not built.
- [ ] Press records (state, action, outcome) and JSON import/export (§6.2).
- [ ] Mask maps derived per form and channel; overlay and map panel in the game UI (§6.3).
- [ ] Active sweep button.
- [ ] Phase 3 exit criteria (§6.5) pass.
- [ ] Research (§6.4): learned co-pilot.
- [ ] Research (§6.4): CARL-like agent that steers by masking, compared with players' maps.

## 8. Suggestions: more game mechanics and the rest of the story (2026-11-09; status 2026-11-13)

A design pass over the game as it stands (mechanics in `game.js`, the script in `script.txt`,
the species in `aquarium233.js`, the lab's tools). The premises for the unused species come
from their names, README notes and the test notes in `script.txt`, not from watching them, so
each wants a short lab session before it is scripted.

### 8.1 Where the script stands

| Chapter | Species | Teaches |
|---|---|---|
| 1 Metamorphosis | 5N7KKM | masking, morphs, the morph book |
| 2 You are not alone | 5N7KKM | dots / energy, the gun, bullet types |
| 3 Brave new world | GDNQYX | clone, recycle |
| 4 Megastructure | Z18A9R, PJ1EF7, XEH4YR | object count, survive a swarm, find all morphs |
| 5 Interstellar Evolution | KZ2BLD, KH97WU, HAESRE, HGMLH8 | find a morph, engineer twin morphs |
| 6 The Rise of Self-Replicators | VT049W | dissolve every object with Paraptera bullets |
| Final: Cosmic Immune System | G6G6CR | make red mass (`[need:mass:r:50]`) by cooperation |

Done 2026-11-13: the story uses all 11 species and has an ending (credits, then sandbox mode).
The typos listed here earlier are fixed.

### 8.2 New mechanics

Roughly from cheap (reusing what exists) to bigger. Status 2026-11-13 in **bold brackets**.

1. **[Half built: `[need:survive_time]`; `hold` not built.]** **Time goals:** `[need:survive:<s>]`, and `[need:hold:<morph>:<s>]` (stay that morph for
   s seconds; hard with unstable morphs, a good test of steering). The first is already wanted
   in the chapter 5 notes. Both are small additions next to the existing conditions.
2. **[Not built.]** **Beacons:** `[beacon:<x>,<y>]` and `[need:reach]`, a target marked in the world view.
   Steering by masking is the core skill, yet no chapter asks the player to go anywhere. A row
   of beacons makes a race or an escort.
3. **[Built as `[need:max_object]` and `[need:object:0]` (chapter 6).]** **A cap on objects:** `[need:object_max:<n>:<s>]` (at most n objects for s seconds), the
   opposite of chapter 4. It suits the replicating species: thin out an outbreak with bullets
   and recycling, which gives the gun a real job.
4. **[Partly: only bullet patterns cost energy (`[new_bullet:<n>]`); firing and cloning are free.]** **Energy as a currency:** firing, cloning and recycling cost energy (the score goes down).
   Dots become fuel rather than points, and the player chooses between shooting, cloning and
   saving.
5. **[Not built.]** **Lab tools as technology:** Draw ("synthesis", adding mass, the partner of recycling),
   Detect: Every ("scanner": names every creature, a help for finding morphs), the mask-radius
   slider ("sensor jammer" upgrade), body-axis view ("gyroscope"). Each unlocks by
   `[need:energy:…]`, as the second bullet type does.
6. **[Not built: radiation is still only the background.]** **Radiation that hurts:** in chapter 4 the gamma bursts are only told and shown (the orange
   background). Let them slowly wear down creatures in the open, but not those near big
   clusters, so the megastructure really shelters the ship. The bullets channel already has the
   wear-down mechanism to reuse.
7. **[Not built.]** **Hostile creatures:** morphs that home on the ship (Carlman's ghosts are the model), or
   creatures firing bullets of their own. Bigger work, but it turns "you are not alone" into a
   real threat.
8. **[Done another way: twins and trios are built by cloning and named as wholes, goals by `[need:discover_morph]` (chapters 3 and 5).]** **Docking:** `[need:merge:<morph>]`, reaching a pair morph by joining another creature
   (KH97WU *impetus+gyrans*, 5N7KKM *lithos+gyrans*); detection already names composite
   morphs. Cooperation as a mechanic, the opposite of shooting.
9. **[Not built.]** **The stats chart as a "spectrometer":** a late chapter sets a target in the lab's phase
   space (e.g. hold this spin rate). Niche, but true to the research side of the project.

### 8.3 Finishing chapter 4

*Status 2026-11-13: chapter 4 is finished more simply: `[need:object:50:min]`, then "Shield in
operation. We survived." The burst-and-survive payoff below is not built; it needs 8.2.6.*

The megastructure line (`[need:object:60:min]`) wants a payoff: build to 60 objects (mostly cloned
*lithos*, the still discs), then a gamma burst (`[clear_world:radiation]` with a stronger pulse,
using 8.2.6), then `[need:survive_time:30]`, then "The shield holds." Build, defend, survive.

### 8.4 Chapters for the unused species

*Status 2026-11-13: superseded. Every species is now in the story (table in 8.1), mostly as
discovery scenes in chapters 4 and 5 rather than as the chapters below. The ideas that need
beacons, hold or hostile creatures are still open.*

| Species | What it's known for | Goal idea | New mechanic |
|---|---|---|---|
| **VT049W** (replicating) | *gyrans 💥* splits into 3–4; *tachus* is fast | **Outbreak:** keep at most 6 objects for 60 s while gyrans keeps splitting; then become *tachus* and outrun it to a beacon. | object cap, beacons |
| **PJ1EF7** (multiplying) | *vagus 💥* multiplies | **Swarm:** survive 45 s in a crowd that keeps growing (as the chapter 5 notes say); Reset becomes a tactical retreat. | time goals |
| **XEH4YR** (exploding) | *auctus 💥*, *cor* | **Core:** find all three morphs in a cramped space (as the notes say), then hold *cor* steady while *auctus* blows up around it. | `[need:discover_morph]`, hold |
| **KH97WU** (dancing) | *impetus*, *impetus+gyrans* | **Dance:** dock with a gyrans to form *impetus+gyrans*; a cooperation chapter where shooting doesn't help. | docking |
| **KZ2BLD** (striped) | *virgatus* | **Corridor:** stripes suggest a long, oriented body; steer through a line of beacons, then become *virgatus*. | beacons, hold |
| **HAESRE** (zigzagging) | *zigzag*, *zigzag+lithos* | **Pilot's test:** collect energy along a path; zigzagging makes precise steering the challenge. | beacons, energy |
| **HGMLH8** (shapeshifting) | *metamorpha*, *paratardus*, *paratachus* | **Metamorphosis II:** become each of the three in order, against the clock. The notes say some of its morphs aren't recognised, so check its detection settings first. | hold, time goals |
| **G6G6CR** (Cilium) | the one non-Aquarium species; *munitus* = "fortified" | **First contact / finale:** a truly alien species; protect a *munitus* from hostile creatures, or recycle the ship into Cilium. | hostile creatures, synthesis |

Suggested order after chapter 4: HAESRE (pure steering), KH97WU (cooperation), VT049W
(containment), XEH4YR, PJ1EF7 (chaos), HGMLH8 (mastery), G6G6CR as the finale. Difficulty
rises and each chapter mostly adds one mechanic; most of the goals need only new conditions
(time goals, the object cap, beacons, docking), not new systems.

## 9. Performance: next steps (2026-11-11; still not built 2026-11-13, FFT mode now the default)

Measured with `?dev=1` (fps / steps per s / ms per step, a step's time including its analysis and
the GPU wait), after the per-frame step cap and FFT mode (`?fft=1`):

| | direct | FFT |
|---|---|---|
| laptop, sandbox | 34 / 30 / 22.3 | 50 / 31 / 18.3 |
| phone, sandbox | 17 / 17 / 51.0 | 23 / 23 / 37.4 |
| phone, chapter 1 | 24 / 24 / 38.4 | 42 / 30 / 24.9 |

FFT mode helps (18-35% per step), but far less than its ~10x fewer texture reads suggest: the
convolution was not the only cost. Plan:

1. **Measure first:** a breakdown in the `?dev=1` readout, per step: the sim pass on the GPU
   (finished: `gl.finish` around it, dev only), the readback wait (`readBlocks`), and the CPU's
   detection work (`Creatures.find`, tracking, classification). Run sandbox on the phone.
2. Then, by what it shows:
   - **`blocks.glsl`'s potential:** each step it convolves again, ~1,000 taps at each of the
     96 x 64 block centres (~31M reads), for detection. In FFT mode the potential is the mean of
     the U's the FFT already made (`P` textures): read it there instead.
   - **The readback after every step:** the CPU waits for the GPU each step (probably much of the
     laptop's 18 ms). Options: analyse once per frame, or read back a frame late (PBO + fence).
   - **CPU detection in a crowded world** (sandbox is slower than chapter 1 on both devices):
     classify fewer tracks per step, or detect less often.
   - **FFT mode's ~100 draw calls a step:** fold `fftmul` into the first inverse pass, and use a
     real-to-complex FFT (half the work).
   - **Dots stepping every other step:** ~4-5% (their 9% of the direct step).

## 10. Evolution in the lab (2026-11-13, not built)

**Why.** The ten species after *A. originale* (5N7KKM) were evolved from it by rounds of rule
mutation, selected by eye for metamorphosis, twins and trios, replication and so on. The
transgyrans tuning (MORPH.md, `output/transgyrans/tune/`) was one such round, done offline:
mutate, run from known morphs, measure, select. The aim here is to do it in the lab, where it can
be watched, and to let it *discover* species whose morphs differ in shape, behaviour and
combinations (twins, trios).

**How far the species are from originale** (`aquarium233.js`, 15 kernels each). Of 60 values
(`m`, `s`, `h`, `r` per kernel), 56-60 differ, by a median of 5-15 % (max 25-76 %). 1-7 kernels
also have changed ring weights `b` (e.g. `1,1/4` → `1,1/3,0`). Channel wiring `c` never changes,
and T stays 1. So a mutation must cover the kernel's shape (`r`, `b`), not only growth (`m`, `s`,
`h`). Steps of 1-3 % (as in the transgyrans tuning) fine-tune a species; reaching a new one took
many rounds of ~5 %.

### 10.1 The engine: a mosaic of experiments

The world view becomes a grid of tiles. Each tile is its own small torus with its own rule, and
the ship view shows one tile (click a tile to pick it).

- **Per-tile torus.** In `sim.glsl` a tap's cell wraps within its tile:
  `q = o + mod(p − o + off, tile)`, with o the tile's origin. Same cost as now. Check with a node
  emulation that one tile matches a single board of that size.
- **Per-tile rule.** Growth (`m`, `s`, `h`) comes from a small parameter texture, fetched once per
  cell, not per tap, so it is nearly free. Kernel shape (`r`, `b`) differs per tile too, so the
  tap table gets 5 rows *per tile*. The tap offsets are the union over tiles (out to the largest
  `r`), and a tile's weights are 0 beyond its own reach. Each tap is still read once. Normalisation
  is per tile.
- **Direct convolution only, at first.** FFT mode (the default since 2026-11-11) convolves the
  whole board with one kernel spectrum. Per tile it would need batched per-tile FFTs and a
  spectrum per tile. Leave that until the direct mosaic is measured.
- **Tile size.** At R = 18 a gyrans is ~40-70 cells across, so 128² holds a single morph and 192²
  a collision. The present 384×256 board gives 6 tiles of 128², and a larger lab board gives more.
  How many tiles a laptop and a phone can run is §9's step 1 question: measure the sim pass alone
  first.
- **Per-tile detection.** `blocks.glsl` already maps the board. `Creatures.find` runs on each
  tile's block sub-grid with the tile's own wrap, and each tile keeps its own tracks and stats.
- **Batch speed.** In experiment mode the lab can skip drawing the tiles it isn't showing, run as
  many steps per frame as fit, and detect only every few steps after the settling time.
- **Offline is still the bulk path.** Large sweeps (the transgyrans search ran ~1,500 runs of
  3,000 steps) remain faster in numpy with multiprocessing. The lab is for watching, for small
  sweeps, and for picking what to keep. Both write the same experiment JSON (10.4).

### 10.2 Experiments

Each experiment is a list of tiles, each with (rule, seed, steps), plus how to score the outcome.

1. **Mutation sweep (morph response).** Take species S and morph M. In each tile, mutate the rule:
   - systematic: one parameter over a grid of factors (e.g. ×0.97 … ×1.03), as `tune.py` did;
   - random: every value × (1 + σ·N(0,1)), with an occasional `b` change.

   Run from M and classify the outcome against M's own stats: unchanged (Δspeed, Δturn, Δmass
   within the detection widths), drifted (same shape, different speed or turn: "slightly changes
   direction or rotation"), became known morph X, became something unknown, split, died, exploded.
   The result is a map: parameter × factor, coloured by outcome. Reproducing the transgyrans
   result (kernel 3 `m` ×0.985 makes the crossing last) is the first test, since its answer is
   known.
2. **Collisions.** Two morphs A and B (B can be A) under one rule. Grid or random over the impact
   offset (perpendicular to A's heading), the relative angle, mirroring, and B's phase in its
   pulse. For phase, run B alone first in its own tile for k·Δ steps, then paste it, using
   `clone.js`'s rotate/mirror paste. Outcome: the objects at the end and their names. Of most
   interest: a new whole (twin or trio, named or not), a new single morph, or annihilation.
3. **Seeds** (what a tile starts from, for either of the above):
   - a morph of the species (the long history's pattern);
   - a morph of *another* species, transplanted (same R; it is often unstable and changes);
   - a random soup: a disc ~1.5 R across with noise per channel, as Lenia's random search did.
     Most soups die or explode, and the survivors show what the rule settles into on its own.

### 10.3 Measuring outcomes, and the trap of a mutated rule

Detection (`detect.js`) names morphs against the long histories, which were recorded under the
species' *unmutated* rule. A mutation of a few percent moves every morph's mass and speed out of
the 3 % widths, so everything reads "???". "Unknown" would then mean "drifted", not "new". The
fix is to rebuild the references per variant:

- **Step 1 under a new rule, always:** run each known morph of the parent alone (experiment 1 with
  no further mutation) and record its settled stats as that variant's reference: "gyrans under
  variant v". Morphs that die or change go in the variant's profile (10.4) too.
- **New morph:** a creature whose stats settle (low variance over a window, as Record expects)
  but that fits none of the variant's references. It becomes a candidate: its pattern (cropped,
  RLE), its stats piece, and a thumbnail. Candidates from different tiles are compared with each
  other (the same kernel density) so that duplicates merge. Wholes work the same way, using the
  twin saddle.
- **The player names it.** A candidate shows in the ship view. Record (Part / Whole) with a name
  adds it to the variant's long history, as now. Automatic finds stay "candidate #n" until
  checked by eye.

### 10.4 Variants, profiles and selection (the evolution loop)

- **A variant** is a parent code plus tweaks, in `tools/rule_tweaks.json`'s format (kernel, key,
  value, was, why), so a kept variant goes through `gen_aquarium.py` like the transgyrans fix.
  It gets a new 6-character code and a "parent" field. The lineage back to 5N7KKM is kept.
- **A profile** per variant, from its experiments:
  - how many distinct steady morphs (known + candidates);
  - the mix of still and moving morphs;
  - spontaneous morph changes per 1,000 steps (metamorphosis);
  - the object count's growth (replication);
  - the wholes found by collisions (twins, trios);
  - the share of seeds that survive (a filter: below ~50 % the variant is dropped, not ranked).
- **Selection**, in increasing automation:
  1. *Interactive* (Picbreeder-like, closest to how the species were found): the mosaic shows 6-12
     variants of the current parent, each with its profile as badges. You keep one or a few, and
     the next generation mutates those.
  2. *Novelty / MAP-Elites:* bin variants by two profile axes (e.g. morph count × replication
     rate, or still/moving × wholes found) and keep the best per bin. Running unattended, it fills
     a map of qualitatively different species instead of climbing one score.
- **Export:** variants, experiments (tiles, seeds, outcomes) and candidates as JSON. They can be
  imported into the offline tools and back.

### 10.5 Milestones (suggested order)

- [ ] §9 step 1 first: the per-step breakdown, which shows how many tiles a laptop and a phone can run.
- [ ] E1 Mosaic engine: per-tile torus and per-tile rule in `sim.glsl` (direct mode), per-tile
  detection and tracks; an Experiment mode in the lab; click a tile to show it in the ship view.
  Checked by node emulation (one tile = a single board).
- [ ] E2 Mutation sweep with outcome classes and the parameter × factor map; reproduces the
  transgyrans result.
- [ ] E3 Per-variant references (10.3) and new-morph candidates, deduplicated, with Record to name
  them.
- [ ] E4 Collisions (offset, angle, mirror, phase) and wholes as outcomes.
- [ ] E5 Seeds from other species and random soups.
- [ ] E6 Variants with lineage and profiles; interactive selection; export to `rule_tweaks.json`
  form and `gen_aquarium.py`.
- [ ] E7 Unattended novelty / MAP-Elites runs.
- [ ] Later: per-tile FFT, so the mosaic can use FFT mode too.

## 11. Challenge mode (2026-11-15, not built)

### 11.1 The brief (paraphrased)

A third mode beside story and sandbox: a list of challenges to choose from at the start, each
marked with a green tick once completed. Some challenges allow any species, others name one. A
run starts from any discovered morph. Morph data syncs with the central inventory, read and
write, as in sandbox mode, but the energy collected does not sync (a run's energy is its own, as
in a story run). Challenges:

- **Courier** (any species).
  - A sequence of beacons, each a pattern in the dots channel (e.g. *Heptadentium scintillans*).
    Eating one puts the next one somewhere else.
  - Goal: after the beacons (zero or more), park in a small goal area as a still morph (e.g.
    *lithos*).
  - The world is much larger than the board. The world view scrolls with the player at its
    centre, and the star field is fixed to the world. Objects are pasted and removed on demand,
    just outside the screen.
  - Variants: a time limit (e.g. 30 s); an energy cost (e.g. 1 per mask press, with dots and
    beacons refilling it); erasing fields; masking fields; hostile morphs.
- **Docking** (named species). Dock the ship on one still morph or a group of them, e.g. as a
  *gyrans* circling a *lithos*. Goal: detected as a fuzzy *circumgyrans*. Candidates: 5N7KKM or
  HAESRE *circumgyrans* (gyrans + lithos), perhaps KH97WU *gyrans + impetus*.
- **Black hole** (any species). The world has erasing fields, as walls or a black hole. Needs
  experiments.
- **Nebula** (any species). The world has masking fields, as nebulae, guardrails or a funnel.
  Needs experiments.
- **Fight** (any species). Against one or a group of hostile morphs ("bots"). Still or rotating
  bots fire bullets at random; moving bots head for the player, bumping and/or firing. Goal:
  survive the attacks, or remove every bot with bullets. A boss fight: one exceptionally hard bot.
  Needs experiments.

### 11.2 What exists to build on

- **The script language.** A challenge is a script block in a `challenges.txt`, using the
  `[on_choose_species]` block mechanism, `[need:…]` conditions with progress, `[dot]`,
  `[min_morph]`, `[put_player]`, `[bullet]`, `[clear_world]`. Most of a challenge is script,
  not code.
- **Sandbox's inventory handling.** Sandbox mode reads and writes the central inventory and
  syncs energy. Challenge mode is the same with energy kept per run, as a story run does.
- **Static fields come almost free.** `sim.glsl` already hides the bullets channel from the
  convolution, per tap (`uBulletMode == 2`, mask), or erases mass under it (mode 1). A static
  *field* texture (erase strength, mask strength, per channel) uses the same lines, read once
  per tap or per cell. This is §5's level mask texture, still open in §7.
- **Bullets** already erase whatever creature they hit, the ship included, so bots can reuse
  `Bullets.fire` from another position and angle.
- **Detection** already names fuzzy morphs and wholes (twins and trios), which is what docking
  needs.

### 11.3 Comments and open questions

1. **The large scrolling world is the hard part, and it isn't needed to try the courier
   challenge.** A bigger board costs GPU time in proportion to its area, and the 384×256 board
   already takes 18-37 ms a step (§9). The brief's on-demand idea fits the torus well: keep the
   board, keep a world offset for the player, and map world positions to board positions modulo
   the board size. Then:
   - the view draws the board centred on the player (an offset in `draw.glsl`), with the stars at
     world positions;
   - a band on the far side of the torus, opposite the player, is the *seam*. Each step it is
     cleared of creatures, dots and bullets, and things that come into range are pasted just
     inside it;
   - the seam must stay off-screen and at least a kernel reach (R) from anything visible,
     because cells within R of it feel the wrap. So the view shows a little less than the board:
     the board minus ~1.5 R on each side, slightly zoomed. The world view's canvas size doesn't
     change.

   Suggestion: build the courier on the present torus first, with beacons and a goal area. Judge
   whether it is fun, and only then add scrolling as a variant ("deep space").
2. **Beacons.** *Heptadentium scintillans* doesn't exist yet. A beacon in the dots channel has to
   be a soliton of the dots' rule (Carlman's, μ 0.38, σ 0.07), and a 7-fold twinkling one has to
   be found first: a lab or numpy search. Fallback: a ring of 7 ordinary dots, which looks
   similar and costs nothing. A beacon should also be marked on the screen (a pointer at the
   view's edge when it is out of sight), since steering is slow.
3. **Parking** is §8.2's `hold` condition with a place added:
   `[need:park:<morph>:<x>,<y>,<r>:<s>]`, meaning inside the circle, named that morph, for s
   seconds. It needs a goal-area overlay. *Lithos* is the obvious morph; *tardus* is a harder
   variant.
4. **Challenges need failure, which the story doesn't have.** The story only waits for success,
   and Dissolved offers Retry. Challenges need: a timer shown as a countdown (`[time_limit:<s>]`),
   a ship lost (as now), and the energy running out (with mask costs). Each ends in "Failed" with
   Retry, and in success with the tick. Store the best time beside the tick, so a done challenge
   can still be improved.
5. **Energy cost of masking.** Charge per second of mask held rather than per press, so a long
   hold isn't free. The energy shows already; it needs a "can't mask" state at 0.
6. **Fields** (black hole, nebula). Experiments for the lab, before the challenges are designed:
   - A **mask field** should repel by agnosiophobia (§6.4's observation): a guardrail or funnel
     that steers without killing. Test full and partial strength (0.3-1) and per-channel fields:
     a red-only nebula may push differently from a green one.
   - An **erase field** kills on contact. A radial gradient (weak at the edge, strong in the
     middle) makes a black hole that drains slowly at its edge. That is a readable danger, and the
     same thing as §8.2.6's radiation. Erasing doesn't *pull*, so a black hole needs another
     mechanism if it should attract: e.g. a mask field on the side away from it, pushing the ship
     inward. Worth testing.
   - The lab needs a way to paint fields (a brush like Draw) to run these experiments.
7. **Docking.** *Circumgyrans* has to be recorded first, as a fuzzy whole for 5N7KKM and HAESRE
   (perhaps KH97WU), with a twin saddle that keeps the pair as one. Then the goal is
   `[need:hold:circumgyrans:<s>]`, so a single pass doesn't count. Check in the lab that it is a
   real attractor and not a passing state, as transgyrans was before its rule tweak (MORPH.md).
8. **Fights.** All creatures share one rule, so a bot is a morph of the *player's* species.
   That's fair, but it means bumping can merge the two into a twin instead of hurting the ship.
   That is itself a mechanic, but not "damage". Suggested order:
   - *turrets*: still or rotating morphs (lithos, gyrans) firing at random, reusing
     `Bullets.fire` from their position. Cheap, and immediately a dodging challenge;
   - *hunters*: moving morphs steered by the game with masks, using a simple feedback controller
     (mask the side away from the target, since a mask pushes). This is a small version of
     §6.4's agent and wants the mask maps' turn zones, so §6.2 pays off here;
   - *boss*: a large morph or a whole (trias) with a fast gun and more mass to wear down.
   Bots need bullets marked by owner, so their shots don't count as the player's and the
   friendly-fire rules are clear. Bot names and their health (mass) show in the world view.
9. **"Any species" challenges vary in difficulty with the species.** Fix the species where a
   challenge depends on a morph (parking needs a still morph; most species have *lithos*, but
   not XEH4YR or PJ1EF7's set). Otherwise list the allowed species per challenge.
10. **The challenge list** is a book like the species book: one slot per challenge, with a tick,
    best time, species, and a short description; locked challenges are greyed (e.g. docking
    needs the morphs discovered).

### 11.4 Decisions (the user's review, 2026-11-15)

These replace the comments in 11.3 where they differ.

- **The script vocabulary** (sample in `challenge.txt`):
  - `[challenge:<name>:<n>]`: the challenge's printed name, stage n;
  - `[allow_species:all|<ids>]`: which species the player may choose;
  - `[clear_world:scroll]`: reserved for the scrolling world;
  - `[dot:<pattern>:<how>]`: put a named pattern (e.g. Heptadentium) in the dots channel instead
    of the round dots, once, with no respawn when eaten;
  - `[need:clear_dot]`: the dots channel is empty (mass below a small epsilon);
  - `[target:<how>]`: a target area, a small square the size of HAESRE *antigyrans*, with a dotted
    orange border;
  - `[need:stay_target:<s>]`: the ship stays inside the target for s seconds, counted as
    `[need:survive_time]` is.
- **Parking isn't checked by morph.** `[need:stay_target]` is enough: any morph that can stay in
  the narrow target long enough will do.
- **Heptadentium** comes from `Lenia/Python/animals.json` (2D7s *H. scintillans*: R 18, T 10, b
  1,5/6, m 0.22, s 0.025, exponential kernel, Gaussian growth). It doesn't have to share the
  dots channel with the round dots: the channel takes one rule at a time, the dots' or the
  pattern's.
- **Docking:** 5N7KKM and HAESRE *circumgyrans* are already recorded as fuzzy patterns.
- **Hijacking is allowed.** Merging with a bot into a twin, or dying in a bump while the focus
  passes to the bot (this already happens in story mode during `[need:survive_time]`), is
  emergent play. The player only has to survive, not to stay the same creature.
- **Bots' guns lock the bullets channel's colour** (e.g. G), so the player can't switch to a
  colour that affects them less. Options to try: only bots fire (Fire off), both fire, or they
  fire in separate channels (costs performance).
- **Hunter steering by a recorded "sweet spot"**, not a learned agent. In the lab, record per
  morph a mask offset and a hold time that turn it, short enough not to change its morph. When the
  bot wants to turn towards the player, the game presses a mask there for that time. It doesn't
  have to be accurate, since the bot keeps correcting.
- **"All species" means the player decides,** even with a species that supposedly can't do it.
  Emergent play is welcome.
- **Order:** C1 framework, C7 Heptadentium, C2 courier, C4 docking, C5 fight, C6 scrolling,
  C3 fields (needs the most lab work).

### 11.5 Milestones (in the order of 11.4)

- [x] C1 Framework (2026-11-15, node syntax checks only, not yet played): `?mode=challenge`,
  `challenge.txt`, the challenge book (ticks, best times), species chosen per `[allow_species]`,
  inventory synced without energy, completion recorded. Then (the user's review, 2026-11-15): titles
  *Courier #1*, card CHALLENGE / name / STAGE n, SUCCESS only at a challenge's end, Hold to clear
  completion & best time, fail states ("Failed": ship dissolved, `[time_limit]`, energy out
  under `[mask_cost]`, with `[energy]`), the retry book ("Choose a morph to retry", from Morphs,
  Retry, Choose; no cloning), and in every mode a retry puts the player's last morph in place of
  the script's. Still open: locks.
- [x] C7 *Heptadentium scintillans* in the dots channel under its own rule (ring weights and the
  exponential kernel core in the dots' taps), `[dot:Heptadentium:<how>]` (2026-11-15; numpy: stays
  put, mass within 1 % at R 18; not yet seen on the GPU).
- [ ] C2 Courier on the present torus: Heptadentium beacons (`[dot:<pattern>]`,
  `[need:clear_dot]`), `[target]` and `[need:stay_target]` (built 2026-11-15, not yet played);
  variants: `[time_limit]`, `[mask_cost]` (built, not yet used in `challenge.txt`).
- [x] C4 Docking (2026-11-15, node checks only, not yet played): `[min_morph:<n>:<morph>=<k>]` puts
  groups (k copies packed with a 3-cell gap, random angles), `[need:stay_morph:circumgyrans:<s>]`.
- [ ] C5 Fight: turrets (gun colour locked), then hunters (lab-recorded sweet spot and hold time), then a boss.
  Built (2026-11-15, node checks only, not yet played): `[put_hostile:<n>:<morphs>:<how>:<bullet>]`
  turrets (bots of the player's species, aimed shots every 1.5-4 s) firing in a separate hostile
  bullets field (white, masking G, bots included: erasing bullets, and masking all three, barely hurt), bots
  tracked by place (what a bot turns into, splits into or merges with stays hostile; hijacked: the
  player's), a red circle and the morph's name on each (twins and trios as one, and firing as one,
  as the player's twin ship), `[need:clear_hostile]`; tokens may span lines. The boss for now: a
  *trias* with a strong bullet (Heptapteryx), splitting into three bots. Still open: hunters, health shown.
- [x] C6 Scrolling world (2026-11-16, node checks only, not yet played): `[clear_world:scroll]`, the
  seam scheme of 11.3.1 without the zoom (the user's choice: the view shows the whole board, following
  the ship; the seam along its edges, a 6-cell band cleared each frame); an object whose body touches
  an edge is taken away and [min_morph] puts its replacement on the opposite edge (round dots too);
  `extended` placement for beacons and the target (kept by world point off the board, yellow arrows).
- [ ] C3 Field texture in `sim.glsl` (mask and erase, per channel, with strength), a field brush
  in the lab, and the experiments of 11.3.6; then the nebula and black-hole challenges.
