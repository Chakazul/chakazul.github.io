# Coolship

Steer a Lenia creature by hiding parts of it from itself. The ship is an Aquarium, a
3-channel Lenia creature from the aquarium collection (by default *5N7KKM limus*, the original
Aquarium), run with LeniaF.py's rule. The main control is **cognitive
masking** (Cool et al. 2026): while you hold the pointer on the magnified creature, that spot of
one colour channel drops out of every kernel that reads the channel,
U = K∗(A·(1−B)) / K∗(1−B). The same circle can instead erase or insert mass, for comparison. See
[PLAN.md](PLAN.md) for the design and the later phases.

**Run:** serve this folder over http and open `index.html`. The page fetches its shaders, so
`file://` won't work.

```
cd Coolship && python3 -m http.server 8000     # then open http://localhost:8000/
```

## Controls

- **Press and hold on the magnified creature** (right panel) to apply the selected tool to a
  circle there. The circle stays at the same place relative to the creature's centre of mass
  until you let go. You can drag while holding. With a mouse, the tool's circle follows the
  crosshair over the magnified view; while pressed it is filled with the channel colour.
- **Tool patches.** The top row picks the channel, **R G B** (keys `1` `2` `3`); the bottom row,
  drawn in that channel's colour, picks the action (keys `4` `5` `6`). G mask is the default: it
  steers and shape-shifts best.
  - **Mask**: a **cognitive mask**, drawn as a dashed outline. The cells are hidden from the
    kernels, but no mass is touched.
  - **Erase**: drawn as a solid outline. The channel inside the circle is set to 0 every step
    while held, like Carlman's interventions. It gives a direct comparison with masking.
  - **Insert**: drawn as a dotted outline. The channel inside the circle is set to 1 every step
    while held.
- **Yellow arrow:** the velocity. Its length is the speed, and it is hidden below 0.02 R per unit
  time (the glider does about 0.072 R, i.e. 1.3 cells at R = 18).
- **Heading up** toggle (`E`, after the Speed slider): turns the ship view's rotation on or off:
  - **board up**: no rotation (the default).
  - **heading up**: the view follows the arrow's own smoothed velocity, so the arrow always
    points straight up. Near stationary, when the arrow is hidden, it holds the last heading.
  - A third mode, **body axis up**, is hidden from the button and only reachable with
    `?headup=axis`: the view follows the line from the B channel's CoM to the G channel's, which
    exists at any speed but is not exactly the direction of travel (about 56° off for *gyrans*
    under the legacy logic).
- Hovering a tool patch shows what it does. Each click or shortcut briefly highlights the controls it used (Space while playing lights both
  Pause and Step).
- **Pause** (`Enter`) · **Step** one step (`Space`) · **Mask radius** in cells.
- **Reset** (`Z`): clear the board and put the pattern in the centre.
- **Add** (`X`): add the pattern, rotated by a random angle, at a random spot, the emptiest of 24 candidates, without clearing
  the board. If nothing is focused, it focuses the new creature; after *Dissolved*, it restarts
  the run.
- **Speed** (`Q` faster, `W` slower): 7 slots, 1/8, 1/4, 1/2, 1 (the default), 2, 4, 8, at a quarter step per frame each:
  speed 1 runs one step every 4 frames and speed 8 runs two steps per frame. Speed doesn't change
  the physics, only how many steps you see per second. The time step dt = 1/T is a separate, fixed
  property of the rule (T = 2, see below): it sets how far the simulation advances per step, and
  changing it changes the dynamics.
- **Zoom**: screen px per cell in the ship view (2× to 16×; the default is the slot nearest
  4 R across, 6× = 64 cells at R = 18).
- **Species** and **Morph** pickers. A species is a rule code, named in the picker as Aquarium
  (5N7KKM), Aquarium replicating (VT049W), Aquarium spinning (Z18A9R), Aquarium dancing (KH97WU),
  Aquarium drifting (GDNQYX) and, last, Cilium (G6G6CR). Choosing one resets the
  world with its first morph. A morph (limus, …) only changes what Reset and Add (`X`) put in, and
  keeps the world. Picking a morph in the list (even the current one) also adds it, like Add;
  ▲ ▼ or `C` `V` only step through the morphs. 42 entries of the aquarium collection from 5N7KKM to
  GDNQYX, in file order. 5N7KKM is the original Aquarium, found by globally
  randomising parameters; the other codes were evolved from it.
  - **5N7KKM limus**, the first entry, is the default.
  - **VT049W rectus** is the steady glider the Phase 1 checks used.
  - **VT049W gyrans (replicating)** is the original self-replicating Aquarium. The other VT049W
    forms were stabilised from it. At R = 18 it divides into 3–4 creatures of about its own mass
    within ~400 steps, and the population then holds at 3–4 (checked to 4,000 steps).
  - **Z18A9R trio, lithos and tardus** were found in Coolship itself (at R = 18). They come after
    the Z18A9R forms and were rebuilt with `tools/from_screenshot.py` from the screenshots in
    `assets/screenshots/`. All three use a brightness correction of 1.8; the size is 1.0 for trio
    and lithos and 1.05 for tardus, with shape similarity 0.96, 0.93 and 0.92.
    - **trio:** three lobes, drifting slowly.
    - **lithos:** a round, still disc. Without the correction it died or became a spinning
      crescent.
    - **tardus:** a round form with a green and red core.
    Extra patterns live in `tools/extra_patterns.json`.
  - **To add another pattern from a screenshot** of the world view (default zoom is fine), run
    `tools/from_screenshot.py`. It needs numpy, scipy and Pillow, e.g. the complexity venv:
    ```
    ~/github/complexity/.venv/bin/python tools/from_screenshot.py shot.png --code Z18A9R --name trio --add
    python3 tools/gen_aquarium.py
    ```
    It reads the cell grid off the image. Screenshots look brighter than the cell values (about
    v^(1/1.8)), so it tries several brightness corrections (`--gammas`) and seed scales
    (`--scales`), settling each pair under the code's rule for 300 steps.
    - **A pair passes** if, over 500 more steps, its mass and number of creatures hold.
    - **The check is short on purpose:** a form can be unstable, or turn into another form in the
      long run, and still be worth having as a pattern.
    - **Of the passes, it keeps the one most like the screenshot:** the correlation of per-channel
      radial profiles, which ignores rotation.
    It writes `seed.png`, `scales.png` and `verify.png` to look at; see the script's docstring for
    the options. A run takes about 2.5 minutes.
  - Left out: *G6G6CR metamorpha beaked*, *KH97WU dot*, and *VT049W gyrotachus*, which is a
    different species (its kernel 0 has μ = 0.302, not VT049W's 0.272). The three unnamed Z18A9R entries are
    *gyrus*, *para* and *antipara*; VT049W's two *para* entries are shown as *paralimus* and *para*.
  - Every pattern runs at the game's R (18 by default), its cells resampled from the entry's own
    radius (20–30).
  - Two legacy test patterns from before the switch to LeniaF's rule are hidden, reachable only
    with `?form=fission` or `?form=fission02` (see `species.js`).
- **Several creatures.** Click a creature in the world view (left) to make it the focus that the
  zoom view, the arrow and the tools follow. The focus is outlined by the zoom frame.
  - If the focused creature splits, the focus stays on one half: the one where it was heading.
  - If the focused creature dies, nothing is focused. The zoom view goes blank until you click
    another one.
  - The run stops as *Dissolved* when no creature is left. There is no "exploded" stop: a board
    full of mass just keeps running.

## URL parameters

- `?form=5N7KKM-limus`: the pattern to start with (the default), by id: code and form joined by `-`, e.g.
  `5N7KKM-limus`, `VT049W-gyrans-replicating`, `5N7KKM-lithos+gyrans`. `aquarium` and
  `replicating` are short for the two VT049W ones. Hidden: `fission`, `fission02`.
- `?R=18`: kernel radius in cells (6–30). The board stays 384×256, so a larger R makes the
  creature larger and the sim slower. Below 16, *gyrans R* stops replicating and *lithos* moves.
- `?dt=0.5`: time step, overriding the rule's own (1/T = 0.5).
- `?speed=1`: the Speed slot at start, snapped to the nearest one (1/8–8; ¼ step per frame per unit).
- `?zoom=8`: ship-view px per cell at start, snapped to the nearest slot (2–16).
- `?radius=4`: mask radius at start, in cells (default 2.5 · R / 12, rounded to 0.5).
- `?headup=vel|axis`: start in heading-up or body-axis-up view (`1` = `vel`).
- `?debug=1`: show GPU ms per step and the tap count in the readout.

## How it works

- `aquarium233.js`: rules and cells from Bert Chan's Lenia "aquarium" collection
  (`Lenia/Python/found/233.json`, lines 2–51), generated by `tools/gen_aquarium.py` (rerun it if
  the source changes; don't edit the output by hand).
- `species.js`: `PATTERNS`, the picker's list (the generated patterns plus the hidden legacy
  ones), and `patternRule()`.
- `rule.js`: RLE decoding, nearest-neighbour resampling to the run radius, LeniaF's `free_h`, and
  the packed kernel tables, with no WebGL dependency.
- `glsim.js` + `shaders/`: the WebGL2 engine.
  - Each step runs `sim.glsl` (convolution, growth, and the mask or erase circle, whose centre the
    CPU passes in).
  - Once per frame, `blocks.glsl` sums 4×4-cell blocks into a 96×64 map of mass and moments,
    plus the potential at each block's centre (one extra convolution at 1/16 of the cells), and
    the map is read back.
  - In body-axis view, `reduce.glsl` and `com.glsl` also compute per-channel CoMs in a window
    around the focus.
  - `draw.glsl` renders the world and the zoom view.
- `creatures.js`: finds creatures with the method of `~/github/complexity`
  (`measures/individuality.py`). `blocks.glsl` also samples the potential K∗A (the mean over
  kernels, unmasked) at each block's centre.
  - **Peaks:** each peak of the potential claims its catchment basin above 0.05.
  - **Merging:** neighbouring basins merge when the saddle between them is at least 0.8× the
    lower peak.
  - **Size:** basins smaller than a disk of radius R/4 are dropped.
  This separates touching creatures, such as the parts of a composite seed or a replicator just
  after it divides, which connected components of the mass counted as one. Each creature's CoM is
  rebuilt exactly from per-block moments. The module also picks a creature for clicks and for
  frame-to-frame tracking.
- `app.js`: layout, input, focus tracking, heading logic, the overlay and the readout.

The rule follows **LeniaF.py**, the program that reads 233.json (`LeniaF.py -c 3 -k 3`):
- polynomial kernel cores (4f(1−f))⁴ over the b rings, each kernel normalised to 1;
- polynomial growth 2·max(0, 1 − (u−m)²/(9s²))⁴ − 1;
- `free_h`: each h becomes round(h / Σh over the same destination channel, 3), then h and T are
  both doubled (T: 1 → 2);
- A_c += (1/T) · Σ h′_k G_k, clipped to [0, 1].

The engine's first version used LeniaNDKC-style exponential kernels and Gaussian growth at
dt = 0.2. That was not what 233.json was found with. Under that rule *gyrans R* didn't replicate
and *lithos* moved. The Phase 1 numbers below are from that version.

## Fidelity checks (Phase 1, legacy rule)

These were run off-line, with a numpy port of LeniaNDKC (FFT convolution) and a node emulation
of `sim.glsl`'s packed maths. The scripts were throwaway and are not in the repo.

- **The packed shader maths matches numpy** to 7·10⁻⁸ after 20 steps, unmasked and with a G
  mask held at a fixed spot (checked against an independent numpy implementation of Eq. 2).
- **The glider is stable.** *rectus* at R = 12 and dt = 0.2 glides for 5,000 steps at about 0.49
  cells per unit time, with mass steady at about 336 and no fission. It also glides at native
  R = 25 and at R = 15, at the same speed measured in units of R. *gyrans* survives at R = 12 as a
  circling form. *lithos*, which should be stationary, ends up moving at R = 12, so it is not
  offered.
- **Masking steers it.** Masks of radius 2.5 cells were held for 50 steps on a 5×5 grid of
  offsets around the CoM, on each channel. The outcomes range from no effect, through turns of
  ±70–160° and slowing down, to death: holding G just ahead of the centre is lethal. R has the
  gentlest map. This is Phase 1's numerical proof that the ship can be controlled; the manual
  check in the browser is PLAN.md §4.4.

## Checks for Phase 1b (legacy rule)

- **No other form of the same rule self-replicates** under the normalised update: *tachus*,
  both *para* forms, *para+sequa*, *rectus+gyrans* (a pair that drifts apart and rejoins) and
  *gyrotachus*, 3,000 steps each. lenia-CE's VT049W doesn't either, at R 12 or 16 and dt 0.25
  or 0.5.
- **Without normalisation, the Aquarium switches sharply** at dt = 0.5. It glides unchanged for
  4,000 steps at dt ≤ 0.4975 and divides at 0.5, so there is no slower-replicating setting
  between the two.
- **A replicator at dt = 0.2** (normalised update, 160² torus, 4,000–6,000 steps). Changing
  one kernel's μ by ±4% or σ by ±6% (60 variants): none divides. Changing μ by ±10%, σ by ±12%,
  or h by 0.8× or 1.25× (90 variants): only kernel 5's μ × 0.90 divides, and it fills the board
  by step 1,600. A fine scan of μ × 0.905–0.955 shows every value dividing, but always by
  swelling to 3–4× the mass first, never by clean Aquarium-sized copies. μ × 0.955 is the
  slowest, hence the test pattern.
- **Focus tracking** was checked by replaying the numpy fission run as block maps, every 3 steps,
  through `creatures.js` in node. The focus keeps one half at every split, follows into merges,
  and is lost once, inside the colliding blob soup.

## Checks for the LeniaF rule

These use a numpy port of LeniaF's update for old-format entries: model from `params[0]`,
`free_h`, polynomial kernel and growth.
- **All 45 entries from 5N7KKM to GDNQYX**, at native R and seed, 2,000 steps on a 192² torus.
  - ***gyrans R* is the only self-replicator.** It splits from step 200, and its mass goes from
    1,618 to about 15,000.
  - **Every other single-creature form stays one stable creature.**
  - **The names fit the behaviour:** every *lithos* is stationary, *tachus* (0.092 R/t) is faster
    than *rectus* (0.072), and the *tardus* forms are slow.
  - The composite seeds (*lithos+gyrans*, *dual*, and so on) come apart into their parts.
- **Radius.** At R = 12, *rectus* is faithful (0.072 R/t), but *gyrans R* becomes a plain glider
  and *lithos* moves. From R = 16 up, *lithos* is still and *gyrans R* replicates:
  - R = 16: slowly, then it fills a third of the board by about step 3,200;
  - R = 18: 3–4 creatures, steady to 4,000 steps;
  - R = 20 or 25: it fills about a third of the board.
  R = 18 is the default.
- **The packed shader maths matches the LeniaF port** to 2·10⁻⁷ after 20 steps at R = 18, for
  *rectus* unmasked and for *gyrans R* with a G mask.

- **The 43 picker patterns from 233.json at R = 18** on the game's 384×256 board, 2,000 steps, against the
  native-R runs above. Most keep their speed, mass and number of parts. The exceptions:
  - *gyrans (replicating)* still divides, to 3 creatures by step 2,000, at a lower total mass
    (4.6× instead of 9.4× the seed);
  - the *para* forms and *GDNQYX paragyrans* are 20–80% faster;
  - G6G6CR (native R = 20) changes most. Its *lithos* now stays still (0.036 → 0.005 R/t).
  - The large composites (*5N7KKM lithos+3gyrans*, *lithos+4gyrans*) end differently.
- **Detection** (node, on numpy block maps of 8 runs every 50 steps):
  - single creatures stay one creature at every frame;
  - *2lithos+gyrans* is found as 3 parts (before: 1) and *para+sequa* as 2 (before: 1);
  - *gyrans (replicating)* gives 3–4 creatures, the same as before or one more right after a
    division.
  A find takes 0.4 ms on desktop node.

Not yet checked on a device: GPU speed on phones (use `?debug=1`) and mobile shader precision.
