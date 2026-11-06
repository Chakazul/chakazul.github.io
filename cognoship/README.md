# Cognoship

Steer a Lenia creature by hiding parts of it from itself. The ship is an Aquarium, a
3-channel Lenia creature from the aquarium collection (by default *5N7KKM limus*, the original
Aquarium), run with LeniaF.py's rule. The main control is **cognitive
masking** (Cool et al. 2026): while you hold the pointer on the magnified creature, that spot of
one colour channel drops out of every kernel that reads the channel,
U = K∗(A·(1−B)) / K∗(1−B). The same circle can instead erase or draw mass, for comparison. See
[PLAN.md](PLAN.md) for the design and the later phases.

Two pages share one engine:
- **`index.html`, the game.** For now a trimmed lab: steer the ship with a smaller set of
  controls. Game logic (goals, levels, story mode, more channels and objects) will be added here.
- **`lab.html`, the lab.** The experiment version, with every control: stats, Record /
  Download, morph detection, tuning. The controls below are the lab's.

**Run:** serve this folder over http and open either page. They fetch their shaders, so
`file://` won't work.

```
cd Cognoship && python3 -m http.server 8000     # game: http://localhost:8000/   lab: http://localhost:8000/lab.html
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
  - **Draw**: drawn as a dotted outline. The channel inside the circle is set to 1 every step
    while held.
- **Yellow arrow:** the velocity, close to the actual one (a ~5-step average, corrected for the
  lag that averaging causes on a turning ship). Its length is the speed, and it is hidden below
  0.005 R per unit time (the glider does about 0.072 R, i.e. 1.3 cells at R = 18; the threshold is
  low enough for the slow *tardus* forms).
- **Heading up** toggle (`E`, after the Speed slider): turns the ship view's rotation on or off:
  - **board up**: no rotation (the default).
  - **heading up**: the view follows the direction of travel averaged over ~30 steps (steady
    through pulses and detection jitter, e.g. the *paragyrans* and *para* forms), predicted
    forward by the measured turn rate so a circling ship is followed with no lag. The arrow therefore points about up, swinging a little with each pulse.
    Near stationary, when the arrow is hidden, the view holds the last heading.
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
  original (5N7KKM), Aquarium replicating (VT049W), Aquarium spinning (Z18A9R), Aquarium dancing (KH97WU),
  Aquarium drifting (GDNQYX), Aquarium shapeshifting (HGMLH8), Aquarium multiplying (PJ1EF7),
  Aquarium striped (KZ2BLD), Aquarium exploding (XEH4YR), Aquarium zigzagging (HAESRE) and, last, Cilium (G6G6CR). Choosing one resets the
  world with its first morph. A morph (limus, …) only changes what Reset and Add (`X`) put in, and
  keeps the world. Picking a morph in the list (even the current one) also adds it, like Add;
  ▲ ▼ or `C` `V` only step through the morphs; `Shift`+`C` / `Shift`+`V` step to the
  previous / next species (and reset, as picking it). 62 entries of the aquarium collection from 5N7KKM to
  HAESRE, in file order (the last five codes added 2026-11-05, not yet checked in the browser or
  numpy; KZ2BLD's unnamed second morph is named *virgatus* here, its *gyrans tachus* is *gyrans*,
  and a repeated name gets a number: XEH4YR *auctus 2*, HAESRE *zigzag 2*). 5N7KKM is the original Aquarium, found by globally
  randomising parameters; the other codes were evolved from it.
  - **5N7KKM limus**, the first entry, is the default.
  - **VT049W rectus** is the steady glider the Phase 1 checks used.
  - The self-replicating morphs carry 💥 after their name (in the file *VT049W gyrans R* and
    *PJ1EF7 vagus R*, and XEH4YR *auctus*); their ids keep `-replicating`, and the long history
    records their name without it.
  - **VT049W gyrans 💥** is the original self-replicating Aquarium. The other VT049W
    forms were stabilised from it. At R = 18 it divides into 3–4 creatures of about its own mass
    within ~400 steps, and the population then holds at 3–4 (checked to 4,000 steps).
  - **Z18A9R trio, lithos and tardus** were found in Cognoship itself (at R = 18). They come after
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
    *gyrus*, *paragyrans* and *antigyrans*; VT049W's two *para* entries are shown as *paralimus* and
    *paragyrans*. Also renamed: G6G6CR *fortified* → *munitus*, KH97WU *antipara* → *antigyrans* and
    *dual* → *impetus+gyrans*, HAESRE *anti* → *antigyrans* and *para* → *paragyrans*.
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
- **Detect** (yellow when not None; after ▲ ▼; also shown in fullscreen, where the Stars button is hidden) cycles
  None, Focused, Every (the default): it names a creature's morph by comparing its stats with the
  species' long history (see *Morph detection* under Stats mode).
  - **Focused:** a name box at the top of the ship view, under the channel patches: the focus's
    morph, or **???** when it is unknown.
  - **Every:** also every creature in the world view: its circle drawn thick in its morph's mean
    hue (in place of the focus's yellow) and the morph's name above it; nothing when unknown.

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
- `statschart.js`: the stats mode's histories and phase-space chart (2D canvas). `assets/stats/<species>.json`: saved long histories.
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
- `core.js`, shared by both pages: layout, canvas input, focus tracking, heading logic, the
  stats histories, morph detection and the overlay. A page reads and sets its state through the
  `Core` object, sets `Core.hooks` (change, record, pattern, species, statsLoaded, frame) and
  `Core.opt` (what the overlay shows), then calls `Core.start()`.
- `page.js`, shared: DOM helpers, theme, fullscreen, music (the Stars button), the error box.
- `lab.js`: the lab's controls, readout, modes, stats chart, Record / Download, tuning, keys.
- `game.js`: the game's controls, HUD and keys (and, later, the game logic).
- `style.css`, shared by both pages.

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

## Music and stars

Background music, royalty free, looped as a playlist with 10 s of silence between tracks, at half
volume, from `assets/music/`: *Long Way* by Rexlambo (`rexlambo-long-way.mp3`), *A Few Jumps Away*
by Arthur Vyncke (`arthur-vyncke-a-few-jumps-away.mp3`), *Out of the Skies, Under the Earth* by
Chris Zabriskie (`chris-zabriskie-out-of-the-skies-under-the-earth.mp3`) and *Billions and
Billions* by Stellardrone (`stellardrone-billions-and-billions.mp3`). In the lab, the **Stars**
button at the end of the second row turns the star field and the music on or off together; in
the game the star field is always on and the **Music** button only plays or pauses the music
(each page load starts with the lab's off and the game's on; the choice is not remembered). Browsers only start
audio after the first click, tap or key.

The **star field**: three layers of small tinted stars, in board coordinates (so both
views show the same ones), drift to the left, swaying slowly up and down, behind the creatures, the near layers brighter
and faster. It is drawn procedurally in `draw.glsl`, so it costs no texture and no CPU, and turns
off with the music.

## Stats mode

Page modes, cycled by the **Mode** button or key `A` (`Shift`+`A` backwards): Action, Stats,
Action+Stats, Tuning; the page always starts in Action. Tuning and stats show their controls
between the screen and the help text. **Stats** shows a phase-space chart of the focused creature (`statschart.js`)
in place of the ship view (the right third of the screen; its controls go below the screen). In
**action+stats** mode the chart sits left of the ship view instead, over the right half of the world
view, so the ship view stays in play. The ◀ button in its top-right corner moves it to the left
edge (the left half of the world view); ▶ moves it back.
Every sim step (the app analyses after each step, not once per frame) records: t, total mass,
mass of each channel (R, G, B, summed over the creature's 4×4-cell blocks by a second target of
`blocks.glsl`), speed and turn rate. Speed and turn rate are the yellow arrow's. Units do not
depend on R, and are a spaceship's: 1 m = R cells, 1 s = one unit time (T steps), 1 g = the mass of
1 m² at cell value 1 (the sum of cell values / R²). So mass in g, speed in m/s, turn rate in
rad/s (its absolute value: either direction counts the same).

- **2D / 3D** button and a stat per axis: time, total mass, channel masses (3 trajectories, R G B),
  total + channel masses (4), hue (the HSV hue of the R G B masses, 0–360°, wrapping at red;
  derived from the channel masses when drawn, not recorded; no point where R = G = B), speed, turn
  rate. Default: 3D, total mass (x), hue (y), speed (z).
- The axis pickers (x, y, z, c) are on the stats panel's second row, with **Lock axes** (on by default): it keeps
  every axis's min and max where they are (a stat changed meanwhile gets a new range, then keeps
  it); in 2D, points outside are clipped. Loading another species while locked sets new ranges:
  locked again once its saved history has loaded, or (no file) once the short history has filled.
- **Colour axis (c):** the same stats, default turn rate (plus none: each trajectory its own colour).
  It sets each sample's hue: for hue, the creature's own hue; for any other stat, its range (outliers
  left out, eased like the others) scaled from blue (240°, min) to red (0°, max). Brightness still
  shows age: the short history fades from bright to dark, the long history is dark. A legend
  ("c: …" with a colour bar) sits top-left (3D: under the axis rows); click it to step the stat.
  Clicking an axis label in the chart (or an axis row of the 3D legend) steps it to the next stat;
  so do keys `D` `F` `G` `H` for x, y, z and colour (`Shift`: the previous stat). Axes
  with several series pair them by colour. Drag the 3D view to rotate it (z is up).
- **Short history:** the last 50 steps, always kept (steps with no focus as empty points, so the
  line runs out after the focus dies), drawn fading from the brightest colour (newest,
  with a dot) to a darker one: white to grey, or bright to darker R/G/B for the channels.
- **Long history:** unlimited, filled by **Record** (button, key `S`): each press records the next
  300 steps of the focus: its short history gathered so far (up to 50 steps, back to its last
  gap) and then new steps, 300 in all; pressing again stops early. A new focus (a click, Add,
  Reset) or the focus dying also stops it early, keeping what it has; **Delete** during a
  recording deletes the piece being recorded; changing species discards it (the species's saved
  file is loaded afresh). The
  readout shows `rec ●` and the progress. Drawn in a dark colour; each sample also stores the morph selected
  in the picker at the time, and each morph is labelled once with its name ("~" before a fuzzy one), at a point on
  its line where the label overlaps no other (from the middle of its pieces outwards; if no
  place is free, at the middle anyway). **Delete** deletes its last piece.
- **Saved long histories:** **Download** saves the long history as `<species>.json` (e.g.
  `5N7KKM.json`: R, dt, and each piece's morph and samples, 4 significant digits; one sample a
  line, each segment's columns aligned on the decimal point; about 23 KB per 300-step piece). Put it in `assets/stats/` and commit it: whenever that species is loaded, the
  page loads the file into the long history, so new recordings add pieces to it and the next
  Download saves them all. A file recorded at another R or dt is skipped (a note in the console).
  In the file, `segs` is a list of segments, each one morph's run of samples: `gap: true` means
  it starts a new piece (a new Record or a new focus), so its line is not joined to the one
  before and its time restarts at 0, and morph detection skips its first 10 samples; `gap: false`
  continues the piece before it under another morph name (the picker changed mid-recording).
  `fuzzy: true` marks a fuzzy morph's segment (see *Morph detection*).
  A species with no file just requests it and gets a 404. Changing species clears both
  histories.
- **Every object** button (on by default): every other detected object keeps a short history too,
  but is drawn as a short tail (its last 5 samples; its whole short history when time is on an
  axis), under the focus's fading line. Clicking one of them
  in the world view focuses it and shows its line at once (its history comes with it; the old
  focus becomes a dot, keeping its own). Hovering the mouse over an object in the world view shows
  its line as well. In this mode (and in tuning mode) the world view draws a circle around every
  object: the focus in yellow, the hovered one brighter. Each object is followed step to step by nearest predicted
  position (within R), with its own copy of the arrow's speed and turn-rate filters; a lost one
  disappears (dropped after 10 steps unmatched). Record still records the focus only.
- A reset (Z) clears all short histories (the long one is kept). A new focus (spawn, reset, click, add) starts a new piece, with a gap in the line; so does each
  Record. So one recording can cover several basins of attraction of one species. t
  is relative: each long-history piece starts at t = 0 (also when loaded from a file), and the short
  history runs from 0 at its oldest sample, so with time on an axis the pieces start together.
- The morph label is the picker's morph, not a detected one: after clicking another creature
  (e.g. a different morph added earlier), the label is still the picker's.
- **Outliers ignored:** a sample far from its time neighbours (5 on each side), more than 8× their
  own spread (robust: median and MAD), is left out of the lines (they join its neighbours) and of
  the axis ranges: a one- or two-step spike, e.g. a detection glitch, is dropped; a lasting jump is
  kept. Axis ranges grow at once but shrink gradually (8% of the way per frame), so the chart
  doesn't jump when an extreme point ages out.
- None of the chart's settings are in the URL.
- From the console: `cognoshipStats.long.get('speed', j)`, j = 0 the oldest, up to `cognoshipStats.long.len - 1`
  (`.short` likewise); `cognoshipStats.long.segs` lists the pieces and morphs.

### Morph detection

`detect.js`, used by the Detect button. Each morph's long-history samples (the saved file plus
anything recorded this session, grouped by the morph recorded with them) are a cloud in a feature
space: log(R, G, B mass + 0.01) (total mass left out, it is their sum; the 0.01 g gives an
empty channel a value, e.g. XEH4YR *cor*, which has no B: before it, cor was never detected), log(speed + 0.002) and
log(turn rate + 0.01), each divided by a fixed kernel width (masses 0.03, i.e. 3%; speed and turn
rate 0.3). 10 samples after each gap are skipped while the speed and turn filters settle, and each
cloud is thinned to 250 points. A creature's newest samples (up to 50 since its last gap, after
the same 10) are scored against each cloud with a Gaussian kernel density; the mean score per
morph goes through a softmax (gain 2) to give a probability. The order of the samples in time is
not used, so where in its cycle a creature is doesn't matter.

- Why logs and fixed widths: the still and the slowly drifting morphs (*lithos*, *tardus*) have
  nearly the same masses and differ in speed (~0.0001 vs ~0.02 m/s), and a copy added at another
  angle has masses 1-2% off the recorded one. The first version scaled each stat by its spread
  over the whole species: speed's spread was set by the fast morphs, so the still / drifting
  difference counted for little, and in a species with a heavy morph (Z18A9R *paragyrans*, twice the
  mass) the mass widths were huge too, so lithos and tardus were too close to name either (???).
  A still morph's cloud is nearly a point however many copies are recorded, so the kernel width,
  not the data, decides how far a copy may be off: a fixed relative width covers that.
- A morph is named at a probability of 0.9 or more, if the samples also fit its cloud (mean
  score ≥ −6, about 3.5 widths off in all). Otherwise it is unknown (**???**), as it is before 25
  usable samples.
- Needs a long history of two morphs or more (a species with no saved file shows only **???**).
- Offline check on the saved files (each morph's first half of samples as the reference, windows
  of its second half as the test; `node` with `statschart.js` and `detect.js`): 91% of windows
  named, 99.3% of those right; the misses were HGMLH8 *metamorpha* (it passes through *tachus*'
  and *tardus*' shapes; then named *para* and *parapedes*). With one morph left out of the reference, its windows were wrongly
  named another morph 10% of the time. Lithos and tardus of 5N7KKM, Z18A9R and G6G6CR stay
  correctly named with every mass 2% off either way, and still at 5% off; two live 5N7KKM lithos
  (red mass 2% under the recorded one) that the first version called ??? and tardus are lithos.
  Both halves come from the same run, so a fresh creature may do worse.
- **Fuzzy morphs.** Some forms vary too much between creatures for the fixed widths, e.g.
  VT049W **cytos**: the nearly still, slowly shape-shifting lithos-like offspring of the
  replicating gyrans (mass 2.5-2.8 against lithos's 3.3, drifting ~0.01 m/s, its turn rate mostly
  noise). A piece marked `"fuzzy": true` in the file is matched with widths from the morph's own
  spread (0.75 x its standard deviation per feature, at least the fixed width), and the softmax
  compares normalised densities, so a broad cloud pays for its breadth: a precise morph that fits
  wins, and the fuzzy one is named only where none is near. Shown as **~cytos** (with a "~"). To
  add one: **Shift+S** (or Shift+click Record) records a fuzzy piece and asks for its name (the
  last fuzzy name is offered; plain Record / `S` also asks, offering the picker's morph; it need not be in the morph picker); Download to keep it. A fuzzy
  piece named like a precise morph is kept apart from it. Check
  (VT049W, each of the six cytos pieces left out of the reference in turn): every window named
  ~cytos; rectus, lithos and tachus still right, and each of them left out of the reference stays
  ??? rather than ~cytos. The other species have no fuzzy pieces and are unchanged.
  HGMLH8 **metamorpha** (the shapeshifter) is fuzzy too. Its two pieces differ a lot (one slow
  and close to *tachus* / *tardus*, one faster and turning), so with one as the reference the
  other is mostly ??? (6 of 25 windows named). With four pieces, each left out in turn: 11, 25,
  20 and 18 of 25 windows named ~metamorpha (most of the rest ???, a few *tachus* / *tardus*). While it passes through
  *tachus*' or *tardus*' shape it is named as those (a precise morph that fits wins).
  Cilium (G6G6CR) *munitus*, *lithos*, *tardus* and *rectus* are fuzzy too (they keep changing
  shape); *limus* stays precise. Check (each morph's first half as the reference, windows of its
  second half as recorded and with every mass ±4%, 150 windows): named right 107 with none fuzzy,
  124 with all five, 128 with all but limus, never a wrong name; munitus alone went from 1 to
  10 of 10 as recorded. With limus fuzzy too, an unrecorded rectus was named ~limus (11 of 65
  windows, against 4 of 65 as limus otherwise).
  VT049W *gyrans* (the replicating one) is fuzzy too. Its pieces are short (26-52 steps: it
  divides), and as a precise morph its 52-step piece, left out of the reference, was unnamed in its
  first 6 of 18 windows; fuzzy, all 18 named, and rectus, lithos, tachus and cytos unchanged.
  Pieces under 35 steps are too short to be named at all (25 samples after the 10 skipped).
  How much to record (cytos references from 1-5 of the six pieces, the rest as tests): 1 piece
  of 300 steps names 83% of windows, 2 pieces 93%, 3 pieces 97%, 5 pieces 100%; pieces of
  150 steps: 53%, 84%, 89%, 95% for 1, 2, 3, 5. So 3+ pieces of 300 steps, each a different
  creature, which matter more than length. In no case was an unrecorded lithos called ~cytos
  with full pieces (up to 7% with 75-step pieces).
- The focus is classified every frame; in Every mode the other creatures' tracks (the
  every-object tracks, kept while Detect is Every even with *Every object* off) 4 per frame in turn.

## Tuning mode

Tuning mode (key `A`) shows an extra row of detection sliders and circles every detected object in the world
view (the focus in yellow): **Merge saddle** (basins merge when the saddle is at least this
fraction of the lower peak; higher = fewer merges), **Threshold** (potential below it is
background; higher cuts thin bridges between creatures) and **Min area** (in units of a disc of
radius R/4). Defaults, chosen by hand to stop the ship view jumping between merging and splitting
objects in a crowded board: threshold 0.1, min area 1 (were 0.05, 1; threshold 0.15 for a while, but it left
1–3% of a creature's mass, at its edges, outside the object, against under 1% at 0.1), and a merge saddle per species, set in
`tools/gen_aquarium.py` (`DETECT`) and written into `aquarium233.js` as `rules[code].detect`: 0.95 for
Aquarium (5N7KKM) and Aquarium replicating (VT049W), 0.9 for Cilium (G6G6CR) and Aquarium spinning
(Z18A9R), 1.0 for Aquarium dancing (KH97WU) and drifting (GDNQYX); it was 0.8 for all. Above 1.0 nothing changes (a saddle never exceeds the lower peak, so no distinct peaks
merge from 1.0 up). The readout adds **in objects**: the share of the board's mass inside detected objects (the
rest is in blocks under the threshold or in dropped small basins). Also
settable as `?saddle=`, `?thr=`, `?area=`.
