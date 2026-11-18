# Cognoship: the full reference

Every detail of the game and the lab: controls, layouts, modes and the inventory, URL parameters,
the script syntax, how it works, the fidelity checks, stats mode, morph detection and tuning. For
the essentials, see [README.md](README.md) (the game, for players) and
[README-LAB.md](README-LAB.md) (the lab).


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

In the game, clicking the selected **R G B** patch again (or pressing its key again) switches the tool from masking to drawing, where `[draw]` allows it (cross-hatched, +R / +G / +B), then to erasing ("recycle": the channel's mass under the circle set to 0, as the lab's Erase; striped, −R / −G / −B), then back to masking: mask → draw → erase → mask. Clicking another channel goes back to masking. While drawing on the world (`[world_draw]`, Nebula or Antimatter in use) the patches instead pick the field's colour, and the selected one again its tool: +R paintbrush, ×R soft brush, −R eraser. The cursor is a hand over the patches. Holding **Shift** or **Alt** shows each control's key in a small badge at its corner (sliders: one at each end; the R G B patches: 1–3).

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
- **Spray** (`X`): add the pattern, rotated by a random angle, at a random spot, the emptiest of 24 candidates, without clearing
  the board. If nothing is focused, it focuses the new creature; after *Dissolved*, it restarts
  the run.
- **Add** (`B`): hold the morph to paste it, as the game's clone mode (`clone.js`), but in both
  views: its ghost follows the mouse over the world view or the ship view (shown in both, at one
  place on the board), and a click in either pastes it there; the simulation pauses while it is
  held, then plays on as before, following what was pasted. Pink round buttons at the top left of
  the world view: ↻ turns it 45° clockwise (`N`), ↔ flips it left to right as seen in the view the
  pointer was last in (`M`), ✗ cancels (`Esc`). A touch places it instead (a drag moves it), and
  ✓ appears to paste it; `Enter` pastes it where it is placed, or under the mouse. Picking another
  morph while holding one holds that one instead (no spray); another species, or Play, drops it.
- **Speed** (`Q` faster, `W` slower): 7 slots, 1/8, 1/4, 1/2, 1 (the default), 2, 4, 8, at 30 steps per real second each,
  by the clock, not per frame, so the pace is the same on any screen (60 or 120 Hz) and at any frame rate: speed 1 runs
  30 steps a second (one every 2 frames at 60 Hz, every 4 at 120 Hz), speed 8 240. A frame runs at most 4 steps, and no more than fit in 28 ms by the measured time of a step (at
  least one), and counts at most 0.1 s, so a device that cannot keep up plays slower, smoothly, one step a frame, instead of
  piling up steps in long frames. Speed doesn't change
  the physics, only how many steps you see per second. The time step dt = 1/T is a separate, fixed
  property of the rule (T = 2, see below): it sets how far the simulation advances per step, and
  changing it changes the dynamics.
- **Zoom**: screen px per cell in the ship view (2× to 16×, the game 3× to 12×; the default is the slot nearest
  4 R across, 6× = 64 cells at R = 18).
- **Species** and **Morph** pickers. A species is a rule code, named in the picker as Aquarium
  originale (5N7KKM), A. novus (GDNQYX), A. vorticosus (Z18A9R), A. multiplicans (PJ1EF7),
  A. explodens (XEH4YR), A. virgatus (KZ2BLD), A. saltans (KH97WU), A. sinuosus (HAESRE),
  A. metamorpha (HGMLH8), A. replicans (VT049W) and Cilium immunitatis (G6G6CR), in `sandbox.txt`'s order. Choosing one resets the
  world with its first morph. A morph (limus, …) only changes what Reset, Spray (`X`) and Add (`B`) put in, and
  keeps the world. Picking a morph in the list (even the current one) also adds it, like Spray;
  ▲ ▼ or `C` `V` only step through the morphs; `Shift`+`C` / `Shift`+`V` step to the
  previous / next species (and reset, as picking it). 63 entries of the aquarium collection from 5N7KKM to
  HAESRE, in file order (the last five codes added 2026-11-05, not yet checked in the browser or
  numpy; KZ2BLD's unnamed second morph is named *virgatus* here, its *gyrans tachus* is *limus*,
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
  - **Z18A9R trias, lithos and tardus** were found in Cognoship itself (at R = 18). They come after
    the Z18A9R forms and were rebuilt with `tools/from_screenshot.py` from the screenshots in
    `assets/screenshots/`. All three use a brightness correction of 1.8; the size is 1.0 for trias
    and lithos and 1.05 for tardus, with shape similarity 0.96, 0.93 and 0.92.
    - **trias** (was *trio*): three lobes, drifting slowly.
    - **lithos:** a round, still disc. Without the correction it died or became a spinning
      crescent.
    - **tardus:** a round form with a green and red core.
  - **PJ1EF7 lithos** was found in Cognoship too (at R = 18) and rebuilt the same way: correction
    1.8, size 0.95, shape similarity 0.95. A round disc with a yellow and green core. In the numpy
    check it keeps its mass for 500 steps, then at about 600-700 steps turns into the crescent
    the other PJ1EF7 morphs have. Every other correction and size became the crescent outright.
    Extra patterns live in `tools/extra_patterns.json`.
  - **Z18A9R trigyrans, currus and flos** (2026-11-18, in the lab only): trigyrans, three gyrans
    rotating round a common centre, made stable by a rule tweak (kernel 4 `m` 0.398 → 0.41,
    `tools/rule_tweaks.json`); currus (a slow oval) and flos (a still flower) use only the red and
    green channels. See MORPH.md.
  - **To add a pattern you see in the lab or the game, capture it** (better than a screenshot):
    press `,` to download the ship (the object followed, or its whole once named: a twin, a trio)
    as it is, or `.` for the whole board. You get a pattern JSON with the exact cell values (only
    the ship's own cells, its object's blocks grown by 2 blocks, every other object's left out;
    cropped) and a PNG of its RGB channels (no stars, circles or arrows). Pausing first makes the
    moment easy to pick. Then:
    ```
    python3 tools/add_pattern.py "~/Downloads/Z18A9R unnamed 2026-11-18-10-20-30.json" --name trigyrans [--check]
    ```
    It adds it to `tools/extra_patterns.json` and regenerates `aquarium233.js`; `--check` (numpy,
    scipy: e.g. the complexity venv's python) first runs it 500 steps and reports the number of
    creatures and the mass drift. The file's name is the detected morph, or "unnamed".
  - **To add another pattern from a screenshot** of the world view (default zoom is fine), run
    `tools/from_screenshot.py`. It needs numpy, scipy and Pillow, e.g. the complexity venv:
    ```
    ~/github/complexity/.venv/bin/python tools/from_screenshot.py shot.png --code Z18A9R --name trias --add
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
  - Left out: *G6G6CR metamorpha beaked* and *VT049W gyrotachus*, which is a
    different species (its kernel 0 has μ = 0.302, not VT049W's 0.272). The three unnamed Z18A9R entries are
    *gyrans*, *paragyrans* and *antigyrans*; VT049W's two *para* entries are shown as *rectogyrans* (until
    2026-11-18 *paralimus*) and *paragyrans*, and its *rectus+gyrans* as *sequerectus* (until then
    *rectogyrans*); `tools/gen_aquarium.py` also lists VT049W's morphs in the lab's picker in the
    morph book's order (`MORPH_ORDER`). Also renamed: G6G6CR *fortified* → *munitus*, KH97WU *antipara* → *antigyrans* and
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

## Layouts (the game)

Side by side (a landscape window, and fullscreen): the world view on the left 2/3 of the canvas, the
ship view on the right 1/3, the message under the world view and Retry / Pause / Speed (one row)
under the ship view. Stacked (a portrait window, not fullscreen, e.g. a phone held upright): the
world view across the full width at the top, the ship view at 1.5x below it, at the right edge, with Fire
and ◀ ▶ left of it at the screen's left edge (level with the ship view's bottom, above Upward / Zoom), then Retry / Pause / Speed, then the message, two lines
high, in a larger font (sized so the script's widest message fits two lines). The chapter card
shrinks to the window, and the books take nearly its full width, their list scrolling sideways when
the slots don't fit (the bullet and species books' slots wide enough for a name on one line), with the
species link, Hold to clear data and Close each on its own row under it. The layout follows the window as it turns or resizes. Underneath, one canvas
either way (`core.js` `setLayout`): stacked, it is 768 × 1296 px, the world view (768 × 512), a
16 px gap, and the ship panel drawn at 1.5× (576 × 768), at the right edge; the overlay is drawn in the same
logical coordinates as side by side, mapped to each part, and the page's buttons sit in a layer
over each view (`.layer-world`, `.layer-ship`), placed relative to it.

## Story mode, sandbox mode, challenge mode and the inventory

The game's title row has three links, centred: **Story** (`index.html`), **Challenge**
(`index.html?mode=challenge`) and **Sandbox** (`index.html?mode=sandbox`); the current one
is underlined, and other URL parameters (e.g. `?dev=1`) are kept when switching. All are the same
game, with the same UI and gameplay; story mode runs `script.txt`, sandbox mode `sandbox.txt`,
challenge mode `challenge.txt`.

The **inventory** is the discovered morphs of every species (each with its colour and picture) and
the bullet patterns the gun has. The **central inventory** is kept in the browser
(`localStorage`, key `cognoship-inventory`) and only grows: whatever a run gains (a discovery; a morph's new picture or colour; a bullet pattern added by `[add_bullet]` or
created) is merged into it at once, during the run, not at its end. Except what is given for testing:
morphs added by `[discover_morph]` (until the ship holds one long enough for its picture, which earns
it) and patterns added by `[add_bullet:all]` stay out of it. A story run always starts from an
empty inventory. A sandbox run starts from the central inventory (loaded before `sandbox.txt` runs:
every species' discoveries, and the bullet list), so it has everything found in every story run and
every sandbox run, and saves into it the same way.

The **energy** (the dots score) is in the inventory too, for sandbox mode only, and it can go
down: creating a bullet pattern spends it (`[new_bullet]`). A sandbox run starts from the central
energy and saves into it what it gains and what it spends (as a change, so two tabs add up). A
story run's energy is its own: it starts from 0 and is never saved.

In sandbox mode only, a **Species** button (`C`) sits at the top of the ship view, centred between
Bullets and Morphs, and the species' name at the morph book's lower right is a link too. Both open
the **species book**, "Choose a species" (as wide as the bullet book): a slot per species added by
`[add_allow_species]`, in order, with a small thumbnail of each of its visible morphs in one row
(up to 8: the centre of its picture if discovered, the whole picture for a twin or trio, i.e. a morph recorded as a whole, and for every Cilium morph (large creatures), listed in `aquarium233.js` `rules[code].wholes` by `tools/gen_aquarium.py` from the stats files; black if discovered without one; the book's own colour if not), above its name; the
current species has a lighter background (so does the player's morph, the one Respawn puts back, in the morph book while choosing a morph to respawn). Every book's title shows its verb in the colour of the slots that do it (the morph book's "respawn" in cyan and "clone" in pink; "Choose" (species) or "use" (bullets) in cyan and "Create" in green in the others). Pickable slots in every book have a darker cyan border than the Respawn and Choose buttons. Keys pick slots, the book's rows on the keyboard's: `Q`–`R`,
`A`–`F`, `Z` `X`; only Escape closes it (as every book; its Close button shows *Esc*). A
pick fades the views out as `[new_scene]` does, does `[choose_species]` with it unseen (its `[on_start]` block sets the world up), fades them back in and plays
on. Retry then redoes the species last chosen, not the script's first `[choose_species]`.

In sandbox and challenge mode, every book (morph, bullet, species) has a **Hold to clear data** link
at its lower left (styled as the species link). Held for 3 seconds (a red bar fills under it; let go
sooner and nothing happens), it empties the central inventory and reloads the page: a fresh run with
nothing discovered. The challenge book has **Hold to clear completion & best time** there instead:
it takes the ticks and the best times off every challenge, without a reload.

**Challenge mode** (PLAN.md §11) reads and writes the central inventory as sandbox mode does (the
morphs and bullet patterns), but its energy is its own: each attempt starts from 0. `challenge.txt`'s
lines before its first `[challenge]` set the game up (their actions run once, at the start); then
each `[challenge:<name>:<n>[:<subtitle>]]` line starts challenge `<name>`, stage `<n>` (a challenge may
have several stages, each its own line), which runs to the next `[challenge]` line (or the file's end).

- **The challenge book** opens at the start, and from the **Challenges** button (`C`, where sandbox
  mode has Species): a slot per challenge, in order, titled `<name>: <subtitle>` (no subtitle: `<name>: <n>`, e.g. *Courier: 1*), after a green ✓
  once it has been completed with every species it allows. Picking one
  opens the species book on the species it allows (even one), each with this challenge's record
  with it, in a row of its own below the species' name: a green ✓ and its best time (*12.3s*,
  *1m 12s*) once completed, a blank row until then. Completion and best times are kept per challenge and species; this species book's
  **Hold to clear completion & best time** clears that challenge's, with every species.
- **Starting:** a pick asks for a species in the species book (only those its `[allow_species]`
  allows, even if only one), then the views fade out, the species is set up (as a
  sandbox pick: `[choose_species]`), and the challenge runs from its `[challenge]` line: a card
  "CHALLENGE", `<name>`, "STAGE `<n>`".
- **Completing:** reaching the challenge's end (its last line's conditions met, or moved on from)
  completes it: SUCCESS (in challenge mode only here, not for each line on the way), a message
  with its time (game seconds from the start of the attempt; paused time doesn't count), the
  record saved (`localStorage`, key `cognoship-challenges`: `{ "<name>:<n>@<species id>": { done, best } }`, so a subtitle may change,
  the best time kept), and 3 s later the species book opens on that challenge (its records with each species, the one just made among them; Back goes to the challenge book).
- **Failing:** the ship dissolving (with `[dissolve:retry]`; by default it is "Dissolved" with
  Respawn as in the other modes), `[time_limit]` running out, or the energy reaching 0 under
  `[mask_cost]` fails the attempt: "Failed" in the ship view (in place of "Dissolved"), a message
  saying why, Respawn and Choose; the script stops waiting.
- **Retrying:** Morphs (`X`), Retry (`Backspace`) and Choose (after a failure) open the morph
  book as "Choose a morph to retry" (but while cloning is allowed, `[clone:on]`, e.g. in the
  header, kept across challenges until `[clone:off]`, and the ship is alive, Clone (`X`) opens the
  book to clone from, as in the other modes; Retry still opens the retry book); a
  pick, or Respawn (`Z`) after a failure, starts the attempt again behind a fade: from the line
  after `[challenge]` (no card), the energy from 0, the time from now, and its first
  `[put_player]` putting that morph (Respawn: the player's morph, the one with the lighter
  background) in place of its own. After a challenge is completed, the same plays it again from
  its card.

Not built yet: locked challenges and the scrolling world.

## URL parameters

- `?form=5N7KKM-limus`: the pattern to start with (the default), by id: code and form joined by `-`, e.g.
  `5N7KKM-limus`, `VT049W-gyrans-replicating`, `5N7KKM-lithos+gyrans`. `aquarium` and
  `replicating` are short for the two VT049W ones. Hidden: `fission`, `fission02`.
- `?R=18`: kernel radius in cells (6–30). The board stays 384×256, so a larger R makes the
  creature larger and the sim slower. Below 16, *gyrans R* stops replicating and *lithos* moves.
- `?dt=0.5`: time step, overriding the rule's own (1/T = 0.5).
- `?speed=1`: the Speed slot at start, snapped to the nearest one (1/8–8; 30 steps per second per unit).
- `?zoom=8`: ship-view px per cell at start, snapped to the nearest slot (2–16).
- `?radius=4`: mask radius at start, in cells (default 2.5 · R / 12, rounded to 0.5).
- `?headup=vel|axis`: start in heading-up or body-axis-up view (`1` = `vel`).
- `?debug=1`: show GPU ms per step and the tap count in the readout.
- `?fft=0` (game and lab): FFT mode off: the creatures' convolution by the direct sum (`sim.glsl`) instead of by FFT on the GPU (`fftsim.js`, on by default, see *FFT mode* below). With `?dev=1` the game's readout ends in "· FFT" while FFT mode is on.
- `?all=1` (the game, either mode): all discovered, for testing. Every visible morph of every species counts as discovered (no picture until the ship holds it), and the gun has every bullet pattern from the start. The central inventory is neither read nor written.
- `?phone=<device>` (the game): the phone test (`1`: the default device). Opened on its own, the page goes to `phone.html`, which shows the game in a phone-sized frame (a phone's CSS viewport, upright, scaled down to fit; picked from a list or by the parameter: `pixel9` (Google Pixel 9, 412 × 923, the default), `iphone14` (390 × 844), `iphone15max` (430 × 932), `iphonese` (375 × 667), `galaxys24` (360 × 780), or a size of one's own in the boxes beside the list (or `&w=` `&h=`); the whole screen, so a little taller than a browser leaves; no scrollbar takes width in the frame, as on a phone. Text may still wrap a little differently: a phone's `system-ui` font (Roboto on Android, SF on iPhone) differs from the desktop's. A phone's real viewport also depends on its display-size setting: `?dev=1` shows it, as width×height @pixel ratio, at the end of the readout), so its phone sizes and the stacked portrait layout apply; the other parameters are kept (e.g. `?phone=iphone14&dev=1&fft=0`). In the frame the game acts as a phone (the small music files) and the music card asks every time, whatever was chosen before.
- `?dev=1` (the game): dev mode. Test aids, as keys and as buttons at the top left of the world view (**Ch**,
  **Ok**, **Do**): `I` goes to the next `[chapter]`, running every action on the lines in between
  (as `[jump_here]` does); `O` completes the waiting line's `[need]` conditions at once; `P` scores
  as if the ship ate a dot (+10 energy). The frame rate, simulation steps per second (over the last
  half second) and the measured time of a step at the top right of the world view. A grey circle around every detected object in the
  world view, besides the ship's yellow one (the one under the pointer brighter).

## The game's script (`script.txt`)

The scripts that drive the game (`script.txt`, `sandbox.txt`, `challenge.txt`): their syntax, every
tag and every condition are in [SCRIPT.md](SCRIPT.md).

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

## FFT mode (on by default; `?fft=0` turns it off)

A plug-in (`fftsim.js`, `shaders/fft*.glsl`) that replaces the creatures' sim pass: `sim.glsl` sums
~1,000 kernel taps per cell, a cost that grows with R²; the FFT's does not. With `?fft=0` the
file does nothing; to remove it, delete it, its `<script>` tags and `shaders/fft*.glsl` (`glsim.js`
keeps two small hooks, `internals.ext.step` / `.rule`, inert then).

- **Per step:** `fftpack` writes the fields as complex pairs, (R, G) and (B, bullets); two forward
  FFTs; `fftmul`, per pair of kernels, multiplies the spectrum of the channel each reads by the
  kernel's spectrum (made once per rule, from the kernel as an image of the board, flipped: the
  direct sum is a correlation); eight inverse FFTs give U for all 16 kernel slots; `fftgrow` does
  the update exactly as `sim.glsl`.
- **The cognitive mask:** its circle is small, so `fftgrow` takes its cells out of U and of the
  normaliser directly, for the pixels within reach of it (the kernels as a small image, `uKimg`):
  the same U = K∗(A·(1−B)) / K∗(1−B).
- **The bullets' mask:** it covers the board, so it is done in the FFT: two more forward FFTs ((R,
  G)·b and B·b, b the bullets field on the channels they act on) taken out of the numerators, and
  eight more inverse ones for each kernel's hidden share (K∗b), so while bullets fly a step costs
  about twice as much. Erase, draw and the bullets' erase act after the update, as in `sim.glsl`.
- **The FFT:** mixed-radix (4, 3, 2) Stockham passes in gather form (`fft.glsl`), along x (384 =
  4·4·4·2·3) then y (256 = 2·2·4·4·4, a 4 split so the count is even and each transform ends in its
  own texture); twiddles from a table computed on the CPU in double precision (`sin` / `cos` are
  coarse on some mobile GPUs).
- **Checked on the CPU** (node, before the GPU code): the passes against a direct 2D DFT (relative
  error 3e-13), the convolution against the taps (6e-17), and the mask and bullets corrections
  against `sim.glsl`'s masked sum (4e-9, the kernels' float32 normalisation); with every operation
  rounded to float32 as on the GPU, U is within 1.4e-7 of the direct sum.
- **Memory:** 29 extra RGBA32F textures of the board (about 46 MB).

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

Background music, royalty free: a playlist played in turn and looped, with 10 s of silence between plays
(the first, until the game's `[music:<list>:<play>]` sets a playlist), at half volume, from `assets/music/`: *Long Way* by Rexlambo (`rexlambo-long-way.mp3`), *A Few Jumps Away*
by Arthur Vyncke (`arthur-vyncke-a-few-jumps-away.mp3`), *Out of the Skies, Under the Earth* by
Chris Zabriskie (`chris-zabriskie-out-of-the-skies-under-the-earth.mp3`) and *Billions and
Billions* by Stellardrone (`stellardrone-billions-and-billions.mp3`). In the lab, the **Stars**
button at the end of the second row turns the star field and the music on or off together; in
the game the star field is always on and the **Music** button (a music note at the top right, left of the light / dark button, highlighted while on; `M`) only plays or pauses the music
(on by default, then as last chosen in this browser; each track is downloaded only when it starts to play, and streamed). On a phone (a touch screen), the tracks come re-encoded small, `assets/music/small/*.m4a` (AAC 64 kbps: 8.2 MB for the four, against 19.5 MB; made with `ffmpeg -i <track>.mp3 -c:a aac -b:a 64k -ar 44100 -movflags +faststart small/<track>.m4a`); and unless the browser reports Wi-Fi (Chrome on Android, without Data Saver; an iPhone reports nothing), a card at the bottom asks once, "Play music? About 1–3 MB per track of mobile data", the answer remembered as the Music button's. Browsers only start
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
  so do keys `G` `H` `J` `K` for x, y, z and colour (`Shift`: the previous stat). Axes
  with several series pair them by colour. Drag the 3D view to rotate it (z is up).
- **Short history:** the last 50 steps, always kept (steps with no focus as empty points, so the
  line runs out after the focus dies), drawn fading from the brightest colour (newest,
  with a dot) to a darker one: white to grey, or bright to darker R/G/B for the channels.
- **Long history:** unlimited, filled by **Record** (button, key `S`): each press records the next
  300 steps of the focus: its short history gathered so far (up to 50 steps, back to its last
  gap) and then new steps, 300 in all; pressing again stops early. A new focus (a click, Spray, Add,
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
  axis), under the focus's fading line. So is every recognized twin or trio (a whole named, as its
  circle in the world view shows it), as one object, besides its parts. Clicking one of them
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
  named, 99.3% of those right; the misses were HGMLH8 *metamorpha* (it passes through *paratachus*'
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
  add one: **Record Fuzzy** (`F`) records a fuzzy piece and asks for its name (the
  last fuzzy name is offered; Record / `S` also asks, offering the picker's morph; it need not be in the morph picker); Download to keep it. A fuzzy
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

**Wholes (the game and the lab).** A twin or trio (paragyrans, rectogyrans, trias…) is one
creature of several parts, which the merge saddle sees as separate objects. A second merge at the
**twin saddle** (`rules[code].detect.wholeSaddle`; the Twin saddle slider; `?wsaddle=`) joins parts
into a *whole* (creatures.js `wholes`). Wholes are tracked and named on their own: the long
history stores each piece's level (Stats mode's **Record Whole**, `D`; `"level": "whole"` in the
file), parts are compared only with morphs recorded as parts, wholes only with those recorded as
wholes. A whole shows, as one circle in its morph's hue with its name (its parts' names hidden),
once named 3 classifications in a row, and goes once not named 3 in a row; a whole of creatures
merely touching is never named, so never shown. Tuning mode draws every whole dashed. While the
followed creature is part of a shown whole, the **ship** is the whole (core.js `ship()`): the ship
view centres on it and turns with its heading, the arrow and the world view's frame and circle are
its, the tool and the page's presses (the lab's Add, the game's clones) are placed relative to it,
and its name is the name box's, the game's discoveries' and the morph book picture's (centred on
it). The part is still what is tracked, followed and recorded (Record Whole records its whole).
A new twin is followed by its part at first, until its whole is named: then the view moves to the
whole's centre. The twin
saddles per species: 0.6, except 0.4 for A. saltans (KH97WU) and 0.2 for A. vorticosus (Z18A9R:
transgyrans's gyrans part between crossings, and at 0.6 it is whole only half the time, in
10-step stretches, too broken to be named; at 0.3, after the trigyrans rule tweak, it still
dropped out for about a third of its cycle; touching Z18A9R creatures join at 0.3 in 90 % of
contacts anyway, so 0.2 costs little), chosen by eye in the lab, starting
from a measurement (`output/whole_saddle/measure.py` and `tradeoff.py`: each twin run long, the
highest twin saddle at which its parts stay one whole all, or nearly all, of the time). A loose twin (gyrans orbiting a lithos) is whole only while
its parts touch, at any value.

Tuning mode (key `A`) shows an extra row of detection sliders and circles every detected object in the world
view (the focus in yellow): **Merge saddle** (basins merge when the saddle is at least this
fraction of the lower peak; higher = fewer merges), **Threshold** (potential below it is
background; higher cuts thin bridges between creatures) and **Min area** (in units of a disc of
radius R/4). Defaults, chosen by hand to stop the ship view jumping between merging and splitting
objects in a crowded board: threshold 0.1, min area 1 (were 0.05, 1; threshold 0.15 for a while, but it left
1–3% of a creature's mass, at its edges, outside the object, against under 1% at 0.1), and a merge saddle per species, set in
`tools/gen_aquarium.py` (`DETECT`) and written into `aquarium233.js` as `rules[code].detect`: 1.0,
except 0.95 for A. virgatus (KZ2BLD) and 0.9 for Cilium (G6G6CR), chosen by eye in the lab
(A. replicans, VT049W, was 0.95 too, until 1.0 was needed for its pararectus); it was 0.8 for all. Above 1.0 nothing changes (a saddle never exceeds the lower peak, so no distinct peaks
merge from 1.0 up). The readout adds **in objects**: the share of the board's mass inside detected objects (the
rest is in blocks under the threshold or in dropped small basins). Also
settable as `?saddle=`, `?thr=`, `?area=`.
