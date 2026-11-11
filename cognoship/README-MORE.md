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

In the game, clicking the selected **R G B** patch again (or pressing its key again) switches the tool from masking to erasing ("recycle": the channel's mass under the circle set to 0, as the lab's Erase); the patch turns striped and reads −R / −G / −B. Clicking it again, or another channel, goes back to masking. The cursor is a hand over the patches. Holding **Shift** or **Alt** shows each control's key in a small badge at its corner (sliders: one at each end; the R G B patches: 1–3).

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
  original (5N7KKM), Aquarium replicating (VT049W), Aquarium spinning (Z18A9R), Aquarium dancing (KH97WU),
  Aquarium drifting (GDNQYX), Aquarium shapeshifting (HGMLH8), Aquarium multiplying (PJ1EF7),
  Aquarium striped (KZ2BLD), Aquarium exploding (XEH4YR), Aquarium zigzagging (HAESRE) and, last, Cilium (G6G6CR). Choosing one resets the
  world with its first morph. A morph (limus, …) only changes what Reset and Add (`X`) put in, and
  keeps the world. Picking a morph in the list (even the current one) also adds it, like Add;
  ▲ ▼ or `C` `V` only step through the morphs; `Shift`+`C` / `Shift`+`V` step to the
  previous / next species (and reset, as picking it). 63 entries of the aquarium collection from 5N7KKM to
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
    *gyrans*, *paragyrans* and *antigyrans*; VT049W's two *para* entries are shown as *paralimus* and
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

## Story mode, sandbox mode and the inventory

The game's title row has two links, centred: **Story Mode** (`index.html`) and **Sandbox Mode**
(`index.html?mode=sandbox`); the current one is underlined, and other URL parameters (e.g.
`?dev=1`) are kept when switching. Both are the same game, with the same UI and gameplay; story
mode runs `script.txt`, sandbox mode `sandbox.txt`.

The **inventory** is the discovered morphs of every species (each with its colour and picture) and
the bullet patterns the gun has. The **central inventory** is kept in the browser
(`localStorage`, key `cognoship-inventory`) and only grows: whatever a run gains (a discovery; a morph's new picture or colour; a bullet pattern added by `[bullet]` or
created) is merged into it at once, during the run, not at its end. Except what is given for testing:
morphs added by `[discover_morph]` (until the ship holds one long enough for its picture, which earns
it) and patterns added by `[bullet:all]` stay out of it. A story run always starts from an
empty inventory. A sandbox run starts from the central inventory (loaded before `sandbox.txt` runs:
every species' discoveries, and the bullet list), so it has everything found in every story run and
every sandbox run, and saves into it the same way.

In sandbox mode only, a **Species** button (`C`) sits at the top of the ship view, centred between
Bullets and Morphs, and the species' name at the morph book's lower right is a link too. Both open
the **species book**, "Choose a species" (as wide as the bullet book): a slot per species added by
`[add_allow_species]`, in order, with a small thumbnail of each of its visible morphs in one row
(up to 8: the centre of its picture if discovered; black if discovered without one; the book's own colour if not), above its name; the
current species has a lighter background (so does the player's morph, the one Respawn puts back, in the morph book while choosing a morph to respawn). Every book's title shows its verb in the colour of the slots that do it (the morph book's "respawn" in cyan and "clone" in pink; "Choose" (species) or "use" (bullets) in cyan and "Create" in green in the others). Pickable slots in every book have a darker cyan border than the Respawn and Choose buttons. Keys pick slots, the book's rows on the keyboard's: `Q`–`R`,
`A`–`F`, `Z` `X`; only Escape closes it (as every book; its Close button shows *Esc*). A
pick fades the views out as `[new_scene]` does, does `[choose_species]` with it unseen (its `[on_choose_species]` block sets the world up), fades them back in and plays
on. Retry then redoes the species last chosen, not the script's first `[choose_species]`.

Also in sandbox mode only, every book (morph, bullet, species) has a **Hold to clear data** link at
its lower left (styled as the species link). Held for 3 seconds (a red bar fills under it; let go
sooner and nothing happens), it empties the central inventory and reloads the page: a fresh sandbox
run with nothing discovered.

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
- `?dev=1` (the game): dev mode. Test aids, as keys and as buttons at the top left of the world view (**Ch**,
  **Ok**, **Do**): `I` goes to the next `[chapter]`, running every action on the lines in between
  (as `[jump_here]` does); `O` completes the waiting line's `[need]` conditions at once; `P` scores
  as if the ship ate a dot (+10 energy). The frame rate, simulation steps per second (over the last
  half second) and the measured time of a step at the top right of the world view. A grey circle around every detected object in the
  world view, besides the ship's yellow one (the one under the pointer brighter).

## The game's script (`script.txt`)

The game (`index.html`) is driven by `script.txt` (sandbox mode: `sandbox.txt`, the same syntax), read when the page loads and run line by line
(`game.js`, `ACTIONS` and `CONDITIONS`). The world starts empty with no species, so the script
sets everything up.

### Syntax

- **A line** is an optional message followed by any number of tokens: `Message text [token][token]`.
  Tokens are only read **at the end of the line**. A `[...]` in the middle of the text stays as text.
- **A token** is `[name]` or `[name:arg:arg...]`. It is either an action, `[name...]`, which runs
  when the line is reached, in order, or a condition, `[need:name...]`, which is checked every frame.
  The line waits until all of its conditions hold.
- **A line with text** shows it in the message bar under the world view. The bar is always one line
  high and never wraps: its font is 22 px (18 on a phone), or smaller if the script's widest message,
  with its ≫ or ⧗ and progress, would not fit the bar's width (measured when the script loads and on
  every resize, fullscreen included). Keep messages short; a long one makes every message smaller. If it has no conditions,
  Space or a click on the message bar moves on. If it has conditions, they move it on; Space and
  clicks do nothing then. If its conditions already hold when the line is reached (e.g. the energy
  is there already), it is shown in the waiting colour but moves on by Space or a click, ending with ≫.
- **A line of tokens only** shows nothing new, so the last message stays. It runs its actions and
  moves on at once, or when its conditions hold.
- **Text styles:** `_italic_`, `*bold*`, `*_italic bold_*`. `{<button>:<key>}` shows `<button>` as a button (in a round-cornered box, in the text's colour: cyan in a message waiting on a condition) followed by `<key>` as a keycap, e.g. `{Fire:A}`; `{<button>}` or `{:<key>}` shows just one. `\n` (a backslash and n) starts a new line in the message box. The same styles work in `[next_msg]`'s tip (a `:` inside `{...}` there is part of it, not a separator).
- **Marks:** each message ends with a grey mark. ⧗ means it is waiting on a condition. The first
  message that moves on by Space or a click gets "Press space or click message to continue ≫" on the row under it, later ones ≫, and the last line none (`[next_msg]` changes them).
- **Skipped:** blank lines, and comment lines starting with `#`.
- **Errors:** an unknown name or a bad value is reported in the browser console and ignored.
- **The end:** the last line stays on screen.
- **Game time:** the game's times are counted in simulation steps, given as *game seconds*: 1 game s = 30 steps, a real second at speed 1 (on any screen; less on a device that cannot keep up). So they are the same against the creatures on any device and at any speed, and none passes while paused: `[need:survive_time]`, a morph's discovery and snapshot (1 game s), `[need:object]` / `[need:mass]` smoothing (1 game s), `[min_morph]`'s spacing (0.5 game s), the gun's cooldown (0.7 game s). The dots drift by game time too, at the simulation's rate while playing (the same pace against the ship at any speed); while paused they keep drifting, at the speed-1 rate in real time. The star field drifts in real time. Music and the UI's animations and fades are in real time.
- **Retry:** the game's Retry button (left of Pause, `Backspace`) rebuilds the scene, e.g. a world grown too crowded: behind a quick fade it empties the bullets channel, puts the dots in afresh (the energy kept), and redoes every action from the last checkpoint, a `[chapter]` or `[new_scene]` line, up to the current line, then the current line itself (its message and conditions start afresh). The checkpoints themselves (no card, no extra fade), `[jump_here]`, `[next_msg]` and `[music]` aren't redone.

Morph ids are a name of the current species (`limus`, `gyrans`, ...) or a full pattern id
(`5N7KKM-limus`). Where a token takes several ids they may be separated by commas or colons:
`[allow_morph:limus,gyrans]` = `[allow_morph:limus:gyrans]`.

### Actions

| Token | What it does |
|---|---|
| `[species:<name>:<code>]` | Load species `<code>` (e.g. `[species:Aquarium originale:5N7KKM]`): its rule and merge saddle. The world is kept. No morph is visible after it. The morph book shows the species' name (genus and species) in italic ("Species: *Aquarium originale*"), at its lower right (`[species:<code>]`: no name, the code is shown). Each species has its own discovered morphs (a name is only a morph of its species: 5N7KKM *tardus* is not GDNQYX *tardus*), so the morph book shows this species' discoveries, kept for when it comes back. |
| `[add_allow_species:<name>:<code>[:<ids>]]` | Add species `<code>` (named `<name>`, as `[species]`) to the visible species, in order, with its visible morphs `<ids>` (separated by `,`, as `[allow_morph]`; none given: `all`), e.g. `[add_allow_species:Aquarium originale:5N7KKM:limus,gyrans,tardus]`. Again: its name and morphs are updated. `[choose_species]` uses it (and the species book, to come). |
| `[on_choose_species[:<code>]]` … `[end]` | A block of lines, taken out of the script's run: its actions run each time a species is chosen (`[choose_species]`, the species book to come), after that species' `[species]` and `[allow_morph]`. Without `<code>` it is for any species; `[on_choose_species:<code>]` is for that species only, run instead of the general one. A block holds actions only: its messages and conditions are ignored (with a warning), and its actions run at once, in order, with no waiting. |
| `[choose_species:<code>]` | As if the player chose species `<code>` in the species book: `[species:<name>:<code>]` and `[allow_morph:<ids>]` from its `[add_allow_species]` line, then the actions of its `[on_choose_species]` block. Sandbox mode starts with one. |
| `[allow_morph:<ids>]` | The visible morphs are exactly these, in this order. A morph that isn't visible is still detected but shows as "???". `all` means every morph of the species; no ids means none. |
| `[discover_morph:<ids>]` | Add these morphs to the morph book's discovered list. Ones the player already found stay. Their pictures are black until the ship is one of them for 1 game s. (A morph found in play gets its picture when found. A picture is then retaken only once after each Respawn or pick in the book, the first time the ship holds a known morph for 1 s, as the ship is then usually alone in a cleared world.) |
| `[morph_book:on\|off]` | Show or hide the Morphs button (default `on`; `X` opens the book, as that button and Choose do; only Escape closes it, the letters picking slots) and the name box of the ship's morph at the bottom of the ship view, above Upward and Zoom. `off` also closes the book, and while it is off no morph is discovered (they are still detected). |
| `[morph_book_action:<type>]` | What a click on a discovered morph in the book does (only while the ship is lost; otherwise the book is read-only). Keys pick slots too, the book's rows on the keyboard's: `A` `S` `D` `F` the first row, `Z` `X` `C` `V` the second (`X` picks its slot rather than closing the book, while that slot can be picked) (holding Shift or Alt shows a pickable slot's letter at its top right, as the other key badges). `clear_put` = `[clear_world][put_player:<id>]` (default), `put` = `[put_player:<id>:random]`, `none` = nothing. |
| `[retry_action:<type>]` | What the Respawn button (and `Z`) does: `clear_put`, `put` or `none` as above, with the player's morph, or `book` to open the morph book. No type: back to the default, `clear_put`. |
| `[clear_world[:<bg>][:empty]]` | Empty the world, and turn `[min_morph]` off (set it again after). The dots go too, and none are put back until a `[dot:<n>]` (the energy, the score, is kept); so do the bullets in flight (the gun and its list are kept). `<bg>`: the background behind the board, `black` (the default) or `radiation` (dark orange, slowly darkening and brightening again, every 8 s). `empty`: a scene of empty space, no player: while there is no ship, the ship view shows the centre of the world instead of "Dissolved" (no Respawn or Choose either), until the next `[clear_world]`. Either option alone or both, e.g. `[clear_world:radiation:empty]`. |
| `[put_player:<id>[:<how>]]` | Put morph `<id>` in as the player's ship, which is then followed. `<id>` becomes the player's morph, which Respawn puts back. `first` is the species' first morph (its first pattern). |
| `[put_morph:<id>[:<how>]]` | Put morph `<id>` in. It is never followed: the ship stays followed, or if the ship is lost, it stays lost (only `[put_player]` and the morph book's pick make a ship). |
| `[min_morph:<n>[:<ids>[:<how>]]]` | Keep at least `n` creatures besides the ship on the board. Whenever there are fewer, put in one more, picked at random from `<ids>` (default, or `all`: the visible list), at most one every 0.5 game s (none while paused). `<how>` defaults to `random`. The creatures it puts in are never followed, as `[put_morph]`. `[min_morph:0]` or `[min_morph:off]` stops it. It keeps running across later lines (and Respawn) until then, or until a `[clear_world]`, which turns it off too. |
| `[dot:<n>\|off]` | The dots channel (`dots.js`): keep `n` dots ("energy") on the board, which drift slowly to the left (or as `[stars]` sets). Any creature eats them; the ship eating one scores 10. The score shows at the top of the world view, centred (in fullscreen too), and starts from 0 when the dots are turned on. `off`: no dots, no score. |
| `[gun[:<options>]]` | The player's gun and its bullets channel (`bullets.js`): the player fires the bullet pattern in use (see `[bullet]`) with `A` or the round Fire button (lower left of the world view, shown while the bullets are on, dimmed while it can't fire), put just in front of the ship (every pattern the same distance: its back edge 8 cells from the ship's edge, its centre of mass on the line ahead) and flying the way the ship moves; a ship standing still (no arrow) fires at a random angle, a new one each shot. Always from the front, except *Asterium* and *Tetralapillium*, which barely move (the ship would drive through them at once): they are dropped behind the ship like bombs, the same distance away. *Helicium*, which stays put too, is still fired ahead, to run into: it changes the ship's shape (a way to find new morphs). Where a bullet overlaps a creature it hides those cells from the creature's own kernels (or, with `erase`, takes their mass out), and the creature's mass wears the bullet down, both gradually, so a bullet sinks some way into a body. A ship faster than its bullets runs into them and is hurt (part of the game: a fast morph is a poor gunship). At most one shot every 0.7 game s (`?bcool=`; the cooldown doesn't pass while paused); if the field holds more than 8 bullets' worth of mass, it is emptied. Once it is (almost) empty, under 1% of a bullet's mass, the rest is wiped and the field goes idle: no bullet steps and no bullet path in the creatures' step until the next shot, so the gun costs nothing while no bullet flies. Options, in any order: `on` (the default) or `off` (no gun); `mask` (hide the creature's cells from its own kernels, softly, as the mask tool; the default) or `erase` (take mass out instead); `color_now` (the channel picked with the R G B patches at the top of the ship view, in its colour; changing the pick switches every bullet already flying too; the default), `rgb` (all channels; white bullets), or `r`, `g` or `b` (that channel only, in its colour). `wear_color` (only the creature channels the bullet acts on wear it down, the default) or `wear_all` (all three channels wear it down, whatever its colour). Changing the pattern empties the bullets channel (the bullets in flight would run under the new rule). Tuning in the URL: `?bhit=0.2` (mass taken per step), `?bwear=0.03` (bullet worn per bullet step), `?bsub=3` (bullet steps per sim step), `?bcool=0.7` (the least time between shots, in seconds). |
| `[bullet:<names>[:<name>]]` | The bullet patterns the gun can fire, a list that starts empty (no list, no shots). `<names>`, in order, separated by `,`: a pattern name from the list below adds it (case and spaces don't matter, e.g. `[bullet:scutium fluens]`); `random` adds one not listed yet, at random; `all` adds every pattern; `clear` empties the list. The gun switches to `<name>` if given (added to the list if not in it, e.g. `[bullet:paraptera:paraptera]` adds and equips *paraptera*), else to the last pattern added; if none was (e.g. only `clear`), the one in use stays if it is still listed, else it is the list's first. The patterns it adds are shown: the bullet book opens and reveals them as new cards (one after another), and the equipped one shows as in use (a lighter background) once all are revealed (not when the line is redone by Retry, `[jump_here]` or dev mode's `I`, nor while a book or chapter card is showing). Its name shows under "Fire". With two or more, ◀ ▶ buttons above Fire (`D` / `F`) switch to the previous / next one, which empties the bullets channel. While the gun is on (and the ship is not lost), a **Bullets** button at the top left of the ship view (mirroring Morphs at its top right) or `S` opens the bullet book (only Escape closes it), "Use a bullet type": every pattern by name (with a small white picture of its stored pattern), 4 to a row (left to right, then top to bottom), those in the list clickable to switch to them (the one in use with a yellow border), the others "?". Keys pick slots too, the book's rows on the keyboard's: `1`–`4` the first row, `Q`–`R` the second, `A`–`F` the third, `Z` `X` the fourth. |
| `[bullet_limit[:<x>]]` | The bullets channel is emptied at once when it holds more than `x` bullets' worth of mass (one bullet: a fresh shot of the pattern in use); bullets that run into each other merge and grow. Default 8 (also with no `<x>`). |
| `[new_bullet:<n>\|off]` | The gun may gain one bullet pattern for every `n` energy gathered this run (the dots score, counted from 0 at the run's start, not from this line), whether the script gave it (`[bullet]`) or the player created it: at 310 energy with `n` = 50, six; if two have been gained this run, four more can be created. Only the patterns the run started with (the sandbox's inventory, `?all=1`) don't count. E.g. in story mode, 30 energy, then *Orbium* given, then `[new_bullet:30]`: the next one at 60. While one can, the **Bullets** button reads **Create**, in green (border and text), with a glint as it turns so. The bullet book is then titled "Create or use a bullet type" ("Create" in green, "use" in cyan, as the slots that do them), and its "?" slots have a green border: a click creates that pattern and switches the gun to it, revealed as a new card, then shown as in use; the book stays open. Slots: cyan border, a pattern to pick, the one in use with a lighter background; green border, a "?" to create. On closing, the button stays **Create** while more can be created. `off`: no more. |
| `[merge:<x>]` | Set the detection's merge saddle to `x`, the lab's Merge saddle slider. Higher means fewer merges and smaller objects; above 1, nothing merges. No `x`: back to the species's own value. `[species]` also resets it. |
| `[chapter:<no>:<title>]` | A title card over the game window: "Chapter `<no>`" above `<title>` (the title may contain colons); `[chapter:final:<title>]` shows "FINAL CHAPTER" instead. It clears the message bar. The world view and the ship view's picture fade out first, then the card fades in; Space or a click anywhere closes it, and that press does nothing else; the card fades out as the views fade back in. The game plays on meanwhile, and the line waits until the card is closed. |
| `[upward:on\|off]` | Turn the ship view's Upward (heading up) on or off, as the Upward button and `E`. |
| `[stars:<angle>:<speed>:<bright>]` | The star field: it drifts at `angle` degrees, clockwise from left (0 = to the left, the default; 90 = up; 180 = right), at `speed` times the default speed and `bright` times the default brightness (1 = default; 0 stops it or hides it). A part left out or empty keeps its default, e.g. `[stars:90]`, `[stars::2]`. Changing the angle or speed moves the stars on from where they are. The dots (`[dot]`) drift the same way, at the same speed factor. |
| `[speed:<x>]` | Set the game speed, as the Speed slider: `x` is `n` or `1/n` (`1/8`, `1/4`, `1/2`, `1`, `2`, `4`, `8`; others snap to the nearest). |
| `[music:<list>[:<play>]]` | Replace the background playlist with `<list>`: track numbers separated by commas (1 = the first in `page.js` TRACKS: 1 *Long Way*, 2 *A Few Jumps Away*, 3 *Out of the Skies, Under the Earth*, 4 *Billions and Billions*), played in turn and looped, 10 s of silence between plays. `<play>`: switch to that track now, the one playing fading out first (1.5 s) (if the music is on; else it is the one that plays when the player turns it on); it isn't restarted if it is already playing. No `<play>`: the track playing now carries on, and the list follows it (from the track after it if it is in the list, else from the list's first). E.g. `[music:2,1:2]`. Track 1 alone until a `[music]` says otherwise. |
| `[clone:on\|off]` | Allow clone mode (`on`, the default; the Morphs button gets a pink border and reads "Clone") or not (`off`: the morph book stays read-only while the ship is alive). While it is allowed and the ship is alive, the morph book is titled "Choose a morph to clone" and its discovered morphs can be picked (pink borders). The picked morph then follows the pointer over the ship view, in place of the tool circle; ↻ and ↔ buttons below Morphs turn it 45° clockwise and flip it left to right (as seen), and a click on the ship view pastes it there. With touch (no pointer to follow), a tap places the copy instead (or a drag, the copy following the finger), shown where it would go; it can still be turned, flipped or moved by another tap, and a ✓ button (shown once it is placed) pastes it. The buttons below Morphs, two to a row: ↻ ↔ / ✓ 🗙. The game pauses from the pick until the paste (or cancel), then plays on if it was playing; the copy keeps following the pointer over ↻ ↔. Pasting, or cancelling with the 🗙 button right of ↔, ends that clone (Morphs meanwhile opens the book to pick another morph to clone instead); clone mode stays allowed (clone again from the book) until `[clone:off]`. While the ship is lost, the book is "Choose a morph to respawn" as usual. The player **boards the clone**: once pasted, it is the followed ship (the old one goes on as another creature; a clone pasted onto the ship makes one object with it, followed as before), and its morph is what Respawn puts back. So a morph can be changed without dissolving, and a twin engineered by becoming one half (a clone away from the ship) and pasting the other half next to it. |
| `[erase:on\|off]` | Allow erase ("recycle": the selected R G B patch clicked again, see the top of this README) or not; `off` also goes back to masking. Allowed by default. |
| `[blank]` | The line is an empty message: the message bar shows only its mark (≫, or ⧗ with a condition) and waits for Space or a click, like any message. |
| `[next_msg[:<tip>[:<next>[:<act>]]]]` | Set the grey tip shown once, after the text on the same line, on the next message that moves on by Space or a click (not when the line is skipped over by `[jump_here]` or dev mode's `I`, so the tip goes only on the script's first such message) (default "Press space or click message to continue"), the mark ending such messages (default ≫), and the mark ending a message that waits on a condition (default ⧗). A part left out or empty goes back to its default. The tip can't contain `]`, nor `:` outside `{...}`. |
| `[new_scene]` | Fade the world view and the ship view's picture to black over 0.5 s (everything drawn in them, "Dissolved" included, and the Respawn and Choose buttons; not the other controls), and wait for it. The actions after it set the scene unseen, and the next message fades the views back in as it shows. The game plays on meanwhile. |
| `[jump_here]` | For writing the script: start at the first line that has it. Every action on the lines before it runs first, in order, so the set-up is as if the game had been played up to there; their messages and conditions are skipped, and so are the chapter cards. |

`<how>` is where a creature goes: `center` (the board centre, unrotated; the default for
`put_player` and `put_morph`), `random` (a random empty spot and angle, as the lab's Add) or `far`
(the empty spot farthest from the ship, at a random angle).

### Bullet patterns

The patterns `[bullet:<name>]` adds to the gun (`bullet_patterns.js`), all single-channel classic Lenia from `~/github/Lenia/Python/animals.json` (the entry with that code; #2, #3 = its second, third entry), T = 10. R is the pattern's native radius; every bullet runs at R 13 (with the default creature R 18; scaled with `?R=`), resampled by 13 / R. Speed: cells per sim step at R 13 and the default `?bsub=3`, measured on a CPU copy (numpy); launched at 0, 90, 200 and 315 degrees, the gliders held their mass and heading except where noted. The ids and names used before (`scutium2`–`scutium5`, `orbium1`, `paraptera1`, `Scutium Mk II`–`IV`) still work.

| Name | Species | Code | Native R | b | Speed |
|---|---|---|---|---|---|
| `Orbium` | Orbium unicaudatus | O2u | 13 | 1 | 3.1; dies at some launch angles |
| `Synorbium` | Synorbium solidus | O4s | 13 | 1 | 3.2 |
| `Scutium` | Scutium valvatus | S1v | 13 | 1 | 2.1 |
| `Scutium fluens` | Scutium serratus fluens | 2S1f | 15 | 1/2,1 | 2.2 |
| `Scutium vagus` | Scutium serratus vagus | 2S1v | 10 | 1/2,1 | 1.9 |
| `Scutium serratus` | Scutium serratus fluens | 2S1f #2 | 18 | 1/6,1/2,1 | 3.0 |
| `Triscutium` | Triscutium valvatus laterale | S3vl | 15 | 1 | 2.2 |
| `Paraptera` | Paraptera arcus pedes | P4ap | 13 | 1 | 2.4 |
| `Heptapteryx` | Heptapteryx arcus labens | P7al | 13 | 1 | 2.3 |
| `Helicium` | Helicium cavus pedes | H3cp | 13 | 1 | spins in place |
| `Geminium` | Hydrogeminium natans | 3GH2n | 18 | 1/2,1,2/3 | 1.5 |
| `Kronium` | Kronium vagus | 3K4v #2 | 27 | 2/3,1/3,1 | 1.1; wanders |
| `Asterium` | Asterium vagus | 3R5v #3 | 27 | 1,1/2,1 | 0.4; wanders, heading unpredictable; dropped behind |
| `Tetralapillium` | Tetralapillium metamorpha | 2L4m | 18 | 1,1/2 | changes shape, creeps; dropped behind |

### Conditions

When a waiting line's conditions are all met, a green **SUCCESS** box fades in in the middle of the game window and slowly fades out (3 s; the game goes on meanwhile), and the script moves on. Also for a line whose conditions already hold when it is reached (it then shows ≫ and moves on by Space / a click).

| Token | Holds when |
|---|---|
| `[need:first_touch]` | The player has pressed on the ship, and its detected morph is now not *limus*. The press is required, since a fresh ship can read as unknown at first. |
| `[need:discover_morph[:<ids>[:any\|all]]]` | One of `<ids>` (`any`, the default), or every one of them (`all`), is in the morph book: discovered earlier (also put in by `[discover_morph:<ids>]`), or discovered now and the book that opened for it closed. With no ids: the player discovers a new morph, any, after the line is reached, and closes the book that opened for it. |
| `[need:create_bullet]` | The player has created a bullet type since the line was reached: clicked **Create** (the Bullets button) and a "?" in the bullet book. `[new_bullet]` must allow one (enough energy), or the line waits until it does. |
| `[need:respawn]` | The ship has dissolved and the player has put it back since the line was reached: by **Respawn** (`Z`), or **Choose** (`X`) and a pick in the morph book. |
| `[need:first_clone]` | The player has placed a clone since the line was reached: picked a morph in the morph book while clone mode is allowed (`[clone:on]`) and pasted it on the ship view (a click, or ✓ after a touch). Cancelling doesn't count. |
| `[need:recycle]` | The player has turned erase (recycle) on, in any channel, and erased some mass on the ship view since the line was reached, and the ship has then dissolved. |
| `[need:object:<n>[:exact\|min\|max]]` | The world holds exactly `n` objects (`exact`, the default), at least `n` (`min`) or at most `n` (`max`): all the creatures detected in it, the ship included. The count is smoothed, the median over the last game second (one sample per step), so a creature flickering between one and two objects doesn't make it jump. The ⧗ mark is followed by the progress in grey, "⧗ morph: 20/60". |
| `[need:mass:<r\|g\|b\|rgb>:<x>[:min\|max]]` | The creatures detected in the world, the ship included, hold at least `x` g of mass (`min`, the default) or at most `x` (`max`) in that channel, or in all three (`rgb`); g as the stats charts count it (the sum of cell values / R²). Smoothed as `[need:object]` (the median over the last second), with its progress after ⧗ to one decimal, "⧗ mass: 12.3/20". |
| `[need:survive_time:<s>]` | The ship has stayed alive for `s` game seconds from when the line is reached (paused time, an open book and a chapter card don't count); dissolving stops the count (it stays shown), and it starts again from 0 once the ship is back (Respawn). Progress after ⧗ (as `[need:object]`) in whole seconds, "⧗ time: 12/30". |
| `[need:upward]` | The ship view's Upward is on. |
| `[need:energy:<n>]` | The energy, the dots score, has reached at least `n`. |

### Example

```
[chapter:1:Metamorphosis]
[species:Aquarium originale:5N7KKM][allow_morph:limus][discover_morph:limus]
[morph_book:on][clear_world][put_player:limus]
Welcome to the Cognoship game.
Click on the ship to change its shape.[need:first_touch]
[dot:3][min_morph:2:gyrans,tardus:far]
Gain 50 energy. [need:energy:50]
```

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
(off until first turned on, then as last chosen in this browser; each track is downloaded only when it starts to play). Browsers only start
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
