# Cognoship lab

`lab.html`: the experiment page behind the game, with every control. Steer an Aquarium (a
3-channel Lenia creature from the aquarium collection, run with LeniaF.py's rule) by **cognitive
masking** (Cool et al. 2026): while you hold the pointer on the magnified creature, that spot of
one channel drops out of every kernel that reads it, U = K∗(A·(1−B)) / K∗(1−B). The same circle
can instead erase or draw mass, for comparison.

- [README.md](README.md): the game, for players.
- [README-MORE.md](README-MORE.md): the full reference (every control in detail, the species and
  morph list, how it works, fidelity checks, stats, morph detection, tuning).

**Run:** `cd Cognoship && python3 -m http.server 8000`, then open http://localhost:8000/lab.html.

## Controls

- **Press and hold on the magnified creature** (right panel): apply the tool to a circle there,
  kept at the same place relative to the creature while you hold; drag to move it.
- **Tool patches:** the top row picks the channel, **R G B** (`1` `2` `3`); the bottom row the
  action (`4` `5` `6`): **Mask** (dashed: hidden from the kernels, no mass touched; G mask
  steers best), **Erase** (solid: set to 0 while held), **Draw** (dotted: set to 1).
- **Pause** (`Enter`) · **Step** (`Space`) · **Mask radius**.
- **Reset** (`Z`): clear the board and put the pattern in the centre. **Spray** (`X`): add it at a
  random empty spot and angle, keeping the board. **Add** (`B`): hold it to paste where you click,
  in the world view or the ship view; ↻ ↔ (`N` `M`) turn and flip it, `Enter` pastes, `Esc` cancels.
- **Speed** (`Q` faster, `W` slower): 1/8 to 8, 30 steps a second at 1, by the clock.
- **Zoom** of the ship view. **Heading up** (`E`): turn the ship view with the ship's morph model's frame when it has one (heading, or body axis for a still morph: its action spots stay put on screen), else with the direction of travel (a still morph, e.g. antigyrans or lithos: its body axis, B's CoM to G's, as the models' frame).
- **Species** and **Morph** pickers (`C` `V` step morphs; `Shift`+`C` / `V` step species). Morphs in `[add_allow_species]` order (core.js ALLOWED_MORPHS, kept as sandbox.txt), the others last.
- Click a creature in the world view to follow it. The run stops as *Dissolved* when none is left.
- **Detect** names a creature's morph: None, Focused (name box), Every (each creature, in its hue).
- **Stars** turns the star field and the music on or off.
- Hold **Shift** or **Alt** to show every control's key.

## Modes

**Mode** (`A`, `Shift`+`A` backwards) cycles Action, Stats, Action+Stats, Model and Tuning.

- **Stats:** a 2D / 3D phase-space chart of the followed creature (mass, R G B masses, hue, speed,
  turn rate) in place of the ship view, or beside it in Action+Stats. Axis keys `G` `H` `J` `K`.
- **Record** (`S`), **Record Whole** (`D`), **Record Fuzzy** (`F`): add the next 300 steps to the
  species' long history: of the followed creature, of the twin or trio it is part of, or of it as a
  fuzzy (loosely matched) morph. **Delete** removes the
  last piece. **Download** saves it as `<species>.json`; put it in `assets/stats/` and it loads
  with the species, and morph detection learns from it.
- Wholes (Record Whole; see Tuning's twin saddle) are only compared with wholes, parts with parts.
- **Model:** the species' morph model (`models/<code>.json`, see [models/README.md](models/README.md))
  drawn on the followed ship while the name box names one of its morphs (nothing while it shows
  ???): at each spot tried, what a G mask held there for the chosen **Hold** (0.5 / 1 / 2 / 3 s,
  the buttons just above Mask Erase Draw) did. An arrow: it stays that morph, pointing to its new
  heading; a square: it becomes the lettered morph (grey: dissolves); blue: it keeps its turning
  direction (or has none), orange: it turns the other way after; green: it becomes a twin or
  triple; magenta: it splits (LL) or replicates (T+); faded: under 3 of 4 trials. A ship turning −
  sees the model mirrored. An oscillator swept over its cycle, where the phase changes the actions
  (HAESRE zigzag, Z18A9R transgyrans, GDNQYX rectus, HGMLH8 paratardus): the actions, the Maneuver
  menu and the details show the actions as they work in its phase now (12 phases, told from a sine
  fitted to its last cycle of turn rates, or speeds for rectus; the info line shows it). The
  Maneuver menu and the graph list what it can do in any phase: an action for one phase ("at phase
  9/12") waits for that phase, then is done (the info line shows the wait). A replication through a
  short-lived morph (VT049W: gyrans) is drawn as the action's arrow to it (blue from a single,
  magenta from a twin: one gyrans or two) and that morph's green arrow to the copies. VT049W
  gyrans itself (short-lived: it splits about 115 steps after forming) is modelled by age (60-110
  steps, from when it is first named): young, many actions turn it into paragyrans or rectogyrans;
  older, it only replicates. A still spinner spinning the other way from the model's (antigyrans)
  sees it mirrored, as a mover turning the other way does. Chaotic
  morphs (HGMLH8 metamorpha, KH97WU impetus) were swept by state, but their actions vary as much
  within a state, so their pooled actions are shown. A model with two mask channels (XEH4YR: G and B) shows the selected
  tool channel's actions (`1` `2` `3`), else its first. Hover one for the details under the screen; click it to
  do it (the mask held there for the hold, as if by the player). Entering Model turns Heading up on.
  **Graph** shows the species' morph graph over the middle of the world view: an arrow for each
  morph a morph can become (a loop: it turns the other way in place; grey: dissolves), by its
  best action (hover the arrow for the action and how often it works; it is highlighted); blue
  keeps the turning direction, orange turns the other way, green turns a single into a twin or
  triple, magenta splits it or replicates it (copies of one morph, e.g. PJ1EF7 rectus into
  several tachus through vagus), dashed under 3 of 4 trials (but at least half). Twins and
  triples are rectangles, dissolved a diamond. Under a name: ↻ / ↺ rotating (a mover's turning; a still one's spin, e.g. antigyrans),
  ∿ oscillating (its heading swings, e.g. HAESRE zigzag, or its mass pulses). The followed ship's morph is ringed yellow; its arrows are thicker, and a
  click on one does it (a fuzzy name, ~rectus, counts as rectus). Morphs the game doesn't detect
  (not in sandbox.txt's `[add_allow_species]`, e.g. lithos+gyrans) are dimmed, with their arrows.
  A click on a morph resets the world with it. On at first in Model mode. With more than one mask
  channel modelled, its key (`mask channels modelled: all · G · R`) picks the graph: all (each
  arrow by its best action of any channel) or one channel's (the same arrows, each by that
  channel's best action: dimmed where another channel's does better, faint and not clickable where
  no action of that channel does it).
  **Maneuver** lists the outcomes the ship can reach now (small / big turn right / left, keeping its
  turning direction or turning the other way after; into another morph; dissolve; only those
  with an action that gives it in 3 of 4 trials); picking one does its best action at once: the
  most reliable, then the name box kept, the shortest hold, the shortest time unnamed
  (`models.js`, for the game's bots and autopilot too). A turn is small under 45°.
- **Tuning:** the detection sliders (merge saddle, threshold, min area, twin saddle), with every
  detected object circled and every whole dashed.

## URL parameters

`?form=` (start pattern, e.g. `5N7KKM-limus`), `?R=` (kernel radius, default 18), `?dt=`,
`?speed=`, `?zoom=`, `?radius=`, `?headup=vel|axis`, `?debug=1` (GPU ms per step),
`?saddle=` `?thr=` `?area=` (detection), `?fft=0` (the direct convolution instead of the FFT, the default). Details in [README-MORE.md](README-MORE.md).
