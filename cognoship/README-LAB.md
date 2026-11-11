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
- **Reset** (`Z`): clear the board and put the pattern in the centre. **Add** (`X`): add it at a
  random empty spot and angle, keeping the board.
- **Speed** (`Q` faster, `W` slower): 1/8 to 8, 30 steps a second at 1, by the clock.
- **Zoom** of the ship view. **Heading up** (`E`): turn the ship view with the direction of travel.
- **Species** and **Morph** pickers (`C` `V` step morphs; `Shift`+`C` / `V` step species).
- Click a creature in the world view to follow it. The run stops as *Dissolved* when none is left.
- **Detect** names a creature's morph: None, Focused (name box), Every (each creature, in its hue).
- **Stars** turns the star field and the music on or off.
- Hold **Shift** or **Alt** to show every control's key.

## Modes

**Mode** (`A`, `Shift`+`A` backwards) cycles Action, Stats, Action+Stats and Tuning.

- **Stats:** a 2D / 3D phase-space chart of the followed creature (mass, R G B masses, hue, speed,
  turn rate) in place of the ship view, or beside it in Action+Stats. Axis keys `D` `F` `G` `H`.
- **Record** (`S`): add the next 300 steps to the species' long history. **Delete** removes the
  last piece. **Download** saves it as `<species>.json`; put it in `assets/stats/` and it loads
  with the species, and morph detection learns from it.
- **Tuning:** the detection sliders (merge saddle, threshold, min area), with every detected
  object circled.

## URL parameters

`?form=` (start pattern, e.g. `5N7KKM-limus`), `?R=` (kernel radius, default 18), `?dt=`,
`?speed=`, `?zoom=`, `?radius=`, `?headup=vel|axis`, `?debug=1` (GPU ms per step),
`?saddle=` `?thr=` `?area=` (detection). Details in [README-MORE.md](README-MORE.md).
