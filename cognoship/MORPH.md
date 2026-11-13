# Morph notes

Findings on individual morphs: new ones, renamed ones, and ones checked and not (yet) added.
Simulations here are the numpy port of the LeniaF rule in `tools/from_screenshot.py`, at R = 18
on Cognoship's 384×256 torus, the same rule and board as the game. Masses are in cells.

## Tools

- `tools/from_screenshot.py` rebuilds a pattern from a screenshot: it reads the cell grid, tries
  20 (gamma, scale) seeds, settles each for 300 steps, checks 500 more (number of creatures
  constant, mass within 1.5 %), and keeps the pass most like the screenshot (radial profiles per
  channel). `--add` appends it to `tools/extra_patterns.json`; then `python3 tools/gen_aquarium.py`.
  For recordings and passing states: a video with `--frame N`, `--crop x0,y0,x1,y1`,
  `--clean-overlays` (paints out the ship view's frame and the detection circles) and `--raw`
  (the grid as read, values ** `--gamma`, not settled or checked). Images and `pattern.json` go to
  `output/from_screenshot_<name>/`.
- To compare a candidate with a species' known morphs, run each for 300 + 600 steps and compare
  mass, speed (R per unit time), turning rate (rad per unit time), number of bodies and radial
  profile. A mirror image has the same mass, speed and profile but the opposite turning rate.

## VT049W rectogyrans (renamed, 2026-11-13)

Formerly *rectus+gyrans* (from 233.json): renamed in `tools/gen_aquarium.py` (`POST_RENAME`),
id `VT049W-rectogyrans` (first named pararectus, then rectogyrans: a rectus and a gyrans, as
GDNQYX's). Its long history is in `assets/stats/VT049W.json`, and it is in VT049W's
allow lists in `script.txt` and `sandbox.txt`.

## GDNQYX pararectus and rectogyrans (new, 2026-11-13)

From `../assets/screenshots/GDNQYX pararectus.png` and `GDNQYX rectogyrans.png`, rebuilt by
`tools/from_screenshot.py` (settled and checked): pararectus, two rectus side by side (gamma 1.4,
64×30 cells, shape 0.97); rectogyrans, a rectus joined to a gyrans (gamma 1.0, 50×54 cells, shape
0.95). Both single bodies moving straight. Against GDNQYX's morphs (R = 18, settled):

| morph | mass | speed | turn |
|---|---|---|---|
| **pararectus** | 1934 | 0.097 | 0 |
| **rectogyrans** | 1863 | 0.110 | 0.003 |
| paragyrans | 1791 | 0.113 | 0 |
| lithos+gyrans | 1897 | 0.056 | 0.021 |
| paratardus | 1981 | 0.009 | 0 |
| rectus | 885 | 0.067 | 0.003 |

Rectogyrans is closest to paragyrans (4 % heavier, as fast): its long history will show whether
the detection tells them apart. In the sandbox's allow list.

## 5N7KKM paralithos (new, 2026-11-13)

From `../assets/screenshots/5K7KKM paralithos.png` ("5K7KKM" in the file name is 5N7KKM; the
`(2)` screenshot has a lithos at its lower right for scale). A still pair of lithos-like discs
joined side by side.

- Rebuilt with gamma 1.8, scale 1.05: 48×48 cells; two bodies throughout, mass within 0.35 % over
  500 steps, shape similarity 0.95. In `tools/extra_patterns.json`, R = 18.
- Distinct from every 5N7KKM morph:

  | morph | mass | speed | turn | bodies |
  |---|---|---|---|---|
  | **paralithos** | 1929 | 0.0002 | −0.064 | 2 |
  | lithos | 913 | 0.0001 | −0.020 | 1 |
  | paratardus | 1976 | 0.0166 | 0.001 | 2 |
  | 2lithos+gyrans | 1938 | 0.0058 | 0.003 | 2 |
  | lithos+gyrans | 1680 | 0.039 | 0.056 | 1 |

  About twice a lithos's mass, and still where the two-body morphs of similar mass drift.
- Not yet in any allow list, and no long history in `assets/stats/5N7KKM.json` (the game shows
  "???" until one is recorded in the lab).

## Z18A9R "transgyrans" (checked, not added)

From `../assets/screenshots/Z18A9R transgyrans.png`, a phone screenshot of the ship view (an
antigyrans at its lower right for reference). The ship's arrow was removed by inpainting along its
line, and the "???" box cropped off.

### One crescent: a mirrored gyrans

The followed crescent alone (the arrow starts at its centre of mass) rebuilds cleanly (16 of 20
seeds pass, shape 0.98), but settles into the existing gyrans mirrored:

| | mass | speed | turn |
|---|---|---|---|
| rebuilt crescent | 868 | 0.1233 | **+0.2745** |
| gyrans | 866 | 0.1233 | **−0.2743** |

The detection does not tell mirror images apart: its features (`detect.js`) are the mass of each
channel, the speed and the turn rate without its sign (`core.js`: `Math.abs(k.omega)`), all the
same for a mirror image (the board's grid is symmetric under a flip, and these runs give equal
values). So a lone mirrored gyrans reads as gyrans. The "???" in the screenshot more likely came
from the encounter itself: the other crescent touching it changes its mass and motion, and after
a gap (the two seen as one object) the detection needs 35 samples again before it names anything.

### Two crescents: gyrans + mirrored gyrans

The followed crescent and the one beside it are mirror images of each other. Rebuilt together
from the screenshot, the pair is not reproducible: over 800 steps, seeds a little apart end as an
S-shaped pair, an arch (paragyrans-like), one gyrans, or several pieces; one in 20 passes, shape
only 0.28.

Built instead from clean parts (a settled gyrans and its mirror image, at distances 1.6, 2.0, 2.5
and 3.0 R, placement angles every 45° and relative rotations every 45°: 256 starts, 3,000 steps
each, classified over the last 1,000):

| outcome | runs |
|---|---|
| two bodies, distance rising and falling | 83 |
| merged into one body (~1.85 gyrans' mass) | 66 |
| grew / multiplied (up to ~20× the mass) | 45 |
| one gyrans left | 37 |
| two, drifting far apart | 12 |
| two at a steady distance | 9 |
| died | 4 |

Two kinds among the 83, followed for 6,000 steps:

- **Mirror pair, fixed** (e.g. 1.6 R, 45°, 135°): each circles on the spot, as mirror images of
  each other; the pair's centre does not move, the bodies never touch. A lone gyrans also circles
  on the spot (radius ~0.45 R), so this may be two neighbours rather than a bound twin.
- **Intercepting pair** (1.6 R, 180°, 135°): the two keep running into each other, merged into one
  body ~14 % of the time, the pair's centre wandering slowly; mass exactly 2 gyrans' throughout.
  This matches the screenshot.

Robustness of the intercepting pair: its state at step 3,000, as it is and with 2 % noise on
every cell (8 copies), run 10,000 more steps each (13,000 in all, 6,500 units of time):

| run | one body | two bodies | mass (gyrans) | end |
|---|---|---|---|---|
| original, noise 1, 2, 3, 7, 8 | 12–13 % | 87–88 % | 1.85–2.17 | 2.00–2.01, intercepting |
| noise 6 | 62 % | 38 % | 1.82–3.88 | 1.83, merged into one S-shaped body |
| noise 4 | 4 % | 26 % | up to 4.0 | 3.0: a third gyrans split off |
| noise 5 | 6 % | 44 % | up to 13.6 | 8.7: multiplied |

The six that hold are frame for frame the same after 2,000 steps: the noise dies out and they
return to one orbit. So the intercepting pair looks like an attractor, a lasting twin morph rather
than a transient, though a fragile one: a third of the nudged copies fall out of it, into the
merged body or into replication.

### Tried as a pattern, removed (2026-11-13)

Rebuilt by `output/transgyrans/make_twin.py` (kept locally, not in git; it writes its images and
states beside it): the intercepting start above, its state at step 3,000 is the pattern (71×66
cells, mass 1742, R = 18), in `tools/extra_patterns.json`. The stored pattern, decoded back from
its RLE, still intercepts over 6,000 steps: one body 13 %, two 87 %, mass 1.85–2.03 gyrans, ending
at 2.01. Added as Z18A9R transgyrans, then removed: watched in the game, it is not what the
screenshot showed. The two only touch a little now and then; they do not run into each other.
So the twin seen in the game is still not reproduced: either a different start, or a transient.

### From the video (2026-11-13)

`../assets/screenshots/Z18A9R transgyrans.mp4` (the game's world view) shows the crossing: two
mirror gyrans rolling over each other again and again, until after a few encounters they become a
gyrans and a lithos. Its frames were read into cell grids (`output/transgyrans/video/extract.py`, since folded into
the tool: `--frame N --clean-overlays --raw`; 4.85 px per cell). Under the rule, the grids reproduce the
crossing, and it is passing there too: over 16 frames, a median of ~500 steps (at most ~3,000),
ending as antigyrans or paragyrans, or with one gyrans gone (never the video's gyrans and lithos
in these runs). Frame 1 (values ** 1.8, the longest-lived, 2,580 steps) is added as
**Z18A9R transgyrans**, unsettled, and in the sandbox's allow list: it crosses, then decays.
A search for small rule changes that make the crossing last, with the other Z18A9R morphs
unchanged, is in `output/transgyrans/tune/`.

### A rule tweak that makes it stable (2026-11-13, applied)

Searched in `output/transgyrans/tune/` (a 192×192 torus, float32 FFTs, the same rule): each of the
45 kernel parameters (`m`, `s`, `h` of 15 kernels) by ±3 %, and smaller time steps, one at a time.
Scored by how long the crossing lasts from 16 video frames (up to 3,000 steps; it ends where the
total mass stops swinging), and by how much gyrans, paragyrans, antigyrans, trias, lithos and
tardus change (settled mass within 2 %, speed and turn rate within 15 %).

- Unchanged rule: the crossing lasts a median of 230 steps (mean 439), none to 3,000.
- Six single changes make it last and keep the other morphs within bounds; the best is
  **kernel 3 (G→G) `m` lowered** (0.123): from ×0.99 the crossing becomes an attractor. At ×0.985
  (`m` 0.1212), 5 of 6 seeds hold it for 9,000 steps, all on the same cycle (mass 1.95 ± 0.056
  gyrans, the two meeting and parting); the sixth lost a gyrans at the start.
- What it costs the other morphs (×0.985): mass gyrans −0.43 %, paragyrans −0.13 %, antigyrans
  −0.21 %, trias −0.16 %, lithos +0.20 %, tardus −0.14 %; speeds within ±1 % (antigyrans and
  lithos, near still, within 0.001 R per unit time); gyrans' turn rate 0.274 → 0.266 (−3 %).
  Within the detection's widths (mass 3 %), so the long histories should still hold.
- Not checked: in the game itself (GPU, float32), and other Z18A9R behaviour than these six.

**Applied** as kernel 3 `m` = 0.121 (×0.984, between the ×0.98 and ×0.985 tested), in
`tools/rule_tweaks.json`, which `tools/gen_aquarium.py` applies when it writes `aquarium233.js` and
`tools/from_screenshot.py` when it simulates. **Z18A9R transgyrans** is now frame 1 of the video run
3,000 steps under it (`output/transgyrans/stable/make.py`): 42×59 cells, mass 1664; decoded back
and run 6,000 more steps it holds the crossing throughout (mass 1.95 ± 0.056 gyrans), repeating
about every 2,000 steps. The long histories in `assets/stats/Z18A9R.json` were recorded under the
old rule (the other morphs moved under 0.5 % in mass).

### Open

- Reproducing the twin as seen: two gyrans that really run into each other, again and again.
