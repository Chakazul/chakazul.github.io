# Morph models

Levels of use: an **action** (a mask: colour, place, hold, phase; `actions`) → a **maneuver** (the best action for an outcome: models.js `maneuvers`, `choose`) → a **composite maneuver** (maneuvers in sequence, tried in experiments: `composites`; tools/models.py `composite`) → a **role** with its moods → a **behavior** (how much of each role and mood: behaviors.js; roles and moods are planned, PLAN-BEHAVIORS.md).

What an action does to a morph, found by experiment: which mask (channel, spot, how long) turns a
ship or changes it into another morph. Twins (two bodies named as one, e.g. paragyrans) are
modelled the same way, as ships of their own (making one by cloning is the player's alone). Meant for game AI (bots steering at the player and firing along their heading, autopilot,
transforming at will) and for reading.

Written by `tools/models.py` (its docstring has the details). It runs the game's rule headless
on the GPU from settled copies of a morph at every phase of its own cycle and 4 start angles,
tries each action on each, and follows the game's name box (`detect.js`, ported) through it.

Actions are kept simple and nameable: **mask one spot for a while**, G unless a species needs
another channel. Spots: the centre, or one of 8 directions in the ship's frame (front,
front-right, right, ..., front-left) on the body (0.45 R from the centre), at its edge (0.9 R) or just outside it (1.3 R).
Holds: 0.5, 1, 2 (and 3) s. A twin (named as one whole, e.g. paragyrans) is about 1.6 times as wide, so
its spots are 1.6 times as far out (`--spot-scale 1.6`). (More spots or holds can be added where a book is too coarse.)

- `<code>.md`: the readable book, per morph a menu of what can be done and how:
  **Turn** (keeps the morph: every action that does, from left to right, so a bot can pick a
  small or a big turn), **Change direction** (keeps the morph, turns the other way: its mirror
  image), **Into X** (another morph, turning + or −), **Dissolve**, **Seen in passing** (morphs
  named for a while on the way), **Unreliable** (actions with no outcome in 3 of 4 trials). Each
  row has the *name box* sequence the game would show (e.g. `limus > ??? > limus`, or *always
  named*), the share of trials `p`, and, only for a morph whose phase matters, the phases it
  works at. Speed is not measured: each morph has its own speed; what counts is its heading and
  which morph it is.
- `<code>.json`: the same as data for the game, and per morph every action tried (`actions`), which
  the lab's Model mode draws on the ship (a click does it). `models.js` reads it: the outcomes a
  ship can reach now (its named morph and turning direction) and the best action for each, for
  the lab's Outcome menu and the game.
- `<code>-<morph>.png`: a panel per hold, the creature in its frame behind (yellow: heading up,
  and its turning direction). Each spot: an **arrow** if it stays the morph, pointing to its new
  heading; a **square** if it becomes something else, lettered with the morph (G = gyrans, ...),
  grey if it dissolves. **Blue**: its turning direction stays (or it has none), **orange**: it
  turns the other way after (a morph's mirror image). Faded: the outcome holds in under 3 of 4
  trials. The lab's Model mode draws the same on the ship.

## Conventions

- **Frame:** a moving morph's heading (core.js `slowHeading`, what Upward shows as up); a still
  one's body axis (B channel's centre of mass to the G channel's).
- **Spot:** as above; as numbers `(fwd, side)` in R from the centre of mass, + fwd ahead, + side to
  the right as seen. The disc (radius 4 cells, the game's) follows the creature and turns with its heading
  while held, as the ship view does with Upward on.
- **Turn:** degrees, + clockwise as seen (right), − anticlockwise; measured once the name box names it
  again (at least 30 steps after the hold), against an untouched copy from the same start (so the
  morph's own turning is taken out).
- **Phase:** steps since the cycle's main peak of the total mass when the action starts (0 to cycle − 1) (the game can
  read it from the ship's stats). Used only where it matters: where an action has no main
  outcome, the book checks whether grouping the trials by phase predicts the outcome better than
  random groups do (by 0.1 or more). For 5N7KKM limus it does not (gain 0.00): its small 9-step
  ripple makes no difference, and the unreliable actions vary from trial to trial. A morph that
  oscillates more (a limit cycle) is expected to need the phase as one more input, in a few bins.
- **Turning direction:** the sign of the turn rate, + clockwise (detection cannot tell a morph
  from its mirror image, the turn rate can). Every morph is started turning +: a pattern that
  settles turning − is mirrored first (the book says so). A ship turning − takes the mirrored
  action: left and right swapped, turn negated.
- **Using it:** a bot decides only what is needed now (the morph, its turning direction, and the
  turn it wants) and picks one row; combinations come by themselves as it decides again later,
  so the book has no plans or chains.
- **Hold:** steps; 30 steps = 1 s at speed 1.
- **p:** the share of trials (start angles and phases) that gave the outcome named.
- **Clone:** `dist` R from the ship's centre at `bearing` (0 ahead, 90 right), the clone's body
  axis turned `turn` degrees from the ship's frame (+ clockwise), mirrored first if `flip`.
  Pasted as the game does (each cell the larger of board and clone).

Outcome names: a morph of the species; `a + b` (bodies that stay apart); `???` (the detection
names nothing: possibly a morph with no long history yet, or a new one); `dissolved`, `exploded`,
`multiplied`.

## Run

```
PY=~/github/evocosmos-drama/.venv/bin/python          # any env with torch (CUDA), scipy, matplotlib
$PY tools/models.py steer --code 5N7KKM --morph limus                 # a few minutes a morph (RTX 5070)
$PY tools/models.py steer --code XEH4YR --morph rectus --channels B
$PY tools/models.py steer --code 5N7KKM --morph paragyrans --dists body,edge,outside --spot-scale 1.6 --size 160   # a twin
$PY tools/models.py book  --code 5N7KKM                               # models/5N7KKM.md .json .png
$PY tools/models.py show  --code 5N7KKM --morph limus --spot 'right edge' --hold 30   # GIF with the name box
```

Raw results (every trial, time series, start states, GIFs) go to `output/models/<code>/`
(gitignored).

## Limits

- 128 × 128 torus (192 × 192 for clones) with nothing else on it; the game's board is bigger and
  busier (dots, bullets, other ships).
- Bodies are found as thresholded connected regions, not the game's potential basins, so a twin
  may be counted as one body or two; both the part and the whole detection are tried.
- An outcome is what it is 400+ steps after the action; a passing state that later changes again is
  named by where it is at that time.
