# Plan: species behaviors and the new story chapter

A plan for the next session. The work so far has built four levels of model use; this plan uses
them to give each Aquarium species its own bot behaviors, then introduces them in a new story
chapter.

| Level | What | Where |
|---|---|---|
| action | one mask: colour, place, hold (and phase / age / side where they matter) | `models/<code>.json` → `morphs.<m>.actions`, `push_table` |
| maneuver | the best action for an outcome (reliability less waiting) | `models.js` `maneuvers`, `choose`, `stepAction`; `behaviors.js` `turnAction`, `reach`, `push` |
| composite maneuver (once strategy) | maneuvers in sequence, tried in experiments | `models/<code>.json` → `composites`; `tools/models.py composite`; `behaviors.js` `s.run`; the lab's Composite menu |
| behavior | when to use which composite (→ roles with moods, and behaviors as their shares, below) | `behaviors.js` `behavior({...})` (`balanced`, `aggressive`: 5N7KKM only so far) |

Tools: `tools/models.py` (sweeps `steer`, `book`, `composite`, `combo`, `props`, `cycles`, `show`);
data in `output/models/<code>/` (gitignored); Python `~/github/evocosmos-drama/.venv/bin/python`.

## Step 0: finish the B and R mask run (run 16)

Started 05:24 on 2026-11-22, about 13 h (estimated end ~18:30). It sweeps every modelled morph of
the 10 Aquarium species (not G6G6CR) with B and R masks (XEH4YR: R only), `steer --like` (the same
starts as each G sweep), and rebuilds each species' model as it goes.

1. Check: `grep -c "^== " output/models/run16.log` (53 sweeps + `== done`), `grep "!!" output/models/run16.log`.
2. Rerun any failed sweep by its line in `output/models/run16.sh`, then `book --code <code>`.
3. Check a model has all three: in the lab's Model mode the graph's line reads `mask channels modelled: G · R · B`.
4. Commit `models/` (and `tools/models.py` if changed).
5. Look for new useful outcomes from B / R masks (e.g. single → twin, replicate, transforms the G
   masks cannot do): list per species, actions with p ≥ 0.75 whose outcome no G action reaches.

## Step 1: rerun the composites that retry

A runner bug failed boards early on a try with no action from what they had become (fixed in
`c63338b`). Rerun the composites with `xN` steps so their numbers are right, and any whose best
action is now B or R (after step 0): the `S ...` lines in `output/models/run12.sh` … `run15.sh`
(names now in words, e.g. `--name 'lithos activates'`). Then `book` the species and commit.

## Levels, redesigned (2026-11-24, the user's)

| Level | Scope | What | Where now → planned |
|---|---|---|---|
| action | condition-specific (morph, phase / age / side) | the muscle control: one mask (colour, place, hold) | `models/<code>.json` actions (unchanged) |
| maneuver | morph-specific | the best action for a goal (into X, turn, split, replicate) | `models.js` `maneuvers` (unchanged) |
| composite maneuver | morph-combo-specific | a series of actions for a goal, e.g. rectus turns via gyrans (once strategies) | `composites`; the lab's **Composite** menu (was Strategy: single maneuvers left out; marks + composite, ↻ repeating, ⏱ long hold) |
| role | morph-combo-specific, with **moods** | e.g. *agile mover*: 😡 attack, 😱 avoid, 😰 escape; for 5N limus, GD rectus/gyrans (attack: turns, or rectus turns via gyrans for big ones; pulses when facing). Not only maneuvers: firing too | `behaviors.js` today's behaviors (their states are moods) → roles |
| behavior | species-wide names, species-specific setup | the share of each role / mood: *aggressive* (more attack), *chill* (rest, sleep), ... | settings' Hostile / Player Behavior → these |

The mood shows as an emoji beside the morph's name (the sandbox already shows 😡 😰 😴).

### Lab: testing roles and moods

- Model mode: the **Composite** menu (done: the Strategy menu renamed and trimmed), run step by step.
- A **Role + Mood** menu: the focused morph (changeable) is the player; every other creature is a bot carrying out the selected role in the selected mood.
- Firing (not in the lab): shown as the bot's yellow / red circle turning white.

### Order

1. After the model re-tests (above: verify, the stack pilot, the redesigned starts and phase sweeps): composite maneuvers re-measured (today's numbers come from the old starts).
2. Lab: the Composite menu (done).
3. The user's list of roles and moods → `behaviors.js` roles (moods as today's states), tried in the lab's Role + Mood menu, changed as tests show.
4. Behaviors as role / mood shares, per species; the game's settings and `[hostile_behavior:...]` use them.
5. The story chapter (as before): each species re-introduced with its roles; challenges for the strongest.

Game gaps still to close for roles (from the old step 3): bots that are twins or triples (use the whole's track and model);
still morphs' frames (body axis / board pushes); step-level timing of composite waits (`WAIT`); firing per species.

## Checks

- `node --check` each changed `.js`; `python3 -m py_compile tools/models.py`.
- No headless browser (CLAUDE.md); the user tries the lab and the game.
- Commit as each step lands; don't push.
