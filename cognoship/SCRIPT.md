# The game's script (`script.txt`)

The reference for the scripts that drive the game: `script.txt` (story mode), `sandbox.txt` and `challenge.txt` (the same syntax). The rest of the game is described in [README-MORE.md](README-MORE.md).

The game (`index.html`) is driven by `script.txt` (sandbox mode: `sandbox.txt`, the same syntax), read when the page loads and run line by line
(`game.js`, `ACTIONS` and `CONDITIONS`). The world starts empty with no species, so the script
sets everything up.

## Syntax

- **A line** is an optional message followed by any number of tokens: `Message text [token][token]`.
  Tokens are only read **at the end of the line**. A `[...]` in the middle of the text stays as text.
- **A token** is `[name]` or `[name:arg:arg...]`. It is either an action, `[name...]`, which runs
  when the line is reached, in order, or a condition, `[need:name...]`, which is checked every frame.
  The line waits until all of its conditions hold. A token left open at a line's end goes on in the
  next line (joined with no space), e.g. a long list of morphs: `[put_hostile:3:` then
  `gyrans,lithos:far:Geminium]`.
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
- **Retry:** the game's Retry button (left of Pause, `Backspace`) first opens the morph book as "Choose a morph to retry" (while the morph book is on and a morph of the species has been discovered; else it goes on at once), and a pick rebuilds the scene with that morph as the player's (e.g. when the world has grown too crowded): behind a quick fade it empties the bullets channel, puts the dots in afresh (the energy kept), and redoes every action from the last checkpoint, a `[chapter]` or `[new_scene]` line, up to the current line, then the current line itself (its message and conditions start afresh). The first `[put_player]` redone puts the player's morph (the last one put or chosen, the one with the lighter background in the morph book) in place of its own, if it is a morph of the species (so does a challenge's retry). The checkpoints themselves (no card, no extra fade), `[jump_here]`, `[next_msg]` and `[music]` aren't redone.

## Conventions

- **Signatures:** `<arg>` is a value to fill in; `[...]` around part of a signature marks it optional; `a|b` is a choice of literal words. A token left with an optional part out takes its default.
- **Lists:** where an argument takes several values, separate them with `,` (morph ids may also be separated by `:`: `[allow_morph:limus,gyrans]` = `[allow_morph:limus:gyrans]`).
- **Morph ids:** a morph name of the current species (`limus`, `gyrans`, …) or a full pattern id (`5N7KKM-limus`). `<morph>=<k>` (where noted) puts a group of `k` of that morph.
- **Game seconds:** times are in game seconds (1 game s = 30 simulation steps), so none pass while paused (see *Game time* above).
- **Persistence:** unless an entry says otherwise, a setting is kept until the script changes it. *Per attempt* means challenge mode resets it at the start of every challenge and retry (before `[on_challenge_default]` runs).
- **Mode:** entries marked *(challenge)* only do something useful in challenge mode.

## Index

| Tag | Summary |
|---|---|
| **Flow and presentation** | |
| [`[chapter]`](#chapter) | Title card; the line waits until it is closed. |
| [`[new_scene]`](#new_scene) | Fade the views to black, set the scene unseen. |
| [`[blank]`](#blank) | An empty message. |
| [`[next_msg]`](#next_msg) | The continue tip and the ≫ / ⧗ marks. |
| [`[jump_here]`](#jump_here) | (Authoring) start the script here. |
| [`[music]`](#music) | The background playlist. |
| [`[speed]`](#speed) | The game speed. |
| [`[stars]`](#stars) | The star field's drift and brightness. |
| [`[upward]`](#upward) | The ship view's heading-up mode. |
| [`[view_damping]`](#view_damping) | Smooth the scrolling world's view. |
| **Species and morphs** | |
| [`[species]`](#species) | Load a species. |
| [`[add_allow_species]`](#add_allow_species) | Add a species to the visible list. |
| [`[choose_species]`](#choose_species) | Choose a species, as the species book does. |
| [`[on_start]` … `[end]`](#on_start--end) | Actions run whenever a species is chosen. |
| [`[allow_morph]`](#allow_morph) | Which morphs are named. |
| [`[hide_morph]`](#hide_morph) | Leave morphs out of the books. |
| [`[discover_morph]`](#discover_morph) | Mark morphs as discovered. |
| [`[merge]`](#merge) | The detection's merge saddle. |
| **World and creatures** | |
| [`[clear_world]`](#clear_world) | Empty the world; background, empty space, scrolling world. |
| [`[put_player]`](#put_player) | Put the player's ship. |
| [`[put_morph]`](#put_morph) | Put other creatures. |
| [`[min_morph]`](#min_morph) | Keep at least `n` other creatures. |
| [`[target]`](#target) | A target area. |
| **Player tools and books** | |
| [`[morph_book]`](#morph_book) | Show or hide the morph book. |
| [`[morph_book_action]`](#morph_book_action) | What a pick in the morph book does. |
| [`[retry_action]`](#retry_action) | What Respawn does. |
| [`[clone]`](#clone) | Allow clone mode. |
| [`[erase]`](#erase) | Allow the erase tool. |
| [`[draw]`](#draw) | Allow the draw tool. |
| [`[settings]`](#settings) | The Settings button and book. |
| **Energy, dots and beacons** | |
| [`[energy]`](#energy) | Set or hide the energy. |
| [`[dot]`](#dot) | Keep `n` energy dots. |
| [`[beacon]`](#beacon) | Put a beacon (a pattern to eat). |
| [`[keep_beacon]`](#keep_beacon) | A beacon that comes back after the ship eats it. |
| **Gun, pulses and world drawing** | |
| [`[gun]`](#gun) | The player's gun and its settings. |
| [`[add_bullet]`](#add_bullet) | Add bullet types (pulses) to the gun. |
| [`[new_bullet]`](#new_bullet) | Let the player create bullet types for energy. |
| [`[bullet_limit]`](#bullet_limit) | The bullets channel's mass limit. |
| [`[bullet_cost]`](#bullet_cost) | Energy per shot. |
| [`[world_draw]`](#world_draw) | Antimatter and Nebula: painting on the world. |
| **Challenge mode** | |
| [`[challenge]`](#challenge) | Start a challenge. |
| [`[allow_species]`](#allow_species) | The species a challenge allows. |
| [`[on_challenge_default]` … `[end]`](#on_challenge_default--end) | Actions run at the start of every challenge and retry. |
| [`[time_limit]`](#time_limit) | A time limit. |
| [`[dissolve]`](#dissolve) | Whether losing the ship fails the challenge. |
| [`[mask_cost]`, `[erase_cost]`, `[draw_cost]`](#mask_cost-erase_cost-draw_cost) | Energy per second of a tool press. |
| [`[world_draw_cost]`](#world_draw_cost) | Energy per second of painting the world. |
| [`[clone_cost]`](#clone_cost) | Energy per clone. |
| [`[success_msg]`](#success_msg) | The completion message. |
| [`[fail_msg]`](#fail_msg) | The failure messages. |
| [`[put_hostile]`](#put_hostile) | Put bots (hostile creatures). |
| [`[min_hostile]`](#min_hostile) | Keep at least `n` bots. |
| [`[hostile_gun]`](#hostile_gun) | The bots' bullets' settings. |
| [`[hostile_bullet_wait]`](#hostile_bullet_wait) | How often bots fire. |
| [`[hostile_bullet_limit]`](#hostile_bullet_limit) | The bots' bullets channel's mass limit. |
| [`[bullet_mutual]`](#bullet_mutual) | The player's and the bots' bullets wearing each other down. |

The conditions (`[need:…]`) are listed in [Conditions](#conditions) below.

## Flow and presentation

### `[chapter]`

```
[chapter:<no>:<title>]
[chapter:final:<title>]
```

A title card over the game window: "Chapter `<no>`" above `<title>` (`final`: "FINAL CHAPTER"). The title may contain colons.

- Clears the message bar. The world view and the ship view's picture fade out, then the card fades in.
- Space or a click anywhere closes it (that press does nothing else); the card fades out as the views fade back in.
- The game plays on meanwhile; the line waits until the card is closed.
- A checkpoint for Retry.

### `[new_scene]`

```
[new_scene]
```

Fade the world view and the ship view's picture (everything drawn in them, "Dissolved", Respawn and Choose included; not the other controls) to black over 0.5 s, and wait for it.

- The actions after it set the scene unseen; the next message fades the views back in as it shows.
- The game plays on meanwhile.
- A checkpoint for Retry.

### `[blank]`

```
[blank]
```

Makes the line an empty message: the bar shows only its mark (≫, or ⧗ with a condition) and waits for Space or a click, as any message.

### `[next_msg]`

```
[next_msg[:<tip>[:<next>[:<act>]]]]
```

The marks and tip of the message bar.

| Argument | Meaning | Default |
|---|---|---|
| `<tip>` | Grey tip shown once, after the text, on the next message that moves on by Space or a click (not on one skipped over by `[jump_here]` or dev mode's `I`). Text styles work; it can't contain `]`, nor `:` outside `{...}`. | "Press space or click message to continue" |
| `<next>` | Mark ending a message that moves on by Space or a click. | ≫ |
| `<act>` | Mark ending a message that waits on a condition. | ⧗ |

A part left out or empty goes back to its default. Not redone by Retry.

### `[jump_here]`

```
[jump_here]
```

For writing the script: start at the first line that has it. Every action on the lines before it runs first, in order, so the set-up is as if the game had been played to there; their messages, conditions and chapter cards are skipped. Not redone by Retry.

### `[music]`

```
[music:<list>[:<play>]]
```

Replace the background playlist.

| Argument | Meaning |
|---|---|
| `<list>` | Track numbers, `,` between, played in turn and looped, 10 s of silence between plays. Tracks (`page.js` TRACKS): 1 *Long Way*, 2 *A Few Jumps Away*, 3 *Out of the Skies, Under the Earth*, 4 *Billions and Billions*. |
| `<play>` | Switch to that track now, the one playing fading out first (1.5 s); not restarted if it is already playing. With the music off, it is the one that plays when the player turns it on. |

- No `<play>`: the track playing carries on, and the list follows it (from the track after it if it is in the list, else from the list's first).
- Default: track 1 alone. Not redone by Retry.

Example: `[music:2,1:2]`

### `[speed]`

```
[speed:<x>]
```

Set the game speed, as the Speed slider. `<x>`: `1/8`, `1/4`, `1/2`, `1`, `2`, `4` or `8` (others snap to the nearest).

### `[stars]`

```
[stars:<angle>:<speed>:<bright>]
[stars:follow[:<k>[:<bright>]]]
```

The star field's drift.

| Form | Meaning |
|---|---|
| `<angle>` | Drift direction in degrees, clockwise from left: 0 to the left (default), 90 up, 180 right. |
| `<speed>` | Times the default speed (1; 0 stops it). |
| `<bright>` | Times the default brightness (1; 0 hides it). |
| `follow:<k>` | Drift opposite the ship's motion, as if it flew past the stars: the middle layer at the ship's speed × `k` (default 1; the near layer faster, the far one slower), turning as the ship turns (smoothed over a few frames); still while the ship stands still or the game is paused. |

- A part left out or empty keeps its default: `[stars:90]`, `[stars::2]`.
- Changing the angle or speed moves the stars on from where they are.
- The round dots drift the same way at the same factor (beacons in challenge mode stay put).
- A plain `[stars:…]` ends `follow`. In a scrolling world the view already follows the ship; `follow` adds to that.

### `[upward]`

```
[upward:on|off]
```

Turn the ship view's Upward (heading up) on or off, as the Upward button and `E`.

### `[view_damping]`

```
[view_damping[:<s>]]
```

In a scrolling world, smooth how the world view follows the ship, over `s` more game seconds. Default 0 (also with no value or `off`): the view follows closely, moving with the ship's own velocity.

- With `s` > 0 the view moves with the ship's velocity averaged over `s` and closes the gap to it `s` more slowly, so a ship that zigzags or wobbles fast (HAESRE *zigzag*) no longer shakes the whole world: the view follows its mean course while the ship swings about the middle.
- Larger values let the ship stray further from the middle when it speeds up or turns; a gap over 3 R is still jumped at once.
- Kept until changed. For one species only, put it in that species' `[on_start:<code>]` block and `[view_damping:0]` in the general `[on_start]` block (which runs first, for every species).

Example: `[view_damping:1]`

## Species and morphs

### `[species]`

```
[species:<name>:<code>]
[species:<code>]
```

Load species `<code>`: its rule and merge saddle. The world is kept; no morph is visible after it.

- The morph book shows `<name>` in italic at its lower right ("Species: *Aquarium originale*"); with no name, the code.
- Each species has its own discovered morphs (5N7KKM *tardus* is not GDNQYX *tardus*), kept for when it comes back.
- Resets `[merge]`.

Example: `[species:Aquarium originale:5N7KKM]`

### `[add_allow_species]`

```
[add_allow_species:<name>:<code>[:<ids>]]
```

Add species `<code>` (named as `[species]`) to the visible species, in order, with its visible morphs `<ids>` (as `[allow_morph]`; none given: `all`). Adding it again updates its name and morphs. Used by `[choose_species]`, the species book and `[allow_species]`.

Example: `[add_allow_species:Aquarium originale:5N7KKM:limus,gyrans,tardus]`

### `[choose_species]`

```
[choose_species:<code>]
```

As if the player chose species `<code>` in the species book: `[species:<name>:<code>]` and `[allow_morph:<ids>]` from its `[add_allow_species]` line, then the general `[on_start]` block and its own `[on_start:<code>]` block. Sandbox mode starts with one.

### `[on_start]` … `[end]`

```
[on_start[:<code>]]
...
[end]
```

A block of lines taken out of the script's run. Its actions run each time a species is chosen (`[choose_species]`, the species book), after that species' `[species]` and `[allow_morph]`, and, in sandbox mode, at each Retry (which chooses the current species again, with the morph picked in the retry book for `[put_player]`, instead of replaying the script).

- No `<code>`: for every species. `[on_start:<code>]`: for that species only, run **after** the general block, so it only needs what differs (and may change what the general block set).
- What a species block sets stays set when another species is chosen; reset it in the general block (e.g. `[view_damping:0]` there, `[view_damping:2]` in `[on_start:HAESRE]`).
- Also run in challenge mode, when the player picks a challenge's species: before `[on_challenge_default]` and the challenge's own lines (which may change what they set); a retry doesn't run them again.
- Actions only, run at once and in order; messages and conditions in it are ignored with a warning.
- Old name `[on_choose_species]`, still read.

### `[allow_morph]`

```
[allow_morph:<ids>|all]
```

The visible morphs are exactly these, in this order. A morph that isn't visible is still detected but shows as "???". `all`: every morph of the species; no ids: none.

### `[hide_morph]`

```
[hide_morph:<morphs>|none[:<code>]]
```

These morphs (names) stay detected and named as usual, but are left out of the morph book (so they can't be picked to respawn, retry or clone, e.g. too powerful) and the species book's thumbnails, and are never discovered.

- `<code>`: for that species only; none: every species. `none`: un-hide all (for `<code>`, or everywhere).
- In challenge mode each challenge (and retry) starts from the hidden morphs `challenge.txt`'s set-up lines left; a challenge's species book hides the thumbnails its own `[hide_morph]` tags hide.

### `[discover_morph]`

```
[discover_morph:<ids>]
```

Add these morphs to the morph book's discovered list (ones already found stay). Their pictures are black until the ship is one of them for 1 game s.

A morph found in play gets its picture when found. A picture is retaken only once after each fresh start, the first time the ship holds a known morph for 1 s: a book pick that clears the world (`[morph_book_action:clear_put]`), a species chosen (the sandbox's start too), a challenge started or retried, or a sandbox Retry. Not after a Respawn or boarding a clone, when the world may be chaotic. Pictures leave out the Antimatter / Nebula paint.

### `[merge]`

```
[merge[:<x>]]
```

Set the detection's merge saddle to `<x>` (the lab's Merge saddle slider). Higher means fewer merges and smaller objects; above 1 nothing merges. No `<x>`: the species' own value (also restored by `[species]`).

## World and creatures

### `[clear_world]`

```
[clear_world[:<bg>][:empty][:scroll]]
```

Empty the world. Also removes: the dots (none put back until `[dot:<n>]`; the energy is kept), the bullets in flight (the gun and its list are kept), the bots and their bullets, the `[target]`, the Antimatter / Nebula paint, and turns `[min_morph]` and `[min_hostile]` off.

| Option | Meaning | Default |
|---|---|---|
| `<bg>` | The background: `black`; `radiation` (dark orange, darkening and brightening again every 8 s); `gray` (static dark gray); `sky` (a gradient from light cyan at the world view's bottom to dark cyan at its top; the ship view shows the same sky round the ship). | `black` |
| `empty` | Empty space, no player: while there is no ship, the ship view shows the world's centre instead of "Dissolved" (no Respawn or Choose), until the next `[clear_world]`. | off |
| `scroll` | A scrolling world (see [The scrolling world](#the-scrolling-world)); any other `[clear_world]` ends it. | off |

The options combine in any order. Example: `[clear_world:sky:empty]`

### `[put_player]`

```
[put_player:<id>[:<how>]]
```

Put morph `<id>` in as the player's ship, which is then followed. `<id>` becomes the player's morph, which Respawn puts back.

| Argument | Meaning | Default |
|---|---|---|
| `<id>` | A morph id; `first`: the species' first morph. A list (`[put_player:paragyrans,gyrans,first]`): the first that is a visible morph of the species (none visible: the first that is a morph of it, with a warning). `<morph>=<k>`: a group of `k`, as `[min_morph]`'s, the ship being the one nearest the group's centre. | |
| `<how>` | Placement, see [`<how>`](#how-placement). | `center` |

On Retry, the first `[put_player]` redone puts the player's last morph instead, if it is a morph of the species.

### `[put_morph]`

```
[put_morph:<id>[:<how>]]
[put_morph:<n>:<ids>[:<how>]]
```

Put creatures in. They are never followed: the ship stays followed, or stays lost (only `[put_player]` and the morph book's pick make a ship).

| Argument | Meaning | Default |
|---|---|---|
| `<id>` | A morph id; `<morph>=<k>`: a group of `k`, as `[min_morph]`'s. | |
| `<n>`, `<ids>` | `n` creatures, each picked at random from `<ids>` (each may be `<morph>=<k>`; `all`: the visible list). With `center`, the first goes to the centre, the rest at random. | |
| `<how>` | Placement. | `center` |

Examples: `[put_morph:lithos=3:center]`, `[put_morph:4:lithos,tardus=2:far]`

### `[min_morph]`

```
[min_morph:<n>[:<ids>[:<how>]]]
[min_morph:0|off]
```

Keep at least `n` creatures besides the ship on the board: whenever there are fewer, put in one more, picked at random from `<ids>`.

| Argument | Meaning | Default |
|---|---|---|
| `<n>` | The least number of other creatures. `0` or `off`: stop. | |
| `<ids>` | Morph ids; `<morph>=<k>` puts a group of `k` (up to 12) at once: each turned at random, packed around the spot one by one at the free place nearest the group's centre, at least 3 cells between their bodies (cells above 0.01). | `all` (the visible list) |
| `<how>` | Placement. In a scrolling world, only the first `n`; later ones come in at its edges. | `random` |

- At most one every 0.5 game s (none while paused). The creatures are never followed.
- Keeps running across later lines and Respawn until stopped or a `[clear_world]`.

Example: `[min_morph:3:lithos=3:far]`

### `[target]`

```
[target[:<size>][:<how>[:<d>]]]
```

A target area, outlined in dotted orange in both views. The ship is in it while its centre of mass is (`[need:stay_target]`).

| Argument | Meaning | Default |
|---|---|---|
| `<size>` | Side of the square, in R. | the larger side of HAESRE *antigyrans*'s pattern at the game's R (65 cells at R 18) |
| `<how>` | `center`, `random` (a spot clear of creatures), `far` (the clear spot farthest from the ship, of 30), `extended[:<d>]` (beyond the view in a scrolling world). Always wholly inside the board. | `random` |

One at a time (a new one replaces it); `[clear_world]` removes it. Example: `[target:3:extended:20]`

## Player tools and books

### `[morph_book]`

```
[morph_book:on|off]
```

Show or hide the Morphs button (`X`) and the name box of the ship's morph at the bottom of the ship view. Default `on`.

- The book closes only by Escape; letter keys pick slots.
- `off` also closes the book; while off, no morph is discovered (they are still detected).

### `[morph_book_action]`

```
[morph_book_action:clear_put|put|none]
```

What a click on a discovered morph in the book does while the ship is lost (otherwise the book is read-only, or picks a clone): `clear_put` = `[clear_world][put_player:<id>]` (default), `put` = `[put_player:<id>:random]`, `none` = nothing.

Keys pick slots, the book's rows on the keyboard's: `A` `S` `D` `F` the first row, `Z` `X` `C` `V` the second (Shift or Alt shows each pickable slot's letter).

### `[retry_action]`

```
[retry_action[:clear_put|put|none|book]]
```

What Respawn (and `Z`) does, with the player's morph: as `[morph_book_action]`, or `book` (open the morph book). No type: the default, `clear_put`.

### `[clone]`

```
[clone:on|off]
```

Allow clone mode. Default `on`.

- While allowed and the ship is alive, the Morphs button has a pink border and reads **Clone**; the book, "Choose a morph to clone", lets discovered morphs be picked (pink borders). `off`: the book is read-only while the ship is alive.
- The picked morph follows the pointer over the ship view; ↻ and ↔ below Morphs turn it 45° clockwise and flip it; a click pastes it. Keys: `Q` turn, `W` flip, `Enter` paste, `Escape` cancel (they replace Q / W's speed and Enter's Play while a clone is held).
- Touch: a tap (or drag) places the copy, which can still be turned, flipped or moved; the ✓ button pastes it. Buttons below Morphs: ↻ ↔ / ✓ 🗙.
- The game pauses from the pick until the paste or cancel.
- The player **boards the clone**: once pasted, it is the followed ship (the old one goes on as another creature; a clone pasted onto the ship makes one object with it), and its morph is what Respawn puts back. So a morph can be changed without dissolving, and a twin engineered by becoming one half and pasting the other next to it.
- See also `[clone_cost]`.

### `[erase]`

```
[erase:on|off]
```

Allow the erase tool ("recycle": the selected R G B patch clicked again; drawn as black disks on the colour, "−R"). `off` also goes back to masking. Default `on`.

### `[draw]`

```
[draw:on|off|if:<code>]
```

Allow the draw tool (add mass: the selected R G B patch clicked again after erase, or after mask if erase is off; drawn as coloured disks on black, "+R"): a press on the ship view sets that channel to 1 inside the tool circle, as the lab's Draw. Default `off`.

- `off` also goes back to masking. `if:<code>`: only while the species is `<code>` (e.g. `[draw:if:G6G6CR]`).
- With erase and draw both on, the patch cycles mask → draw → erase → mask.

### `[settings]`

```
[settings:on|off]
```

A **Settings** button (`V`) at the lower centre of the world view, opening the settings book. Default `off`. Each setting is remembered in the browser.

| Setting | Meaning | Until first set |
|---|---|---|
| **Autospawn** | Once the ship has dissolved, Respawn by itself a second later, with no "Dissolved" panel; the ship view stays where the ship was last seen. | off |
| **Scroll** | Every `[clear_world]` makes a scrolling world, with `[stars:follow:0.1]`. Changed in the book: the current species is chosen again as the book closes. | off |
| **Objects** | Off, 1, 2, 5, 10, 20 or 40: `[min_morph:<n>]`, run after every `[on_start]` block (overriding its count; its morphs and placing stay) and at once when moved. Lowering it removes nobody. | 2 |

## Energy, dots and beacons

### `[energy]`

```
[energy:<n>|on|off]
```

Set the energy (the dots score) to `<n>` and show it at the top of the world view. `on`: show it, unchanged. `off`: hide it (dots and beacons still add to it); `on`, `[energy:<n>]` or a new challenge attempt shows it again.

### `[dot]`

```
[dot:<n>[:<how>]]
[dot:off]
```

The dots channel (`dots.js`): keep `n` energy dots on the board, which drift slowly (as `[stars]` sets). Any creature eats them; the ship eating one scores 10.

| Argument | Meaning | Default |
|---|---|---|
| `<n>` | The number of dots kept. `off`: no dots, no scoring. | |
| `<how>` | Where each is put and replaced. A scrolling world puts only the first `n` so and replaces them at its edges. | `random` |

- The score shows at the top, centred, and is kept while the dots are off (a story run starts from 0, a sandbox run from the central inventory's energy).
- If the channel grows past `n + 2` dots' worth of mass, it is emptied and the dots put in afresh (no score). For beacons: 3 beacons' worth, stamped afresh where they were.

### `[beacon]`

```
[beacon:<patterns>[:<n>][:<how>[:<d>]]]
```

Put one beacon, a pattern to be eaten, in the dots channel instead of the round dots. It is not replaced; several lines put several.

| Argument | Meaning | Default |
|---|---|---|
| `<patterns>` | Pattern names (`,` between; first word, full name or key, case ignored), one picked at random; `all`: any. If beacons are already on the board and their pattern is listed, that one again (the channel runs one rule at a time: another would empty it). | |
| `<n>` | Energy the ship gains eating it. | 10 |
| `<how>` | `center`, `random` (clear of creatures), `far` (the clear spot farthest from the ship, of 30), or `extended[:<d>]` (beyond the view in a scrolling world; `d` R from the ship). | `random` |

| Pattern | From `animals.json` | Notes |
|---|---|---|
| `Heptadentium` | *H. scintillans* (2D7s) | |
| `Helicium` | *H. serratus turbulentus* (4H3t) | |
| `Tetrahelicium` | *T. serratus turbulentus* (2H4t) | runs at R 12 |
| `Circoechinium` | *C. ventilans* (3ECv) | wanders slowly |
| `Gyrogeminium` | *G. serratus* (3GG2s) | wanders slowly |
| `Decadentium` | *D. rotans bioculi* (2D10rb) | |
| `Decadentium-3` | *D. rotans trioculi* (2D10rt) | |
| `Nivium` | *N. incarceratus* (3R6i) | |
| `Trilapillium` | *T. inversus* (2L3i) | |

- Each runs under its own rule, at its own R (at the game's R 18) unless noted. The wanderers move about 0.25 cells per game second; a beacon's site follows its pattern's centre of mass, so it is still found, eaten and counted where it has gone.
- The round dots go (a later `[dot:<n>]` brings them back, emptying the channel).
- It drifts with the dots, except in challenge mode, where beacons stay put.
- Any creature eats it. One eaten by another creature, or worn away, is put back nearby (2 to 4 R away, clear of creatures, the same energy); see `[need:clear_beacon]`.
- `[dot:<pattern>…]`, the old form, still works with a warning.

Example: `[beacon:Heptadentium:30:far]`, `[beacon:Helicium,Tetrahelicium]`

### `[keep_beacon]`

```
[keep_beacon:<patterns>[:<n>][:<how>[:<d>]]]
```

As `[beacon]`, but once the ship eats it, another is put in at a new spot by the same `<how>` and `<d>`, with the same `<n>`, again and again (no clear spot: tried again the next frame).

- The new one is a pattern picked afresh from `<patterns>`, not the one just eaten (if the list has another), while no other beacon is on the board; with others about, the same pattern.
- One eaten by another creature, or worn away, is put back nearby as `[beacon]`'s.
- `[need:clear_beacon]` never completes while one is kept; `[clear_world]` ends them.

Example: `[keep_beacon:all:20:extended:20]`

## Gun, pulses and world drawing

The UI calls the bullets *pulses* (the **Pulses** button, the pulse book); the tags keep `bullet`.

### `[gun]`

```
[gun[:<options>]]
```

The player's gun and its bullets channel (`bullets.js`). Options in any order, `,` or `:` between:

| Option | Meaning | Default |
|---|---|---|
| `on` / `off` | The gun on, or no gun. | `on` |
| `mask` / `erase` | A bullet hides the cells it overlaps from the creature's own kernels, softly, as the mask tool; or takes their mass out. | `mask` |
| `color_now` / `rgb` / `r` / `g` / `b` | The channels it acts on: the R G B patch selected (live: bullets in flight switch too), all three (white bullets), or one. | `color_now` |
| `wear_color` / `wear_all` | Which creature channels wear a bullet down: those it acts on, or all three. | `wear_color` |

- **Firing:** `A` or the round Fire button (lower left of the world view; dimmed while it can't fire). A bullet is put just in front of the ship (its back edge 8 cells from the ship's edge; a twin or trio ship fires from the whole) and flies the way the ship moves; a ship standing still fires at a random angle. *Asterium* and *Tetralapillium*, which barely move, are dropped behind like bombs. *Helicium* is still fired ahead, to run into: it changes the ship's shape.
- **Effect:** where a bullet overlaps a creature, it masks (or erases) gradually while the creature's mass wears it down, so a bullet sinks some way into a body. A ship faster than its bullets runs into them.
- **Limits:** at most one shot every 0.7 game s. Past `[bullet_limit]` bullets' worth of mass, the channel is emptied; under 1% of a bullet's mass, it is wiped and goes idle (no cost) until the next shot.
- Changing the pattern empties the channel.
- Tuning in the URL: `?bhit=0.2` (mass taken per step), `?bwear=0.03` (bullet worn per bullet step), `?bsub=3` (bullet steps per sim step), `?bcool=0.7` (seconds between shots).

### `[add_bullet]`

```
[add_bullet:<names>[:<equip>]]
```

Add bullet types to the gun's list, which starts empty (no list, no shots).

| Argument | Meaning |
|---|---|
| `<names>` | In order, `,` between: a pattern name (see [Bullet patterns](#bullet-patterns); case and spaces ignored); `random` (one not listed yet); `all` (every pattern; not saved to the inventory); `clear` (empty the list of patterns); `antimatter` / `nebula` (that world-draw type, see `[world_draw]`; free, never saved, always revealed when named). |
| `<equip>` | Switch to this type (added if not listed). None: the last added; none added: the one in use if still listed, else the list's first. |

- The types added are shown: the pulse book opens and reveals them as new cards one after another; the equipped one shows as in use once all are revealed. Not when redone by Retry, `[jump_here]` or dev `I`, nor while a book or chapter card shows.
- Each costs `[new_bullet]`'s price while that is on.
- The type in use is named under Fire. With two or more, ◀ ▶ above Fire (`D` / `F`) switch type (emptying the bullets channel).
- **The pulse book:** the **Pulses** button (top left of the ship view) or `S`; only Escape closes it. Every pattern by name with a small picture, 4 to a row; listed ones clickable (the one in use with a yellow border), the others "?". Keys: `1`–`4`, `Q`–`R`, `A`–`F`, `Z` `X` by row.

Examples: `[add_bullet:Orbium:Orbium]`, `[add_bullet:antimatter:antimatter]`

### `[new_bullet]`

```
[new_bullet:<n>|off]
```

Let the player create bullet types for `n` energy each (`off`: none may be created).

- While the energy is at least `n`, the **Pulses** button reads **Create** in green, with a glint as it turns so; the pulse book ("Create or use a bullet type") shows its "?" slots with a green border, and a click creates that pattern, switches to it and reveals it.
- Types the script gives (`[add_bullet]`) cost `n` too while it is on (the energy never goes below 0).

### `[bullet_limit]`

```
[bullet_limit[:<x>]]
```

Empty the bullets channel when it holds more than `x` bullets' worth of mass (bullets that run into each other merge and grow). Default 8, also with no `<x>`.

### `[bullet_cost]`

```
[bullet_cost:<n>|off]
```

Each shot costs `n` energy (the energy then shows at the top); with less, the gun doesn't fire. The bots' bullets are free. `0` / `off`: free (default). Per attempt.

### `[world_draw]`

```
[world_draw:<types>]
```

Put the world-draw types **Antimatter** and / or **Nebula** into the pulse book, quietly, after the patterns. Neither is there by default.

| `<types>` (`,` or `:` between) | Meaning |
|---|---|
| `antimatter`, `nebula` | Add that type. Each tag adds to what the book has. |
| `on` | Both. |
| `off` | Take both out, and empty the fields (one in use: back to the list's last pattern). |

`[add_bullet:antimatter]` / `[add_bullet:nebula]` adds one with the card reveal instead (the story).

- **Painting:** in use, they fire nothing: the player paints three static fields, R, G and B, onto the world with a press in the world view or the ship view (`field.js`). A field never grows, moves or wears down. **Antimatter** erases its channel of every creature (caps it at 1 − field each step, at once; drawn dark). **Nebula** masks it (hides it from the creature's own kernels; drawn as a translucent glow).
- **Tools:** the R G B patches pick the colour; the selected one again cycles **Paintbrush** (+, raises the field to full), **Watercolor** (×, a soft disc fading to its edge; going over it adds nothing), **Eraser** (−). The tool's name shows at the world view's lower right.
- **Controls:** Zoom becomes **Width** (2–32 cells; `R` `T`; watercolor and eraser twice as wide), Fire becomes **Clear** (the selected colour's field; `A`), Upward becomes **Undo** / **Redo** (`E`, one level).
- **With other channels:** the player's bullets are replaced (those in flight emptied). The bots keep firing: with Antimatter their bullets act as usual; with Nebula the two masks merge (per channel the stronger; 2 more FFTs per step), but with erasing bots (`[hostile_gun:erase]`, testing only) Nebula alone acts. The dots and beacons stay.
- The fields are kept (hidden, no effect) while another type is in use, and emptied by `[clear_world]` and `[world_draw:off]`.
- See also `[world_draw_cost]`, `[need:draw]`.

### Bullet patterns

The patterns `[add_bullet:<name>]` adds (`bullet_patterns.js`), all single-channel classic Lenia from `~/github/Lenia/Python/animals.json` (the entry with that code; #2, #3 = its second, third entry), T = 10. Native R is the pattern's own radius; every bullet runs at R 13 (with the default creature R 18; scaled with `?R=`). Speed: cells per sim step at R 13 and `?bsub=3`, measured in numpy at launch angles 0, 90, 200 and 315°; the gliders held their mass and heading except where noted. The old ids (`scutium2`–`scutium5`, `orbium1`, `paraptera1`, `Scutium Mk II`–`IV`) still work.

| Name | Species | Code | Native R | b | Speed |
|---|---|---|---|---|---|
| `Orbium` | Orbium unicaudatus | O2u | 13 | 1 | 3.1; dies at some launch angles |
| `Synorbium` | Synorbium solidus | O4s | 13 | 1 | 3.2 |
| `Scutium` | Scutium valvatus | S1v | 13 | 1 | 2.1 |
| `Scutium fluens` | Scutium serratus fluens | 2S1f | 15 | 1/2,1 | 2.2 |
| `Scutium vagus` | Scutium serratus vagus | 2S1v | 10 | 1/2,1 | 1.9 |
| `Scutium serratus` | Scutium serratus fluens | 2S1f #2 | 18 | 1/6,1/2,1 | 3.0 |
| `Triscutium` | Triscutium valvatus laterale | S3vl | 15 | 1 | 2.2 |
| `Helicium` | Helicium cavus pedes | H3cp | 13 | 1 | spins in place |
| `Geminium` | Hydrogeminium natans | 3GH2n | 18 | 1/2,1,2/3 | 1.5 |
| `Kronium` | Kronium vagus | 3K4v #2 | 27 | 2/3,1/3,1 | 1.1; wanders |
| `Asterium` | Asterium vagus | 3R5v #3 | 27 | 1,1/2,1 | 0.4; wanders, heading unpredictable; dropped behind |
| `Tetralapillium` | Tetralapillium metamorpha | 2L4m | 18 | 1,1/2 | changes shape, creeps; dropped behind |
| `Paraptera` | Paraptera arcus pedes | P4ap | 13 | 1 | 2.4 |
| `Heptapteryx` | Heptapteryx arcus labens | P7al | 13 | 1 | 2.3 |

## Challenge mode

### `[challenge]`

```
[challenge:<name>:<n>[:<subtitle>]]
```

Start challenge `<name>`, stage `<n>`: a title card as `[chapter]`'s, "`<NAME>` CHALLENGE" small above, `<subtitle>` big below (none: "Stage `<n>`"). The line waits until it is closed.

- In the challenge book and messages: *`<name>`: `<subtitle>`* (no subtitle: *`<name>`: `<n>`*). The record is kept as `<name>:<n>`.
- A retry starts from the line after it.

Example: `[challenge:Courier:1:The Delivery]`

### `[allow_species]`

```
[allow_species:all|<codes>]
```

On a challenge's lines: the species the player may choose for it. `all` (default): every species added by `[add_allow_species]`. Does nothing when run.

Example: `[allow_species:5N7KKM,HAESRE]`

### `[on_challenge_default]` … `[end]`

```
[on_challenge_default]
...
[end]
```

A block of lines taken out of the script's run. Its actions run, in order and at once, at the start of every challenge and retry, after the attempt's reset (no costs, energy 0, no time limit, default messages) and before the challenge's own lines, which may change them. Actions only (messages and conditions ignored, with a warning).

### `[time_limit]`

```
[time_limit:<s>|off]
```

The challenge must be completed within `s` game seconds from here, or it fails (`time`). The time left shows at the top left (⧖ 23.4), red in the last 5 s. `off` (or no value): no limit. Per attempt.

### `[dissolve]`

```
[dissolve:retry|respawn]
```

What the ship dissolving does. `retry`: the challenge fails (`lost`); Respawn or Choose start the attempt again. `respawn` (default): as in the other modes ("Dissolved"; Respawn and Choose put a ship back by `[retry_action]`, and the attempt goes on, its time counted from the respawn).

### `[mask_cost]`, `[erase_cost]`, `[draw_cost]`

```
[mask_cost:<n>|off]
[erase_cost:<n>|off]
[draw_cost:<n>|off]
```

A press on the ship view with that tool costs `n` energy per game second held. The energy shows at the top; reaching 0 fails the challenge (`energy`), except while a `[need:survive]` line waits (then the press is free at 0). Give some first with `[energy:<n>]`. `0` / `off`: free. Per attempt.

### `[world_draw_cost]`

```
[world_draw_cost:<n>|off]
```

Painting the world with Antimatter or Nebula (the paintbrush or watercolor held; the eraser, Clear and Undo are free) costs `n` energy per real second (painting works while paused). Reaching 0 fails as `[mask_cost]`. Per attempt.

### `[clone_cost]`

```
[clone_cost:<n>|off]
```

Each clone pasted costs `n` energy (the energy then shows at the top); with less, it can't be pasted (it stays held, the Clone button flashes). `0` / `off`: free (default). Per attempt.

### `[success_msg]`

```
[success_msg[:<text>]]
```

The message when the challenge is completed. Placeholders: `{title}` (the challenge's subtitle), `{time}` (*12.3s*, *1m 12s*), `{best}` (*(best 1m 12s)*, *(new completion!)* or *(new record!)*). A `:` in the text is kept (as `: `). No text: the default. Per attempt.

Default: `_{title}_ completed in {time} {best}` (`survived for` when the challenge waits on `[need:survive]`).

Example: `[success_msg:You held out for {time}! {best}]`

### `[fail_msg]`

```
[fail_msg:<cause>[:<text>]]
```

The message when the challenge fails by `<cause>`. Placeholders: `{title}`, `{time}` (the attempt so far). No text: the default. Per attempt.

| `<cause>` | Fails when | Enabled by | Default reason |
|---|---|---|---|
| `lost` | the ship dissolves | `[dissolve:retry]` (not while a `[need:survive]` line waits) | the ship dissolved |
| `short` | the ship dissolves before `n` s | `[need:survive:<n>]` | the ship dissolved too soon |
| `time` | the time runs out | `[time_limit]` | time's up |
| `energy` | the energy reaches 0 | `[mask_cost]`, `[erase_cost]`, `[draw_cost]` or `[world_draw_cost]` above 0 (not while a `[need:survive]` line waits) | out of energy |

Default message: `_{title}_ failed: <reason>.` A failed challenge stays failed (Space and clicks do nothing) until a retry.

Example: `[fail_msg:short:_{title}_: the swarm got you after {time}.]`

### `[put_hostile]`

```
[put_hostile:[<n>:]<morphs>[:<how>][:<bullet>]]
```

Put `n` bots into the world: hostile creatures that fire at the ship.

| Argument | Meaning | Default |
|---|---|---|
| `<n>` | How many. | 1 |
| `<morphs>` | Morph names of the current species, one picked at random for each (those it lacks are left out). | |
| `<how>` | Placement, as `[put_morph]`; `center`: the first only, the rest random. | `random` |
| `<bullet>` | The bullet pattern the bots fire (one for all the bots: the last `[put_hostile]` or `[min_hostile]`'s). | `Scutium` |

- Each is marked by a red circle with its morph's name above it (a twin or trio named as one shows as one). Never followed.
- **Firing:** every `[hostile_bullet_wait]` seconds a bot fires at the ship, give or take 20° (no ship: anywhere); a named twin or trio fires as one. Their bullets have their own field, drawn white, acting on whatever they hit, bots included, with `[hostile_gun]`'s settings (by default the player's `[gun]`'s).
- **Tracking:** bots are tracked by place: whatever a bot turns into, both halves when it splits, and what it merges with all stay bots; a bot the ship hijacks (merging, then followed) is the player's; one missing from detection is kept 1 game s.
- `[clear_world]` removes the bots and their bullets.

Examples: `[put_hostile:3:gyrans,lithos:far:Scutium serratus]`, `[put_hostile:trias:far:Heptapteryx]`

### `[min_hostile]`

```
[min_hostile:<n>[:<morphs>[:<how>]][:<bullet>]]
[min_hostile:0|off]
```

Keep at least `n` bots on the board, as `[min_morph]` keeps creatures: whenever there are fewer (a bot dissolved, shot apart or hijacked), put in one more as `[put_hostile]` does.

| Argument | Meaning | Default |
|---|---|---|
| `<n>` | The least number of bots. `0` or `off`: stop. | |
| `<morphs>` | Morph names; `all`: the visible list. | `all` |
| `<how>` | Placement; in a scrolling world only the first `n`, later ones at its edges. | `random` |
| `<bullet>` | As `[put_hostile]`'s. | `Scutium` |

At most one every 0.5 game s (none while paused). `[clear_world]` turns it off.

Example: `[min_hostile:5:all:random:Geminium]`

### `[hostile_gun]`

```
[hostile_gun[:<options>]]
```

The bots' bullets' own settings, instead of following the player's `[gun]` (the default). `[gun]`'s options, over the defaults `mask`, `color_now`, `wear_color`:

| Option | Meaning |
|---|---|
| `mask` / `erase` | As `[gun]`. |
| `color_now` / `rgb` / `r` / `g` / `b` | As `[gun]` (`color_now`: the player's R G B pick, live). |
| `wear_color` / `wear_all` | As `[gun]`. |
| `off` / `on` | The bots hold fire / fire again (alone, the other settings are kept). |
| `follow` (or no options) | Back to following the player's `[gun]`. |

- Kept until changed, across challenges too (put it in `[on_challenge_default]` to reset it).
- The bots' bullets are drawn white whatever their channels.
- While both the player's and the bots' bullets fly, the creatures' step takes them as one field: the same settings, added; both masking on different channels, one per-channel mask (each channel the stronger; the FFT path then transforms it every step); otherwise (testing only) both act as the player's, with a console warning.

Example: `[hostile_gun:on,mask,rgb,wear_all]`

### `[hostile_bullet_wait]`

```
[hostile_bullet_wait[:<a>:<b>:<c>]]
```

Bots fire every `a` to `b` game s at random, a new bot's first shot after `c` more. Default `1.5:4:2`, also with no values.

### `[hostile_bullet_limit]`

```
[hostile_bullet_limit[:<x>]]
```

As `[bullet_limit]`, for the bots' bullets field. Default 8.

### `[bullet_mutual]`

```
[bullet_mutual:none|wear[:<k>]|masking]
```

The player's and the bots' bullets together.

| Value | Meaning |
|---|---|
| `none` | Apart (default). |
| `wear[:<k>]` | Where they overlap, each wears the other down at its own wear rate (0.03 per bullet step, `?bwear=`) per unit of the other's mass, times `k` (default 1; `2` twice as fast, `0.5` half), so shots that meet cancel out. |
| `masking` | Not built yet (as `none`, with a warning). |

## Shared values

### `<how>`: placement

| Value | Meaning |
|---|---|
| `center` | The board's centre (a scrolling world: the view's middle), unrotated. |
| `random` | A random empty spot and angle, as the lab's Add. |
| `far` | The empty spot farthest from the ship, at a random angle. |
| `extended[:<d>]` | Beyond the world view, for `[beacon]`, `[keep_beacon]` and `[target]` in a scrolling world (`d` R from the ship); otherwise as `far`. |

A creature is put at least 5 R (90 cells) from the beacons and the target area (so it doesn't eat the one or sit in the other; in a scrolling world those off the board too), except with `center`.

### The scrolling world

`[clear_world:scroll]`: the board (still a 384 × 256 torus) becomes a window on an unbounded world that moves with the ship.

- **The view** (the whole board) follows the ship smoothly, keeping it at its middle: it moves with the ship's smoothed velocity and eases out the rest (over ~20 game steps; 1 s while paused), so the detected centre's wobble doesn't shake it; a jump of over 3 R (a respawn elsewhere) is taken at once. `[view_damping]` smooths it more.
- **The seam:** the board's wrap runs along the view's edges, under a black stripe 2 R (36 cells) wide that is only drawn (the world under it runs and is detected). A band 6 cells deep inside each edge is cleared every frame (creatures, dots, bullets), so nothing wraps round to the ship's front.
- **Objects taken away:** an object (a creature part, not the ship's) is taken away as soon as its body (its circle, twice its radius) touches an edge's band. While `[min_morph]` runs, its replacement goes on the opposite edge, 0.1 R inside the band (moved along the edge if the spot is taken; not taken away again for 1.5 game s). Every later put (one eaten, shot or dissolved) comes in on the edge the ship is heading to (give or take 40°; standing still, any edge). Round dots the same. Other creatures are put at least 25.8 cells inside the edges.
- **Things beyond the view:** a beacon or target put `extended` is kept at a world point beyond the view (at a random angle, 1.25 to 2 times the distance to the view's edge that way; `extended:<d>`: exactly `d` R from the ship), with a big blinking yellow arrow at the view's edge pointing to it (none after a failed attempt). It is put on the board as the ship nears it (about 31 cells inside the edge, partly under the stripe) and taken off again, still kept, once the ship has left it behind. So `[need:clear_dot]` waits for beacons off the board too.

## Conditions

```
<message>[need:<name>[:<args>]]
```

A line waits until all of its conditions hold. When they do, a green **SUCCESS** box fades in in the middle of the game window and out again (3 s; the game goes on), and the script moves on. A line whose conditions already hold when reached shows ≫ and moves on by Space or a click (with SUCCESS too).

Where noted, the ⧗ mark is followed by the progress in grey. Counts in game seconds don't pass while paused, nor while a book or a chapter card is open.

| Condition | Holds when | Progress |
|---|---|---|
| [`first_touch`](#needfirst_touch) | The player pressed on the ship and it changed morph. | |
| [`upward`](#needupward) | Upward is on. | |
| [`discover_morph`](#needdiscover_morph) | Morphs are discovered. | |
| [`respawn`](#needrespawn) | The ship was put back after dissolving. | |
| [`recycle`](#needrecycle) | The player erased mass, then the ship dissolved. | |
| [`first_clone`](#needfirst_clone) | A clone was pasted. | |
| [`create_bullet`](#needcreate_bullet) | A bullet type was created. | |
| [`draw`](#needdraw) | Enough paint put on the world. | draw: 4.2/10 |
| [`energy`](#needenergy) | The energy reached `n`. | |
| [`object`](#needobject) | The number of objects. | morph: 20/60 |
| [`mass`](#needmass) | The creatures' mass in a channel. | mass: 12.3/20 |
| [`survive_time`](#needsurvive_time) | The ship stayed alive `s` s. | time: 12/30 |
| [`survive`](#needsurvive) | The ship dissolved after lasting `n` s. | time: 42/60 |
| [`stay_morph`](#needstay_morph) | The ship stayed one morph `s` s. | time: 2/3 |
| [`stay_target`](#needstay_target) | The ship stayed in the target `s` s. | time: 4/10 |
| [`clear_beacon`](#needclear_beacon) | Every beacon eaten by the ship. | |
| [`clear_dot`](#needclear_dot) | The dots channel is empty. | left: 1 |
| [`clear_hostile`](#needclear_hostile) | No bots left. | bots: 2 |

### `[need:first_touch]`

The player has pressed on the ship, and its detected morph is now not *limus*. The press is required, since a fresh ship can read as unknown at first.

### `[need:upward]`

The ship view's Upward is on.

### `[need:discover_morph]`

```
[need:discover_morph[:<ids>[:any|all]]]
```

One of `<ids>` (`any`, default) or every one (`all`) is in the morph book: discovered earlier (or put in by `[discover_morph]`), or discovered now and the book that opened for it closed. No ids: the player discovers any new morph after the line is reached, and closes the book that opened for it.

### `[need:respawn]`

The ship has dissolved and the player has put it back since the line was reached: by Respawn (`Z`), or Choose (`X`) and a pick in the morph book.

### `[need:recycle]`

The player has erased some mass on the ship view (any channel) since the line was reached, and the ship has then dissolved.

### `[need:first_clone]`

The player has pasted a clone since the line was reached (`[clone:on]`; picked in the morph book and pasted by a click, or ✓ after a touch). Cancelling doesn't count.

### `[need:create_bullet]`

The player has created a bullet type since the line was reached: **Create** (the Pulses button) and a "?" in the pulse book. `[new_bullet]` must allow one, or the line waits until it does.

### `[need:draw]`

```
[need:draw:<x>]
```

The player has painted at least `x` g into the world's R G B fields (Antimatter or Nebula; the three together) since the line was reached: what the paintbrush and watercolor add, in g as `[need:mass]` (10 g is 3,240 full cells at R 18). The eraser, Clear, Undo and the scrolling world's seam don't take it back.

### `[need:energy]`

```
[need:energy:<n>]
```

The energy (the dots score) has reached at least `n`.

### `[need:object]`

```
[need:object:<n>[:exact|min|max]]
```

The world holds exactly `n` objects (`exact`, default), at least `n` (`min`) or at most `n` (`max`): every creature detected, the ship included. The count is smoothed (the median over the last game second), so a creature flickering between one and two objects doesn't make it jump.

### `[need:mass]`

```
[need:mass:<r|g|b|rgb>:<x>[:min|max]]
```

The creatures detected, the ship included, hold at least `x` g (`min`, default) or at most `x` (`max`) of mass in that channel, or all three (`rgb`); g as the stats charts count it (cell values summed / R²). Smoothed as `[need:object]`.

### `[need:survive_time]`

```
[need:survive_time:<s>]
```

The ship has stayed alive for `s` game seconds since the line was reached. Dissolving stops the count (still shown); it starts again from 0 once the ship is back.

### `[need:survive]`

```
[need:survive[:<n>]]
```

The ship dissolves (or is otherwise lost), having lasted at least `n` game seconds (default 0): a success, with SUCCESS, the green ✓ and the time recorded. Meant for a challenge's last line, a survival round (e.g. with `[min_hostile]`).

- **Too soon:** in challenge mode the challenge fails (`short`); elsewhere the line waits on.
- **The time:** the challenge's (from the attempt's start, or the last respawn); elsewhere from when the line was reached. The count goes on past `n` ("time: 75/60").
- **While it waits:** `[dissolve:retry]` doesn't fail the challenge, and neither does running out of energy (the energy stays at 0; firing and cloning wait for more).
- **Best time:** the challenge's best is its **longest** ("survived for 1m 12s").

### `[need:stay_morph]`

```
[need:stay_morph:<morph>:<s>]
```

The ship has been detected as `<morph>` for `s` game seconds in a row; any other name, unknown, or dissolving starts the count again. The name needn't be visible. Docking: a gyrans circling a lithos is detected as one object, the fuzzy *circumgyrans* (5N7KKM, HAESRE).

Example: `[need:stay_morph:circumgyrans:3]`

### `[need:stay_target]`

```
[need:stay_target:<s>]
```

The ship has stayed in the `[target]` area for `s` game seconds in a row; leaving it, or dissolving, starts the count again.

### `[need:clear_beacon]`

Every beacon has been eaten by the player's ship (a twin or trio ship: any of its parts). A beacon eaten by another creature, or worn away, is put back nearby, so the line waits for the ship to eat it. Never holds while a `[keep_beacon]` is kept.

### `[need:clear_dot]`

The dots channel holds nothing: every dot or beacon eaten (beacons off the board in a scrolling world, and put-back ones, included), or the dots off.

### `[need:clear_hostile]`

(Challenge) No bots left: every one destroyed or hijacked. Other creatures (e.g. `[put_morph]`'s) may stay.

## Example

```
[chapter:1:Metamorphosis]
[species:Aquarium originale:5N7KKM][allow_morph:limus][discover_morph:limus]
[morph_book:on][clear_world][put_player:limus]
Welcome to the Cognoship game.
Click on the ship to change its shape.[need:first_touch]
[dot:3][min_morph:2:gyrans,tardus:far]
Gain 50 energy. [need:energy:50]
```
