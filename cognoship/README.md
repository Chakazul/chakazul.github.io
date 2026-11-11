# Cognoship

A browser game about steering a living spaceship. Your ship is made of *Lenium*, a Lenia
creature (a continuous cellular automaton) that senses and reacts to its surroundings. You
can't push it; you steer it by **cognitive masking**: hiding parts of its senses from itself
(Cool, Hartl, Levin and Petti), so it turns, speeds up or changes shape on its own.

**Play:** https://chakazul.github.io/cognoship/ (or locally, below). It runs in a desktop or
mobile browser. A phone held upright gets a stacked layout; turn it sideways or go fullscreen for
the side-by-side one.

- [README-LAB.md](README-LAB.md): the lab, the experiment page with every control.
- [README-MORE.md](README-MORE.md): the full reference, including how to write the game's script.

## Modes

- **Story Mode:** a journey across the galaxy in chapters, told in the message bar under the
  world view. Each run starts with nothing discovered.
- **Sandbox Mode:** free play with everything you have discovered in every run, in any species
  you have reached (the **Species** button).

Switch with the links at the top of the page. What you discover (morphs, with their pictures,
and bullet types) is saved in your browser as you play and shows up in Sandbox Mode. In Sandbox
Mode, **Hold to clear data** in any book (held for 3 seconds) erases it.

## The screen

- **World view** (left, or top): the whole world. Your ship is in the white frame, circled in
  yellow. Your **energy** shows at the top.
- **Ship view** (right, or below): a close-up of your ship. The yellow arrow is its velocity; the
  name box shows its current **morph** (shape and behaviour), or **???** if you haven't been told
  about it.
- **Message bar:** the story. A message ending in **≫** moves on with a click or Space; one ending
  in **⧗** waits for you to do something (its progress may show after it).

## How to play

- **Steer:** press and hold on the ship in the ship view. The circle you hold hides that spot from
  the ship's senses, and it reacts. Gentle pats on a side or the tail change its shape; holding
  ahead or to one side changes its course.
- **R G B** (top of the ship view, keys `1` `2` `3`): which of the ship's three channels you hide.
  Clicking the selected one again switches to *recycling* (erasing that channel's mass); again
  goes back.
- **Upward** (`E`): turns the ship view so the ship always heads up. **Zoom** (`R` `T`).
- **Morphs:** hold a new morph steady for about a second to discover it. The **Morphs** button
  (`X`) opens the morph book.
- **Dissolved:** if your ship falls apart, **Respawn** (`Z`) brings it back, or **Choose** (`X`)
  another morph you've discovered.
- **Energy:** steer into the glowing dots to collect energy. Other ships eat them too.
- **Fire** (`A`, the round button): shoots the bullet type in use the way the ship is heading
  (a standing ship fires at random). Bullets confuse the senses of what they hit, including your
  own ship if you fly into them. ◀ ▶ (`D` `F`) or the **Bullets** button (`S`) switch types; when
  it reads **Create**, you have enough energy for a new type: open it and pick a **?**.
- **Clone:** when the Morphs button reads **Clone**, pick a discovered morph and place a copy of
  it near your ship (↻ turns it, ↔ flips it). Two morphs together can become a new one.
- **Species** (`C`, Sandbox Mode): switch to another species you have reached.

## Controls

| Key | Does |
|---|---|
| Space or a click on the message | Next message, close a chapter card |
| `Enter` | Pause / play |
| `Backspace` | Retry: start the current scene again |
| `Q` / `W` | Slower / faster |
| `1` `2` `3` | Hide the R, G or B channel (the selected one again: recycle) |
| `E` | Upward |
| `R` / `T` | Zoom out / in |
| `A` | Fire |
| `D` / `F` | Previous / next bullet type |
| `S` | Bullet book |
| `X` | Morph book; Choose when dissolved |
| `C` | Species book (Sandbox Mode) |
| `Z` | Respawn when dissolved |
| Letters, in a book | Pick that slot (hold Shift to see them) |
| `Escape` | Close a book |
| `M` | Music on / off |
| Shift (held) | Show every control's key |

The buttons at the top right: music, light / dark page, fullscreen.

## Run it locally

Serve this folder over http (it fetches its shaders and script, so `file://` won't work):

```
cd Cognoship && python3 -m http.server 8000     # game: http://localhost:8000/   lab: http://localhost:8000/lab.html
```

## Credits

*Lenia* by Bert Chan; *cognitive masking* by Cool, Hartl, Levin and Petti. Music: *Long Way* by
Rexlambo, *A Few Jumps Away* by Arthur Vyncke, *Out of the Skies, Under the Earth* by Chris
Zabriskie, *Billions and Billions* by Stellardrone.
