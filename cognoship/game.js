"use strict";
// ============================================================================================
//  Cognoship, the game (index.html), on top of the shared core (core.js). For now a trimmed
//  lab: steer the default species's ship (or ?form=) by masking. Game logic (goals, levels, story
//  mode, more channels and objects) goes here; experiments stay in the lab (lab.html, lab.js).
// ============================================================================================

(function () {
  const C = Core;
  const { $, flash, fillRange } = Page;

  // The ship view shows only the channel patches: the action is always Mask (Erase and Draw
  // stay in the core, for the lab). The star field is always on; the Music button only plays or
  // pauses the music.
  C.opt.stars = () => true;
  C.opt.viewLabels = () => false;
  C.opt.actionPatches = () => false;
  C.action = 0;
  // the R G B patch again (or its key): erase ("recycle"), striped, -R / -G / -B; allowed unless [erase:off]
  let eraseOn = true;
  C.opt.patchErase = () => eraseOn;
  C.opt.morphItalic = () => true;
  C.setDetectMode('focus');   // Detect: Focused, the focus's morph at the top of the ship view
  C.opt.nameBox = () => book.on;   // the name box only with the morph book on ([morph_book])
  C.opt.nameDot = () => true;
  C.opt.nameBoxLow = () => true;   // the name box at the bottom of the ship view, above Upward / Zoom
  // ?dev=1: dev mode: test keys and buttons (I, O, P: see devSkip), the frame rate and steps per
  // second at the top right of the world view, and a grey circle around every other detected object
  // in the world view too (the core's circles, as the lab's: the ship's yellow, the one under the
  // pointer brighter)
  const DEV = new URLSearchParams(location.search).get('dev') === '1';
  // ?mode=sandbox: sandbox mode, the same game run by sandbox.txt (not script.txt), starting from the
  // central inventory (see the inventory, below); else story mode (script.txt), from an empty one.
  // The title row's two links switch between them (other URL parameters kept).
  const MODE = new URLSearchParams(location.search).get('mode') === 'sandbox' ? 'sandbox' : 'story';
  const SCRIPT_FILE = MODE === 'sandbox' ? 'sandbox.txt' : 'script.txt';
  // ?all=1: all discovered (for testing): every visible morph of every species counts as discovered
  // (allDiscovered) and the gun has every bullet pattern from the start; the central inventory is
  // neither read nor written (in either mode).
  const ALL = new URLSearchParams(location.search).get('all') === '1';
  for (const [id, m] of [['mode-story', 'story'], ['mode-sandbox', 'sandbox']]) {
    const q = new URLSearchParams(location.search);
    if (m === 'sandbox') q.set('mode', 'sandbox'); else q.delete('mode');
    const s = q.toString();
    $(id).href = s ? `?${s}` : location.pathname;
    $(id).classList.toggle('on', m === MODE);
    if (m === MODE) $(id).setAttribute('aria-current', 'page');
  }
  C.opt.circlesAll = () => DEV;
  C.opt.worldFrame = () => true;   // the world view framed as the ship view, its corners rounded alike
  C.opt.focusCircle = () => true;   // a yellow circle around the ship's detected object, inside the frame
  C.opt.frameColor = () => '#fff';   // the world view's frame around the ship: white      // with the morph's colour dot, as in the morph book
  // Only the player's ship ([put_player]) is ever followed: no creature is focused by being put
  // in ([put_morph], [min_morph]) or by being near the centre after a cleared world.
  C.opt.autoFocus = () => false;
  C.opt.overFill = () => '#000';   // the ship view goes black while the ship is lost
  // "Dissolved" alone, centred in the ship view (the Respawn button below it), at 1.8 x the message
  // area's former text size (18 CSS px, 15 on a phone), at least 44 canvas px.
  C.opt.banner = () => {
    const css = matchMedia('(max-width:560px)').matches && !stacked ? 15 : 18;
    const s = css * C.logicalPerCss(true);
    return { title: Math.max(44, 1.8 * s), line: 0, inShip: true, dy: -40 };   // a little above the middle
  };

  // Text with _italic_, *bold* and *_italic bold_* (either nesting), {<button>:<key>} (the button's
  // text in a round-cornered box, then the key as a keycap; {<button>} or {:<key>}: one of them) and
  // \n (a line break) as nodes: plain text, <i>, <b>, <span>, <br>.
  function styled(text) {
    const out = [];
    let last = 0;
    for (const m of text.matchAll(/\{([^{}:]*)(?::([^{}]*))?\}|\*_([^*_]+)_\*|_\*([^*_]+)\*_|\*([^*]+)\*|_([^_]+)_/g)) {
      if (m.index > last) out.push(text.slice(last, m.index));
      last = m.index + m[0].length;
      if (m[0][0] === '{') {
        const b = m[1].trim(), k = (m[2] || '').trim();
        if (b) { const el = document.createElement('span'); el.className = 'msg-btn'; el.textContent = b; out.push(el); }
        if (k) { const el = document.createElement('span'); el.className = 'msg-key'; el.textContent = k; out.push(el); }
        continue;
      }
      m.splice(1, 2);   // the groups below as before: 1-2 both, 3 bold, 4 italic
      const both = m[1] ?? m[2], el = document.createElement(both !== undefined || m[3] !== undefined ? 'b' : 'i');
      if (both !== undefined) { const i = document.createElement('i'); i.textContent = both; el.append(i); }
      else el.textContent = m[3] ?? m[4];
      out.push(el);
    }
    if (last < text.length) out.push(text.slice(last));
    // \n (a backslash and n, as written in script.txt): a line break
    return out.flatMap(x => typeof x !== 'string' ? [x]
      : x.split('\\n').flatMap((t, k) => k ? [document.createElement('br'), t] : [t]));
  }

  // ---- messages -----------------------------------------------------------------------------
  // The message area under the world view.
  // hint: a grey note after the text, e.g. what key goes on.
  // One inner span holds the text and the hint, so they wrap as one paragraph.
  // The text's markup: _italic_, *bold*, and *_both_* (or _*both*_).
  // wait: shown in the waiting colour (the player has something to do).
  // A change fades: the old message fades out (if there is one), then the new one fades in
  // (MSG_FADE_MS each, as .msg's transition in style.css). A message given during the fade-out
  // replaces the one waiting to come in.
  const MSG_FADE_MS = 300;
  let msgNext = null, msgTimer = null;
  // spc: the hint ends with the mark that Space answers (marks.next), which gets its key badge.
  // extra: an element put after the hint (a condition's progress, after ⧗)
  function say(text, hint = '', wait = false, spc = false, extra = null) {
    const p = document.createElement('span');
    p.append(...styled(text));
    // a hint starting with a line break goes on its own row, under the text
    if (hint) {
      const h = document.createElement('span'); h.className = 'hint';
      let s = hint.startsWith('\n') ? hint.slice(1) : ' ' + hint;
      if (hint.startsWith('\n')) p.append(document.createElement('br'));
      if (spc && s.endsWith(marks.next)) {
        h.append(...styled(s.slice(0, -marks.next.length)));
        const m = document.createElement('span'); m.className = 'mark'; m.dataset.key = 'spc'; m.textContent = marks.next;
        h.append(m);
      } else h.append(...styled(s));
      if (extra) h.append(' ', extra);
      p.append(h);
    }
    msgNext = { p, wait, empty: !text && !hint };
    if (msgTimer) return;                       // fading out: the new one comes in after
    const el = $('msg');
    if (!el.textContent) { showMsg(); return; }  // nothing to fade out
    el.classList.add('fade');
    msgTimer = setTimeout(showMsg, MSG_FADE_MS);
  }
  // The message bar is one line of a fixed height (style.css): its font is --msg-base, or smaller if
  // the script's widest message, with what follows it (≫, or ⧗ and a condition's progress), would
  // not fit its width at that size. Measured on the page (the same fonts, bold and italic), once the
  // script is loaded and on every resize (a window resized or rotated, fullscreen in or out).
  // The progress a line's conditions show after ⧗ at its widest (CONDITIONS' t.progress), for fitMsg
  function progressSample(conds) {
    const PROG = { object: ([n]) => `morph: ${n}/${n}`, mass: ([, x]) => `mass: ${(+x).toFixed(1)}/${x}`,
                   survive_time: ([s]) => `time: ${s}/${s}` };
    return conds.filter(c => PROG[c.name]).map(c => ' ' + PROG[c.name](c.args)).join('');
  }
  function fitMsg() {
    const el = $('msg');
    if (!SCRIPT.length || !el.clientWidth) return;
    el.style.removeProperty('--msg-font');
    const cs = getComputedStyle(el), base = parseFloat(cs.fontSize), avail = el.clientWidth - parseFloat(cs.paddingRight) - 2;
    const probe = document.createElement('div');
    probe.style.cssText = 'position:absolute;left:0;top:0;visibility:hidden;font-size:100px;white-space:nowrap;line-height:1';
    let tip = null;   // a [next_msg]'s tip, for the next message that moves on by Space
    for (const L of SCRIPT) {
      const nm = L.actions.find(a => a.name === 'next_msg');
      if (nm) tip = nm.args[0] || MARKS.tip;
      if (!L.text) continue;
      const d = document.createElement('div');
      d.style.width = 'max-content';
      const t = !L.conds.length && tip ? ` ${tip}` : '';
      if (!L.conds.length) tip = null;
      d.append(...styled(L.text), ...styled(L.conds.length ? ` ${marks.act}${progressSample(L.conds)}` : `${t} ${marks.next}`));
      probe.append(d);
    }
    el.parentNode.append(probe);
    const widest = Math.max(0, ...Array.from(probe.children, d => d.getBoundingClientRect().width));
    probe.remove();
    // stacked: two lines (wrapping at words, so allow some slack)
    const room = stacked ? avail * 2 * 0.85 : avail;
    if (widest > 0 && widest * base / 100 > room) el.style.setProperty('--msg-font', `${Math.floor(10 * room * 100 / widest) / 10}px`);
  }
  window.addEventListener('resize', fitMsg);

  // ---- layout ---------------------------------------------------------------------------------
  // A portrait window, not fullscreen (a phone held upright): the stacked layout (core setLayout,
  // style.css :root.stacked): the world view at the top, the ship view at 1.5x below it, at the right
  // edge, with Fire and < > left of it at the screen's left edge (moved into its layer, placed
  // outside it), then the controls, then the message (two lines). Otherwise side by side. Checked
  // on every resize and fullscreen change.
  const BAND = 16, SHIP_SCALE = 1.5;   // (style.css :root.stacked has the same numbers: 768 x 1296)
  let stacked = false;
  function applyLayout() {
    const want = innerHeight > innerWidth && !document.documentElement.classList.contains('fs');
    if (want === stacked) return;
    stacked = want;
    document.documentElement.classList.toggle('stacked', stacked);
    C.setLayout(stacked, BAND, SHIP_SCALE, 1);   // the ship view at the right edge
    const layer = document.querySelector(stacked ? '.layer-ship' : '.layer-world');
    layer.append($('bullet-nav'), $('fire'));
    let col = $('dev-col');
    if (!col) {   // (made here if the page lacks it, e.g. an index.html cached from before it)
      col = document.createElement('div'); col.id = 'dev-col'; col.className = 'dev-col';
      document.querySelector('.layer-ship').append(col);
    }
    (stacked ? col : layer).append($('dev'), $('dev-fps'));   // ?dev=1: left of the ship view, a column
    fitMsg();
  }
  window.addEventListener('resize', applyLayout);
  new MutationObserver(applyLayout).observe(document.documentElement, { attributes: true, attributeFilter: ['class'] });   // (the first call: at the start, below)

  function showMsg() {
    msgTimer = null;
    const el = $('msg'), m = msgNext;
    msgNext = null;
    el.classList.add('fade');
    el.replaceChildren(m.p);
    el.classList.toggle('wait', m.wait);
    if (m.empty) return;
    void el.offsetWidth;   // start the fade-in from opacity 0
    el.classList.remove('fade');
  }

  // The script (script.txt), run line by line. A line is a message, then any number of [tokens]:
  // [name] or [name:arg:arg...], each an action (ACTIONS, run when the line is reached, in order)
  // or a condition, [need:name] or [need:name:arg...] (CONDITIONS, checked every frame; the line
  // waits until all of them hold).
  //   - A line with text shows it. With no condition it moves on by Space or a click on the
  //     message bar; with conditions, when they hold (Space and clicks then do nothing).
  //   - A line of tokens only shows nothing (the last message stays): it runs its actions, then
  //     moves on at once, or when its conditions hold.
  // Message text may use _italic_, *bold* and *_italic bold_* (say / styled).
  // Blank lines and comment lines (starting with #) are skipped. Unknown names are reported in the console
  // and ignored. The last line stays. A message ends with a grey mark: one waiting on a condition
  // with marks.act (the player has something to do); else marks.next (Space / click moves on;
  // none on the last line), the first such message with marks.tip before it ([next_msg] sets
  // them, and shows the tip again on the next such message).
  const MARKS = { tip: 'Press space or click message to continue', next: '\u226b', act: '\u29d7' };   // ≫ ⧗
  const marks = { ...MARKS };
  let tipPending = true;   // the tip goes on the next message that moves on by Space / a click

  // Visible morphs: the names the player may see. Every morph is still detected, but one not
  // visible shows as "???" (name box, labels). [species] empties the list, [allow_morph:<ids>]
  // sets it; "all" stands for every morph of the species (also ones detected from the recorded
  // stats only, as paralimus).
  const vis = { all: false, allow: new Set() };
  const isVisible = name => vis.all || vis.allow.has(name);
  C.opt.morphLabel = det => isVisible(det.name) ? det.name : null;   // a fuzzy morph too: no "~"
  // The ship is the one followed: a click on the world view does not switch to another creature.
  C.opt.clickFocus = () => false;
  C.opt.noFocusText = () => [];   // nothing followed: the ship view stays empty
  // The player is dead when the followed ship is gone (other creatures may live on): "Dissolved"
  // and Respawn. (steps > 0: not in the moment between a cleared world and the ship put in.)
  const shipLost = () => C.status !== 'alive' || (!C.focus && C.steps > 0);
  // [clear_world:empty]: a scene of empty space, no player: while the ship is lost the ship view
  // shows the world's centre, with no "Dissolved", Respawn or Choose; until the next [clear_world]
  let emptySpace = false;
  const lostShown = () => shipLost() && !emptySpace;   // "Dissolved" and Respawn / Choose show
  C.opt.over = lostShown;
  C.opt.emptyView = () => emptySpace;
  C.opt.pauseWhenEmpty = () => false;   // "Dissolved" does not pause the game
  // morph ids from a token's args: [allow_morph:limus:gyrans] or [allow_morph:limus,gyrans]
  const morphIds = args => args.flatMap(a => a.split(/[,\s]+/)).filter(Boolean);

  // ---- morph book -----------------------------------------------------------------------------
  // The morphs the player has discovered, empty at first. [discover_morph:<ids>] adds to it;
  // while the book is on, a morph is also added by being found: a visible morph not yet
  // discovered that the followed ship has been detected as for DISCOVER_STEPS (a game second) of play in a row; that
  // also opens the book, titled "New morph discovered!", the new morph marked.
  // [morph_book:on] shows the Morphs button (it opens the book). The book has a slot per visible
  // morph, in the visible list's order: the name if discovered, else "?"; a click on a discovered
  // one does the [morph_book_action], but only while the ship is lost (else read-only). The game pauses while it is open. Close, a click outside the
  // box, or Escape closes it.
  // per species: a morph name is only a morph of its species (5N7KKM tardus is not GDNQYX tardus),
  // so each species keeps its own list; [species] switches to its list (kept for when it comes back)
  const discoveredBy = new Map();
  let speciesName = '';   // the species' name (genus and species), as [species:<name>:<id>] gives it (the morph book shows it)
  let discovered = [];   // { name, hue }
  const book = { mode: 'morph', on: false, open: false, wasRunning: false, newName: null, fromDiscovery: false, action: 'clear_put',
                 found: [] };   // the morphs found in play, each added when the book its discovery opened closes ([need:discover_morph])
  // The visible morphs, in order: the allowed ones as allowed; with "all", then every other morph
  // of the species (its patterns, then the names in the recorded stats).
  function visibleList() {
    const out = [...vis.allow];
    if (vis.all) {
      const sp = C.speciesOf(C.pattern);
      const names = PATTERNS.filter(p => C.speciesOf(p) === sp && !p.hidden).map(p => C.morphName(p).replace(C.BOOM, ''))
        .concat(C.stats.long.segs.map(sg => sg.morph));
      for (const n of names) if (n && !out.includes(n)) out.push(n);
    }
    return out;
  }
  // One slot per visible morph, in the visible list's order: if discovered, its picture (when it
  // was found in play) above its name, else "?".
  // Title: "New morph discovered!" when a discovery opened the book.
  function renderBook() {
    $('book-grid').classList.toggle('bullets', book.mode === 'bullet');
    $('book-grid').classList.toggle('species', book.mode === 'species');
    document.querySelector('.book-box').classList.toggle('wide', book.mode !== 'morph');   // the bullet and species books: wider
    $('book-close').dataset.key = 'Esc';   // only Escape closes a book (its letters pick slots)   // its key badge: the key that closes this book
    // the morph book's species, at its lower right (a morph is only a morph of its species)
    $('book-species').hidden = book.mode !== 'morph' || !C.pattern;
    if (C.pattern) {
      const code = document.createElement('i');
      code.textContent = speciesName || C.speciesOf(C.pattern);   // no name: the id
      $('book-species').replaceChildren('Species: ', code);
    }
    // sandbox mode: the species' name goes to the species book
    $('book-species').classList.toggle('link', MODE === 'sandbox');
    $('book-species').title = MODE === 'sandbox' ? 'Choose another species (C)' : '';
    $('book-clear').hidden = MODE !== 'sandbox';
    if (book.mode === 'bullet') return renderBullets();
    if (book.mode === 'species') return renderSpecies();
    allDiscovered();
    const canPick = book.action !== 'none' && shipLost();   // a click on a discovered morph does something
    const canClone = clone.on && !shipLost() && !book.fromDiscovery;   // [clone:on]: pick one to clone (not in the discovery's book)
    setTitle(book.fromDiscovery ? ['New morph discovered!'] : canPick ? ['Choose a morph to ', ['respawn', 'choose']]
      : canClone ? ['Choose a morph to ', ['clone', 'clone']] : ['Dissolve to respawn']);   // read-only: how to get a pick
    const grid = $('book-grid');
    grid.replaceChildren(...visibleList().map((name, slot) => {
      const d = discovered.find(x => x.name === name), el = document.createElement('div');
      el.className = 'book-slot' + (d ? '' : ' empty') + (d && d.name === book.newName ? ' new' : '');
      if (d && d.name === book.newName && book.reveal) revealSlot(el, 0);
      // Every slot is the same size: a picture box (SNAP_CELLS across, at the world view's scale)
      // above a label line.
      const scale = $('gl').getBoundingClientRect().width / $('gl').width;   // CSS px per canvas px
      const pic = document.createElement('div'), label = document.createElement('span');
      pic.className = 'pic';
      pic.style.width = Math.round(SNAP_CELLS * C.WS) * scale + 'px';   // as the pictures (snapshot), whatever the zoom
      label.className = 'label';
      if (d) {
        const dot = document.createElement('span');
        dot.className = 'dot';
        dot.style.background = Number.isFinite(d.hue) ? `hsl(${d.hue.toFixed(0)},100%,55%)` : '#ccc';
        // the picture, taken at its discovery, at the world view's scale on screen; none (e.g. put
        // in by [discover_morph:<ids>]): a black box
        if (d.img) {
          const img = new Image();
          img.src = d.img;
          img.alt = d.name;
          img.onload = () => { img.style.width = img.naturalWidth * scale + 'px'; };
          pic.append(img);
        } else pic.classList.add('blank');
        label.append(dot, d.name);
        // a pick only while the ship is lost (else the book is read-only)
        if (canPick || canClone) {
          el.classList.add('pick');
          if (canPick) {
            if (d.name === playerMorph) el.classList.add('current');   // what Respawn would put back: a lighter background
            el.title = book.action === 'clear_put' ? `Start again with ${d.name}` : `Put a ${d.name} into the world`;
            el.onclick = () => pickMorph(d.name);
          } else {   // clone mode: pink, and the morph follows the pointer to be pasted
            el.classList.add('clone');
            el.title = `Clone ${d.name}: click on the ship view to paste it`;
            el.onclick = () => startClone(d.name);
          }
          if (SLOT_KEYS.morph[slot]) el.dataset.key = SLOT_KEYS.morph[slot];   // its key picks it too (badge while Shift / Alt is held)
        }
      } else label.textContent = '\u00a0';   // an empty label line; the big ? is drawn over the slot (style.css)
      el.append(pic, label);
      return el;
    }));
    fitBook();
  }
  // The book never scrolls: if it is taller than the game window (two rows or more of pictures),
  // shrink every picture box alike until it fits. Only measurable while it is shown.
  // (The species book's thumbnail grids alike, each narrowed by its own height / width.)
  function fitBook() {
    const box = document.querySelector('.book-box'), pics = $('book-grid').querySelectorAll('.pic, .thumbs');
    if ($('book').hidden || !pics.length) return;
    const over = box.scrollHeight - $('book').clientHeight * 0.98;   // scrollHeight: its full content (it is clipped at 98%)
    if (over <= 0) return;
    const rows = Math.ceil(pics.length / 4), r = pics[0].getBoundingClientRect(), aspect = r.height / r.width || 1;
    const fit = Math.max(16, r.width - (over / rows) / aspect - 1);
    pics.forEach(p => { p.style.width = fit + 'px'; });
  }
  // A click on a discovered morph in the book: close it, do the morph_book_action, and play (even
  // if the game was paused before the book opened).
  let respawns = 0;   // ships put back after dissolving, by Respawn or Choose ([need:respawn])
  function pickMorph(name) {
    closeBook();
    // either way the new creature is followed
    if (book.action === 'clear_put') { C.clearWorld(); ACTIONS.put_player(name, 'center'); }   // (not the script's [clear_world]: background and [min_morph] kept)
    else if (book.action === 'put') ACTIONS.put_player(name, 'random');
    if (book.action !== 'none') { retake = true; cand = null; respawns++; }   // its picture may be retaken; [need:respawn]
    C.setRunning(true);
  }
  // ---- clone mode ([clone:on]) ------------------------------------------------------------
  // While it is on and the ship is alive, the morph book offers the discovered morphs to clone.
  // Picking one closes the book; the morph then follows the pointer over the ship view (instead
  // of the tool circle), turned by ↻ (45 degrees clockwise) and flipped by ↔ (left to right, as
  // seen), until a click on the ship view pastes it there, or the Morphs button cancels it. Clone
  // mode stays allowed until [clone:off].
  const newBullets = new Set();   // the patterns created in the bullet book, green in it until it is closed
  let bulletNew = false;   // a new pattern can be created ([new_bullet]): the Bullets button reads "Create", green
  let revealName = null;   // the pattern just created: revealed as a new card as the book re-renders
  let clonesPlaced = 0;   // clones pasted into the world so far ([need:first_clone] counts on from its line)
  const clone = { on: true, held: null };   // on: allowed (the default, until [clone:off])   // held: { name, M (board-frame 2x2), src form, form, img }
  function startClone(name) {
    const p = findMorph(name);
    if (!p) return;
    // as before the book opened; picking again while holding one: as before the first pick
    const wasRunning = clone.held ? clone.held.wasRunning : book.wasRunning;
    closeBook();
    C.setRunning(false);   // paused until it is placed (or cancelled), then as it was
    clone.held = Object.assign(Clone.hold(p, C.R, name), { wasRunning });
    syncButtons();
  }
  // Drop the held clone (pasted or cancelled); off: also disallow clone mode ([clone:off]).
  function endClone(off = false) {
    const was = clone.held && clone.held.wasRunning;
    if (off) clone.on = false;
    clone.held = null; syncButtons();
    if (book.open && book.mode === 'morph') renderBook();
    if (was) C.setRunning(true);
  }
  // the ghost keeps following the pointer over ↻ / ↔ (they sit on the ship view)
  C.opt.keepHover = el => !!clone.held && !!el && $('clone-tools').contains(el);
  $('clone-tools').addEventListener('pointermove', e => { if (clone.held) C.hoverAt(e); });
  // ↻: 45 degrees clockwise; ↔: left to right as seen in the ship view (clone.js)
  function rotateClone() { if (clone.held) Clone.rotate(clone.held); }
  function flipClone() { if (clone.held) Clone.flip(clone.held, C.viewAngle()); }
  // the ghost over the ship view, drawn as the view shows the board (phi: board = rot(view, phi)):
  // where it is placed (a touch, awaiting ✓), else under the hovering mouse
  C.opt.drawGhost = (ctx, hoverAt, centre, cellPx, phi, focus) => {
    const H = clone.held;
    if (!H) return false;
    let x, y;
    if (H.at) {
      const b = [C.wrapD(H.at[0] - focus.col, C.BW), C.wrapD(H.at[1] - focus.row, C.BH)];
      const c = Math.cos(-phi), s = Math.sin(-phi);   // board -> view
      x = centre[0] + (b[0] * c - b[1] * s) * cellPx; y = centre[1] + (b[0] * s + b[1] * c) * cellPx;
    } else if (hoverAt) [x, y] = hoverAt;
    else return true;   // nothing to show, and no tool circle either
    Clone.draw(ctx, H, x, y, phi, cellPx);
    return true;
  };
  // A click on the ship view pastes the held clone there (centred on the pointer). A touch (no
  // mouse, so no preview) places it there instead, shown, to be turned / flipped / moved by another
  // touch, and pasted by ✓.
  // (a touch dragged: the copy follows it, placed all along)
  C.opt.zoomDrag = (off, e) => {
    const H = clone.held;
    if (!H || !H.at || !C.focus || (e && e.pointerType === 'mouse')) return;
    H.at = [C.focus.col + off[0], C.focus.row + off[1]];
  };
  C.opt.zoomPress = (off, e) => {
    const H = clone.held;
    if (!H || !C.focus) return false;
    const at = [C.focus.col + off[0], C.focus.row + off[1]];
    if (e && e.pointerType !== 'mouse') { H.at = at; syncButtons(); return true; }
    pasteClone(at);
    return true;
  };
  // The player boards the clone: the game follows it from now on (the old ship goes on as another
  // creature; pasted onto it, the two are one object, followed as before), and its morph is what
  // Respawn puts back. So a morph can be changed without dissolving, and a twin engineered by
  // becoming one half and pasting the other next to it.
  function pasteClone(at) {
    const H = clone.held;
    Clone.paste(H, at);
    C.analyze(0);
    const i = C.found ? Creatures.pick(C.found, at[0], at[1], 2 * C.R, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
    if (i >= 0) { C.setFocus(i); playerMorph = H.name; retake = true; cand = null; }
    clonesPlaced++;   // [need:first_clone]
    endClone();
  }
  $('clone-rot').onclick = rotateClone;
  $('clone-flip').onclick = flipClone;
  $('clone-cancel').onclick = () => endClone();
  $('clone-ok').onclick = () => { if (clone.held && clone.held.at) pasteClone(clone.held.at); };

  // The bullet book shares the morph book's window (book.mode 'bullet'): every bullet pattern by
  // name, 4 to a row (as the morph book); the ones the gun has ([bullet]) can be
  // clicked to switch to them (the one in use with a yellow border), the others are "?".
  // Keys picking a slot, by its place: the books' rows of 4 on the keyboard's rows (morphs: ASDF,
  // ZXCV; bullets: 1234, QWER, ASDF, ZX..).
  const SLOT_KEYS = { morph: [...'ASDFZXCV'], bullet: [...'1234QWERASDFZXCV'], species: [...'QWERASDFZX'] };
  // [new_bullet:<n>]: creating a pattern costs n energy (the score): how many more can be created
  // now, with the energy there is (none while it is below n), at most as many as are not listed
  let bulletsCreated = 0;   // patterns created this run ([need:create_bullet])
  function bulletCredits() {
    if (!newBullet || !Bullets.on) return 0;
    const have = Bullets.list.length, all = Object.keys(Bullets.KINDS).length;
    return Math.max(0, Math.min(all - have, Math.floor(energy / newBullet.every)));
  }
  // The Bullets button reads "Create" (green) while one can be created; turning so, with a glint
  function syncNew() {
    const was = bulletNew;
    bulletNew = bulletCredits() > 0;
    if (bulletNew === was) return;
    syncButtons();
    if (bulletNew) revealButton($('bullets'));
  }
  // A "?" slot clicked while one can be created: that pattern is added to the list and the gun
  // switches to it, revealed as a new card (green), turning yellow (in use) once revealed; the
  // book stays open
  function createBullet(name) {
    if (bulletCredits() <= 0 || !Bullets.addRandom(name)) return;
    bulletsCreated++;
    payBullets(1);   // its price, from the energy
    newBullets.add(name);
    Bullets.select(name);
    revealing.add(name);
    revealName = name; renderBook(); revealName = null;
    scrollToReveal();
    syncNew();
    setTimeout(() => {
      revealing.delete(name);
      if (book.open && book.mode === 'bullet') renderBook();
    }, REVEAL_FLIP_MS);
  }
  // Each slot: yellow border, the one in use; cyan, another listed (a click fires it); green, one
  // created since the book opened; "?", not listed (cyan, a click creates it, while one can be:
  // "Create or choose bullet").
  // [bullet:<list>:<name>]: open the bullet book on the patterns just added, revealed as new cards
  // (green); the one equipped turns yellow once they are all revealed
  const REVEAL_FLIP_MS = 750;   // as .book-slot.reveal's card-flip (style.css)
  const revealing = new Set();
  function showNewBullets(names) {
    for (const n of names) newBullets.add(n);
    revealing.clear(); names.forEach(n => revealing.add(n));
    openBulletBook();
    setTimeout(() => {
      revealing.clear();
      if (book.open && book.mode === 'bullet') renderBook();
    }, ZAP.length * ZAP_MS + (names.length - 1) * REVEAL_STAGGER_MS + REVEAL_FLIP_MS);
  }
  function renderBullets() {
    const create = bulletCredits() > 0;
    setTitle(create ? [['Create', 'create'], ' or ', ['use', 'choose'], ' a bullet type'] : [['Use', 'choose'], ' a bullet type']);
    const have = Bullets.list;
    let nNew = 0;
    $('book-grid').replaceChildren(...Object.keys(Bullets.KINDS).map((name, slot) => {
      const el = document.createElement('div'), ok = have.includes(name);
      const current = name === Bullets.kind && !revealing.has(name);   // (yellow once revealed)
      el.className = 'book-slot' + (ok ? ' pick' : ' empty') + (current ? ' current' : '');   // (a new one: revealed, no colour of its own)
      if (newBullets.has(name) && book.reveal) revealSlot(el, nNew++);
      if (name === revealName) revealSlot(el, 0, 0);   // just created: at once
      const label = document.createElement('span');
      label.className = 'label';
      if (ok) label.append(bulletPic(name), name);   // its stored pattern, small, in white, before the name
      else label.textContent = '\u00a0';
      el.append(label);
      if (ok) {
        el.title = `Fire ${name}`; el.onclick = () => { Bullets.select(name); closeBook(); };
        if (SLOT_KEYS.bullet[slot]) el.dataset.key = SLOT_KEYS.bullet[slot];
      } else if (create) {
        el.classList.add('pick', 'create');
        el.title = 'Create this new bullet'; el.onclick = () => createBullet(name);
        if (SLOT_KEYS.bullet[slot]) el.dataset.key = SLOT_KEYS.bullet[slot];
      }
      return el;
    }));
  }
  // A small picture of a bullet pattern as stored (bullet_patterns.js): white, its cell values as
  // opacity, on nothing; the same scale for all (BULLET_PIC_PX css px per cell, then capped to
  // BULLET_PIC_MAX across), made once per pattern.
  const BULLET_PIC_PX = 0.9, BULLET_PIC_MAX = 34, bulletPics = new Map();
  function bulletPic(name) {
    if (!bulletPics.has(name)) {
      const rows = Rule.rle2cells(Bullets.KINDS[name].rle), h = rows.length, w = Math.max(...rows.map(r => r.length));
      const cv = document.createElement('canvas');
      cv.width = w; cv.height = h;
      const g = cv.getContext('2d'), im = g.createImageData(w, h);
      rows.forEach((r, y) => r.forEach((v, x) => { const i = (y * w + x) * 4; im.data.set([255, 255, 255, Math.round(255 * v)], i); }));
      g.putImageData(im, 0, 0);
      const k = Math.min(BULLET_PIC_PX, BULLET_PIC_MAX / Math.max(w, h));
      bulletPics.set(name, { src: cv.toDataURL(), w: w * k, h: h * k });
    }
    const p = bulletPics.get(name), img = new Image();
    img.src = p.src; img.alt = ''; img.className = 'bpic';
    img.style.width = p.w + 'px'; img.style.height = p.h + 'px';
    return img;
  }
  // ---- species book (sandbox mode) ------------------------------------------------------------
  // The Species button (C), or the species' name in the morph book: a slot per visible species
  // ([add_allow_species]), in order: a small thumbnail of each of its visible morphs in a row (the
  // centre of its picture if discovered, else black), and its name. The current one has a
  // yellow border. A click chooses it (chooseSpecies: its [on_choose_species] block sets the world up).
  // A species' visible morphs: its [add_allow_species] list, or for all, every morph of it (its patterns)
  const speciesMorphs = e => !e.morphs.includes('all') ? e.morphs
    : [...new Set(PATTERNS.filter(p => C.speciesOf(p) === e.id && !p.hidden).map(p => C.morphName(p).replace(C.BOOM, '')))];
  // A book's title: plain text, and its verbs ([text, kind]) in the colour of the slots that do
  // that (style.css .verb-*): choose, cyan; create, green; clone, pink
  function setTitle(parts) {
    $('book-title').replaceChildren(...parts.map(p => {
      if (typeof p === 'string') return p;
      const s = document.createElement('span'); s.className = 'verb-' + p[1]; s.textContent = p[0]; return s;
    }));
  }
  function openSpeciesBook() {
    if (MODE !== 'sandbox') return;
    if (book.open && book.mode !== 'species') closeBook();
    book.mode = 'species';
    showBook();
  }
  function renderSpecies() {
    setTitle([['Choose', 'choose'], ' a species']);
    const cur = C.pattern && C.speciesOf(C.pattern);
    $('book-grid').replaceChildren(...visSpecies.map((e, slot) => {
      const found = discoveredBy.get(e.id) || [];
      const el = document.createElement('div'), thumbs = document.createElement('div'), label = document.createElement('span');
      el.className = 'book-slot pick species' + (e.id === cur ? ' current' : '');
      thumbs.className = 'thumbs';
      for (const m of speciesMorphs(e)) {
        const d = found.find(x => x.name === m), t = document.createElement('div');
        t.className = 'thumb' + (d ? '' : ' unknown');   // not discovered: the box's colour (discovered with no picture: black)
        if (d) { t.title = m; if (d.img) { const img = new Image(); img.src = d.img; img.alt = m; t.append(img); } }
        thumbs.append(t);
      }
      label.className = 'label';
      label.textContent = e.name || e.id;
      el.append(thumbs, label);
      el.title = `Choose ${e.name || e.id}`;
      el.onclick = () => pickSpecies(e.id);
      if (SLOT_KEYS.species[slot]) el.dataset.key = SLOT_KEYS.species[slot];
      return el;
    }));
    fitBook();
  }
  let chosenSpecies = null;   // the species last chosen in the book: Retry's replay chooses it again
  // As a [new_scene]: the views fade to black (WORLD_FADE_MS), the species is set up unseen, then
  // they fade back in (as Retry does)
  let speciesFading = false;
  function pickSpecies(id) {
    closeBook();
    if (clone.held) endClone();
    if (speciesFading) return;
    speciesFading = true;
    fadeWorld(1);
    setTimeout(() => {
      chosenSpecies = id;
      chooseSpecies(id);
      retake = true; cand = null;
      sceneHidden = false;
      fadeWorld(0);
      C.setRunning(true);
      speciesFading = false;
    }, WORLD_FADE_MS);
  }
  function openBulletBook() {
    if (!Bullets.on) return;
    if (book.open && book.mode === 'morph') closeBook();
    book.mode = 'bullet';
    showBook();
  }
  function openBook(newName = null) {
    if (chapterShowing()) return;   // never with a chapter card (coming, up or fading out)
    book.mode = 'morph';
    book.newName = newName; book.fromDiscovery = book.fromDiscovery || !!newName;
    showBook();
  }
  // A new morph or bullet, in the book opened to show it: revealed once as a new card (style.css
  // .reveal: flipped over, a glint across it, a green glow), after the zap in; several in turn.
  const REVEAL_STAGGER_MS = 180;
  function revealSlot(el, k, base = ZAP.length * ZAP_MS) {   // base: the zap in, by default
    el.classList.add('reveal');
    el.style.setProperty('--reveal-delay', base + k * REVEAL_STAGGER_MS + 'ms');
  }
  // The book's list scrolling sideways (the stacked layout): brings the first card being revealed
  // to its middle, so the reveal is seen
  function scrollToReveal() {
    const grid = $('book-grid'), el = grid.querySelector('.book-slot.reveal');
    if (!el || grid.scrollWidth <= grid.clientWidth) return;
    const g = grid.getBoundingClientRect(), r = el.getBoundingClientRect();
    grid.scrollLeft += r.left - g.left - (g.width - r.width) / 2;
  }
  // The same reveal on a button (the Bullets button turning "Create"), once, from now
  function revealButton(el) {
    el.classList.remove('reveal');
    void el.offsetWidth;   // restart the animation
    el.classList.add('reveal');
    el.addEventListener('animationend', e => { if (e.animationName === 'card-glow') el.classList.remove('reveal'); });
  }
  function showBook() {
    if (chapterShowing()) return;
    book.reveal = !book.open;   // only as it opens: a later re-render shows the slots as they are
    renderBook();
    book.reveal = false;
    if (book.open) return;
    book.open = true; book.wasRunning = C.running;
    C.setRunning(false);
    $('book').hidden = false;
    fitBook();
    scrollToReveal();
    zap(true);
    $('book-close').focus();
  }
  function closeBook() {
    if (!book.open) return;
    book.open = false;
    if (book.mode === 'bullet') { newBullets.clear(); revealing.clear(); }   // seen: no longer new (the button stays "Create" while more can be created)
    zap(false);
    if (book.fromDiscovery && book.newName) book.found.push(book.newName);
    book.fromDiscovery = false; book.newName = null;
    if (book.wasRunning) C.setRunning(true);
  }
  // "Zap in": ZAP.length frames, one every ZAP_MS, growing and firming up to the box's own size,
  // stack up over the game window, then the box shows, the window dims behind it (not before),
  // and the frames go. "Zap out" is the same backwards: the box goes and the window undims, the
  // frames show stacked, then they go from the largest down, and the book is hidden.
  const ZAP = [[4, 0.25], [3, 0.45], [2, 0.65], [1, 0.85], [0, 1]];   // [steps of ZAP_STEP px smaller than the box, opacity]
  const ZAP_STEP = 26;   // px, in width and height alike
  const ZAP_MS = 45;
  let zapTimers = [];
  function zap(opening) {
    zapTimers.forEach(clearTimeout); zapTimers = [];
    const bookEl = $('book'), box = bookEl.querySelector('.book-box');
    bookEl.querySelectorAll('.zap').forEach(z => z.remove());
    const r = box.getBoundingClientRect();
    if (opening) box.style.visibility = 'hidden';
    bookEl.classList.add('zapping');   // no dimming while the frames show
    const frames = ZAP.map(([k, a]) => {
      const z = document.createElement('div');
      z.className = 'zap';
      z.style.width = Math.max(0, r.width - k * ZAP_STEP) + 'px'; z.style.height = Math.max(0, r.height - k * ZAP_STEP) + 'px';
      z.style.opacity = a;
      if (opening) z.style.visibility = 'hidden';
      bookEl.append(z);
      return z;
    });
    const at = (k, f) => zapTimers.push(setTimeout(f, k * ZAP_MS));
    if (opening) {
      frames.forEach((z, k) => at(k, () => { z.style.visibility = ''; }));
      at(frames.length, () => { box.style.visibility = ''; bookEl.classList.remove('zapping'); frames.forEach(z => z.remove()); });
    } else {
      box.style.visibility = 'hidden';
      frames.slice().reverse().forEach((z, k) => at(k + 1, () => z.remove()));
      at(frames.length + 1, () => { bookEl.hidden = true; bookEl.classList.remove('zapping'); box.style.visibility = ''; });
    }
  }
  function setBook(on) { book.on = on; syncButtons(); if (!on) closeBook(); }
  $('morphs').onclick = () => openBook();   // (holding a clone: the book to pick another one to clone)
  $('bullets').onclick = () => openBulletBook();
  $('species').onclick = () => openSpeciesBook();
  $('book-species').onclick = () => { if (book.mode === 'morph') openSpeciesBook(); };
  $('book-close').onclick = closeBook;
  $('book').addEventListener('pointerdown', e => { if (e.target === $('book')) closeBook(); });   // outside the box
  // Each frame: has the followed ship been a visible morph not discovered yet, without a break, for
  // DISCOVER_STEPS (a game second, in steps: paused time does not count)? A passing detection is
  // not a discovery. A discovered one gets a new picture (after SNAPSHOT_STEPS) only once after each
  // Respawn or book pick (retake), when the ship is usually alone in a cleared world.
  const DISCOVER_STEPS = C.STEPS_PER_S, SNAPSHOT_STEPS = C.STEPS_PER_S;
  let cand = null, lastStep = 0;   // cand: { name, n: steps } the morph the ship is in now, if it would be new
  let retake = false;   // set by Respawn and a book pick: the next known morph held gets a new picture
  function checkDiscovery() {
    const dn = C.simSteps - lastStep;
    lastStep = C.simSteps;
    if (chapterShowing()) return;   // no discovery (nor its book) under a chapter card; it waits
    const d = C.focus && C.stats.short.det;
    const known = d && d.name && discovered.find(x => x.name === d.name);
    if (known && !Number.isFinite(known.hue) && Number.isFinite(d.hue)) { known.hue = d.hue; saveInventory(); }   // put in unseen: its colour now
    // a new discovery; or a known morph, whose picture is retaken quietly after SNAPSHOT_STEPS (once
    // the detection has settled), replacing the old one: only once after a Respawn or a book pick
    if (!book.on || !d || !d.name || !isVisible(d.name)) { cand = null; return; }
    if (known && known.img && !retake) return;   // (one with no picture yet, from [discover_morph], gets it any time)
    if (!cand || cand.name !== d.name) { cand = { name: d.name, n: 0, taken: false }; return; }
    if (cand.taken) return;
    cand.n += dn;
    if (cand.n < (known ? SNAPSHOT_STEPS : DISCOVER_STEPS)) return;
    cand.taken = true;   // kept, so the same stretch as this morph takes no more pictures
    retake = false;
    if (known) { known.img = snapshot(); known.given = false; saveInventory(); if (book.open) renderBook(); return; }   // (a given one held: earned now)
    discover(d.name, d.hue, snapshot());
    openBook(d.name);
  }
  // Add a morph to the discovered list (its hue, the book's dot colour, once it has been detected;
  // img: its picture in the book, if it was found in play).
  // ?all=1: the visible morphs not discovered yet are, quietly (no picture yet: one is taken the
  // first time the ship holds each)
  function allDiscovered() {
    if (!ALL || !C.pattern) return;
    for (const n of visibleList()) if (findMorph(n) && !discovered.some(x => x.name === n)) discovered.push({ name: n, hue: NaN, img: null });
  }
  function discover(name, hue = NaN, img = null) {
    if (!discovered.some(x => x.name === name)) { discovered.push({ name, hue, img }); saveInventory(); }
  }

  // ---- inventory ------------------------------------------------------------------------------
  // The inventory: the discovered morphs of every species (discoveredBy: name, colour, picture)
  // and the bullet patterns the gun has (Bullets.list), and the energy (the dots score; sandbox
  // only: a sandbox run starts from the central energy and saves what it gains and spends; a story
  // run's energy, from 0, is its own). The central inventory, in localStorage
  // (INV_KEY), gathers them from every run: whatever a run gains is merged into it as it comes
  // (a discovery, a new picture or colour, a bullet added or created), never removed. A story run
  // starts from an empty inventory; a sandbox run starts from the central one (loaded at the start,
  // before sandbox.txt runs), and saves into it as a story run does.
  const INV_KEY = 'cognoship-inventory';
  function loadInventory() {
    try {
      const o = JSON.parse(localStorage.getItem(INV_KEY));
      if (o && typeof o === 'object') return { morphs: o.morphs && typeof o.morphs === 'object' ? o.morphs : {},
        bullets: Array.isArray(o.bullets) ? o.bullets : [], energy: Number.isFinite(o.energy) ? Math.max(0, o.energy) : 0 };
    } catch (e) {}
    return { morphs: {}, bullets: [], energy: 0 };
  }
  // The energy: what the ship has eaten (dots.js counts it: its score's rises), less what bullet
  // patterns cost. energyUnsaved (sandbox only): its change since it was last saved, added to the
  // central energy (so two tabs add up).
  let energy = 0, dotsSeen = 0, energyUnsaved = 0;
  function changeEnergy(d) {
    d = Math.max(-energy, d);
    energy += d;
    if (d && MODE === 'sandbox') { energyUnsaved += d; saveInventory(); }   // (a story run's energy is its own)
    showEnergy();
  }
  // n bullet patterns gained (created, or given by [bullet]): [new_bullet:<n>]'s price each, from
  // the energy (none while [new_bullet] is off; the energy never below 0)
  function payBullets(n) { if (newBullet && n > 0) changeEnergy(-n * newBullet.every); }
  let invTimer = null;
  const unsavedBullets = new Set();   // patterns [bullet:all] gave: not saved
  // Merge this run's inventory into the central one, soon (several changes in a frame: one write)
  let invCleared = false;   // the central inventory was just cleared: nothing more is saved (the page reloads)
  function saveInventory() { if (!invTimer && !invCleared && !ALL) invTimer = setTimeout(flushInventory, 100); }
  function flushInventory() {
    clearTimeout(invTimer); invTimer = null;
    if (invCleared) return;
    const inv = loadInventory();
    for (const [sp, list] of discoveredBy) {
      const c = Array.isArray(inv.morphs[sp]) ? inv.morphs[sp] : (inv.morphs[sp] = []);
      for (const d of list) {
        if (d.given) continue;   // [discover_morph]: not earned
        const e = c.find(x => x.name === d.name);
        if (!e) c.push({ name: d.name, hue: Number.isFinite(d.hue) ? d.hue : null, img: d.img || null });
        else {
          if (d.img) e.img = d.img;   // the latest picture
          if (Number.isFinite(d.hue)) e.hue = d.hue;
        }
      }
    }
    for (const b of Bullets.list) if (!unsavedBullets.has(b) && !inv.bullets.includes(b)) inv.bullets.push(b);
    inv.energy = Math.max(0, inv.energy + energyUnsaved); energyUnsaved = 0;
    try { localStorage.setItem(INV_KEY, JSON.stringify(inv)); }
    catch (e) { console.warn('Cognoship: the inventory could not be saved', e); }
  }
  window.addEventListener('pagehide', () => { if (invTimer) flushInventory(); });
  // Sandbox mode: "Hold to clear data" in the books, held for CLEAR_HOLD_MS, empties the central
  // inventory, then reloads the page: a fresh sandbox run from nothing. Let go sooner: nothing happens.
  const CLEAR_HOLD_MS = 3000;
  let clearTimer = null;
  $('book-clear').addEventListener('pointerdown', e => {
    e.preventDefault();
    $('book-clear').setPointerCapture?.(e.pointerId);
    $('book-clear').classList.add('holding');
    clearTimer = setTimeout(() => {
      clearTimer = null;
      invCleared = true; clearTimeout(invTimer); invTimer = null;
      try { localStorage.removeItem(INV_KEY); } catch (err) {}
      $('book-clear').classList.remove('holding'); $('book-clear').classList.add('done');
      $('book-clear').textContent = 'Data cleared';
      setTimeout(() => location.reload(), 600);
    }, CLEAR_HOLD_MS);
  });
  const clearCancel = () => { if (clearTimer) { clearTimeout(clearTimer); clearTimer = null; } $('book-clear').classList.remove('holding'); };
  for (const ev of ['pointerup', 'pointercancel', 'lostpointercapture']) $('book-clear').addEventListener(ev, clearCancel);
  // Sandbox: start from the central inventory
  function loadCentralInventory() {
    const inv = loadInventory();
    for (const [sp, list] of Object.entries(inv.morphs))
      if (Array.isArray(list))
        discoveredBy.set(sp, list.filter(d => d && d.name).map(d => ({ name: d.name, hue: Number.isFinite(d.hue) ? d.hue : NaN, img: d.img || null })));
    if (inv.bullets.length) Bullets.setList(inv.bullets.join(','));
    energy = inv.energy; showEnergy();   // (not a gain)
  }
  // A picture of the ship view as it is now (as drawn this frame, so called from the frame hook,
  // right after the render), scaled down to the world view's scale (C.WS px per cell): a data URL.
  // Whatever the zoom: the SNAP_CELLS-wide square the ship view shows at its middle zoom (6 px per
  // cell), taken from the world view (so unrotated), round the board's edges, at its scale.
  const SNAP_CELLS = C.ZOOM.edge / 6;   // 64 cells
  function snapshot() {
    try {
      const f = C.focus, WS = C.WS, n = Math.round(SNAP_CELLS * WS), W = C.BW * WS, H = C.BH * WS;
      if (!f) return null;
      const x0 = Math.round((f.col + 0.5) * WS - n / 2), y0 = Math.round((f.row + 0.5) * WS - n / 2);
      const cv = document.createElement('canvas');
      cv.width = cv.height = n;
      const g = cv.getContext('2d');
      for (const dx of [-W, 0, W]) for (const dy of [-H, 0, H]) {   // the world view's copies round the torus
        const sx = x0 + dx, sy = y0 + dy;
        const ax = Math.max(0, sx), ay = Math.max(0, sy), bx = Math.min(W, sx + n), by = Math.min(H, sy + n);
        if (bx > ax && by > ay) g.drawImage($('gl'), ax, ay, bx - ax, by - ay, ax - sx, ay - sy, bx - ax, by - ay);
      }
      return cv.toDataURL();
    } catch (e) { return null; }
  }

  // The player's morph, set by [put_player] and a pick in the morph book: what Respawn puts back.
  let playerMorph = null;
  // Put morph <id> in by `how` (center, random, far); follow it if focusNew. False if unknown.
  function putMorph(token, id, how, focusNew) {
    const p = findMorph(id);
    if (!p) { console.warn(`script.txt: unknown morph [${token}:${id}]`); return false; }
    if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [${token}:${id}:${how}], using center`); how = 'center'; }
    C.selectMorph(p.id); C.putPattern(how, focusNew);
    return true;
  }
  // A morph of the current species by its name (limus, gyrans, ...), or a pattern id (5N7KKM-limus).
  const findMorph = id => PATTERNS.find(p => C.speciesOf(p) === C.speciesOf(C.pattern) &&
                                             C.morphName(p).replace(C.BOOM, '') === id) || PATTERNS.find(p => p.id === id);
  // The current species' first morph (its first pattern, as listed), by name
  function firstMorph() {
    const p = C.pattern && PATTERNS.find(p => C.speciesOf(p) === C.speciesOf(C.pattern) && !p.hidden);
    return p ? C.morphName(p).replace(C.BOOM, '') : null;
  }
  // The visible species ([add_allow_species]), in order: { name, id, morphs: its visible morphs }.
  // For the species chooser to come; nothing reads it yet ([species] loads any species).
  const visSpecies = [];
  const speciesBlocks = new Map();   // [on_choose_species[:<id>]] blocks: id ('' = any species) -> its lines
  // A species chosen (the species book, to come; [choose_species:<id>]): load it ([species], with
  // its [add_allow_species] name), make its morphs visible ([allow_morph], its [add_allow_species]
  // list), then run the actions of its [on_choose_species:<id>] block, or else of the general one, in
  // order and at once (no waiting: a block holds actions only).
  function chooseSpecies(id) {
    const e = visSpecies.find(x => x.id === id);
    if (!e) return console.warn(`script.txt: [choose_species:${id}]: not added by [add_allow_species]`);
    ACTIONS.species(e.name, id);
    ACTIONS.allow_morph(e.morphs.join(','));
    for (const L of speciesBlocks.get(id) ?? speciesBlocks.get('') ?? [])
      for (const a of L.actions) ACTIONS[a.name](...a.args);
  }
  const PUT_HOW = ['center', 'random', 'far'];
  const ACTIONS = {
    // add_allow_species:<name>:<id>[:<ids>]: add species <id> (named <name>, as [species]) to the
    // visible species, with its visible morphs <ids> (as [allow_morph]; none given: all). Again:
    // its name and morphs updated.
    // on_choose_species[:<id>] ... end: a block of lines, taken out of the script's run
    // (parseScript), whose actions chooseSpecies runs. No id: for any species; with one: for that
    // species only, instead of the general block.
    on_choose_species: () => {},
    end: () => {},
    // choose_species:<id>: as if the player chose species <id> in the species book (chooseSpecies)
    // (redone by Retry, in its replay and in the current line run afresh: the species the player
    // chose since, if any)
    choose_species: (id) => chooseSpecies(resetting && chosenSpecies ? chosenSpecies : id),
    add_allow_species: (name, id, ids) => {
      if (!id || !PATTERNS.some(p => C.speciesOf(p) === id)) return console.warn(`script.txt: unknown species [add_allow_species:${name}:${id}]`);
      const morphs = ids === undefined ? ['all'] : morphIds([ids]);
      const e = visSpecies.find(x => x.id === id);
      if (e) Object.assign(e, { name: name || e.name, morphs }); else visSpecies.push({ name, id, morphs });
      syncButtons();   // (sandbox: the Species button)
    },
    // species:<name>:<id>: load species <id> (its rule; the world is kept), with no visible
    // morphs; the morph book names it "Species: <name>" (genus and species, e.g. Aquarium
    // originale). ([species:<id>]: no name, the book shows the id.)
    species: (name, id) => {
      if (id === undefined) [name, id] = ['', name];
      const p = PATTERNS.find(p => C.speciesOf(p) === id);
      if (!p) return console.warn(`script.txt: unknown species [species:${name}:${id}]`);
      C.setPattern(p.id);
      speciesName = name;
      vis.all = false; vis.allow.clear();
      if (!discoveredBy.has(id)) discoveredBy.set(id, []);
      discovered = discoveredBy.get(id); cand = null;
      if (book.open && book.mode === 'morph') renderBook();
    },
    // allow_morph:<ids>: the visible morphs are exactly these, in this order ("all": every morph of
    // the species; no ids: none)
    allow_morph: (...args) => {
      const ids = morphIds(args);
      vis.all = ids.includes('all');
      vis.allow = new Set(ids.filter(id => id !== 'all'));
      allDiscovered();
    },
    // morph_book:on|off: show or hide the Morphs button (off also closes the book)
    morph_book: (state = 'on') => setBook(state !== 'off'),
    // clear_world[:<bg>][:empty]: empty the world, with no dots (until [dot:<n>]), no bullets and no
    // [min_morph]; empty: empty space, no player (emptySpace), until the next [clear_world]
    clear_world: (...args) => {
      emptySpace = args.includes('empty');
      let bg = args.find(a => a !== 'empty') || 'black';
      if (!BACKGROUNDS.includes(bg)) { console.warn(`script.txt: unknown background [clear_world:${bg}], using black`); bg = 'black'; }
      background = bg;
      minMorph = null;   // a new scene: no [min_morph] until the script sets one again
      if (Dots.on) { Dots.set(0); Dots.clear(); }   // nor dots (none put back) until a [dot:<n>]; the energy kept
      Bullets.clear();   // nor bullets in flight
      C.clearWorld();
    },
    // put_morph:<id>[:center|random|far]: put morph <id> of the current species into the world:
    // at the centre (default, unrotated), at a random empty spot and angle (the lab's Add), or at
    // a random angle in the empty spot farthest from the followed ship. Never followed: the ship
    // stays followed (if none, none is).
    put_morph: (id, how = 'center') => putMorph('put_morph', id, how, 'never'),
    // put_player:<id>[:center|random|far]: as put_morph, but the new creature is the player's
    // ship (followed), and <id> becomes the player's morph: what Respawn puts back.
    // put_player:first: the species' first morph (its first pattern).
    put_player: (id, how = 'center') => {
      if (id === 'first') id = firstMorph() ?? id;
      if (putMorph('put_player', id, how, true)) playerMorph = id;
    },
    // discover_morph:<ids>: add these morphs to the discovered list (the player may have found
    // others already: they stay) (the condition is [need:discover_morph])
    // Given, not earned: kept out of the central inventory (for testing), until the ship holds it
    // long enough for its picture.
    discover_morph: (...args) => {
      for (const id of morphIds(args)) if (!discovered.some(x => x.name === id)) { discover(id); discovered.find(x => x.name === id).given = true; }
    },
    // morph_book_action:<type>: what a click on a discovered morph in the book does (it also closes
    // the book): clear_put = [clear_world][put_player:<id>]; put = [put_player:<id>:random];
    // none = nothing. The default is clear_put.
    // jump_here: the script starts at this line, after running every action before it, with no
    // messages, conditions or chapter cards (see the start, at the end of this file)
    jump_here: () => {},
    // chapter:<no>:<title>: a title card over the game window, "Chapter <no>" (<no> final: "FINAL
    // CHAPTER") above <title>, faded in; Space or a
    // click anywhere closes it (fading out) and is used up. The game plays on meanwhile, and the
    // line waits until it is closed.
    chapter: (no = '', ...title) => { showChapter(no, title.join(':')); return () => !chapterOpen; },
    // merge:<x>: set the detection's merge saddle to <x> (basins merge across saddles >= x times
    // the lower peak; the lab's Merge saddle slider). No x: back to the species's own value.
    // [species] also resets it to that.
    merge: (x) => {
      const v = x === undefined ? AQUARIUM233.rules[C.pattern.rule]?.detect?.saddle ?? 0.9 : parseFloat(x);
      if (!Number.isFinite(v) || v <= 0) return console.warn(`script.txt: bad value [merge:${x}], ignored`);
      C.DETECT.saddle = v; C.analyze(0);
    },
    // min_morph:<n>[:<ids>[:center|random|far]]: keep at least n creatures besides the player's
    // ship on the board: whenever there are fewer, put in one more (as [put_morph]), a morph
    // picked at random from <ids> (default, or all: the visible list), placed as given (default random).
    // n = 0 or off: stop.
    min_morph: (n, ids = '', how = 'random') => {
      const k = parseInt(n, 10);
      if (n === 'off' || k === 0) { minMorph = null; return; }
      if (!(k > 0)) return console.warn(`script.txt: bad value [min_morph:${n}], ignored`);
      if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [min_morph:${n}:${ids}:${how}], using random`); how = 'random'; }
      minMorph = { n: k, ids: ids === 'all' ? [] : morphIds([ids]), how, next: 0 };
    },
    // dot:<n>|off: the dots channel (dots.js): keep <n> dots on the board, which any creature
    // eats; the player's ship eating one scores 10 (the score shows at the top, from 0 when the
    // dots are turned on); off: no dots, no score.
    dot: (v) => {
      if (v === 'off') return Dots.set('off');
      const n = parseInt(v, 10);
      if (!(n >= 0)) return console.warn(`script.txt: bad value [dot:${v}], ignored`);
      Dots.set(n);
    },
    // gun[:<options>]: the player's gun (bullets.js): C or Fire shoots the pattern in use, which
    // masks (or erases) the creatures it hits while they wear it down. Options, any order:
    // off, erase|mask, rgb|r|g|b|color_now, wear_color|wear_all (defaults mask, color_now, wear_color).
    gun: (...options) => Bullets.setGun(options),
    // bullet:<names>[:<name>]: the patterns the gun can fire: add names (separated by commas),
    // random (one more at random), all, or clear (none); then equip <name> (added too if not
    // listed). With two or more, ◀ ▶ (D / F) switch between them. The patterns it adds are shown:
    // the bullet book opens and reveals them (not when replayed: [jump_here], Reset, I).
    bullet: (names = '', equip) => {
      const before = Bullets.list;
      Bullets.setList(names, equip);
      const added = Bullets.list.filter(n => !before.includes(n));
      // all: for testing, kept out of the central inventory (patterns named, or random, are saved)
      if (names.split(',').some(t => t.trim().toLowerCase() === 'all')) {
        const named = names.split(',').map(t => t.trim()).concat(equip ? [equip.trim()] : []).map(Bullets.kindOf);
        for (const n of added) if (!named.includes(n)) unsavedBullets.add(n);
      }
      payBullets(added.filter(n => !unsavedBullets.has(n)).length);   // each costs [new_bullet]'s price, as a created one ([bullet:all]'s extras free)
      if (!replaying && added.length && Bullets.on && !book.open && !chapterShowing()) showNewBullets(added);
    },
    // bullet_limit:<x>: the bullets channel is emptied once it holds more than x bullets' worth of
    // mass (bullets run into each other and grow); default 8
    bullet_limit: (x = '8') => {
      const v = parseFloat(x);
      if (!(v > 0)) { console.warn(`script.txt: bad value [bullet_limit:${x}], ignored`); return; }
      Bullets.setLimit(v);
    },
    // new_bullet:<n>|off: a bullet pattern may be created for n energy (the score): while the
    // energy is at least n, the Bullets button reads "Create" and a "?" in the bullet book creates
    // that pattern, taking n from the energy (bulletCredits); off: none may
    new_bullet: (n) => {
      const v = parseFloat(n);
      newBullet = n === 'off' ? null : v > 0 ? { every: v } : null;
      if (n !== 'off' && !(v > 0)) console.warn(`script.txt: bad value [new_bullet:${n}], ignored`);
      syncNew();
    },
    // new_scene: fade the world view and the ship view's picture to black (WORLD_FADE_MS; the
    // line waits for it), so the actions after it set the scene unseen; the next message fades
    // them back in as it shows. The game plays on.
    new_scene: () => { fadeWorld(1); sceneHidden = true; return fadeDone; },
    // next_msg[:<tip>[:<next>[:<act>]]]: the grey tip shown once on the next message that moves on
    // by Space / a click (default "Press space or click message to continue"), the mark ending
    // such messages (default ≫), and the mark ending a message waiting on a condition (default
    // ⧗). A part left out or empty is back to its default.
    next_msg: (tip = '', next = '', act = '') => {
      marks.tip = tip || MARKS.tip; marks.next = next || MARKS.next; marks.act = act || MARKS.act;
      tipPending = !replaying;   // not when skipped to ([jump_here], dev I): the tip is for the script's first message
    },
    // speed:<x>: the game speed (the Speed slider), x = n or 1/n (1/8 .. 8), snapped to the
    // nearest slot
    speed: (x = '1') => {
      const [a, b] = String(x).split('/'), v = b !== undefined ? +a / +b : +a;
      if (!(v > 0)) return console.warn(`script.txt: bad value [speed:${x}], ignored`);
      C.setSpeed(C.SPEEDS.reduce((best, s, i) => Math.abs(Math.log(s / v)) < Math.abs(Math.log(C.SPEEDS[best] / v)) ? i : best, 0));
    },
    // upward:on|off: the ship view's Upward (heading up) on or off, as the Upward button / E
    upward: (v = 'on') => {
      if (v !== 'on' && v !== 'off') return console.warn(`script.txt: bad value [upward:${v}], ignored`);
      if ((C.viewMode !== 'board') !== (v === 'on')) C.toggleHeadUp();
    },
    // stars:<angle>:<speed>:<bright>: the star field drifts at <angle> degrees, clockwise from left
    // (0: to the left, the default; 90: up), at <speed> and <bright> times the default (1); a part
    // left out or empty is the default
    stars: (angle = '', speed = '', bright = '') => {
      const num = (x, d) => x === '' ? d : parseFloat(x);
      const [a, s, b] = [num(angle, 0), num(speed, 1), num(bright, 1)];
      if (![a, s, b].every(Number.isFinite) || s < 0 || b < 0)
        return console.warn(`script.txt: bad value [stars:${angle}:${speed}:${bright}], ignored`);
      C.setStars(a, s, b);
    },
    // music:<list>[:<play>]: the background playlist becomes <list> (track numbers separated by
    // commas, 1 = the first of page.js TRACKS), played in turn and looped; <play>: switch to that
    // track now; none: the track playing carries on, and the list follows it
    music: (listArg = '', playArg) => {
      const ok = n => { const i = parseInt(n, 10) - 1; return i >= 0 && i < Page.TRACKS.length ? i : -1; };
      const ns = listArg.split(',').map(s => s.trim()).filter(Boolean);
      const list = ns.map(ok).filter((i, k) => {
        if (i < 0) console.warn(`script.txt: no such track [music:...${ns[k]}...] (1-${Page.TRACKS.length}), ignored`);
        return i >= 0;
      });
      let play;
      if (playArg !== undefined && playArg.trim()) {
        play = ok(playArg);
        if (play < 0) { console.warn(`script.txt: no such track to play [music:${listArg}:${playArg}], ignored`); play = undefined; }
      }
      if (list.length) Page.setPlaylist(list, play);
      else if (play !== undefined) Page.setPlaylist([play], play);   // no list: that track alone
    },
    // clone:on|off: allow clone mode or not (see startClone): while allowed, the morph book offers the
    // discovered morphs to clone into the world while the ship is alive, as often as wanted, until
    // [clone:off]
    clone: (v = 'on') => { if (v === 'off') { endClone(true); syncButtons(); } else { clone.on = true; syncButtons(); if (book.open) renderBook(); } },
    // erase:on|off: allow erase ("recycle", the selected R G B patch again) or not (back to masking)
    erase: (v = 'on') => { eraseOn = v !== 'off'; if (!eraseOn) C.action = 0; },
    // blank: the line is an empty message: it shows nothing but its mark (≫, or ⧗ with a
    // condition) and waits for Space / a click like any message (see parseScript: L.msg)
    blank: () => {},
    // retry_action:<type>: what Respawn (and Z) does: clear_put, put or none as morph_book_action,
    // but with the player's morph; or book (open the morph book). No type: back to the default
    // (clear_put).
    retry_action: (type) => {
      if (type !== undefined && !['clear_put', 'put', 'none', 'book'].includes(type)) {
        console.warn(`script.txt: unknown [retry_action:${type}]`); return;
      }
      retryAction = type ?? null;
    },
    morph_book_action: (type = 'clear_put') => {
      if (!['clear_put', 'put', 'none'].includes(type)) { console.warn(`script.txt: unknown [morph_book_action:${type}]`); type = 'none'; }
      book.action = type; renderBook();
    },
  };
  const CONDITIONS = {
    // name: (...args) => test, made when its line is reached; test() is true once the event has
    // happened (called every frame while the line waits).
    // first_touch: the player pressed on the ship, and its detected morph is now not limus (any
    // other, or unknown); the press is needed, since a fresh ship can read as unknown at first
    // discover_morph[:<ids>[:any|all]]: one of <ids> (any, the default) or every one of them
    // (all) is in the morph book, found earlier or now (then once the book its discovery opened is
    // closed); no ids: the player discovers a new morph (any, after the line is reached) and
    // closes the morph book it opened
    // first_clone: the player has placed a clone (picked in the morph book while clone mode is
    // allowed, then pasted on the ship view) since the line was reached
    first_clone: () => { const n0 = clonesPlaced; return () => clonesPlaced > n0; },
    // respawn: the ship has dissolved and the player has put it back, by Respawn (Z) or Choose (a
    // pick in the morph book), since the line was reached
    respawn: () => { const n0 = respawns; return () => respawns > n0; },
    // create_bullet: the player has created a bullet type (Create, then a "?" in the bullet book)
    // since the line was reached ([new_bullet] must allow one)
    create_bullet: () => { const n0 = bulletsCreated; return () => bulletsCreated > n0; },
    // recycle: the player has erased some mass (erase on, any channel, pressed on the ship view)
    // since the line was reached, and the ship has then dissolved
    recycle: () => {
      let erased = false;
      return () => {
        if (C.masking && C.action === 1) erased = true;
        return erased && shipLost();
      };
    },
    // object:<n>[:exact|min|max]: the world holds exactly (the default) / at least / at most n
    // objects (creatures detected), the ship included; the count smoothed (smoothed).
    // Progress: "morph: <count>/<n>".
    object: (n, how = 'exact') => {
      const v = parseInt(n, 10);
      if (!(v >= 0)) { console.warn(`script.txt: bad value [need:object:${n}], ignored`); return () => true; }
      const CMP = { exact: k => k === v, min: k => k >= v, max: k => k <= v };
      if (!CMP[how]) { console.warn(`script.txt: unknown option [need:object:${n}:${how}], using exact`); how = 'exact'; }
      const count = smoothed(allObjects);
      const t = () => { const k = count(); return k !== null && CMP[how](k); };
      t.progress = () => `morph: ${count() ?? 0}/${v}`;
      return t;
    },
    // mass:<r|g|b|rgb>:<x>[:min|max]: the objects in the world (the ship included) hold at least
    // (the default) / at most x g of mass in that channel (rgb: all three), as the stats count it;
    // smoothed as [need:object]. Progress: "mass: <mass>/<x>", one decimal.
    mass: (ch, x, how = 'min') => {
      const CH = { r: [0], g: [1], b: [2], rgb: [0, 1, 2] }, v = parseFloat(x);
      if (!CH[ch] || !(v >= 0)) { console.warn(`script.txt: bad [need:mass:${ch}:${x}], ignored`); return () => true; }
      const CMP = { min: m => m >= v, max: m => m <= v };
      if (!CMP[how]) { console.warn(`script.txt: unknown option [need:mass:${ch}:${x}:${how}], using min`); how = 'min'; }
      const mass = smoothed(() => { const m = C.worldMass(); return m && CH[ch].reduce((s, c) => s + m[c], 0); });
      const t = () => { const m = mass(); return m !== null && CMP[how](m); };
      t.progress = () => `mass: ${(mass() ?? 0).toFixed(1)}/${v}`;
      return t;
    },
    // survive_time:<s>: the ship has stayed alive for s game seconds (C.STEPS_PER_S steps each:
    // paused time does not count, and the same against the creatures on any device or speed), from
    // when the line is reached, not under a chapter card; dissolving stops the count, and the next
    // ship (Respawn) starts it again from 0. Progress: "time: <t>/<s>".
    survive_time: (s) => {
      const v = parseFloat(s);
      if (!(v >= 0)) { console.warn(`script.txt: bad value [need:survive_time:${s}], ignored`); return () => true; }
      let n = 0, last = C.simSteps, lost = false;
      const t = () => {
        const dn = C.simSteps - last;
        last = C.simSteps;
        if (shipLost()) lost = true;   // the count stops (still shown) ...
        else if (lost) { lost = false; n = 0; }   // ... and starts again from 0 once a ship is back (Respawn)
        else if (!book.open && !chapterShowing()) n += dn;
        return n >= v * C.STEPS_PER_S;
      };
      t.progress = () => `time: ${Math.min(v, Math.floor(n / C.STEPS_PER_S))}/${v}`;
      return t;
    },
    // upward: the ship view's Upward is on
    upward: () => () => C.viewMode !== 'board',
    // energy:<n>: the energy (the score from eating dots, [dot]) has reached at least n
    energy: (n) => {
      const v = parseFloat(n);
      if (!Number.isFinite(v)) { console.warn(`script.txt: bad value [need:energy:${n}], ignored`); return () => true; }
      return () => energy >= v;
    },
    discover_morph: (...args) => {
      let mode = 'any';
      if (args.length > 1 && ['any', 'all'].includes(args[args.length - 1])) mode = args.pop();
      const ids = morphIds(args), c0 = book.found.length;
      if (!ids.length) return () => book.found.length > c0;
      const has = id => discovered.some(x => x.name === id);
      return () => !book.open && (mode === 'all' ? ids.every(has) : ids.some(has));
    },
    first_touch: () => {
      let touched = false;
      return () => {
        if (C.masking) touched = true;
        const d = C.stats.short.det;
        return touched && !!C.focus && !!d && d.name !== 'limus';
      };
    },
  };

  // the objects in the world, the player's ship included (null: not detected yet)
  const allObjects = () => C.found ? C.found.comps.length : null;
  // A value smoothed (the object count, the world's mass): the detection flickers (a creature
  // splits or merges for a step or two), so the value is the median of fn()'s values over the last
  // OBJECT_SMOOTH_STEPS steps, one sample per step (the latest value at that step: it may change
  // while paused, e.g. a clone pasted) (the same value for the test and the progress shown).
  // fn() null: not yet.
  const OBJECT_SMOOTH_STEPS = C.STEPS_PER_S;   // a game second
  function smoothed(fn) {
    const samples = [];   // [step, value]
    let med = null;
    return () => {
      const now = C.simSteps, k = fn();
      if (k === null) return med;
      if (samples.length && samples[samples.length - 1][0] === now) samples[samples.length - 1][1] = k;
      else samples.push([now, k]);
      while (samples.length > 1 && now - samples[0][0] >= OBJECT_SMOOTH_STEPS) samples.shift();
      const sorted = samples.map(x => x[1]).sort((a, b) => a - b);
      med = sorted[Math.floor(sorted.length / 2)];
      return med;
    };
  }

  let SCRIPT = [], line = -1, tests = [];
  // the waiting message's progress (after ⧗), from its conditions' progress(), updated every frame
  let progress = null;   // { el, fns }
  function showProgress() {
    if (progress) progress.el.textContent = progress.fns.map(f => f()).join(' ');
  }
  // Then the [on_choose_species[:<id>]] ... [end] blocks are taken out of the script's run (into
  // speciesBlocks, see chooseSpecies).
  function parseScript(txt) {
    const lines = txt.split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#')).map(l => {
      const toks = [];
      let m;
      while ((m = l.match(/\s*\[([^\]]*)\]\s*$/))) { toks.unshift(m[1].trim()); l = l.slice(0, m.index); }
      const actions = [], conds = [];
      for (const t of toks) {
        const [name, ...args] = t.split(/:(?![^{}]*\})/).map(x => x.trim());   // a ':' inside {...} stays
        if (name === 'need') {
          const [cond, ...cargs] = args;
          if (CONDITIONS[cond]) conds.push({ name: cond, args: cargs });
          else console.warn(`script.txt: unknown condition [${t}], ignored`);
        }
        else if (ACTIONS[name]) actions.push({ name, args });
        else console.warn(`script.txt: unknown action [${t}], ignored`);
      }
      // msg: the line is a message (it has text, or [blank]: an empty one)
      return { text: l.trim(), actions, conds, msg: !!l.trim() || actions.some(a => a.name === 'blank') };
    });
    const out = [];
    let block = null;   // the lines of the block being read
    for (const L of lines) {
      const open = L.actions.find(a => a.name === 'on_choose_species');
      if (open) {
        if (block) console.warn('script.txt: [on_choose_species] inside another block, ignored');
        else { block = []; speciesBlocks.set(open.args[0] || '', block); }
        continue;
      }
      if (L.actions.some(a => a.name === 'end')) {
        if (!block) console.warn('script.txt: [end] with no [on_choose_species], ignored');
        block = null;
        continue;
      }
      if (!block) { out.push(L); continue; }
      if (L.text || L.conds.length) console.warn(`script.txt: a message or condition in an [on_choose_species] block, ignored: ${L.text}`);
      block.push(L);
    }
    if (block) console.warn('script.txt: [on_choose_species] with no [end]: the block runs to the end of the file');
    return out;
  }
  // A line with [new_scene] or [chapter] waits while the morph book shows (and zaps out): the
  // views must not fade, nor the card come up, behind it. checkScript runs it once it is gone.
  let deferred = -1;
  const bookShowing = () => !$('book').hidden;
  // Reset: rebuild the scene (e.g. a world grown too crowded) by redoing, behind a fade, every
  // action from the last checkpoint ([chapter] or [new_scene]) up to the current line, then the
  // current line itself (its message and conditions afresh). Not the checkpoints themselves (no
  // card, no extra fade), nor what only steers the script or the music.
  const CHECKPOINTS = ['chapter', 'new_scene'], REPLAY_SKIP = new Set([...CHECKPOINTS, 'jump_here', 'next_msg', 'music']);
  let resetting = false;
  let replaying = false;   // actions redone (Reset, [jump_here], dev I): nothing pops up ([bullet]'s book)
  function resetScene() {
    if (line < 0 || resetting || chapterShowing() || bookShowing()) return;
    if (clone.held) endClone();
    const cur = deferred >= 0 ? deferred : line;
    let cp = cur;
    while (cp > 0 && !SCRIPT[cp].actions.some(a => CHECKPOINTS.includes(a.name))) cp--;
    resetting = true;
    fadeWorld(1);
    setTimeout(() => {
      Dots.clear(); Bullets.clear();   // the dots put in afresh, no bullets in flight
      replaying = true;
      for (let i = cp; i < cur; i++)
        for (const a of SCRIPT[i].actions) if (!REPLAY_SKIP.has(a.name)) ACTIONS[a.name](...a.args);
      replaying = false;
      sceneHidden = false;
      run(cur);
      fadeWorld(0);
      resetting = false;
    }, WORLD_FADE_MS);
  }
  function run(i) {
    const L = SCRIPT[i];
    if (bookShowing() && L.actions.some(a => a.name === 'new_scene' || a.name === 'chapter')) {
      line = i; deferred = i; tests = []; return;
    }
    deferred = -1;
    line = i;
    progress = null;   // (a line of tokens only keeps the message, not its progress)
    // an action may return a test too: the line also waits for it (chapter: until the title closes)
    const waits = L.actions.map(a => ACTIONS[a.name](...a.args)).filter(f => typeof f === 'function');
    tests = L.conds.map(c => CONDITIONS[c.name](...c.args)).concat(waits);
    if (L.msg && sceneHidden) { sceneHidden = false; fadeWorld(0); }   // the scene is set: show it
    if (L.msg) {
      // conditions already met on arrival (e.g. the energy is there already): shown in the waiting
      // colour, but it moves on by Space / a click, as a message with no condition
      const met = tests.length > 0 && tests.map(t => t()).every(Boolean);
      const wait = tests.length > 0;
      if (met) { tests = []; if (L.conds.length) showSuccess(); }   // "SUCCESS" too
      const last = i + 1 >= SCRIPT.length, tip = !tests.length && !last && tipPending;
      if (tip) tipPending = false;
      const fns = tests.map(t => t.progress).filter(Boolean);
      progress = fns.length ? { el: document.createElement('span'), fns } : null;
      showProgress();
      say(L.text, tests.length ? marks.act : last ? '' : tip ? `${marks.tip} ${marks.next}` : marks.next,
          wait, !tests.length && !last, progress?.el);   // cyan: the player has (or had) something to do
    }
    else if (!tests.length) next();   // tokens only, nothing to wait for
  }
  function next() { if (line + 1 < SCRIPT.length) run(line + 1); }
  // Space / a click on the message bar moves on: a message with no condition
  const byContinue = () => deferred < 0 && line >= 0 && SCRIPT[line].msg && !tests.length;
  // ?dev=1: O completes the waiting line's [need] conditions at once (a line waiting only on an
  // action, a chapter card, is left alone); P scores as if the ship ate a dot (+10 energy); I goes
  // to the next chapter
  function devSkip() { if (DEV && deferred < 0 && line >= 0 && tests.length && SCRIPT[line].conds.length) { showSuccess(); next(); } }
  // I: on to the next [chapter] line. Every action on the lines in between runs first, in order
  // (as [jump_here] does), with their messages, conditions and chapter cards skipped.
  function devChapter() {
    if (!DEV || line < 0 || resetting) return;
    const cur = deferred >= 0 ? deferred : line + 1;   // a deferred line has not run its actions yet
    const to = SCRIPT.findIndex((L, i) => i > line && L.actions.some(a => a.name === 'chapter'));
    if (to < 0) return;
    if (clone.held) endClone();
    if (chapterOpen) closeChapter();
    replaying = true;
    for (let i = cur; i < to; i++)
      for (const a of SCRIPT[i].actions) if (a.name !== 'chapter') ACTIONS[a.name](...a.args);
    replaying = false;
    run(to);
  }
  // its buttons, at the top left of the world view: Ch (I), Ok (O), Do (P)
  $('dev').hidden = !DEV;
  $('dev-ch').onclick = devChapter; $('dev-ok').onclick = devSkip; $('dev-do').onclick = () => { if (DEV) Dots.eat(); };
  // the frame rate and the simulation steps per second, over the last half second
  $('dev-fps').hidden = !DEV;
  const fps = { t: 0, frames: 0, steps: 0 };
  function devFps() {
    if (!DEV) return;
    const now = performance.now();
    if (!fps.t) { Object.assign(fps, { t: now, frames: 0, steps: C.simSteps }); return; }
    fps.frames++;
    if (now - fps.t < 500) return;
    const s = (now - fps.t) / 1000;
    const parts = [`${Math.round(fps.frames / s)} fps`, `${Math.round((C.simSteps - fps.steps) / s)} steps/s`, `${C.stepMs.toFixed(1)} ms/step`,
                   ...(FFTSim.on ? ['FFT'] : []), `${innerWidth}×${innerHeight}`, `@${devicePixelRatio.toFixed(2)}`];   // (the CSS viewport and pixel ratio: for the phone test's sizes)
    $('dev-fps').textContent = parts.join(stacked ? '\n' : ' · ');   // stacked: one per row
    Object.assign(fps, { t: now, frames: 0, steps: C.simSteps });
  }
  // "SUCCESS" over the game window when a waiting line's [need] conditions are met: fades in, holds,
  // then slowly fades out (the .success animation, SUCCESS_MS); the game goes on meanwhile
  const SUCCESS_MS = 3000;
  let successTimer = null;
  function showSuccess() {
    const el = $('success');
    clearTimeout(successTimer);
    el.hidden = false; el.classList.remove('show');
    void el.offsetWidth;   // restart the animation
    el.classList.add('show');
    successTimer = setTimeout(() => { el.hidden = true; el.classList.remove('show'); }, SUCCESS_MS);
  }
  function checkScript() {
    if (deferred >= 0) { if (!bookShowing()) run(deferred); return; }
    if (line < 0 || !tests.length) return;
    const ok = tests.map(t => t());   // every test, every frame (one may be watching for an action)
    showProgress();
    if (ok.every(Boolean)) { if (SCRIPT[line].conds.length) showSuccess(); next(); }
  }
  // a click on the message bar (a hand cursor over it while that moves on)
  // At the start of story mode the message bar blinks slowly (style.css .msg.attract), until the player first
  // moves a message on
  const playerNext = () => { if (!byContinue()) return; $('msg').classList.remove('attract'); next(); };
  $('msg').addEventListener('pointerdown', playerNext);

  // ---- background ([clear_world:<bg>]) --------------------------------------------------------
  // black (the default), or radiation: a dark orange behind the board, slowly darkening and
  // brightening again (a period of BG_PERIOD_S)
  const BACKGROUNDS = ['black', 'radiation'], BG_PERIOD_S = 8;
  let background = 'black';
  function updateBackground() {
    if (background !== 'radiation') return SimGL.setBackground([0, 0, 0]);
    const k = 0.45 + 0.15 * Math.sin(2 * Math.PI * performance.now() / 1000 / BG_PERIOD_S);   // 0.3 .. 0.6 of the full orange
    SimGL.setBackground([0.40 * k, 0.15 * k, 0.0]);
  }

  // ---- scene fades ([new_scene], [chapter]) -----------------------------------------------------
  const WORLD_FADE_MS = 500;
  let worldFade = { from: 0, to: 0, t0: 0, ms: WORLD_FADE_MS };
  function worldShade() {
    const f = worldFade.ms > 0 ? Math.min(1, (performance.now() - worldFade.t0) / worldFade.ms) : 1;
    return worldFade.from + (worldFade.to - worldFade.from) * f;
  }
  const fadeDone = () => performance.now() - worldFade.t0 >= worldFade.ms;
  function fadeWorld(to, ms = WORLD_FADE_MS) { worldFade = { from: worldShade(), to, t0: performance.now(), ms }; }
  let sceneHidden = false;   // faded out by [new_scene]: the next message fades the views back in
  C.opt.worldShade = worldShade;

  // ---- [min_morph] ------------------------------------------------------------------------
  // Each frame: fewer creatures than minMorph.n besides the ship -> put one more in. At most one
  // every MIN_MORPH_STEPS (half a game second), so a creature just put in is counted (detected)
  // before the next. (In steps: none is put in while paused.)
  const MIN_MORPH_STEPS = Math.round(C.STEPS_PER_S / 2);
  let minMorph = null;   // { n, ids, how, next: the step (C.simSteps) of the next allowed put }
  function checkMinMorph() {
    if (!minMorph || !C.found || C.simSteps < minMorph.next) return;
    const others = C.found.comps.length - (C.focusIdx >= 0 ? 1 : 0);
    if (others >= minMorph.n) return;
    const pool = (minMorph.ids.length ? minMorph.ids : visibleList()).filter(id => findMorph(id));
    if (!pool.length) return;
    putMorph('min_morph', pool[Math.floor(Math.random() * pool.length)], minMorph.how, 'never');   // never the new focus
    minMorph.next = C.simSteps + MIN_MORPH_STEPS;
  }

  // ---- chapter title card ([chapter]) ----------------------------------------------------------
  // The views fade out first (WORLD_FADE_MS), then the card fades in; closing it fades the card
  // out and the views back in together.
  const FADE_MS = 1500;   // as .chapter's transition in style.css
  let chapterOpen = false, cardShown = false, chapterTimer = null;
  const chapterShowing = () => chapterOpen || !$('chapter').hidden;
  function showChapter(no, title) {
    $('chapter-no').textContent = no.toLowerCase() === 'final' ? 'FINAL CHAPTER' : `Chapter ${no}`; $('chapter-title').textContent = title;
    say('');   // a new chapter: the message bar starts empty
    clearTimeout(chapterTimer);
    chapterOpen = true; cardShown = false;
    fadeWorld(1);
    chapterTimer = setTimeout(() => {
      const el = $('chapter');
      el.hidden = false;
      void el.offsetWidth;   // start the fade from opacity 0
      el.classList.add('show');
      cardShown = true;
    }, WORLD_FADE_MS);
  }
  function closeChapter() {
    chapterOpen = false; cardShown = false;
    const el = $('chapter');
    el.classList.remove('show');
    if (!sceneHidden) fadeWorld(0, FADE_MS);   // the views come back as the card goes
    clearTimeout(chapterTimer);
    chapterTimer = setTimeout(() => { el.hidden = true; }, FADE_MS);
  }
  // Space or a click anywhere closes it, before anything else sees the event (so it does nothing
  // else). Other keys pass through to the game as usual.
  for (const ev of ['keydown', 'pointerdown'])
    window.addEventListener(ev, e => {
      if (!cardShown || (ev === 'keydown' && e.key !== ' ')) return;   // (not before the card shows)
      e.preventDefault(); e.stopImmediatePropagation();
      closeChapter();
    }, { capture: true });

  // ---- key badges ---------------------------------------------------------------------------
  // Holding Shift or Alt shows each control's key in a small badge at its top right corner (from
  // its data-key; a slider gets one at each end, from data-key-lo / data-key-hi, over the track),
  // and the channel patches' 1-3 (core opt.showKeys). Space is left out (it moves the messages on).
  const KEYS = { play: '\u21b5', sound: 'M', headup: 'E', morphs: 'X', reset: 'Z', 'morphs-dead': 'X', fire: 'A',
                 'book-close': 'Esc', 'bullet-prev': 'D', 'bullet-next': 'F', bullets: 'S', species: 'C', restart: '\u2190',
                 'dev-ch': 'I', 'dev-ok': 'O', 'dev-do': 'P' };
  const SLIDER_KEYS = { speed: ['Q', 'W'], zoom: ['R', 'T'] };   // [left end (less), right end (more)]
  const SLIDER_INSET = 10;   // px in from each end of the track
  // (a button with a .kt text span gets its badge there, at the text's top right)
  for (const [id, k] of Object.entries(KEYS)) ($(id).querySelector('.kt') || $(id)).dataset.key = k;
  let keysOn = false;
  C.opt.showKeys = () => keysOn;
  function showKeys(on) {
    if (on === keysOn) return;
    keysOn = on;
    document.documentElement.classList.toggle('keys', on);
    document.querySelectorAll('.kbd-slider').forEach(b => b.remove());
    if (!on) return;
    for (const el of document.querySelectorAll('[data-key]'))
      if (getComputedStyle(el).position === 'static') el.style.position = 'relative';   // the badge's anchor
    for (const [id, [lo, hi]] of Object.entries(SLIDER_KEYS)) {
      const input = $(id), box = input.parentElement;
      if (!input.offsetWidth) continue;   // hidden
      if (getComputedStyle(box).position === 'static') box.style.position = 'relative';
      for (const [k, end] of [[lo, 0], [hi, 1]]) {
        const b = document.createElement('span');
        b.className = 'kbd kbd-slider'; b.textContent = k;
        b.style.left = input.offsetLeft + (end ? input.offsetWidth - SLIDER_INSET : SLIDER_INSET) + 'px';
        b.style.top = input.offsetTop + 3 + 'px';
        box.append(b);
      }
    }
  }
  document.addEventListener('keydown', e => { if (e.key === 'Shift' || e.key === 'Alt') { if (e.key === 'Alt') e.preventDefault(); showKeys(true); } });
  document.addEventListener('keyup', e => { if (e.key === 'Shift' || e.key === 'Alt') showKeys(e.shiftKey || e.altKey); });
  window.addEventListener('blur', () => showKeys(false));

  // ---- core hooks -------------------------------------------------------------------------
  C.hooks.change = () => syncButtons();
  // the score, at the top of the world view (centred) while the dots are on; in fullscreen too
  Dots.hooks.score = () => {
    const d = Dots.score - dotsSeen;   // eaten since (not a fall: dots.js may restart its count)
    dotsSeen = Dots.score;
    if (d > 0) changeEnergy(d); else showEnergy();
  };
  function showEnergy() {
    $('score').hidden = !Dots.on; $('score-v').textContent = energy;
    const was = bulletNew;
    syncNew();   // [new_bullet:<n>]: one more may be creatable
    if (bulletNew !== was && book.open && book.mode === 'bullet') renderBook();   // "Create or choose bullet" now, or no more
  };
  let newBullet = null;   // [new_bullet:<n>]: { every: n }
  // in the dots' colour (over the world view, dark on either page)
  Dots.hooks.color = (hue, rgb) => {
    $('score').style.color = `rgb(${rgb.map(v => Math.round(v * 255)).join(',')})`;
  };
  // Fire (and X): dimmed while it can't fire (no ship, or just fired)
  function doFire() { if (Bullets.fire()) flash('fire'); }
  $('fire').addEventListener('pointerdown', e => { e.preventDefault(); doFire(); });
  $('bullet-prev').addEventListener('pointerdown', e => { e.preventDefault(); if (Bullets.cycle(-1)) flash('bullet-prev'); });
  $('bullet-next').addEventListener('pointerdown', e => { e.preventDefault(); if (Bullets.cycle(1)) flash('bullet-next'); });
  // the pattern's name under "Fire"; < > above it while the list holds two or more (D / F)
  Bullets.hooks.change = () => {
    const many = Bullets.on && Bullets.list.length > 1;
    $('fire').hidden = !Bullets.on;
    $('fire-id').hidden = !Bullets.kind;
    $('fire-id').textContent = (Bullets.kind || '').replace(' ', '\n');   // two words: two rows
    $('bullet-nav').hidden = !many;
    syncButtons();   // the Bullets button
    syncNew();   // the list changed: as many more creatable?
    saveInventory();   // (a pattern added: into the central inventory)
    if (book.open && book.mode === 'bullet') { if (Bullets.on) renderBook(); else closeBook(); }
    $('fire').title = 'Fire a bullet the way the ship moves (A)' + (many ? '; D / F: the previous / next pattern' : '');
  };
  C.hooks.frame = () => {
    devFps();
    updateBackground();
    // Respawn / Choose fade with the views ([new_scene], [chapter]), as "Dissolved" (core)
    const shade = worldShade(), dead = document.querySelector('.ship-dead');
    dead.style.opacity = 1 - shade;
    dead.style.visibility = shade > 0.99 ? 'hidden' : '';
    $('msg').classList.toggle('next', byContinue());
    if (Bullets.on) {   // Fire: its ring in the bullets' colour
      $('fire').classList.toggle('ready', Bullets.canFire());
      const rgb = `rgb(${Bullets.color.map(v => Math.round(v * 255)).join(',')})`;
      $('fire').style.setProperty('--fire', rgb); $('bullet-nav').style.setProperty('--fire', rgb);   // (< > too)
    }
    Dots.frame(); Bullets.frame(); checkMinMorph(); checkDiscovery(); checkScript(); if ($('reset').hidden === lostShown()) syncButtons(); };

  // ---- controls ---------------------------------------------------------------------------
  function syncButtons() {
    $('play').setAttribute('aria-pressed', !C.running);   // "Pause", highlighted while paused
    $('reset').hidden = !lostShown();
    $('morphs-dead').hidden = !lostShown() || !book.on;
    $('morphs').hidden = !book.on || shipLost();   // while the ship is lost, Choose takes its place
    if (clone.held && shipLost()) clone.held = null;   // lost while holding: the clone goes
    $('clone-tools').hidden = !clone.held;
    $('morphs').classList.toggle('clone-on', clone.on);   // pink, and "Clone", while clone mode is allowed
    $('morphs').textContent = clone.on ? 'Clone' : 'Morphs';
    $('morphs').title = clone.on ? 'Clone a discovered morph into the world: the morph book (X)'
                                 : 'The morph book: the morphs you have discovered (X)';
    $('clone-ok').hidden = !(clone.held && clone.held.at);   // ✓ once a touch has placed it
    $('bullets').classList.toggle('new', bulletNew);
    $('bullets').textContent = bulletNew ? 'Create' : 'Bullets';
    $('bullets').hidden = !Bullets.on || shipLost();   // the bullet book's button: while the gun is on, not while lost
    $('species').hidden = MODE !== 'sandbox' || !visSpecies.length;   // sandbox mode: the species book's button
    $('headup').setAttribute('aria-pressed', C.viewMode !== 'board');
    $('speed').value = C.speedIdx + 1;
    $('zoom').value = C.zoomIdx + 1;
    for (const el of document.querySelectorAll('input[type=range]')) fillRange(el);
  }

  // Play while holding a clone (paused): leaves clone mode (the clone is dropped) and plays
  const doPlay = () => {
    flash('play'); C.togglePlay();
    if (clone.held && C.running) endClone();
  };
  // Respawn (and Z), by [retry_action:<type>]: clear_put (clear the world, put the player's morph at
  // the centre; the default), put (put it at a random spot), none (nothing), or book (open the
  // morph book to pick the morph, by its morph_book_action). The Morphs button under Respawn
  // (with the book on) opens the book, where a pick then does the morph_book_action. Either put makes the new creature the player's ship and plays. The player's morph:
  // the last [put_player] or book pick (none yet: the last morph put in).
  let retryAction = null;
  function doReset() {
    const type = retryAction ?? 'clear_put';
    const name = playerMorph ?? C.morphName(C.pattern).replace(C.BOOM, '');
    if (type === 'book') openBook();
    else if (type === 'clear_put') { C.clearWorld(); ACTIONS.put_player(name, 'center'); C.setRunning(true); }
    else if (type === 'put') { ACTIONS.put_player(name, 'random'); C.setRunning(true); }
    if (type === 'clear_put' || type === 'put') { retake = true; cand = null; respawns++; }   // its picture may be retaken; [need:respawn]
  }
  const doView = () => { flash('headup'); C.toggleHeadUp(); };
  const doSpeed = d => { flash('speed'); C.changeSpeed(d); };
  // The game's zooms: the core's but the outermost two (2 and 16 px per cell): slots 1 .. 5
  const ZOOM_MIN = 1, ZOOM_MAX = C.ZOOMS.length - 2;
  const clampZoom = i => Math.max(ZOOM_MIN, Math.min(ZOOM_MAX, i));
  if (C.zoomIdx !== clampZoom(C.zoomIdx)) C.setZoom(clampZoom(C.zoomIdx));   // e.g. ?zoom=2
  const doZoom = d => { flash('zoom'); C.setZoom(clampZoom(C.zoomIdx + d)); };   // R out, T in
  $('play').onclick = doPlay;
  const doRestart = () => { flash('restart'); resetScene(); };   // Retry (and Backspace)
  $('restart').onclick = doRestart;
  $('reset').onclick = doReset;
  $('morphs-dead').onclick = () => openBook();
  $('headup').onclick = doView;
  $('speed').oninput = e => C.setSpeed(+e.target.value - 1);
  $('zoom').value = C.zoomIdx + 1;
  $('zoom').oninput = e => C.setZoom(clampZoom(+e.target.value - 1));

  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (book.open) {   // the book takes the keys while open; Space does nothing
      // a slot's letter (SLOT_KEYS) picks it, if it can be picked (which closes the book too)
      const el = k.length === 1 && $('book-grid').querySelector(`.book-slot.pick[data-key="${k.toUpperCase()}"]`);
      if (el) el.click();
      // else Escape (only: the letters that open the books pick slots in them) closes it
      else if (k === 'Escape') closeBook();
      else if (e.code === 'Space') e.preventDefault();   // not even the focused Close button
      return;
    }
    if (e.code === 'Space') { e.preventDefault(); playerNext(); }
    else if (k === 'o') { if (DEV) flash('dev-ok'); devSkip(); }
    else if (k === 'p') { if (DEV) { flash('dev-do'); Dots.eat(); } }
    else if (k === 'i') { if (DEV) flash('dev-ch'); devChapter(); }
    else if (k === 'Enter') { e.preventDefault(); doPlay(); }
    else if (k === 'z') { if (lostShown()) doReset(); }   // only while Respawn shows
    else if (k === 'a') doFire();
    else if (k === 'x') { if (book.on) openBook(); }   // the morph book, as the Morphs / Choose buttons
    else if (k === 's') { if (!$('bullets').hidden) openBulletBook(); }   // the bullet book, as the Bullets button
    else if (k === 'c') { if (!$('species').hidden) openSpeciesBook(); }   // the species book (sandbox mode)
    else if (k === 'd') { if (Bullets.cycle(-1)) flash('bullet-prev'); }
    else if (k === 'f') { if (Bullets.cycle(1)) flash('bullet-next'); }
    else if (k === 'q') doSpeed(-1);   // slower
    else if (k === 'w') doSpeed(1);    // faster
    else if (k === 'Backspace') { e.preventDefault(); doRestart(); }
    else if (k === 'e') doView();
    else if (k === 'r') doZoom(-1);
    else if (k === 't') doZoom(1);
    else if (k === 'm') Page.toggleSound();
    else if (e.key >= '1' && e.key <= '3') C.pickChannel(+e.key - 1);   // the same one again: erase on / off
  });

  // No keyboard focus cycling (Tab, Shift+Tab) anywhere on the page, whatever has the focus (a
  // button, a slider), and no control takes part in it; nor a focus ring on a control after a click
  // or a change (style.css .game-page)
  window.addEventListener('keydown', e => { if (e.key === 'Tab') e.preventDefault(); }, { capture: true });
  for (const el of document.querySelectorAll('button, input, select, a[href]')) el.tabIndex = -1;
  document.documentElement.classList.add('game-page');

  // ---- start ------------------------------------------------------------------------------
  // The world starts empty with no species; the script, run once the simulation has started,
  // loads the species and puts the creatures in.
  applyLayout();
  if (MODE === 'story') $('msg').classList.add('attract');   // (story mode only)
  const scriptText = fetch(SCRIPT_FILE, { cache: 'no-cache' })
    .then(r => { if (!r.ok) throw new Error(`${SCRIPT_FILE}: HTTP ${r.status}`); return r.text(); });
  C.start({ empty: true }).then(() => scriptText.then(txt => {
    if (ALL) Bullets.setList('all');
    else if (MODE === 'sandbox') loadCentralInventory();
    SCRIPT = parseScript(txt);
    fitMsg();
    // [jump_here] (for writing the script): start at the first line that has it. Every action on
    // the lines before it runs first, in order, so the set-up is as if played up to there; their
    // messages and conditions are skipped, and so are the chapter cards.
    const from = Math.max(0, SCRIPT.findIndex(L => L.actions.some(a => a.name === 'jump_here')));
    replaying = true;
    for (const L of SCRIPT.slice(0, from))
      for (const a of L.actions) if (a.name !== 'chapter') ACTIONS[a.name](...a.args);
    replaying = false;
    if (from < SCRIPT.length) run(from);
  }, err => { console.error(err); say(`Could not load the script (${SCRIPT_FILE}).`); }), Page.showError);
})();
