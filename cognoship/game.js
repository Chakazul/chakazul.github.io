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
  C.opt.morphItalic = () => true;
  C.setDetectMode('focus');   // Detect: Focused, the focus's morph at the top of the ship view
  // "Dissolved" alone, centred in the ship view (the Retry button below it), at 1.8 x the message
  // area's text size (18 CSS px, 15 on a phone), at least 44 canvas px.
  C.opt.banner = () => {
    const css = parseFloat(getComputedStyle($('msg')).fontSize) || 18;
    const s = css * 1152 / ($('ui').getBoundingClientRect().width || 1152);
    return { title: Math.max(44, 1.8 * s), line: 0, inShip: true };
  };

  // Text with _italic_, *bold* and *_italic bold_* (either nesting) as nodes: plain text, <i>, <b>.
  function styled(text) {
    const out = [];
    let last = 0;
    for (const m of text.matchAll(/\*_([^*_]+)_\*|_\*([^*_]+)\*_|\*([^*]+)\*|_([^_]+)_/g)) {
      if (m.index > last) out.push(text.slice(last, m.index));
      const both = m[1] ?? m[2], el = document.createElement(both !== undefined || m[3] !== undefined ? 'b' : 'i');
      if (both !== undefined) { const i = document.createElement('i'); i.textContent = both; el.append(i); }
      else el.textContent = m[3] ?? m[4];
      out.push(el);
      last = m.index + m[0].length;
    }
    if (last < text.length) out.push(text.slice(last));
    return out;
  }

  // ---- messages -----------------------------------------------------------------------------
  // The message area under the world view.
  // hint: a grey note after the text, e.g. what key goes on.
  // One inner span holds the text and the hint, so they wrap as one paragraph.
  // The text's markup: _italic_, *bold*, and *_both_* (or _*both*_).
  function say(text, hint = '') {
    const p = document.createElement('span');
    p.append(...styled(text));
    if (hint) { const h = document.createElement('span'); h.className = 'hint'; h.textContent = ' ' + hint; p.append(h); }
    $('msg').replaceChildren(p);
  }

  // The script (script.txt), run line by line. A line is a message, then any number of [tokens]:
  // [name] or [name:arg:arg...], each an action (ACTIONS, run when the line is reached, in order)
  // or a condition (CONDITIONS, checked every frame; the line waits until all of them hold).
  //   - A line with text shows it. With no condition it moves on by Space or a click on the
  //     world view; with conditions, when they hold (Space and clicks then do nothing).
  //   - A line of tokens only shows nothing (the last message stays): it runs its actions, then
  //     moves on at once, or when its conditions hold.
  // Message text may use _italic_, *bold* and *_italic bold_* (say / styled).
  // Blank lines and comment lines (starting with # or !) are skipped. Unknown names are reported in the console
  // and ignored. The last line stays. A message ends with a grey mark: one waiting on a condition
  // with ACT_MARK (the player has something to do); else the first message with FIRST_HINT, later
  // ones with NEXT_MARK (Space / click moves on; none on the last line).
  const FIRST_HINT = 'Press space or click world to continue \u226b';   // ≫
  const NEXT_MARK = '\u226b';         // ≫
  const ACT_MARK = '\u29d7';         // ⧗

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
  // and Retry. (steps > 0: not in the moment between a cleared world and the ship put in.)
  const shipLost = () => C.status !== 'alive' || (!C.focus && C.steps > 0);
  C.opt.over = shipLost;
  C.opt.pauseWhenEmpty = () => false;   // "Dissolved" does not pause the game
  // morph ids from a token's args: [allow_morph:limus:gyrans] or [allow_morph:limus,gyrans]
  const morphIds = args => args.flatMap(a => a.split(/[,\s]+/)).filter(Boolean);

  // ---- morph book -----------------------------------------------------------------------------
  // The morphs the player has discovered, empty at first. [discover_morph:<ids>] adds to it;
  // while the book is on, a morph is also added by being found: a visible morph not yet
  // discovered that the followed ship has been detected as for DISCOVER_MS of play in a row; that
  // also opens the book, titled "New morph discovered!", the new morph marked.
  // [morph_book:on] shows the Morphs button (it opens the book). The book has a slot per visible
  // morph, in the visible list's order: the name if discovered, else "?"; a click on a discovered
  // one does the [morph_book_action], but only while the ship is lost (else read-only). The game pauses while it is open. Close, a click outside the
  // box, or Escape closes it.
  const discovered = [];   // { name, hue }
  const book = { on: false, open: false, wasRunning: false, newName: null, fromDiscovery: false, action: 'clear_put',
                 discoveryCloses: 0 };   // book closes after a discovery opened it ([discover_morph])
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
  // One slot per visible morph, in the visible list's order: its name if discovered, else "?".
  // Title: "New morph discovered!" when a discovery opened the book.
  function renderBook() {
    const canPick = book.action !== 'none' && shipLost();   // a click on a discovered morph does something
    $('book-title').textContent = book.fromDiscovery ? 'New morph discovered!' : canPick ? 'Choose a morph' : 'Morph book';
    const grid = $('book-grid');
    grid.replaceChildren(...visibleList().map(name => {
      const d = discovered.find(x => x.name === name), el = document.createElement('div');
      el.className = 'book-slot' + (d ? '' : ' empty') + (d && d.name === book.newName ? ' new' : '');
      if (d) {
        const dot = document.createElement('span');
        dot.className = 'dot';
        dot.style.background = Number.isFinite(d.hue) ? `hsl(${d.hue.toFixed(0)},100%,55%)` : '#ccc';
        el.append(dot, d.name);
        // a pick only while the ship is lost (else the book is read-only)
        if (canPick) {
          el.classList.add('pick');
          el.title = book.action === 'clear_put' ? `Start again with ${d.name}` : `Put a ${d.name} into the world`;
          el.onclick = () => pickMorph(d.name);
        }
      } else el.textContent = '?';
      return el;
    }));
  }
  // A click on a discovered morph in the book: close it, do the morph_book_action, and play (even
  // if the game was paused before the book opened).
  function pickMorph(name) {
    closeBook();
    // either way the new creature is followed
    if (book.action === 'clear_put') { ACTIONS.clear_world(); ACTIONS.put_player(name, 'center'); }
    else if (book.action === 'put') ACTIONS.put_player(name, 'random');
    C.setRunning(true);
  }
  function openBook(newName = null) {
    book.newName = newName; book.fromDiscovery = book.fromDiscovery || !!newName;
    renderBook();
    if (book.open) return;
    book.open = true; book.wasRunning = C.running;
    C.setRunning(false);
    $('book').hidden = false;
    $('book-close').focus();
  }
  function closeBook() {
    if (!book.open) return;
    book.open = false;
    $('book').hidden = true;
    if (book.fromDiscovery) book.discoveryCloses++;
    book.fromDiscovery = false; book.newName = null;
    if (book.wasRunning) C.setRunning(true);
  }
  function setBook(on) { book.on = on; $('morphs').hidden = !on; if (!on) closeBook(); }
  $('morphs').onclick = () => openBook();
  $('book-close').onclick = closeBook;
  $('book').addEventListener('pointerdown', e => { if (e.target === $('book')) closeBook(); });   // outside the box
  // Each frame: has the followed ship been a visible morph not discovered yet, without a break,
  // for DISCOVER_MS of play (paused time does not count, whatever the speed)? A passing
  // detection is not a discovery.
  const DISCOVER_MS = 2000;
  let cand = null, lastT = 0;   // cand: { name, ms } the morph the ship is in now, if it would be new
  function checkDiscovery() {
    const now = performance.now(), dt = Math.min(100, now - lastT);   // capped: a hidden tab's gap
    lastT = now;
    const d = C.focus && C.stats.short.det;
    const known = d && d.name && discovered.find(x => x.name === d.name);
    if (known && !Number.isFinite(known.hue) && Number.isFinite(d.hue)) known.hue = d.hue;   // put in unseen: its colour now
    if (!book.on || !d || !d.name || !isVisible(d.name) || known) { cand = null; return; }
    if (!cand || cand.name !== d.name) { cand = { name: d.name, ms: 0 }; return; }
    if (C.running) cand.ms += dt;
    if (cand.ms < DISCOVER_MS) return;
    cand = null;
    discover(d.name, d.hue);
    openBook(d.name);
  }
  // Add a morph to the discovered list (its hue, the book's dot colour, once it has been detected).
  function discover(name, hue = NaN) {
    if (!discovered.some(x => x.name === name)) discovered.push({ name, hue });
  }

  // The player's morph, set by [put_player] and a pick in the morph book: what Retry puts back.
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
  const PUT_HOW = ['center', 'random', 'far'];
  const ACTIONS = {
    // species:<id>: load species <id> (its rule; the world is kept), with no visible morphs
    species: id => {
      const p = PATTERNS.find(p => C.speciesOf(p) === id);
      if (!p) return console.warn(`script.txt: unknown species [species:${id}]`);
      C.setPattern(p.id);
      vis.all = false; vis.allow.clear();
    },
    // allow_morph:<ids>: the visible morphs are exactly these, in this order ("all": every morph of
    // the species; no ids: none)
    allow_morph: (...args) => {
      const ids = morphIds(args);
      vis.all = ids.includes('all');
      vis.allow = new Set(ids.filter(id => id !== 'all'));
    },
    // morph_book:on|off: show or hide the Morphs button (off also closes the book)
    morph_book: (state = 'on') => setBook(state !== 'off'),
    // clear_world: empty the world
    clear_world: () => C.clearWorld(),
    // put_morph:<id>[:center|random|far]: put morph <id> of the current species into the world:
    // at the centre (default, unrotated), at a random empty spot and angle (the lab's Add), or at
    // a random angle in the empty spot farthest from the followed ship. The ship stays followed
    // (if none, the new creature is).
    put_morph: (id, how = 'center') => putMorph('put_morph', id, how, false),
    // put_player:<id>[:center|random|far]: as put_morph, but the new creature is the player's
    // ship (followed), and <id> becomes the player's morph: what Retry puts back.
    put_player: (id, how = 'center') => { if (putMorph('put_player', id, how, true)) playerMorph = id; },
    // discover_morph:<ids>: add these morphs to the discovered list (the player may have found
    // others already: they stay); with no ids it is the condition
    discover_morph: (...args) => { for (const id of morphIds(args)) discover(id); },
    // morph_book_action:<type>: what a click on a discovered morph in the book does (it also closes
    // the book): clear_put = [clear_world][put_player:<id>]; put = [put_player:<id>:random];
    // none = nothing. The default is clear_put.
    // start_here: the script starts at this line (see the start, at the end of this file)
    start_here: () => {},
    // chapter:<no>:<title>: a title card over the game window, "Chapter <no>" above <title>, faded in; any
    // key or click closes it (fading out) and is used up. The game plays on meanwhile, and the
    // line waits until it is closed.
    chapter: (no = '', ...title) => { showChapter(no, title.join(':')); return () => !chapterOpen; },
    // retry_action:<type>: what Retry (and Z) does: clear_put, put or none as morph_book_action,
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
    // discover_morph: the player discovers a new morph, and closes the morph book it opened
    discover_morph: () => {
      const c0 = book.discoveryCloses;
      return () => book.discoveryCloses > c0;
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

  let SCRIPT = [], line = -1, tests = [], firstText = -1;
  function parseScript(txt) {
    return txt.split(/\r?\n/).map(l => l.trim()).filter(l => l && !l.startsWith('#') && !l.startsWith('!')).map(l => {
      const toks = [];
      let m;
      while ((m = l.match(/\s*\[([^\]]*)\]\s*$/))) { toks.unshift(m[1].trim()); l = l.slice(0, m.index); }
      const actions = [], conds = [];
      for (const t of toks) {
        const [name, ...args] = t.split(':').map(x => x.trim());
        // a name that is both (discover_morph): the action with args, the condition without
        if (ACTIONS[name] && (args.length || !CONDITIONS[name])) actions.push({ name, args });
        else if (CONDITIONS[name]) conds.push({ name, args });
        else console.warn(`script.txt: unknown action or condition [${t}], ignored`);
      }
      return { text: l.trim(), actions, conds };
    });
  }
  function run(i) {
    line = i;
    const L = SCRIPT[i];
    // an action may return a test too: the line also waits for it (chapter: until the title closes)
    const waits = L.actions.map(a => ACTIONS[a.name](...a.args)).filter(f => typeof f === 'function');
    tests = L.conds.map(c => CONDITIONS[c.name](...c.args)).concat(waits);
    if (L.text) {
      say(L.text, tests.length ? ACT_MARK : i === firstText ? FIRST_HINT : i + 1 < SCRIPT.length ? NEXT_MARK : '');
      $('msg').classList.toggle('wait', tests.length > 0);   // yellow: the player has something to do
    }
    else if (!tests.length) next();   // tokens only, nothing to wait for
  }
  function next() { if (line + 1 < SCRIPT.length) run(line + 1); }
  // Space / a world-view click moves on: a message with no condition
  const byContinue = () => line >= 0 && !!SCRIPT[line].text && !tests.length;
  function checkScript() {
    if (line < 0 || !tests.length) return;
    const ok = tests.map(t => t());   // every test, every frame (one may be watching for an action)
    if (ok.every(Boolean)) next();
  }
  // a click on the world view (the canvas's left 2/3), besides focusing the creature there
  $('ui').addEventListener('pointerdown', e => {
    const r = $('ui').getBoundingClientRect();
    if (e.clientX - r.left < r.width * 2 / 3 && byContinue()) next();
  });

  // ---- chapter title card ([chapter]) ----------------------------------------------------------
  const FADE_MS = 1500;   // as .chapter's transition in style.css
  let chapterOpen = false, chapterTimer = null;
  function showChapter(no, title) {
    $('chapter-no').textContent = `Chapter ${no}`; $('chapter-title').textContent = title;
    say(''); $('msg').classList.remove('wait');   // a new chapter: the message bar starts empty
    clearTimeout(chapterTimer);
    const el = $('chapter');
    el.hidden = false;
    void el.offsetWidth;   // start the fade from opacity 0
    el.classList.add('show');
    chapterOpen = true;
  }
  function closeChapter() {
    chapterOpen = false;
    const el = $('chapter');
    el.classList.remove('show');
    chapterTimer = setTimeout(() => { el.hidden = true; }, FADE_MS);
  }
  // Any key or click anywhere closes it, before anything else sees the event (so it does nothing else).
  for (const ev of ['keydown', 'pointerdown'])
    window.addEventListener(ev, e => {
      if (!chapterOpen || (ev === 'keydown' && ['Shift', 'Control', 'Alt', 'Meta'].includes(e.key))) return;
      e.preventDefault(); e.stopImmediatePropagation();
      closeChapter();
    }, { capture: true });

  // ---- core hooks -------------------------------------------------------------------------
  C.hooks.change = () => syncButtons();
  C.hooks.frame = () => { checkDiscovery(); checkScript(); if ($('reset').hidden === shipLost()) syncButtons(); };

  // ---- controls ---------------------------------------------------------------------------
  function syncButtons() {
    $('play').setAttribute('aria-pressed', C.running);   // "Play", highlighted while playing
    $('reset').hidden = !shipLost();
    $('morphs-dead').hidden = !shipLost() || !book.on;
    $('headup').setAttribute('aria-pressed', C.viewMode !== 'board');
    $('speed').value = C.speedIdx + 1;
    $('zoom').value = C.zoomIdx + 1;
    for (const el of document.querySelectorAll('input[type=range]')) fillRange(el);
  }

  const doPlay = () => { flash('play'); C.togglePlay(); };
  // Retry (and Z), by [retry_action:<type>]: clear_put (clear the world, put the player's morph at
  // the centre; the default), put (put it at a random spot), none (nothing), or book (open the
  // morph book to pick the morph, by its morph_book_action). The Morphs button under Retry
  // (with the book on) opens the book, where a pick then does the morph_book_action. Either put makes the new creature the player's ship and plays. The player's morph:
  // the last [put_player] or book pick (none yet: the last morph put in).
  let retryAction = null;
  function doReset() {
    const type = retryAction ?? 'clear_put';
    const name = playerMorph ?? C.morphName(C.pattern).replace(C.BOOM, '');
    if (type === 'book') openBook();
    else if (type === 'clear_put') { ACTIONS.clear_world(); ACTIONS.put_player(name, 'center'); C.setRunning(true); }
    else if (type === 'put') { ACTIONS.put_player(name, 'random'); C.setRunning(true); }
  }
  const doView = () => { flash('headup'); C.toggleHeadUp(); };
  const doSpeed = d => { flash('speed'); C.changeSpeed(d); };
  const doZoom = d => { flash('zoom'); C.setZoom(Math.max(0, Math.min(C.ZOOMS.length - 1, C.zoomIdx + d))); };   // R out, T in
  $('play').onclick = doPlay;
  $('reset').onclick = doReset;
  $('morphs-dead').onclick = () => openBook();
  $('headup').onclick = doView;
  $('speed').oninput = e => C.setSpeed(+e.target.value - 1);
  $('zoom').value = C.zoomIdx + 1;
  $('zoom').oninput = e => C.setZoom(+e.target.value - 1);

  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (book.open) { if (k === 'Escape') closeBook(); return; }   // the book takes the keys while open
    if (e.code === 'Space') { e.preventDefault(); if (byContinue()) next(); }
    else if (k === 'Enter') { e.preventDefault(); doPlay(); }
    else if (k === 'z') doReset();
    else if (k === 'q') doSpeed(1);
    else if (k === 'w') doSpeed(-1);
    else if (k === 'e') doView();
    else if (k === 'r') doZoom(-1);
    else if (k === 't') doZoom(1);
    else if (k === 'm') Page.toggleSound();
    else if (k === 'n') { flash('sound'); Page.nextTrack(); }
    else if (e.key >= '1' && e.key <= '3') C.channel = +e.key - 1;
  });

  // ---- start ------------------------------------------------------------------------------
  // The world starts empty with no species; the script, run once the simulation has started,
  // loads the species and puts the creatures in.
  const scriptText = fetch('script.txt', { cache: 'no-cache' })
    .then(r => { if (!r.ok) throw new Error(`script.txt: HTTP ${r.status}`); return r.text(); });
  C.start({ empty: true }).then(() => scriptText.then(txt => {
    SCRIPT = parseScript(txt);
    // [start_here] (for writing the script): start at the first line that has it, skipping
    // everything before it, actions included (so that line needs its own set-up, e.g. [species])
    const from = Math.max(0, SCRIPT.findIndex(L => L.actions.some(a => a.name === 'start_here')));
    firstText = SCRIPT.findIndex((L, i) => i >= from && L.text);
    if (from < SCRIPT.length) run(from);
  }, err => { console.error(err); say('Could not load the script (script.txt).'); }), Page.showError);
})();
