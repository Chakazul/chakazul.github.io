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
  // the R G B patch again (or its key): draw if allowed (below), then erase ("recycle"), striped,
  // -R / -G / -B, allowed unless [erase:off]: mask -> draw -> erase -> mask
  let eraseOn = true;
  C.opt.patchErase = () => eraseOn;
  // ... and (before erase, if on) draw: add mass, cross-hatched, +R / +G / +B; allowed by [draw]
  // (off by default): on, or if the current species is a given one (checked as it changes)
  let drawRule = { on: false, species: null };
  const drawOn = () => drawRule.on && (!drawRule.species || (!!C.pattern && C.speciesOf(C.pattern) === drawRule.species));
  C.opt.patchDraw = drawOn;
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
  // central inventory (see the inventory, below); ?mode=challenge: challenge mode (challenge.txt, see
  // the challenges, below), from the central inventory too but with its own energy; else story mode
  // (script.txt), from an empty one. The title row's links switch between them (other URL
  // parameters kept).
  const MODE = { sandbox: 'sandbox', challenge: 'challenge' }[new URLSearchParams(location.search).get('mode')] ?? 'story';
  const SCRIPT_FILE = { sandbox: 'sandbox.txt', challenge: 'challenge.txt' }[MODE] ?? 'script.txt';
  // ?all=1: all discovered (for testing): every visible morph of every species counts as discovered
  // (allDiscovered) and the gun has every bullet pattern from the start; the central inventory is
  // neither read nor written (in either mode).
  const ALL = new URLSearchParams(location.search).get('all') === '1';
  for (const [id, m] of [['mode-story', 'story'], ['mode-sandbox', 'sandbox'], ['mode-challenge', 'challenge']]) {
    const q = new URLSearchParams(location.search);
    if (m !== 'story') q.set('mode', m); else q.delete('mode');
    const s = q.toString();
    $(id).href = s ? `?${s}` : location.pathname;
    $(id).classList.toggle('on', m === MODE);
    if (m === MODE) $(id).setAttribute('aria-current', 'page');
  }
  C.opt.circlesAll = () => DEV;
  C.opt.worldFrame = () => true;   // the world view framed as the ship view, its corners rounded alike
  C.opt.wholeNames = () => false;   // (?dev=1's circles: a twin's without its name; a hostile's named by its red mark)
  C.opt.wholes = () => true;   // twins and trios found, named and followed as one (core.js wholes, ship)
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
    return { title: Math.max(44, 1.8 * s), line: 0, inShip: true, dy: -40,   // a little above the middle
             text: lostFails() ? 'Failed' : 'Dissolved' };   // challenge mode: a lost ship ([dissolve:retry]), time's up or no energy fails it
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
        const m = document.createElement('span'); m.className = 'mark'; m.dataset.key = 'Space'; m.textContent = marks.next;
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
                   survive_time: ([s]) => `time: ${s}/${s}`, survive: ([n]) => n ? `time: 999/${n}` : 'time: 999' };
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
    layer.append($('bullet-nav'), $('fire'), $('bullets'));   // (Pulses with them: right of Fire, < > above it; stacked, top down: Fire, < >, Pulses)
    let col = $('dev-col');
    if (!col) {   // (made here if the page lacks it, e.g. an index.html cached from before it)
      col = document.createElement('div'); col.id = 'dev-col'; col.className = 'dev-col';
      document.querySelector('.layer-ship').append(col);
    }
    (stacked ? col : layer).append($('dev'), $('dev-fps'));   // ?dev=1: left of the ship view, a column
    fitMsg();
  }
  // Pulses: beside Fire in the world view (also with an index.html cached from when it was in the
  // ship view: moved there, its class set)
  if (!$('bullets').classList.contains('world-bullets')) {
    $('bullets').className = 'world-bullets';
    document.querySelector('.layer-world').append($('bullets'));
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
  // stats only).
  const vis = { all: false, allow: new Set() };
  const isVisible = name => vis.all || vis.allow.has(name);
  C.opt.morphLabel = det => isVisible(det.name) ? det.name : null;   // a fuzzy morph too: no "~"
  // Hidden morphs ([hide_morph]): still detected and named where visible (name box, world view), but
  // left out of the morph book (so of respawn, retry and clone picks) and the species book's
  // thumbnails, and never discovered. Keys "<code>:<name>", or "*:<name>" for every species.
  const hiddenMorphs = new Set();
  const curCode = () => C.pattern ? C.speciesOf(C.pattern) : '';
  const isHidden = (name, code = curCode()) => hiddenMorphs.has(`*:${name}`) || hiddenMorphs.has(`${code}:${name}`);
  function hideMorphs(morphs = '', code) {   // [hide_morph]'s change to hiddenMorphs
    const names = morphIds([morphs]);
    if (names.includes('none')) {
      for (const k of [...hiddenMorphs]) if (!code || k.startsWith(code + ':')) hiddenMorphs.delete(k);
    } else for (const n of names) hiddenMorphs.add(`${code || '*'}:${n}`);
  }
  let baseHidden = new Set();   // challenge mode: the hidden morphs after its set-up lines (each challenge starts from them)
  // The ship is the one followed: a click on the world view does not switch to another creature.
  C.opt.clickFocus = () => false;
  C.opt.noFocusText = () => [];   // nothing followed: the ship view stays empty
  // The player is dead when the followed ship is gone (other creatures may live on): "Dissolved"
  // and Respawn. (steps > 0: not in the moment between a cleared world and the ship put in.)
  const shipLost = () => C.status !== 'alive' || (!C.focus && C.steps > 0);
  // [clear_world:empty]: a scene of empty space, no player: while the ship is lost the ship view
  // shows the world's centre, with no "Dissolved", Respawn or Choose; until the next [clear_world]
  let emptySpace = false;
  // The banner reads "Failed" (not "Dissolved"): the challenge has failed, or the ship is lost
  // where that fails it ([dissolve:retry]; from the first frame, before challengeFrame calls it)
  const lostFails = () => MODE === 'challenge' && (!!challenge.failed || (!!challenge.cur && challenge.dissolveFails));
  const lostShown = () => (shipLost() && !emptySpace) || !!challenge.failed;   // "Dissolved" (challenge mode: "Failed") and Respawn / Choose show
  // Autospawn on ([settings]): no "Dissolved", Respawn or Choose while it waits; the ship view stays
  // on the spot where the ship was last seen, as it was turned (lastShip)
  const autoLost = () => settingsOn && autospawn && !lostFails();
  const lostPanel = () => lostShown() && !autoLost();
  let lastShip = null;   // [col, row, phi] of the ship, each frame it is followed
  C.opt.over = lostPanel;
  C.opt.emptyView = () => emptySpace || (autoLost() && shipLost() && lastShip) || false;
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
    $('book-grid').classList.toggle('bullets', book.mode === 'bullet' || book.mode === 'hostile');
    $('book-grid').classList.toggle('species', book.mode === 'species' || book.mode === 'challenge');
    document.querySelector('.book-box').classList.toggle('wide', book.mode !== 'morph');   // the bullet and species books: wider
    document.querySelector('.book-box').classList.toggle('challenges', book.mode === 'challenge');   // the challenge book: a little wider still (the ticks)
    $('book-close').textContent = bookBack() ? 'Back' : 'Close';   // the species book for a challenge: back to the challenge book
    $('book-close').title = !bookBack() ? '' : book.mode === 'hostile' ? 'Back to the settings (Esc)' : 'Back to the challenge book (Esc)';
    $('book-close').dataset.key = 'Esc';   // only Escape closes a book (its letters pick slots)   // its key badge: the key that closes this book
    // the morph book's species, at its lower right (a morph is only a morph of its species)
    $('book-species').hidden = book.mode !== 'morph' || !C.pattern;
    $('book-settings').hidden = book.mode !== 'settings';   // ([settings:on]: the settings book, in place of the slots)
    $('book-grid').hidden = book.mode === 'settings';
    syncObjects();
    if (C.pattern) {
      const code = document.createElement('i');
      code.textContent = speciesName || C.speciesOf(C.pattern);   // no name: the id
      $('book-species').replaceChildren('Species: ', code);
    }
    // sandbox mode: the species' name goes to the species book
    $('book-species').classList.toggle('link', MODE === 'sandbox');
    $('book-species').title = MODE === 'sandbox' ? 'Choose another species (D)' : '';
    $('book-clear').hidden = MODE === 'story' || book.mode === 'settings' || book.mode === 'hostile';
    // the challenge book's: clears the completion (the ticks) and the best times
    // (the species book for a challenge: that challenge's, with every species)
    const recs = book.mode === 'challenge' || (book.mode === 'species' && challenge.pending);
    const clearWhat = recs ? 'completion & best time' : 'data';
    if (!$('book-clear').classList.contains('done')) $('book-clear').textContent = `Hold to clear ${clearWhat}`;
    $('book-clear').title = clearWhat === 'data' ? 'Hold for 3 seconds to clear every discovered morph and pulse, and the energy, from every run'
      : book.mode === 'species' ? `Hold for 3 seconds to clear the completion and best times of ${challenge.pending.title}, with every species`
      : 'Hold for 3 seconds to clear the completion and best time of every challenge';
    if (book.mode === 'bullet') return renderBullets();
    if (book.mode === 'species') return renderSpecies();
    if (book.mode === 'settings') return setTitle(['Settings']);
    if (book.mode === 'hostile') return renderHostilePulses();
    if (book.mode === 'challenge') return renderChallenges();
    allDiscovered();
    // a click on a discovered morph does something: (challenge mode: book.retry, the retry book) restart the challenge with it
    const canPick = book.retry || (book.action !== 'none' && shipLost());
    const canClone = clone.on && !shipLost() && !book.fromDiscovery;   // [clone:on]: pick one to clone (not in the discovery's book)
    setTitle(book.fromDiscovery ? ['New morph discovered!'] : book.retry ? ['Choose a morph to ', ['retry', 'choose']]
      : canPick ? ['Choose a morph to ', ['respawn', 'choose']]
      : canClone ? ['Choose a morph to ', ['clone', 'clone']] : ['Dissolve to respawn']);   // read-only: how to get a pick
    const grid = $('book-grid');
    grid.replaceChildren(...visibleList().filter(n => !isHidden(n)).map((name, slot) => {
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
            el.title = book.retry ? `Retry the ${MODE === 'challenge' ? 'challenge' : 'scene'} with ${d.name}` : book.action === 'clear_put' ? `Start again with ${d.name}` : `Put a ${d.name} into the world`;
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
    if (book.retry) {
      closeBook(); playerMorph = name;
      if (MODE === 'challenge') return retryChallenge();
      // story, sandbox: the scene again, once the book has zapped out (resetScene waits for none)
      const go = () => bookShowing() ? requestAnimationFrame(go) : retryScene();
      return go();
    }
    closeBook();
    // either way the new creature is followed
    if (book.action === 'clear_put') { C.clearWorld(); ACTIONS.put_player(name, 'center'); }   // (not the script's [clear_world]: background and [min_morph] kept)
    else if (book.action === 'put') ACTIONS.put_player(name, 'random');
    if (book.action !== 'none') { retake = true; cand = null; respawns++; respawned(); }   // [need:respawn]; its picture may be retaken
    C.setRunning(true);
  }
  // ---- clone mode ([clone:on]) ------------------------------------------------------------
  // While it is on and the ship is alive, the morph book offers the discovered morphs to clone.
  // Picking one closes the book; the morph then follows the pointer over the ship view (instead
  // of the tool circle), turned by ↻ (45 degrees clockwise) and flipped by ↔ (left to right, as
  // seen), until a click on the ship view pastes it there, or the Morphs button cancels it. Keys:
  // F flips, G turns, Enter pastes (where it is shown), Escape cancels. Clone mode stays allowed
  // until [clone:off].
  const newBullets = new Set();   // the patterns created in the bullet book, green in it until it is closed
  let bulletNew = false;   // a new pattern can be created ([new_bullet]): the Bullets button reads "Create", green
  let revealName = null;   // the pattern just created: revealed as a new card as the book re-renders
  let clonesPlaced = 0;   // clones pasted into the world so far ([need:first_clone] counts on from its line)
  const clone = { on: MODE !== 'challenge', held: null };   // on: allowed (the default, until [clone:off])   // held: { name, M (board-frame 2x2), src form, form, img }
  function startClone(name) {
    const p = findMorph(name);
    if (!p) return;
    // as before the book opened; picking again while holding one: as before the first pick
    const wasRunning = clone.held ? clone.held.wasRunning : book.wasRunning;
    closeBook();
    C.setRunning(false);   // paused until it is placed (or cancelled), then as it was
    clone.held = Object.assign(Clone.hold(p, C.R, name), { wasRunning });
    if (lastMouse) C.hoverAt(lastMouse);   // shown at the mouse at once, before it moves
    syncButtons();
  }
  // The mouse's last position (a pointer event's clientX / clientY), so a clone picked in the book
  // (by a click or a key) shows under the mouse before it moves; null after a touch
  let lastMouse = null;
  for (const ev of ['pointermove', 'pointerdown'])
    document.addEventListener(ev, e => {
      lastMouse = e.pointerType === 'mouse' ? { clientX: e.clientX, clientY: e.clientY, pointerType: 'mouse' } : null;
    }, { capture: true });
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
    if (!H && fieldMode()) { if (hoverAt) drawBrush(ctx, hoverAt[0], hoverAt[1], Field.widthOf(C.fieldTool) / 2 * cellPx); return true; }   // the brush, not the ship tool's circle
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
    if (!H && Field.active && C.focus) { const sh = C.ship(); Field.move([sh.col + off[0], sh.row + off[1]]); return; }   // painting the field
    if (!H || !H.at || !C.focus || (e && e.pointerType === 'mouse')) return;
    const sh = C.ship();   // (off is from the ship view's centre: the ship, a twin's whole)
    H.at = [sh.col + off[0], sh.row + off[1]];
  };
  C.opt.zoomPress = (off, e) => {
    const H = clone.held;
    if (!H && fieldMode() && C.focus) { const sh = C.ship(); Field.begin([sh.col + off[0], sh.row + off[1]], C.channel, C.fieldTool); return true; }
    if (!H || !C.focus) return false;
    const sh = C.ship(), at = [sh.col + off[0], sh.row + off[1]];
    if (e && e.pointerType !== 'mouse') { H.at = at; syncButtons(); return true; }
    pasteClone(at);
    return true;
  };
  // The player boards the clone: the game follows it from now on (the old ship goes on as another
  // creature; pasted onto it, the two are one object, followed as before), and its morph is what
  // Respawn puts back. So a morph can be changed without dissolving, and a twin engineered by
  // becoming one half and pasting the other next to it.
  function pasteClone(at) {
    if (cloneCost && energy < cloneCost) { flash('morphs'); return; }   // [clone_cost]: not enough energy (the clone stays held)
    const H = clone.held;
    Clone.paste(H, at);
    if (cloneCost) changeEnergy(-cloneCost);
    C.analyze(0);
    const i = C.found ? Creatures.pick(C.found, at[0], at[1], 2 * C.R, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
    if (i >= 0) { C.setFocus(i); playerMorph = H.name; cand = null; }   // (no retake: a world with the old ship in it)
    clonesPlaced++;   // [need:first_clone]
    endClone();
  }
  $('clone-rot').onclick = rotateClone;
  $('clone-flip').onclick = flipClone;
  $('clone-cancel').onclick = () => { flashGo('clone-cancel'); endClone(); };
  $('clone-ok').onclick = () => { if (clone.held && clone.held.at) { flashGo('clone-ok'); pasteClone(clone.held.at); } };
  // Enter: paste where the clone is shown (placed by a touch, else under the mouse on the ship view)
  function confirmClone() {
    const H = clone.held;
    if (!H || !C.focus) return;
    const at = H.at || C.pointerCell(true);
    if (at) pasteClone(at);
  }

  // ---- world draw ([world_draw], [add_bullet:nebula|antimatter]: the bullet types Nebula and Antimatter) -------------------
  // While one of them is in use (fieldMode) the player paints the R G B fields (field.js) with a
  // press in either view: the R G B patches pick the colour, the selected one again the tool (+
  // paintbrush, × soft brush, − eraser); the Zoom slider is Width (2 to 32 cells; the soft brush and the
  // eraser twice that), Fire is Clear (the
  // selected colour's field emptied), Upward is Undo (the last stroke or Clear; again: redo it).
  // The bots don't fire meanwhile (their bullets' channel is a field).
  const fieldMode = () => !!Bullets.field;
  SimGL.setFieldKeepsDots(true);   // the dots (and beacons, to be eaten) stay while drawing on the world
  C.opt.fieldTools = fieldMode;
  C.opt.worldCursor = fieldMode;   // the crosshair over the world view only while drawing on it (else nothing to do there)
  C.opt.worldPress = at => { if (!fieldMode()) return false; Field.begin(at, C.channel, C.fieldTool); return true; };
  C.opt.worldDrag = at => Field.move(at);
  C.opt.pressEnd = () => Field.end();
  // the brush's outline: in the colour being painted; solid (paintbrush), dotted (soft brush), dashed (eraser)
  function drawBrush(ctx, x, y, r) {
    ctx.save();
    ctx.strokeStyle = C.CH_COLORS[C.channel]; ctx.lineWidth = 2; ctx.setLineDash([[], [2, 4], [6, 4]][C.fieldTool] || []);
    ctx.beginPath(); ctx.arc(x, y, Math.max(2, r), 0, 2 * Math.PI); ctx.stroke();
    ctx.restore();
  }
  C.opt.worldGhost = hover => fieldMode() && hover
    ? { at: hover, r: Field.widthOf(C.fieldTool) / 2 * C.WS + 2, draw: (ctx, x, y) => drawBrush(ctx, x, y, Field.widthOf(C.fieldTool) / 2 * C.WS) } : null;
  Field.hooks.change = () => syncButtons();
  const doClear = () => { flash('fire'); Field.clear(C.channel); };
  const doUndo = () => { if (Field.undo()) flash('headup'); };
  // A small picture of a field type for the bullet book: a soft cloud of the three colours (Nebula),
  // or the same in cyan, magenta and yellow (Antimatter)
  const fieldPics = new Map();
  function fieldPic(name) {
    if (!fieldPics.has(name)) {
      const n = 32, cv = document.createElement('canvas');
      cv.width = cv.height = n;
      const g = cv.getContext('2d');
      // three soft blobs: Nebula's red, green, blue; Antimatter's cyan, magenta, yellow
      const cols = name === 'Antimatter' ? ['0,255,255', '255,59,255', '255,255,59'] : ['255,59,59', '59,255,90', '59,123,255'];
      [[11, 13], [21, 12], [16, 21]].forEach(([x, y], i) => {
        const gr = g.createRadialGradient(x, y, 0, x, y, 12);
        gr.addColorStop(0, `rgba(${cols[i]},.9)`); gr.addColorStop(1, `rgba(${cols[i]},0)`);
        g.fillStyle = gr; g.fillRect(0, 0, n, n);
      });
      fieldPics.set(name, cv.toDataURL());
    }
    const img = new Image();
    img.src = fieldPics.get(name); img.alt = ''; img.className = 'bpic';
    img.style.width = img.style.height = '28px';
    return img;
  }

  // The bullet book shares the morph book's window (book.mode 'bullet'): every bullet pattern by
  // name, 4 to a row (as the morph book); the ones the gun has ([bullet]) can be
  // clicked to switch to them (the one in use with a yellow border), the others are "?".
  // Keys picking a slot, by its place: the books' rows of 4 on the keyboard's rows (morphs: ASDF,
  // ZXCV; bullets: 1234, QWER, ASDF, ZX..).
  const SLOT_KEYS = { morph: [...'ASDFZXCV'], bullet: [...'1234QWERASDFZXCV'], species: [...'QWERASDFZXCV'] };
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
    setTitle(create ? [['Create', 'create'], ' or ', ['use', 'choose'], ' a pulse type'] : [['Use', 'choose'], ' a pulse type']);
    const have = Bullets.list, fields = Bullets.fieldKinds;   // (world draw: Nebula, Antimatter after the patterns, always usable)
    let nNew = 0;
    $('book-grid').replaceChildren(...Object.keys(Bullets.KINDS).concat(fields).map((name, slot) => {
      const el = document.createElement('div'), isField = fields.includes(name), ok = isField || have.includes(name);
      const current = name === Bullets.kind && !revealing.has(name);   // (yellow once revealed)
      el.className = 'book-slot' + (ok ? ' pick' : ' empty') + (current ? ' current' : '');   // (a new one: revealed, no colour of its own)
      if (newBullets.has(name) && book.reveal) revealSlot(el, nNew++);
      if (name === revealName) revealSlot(el, 0, 0);   // just created: at once
      const label = document.createElement('span');
      label.className = 'label';
      if (ok) label.append(isField ? fieldPic(name) : bulletPic(name), name);   // its stored pattern, small, in white, before the name
      else label.textContent = '\u00a0';
      el.append(label);
      if (ok) {
        el.title = isField ? `Draw ${name} on the world (it ${name === 'Nebula' ? 'masks' : 'erases'} the creatures' channels)` : `Fire ${name}`;
        el.onclick = () => { Bullets.select(name); closeBook(); };
        if (SLOT_KEYS.bullet[slot]) el.dataset.key = SLOT_KEYS.bullet[slot];
      } else if (create) {
        el.classList.add('pick', 'create');
        el.title = 'Create this new pulse'; el.onclick = () => createBullet(name);
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
  // The Species button (D), or the species' name in the morph book: a slot per visible species
  // ([add_allow_species]), in order: a small thumbnail of each of its visible morphs in a row (the
  // centre of its picture if discovered, else black), and its name. The current one has a
  // yellow border. A click chooses it (chooseSpecies: its [on_start] block sets the world up).
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
    if (MODE !== 'sandbox') return;   // (challenge mode: pickChallenge shows it)
    if (book.open && book.mode !== 'species') closeBook();
    book.mode = 'species';
    showBook();
  }
  const WIDE_THUMBS = ['G6G6CR'];   // species whose morphs are large (Cilium): every thumbnail the whole picture
  function renderSpecies() {
    const ch = challenge.pending;   // challenge mode: the species for this challenge
    setTitle(ch ? [['Choose', 'choose'], ` a species for ${ch.title}`] : [['Choose', 'choose'], ' a species']);
    const cur = C.pattern && C.speciesOf(C.pattern);
    const hid = ch ? challengeHidden(ch) : null;   // challenge mode: that challenge's hidden morphs (its [hide_morph] tags), not those in force now
    const hidden = (m, code) => hid ? hid.has(`*:${m}`) || hid.has(`${code}:${m}`) : isHidden(m, code);
    $('book-grid').replaceChildren(...(ch ? challengeSpecies(ch) : visSpecies).map((e, slot) => {
      const found = discoveredBy.get(e.id) || [];
      const el = document.createElement('div'), thumbs = document.createElement('div'), label = document.createElement('span');
      el.className = 'book-slot pick species' + (e.id === cur ? ' current' : '');
      thumbs.className = 'thumbs';
      for (const m of speciesMorphs(e).filter(m => !hidden(m, e.id))) {
        const d = found.find(x => x.name === m), t = document.createElement('div');
        const twin = (AQUARIUM233.rules[e.id]?.wholes || []).includes(m) || WIDE_THUMBS.includes(e.id);   // recorded as a whole (a twin, a trio), or a species of large creatures
        t.className = 'thumb' + (d ? '' : ' unknown') + (twin ? ' whole' : '');   // not discovered: the box's colour (discovered with no picture: black)
        if (d) { t.title = m; if (d.img) { const img = new Image(); img.src = d.img; img.alt = m; t.append(img); } }
        thumbs.append(t);
      }
      label.className = 'label';
      label.textContent = e.name || e.id;
      if (ch) {   // challenge mode: this challenge's record with this species, in a row of its own below the name: a green tick and the best time (blank: not completed yet)
        const r = loadRecords()[recordKey(ch.id, e.id)], row = document.createElement('div');
        row.className = 'crow';
        if (r && r.done) {
          const tick = document.createElement('span'), best = document.createElement('span');
          tick.className = 'tick'; tick.textContent = '✓'; tick.title = 'Completed';
          best.className = 'best'; best.textContent = Number.isFinite(r.best) ? fmtTime(r.best) : '';
          row.append(tick, best);
        } else row.textContent = '\u00a0';   // (the row kept: every slot the same height)
        el.append(thumbs, label, row);
        if (book.reveal && newRecord && newRecord.id === ch.id && newRecord.species === e.id) revealSlot(el, 0);   // the new record: revealed as it opens
      } else el.append(thumbs, label);
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
    if (challenge.pending) return startChallenge(challenge.pending, id);
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

  // ---- challenges (challenge mode) --------------------------------------------------------------
  // challenge.txt: the lines before the first [challenge] set the game up (their actions run once,
  // at the start); then each [challenge:<name>:<n>] line starts a challenge, up to the next one.
  // The challenge book (shown at the start; the Challenges button, C) lists them, a green tick on
  // each one completed, with its best time. A pick asks for a species (one its [allow_species]
  // allows; if only one, that one), then runs the challenge from its [challenge] line. Reaching its
  // end completes it: its time (game seconds from its start) is recorded in localStorage
  // (CHALLENGE_KEY, the best kept) and the book opens again. Retry after that plays it again.
  const CHALLENGE_KEY = 'cognoship-challenges';
  // list: { id, name, stage, title, start, end (line indices), allow (species ids, or ['all']) };
  // cur: the one being played; pending: picked, its species being chosen; last: the last one played;
  // failed: why the one being played failed (null: it hasn't): 'lost' (the ship dissolved), 'time'
  // ([time_limit] ran out), 'energy' ([mask_cost] took it to 0); deadline: the step [time_limit]
  // ends at (null: none); maskCost: energy per game second of a press ([mask_cost])
  // dissolveFails: [dissolve:retry], a lost ship fails the challenge
  // longest: the attempt waits on [need:survive], so its best time is its longest
  const challenge = { list: [], cur: null, pending: null, last: null, t0: 0, failed: null, deadline: null, maskCost: 0, eraseCost: 0, drawCost: 0, worldDrawCost: 0, dissolveFails: false, longest: false };
  // The challenge book's "Hold to clear completion & best time": every record goes (ticks, best times)
  // id: that challenge's records only (with every species); none: every record
  function clearCompletion(id = null) {
    if (!id) { try { localStorage.removeItem(CHALLENGE_KEY); } catch (e) {} return; }
    const r = loadRecords();
    for (const k of Object.keys(r)) if (k.startsWith(`${id}@`)) delete r[k];
    try { localStorage.setItem(CHALLENGE_KEY, JSON.stringify(r)); } catch (e) {}
  }
  function loadRecords() {
    try { return JSON.parse(localStorage.getItem(CHALLENGE_KEY)) || {}; } catch (e) { return {}; }
  }
  // The records are per challenge and species: "<name>:<n>@<species id>" -> { done, best }
  const recordKey = (id, species) => `${id}@${species}`;
  // Challenge id completed with species in sec seconds: its record, saved (not with ?all=1); longest:
  // the best is the longest time ([need:survive]), else the shortest
  function saveRecord(id, sec, species, longest = false) {
    const r = loadRecords(), key = recordKey(id, species), o = r[key] || {}, first = !o.done;
    o.done = true;
    if (!Number.isFinite(o.best) || (longest ? sec > o.best : sec < o.best)) o.best = sec;
    r[key] = o;
    if (!ALL) {
      try { localStorage.setItem(CHALLENGE_KEY, JSON.stringify(r)); }
      catch (e) { console.warn('Cognoship: the challenge record could not be saved', e); }
    }
    return { ...o, first };   // (first: not completed with this species before)
  }
  function parseChallenges() {
    challenge.list = [];
    SCRIPT.forEach((L, i) => {
      const a = L.actions.find(a => a.name === 'challenge');
      if (!a) return;
      const [name = 'Challenge', n = '', sub = ''] = a.args;
      // id (the record's key): "<name>:<n>" (not the subtitle, which may change); title: "<name>: <subtitle>" (none: "<name>: <n>")
      challenge.list.push({ id: `${name}:${n}`, name, stage: n, sub, title: sub ? `${name}: ${sub}` : n ? `${name}: ${n}` : name,
                            start: i, end: SCRIPT.length, allow: ['all'] });
    });
    challenge.list.forEach((c, k) => {
      if (k + 1 < challenge.list.length) c.end = challenge.list[k + 1].start;
      for (let i = c.start; i < c.end; i++) {
        const a = SCRIPT[i].actions.find(a => a.name === 'allow_species');
        if (a) { c.allow = a.args.flatMap(x => x.split(',')).map(x => x.trim()).filter(Boolean); break; }
      }
      // its [hide_morph] tags, in order (hide: their args), for its species book's thumbnails
      c.hide = [];
      for (let i = c.start; i < c.end; i++) for (const a of SCRIPT[i].actions) if (a.name === 'hide_morph') c.hide.push(a.args);
    });
  }
  // The hidden morphs while challenge c is played: the set-up's (baseHidden), then its own [hide_morph] tags
  function challengeHidden(c) {
    const saved = new Set(hiddenMorphs);
    hiddenMorphs.clear(); for (const k of baseHidden) hiddenMorphs.add(k);
    for (const args of c.hide) hideMorphs(...args);
    const out = new Set(hiddenMorphs);
    hiddenMorphs.clear(); for (const k of saved) hiddenMorphs.add(k);
    return out;
  }
  // The species challenge c allows, of the visible ones ([add_allow_species]), in their order
  const challengeSpecies = c => c.allow.includes('all') ? visSpecies : visSpecies.filter(e => c.allow.includes(e.id));
  // "12.3s" under a minute, else "1m 12s" (whole seconds)
  const fmtTime = sec => { if (Math.round(sec * 10) < 600) return `${sec.toFixed(1)}s`; const t = Math.round(sec); return `${Math.floor(t / 60)}m ${t % 60}s`; };
  // a new best time just made ({ id, species }): the challenge's species book, opening, reveals that species' slot
  let newRecord = null;
  // c: open on its species book instead (after completing it: its records, another species to try; Back: the challenge book)
  function openChallengeBook(c = null) {
    if (MODE !== 'challenge') return;
    if (book.open && book.mode !== 'challenge') closeBook();
    challenge.pending = c;
    book.mode = c ? 'species' : 'challenge';
    showBook();
  }
  // A slot per challenge: a green tick if completed, its title, then its best time and the species
  // it was made with
  // A slot per challenge, its title, after a green tick once completed with every species it allows
  // (completion and best times per species: in the species book)
  function renderChallenges() {
    setTitle([['Choose', 'choose'], ' a challenge']);
    const rec = loadRecords();
    $('book-grid').replaceChildren(...challenge.list.map((c, slot) => {
      const el = document.createElement('div'), label = document.createElement('span');
      el.className = 'book-slot pick challenge' + (c === (challenge.cur || challenge.last) ? ' current' : '');   // the one being played, or last played: a lighter background
      label.className = 'label';
      label.textContent = c.title;
      const sp = challengeSpecies(c);
      if (sp.length && sp.every(e => rec[recordKey(c.id, e.id)]?.done)) {
        const t = document.createElement('span'); t.className = 'tick'; t.textContent = '✓'; t.title = 'Completed with every species'; label.prepend(t);
      }
      el.append(label);
      el.title = `Play ${c.title}`;
      el.onclick = () => pickChallenge(c);
      if (SLOT_KEYS.species[slot]) el.dataset.key = SLOT_KEYS.species[slot];
      return el;
    }));
    fitBook();
  }
  // A challenge picked: its species next (the species book turns to the ones it allows, each with its
  // completion and best time for this challenge; even only one)
  function pickChallenge(c) {
    const sp = challengeSpecies(c);
    if (!sp.length) return console.warn(`challenge.txt: no visible species allowed for ${c.title}`);
    challenge.pending = c;
    book.mode = 'species';
    renderBook();
  }
  // Start challenge c with species id: as a species pick (the views fade, it is set up unseen), then
  // its lines from its [challenge] line (whose card fades them back in)
  // (morph: a retry's, for its first [put_player]; null: the challenge's own)
  function startChallenge(c, id, morph = null) {
    challenge.pending = null;
    closeBook();
    if (clone.held) endClone();
    if (speciesFading) return;
    speciesFading = true;
    if (chapterOpen) closeChapter();
    fadeWorld(1);
    setTimeout(() => {
      chosenSpecies = id;
      chooseSpecies(id);
      retake = true; cand = null;
      challenge.cur = challenge.last = c;
      freshAttempt();
      challengeMsgs = { success: null, fail: {} };   // (each challenge its own: the defaults, then its header)
      runChallengeDefaults();
      useMorph = morph;
      speciesFading = false;
      C.setRunning(true);
      run(c.start);
    }, WORLD_FADE_MS);
  }
  // [on_challenge_default]'s lines (actions only), in order and at once
  let challengeDefaults = [];
  function runChallengeDefaults() {
    for (const L of challengeDefaults) for (const a of L.actions) ACTIONS[a.name](...a.args);
  }
  // A new attempt: no failure, no time limit, mask or bullet cost (until the script sets them), the energy
  // from 0, the time from now
  function freshAttempt() {
    hiddenMorphs.clear(); for (const k of baseHidden) hiddenMorphs.add(k);   // (the hidden morphs as the set-up left them: a challenge's [hide_morph] is its own)
    Object.assign(challenge, { failed: null, deadline: null, maskCost: 0, eraseCost: 0, drawCost: 0, worldDrawCost: 0, longest: false, t0: C.simSteps });
    bulletCost = 0; cloneCost = 0;   // ([bullet_cost], [clone_cost] too: each challenge sets its own costs)
    energyHidden = false;   // (and shows the energy, unless [energy:off] again)
    energy = 0; showEnergy();
  }
  // Retry (the retry book's pick, Respawn after a failure): the challenge being played (or the last
  // one, if it was completed) again, from the line after its [challenge] (no card), behind a fade;
  // [put_player] puts the player's morph (the last chosen) in place of its own (useMorph)
  function retryChallenge() {
    const c = challenge.cur || challenge.last;
    if (!c) return;
    if (!challenge.cur) return startChallenge(c, chosenSpecies, playerMorph);   // completed: from its card
    if (clone.held) endClone();
    if (speciesFading) return;
    speciesFading = true;
    fadeWorld(1);
    setTimeout(() => {
      Dots.clear(); Bullets.clear(); Hostile.clear();
      freshAttempt();
      challengeMsgs = { success: null, fail: {} };   // (each challenge its own: the defaults, then its header)
      runChallengeDefaults();
      speciesFading = false;
      useMorph = playerMorph;
      retake = true; cand = null;   // the morph's picture may be retaken, as at a challenge's start
      if (c.start + 1 < c.end) run(c.start + 1); else completeChallenge();
      sceneHidden = false;
      fadeWorld(0);
      C.setRunning(true);
    }, WORLD_FADE_MS);
  }
  // The challenge being played fails (why: see challenge.failed): "Failed" in the ship view, with
  // Respawn and Choose (lostShown); the script stops waiting
  const FAIL_TEXT = { lost: 'the ship dissolved', time: "time's up", energy: 'out of energy', short: 'the ship dissolved too soon' };
  // [success_msg], [fail_msg]: the challenge's own messages (null: the defaults), set afresh for each
  // attempt (before [on_challenge_default] and the challenge's lines run again); {title}, {time},
  // {best} filled in (fillMsg)
  let challengeMsgs = { success: null, fail: {} };
  const fillMsg = (tpl, v) => tpl.replace(/\{(title|time|best)\}/g, (_, k) => v[k] ?? '').replace(/\s+$/, '');
  function failChallenge(why) {
    if (!challenge.cur || challenge.failed) return;
    challenge.failed = why;
    tests = []; progress = null; deferred = -1;
    const tpl = challengeMsgs.fail[why] ?? `_{title}_ failed: ${FAIL_TEXT[why]}.`;
    say(fillMsg(tpl, { title: challenge.cur.title, time: fmtTime((C.simSteps - challenge.t0) / C.STEPS_PER_S) }));
    syncButtons();
  }
  // Each frame, during a challenge: a lost ship fails it ([dissolve:retry]); [time_limit] runs out; [mask_cost] spends
  // the energy while the player presses (mask or erase), and fails it at 0
  let costSteps = 0;   // C.simSteps when the mask cost was last taken
  let drawCostAt = 0;   // performance.now() when the world-draw cost was last taken
  function challengeFrame() {
    const c = challenge.cur, now = C.simSteps;
    const dn = now - costSteps; costSteps = now;
    const t = performance.now(), dtS = drawCostAt ? Math.min(0.25, (t - drawCostAt) / 1000) : 0; drawCostAt = t;
    showTimer();
    if (!c || challenge.failed || chapterShowing() || speciesFading) return;
    const surv = tests.find(t => t.survive)?.survive;   // [need:survive]: the ship lost is the goal, if it lasted long enough
    if (surv && shipLost() && !emptySpace && surv.lasted() < surv.need) return failChallenge('short');
    if (challenge.dissolveFails && shipLost() && !emptySpace && !surv) return failChallenge('lost');
    if (challenge.deadline !== null && now >= challenge.deadline) return failChallenge('time');
    // the press's cost by its action: mask, erase or draw ([mask_cost], [erase_cost], [draw_cost])
    const costs = [challenge.maskCost, challenge.eraseCost, challenge.drawCost];
    if (costs.some(x => x > 0)) {
      const cost = costs[C.action] || 0;
      if (C.masking && dn > 0 && cost > 0) changeEnergy(-cost * dn / C.STEPS_PER_S);
      if (energy <= 0 && !surv) return failChallenge('energy');   // ([need:survive]: play on at 0)
    }
    // [world_draw_cost]: a stroke of Nebula / Antimatter's paintbrush or soft brush held (not the
    // eraser), per real second (drawing works while paused too)
    if (challenge.worldDrawCost > 0) {
      if (fieldMode() && Field.active && Field.tool !== 2 && dtS > 0) changeEnergy(-challenge.worldDrawCost * dtS);
      if (energy <= 0 && !surv) return failChallenge('energy');   // ([need:survive]: play on at 0)
    }
  }
  // [mask_cost], [erase_cost], [draw_cost]: energy per game second of a press with that action
  function setPressCost(what, v) {
    const n = v === 'off' ? 0 : parseFloat(v);
    if (!(n >= 0)) return console.warn(`script.txt: bad value [${what.replace('worldDraw', 'world_draw')}_cost:${v}], ignored`);
    challenge[what + 'Cost'] = n;
    showEnergy();
  }
  // the time left ([time_limit]) at the top left of the world view: red in the last 5 s
  function showTimer() {
    const on = !!challenge.cur && challenge.deadline !== null;
    $('timer').hidden = !on;
    if (!on) return;
    const left = Math.max(0, (challenge.deadline - C.simSteps) / C.STEPS_PER_S);
    $('timer').textContent = '⧖ ' + left.toFixed(1);
    $('timer').classList.toggle('low', left <= 5);
  }
  // A ship put back after dissolving ([dissolve:respawn], the challenge not failed): its completion
  // time counts from now
  function respawned() { if (challenge.cur && !challenge.failed) challenge.t0 = C.simSteps; }
  // The end of the challenge being played: recorded, said, and the book opens again after the SUCCESS box
  function completeChallenge() {
    const c = challenge.cur, sec = (C.simSteps - challenge.t0) / C.STEPS_PER_S;
    challenge.cur = null; line = -1; tests = []; progress = null;
    challenge.deadline = null; challenge.maskCost = challenge.eraseCost = challenge.drawCost = challenge.worldDrawCost = 0;
    showSuccess(true);   // (challenge mode: only here, at its end)
    const longest = challenge.longest;   // ([need:survive]: the longer the better)
    const r = saveRecord(c.id, sec, chosenSpecies, longest);
    newRecord = r.best === sec ? { id: c.id, species: chosenSpecies } : null;   // its species book reveals the slot
    const best = r.best !== sec ? `(best ${fmtTime(r.best)})` : r.first ? '(new completion!)' : '(new record!)';
    const tpl = challengeMsgs.success ?? `_{title}_ ${longest ? 'survived for' : 'completed in'} {time} {best}`;
    say(fillMsg(tpl, { title: c.title, time: fmtTime(sec), best }));
    setTimeout(() => { if (!challenge.cur && !book.open && !speciesFading) openChallengeBook(c); }, SUCCESS_MS);   // its species book (the record just made; Back: the challenge book)
  }
  // the line after a challenge's last, or the script's end
  const scriptEnd = () => challenge.cur ? challenge.cur.end : SCRIPT.length;

  // ---- target area ([target], [need:stay_target]) ----------------------------------------------
  // A square the size of HAESRE antigyrans (its pattern's larger side at the game's R), outlined in
  // dotted orange in both views (core opt.marks). The ship is in it while its CoM is.
  const TARGET_COLOR = '#ff9a1f', TARGET_DASH = [2, 3];
  let target = null;   // { col, row, size } (board cells)
  let targetSize = 0;
  function targetSide() {
    if (targetSize) return targetSize;
    const p = AQUARIUM233.patterns.find(p => p.id === 'HAESRE-antigyrans');
    const g = p ? Rule.rle2cells(p.cells[0]) : null;
    return (targetSize = g ? Math.round(Math.max(g.length, g[0].length) * C.R / p.R) : 3 * C.R);
  }
  // Place it by how: center (the board's centre), random or far (a spot clear of every creature;
  // far: of 30 such, the one farthest from the ship), always wholly inside the board (not across its edge)
  // (the scrolling world: center is the view's middle; extended a world point out of the view, the
  // target then kept by its world point, wx, wy: scrollFrame puts it on the board, or off: out of reach)
  // side: its side, cells (default targetSide())
  function placeTarget(how, dist, side) {
    const size = side || targetSide(), W = C.BW, H = C.BH;
    if (how === 'extended') {
      if (scroll) { const [wx, wy] = extendedPoint(size, dist); return (target = { wx, wy, size, col: 0, row: 0, off: true }); }
      how = 'far';
    }
    if (how === 'center') { const [col, row] = centerPt(); return (target = withWorld({ col, row, size })); }
    const comps = C.found ? C.found.comps : [], ship = how === 'far' ? C.ship() : null;
    let best = null, bestD = -1, found = 0;
    for (let k = 0; k < 400 && found < (ship ? 30 : 1); k++) {
      const col = size / 2 + Math.random() * (W - size), row = size / 2 + Math.random() * (H - size);   // wholly inside the board, not across its edge
      if (!spotOk(col, row, size)) continue;
      if (comps.some(c => Math.hypot(C.wrapD(c.col - col, W), C.wrapD(c.row - row, H)) < (c.radius || C.R) + size * 0.75)) continue;
      found++;
      const d = ship ? Math.hypot(C.wrapD(ship.col - col, W), C.wrapD(ship.row - row, H)) : 0;
      if (d > bestD) { best = { col, row }; bestD = d; }
    }
    if (!best) { console.warn('challenge: no clear spot for [target], placed at random'); best = { col: size / 2 + Math.random() * (W - size), row: size / 2 + Math.random() * (H - size) }; }
    return (target = withWorld({ ...best, size }));
  }
  const withWorld = t => { if (scroll) [t.wx, t.wy] = toWorld(t.col, t.row); return t; };
  function inTarget() {
    const sh = C.ship();
    if (!target || target.off || !sh || shipLost()) return false;
    const h = target.size / 2;
    return Math.abs(C.wrapD(sh.col - target.col, C.BW)) <= h && Math.abs(C.wrapD(sh.row - target.row, C.BH)) <= h;
  }
  C.opt.marks = () => (target && !target.off ? [{ col: target.col, row: target.row, w: target.size, h: target.size, color: TARGET_COLOR, dash: TARGET_DASH }] : []).concat(botMarks());

  function openBulletBook() {
    if (!Bullets.on) return;
    if (book.open && book.mode === 'morph') closeBook();
    book.mode = 'bullet';
    showBook();
  }
  // Story and sandbox's retry book (Retry): the morph book, a pick redoes the scene with it (pickMorph)
  function openRetryBook() {
    if (chapterShowing() || bookShowing()) return;
    book.mode = 'morph'; book.retry = true;
    showBook();
  }
  // retry: (challenge mode) the retry book, even with cloning on (Retry)
  function openBook(newName = null, retry = false) {
    if (chapterShowing()) return;   // never with a chapter card (coming, up or fading out)
    // challenge mode, during a challenge or after one: the retry book; not while a ship lost
    // without failing ([dissolve:respawn]) waits for a respawn, nor (but from Retry) while the
    // challenge allows cloning ([clone:on]) and the ship is alive: the book to clone from
    const cloning = clone.on && challenge.cur && !challenge.failed && !shipLost();
    if (MODE === 'challenge' && !newName && (challenge.cur || challenge.last) && (challenge.failed || !challenge.cur || !shipLost()) && (retry || !cloning)) {
      if (book.open && book.mode !== 'morph') closeBook();
      book.mode = 'morph'; book.retry = true;
      return showBook();
    }
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
    book.open = false; book.retry = false; newRecord = null;   // (revealed once)
    if (book.mode === 'bullet') { newBullets.clear(); revealing.clear(); }   // seen: no longer new (the button stays "Create" while more can be created)
    zap(false);
    if (book.fromDiscovery && book.newName) book.found.push(book.newName);
    book.fromDiscovery = false; book.newName = null;
    if (book.wasRunning) C.setRunning(true);
    if (book.mode === 'settings' || book.mode === 'hostile') settingsClosed();
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
  // the book buttons (and their keys, X S C V) flash yellow as they open their book, as the other controls
  $('morphs').onclick = () => { flash('morphs'); openBook(); };   // (holding a clone: the book to pick another one to clone)
  $('bullets').onclick = () => { flash('bullets'); openBulletBook(); };
  $('species').onclick = () => { flash('species'); MODE === 'challenge' ? openChallengeBook() : openSpeciesBook(); };
  if (MODE === 'challenge') { $('species').textContent = 'Challenges'; $('species').title = 'The challenge book: choose a challenge (D)'; }
  $('book-species').onclick = () => { if (book.mode === 'morph') openSpeciesBook(); };
  // The species book for a challenge goes back to the challenge book (its Back button, Escape); any other closes
  const bookBack = () => (book.mode === 'species' && !!challenge.pending) || book.mode === 'hostile';
  function closeOrBack() {
    if (!bookBack()) return closeBook();
    if (book.mode === 'hostile') { book.mode = 'settings'; return renderBook(); }   // (the hostile pulse book: back to the settings book)
    challenge.pending = null;
    book.mode = 'challenge';
    renderBook();
  }
  $('book-close').onclick = closeOrBack;
  $('book').addEventListener('pointerdown', e => { if (e.target === $('book')) closeBook(); });   // outside the box
  // Each frame: has the followed ship been a visible morph not discovered yet, without a break, for
  // DISCOVER_STEPS (a game second, in steps: paused time does not count)? A passing detection is
  // not a discovery. A discovered one gets a new picture (after SNAPSHOT_STEPS) only once after each
  // retake: a Respawn, a book pick (Choose), a species chosen (the sandbox's start too), a challenge
  // started or retried, a sandbox Retry (not boarding a clone: the old ship is still beside it).
  const DISCOVER_STEPS = C.STEPS_PER_S, SNAPSHOT_STEPS = C.STEPS_PER_S;
  let cand = null, lastStep = 0;   // cand: { name, n: steps } the morph the ship is in now, if it would be new
  let retake = false;   // set by a Respawn, a book pick, a species chosen, a challenge (re)started, a sandbox Retry: the next known morph held gets a new picture
  function checkDiscovery() {
    const dn = C.simSteps - lastStep;
    lastStep = C.simSteps;
    if (chapterShowing()) return;   // no discovery (nor its book) under a chapter card; it waits
    const d = C.focus && C.shipDet();   // (a twin ship: its whole's morph)
    const known = d && d.name && discovered.find(x => x.name === d.name);
    if (known && !Number.isFinite(known.hue) && Number.isFinite(d.hue)) { known.hue = d.hue; saveInventory(); }   // put in unseen: its colour now
    // a new discovery; or a known morph, whose picture is retaken quietly after SNAPSHOT_STEPS (once
    // the detection has settled), replacing the old one: only once after a Respawn or a book pick
    if (!book.on || !d || !d.name || !isVisible(d.name) || isHidden(d.name)) { cand = null; return; }   // (a hidden morph: never discovered)
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
  // Morphs renamed since, in versions (species -> old name -> new): an inventory keeps them under
  // the new names. Each version applied once (the inventory's rv: the last one applied; none: 0),
  // its renames all at once (so a name can pass to another morph: VT049W's rectogyrans became
  // sequerectus as paralimus became rectogyrans); a name already taken is left as it is.
  const MORPH_RENAMES = [
    { KZ2BLD: { gyrans: 'limus' } },                                   // rv 1
    { VT049W: { rectogyrans: 'sequerectus', paralimus: 'rectogyrans' } },   // rv 2 (2026-11-18)
  ];
  function renameMorphs(morphs, rv = 0) {
    MORPH_RENAMES.slice(rv).forEach(step => {
      for (const [sp, ren] of Object.entries(step)) {
        const list = morphs[sp];
        if (!Array.isArray(list)) continue;
        const was = list.map(d => d && d.name), now = was.map(n => ren[n] ?? n);
        list.forEach((d, i) => { if (d && ren[was[i]] && !now.some((n, j) => j !== i && n === now[i])) d.name = now[i]; });
      }
    });
    return morphs;
  }
  function loadInventory() {
    try {
      const o = JSON.parse(localStorage.getItem(INV_KEY));
      if (o && typeof o === 'object') return { morphs: o.morphs && typeof o.morphs === 'object' ? renameMorphs(o.morphs, Number.isInteger(o.rv) ? o.rv : 0) : {},
        bullets: Array.isArray(o.bullets) ? o.bullets : [], energy: Number.isFinite(o.energy) ? Math.max(0, o.energy) : 0, rv: MORPH_RENAMES.length };
    } catch (e) {}
    return { morphs: {}, bullets: [], energy: 0, rv: MORPH_RENAMES.length };
  }
  // The energy: what the ship has eaten (dots.js counts it: its score's rises), less what bullet
  // patterns cost. energyUnsaved (sandbox only): its change since it was last saved, added to the
  // central energy (so two tabs add up).
  let energy = 0, dotsSeen = 0, energyUnsaved = 0;
  let bulletCost = 0;   // [bullet_cost:<n>]: energy per shot of the player's gun
  let cloneCost = 0;    // [clone_cost:<n>]: energy per clone pasted
  let energyHidden = false;   // [energy:off]: not shown
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
      if (book.mode === 'challenge' || (book.mode === 'species' && challenge.pending)) {   // the ticks and the best times go
        clearCompletion(book.mode === 'species' ? challenge.pending.id : null);
        $('book-clear').classList.remove('holding'); $('book-clear').classList.add('done');
        $('book-clear').textContent = 'Completion & best time cleared';
        renderBook();
        setTimeout(() => { $('book-clear').classList.remove('done'); if (book.open) renderBook(); }, 1500);
        return;
      }
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
    if (inv.bullets.length) {   // the gun starts with the list's last (in the book's order)
      const ids = Object.keys(Bullets.KINDS), names = inv.bullets.filter(b => ids.includes(b));
      Bullets.setList(inv.bullets.join(','), names.sort((a, b) => ids.indexOf(a) - ids.indexOf(b)).pop());
    }
    if (MODE === 'sandbox') { energy = inv.energy; showEnergy(); }   // (not a gain) (challenge mode: its own, from 0)
  }
  // A picture of the ship view as it is now (as drawn this frame, so called from the frame hook,
  // right after the render), scaled down to the world view's scale (C.WS px per cell): a data URL.
  // Whatever the zoom: the SNAP_CELLS-wide square the ship view shows at its middle zoom (6 px per
  // cell), taken from the world view (so unrotated), round the board's edges, at its scale.
  const SNAP_CELLS = C.ZOOM.edge / 6;   // 64 cells
  function snapshot() {
    try {
      const f = C.ship(), WS = C.WS, n = Math.round(SNAP_CELLS * WS), W = C.BW * WS, H = C.BH * WS;
      if (!f) return null;
      const cv = document.createElement('canvas');
      cv.width = cv.height = n;
      const g = cv.getContext('2d');
      return SimGL.withoutField(() => {
        const x0 = Math.round(C.wxp(f.col) - n / 2), y0 = Math.round(C.wyp(f.row) - n / 2);   // (the world view's place: shifted in a scrolling world)
        for (const dx of [-W, 0, W]) for (const dy of [-H, 0, H]) {   // the world view's copies round the torus
          const sx = x0 + dx, sy = y0 + dy;
          const ax = Math.max(0, sx), ay = Math.max(0, sy), bx = Math.min(W, sx + n), by = Math.min(H, sy + n);
          if (bx > ax && by > ay) g.drawImage($('gl'), ax, ay, bx - ax, by - ay, ax - sx, ay - sy, bx - ax, by - ay);
        }
        return cv.toDataURL();   // (taken with the world-draw field left out: SimGL.withoutField)
      });
    } catch (e) { return null; }
  }

  // The player's morph, set by [put_player] and a pick in the morph book: what Respawn puts back.
  let playerMorph = null;
  // Put morph <id> in by `how` (center, random, far); follow it if focusNew. False if unknown.
  function putMorph(token, id, how, focusNew) {
    const p = findMorph(id);
    if (!p) { console.warn(`script.txt: unknown morph [${token}:${id}]`); return false; }
    if (typeof how === 'string' && !PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [${token}:${id}:${how}], using center`); how = 'center'; }
    C.selectMorph(p.id); C.putPattern(creatureHow(how), focusNew);
    return true;
  }
  // A morph of the current species by its name (limus, gyrans, ...), or a pattern id (5N7KKM-limus).
  const findMorph = id => PATTERNS.find(p => C.speciesOf(p) === C.speciesOf(C.pattern) &&
                                             C.morphName(p).replace(C.BOOM, '') === id) || PATTERNS.find(p => p.id === id);
  // Whether name is a morph of the current species (one of its patterns)
  const inSpecies = name => !!C.pattern && PATTERNS.some(p => C.speciesOf(p) === C.speciesOf(C.pattern) && C.morphName(p).replace(C.BOOM, '') === name);
  // The current species' first morph (its first pattern, as listed), by name
  function firstMorph() {
    const p = C.pattern && PATTERNS.find(p => C.speciesOf(p) === C.speciesOf(C.pattern) && !p.hidden);
    return p ? C.morphName(p).replace(C.BOOM, '') : null;
  }
  // The visible species ([add_allow_species]), in order: { name, id, morphs: its visible morphs }.
  // For the species chooser to come; nothing reads it yet ([species] loads any species).
  const visSpecies = [];
  const speciesBlocks = new Map();   // [on_start[:<id>]] blocks: id ('' = any species) -> its lines
  // A species chosen (the species book, to come; [choose_species:<id>]): load it ([species], with
  // its [add_allow_species] name), make its morphs visible ([allow_morph], its [add_allow_species]
  // list), then run the actions of the general [on_start] block, then of its own [on_start:<id>]
  // block (which may change what the general one set), in order and at once (no waiting: a block
  // holds actions only).
  function chooseSpecies(id) {
    const e = visSpecies.find(x => x.id === id);
    if (!e) return console.warn(`script.txt: [choose_species:${id}]: not added by [add_allow_species]`);
    ACTIONS.species(e.name, id);
    ACTIONS.allow_morph(e.morphs.join(','));
    for (const L of (speciesBlocks.get('') ?? []).concat(speciesBlocks.get(id) ?? []))
      for (const a of L.actions) ACTIONS[a.name](...a.args);
    if (settingsOn) applyObjects();   // the Objects slider's [min_morph:<n>], over the block's
    retake = true; cand = null;   // a fresh world: the morph's picture may be retaken (the sandbox's start too)
  }
  // ---- settings ([settings:on]) -------------------------------------------------------------
  // The Settings button (no key; the world view's lower right) opens the settings book: World, Autospawn,
  // Objects and Hostile, each kept in the browser across runs.
  function openSettingsBook() {
    if (!settingsOn) return;
    if (book.open && book.mode !== 'settings') closeBook();
    book.mode = 'settings'; worldWas = worldMode;
    showBook();
  }
  $('settings').onclick = () => { flash('settings'); openSettingsBook(); };
  // World, a button cycling three modes. Normal: the world as the script makes it. Centered: the same
  // (ordinary) world, the world view rolled round the torus so the ship stays in its middle (roll,
  // centerFrame), at once. Scroll: every [clear_world] makes a scrolling world, as
  // [clear_world:scroll]; changed to or from it in the book: as it closes, the species is chosen
  // again (the world afresh, with it). Centered and Scroll: the stars follow the ship ([stars:follow:0.1]).
  const WORLD_KEY = 'cognoship-world', WORLDS = ['normal', 'centered', 'scroll'];
  let worldMode = 'normal', worldWas = 'normal', scrollStars = false;   // scrollStars: the setting's [stars:follow] is in force
  try {
    const v = localStorage.getItem(WORLD_KEY);
    if (WORLDS.includes(v)) worldMode = v;
    else if (localStorage.getItem('cognoship-scroll') === '1') worldMode = 'scroll';   // (the old Scroll toggle's)
  } catch (e) {}
  $('world').onclick = () => {
    worldMode = WORLDS[(WORLDS.indexOf(worldMode) + 1) % WORLDS.length];
    try { localStorage.setItem(WORLD_KEY, worldMode); } catch (e) {}
    if (!scroll) { setRoll(settingsOn && worldMode === 'centered'); followStars(); }   // (Centered on or off: at once)
    syncObjects();
  };
  // the setting's star drift: on with Centered or Scroll, back off (as the script set them) without
  function followStars() {
    const on = settingsOn && worldMode !== 'normal';
    if (on) { C.setStarsFollow(0.1); scrollStars = true; }   // ([stars:follow:0.1], the brightness and twinkle kept)
    else if (scrollStars) { C.setStarsFollow(0); scrollStars = false; }
  }
  function settingsClosed() {
    const was = worldWas; worldWas = worldMode;
    if ((was === 'scroll') === (worldMode === 'scroll')) return;
    const id = chosenSpecies ?? (C.pattern && C.speciesOf(C.pattern));
    if (id && visSpecies.some(e => e.id === id)) pickSpecies(id);
  }
  // The Objects slider: [min_morph:<n>], n one of OBJECTS (0: off), run after each
  // [on_start] block (and at once, on a change), so it overrides the block's count. Only
  // the count: the morphs and the placing stay the block's (the last [min_morph]'s; none yet: all,
  // random). Kept in the browser (OBJ_KEY) across runs; until first set, 2.
  const OBJECTS = [0, 1, 2, 5, 10, 20, 40], OBJ_KEY = 'cognoship-objects';
  let settingsOn = false, minSpec = { ids: [], how: 'random' };   // minSpec: the last [min_morph]'s morphs and placing
  let objectsIdx = 2;
  try { const v = parseInt(localStorage.getItem(OBJ_KEY), 10); if (v >= 0 && v < OBJECTS.length) objectsIdx = v; } catch (e) {}
  function applyObjects() {
    const n = OBJECTS[objectsIdx];
    if (!n) minMorph = null;
    else if (minMorph) minMorph.n = n;
    else { minMorph = { n, ids: minSpec.ids, how: minSpec.how, next: 0 }; edgeQueue = []; }
    applyHostiles();
  }
  // The Hostile slider: [min_hostile:<n>] as Objects is [min_morph:<n>], the same counts (its morphs
  // and placing the last [min_hostile]'s; none yet: all, random). Until first set, off. They fire the
  // Hostile Pulses button's pattern (beside it: the hostile pulse book), whatever the script's. Both
  // kept in the browser.
  const HOSTILES = OBJECTS, HOST_KEY = 'cognoship-hostile-n', HOST_BULLET_KEY = 'cognoship-hostile-pulse';
  let hostSpec = { ids: [], how: 'random', bullet: 'Scutium' }, hostilesIdx = 0, hostBullet = 'Scutium serratus';
  try { const v = parseInt(localStorage.getItem(HOST_KEY), 10); if (v >= 0 && v < HOSTILES.length) hostilesIdx = v; } catch (e) {}
  try { const v = localStorage.getItem(HOST_BULLET_KEY); if (v && Bullets.kindOf(v)) hostBullet = Bullets.kindOf(v); } catch (e) {}
  function applyHostiles() {
    const n = HOSTILES[hostilesIdx];
    if (!n) { minHostile = null; return; }
    Hostile.set(hostPulse());
    if (minHostile) minHostile.n = n;
    else minHostile = { n, ids: hostSpec.ids, how: hostSpec.how, next: 0, puts: 0 };
  }
  // The pattern the bots fire: the pick, if the player has it (created: Bullets.list), else the
  // player's first (none: the pick anyway)
  const hostPulse = () => Bullets.list.includes(hostBullet) ? hostBullet : Bullets.list[0] || hostBullet;
  // The hostile pulse book (Hostile Pulses, in the settings book): as the pulse book, the patterns the
  // player has (no world-draw types), the one in use yellow, the others' slots empty; a pick arms the
  // bots with it and goes back to the settings book (as does Back, Escape)
  $('hostile-pulses').onclick = () => { flash('hostile-pulses'); book.mode = 'hostile'; renderBook(); };
  function renderHostilePulses() {
    setTitle([['Choose', 'choose'], ' the hostile pulse type']);
    const have = Bullets.list, cur = hostPulse();
    $('book-grid').replaceChildren(...Object.keys(Bullets.KINDS).map((name, slot) => {
      const el = document.createElement('div'), label = document.createElement('span'), ok = have.includes(name);
      el.className = 'book-slot' + (ok ? ' pick' : ' empty') + (name === cur ? ' current' : '');
      label.className = 'label';
      if (ok) label.append(bulletPic(name), name); else label.textContent = '\u00a0';
      el.append(label);
      if (!ok) return el;
      el.title = `The hostile creatures fire ${name}`;
      el.onclick = () => {
        hostBullet = name;
        try { localStorage.setItem(HOST_BULLET_KEY, name); } catch (e) {}
        if (minHostile) Hostile.set(name);   // (the bots now: their bullets in flight cleared)
        book.mode = 'settings'; renderBook();
      };
      if (SLOT_KEYS.bullet[slot]) el.dataset.key = SLOT_KEYS.bullet[slot];
      return el;
    }));
  }
  // Autospawn (beside it): once the ship is lost (Respawn shows), Respawn by itself after
  // AUTOSPAWN_MS, as the button would (not while a book, a card or a fade is up). Kept as Objects.
  const AUTO_KEY = 'cognoship-autospawn', AUTOSPAWN_MS = 1000;
  let autospawn = false, lostAt = null;
  try { autospawn = localStorage.getItem(AUTO_KEY) === '1'; } catch (e) {}
  function autospawnFrame() {
    if (!settingsOn || !autospawn || !lostShown() || lostFails()) { lostAt = null; return; }
    if (book.open || chapterShowing() || speciesFading || resetting) return;
    const now = performance.now();
    if (lostAt === null) lostAt = now;
    if (now - lostAt >= AUTOSPAWN_MS) { lostAt = null; doReset(); }
  }
  // Autofire (after it): while on, with the gun on (not a world-draw type) and a ship, the ship fires
  // by itself every [bullet_wait] a to b game s (at random), the first c after it can start (Autofire
  // on, the ship back, a scene); at the nearest hostile object (hostileObjects), else the nearest
  // other one, the ship's excluded; none: no shot. Not under a book, a card or a fade. Kept as Objects.
  const AUTOFIRE_KEY = 'cognoship-autofire';
  let autofire = false, autoNext = null;
  try { autofire = localStorage.getItem(AUTOFIRE_KEY) === '1'; } catch (e) {}
  function autofireFrame() {
    const ship = C.ship();
    if (!settingsOn || !autofire || !Bullets.on || fieldMode() || !ship || shipLost() || speciesFading || resetting || challenge.failed) { autoNext = null; return; }
    if (book.open || chapterShowing() || !C.found) return;
    const now = C.simSteps;
    if (autoNext === null) { autoNext = now + Math.round(fireWait.first * C.STEPS_PER_S); return; }
    if (now < autoNext || !canPay() || !Bullets.canFire()) return;
    const dist = o => Math.hypot(C.wrapD(o.col - ship.col, C.BW), C.wrapD(o.row - ship.row, C.BH));
    const nearest = list => list.reduce((a, o) => !a || dist(o) < dist(a) ? o : a, null);
    const hostile = hostileObjects();
    const t = nearest(hostile) || nearest(C.objects().filter(o => !o.parts.includes(C.focusIdx)));
    if (!t) return;
    if (Bullets.fire(Math.atan2(C.wrapD(t.row - ship.row, C.BH), C.wrapD(t.col - ship.col, C.BW)))) {
      flash('fire'); if (bulletCost) changeEnergy(-bulletCost);
      autoNext = now + Math.round((fireWait.min + Math.random() * (fireWait.max - fireWait.min)) * C.STEPS_PER_S);
    }
  }
  // [bullet_wait]: Autofire's a to b game s between shots, and c before the first
  const FIRE_WAIT = { min: 1.5, max: 4, first: 2 };
  let fireWait = { ...FIRE_WAIT };
  $('autofire').onclick = () => {
    autofire = !autofire; autoNext = null;
    try { localStorage.setItem(AUTOFIRE_KEY, autofire ? '1' : '0'); } catch (e) {}
    syncObjects();
  };
  $('autospawn').onclick = () => {
    autospawn = !autospawn;
    try { localStorage.setItem(AUTO_KEY, autospawn ? '1' : '0'); } catch (e) {}
    syncObjects();
  };
  function syncObjects() {
    $('autospawn').setAttribute('aria-pressed', autospawn);
    $('autofire').setAttribute('aria-pressed', autofire);
    $('world').setAttribute('aria-pressed', worldMode !== 'normal');
    $('world').textContent = 'World: ' + worldMode[0].toUpperCase() + worldMode.slice(1);
    $('objects').value = objectsIdx;
    $('objects-v').textContent = OBJECTS[objectsIdx] || 'Off';
    fillRange($('objects'));
    $('hostiles').value = hostilesIdx;
    $('hostiles-v').textContent = HOSTILES[hostilesIdx] || 'Off';
    $('hostile-pulses').title = `Hostile Pulses: the pulse type the hostile creatures fire (now ${hostPulse()})`;
    fillRange($('hostiles'));
  }
  $('hostiles').oninput = e => {
    hostilesIdx = +e.target.value;
    try { localStorage.setItem(HOST_KEY, String(hostilesIdx)); } catch (err) {}
    syncObjects();
    applyHostiles();   // (the world now too)
  };
  $('objects').oninput = e => {
    objectsIdx = +e.target.value;
    try { localStorage.setItem(OBJ_KEY, String(objectsIdx)); } catch (err) {}
    syncObjects();
    applyObjects();   // (the world now too, not only the next species chosen)
  };
  // extended: out of the world view (the scrolling world's beacons and target; creatures, which
  // live on the board, and anything outside a scrolling world: as far)
  const PUT_HOW = ['center', 'random', 'far', 'extended'];
  // how as core putPattern takes it: extended as far (creatures live on the board); a spot [col, row] as it is
  const creatureHow = how => how === 'extended' ? 'far' : how;   // (or a placer: the scrolling world's edges)
  const ACTIONS = {
    // add_allow_species:<name>:<id>[:<ids>]: add species <id> (named <name>, as [species]) to the
    // visible species, with its visible morphs <ids> (as [allow_morph]; none given: all). Again:
    // its name and morphs updated.
    // on_start[:<id>] ... end (old name: on_choose_species): a block of lines, taken out of the script's run
    // (parseScript), whose actions chooseSpecies runs. No id: for any species; with one: for that
    // species only, instead of the general block.
    on_start: () => {},
    on_choose_species: () => {},   // (its old name, still read)
    // on_challenge_default ... end: (challenge mode) a block of lines (actions only), taken out of
    // the script's run (parseScript), run at the start of every challenge and of every retry, after
    // its fresh attempt (no costs, energy 0) and before its own lines: e.g. defaults that each
    // challenge may then change
    on_challenge_default: () => {},
    end: () => {},
    // choose_species:<id>: as if the player chose species <id> in the species book (chooseSpecies)
    // (redone by Retry, in its replay and in the current line run afresh: the species the player
    // chose since, if any)
    choose_species: (id) => chooseSpecies(resetting && chosenSpecies ? chosenSpecies : id),
    // challenge:<name>[:<n>]: (challenge mode) starts challenge <name>, stage <n>: a title card
    // ("Challenge · Stage <n>" above <name>), as [chapter]'s; the line waits until it closes.
    // Also a checkpoint, as [chapter]: Retry replays the challenge's actions from here.
    challenge: (name = '', n = '', sub = '') => { showChapter(sub || n, name, 'challenge'); return () => !chapterOpen; },
    // time_limit:<s>|off: (challenge mode) the challenge must be completed within s game seconds
    // from here (shown at the top left of the world view), or it fails; off (or none): no limit
    time_limit: (v) => {
      const t = parseFloat(v);
      if (v === undefined || v === 'off') { challenge.deadline = null; return; }
      if (!(t > 0)) return console.warn(`script.txt: bad value [time_limit:${v}], ignored`);
      challenge.deadline = C.simSteps + Math.round(t * C.STEPS_PER_S);
    },
    // mask_cost:<n>: (challenge mode) a press on the ship view (mask or erase) costs n energy per
    // game second it is held; the energy at 0 or below fails the challenge (so set some first, by
    // [energy:<n>]); 0 or off: free
    mask_cost: (v) => setPressCost('mask', v),
    // erase_cost:<n> / draw_cost:<n>: (challenge mode) as mask_cost, for a press that erases / draws
    // (each action its own cost; mask_cost is only the mask's)
    erase_cost: (v) => setPressCost('erase', v),
    draw_cost: (v) => setPressCost('draw', v),
    // world_draw_cost:<n>|off: (challenge mode) drawing in the world with Nebula or Antimatter
    // ([world_draw]: the paintbrush or the soft brush, a stroke held; not the eraser, Clear or Undo)
    // costs n energy per second; 0 or below fails the challenge, as mask_cost; 0 or off: free
    world_draw_cost: (v) => setPressCost('worldDraw', v),
    // bullet_cost:<n>|off: each shot of the player's gun costs n energy (the energy then shows at the
    // top); with less than n the gun doesn't fire (Fire dimmed); 0 or off: free. Kept until changed
    // (challenge mode: a challenge's start or retry sets it free again, as [mask_cost]).
    bullet_cost: (v) => {
      const n = v === 'off' ? 0 : parseFloat(v);
      if (!(n >= 0)) return console.warn(`script.txt: bad value [bullet_cost:${v}], ignored`);
      bulletCost = n;
      showEnergy();
    },
    // clone_cost:<n>|off: each clone pasted costs n energy (shown at the top); with less than n it
    // can't be pasted (it stays held); 0 or off: free. Challenge mode: reset with each attempt.
    clone_cost: (v) => {
      const n = v === 'off' ? 0 : parseFloat(v);
      if (!(n >= 0)) return console.warn(`script.txt: bad value [clone_cost:${v}], ignored`);
      cloneCost = n;
      showEnergy();
    },
    // dissolve:retry|respawn: (challenge mode) what the ship dissolving does: retry: it fails the
    // challenge ("Failed", then Respawn or Choose start it again); respawn (the default): as in the
    // other modes ("Dissolved", Respawn / Choose put a ship back, by [retry_action])
    dissolve: (v = 'respawn') => {
      if (!['retry', 'respawn'].includes(v)) return console.warn(`script.txt: unknown [dissolve:${v}], ignored`);
      challenge.dissolveFails = v === 'retry';
    },
    // energy:<n>: set the energy (the dots score) to n (a challenge starts from 0)
    // energy:off: the energy not shown (the dots and beacons eaten still add to it); [energy:<n>]
    // shows it again (challenge mode: a new attempt too)
    energy: (v) => {
      if (v === 'off') { energyHidden = true; showEnergy(); return; }
      if (v === 'on') { energyHidden = false; showEnergy(); return; }   // (shown, unchanged)
      const n = parseFloat(v);
      if (!Number.isFinite(n)) return console.warn(`script.txt: bad value [energy:${v}], ignored`);
      energyHidden = false;
      changeEnergy(n - energy);
    },
    // allow_species:all|<ids>: (challenge mode, on a challenge's lines) the species the player may
    // choose for it (parseChallenges reads it; nothing to do when run)
    allow_species: () => {},
    // target[:center|random|far]: put the target area (a dotted orange square, see placeTarget) at
    // the board's centre, at a random spot clear of creatures (the default), or at the clear spot
    // farthest from the ship; [clear_world] removes it
    // target[:<size>][:<how>[:<d>]]: size: its side, in R (default: HAESRE antigyrans's size);
    // extended:<d>: (a scrolling world) d R from the ship, at a random angle
    target: (...args) => {
      const size = /^\d+(\.\d+)?$/.test(args[0] ?? '') ? parseFloat(args.shift()) : undefined;
      let [how = 'random', dist] = args;
      if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [target:${how}], using random`); how = 'random'; }
      const d = parseFloat(dist);
      if (dist !== undefined && !(d > 0)) console.warn(`script.txt: bad distance [target:${how}:${dist}], ignored`);
      placeTarget(how, d > 0 ? d : undefined, size > 0 ? size * C.R : undefined);
    },
    add_allow_species: (name, id, ids) => {
      if (!id || !PATTERNS.some(p => C.speciesOf(p) === id)) return console.warn(`script.txt: unknown species [add_allow_species:${name}:${id}]`);
      const morphs = ids === undefined ? ['all'] : morphIds([ids]);
      const e = visSpecies.find(x => x.id === id);
      if (e) Object.assign(e, { name: name || e.name, morphs }); else visSpecies.push({ name, id, morphs });
      syncButtons();   // (sandbox: the Species button)
    },
    // species:<name>:<id>: load species <id> (its rule; the world is kept), with no visible
    // morphs; the morph book names it "Species: <name>" (genus and species, e.g. Aquarium
    // originale). ([species:<id>]: no name, the book shows the id.) Resets [hide_morph] (none hidden).
    species: (name, id) => {
      if (id === undefined) [name, id] = ['', name];
      const p = PATTERNS.find(p => C.speciesOf(p) === id);
      if (!p) return console.warn(`script.txt: unknown species [species:${name}:${id}]`);
      C.setPattern(p.id);
      speciesName = name;
      vis.all = false; vis.allow.clear();
      hiddenMorphs.clear();
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
    // hide_morph:<morphs>[:<code>]: these morphs (names, ',' between) stay detected and named (if
    // visible: [allow_morph], [add_allow_species]), but are left out of the morph book (no respawn,
    // retry or clone pick of them; e.g. too powerful) and the species book, and are never
    // discovered; for species <code> only, or with no code every species. none: none hidden (for
    // <code>, or with no code at all). Kept until changed.
    hide_morph: (morphs = '', code) => { hideMorphs(morphs, code); if (book.open) renderBook(); },
    // morph_book:on|off: show or hide the Morphs button (off also closes the book)
    morph_book: (state = 'on') => setBook(state !== 'off'),
    // clear_world[:<bg>][:empty][:scroll]: empty the world, with no dots (until [dot:<n>]), no bullets,
    // no [min_morph] and no [target]; empty: empty space, no player (emptySpace), until the next
    // [clear_world]; scroll: reserved for a scrolling world (not built: an ordinary one)
    clear_world: (...args) => {
      emptySpace = args.includes('empty');
      const setScroll = settingsOn && worldMode === 'scroll';   // the World setting's Scroll: a scrolling world
      startScroll(args.includes('scroll') || setScroll);   // a scrolling world, or the ordinary torus
      setRoll(!scroll && settingsOn && worldMode === 'centered');   // its Centered: the ordinary one, rolled to keep the ship in the middle
      followStars();   // (Centered, Scroll: the stars follow the ship, [stars:follow:0.1]; off again: they drift as the script set them)
      let bg = args.find(a => a !== 'empty' && a !== 'scroll') || 'black';
      target = null;
      if (!BACKGROUNDS.includes(bg)) { console.warn(`script.txt: unknown background [clear_world:${bg}], using black`); bg = 'black'; }
      background = bg;
      minMorph = null; minHostile = null;   // a new scene: no [min_morph] nor [min_hostile] until the script sets one again
      if (Dots.on) { Dots.set(0); Dots.clear(); }   // nor dots (none put back) until a [dot:<n>]; the energy kept
      Bullets.clear();   // nor bullets in flight
      bots = []; Hostile.clear();   // nor bots, nor their bullets
      Field.reset();   // nor anything drawn ([world_draw])
      C.clearWorld();
    },
    // world_draw:on|off|nebula|antimatter: put the bullet types Nebula and / or Antimatter (on: both;
    // ',' or ':' between) into the bullet book, quietly (after the patterns; the gun must be on),
    // or take both out (off, which also empties the fields; neither is in the book by default). In
    // use, the player paints the world's R G B fields instead of firing (see fieldMode).
    // [add_bullet:nebula|antimatter] puts one in with the book's reveal instead (the story).
    world_draw: (...args) => {
      const ts = args.flatMap(a => a.split(',')).map(t => t.trim().toLowerCase()).filter(Boolean);
      if (!ts.length) ts.push('on');
      if (ts.includes('off')) { Field.reset(); Bullets.setWorldDraw(null); return; }
      const names = ts.flatMap(t => t === 'on' ? ['nebula', 'antimatter'] : [t]);
      if (!Bullets.setWorldDraw(names)) console.warn(`script.txt: unknown option [world_draw:${args.join(':')}], ignored`);
    },
    // put_morph:<id>[:center|random|far]: put morph <id> of the current species into the world:
    // at the centre (default, unrotated), at a random empty spot and angle (the lab's Add), or at
    // a random angle in the empty spot farthest from the followed ship. Never followed: the ship
    // stays followed (if none, none is).
    // <id> written <morph>=<k> (e.g. lithos=3): a group of k, as [min_morph]'s (putGroup).
    // put_morph:<n>:<ids>[:<how>]: n creatures, each picked at random from <ids> (names, ',' between,
    // each may be <morph>=<k>; all: the visible list), as [min_morph] picks; center: the first only,
    // the rest random.
    put_morph: (...args) => {
      const n = /^\d+$/.test(args[0] ?? '') ? parseInt(args.shift(), 10) : 1;
      const [ids = '', how = 'center'] = args;
      const pool = ids === 'all' ? visibleList().map(id => ({ id, k: 1 })) : morphIds([ids]).map(x => groupArg('put_morph', x));
      if (!pool.length) return console.warn(`script.txt: [put_morph:${ids}]: no morph`);
      for (let i = 0; i < n; i++) {
        const g = pool[Math.floor(Math.random() * pool.length)], h = i > 0 && how === 'center' ? 'random' : how;
        if (g.k > 1) putGroup(g.id, g.k, h, 'put_morph');
        else putMorph('put_morph', g.id, h, 'never');
      }
    },
    // put_player:<id>[:center|random|far]: as put_morph, but the new creature is the player's
    // ship (followed), and <id> becomes the player's morph: what Respawn puts back.
    // put_player:first: the species' first morph (its first pattern).
    // <morph>=<k>: a group of k, the ship the one nearest the group's centre (the rest are not
    // the player's: Respawn puts back one)
    // <morphs> a list (',' between): the first of them that is a visible morph of the species
    // (else the first that is a morph of it, else the first, warned)
    put_player: (ids, how = 'center') => {
      const cands = String(ids ?? '').split(',').map(x => x.trim()).filter(Boolean).map(x => groupArg('put_player', x));
      if (!cands.length) return console.warn('script.txt: [put_player] with no morph, ignored');
      const name = e => e.id === 'first' ? firstMorph() ?? e.id : e.id;
      const g = cands.find(e => inSpecies(name(e)) && isVisible(name(e))) || cands.find(e => inSpecies(name(e))) || cands[0];
      if (cands.length > 1 && !(inSpecies(name(g)) && isVisible(name(g)))) console.warn(`script.txt: [put_player:${ids}]: none of these is a visible morph of the species; ${name(g)}`);
      let id = name(g);
      if (useMorph && inSpecies(useMorph)) id = useMorph;   // a retry: the player's morph (the group's size kept)
      useMorph = null;
      if (g.k > 1 ? putGroup(id, g.k, how, 'put_player', true) : putMorph('put_player', id, how, true)) playerMorph = id;
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
    // sub_chapter:<name>: (authoring) a place dev mode's I (devChapter) also stops at; nothing else
    // (<name> for reference only)
    sub_chapter: () => {},
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
    // An id written <morph>=<k> (e.g. lithos=3) puts a group of k of that morph instead: packed
    // together around the spot, each at a random angle, a thin gap between them (putGroup).
    // n = 0 or off: stop.
    min_morph: (n, ids = '', how = 'random') => {
      const k = parseInt(n, 10);
      if (n === 'off' || k === 0) { minMorph = null; return; }
      if (!(k > 0)) return console.warn(`script.txt: bad value [min_morph:${n}], ignored`);
      if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [min_morph:${n}:${ids}:${how}], using random`); how = 'random'; }
      const entries = ids === 'all' ? [] : morphIds([ids]).map(x => { const g = groupArg('min_morph', x); return { id: g.id, n: g.k }; });
      minMorph = { n: k, ids: entries, how, next: 0 };
      minSpec = { ids: entries, how };   // (for the Objects slider, [settings])
      edgeQueue = [];
    },
    // settings:on|off: the species book's Objects slider (sandbox mode), or not (the default); see applyObjects
    settings: (v = 'on') => {
      if (v !== 'on' && v !== 'off') return console.warn(`script.txt: unknown option [settings:${v}], ignored`);
      settingsOn = v === 'on';
      if (!settingsOn && !scroll) setRoll(false);
      if (!settingsOn) followStars();
      if (!settingsOn && book.open && book.mode === 'settings') closeBook();
      syncButtons();
    },
    // success_msg:<text>: (challenge mode) the message when this challenge is completed, in place of
    // "_{title}_ completed in {time} {best}" ("survived for", waiting on [need:survive]); {title}, {time}
    // (as 1m 12s) and {best} ("(best …)", "(new completion!)" or "(new record!)") are filled in. A ':'
    // in the text is kept (as ': '). No text: back to the default. Set afresh for each attempt.
    success_msg: (...text) => { challengeMsgs.success = text.join(': ') || null; },
    // fail_msg:<cause>:<text>: (challenge mode) the message when this challenge fails by <cause> (lost,
    // short, time, energy: see FAIL_TEXT), in place of "_{title}_ failed: <reason>."; {title} and {time}
    // (so far) filled in. No text: back to the default.
    fail_msg: (cause, ...text) => {
      if (!FAIL_TEXT[cause]) return console.warn(`script.txt: unknown cause [fail_msg:${cause}] (${Object.keys(FAIL_TEXT).join(', ')}), ignored`);
      if (text.join('')) challengeMsgs.fail[cause] = text.join(': '); else delete challengeMsgs.fail[cause];
    },
    // view_damping:<s>: (the scrolling world) smooth the world view's following of the ship over s game
    // seconds more (0 or no value: the default, close following); see moveCamera
    view_damping: (s) => {
      const v = s === undefined || s === 'off' ? 0 : parseFloat(s);
      if (!(v >= 0)) return console.warn(`script.txt: bad value [view_damping:${s}], ignored`);
      camDamp = v;
    },
    // put_hostile:[<n>:]<morphs>[:center|random|far][:<bullet>]: put n bots (default 1) into the
    // world, each a morph picked at random from <morphs> (names, ',' between) of the current
    // species (those it lacks are left out), placed as put_morph's (default random; center: the
    // first only, the rest random), never followed. Bots fire <bullet> (a bullet pattern, default
    // Scutium; one for all bots: the last [put_hostile]'s) now and then at the ship, in their own
    // bullets field: white, masking the G channel of whatever it hits, bots included (Hostile).
    // What they split into or merge with stays hostile (trackBots); [clear_world] removes them.
    put_hostile: (...args) => {
      const n = /^\d+$/.test(args[0] ?? '') ? parseInt(args.shift(), 10) : 1;
      let [ids = '', how = 'random', bullet = 'Scutium'] = args;
      if (!PUT_HOW.includes(how)) {
        if (args.length === 2 && Bullets.kindOf(how)) { bullet = how; how = 'random'; }   // [put_hostile:<morphs>:<bullet>]
        else { console.warn(`script.txt: unknown option [put_hostile:...:${how}], using random`); how = 'random'; }
      }
      if (!Hostile.set(bullet)) { console.warn(`script.txt: unknown bullet [put_hostile:...:${bullet}], using Scutium`); Hostile.set('Scutium'); }
      const pool = morphIds([ids]).filter(id => findMorph(id));
      if (!pool.length) return console.warn(`script.txt: [put_hostile:${ids}]: none of these morphs is in the species`);
      for (let k = 0; k < n; k++) putBot(pool[Math.floor(Math.random() * pool.length)], k > 0 && how === 'center' ? 'random' : creatureHow(how));
    },
    // min_hostile:<n>[:<morphs>[:center|random|far]][:<bullet>]: keep at least n bots on the board, as
    // [min_morph] keeps creatures: whenever there are fewer (bots dissolved, shot, or dropped), put in
    // one more, as [put_hostile] (a morph picked at random from <morphs>, names with ',' between,
    // default or all: the visible list; placed by how, default random; at most one every half game
    // second), firing <bullet> (default Scutium; one for all bots, as [put_hostile]'s). n = 0 or off:
    // stop. [clear_world] turns it off too.
    min_hostile: (n, ...args) => {
      const k = parseInt(n, 10);
      if (n === 'off' || k === 0) { minHostile = null; return; }
      if (!(k > 0)) return console.warn(`script.txt: bad value [min_hostile:${n}], ignored`);
      let [ids = '', how = 'random', bullet = 'Scutium'] = args;
      if (!PUT_HOW.includes(how)) {
        if (args.length === 2 && Bullets.kindOf(how)) { bullet = how; how = 'random'; }   // [min_hostile:<n>:<morphs>:<bullet>]
        else { console.warn(`script.txt: unknown option [min_hostile:...:${how}], using random`); how = 'random'; }
      }
      if (!Hostile.set(bullet)) { console.warn(`script.txt: unknown bullet [min_hostile:...:${bullet}], using Scutium`); Hostile.set('Scutium'); }
      const ids2 = ids === 'all' || !ids ? [] : morphIds([ids]);
      if (ids2.length && !ids2.some(id => findMorph(id))) console.warn(`script.txt: [min_hostile:${n}:${ids}]: none of these morphs is in the species`);
      minHostile = { n: k, ids: ids2, how, next: 0, puts: 0 };
      hostSpec = { ids: ids2, how, bullet };   // (for the Hostile slider, [settings])
    },
    // dot:<n>|off: the dots channel (dots.js): keep <n> dots on the board, which any creature
    // eats; the player's ship eating one scores 10 (the score shows at the top, from 0 when the
    // dots are turned on); off: no dots, no score.
    // dot:<pattern>[:<n>][:center|random|far]: put one <pattern> (Heptadentium) in the dots channel
    // instead of the round dots (the channel takes its rule: the round dots go), at the board's
    // centre, at a random clear spot (the default) or at the clear spot farthest from the ship.
    // Eaten by the ship, it gives <n> energy (default 10, as a dot); it is not replaced
    // ([need:clear_dot] waits for that).
    dot: (v, ...rest) => {
      if (v === 'off') return Dots.set('off');
      if (Dots.patternKey(v)) { console.warn(`script.txt: [dot:${v}] is now [beacon:${v}]`); return ACTIONS.beacon(v, ...rest); }
      const n = parseInt(v, 10);
      if (!(n >= 0)) return console.warn(`script.txt: bad value [dot:${v}], ignored`);
      let how = rest[0] ?? 'random';   // [dot:<n>:<how>]: where the eaten (or, scrolling, lost) dots are replaced
      if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [dot:${v}:${how}], using random`); how = 'random'; }
      Dots.set(n, how);
    },
    // beacon:<patterns>[:<n>][:center|random|far|extended[:<d>]]: put one beacon in the dots channel
    // (instead of the round dots: the channel takes its rule): a pattern picked at random from
    // <patterns> (names, ',' between; all: every one, dots.js PATTERNS), or the one the beacons on
    // the board already are, if listed (the channel runs one rule at a time: another would empty
    // it), placed by how (default random; extended: out of the view, in a scrolling world; <d>: d R
    // from the ship). Eaten by the ship, it gives n energy (default 10); not replaced
    // ([need:clear_dot] waits for it). In challenge mode beacons stay put (the round dots drift with
    // the star field).
    beacon: (list, ...rest) => putBeacon('beacon', list, rest, false),
    // keep_beacon:…: as beacon, but once the ship eats it another is put in at a new spot (by the same
    // how, the same energy), a pattern picked afresh from the list (dots.js put), again and again
    // ([need:clear_beacon] never completes)
    keep_beacon: (list, ...rest) => putBeacon('keep_beacon', list, rest, true),
    // gun[:<options>]: the player's gun (bullets.js): C or Fire shoots the pattern in use, which
    // masks (or erases) the creatures it hits while they wear it down. Options, any order:
    // off, erase|mask, rgb|r|g|b|color_now, wear_color|wear_all (defaults mask, color_now, wear_color).
    gun: (...options) => Bullets.setGun(options),
    // hostile_gun[:<options>]: the bots' bullets' own settings ([put_hostile], [min_hostile]): [gun]'s
    // options (mask | erase, rgb | r | g | b | color_now, wear_color | wear_all; off: the bots hold
    // fire, on: they fire again); follow, or none: as the player's [gun] (the default). Kept until changed.
    hostile_gun: (...options) => Hostile.setGun(options),
    // bullet_mutual:none|wear|masking: the player's and the bots' bullets together: apart (none, the
    // default), or each wearing the other down where they overlap (wear: shots meeting cancel out);
    // masking: not built yet (as none)
    // wear:<k>: the rate, a multiple of the bullets' own wear rate (0.03, ?bwear=; default 1)
    bullet_mutual: (v = 'none', k) => {
      if (v === 'masking') { console.warn('script.txt: [bullet_mutual:masking] is not built yet: none'); v = 'none'; }
      if (v !== 'none' && v !== 'wear') return console.warn(`script.txt: unknown option [bullet_mutual:${v}], ignored`);
      let rate = k === undefined ? 1 : parseFloat(k);
      if (!(rate >= 0)) { console.warn(`script.txt: bad rate [bullet_mutual:${v}:${k}], using 1`); rate = 1; }
      SimGL.setBulletMutual(v, rate);
    },
    // bullet:… (the old name): as add_bullet
    bullet: (...args) => { console.warn('script.txt: [bullet] is now [add_bullet]'); return ACTIONS.add_bullet(...args); },
    // add_bullet:<names>[:<name>]: the patterns the gun can fire: add names (separated by commas),
    // random (one more at random), all, or clear (none); nebula / antimatter put that world-draw
    // type in the book ([world_draw]); then equip <name> (added too if not listed). With two or more, ◀ ▶ (C / V) switch between them. The patterns it adds are shown:
    // the bullet book opens and reveals them (not when replayed: [jump_here], Reset, I).
    add_bullet: (names = '', equip) => {
      const before = Bullets.choices;   // (the world-draw types too, revealed as the patterns)
      Bullets.setList(names, equip);
      const added = Bullets.choices.filter(n => !before.includes(n));
      const fieldsAdded = added.filter(n => Bullets.fieldKinds.includes(n));   // (free, never saved)
      // all: for testing, kept out of the central inventory (patterns named, or random, are saved)
      if (names.split(',').some(t => t.trim().toLowerCase() === 'all')) {
        const named = names.split(',').map(t => t.trim()).concat(equip ? [equip.trim()] : []).map(Bullets.kindOf);
        for (const n of added) if (!named.includes(n)) unsavedBullets.add(n);
      }
      payBullets(added.filter(n => !unsavedBullets.has(n) && !fieldsAdded.includes(n)).length);   // each costs [new_bullet]'s price, as a created one ([bullet:all]'s extras free)
      // every pulse named in the list is revealed, even if the book had it already (created by the
      // player, or Antimatter / Nebula by [world_draw]); only the new ones are paid for. The one to
      // equip only if it is new (added): [add_bullet::<name>] just equipping one reveals nothing
      const named = names.split(',').map(t => t.trim().toLowerCase()).filter(Boolean);
      const shown = added.concat(Bullets.choices.filter(n => !added.includes(n) && named.some(t => t === n.toLowerCase() || Bullets.kindOf(t) === n)));
      if (!replaying && shown.length && Bullets.on && !book.open && !chapterShowing()) showNewBullets(shown);
    },
    // bullet_limit:<x>: the bullets channel is emptied once it holds more than x bullets' worth of
    // mass (bullets run into each other and grow); default 8
    bullet_limit: (x = '8') => {
      const v = parseFloat(x);
      if (!(v > 0)) { console.warn(`script.txt: bad value [bullet_limit:${x}], ignored`); return; }
      Bullets.setLimit(v);
    },
    // hostile_bullet_wait:<a>:<b>:<c>: bots fire every a to b game s (at random), a new bot's
    // first shot after c more; default 1.5:4:2 (also with no values)
    // bullet_wait:<a>:<b>:<c>: Autofire ([settings]) fires every a to b game s (at random), the first
    // shot after c; default 1.5:4:2 (also with no values)
    bullet_wait: (a, b, c) => {
      if (a === undefined) { fireWait = { ...FIRE_WAIT }; return; }
      const v = [a, b, c].map((x, i) => x === undefined ? Object.values(FIRE_WAIT)[i] : parseFloat(x));
      if (!v.every(x => x >= 0) || v[1] < v[0]) return console.warn(`script.txt: bad value [bullet_wait:${[a, b, c].filter(x => x !== undefined).join(':')}], ignored`);
      fireWait = { min: v[0], max: v[1], first: v[2] };
    },
    hostile_bullet_wait: (a, b, c) => {
      if (a === undefined) { botWait = { ...BOT_WAIT }; return; }
      const v = [a, b, c].map((x, i) => x === undefined ? Object.values(BOT_WAIT)[i] : parseFloat(x));
      if (!v.every(x => x >= 0) || v[1] < v[0]) return console.warn(`script.txt: bad value [hostile_bullet_wait:${[a, b, c].filter(x => x !== undefined).join(':')}], ignored`);
      botWait = { min: v[0], max: v[1], first: v[2] };
    },
    // hostile_bullet_limit:<x>: as bullet_limit, for the bots' bullets field ([put_hostile]); default 8
    hostile_bullet_limit: (x = '8') => {
      const v = parseFloat(x);
      if (!(v > 0)) { console.warn(`script.txt: bad value [hostile_bullet_limit:${x}], ignored`); return; }
      Hostile.setLimit(v);
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
    // stars:follow:<k>[:<bright>]: the star field (and the round dots) drift opposite the ship's
    // motion, as if it flew past them: the stars' middle layer and the dots at the ship's speed times k
    // <rate> (the fourth part, also after follow:<k>:<bright>): each star's brightness flickers on
    // its own, faster the higher (1: about one flicker every 2 to 7 s); 0 or left out: steady (draw.glsl)
    stars: (angle = '', speed = '', bright = '', rate = '') => {
      const num = (x, d) => x === '' ? d : parseFloat(x);
      const b = num(bright, 1), tw = num(rate, 0);
      if (!(tw >= 0)) return console.warn(`script.txt: bad twinkle rate [stars:…:${rate}], ignored`);
      if (angle === 'follow') {
        const k = num(speed, 1);
        if (!(k > 0) || !(b >= 0)) return console.warn(`script.txt: bad value [stars:follow:${speed}:${bright}], ignored`);
        return C.setStarsFollow(k, b, tw);
      }
      const [a, s] = [num(angle, 0), num(speed, 1)];
      if (![a, s, b].every(Number.isFinite) || s < 0 || b < 0)
        return console.warn(`script.txt: bad value [stars:${angle}:${speed}:${bright}], ignored`);
      C.setStars(a, s, b, tw);
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
    // (challenge mode: never on; the Morphs button opens the retry book there)
    // (challenge mode too, kept across challenges: e.g. set once in challenge.txt's header)
    clone: (v = 'on') => { if (v === 'off') { endClone(true); syncButtons(); } else { clone.on = true; syncButtons(); if (book.open) renderBook(); } },
    // erase:on|off: allow erase ("recycle", the selected R G B patch again) or not (back to masking)
    erase: (v = 'on') => { eraseOn = v !== 'off'; if (!eraseOn && C.action === 1) C.action = 0; },
    // draw:on|off|if:<species>: allow draw (add mass; the selected R G B patch again, before erase if
    // that is allowed: mask -> draw -> erase -> mask) or not (the default); if: only while the
    // species is <species> (its code, e.g. G6G6CR)
    draw: (v = 'on', sp) => {
      if (v === 'if') {
        if (!sp) return console.warn('script.txt: [draw:if] needs a species, e.g. [draw:if:G6G6CR]; ignored');
        drawRule = { on: true, species: sp };
      } else if (v === 'on' || v === 'off') drawRule = { on: v === 'on', species: null };
      else return console.warn(`script.txt: unknown option [draw:${v}], ignored`);
      if (!drawOn() && C.action === 2) C.action = 0;
    },
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
    // respawn: the ship has dissolved and the player has put it back, by Respawn (A) or Choose (a
    // pick in the morph book), since the line was reached
    respawn: () => { const n0 = respawns; return () => respawns > n0; },
    // draw:<x>: the player has painted at least x g of mass into the world's R G B fields (Nebula
    // or Antimatter, [world_draw]; all three together) since the line was reached: what the
    // paintbrush and soft brush add (cell values summed / R², as [need:mass]); the eraser, Clear
    // and Undo don't take it back. Progress: "draw: <mass>/<x>", one decimal.
    draw: (x) => {
      const v = parseFloat(x);
      if (!(v >= 0)) { console.warn(`script.txt: bad value [need:draw:${x}], ignored`); return () => true; }
      const d0 = Field.drawn, m = () => (Field.drawn - d0) / (C.R * C.R);
      const t = () => m() >= v;
      t.progress = () => `draw: ${Math.min(v, m()).toFixed(1)}/${v}`;
      return t;
    },
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
    // clear_hostile: no bots left ([put_hostile]; one missing from detection counts until dropped)
    // (counted as the red circles show: a twin or composite bot is one, whatever its parts)
    clear_hostile: () => { const t = () => bots.length === 0; t.progress = () => `hostile: ${hostileCount()}`; return t; },
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
    // survive[:<n>]: the ship dissolves (or is otherwise lost) after the line is reached, having
    // lasted at least n game seconds (default 0): a success, not a failure (in challenge mode,
    // [dissolve:retry] doesn't fail it meanwhile); lost sooner: in challenge mode it fails the
    // challenge (challengeFrame), elsewhere the line waits on. The time: the challenge's (from the
    // attempt's start, or the last respawn), elsewhere from when the line was reached. A challenge
    // that waits on it keeps its longest time as its best (challenge.longest). Progress: "time:
    // <t>/<n>" (no n: "time: <t>"), game seconds (paused time not counting).
    survive: (x) => {
      const v = x === undefined ? 0 : parseFloat(x);
      if (!(v >= 0)) { console.warn(`script.txt: bad value [need:survive:${x}], using 0`); }
      const need = v >= 0 ? v : 0;
      if (challenge.cur) challenge.longest = true;
      const n0 = C.simSteps;
      const lasted = () => (C.simSteps - (challenge.cur ? challenge.t0 : n0)) / C.STEPS_PER_S;
      const t = () => shipLost() && !emptySpace && lasted() >= need;
      t.survive = { need, lasted };
      t.progress = () => need > 0 ? `time: ${Math.floor(lasted())}/${need}` : `time: ${Math.floor(lasted())}`;   // (counting on past n)
      return t;
    },
    // clear_beacon: every beacon ([beacon]) eaten by the player's ship (one eaten by another creature,
    // or worn away, is put back nearby: dots.js respawnNear). Progress: "left: <n>".
    clear_beacon: () => {
      const t = () => Dots.beaconCount() === 0;
      return t;
    },
    // clear_dot: the dots channel holds nothing (every dot or pattern eaten; dots off counts too).
    // Progress: "left: <n>".
    clear_dot: () => {
      const t = () => !Dots.on || Dots.count === 0;
      t.progress = () => `left: ${Dots.on ? Dots.count : 0}`;
      return t;
    },
    // stay_target:<s>: the ship has stayed in the target area ([target]) for s game seconds in a
    // row (counted as [need:survive_time], paused time not counting); leaving it, or dissolving,
    // starts the count again from 0. Progress: "time: <t>/<s>".
    stay_target: (s) => {
      const v = parseFloat(s);
      if (!(v >= 0)) { console.warn(`script.txt: bad value [need:stay_target:${s}], ignored`); return () => true; }
      let n = 0, last = C.simSteps;
      const t = () => {
        const dn = C.simSteps - last;
        last = C.simSteps;
        if (!inTarget()) n = 0;
        else if (!book.open && !chapterShowing()) n += dn;
        return n >= v * C.STEPS_PER_S;
      };
      t.progress = () => `time: ${Math.min(v, Math.floor(n / C.STEPS_PER_S))}/${v}`;
      return t;
    },
    // stay_morph:<morph>:<s>: the ship has been detected as <morph> (e.g. circumgyrans, docked:
    // a gyrans circling a lithos, one object) for s game seconds in a row (counted as
    // [need:survive_time]); anything else, or dissolving, starts the count again from 0.
    // Progress: "time: <t>/<s>".
    stay_morph: (name, s) => {
      const v = parseFloat(s);
      if (!name || !(v >= 0)) { console.warn(`script.txt: bad [need:stay_morph:${name}:${s}], ignored`); return () => true; }
      let n = 0, last = C.simSteps;
      const t = () => {
        const dn = C.simSteps - last;
        last = C.simSteps;
        const d = C.focus && !shipLost() && C.shipDet();
        if (!d || d.name !== name) n = 0;
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
        const d = C.shipDet();
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
  // Then the [on_start[:<id>]] ... [end] blocks are taken out of the script's run (into
  // speciesBlocks, see chooseSpecies).
  function parseScript(txt) {
    // a tag left open at a line's end ([put_hostile:3:<newline>gyrans,...]) goes on in the next line
    const raw = [];
    for (const l of txt.split(/\r?\n/).map(l => l.trim())) {
      const open = raw.length && (raw[raw.length - 1].match(/\[/g) || []).length > (raw[raw.length - 1].match(/\]/g) || []).length;
      if (open && !l.startsWith('#')) raw[raw.length - 1] += l;
      else raw.push(l);
    }
    const lines = raw.filter(l => l && !l.startsWith('#')).map(l => {
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
      const open = L.actions.find(a => a.name === 'on_start' || a.name === 'on_choose_species' || a.name === 'on_challenge_default');
      if (open) {
        if (block) console.warn(`script.txt: [${open.name}] inside another block, ignored`);
        else if (open.name === 'on_challenge_default') { block = []; challengeDefaults = block; }
        else { block = []; speciesBlocks.set(open.args[0] || '', block); }
        continue;
      }
      if (L.actions.some(a => a.name === 'end')) {
        if (!block) console.warn('script.txt: [end] with no [on_start] or [on_challenge_default], ignored');
        block = null;
        continue;
      }
      if (!block) { out.push(L); continue; }
      if (L.text || L.conds.length) console.warn(`script.txt: a message or condition in a block ([on_start], [on_challenge_default]), ignored: ${L.text}`);
      block.push(L);
    }
    if (block) console.warn('script.txt: a block ([on_start], [on_challenge_default]) with no [end]: it runs to the end of the file');
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
  const CHECKPOINTS = ['chapter', 'new_scene', 'challenge'], REPLAY_SKIP = new Set([...CHECKPOINTS, 'jump_here', 'next_msg', 'music']);
  let resetting = false;
  // A retry (Retry in any mode; a challenge's retry) sets useMorph to the player's morph, the last
  // one put or chosen (the morph book's lighter one): the next [put_player] puts it in place of its
  // own (the initial morph), if it is a morph of the species, and clears it
  let useMorph = null;
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
      Dots.clear(); Bullets.clear(); Hostile.clear();   // the dots put in afresh, no bullets in flight
      replaying = true; useMorph = playerMorph;
      for (let i = cp; i < cur; i++)
        for (const a of SCRIPT[i].actions) if (!REPLAY_SKIP.has(a.name)) ACTIONS[a.name](...a.args);
      replaying = false;
      if (MODE === 'sandbox') { retake = true; cand = null; }   // a retried scene: the morph's picture may be retaken
      sceneHidden = false;
      run(cur);
      useMorph = null;
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
      const last = i + 1 >= scriptEnd(), tip = !tests.length && !last && tipPending;
      if (tip) tipPending = false;
      const fns = tests.map(t => t.progress).filter(Boolean);
      progress = fns.length ? { el: document.createElement('span'), fns } : null;
      showProgress();
      say(L.text, tests.length ? marks.act : last ? '' : tip ? `${marks.tip} ${marks.next}` : marks.next,
          wait, !tests.length && !last, progress?.el);   // cyan: the player has (or had) something to do
    }
    else if (!tests.length) next();   // tokens only, nothing to wait for
  }
  // (challenge mode: the challenge's end completes it)
  function next() {
    if (challenge.failed) return;   // (a failed challenge goes nowhere until a retry: Space / a click must not move it on, nor complete it)
    if (line + 1 < scriptEnd()) run(line + 1);
    else if (challenge.cur) completeChallenge();
  }
  // Space / a click on the message bar moves on: a message with no condition
  const byContinue = () => deferred < 0 && line >= 0 && SCRIPT[line].msg && !tests.length && !challenge.failed;
  // ?dev=1: O completes the waiting line's [need] conditions at once (a line waiting only on an
  // action, a chapter card, is left alone); P scores as if the ship ate a dot (+10 energy); I goes
  // to the next chapter; L records the frame times (devRecord)
  function devSkip() { if (DEV && deferred < 0 && line >= 0 && tests.length && SCRIPT[line].conds.length) { showSuccess(); next(); } }
  // I: on to the next [chapter] line. Every action on the lines in between runs first, in order
  // (as [jump_here] does), with their messages, conditions and chapter cards skipped.
  function devChapter() {
    if (!DEV || line < 0 || resetting) return;
    const cur = deferred >= 0 ? deferred : line + 1;   // a deferred line has not run its actions yet
    const to = SCRIPT.findIndex((L, i) => i > line && L.actions.some(a => a.name === 'chapter' || a.name === 'sub_chapter'));   // the next [chapter] or [sub_chapter]
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
  $('dev-rc').onclick = () => devRecord();
  // L (Rc): record every frame's times (core prof: the steps' GPU time, the analyses' readback and
  // CPU parts, ...) and what flies (the ship's pulses, the bots', both, none), up to REC_MAX_S; L again
  // stops and downloads a summary by those four cases (devRecSummary), to paste back
  const REC_MAX_S = 300;
  let rec = null;
  function devRecord() {
    if (!DEV) return;
    if (!rec) { rec = { t0: performance.now(), last: 0, rows: [] }; C.prof.on = true; $('dev-rc').classList.add('rec'); return; }
    const text = devRecSummary(rec);
    rec = null; C.prof.on = false; $('dev-rc').classList.remove('rec'); $('dev-rc').textContent = 'Rc';
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
    a.download = `cognoship-perf-${new Date().toISOString().slice(0, 19).replace(/[-:T]/g, '')}.txt`; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
  }
  // a frame's row (at the end of C.hooks.frame: hook = the game's frame work)
  function devRecFrame(hookT0) {
    if (!rec) return;
    const now = performance.now(), P = C.prof;
    rec.rows.push({ t: (now - rec.t0) / 1000, dt: rec.last ? hookT0 - rec.last : NaN, steps: P.steps,
      loop: P.loop, gpu: P.gpu, read: P.read, cpu: P.cpu, blocks: P.blocks, detect: P.detect, render: P.render, hook: now - hookT0,
      me: !SimGL.bulletsIdle, bot: !SimGL.hostileIdle, mask: SimGL.maskOn, vec: SimGL.vecOn, ...SimGL.regionStats,
      comps: C.found ? C.found.comps.length : 0, hostiles: hostileCount(), speed: C.SPEEDS[C.speedIdx] });
    rec.last = hookT0;
    const s = Math.floor((now - rec.t0) / 1000);
    $('dev-rc').textContent = `${s}s`;
    if (s >= REC_MAX_S) devRecord();
  }
  function devRecSummary(rec) {
    const rows = rec.rows.filter(r => !isNaN(r.dt));
    const f = (v, d = 1) => isFinite(v) ? v.toFixed(d) : '-';
    const mean = a => a.length ? a.reduce((x, y) => x + y, 0) / a.length : NaN;
    const pct = (a, p) => { if (!a.length) return NaN; const b = [...a].sort((x, y) => x - y); return b[Math.min(b.length - 1, Math.floor(p * b.length))]; };
    const sum = (a, k) => a.reduce((x, r) => x + r[k], 0);
    const caseOf = r => r.steps === 0 ? 'paused' : r.me && r.bot ? 'both' : r.me ? 'me' : r.bot ? 'bot' : 'none';
    const out = [];
    out.push('# Cognoship frame-time recording (dev L)');
    out.push(`url: ${location.search || '(none)'}  batch=${C.BATCH ? 1 : 0}  FFT=${FFTSim.on ? 1 : 0}  board ${SimGL.W}x${SimGL.H}  R=${C.R}`);
    out.push(`screen ${innerWidth}x${innerHeight} @${devicePixelRatio.toFixed(2)}  ${navigator.userAgent}`);
    const dbg = document.createElement('canvas').getContext('webgl2'), ext = dbg && dbg.getExtension('WEBGL_debug_renderer_info');
    if (ext) out.push(`gpu: ${dbg.getParameter(ext.UNMASKED_RENDERER_WEBGL)}`);
    out.push(`recorded ${f(rows.length ? rows[rows.length - 1].t : 0)} s, ${rows.length} frames, ${sum(rows, 'steps')} steps`);
    out.push('');
    out.push('Per case: frame = per frame, ms; step = per step, ms (gpu: the GPU work, read: the blocks readback, cpu: detection + tracking).');
    out.push('other = frame interval - (loop + blocks + detect + render + hook): the browser, compositing, the wait for vsync.');
    out.push('');
    const H = ['case', 'sec', 'frames', 'fps', 'steps/s', 'st/fr', 'n=1/2/3/4 %', 'int p50/p90', 'loop', 'gpu/st', 'read/st', 'cpu/st', 'blocks', 'detect', 'render', 'hook', 'other', 'mask%', 'vec%', 'stepped me/bot/dot', 'comps', 'hostile', 'speed'];
    const T = [H];
    for (const c of ['none', 'me', 'bot', 'both', 'paused']) {
      const a = rows.filter(r => caseOf(r) === c);
      if (!a.length) continue;
      const sec = sum(a, 'dt') / 1000, st = sum(a, 'steps'), n = a.length;
      const hist = [1, 2, 3, 4].map(k => Math.round(100 * a.filter(r => r.steps === k).length / n)).join('/');
      const pf = k => sum(a, k) / n, ps = k => st ? sum(a, k) / st : NaN;
      const other = mean(a.map(r => r.dt - r.loop - r.blocks - r.detect - r.render - r.hook));
      const ints = a.map(r => r.dt);
      T.push([c, f(sec), n, f(n / sec, 0), f(st / sec, 0), f(st / n, 2), hist, `${f(pct(ints, 0.5))}/${f(pct(ints, 0.9))}`,
        f(pf('loop')), f(ps('gpu'), 2), f(ps('read'), 2), f(ps('cpu'), 2), f(pf('blocks'), 2), f(pf('detect'), 2), f(pf('render'), 2), f(pf('hook'), 2), f(other),
        Math.round(100 * a.filter(r => r.mask).length / n), Math.round(100 * a.filter(r => r.vec).length / n),
        `${f(pf('bul'), 0)}/${f(pf('hos'), 0)}/${f(pf('dot'), 0)}%`, f(pf('comps'), 0), f(pf('hostiles'), 1), [...new Set(a.map(r => r.speed))].join(',')]);
    }
    const w = H.map((_, i) => Math.max(...T.map(r => String(r[i]).length)));
    for (const r of T) out.push(r.map((v, i) => String(v).padStart(w[i])).join('  '));
    out.push('');
    out.push('Per second: t, case (most frames), fps, steps/s, gpu/st, read/st, cpu/st, frame hook+detect+render ms, % of the board the pulses stepped me/bot');
    for (let s = 0; ; s++) {
      const a = rows.filter(r => Math.floor(r.t) === s);
      if (!a.length) { if (s > (rows.length ? rows[rows.length - 1].t : 0)) break; continue; }
      const cnt = {}; for (const r of a) cnt[caseOf(r)] = (cnt[caseOf(r)] || 0) + 1;
      const c = Object.entries(cnt).sort((x, y) => y[1] - x[1])[0][0];
      const sec = sum(a, 'dt') / 1000, st = sum(a, 'steps');
      out.push(`${String(s).padStart(4)}  ${c.padEnd(6)} ${f(a.length / sec, 0).padStart(3)} fps ${f(st / sec, 0).padStart(4)} st/s  gpu ${f(st ? sum(a, 'gpu') / st : NaN, 2)}  read ${f(st ? sum(a, 'read') / st : NaN, 2)}  cpu ${f(st ? sum(a, 'cpu') / st : NaN, 2)}  rest ${f((sum(a, 'hook') + sum(a, 'detect') + sum(a, 'render')) / a.length)}  stepped ${f(sum(a, 'bul') / a.length, 0)}/${f(sum(a, 'hos') / a.length, 0)}%`);
    }
    return out.join('\n') + '\n';
  }
  // the frame rate and the simulation steps per second, over the last half second
  $('dev-fps').hidden = !DEV;
  const fps = { t: 0, frames: 0, steps: 0 };
  // the tiles the bullets, the bots' bullets and the dots last stepped, of all (glsim.js regions)
  function regionText() { const r = SimGL.regionStats; return `stepped ${r.bul}/${r.hos}/${r.dot}%`; }
  function devFps() {
    if (!DEV) return;
    const now = performance.now();
    if (!fps.t) { Object.assign(fps, { t: now, frames: 0, steps: C.simSteps }); return; }
    fps.frames++;
    if (now - fps.t < 500) return;
    const s = (now - fps.t) / 1000;
    const parts = [`${Math.round(fps.frames / s)} fps`, `${Math.round((C.simSteps - fps.steps) / s)} steps/s`, `${C.stepMs.toFixed(1)} ms/step`, ...(C.DEBUG ? [`gpu ${C.msPerStep.toFixed(1)} · read ${C.msRead.toFixed(1)} · cpu ${C.msCpu.toFixed(1)} ms`] : []),
                   ...(FFTSim.on ? ['FFT'] : []), regionText(), `${innerWidth}×${innerHeight}`, `@${devicePixelRatio.toFixed(2)}`];   // (the CSS viewport and pixel ratio: for the phone test's sizes)
    $('dev-fps').textContent = parts.join(stacked ? '\n' : ' · ');   // stacked: one per row
    Object.assign(fps, { t: now, frames: 0, steps: C.simSteps });
  }
  // "SUCCESS" over the game window when a waiting line's [need] conditions are met: fades in, holds,
  // then slowly fades out (the .success animation, SUCCESS_MS); the game goes on meanwhile
  const SUCCESS_MS = 3000;
  let successTimer = null;
  // (challenge mode: only at a challenge's end, final, not for each line on the way)
  function showSuccess(final = false) {
    if (MODE === 'challenge' && !final) return;
    const el = $('success');
    clearTimeout(successTimer);
    el.hidden = false; el.classList.remove('show');
    void el.offsetWidth;   // restart the animation
    el.classList.add('show');
    successTimer = setTimeout(() => { el.hidden = true; el.classList.remove('show'); }, SUCCESS_MS);
  }
  function checkScript() {
    if (challenge.failed) return;   // (challenge mode: failed, waiting for a retry)
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
  // (a book open: a press on the message bar does nothing: the book is closed first, by Close, Esc
  // or a click beside it, and then the story goes on)
  $('msg').addEventListener('pointerdown', () => { if (!book.open) playerNext(); });

  // ---- background ([clear_world:<bg>]) --------------------------------------------------------
  // black (the default), or radiation: a dark orange behind the board, slowly darkening and
  // brightening again (a period of BG_PERIOD_S)
  // gray: a static dark gray
  // sky: a circular gradient, light cyan at the world view's bottom to dark cyan at its top (draw.glsl sky)
  const BACKGROUNDS = ['black', 'radiation', 'gray', 'sky'], BG_PERIOD_S = 8, GRAY = [0.13, 0.13, 0.14];
  let background = 'black';
  function updateBackground() {
    if (background === 'sky') return SimGL.setSky();
    if (background === 'gray') return SimGL.setBackground(GRAY);
    if (background !== 'radiation') return SimGL.setBackground([0, 0, 0]);
    const k = 0.45 + 0.15 * Math.sin(2 * Math.PI * performance.now() / 1000 / BG_PERIOD_S);   // 0.3 .. 0.6 of the full orange
    SimGL.setBackground([0.40 * k, 0.15 * k, 0.0]);
  }

  // ---- scene fades ([new_scene], [chapter]) -----------------------------------------------------
  const WORLD_FADE_MS = 500;
  // (the start, every mode: the views black (START_FADE) until the script's set-up has run, then faded
  // in as after a [new_scene]; unless a fade of the set-up's own has taken over: a chapter card, a species pick)
  const START_FADE = { from: 1, to: 1, t0: 0, ms: WORLD_FADE_MS };
  let worldFade = START_FADE;
  function worldShade() {
    const f = worldFade.ms > 0 ? Math.min(1, (performance.now() - worldFade.t0) / worldFade.ms) : 1;
    return worldFade.from + (worldFade.to - worldFade.from) * f;
  }
  const fadeDone = () => performance.now() - worldFade.t0 >= worldFade.ms;
  function fadeWorld(to, ms = WORLD_FADE_MS) { worldFade = { from: worldShade(), to, t0: performance.now(), ms }; }
  let sceneHidden = false;   // faded out by [new_scene]: the next message fades the views back in
  C.opt.worldShade = worldShade;

  // [beacon] / [keep_beacon] (tag; keep: put in afresh once the ship eats it): see beacon
  function putBeacon(tag, list, rest, keep) {
    const names = list.split(',').map(t => t.trim()).filter(Boolean);
    const all = names.some(t => t.toLowerCase() === 'all'), named = names.filter(t => t.toLowerCase() !== 'all');
    const bad = named.filter(t => !Dots.patternKey(t));
    if (bad.length) console.warn(`script.txt: unknown beacon pattern [${tag}:${bad.join(',')}], ignored`);
    const ok = all ? Dots.PATTERN_KEYS : named.map(t => Dots.patternKey(t)).filter(Boolean);
    if (!ok.length) return;
    const v = Dots.beaconCount() > 0 && ok.includes(Dots.kind) ? Dots.kind : ok[Math.floor(Math.random() * ok.length)];
    const points = rest.length && /^\d+(\.\d+)?$/.test(rest[0]) ? parseFloat(rest.shift()) : undefined;
    let how = rest[0] ?? 'random';
    if (!PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [${tag}:${v}:${how}], using random`); how = 'random'; }
    const dist = rest[1] !== undefined ? parseFloat(rest[1]) : undefined;
    if (rest[1] !== undefined && !(dist > 0)) console.warn(`script.txt: bad distance [${tag}:${v}:…:${rest[1]}], ignored`);
    if (!Dots.put(v, how, points, dist > 0 ? dist : undefined, keep ? ok : null)) console.warn(`script.txt: [${tag}:${v}:${how}]: no clear spot found`);
  }
  // A bot ([put_hostile], [min_hostile]): morph id put in by how (creatureHow's), never followed,
  // firing after the first wait
  function putBot(id, how) {
    C.selectMorph(findMorph(id).id);
    const at = C.putPattern(how, 'never');
    bots.push({ col: at[0], row: at[1], radius: C.R, next: C.simSteps + Math.round(botWait.first * C.STEPS_PER_S) + botDelay(), seen: C.simSteps, live: true });
    botsFound = null;
  }
  // ---- [min_hostile] ------------------------------------------------------------------------
  // Each frame: fewer bots than minHostile.n -> put one more in (as checkMinMorph, at most one every
  // MIN_MORPH_STEPS; its how places the first n, a scrolling world's later ones come in at its edges)
  let minHostile = null;   // { n, ids (none: the visible list), how, next, puts }
  function checkMinHostile() {
    if (!minHostile || !C.found || C.simSteps < minHostile.next || bots.length >= minHostile.n) return;
    const pool = (minHostile.ids.length ? minHostile.ids : visibleList()).filter(id => findMorph(id));
    if (!pool.length) return;
    const first = minHostile.puts < minHostile.n;
    const how = scroll && !first ? scrollRefill(C.R) || minHostile.how : minHostile.how;
    putBot(pool[Math.floor(Math.random() * pool.length)], creatureHow(how));
    minHostile.puts++;
    minHostile.next = C.simSteps + MIN_MORPH_STEPS;
  }

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
    const pool = (minMorph.ids.length ? minMorph.ids : visibleList().map(id => ({ id, n: 1 }))).filter(e => findMorph(e.id));
    if (!pool.length) return;
    const e = pool[Math.floor(Math.random() * pool.length)];
    // its how places the first n; in a scrolling world every later one is a replacement: of one taken
    // away at an edge on the opposite edge, else on the edge the ship heads to (scrollRefill); by
    // count of puts, not of creatures (while the ship moves, objects may leave at the edges as fast
    // as they come, and the count never reach n)
    const first = (minMorph.puts || 0) < minMorph.n;
    const how = scroll && !first ? scrollRefill(e.n > 1 ? 2 * C.R : C.R) || minMorph.how : minMorph.how;
    if (e.n > 1) putGroup(e.id, e.n, how, 'min_morph');
    else putMorph('min_morph', e.id, how, 'never');   // never the new focus
    minMorph.puts = (minMorph.puts || 0) + 1;
    minMorph.next = C.simSteps + MIN_MORPH_STEPS;
  }

  // ---- the scrolling world ([clear_world:scroll]) -------------------------------------------
  // The board stays a torus, a window on an unbounded world that moves with the ship: the world
  // view (the whole board, as ever) follows the ship, keeping it at its middle (core
  // opt.viewShift; smoothly: moveCamera), so the board's wrap, the seam, runs along the view's edges, opposite the ship;
  // a band there is cleared each frame (glsim clearSeam: creatures, dots, bullets; the bullets also
  // every bullet step, wider: glsim setBulletSeam; the world-draw fields: field.js), so nothing
  // wraps round to the ship's front. The view has no outside (the board is all there is), so an
  // object (a part, not the ship's) is taken away (a disc round it cleared) as soon as its body
  // (its circle, twice its radius) touches an edge's band, and [min_morph] puts its replacement on
  // the opposite edge, its centre EDGE_IN inside: the band, a body (R) and a margin of R / 10
  // (edgeQueue); [min_morph]'s other replacements (an object eaten, shot or dissolved) come in on
  // the edge the ship heads to (standing still: any), as far inside (edgeSpot). Its how places
  // only its first n.
  // Round dots the same (dots.js). World points: the camera's (cx, cy) moves as the ship does on
  // the board (bx, by); a board cell is at cx + its offset from bx (round the torus). What lies
  // beyond the board (extended: the beacons, the target) is kept by its world point and put on the
  // board as the ship nears it; a yellow arrow at the view's edge points to each one off the board.
  const SEAM_BAND = 6, EDGE_IN = SEAM_BAND + 1.1 * 18;   // 25.8 cells: creatures put elsewhere stay this far inside the edges
  const VX = 384 / 2, VY = 256 / 2;   // the view's half-size (the board's), cells
  const POINTER_COLOR = '#fff200';   // bright yellow (core draws it big and blinking: blink)
  let scroll = null;   // { cx, cy, bx, by, ux, uy (the heading, or 0, 0) } (cells), or null: the ordinary torus
  let edgeQueue = [];   // board spots for [min_morph]'s next replacements (opposite an object taken away at an edge)
  const toWorld = (col, row) => [scroll.cx + C.wrapD(col - scroll.bx, C.BW), scroll.cy + C.wrapD(row - scroll.by, C.BH)];
  // a world point's board cell, or null if farther than half the board less margin from the ship (on either axis)
  function toBoard(wx, wy, margin = 0) {
    const dx = wx - scroll.cx, dy = wy - scroll.cy;
    if (Math.abs(dx) > C.BW / 2 - margin || Math.abs(dy) > C.BH / 2 - margin) return null;
    return [(((scroll.bx + dx) % C.BW) + C.BW) % C.BW, (((scroll.by + dy) % C.BH) + C.BH) % C.BH];
  }
  // whether a world point (offset dx, dy from the ship) is in the view, give or take margin
  const inView = (dx, dy, margin = 0) => Math.abs(dx) <= VX + margin && Math.abs(dy) <= VY + margin;
  const off = (col, row) => [C.wrapD(col - scroll.bx, C.BW), C.wrapD(row - scroll.by, C.BH)];
  const toCell = (dx, dy) => [(((scroll.bx + dx) % C.BW) + C.BW) % C.BW, (((scroll.by + dy) % C.BH) + C.BH) % C.BH];
  // A board spot where a thing of that reach (cells) may go: anywhere on the torus; in a scrolling
  // world, at least EDGE_IN (or half its reach) inside the view's edges
  // world, at least EDGE_IN (or half its reach) inside the view's edges. Either way, away from the
  // beacons and the target (keepClear).
  const spotOk = (col, row, reach) => {
    if (!keepClear(col, row, reach)) return false;
    if (!scroll) return true;
    const [dx, dy] = off(col, row), m = Math.max(EDGE_IN, reach / 2);
    return Math.abs(dx) < VX - m && Math.abs(dy) < VY - m;
  };
  // A creature of that reach (cells) at (col, row) would be far from every beacon (it would eat it)
  // and from the target area (it would sit in it): KEEP_AWAY between them. In a scrolling world by
  // world points, so the beacons and the target off the board (parked) are kept clear too.
  const KEEP_AWAY = 5 * 18;   // cells (5 R)
  function keepClear(col, row, reach) {
    const gap = reach + KEEP_AWAY;
    if (scroll) {
      const [wx, wy] = toWorld(col, row);
      for (const [bx, by] of Dots.beacons()) if (Math.hypot(bx - wx, by - wy) < Dots.beaconRadius() + gap) return false;
      if (target && target.wx !== undefined) {
        const h = target.size / 2 + gap;
        if (Math.abs(target.wx - wx) < h && Math.abs(target.wy - wy) < h) return false;
      }
      return true;
    }
    for (const [c, r, rad] of Dots.beaconCells()) if (Math.hypot(C.wrapD(c - col, C.BW), C.wrapD(r - row, C.BH)) < rad + gap) return false;
    if (target) {
      const h = target.size / 2 + gap;
      if (Math.abs(C.wrapD(target.col - col, C.BW)) < h && Math.abs(C.wrapD(target.row - row, C.BH)) < h) return false;
    }
    return true;
  }
  C.opt.spotOk = spotOk;
  C.opt.center = () => scroll ? [scroll.bx, scroll.by] : null;
  // center: the board's centre; in a scrolling world the view's middle (the ship's place)
  const centerPt = () => scroll ? [scroll.bx, scroll.by] : [C.BW / 2, C.BH / 2];
  // Edge placements: a placer { place(w, h) } puts a thing with a w x h bounding box (the rotated
  // form's, known only when it is put: core putPattern, putGroup, dots.js) against the inner side of
  // an edge's band, EDGE_GAP from it: sx, sy = +1 / -1 for the right / bottom, left / top edge (0:
  // not that axis), ax, ay = the other coordinate (offset from the ship), kept inside the band.
  // place returns the board cell for its centre, and marks it fresh (FRESH_S: not taken away again
  // meanwhile, as its detection circle may touch the band).
  const EDGE_GAP = 0.1 * 18, FRESH_S = 1.5;
  let fresh = [];   // [col, row, step] of things just put at an edge
  function placer(sx, sy, ax, ay) {
    const lx = VX - SEAM_BAND - EDGE_GAP, ly = VY - SEAM_BAND - EDGE_GAP;
    return { sx, sy, ax, ay, place(w, h) {
      const dx = sx ? sx * (lx - w / 2) : Math.max(-(lx - w / 2), Math.min(lx - w / 2, ax));
      const dy = sy ? sy * (ly - h / 2) : Math.max(-(ly - h / 2), Math.min(ly - h / 2, ay));
      const c = toCell(dx, dy);
      fresh.push([c[0], c[1], C.simSteps]);
      return c;
    } };
  }
  // where a placer would put a thing of half-size body, without marking it (to check it is clear)
  const placeAt = (p, body) => { const c = p.place(2 * body, 2 * body); fresh.pop(); return c; };
  // An object at (col, row) whose body (half-size `body`, cells) touches an edge's band: a placer on
  // the opposite edge (the other coordinate kept) for its replacement; null if it touches none
  function edgeOut(col, row, body) {
    if (!scroll) return null;
    const [dx, dy] = off(col, row);
    const outX = Math.abs(dx) + body >= VX - SEAM_BAND, outY = Math.abs(dy) + body >= VY - SEAM_BAND;
    if (!outX && !outY) return null;
    return placer(outX ? -Math.sign(dx) : 0, outY ? -Math.sign(dy) : 0, dx, dy);
  }
  // A placer on the edge the ship heads to (standing still: any), give or take 40°, clear of
  // creatures (body: the thing's half-size); null if not a scrolling world or none was found
  function edgeSpot(body = C.R) {
    if (!scroll) return null;
    const moving = scroll.ux || scroll.uy, a0 = moving ? Math.atan2(scroll.uy, scroll.ux) : Math.random() * 2 * Math.PI;
    const comps = C.found ? C.found.comps : [], lx = VX - SEAM_BAND, ly = VY - SEAM_BAND;
    for (let tries = 0; tries < 40; tries++) {
      const a = a0 + (2 * Math.random() - 1) * 0.7, ux = Math.cos(a), uy = Math.sin(a);
      const tx = Math.abs(ux) > 1e-6 ? lx / Math.abs(ux) : Infinity, ty = Math.abs(uy) > 1e-6 ? ly / Math.abs(uy) : Infinity;
      const p = tx < ty ? placer(Math.sign(ux), 0, 0, tx * uy) : placer(0, Math.sign(uy), ty * ux, 0);   // the edge the ray meets
      const [col, row] = placeAt(p, body);
      if (clearAt(col, row, body, comps)) return p;
    }
    return null;
  }
  const clearAt = (col, row, body, comps = C.found ? C.found.comps : []) => keepClear(col, row, body) &&
    !comps.some(c => Math.hypot(C.wrapD(c.col - col, C.BW), C.wrapD(c.row - row, C.BH)) < (c.radius || C.R) + body);
  // A queued placer (edgeQueue), moved along its edge (up to 40 cells either way) until clear of creatures
  function queuedSpot(p, body = C.R) {
    for (let k = 0; k < 9; k++) {
      const s = k === 0 ? 0 : (k % 2 ? 1 : -1) * 10 * Math.ceil(k / 2);
      const q = placer(p.sx, p.sy, p.ax + (p.sx ? 0 : s), p.ay + (p.sx ? s : 0));
      const [col, row] = placeAt(q, body);
      if (clearAt(col, row, body)) return q;
    }
    return p;
  }
  // [min_morph]'s placer for a creature in a scrolling world: the next queued one (opposite an
  // object taken away at an edge), else (edge: its first n are in) on the edge the ship heads to; null: as its how
  function scrollRefill(body, edge = true) {
    if (!scroll) return null;
    if (edgeQueue.length) return queuedSpot(edgeQueue.shift(), body);
    return edge ? edgeSpot(body) : null;
  }
  // A world point out of the world view: at a random angle, 1.25 to 2 times the distance from the
  // ship to the view's edge that way, plus reach; or (dist, in R) exactly dist R from the ship, at a
  // random angle (a run's distances, comparable between runs)
  function extendedPoint(reach, dist) {
    const a = Math.random() * 2 * Math.PI, ux = Math.cos(a), uy = Math.sin(a);
    if (dist > 0) return [scroll.cx + dist * C.R * ux, scroll.cy + dist * C.R * uy];
    const t = Math.min(Math.abs(ux) > 1e-6 ? VX / Math.abs(ux) : Infinity, Math.abs(uy) > 1e-6 ? VY / Math.abs(uy) : Infinity);
    const d = t * (1.25 + 0.75 * Math.random()) + reach;
    return [scroll.cx + d * ux, scroll.cy + d * uy];
  }
  function startScroll(on) {
    scroll = on ? { cx: 0, cy: 0, bx: C.BW / 2, by: C.BH / 2, ux: 0, uy: 0 } : null;
    edgeQueue = []; fresh = [];
    if (!on) { SimGL.setBulletSeam(null); Field.clearSeam(null); }   // (the ordinary torus: no seam)
    Dots.setWorld(on ? { toBoard, toWorld, extended: extendedPoint, ok: spotOk, center: () => [scroll.bx, scroll.by], edgeOut, edgeSpot, queuedSpot } : null);
  }
  // Each frame: the camera follows the ship; the objects at an edge are taken away (their
  // replacements queued while [min_morph] runs); the seam is cleared; the target is put on the
  // board (or off it)
  // The camera (the view's centre, bx, by) follows the ship smoothly, not snapped to it each frame:
  // the ship's detected centre wobbles (a twin's parts not yet named as one, a part's detection
  // jumping), and the whole moves the view. Each frame it moves with the ship's smoothed velocity
  // (core kin.vel, as the arrow: no lag behind a moving ship), then closes CAM_TAU_STEPS of the gap
  // (game steps: as fast at any speed; also CAM_TAU_S real seconds, so it settles while paused); a
  // gap over CAM_SNAP R (a respawn elsewhere, a new scene) is jumped at once.
  // [view_damping:<s>] (camDamp, game seconds; 0, the default: as above): the velocity the view moves
  // with is itself smoothed over s (camV, an exponential average of the ship's), and the gap closes s
  // more slowly, so a ship that zigzags or wobbles fast (HAESRE zigzag) no longer shakes the world: the
  // view follows its mean course, the ship swinging about the middle.
  const CAM_TAU_STEPS = 20, CAM_TAU_S = 1, CAM_SNAP = 3;
  let camT = 0, camSteps = 0, camDamp = 0, camV = null;
  function moveCamera(cam, sh) {   // cam: scroll, or roll (no world point: cx, cy)
    const now = performance.now(), dtS = camT ? Math.min(0.25, (now - camT) / 1000) : 0, dn = C.simSteps - camSteps;
    camT = now; camSteps = C.simSteps;
    let gx = C.wrapD(sh.col - cam.bx, C.BW), gy = C.wrapD(sh.row - cam.by, C.BH);
    let mx, my;
    const v0 = sh.kin ? sh.kin.vel : [0, 0];
    if (Math.hypot(gx, gy) > CAM_SNAP * C.R || !cam.cam) { mx = gx; my = gy; camV = v0.slice(); }   // jump (and the first frame)
    else {
      const ds = camDamp * C.STEPS_PER_S;   // (in steps)
      if (!camV || !ds) camV = v0.slice();
      else { const b = 1 - Math.exp(-dn / ds); camV[0] += b * (v0[0] - camV[0]); camV[1] += b * (v0[1] - camV[1]); }
      const v = camV, px = v[0] * C.dt * dn, py = v[1] * C.dt * dn;   // the ship's own motion since the last frame (damped: its mean)
      gx -= px; gy -= py;
      const a = 1 - Math.exp(-(dn / (CAM_TAU_STEPS + ds) + dtS / (CAM_TAU_S + camDamp)));
      mx = px + a * gx; my = py + a * gy;
    }
    cam.cam = true;
    if (cam.cx !== undefined) { cam.cx += mx; cam.cy += my; }
    cam.bx = (((cam.bx + mx) % C.BW) + C.BW) % C.BW; cam.by = (((cam.by + my) % C.BH) + C.BH) % C.BH;
  }
  // World: Centered ([settings]): the ordinary torus, the world view rolled round it (core
  // opt.viewShift) so the ship stays in its middle; the camera as the scrolling world's (moveCamera),
  // with nothing else of it (no seam, no stripe, no arrows)
  let roll = null;   // { bx, by (the view's centre, cells), cam }, or null: the view unrolled
  function setRoll(on) {
    if (!on) roll = null;
    else if (!roll) { roll = { bx: C.BW / 2, by: C.BH / 2, cam: false }; camT = 0; }
  }
  function centerFrame() {
    if (!roll) return;
    const sh = C.ship();
    if (sh && !shipLost()) moveCamera(roll, sh); else camT = 0;
  }
  function scrollFrame() {
    if (!scroll) return;
    const sh = C.ship();
    if (sh && !shipLost()) {
      moveCamera(scroll, sh);
      const k = sh.kin, moving = k && k.speed >= C.SPEED_MIN;
      scroll.ux = moving ? Math.cos(k.velAng) : 0; scroll.uy = moving ? Math.sin(k.velAng) : 0;
    } else { scroll.ux = scroll.uy = 0; camT = 0; }
    const discs = [];
    fresh = fresh.filter(f => C.simSteps - f[2] < FRESH_S * C.STEPS_PER_S);
    if (C.found) {
      const ship = C.focusIdx >= 0 ? (C.objects().find(o => o.parts.includes(C.focusIdx))?.parts ?? [C.focusIdx]) : [];
      C.found.comps.forEach((c, i) => {
        if (ship.includes(i) || discs.length >= 8) return;
        if (fresh.some(([fc, fr]) => Math.hypot(C.wrapD(c.col - fc, C.BW), C.wrapD(c.row - fr, C.BH)) < C.R)) return;   // just put at an edge
        const body = 2 * (c.radius || C.R);   // (its circle, as the lab draws it)
        const p = edgeOut(c.col, c.row, body);
        if (!p) return;
        discs.push([c.col, c.row, body + 4]);
        if (minMorph) edgeQueue.push(p);   // (its replacement: the opposite edge)
      });
    }
    SimGL.clearSeam(scroll.bx + C.BW / 2, scroll.by + C.BH / 2, SEAM_BAND, discs);
    SimGL.setBulletSeam(scroll.bx + C.BW / 2, scroll.by + C.BH / 2);   // the bullets: cleared there every bullet step too, wider (bulletsim.glsl)
    Field.clearSeam(scroll.bx + C.BW / 2, scroll.by + C.BH / 2, SEAM_BAND);   // the world-draw fields (field.js)
    if (discs.length) C.analyze(0);   // (gone now: not counted, nor taken again)
    if (target && target.wx !== undefined) {
      const b = toBoard(target.wx, target.wy, target.size / 2 + SEAM_BAND);
      target.off = !b;
      if (b) [target.col, target.row] = b;
    }
  }
  C.opt.viewShift = () => scroll ? [scroll.bx - C.BW / 2, scroll.by - C.BH / 2] : roll ? [roll.bx - C.BW / 2, roll.by - C.BH / 2] : null;
  // a black stripe 2 R wide inside the world view's edges (the seam and the objects being taken away
  // out of sight); only drawn: nothing else reads it
  const STRIPE = () => 2 * C.R;
  C.opt.viewInset = () => scroll ? STRIPE() : 0;
  // the arrows: to the beacons and the target out of the view
  C.opt.pointers = () => {
    if (!scroll || challenge.failed) return [];   // (a failed attempt: none)
    const pts = Dots.beacons();
    if (target && target.wx !== undefined) pts.push([target.wx, target.wy]);
    return pts.map(([wx, wy]) => ({ dx: wx - scroll.cx, dy: wy - scroll.cy, color: POINTER_COLOR, blink: true })).filter(p => !inView(p.dx, p.dy, -STRIPE()));   // (out of the view's visible part: behind its black stripe)
  };

  // ---- bots ([put_hostile]) ---------------------------------------------------------------
  // Detection names objects afresh each frame, with no identities, so bots are tracked by place:
  // at each new detection, every object besides the ship within a bot's reach (the larger of the
  // two radii, plus TRACK_GAP) of a bot's last place is a bot. So a bot stays hostile whatever it
  // turns into (another morph, ???), both halves of a bot that splits are bots, and an object a bot
  // merges with is one; a bot the ship hijacks (now the followed ship) is the player's. A bot
  // missing from a detection (detection flickers in a collision) is kept at its last place for
  // BOT_GRACE seconds, not firing, then dropped.
  // game s between a bot's shots (min to max, at random), and before a new bot's first (first more): [hostile_bullet_wait]
  const BOT_WAIT = { min: 1.5, max: 4, first: 2 };
  let botWait = { ...BOT_WAIT };
  const BOT_SPREAD = 0.35;                     // aim: towards the ship, give or take this (radians)
  const BOT_GRACE = 1;                         // game s
  const BOT_COLOR = '#ff4d4d';
  let bots = [];   // { col, row, radius (cells), next: the step of its next shot, seen: the step last detected, live: in the last detection, k: its index in that detection's comps }
  let botsFound = null;   // the detection last tracked (C.found)
  const botDelay = () => Math.round((botWait.min + Math.random() * (botWait.max - botWait.min)) * C.STEPS_PER_S);
  function trackBots() {
    if (!bots.length || !C.found || C.found === botsFound) return;
    botsFound = C.found;
    const now = C.simSteps, gap = C.R / 2, used = new Set(), out = [];
    C.found.comps.forEach((c, i) => {
      if (i === C.focusIdx) return;
      const r = c.radius || C.R;
      let best = null, bd = Infinity;
      for (const b of bots) {
        const d = Math.hypot(C.wrapD(c.col - b.col, C.BW), C.wrapD(c.row - b.row, C.BH));
        if (d <= Math.max(b.radius, r) + gap && d < bd) { best = b; bd = d; }
      }
      if (!best) return;
      out.push({ col: c.col, row: c.row, radius: r, next: used.has(best) ? now + botDelay() : best.next, seen: now, live: true, k: i });   // a second half: its own timing
      used.add(best);
    });
    for (const b of bots) if (!used.has(b) && now - b.seen < BOT_GRACE * C.STEPS_PER_S) { b.live = false; out.push(b); }
    bots = out;
  }
  // The hostile objects as shown (core objects(): a twin or trio named as one is one object): those
  // with a bot among their parts, not the ship's; each with its bots ({ ...object, bots })
  function hostileObjects() {
    if (!bots.length || !C.found) return [];
    if (botsFound !== C.found) trackBots();   // (the same detection as the bots' indices)
    const byK = new Map();
    for (const b of bots) if (b.live) byK.set(b.k, b);
    const out = [];
    for (const o of C.objects()) {
      if (o.parts.includes(C.focusIdx)) continue;
      const bs = o.parts.map(k => byK.get(k)).filter(Boolean);
      if (bs.length) out.push({ ...o, bots: bs });
    }
    return out;
  }
  // The bots as shown, for [need:clear_hostile]'s count: the red circles (botMarks), one per hostile
  // object (a twin or composite as one). Bots missing from detection a moment (BOT_GRACE) have no
  // circle and aren't counted, but the count shows at least 1 until every bot is gone (the test).
  const hostileCount = () => bots.length ? Math.max(1, C.found ? hostileObjects().length : 0) : 0;
  // Each frame: every hostile object whose time has come (its bots' earliest) fires at the ship
  // (none: anywhere), from the object as a whole (its centre and reach: a twin fires as one); not
  // under a card, in a species fade, or after a failure
  function botsFire() {
    if (!bots.length || !Hostile.on || chapterShowing() || speciesFading || challenge.failed) return;
    const now = C.simSteps, ship = C.ship();
    for (const o of hostileObjects()) {
      if (now < Math.min(...o.bots.map(b => b.next))) continue;
      const next = now + botDelay();
      for (const b of o.bots) b.next = next;
      const ang = ship ? Math.atan2(C.wrapD(ship.row - o.row, C.BH), C.wrapD(ship.col - o.col, C.BW)) + (2 * Math.random() - 1) * BOT_SPREAD
                       : Math.random() * 2 * Math.PI;
      Hostile.fire(o.col, o.row, o.reach, ang);
    }
  }
  // a red circle round each hostile object, its morph's name above it (??? if unknown)
  // (its name in white; an unknown morph: none)
  const botMarks = () => hostileObjects().map(o => ({ col: o.col, row: o.row, r: o.radius + 2, color: BOT_COLOR, label: o.label || null, labelColor: '#fff' }));
  C.opt.detectNear = () => bots.map(b => [b.col, b.row]);   // the bots' morphs are named (Detect: Focused names the ship's only)

  // ---- groups ([min_morph:<n>:<morph>=<k>], [put_morph:<morph>=<k>], [put_player:<morph>=<k>]) --
  // k copies of a morph packed together: each turned at random, then placed in turn at the free
  // spot nearest the group's centre (rings outward, each ring from a random angle), its body
  // (cells above GROUP_BODY in any channel) at least GROUP_GAP cells from the bodies already
  // placed. The group's centre: a spot clear of creatures, chosen by how (center, random, far: of
  // 30 such, the one farthest from the ship). focusNew: the first placed (nearest the centre)
  // becomes the followed ship ([put_player]); otherwise none is followed. False if the morph is unknown.
  const GROUP_MAX = 12, GROUP_GAP = 3, GROUP_BODY = 0.01;
  // <morph>[=<k>] -> { id, k } (k: 1 to GROUP_MAX)
  function groupArg(token, x) {
    const [id, m] = String(x).split('='), k = m === undefined ? 1 : parseInt(m, 10);
    if (!(k >= 1)) console.warn(`script.txt: bad group size [${token}:${x}], using 1`);
    return { id, k: k >= 1 ? Math.min(k, GROUP_MAX) : 1 };
  }
  function putGroup(id, k, how, token, focusNew = false) {
    const p = findMorph(id);
    if (!p) { console.warn(`script.txt: unknown morph [${token}:${id}=${k}]`); return false; }
    if (typeof how === 'string' && !PUT_HOW.includes(how)) { console.warn(`script.txt: unknown option [${token}:${id}=${k}:${how}], using center`); how = 'center'; }
    if (how === 'extended') how = 'far';
    const f0 = Rule.formCells(p.cells, p.R, C.R);
    const forms = Array.from({ length: k }, () => Rule.rotate(f0, Math.random() * 2 * Math.PI));
    const D = Math.max(...forms.map(f => Math.max(f.w, f.h)));
    const G = Math.ceil(D * (k + 1) + 2 * GROUP_GAP * (k + 1)), c0 = G >> 1;
    const blocked = new Uint8Array(G * G);
    // a form's body: the cells (x, y) above GROUP_BODY in any channel
    const body = f => {
      const out = [];
      for (let y = 0; y < f.h; y++) for (let x = 0; x < f.w; x++)
        if (f.cells[0][y][x] > GROUP_BODY || f.cells[1][y][x] > GROUP_BODY || f.cells[2][y][x] > GROUP_BODY) out.push(x, y);
      return out;
    };
    const fits = (b, ox, oy) => {
      for (let i = 0; i < b.length; i += 2) {
        const x = ox + b[i], y = oy + b[i + 1];
        if (x < 0 || y < 0 || x >= G || y >= G || blocked[y * G + x]) return false;
      }
      return true;
    };
    const disc = [];   // the gap's disc of offsets
    for (let dy = -GROUP_GAP; dy <= GROUP_GAP; dy++) for (let dx = -GROUP_GAP; dx <= GROUP_GAP; dx++) if (dx * dx + dy * dy <= GROUP_GAP * GROUP_GAP) disc.push(dx, dy);
    const placed = [];   // { f, ox, oy }: top-left in the group's grid
    for (const f of forms) {
      const b = body(f);
      let at = null;
      for (let r = 0; r < c0 && !at; r++) {
        const steps = Math.max(1, Math.ceil(2 * Math.PI * r)), a0 = Math.random() * 2 * Math.PI;
        for (let j = 0; j < steps && !at; j++) {
          const a = a0 + 2 * Math.PI * j / steps;
          const ox = Math.round(c0 + r * Math.cos(a) - f.w / 2), oy = Math.round(c0 + r * Math.sin(a) - f.h / 2);
          if (fits(b, ox, oy)) at = [ox, oy];
        }
      }
      if (!at) { console.warn(`[${token}]: no room for ${id} #${placed.length + 1} in its group`); break; }
      for (let i = 0; i < b.length; i += 2)
        for (let d = 0; d < disc.length; d += 2) {
          const x = at[0] + b[i] + disc[d], y = at[1] + b[i + 1] + disc[d + 1];
          if (x >= 0 && y >= 0 && x < G && y < G) blocked[y * G + x] = 1;
        }
      placed.push({ f, ox: at[0], oy: at[1] });
    }
    if (!placed.length) return true;
    // the group's extent from its centre, then a spot for it
    let ext = 0;
    for (const q of placed) ext = Math.max(ext, Math.hypot(q.ox - c0, q.oy - c0), Math.hypot(q.ox + q.f.w - c0, q.oy + q.f.h - c0));
    const at = groupSpot(ext, how);
    for (const q of placed) SimGL.add(q.f, Math.round(at[1] - c0 + q.oy), Math.round(at[0] - c0 + q.ox));
    C.revive(); C.analyze(0);
    if (focusNew) {   // the ship: the first placed, at its form's middle
      const q = placed[0];
      const x = ((at[0] - c0 + q.ox + q.f.w / 2) % C.BW + C.BW) % C.BW, y = ((at[1] - c0 + q.oy + q.f.h / 2) % C.BH + C.BH) % C.BH;
      const i = C.found ? Creatures.pick(C.found, x, y, 2 * C.R, SimGL.BW, SimGL.BH, SimGL.BLOCK, C.BW, C.BH) : -1;
      if (i >= 0) C.setFocus(i);
    }
    C.hooks.change();
    return true;
  }
  // A spot for a group of radius ext (cells): the board's centre, or one clear of every creature
  // (random; far: of 30 such, the one farthest from the ship); the emptiest tried if none is clear
  function groupSpot(ext, how) {
    const W = C.BW, H = C.BH;
    if (how === 'center') return centerPt();
    if (how && how.place) return how.place(2 * ext, 2 * ext);   // (a placer: the scrolling world's edge)
    const comps = C.found ? C.found.comps : [], ship = how === 'far' ? C.ship() : null;
    let best = null, bestD = -1, found = 0, fallback = null, fallbackD = -1;
    for (let k = 0; k < 400 && found < (ship ? 30 : 1); k++) {
      const col = Math.random() * W, row = Math.random() * H;
      if (!spotOk(col, row, ext + C.R)) continue;   // (the scrolling world: away from the seam)
      const near = Math.min(Infinity, ...comps.map(c => Math.hypot(C.wrapD(c.col - col, W), C.wrapD(c.row - row, H)) - (c.radius || C.R)));
      if (near > fallbackD) { fallback = [col, row]; fallbackD = near; }
      if (near < ext + C.R) continue;
      found++;
      const d = ship ? Math.hypot(C.wrapD(ship.col - col, W), C.wrapD(ship.row - row, H)) : 0;
      if (d > bestD) { best = [col, row]; bestD = d; }
    }
    return best || fallback || [Math.random() * W, Math.random() * H];
  }

  // ---- chapter title card ([chapter]) ----------------------------------------------------------
  // The views fade out first (WORLD_FADE_MS), then the card fades in; closing it fades the card
  // out and the views back in together.
  const FADE_MS = 1500;   // as .chapter's transition in style.css
  let chapterOpen = false, cardShown = false, chapterTimer = null;
  const chapterShowing = () => chapterOpen || !$('chapter').hidden;
  // kind: 'chapter' ([chapter]) or 'challenge' ([challenge]: "<title> Challenge" above, <no> (its subtitle) as the title)
  function showChapter(no, title, kind = 'chapter') {
    // a challenge: "<name> Challenge" above (small), its subtitle (a number: "Stage <n>"; none: its name) below
    const stage = /^\d+$/.test(no) ? `Stage ${no}` : no;
    $('chapter-no').textContent = kind === 'challenge' ? `${title} Challenge` : no.toLowerCase() === 'final' ? 'FINAL CHAPTER' : `Chapter ${no}`;
    $('chapter-title').textContent = kind === 'challenge' ? stage || title : title;
    $('chapter-stage').hidden = true;
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
  // Q W E R T: Speed down, up, Upward (Undo), Zoom (Width) out, in; A S D: Respawn, Choose / Morphs /
  // Clone, Species / Challenges, and holding a clone F flips, G turns; Z X C V: Fire (Clear), Pulses,
  // the previous / next pulse (Settings: no key)
  const KEYS = { play: '\u21b5', sound: 'M', headup: 'E', morphs: 'S', reset: 'A', 'morphs-dead': 'S', fire: 'Z',
                 'book-close': 'Esc', 'bullet-prev': 'C', 'bullet-next': 'V', bullets: 'X', species: 'D', restart: '\u2190',
                 'dev-ch': 'I', 'dev-ok': 'O', 'dev-do': 'P',
                 'clone-rot': 'G', 'clone-flip': 'F', 'clone-ok': '\u21b5', 'clone-cancel': 'Esc' };
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
  Dots.beaconsDrift = MODE !== 'challenge';   // challenge mode: beacons stay put (the round dots still drift)
  Dots.hooks.score = () => {
    const d = Dots.score - dotsSeen;   // eaten since (not a fall: dots.js may restart its count)
    dotsSeen = Dots.score;
    if (d > 0) changeEnergy(d); else showEnergy();
  };
  function showEnergy() {
    $('score').hidden = energyHidden || (!Dots.on && !challenge.maskCost && !challenge.eraseCost && !challenge.drawCost && !challenge.worldDrawCost && !bulletCost && !cloneCost); $('score-v').textContent = Math.round(energy);
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
  // [bullet_cost:<n>]: each shot costs n energy (bulletCost, with the energy); with less, the gun doesn't fire
  const canPay = () => energy >= bulletCost;
  function doFire() {
    if (fieldMode()) return doClear();   // (Fire reads Clear)
    if (!canPay()) return;
    if (Bullets.fire()) { flash('fire'); if (bulletCost) changeEnergy(-bulletCost); }
  }
  $('fire').addEventListener('pointerdown', e => { e.preventDefault(); doFire(); });
  $('bullet-prev').addEventListener('pointerdown', e => { e.preventDefault(); if (Bullets.cycle(-1)) flash('bullet-prev'); });
  $('bullet-next').addEventListener('pointerdown', e => { e.preventDefault(); if (Bullets.cycle(1)) flash('bullet-next'); });
  // the pattern's name under "Fire"; < > above it while the list holds two or more (C / V)
  Bullets.hooks.change = () => {
    const many = Bullets.on && Bullets.count > 1, fm = fieldMode();
    $('fire').querySelector('.kt').textContent = fm ? 'Clear' : 'Fire';
    $('fire').hidden = !Bullets.on;
    $('fire-id').hidden = !Bullets.kind;
    $('fire-id').textContent = (Bullets.kind || '').replace(' ', '\n');   // two words: two rows
    $('bullet-nav').hidden = !many;
    syncButtons();   // the Bullets button
    syncNew();   // the list changed: as many more creatable?
    saveInventory();   // (a pattern added: into the central inventory)
    if (book.open && book.mode === 'bullet') { if (Bullets.on) renderBook(); else closeBook(); }
    $('fire').title = (fm ? 'Clear the selected colour\'s field (Z)' : 'Fire a pulse the way the ship moves (Z)') + (many ? '; C / V: the previous / next type' : '');
    if (!fm) Field.end();
  };
  C.hooks.frame = () => {
    const hookT0 = performance.now();
    devFps();
    updateBackground();
    // Respawn / Choose fade with the views ([new_scene], [chapter]), as "Dissolved" (core)
    const shade = worldShade(), dead = document.querySelector('.ship-dead');
    dead.style.opacity = 1 - shade;
    dead.style.visibility = shade > 0.99 ? 'hidden' : '';
    $('msg').classList.toggle('next', byContinue());
    if (Bullets.on) {   // Fire: its ring in the bullets' colour
      $('fire').classList.toggle('ready', fieldMode() || (Bullets.canFire() && canPay()));
      const rgb = `rgb(${Bullets.color.map(v => Math.round(v * 255)).join(',')})`;
      $('fire').style.setProperty('--fire', rgb); $('bullet-nav').style.setProperty('--fire', rgb);   // (< > too)
    }
    if (C.action === 2 && !drawOn()) C.action = 0;   // (draw no longer allowed: the species changed)
    syncToolNames();
    Field.tick();
    autospawnFrame();
    scrollFrame(); centerFrame(); Dots.frame(); Bullets.frame(); Hostile.frame(); trackBots(); botsFire(); autofireFrame(); checkMinMorph(); checkMinHostile(); checkDiscovery(); challengeFrame(); checkScript(); if ($('reset').hidden === (lostPanel() || flashing('reset'))) syncButtons();
    const sh = C.ship(); if (sh) lastShip = [sh.col, sh.row, C.viewAngle()];
    devRecFrame(hookT0); };

  // ---- controls ---------------------------------------------------------------------------
  // The tool in use, by name, in italic: the ship's (C.action) at the ship view's lower left while
  // a ship is there; drawing on the world (fieldMode), the brush (C.fieldTool) at the world view's
  // lower right instead; both in the selected channel's colour. Holding a clone: "Cloning", pink
  const SHIP_TOOLS = ['Masking', 'Erasing', 'Drawing'], FIELD_TOOLS = ['Paintbrush', 'Watercolor', 'Eraser'];
  function syncToolNames() {
    const fm = fieldMode(), st = $('ship-tool'), ft = $('field-tool');
    const cloning = !!clone.held && !fm;   // holding a clone (↻ ↔ ✓ ✗ shown): "Cloning", in the clone pink
    const shipText = cloning ? 'Cloning' : !fm && !shipLost() && !emptySpace ? SHIP_TOOLS[C.action] ?? '' : '', fieldText = fm ? FIELD_TOOLS[C.fieldTool] ?? '' : '';
    if (st.textContent !== shipText) st.textContent = shipText;
    if (ft.textContent !== fieldText) ft.textContent = fieldText;
    st.hidden = !shipText; ft.hidden = !fieldText;
    ft.style.color = C.CH_COLORS[C.channel] || '';
    st.style.color = cloning ? 'var(--clone)' : ft.style.color;
    st.style.opacity = ft.style.opacity = 1 - worldShade();   // (faded with the views: [new_scene], [chapter])
  }
  // A button that goes as it is used (Respawn, ✓, ✗) stays until its flash is over (Page.flash: 180 ms)
  const flashing = id => $(id).classList.contains('flash');
  function flashGo(id) { flash(id); setTimeout(syncButtons, 200); }
  function syncButtons() {
    $('play').setAttribute('aria-pressed', !C.running);   // "Pause", highlighted while paused
    $('reset').hidden = !lostPanel() && !flashing('reset');   // (flashing: kept until the flash is over, so it is seen)
    $('morphs-dead').hidden = !lostPanel() || !book.on;
    $('morphs').hidden = !book.on || shipLost();   // while the ship is lost, Choose takes its place
    if (clone.held && shipLost()) clone.held = null;   // lost while holding: the clone goes
    $('clone-tools').hidden = !clone.held && !flashing('clone-ok') && !flashing('clone-cancel');
    $('morphs').classList.toggle('clone-on', clone.on);   // pink, and "Clone", while clone mode is allowed
    $('morphs').textContent = clone.on ? 'Clone' : 'Morphs';
    $('morphs').title = clone.on ? 'Clone a discovered morph into the world: the morph book (S)'
                                 : 'The morph book: the morphs you have discovered (S)';
    $('clone-ok').hidden = !(clone.held && clone.held.at) && !flashing('clone-ok');   // ✓ once a touch has placed it
    $('bullets').classList.toggle('new', bulletNew);
    $('bullets').textContent = bulletNew ? 'Create' : 'Pulses';   // (the UI's name for the bullets)
    $('bullets').hidden = !Bullets.on;   // the bullet book's button: while the gun is on (the ship lost too)
    $('settings').hidden = !settingsOn;   // ([settings:on])
    $('species').hidden = MODE === 'story' || (MODE === 'sandbox' && !visSpecies.length);   // sandbox mode: the species book's button; challenge mode: the challenge book's
    // world draw (fieldMode): Upward is Undo (Redo once undone), Zoom is Width
    const fm = fieldMode();
    $('headup').textContent = fm ? (Field.undone ? 'Redo' : 'Undo') : 'Upward';
    $('headup').classList.toggle('dim', fm && !Field.canUndo);
    $('headup').title = fm ? 'Undo the last stroke or Clear; again: redo it (E)' : 'Upward: rotate the ship view so the creature\'s heading points up (E toggles)';
    $('headup').setAttribute('aria-pressed', !fm && C.viewMode !== 'board');
    $('zoom-name').textContent = fm ? 'Width' : 'Zoom';
    $('zoom').parentElement.title = fm ? `Brush width: ${Field.WIDTHS.join(' / ')} cells; the soft brush and the eraser twice that (R narrower, T wider)` : 'Zoom of the ship view, in screen px per cell (R out, T in)';
    $('speed').value = C.speedIdx + 1;
    $('zoom').value = fm ? Field.widthIdx + 2 : C.zoomIdx + 1;
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
    // challenge mode, failed (or completed): the challenge again, with the player's morph
    if (MODE === 'challenge' && (challenge.failed || (!challenge.cur && challenge.last))) return retryChallenge();
    const type = retryAction ?? 'clear_put';
    const name = playerMorph ?? C.morphName(C.pattern).replace(C.BOOM, '');
    if (type === 'book') openBook();
    else if (type === 'clear_put') { C.clearWorld(); ACTIONS.put_player(name, 'center'); C.setRunning(true); }
    else if (type === 'put') { ACTIONS.put_player(name, 'random'); C.setRunning(true); }
    if (type === 'clear_put' || type === 'put') { retake = true; cand = null; respawns++; respawned(); }   // [need:respawn]; the morph's picture may be retaken
  }
  const doView = () => { if (fieldMode()) return doUndo(); flash('headup'); C.toggleHeadUp(); };
  const doSpeed = d => { flash('speed'); C.changeSpeed(d); };
  // The game's zooms: the core's but the outermost two (2 and 16 px per cell): slots 1 .. 5
  const ZOOM_MIN = 1, ZOOM_MAX = C.ZOOMS.length - 2;
  const clampZoom = i => Math.max(ZOOM_MIN, Math.min(ZOOM_MAX, i));
  if (C.zoomIdx !== clampZoom(C.zoomIdx)) C.setZoom(clampZoom(C.zoomIdx));   // e.g. ?zoom=2
  const doZoom = d => {   // R out, T in (world draw: the brush's width)
    flash('zoom');
    if (fieldMode()) { Field.setWidth(Field.widthIdx + d); syncButtons(); } else C.setZoom(clampZoom(C.zoomIdx + d));
  };
  $('play').onclick = doPlay;
  // Retry (and Backspace); challenge mode: the retry book (a pick plays the challenge again)
  // Sandbox mode: Retry is the species chosen again (its [on_start] block, with the player's morph
  // for [put_player], useMorph), behind a fade, not the whole script replayed
  function retrySandbox() {
    if (resetting || speciesFading || chapterShowing() || bookShowing()) return;
    const id = chosenSpecies ?? (C.pattern && C.speciesOf(C.pattern));
    if (!id || !visSpecies.some(e => e.id === id)) return resetScene();
    if (clone.held) endClone();
    resetting = true;
    fadeWorld(1);
    setTimeout(() => {
      Dots.clear(); Bullets.clear(); Hostile.clear();
      useMorph = playerMorph;
      chooseSpecies(id);
      useMorph = null;
      fadeWorld(0);
      resetting = false;
    }, WORLD_FADE_MS);
  }
  const retryScene = () => MODE === 'sandbox' ? retrySandbox() : resetScene();
  const doRestart = () => {
    flash('restart');
    if (MODE === 'challenge' && (challenge.cur || challenge.last)) return openBook(null, true);   // the retry book
    // story, sandbox: the retry book too (a pick redoes the scene with that morph), if the morph
    // book is on and a morph of this species is discovered; else the scene again at once
    if (book.on && discovered.some(d => isVisible(d.name) && !isHidden(d.name)) && line >= 0 && !resetting) return openRetryBook();
    retryScene();
  };
  $('restart').onclick = doRestart;
  // Respawn (A): the button flashes first, then the respawn (RESPAWN_DELAY_MS later; a second press meanwhile: nothing)
  const RESPAWN_DELAY_MS = 180;
  let respawnPending = false;
  $('reset').onclick = () => {
    if (respawnPending) return;
    respawnPending = true; flashGo('reset');
    setTimeout(() => { respawnPending = false; if (lostShown()) doReset(); }, RESPAWN_DELAY_MS);
  };
  $('morphs-dead').onclick = () => { flash('morphs-dead'); openBook(); };
  $('headup').onclick = doView;
  $('speed').oninput = e => C.setSpeed(+e.target.value - 1);
  $('zoom').value = C.zoomIdx + 1;
  $('zoom').oninput = e => { if (fieldMode()) { Field.setWidth(+e.target.value - 2); syncButtons(); } else C.setZoom(clampZoom(+e.target.value - 1)); };

  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'SELECT' || (e.target.tagName === 'INPUT' && e.target.type !== 'range')) return;
    const k = e.key.length === 1 ? e.key.toLowerCase() : e.key;
    if (book.open) {   // the book takes the keys while open; Space does nothing
      // a slot's letter (SLOT_KEYS) picks it, if it can be picked (which closes the book too)
      const el = k.length === 1 && $('book-grid').querySelector(`.book-slot.pick[data-key="${k.toUpperCase()}"]`);
      if (el) el.click();
      // Enter: the current slot (the lighter one: the player's morph, the bullet in use, the current
      // species, the challenge being played), as a click on it; none: nothing (not even the focused Close)
      else if (k === 'Enter') { e.preventDefault(); $('book-grid').querySelector('.book-slot.pick.current')?.click(); }
      // else Escape (only: the letters that open the books pick slots in them) closes it
      else if (k === 'Escape') closeOrBack();
      else if (e.code === 'Space') e.preventDefault();   // not even the focused Close button
      return;
    }
    if (clone.held) {   // holding a clone: F flips, G turns, Enter pastes, Escape cancels
      if (k === 'f') { flash('clone-flip'); return flipClone(); }
      if (k === 'g') { flash('clone-rot'); return rotateClone(); }
      if (k === 'Enter') { e.preventDefault(); if (!$('clone-ok').hidden) flashGo('clone-ok'); return confirmClone(); }
      if (k === 'Escape') { flashGo('clone-cancel'); return endClone(); }
    }
    if (e.code === 'Space') { e.preventDefault(); playerNext(); }
    else if (k === 'o') { if (DEV) flash('dev-ok'); devSkip(); }
    else if (k === 'p') { if (DEV) { flash('dev-do'); Dots.eat(); } }
    else if (k === 'i') { if (DEV) flash('dev-ch'); devChapter(); }
    else if (k === 'l') { if (DEV) devRecord(); }
    else if (k === 'Enter') { e.preventDefault(); doPlay(); }
    // (the keys: see KEYS)
    else if (k === 'a') { if (lostShown()) $('reset').onclick(); }   // only while Respawn shows
    else if (k === 's') { if (book.on) { const b = ['morphs', 'morphs-dead'].find(id => !$(id).hidden); b ? $(b).onclick() : openBook(); } }   // the morph book, as the Morphs / Choose buttons
    else if (k === 'd') { if (!$('species').hidden) $('species').onclick(); }   // the species book (sandbox mode), the challenge book (challenge mode)
    else if (k === 'q') doSpeed(-1);   // slower
    else if (k === 'w') doSpeed(1);    // faster
    else if (k === 'e') doView();
    else if (k === 'r') doZoom(-1);
    else if (k === 't') doZoom(1);
    else if (k === 'z') doFire();
    else if (k === 'x') { if (!$('bullets').hidden) $('bullets').onclick(); }   // the bullet book, as the Bullets button
    else if (k === 'c') { if (Bullets.cycle(-1)) flash('bullet-prev'); }
    else if (k === 'v') { if (Bullets.cycle(1)) flash('bullet-next'); }
    else if (k === 'Backspace') { e.preventDefault(); doRestart(); }
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
    else if (MODE !== 'story') loadCentralInventory();
    SCRIPT = parseScript(txt);
    fitMsg();
    // challenge mode: the set-up lines' actions (before the first challenge), then the challenge book
    if (MODE === 'challenge') {
      parseChallenges();
      const first = challenge.list.length ? challenge.list[0].start : SCRIPT.length;
      for (const L of SCRIPT.slice(0, first)) for (const a of L.actions) ACTIONS[a.name](...a.args);
      baseHidden = new Set(hiddenMorphs);
      if (!challenge.list.length) return say('No challenges in challenge.txt.');
      return openChallengeBook();
    }
    // [jump_here] (for writing the script): start at the first line that has it. Every action on
    // the lines before it runs first, in order, so the set-up is as if played up to there; their
    // messages and conditions are skipped, and so are the chapter cards.
    const from = Math.max(0, SCRIPT.findIndex(L => L.actions.some(a => a.name === 'jump_here')));
    replaying = true;
    for (const L of SCRIPT.slice(0, from))
      for (const a of L.actions) if (a.name !== 'chapter') ACTIONS[a.name](...a.args);
    replaying = false;
    if (from < SCRIPT.length) run(from);
  }, err => { console.error(err); say(`Could not load the script (${SCRIPT_FILE}).`); }), Page.showError)
    .finally(() => { Page.releaseMusic(); if (worldFade === START_FADE) fadeWorld(0); });   // (the views fade in: START_FADE)   // (the music was held until the script set its playlist: page.js; with no [music], the default)
})();
