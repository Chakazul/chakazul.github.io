"use strict";
// ============================================================================================
//  Cognoship page chrome, shared by the game (index.html) and the lab (lab.html): DOM helpers,
//  the hand-drawn slider fill, the light / dark theme, fullscreen, and the background music
//  (the Stars button). It wires whichever of these elements the page has, by id: theme,
//  fullscreen, fs-exit, sound.
// ============================================================================================

const Page = (function () {
  const $ = id => document.getElementById(id);

  // Briefly highlight the controls an action used (a click or its shortcut key).
  function flash(...ids) {
    for (const id of ids) {
      const el = $(id).tagName === 'INPUT' ? $(id).parentElement : $(id);   // a slider: its label
      el.classList.remove('flash'); void el.offsetWidth; el.classList.add('flash');
      clearTimeout(el._flash); el._flash = setTimeout(() => el.classList.remove('flash'), 180);
    }
  }
  // The filled part of a slider bar (CSS --pct), since the bar is drawn by hand for the theme.
  function fillRange(el) { el.style.setProperty('--pct', 100 * (el.value - el.min) / (el.max - el.min) + '%'); }
  document.addEventListener('input', e => { if (e.target.type === 'range') fillRange(e.target); });

  // Background music, royalty free: a playlist looped with 10 s of silence between tracks. The
  // sound button (the lab's Stars, the game's Music) turns the music on or off; in the lab the
  // star field follows it too. Off by default; browsers
  // only let audio start after a user gesture, so it begins on the first click, tap or key.
  const TRACKS = ['assets/music/rexlambo-long-way.mp3',                                 // "Long Way" by Rexlambo
                  'assets/music/arthur-vyncke-a-few-jumps-away.mp3',                     // "A Few Jumps Away" by Arthur Vyncke
                  'assets/music/chris-zabriskie-out-of-the-skies-under-the-earth.mp3',  // "Out of the Skies, Under the Earth" by Chris Zabriskie
                  'assets/music/stellardrone-billions-and-billions.mp3'];                // "Billions and Billions" by Stellardrone
  const TRACK_GAP_MS = 10000;
  const music = new Audio();
  music.volume = 0.5;
  let track = 0, gapTimer = null;
  music.src = TRACKS[track];
  music.addEventListener('ended', () => {
    track = (track + 1) % TRACKS.length;
    music.src = TRACKS[track];
    gapTimer = setTimeout(() => { gapTimer = null; if (soundOn) music.play().catch(() => {}); }, TRACK_GAP_MS);
  });
  // On or off at each page load: the button's data-default ("on" in the game), else off (the lab).
  // The choice is not remembered.
  let soundOn = $('sound')?.dataset.default === 'on';
  // The button's label: data-label (the game: Music), else Stars.
  function syncSound() {
    if ($('sound')) {
      // data-label (the game: Music): the label alone, highlighted when on; else "Stars: ON / OFF"
      const label = $('sound').dataset.label;
      $('sound').textContent = label || `Stars: ${soundOn ? 'ON' : 'OFF'}`;
      $('sound').setAttribute('aria-pressed', soundOn);
    }
    if (!soundOn) music.pause();
    else if (!gapTimer) music.play().catch(() => {});   // in a gap between tracks, the timer starts the next
  }
  function toggleSound() { if ($('sound')) flash('sound'); soundOn = !soundOn; syncSound(); }
  if ($('sound')) $('sound').onclick = toggleSound;
  // Skip to the next track now (also from a gap between tracks); while off, it plays when turned on.
  function nextTrack() {
    clearTimeout(gapTimer); gapTimer = null;
    track = (track + 1) % TRACKS.length;
    music.src = TRACKS[track];
    if (soundOn) music.play().catch(() => {});
  }
  const startMusic = () => { if (soundOn && music.paused && !gapTimer) music.play().catch(() => {}); };
  for (const ev of ['pointerdown', 'keydown']) document.addEventListener(ev, startMusic, { capture: true });   // retried until it plays
  syncSound();

  // theme and fullscreen
  if ($('theme')) $('theme').onclick = () => {
    const light = document.documentElement.getAttribute('data-theme') !== 'light';
    if (light) document.documentElement.setAttribute('data-theme', 'light');
    else document.documentElement.removeAttribute('data-theme');
    try { localStorage.setItem('cognoship-theme', light ? 'light' : 'dark'); } catch (e) {}
  };
  // In fullscreen the two control rows stay below the canvas; --rows is their height, so the
  // canvas takes what is left.
  // the game: one row (.below: messages and controls); the lab: the two control rows
  const FS_ROWS = document.querySelector('.below') ? ['.below'] : ['.controls', '.row2'];
  function fitFs() {
    if (!document.documentElement.classList.contains('fs')) return;
    const rows = FS_ROWS.map(q => document.querySelector(q)).filter(Boolean);
    const h = rows.reduce((a, el) => a + el.offsetHeight + parseFloat(getComputedStyle(el).marginTop), 0) + 8;
    document.documentElement.style.setProperty('--rows', h + 'px');
  }
  function setFs(on) { document.documentElement.classList.toggle('fs', on); fitFs(); }
  window.addEventListener('resize', fitFs);
  if (window.ResizeObserver) {   // the rows can re-wrap when the canvas width changes
    const ro = new ResizeObserver(fitFs);
    FS_ROWS.map(q => document.querySelector(q)).filter(Boolean).forEach(el => ro.observe(el));
  }
  if ($('fullscreen')) $('fullscreen').onclick = async () => {
    setFs(true);
    try { await document.documentElement.requestFullscreen?.(); } catch (e) {}
    try { await screen.orientation?.lock?.('landscape'); } catch (e) {}
  };
  if ($('fs-exit')) $('fs-exit').onclick = () => { setFs(false); if (document.fullscreenElement) document.exitFullscreen(); };
  document.addEventListener('fullscreenchange', () => { if (!document.fullscreenElement) setFs(false); });

  // A start-up failure, shown over the canvas (the #error box).
  function showError(err) {
    const box = $('error');
    box.hidden = false;
    box.textContent = 'Could not start: ' + err.message +
      (location.protocol === 'file:' ? ' (serve this folder over http, e.g. python3 -m http.server)' : '');
    console.error(err);
  }

  return { $, flash, fillRange, showError, toggleSound, nextTrack, get soundOn() { return soundOn; } };
})();
