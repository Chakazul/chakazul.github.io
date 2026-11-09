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

  // Background music, royalty free: a playlist of TRACKS (setPlaylist; the first alone at first),
  // played in turn and looped, with 10 s of silence between plays. The
  // sound button (the lab's Stars, the game's Music) turns the music on or off; in the lab the
  // star field follows it too. Off by default (then as last chosen); browsers
  // only let audio start after a user gesture, so it begins on the first click, tap or key.
  // Each track is fetched only when it starts to play (no src, preload none, until then).
  const TRACKS = ['assets/music/rexlambo-long-way.mp3',                                 // "Long Way" by Rexlambo
                  'assets/music/arthur-vyncke-a-few-jumps-away.mp3',                     // "A Few Jumps Away" by Arthur Vyncke
                  'assets/music/chris-zabriskie-out-of-the-skies-under-the-earth.mp3',  // "Out of the Skies, Under the Earth" by Chris Zabriskie
                  'assets/music/stellardrone-billions-and-billions.mp3'];                // "Billions and Billions" by Stellardrone
  const TRACK_GAP_MS = 10000;
  const music = new Audio();
  const VOLUME = 0.5, MUSIC_FADE_MS = 1500;   // a switch to another track fades the old one out first
  music.volume = VOLUME;
  music.preload = 'none';
  let track = 0, gapTimer = null, loaded = -1;   // loaded: the track in music.src
  let playlist = [0], pos = 0;   // the tracks played in turn, looped (setPlaylist); pos: track's place in it
  // Play the current track, giving it its file first if it has not got it yet.
  let fadeTimer = null;
  function play() {
    clearInterval(fadeTimer); fadeTimer = null; music.volume = VOLUME;
    if (loaded !== track) { music.src = TRACKS[track]; loaded = track; }
    music.play().catch(() => {});
  }
  music.addEventListener('ended', () => {   // the playlist's next track (or the same one), after the gap
    pos = (pos + 1) % playlist.length; track = playlist[pos];   // its file is fetched when it starts (pos -1: the first)
    gapTimer = setTimeout(() => { gapTimer = null; if (soundOn) play(); }, TRACK_GAP_MS);
  });
  // On or off as last chosen in this browser (localStorage, one key per page: the game's Music,
  // the lab's Stars); never chosen: the button's data-default ("on"), else off (both pages now).
  const SOUND_KEY = `cognoship-sound-${$('sound')?.dataset.label ? 'game' : 'lab'}`;
  let soundOn = $('sound')?.dataset.default === 'on';
  try { const v = localStorage.getItem(SOUND_KEY); if (v) soundOn = v === 'on'; } catch (e) {}
  // The button: the game's Music (data-label, a note icon in the title row), else the lab's Stars.
  function syncSound() {
    if ($('sound')) {
      // data-label (the game: Music, an icon button): unchanged, highlighted when on; else "Stars: ON / OFF"
      if (!$('sound').dataset.label) $('sound').textContent = `Stars: ${soundOn ? 'ON' : 'OFF'}`;
      $('sound').setAttribute('aria-pressed', soundOn);
    }
    if (!soundOn) music.pause();
    else if (!gapTimer) play();   // in a gap between tracks, the timer starts the next
  }
  function toggleSound() {
    if ($('sound')) flash('sound');
    soundOn = !soundOn; syncSound();
    try { localStorage.setItem(SOUND_KEY, soundOn ? 'on' : 'off'); } catch (e) {}
  }
  if ($('sound')) $('sound').onclick = toggleSound;
  // Skip to the next track now (also from a gap between tracks); while off, it plays when turned on.
  function nextTrack() { switchTo((track + 1) % TRACKS.length); }
  // Switch to track i now: the one playing fades out (MUSIC_FADE_MS), then i starts; while off, i is
  // the one that plays when turned on.
  function switchTo(i) {
    clearTimeout(gapTimer); gapTimer = null;
    track = i;
    if (!soundOn) { music.pause(); return; }
    if (music.paused) { play(); return; }
    if (fadeTimer) return;   // already fading out: the new track starts when it ends
    const v0 = music.volume, t0 = performance.now();
    fadeTimer = setInterval(() => {
      const k = (performance.now() - t0) / MUSIC_FADE_MS;
      if (k < 1) { music.volume = v0 * (1 - k); return; }
      clearInterval(fadeTimer); fadeTimer = null;
      music.pause();
      if (soundOn) play();
    }, 40);
  }
  const startMusic = () => { if (soundOn && music.paused && !gapTimer) play(); };
  // retried on every press until it plays; on window, capturing, and added before the pages' own
  // handlers, so it sees every press first, even one a handler swallows (the game's chapter card,
  // closed by Space or a click, stops the event there)
  for (const ev of ['pointerdown', 'keydown']) window.addEventListener(ev, startMusic, { capture: true });
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

  // The playlist becomes `list` (0-based track numbers; the game's [music:<list>:<play>] gives
  // n - 1), played in turn and looped. start: switch to that track now (from its beginning, unless
  // it is the one playing; while the music is off, it is the one that plays when turned on). No
  // start: the track playing now carries on, and the list follows it (from the track after it if
  // it is in the list, else from the list's first).
  function setPlaylist(list, start) {
    list = list.filter(i => i >= 0 && i < TRACKS.length);
    if (!list.length) return false;
    playlist = list;
    if (start !== undefined && start >= 0 && start < TRACKS.length && start !== track) switchTo(start);
    pos = list.indexOf(track);   // -1 (not in it): the next is the list's first
    return true;
  }
  const setTrack = i => setPlaylist([i], i);
  return { $, flash, fillRange, showError, toggleSound, nextTrack, setTrack, setPlaylist, TRACKS, get soundOn() { return soundOn; } };
})();
