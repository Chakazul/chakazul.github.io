"use strict";
// ============================================================================================
//  Generative ambient music (a draft, tried on ambient.html; not wired into the game yet): an
//  endless piece made live with the Web Audio API, nothing downloaded. Four moods, each set up from
//  an analysis of one of the game's tracks (key, mode, tempo, brightness, the chords it dwells on;
//  ffmpeg + numpy), so it sounds as if it belongs with them without using their tunes.
//
//  Layers, each into a shared reverb (a convolver on a few seconds of decaying noise) and a delay:
//    pad    the chord, each note two detuned sawtooths and a triangle an octave down, through a
//           low-pass filter that opens and closes over the chord, swelling in and out over seconds
//    drone  the chord's root two octaves down, a soft sine
//    bells  now and then a note of the mode's pentatonic subset, an FM bell (sine carrier, sine
//           modulator at a ratio, its depth decaying), echoed by the delay
//  The chords walk (seeded random) along the mood's progression, each held for `bars` beats; a
//  scheduler queues the notes ahead on the AudioContext's clock, so frame drops don't matter.
// ============================================================================================

const Ambient = (function () {
  const MINOR = [0, 2, 3, 5, 7, 8, 10], MAJOR = [0, 2, 4, 5, 7, 9, 11];
  // degrees: the chords (scale degrees, 0 = the tonic), the walk going from one to any listed after
  // it, weighted. root: the tonic, MIDI (C#3 = 49). bpm: beats per minute. bright: the pads' filter
  // at its widest, Hz. bells: chance of a bell per beat. bars: beats per chord.
  const MOODS = [
    { name: 'Long Way (C# minor)', root: 49, scale: MINOR, bpm: 46, bright: 1300, bells: 0.28, bars: 8,
      degrees: { 0: [5, 3, 2], 5: [3, 0, 2], 3: [0, 5, 6], 2: [5, 0], 6: [0, 3] } },
    { name: 'A Few Jumps Away (C minor)', root: 48, scale: MINOR, bpm: 46, bright: 950, bells: 0.2, bars: 8,
      degrees: { 0: [2, 5, 6], 2: [5, 0], 5: [3, 0, 6], 3: [0, 6], 6: [0, 2] } },
    { name: 'Out of the Skies (Db major)', root: 49, scale: MAJOR, bpm: 54, bright: 900, bells: 0.22, bars: 8,
      degrees: { 0: [5, 3, 4], 5: [3, 0], 3: [4, 0], 4: [0, 5] } },
    { name: 'Billions and Billions (Bb minor)', root: 46, scale: MINOR, bpm: 65, bright: 2600, bells: 0.42, bars: 12,
      degrees: { 0: [6, 2, 3], 6: [0, 2], 2: [3, 0], 3: [0, 6] } },
  ];

  let ctx = null, out = null, padBus = null, bellBus = null, reverb = null, delay = null;
  let mood = null, rand = Math.random, timer = null, nextChordAt = 0, nextBeat = 0, degree = 0, playing = false;
  const LOOKAHEAD = 0.6;   // s queued ahead
  const hz = m => 440 * Math.pow(2, (m - 69) / 12);

  // A small seeded random (mulberry32), so a mood plays the same way each time it starts
  function seeded(a) {
    return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  }

  function setup() {
    ctx = new (window.AudioContext || window.webkitAudioContext)();
    const comp = ctx.createDynamicsCompressor();
    comp.threshold.value = -18; comp.ratio.value = 3;
    out = ctx.createGain(); out.gain.value = 0;
    out.connect(comp).connect(ctx.destination);
    // reverb: 6 s of stereo noise, decaying
    reverb = ctx.createConvolver();
    const len = Math.round(ctx.sampleRate * 6), ir = ctx.createBuffer(2, len, ctx.sampleRate);
    for (let c = 0; c < 2; c++) {
      const d = ir.getChannelData(c);
      for (let i = 0; i < len; i++) d[i] = (Math.random() * 2 - 1) * Math.pow(1 - i / len, 3.2);
    }
    reverb.buffer = ir;
    const wet = ctx.createGain(); wet.gain.value = 0.55;
    reverb.connect(wet).connect(out);
    // delay: echoes, darker each time
    delay = ctx.createDelay(4);
    const fb = ctx.createGain(), tone = ctx.createBiquadFilter(), dsend = ctx.createGain();
    fb.gain.value = 0.38; tone.type = 'lowpass'; tone.frequency.value = 2200; dsend.gain.value = 0.5;
    delay.connect(tone).connect(fb).connect(delay);
    tone.connect(dsend).connect(out); tone.connect(reverb);
    padBus = ctx.createGain(); padBus.gain.value = 0.5;
    padBus.connect(out); padBus.connect(reverb);
    bellBus = ctx.createGain(); bellBus.gain.value = 0.35;
    bellBus.connect(out); bellBus.connect(reverb); bellBus.connect(delay);
  }

  // The chord on scale degree d, voiced near MIDI `centre`: four notes (a triad and the ninth)
  function chord(d, centre = 57) {
    const sc = mood.scale, note = k => mood.root + sc[(d + k) % 7] + 12 * Math.floor((d + k) / 7);
    return [0, 2, 4, 8].map(k => { let m = note(k); while (m < centre - 8) m += 12; while (m > centre + 6) m -= 12; return m; });
  }

  function pad(m, t, dur) {
    const att = Math.min(4, dur * 0.35), rel = 5, level = 0.07;
    const g = ctx.createGain(), f = ctx.createBiquadFilter(), pan = ctx.createStereoPanner();
    g.gain.setValueAtTime(0, t);
    g.gain.linearRampToValueAtTime(level, t + att);
    g.gain.setValueAtTime(level, t + dur);
    g.gain.linearRampToValueAtTime(0, t + dur + rel);
    f.type = 'lowpass'; f.Q.value = 0.8;
    f.frequency.setValueAtTime(mood.bright * 0.35, t);
    f.frequency.linearRampToValueAtTime(mood.bright * (0.7 + 0.3 * rand()), t + dur * 0.5);
    f.frequency.linearRampToValueAtTime(mood.bright * 0.3, t + dur + rel);
    pan.pan.value = (rand() * 2 - 1) * 0.45;
    f.connect(g).connect(pan).connect(padBus);
    const oscs = [['sawtooth', 0, -7], ['sawtooth', 0, 7], ['triangle', -12, 0]].map(([type, oct, cents]) => {
      const o = ctx.createOscillator();
      o.type = type; o.frequency.value = hz(m + oct); o.detune.value = cents + (rand() - 0.5) * 4;
      o.connect(f); o.start(t); o.stop(t + dur + rel + 0.1);
      return o;
    });
    return oscs;
  }

  function drone(m, t, dur) {
    const g = ctx.createGain(), o = ctx.createOscillator();
    o.type = 'sine'; o.frequency.value = hz(m);
    g.gain.setValueAtTime(0, t);
    g.gain.linearRampToValueAtTime(0.12, t + 4);
    g.gain.setValueAtTime(0.12, t + dur);
    g.gain.linearRampToValueAtTime(0, t + dur + 4);
    o.connect(g).connect(padBus);
    o.start(t); o.stop(t + dur + 4.1);
  }

  function bell(m, t) {
    const car = ctx.createOscillator(), mod = ctx.createOscillator(), mg = ctx.createGain(), g = ctx.createGain(), pan = ctx.createStereoPanner();
    const f = hz(m), ratio = rand() < 0.5 ? 2 : 3.5, decay = 3 + rand() * 3;
    car.type = mod.type = 'sine';
    car.frequency.value = f; mod.frequency.value = f * ratio;
    mg.gain.setValueAtTime(f * 1.6, t);
    mg.gain.exponentialRampToValueAtTime(f * 0.05, t + decay * 0.6);
    g.gain.setValueAtTime(0.0001, t);
    g.gain.exponentialRampToValueAtTime(0.16, t + 0.015);
    g.gain.exponentialRampToValueAtTime(0.0001, t + decay);
    pan.pan.value = (rand() * 2 - 1) * 0.7;
    mod.connect(mg).connect(car.frequency);
    car.connect(g).connect(pan).connect(bellBus);
    car.start(t); mod.start(t); car.stop(t + decay + 0.1); mod.stop(t + decay + 0.1);
  }

  // The next chord on the walk: one of those listed after the current one, the first ones likelier
  function step() {
    const next = mood.degrees[degree] || [0];
    const w = next.map((_, i) => 1 / (i + 1)), sum = w.reduce((a, b) => a + b, 0);
    let r = rand() * sum, i = 0;
    while (r > w[i]) { r -= w[i]; i++; }
    degree = next[Math.min(i, next.length - 1)];
  }

  function schedule() {
    const beat = 60 / mood.bpm, until = ctx.currentTime + LOOKAHEAD;
    while (nextChordAt < until) {
      const dur = beat * mood.bars, notes = chord(degree);
      notes.forEach(m => pad(m, nextChordAt, dur));
      drone(mood.root + mood.scale[degree] - 12, nextChordAt, dur);   // the chord's root, an octave below the tonic's octave
      nextChordAt += dur;
      step();
    }
    while (nextBeat < until) {
      if (rand() < mood.bells) {
        // a pentatonic note (degrees 0 1 2 4 5 of the mode) in the upper octaves
        const pent = [0, 1, 2, 4, 5], k = pent[Math.floor(rand() * pent.length)];
        bell(mood.root + mood.scale[k] + (rand() < 0.6 ? 24 : 36), nextBeat + (rand() < 0.3 ? beat / 2 : 0));
      }
      nextBeat += beat;
    }
  }

  // Start mood i (stopping any other); call from a click or key (browsers only start audio then)
  function start(i = 0, seed = 1) {
    if (!ctx) setup();
    ctx.resume();
    stop(true);
    mood = MOODS[i]; rand = seeded(seed * 7919 + i); degree = 0;
    nextChordAt = nextBeat = ctx.currentTime + 0.1;
    out.gain.cancelScheduledValues(ctx.currentTime);
    out.gain.setValueAtTime(out.gain.value, ctx.currentTime);
    out.gain.linearRampToValueAtTime(0.8, ctx.currentTime + 3);
    playing = true;
    schedule();
    timer = setInterval(schedule, 150);
  }
  // Fade out over 3 s (quick: at once, for a switch; the notes already queued ring out)
  function stop(quick = false) {
    clearInterval(timer); timer = null;
    if (!ctx || !playing) return;
    playing = false;
    const now = ctx.currentTime;
    out.gain.cancelScheduledValues(now);
    out.gain.setValueAtTime(out.gain.value, now);
    out.gain.linearRampToValueAtTime(0, now + (quick ? 0.3 : 3));
  }
  function volume(v) { if (out) out.gain.value = v; }

  return { MOODS, start, stop, volume, get playing() { return playing; } };
})();
