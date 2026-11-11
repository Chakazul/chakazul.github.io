"use strict";
// ============================================================================================
//  FFT mode (?fft=1): the creatures' convolution by FFT on the GPU, in place of sim.glsl's direct
//  sum over ~1,000 kernel taps per cell, whose cost grows with R^2. A plug-in: with ?fft=1 it
//  wraps SimGL.init and, once the engine is up, sets SimGL.internals.ext.step / .rule; without
//  it, this file does nothing. To remove FFT mode, delete this file, its <script> tags and
//  shaders/fft*.glsl (glsim.js keeps two small hooks, ext and internals, which do nothing then).
//
//  Per step (step):
//    fftpack   the fields as complex pairs: (R, G), (B, b), and where the bullets mask, (R b, G b),
//              (B b) (b the bullets field, times the channels they act on)
//    fft x 2   forward FFTs of those (4 while the bullets mask)
//    fftmul    per pair of kernels: the spectrum of the channel each reads (less its bullets
//              part) times the kernel's spectrum; and the bullets' spectrum times the kernel's
//    fft x 8   inverse FFTs: U for every kernel (8 more for the bullets' hidden shares; gathered
//              by fftgather into 4 textures)
//    fftgrow   the update, as sim.glsl, with the cognitive mask's small circle taken out of U
//              directly (its cells only), and erase / draw / the bullets' erase after it
//  The kernels' spectra are made once per rule (ext.rule): each kernel as an image of the board,
//  flipped (sim.glsl's U(p) = sum w(o) A(p + o) is a correlation), forward FFT'd.
//
//  The FFT: mixed radix (4, 3, 2) Stockham passes, gather form (fft.glsl), along x (384 =
//  4 4 4 2 3) then y (256 = 4 4 4 4, one 4 split into 2 2 so the pass count is even and a
//  transform ends where it started). Twiddles from a table made on the CPU in double precision.
//  Checked on the CPU first: the passes against a direct 2D DFT (relative error 3e-13) and the
//  convolution against the taps (error 6e-17), in double precision; the GPU runs float32.
// ============================================================================================

const FFTSim = (function () {
  const ON = new URLSearchParams(location.search).get('fft') === '1';
  if (!ON) return { on: false };

  const NAMES = ['fft', 'fftpack', 'fftmul', 'fftgather', 'fftgrow'];
  const PAIRS = 8;   // kernel slots (16) in pairs, two complex numbers per texel
  let I = null, gl = null, W = 0, H = 0;
  const prog = {};
  let passes = [];   // { axis, N, R, Ns }
  let twTex = null, scratch = null;
  let F = [], Kspec = [], P = [], D = [];   // fields' spectra (4), kernels' spectra, U, bullets' shares (8 each)
  const fbo = new Map();   // texture -> its own framebuffer (an FFT pass writes one texture)
  let packFbo = null, mulFbo = [], kimgTex = null, KR = 0, pairsUsed = PAIRS, ready = false;

  // Radices for a length: 4s, then 2s, then 3s (384 = 4 4 4 2 3, 256 = 4 4 4 4)
  function radices(N) {
    const r = [];
    while (N % 4 === 0) { r.push(4); N /= 4; }
    while (N % 2 === 0) { r.push(2); N /= 2; }
    while (N % 3 === 0) { r.push(3); N /= 3; }
    if (N !== 1) throw new Error('FFT mode: board side not a product of 2s and 3s');
    return r;
  }
  function plan() {
    const px = radices(W), py = radices(H);
    if ((px.length + py.length) % 2) {   // an even count: a transform ends in its own texture
      const i = py.indexOf(4), j = px.indexOf(4);
      if (i >= 0) py.splice(i, 1, 2, 2); else if (j >= 0) px.splice(j, 1, 2, 2); else py.push(1);
    }
    passes = [];
    let Ns = 1;
    for (const R of px) { passes.push({ axis: 0, N: W, R, Ns }); Ns *= R; }
    Ns = 1;
    for (const R of py) { passes.push({ axis: 1, N: H, R, Ns }); Ns *= R; }
  }

  const tex = () => { const t = I.makeTex(W, H); fbo.set(t, I.makeFbo([t])); return t; };

  async function attach() {
    I = SimGL.internals; gl = I.gl; W = I.W; H = I.H;
    const src = await Promise.all(NAMES.map(async n => {
      const r = await fetch('./shaders/' + n + '.glsl');
      if (!r.ok) throw new Error('shader ' + n + ': HTTP ' + r.status);
      return r.text();
    }));
    NAMES.forEach((n, i) => { prog[n] = I.link(src[i]); });
    plan();
    // twiddles: row 0 exp(-2 pi i m / W) for m < W, row 1 the same over H
    const tw = new Float32Array(Math.max(W, H) * 2 * 4);
    for (let m = 0; m < W; m++) { tw[m * 4] = Math.cos(-2 * Math.PI * m / W); tw[m * 4 + 1] = Math.sin(-2 * Math.PI * m / W); }
    for (let m = 0; m < H; m++) { const i = (Math.max(W, H) + m) * 4; tw[i] = Math.cos(-2 * Math.PI * m / H); tw[i + 1] = Math.sin(-2 * Math.PI * m / H); }
    twTex = I.makeTex(Math.max(W, H), 2, tw);
    scratch = tex();
    F = [0, 1, 2, 3].map(tex);
    Kspec = Array.from({ length: PAIRS }, tex);
    P = Array.from({ length: PAIRS }, tex);
    D = Array.from({ length: PAIRS }, tex);
    packFbo = I.makeFbo(F);
    mulFbo = P.map((t, j) => I.makeFbo([t, D[j]]));
    I.ext.rule = buildKernels;
    if (I.taps) buildKernels(I.taps);
    I.ext.step = step;
    ready = true;
    console.info(`Cognoship: FFT mode on (${passes.map(p => p.R).join(' ')} passes, ${W} x ${H})`);
  }

  // One FFT of t (both complex pairs), in place by way of the scratch texture. sign -1 forward,
  // +1 inverse (unscaled: fftmul scales by 1 / (W H)).
  function fft(t, sign) {
    const p = prog.fft;
    gl.useProgram(p);
    gl.viewport(0, 0, W, H);
    gl.uniform1f(I.loc(p, 'uSign'), sign);
    I.bindTex(1, twTex, p, 'uTw');
    let a = t, b = scratch;
    for (const ps of passes) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.get(b));
      I.bindTex(0, a, p, 'uIn');
      gl.uniform1i(I.loc(p, 'uAxis'), ps.axis);
      gl.uniform1i(I.loc(p, 'uN'), ps.N);
      gl.uniform1i(I.loc(p, 'uR'), ps.R);
      gl.uniform1i(I.loc(p, 'uNs'), ps.Ns);
      I.draw4();
      [a, b] = [b, a];
    }
  }

  // The rule's kernels: their spectra (pairs), and the kernels as one small image for fftgrow
  function buildKernels(taps) {
    const SL = Rule.SLOTS;
    KR = 0;
    for (const [dx, dy] of taps.offs) KR = Math.max(KR, Math.abs(dx), Math.abs(dy));
    const side = 2 * KR + 1, kw = 4 * side;
    const kimg = new Float32Array(kw * side * 4);
    const imgs = Array.from({ length: PAIRS }, () => new Float32Array(W * H * 4));
    let used = 0;
    taps.offs.forEach(([dx, dy], t) => {
      const x = ((-dx) % W + W) % W, y = ((-dy) % H + H) % H;
      for (let k = 0; k < SL; k++) {
        const w = taps.weights[t * SL + k];
        if (!w) continue;
        used = Math.max(used, k + 1);
        imgs[k >> 1][(y * W + x) * 4 + (k & 1) * 2] += w;
        const g = k >> 2, i = k & 3;
        kimg[((dy + KR) * kw + g * side + dx + KR) * 4 + i] = w;
      }
    });
    pairsUsed = Math.ceil(used / 2);
    if (kimgTex) gl.deleteTexture(kimgTex);
    kimgTex = I.makeTex(kw, side, kimg);
    imgs.forEach((img, j) => {
      gl.bindTexture(gl.TEXTURE_2D, Kspec[j]);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, W, H, gl.RGBA, gl.FLOAT, img);
      fft(Kspec[j], -1);
    });
    for (const t of P.concat(D)) {   // pairs left unused stay 0
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo.get(t));
      gl.clearColor(0, 0, 0, 0); gl.clear(gl.COLOR_BUFFER_BIT);
    }
  }

  // kernel k's source channel, one-hot (params.src[g]: column-major mat3x4, column c, row i)
  function srcOf(params, k) {
    const S = params.src[k >> 2], i = k & 3;
    return [S[i], S[4 + i], S[8 + i]];
  }

  // The creatures' step, into stateFbo[1 - cur] (glsim.js flips cur after)
  function step(tool) {
    const params = I.params, bul = I.bul, bActive = I.bulletsActive;
    const bmode = bActive ? bul.mode : 0, bmask = bmode === 2;
    const state = I.state[I.cur];
    gl.viewport(0, 0, W, H);

    let p = prog.fftpack;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, packFbo);
    I.bindTex(0, state, p, 'uState');
    I.bindTex(1, I.bulTex, p, 'uBullets');
    gl.uniform1i(I.loc(p, 'uBMask'), bmask ? 1 : 0);
    gl.uniform3fv(I.loc(p, 'uBulletCh'), bul.ch);
    I.draw4();
    for (let i = 0; i < (bmask ? 4 : 2); i++) fft(F[i], -1);

    p = prog.fftmul;
    for (let j = 0; j < pairsUsed; j++) {
      gl.useProgram(p);
      gl.viewport(0, 0, W, H);
      gl.bindFramebuffer(gl.FRAMEBUFFER, mulFbo[j]);
      F.forEach((t, i) => I.bindTex(i, t, p, 'uF' + i));
      I.bindTex(4, Kspec[j], p, 'uK');
      gl.uniform3fv(I.loc(p, 'uSrc0'), srcOf(params, 2 * j));
      gl.uniform3fv(I.loc(p, 'uSrc1'), srcOf(params, 2 * j + 1));
      gl.uniform1i(I.loc(p, 'uBMask'), bmask ? 1 : 0);
      gl.uniform3fv(I.loc(p, 'uBulletCh'), bul.ch);
      gl.uniform1f(I.loc(p, 'uScale'), 1 / (W * H));
      I.draw4();
    }
    for (let j = 0; j < pairsUsed; j++) { fft(P[j], 1); if (bmask) fft(D[j], 1); }

    if (bmask) {   // the bullets' hidden shares, by kernel group, into F (free now)
      p = prog.fftgather;
      gl.useProgram(p);
      gl.viewport(0, 0, W, H);
      gl.bindFramebuffer(gl.FRAMEBUFFER, packFbo);
      D.forEach((t, j) => I.bindTex(j, t, p, 'uP' + j));
      I.draw4();
    }

    p = prog.fftgrow;
    gl.useProgram(p);
    gl.viewport(0, 0, W, H);
    gl.bindFramebuffer(gl.FRAMEBUFFER, I.stateFbo[1 - I.cur]);
    P.forEach((t, j) => I.bindTex(j, t, p, 'uP' + j));
    F.forEach((t, g) => I.bindTex(8 + g, t, p, 'uD' + g));
    I.bindTex(12, state, p, 'uState');
    I.bindTex(13, I.bulTex, p, 'uBullets');
    I.bindTex(14, kimgTex, p, 'uKimg');
    gl.uniform1i(I.loc(p, 'uKR'), KR);
    gl.uniform2i(I.loc(p, 'uSize'), W, H);
    for (let g = 0; g < 4; g++) {
      gl.uniformMatrix3x4fv(I.loc(p, `uSrc[${g}]`), false, params.src[g]);
      gl.uniformMatrix4x3fv(I.loc(p, `uDst[${g}]`), false, params.dst[g]);
      gl.uniform4fv(I.loc(p, `uM[${g}]`), params.m[g]);
      gl.uniform4fv(I.loc(p, `uS[${g}]`), params.s[g]);
      gl.uniform4fv(I.loc(p, `uH[${g}]`), params.h[g]);
    }
    gl.uniform3fv(I.loc(p, 'uGain'), I.gain);
    gl.uniform1f(I.loc(p, 'uDt'), I.dt);
    gl.uniform1i(I.loc(p, 'uGrowth'), I.growthFn);
    gl.uniform1i(I.loc(p, 'uMaskOn'), tool ? tool.action + 1 : 0);
    gl.uniform2f(I.loc(p, 'uMaskC'), tool ? tool.col : 0, tool ? tool.row : 0);
    gl.uniform1f(I.loc(p, 'uMaskR'), tool ? tool.r : 0);
    gl.uniform3f(I.loc(p, 'uMaskCh'), tool && tool.ch === 0 ? 1 : 0, tool && tool.ch === 1 ? 1 : 0, tool && tool.ch === 2 ? 1 : 0);
    gl.uniform1f(I.loc(p, 'uReach'), I.reach);
    gl.uniform1i(I.loc(p, 'uBulletMode'), bmode);
    gl.uniform3fv(I.loc(p, 'uBulletCh'), bul.ch);
    gl.uniform1f(I.loc(p, 'uBulletHit'), bul.hit);
    I.draw4();
  }

  // Plug in once the engine is up
  const init0 = SimGL.init;
  SimGL.init = async function (canvas, opts) { await init0(canvas, opts); await attach(); };

  return { on: true, get ready() { return ready; } };
})();
