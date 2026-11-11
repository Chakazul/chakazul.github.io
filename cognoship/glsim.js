"use strict";
// ============================================================================================
//  WebGL2 engine for Cognoship: 3-channel Lenia on a torus, with cognitive masking.
//
//  Per step: one sim pass (convolution + growth + mask/erase). The mask's absolute centre comes
//  from the CPU each step (core.js extrapolates the focused creature's CoM within a frame).
//
//  Per rendered frame, two small readbacks:
//    readBlocks()      coarse 4x4-block mass map, so core.js can find every creature
//                      (connected components) and follow the focused one;
//    readChannelComs() per-channel CoMs inside a window around the focused creature (reduce ->
//                      com), for the body-axis view mode.
//
//  Dots channel (setDots, ported from Carlman): a separate single-channel Lenia field of small
//  stationary solitons under its own rule, stepped after the creatures and eaten wherever their
//  mass is (dotsim.glsl); readDotSites() gives the mass left around each dot, for scoring.
//  Bullets channel (setBullets): another single-channel field, of gliders (bulletsim.glsl), which
//  takes mass out of the creatures (or masks them, sim.glsl) while they wear it down;
//  readBulletBlocks() gives its coarse mass map, so the page can count and find the bullets.
//
//  Board state is a ping-pong pair of RGBA32F textures (rgb = channels). Needs
//  EXT_color_buffer_float, as Carlman does.
// ============================================================================================

const SimGL = (function () {
  const SHADER_DIR = './shaders/';
  const NAMES = ['vertex', 'sim', 'reduce', 'com', 'blocks', 'draw', 'dotsim', 'dotsites', 'bulletsim', 'blocksum'];
  const GRID = 16;    // reduce.glsl output is GRID x GRID
  const BLOCK = 4;    // blocks.glsl block edge

  let gl = null;
  let prog = {}, uCache = new Map();
  let W = 0, H = 0, dt = 0.5, gain = [1, 1, 1], growthFn = 1, reach = 0;
  let state = [], stateFbo = [], cur = 0;
  let tapTex = null, tapCount = 0, params = null, kernelCount = 1;
  let redTex = [], redFbo = null, redBlock = 1;
  let comTex = null, comFbo = null;
  let blkTex = null, blkChanTex = null, blkFbo = null, BW = 0, BH = 0, blkBuf = null, blkChanBuf = null;
  const comBuf = new Float32Array(16);
  // dots channel: ping-pong pair (R used), its rule's tap table, and the per-dot site pass
  const MAX_SITES = 256;
  const DOT_EAT = 0.1;   // creature mass above this erases the dot there (Carlman's EAT_THRESHOLD)
  let dotTex = [], dotFbo = [], dcur = 0, dotsOn = false, dotRGB = [0.8, 0.86, 1];
  let dotOff = [0, 0];   // the dots' drift, cells in [0, W) x [0, H): dot cell (x, y) is at board cell (x, y) - dotOff
  let dotTapTex = null, dotTapCount = 0, dotRule = null;
  let siteTex = null, siteOutTex = null, siteOutFbo = null;
  const siteIn = new Float32Array(MAX_SITES * 4), siteOut = new Float32Array(MAX_SITES * 4);
  // bullets channel: ping-pong pair (R used), its rule's tap table, and how it acts (setBullets)
  const BUL_BLOCK = 8;   // readBulletBlocks block edge
  let bulTex = [], bulFbo = [], bcur = 0, bulletsOn = false, bulRGB = [1, 0.42, 0.24];
  let bulTapTex = null, bulTapCount = 0, bulRule = null;
  // idle: the field is empty (cleared, or nothing fired since): its steps and its effect on the
  // creatures' step (sim.glsl's bullet path) are skipped until putBullets puts some in.
  let bulletsIdle = true;
  let bul = { mode: 1, ch: [1, 1, 1], wearCh: [1, 1, 1], hit: 0.2, wear: 0.05, sub: 2 };
  let bblkTex = null, bblkFbo = null, BBW = 0, BBH = 0, bblkBuf = null;
  // Plug-ins (fftsim.js, FFT mode, on unless ?fft=0): ext.step(tool) replaces the creatures' sim pass (it must write the
  // next board into stateFbo[1 - cur]; step() then flips cur, steps the dots and bullets);
  // ext.rule(taps) follows setRule. Both null: the direct convolution (sim.glsl) as ever.
  const ext = { step: null, rule: null };
  let lastTaps = null;   // setRule's tap table (offs, weights), for ext.rule

  // ---- setup ----------------------------------------------------------------------------
  async function loadSources() {
    const out = {};
    await Promise.all(NAMES.map(async n => {
      const r = await fetch(SHADER_DIR + n + '.glsl');
      if (!r.ok) throw new Error('shader ' + n + ': HTTP ' + r.status);
      out[n] = await r.text();
    }));
    return out;
  }

  function compile(type, src, name) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
      throw new Error(name + ': ' + gl.getShaderInfoLog(s));
    return s;
  }

  function link(vs, fsSrc, name) {
    const p = gl.createProgram();
    gl.attachShader(p, vs);
    gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fsSrc, name));
    gl.bindAttribLocation(p, 0, 'a_position');
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS))
      throw new Error(name + ': ' + gl.getProgramInfoLog(p));
    return p;
  }

  function loc(p, name) {
    const key = p.__id + ':' + name;
    if (!uCache.has(key)) uCache.set(key, gl.getUniformLocation(p, name));
    return uCache.get(key);
  }

  function makeTex(w, h, data = null) {
    const t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, w, h, 0, gl.RGBA, gl.FLOAT, data);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return t;
  }

  function makeFbo(texs) {
    const f = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, f);
    texs.forEach((t, i) => gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + i, gl.TEXTURE_2D, t, 0));
    if (texs.length > 1) gl.drawBuffers(texs.map((_, i) => gl.COLOR_ATTACHMENT0 + i));
    const st = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (st !== gl.FRAMEBUFFER_COMPLETE) throw new Error('framebuffer incomplete: 0x' + st.toString(16));
    return f;
  }

  // opts: { W, H, R, kernels }
  async function init(canvas, opts) {
    gl = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: false });
    if (!gl) throw new Error('WebGL2 is not available in this browser.');
    if (!gl.getExtension('EXT_color_buffer_float'))
      throw new Error('This GPU cannot render to float textures (EXT_color_buffer_float).');
    if (gl.getParameter(gl.MAX_DRAW_BUFFERS) < 4)
      throw new Error('This GPU supports fewer than 4 draw buffers.');

    const src = await loadSources();
    VERTEX_SRC = src.vertex;
    const vs = compile(gl.VERTEX_SHADER, src.vertex, 'vertex');
    let id = 0;
    for (const n of NAMES) if (n !== 'vertex') { prog[n] = link(vs, src[n], n); prog[n].__id = id++; }

    const quad = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, quad);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);

    W = opts.W; H = opts.H;
    if (Math.max(W, H) > GRID * 64) throw new Error('board too large for the reduction');
    state = [makeTex(W, H), makeTex(W, H)];
    stateFbo = state.map(t => makeFbo([t]));
    redBlock = Math.ceil(Math.max(W, H) / GRID);
    redTex = [0, 1, 2, 3].map(() => makeTex(GRID, GRID));
    redFbo = makeFbo(redTex);
    comTex = makeTex(4, 1);
    comFbo = makeFbo([comTex]);
    BW = Math.ceil(W / BLOCK); BH = Math.ceil(H / BLOCK);
    blkTex = makeTex(BW, BH);
    blkChanTex = makeTex(BW, BH);
    blkFbo = makeFbo([blkTex, blkChanTex]);
    blkBuf = new Float32Array(BW * BH * 4);
    blkChanBuf = new Float32Array(BW * BH * 4);
    dotTex = [makeTex(W, H), makeTex(W, H)];
    dotFbo = dotTex.map(t => makeFbo([t]));
    bulTex = [makeTex(W, H), makeTex(W, H)];
    bulFbo = bulTex.map(t => makeFbo([t]));
    BBW = Math.ceil(W / BUL_BLOCK); BBH = Math.ceil(H / BUL_BLOCK);
    bblkTex = makeTex(BBW, BBH);
    bblkFbo = makeFbo([bblkTex]);
    bblkBuf = new Float32Array(BBW * BBH * 4);
    siteTex = makeTex(MAX_SITES, 1);
    siteOutTex = makeTex(MAX_SITES, 1);
    siteOutFbo = makeFbo([siteOutTex]);
    setRule(opts.kernels, opts.R, opts.kn ?? 1);
    clear();
  }

  // kn: kernel core (LeniaF numbering, see rule.js shell).
  function setRule(kernels, R, kn = 1) {
    const taps = Rule.buildTaps(kernels, R, kn);
    tapCount = taps.count;
    if (tapCount > gl.getParameter(gl.MAX_TEXTURE_SIZE)) throw new Error('kernel too large');
    const data = new Float32Array(tapCount * 5 * 4);
    let rmax = 0;
    taps.offs.forEach(([dx, dy], t) => {
      data[t * 4] = dx; data[t * 4 + 1] = dy;
      rmax = Math.max(rmax, Math.hypot(dx, dy));
      for (let g = 0; g < 4; g++)
        for (let j = 0; j < 4; j++)
          data[((g + 1) * tapCount + t) * 4 + j] = taps.weights[t * Rule.SLOTS + g * 4 + j];
    });
    if (tapTex) gl.deleteTexture(tapTex);
    tapTex = makeTex(tapCount, 5, data);
    reach = rmax;
    params = Rule.buildParams(kernels);
    kernelCount = kernels.length;
    lastTaps = taps;
    if (ext.rule) ext.rule(taps);
  }

  // The update: A += dt * gain * sum h G(U), G = growth function gn (1 polynomial, 2 Gaussian).
  // normalised = divide by sum h per destination channel (legacy LeniaNDKC style); LeniaF's
  // kernels arrive with h already normalised (rule.js leniaF), so they use normalised = false.
  function setUpdate(normalised, newDt, gn = 1) {
    dt = newDt; growthFn = gn;
    gain = normalised ? params.hsum.map(h => h > 0 ? 1 / h : 0) : [1, 1, 1];
  }

  // ---- board ----------------------------------------------------------------------------
  function clear() {
    bulletsIdle = true;
    for (const f of stateFbo.concat(dotFbo, bulFbo)) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, f);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
  }

  // Place a pattern (formCells output) with its top-left at (row0, col0), on a cleared board.
  function seed(form, row0, col0) {
    const buf = new Float32Array(W * H * 4);
    for (let c = 0; c < 3; c++)
      for (let y = 0; y < form.h; y++)
        for (let x = 0; x < form.w; x++) {
          const r = ((row0 + y) % H + H) % H, q = ((col0 + x) % W + W) % W;
          buf[(r * W + q) * 4 + c] = form.cells[c][y][x];
        }
    cur = 0;
    for (const t of state) {
      gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
    }
  }

  // Add a form onto the current board (cell-wise max with what is there), top-left at
  // (row0, col0), wrapping round the torus. Reads the board back once, so only for clicks.
  function add(form, row0, col0) {
    const buf = new Float32Array(W * H * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, stateFbo[cur]);
    gl.readPixels(0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
    for (let c = 0; c < 3; c++)
      for (let y = 0; y < form.h; y++)
        for (let x = 0; x < form.w; x++) {
          const r = ((row0 + y) % H + H) % H, q = ((col0 + x) % W + W) % W, i = (r * W + q) * 4 + c;
          buf[i] = Math.max(buf[i], form.cells[c][y][x]);
        }
    gl.bindTexture(gl.TEXTURE_2D, state[cur]);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
  }

  // ---- dots channel ---------------------------------------------------------------------
  // rule: { R, mu, sigma, dt }: one kernel of the standard Lenia shell, quad4 core, beta [1],
  // normalised, taps below 1e-7 dropped (as Carlman's kernelData), and a Gaussian growth.
  // The classic Lenia kernel: rings b (betas) over r = distance / R, each a quad4 bump
  // (Lenia/Python Board kernel, kn = 1), normalised, taps below 1e-7 dropped.
  function shellTaps(R, betas = [1]) {
    const quad4 = r => { const q = 4 * r * (1 - r); return q * q * q * q; };
    const B = betas.length;
    const KS = 2 * R + 1, taps = [];
    let sum = 0;
    for (let iy = 0; iy < KS; iy++)
      for (let ix = 0; ix < KS; ix++) {
        const yy = -1 + 2 * iy / (KS - 1), xx = -1 + 2 * ix / (KS - 1);
        const r = Math.hypot(xx, yy), Br = B * r;
        const v = r > 1 || Br >= B ? 0 : betas[Math.min(B - 1, Math.floor(Br))] * quad4(Br - Math.floor(Br));
        if (v > 0) { taps.push([ix - R, iy - R, v]); sum += v; }
      }
    const kept = taps.filter(t => t[2] / sum > 1e-7);
    const data = new Float32Array(kept.length * 4);
    kept.forEach(([dx, dy, v], i) => data.set([dx, dy, v / sum, 0], i * 4));
    return { tex: makeTex(kept.length, 1, data), count: kept.length };
  }
  function setDotRule(rule) {
    if (dotTapTex) gl.deleteTexture(dotTapTex);
    ({ tex: dotTapTex, count: dotTapCount } = shellTaps(rule.R));
    dotRule = rule;
  }
  function clearFbos(fbos) {
    for (const f of fbos) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, f);
      gl.clearColor(0, 0, 0, 0);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
  }
  // Write an h x w block (row-major `flat`, or zeros if null) into the R of tex, top-left at
  // (row0, col0), replacing what is there. The block must lie inside the board (no wrap).
  function putField(tex, flat, h, w, row0, col0) {
    const buf = new Float32Array(w * h * 4);
    if (flat) for (let i = 0; i < w * h; i++) buf[i * 4] = flat[i];
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, col0, row0, w, h, gl.RGBA, gl.FLOAT, buf);
  }
  // Dots on or off. Either way the field starts empty.
  function setDots(on) { dotsOn = on; clearFbos(dotFbo); }
  // Write an h x w block of dot values (row-major `flat`, or zeros if null), top-left at (row0, col0).
  function putDots(flat, h, w, row0, col0) { putField(dotTex[dcur], flat, h, w, row0, col0); }

  // ---- bullets channel ------------------------------------------------------------------
  // rule: { R, betas, mu, sigma, dt }: the classic Lenia rule, one kernel of quad4 rings (betas,
  // default [1]), polynomial growth.
  function setBulletRule(rule) {
    if (bulTapTex) gl.deleteTexture(bulTapTex);
    ({ tex: bulTapTex, count: bulTapCount } = shellTaps(rule.R, rule.betas));
    bulRule = rule;
  }
  // Bullets on (opts: { mode: 1 erase | 2 mask, ch: [r, g, b] 0/1 (acted on), wearCh (wearing them),
  // hit, wear, sub: bullet steps per
  // step }) or off (null). Turning them on or off empties the field; a change of opts while on keeps it.
  function setBullets(opts) {
    const was = bulletsOn;
    bulletsOn = !!opts;
    if (opts) bul = { ...bul, ...opts };
    if (was !== bulletsOn) { clearFbos(bulFbo); bulletsIdle = true; }
  }
  function putBullets(flat, h, w, row0, col0) { putField(bulTex[bcur], flat, h, w, row0, col0); bulletsIdle = false; }
  function clearBullets() { clearFbos(bulFbo); bulletsIdle = true; }
  function stepBullets() {
    const p = prog.bulletsim;
    gl.useProgram(p);
    gl.viewport(0, 0, W, H);
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uTapCount'), bulTapCount);
    gl.uniform1f(loc(p, 'uMu'), bulRule.mu);
    gl.uniform1f(loc(p, 'uSigma'), bulRule.sigma);
    gl.uniform1f(loc(p, 'uDt'), bulRule.dt);
    gl.uniform1f(loc(p, 'uWear'), bul.wear);
    gl.uniform3fv(loc(p, 'uWearCh'), bul.wearCh);
    for (let i = 0; i < bul.sub; i++) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, bulFbo[1 - bcur]);
      bindTex(0, bulTex[bcur], p, 'uField');
      bindTex(1, state[cur], p, 'uState');
      bindTex(2, bulTapTex, p, 'uTaps');
      draw4();
      bcur = 1 - bcur;
    }
  }
  // The bullets' coarse mass map: BBW x BBH sums over BUL_BLOCK x BUL_BLOCK cells, row-major
  // (row 0 = board row 0), one float per block.
  function readBulletBlocks() {
    const p = prog.blocksum;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, bblkFbo);
    gl.viewport(0, 0, BBW, BBH);
    bindTex(0, bulTex[bcur], p, 'uField');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uBlock'), BUL_BLOCK);
    draw4();
    gl.readPixels(0, 0, BBW, BBH, gl.RGBA, gl.FLOAT, bblkBuf);
    const out = new Float32Array(BBW * BBH);
    for (let i = 0; i < out.length; i++) out[i] = bblkBuf[i * 4];
    return out;
  }
  function setBulletColor(rgb) { bulRGB = rgb; }
  // The dots step, on the board the creatures' step just wrote (state[cur]).
  function stepDots() {
    const p = prog.dotsim;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, dotFbo[1 - dcur]);
    gl.viewport(0, 0, W, H);
    bindTex(0, dotTex[dcur], p, 'uDots');
    bindTex(1, state[cur], p, 'uState');
    bindTex(2, dotTapTex, p, 'uTaps');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uTapCount'), dotTapCount);
    gl.uniform1f(loc(p, 'uMu'), dotRule.mu);
    gl.uniform1f(loc(p, 'uSigma'), dotRule.sigma);
    gl.uniform1f(loc(p, 'uDt'), dotRule.dt);
    gl.uniform1f(loc(p, 'uEat'), DOT_EAT);
    gl.uniform2i(loc(p, 'uOff'), Math.round(dotOff[0]) % W, Math.round(dotOff[1]) % H);
    draw4();
    dcur = 1 - dcur;
  }
  // The dot mass within `radius` cells (a square window) of each site [col, row], as an array
  // (at most MAX_SITES sites).
  function readDotSites(sites, radius) {
    const n = Math.min(sites.length, MAX_SITES);
    if (!n) return [];
    siteIn.fill(0);
    for (let i = 0; i < n; i++) { siteIn[i * 4] = sites[i][0]; siteIn[i * 4 + 1] = sites[i][1]; }
    gl.bindTexture(gl.TEXTURE_2D, siteTex);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, MAX_SITES, 1, gl.RGBA, gl.FLOAT, siteIn);
    const p = prog.dotsites;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, siteOutFbo);
    gl.viewport(0, 0, n, 1);
    bindTex(0, dotTex[dcur], p, 'uDots');
    bindTex(1, siteTex, p, 'uSites');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uRadius'), radius);
    draw4();
    gl.readPixels(0, 0, n, 1, gl.RGBA, gl.FLOAT, siteOut);
    return Array.from({ length: n }, (_, i) => siteOut[i * 4]);
  }
  function setDotColor(rgb) { dotRGB = rgb; }
  let bgRGB = [0, 0, 0];
  function setBackground(rgb) { bgRGB = rgb; }   // behind the board, where it is empty (draw.glsl uBg)
  function setDotOffset(x, y = 0) { dotOff = [((x % W) + W) % W, ((y % H) + H) % H]; }

  // Erase (action 1: set to 0) or draw (2: set to 1) channel ch inside a circle (centre col, row,
  // radius r, cells; round the torus) on the board as it is, with no step: the tools while paused.
  // Reads the board back, so only while a tool is held.
  function paint(tool) {
    const buf = new Float32Array(W * H * 4);
    gl.bindFramebuffer(gl.FRAMEBUFFER, stateFbo[cur]);
    gl.readPixels(0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
    const v = tool.action === 2 ? 1 : 0, R = Math.ceil(tool.r), r2 = tool.r * tool.r;
    const c0 = Math.round(tool.col), r0 = Math.round(tool.row);
    for (let dy = -R; dy <= R; dy++)
      for (let dx = -R; dx <= R; dx++) {
        const ddx = c0 + dx - tool.col, ddy = r0 + dy - tool.row;
        if (ddx * ddx + ddy * ddy > r2) continue;
        const y = ((r0 + dy) % H + H) % H, x = ((c0 + dx) % W + W) % W;
        buf[(y * W + x) * 4 + tool.ch] = v;
      }
    gl.bindTexture(gl.TEXTURE_2D, state[cur]);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, W, H, gl.RGBA, gl.FLOAT, buf);
  }

  function bindTex(unit, tex, p, name) {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.uniform1i(loc(p, name), unit);
  }

  function draw4() { gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4); }

  // tool: null, or { col, row, r, ch, action }: absolute circle centre in board cells;
  // action 0 = occlude (mask), 1 = erase (set to 0), 2 = draw (set to 1) inside the circle.
  function step(tool) {
    if (ext.step) {   // a plug-in's sim pass (fftsim.js)
      ext.step(tool);
      cur = 1 - cur;
      if (dotsOn && dotTapTex) stepDots();
      if (bulletsOn && bulTapTex && !bulletsIdle) stepBullets();
      return;
    }
    const p = prog.sim;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, stateFbo[1 - cur]);
    gl.viewport(0, 0, W, H);
    bindTex(0, state[cur], p, 'uState');
    bindTex(1, tapTex, p, 'uTaps');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uTapCount'), tapCount);
    for (let g = 0; g < 4; g++) {
      gl.uniformMatrix3x4fv(loc(p, `uSrc[${g}]`), false, params.src[g]);
      gl.uniformMatrix4x3fv(loc(p, `uDst[${g}]`), false, params.dst[g]);
      gl.uniform4fv(loc(p, `uM[${g}]`), params.m[g]);
      gl.uniform4fv(loc(p, `uS[${g}]`), params.s[g]);
      gl.uniform4fv(loc(p, `uH[${g}]`), params.h[g]);
    }
    gl.uniform3fv(loc(p, 'uGain'), gain);
    gl.uniform1f(loc(p, 'uDt'), dt);
    gl.uniform1i(loc(p, 'uGrowth'), growthFn);
    gl.uniform1i(loc(p, 'uMaskOn'), tool ? tool.action + 1 : 0);
    gl.uniform2f(loc(p, 'uMaskC'), tool ? tool.col : 0, tool ? tool.row : 0);
    gl.uniform1f(loc(p, 'uMaskR'), tool ? tool.r : 0);
    gl.uniform3f(loc(p, 'uMaskCh'), tool && tool.ch === 0 ? 1 : 0, tool && tool.ch === 1 ? 1 : 0, tool && tool.ch === 2 ? 1 : 0);
    gl.uniform1f(loc(p, 'uReach'), reach);
    bindTex(2, bulTex[bcur], p, 'uBullets');
    gl.uniform1i(loc(p, 'uBulletMode'), bulletsOn && bulTapTex && !bulletsIdle ? bul.mode : 0);
    gl.uniform3fv(loc(p, 'uBulletCh'), bul.ch);
    gl.uniform1f(loc(p, 'uBulletHit'), bul.hit);
    draw4();
    cur = 1 - cur;
    if (dotsOn && dotTapTex) stepDots();
    if (bulletsOn && bulTapTex && !bulletsIdle) stepBullets();
  }

  // Coarse map: BW x BH texels of (mass, sum m*dx, sum m*dy, potential) per BLOCK x BLOCK
  // patch (blocks.glsl). The same call fills blockChannels: per block (mass R, G, B, 0).
  function readBlocks() {
    const p = prog.blocks;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, blkFbo);
    gl.viewport(0, 0, BW, BH);
    bindTex(0, state[cur], p, 'uState');
    bindTex(1, tapTex, p, 'uTaps');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uTapCount'), tapCount);
    for (let g = 0; g < 4; g++) gl.uniformMatrix3x4fv(loc(p, `uSrc[${g}]`), false, params.src[g]);
    gl.uniform1f(loc(p, 'uInvK'), 1 / kernelCount);
    draw4();
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    gl.readPixels(0, 0, BW, BH, gl.RGBA, gl.FLOAT, blkBuf);
    gl.readBuffer(gl.COLOR_ATTACHMENT1);
    gl.readPixels(0, 0, BW, BH, gl.RGBA, gl.FLOAT, blkChanBuf);
    gl.readBuffer(gl.COLOR_ATTACHMENT0);
    return blkBuf;
  }

  // Per-channel CoMs within `radius` cells of (col, row): 4 x (row, col, mass, valid) for
  // channels R, G, B, then all three together.
  function readChannelComs(col, row, radius) {
    let p = prog.reduce;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, redFbo);
    gl.viewport(0, 0, GRID, GRID);
    bindTex(0, state[cur], p, 'uState');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform1i(loc(p, 'uBlock'), redBlock);
    gl.uniform2f(loc(p, 'uCenter'), col, row);
    gl.uniform1f(loc(p, 'uRadius'), radius);
    draw4();

    p = prog.com;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, comFbo);
    gl.viewport(0, 0, 4, 1);
    for (let i = 0; i < 4; i++) bindTex(i, redTex[i], p, 'uRed' + i);
    gl.uniform2i(loc(p, 'uSize'), W, H);
    draw4();
    gl.readPixels(0, 0, 4, 1, gl.RGBA, gl.FLOAT, comBuf);
    return comBuf;
  }

  // view: { canvasW, canvasH, worldScale, zoom: [x0, y0, edge, pxPerCell], rot: mat2 (column-
  // major, view (right, down) -> board (col, row)), focus: [col, row, on],
  //   stars: star field on, time: seconds (star sway), starOff: [x, y] drift (cells, at unit
  //   layer speed), starDir: [x, y] unit drift direction, starSpeed, starBright: 1 = default }
  function render(view) {
    const p = prog.draw;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, view.canvasW, view.canvasH);
    bindTex(0, state[cur], p, 'uState');
    bindTex(1, dotTex[dcur], p, 'uDots');
    bindTex(2, bulTex[bcur], p, 'uBullets');
    gl.uniform1f(loc(p, 'uBulletsOn'), bulletsOn ? 1 : 0);
    gl.uniform3fv(loc(p, 'uBulletRGB'), bulRGB);
    gl.uniform1f(loc(p, 'uDotsOn'), dotsOn ? 1 : 0);
    gl.uniform3fv(loc(p, 'uDotRGB'), dotRGB);
    gl.uniform2fv(loc(p, 'uDotOff'), dotOff);
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform2f(loc(p, 'uCanvas'), view.canvasW, view.canvasH);
    gl.uniform1f(loc(p, 'uWorldScale'), view.worldScale);
    gl.uniform4fv(loc(p, 'uZoom'), view.zoom);
    gl.uniformMatrix2fv(loc(p, 'uView'), false, view.rot);
    gl.uniform3fv(loc(p, 'uFocus'), view.focus);
    gl.uniform1f(loc(p, 'uStars'), view.stars ? 1 : 0);
    gl.uniform1f(loc(p, 'uTime'), view.time || 0);
    gl.uniform2fv(loc(p, 'uStarOff'), view.starOff || [-(view.time || 0), 0]);
    gl.uniform2fv(loc(p, 'uStarDir'), view.starDir || [-1, 0]);
    gl.uniform1f(loc(p, 'uStarSpeed'), view.starSpeed ?? 1);
    gl.uniform1f(loc(p, 'uStarBright'), view.starBright ?? 1);
    gl.uniform3fv(loc(p, 'uBg'), bgRGB);
    draw4();
  }

  function finish() { gl.finish(); }

  // For plug-ins (fftsim.js): the engine's context, helpers and live state. Read the getters each
  // time: the values change (cur, the rule, the bullets).
  const internals = {
    ext, makeTex, makeFbo, loc, bindTex, draw4,
    // a fragment shader linked with the engine's vertex shader
    link: fsSrc => { const vs = compile(gl.VERTEX_SHADER, VERTEX_SRC, 'vertex'), p = link(vs, fsSrc, 'plugin'); p.__id = 'x' + (++pluginIds); return p; },
    get gl() { return gl; }, get W() { return W; }, get H() { return H; },
    get state() { return state; }, get stateFbo() { return stateFbo; }, get cur() { return cur; },
    get params() { return params; }, get dt() { return dt; }, get gain() { return gain; }, get growthFn() { return growthFn; },
    get reach() { return reach; }, get taps() { return lastTaps; },
    get bulletsActive() { return bulletsOn && !!bulTapTex && !bulletsIdle; }, get bul() { return bul; }, get bulTex() { return bulTex[bcur]; },
  };
  let VERTEX_SRC = '', pluginIds = 0;

  return { init, setRule, setUpdate, seed, add, paint, step, readBlocks, readChannelComs, render, finish, BLOCK,
           setDotRule, setDots, putDots, readDotSites, setDotColor, setDotOffset, setBackground, MAX_SITES,
           setBulletRule, setBullets, putBullets, clearBullets, readBulletBlocks, setBulletColor, BUL_BLOCK,
           get bulletsIdle() { return bulletsIdle; }, get BBW() { return BBW; }, get BBH() { return BBH; },
           get W() { return W; }, get H() { return H; }, get BW() { return BW; }, get BH() { return BH; },
           get blockChannels() { return blkChanBuf; },
           get tapCount() { return tapCount; }, internals };
})();
