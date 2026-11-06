"use strict";
// ============================================================================================
//  WebGL2 engine for Cognoship: 3-channel Lenia on a torus, with cognitive masking.
//
//  Per step: one sim pass (convolution + growth + mask/erase). The mask's absolute centre comes
//  from the CPU each step (app.js extrapolates the focused creature's CoM within a frame).
//
//  Per rendered frame, two small readbacks:
//    readBlocks()      coarse 4x4-block mass map, so app.js can find every creature
//                      (connected components) and follow the focused one;
//    readChannelComs() per-channel CoMs inside a window around the focused creature (reduce ->
//                      com), for the body-axis view mode.
//
//  Board state is a ping-pong pair of RGBA32F textures (rgb = channels). Needs
//  EXT_color_buffer_float, as Carlman does.
// ============================================================================================

const SimGL = (function () {
  const SHADER_DIR = './shaders/';
  const NAMES = ['vertex', 'sim', 'reduce', 'com', 'blocks', 'draw'];
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
    for (const f of stateFbo) {
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

  function bindTex(unit, tex, p, name) {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.uniform1i(loc(p, name), unit);
  }

  function draw4() { gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4); }

  // tool: null, or { col, row, r, ch, action }: absolute circle centre in board cells;
  // action 0 = occlude (mask), 1 = erase (set to 0), 2 = draw (set to 1) inside the circle.
  function step(tool) {
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
    draw4();
    cur = 1 - cur;
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
  //   stars: star field on, time: seconds (star drift) }
  function render(view) {
    const p = prog.draw;
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    gl.viewport(0, 0, view.canvasW, view.canvasH);
    bindTex(0, state[cur], p, 'uState');
    gl.uniform2i(loc(p, 'uSize'), W, H);
    gl.uniform2f(loc(p, 'uCanvas'), view.canvasW, view.canvasH);
    gl.uniform1f(loc(p, 'uWorldScale'), view.worldScale);
    gl.uniform4fv(loc(p, 'uZoom'), view.zoom);
    gl.uniformMatrix2fv(loc(p, 'uView'), false, view.rot);
    gl.uniform3fv(loc(p, 'uFocus'), view.focus);
    gl.uniform1f(loc(p, 'uStars'), view.stars ? 1 : 0);
    gl.uniform1f(loc(p, 'uTime'), view.time || 0);
    draw4();
  }

  function finish() { gl.finish(); }

  return { init, setRule, setUpdate, seed, add, step, readBlocks, readChannelComs, render, finish, BLOCK,
           get W() { return W; }, get H() { return H; }, get BW() { return BW; }, get BH() { return BH; },
           get blockChannels() { return blkChanBuf; },
           get tapCount() { return tapCount; } };
})();
