"use strict";
// ============================================================================================
//  GPU Lenia engine for the CARL maze demo.
//
//  The CPU demo (CARL/maze_playground.html) ran the whole simulation in JS: an FFT-accelerated
//  convolution per step, a full-board sweep for the center of mass, and a per-pixel loop to
//  paint the canvas. All three are per-cell-independent work, so all three live here as
//  fragment shader passes instead. What stays on the CPU is the part that genuinely cannot
//  move: the policy network, which runs in onnxruntime-web's WASM backend.
//
//  The design constraint that shapes this file is that the CPU needs exactly two things back
//  from the GPU each step -- where the soliton is, and the 96x96x4 window the policy reads --
//  and every readback is a pipeline stall. So a step queues five passes with no synchronization
//  between them (action, sim, reduce, com, crop) and the caller then picks up both results in
//  one readback() call. The crop pass reads the CoM out of a texture rather than a uniform
//  precisely so that it does not have to wait for the CPU to be told where to look.
//
//  Board state lives in a ring of four R32F textures, not the usual ping-pong pair, because the
//  policy is fed a 4-frame stack and all four frames have to be croppable at a shared origin.
//  Stepping writes over the four-steps-ago frame, which is exactly the one falling out of the
//  stack, so four slots is the whole requirement.
// ============================================================================================

const SimGL = (function () {
  const SHADER_DIR = './shaders/';
  const NAMES = ['vertex', 'sim', 'action', 'reduce', 'com', 'crop', 'draw'];
  const RING = 4;          // frames kept for the policy's frame stack
  const GRID = 16;         // stage-1 reduction output is GRID x GRID (see reduce.glsl)

  let gl = null, canvas = null;
  let prog = {}, uCache = new Map(), progSeq = 0;
  let W = 0, H = 0, net = 96, kR = 18;
  let ruleMu = 0.24, ruleSigma = 0.024, ruleDt = 0.1;

  let ringTex = [], ringFbo = [], ringIdx = 0;
  let scratchTex = null, scratchFbo = null;
  let wallTex = null, kernelTex = null;
  let redTex0 = null, redTex1 = null, redFbo = null, redBlock = 1;
  let comTex = null, comFbo = null;
  let cropTex = null, cropFbo = null;
  let vao = null;

  const comPix = new Float32Array(4);
  let cropPix = null;

  // ------------------------------------------------------------------------------------------
  //  Boilerplate
  // ------------------------------------------------------------------------------------------
  function compile(srcText, type, label) {
    const sh = gl.createShader(type);
    gl.shaderSource(sh, srcText);
    gl.compileShader(sh);
    if (!gl.getShaderParameter(sh, gl.COMPILE_STATUS))
      throw new Error(label + ': ' + gl.getShaderInfoLog(sh));
    return sh;
  }

  function link(vsSrc, fsSrc, label) {
    const p = gl.createProgram();
    gl.attachShader(p, compile(vsSrc, gl.VERTEX_SHADER, label + ' (vertex)'));
    gl.attachShader(p, compile(fsSrc, gl.FRAGMENT_SHADER, label + ' (fragment)'));
    // Pinning a_position to slot 0 in every program lets one VAO drive all of them.
    gl.bindAttribLocation(p, 0, 'a_position');
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS))
      throw new Error(label + ': ' + gl.getProgramInfoLog(p));
    return p;
  }

  function u(p, name) {
    if (p.__id === undefined) p.__id = ++progSeq;
    const key = p.__id + ':' + name;
    if (!uCache.has(key)) uCache.set(key, gl.getUniformLocation(p, name));
    return uCache.get(key);
  }

  function tex(internal, format, type, w, h, data) {
    const t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texImage2D(gl.TEXTURE_2D, 0, internal, w, h, 0, format, type, data || null);
    // Float targets are not filterable without OES_texture_float_linear, and every read in
    // these shaders is a texelFetch anyway -- wrapping is done by hand so the sampler never
    // needs to know about the torus.
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    return t;
  }

  function fbo(attachments) {
    const f = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, f);
    attachments.forEach((t, i) =>
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + i, gl.TEXTURE_2D, t, 0));
    if (attachments.length > 1)
      gl.drawBuffers(attachments.map((_, i) => gl.COLOR_ATTACHMENT0 + i));
    const st = gl.checkFramebufferStatus(gl.FRAMEBUFFER);
    if (st !== gl.FRAMEBUFFER_COMPLETE) throw new Error('incomplete framebuffer: 0x' + st.toString(16));
    return f;
  }

  function bind(p, name, unit, t) {
    gl.activeTexture(gl.TEXTURE0 + unit);
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.uniform1i(u(p, name), unit);
  }

  function pass(p, target, w, h) {
    gl.useProgram(p);
    gl.bindFramebuffer(gl.FRAMEBUFFER, target);
    gl.viewport(0, 0, w, h);
    return p;
  }

  function drawQuad() { gl.drawArrays(gl.TRIANGLES, 0, 6); }

  function del(t, isFbo) { if (t) (isFbo ? gl.deleteFramebuffer(t) : gl.deleteTexture(t)); }

  // ------------------------------------------------------------------------------------------
  //  Kernel -- built with the same arithmetic as buildKernel() in the CPU demo, including the
  //  1e-7 tap threshold and the normalization by the *unthresholded* sum, so the weights the
  //  shader multiplies by are the ones the CPU version used.
  // ------------------------------------------------------------------------------------------
  function quad4(r) { const q = 4 * r * (1 - r); return q * q * q * q; }

  function kernelData(R, betas) {
    const KS = 2 * R + 1, raw = new Float32Array(KS * KS);
    let sum = 0;
    for (let iy = 0; iy < KS; iy++) {
      const yy = -1 + 2 * iy / (KS - 1);
      for (let ix = 0; ix < KS; ix++) {
        const xx = -1 + 2 * ix / (KS - 1);
        let r = Math.sqrt(xx * xx + yy * yy); if (r > 1) r = 0;
        let v; const b = betas.length;
        if (b > 1) {
          const Br = b * r; let idx = Math.floor(Br); if (idx > b - 1) idx = b - 1;
          v = betas[idx] * quad4(Br % 1);
        } else v = betas[0] * quad4(r);
        raw[iy * KS + ix] = v; sum += v;
      }
    }
    if (sum < 1e-6) sum = 1;
    const out = new Float32Array(KS * KS);
    for (let i = 0; i < raw.length; i++) { const w = raw[i] / sum; out[i] = w > 1e-7 ? w : 0; }
    return out;
  }

  // ------------------------------------------------------------------------------------------
  //  Public API
  // ------------------------------------------------------------------------------------------
  async function init(cv, netSize) {
    canvas = cv; net = netSize || 96;
    // preserveDrawingBuffer, because the board is not redrawn every animation frame: at low sim
    // speeds most frames step nothing, and while paused the loop stops entirely. Without it the
    // drawing buffer is cleared after compositing and the canvas would blank out between draws.
    gl = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL 2 unavailable — try another browser (see caniuse.com/webgl2).');
    // Without this the float textures the sim, reduction and crop all render into are not
    // color-renderable, and none of the passes below can exist.
    if (!gl.getExtension('EXT_color_buffer_float'))
      throw new Error('EXT_color_buffer_float unavailable — this GPU/browser cannot render to float textures.');

    const srcs = await Promise.all(NAMES.map(n =>
      fetch(SHADER_DIR + n + '.glsl').then(r => {
        if (!r.ok) throw new Error('failed to load ' + n + '.glsl (' + r.status + ') — serve over http(s), not file://');
        return r.text();
      })));
    const src = {}; NAMES.forEach((n, i) => src[n] = srcs[i]);
    ['sim', 'action', 'reduce', 'com', 'crop', 'draw'].forEach(n => prog[n] = link(src.vertex, src[n], n));

    vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER,
      new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);

    gl.disable(gl.DEPTH_TEST);
    gl.disable(gl.BLEND);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);

    cropPix = new Float32Array(net * net * 4);
    redTex0 = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID, GRID);
    redTex1 = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID, GRID);
    redFbo = fbo([redTex0, redTex1]);
    comTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, 1, 1);
    comFbo = fbo([comTex]);
    cropTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, net, net);
    cropFbo = fbo([cropTex]);
  }

  function setBoard(w, h) {
    W = w; H = h;
    canvas.width = W; canvas.height = H;

    ringTex.forEach(t => del(t, false));
    ringFbo.forEach(f => del(f, true));
    del(scratchTex, false); del(scratchFbo, true); del(wallTex, false);

    ringTex = []; ringFbo = [];
    for (let i = 0; i < RING; i++) {
      const t = tex(gl.R32F, gl.RED, gl.FLOAT, W, H);
      ringTex.push(t); ringFbo.push(fbo([t]));
    }
    ringIdx = RING - 1;
    scratchTex = tex(gl.R32F, gl.RED, gl.FLOAT, W, H);
    scratchFbo = fbo([scratchTex]);
    wallTex = tex(gl.R8, gl.RED, gl.UNSIGNED_BYTE, W, H, new Uint8Array(W * H));

    // Stage-1 blocks tile the board across a fixed GRID x GRID output. reduce.glsl caps its
    // per-block loops at 16 so the shader has a constant bound, which puts a ceiling of
    // GRID*16 = 256 on the board edge -- past that the reduction would silently miss cells.
    redBlock = Math.ceil(Math.max(W, H) / GRID);
    if (redBlock > 16) throw new Error('board edge ' + Math.max(W, H) + ' exceeds the 256 the reduction supports');
  }

  // The maze mask arrives as 0/1 bytes, which is what the JS side wants for its own masking.
  // R8 is a *normalized* format though, so a stored 1 samples as 1/255 and every "is this a
  // wall" test in the shaders would read false. Expand to 0/255 so it samples as 0.0/1.0.
  function setWall(mask) {
    const t = new Uint8Array(mask.length);
    for (let i = 0; i < mask.length; i++) t[i] = mask[i] ? 255 : 0;
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, W, H, 0, gl.RED, gl.UNSIGNED_BYTE, t);
  }

  function setRule(rule) {
    kR = rule.R;
    del(kernelTex, false);
    const KS = 2 * kR + 1;
    kernelTex = tex(gl.R32F, gl.RED, gl.FLOAT, KS, KS, kernelData(kR, rule.betas));
    ruleMu = rule.mu; ruleSigma = rule.sigma; ruleDt = rule.dt;
  }

  function uploadState(arr) {
    // Fills every frame of the ring with the same state — the equivalent of fillStack() at
    // spawn, so the policy's first input is four copies of the starting board.
    for (let i = 0; i < RING; i++) {
      gl.bindTexture(gl.TEXTURE_2D, ringTex[i]);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, W, H, 0, gl.RED, gl.FLOAT, arr);
    }
    ringIdx = RING - 1;
  }

  function runAction(srcTex, a) {
    const p = pass(prog.action, scratchFbo, W, H);
    bind(p, 'uState', 0, srcTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform2i(u(p, 'uCenter'), a.x, a.y);
    gl.uniform1f(u(p, 'uDelta'), a.delta);
    gl.uniform1i(u(p, 'uRadius'), a.radius);
    drawQuad();
  }

  function runSim(srcTex, dstFbo) {
    const p = pass(prog.sim, dstFbo, W, H);
    bind(p, 'uState', 0, srcTex);
    bind(p, 'uWall', 1, wallTex);
    bind(p, 'uKernel', 2, kernelTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uR'), kR);
    gl.uniform1f(u(p, 'uMu'), ruleMu);
    gl.uniform1f(u(p, 'uSigma'), ruleSigma);
    gl.uniform1f(u(p, 'uDt'), ruleDt);
    drawQuad();
  }

  // reduce -> com -> crop. Queued together with no readback in between; the caller collects
  // both results afterwards in one go.
  function runAnalysis() {
    let p = pass(prog.reduce, redFbo, GRID, GRID);
    bind(p, 'uState', 0, ringTex[ringIdx]);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uBlock'), redBlock);
    drawQuad();

    p = pass(prog.com, comFbo, 1, 1);
    bind(p, 'uRed0', 0, redTex0);
    bind(p, 'uRed1', 1, redTex1);
    gl.uniform2i(u(p, 'uSize'), W, H);
    drawQuad();

    p = pass(prog.crop, cropFbo, net, net);
    for (let k = 0; k < RING; k++)                        // oldest -> newest
      bind(p, 'uF' + k, k, ringTex[(ringIdx + 1 + k) % RING]);
    bind(p, 'uCom', RING, comTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uNet'), net);
    drawQuad();
  }

  function prime() { runAnalysis(); }

  function step(action) {
    let input = ringTex[ringIdx];
    if (action) { runAction(input, action); input = scratchTex; }
    const dst = (ringIdx + 1) % RING;                     // overwrites the frame aging out
    runSim(input, ringFbo[dst]);
    ringIdx = dst;
    runAnalysis();
  }

  // The one synchronization point per step. The CoM read is what actually stalls on the queued
  // passes; the crop read that follows is already resident by then.
  function readback() {
    gl.bindFramebuffer(gl.FRAMEBUFFER, comFbo);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, comPix);
    gl.bindFramebuffer(gl.FRAMEBUFFER, cropFbo);
    gl.readPixels(0, 0, net, net, gl.RGBA, gl.FLOAT, cropPix);
    return { valid: comPix[3] > 0.5, row: comPix[0], col: comPix[1], mass: comPix[2], crop: cropPix };
  }

  function draw(o) {
    const p = pass(prog.draw, null, W, H);
    bind(p, 'uState', 0, ringTex[ringIdx]);
    bind(p, 'uWall', 1, wallTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform2f(u(p, 'uTarget'), o.targetRow, o.targetCol);
    gl.uniform1f(u(p, 'uTargetSigma'), o.targetSigma);
    gl.uniform3f(u(p, 'uWallRGB'), o.wallRGB[0] / 255, o.wallRGB[1] / 255, o.wallRGB[2] / 255);
    gl.uniform3f(u(p, 'uTargetRGB'), o.targetRGB[0] / 255, o.targetRGB[1] / 255, o.targetRGB[2] / 255);
    drawQuad();
  }

  return { init, setBoard, setWall, setRule, uploadState, prime, step, readback, draw };
})();
