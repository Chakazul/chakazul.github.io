"use strict";
// ============================================================================================
//  GPU Lenia engine for the CARL maze demo.
//
//  The CPU demo (CARL/maze_playground.html) ran the whole sim in JS: FFT convolution per step,
//  a full-board sweep for the centre of mass, a per-pixel canvas paint. All three are
//  per-cell-independent, so all three are fragment shader passes here; only the policy network
//  (onnxruntime-web/WASM) stays on the CPU.
//
//  Every readback is a pipeline stall, and the CPU needs exactly two things back each step:
//  where the soliton is, and the 96x96x4 window the policy reads. So a step queues five passes
//  with no synchronization between them (action, sim, reduce, com, crop), and readback()
//  collects both results in one call. crop.glsl reads the CoM out of a texture rather than a
//  uniform for the same reason -- it doesn't have to wait to be told where to look.
//
//  Board state is a ring of four R32F textures, not a ping-pong pair, because the policy needs
//  a 4-frame stack all croppable at a shared origin. Stepping overwrites the four-steps-ago
//  frame -- exactly the one falling out of the stack.
// ============================================================================================

const SimGL = (function () {
  const SHADER_DIR = './shaders/';
  const NAMES = ['vertex', 'sim', 'action', 'reduce', 'com', 'crop', 'draw', 'eatreduce', 'eatsum',
                 'rotate', 'ghostsim', 'ghostblit', 'tilered', 'tilecom', 'dotsites'];
  const RING = 4;          // frames kept for the policy's frame stack
  const GRID = 16;         // stage-1 reduction output is GRID x GRID (see reduce.glsl)
  const EAT_THRESHOLD = 0.1;  // Pac-Man mass above this erases the dots channel at that cell
  // Frightened-window speed multipliers, applied to `dt` since a Lenia soliton's speed is an
  // emergent property of its rule and dt is the knob that scales growth without touching the
  // rule itself. A frightened ghost's tile runs at GHOST_FRIGHTENED_SPEED of normal (see
  // runGhostSim()); Pac-Man runs at PACMAN_FRIGHTENED_SPEED whenever *any* ghost is frightened.
  const GHOST_FRIGHTENED_SPEED = 0.6;
  const PACMAN_FRIGHTENED_SPEED = 1.0;
  // Edge of each ghost's private tile. The soliton's mass reaches ~23px from centre plus a kernel
  // radius of growth plus another radius of taps reading that -- 96 has the room; 48 would clip.
  // Also bounds cost: 9216 cells per tile against the board's 137500.
  const GHOST_WIN = 96;
  const MAX_GHOSTS = 9;   // must match the #define in draw.glsl, ghostblit.glsl and ghostsim.glsl
  // Pac-Man's own window, wider than a ghost's because he takes interventions: CARL can act up to
  // 48px off centre plus a 7px action radius, so real mass can arrive 55px out -- a 96 window
  // (half-width 48) would clip that. Unlike the ghosts he keeps a board-sized texture; only the
  // simulation is windowed, so the frame stack, CoM and eat checks need no changes.
  const PAC_WIN = 128;

  let gl = null, canvas = null;
  let prog = {}, uCache = new Map(), progSeq = 0;
  let W = 0, H = 0, net = 96, kR = 18;
  let ruleMu = 0.24, ruleSigma = 0.024, ruleDt = 0.1;

  let ringTex = [], ringFbo = [], ringIdx = 0;
  // Board position of Pac-Man's simulation window, tracked from the last readback's CoM.
  let pacWin = null;
  // Channel 2 (dots): a second Lenia field, own rule, plain ping-pong pair rather than a ring --
  // the policy never sees it, so there's no frame stack to crop and no CoM to find.
  let ch2Tex = [], ch2Fbo = [], ch2Idx = 0, ch2Kernel = null;
  let ch2R = 18, ch2Mu = 0.24, ch2Sigma = 0.024, ch2Dt = 0.1;
  // Channel 3 (ghosts): not a board-sized field like the other two -- each ghost is simulated in
  // its own private GHOST_WIN-square tile of one atlas, composited to board space afterwards.
  // Tiles are what stop two ghosts being two solitons that destroy each other on contact, and
  // what makes cost scale with ghost count rather than board size.
  let gAtlas = [], gAtlasFbo = [], gAtlasIdx = 0, ch3Kernel = null;
  let ch3R = 18, ch3Mu = 0.24, ch3Sigma = 0.024, ch3Dt = 0.1;
  let gBoardTex = null, gBoardFbo = null;              // the composited board-space result
  // Same composite with a currently-frightened ghost's mass excluded -- what the ghosts-eat-
  // Pac-Man check (step(), below) reads instead, so a frightened ghost is harmless per-tile.
  // Second render target of runGhostBlit(), not a separate pass.
  let gDangerTex = null;
  let gRedTex = null, gRedFbo = null;                  // per-tile CoM, stage 1
  let gComTex = null, gComFbo = null;                  // per-tile CoM, stage 2 -- one texel each
  // Per-ghost, owned by the caller: board position of each tile's (0,0), how far it slides this
  // step to stay centred, and whether Pac-Man may currently eat it (frightened).
  let gOrigin = [], gShift = [], gFrightened = [], gCount = 0;
  let scratchTex = null, scratchFbo = null;
  let wallTex = null, powerTex = null, kernelTex = null;
  let redTex0 = null, redTex1 = null, redFbo = null, redBlock = 1;
  let comTex = null, comFbo = null;
  // Total remaining dots-channel mass, same reduce.glsl -> com.glsl shape as the Pac-Man CoM but
  // pointed at ch2Tex. Reused each time rather than accumulated from "eaten" mass on the CPU,
  // since the dots regrow under their own rule and mass isn't conserved. Only com.glsl's mass
  // component (.z) is used; its position (.x/.y) is not.
  let dotsRedTex0 = null, dotsRedTex1 = null, dotsRedFbo = null, dotsComTex = null, dotsComFbo = null;
  let hasDots = false;      // false whenever there is no dots channel to have a total mass
  // "How much dot mass did Pac-Man just eat", same two-stage shape, see eatreduce.glsl.
  let eatRedTex = null, eatRedFbo = null, eatSumTex = null, eatSumFbo = null;
  // Per-dot presence (dotsites.glsl): one texel per dot, channel-2 mass left around its stamp.
  // Positions arrive as a texture, not a uniform array, since their count isn't a compile-time constant.
  let siteTex = null, siteOutTex = null, siteOutFbo = null, siteCount = 0, siteRadius = 0;
  let sitePix = null, siteMass = null;
  let cropTex = null, cropFbo = null;
  let vao = null;

  const comPix = new Float32Array(4);
  const gComPix = new Float32Array(4 * MAX_GHOSTS);
  const eatPix = new Float32Array(4);
  const dotsComPix = new Float32Array(4);
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
    // Float targets aren't filterable without OES_texture_float_linear, and every read here is a
    // texelFetch anyway, so wrapping is done by hand rather than via the sampler.
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
  //  Kernel -- same arithmetic as the CPU demo's buildKernel() (1e-7 tap threshold, normalized
  //  by the unthresholded sum), so the weights match.
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
    // preserveDrawingBuffer: the board isn't redrawn every frame (paused, or low sim speed), and
    // without this the drawing buffer clears after compositing and the canvas blanks between draws.
    gl = canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: true });
    if (!gl) throw new Error('WebGL 2 unavailable — try another browser (see caniuse.com/webgl2).');
    // Needed for the sim/reduction/crop float textures to be colour-renderable at all.
    if (!gl.getExtension('EXT_color_buffer_float'))
      throw new Error('EXT_color_buffer_float unavailable — this GPU/browser cannot render to float textures.');

    const srcs = await Promise.all(NAMES.map(n =>
      fetch(SHADER_DIR + n + '.glsl').then(r => {
        if (!r.ok) throw new Error('failed to load ' + n + '.glsl (' + r.status + ') — serve over http(s), not file://');
        return r.text();
      })));
    const src = {}; NAMES.forEach((n, i) => src[n] = srcs[i]);
    ['sim', 'action', 'reduce', 'com', 'crop', 'draw', 'eatreduce', 'eatsum', 'rotate',
     'ghostsim', 'ghostblit', 'tilered', 'tilecom', 'dotsites']
      .forEach(n => prog[n] = link(src.vertex, src[n], n));

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
    dotsRedTex0 = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID, GRID);
    dotsRedTex1 = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID, GRID);
    dotsRedFbo = fbo([dotsRedTex0, dotsRedTex1]);
    dotsComTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, 1, 1);
    dotsComFbo = fbo([dotsComTex]);
    // Ghost atlas: MAX_GHOSTS tiles side by side, ping-ponged like any Lenia field.
    for (let i = 0; i < 2; i++) {
      const t = tex(gl.R32F, gl.RED, gl.FLOAT, GHOST_WIN * MAX_GHOSTS, GHOST_WIN,
                    new Float32Array(GHOST_WIN * MAX_GHOSTS * GHOST_WIN));
      gAtlas.push(t); gAtlasFbo.push(fbo([t]));
    }
    gRedTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID * MAX_GHOSTS, GRID);
    gRedFbo = fbo([gRedTex]);
    gComTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, MAX_GHOSTS, 1);
    gComFbo = fbo([gComTex]);
    eatRedTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, GRID, GRID);
    eatRedFbo = fbo([eatRedTex]);
    eatSumTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, 1, 1);
    eatSumFbo = fbo([eatSumTex]);
    cropTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, net, net);
    cropFbo = fbo([cropTex]);
  }

  function setBoard(w, h) {
    W = w; H = h;
    canvas.width = W; canvas.height = H;

    ringTex.forEach(t => del(t, false));
    ringFbo.forEach(f => del(f, true));
    ch2Tex.forEach(t => del(t, false));
    ch2Fbo.forEach(f => del(f, true));
    del(gBoardTex, false); del(gDangerTex, false); del(gBoardFbo, true);
    del(scratchTex, false); del(scratchFbo, true); del(wallTex, false); del(powerTex, false);

    ringTex = []; ringFbo = [];
    for (let i = 0; i < RING; i++) {
      const t = tex(gl.R32F, gl.RED, gl.FLOAT, W, H);
      ringTex.push(t); ringFbo.push(fbo([t]));
    }
    ringIdx = RING - 1;
    // Zero-filled so the draw pass has something to sample before the first upload.
    ch2Tex = []; ch2Fbo = [];
    const empty = new Float32Array(W * H);
    for (let i = 0; i < 2; i++) {
      const t2 = tex(gl.R32F, gl.RED, gl.FLOAT, W, H, empty);
      ch2Tex.push(t2); ch2Fbo.push(fbo([t2]));
    }
    ch2Idx = 0;
    // Two channels: mass, and which ghost owns the pixel. See ghostblit.glsl.
    gBoardTex = tex(gl.RG32F, gl.RG, gl.FLOAT, W, H, new Float32Array(W * H * 2));
    // Second render target of the same pass: mass with any frightened ghost excluded.
    gDangerTex = tex(gl.R32F, gl.RED, gl.FLOAT, W, H, new Float32Array(W * H));
    gBoardFbo = fbo([gBoardTex, gDangerTex]);
    scratchTex = tex(gl.R32F, gl.RED, gl.FLOAT, W, H);
    scratchFbo = fbo([scratchTex]);
    wallTex = tex(gl.R8, gl.RED, gl.UNSIGNED_BYTE, W, H, new Uint8Array(W * H));
    // 1 where a power pellet was stamped rather than a plain dot, 0 elsewhere -- dots and pellets
    // are identical mass under the shared channel-2 rule, so this is what draw.glsl colours from.
    // Static once written: the dots channel never moves between placeDots() calls.
    powerTex = tex(gl.R8, gl.RED, gl.UNSIGNED_BYTE, W, H, new Uint8Array(W * H));

    // Stage-1 blocks tile the board across a fixed GRID x GRID output, whatever the board size --
    // reduce.glsl loops on this value rather than a constant 16, which used to cap the edge at 256.
    redBlock = Math.ceil(Math.max(W, H) / GRID);
  }

  // R8 is a normalized format, so a stored 1 would sample as 1/255 (false in the shaders' "is
  // this a wall" tests). Expand the incoming 0/1 mask to 0/255 so it samples as 0.0/1.0.
  function setWall(mask) {
    const t = new Uint8Array(mask.length);
    for (let i = 0; i < mask.length; i++) t[i] = mask[i] ? 255 : 0;
    gl.bindTexture(gl.TEXTURE_2D, wallTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, W, H, 0, gl.RED, gl.UNSIGNED_BYTE, t);
  }

  // Same 0/255 expansion as setWall(), for the power-pellet colour mask.
  function setPowerMask(mask) {
    const t = new Uint8Array(mask.length);
    for (let i = 0; i < mask.length; i++) t[i] = mask[i] ? 255 : 0;
    gl.bindTexture(gl.TEXTURE_2D, powerTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, W, H, 0, gl.RED, gl.UNSIGNED_BYTE, t);
  }

  // Where each dot was stamped ([row, col]) plus the half-width of the window summed around each.
  // readback()'s dotSites comes back in this same order. An empty list turns the pass off.
  function setDotSites(sites, radius) {
    del(siteTex, false); del(siteOutTex, false); del(siteOutFbo, true);
    siteTex = siteOutTex = siteOutFbo = null;
    siteCount = sites.length; siteRadius = radius;
    if (!siteCount) return;
    if (siteCount > gl.getParameter(gl.MAX_TEXTURE_SIZE))
      throw new Error(`${siteCount} dots exceed this GPU's texture width limit`);
    const pos = new Float32Array(siteCount * 4);
    sites.forEach(([r, c], i) => { pos[i * 4] = ((c % W) + W) % W; pos[i * 4 + 1] = ((r % H) + H) % H; });
    siteTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, siteCount, 1, pos);
    siteOutTex = tex(gl.RGBA32F, gl.RGBA, gl.FLOAT, siteCount, 1);
    siteOutFbo = fbo([siteOutTex]);
    sitePix = new Float32Array(siteCount * 4);
    siteMass = new Float32Array(siteCount);
  }

  function setRule(rule) {
    kR = rule.R;
    del(kernelTex, false);
    const KS = 2 * kR + 1;
    kernelTex = tex(gl.R32F, gl.RED, gl.FLOAT, KS, KS, kernelData(kR, rule.betas));
    ruleMu = rule.mu; ruleSigma = rule.sigma; ruleDt = rule.dt;
  }

  function setRule2(rule) {
    ch2R = rule.R;
    del(ch2Kernel, false);
    const KS = 2 * ch2R + 1;
    ch2Kernel = tex(gl.R32F, gl.RED, gl.FLOAT, KS, KS, kernelData(ch2R, rule.betas));
    ch2Mu = rule.mu; ch2Sigma = rule.sigma; ch2Dt = rule.dt;
  }

  function setRule3(rule) {
    ch3R = rule.R;
    del(ch3Kernel, false);
    const KS = 2 * ch3R + 1;
    ch3Kernel = tex(gl.R32F, gl.RED, gl.FLOAT, KS, KS, kernelData(ch3R, rule.betas));
    ch3Mu = rule.mu; ch3Sigma = rule.sigma; ch3Dt = rule.dt;
  }

  function uploadState(arr) {
    // Fills every ring frame with the same state, so the policy's first input is four copies of
    // the starting board.
    for (let i = 0; i < RING; i++) {
      gl.bindTexture(gl.TEXTURE_2D, ringTex[i]);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, W, H, 0, gl.RED, gl.FLOAT, arr);
    }
    ringIdx = RING - 1;
  }

  function uploadState2(arr) {
    gl.bindTexture(gl.TEXTURE_2D, ch2Tex[ch2Idx]);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R32F, W, H, 0, gl.RED, gl.FLOAT, arr);
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

  // Both channels are the same Lenia step with different weights, so one program serves both.
  // `eat`, when given, is {tex, threshold}: another channel's state that erases this one wherever
  // it exceeds threshold -- one fetch-and-compare rather than a second convolution to detect
  // overlap. `wallEnabled` skips the wall fetch for a channel whose solitons never move (the dots
  // can never reach a wall). `win`, when given, confines the step to a PAC_WIN-square patch
  // around that board position, clearing the rest of the destination -- every other pass keeps
  // reading an ordinary full board, since what's dropped is always empty (soliton ~46px, window
  // ~128px, interventions reach at most 55px out).
  function runSim(srcTex, dstFbo, kern, R, mu, sigma, dt, eat, wallEnabled, win) {
    if (win) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, dstFbo);
      gl.viewport(0, 0, W, H);
      gl.clearBufferfv(gl.COLOR, 0, [0, 0, 0, 1]);
    }
    const p = pass(prog.sim, dstFbo, W, H);
    bind(p, 'uState', 0, srcTex);
    bind(p, 'uWall', 1, wallTex);
    bind(p, 'uKernel', 2, kern);
    bind(p, 'uEat', 3, eat ? eat.tex : srcTex);   // unit needs a valid texture even when disabled
    gl.uniform1i(u(p, 'uEatEnabled'), eat ? 1 : 0);
    gl.uniform1f(u(p, 'uEatThreshold'), eat ? eat.threshold : 0.0);
    gl.uniform1i(u(p, 'uWallEnabled'), wallEnabled ? 1 : 0);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uR'), R);
    gl.uniform1f(u(p, 'uMu'), mu);
    gl.uniform1f(u(p, 'uSigma'), sigma);
    gl.uniform1f(u(p, 'uDt'), dt);
    if (!win) { drawQuad(); return; }
    gl.enable(gl.SCISSOR_TEST);
    for (const [x, y, w, h] of windowRects(win[0], win[1], PAC_WIN)) { gl.scissor(x, y, w, h); drawQuad(); }
    gl.disable(gl.SCISSOR_TEST);
  }

  // reduce -> com -> crop, queued with no readback in between; the caller collects both results
  // afterwards in one go.
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

  // One texel per dot site -- see dotsites.glsl.
  function runDotSites() {
    const p = pass(prog.dotsites, siteOutFbo, siteCount, 1);
    bind(p, 'uState', 0, ch2Tex[ch2Idx]);
    bind(p, 'uSites', 1, siteTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uRadius'), siteRadius);
    drawQuad();
  }

  // reduce -> com again, over ch2Tex instead of the ring: total dots-channel mass left, plus
  // per-dot-site mass. Queued alongside runAnalysis() with no readback in between.
  function runDotsAnalysis() {
    let p = pass(prog.reduce, dotsRedFbo, GRID, GRID);
    bind(p, 'uState', 0, ch2Tex[ch2Idx]);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uBlock'), redBlock);
    drawQuad();

    p = pass(prog.com, dotsComFbo, 1, 1);
    bind(p, 'uRed0', 0, dotsRedTex0);
    bind(p, 'uRed1', 1, dotsRedTex1);
    gl.uniform2i(u(p, 'uSize'), W, H);
    drawQuad();

    if (siteCount) runDotSites();
  }

  // A window near a board edge wraps toroidally, but a scissor rectangle can't -- cut it into the
  // one to four pieces that do lie on the board.
  function spans(start, len, size) {
    const s0 = ((start % size) + size) % size;
    return s0 + len <= size ? [[s0, len]] : [[s0, size - s0], [0, len - (size - s0)]];
  }
  function windowRects(cy, cx, win) {
    const half = win >> 1, out = [];
    for (const [x, w] of spans(Math.round(cx) - half, win, W))
      for (const [y, h] of spans(Math.round(cy) - half, win, H)) out.push([x, y, w, h]);
    return out;
  }

  // `origins` are each tile's board (0,0), `shifts` how far it slides this step (whole cells,
  // since a fractional slide would resample the soliton and smear it), `frightened` whether
  // Pac-Man currently eats that tile instead of the reverse. Rebuilt from the caller's ghost list
  // on every call (steerGhosts() calls this every step), so frightened state never goes stale.
  function setGhostTiles(origins, shifts, frightened) {
    gCount = Math.min(origins.length, MAX_GHOSTS);
    gOrigin = origins.slice(0, gCount).map(o => [Math.round(o[0]), Math.round(o[1])]);
    gShift = shifts ? shifts.slice(0, gCount).map(v => [Math.round(v[0]), Math.round(v[1])])
                    : gOrigin.map(() => [0, 0]);
    gFrightened = frightened ? frightened.slice(0, gCount).map(v => v ? 1 : 0) : gOrigin.map(() => 0);
  }

  function setTileUniforms(p, name, list) {
    for (let i = 0; i < gCount; i++) gl.uniform2i(u(p, `${name}[${i}]`), list[i][1], list[i][0]);
  }

  function setTileFlagUniforms(p, name, list) {
    for (let i = 0; i < gCount; i++) gl.uniform1i(u(p, `${name}[${i}]`), list[i]);
  }

  function setTileFloatUniforms(p, name, list) {
    for (let i = 0; i < gCount; i++) gl.uniform1f(u(p, `${name}[${i}]`), list[i]);
  }

  // texSubImage2D rather than a whole-atlas upload, so placing/respawning one ghost leaves every
  // other tile untouched.
  function uploadGhostTile(i, data) {
    if (i >= MAX_GHOSTS) return;
    gl.bindTexture(gl.TEXTURE_2D, gAtlas[gAtlasIdx]);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, i * GHOST_WIN, 0, GHOST_WIN, GHOST_WIN,
                     gl.RED, gl.FLOAT, data);
  }

  // One Lenia step for every ghost, each confined to its own tile. A single draw covers the whole
  // pack -- tiles are contiguous from 0, so one scissor bounds the used ones. A tile whose ghost
  // is frightened (gFrightened, per-tile) is erased wherever `pacmanTex` exceeds EAT_THRESHOLD --
  // the reverse of the ordinary ghost-eats-Pac-Man check -- so a respawned ghost mid-window comes
  // back unfrightened while its packmates stay frightened.
  function runGhostSim(pacmanTex) {
    const dst = 1 - gAtlasIdx;
    const p = pass(prog.ghostsim, gAtlasFbo[dst], GHOST_WIN * MAX_GHOSTS, GHOST_WIN);
    bind(p, 'uAtlas', 0, gAtlas[gAtlasIdx]);
    bind(p, 'uWall', 1, wallTex);
    bind(p, 'uKernel', 2, ch3Kernel);
    bind(p, 'uPacman', 3, pacmanTex);
    gl.uniform1f(u(p, 'uEatThreshold'), EAT_THRESHOLD);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uWin'), GHOST_WIN);
    gl.uniform1i(u(p, 'uR'), ch3R);
    gl.uniform1f(u(p, 'uMu'), ch3Mu);
    gl.uniform1f(u(p, 'uSigma'), ch3Sigma);
    setTileUniforms(p, 'uOrigin', gOrigin);
    setTileUniforms(p, 'uShift', gShift);
    setTileFlagUniforms(p, 'uFrightened', gFrightened);
    // Per-tile, not a shared uDt: a frightened ghost runs slower than its still-normal packmates.
    setTileFloatUniforms(p, 'uDt', gFrightened.map(f => f ? ch3Dt * GHOST_FRIGHTENED_SPEED : ch3Dt));
    gl.enable(gl.SCISSOR_TEST);
    gl.scissor(0, 0, GHOST_WIN * gCount, GHOST_WIN);
    drawQuad();
    gl.disable(gl.SCISSOR_TEST);
    gAtlasIdx = dst;
  }

  // Tiles -> board space, for draw.glsl and for the eat check in Pac-Man's own sim pass. Two
  // render targets from one pass: gBoardTex (total mass + owner) and gDangerTex (frightened
  // tiles excluded).
  function runGhostBlit() {
    const p = pass(prog.ghostblit, gBoardFbo, W, H);
    bind(p, 'uAtlas', 0, gAtlas[gAtlasIdx]);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uWin'), GHOST_WIN);
    gl.uniform1i(u(p, 'uCount'), gCount);
    setTileUniforms(p, 'uOrigin', gOrigin);
    setTileFlagUniforms(p, 'uFrightened', gFrightened);
    drawQuad();
  }

  // Where each ghost sits inside its own tile. Two draws for the whole pack, not two per ghost:
  // the tiles' partials sit side by side, so one reduction covers all of them.
  function runGhostCom() {
    let p = pass(prog.tilered, gRedFbo, GRID * MAX_GHOSTS, GRID);
    bind(p, 'uAtlas', 0, gAtlas[gAtlasIdx]);
    gl.uniform1i(u(p, 'uWin'), GHOST_WIN);
    gl.uniform1i(u(p, 'uBlock'), Math.ceil(GHOST_WIN / GRID));
    drawQuad();

    p = pass(prog.tilecom, gComFbo, MAX_GHOSTS, 1);
    bind(p, 'uRed', 0, gRedTex);
    drawQuad();
  }

  // Exact quarter-turn permutation of the ghost's live tile, not a fresh stamp -- see rotate.glsl.
  // Pivot is tile-local and integer; a half-pixel pivot would break the exact texel-to-texel mapping.
  function rotateGhost(i, lx, ly, turns) {
    if (!ch3Kernel || i >= gCount || !(turns % 4)) return;
    const dst = 1 - gAtlasIdx;
    const p = pass(prog.rotate, gAtlasFbo[dst], GHOST_WIN * MAX_GHOSTS, GHOST_WIN);
    bind(p, 'uAtlas', 0, gAtlas[gAtlasIdx]);
    gl.uniform1i(u(p, 'uWin'), GHOST_WIN);
    gl.uniform1i(u(p, 'uTile'), i);
    gl.uniform2i(u(p, 'uPivot'), lx, ly);
    gl.uniform1i(u(p, 'uTurns'), turns % 4);
    drawQuad();
    gAtlasIdx = dst;
  }

  let hasGhost = false;       // false whenever there is no ghost to locate

  function prime() {
    runAnalysis();
    hasGhost = !!ch3Kernel && gCount > 0;
    if (hasGhost) { runGhostBlit(); runGhostCom(); }
    hasDots = !!ch2Kernel;
    if (hasDots) runDotsAnalysis();
  }

  // How much dots-channel mass is about to be erased by the channel-2 sim pass's eat check below
  // -- same two textures, same threshold, just measured rather than acted on. dotsTex is channel
  // 2's pre-step state, still resident until the sim pass after this one overwrites it.
  function runEatDetect(dotsTex, pacmanTex) {
    let p = pass(prog.eatreduce, eatRedFbo, GRID, GRID);
    bind(p, 'uDots', 0, dotsTex);
    bind(p, 'uPacman', 1, pacmanTex);
    bind(p, 'uPower', 2, powerTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uBlock'), redBlock);
    gl.uniform1f(u(p, 'uThreshold'), EAT_THRESHOLD);
    drawQuad();

    p = pass(prog.eatsum, eatSumFbo, 1, 1);
    bind(p, 'uRed', 0, eatRedTex);
    drawQuad();
  }

  let hasEatSignal = false;   // false whenever there is no dots channel to have eaten anything

  function step(action) {
    let input = ringTex[ringIdx];
    if (action) { runAction(input, action); input = scratchTex; }
    const dst = (ringIdx + 1) % RING;                     // overwrites the frame aging out
    // Ghosts eat Pac-Man the way he eats the dots, just reversed -- a frightened ghost contributes
    // nothing to gDangerTex (runGhostBlit()), so it's harmless without special-casing here. Reads
    // the ghosts' pre-step state (their own sim pass is queued below); one Lenia step staler than
    // the dots' own eat check, which doesn't matter at this drift rate.
    const ghostEat = (ch3Kernel && gCount) ? { tex: gDangerTex, threshold: EAT_THRESHOLD } : null;
    // Sped up while any ghost is frightened, not just while he's near one -- the speed change is a
    // property of the window, matching the ghosts' own.
    const pacDt = gFrightened.some(f => f) ? ruleDt * PACMAN_FRIGHTENED_SPEED : ruleDt;
    runSim(input, ringFbo[dst], kernelTex, kR, ruleMu, ruleSigma, pacDt, ghostEat, true, pacWin);
    ringIdx = dst;
    // Channel 2 (dots) steps on the same schedule with no action pass and no analysis -- CARL
    // never sees or steers it. It does eat: ringTex[dst] is Pac-Man's just-stepped state, so a dot
    // is erased the instant his mass overlaps it. Wall masking is skipped (dots never move, so
    // they can never drift into a wall).
    hasEatSignal = !!ch2Kernel;
    hasDots = !!ch2Kernel;
    if (ch2Kernel) {
      const d2 = 1 - ch2Idx;
      runEatDetect(ch2Tex[ch2Idx], ringTex[dst]);
      runSim(ch2Tex[ch2Idx], ch2Fbo[d2], ch2Kernel, ch2R, ch2Mu, ch2Sigma, ch2Dt,
        { tex: ringTex[dst], threshold: EAT_THRESHOLD }, false);
      ch2Idx = d2;
      runDotsAnalysis();     // measures the just-erased state, so a win reads the same step it happens
    }
    // Ghosts: free-running like the dots but mobile, so wall-masked unlike them. Ordinarily
    // nothing eats them (the ghostEat above reads the board-space composite they leave behind); a
    // frightened one is instead erased by Pac-Man's just-stepped state, symmetric with a dot.
    hasGhost = !!ch3Kernel && gCount > 0;
    if (hasGhost) { runGhostSim(ringTex[dst]); runGhostBlit(); runGhostCom(); }
    runAnalysis();
  }

  // The one synchronization point per step. The CoM read is what actually stalls on the queued
  // passes; the crop read that follows is already resident by then.
  function readback() {
    gl.bindFramebuffer(gl.FRAMEBUFFER, comFbo);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, comPix);
    // The next step's window follows him; held at its last position once he's dissolved (episode over).
    if (comPix[3] > 0.5) pacWin = [comPix[0], comPix[1]];
    let eaten = 0, pelletEaten = 0;
    if (hasEatSignal) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, eatSumFbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, eatPix);
      eaten = eatPix[0];
      pelletEaten = eatPix[1];
    }
    // Every ghost in one read (neighbouring texels of one row), same cost as reading a single
    // ghost. Positions come back tile-local; the caller adds the tile origin.
    const ghosts = [];
    if (hasGhost) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, gComFbo);
      gl.readPixels(0, 0, MAX_GHOSTS, 1, gl.RGBA, gl.FLOAT, gComPix);
      for (let i = 0; i < gCount; i++)
        ghosts.push({ valid: gComPix[i * 4 + 3] > 0.5, localRow: gComPix[i * 4],
                      localCol: gComPix[i * 4 + 1], mass: gComPix[i * 4 + 2] });
    }
    let dotsMass = 0;
    if (hasDots) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, dotsComFbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, dotsComPix);
      dotsMass = dotsComPix[2];
    }
    // Every dot site in one row, same shape as the ghosts above. Reused buffer: the caller reads
    // it before the next readback() overwrites it.
    let dotSites = null;
    if (hasDots && siteCount) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, siteOutFbo);
      gl.readPixels(0, 0, siteCount, 1, gl.RGBA, gl.FLOAT, sitePix);
      for (let i = 0; i < siteCount; i++) siteMass[i] = sitePix[i * 4];
      dotSites = siteMass;
    }
    gl.bindFramebuffer(gl.FRAMEBUFFER, cropFbo);
    gl.readPixels(0, 0, net, net, gl.RGBA, gl.FLOAT, cropPix);
    return {
      valid: comPix[3] > 0.5, row: comPix[0], col: comPix[1], mass: comPix[2], crop: cropPix, eaten,
      pelletEaten, ghosts, hasDots, dotsMass, dotSites,
    };
  }

  function draw(o) {
    const p = pass(prog.draw, null, W, H);
    bind(p, 'uState', 0, ringTex[ringIdx]);
    bind(p, 'uWall', 1, wallTex);
    bind(p, 'uState2', 2, ch2Tex[ch2Idx]);
    bind(p, 'uState3', 3, gBoardTex);
    bind(p, 'uPower', 4, powerTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    const rgb = (name, c) => gl.uniform3f(u(p, name), c[0] / 255, c[1] / 255, c[2] / 255);
    rgb('uWallRGB', o.wall);
    rgb('uSolitonRGB', o.soliton);
    rgb('uSoliton2RGB', o.soliton2);
    rgb('uPelletRGB', o.pellet || o.soliton2);
    rgb('uSoliton3RGB', o.soliton3);
    rgb('uBackRGB', o.background);
    const cols = o.ghosts || [];
    for (let i = 0; i < MAX_GHOSTS; i++) {
      const c = cols.length ? cols[i % cols.length] : o.soliton3;   // cycled: every ghost gets one
      gl.uniform3f(u(p, `uGhostRGB[${i}]`), c[0] / 255, c[1] / 255, c[2] / 255);
    }
    drawQuad();
  }

  return {
    init, setBoard, setWall, setPowerMask, setDotSites, setRule, setRule2, setRule3,
    uploadState, uploadState2, setGhostTiles, uploadGhostTile, rotateGhost,
    prime, step, readback, draw, GHOST_WIN,
  };
})();
