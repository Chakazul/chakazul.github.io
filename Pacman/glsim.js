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
  const NAMES = ['vertex', 'sim', 'action', 'reduce', 'com', 'crop', 'draw', 'eatreduce', 'eatsum',
                 'rotate', 'ghostsim', 'ghostblit', 'tilered', 'tilecom', 'dotsites'];
  const RING = 4;          // frames kept for the policy's frame stack
  const GRID = 16;         // stage-1 reduction output is GRID x GRID (see reduce.glsl)
  const EAT_THRESHOLD = 0.1;  // Pac-Man mass above this erases the dots channel at that cell
  // Frightened-window speed multipliers, applied to `dt` (the growth function's per-step time
  // increment) rather than to any explicit velocity -- nothing here tracks one, a Lenia soliton's
  // speed is an emergent property of its rule, and dt is the one knob that scales how much of a
  // step's growth is applied without touching the rule (mu/sigma/betas) itself, which is what
  // actually shapes the pattern. A frightened ghost's own tile runs at GHOST_FRIGHTENED_SPEED of
  // its ordinary dt (see runGhostSim()); Pac-Man runs at PACMAN_FRIGHTENED_SPEED of his whenever
  // *any* ghost is currently frightened (see step()), not just while he's actually near one.
  const GHOST_FRIGHTENED_SPEED = 0.6;
  const PACMAN_FRIGHTENED_SPEED = 1.0;
  // Edge of the private world each ghost is simulated in. The soliton's mass reaches ~23px from
  // its centre, growth can appear a kernel radius beyond that, and computing those cells reads
  // another radius further again -- so 96 has the room it needs, where the obvious 48 would clip
  // the pattern. It is also what bounds the cost: a tile is 9216 cells against the board's 137500,
  // so several ghosts together come to a fraction of one whole-board pass.
  const GHOST_WIN = 96;
  const MAX_GHOSTS = 9;   // must match the #define in draw.glsl, ghostblit.glsl and ghostsim.glsl
  // Pac-Man's own window. Wider than a ghost's, and for a reason a ghost does not have: he takes
  // interventions. CARL picks a cell anywhere in its netSize view -- 48px off his centre of mass --
  // and lays down a disc of actionRadius on top, so real mass can arrive 55px out. A 96 window
  // (half-width 48) would silently clip any intervention past 41px; 64 covers the whole reach with
  // room for the growth that follows it. Unlike the ghosts he keeps a board-sized texture: only the
  // *simulation* is windowed, so the policy's frame stack, the CoM and the eat checks all still see
  // an ordinary full board and need no changes at all.
  const PAC_WIN = 128;

  let gl = null, canvas = null;
  let prog = {}, uCache = new Map(), progSeq = 0;
  let W = 0, H = 0, net = 96, kR = 18;
  let ruleMu = 0.24, ruleSigma = 0.024, ruleDt = 0.1;

  let ringTex = [], ringFbo = [], ringIdx = 0;
  // Board position of Pac-Man's simulation window, tracked from the CoM the last readback returned
  // -- the engine already has it, so nothing has to be plumbed in from the caller.
  let pacWin = null;
  // Channel 2: a second Lenia field on the same board, stepped under its own rule. It gets a
  // plain ping-pong pair rather than a ring, because nothing here needs its history -- the policy
  // never sees it, so there is no frame stack to crop and no CoM to find.
  let ch2Tex = [], ch2Fbo = [], ch2Idx = 0, ch2Kernel = null;
  let ch2R = 18, ch2Mu = 0.24, ch2Sigma = 0.024, ch2Dt = 0.1;
  // Channel 3 (the ghosts). Unlike the other two this is not a board-sized field: the ghosts are
  // simulated in an atlas of private GHOST_WIN-square tiles, one each, and composited back to a
  // board-sized texture afterwards for everything downstream to read. Tiles are what stop two
  // ghosts being two solitons in one field -- which is to say, what stops them destroying each
  // other on contact -- and what makes their cost scale with the number of ghosts, not the board.
  let gAtlas = [], gAtlasFbo = [], gAtlasIdx = 0, ch3Kernel = null;
  let ch3R = 18, ch3Mu = 0.24, ch3Sigma = 0.024, ch3Dt = 0.1;
  let gBoardTex = null, gBoardFbo = null;              // the composited board-space result
  // Same composite, but mass from a currently-frightened ghost's tile excluded -- what the
  // ordinary ghosts-eat-Pac-Man check (in step(), below) reads instead of gBoardTex, so a
  // frightened ghost simply contributes nothing to it and is harmless per-tile rather than only
  // globally. A second render target of the same runGhostBlit() pass, not a separate pass.
  let gDangerTex = null;
  let gRedTex = null, gRedFbo = null;                  // per-tile CoM, stage 1
  let gComTex = null, gComFbo = null;                  // per-tile CoM, stage 2 -- one texel each
  // Board position of each tile's (0,0), how far its contents slide this step to keep the soliton
  // centred, and whether Pac-Man is currently allowed to eat it (frightened) -- all three owned by
  // the caller, which is the side that tracks which ghost is which and its frightened state.
  let gOrigin = [], gShift = [], gFrightened = [], gCount = 0;
  let scratchTex = null, scratchFbo = null;
  let wallTex = null, powerTex = null, kernelTex = null;
  let redTex0 = null, redTex1 = null, redFbo = null, redBlock = 1;
  let comTex = null, comFbo = null;
  // Total remaining dots-channel mass -- same two-stage reduce.glsl -> com.glsl shape as the
  // Pac-Man CoM above, just pointed at ch2Tex. Reused rather than measured on the CPU from
  // cumulative "eaten" mass, because the dots are themselves free-running Lenia solitons: eating
  // only part of one can leave the rest to regrow under its own growth rule, so mass is not
  // conserved and a running subtraction drifts from the field's actual total. com.glsl's mass
  // component (.z) is exactly the sum this needs; its circular-mean position (.x/.y) is unused.
  let dotsRedTex0 = null, dotsRedTex1 = null, dotsRedFbo = null, dotsComTex = null, dotsComFbo = null;
  let hasDots = false;      // false whenever there is no dots channel to have a total mass
  // The "how much dot mass did Pac-Man just eat" reduction -- same two-stage shape as
  // reduce.glsl -> com.glsl, just over a different pair of textures (see eatreduce.glsl).
  let eatRedTex = null, eatRedFbo = null, eatSumTex = null, eatSumFbo = null;
  // Per-dot presence (see dotsites.glsl): one output texel per dot, the channel-2 mass left in a
  // window around where it was stamped. The positions arrive as a texture rather than a uniform
  // array because their count is whatever the layout holds, not a compile-time constant.
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
    // The ghost atlas: MAX_GHOSTS tiles side by side, ping-ponged like any Lenia field.
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
    // Zero-filled, so the draw pass has something to sample before the first upload.
    ch2Tex = []; ch2Fbo = [];
    const empty = new Float32Array(W * H);
    for (let i = 0; i < 2; i++) {
      const t2 = tex(gl.R32F, gl.RED, gl.FLOAT, W, H, empty);
      ch2Tex.push(t2); ch2Fbo.push(fbo([t2]));
    }
    ch2Idx = 0;
    // Two channels: mass, and which ghost owns the pixel. See ghostblit.glsl.
    gBoardTex = tex(gl.RG32F, gl.RG, gl.FLOAT, W, H, new Float32Array(W * H * 2));
    // Second render target of the same pass: mass with any currently-frightened ghost's
    // contribution excluded. See gDangerTex's own comment above.
    gDangerTex = tex(gl.R32F, gl.RED, gl.FLOAT, W, H, new Float32Array(W * H));
    gBoardFbo = fbo([gBoardTex, gDangerTex]);
    scratchTex = tex(gl.R32F, gl.RED, gl.FLOAT, W, H);
    scratchFbo = fbo([scratchTex]);
    wallTex = tex(gl.R8, gl.RED, gl.UNSIGNED_BYTE, W, H, new Uint8Array(W * H));
    // 1 where placeDots() stamped a power pellet rather than an ordinary dot, 0 elsewhere -- same
    // 0/255-expansion trick as wallTex, since R8 is normalized. Dots and pellets are otherwise
    // identical mass under the one shared channel-2 rule (a Lenia pattern's size is fixed by the
    // kernel radius it runs at, so stamping one bigger just has it relax back to the same
    // equilibrium size, not stay distinct), so colour is what draw.glsl uses to tell them apart;
    // this mask is what tells draw.glsl which colour a given pixel of channel 2 gets. Static once
    // written -- the dots channel never moves, so it never needs updating between placeDots() calls.
    powerTex = tex(gl.R8, gl.RED, gl.UNSIGNED_BYTE, W, H, new Uint8Array(W * H));

    // Stage-1 blocks tile the board across a fixed GRID x GRID output, whatever the board size:
    // reduce.glsl loops on this value rather than a constant, so there is no board-edge ceiling
    // any more (it used to cap the loops at 16, limiting the edge to GRID*16 = 256).
    redBlock = Math.ceil(Math.max(W, H) / GRID);
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

  // Same 0/255 expansion as setWall(), for the power-pellet colour mask. Called by placeDots()
  // alongside uploadState2(), whenever the dots channel is (re)stamped.
  function setPowerMask(mask) {
    const t = new Uint8Array(mask.length);
    for (let i = 0; i < mask.length; i++) t[i] = mask[i] ? 255 : 0;
    gl.bindTexture(gl.TEXTURE_2D, powerTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.R8, W, H, 0, gl.RED, gl.UNSIGNED_BYTE, t);
  }

  // Where each dot was stamped, as [row, col] board positions, plus the half-width of the square
  // window summed around each. Called by placeDots() whenever the dots channel is (re)stamped;
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
    // Fills every frame of the ring with the same state — the equivalent of fillStack() at
    // spawn, so the policy's first input is four copies of the starting board.
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

  // Both channels are the same Lenia step with different weights, so one program serves both --
  // a channel is just its own kernel plus its own growth parameters. `eat`, when given, is
  // {tex, threshold}: another channel's state texture that erases this one wherever it exceeds
  // threshold -- one extra fetch-and-compare, the same shape as the wall check, rather than a
  // second convolution pass to detect overlap. `wallEnabled` skips the wall fetch entirely for a
  // channel whose solitons never move -- true stationary spots (the dots) can never reach a wall,
  // so masking against one every step is pure waste for that channel.
  // `win`, when given, is a board position to confine the step to: the destination is cleared and
  // only a PAC_WIN-square patch around that point is computed. The texture stays board-sized, so
  // every other pass -- the crop, the reduction, the eat checks, the draw -- carries on reading an
  // ordinary full board and none of them had to change. What is dropped is only ever empty: the
  // soliton spans ~46px and the window ~128, and interventions reach 55px at the very most.
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

  // One texel per dot site -- see dotsites.glsl.
  function runDotSites() {
    const p = pass(prog.dotsites, siteOutFbo, siteCount, 1);
    bind(p, 'uState', 0, ch2Tex[ch2Idx]);
    bind(p, 'uSites', 1, siteTex);
    gl.uniform2i(u(p, 'uSize'), W, H);
    gl.uniform1i(u(p, 'uRadius'), siteRadius);
    drawQuad();
  }

  // reduce -> com again, over ch2Tex instead of the ring: how much dots-channel mass is left on
  // the whole board right now -- and, per dot site, around each dot. Queued alongside
  // runAnalysis() with no readback in between.
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

  // A window near a board edge runs off it, and the board is a torus, so it comes back on the far
  // side -- while a scissor rectangle cannot wrap. Cut it into the one to four pieces that do lie
  // on the board.
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

  // The caller owns which ghost is where. `origins` are the board positions of each tile's (0,0),
  // `shifts` how far each tile's contents slide this step -- whole cells, because a fractional
  // slide would mean resampling the soliton, and resampling it every step smears it away -- and
  // `frightened` whether Pac-Man currently eats that tile on overlap instead of the other way
  // round. All three are rebuilt from the caller's own ghost list on every call (steerGhosts()
  // calls this every step), so per-tile frightened state never goes stale between calls the way a
  // separate setter touched only sometimes would.
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

  // Writes one soliton into one tile. texSubImage2D rather than a whole-atlas upload, so placing
  // or respawning a single ghost leaves every other tile running untouched.
  function uploadGhostTile(i, data) {
    if (i >= MAX_GHOSTS) return;
    gl.bindTexture(gl.TEXTURE_2D, gAtlas[gAtlasIdx]);
    gl.texSubImage2D(gl.TEXTURE_2D, 0, i * GHOST_WIN, 0, GHOST_WIN, GHOST_WIN,
                     gl.RED, gl.FLOAT, data);
  }

  // One Lenia step for every ghost, each confined to its own tile. A single draw covers the whole
  // pack -- the tiles are contiguous from 0, so one scissor bounds the used ones. `pacmanTex` is
  // Pac-Man's board-space state; a tile whose ghost is currently frightened (gFrightened, per-tile)
  // gets erased wherever it exceeds EAT_THRESHOLD there -- the frightened-only reverse of the
  // ordinary ghost-eats-Pac-Man check, and per-ghost rather than all-or-nothing so an already-eaten
  // ghost that respawned mid-window comes back not frightened while its packmates still are.
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
    // Per-tile, not a shared uDt: a frightened ghost runs slower than its (possibly still normal)
    // packmates, which a single uniform for the whole atlas couldn't express.
    setTileFloatUniforms(p, 'uDt', gFrightened.map(f => f ? ch3Dt * GHOST_FRIGHTENED_SPEED : ch3Dt));
    gl.enable(gl.SCISSOR_TEST);
    gl.scissor(0, 0, GHOST_WIN * gCount, GHOST_WIN);
    drawQuad();
    gl.disable(gl.SCISSOR_TEST);
    gAtlasIdx = dst;
  }

  // Tiles -> board space, for draw.glsl and for the eat check in Pac-Man's own sim pass. Two
  // render targets: gBoardTex (total mass + owner, for rendering) and gDangerTex (mass with any
  // frightened tile's contribution left out, for the ghosts-eat-Pac-Man check -- see its own
  // comment). Both come out of this one pass since it already sums every tile per pixel.
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

  // Turn one ghost without disturbing what it is: an exact quarter-turn permutation of its live
  // tile, rather than a fresh stamp of the canonical pattern. Pivot is tile-local and integer -- a
  // half-pixel pivot would not land texel centres on texel centres and the rotation would stop
  // being exact. See rotate.glsl.
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

  // How much dots-channel mass sits under Pac-Man's just-stepped state, i.e. how much is about
  // to be erased by the eat check inside the channel-2 sim pass below -- same two textures, same
  // threshold, just measured rather than acted on. dotsTex is channel 2's *pre*-step state (still
  // resident: ch2Tex's ping-pong means it isn't overwritten until the sim pass after this one).
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
    // The ghosts eat Pac-Man exactly the way Pac-Man eats the dots, just pointed the other way --
    // except a currently-frightened ghost contributes nothing to gDangerTex (see runGhostBlit()),
    // so it is individually harmless without anything needing to be suppressed here. This one
    // reads the ghosts' *pre*-step state, since their own sim pass is queued below -- a frame
    // staler than the dots' eat check, and unnoticeable at one Lenia step of ghost drift. It stays
    // readable through that pass regardless: ping-pong means the ghosts write the other slot of
    // their pair, not the one bound here.
    const ghostEat = (ch3Kernel && gCount) ? { tex: gDangerTex, threshold: EAT_THRESHOLD } : null;
    // Sped up while any ghost is currently frightened (gFrightened, the same per-tile flags the
    // eat checks use) -- not just while he's near one, matching the ghosts' own speed change being
    // a property of the window rather than of proximity.
    const pacDt = gFrightened.some(f => f) ? ruleDt * PACMAN_FRIGHTENED_SPEED : ruleDt;
    runSim(input, ringFbo[dst], kernelTex, kR, ruleMu, ruleSigma, pacDt, ghostEat, true, pacWin);
    ringIdx = dst;
    // Channel 2 steps on the same schedule but takes no action pass and no analysis: CARL
    // neither sees it nor steers it, so there is nothing to intervene on and nothing to read back.
    // It does eat, though: ringTex[dst] is Pac-Man's just-stepped state, so a dot is erased the
    // instant Pac-Man's mass overlaps it this same tick, not one frame late. Wall masking is
    // skipped (wallEnabled=false): the dots are stationary spots stamped only onto open floor,
    // so they can never drift into a wall, and there is nothing to check for.
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
    // The ghosts: free-running like the dots, but mobile, so unlike the dots they are wall-masked.
    // Ordinarily nothing eats them -- the eating they take part in is the ghostEat above, which
    // reads the board-space composite this leaves behind. A frightened one (per-tile, gFrightened)
    // is the other way round: Pac-Man's just-stepped state (ringTex[dst], same frame the dots' own
    // eat check reads) erases it instead, symmetric with how a dot gets erased under him.
    hasGhost = !!ch3Kernel && gCount > 0;
    if (hasGhost) { runGhostSim(ringTex[dst]); runGhostBlit(); runGhostCom(); }
    runAnalysis();
  }

  // The one synchronization point per step. The CoM read is what actually stalls on the queued
  // passes; the crop read that follows is already resident by then.
  function readback() {
    gl.bindFramebuffer(gl.FRAMEBUFFER, comFbo);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, comPix);
    // The next step's window follows him. Held at its last position when the CoM goes invalid,
    // which only happens once he has dissolved and the episode is over anyway.
    if (comPix[3] > 0.5) pacWin = [comPix[0], comPix[1]];
    let eaten = 0, pelletEaten = 0;
    if (hasEatSignal) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, eatSumFbo);
      gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, eatPix);
      eaten = eatPix[0];
      pelletEaten = eatPix[1];
    }
    // Another 1x1 read, and the pipeline is already flushed by the CoM read above, so it costs
    // essentially nothing on top of the stall that was happening anyway.
    // Every ghost in one read: their results are neighbouring texels of a single row, so the whole
    // pack costs the same one call a single ghost did. Positions come back tile-local; the caller
    // adds the tile origin, since it is the side that knows where each tile is.
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
    // Every dot site in one row, same one-call shape as the ghosts above. Reused buffer: the caller
    // reads it before the next readback() overwrites it.
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
