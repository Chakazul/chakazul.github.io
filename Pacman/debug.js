"use strict";
// ============================================================================================
//  ?debug=1 -- on-device performance instrumentation (see PERF.md, "Measuring on a device").
//
//  Inert without the param. With it, wraps app.js/glsim.js functions from the outside (no
//  change to their code paths) to time every step, and adds a small panel with:
//   - a live readout: inference vs. sim+readback vs. frame time, split by whether CARL inferred;
//   - Bench: pauses the game and times the policy alone, per execution provider (webgpu, wasm)
//     and per crop size (96/64/48), plus whether inference blocks the main thread;
//   - Sweep threads: reloads the page with ?threads=1,2,4,max and benches each, since onnxruntime's
//     WASM thread count is fixed for the life of a page;
//   - Copy report: everything above as one JSON blob, kept across reloads in localStorage.
// ============================================================================================
(function () {
  const params = new URLSearchParams(location.search);
  if (!params.has('debug') || params.get('debug') === '0') return;

  const BENCH_KEY = 'carl-debug-bench';   // localStorage: bench results, accumulated across reloads
  const SWEEP_KEY = 'carl-debug-sweep';   // sessionStorage: thread counts still to visit
  const now = () => performance.now();
  const store = (s, k, v) => { try { v === undefined ? s.removeItem(k) : s.setItem(k, JSON.stringify(v)); } catch (e) {} };
  const load = (s, k, d) => { try { return JSON.parse(s.getItem(k)) ?? d; } catch (e) { return d; } };

  // ------------------------------------------------------------------------------------------
  //  Live timing series
  // ------------------------------------------------------------------------------------------
  const MAX_SAMPLES = 4000;
  const series = {};
  function rec(name, v) {
    const a = series[name] || (series[name] = []);
    a.push(v);
    if (a.length > 2 * MAX_SAMPLES) a.splice(0, a.length - MAX_SAMPLES);
  }
  function summarize(a) {
    if (!a || !a.length) return null;
    const s = a.slice(-MAX_SAMPLES).sort((x, y) => x - y), q = p => s[Math.min(s.length - 1, Math.floor(p * s.length))];
    const r = x => Math.round(x * 100) / 100;
    return { n: s.length, mean: r(s.reduce((x, y) => x + y, 0) / s.length), p50: r(q(0.5)), p90: r(q(0.9)), p99: r(q(0.99)), max: r(s[s.length - 1]) };
  }

  let benching = false, inferred = false, stepsThisFrame = 0;

  // Top-level function declarations in app.js are global bindings, so reassigning them here is
  // what app.js's own call sites pick up.
  const origStep = agentStep;
  agentStep = async function () {
    inferred = false;
    const t = now();
    await origStep();
    rec(inferred ? 'step.carl' : 'step.idle', now() - t);
    stepsThisFrame++;
  };
  const origAct = agentAct;
  agentAct = async function () {
    inferred = true;
    const t = now();
    const r = await origAct();
    rec('infer', now() - t);   // de-interleave + session.run + argmax
    return r;
  };
  const origRender = render;
  render = function () { const t = now(); origRender(); rec('render', now() - t); };

  const origSimStep = SimGL.step, origReadback = SimGL.readback;
  SimGL.step = function (a) { const t = now(); const r = origSimStep(a); rec('gpu.submit', now() - t); return r; };
  SimGL.readback = function () { const t = now(); const r = origReadback(); rec('gpu.readback', now() - t); return r; };

  const origRun = ort.InferenceSession.prototype.run;
  ort.InferenceSession.prototype.run = async function (...a) {
    const t = now();
    const r = await origRun.apply(this, a);
    if (!benching) rec('infer.run', now() - t);
    return r;
  };

  // Independent rAF ticker: frame gaps show main-thread blocking however it's caused.
  let lastFrame = 0;
  (function tick(t) {
    if (running && lastFrame) { rec('frame', t - lastFrame); rec('stepsPerFrame', stepsThisFrame); }
    lastFrame = t; stepsThisFrame = 0;
    requestAnimationFrame(tick);
  })(0);

  // ------------------------------------------------------------------------------------------
  //  Environment
  // ------------------------------------------------------------------------------------------
  async function environment() {
    const env = {
      ua: navigator.userAgent, platform: navigator.userAgentData?.platform ?? navigator.platform,
      cores: navigator.hardwareConcurrency, memoryGB: navigator.deviceMemory,
      dpr: devicePixelRatio, screen: `${screen.width}x${screen.height}`, viewport: `${innerWidth}x${innerHeight}`,
      crossOriginIsolated: window.crossOriginIsolated, sharedArrayBuffer: typeof SharedArrayBuffer !== 'undefined',
      swControlled: !!navigator.serviceWorker?.controller,
      ortVersion: ort.env.versions?.web, wasmThreads: ort.env.wasm.numThreads, wasmSimd: ort.env.wasm.simd,
      requestedEPs: EXECUTION_PROVIDERS, webgpuApi: !!navigator.gpu,
    };
    // Whether the game's own session really got the WebGPU backend, or fell back to WASM quietly.
    try { env.webgpuBackendUp = !!ort.env.webgpu?.adapter; } catch (e) { env.webgpuBackendUp = 'unknown'; }
    if (navigator.gpu) {
      try {
        const ad = await navigator.gpu.requestAdapter();
        if (ad) {
          const info = ad.info || (ad.requestAdapterInfo ? await ad.requestAdapterInfo() : {});
          env.webgpuAdapter = { vendor: info.vendor, architecture: info.architecture, device: info.device, description: info.description };
          env.webgpuF16 = ad.features.has('shader-f16');
        } else env.webgpuAdapter = null;
      } catch (e) { env.webgpuAdapter = 'error: ' + e.message; }
    }
    try {
      const gl = document.createElement('canvas').getContext('webgl2');
      const ext = gl.getExtension('WEBGL_debug_renderer_info');
      env.glRenderer = ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
      env.glVendor = ext ? gl.getParameter(ext.UNMASKED_VENDOR_WEBGL) : gl.getParameter(gl.VENDOR);
      const hp = gl.getShaderPrecisionFormat(gl.FRAGMENT_SHADER, gl.HIGH_FLOAT);
      env.fragHighpFloat = hp.precision;
      gl.getExtension('WEBGL_lose_context')?.loseContext();
    } catch (e) { env.glRenderer = 'error: ' + e.message; }
    let d = Infinity, t = now();
    for (let i = 0; i < 1e5 && d === Infinity; i++) { const u = now(); if (u > t) d = u - t; }
    env.timerResolutionMs = d;
    return env;
  }

  // ------------------------------------------------------------------------------------------
  //  Bench: the policy alone, per execution provider and crop size
  // ------------------------------------------------------------------------------------------
  const ready = () => new Promise(res => { (function w() { session && lastCrop ? res() : setTimeout(w, 200); })(); });

  // The live crop's centre SxS, so the bench sees a real soliton rather than an empty window.
  function feeds(S) {
    const net = CFG.netSize, off = (net - S) >> 1, SS = S * S;
    const data = new Float32Array(CFG.K * SS);
    for (let k = 0; k < CFG.K; k++)
      for (let y = 0; y < S; y++)
        for (let x = 0; x < S; x++) data[k * SS + y * S + x] = lastCrop[((y + off) * net + x + off) * 4 + k];
    return { state: new ort.Tensor('float32', data, [1, CFG.K, S, S]), context: new ort.Tensor('float32', buildContext(), [1, 4]) };
  }

  // A MessageChannel ping-pong runs alongside inference: its longest gap says how long the main
  // thread was held. ~runMs means inference blocks the page (pipelining can't overlap it);
  // a few ms means it runs off-thread.
  async function blocking(s, f) {
    const ch = new MessageChannel();
    let last = now(), maxGap = 0, ticks = 0, on = true;
    ch.port1.onmessage = () => { const t = now(); maxGap = Math.max(maxGap, t - last); last = t; ticks++; if (on) ch.port2.postMessage(0); };
    ch.port2.postMessage(0);
    const RUNS = 5, t0 = now();
    for (let i = 0; i < RUNS; i++) await s.run(f);
    const runMs = (now() - t0) / RUNS;
    on = false; ch.port1.close();
    return { runMs: Math.round(runMs * 100) / 100, maxMainThreadGapMs: Math.round(maxGap * 100) / 100, ticksPerRun: Math.round(ticks / RUNS) };
  }

  async function bench() {
    await ready();
    const wasRunning = running;
    setRunning(false);
    while (busy) await new Promise(r => setTimeout(r, 20));   // let an in-flight game step finish
    benching = true;
    const out = { when: new Date().toISOString(), threads: ort.env.wasm.numThreads, crossOriginIsolated: window.crossOriginIsolated, results: [] };
    const eps = navigator.gpu ? ['webgpu', 'wasm'] : ['wasm'];
    try {
      for (const ep of eps) {
        const r = { ep };
        status(`bench: ${ep} …`);
        try {
          let t = now();
          const s = await ort.InferenceSession.create(CFG.modelUrl, { executionProviders: [ep] });
          r.createMs = Math.round(now() - t);
          for (const S of [96, 64, 48]) {
            if (S > CFG.netSize) continue;
            status(`bench: ${ep} ${S}² …`);
            const f = feeds(S);
            t = now(); await s.run(f); r[`firstRun${S}`] = Math.round(now() - t);   // JIT / shader compile
            await s.run(f); await s.run(f);
            const xs = [], RUNS = S === 96 ? 20 : 10;
            for (let i = 0; i < RUNS; i++) { t = now(); await s.run(f); xs.push(now() - t); }
            r[`run${S}`] = summarize(xs);
          }
          status(`bench: ${ep} main-thread test …`);
          r.blocking = await blocking(s, feeds(Math.min(96, CFG.netSize)));
          await s.release?.();
        } catch (e) { r.error = String(e?.message || e); }
        out.results.push(r);
      }
    } finally {
      benching = false;
      if (wasRunning) setRunning(true);
    }
    const all = load(localStorage, BENCH_KEY, []);
    all.push(out);
    store(localStorage, BENCH_KEY, all);
    status(`bench done (${all.length} stored)`);
    return out;
  }

  // ------------------------------------------------------------------------------------------
  //  Thread sweep: one bench per page load, since numThreads can't change after WASM init
  // ------------------------------------------------------------------------------------------
  function gotoThreads(n) {
    const u = new URL(location.href);
    u.searchParams.set('threads', n);
    location.href = u.toString();
  }
  function startSweep() {
    const hc = navigator.hardwareConcurrency || 4;
    const list = [...new Set([1, 2, 4, Math.min(hc, 8)])].filter(n => n <= Math.max(hc, 1));
    store(sessionStorage, SWEEP_KEY, list);
    gotoThreads(list[0]);
  }
  async function continueSweep() {
    const list = load(sessionStorage, SWEEP_KEY, null);
    if (!list || !list.length) return;
    if (+params.get('threads') !== list[0]) { gotoThreads(list[0]); return; }
    await bench();
    list.shift();
    if (list.length) { store(sessionStorage, SWEEP_KEY, list); gotoThreads(list[0]); }
    else { store(sessionStorage, SWEEP_KEY); status('thread sweep done -- Copy report'); }
  }

  // ------------------------------------------------------------------------------------------
  //  Report + panel
  // ------------------------------------------------------------------------------------------
  // Gathered once up front: Copy has to stay synchronous, since iOS Safari only allows clipboard
  // writes inside the tap's own task.
  let envCache = null;
  ready().then(environment).then(e => { envCache = e; });

  function report() {
    const live = {};
    for (const k of Object.keys(series).sort()) live[k] = summarize(series[k]);
    return {
      reportVersion: 1, when: new Date().toISOString(), url: location.href,
      env: envCache,
      game: { actorMode, sps, measSps: Math.round(measSps), level, stride: inferenceStride, sometimesWindow,
              frameBudgetMs: FRAME_BUDGET_MS, netSize: CFG.netSize, board: `${W}x${H}` },
      live, bench: load(localStorage, BENCH_KEY, []),
    };
  }

  const panel = document.createElement('div');
  panel.style.cssText = 'position:fixed;left:8px;right:8px;bottom:8px;max-width:560px;z-index:1000;' +
    'background:var(--raised);color:var(--ink);border:1px solid var(--line2);border-radius:10px;' +
    'padding:8px 10px;font:11px/1.45 "IBM Plex Mono",ui-monospace,monospace;box-shadow:0 4px 18px rgba(0,0,0,.25)';
  panel.innerHTML =
    '<div id="dbg-live" style="white-space:pre-wrap">debug: waiting for model…</div>' +
    '<div id="dbg-status" style="color:var(--muted)"></div>' +
    '<div style="display:flex;flex-wrap:wrap;gap:6px;margin-top:6px">' +
    ['bench:Bench', 'sweep:Sweep threads', 'copy:Copy report', 'reset:Reset live', 'clear:Clear benches', 'hide:Hide']
      .map(s => { const [id, label] = s.split(':'); return `<button id="dbg-${id}" style="font:inherit;padding:4px 8px;border-radius:6px;border:1px solid var(--line2);background:var(--bg);color:var(--ink)">${label}</button>`; }).join('') +
    '</div><textarea id="dbg-out" readonly hidden style="width:100%;height:140px;margin-top:6px;font:10px ui-monospace,monospace;box-sizing:border-box"></textarea>';
  document.body.appendChild(panel);
  const el = id => panel.querySelector('#dbg-' + id);
  function status(msg) { el('status').textContent = msg; }

  el('bench').onclick = () => bench();
  el('sweep').onclick = () => startSweep();
  el('reset').onclick = () => { for (const k in series) delete series[k]; status('live stats reset'); };
  el('clear').onclick = () => { store(localStorage, BENCH_KEY); status('stored benches cleared'); };
  el('hide').onclick = () => {
    const hid = el('live').hidden = !el('live').hidden;
    el('status').hidden = hid; el('out').hidden = true;
    el('hide').textContent = hid ? 'Show' : 'Hide';
  };
  el('copy').onclick = async () => {
    const txt = JSON.stringify(report());
    const out = el('out');
    out.value = txt; out.hidden = false; out.select();
    try { await navigator.clipboard.writeText(txt); status(`report copied (${txt.length} chars)`); }
    catch (e) { status('clipboard blocked -- select the text below and copy it'); }
  };

  setInterval(() => {
    if (el('live').hidden) return;
    const p = name => { const s = summarize(series[name]); return s ? `${s.p50}/${s.p90}` : '–'; };
    let ep = EXECUTION_PROVIDERS.join('>');
    try { ep += ort.env.webgpu?.adapter ? ' (webgpu up)' : ''; } catch (e) {}
    el('live').textContent =
      `${ep}  thr ${ort.env.wasm.numThreads}  COI ${window.crossOriginIsolated ? 'yes' : 'NO'}  ${actorMode}  L${level}\n` +
      `p50/p90 ms  infer ${p('infer')}  run ${p('infer.run')}  readback ${p('gpu.readback')}  submit ${p('gpu.submit')}\n` +
      `step carl ${p('step.carl')}  idle ${p('step.idle')}  frame ${p('frame')}  render ${p('render')}\n` +
      `rate ${Math.round(measSps)}/${sps} steps/s`;
  }, 500);

  continueSweep();
})();
