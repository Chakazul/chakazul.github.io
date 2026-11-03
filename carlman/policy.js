"use strict";
// ============================================================================================
//  CARL policy runtime -- loaded by the page (main-thread fallback, ?debug=1's bench) and by
//  policy-worker.js (the default), so both run the same code.
//
//  The policy is the export split by tools/split_model.py into a FiLM model (context -> 28 γ/β
//  tensors) and a conv core, plus the export itself ('orig'). Which core and which backend run
//  is decided per device at load by timing (pick()): neither backend wins everywhere -- see
//  PERF.md, "Device measurements".
//
//  Agents other than the paper model (see AGENTS in app.js) ship only the export itself, so they
//  always run 'orig'.
//
//  State goes in as a frame-major Float32Array [K*S*S] and context as Float32Array [ctxDim] (4
//  with the time input, 3 without); infer() returns the flat argmax index into the [3,S,S] Q map.
//  The caller turns that into a board position, since only it knows where the crop was taken.
// ============================================================================================

// `cfg`: { filmUrl, modelUrl, coreUrl: {slice, gather} | null, netSize, K, ctxDim, epParam, modelParam }.
// Output channel order is add / remove / no-op (CFG.actionValue in app.js).
const NOOP_CHANNEL = 2;

function createPolicyRuntime(ort, cfg) {
  const S = cfg.netSize, now = () => performance.now();

  // FiLM γ/β for a context. The context only changes on a steer or a cost-slider move, so this
  // runs on those, not per inference. CPU (WASM) on purpose: 112 tiny ops, whose outputs feed the
  // core as ordinary CPU tensors whatever backend the core runs on.
  let filmSession = null, filmCache = { key: null, feeds: null };
  async function filmFeeds(ctx) {
    const key = Array.prototype.join.call(ctx, ',');
    if (key !== filmCache.key) {
      if (!filmSession) filmSession = await ort.InferenceSession.create(cfg.filmUrl, { executionProviders: ['wasm'] });
      filmCache = { key, feeds: await filmSession.run({ context: new ort.Tensor('float32', ctx, [1, cfg.ctxDim]) }) };
    }
    return filmCache.feeds;
  }

  async function feedsFor(model, state, ctx) {
    const st = new ort.Tensor('float32', state, [1, cfg.K, S, S]);
    if (model === 'orig') return { state: st, context: new ort.Tensor('float32', ctx, [1, cfg.ctxDim]) };
    return { state: st, ...(await filmFeeds(ctx)) };
  }

  const urlFor = model => (model === 'orig' ? cfg.modelUrl : cfg.coreUrl[model]);
  // The graph that measured fastest on each backend: gather on WebGPU (but it exists only at a
  // 96 crop), slice on WASM, whose Gathers are slow on CPU. ?model= forces one for both.
  function graphFor(ep) {
    if (!cfg.coreUrl) return 'orig';
    const m = cfg.modelParam || (ep === 'webgpu' ? 'gather' : 'slice');
    return m === 'gather' && S !== 96 ? 'slice' : m;
  }

  async function hasGpu() {
    return !!(typeof navigator !== 'undefined' && navigator.gpu
      && await navigator.gpu.requestAdapter().catch(() => null));
  }

  // Times each candidate backend with its graph on (state, ctx) and keeps the faster, warm.
  // `eps` restricts the candidates (default: every backend available here, minus ?ep= exclusions).
  // Cost: each candidate's session creation, a compile/JIT run, then timed windows of RUNS runs
  // until the speed settles (at most ~2 s per candidate).
  let session = null, model = null, ep = null;
  async function pick(state, ctx, eps) {
    const RUNS = 5, gpu = cfg.epParam !== 'wasm' && await hasGpu();
    if (!eps) {
      eps = [];
      if (gpu) eps.push('webgpu');
      if (cfg.epParam !== 'webgpu' || !gpu) eps.push('wasm');
    } else if (!gpu) eps = eps.filter(e => e !== 'webgpu');
    const cands = eps.map(e => ({ ep: e, model: graphFor(e) }));
    for (const c of cands) {
      try {
        let t = now();
        c.session = await ort.InferenceSession.create(urlFor(c.model), { executionProviders: [c.ep] });
        c.createMs = Math.round(now() - t);
        const feeds = await feedsFor(c.model, state, ctx);
        t = now(); await c.session.run(feeds); c.firstRunMs = Math.round(now() - t);
        // Warm until the speed settles, not for a fixed count. On a cold load (first visit, or
        // any incognito session before the browser has cached compiled code) WASM starts on V8's
        // baseline tier and speeds up 2-3x as hot functions are recompiled optimized in the
        // background, which can take longer than a couple of runs on a phone. Timed that early,
        // WASM loses to WebGPU on devices where it is really the faster one, and the pick then
        // holds for the whole visit (~30 instead of ~55 steps/s on a phone). So: windows of RUNS
        // runs until a window's median is no more than 5% better than the last, capped by time.
        const median = xs => xs.slice().sort((a, b) => a - b)[xs.length >> 1];
        const WARM_CAP_MS = 2000, t0 = now();
        let prev = Infinity, windows = 0;
        for (;;) {
          const xs = [];
          for (let i = 0; i < RUNS; i++) { t = now(); await c.session.run(feeds); xs.push(now() - t); }
          const m = median(xs); windows++;
          c.ms = Math.min(prev, m);
          if (m > prev * 0.95 || now() - t0 > WARM_CAP_MS) break;
          prev = m;
        }
        c.warmRuns = windows * RUNS; c.warmMs = Math.round(now() - t0);
      } catch (e) {
        c.error = String(e?.message || e);
        try { await c.session?.release(); } catch (e2) { /* already gone */ }
        c.session = null;
      }
    }
    const ok = cands.filter(c => c.session);
    if (!ok.length) throw new Error(cands.map(c => `${c.ep}: ${c.error}`).join('; ') || 'no backend available');
    const best = ok.reduce((a, b) => (faster(b.ep, b.ms, a.ep, a.ms) ? b : a));
    // The runner-up is kept, warm, for the re-check during play (see infer()); anything else goes.
    const runnerUp = ok.find(c => c !== best) || null;
    for (const c of ok) if (c !== best && c !== runnerUp) try { await c.session.release(); } catch (e) { /* ignore */ }
    session = best.session; model = best.model; ep = best.ep;
    alt = runnerUp && { session: runnerUp.session, model: runnerUp.model, ep: runnerUp.ep, xs: [] };
    live = []; recheck = { started: now(), lastSample: 0 };
    return { ep, model, ms: best.ms, candidates: cands.map(({ session: _, ...rest }) => rest) };
  }

  // Is backend `a` at `aMs` the better pick than `b` at `bMs`? WebGPU has to be clearly faster to
  // win. WASM wins ties and near-ties: phone timings are noisy, and with the worker pipelining
  // inference against the ~18 ms phone sim, a few ms saved buys nothing. On every device in
  // PERF.md's table this picks the measured-fastest backend (iPhone and M2 WebGPU, the rest WASM).
  const WEBGPU_MARGIN = 0.8;
  function faster(a, aMs, b, bMs) {
    if (a === b) return aMs < bMs;
    return a === 'webgpu' ? aMs < bMs * WEBGPU_MARGIN : !(bMs < aMs * WEBGPU_MARGIN);
  }

  // Re-check during play. The load-time pick runs while the page is still busy loading and the
  // browser may not have optimized its code yet, so it can be wrong in a way no warm-up fully
  // rules out. Once the game is running, every inference's time on the chosen backend is
  // recorded, and about once a second one extra inference runs on the runner-up, under the same
  // load. After RECHECK_SAMPLES of those the medians are compared with the same margin, the
  // backend is switched live if the runner-up wins, and the loser is released. One extra
  // inference a second costs a step or so for those ~10 s; it happens once per page load.
  const RECHECK_DELAY_MS = 3000, RECHECK_EVERY_MS = 1000, RECHECK_SAMPLES = 8;
  let alt = null, live = [], recheck = null, switched = null;
  async function recheckStep(state, ctx) {
    const t = now();
    if (t - recheck.started < RECHECK_DELAY_MS || t - recheck.lastSample < RECHECK_EVERY_MS) return;
    recheck.lastSample = t;
    const feeds = await feedsFor(alt.model, state, ctx);
    const t1 = now(); await alt.session.run(feeds); alt.xs.push(now() - t1);
    if (alt.xs.length < RECHECK_SAMPLES) return;
    const median = xs => xs.slice().sort((a, b) => a - b)[xs.length >> 1];
    const curMs = median(live.slice(-30)), altMs = median(alt.xs);
    let loser = alt;
    if (faster(alt.ep, altMs, ep, curMs)) {
      switched = { from: ep, to: alt.ep, fromMs: Math.round(curMs * 10) / 10, toMs: Math.round(altMs * 10) / 10 };
      loser = { session };
      session = alt.session; model = alt.model; ep = alt.ep;
    }
    alt = null; recheck = null; live = [];
    try { await loser.session.release(); } catch (e) { /* ignore */ }
  }

  // `thr`: the action threshold. An add/remove only goes ahead if its Q beats the best no-op by at
  // least this much; otherwise the step becomes a no-op at the same cell. 0 is the plain argmax.
  // Damps reactions to a soliton's natural wobble (see the steering_small report in rlca).
  async function infer(state, ctx, thr = 0) {
    const feeds = await feedsFor(model, state, ctx);
    const t = now();
    const q = (await session.run(feeds)).q.data;
    const runMs = now() - t;
    let best = 0, bv = q[0];
    for (let i = 1; i < q.length; i++) if (q[i] > bv) { bv = q[i]; best = i; }
    const SS = S * S, noop = NOOP_CHANNEL;
    let suppressed = false;
    if (thr > 0 && Math.floor(best / SS) !== noop) {
      let bn = -Infinity;
      for (let i = noop * SS; i < (noop + 1) * SS; i++) if (q[i] > bn) bn = q[i];
      if (bv - bn < thr) { best = noop * SS + (best % SS); suppressed = true; }
    }
    if (recheck) {
      live.push(runMs);
      try { await recheckStep(state, ctx); }
      catch (e) { try { await alt?.session.release(); } catch (e2) { /* ignore */ } alt = null; recheck = null; }
    }
    // ep/model ride along so the caller sees a live switch; `switched` reports it (once known).
    return { best, runMs, suppressed, ep, model, switched };
  }

  async function release() {
    try { await session?.release(); } catch (e) { /* ignore */ }
    try { await alt?.session.release(); } catch (e) { /* ignore */ }
    session = null; alt = null; recheck = null;
  }

  return { pick, infer, feedsFor, release, get model() { return model; }, get ep() { return ep; } };
}
