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
//  State goes in as a frame-major Float32Array [K*S*S] and context as Float32Array [4]; infer()
//  returns the flat argmax index into the [3,S,S] Q map. The caller turns that into a board
//  position, since only it knows where the crop was taken.
// ============================================================================================

// `cfg`: { filmUrl, modelUrl, coreUrl: {slice, gather}, netSize, K, epParam, modelParam }.
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
      filmCache = { key, feeds: await filmSession.run({ context: new ort.Tensor('float32', ctx, [1, 4]) }) };
    }
    return filmCache.feeds;
  }

  async function feedsFor(model, state, ctx) {
    const st = new ort.Tensor('float32', state, [1, cfg.K, S, S]);
    if (model === 'orig') return { state: st, context: new ort.Tensor('float32', ctx, [1, 4]) };
    return { state: st, ...(await filmFeeds(ctx)) };
  }

  const urlFor = model => (model === 'orig' ? cfg.modelUrl : cfg.coreUrl[model]);
  // The graph that measured fastest on each backend: gather on WebGPU (but it exists only at a
  // 96 crop), slice on WASM, whose Gathers are slow on CPU. ?model= forces one for both.
  function graphFor(ep) {
    const m = cfg.modelParam || (ep === 'webgpu' ? 'gather' : 'slice');
    return m === 'gather' && S !== 96 ? 'slice' : m;
  }

  async function hasGpu() {
    return !!(typeof navigator !== 'undefined' && navigator.gpu
      && await navigator.gpu.requestAdapter().catch(() => null));
  }

  // Times each candidate backend with its graph on (state, ctx) and keeps the faster, warm.
  // `eps` restricts the candidates (default: every backend available here, minus ?ep= exclusions).
  // Cost: each candidate's session creation, a compile/JIT run, a warm run, then RUNS timed runs.
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
        await c.session.run(feeds);
        const xs = [];
        for (let i = 0; i < RUNS; i++) { t = now(); await c.session.run(feeds); xs.push(now() - t); }
        c.ms = xs.sort((a, b) => a - b)[RUNS >> 1];
      } catch (e) {
        c.error = String(e?.message || e);
        try { await c.session?.release(); } catch (e2) { /* already gone */ }
        c.session = null;
      }
    }
    const ok = cands.filter(c => c.session);
    if (!ok.length) throw new Error(cands.map(c => `${c.ep}: ${c.error}`).join('; ') || 'no backend available');
    const best = ok.reduce((a, b) => (b.ms < a.ms ? b : a));
    for (const c of ok) if (c !== best) try { await c.session.release(); } catch (e) { /* ignore */ }
    session = best.session; model = best.model; ep = best.ep;
    return { ep, model, ms: best.ms, candidates: cands.map(({ session: _, ...rest }) => rest) };
  }

  async function infer(state, ctx) {
    const feeds = await feedsFor(model, state, ctx);
    const t = now();
    const q = (await session.run(feeds)).q.data;
    const runMs = now() - t;
    let best = 0, bv = q[0];
    for (let i = 1; i < q.length; i++) if (q[i] > bv) { bv = q[i]; best = i; }
    return { best, runMs };
  }

  async function release() {
    try { await session?.release(); } catch (e) { /* ignore */ }
    session = null;
  }

  return { pick, infer, feedsFor, release, get model() { return model; }, get ep() { return ep; } };
}
