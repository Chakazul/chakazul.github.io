"use strict";
// CARL policy inference in a dedicated worker, so it runs alongside the main thread's GPU step
// instead of before it (see agentStep()'s pipelining in app.js). WASM inference otherwise holds
// the main thread for its whole run, and WebGPU for most of it (PERF.md, main-thread test).
//
// Protocol (main -> worker -> main):
//   {type:'init', cfg, threads, wasmPaths, state, ctx}  ->  {type:'ready', choice} | {type:'error', message}
//   {type:'infer', id, state, ctx}                       ->  {type:'result', id, best, runMs} | {type:'error', id, message}
// `state` is a transferred frame-major Float32Array; `best` the flat argmax index into the Q map.
importScripts('ort.webgpu.min.js', 'policy.js');

let policy = null;

self.onmessage = async ({ data: m }) => {
  if (m.type === 'init') {
    try {
      ort.env.wasm.wasmPaths = m.wasmPaths;
      // Multi-threaded WASM inside a worker needs the worker itself to be cross-origin isolated.
      ort.env.wasm.numThreads = self.crossOriginIsolated ? m.threads : 1;
      policy = createPolicyRuntime(ort, m.cfg);
      const choice = await policy.pick(m.state, m.ctx);
      choice.threads = ort.env.wasm.numThreads;
      choice.crossOriginIsolated = self.crossOriginIsolated;
      choice.gpuInWorker = !!self.navigator?.gpu;
      self.postMessage({ type: 'ready', choice });
    } catch (e) {
      self.postMessage({ type: 'error', message: String(e?.message || e) });
    }
  } else if (m.type === 'infer') {
    try {
      const { best, runMs } = await policy.infer(m.state, m.ctx);
      self.postMessage({ type: 'result', id: m.id, best, runMs });
    } catch (e) {
      self.postMessage({ type: 'error', id: m.id, message: String(e?.message || e) });
    }
  }
};
