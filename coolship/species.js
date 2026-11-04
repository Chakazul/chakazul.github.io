"use strict";
// ============================================================================================
//  Species data for Coolship: the picker's patterns and their rules.
//
//  Source: Bert Chan's Lenia "aquarium" collection, Lenia/Python/found/233.json, the entries from
//  5N7KKM to GDNQYX, generated into aquarium233.js by tools/gen_aquarium.py. 5N7KKM is the
//  original Aquarium, found by globally randomising parameters; the other codes were evolved from
//  it. Forms of one code share its rule and differ in their cells (VT049W gyrotachus, whose kernel 0
//  differs, is a different species and left out). 5N7KKM limus is the default ship. VT049W rectus is the steady glider used for
//  the Phase 1 checks; VT049W gyrans R is the
//  original self-replicating Aquarium ("R" = replicating), and the other VT049W forms were
//  stabilised from it.
//
//  The rule follows LeniaF.py, the program that reads 233.json (run as LeniaF.py -c 3 -k 3):
//    kernel: shell rings b over r*R, polynomial core (4f(1-f))^4 (kn = 1), each normalised to 1;
//    growth: G = 2 max(0, 1 - (u-m)^2 / (9 s^2))^4 - 1 (gn = 1);
//    h: free_h, h' = round(h / sum h over the same destination, 3), doubled with T (T: 1 -> 2);
//    update: A_c += (1/T) * sum_{k: dst=c} h'_k G_k, clipped to [0, 1].
//  c = [source channel, destination channel]. glsim.js rebuilds the kernels at the run radius
//  and the cells are resampled nearest-neighbour by R / (the entry's native R). Checks (numpy port of LeniaF,
//  README "Checks for Phase 1b"): at native R the forms match their names (lithos still, tachus
//  fast, only gyrans R replicates); at R = 18 rectus glides and gyrans R replicates to 3-4
//  creatures. Below R = 16, gyrans R stops replicating and lithos moves.
// ============================================================================================

// Patterns offered by the picker (app.js), in file order; DEFAULT_PATTERN is selected first.
//   rule: a code in AQUARIUM233.rules; R: the entry's native radius; cells: RLE per channel.
//   logic: 'leniaF' (the rule above), or 'legacy', the engine's first version: exponential bump
//   kernel and Gaussian growth (kn = gn = 2), original h, optionally normalised per channel by
//   sum h (norm), at its own dt. The legacy patterns are hidden: kept for reference, reachable
//   only with ?form=.
const DEFAULT_PATTERN = AQUARIUM233.patterns[0].id;   // 5N7KKM limus, the original Aquarium
const PATTERNS = (function () {
  const list = AQUARIUM233.patterns.map(p => ({ ...p, logic: 'leniaF' }));
  const rectus = list.find(p => p.id === 'VT049W-rectus');
  const legacy = { rule: 'VT049W', R: rectus.R, cells: rectus.cells, logic: 'legacy', hidden: true };
  // Legacy logic, R = 12 era. Unnormalised at dt = 0.5 (Lenia/WebGL index4.html "Fission
  // cells"): divides from ~350 steps, fills the board by ~1250.
  list.push({ ...legacy, id: 'fission', name: 'Fission Aquarium (legacy)', norm: false, dt: 0.5 });
  // Normalised at dt = 0.2 with kernel 5 (G->G) mu x 0.955: swells to ~3x, then splits unevenly.
  list.push({ ...legacy, id: 'fission02', name: 'Fission Aquarium, dt 0.2 (legacy test)', norm: true, dt: 0.2,
    tweak: { k: 5, m: 0.955 } });
  return list;
})();

// The pattern's rule: { kernels, kn, gn, norm, dt }, ready for SimGL.setRule / setUpdate.
function patternRule(p) {
  const base = AQUARIUM233.rules[p.rule];
  let ks = (p.kernels || base.kernels).map(k => ({ ...k }));
  if (p.tweak) for (const key of ['m', 's', 'h']) if (p.tweak[key]) ks[p.tweak.k][key] *= p.tweak[key];
  if (p.logic === 'legacy') return { kernels: ks, kn: 2, gn: 2, norm: p.norm, dt: p.dt };
  const f = Rule.leniaF(ks, base.T);
  return { kernels: f.kernels, kn: 1, gn: 1, norm: false, dt: 1 / f.T };
}
