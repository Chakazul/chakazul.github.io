"use strict";
// ============================================================================================
//  Clone: a morph held to be pasted into the world (the game's clone mode, the lab's Add).
//  A held clone is { name, M (board-frame 2x2, row-major), src (the form at the run's R), form
//  (src turned and flipped by M: { w, h, cells }), img (its picture: rgb = channels, alpha = their
//  largest) }. The page decides where it goes; this turns, flips, draws and pastes it.
// ============================================================================================

const Clone = (function () {
  // Hold pattern p (cells at its native p.R) at the run's radius R
  function hold(p, R, name = p.id) {
    const H = { name, M: [1, 0, 0, 1], src: Rule.formCells(p.cells, p.R, R) };
    build(H);
    return H;
  }
  // Another pattern into the same hold, its turn and flip kept
  function swap(H, p, R, name = p.id) {
    H.name = name; H.src = Rule.formCells(p.cells, p.R, R);
    build(H);
  }
  // The held form turned and flipped by M (board frame), resampled bilinearly, and its picture
  function build(H) {
    const s = H.src, [a, b, c, d] = H.M, det = a * d - b * c;
    const inv = [d / det, -b / det, -c / det, a / det];   // M^-1, row-major
    const S = Math.ceil(Math.hypot(s.w, s.h)) + 2, m = (S - 1) / 2, sx0 = (s.w - 1) / 2, sy0 = (s.h - 1) / 2;
    const cells = [0, 1, 2].map(() => Array.from({ length: S }, () => new Array(S).fill(0)));
    const at = (k, y, x) => (y >= 0 && y < s.h && x >= 0 && x < s.w ? s.cells[k][y][x] : 0);
    for (let y = 0; y < S; y++)
      for (let x = 0; x < S; x++) {
        const dx = x - m, dy = y - m, sx = inv[0] * dx + inv[1] * dy + sx0, sy = inv[2] * dx + inv[3] * dy + sy0;
        const x0 = Math.floor(sx), y0 = Math.floor(sy), fx = sx - x0, fy = sy - y0;
        for (let k = 0; k < 3; k++)
          cells[k][y][x] = (at(k, y0, x0) * (1 - fx) + at(k, y0, x0 + 1) * fx) * (1 - fy) + (at(k, y0 + 1, x0) * (1 - fx) + at(k, y0 + 1, x0 + 1) * fx) * fy;
      }
    H.form = { w: S, h: S, cells };
    const cv = document.createElement('canvas');
    cv.width = cv.height = S;
    const g = cv.getContext('2d'), im = g.createImageData(S, S);
    for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
      const v = [0, 1, 2].map(k => cells[k][y][x]), i = (y * S + x) * 4;
      im.data.set([...v.map(q => Math.round(255 * Math.min(1, q))), Math.round(255 * Math.min(1, Math.max(...v)))], i);
    }
    g.putImageData(im, 0, 0);
    H.img = cv;
  }
  const mul = (P, Q) => [P[0] * Q[0] + P[1] * Q[2], P[0] * Q[1] + P[1] * Q[3], P[2] * Q[0] + P[3] * Q[2], P[2] * Q[1] + P[3] * Q[3]];
  // 45 degrees clockwise as seen (a rotation looks the same in any view)
  function rotate(H) {
    const r = Math.PI / 4, c = Math.cos(r), s = Math.sin(r);
    H.M = mul([c, -s, s, c], H.M); build(H);
  }
  // Left to right as seen in a view turned by phi (board = rot(view, phi); the world view: 0)
  function flip(H, phi = 0) {
    const c = Math.cos(phi), s = Math.sin(phi);
    H.M = mul(mul(mul([c, -s, s, c], [-1, 0, 0, 1]), [c, s, -s, c]), H.M); build(H);
  }
  // The ghost centred at (x, y) on a view turned by phi, cellPx px per cell
  function draw(ctx, H, x, y, phi, cellPx) {
    ctx.save();
    ctx.translate(x, y); ctx.rotate(-phi); ctx.scale(cellPx, cellPx);
    ctx.imageSmoothingEnabled = false; ctx.globalAlpha = 0.75;
    ctx.drawImage(H.img, -H.form.w / 2, -H.form.h / 2);
    ctx.restore();
  }
  // Add the held form into the world, centred on board cell at = [col, row]
  function paste(H, at) {
    SimGL.add(H.form, Math.round(at[1] - H.form.h / 2), Math.round(at[0] - H.form.w / 2));
  }
  return { hold, swap, build, rotate, flip, draw, paste };
})();
