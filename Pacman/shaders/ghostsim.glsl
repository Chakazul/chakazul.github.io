#version 300 es
// One Lenia step per ghost, each confined to its own 96x96 tile of one atlas texture (ghost i at
// x in [i*uWin, (i+1)*uWin)); the convolution never reads past its own tile's edge. Two ghosts
// sharing a field would be two solitons that annihilate or blow up on contact -- tiles let them
// pass through each other instead, and make cost scale with ghost count rather than board size.
// uOrigin anchors each tile to real board coordinates for wall lookups; uShift re-centres the
// tile by whole cells each step (a fractional shift would resample the soliton and smear it away).
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

#define MAX_GHOSTS 9   // must match glsim.js's MAX_GHOSTS

uniform sampler2D uAtlas;    // source atlas, R = cell value
uniform sampler2D uWall;     // board wall mask
#ifdef KERNEL_TEX
uniform sampler2D uKernel;   // (2R+1)^2 normalized growth kernel
#else
// Per-row nonzero runs, same layout as sim.glsl's -- see makeKernel() in glsim.js.
layout(std140) uniform Kernel {
    ivec4 uRows[37];
    vec4  uK[370];
};
#endif
uniform sampler2D uPacman;   // Pac-Man's board-space state, post-step -- read for a tile only when
                             // uFrightened[i], to erase that ghost on overlap (per-tile, so one
                             // frightened ghost doesn't affect its packmates)
uniform int   uFrightened[MAX_GHOSTS];
uniform float uEatThreshold;
uniform ivec2 uSize;         // board (W, H)
uniform int   uWin;          // tile edge
uniform int   uR;            // kernel radius
uniform float uMu;
uniform float uSigma;
uniform float uDt[MAX_GHOSTS];       // per-tile: a frightened ghost runs at a fraction of this
uniform ivec2 uOrigin[MAX_GHOSTS];   // board position of each tile's (0,0), after this step's shift
uniform ivec2 uShift[MAX_GHOSTS];    // whole cells the tile contents move by this step

out vec4 fragColor;

void main() {
    ivec2 a = ivec2(gl_FragCoord.xy);
    int i = a.x / uWin;                       // which ghost's tile this fragment is in
    ivec2 l = ivec2(a.x - i * uWin, a.y);     // destination cell, tile-local
    ivec2 sl = l + uShift[i];                 // the source cell it is fed by
    int base = i * uWin;

    // Taps outside this tile read as empty rather than wrapping into the neighbouring ghost's tile.
    float conv = 0.0;
#ifdef KERNEL_TEX
    for (int dy = -uR; dy <= uR; dy++) {
        int y = sl.y + dy;
        if (y < 0 || y >= uWin) continue;
        for (int dx = -uR; dx <= uR; dx++) {
            float w = texelFetch(uKernel, ivec2(dx + uR, dy + uR), 0).r;
            if (w == 0.0) continue;
            int x = sl.x + dx;
            if (x < 0 || x >= uWin) continue;
            conv += w * texelFetch(uAtlas, ivec2(base + x, y), 0).r;
        }
    }
#else
    for (int r = 0; r <= 2 * uR; r++) {
        int y = sl.y + r - uR;
        if (y < 0 || y >= uWin) continue;
        ivec4 row = uRows[r];
        for (int c = 0; c < row.y; c++) {
            vec4 w = uK[row.z + c];
            int x = sl.x + row.x + 4 * c;
            // Out-of-tile taps contribute 0, the same as being skipped, in tap order.
            conv += w.x * (x     >= 0 && x     < uWin ? texelFetch(uAtlas, ivec2(base + x,     y), 0).r : 0.0);
            conv += w.y * (x + 1 >= 0 && x + 1 < uWin ? texelFetch(uAtlas, ivec2(base + x + 1, y), 0).r : 0.0);
            conv += w.z * (x + 2 >= 0 && x + 2 < uWin ? texelFetch(uAtlas, ivec2(base + x + 2, y), 0).r : 0.0);
            conv += w.w * (x + 3 >= 0 && x + 3 < uWin ? texelFetch(uAtlas, ivec2(base + x + 3, y), 0).r : 0.0);
        }
    }
#endif

    float s = (sl.x >= 0 && sl.x < uWin && sl.y >= 0 && sl.y < uWin)
            ? texelFetch(uAtlas, ivec2(base + sl.x, sl.y), 0).r : 0.0;
    float d = uMu - conv;
    float v = clamp(s + (2.0 * exp(-(d * d) / (2.0 * uSigma * uSigma)) - 1.0) * uDt[i], 0.0, 1.0);

    // Wall collision at the board cell this tile cell currently maps to.
    ivec2 b = uOrigin[i] + l;
    b.x -= uSize.x * int(floor(float(b.x) / float(uSize.x)));
    b.y -= uSize.y * int(floor(float(b.y) / float(uSize.y)));
    if (texelFetch(uWall, b, 0).r > 0.5) v = 0.0;
    else if (uFrightened[i] != 0 && texelFetch(uPacman, b, 0).r > uEatThreshold) v = 0.0;

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
