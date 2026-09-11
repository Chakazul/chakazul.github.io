#version 300 es
// One Lenia step for every ghost, each inside its own private 96x96 world.
//
// The ghosts do not share a field. They live in tiles of one atlas texture -- ghost i occupies
// x in [i*uWin, (i+1)*uWin) -- and the convolution below refuses to read past its own tile's edge.
// That is the whole point: two ghosts in a shared field are two Lenia solitons, and two Lenia
// solitons that meet annihilate or blow up. Tiles make them pass through each other instead, which
// is what an arcade ghost does. It also means the cost is per ghost rather than per board: a tile
// is 9216 cells against the board's 137500, so several ghosts together come to a fraction of the
// single whole-board pass this replaces.
//
// A tile is a window onto the board, not a world of its own in every respect: uOrigin[i] says
// where its top-left sits, so the wall lookup at the bottom is a real maze lookup and the ghost
// still runs the corridors. uShift[i] slides the tile's contents by whole cells to keep the
// soliton centred as it travels -- whole cells because a fractional shift would mean resampling,
// and resampling a soliton every step smears it away.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

#define MAX_GHOSTS 9   // must match glsim.js's MAX_GHOSTS

uniform sampler2D uAtlas;    // source atlas, R = cell value
uniform sampler2D uWall;     // board wall mask
uniform sampler2D uKernel;   // (2R+1)^2 normalized growth kernel
uniform sampler2D uPacman;   // Pac-Man's board-space state, post-step -- read for a tile only when
                             // uFrightened[i] is set (Pac-Man eats that ghost on overlap, the same
                             // erase-on-threshold shape sim.glsl uses the other way round for the
                             // ordinary ghosts-eat-Pac-Man case). Per-tile, not all-or-nothing, so
                             // one ghost being frightened doesn't touch its packmates.
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

    // Toroidal nowhere: taps that fall outside this tile read as empty rather than wrapping to the
    // far side or reaching into the neighbouring ghost. The soliton is ~46px wide and kept centred,
    // so it never comes near an edge and never notices the difference.
    float conv = 0.0;
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

    float s = (sl.x >= 0 && sl.x < uWin && sl.y >= 0 && sl.y < uWin)
            ? texelFetch(uAtlas, ivec2(base + sl.x, sl.y), 0).r : 0.0;
    float d = uMu - conv;
    float v = clamp(s + (2.0 * exp(-(d * d) / (2.0 * uSigma * uSigma)) - 1.0) * uDt[i], 0.0, 1.0);

    // Wall collision, applied after the step exactly as the board-space channels do it -- at the
    // board cell this tile cell currently maps to.
    ivec2 b = uOrigin[i] + l;
    b.x -= uSize.x * int(floor(float(b.x) / float(uSize.x)));
    b.y -= uSize.y * int(floor(float(b.y) / float(uSize.y)));
    if (texelFetch(uWall, b, 0).r > 0.5) v = 0.0;
    else if (uFrightened[i] != 0 && texelFetch(uPacman, b, 0).r > uEatThreshold) v = 0.0;

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
