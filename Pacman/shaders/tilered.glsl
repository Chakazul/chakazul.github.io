#version 300 es
// Stage 1 of the per-ghost CoM: sums one block of one tile. Plain (not circular) sums -- a tile
// has hard edges and the soliton is kept off them, so there's no wrap to average around. Output
// is 16x16 per tile, laid side by side so one draw covers every ghost (tile i at x in
// [i*16, (i+1)*16)).
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uAtlas;
uniform int uWin;     // tile edge
uniform int uBlock;   // block edge = ceil(uWin / 16)

out vec4 fragColor;   // (mass, sum v*x, sum v*y, -), all tile-local

void main() {
    ivec2 b = ivec2(gl_FragCoord.xy);
    int i = b.x / 16;                  // tile
    int bx = b.x - i * 16, by = b.y;   // block within the tile
    int x0 = bx * uBlock, y0 = by * uBlock;

    float mass = 0.0, mx = 0.0, my = 0.0;
    for (int j = 0; j < uBlock; j++) {
        int ly = y0 + j;
        if (ly >= uWin) break;
        for (int k = 0; k < uBlock; k++) {
            int lx = x0 + k;
            if (lx >= uWin) break;
            float v = texelFetch(uAtlas, ivec2(i * uWin + lx, ly), 0).r;
            if (v <= 0.0) continue;
            mass += v;
            mx += v * float(lx);
            my += v * float(ly);
        }
    }
    fragColor = vec4(mass, mx, my, 0.0);
}
