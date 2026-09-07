#version 300 es
// Paints the ghost tiles back onto board-sized textures, which is what the rest of the demo
// actually consumes: draw.glsl colours from o0 and sim.glsl reads o1 to erase Pac-Man where a
// non-frightened ghost overlaps him. Nothing downstream needs to know the ghosts live in tiles.
//
// One pass over the board rather than one small draw per tile, because tiles can overlap in board
// space once two ghosts are close, and separate draws would have the later one *replace* the
// earlier instead of adding to it -- quietly deleting a ghost. Summing here also needs no blending,
// which on a float target would mean depending on EXT_float_blend.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

#define MAX_GHOSTS 8

uniform sampler2D uAtlas;
uniform ivec2 uSize;                 // board (W, H)
uniform int   uWin;                  // tile edge
uniform int   uCount;                // ghosts actually present
uniform ivec2 uOrigin[MAX_GHOSTS];
uniform int   uFrightened[MAX_GHOSTS];   // per-ghost: excluded from o1 (the "dangerous" total)

layout(location = 0) out vec4 o0;   // (summed mass, 1-based index of the ghost contributing most, -, -)
layout(location = 1) out vec4 o1;   // (mass from non-frightened ghosts only, -, -, -)

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);

    float mass = 0.0, danger = 0.0, best = 0.0, owner = 0.0;
    for (int i = 0; i < uCount; i++) {
        ivec2 d = p - uOrigin[i];
        // Nearest image: a tile whose origin is near a board edge covers cells on the far side.
        if (d.x >  uSize.x / 2) d.x -= uSize.x;
        if (d.x < -uSize.x / 2) d.x += uSize.x;
        if (d.y >  uSize.y / 2) d.y -= uSize.y;
        if (d.y < -uSize.y / 2) d.y += uSize.y;
        if (d.x < 0 || d.x >= uWin || d.y < 0 || d.y >= uWin) continue;

        float v = texelFetch(uAtlas, ivec2(i * uWin + d.x, d.y), 0).r;
        mass += v;
        if (uFrightened[i] == 0) danger += v;
        // Whoever has the most mass here owns the pixel's colour. Carried in the texture rather
        // than re-derived in draw.glsl from ghost positions, so it is exact even where tiles
        // overlap, and costs nothing extra -- this loop already knows the answer.
        if (v > best) { best = v; owner = float(i) + 1.0; }
    }

    o0 = vec4(mass, owner, 0.0, 1.0);
    o1 = vec4(danger, 0.0, 0.0, 1.0);
}
