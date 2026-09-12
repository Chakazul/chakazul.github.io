#version 300 es
// Composites the ghost tiles onto board-sized textures: o0 for draw.glsl's colouring, o1 (mass
// with frightened ghosts excluded) for sim.glsl's ghost-eats-Pac-Man check -- nothing downstream
// needs to know the ghosts live in tiles. One pass over the whole board rather than a draw per
// tile: tiles can overlap once two ghosts are close, and separate draws would have the later one
// *replace* the earlier instead of summing, silently deleting a ghost.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

#define MAX_GHOSTS 9   // must match glsim.js's MAX_GHOSTS

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
        // Whoever has the most mass here owns the pixel's colour -- exact even where tiles overlap.
        if (v > best) { best = v; owner = float(i) + 1.0; }
    }

    o0 = vec4(mass, owner, 0.0, 1.0);
    o1 = vec4(danger, 0.0, 0.0, 1.0);
}
