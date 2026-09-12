#version 300 es
// Stage 2 of the per-ghost CoM: folds one tile's 16x16 partials into one texel, one per ghost in
// a single output row, so the whole pack comes back in the same one-stall readback.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uRed;

out vec4 fragColor;   // (local row, local col, mass, valid)

void main() {
    int i = int(gl_FragCoord.x);      // ghost index

    float mass = 0.0, mx = 0.0, my = 0.0;
    for (int j = 0; j < 16; j++) {
        for (int k = 0; k < 16; k++) {
            vec4 a = texelFetch(uRed, ivec2(i * 16 + k, j), 0);
            mass += a.x; mx += a.y; my += a.z;
        }
    }

    // Nothing left in this tile: the ghost dissolved -- JS treats that as a death.
    if (mass < 1e-6) { fragColor = vec4(0.0, 0.0, 0.0, 0.0); return; }
    fragColor = vec4(my / mass, mx / mass, mass, 1.0);
}
