#version 300 es
// Coarse map for finding creatures (creatures.js segments it on the CPU). Each output texel
// covers one BLOCK x BLOCK patch of the board:
//   rgb = (mass, sum m*dx, sum m*dy) over the patch and all three channels, with dx, dy the
//         cell's offset inside the block, so a creature's exact CoM can be rebuilt from its blocks;
//   a   = the potential at the patch's centre cell: the mean over kernels of K_k * A_src(k)
//         (each kernel normalised, no mask), how the rule "sees" the world. Creatures are the
//         catchment basins of its peaks (the method of ~/github/complexity, individuality.py).
// The potential is smooth on the scale of R, so one sample per block is enough.
// Second target: (mass R, mass G, mass B, 0) over the patch, for the per-channel mass stats.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uState;
uniform sampler2D uTaps;    // as in sim.glsl
uniform ivec2 uSize;        // board (W, H)
uniform int   uTapCount;
uniform mat3x4 uSrc[4];
uniform float uInvK;        // 1 / kernel count

layout(location = 0) out vec4 fragColor;
layout(location = 1) out vec4 chanColor;

const int BLOCK = 4;    // must match glsim.js

void main() {
    ivec2 b = ivec2(gl_FragCoord.xy) * BLOCK;
    float m = 0.0, mx = 0.0, my = 0.0;
    vec3 mc = vec3(0.0);
    for (int j = 0; j < BLOCK; j++) {
        for (int i = 0; i < BLOCK; i++) {
            ivec2 p = b + ivec2(i, j);
            if (p.x >= uSize.x || p.y >= uSize.y) continue;
            vec3 v = texelFetch(uState, p, 0).rgb;
            float s = v.r + v.g + v.b;
            m += s;  mx += s * float(i);  my += s * float(j);  mc += v;
        }
    }

    ivec2 c = min(b + BLOCK / 2, uSize - 1);
    vec4 U = vec4(0.0);
    for (int t = 0; t < uTapCount; t++) {
        ivec2 off = ivec2(texelFetch(uTaps, ivec2(t, 0), 0).xy);
        vec3 a = texelFetch(uState, (c + off + uSize) % uSize, 0).rgb;
        U += texelFetch(uTaps, ivec2(t, 1), 0) * (uSrc[0] * a);
        U += texelFetch(uTaps, ivec2(t, 2), 0) * (uSrc[1] * a);
        U += texelFetch(uTaps, ivec2(t, 3), 0) * (uSrc[2] * a);
        U += texelFetch(uTaps, ivec2(t, 4), 0) * (uSrc[3] * a);
    }
    fragColor = vec4(m, mx, my, dot(U, vec4(1.0)) * uInvK);
    chanColor = vec4(mc, 0.0);
}
