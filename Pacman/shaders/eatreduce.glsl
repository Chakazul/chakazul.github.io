#version 300 es
// Stage 1 of "how much dot mass did Pac-Man just eat": same 16x16 block-tiling as reduce.glsl,
// summing dots-channel mass (pre-step, uDots) at cells sim.glsl's own erase check will zero this
// step (uPacman post-step vs uEatThreshold) -- the real eat condition, not a before/after mass
// delta, which would also pick up channel 2's own unrelated growth noise. Also splits out
// power-pellet mass (.g, via uPower) to trigger the frightened window.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uDots;      // channel 2 (dots + pellets), pre-step
uniform sampler2D uPacman;    // channel 1 (Pac-Man), post-step
uniform sampler2D uPower;     // 1 at a power-pellet cell, 0 at a plain dot
uniform ivec2 uSize;          // board (W, H)
uniform int   uBlock;         // block edge = ceil(max(W,H)/16), same as reduce.glsl
uniform float uThreshold;     // must match glsim.js's EAT_THRESHOLD

out vec4 fragColor;   // (total eaten this block, power-pellet eaten this block, -, -)

void main() {
    ivec2 b = ivec2(gl_FragCoord.xy);
    int x0 = b.x * uBlock, y0 = b.y * uBlock;

    float sum = 0.0, pelletSum = 0.0;
    for (int j = 0; j < uBlock; j++) {
        int y = y0 + j;
        if (y >= uSize.y) break;
        for (int i = 0; i < uBlock; i++) {
            int x = x0 + i;
            if (x >= uSize.x) break;
            float dot = texelFetch(uDots, ivec2(x, y), 0).r;
            if (dot <= 0.0) continue;
            if (texelFetch(uPacman, ivec2(x, y), 0).r > uThreshold) {
                sum += dot;
                if (texelFetch(uPower, ivec2(x, y), 0).r > 0.5) pelletSum += dot;
            }
        }
    }

    fragColor = vec4(sum, pelletSum, 0.0, 0.0);
}
