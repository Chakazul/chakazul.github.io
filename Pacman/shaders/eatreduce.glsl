#version 300 es
// Stage 1 of the "how much dot mass did Pac-Man just eat" reduction -- same 16x16
// block-tiling as reduce.glsl's CoM reduction, but summing a different quantity: the
// dots-channel mass sitting at cells sim.glsl is about to erase this step (see its
// `uEatEnabled`/`uEatThreshold` check). uDots is channel 2's state *before* this
// step's growth/erasure and uPacman is channel 1's state *after* its own step, which
// is exactly the pairing sim.glsl's own erase check uses -- so this measures the real
// eat condition rather than inferring it from a before/after mass delta (which would
// also pick up channel 2's own, unrelated, growth fluctuation).
//
// Also splits off how much of that was power-pellet mass specifically (uPower, see
// glsim.js's setPowerMask()), which is what tells the CPU side to start the frightened
// window -- into fragColor.g, alongside the total in .r, so both come back from one pass.
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
