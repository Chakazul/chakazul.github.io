#version 300 es
// One Lenia step: toroidal convolution, growth function, then eat/wall masking. Direct-sum
// convolution is the GPU's fast path (every cell an independent thread), replacing the CPU
// demo's FFT, which was itself forced onto the slow Bluestein path by non-power-of-2 board
// sizes. The kernel texture is built by the same arithmetic as the CPU's buildKernel()
// (normalized, 1e-7-thresholded), so weights match bit for bit.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uState;    // R = cell value in [0,1]
uniform sampler2D uWall;     // R > 0.5 where the maze wall is
uniform sampler2D uKernel;   // R = normalized kernel weight, (2R+1)^2, center at (uR,uR)
uniform ivec2 uSize;         // board (W, H)
uniform int   uR;            // kernel radius
uniform float uMu;           // growth center
uniform float uSigma;        // growth width
uniform float uDt;           // time step
uniform sampler2D uEat;      // another channel's state; where it exceeds uEatThreshold, erase
                              // this cell (Pac-Man eating a dot). Ignored unless uEatEnabled != 0.
uniform int   uEatEnabled;
uniform float uEatThreshold;
uniform int   uWallEnabled;  // 0 skips the uWall fetch below entirely -- for a channel whose
                              // solitons never move (the dots), wall collision can never trigger,
                              // so there is nothing worth spending the fetch on.

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    int W = uSize.x, H = uSize.y;

    // Toroidal tap sum -- taps never exceed the kernel radius, so wrapping is one add/subtract.
    float conv = 0.0;
    for (int dy = -uR; dy <= uR; dy++) {
        int y = p.y + dy;
        if (y < 0) y += H; else if (y >= H) y -= H;
        for (int dx = -uR; dx <= uR; dx++) {
            float w = texelFetch(uKernel, ivec2(dx + uR, dy + uR), 0).r;
            if (w == 0.0) continue;          // outside the disc, or below the tap threshold
            int x = p.x + dx;
            if (x < 0) x += W; else if (x >= W) x -= W;
            conv += w * texelFetch(uState, ivec2(x, y), 0).r;
        }
    }

    float s = texelFetch(uState, p, 0).r;
    float d = uMu - conv;
    float v = s + (2.0 * exp(-(d * d) / (2.0 * uSigma * uSigma)) - 1.0) * uDt;
    v = clamp(v, 0.0, 1.0);

    // Eating: wherever the other channel (Pac-Man) has meaningful mass, this cell's mass is
    // erased -- the same one-fetch-and-zero shape as the wall check below, just against a
    // second channel's state texture instead of the fixed wall mask.
    if (uEatEnabled != 0 && texelFetch(uEat, p, 0).r > uEatThreshold) v = 0.0;

    // Wall collision, applied after the step so mass an action dropped into a wall is still
    // visible to this step's convolution, and only then removed.
    if (uWallEnabled != 0 && texelFetch(uWall, p, 0).r > 0.5) v = 0.0;

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
