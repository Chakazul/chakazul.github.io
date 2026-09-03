#version 300 es
// One Lenia step: toroidal convolution with the growth kernel, then the growth
// function, then wall masking. This is the pass that replaces the CPU FFT in
// CARL/maze_playground.html -- the math below is the direct-sum form of the same
// convolution, which is what the FFT there was accelerating. On the GPU the direct
// sum is the fast path: every cell is an independent thread, so the O(knn) tap loop
// costs nothing that matters, and it avoids the Bluestein path the CPU version was
// forced onto by the non-power-of-2 board sizes (100/150/200/250).
//
// The kernel arrives as a (2R+1)^2 texture built by the JS side with the exact same
// code as the CPU demo's buildKernel(), already normalized by its sum and already
// zeroed below the 1e-7 tap threshold. Keeping it a texture rather than recomputing
// quad4() per tap means the weights are bit-identical to the CPU version's, and the
// tiny kernel stays resident in texture cache.
precision highp float;
precision highp sampler2D;

uniform sampler2D uState;    // R = cell value in [0,1]
uniform sampler2D uWall;     // R > 0.5 where the maze wall is
uniform sampler2D uKernel;   // R = normalized kernel weight, (2R+1)^2, center at (uR,uR)
uniform ivec2 uSize;         // board (W, H)
uniform int   uR;            // kernel radius
uniform float uMu;           // growth center
uniform float uSigma;        // growth width
uniform float uDt;           // time step

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    int W = uSize.x, H = uSize.y;

    // Toroidal tap sum. dy/dx never exceed the kernel radius (18) and the board is
    // always >= 100 wide, so wrapping is a single add/subtract, not a modulo.
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

    // Wall collision, applied after the step exactly as applyWallCollision() does --
    // so mass an action dropped into a wall is still visible to this step's
    // convolution, and only then removed.
    if (texelFetch(uWall, p, 0).r > 0.5) v = 0.0;

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
