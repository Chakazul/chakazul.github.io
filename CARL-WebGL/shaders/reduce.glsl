#version 300 es
// Stage 1 of the center-of-mass reduction: each output texel sums one block of the
// board into the five accumulators the toroidal circular mean needs.
//
// The CPU demo walked the whole board every step to find the soliton. Reading the
// board back to do that here would defeat the point of running the sim on the GPU, so
// the sum happens on the GPU and only the 1x1 result of stage 2 ever crosses back.
// Output is a fixed 16x16 grid regardless of board size (uBlock = ceil(N/16)), which
// keeps stage 2's serial sum at a constant 256 texels; board sizes here top out at
// 250, so a block never exceeds 16x16 either.
//
// Five accumulators do not fit one RGBA target, hence the two attachments.
precision highp float;
precision highp sampler2D;

uniform sampler2D uState;
uniform ivec2 uSize;    // board (W, H)
uniform int   uBlock;   // block edge = ceil(max(W,H)/16)

layout(location = 0) out vec4 o0;   // (mass, sum v*sinY, sum v*cosY, sum v*sinX)
layout(location = 1) out vec4 o1;   // (sum v*cosX, -, -, -)

const float TWO_PI = 6.283185307179586;

void main() {
    ivec2 b = ivec2(gl_FragCoord.xy);          // block index, 0..15
    int x0 = b.x * uBlock, y0 = b.y * uBlock;

    float twoPiH = TWO_PI / float(uSize.y);
    float twoPiW = TWO_PI / float(uSize.x);

    float mass = 0.0, sy = 0.0, cy = 0.0, sx = 0.0, cx = 0.0;
    for (int j = 0; j < 16; j++) {
        if (j >= uBlock) break;
        int y = y0 + j;
        if (y >= uSize.y) break;
        float ay = float(y) * twoPiH;
        float sinY = sin(ay), cosY = cos(ay);
        for (int i = 0; i < 16; i++) {
            if (i >= uBlock) break;
            int x = x0 + i;
            if (x >= uSize.x) break;
            float v = texelFetch(uState, ivec2(x, y), 0).r;
            if (v <= 0.0) continue;
            float ax = float(x) * twoPiW;
            mass += v;
            sy += v * sinY;  cy += v * cosY;
            sx += v * sin(ax);  cx += v * cos(ax);
        }
    }

    o0 = vec4(mass, sy, cy, sx);
    o1 = vec4(cx, 0.0, 0.0, 0.0);
}
