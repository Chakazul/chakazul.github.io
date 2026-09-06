#version 300 es
// Applies one intervention -- CARL's or the user's -- as a disc of added or removed
// mass, writing the result to a scratch texture that the sim pass then steps from.
//
// This is a separate pass rather than something folded into sim.glsl on purpose: the
// sim pass reads ~1000 neighbours per cell, so having it apply the action inline would
// mean re-testing the action disc at every one of those taps. A standalone pass costs
// one fetch per cell instead, which is nothing next to the convolution. It is skipped
// entirely on no-op steps.
precision highp float;
precision highp sampler2D;

uniform sampler2D uState;
uniform ivec2 uSize;      // board (W, H)
uniform ivec2 uCenter;    // action center in board coords (x = col, y = row)
uniform float uDelta;     // +/- action magnitude
uniform int   uRadius;    // action radius in cells

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    float v = texelFetch(uState, p, 0).r;

    // Shortest toroidal offset to the action center. The radius (7) is far smaller
    // than half the board, so the nearest-image offset is the right one.
    int dx = p.x - uCenter.x;
    if (dx >  uSize.x / 2) dx -= uSize.x;
    if (dx < -uSize.x / 2) dx += uSize.x;
    int dy = p.y - uCenter.y;
    if (dy >  uSize.y / 2) dy -= uSize.y;
    if (dy < -uSize.y / 2) dy += uSize.y;

    if (dx * dx + dy * dy <= uRadius * uRadius) v = clamp(v + uDelta, 0.0, 1.0);

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
