#version 300 es
// Applies one intervention (add/remove mass) as a disc into a scratch texture the sim pass then
// steps from. A separate pass rather than folded into sim.glsl, so the ~1000-tap convolution
// isn't re-testing the action disc at every neighbour; skipped entirely on no-op steps.
precision highp float;
precision highp sampler2D;
// Required, not cosmetic: GLSL ES 3.00 defaults fragment-stage int to mediump (only guaranteed
// 16-bit), which most mobile GPUs implement literally (desktop drivers quietly hand out 32 bits).
// The squared distance below reaches ~75625, which wraps mod 65536 under 16-bit int, turning the
// action radius into a board-spanning ring of injected mass -- exploding the soliton on the first
// action on affected devices.
precision highp int;

uniform sampler2D uState;
uniform ivec2 uSize;      // board (W, H)
uniform ivec2 uCenter;    // action center in board coords (x = col, y = row)
uniform float uDelta;     // +/- action magnitude
uniform int   uRadius;    // action radius in cells

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    float v = texelFetch(uState, p, 0).r;

    // Shortest toroidal offset to the action center.
    int dx = p.x - uCenter.x;
    if (dx >  uSize.x / 2) dx -= uSize.x;
    if (dx < -uSize.x / 2) dx += uSize.x;
    int dy = p.y - uCenter.y;
    if (dy >  uSize.y / 2) dy -= uSize.y;
    if (dy < -uSize.y / 2) dy += uSize.y;

    if (dx * dx + dy * dy <= uRadius * uRadius) v = clamp(v + uDelta, 0.0, 1.0);

    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
