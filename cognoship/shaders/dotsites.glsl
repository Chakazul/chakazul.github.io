#version 300 es
// Per-dot presence, for scoring (ported from Carlman): texel i sums the dots channel in a
// (2 uRadius + 1)^2 window around dot i's centre. A dot is bistable (bitten deep enough, it
// dissolves to ~0 within a few steps; grazed or missed, it stays near full mass), so the CPU
// needs only one threshold on this sum.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uDots;    // dots channel, R = cell value
uniform sampler2D uSites;   // texel (i, 0).xy = dot i's centre (col, row), whole numbers
uniform ivec2 uSize;        // board (W, H)
uniform int   uRadius;      // window half-width, far below the board edge

out vec4 fragColor;         // (mass around dot i, -, -, -)

void main() {
    int i = int(gl_FragCoord.x);
    ivec2 c = ivec2(floor(texelFetch(uSites, ivec2(i, 0), 0).xy + 0.5));
    float sum = 0.0;
    for (int dy = -uRadius; dy <= uRadius; dy++) {
        int y = (c.y + dy + uSize.y) % uSize.y;
        for (int dx = -uRadius; dx <= uRadius; dx++)
            sum += texelFetch(uDots, ivec2((c.x + dx + uSize.x) % uSize.x, y), 0).r;
    }
    fragColor = vec4(sum, 0.0, 0.0, 0.0);
}
