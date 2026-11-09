#version 300 es
// Coarse mass map of a single-channel field (the bullets): texel (i, j) = the sum of R over the
// uBlock x uBlock cells from (i, j) * uBlock (cells past the board edge left out).
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uField;
uniform ivec2 uSize;        // board (W, H)
uniform int   uBlock;

out vec4 fragColor;

void main() {
    ivec2 b0 = ivec2(gl_FragCoord.xy) * uBlock;
    float sum = 0.0;
    for (int y = 0; y < uBlock; y++)
        for (int x = 0; x < uBlock; x++) {
            ivec2 q = b0 + ivec2(x, y);
            if (q.x < uSize.x && q.y < uSize.y) sum += texelFetch(uField, q, 0).r;
        }
    fragColor = vec4(sum, 0.0, 0.0, 0.0);
}
