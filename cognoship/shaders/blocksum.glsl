#version 300 es
// Coarse mass maps of the three single-channel fields at once (one readback for all): texel (i, j)
// = the sums of R over the uBlock x uBlock cells from (i, j) * uBlock (cells past the board edge
// left out), of the bullets (r), the bots' bullets (g) and the dots (b; in their own frame).
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uA;       // bullets channel
uniform sampler2D uB;       // hostile bullets channel
uniform sampler2D uC;       // dots channel
uniform ivec2 uSize;        // board (W, H)
uniform int   uBlock;

out vec4 fragColor;

void main() {
    ivec2 b0 = ivec2(gl_FragCoord.xy) * uBlock;
    vec3 sum = vec3(0.0);
    for (int y = 0; y < uBlock; y++)
        for (int x = 0; x < uBlock; x++) {
            ivec2 q = b0 + ivec2(x, y);
            if (q.x < uSize.x && q.y < uSize.y)
                sum += vec3(texelFetch(uA, q, 0).r, texelFetch(uB, q, 0).r, texelFetch(uC, q, 0).r);
        }
    fragColor = vec4(sum, 0.0);
}
