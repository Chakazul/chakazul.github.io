#version 300 es
// Stage 1 of the CoM reduction (after Carlman's reduce.glsl), per channel. Each output texel
// sums one board block into toroidal-mean accumulators: mass, and mass-weighted sin/cos of the
// row and column angles, for each of the three channels (15 numbers, packed in 4 targets).
// Per-channel CoMs give the body axis (G leads B along the heading), see app.js. Only cells within
// uRadius of uCenter (torus distance) count, so the CoMs are the focused creature's, not the
// whole board's (uRadius <= 0 counts everything).
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uState;
uniform ivec2 uSize;    // board (W, H)
uniform int   uBlock;   // block edge = ceil(max(W, H) / 16)
uniform vec2  uCenter;  // window centre, cells (col, row)
uniform float uRadius;  // window radius, cells

layout(location = 0) out vec4 o0;   // (m.r, m.g, m.b, sinY.r)
layout(location = 1) out vec4 o1;   // (sinY.g, sinY.b, cosY.r, cosY.g)
layout(location = 2) out vec4 o2;   // (cosY.b, sinX.r, sinX.g, sinX.b)
layout(location = 3) out vec4 o3;   // (cosX.r, cosX.g, cosX.b, -)

const float TWO_PI = 6.283185307179586;

void main() {
    ivec2 b = ivec2(gl_FragCoord.xy);
    int x0 = b.x * uBlock, y0 = b.y * uBlock;
    float twoPiH = TWO_PI / float(uSize.y);
    float twoPiW = TWO_PI / float(uSize.x);

    vec3 m = vec3(0.0), sy = vec3(0.0), cy = vec3(0.0), sx = vec3(0.0), cx = vec3(0.0);
    for (int j = 0; j < uBlock; j++) {
        int y = y0 + j;
        if (y >= uSize.y) break;
        float ay = float(y) * twoPiH;
        float sinY = sin(ay), cosY = cos(ay);
        for (int i = 0; i < uBlock; i++) {
            int x = x0 + i;
            if (x >= uSize.x) break;
            vec3 v = texelFetch(uState, ivec2(x, y), 0).rgb;
            if (v.r + v.g + v.b <= 0.0) continue;
            if (uRadius > 0.0) {
                vec2 s = vec2(uSize);
                vec2 d = vec2(x, y) - uCenter;
                d -= s * floor(d / s + 0.5);
                if (dot(d, d) > uRadius * uRadius) continue;
            }
            float ax = float(x) * twoPiW;
            m += v;
            sy += v * sinY;  cy += v * cosY;
            sx += v * sin(ax);  cx += v * cos(ax);
        }
    }

    o0 = vec4(m, sy.r);
    o1 = vec4(sy.g, sy.b, cy.r, cy.g);
    o2 = vec4(cy.b, sx);
    o3 = vec4(cx, 0.0);
}
