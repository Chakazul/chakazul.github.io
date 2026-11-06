#version 300 es
// Stage 2 of the CoM reduction: folds reduce.glsl's 16x16 partials into a 4x1 texture.
// Texel i < 3 is channel i; texel 3 is the whole creature (all channels summed).
// Each texel = (comRow, comCol, mass, valid). Kept on the GPU so sim.glsl can place the mask
// from it with no CPU round trip; the CPU reads the same 4 texels once per frame.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uRed0;
uniform sampler2D uRed1;
uniform sampler2D uRed2;
uniform sampler2D uRed3;
uniform ivec2 uSize;

out vec4 fragColor;

const float TWO_PI = 6.283185307179586;

void main() {
    int which = int(gl_FragCoord.x);   // 0, 1, 2 = channel; 3 = total
    vec3 m = vec3(0.0), sy = vec3(0.0), cy = vec3(0.0), sx = vec3(0.0), cx = vec3(0.0);
    for (int j = 0; j < 16; j++) {
        for (int i = 0; i < 16; i++) {
            ivec2 t = ivec2(i, j);
            vec4 a = texelFetch(uRed0, t, 0);
            vec4 b = texelFetch(uRed1, t, 0);
            vec4 c = texelFetch(uRed2, t, 0);
            vec4 d = texelFetch(uRed3, t, 0);
            m += a.xyz;
            sy += vec3(a.w, b.x, b.y);
            cy += vec3(b.z, b.w, c.x);
            sx += c.yzw;
            cx += d.xyz;
        }
    }

    float M, SY, CY, SX, CX;
    if (which == 3) {
        M = m.r + m.g + m.b;  SY = sy.r + sy.g + sy.b;  CY = cy.r + cy.g + cy.b;
        SX = sx.r + sx.g + sx.b;  CX = cx.r + cx.g + cx.b;
    } else {
        M = m[which];  SY = sy[which];  CY = cy[which];  SX = sx[which];  CX = cx[which];
    }

    if (M < 1e-6) { fragColor = vec4(0.0); return; }

    float H = float(uSize.y), W = float(uSize.x);
    float my = atan(SY, CY) / (TWO_PI / H);  if (my < 0.0) my += H;
    float mx = atan(SX, CX) / (TWO_PI / W);  if (mx < 0.0) mx += W;
    fragColor = vec4(my, mx, M, 1.0);
}
