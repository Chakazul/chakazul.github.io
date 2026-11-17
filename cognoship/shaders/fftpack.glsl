#version 300 es
// FFT mode (fftsim.js): the fields to transform, as complex pairs (real parts; imaginary 0):
//   out0 = (R, G)   out1 = (B, b)   out2 = (R b chR, G b chG)   out3 = (B b chB, 0)
// with b the bullets field where the bullets mask (uBMask; else 0) and ch the channels they act on.
// uBVec = 1 (the world-draw field, its rgb one per channel): out2 / out3 take b's own channel
// (R bR chR, ...) and out1's b is 0 (fftmul reads the field's spectra from uFB0 / uFB1 instead).
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uState;
uniform sampler2D uBullets;
uniform int uBMask;
uniform vec3 uBulletCh;
uniform int uBVec;

layout(location = 0) out vec4 out0;
layout(location = 1) out vec4 out1;
layout(location = 2) out vec4 out2;
layout(location = 3) out vec4 out3;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec3 a = texelFetch(uState, p, 0).rgb;
    vec4 t = uBMask == 1 ? texelFetch(uBullets, p, 0) : vec4(0.0);
    float b = uBVec == 1 ? 0.0 : t.r;
    vec3 ab = a * (uBVec == 1 ? t.rgb : vec3(t.r)) * uBulletCh;
    out0 = vec4(a.r, 0.0, a.g, 0.0);
    out1 = vec4(a.b, 0.0, b, 0.0);
    out2 = vec4(ab.r, 0.0, ab.g, 0.0);
    out3 = vec4(ab.b, 0.0, 0.0, 0.0);
}
