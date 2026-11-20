#version 300 es
// FFT mode (fftsim.js): the fields to transform, two real fields to a complex number (one as the
// real part, one as the imaginary: fftmul.glsl separates their spectra again):
//   out0 = (R + i G, B + i b)   out1 = (R b chR + i G b chG, B b chB)
// with b the bullets field where the bullets mask (uBMask; else 0) and ch the channels they act on.
// uBVec = 1 (the world-draw field, its rgb one per channel): out1 takes b's own channel
// (R bR chR, ...) and out0's b is 0 (fftmul reads the field's spectra from uFB0 instead, this
// pass's out0 of the field itself).
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

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec3 a = texelFetch(uState, p, 0).rgb;
    vec4 t = uBMask == 1 ? texelFetch(uBullets, p, 0) : vec4(0.0);
    float b = uBVec == 1 ? 0.0 : t.r;
    vec3 ab = a * (uBVec == 1 ? t.rgb : vec3(t.r)) * uBulletCh;
    out0 = vec4(a.r, a.g, a.b, b);
    out1 = vec4(ab.r, ab.g, ab.b, 0.0);
}
