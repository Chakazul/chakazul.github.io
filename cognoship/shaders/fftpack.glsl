#version 300 es
// FFT mode (fftsim.js): the fields to transform, as complex pairs (real parts; imaginary 0):
//   out0 = (R, G)   out1 = (B, b)   out2 = (R b chR, G b chG)   out3 = (B b chB, 0)
// with b the bullets field where the bullets mask (uBMask; else 0) and ch the channels they act on.
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uState;
uniform sampler2D uBullets;
uniform int uBMask;
uniform vec3 uBulletCh;

layout(location = 0) out vec4 out0;
layout(location = 1) out vec4 out1;
layout(location = 2) out vec4 out2;
layout(location = 3) out vec4 out3;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec3 a = texelFetch(uState, p, 0).rgb;
    float b = uBMask == 1 ? texelFetch(uBullets, p, 0).r : 0.0;
    vec3 ab = a * b * uBulletCh;
    out0 = vec4(a.r, 0.0, a.g, 0.0);
    out1 = vec4(a.b, 0.0, b, 0.0);
    out2 = vec4(ab.r, 0.0, ab.g, 0.0);
    out3 = vec4(ab.b, 0.0, 0.0, 0.0);
}
