#version 300 es
// FFT mode (fftsim.js): for one pair of kernels (2j, 2j+1), multiply each kernel's spectrum by
// the spectrum of the channel it reads (uSrc0 / uSrc1, one-hot), less, where the bullets mask
// (uBMask), that of the channel times the bullets field; scaled by 1 / (W H) for the inverse.
//   out0 = the numerators' spectra (pair)   out1 = the bullets' spectra times src . bulletCh (pair),
// whose inverse is how much of each kernel the bullets hide.
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uF0, uF1, uF2, uF3;   // spectra of (R, G), (B, b), (R b chR, G b chG), (B b chB, -)
uniform sampler2D uK;                   // the pair's kernel spectra
uniform vec3 uSrc0, uSrc1;
uniform int uBMask;
uniform vec3 uBulletCh;
uniform float uScale;

layout(location = 0) out vec4 out0;
layout(location = 1) out vec4 out1;

vec2 cmul(vec2 a, vec2 b) { return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec4 f0 = texelFetch(uF0, p, 0), f1 = texelFetch(uF1, p, 0), k = texelFetch(uK, p, 0);
    vec2 a0 = uSrc0.x * f0.xy + uSrc0.y * f0.zw + uSrc0.z * f1.xy;
    vec2 a1 = uSrc1.x * f0.xy + uSrc1.y * f0.zw + uSrc1.z * f1.xy;
    vec4 d = vec4(0.0);
    if (uBMask == 1) {
        vec4 f2 = texelFetch(uF2, p, 0), f3 = texelFetch(uF3, p, 0);
        a0 -= uSrc0.x * f2.xy + uSrc0.y * f2.zw + uSrc0.z * f3.xy;
        a1 -= uSrc1.x * f2.xy + uSrc1.y * f2.zw + uSrc1.z * f3.xy;
        d = vec4(cmul(k.xy, f1.zw) * dot(uSrc0, uBulletCh), cmul(k.zw, f1.zw) * dot(uSrc1, uBulletCh)) * uScale;
    }
    out0 = vec4(cmul(k.xy, a0), cmul(k.zw, a1)) * uScale;
    out1 = d;
}
