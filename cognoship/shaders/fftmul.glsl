#version 300 es
// FFT mode (fftsim.js): for one group of 4 kernels (4g .. 4g + 3), multiply each kernel's spectrum
// by the spectrum of the channel it reads (uSrc[i], one-hot), less, where the bullets mask (uBMask),
// that of the channel times the bullets field; scaled by 1 / (W H) for the inverse.
//   out0 = the numerators' spectra   out1 = the bullets' spectra times src . bulletCh, whose
// inverse is how much of each kernel the bullets hide.
// The fields came in two to a complex number (fftpack.glsl): z = FFT(x + i y), so with zm = z at
// the mirrored frequency (-k), X = (z + conj(zm)) / 2 and Y = (z - conj(zm)) / (2 i) (x, y real).
// The results go out the same way, two kernels to a complex number, S0 + i S1 and S2 + i S3: their
// inverse FFT is U0 + i U1, U2 + i U3 (each U real), the group's 4 kernels in x, y, z, w.
// uBVec = 1 (the world-draw field, its rgb one per channel): the bullets' spectrum is that of the
// kernel's own channel of the field, uFB0 = FFT(bR + i bG, bB) (fftsim.js makes it once per change
// of the static field).
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uF0, uF1;    // FFT(R + i G, B + i b), FFT(R b chR + i G b chG, B b chB)
uniform sampler2D uK0, uK1;    // the kernels' spectra: 4g, 4g + 1 in uK0's xy, zw; 4g + 2, 4g + 3 in uK1's
uniform vec3 uSrc[4];
uniform int uBMask;
uniform int uBVec;
uniform sampler2D uFB0;
uniform vec3 uBulletCh;
uniform float uScale;
uniform ivec2 uSize;

layout(location = 0) out vec4 out0;
layout(location = 1) out vec4 out1;

vec2 cmul(vec2 a, vec2 b) { return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
vec2 conj2(vec2 a) { return vec2(a.x, -a.y); }
vec2 times_i(vec2 a) { return vec2(-a.y, a.x); }

// the spectra of the 4 real fields packed in t (xy: x0 + i y0, zw: x1 + i y1), at p (z) and -p (zm)
void split(vec4 z, vec4 zm, out vec2 x0, out vec2 y0, out vec2 x1, out vec2 y1) {
    vec2 s = z.xy + conj2(zm.xy), d = z.xy - conj2(zm.xy);
    x0 = 0.5 * s; y0 = 0.5 * vec2(d.y, -d.x);   // d / (2 i)
    s = z.zw + conj2(zm.zw); d = z.zw - conj2(zm.zw);
    x1 = 0.5 * s; y1 = 0.5 * vec2(d.y, -d.x);
}

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    ivec2 pm = (uSize - p) % uSize;
    vec2 R, G, B, b;
    split(texelFetch(uF0, p, 0), texelFetch(uF0, pm, 0), R, G, B, b);
    vec2 Rb = vec2(0.0), Gb = vec2(0.0), Bb = vec2(0.0), bR = b, bG = b, bB = b, unused;
    bool bmask = uBMask == 1;
    if (bmask) {
        split(texelFetch(uF1, p, 0), texelFetch(uF1, pm, 0), Rb, Gb, Bb, unused);
        if (uBVec == 1) split(texelFetch(uFB0, p, 0), texelFetch(uFB0, pm, 0), bR, bG, bB, unused);
    }
    vec4 k0 = texelFetch(uK0, p, 0), k1 = texelFetch(uK1, p, 0);
    vec2 K[4] = vec2[4](k0.xy, k0.zw, k1.xy, k1.zw);
    vec2 S[4], Dk[4];
    for (int i = 0; i < 4; i++) {
        vec3 s = uSrc[i];
        vec2 a = s.x * (R - Rb) + s.y * (G - Gb) + s.z * (B - Bb);
        S[i] = cmul(K[i], a);
        Dk[i] = bmask ? cmul(K[i], s.x * bR + s.y * bG + s.z * bB) * dot(s, uBulletCh) : vec2(0.0);
    }
    out0 = vec4(S[0] + times_i(S[1]), S[2] + times_i(S[3])) * uScale;
    out1 = vec4(Dk[0] + times_i(Dk[1]), Dk[2] + times_i(Dk[3])) * uScale;
}
