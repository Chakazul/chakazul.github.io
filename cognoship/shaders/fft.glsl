#version 300 es
// FFT mode (?fft=1, fftsim.js): one pass of a mixed-radix Stockham FFT along x or y, in gather form.
// Each texel holds two complex numbers (xy, zw), transformed alike. A pass of radix R after
// passes whose radices multiply to Ns: output d = sum over q < R of
//   in[j + q N/R] * exp(sign * 2 pi i q (k + r Ns) / (Ns R)),
// with r = (d / Ns) % R, k = d % Ns, j = (d / (Ns R)) Ns + k. Passes of radix 4, 3, 2 (or 1, a
// copy) in any order give the DFT in natural order (checked against a direct DFT on the CPU).
// The twiddles come from uTw (row 0 for x, row 1 for y: exp(-2 pi i m / N), computed on the CPU in
// double precision), not from sin / cos, which are coarse on some mobile GPUs.
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uIn;
uniform sampler2D uTw;
uniform int uAxis;     // 0: along x, 1: along y
uniform int uN;        // the length along that axis
uniform int uR;        // this pass's radix (1-4)
uniform int uNs;       // the product of the radices before it
uniform float uSign;   // -1 forward, +1 inverse (unscaled)

out vec4 fragColor;

vec2 cmul(vec2 a, vec2 b) { return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    int d = uAxis == 0 ? p.x : p.y;
    int r = (d / uNs) % uR, k = d % uNs, j = (d / (uNs * uR)) * uNs + k;
    int stride = uN / uR, period = uNs * uR, scale = uN / period;
    vec4 acc = vec4(0.0);
    for (int q = 0; q < 4; q++) {
        if (q >= uR) break;
        int idx = j + q * stride;
        vec4 v = texelFetch(uIn, uAxis == 0 ? ivec2(idx, p.y) : ivec2(p.x, idx), 0);
        vec2 w = texelFetch(uTw, ivec2(((q * (k + r * uNs)) % period) * scale, uAxis), 0).xy;
        if (uSign > 0.0) w.y = -w.y;   // the inverse: the conjugate
        acc.xy += cmul(v.xy, w);
        acc.zw += cmul(v.zw, w);
    }
    fragColor = acc;
}
