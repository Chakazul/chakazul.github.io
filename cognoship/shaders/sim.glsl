#version 300 es
// One Lenia step for the 3-channel board, with cognitive masking (Cool et al. 2026, Eq. 2).
//
// Direct convolution over a tap table (rule.js buildTaps): per tap, ONE board fetch gives all
// three channels, and each vec4 group of kernels picks its source channels with a 3x4 selection
// matrix, so 15 kernels cost one state fetch plus four weight fetches per tap.
//
// Masking hides channel c inside a circle that follows the focused creature (core.js passes the
// circle's absolute centre each step, extrapolating the creature's CoM). Every
// kernel reading channel c then sees U = K*(A(1-B)) / K*(1-B): the occluded cells drop out of
// the numerator, and the denominator renormalises by the visible kernel mass. The mask hides
// information only; it never removes mass. Pixels whose kernel can't reach the circle take the
// plain path, where the denominator is exactly 1 (kernels are normalised).
//
// The same circle can instead erase: uMaskOn = 2 sets the chosen channel to 0 inside it after
// the update (removing mass, like Carlman's interventions), with no occlusion; or draw:
// uMaskOn = 3 sets it to 1 there.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump on mobile (see Carlman action.glsl)

uniform sampler2D uState;   // RGBA32F, rgb = channels R, G, B
uniform sampler2D uTaps;    // width = tap count; row 0 = (dx, dy), rows 1..4 = weights of group 0..3
uniform ivec2 uSize;        // board (W, H)
uniform int   uTapCount;

uniform mat3x4 uSrc[4];     // channel -> kernel-slot selection, per group
uniform mat4x3 uDst[4];     // kernel-slot -> channel selection, per group
uniform vec4  uM[4];
uniform vec4  uS[4];
uniform vec4  uH[4];
uniform vec3  uGain;        // per destination channel: 1 for LeniaF (h already normalised by
                            // free_h) and the unnormalised legacy update, or 1 / sum h for the
                            // legacy LeniaNDKC-style normalised update
uniform float uDt;
uniform int   uGrowth;      // growth function (LeniaF numbering): 1 = polynomial quad4, 2 = Gaussian

uniform int   uMaskOn;      // 0 off, 1 cognitive mask, 2 erase mass, 3 draw mass
uniform vec2  uMaskC;       // mask centre, board cells (col, row)
uniform float uMaskR;       // mask radius, cells
uniform vec3  uMaskCh;      // one-hot: which channel the mask hides
uniform float uReach;       // kernel reach, cells (max tap distance)

out vec4 fragColor;

vec2 wrapDelta(vec2 d) {
    vec2 s = vec2(uSize);
    return d - s * floor(d / s + 0.5);
}

// Growth in [-1, 1] for four kernels at once.
vec4 growth(vec4 u, vec4 m, vec4 s) {
    vec4 d = u - m;
    if (uGrowth == 1) {
        vec4 q = max(vec4(0.0), 1.0 - d * d / (9.0 * s * s));
        q *= q;
        return q * q * 2.0 - 1.0;
    }
    return exp(-d * d / (2.0 * s * s)) * 2.0 - 1.0;
}

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);   // (col, row)

    // Where is the mask, and can this pixel's kernel see it at all?
    bool near = false;
    vec2 mc = vec2(0.0);
    bool erase = false;
    if (uMaskOn != 0) {
        mc = uMaskC;
        vec2 d = wrapDelta(vec2(p) - mc);
        if (uMaskOn == 1) {
            float lim = uReach + uMaskR + 1.0;
            near = dot(d, d) <= lim * lim;
        } else {
            erase = dot(d, d) <= uMaskR * uMaskR;
        }
    }

    vec4 U0 = vec4(0.0), U1 = vec4(0.0), U2 = vec4(0.0), U3 = vec4(0.0);
    vec4 N0 = vec4(0.0), N1 = vec4(0.0), N2 = vec4(0.0), N3 = vec4(0.0);
    float r2 = uMaskR * uMaskR;

    for (int t = 0; t < uTapCount; t++) {
        ivec2 off = ivec2(texelFetch(uTaps, ivec2(t, 0), 0).xy);
        ivec2 q = (p + off + uSize) % uSize;
        vec3 a = texelFetch(uState, q, 0).rgb;
        vec4 w0 = texelFetch(uTaps, ivec2(t, 1), 0);
        vec4 w1 = texelFetch(uTaps, ivec2(t, 2), 0);
        vec4 w2 = texelFetch(uTaps, ivec2(t, 3), 0);
        vec4 w3 = texelFetch(uTaps, ivec2(t, 4), 0);
        if (near) {
            vec2 d = wrapDelta(vec2(q) - mc);
            vec3 vis = dot(d, d) <= r2 ? vec3(1.0) - uMaskCh : vec3(1.0);
            a *= vis;
            N0 += w0 * (uSrc[0] * vis);
            N1 += w1 * (uSrc[1] * vis);
            N2 += w2 * (uSrc[2] * vis);
            N3 += w3 * (uSrc[3] * vis);
        }
        U0 += w0 * (uSrc[0] * a);
        U1 += w1 * (uSrc[1] * a);
        U2 += w2 * (uSrc[2] * a);
        U3 += w3 * (uSrc[3] * a);
    }

    if (near) {
        // A kernel whose whole reach is occluded sees nothing: U = 0.
        U0 = step(1e-6, N0) * U0 / max(N0, 1e-6);
        U1 = step(1e-6, N1) * U1 / max(N1, 1e-6);
        U2 = step(1e-6, N2) * U2 / max(N2, 1e-6);
        U3 = step(1e-6, N3) * U3 / max(N3, 1e-6);
    }

    vec3 delta = vec3(0.0);
    delta += uDst[0] * (uH[0] * growth(U0, uM[0], uS[0]));
    delta += uDst[1] * (uH[1] * growth(U1, uM[1], uS[1]));
    delta += uDst[2] * (uH[2] * growth(U2, uM[2], uS[2]));
    delta += uDst[3] * (uH[3] * growth(U3, uM[3], uS[3]));

    vec3 A = texelFetch(uState, p, 0).rgb;
    vec3 Anew = clamp(A + uDt * delta * uGain, 0.0, 1.0);
    if (erase) Anew = uMaskOn == 3 ? max(Anew, uMaskCh) : Anew * (vec3(1.0) - uMaskCh);
    fragColor = vec4(Anew, 1.0);
}
