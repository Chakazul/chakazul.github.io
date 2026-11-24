#version 300 es
// FFT mode (?fft=1, fftsim.js): the creatures' step from the convolutions the FFT made, the same
// update as sim.glsl (read that first): U per kernel from the inverse FFTs uP0..uP3 (kernel 4g + i
// in component i of uPg: fftmul.glsl), already K*(A (1 - bullets)) where the bullets mask; uD0..uD3
// (laid out alike) how much of each kernel the bullets hide (K*b, times src . bulletCh): visible 1 - D.
// The cognitive mask's circle is small, so instead of more FFTs its cells are taken out directly:
// for a pixel within reach of it, each circle cell q takes w(q - p) of its hidden channel out of the
// numerator and the denominator (w from uKimg, the kernels as an image). Erase / draw and the
// bullets' erase act after the update, as in sim.glsl.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uP0, uP1, uP2, uP3;
uniform sampler2D uD0, uD1, uD2, uD3;
uniform sampler2D uState;
uniform sampler2D uBullets;
uniform sampler2D uKimg;    // (2 KR + 1) x (2 KR + 1) per group, groups side by side: w at offset (dx, dy)
uniform int   uKR;
uniform ivec2 uSize;

uniform mat3x4 uSrc[4];
uniform mat4x3 uDst[4];
uniform vec4  uM[4];
uniform vec4  uS[4];
uniform vec4  uH[4];
uniform vec3  uGain;
uniform float uDt;
uniform int   uGrowth;

uniform int   uMaskOn;
uniform vec2  uMaskC;
uniform float uMaskR;
uniform vec3  uMaskCh;
const int XMAX = 8;          // more cognitive masks (the game's bots): uXN circles (col, row, radius, channel); apart from each other
uniform int   uXN;
uniform vec4  uX[XMAX];
uniform float uReach;
uniform int   uBulletMode;
uniform vec3  uBulletCh;
uniform float uBulletHit;
uniform sampler2D uCap;     // Antimatter (the world-draw field, erasing) beside the bullets: each channel capped at 1 - its field (glsim.js capOn)
uniform int   uCapOn;
uniform int   uBulletVec;   // 1: the bullets texture is the world-draw field, its rgb one per channel (sim.glsl)
uniform float uInvK;        // 1 / kernel count

layout(location = 0) out vec4 fragColor;
// the potential, for the creature finder (blocks.glsl): the mean over kernels of U, as the FFT gave
// it (the cognitive mask's circle not taken out; where the bullets mask, what they leave visible)
layout(location = 1) out vec4 potColor;

vec2 wrapDelta(vec2 d) {
    vec2 s = vec2(uSize);
    return d - s * floor(d / s + 0.5);
}

vec3 bulletAt(ivec2 q) {
    vec4 t = texelFetch(uBullets, q, 0);
    return (uBulletVec == 1 ? t.rgb : vec3(t.r)) * uBulletCh;
}

// Take a mask circle's cells (centre mc, radius mr, hiding channel ch) out of U and N at pixel p:
// each circle cell q takes w(q - p) of its hidden channel out of the numerator and the denominator
void hideCircle(ivec2 p, vec2 mc, float rad, vec3 ch, bool bmask,
                inout vec4 U0, inout vec4 U1, inout vec4 U2, inout vec4 U3,
                inout vec4 N0, inout vec4 N1, inout vec4 N2, inout vec4 N3) {
        int mr = int(ceil(rad)) + 1, side = 2 * uKR + 1;
        ivec2 c0 = ivec2(floor(mod(mc, vec2(uSize)) + 0.5));   // in [0, size]: the sums below stay >= 0 (% of a negative int is undefined in GLSL)
        float r2 = rad * rad;
        for (int dy = -mr; dy <= mr; dy++) {
            for (int dx = -mr; dx <= mr; dx++) {
                ivec2 q = (c0 + ivec2(dx, dy) + 2 * uSize) % uSize;
                vec2 dq = wrapDelta(vec2(q) - mc);
                if (dot(dq, dq) > r2) continue;
                ivec2 o = ivec2(wrapDelta(vec2(q - p)));   // the kernel offset of q from p
                if (abs(o.x) > uKR || abs(o.y) > uKR) continue;
                vec3 vb = bmask ? vec3(1.0) - bulletAt(q) : vec3(1.0);
                vec3 hid = vb * ch;                           // what the circle hides at q
                vec3 ah = texelFetch(uState, q, 0).rgb * hid;
                ivec2 t = o + uKR;
                vec4 w0 = texelFetch(uKimg, t, 0), w1 = texelFetch(uKimg, t + ivec2(side, 0), 0);
                vec4 w2 = texelFetch(uKimg, t + ivec2(2 * side, 0), 0), w3 = texelFetch(uKimg, t + ivec2(3 * side, 0), 0);
                U0 -= w0 * (uSrc[0] * ah); N0 -= w0 * (uSrc[0] * hid);
                U1 -= w1 * (uSrc[1] * ah); N1 -= w1 * (uSrc[1] * hid);
                U2 -= w2 * (uSrc[2] * ah); N2 -= w2 * (uSrc[2] * hid);
                U3 -= w3 * (uSrc[3] * ah); N3 -= w3 * (uSrc[3] * hid);
            }
        }
}

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
    ivec2 p = ivec2(gl_FragCoord.xy);

    bool near = false, erase = false;
    vec2 mc = uMaskC;
    if (uMaskOn != 0) {
        vec2 d = wrapDelta(vec2(p) - mc);
        if (uMaskOn == 1) {
            float lim = uReach + uMaskR + 1.0;
            near = dot(d, d) <= lim * lim;
        } else {
            erase = dot(d, d) <= uMaskR * uMaskR;
        }
    }
    bool bmask = uBulletMode == 2;
    bool masked = near || bmask;   // (and any extra mask within reach: below)

    vec4 U0 = texelFetch(uP0, p, 0), U1 = texelFetch(uP1, p, 0);
    vec4 U2 = texelFetch(uP2, p, 0), U3 = texelFetch(uP3, p, 0);
    potColor = vec4(dot(U0 + U1 + U2 + U3, vec4(1.0)) * uInvK, 0.0, 0.0, 0.0);
    vec4 N0 = vec4(1.0), N1 = vec4(1.0), N2 = vec4(1.0), N3 = vec4(1.0);
    if (bmask) {
        N0 -= texelFetch(uD0, p, 0); N1 -= texelFetch(uD1, p, 0);
        N2 -= texelFetch(uD2, p, 0); N3 -= texelFetch(uD3, p, 0);
    }

    if (near) hideCircle(p, mc, uMaskR, uMaskCh, bmask, U0, U1, U2, U3, N0, N1, N2, N3);
    for (int i = 0; i < XMAX; i++) {   // the extra masks (the game's bots)
        if (i >= uXN) break;
        vec2 d = wrapDelta(vec2(p) - uX[i].xy);
        float lim = uReach + uX[i].z + 1.0;
        if (dot(d, d) > lim * lim) continue;
        masked = true;
        hideCircle(p, uX[i].xy, uX[i].z, vec3(equal(ivec3(int(uX[i].w)), ivec3(0, 1, 2))), bmask, U0, U1, U2, U3, N0, N1, N2, N3);
    }

    if (masked) {
        // a kernel whose whole reach is hidden sees nothing: U = 0
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
    // the bullets' erase takes mass out gradually; the world-draw field's (Antimatter) caps it at once: at most 1 - field
    if (uBulletMode == 1) Anew = uBulletVec == 1 ? min(Anew, vec3(1.0) - bulletAt(p)) : max(Anew - uBulletHit * bulletAt(p), 0.0);
    if (uCapOn == 1) Anew = min(Anew, vec3(1.0) - texelFetch(uCap, p, 0).rgb);
    fragColor = vec4(Anew, 1.0);
}
