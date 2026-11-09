#version 300 es
// One Lenia step of the bullets channel: a single-channel field of gliders (Scutium solidus by
// default) under the classic Lenia rule (polynomial kernel core and growth, as Lenia's animals.json
// kn = gn = 1). Run after the creatures' step, possibly several times per step (glsim.js
// setBullets sub). The creatures wear the bullets down: the bullet mass at a cell drops by uWear
// times the creature mass there (the largest of the uWearCh channels), every bullet step. Slower than
// the dots, which are erased outright, so a bullet sinks some way into a body before it is gone;
// meanwhile it takes mass out of the body (sim.glsl uBullet*).
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump on mobile (see Carlman action.glsl)

uniform sampler2D uField;   // bullets channel, R = cell value
uniform sampler2D uState;   // the creatures' board, just stepped (rgb = channels)
uniform sampler2D uTaps;    // width = tap count, row 0: (dx, dy, weight, -); weights sum to 1
uniform ivec2 uSize;        // board (W, H)
uniform int   uTapCount;
uniform float uMu;          // growth centre
uniform float uSigma;       // growth width (polynomial)
uniform float uDt;
uniform float uWear;        // bullet mass removed per step per unit of creature mass
uniform vec3  uWearCh;      // which creature channels wear the bullets (one-hot or all ones)

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    float conv = 0.0;
    for (int t = 0; t < uTapCount; t++) {
        vec4 tap = texelFetch(uTaps, ivec2(t, 0), 0);
        ivec2 q = (p + ivec2(tap.xy) + uSize) % uSize;
        conv += tap.z * texelFetch(uField, q, 0).r;
    }
    float d = conv - uMu;
    float g = max(0.0, 1.0 - d * d / (9.0 * uSigma * uSigma));
    g *= g; g *= g;
    float v = texelFetch(uField, p, 0).r + uDt * (2.0 * g - 1.0);
    vec3 a = texelFetch(uState, p, 0).rgb * uWearCh;
    v = clamp(v - uWear * max(a.r, max(a.g, a.b)), 0.0, 1.0);
    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
