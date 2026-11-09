#version 300 es
// One Lenia step of the dots channel (ported from Carlman): a single-channel field of small
// stationary solitons under their own rule, which never moves them, only makes them breathe. It is
// stepped right after the creatures' step, and eating is applied to its result: wherever any
// creature channel of the just-stepped board is above uEat, the dot mass there is erased. The
// dots do not act back on the creatures.
//
// The field drifts to the left (as the star field's near layer): this texture is the dots' own
// frame, and dot cell x is at board cell x - uOffX, so eating reads the board there.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump on mobile (see Carlman action.glsl)

uniform sampler2D uDots;    // dots channel, R = cell value
uniform sampler2D uState;   // the creatures' board, just stepped (rgb = channels)
uniform sampler2D uTaps;    // width = tap count, row 0: (dx, dy, weight, -); weights sum to 1
uniform ivec2 uSize;        // board (W, H)
uniform int   uTapCount;
uniform float uMu;          // growth centre
uniform float uSigma;       // growth width (Gaussian)
uniform float uDt;
uniform float uEat;         // creature mass above this erases the dot mass at that cell
uniform int   uOffX;        // the drift, whole cells in [0, W)

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    float conv = 0.0;
    for (int t = 0; t < uTapCount; t++) {
        vec4 tap = texelFetch(uTaps, ivec2(t, 0), 0);
        ivec2 q = (p + ivec2(tap.xy) + uSize) % uSize;
        conv += tap.z * texelFetch(uDots, q, 0).r;
    }
    float d = conv - uMu;
    float v = texelFetch(uDots, p, 0).r + uDt * (2.0 * exp(-d * d / (2.0 * uSigma * uSigma)) - 1.0);
    v = clamp(v, 0.0, 1.0);
    vec3 a = texelFetch(uState, ivec2((p.x - uOffX + uSize.x) % uSize.x, p.y), 0).rgb;
    if (max(a.r, max(a.g, a.b)) > uEat) v = 0.0;
    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
