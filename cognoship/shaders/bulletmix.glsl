#version 300 es
// The two bullets fields as one, for the creatures' step (glsim.js effBullets): the player's and
// the bots' (hostile) bullets added, at most 1, so both act on the creatures exactly as the
// player's bullets do (sim.glsl / the FFT path's uBullets).
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uA;   // the bullets field (R); uVec 1: the world-draw field (Nebula), rgb one per channel
uniform sampler2D uB;   // the hostile bullets field (R)
uniform int  uVec;      // 1: two masks as one per-channel mask: per channel the stronger of A's and the bots' bullets' (on their channels uBCh)
uniform int  uAVec;     // (uVec) 1: A is the world-draw field (Nebula), rgb one per channel; 0: the player's bullets, R on their channels uACh
uniform vec3 uACh;
uniform vec3 uBCh;

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    if (uVec == 1) {
        vec4 a = texelFetch(uA, p, 0);
        fragColor = vec4(max(uAVec == 1 ? a.rgb : a.r * uACh, texelFetch(uB, p, 0).r * uBCh), 1.0);
        return;
    }
    fragColor = vec4(min(texelFetch(uA, p, 0).r + texelFetch(uB, p, 0).r, 1.0), 0.0, 0.0, 1.0);
}
