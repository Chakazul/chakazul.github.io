#version 300 es
// The two bullets fields as one, for the creatures' step (glsim.js effBullets): the player's and
// the bots' (hostile) bullets added, at most 1, so both act on the creatures exactly as the
// player's bullets do (sim.glsl / the FFT path's uBullets).
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uA;   // the bullets field (R)
uniform sampler2D uB;   // the hostile bullets field (R)

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    fragColor = vec4(min(texelFetch(uA, p, 0).r + texelFetch(uB, p, 0).r, 1.0), 0.0, 0.0, 1.0);
}
