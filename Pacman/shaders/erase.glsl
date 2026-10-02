#version 300 es
// The dots channel's off-step (see setDotsEvery() in glsim.js): no Lenia step, just the eat
// check sim.glsl applies after one -- wherever Pac-Man's mass exceeds uEatThreshold, the dot mass
// there is erased. Keeps eating instant while the dots' own dynamics run at a fraction of the
// step rate. Three fetches a cell, against hundreds for a full step.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uState;    // dots channel, R = cell value
uniform sampler2D uEat;      // Pac-Man's just-stepped state
uniform float uEatThreshold;

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    float v = texelFetch(uState, p, 0).r;
    if (texelFetch(uEat, p, 0).r > uEatThreshold) v = 0.0;
    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
