#version 300 es
// FFT mode (fftsim.js): the real parts of 8 pair textures (kernels 2j in .x, 2j+1 in .z) gathered
// into 4 group textures, kernel 4g + i in component i of out g: the bullets' hidden shares, for
// fftgrow.glsl's uD0..uD3 (so it needs 15 samplers, not 19).
precision highp float;
precision highp int;
precision highp sampler2D;

uniform sampler2D uP0, uP1, uP2, uP3, uP4, uP5, uP6, uP7;

layout(location = 0) out vec4 out0;
layout(location = 1) out vec4 out1;
layout(location = 2) out vec4 out2;
layout(location = 3) out vec4 out3;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    out0 = vec4(texelFetch(uP0, p, 0).xz, texelFetch(uP1, p, 0).xz);
    out1 = vec4(texelFetch(uP2, p, 0).xz, texelFetch(uP3, p, 0).xz);
    out2 = vec4(texelFetch(uP4, p, 0).xz, texelFetch(uP5, p, 0).xz);
    out3 = vec4(texelFetch(uP6, p, 0).xz, texelFetch(uP7, p, 0).xz);
}
