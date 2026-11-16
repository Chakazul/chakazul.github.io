#version 300 es
// The scrolling world's seam (glsim.js clearSeam): a field copied with the cells within uBand of
// the seam's column or row (round the torus) set to 0, so what drifts away from the ship is
// gone instead of wrapping round to its front; and (the creatures' field) the cells within the
// discs uDisc (col, row, radius; uDiscN of them, round the torus) set to 0: objects taken away.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uSrc;
uniform ivec2 uSize;
uniform vec2  uSeam;    // (col, row) of the seam lines
uniform float uBand;    // half-width, cells
uniform vec3  uDisc[8];
uniform int   uDiscN;

out vec4 fragColor;

void main() {
    ivec2 p = ivec2(gl_FragCoord.xy);
    vec2 s = vec2(uSize), d = vec2(p) - uSeam;
    d = abs(d - s * floor(d / s + 0.5));
    bool gone = d.x < uBand || d.y < uBand;
    for (int i = 0; i < 8; i++) {
        if (i >= uDiscN) break;
        vec2 e = vec2(p) - uDisc[i].xy;
        e -= s * floor(e / s + 0.5);
        if (dot(e, e) <= uDisc[i].z * uDisc[i].z) gone = true;
    }
    fragColor = gone ? vec4(0.0) : texelFetch(uSrc, p, 0);
}
