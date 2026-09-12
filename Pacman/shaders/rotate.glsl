#version 300 es
// Turns one ghost by uTurns quarter-turns about a pivot, leaving every other tile untouched.
// Re-stamping the canonical pattern at the new heading would reset whatever the soliton had
// evolved into (a visible blink); rotating the texture keeps it exactly. A 90deg rotation about
// an *integer* pivot maps texel centres onto texel centres, so it's a permutation -- one
// texelFetch each, no filtering, mass preserved to the bit.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uAtlas;
uniform int   uWin;       // tile edge
uniform int   uTile;      // which ghost is turning
uniform ivec2 uPivot;     // rotation centre, tile-local
uniform int   uTurns;     // clockwise quarter turns to apply, 1..3

out vec4 fragColor;

void main() {
    ivec2 a = ivec2(gl_FragCoord.xy);
    int i = a.x / uWin;
    ivec2 l = ivec2(a.x - i * uWin, a.y);

    // Every other ghost passes through untouched -- this pass covers the whole atlas.
    if (i != uTile) {
        fragColor = vec4(texelFetch(uAtlas, a, 0).r, 0.0, 0.0, 1.0);
        return;
    }

    // Inverse of uTurns clockwise quarter-turns: forward sends (x,y)->(-y,x), so backward is
    // (x,y)->(y,-x) -- same convention as the CPU's rotate90(), so a heading index agrees on both.
    ivec2 d = l - uPivot;
    for (int t = 0; t < uTurns; t++) d = ivec2(d.y, -d.x);
    ivec2 s = d + uPivot;

    // The tile has hard edges, so anything rotated in from outside it is empty.
    float v = (s.x >= 0 && s.x < uWin && s.y >= 0 && s.y < uWin)
            ? texelFetch(uAtlas, ivec2(i * uWin + s.x, s.y), 0).r : 0.0;
    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
