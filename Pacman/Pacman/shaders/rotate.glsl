#version 300 es
// Turns one ghost in place: rotates its tile by a quarter turn (or two, or three) about a pivot,
// without touching what the pattern actually is, and without touching any other ghost's tile.
//
// This exists because the obvious way to turn a soliton -- re-stamping the canonical pattern at
// the new heading -- visibly resets it: whatever the running pattern had evolved into is thrown
// away and replaced with the pristine one from the bank, which reads as a blink rather than a
// turn. Rotating the texture instead keeps the exact mass the channel already had.
//
// It is exact, not resampled. A 90deg rotation about an *integer* pivot maps texel centres
// precisely onto texel centres, so every output cell has one source cell and the whole thing is a
// permutation: one texelFetch each, no filtering, no interpolation blur, and mass preserved to the
// bit. That is only true at exact quarter turns -- any other angle would need real resampling.
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

    // Every other ghost passes through untouched. One ghost turning must leave the rest of the
    // pack exactly as it was -- this pass covers the whole atlas, so saying nothing about the
    // other tiles means copying them, not clearing them.
    if (i != uTile) {
        fragColor = vec4(texelFetch(uAtlas, a, 0).r, 0.0, 0.0, 1.0);
        return;
    }

    // Walk the offset *backwards* uTurns quarter turns to find which cell fed this one. One
    // clockwise quarter turn sends offset (x,y) to (-y,x), so its inverse sends (x,y) to (y,-x) --
    // the same convention rotate90()'s k uses on the CPU, so a heading index means the same thing
    // to both.
    ivec2 d = l - uPivot;
    for (int t = 0; t < uTurns; t++) d = ivec2(d.y, -d.x);
    ivec2 s = d + uPivot;

    // The tile has hard edges, so anything rotated in from outside it is empty. The soliton is
    // kept centred and is far narrower than the tile, so nothing real is ever lost here.
    float v = (s.x >= 0 && s.x < uWin && s.y >= 0 && s.y < uWin)
            ? texelFetch(uAtlas, ivec2(i * uWin + s.x, s.y), 0).r : 0.0;
    fragColor = vec4(v, 0.0, 0.0, 1.0);
}
