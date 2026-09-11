#version 300 es
// Per-dot presence, for scoring: output texel i sums the dots-channel mass left in a square window
// centred on where dot i was stamped. The dots never move, so a fixed window is all it takes, and a
// dot is bistable -- one Pac-Man bites into dissolves all the way to 0 within a few steps, one he
// misses stays near its full mass -- so the CPU side only has to threshold each sum once.
//
// Not the eaten-mass reduction (eatreduce.glsl) with a different threshold: that measures only the
// mass under Pac-Man's pixels, and the part of a bitten dot he never covered dissolves on its own
// afterwards, uncounted. Not the channel's total mass (reduce.glsl -> com.glsl) either: every dot
// breathes in phase with every other, so the total swings by many dots' worth.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uState;   // channel 2 (dots + pellets), post-step
uniform sampler2D uSites;   // texel (i,0).xy = dot i's board position (col, row)
uniform ivec2 uSize;        // board (W, H)
uniform int   uRadius;      // window half-width; always far below the board edge

out vec4 fragColor;   // (mass left around dot i, -, -, -)

void main() {
    int i = int(gl_FragCoord.x);
    // Whole numbers stored as floats; rounded rather than truncated so nothing rides on exactness.
    ivec2 c = ivec2(floor(texelFetch(uSites, ivec2(i, 0), 0).xy + 0.5));

    float sum = 0.0;
    for (int dy = -uRadius; dy <= uRadius; dy++) {
        int y = c.y + dy;
        if (y < 0) y += uSize.y; else if (y >= uSize.y) y -= uSize.y;
        for (int dx = -uRadius; dx <= uRadius; dx++) {
            int x = c.x + dx;
            if (x < 0) x += uSize.x; else if (x >= uSize.x) x -= uSize.x;
            sum += texelFetch(uState, ivec2(x, y), 0).r;
        }
    }

    fragColor = vec4(sum, 0.0, 0.0, 0.0);
}
