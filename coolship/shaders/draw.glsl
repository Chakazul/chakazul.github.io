#version 300 es
// Renders the whole canvas in one pass: the world (left 2/3, the full torus at uWorldScale px
// per cell) and the zoom view (a square in the right panel, centred on the focused creature's CoM,
// nearest-neighbour so cells show as big pixels, optionally rotated heading-up by uView).
// Everything else is black; the 2D overlay canvas draws the controls, arrow and mask outline.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uState;
uniform vec3  uFocus;       // (col, row, on): the focused creature's CoM; off = blank zoom view
uniform ivec2 uSize;        // board (W, H)
uniform vec2  uCanvas;      // canvas (w, h), px
uniform float uWorldScale;  // px per cell in the world view
uniform vec4  uZoom;        // (x0, y0, edge, px per cell), top-left origin
uniform mat2  uView;        // view vector (right, down) -> board vector (col, row)

out vec4 fragColor;

// Integer wrap. Float mod() is x - y*floor(x/y), and with an approximate reciprocal it can return
// y itself for a cell exactly at the edge: texelFetch then reads past the texture and returns
// black, a dark line along the board edge in the zoom view. Shifting by a multiple of the board
// first keeps the % operands non-negative (negative % is undefined in GLSL ES 3.00).
vec3 cellColor(vec2 cell) {
    ivec2 q = (ivec2(floor(cell)) + uSize * 8) % uSize;
    return texelFetch(uState, q, 0).rgb;
}

void main() {
    vec2 s = vec2(gl_FragCoord.x, uCanvas.y - gl_FragCoord.y);   // top-left origin
    vec3 col = vec3(0.0);

    float worldW = float(uSize.x) * uWorldScale;
    if (s.x < worldW) {
        col = cellColor(s / uWorldScale);
    } else if (s.x >= uZoom.x && s.x < uZoom.x + uZoom.z && s.y >= uZoom.y && s.y < uZoom.y + uZoom.z) {
        if (uFocus.z > 0.5) {
            vec2 v = (s - uZoom.xy - 0.5 * uZoom.z) / uZoom.w;   // cells, view frame
            vec2 b = uView * v;                                   // cells, board frame (col, row)
            col = cellColor(uFocus.xy + b + 0.5);
        }
    }
    fragColor = vec4(col, 1.0);
}
