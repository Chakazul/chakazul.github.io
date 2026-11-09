#version 300 es
// Renders the whole canvas in one pass: the world (left 2/3, the full torus at uWorldScale px
// per cell) and the zoom view (a square in the right panel, centred on the focused creature's CoM,
// nearest-neighbour so cells show as big pixels, optionally rotated heading-up by uView).
// Everything else is black; the 2D overlay canvas draws the controls, arrow and mask outline.
//
// Star field (uStars on): three layers of small tinted points that live in board coordinates, so
// both views show the same stars, drifting to the left with a slow up-and-down sway and wrapping round the torus.
// Near layers are brighter and faster, far ones dimmer and slower. A star is a few screen px in
// either view (its size is in px, not cells), and it is drawn only behind the creatures.
precision highp float;
precision highp sampler2D;
precision highp int;

uniform sampler2D uState;
uniform sampler2D uDots;    // dots channel (R), shown in uDotRGB over the creatures where uDotsOn
uniform float uDotsOn;
uniform vec3  uDotRGB;
uniform sampler2D uBullets; // bullets channel (R), shown in uBulletRGB over everything where uBulletsOn
uniform float uBulletsOn;
uniform vec3  uBulletRGB;
uniform float uDotOff;      // the dots' drift to the left, cells: board x shows dot x + uDotOff
uniform vec3  uFocus;       // (col, row, on): the focused creature's CoM; off = blank zoom view
uniform ivec2 uSize;        // board (W, H)
uniform vec2  uCanvas;      // canvas (w, h), px
uniform float uWorldScale;  // px per cell in the world view
uniform vec4  uZoom;        // (x0, y0, edge, px per cell), top-left origin
uniform mat2  uView;        // view vector (right, down) -> board vector (col, row)
uniform float uStars;       // 1 = star field on
uniform float uTime;        // seconds, for the star drift
uniform vec3  uBg;          // the background behind the board (glsim.js setBackground), black by default

out vec4 fragColor;

// Integer wrap. Float mod() is x - y*floor(x/y), and with an approximate reciprocal it can return
// y itself for a cell exactly at the edge: texelFetch then reads past the texture and returns
// black, a dark line along the board edge in the zoom view. Shifting by a multiple of the board
// first keeps the % operands non-negative (negative % is undefined in GLSL ES 3.00).
vec3 cellColor(vec2 cell) {
    ivec2 q = (ivec2(floor(cell)) + uSize * 8) % uSize;
    vec3 c = texelFetch(uState, q, 0).rgb;
    if (uDotsOn > 0.5) {
        // linear across x, so the drift is smooth rather than whole cells at a time
        float x = cell.x + uDotOff - 0.5, f = fract(x);
        int x0 = (int(floor(x)) + uSize.x * 8) % uSize.x, x1 = (x0 + 1) % uSize.x;
        float v = mix(texelFetch(uDots, ivec2(x0, q.y), 0).r, texelFetch(uDots, ivec2(x1, q.y), 0).r, f);
        c = mix(c, uDotRGB, v);
    }
    if (uBulletsOn > 0.5) c = mix(c, uBulletRGB, texelFetch(uBullets, q, 0).r);
    return c;
}

uint hash(uvec3 v) {
    uint h = v.x * 0x8da6b343u ^ v.y * 0xd8163841u ^ v.z * 0xcb1ab31fu;
    h ^= h >> 15; h *= 0x2c1b3c6du; h ^= h >> 12; h *= 0x297a2d39u; h ^= h >> 15;
    return h;
}
float rnd(uvec3 v) { return float(hash(v) & 0xffffffu) / 16777216.0; }

// One layer: a grid of g x g cells (g divides the board, so it wraps), at most one star per grid
// cell, placed away from the cell edge so only this cell needs testing. pxPerCell = view scale.
// Drift: steadily to the left, with a slow sine sway up and down. Both scale with the layer's
// speed, so near layers sway further (parallax). Period 2 pi / STAR_SWAY_W seconds.
const float STAR_SWAY_A = 1.6;    // sway amplitude, cells per unit of layer speed
const float STAR_SWAY_W = 0.3;    // sway angular frequency, rad/s (about 21 s per swing)

vec3 starLayer(vec2 b, float pxPerCell, float g, float speed, float bright, float density, uint layer) {
    vec2 sz = vec2(uSize);
    vec2 shift = vec2(-uTime, STAR_SWAY_A * sin(STAR_SWAY_W * uTime)) * speed;
    vec2 q = mod(b - shift, sz);
    vec2 gi = floor(q / g);
    uvec3 key = uvec3(uvec2(gi), layer);
    if (rnd(key) > density) return vec3(0.0);
    vec2 pos = (gi + 0.15 + 0.7 * vec2(rnd(key + uvec3(0u, 0u, 101u)), rnd(key + uvec3(0u, 0u, 202u)))) * g;
    vec2 d = q - pos;
    d -= sz * round(d / sz);                       // torus
    float dpx = length(d) * pxPerCell;
    float r = 0.8 + 0.85 * bright;                 // radius in screen px: ~1.05, 1.25, 1.6 by layer
    float a = clamp(r + 0.5 - dpx, 0.0, 1.0);
    if (a <= 0.0) return vec3(0.0);
    // mostly white, with a faint cyan, red or yellow tint
    float t = rnd(key + uvec3(0u, 0u, 303u));
    vec3 tint = t < 0.25 ? vec3(0.7, 1.0, 1.0) : t < 0.4 ? vec3(1.0, 0.75, 0.7)
              : t < 0.6 ? vec3(1.0, 0.95, 0.7) : vec3(1.0);
    float level = bright * (0.6 + 0.4 * rnd(key + uvec3(0u, 0u, 404u)));
    return tint * level * a;
}

vec3 stars(vec2 b, float pxPerCell) {
    return starLayer(b, pxPerCell, 8.0, 0.6, 0.30, 0.10, 1u)     // far: dim, slow
         + starLayer(b, pxPerCell, 16.0, 1.5, 0.55, 0.12, 2u)    // middle
         + starLayer(b, pxPerCell, 32.0, 3.5, 0.95, 0.15, 3u);   // near: bright, fast
}

// Stars only show through where the board is (nearly) empty.
// Stars and the background colour (uBg, black by default) only show where the board is (nearly) empty.
vec3 withStars(vec3 col, vec2 b, float pxPerCell) {
    float cover = clamp(max(col.r, max(col.g, col.b)) * 4.0, 0.0, 1.0);
    vec3 under = uBg + (uStars < 0.5 ? vec3(0.0) : stars(b, pxPerCell));
    return col + under * (1.0 - cover);
}

void main() {
    vec2 s = vec2(gl_FragCoord.x, uCanvas.y - gl_FragCoord.y);   // top-left origin
    vec3 col = vec3(0.0);

    float worldW = float(uSize.x) * uWorldScale;
    if (s.x < worldW) {
        col = withStars(cellColor(s / uWorldScale), s / uWorldScale, uWorldScale);
    } else if (s.x >= uZoom.x && s.x < uZoom.x + uZoom.z && s.y >= uZoom.y && s.y < uZoom.y + uZoom.z) {
        if (uFocus.z > 0.5) {
            vec2 v = (s - uZoom.xy - 0.5 * uZoom.z) / uZoom.w;   // cells, view frame
            vec2 b = uView * v;                                   // cells, board frame (col, row)
            col = withStars(cellColor(uFocus.xy + b + 0.5), uFocus.xy + b + 0.5, uZoom.w);
        }
    }
    fragColor = vec4(col, 1.0);
}
