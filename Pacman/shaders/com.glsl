#version 300 es
// Stage 2 of the CoM reduction: folds reduce.glsl's 16x16 partials into one texel with the
// soliton's toroidal centre of mass and total mass (same formula as the CPU's computeCoM()).
// Written to a texture rather than read back here, so crop.glsl can queue off it in the same
// batch and the CPU collects CoM + crop together in one stall.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uRed0;
uniform sampler2D uRed1;
uniform ivec2 uSize;    // board (W, H)

out vec4 fragColor;     // (comRow, comCol, mass, valid)

const float TWO_PI = 6.283185307179586;

void main() {
    float mass = 0.0, sy = 0.0, cy = 0.0, sx = 0.0, cx = 0.0;
    for (int j = 0; j < 16; j++) {
        for (int i = 0; i < 16; i++) {
            vec4 a = texelFetch(uRed0, ivec2(i, j), 0);
            mass += a.x;  sy += a.y;  cy += a.z;  sx += a.w;
            cx += texelFetch(uRed1, ivec2(i, j), 0).r;
        }
    }

    // Empty board: no soliton to report. JS reads .w=0 as a dissolve.
    if (mass < 1e-6) {
        fragColor = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }

    float H = float(uSize.y), W = float(uSize.x);
    float twoPiH = TWO_PI / H, twoPiW = TWO_PI / W;

    float my = atan(sy, cy) / twoPiH;  if (my < 0.0) my += H;
    float mx = atan(sx, cx) / twoPiW;  if (mx < 0.0) mx += W;

    fragColor = vec4(my, mx, mass, 1.0);
}
