#version 300 es
// Stage 2 of the reduction: fold the 16x16 partial sums into one texel holding the
// soliton's toroidal center of mass and total mass. Same circular-mean formula as
// computeCoM() in the CPU demo.
//
// Writing this to a texture rather than reading the partials back is what keeps the
// step down to a single GPU->CPU stall: crop.glsl reads the CoM straight out of this
// texture, so the crop can be queued in the same batch as the step that produced it,
// and the CPU picks up the CoM and the crop together in one readback.
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

    // Below this the board is empty and there is no soliton left to locate -- the JS
    // side reads .w and treats the episode as a dissolve, matching computeCoM()
    // returning null.
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
