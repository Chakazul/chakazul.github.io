#version 300 es
// Display pass: grayscale state, wall color, and the goal blob tinted over the top.
// Replaces the per-pixel JS loop + putImageData in the CPU demo's render().
//
// The goal blob is evaluated analytically instead of being kept as a board-sized
// texture the way buildTargetBlob() built one -- it is only ever used for display, and
// this is the same Gaussian, non-toroidal like the original.
//
// Row 0 of the board is drawn at the top of the canvas, matching the ImageData
// orientation the overlay canvas (arrows, action discs) and the mouse handlers assume.
precision highp float;
precision highp sampler2D;

uniform sampler2D uState;
uniform sampler2D uWall;
uniform ivec2 uSize;         // board (W, H)
uniform vec2  uTarget;       // goal (row, col)
uniform float uTargetSigma;
uniform vec3  uWallRGB;      // 0..1
uniform vec3  uTargetRGB;    // 0..1

out vec4 fragColor;

void main() {
    int col = int(gl_FragCoord.x);
    int row = uSize.y - 1 - int(gl_FragCoord.y);   // GL y is up, board row 0 is at top
    ivec2 p = ivec2(col, row);

    if (texelFetch(uWall, p, 0).r > 0.5) {
        fragColor = vec4(uWallRGB, 1.0);
        return;
    }

    vec3 c = vec3(clamp(texelFetch(uState, p, 0).r, 0.0, 1.0));

    float dr = float(row) - uTarget.x;
    float dc = float(col) - uTarget.y;
    float tb = exp(-(dr * dr + dc * dc) / (2.0 * uTargetSigma * uTargetSigma));
    if (tb > 0.02) c = mix(c, uTargetRGB, tb);

    fragColor = vec4(c, 1.0);
}
