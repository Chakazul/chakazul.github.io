#version 300 es
// Display pass: the state ramped between two configurable colours, with walls over the top.
// Replaces the per-pixel JS loop + putImageData in the CPU demo's render().
//
// Row 0 of the board is drawn at the top of the canvas, matching the ImageData
// orientation the overlay canvas (arrows, action discs) and the mouse handlers assume.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uState;
uniform sampler2D uState2;   // second Lenia channel, stepped but never steered
uniform sampler2D uWall;
uniform ivec2 uSize;         // board (W, H)
uniform vec3  uWallRGB;      // 0..1, flat fill for wall cells
uniform vec3  uSolitonRGB;   // 0..1, the colour mass 1 reaches
uniform vec3  uSoliton2RGB;  // 0..1, the colour channel 2 reaches
uniform vec3  uBackRGB;      // 0..1, the colour mass 0 sits at

out vec4 fragColor;

void main() {
    int col = int(gl_FragCoord.x);
    int row = uSize.y - 1 - int(gl_FragCoord.y);   // GL y is up, board row 0 is at top
    ivec2 p = ivec2(col, row);

    if (texelFetch(uWall, p, 0).r > 0.5) {
        fragColor = vec4(uWallRGB, 1.0);
        return;
    }

    // Mass ramps the board from the background colour up to the soliton colour, replacing the
    // plain black->white ramp the CPU demo used. Channel 2 is laid over that ramp by its own
    // mass, so where the two overlap the denser channel is the one that shows.
    float v  = clamp(texelFetch(uState,  p, 0).r, 0.0, 1.0);
    float v2 = clamp(texelFetch(uState2, p, 0).r, 0.0, 1.0);
    vec3 c = mix(uBackRGB, uSolitonRGB, v);
    fragColor = vec4(mix(c, uSoliton2RGB, v2), 1.0);
}
