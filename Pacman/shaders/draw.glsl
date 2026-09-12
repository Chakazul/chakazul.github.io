#version 300 es
// Display pass: state ramped between configurable colours, walls on top. Replaces the CPU
// demo's per-pixel putImageData loop. Row 0 is drawn at the canvas top, matching the overlay
// canvas and mouse handlers.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uState;
uniform sampler2D uState2;   // second Lenia channel (dots + power pellets), stepped but never steered
uniform sampler2D uState3;   // ghosts composited to board space by ghostblit.glsl:
                             // R = summed mass, G = 1-based index of the owning ghost
uniform sampler2D uWall;
uniform sampler2D uPower;    // 1 where a power pellet was stamped rather than a plain dot -- same
                             // mass, same rule, coloured differently (see glsim.js's setPowerMask())
uniform ivec2 uSize;         // board (W, H)
uniform vec3  uWallRGB;      // 0..1, flat fill for wall cells
uniform vec3  uSolitonRGB;   // 0..1, the colour mass 1 reaches
uniform vec3  uSoliton2RGB;  // 0..1, the colour an ordinary dot reaches
uniform vec3  uPelletRGB;    // 0..1, the colour a power pellet reaches
// One colour per ghost. Which one applies is read straight out of the composited texture rather
// than re-derived from ghost positions here, so it stays exact even where two tiles overlap.
#define MAX_GHOSTS 9   // must match glsim.js's MAX_GHOSTS
uniform vec3  uGhostRGB[MAX_GHOSTS];
uniform vec3  uSoliton3RGB;  // 0..1, fallback for mass with no owner recorded
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

    // Mass ramps background -> soliton colour; the free-running channels layer over that by their
    // own mass. Ghosts composite last, over Pac-Man, since they're the channel that erases him.
    float v  = clamp(texelFetch(uState,  p, 0).r, 0.0, 1.0);
    float v2 = clamp(texelFetch(uState2, p, 0).r, 0.0, 1.0);
    vec2 g3 = texelFetch(uState3, p, 0).rg;
    float v3 = clamp(g3.r, 0.0, 1.0);
    vec3 c = mix(uBackRGB, uSolitonRGB, v);
    vec3 dotColor = texelFetch(uPower, p, 0).r > 0.5 ? uPelletRGB : uSoliton2RGB;
    c = mix(c, dotColor, v2);

    int owner = int(g3.g) - 1;
    vec3 ghost = (owner >= 0 && owner < MAX_GHOSTS) ? uGhostRGB[owner] : uSoliton3RGB;
    fragColor = vec4(mix(c, ghost, v3), 1.0);
}
