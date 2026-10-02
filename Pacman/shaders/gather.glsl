#version 300 es
// Copies every small per-step result into one row, so readback() is a single readPixels instead
// of one per result. Layout (texel x), must match G_* in glsim.js:
//   0                      Pac-Man CoM        (row, col, mass, valid)          com.glsl
//   1                      eaten mass         (total, power pellet, -, -)      eatsum.glsl
//   2                      dots channel CoM   (-, -, total mass, -)            com.glsl
//   3 .. 3+MAX_GHOSTS-1    per-ghost CoM      (local row, local col, mass, valid)  tilecom.glsl
//   3+MAX_GHOSTS ..        per-dot mass       (mass, -, -, -)                  dotsites.glsl
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

#define MAX_GHOSTS 9   // must match glsim.js's MAX_GHOSTS

uniform sampler2D uCom;
uniform sampler2D uEat;
uniform sampler2D uDotsCom;
uniform sampler2D uGhostCom;
uniform sampler2D uSites;

out vec4 fragColor;

void main() {
    int i = int(gl_FragCoord.x);
    if (i == 0)                   fragColor = texelFetch(uCom, ivec2(0, 0), 0);
    else if (i == 1)              fragColor = texelFetch(uEat, ivec2(0, 0), 0);
    else if (i == 2)              fragColor = texelFetch(uDotsCom, ivec2(0, 0), 0);
    else if (i < 3 + MAX_GHOSTS)  fragColor = texelFetch(uGhostCom, ivec2(i - 3, 0), 0);
    else                          fragColor = texelFetch(uSites, ivec2(i - 3 - MAX_GHOSTS, 0), 0);
}
