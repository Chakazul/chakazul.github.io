#version 300 es
// Stage 2: folds eatreduce.glsl's partials into one texel. Both totals are read on the CPU as
// booleans ("was anything/any pellet eaten"), not precise counts.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uRed;

out vec4 fragColor;   // (total eaten, power-pellet eaten, -, -)

void main() {
    float sum = 0.0, pelletSum = 0.0;
    for (int j = 0; j < 16; j++)
        for (int i = 0; i < 16; i++) {
            vec4 a = texelFetch(uRed, ivec2(i, j), 0);
            sum += a.x; pelletSum += a.y;
        }

    fragColor = vec4(sum, pelletSum, 0.0, 0.0);
}
