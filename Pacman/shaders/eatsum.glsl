#version 300 es
// Stage 2 of the eaten-mass reduction: fold eatreduce.glsl's 16x16 partial sums into
// one texel, the same shape as com.glsl folding reduce.glsl's output. The total (.r)
// only has to answer "was anything eaten this step", so the CPU side treats it as a
// boolean threshold rather than a precise mass count; the pellet-only total (.g) is
// read the same way, to decide whether to (re)start the frightened window.
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
