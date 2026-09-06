#version 300 es
// Stage 2 of the eaten-mass reduction: fold eatreduce.glsl's 16x16 partial sums into
// one texel, the same shape as com.glsl folding reduce.glsl's output. The result only
// has to answer "was anything eaten this step", so the CPU side treats it as a
// boolean threshold rather than a precise mass count.
precision highp float;
precision highp sampler2D;

uniform sampler2D uRed;

out vec4 fragColor;

void main() {
    float sum = 0.0;
    for (int j = 0; j < 16; j++)
        for (int i = 0; i < 16; i++)
            sum += texelFetch(uRed, ivec2(i, j), 0).r;

    fragColor = vec4(sum, 0.0, 0.0, 0.0);
}
