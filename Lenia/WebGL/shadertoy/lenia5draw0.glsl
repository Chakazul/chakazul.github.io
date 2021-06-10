#version 300 es
// Image Drawing Fragment shader template ("Image" tab in Shadertoy)
precision highp float;

uniform sampler2D iChannel0;
out vec4 fragColor;

// high precision = species, low precision = value
const float highSize = 4.;  // 2 bits species = none + max 3 species
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(floor(texel * highSize));
}
vec3 unpackValue(in vec3 texel) {
    return (fract(texel * highSize) - 0.1) / 0.8;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (vec3(species) + val * 0.8 + 0.1) / highSize;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    fragColor = vec4(unpackValue(texelFetch(iChannel0, ivec2(fragCoord.xy), 0).rgb), 1.);
}

void main() {
    mainImage(fragColor, gl_FragCoord.xy);
}
