#version 300 es
// Image Drawing Fragment shader template ("Image" tab in Shadertoy)
precision highp float;

#define useHigherDigitsForSpecies 1

const int speciesNum = 2;

uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform sampler2D iChannel0;             // input channel 0

out vec4 fragColor;

#if useHigherDigitsForSpecies == 1
// higher digits = species, lower digits = value
const float highSize = 8.;  // 2 bits species = none + max 3 species
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(floor(texel * highSize));
}
vec3 unpackValue(in vec3 texel) {
    return (fract(texel * highSize) - 0.1) / 0.8;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (vec3(species) + val * 0.8 + 0.1) / highSize;
}
#else
// higher digits = value, lower digits = species
const float highSize = 64.;  // 6 bits value
const float lowSize = 4.;  // 2 bits species = none + max 3 species
const float valMargin = 0.01;  // value 0..1 pack into 0.01..0.99
const float valBand = 1. - 2. * valMargin;
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(fract(texel * highSize) * lowSize + 0.5);
}
vec3 unpackValue(in vec3 texel) {
    return (floor(texel * highSize) / highSize -valMargin)/valBand;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (floor((val*valBand+valMargin) * highSize) + vec3(species) / lowSize) / highSize;
}
#endif

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec3 texel = texelFetch(iChannel0, ivec2(fragCoord.xy), 0).rgb;
    if (iMouse.z > 0.)
        texel = vec3(unpackSpecies(texel)) / float(speciesNum);
    else
        texel = unpackValue(texel);
    fragColor = vec4(texel, 1.);
}

void main() {
    mainImage(fragColor, gl_FragCoord.xy);
}
