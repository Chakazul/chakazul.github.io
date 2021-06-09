#version 300 es
// Image Drawing Fragment shader template ("Image" tab in Shadertoy)
precision highp float;

#define speciesNum 3

uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform sampler2D iChannel0;             // input channel 0

out vec4 fragColor;

// high precision = species, low precision = value
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
