#version 300 es
// Image Drawing Fragment shader template ("Image" tab in Shadertoy)
precision highp float;

uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform sampler2D iChannel0;             // input channel 0

out vec4 fragColor;

void main() {
    fragColor = texelFetch(iChannel0, ivec2(gl_FragCoord.xy), 0);
}
