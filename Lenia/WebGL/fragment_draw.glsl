#version 300 es
// Image Drawing Fragment shader template ("Image" tab in Shadertoy)
precision highp float;

uniform sampler2D iChannel0;
out vec4 fragColor;

void main()
{
    fragColor = texelFetch(iChannel0, ivec2(gl_FragCoord.xy), 0);
}
