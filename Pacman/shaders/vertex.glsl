#version 300 es
// Shared vertex shader: every pass is a full-target quad. Each pass reads gl_FragCoord rather
// than a varying, since targets differ in size and all of them want integer texel addressing.

in vec2 a_position;

void main() {
    gl_Position = vec4(a_position, 0.0, 1.0);
}
