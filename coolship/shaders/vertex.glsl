#version 300 es
// Shared vertex shader: every pass is a full-target quad and addresses texels by gl_FragCoord.

in vec2 a_position;

void main() {
    gl_Position = vec4(a_position, 0.0, 1.0);
}
