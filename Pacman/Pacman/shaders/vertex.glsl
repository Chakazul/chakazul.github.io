#version 300 es
// Shared vertex shader: every pass is a full-target quad, so the vertex stage never
// does anything except pass the two triangles through. Each pass decides what it is
// looking at from gl_FragCoord, not from varyings -- targets differ in size (board,
// 16x16 reduction, 1x1 CoM, 96x96 crop) and integer texel addressing is what all of
// them want anyway.

in vec2 a_position;

void main() {
    gl_Position = vec4(a_position, 0.0, 1.0);
}
