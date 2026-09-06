#version 300 es
// Builds the agent's input window: a netSize x netSize toroidal crop centered on the
// soliton, taken from the last K=4 board states and packed one frame per channel.
//
// Two things matter here. First, all four frames are cropped at the *same* origin --
// the one derived from the current CoM -- exactly as buildCroppedStateInput() does. If
// each frame were instead cropped at the CoM it had when it was captured, the stack
// would show a stationary soliton and the velocity information the policy reads would
// be gone. Keeping the four board states on the GPU and cropping them together here is
// what makes that possible without reading whole boards back.
//
// Second, the origin comes from the CoM texture rather than a uniform, so this pass
// does not have to wait for the CPU to learn where the soliton is.
precision highp float;
precision highp sampler2D;
precision highp int;      // fragment-stage int defaults to mediump -- see action.glsl

uniform sampler2D uF0;    // oldest state  (model channel k=0)
uniform sampler2D uF1;
uniform sampler2D uF2;
uniform sampler2D uF3;    // newest state  (model channel k=3)
uniform sampler2D uCom;   // 1x1, (comRow, comCol, mass, valid)
uniform ivec2 uSize;      // board (W, H)
uniform int   uNet;       // window edge

out vec4 fragColor;

void main() {
    ivec2 q = ivec2(gl_FragCoord.xy);          // window cell, 0..uNet-1
    vec4 com = texelFetch(uCom, ivec2(0, 0), 0);

    // floor(x+0.5) is JS Math.round, so this origin matches agentWindowOrigin() on the
    // CPU side -- which the JS uses to map the policy's argmax back to board coords.
    int halfNet = uNet / 2;
    int oy = int(floor(com.x + 0.5)) - halfNet;
    int ox = int(floor(com.y + 0.5)) - halfNet;

    int y = oy + q.y;  y -= uSize.y * int(floor(float(y) / float(uSize.y)));
    int x = ox + q.x;  x -= uSize.x * int(floor(float(x) / float(uSize.x)));
    ivec2 p = ivec2(x, y);

    fragColor = vec4(texelFetch(uF0, p, 0).r,
                     texelFetch(uF1, p, 0).r,
                     texelFetch(uF2, p, 0).r,
                     texelFetch(uF3, p, 0).r);
}
