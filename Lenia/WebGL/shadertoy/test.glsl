// from https://smoothstep.io/anim/534226b1d5ed
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  // Compute coordinates that range from x = 0 to 1 and y = 0 to 1.
  vec2 uv = fragCoord / iResolution.xy;

  // Compute a color that varies with uv.x, uv.y and iTime.
  vec3 col = 0.5 + 0.5 * cos(iTime + uv.xyx + vec3(0, 2, 4));

  // Output the colour for this pixel.
  fragColor = vec4(col, 1.0);
}
