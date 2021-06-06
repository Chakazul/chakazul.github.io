// from https://smoothstep.io/anim/26d46785c1f6
// iPrevFrame is a texture with the contents of the last rendered frame.
// We can create stateful shaders by using its contents.

int getCell(in ivec2 p) {
  // Get the value of the cell in the previous frame at position p. 0 or 1.
  return (texelFetch(iPrevFrame, p, 0 ).x > 0.1 ) ? 1 : 0;
}

float getCellAge(in ivec2 p) {
  // Get the y value of the cell, used to simulate the age of the cell, for coloring.
  return texelFetch(iPrevFrame, p, 0 ).y;
}

int getCellWrapped(in ivec2 p) {
  // getCell but with texture wrapping.
  ivec2 r = ivec2(textureSize(iPrevFrame, 0));
  p = (p + r) % r;
  return (texelFetch(iPrevFrame, p, 0 ).x > 0.1 ) ? 1 : 0;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  ivec2 px = ivec2(fragCoord);
  float newCell = 0.;
  float age = 0.;

  ivec3 c = ivec3(1, -1, 0);


  if (px.y > int(iResolution.y) / 4) {
    // Conway's Game of Life.
    int numNeighbours =   (
      getCell(px+c.yy) + getCell(px+c.zy) + getCell(px+c.xy)
      + getCell(px+c.yz)				    + getCell(px+c.xz)
      + getCell(px+c.yx) + getCell(px+c.zx) + getCell(px+c.xx)
    );
    int thisCell = getCell(px);
    newCell = (
      ((numNeighbours == 2) && (thisCell == 1)) || (numNeighbours == 3)
    ) ? 1.0 : 0.0;
    age = getCellAge(px);
    // Decay the age from 1.
    age = newCell + 0.97 * (1. - newCell) * age;
  } else if (px.y != 0) {
    // Copy the row below.
    newCell = float(getCell(px + c.zy));
    age = newCell * 1.0 / (0.8 + 0.01 * fragCoord.y);
  } else {
    // Rule 30 automaton.
    int left = getCellWrapped(px + c.yz);
    int mid = getCellWrapped(px + c.zz);
    int right = getCellWrapped(px + c.xz);
    int val = 4 * left + 2 * mid + right;

    if (1 <= val && val <= 4) {
      newCell = 1.;
    }
  }

  // Initialise.
  if (iFrame < 2) {
    newCell = 0.;
    age = 0.;
    if (px.y == 0 && px.x == int(iResolution.x) / 2) {
      newCell = 1.;
    }
  }
  
    
	float blue = 0.5 - 0.5 * cos(2.0 * newCell + 3.0 * age);
	fragColor = vec4(vec3(newCell, age, blue), 1.);
    
}