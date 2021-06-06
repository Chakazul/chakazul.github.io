// from https://smoothstep.io/anim/68c8e9a141e5
/*
transition of 3 CAs: Wolfram > Conway > Lenia
  1st stage: Wolfram's elementary cellular automata, rule 30
  2nd stage: Conway's Game of Life
  3rd stage: Lenia (multi-kernel version)
based on example - https://smoothstep.io/anim/26d46785c1f6
https://chakazul.github.io/lenia.html
*/

const float R = 10.;       // space resolution = kernel radius
const float T = 10.;       // time resolution = number of divisions per unit time
const float dt = 1./T;     // time step
const float ageRatio = 0.4;

// choose a species by uncommenting block (add/remove star at top)

/*/
// species: glider
const vec3 B = vec3( 3., 1., 1. );    // kernel ring number
const mat3 b = mat3( 0.66, 1., 1.,  1., 0., 0.,  0., 0., 0. );  // kernel ring heights
const vec3 m = vec3( 0.3, 0.23, 0.19 );    // growth center
const vec3 s = vec3( 0.03, 0.038, 0.027 );    // growth width
const vec3 h = vec3( 0.666, 0.666, 0.666 );     // growth strength
/**/

/**/
// species: ghost
const vec3 B = vec3( 2., 3., 1. );    // kernel ring number
const mat3 b = mat3( 0.25, 1., 1.,  1., 0.75, 0.,  0., 0.75, 0. );  // kernel ring heights
const vec3 m = vec3( 0.16, 0.22, 0.28 );    // growth center
const vec3 s = vec3( 0.025, 0.042, 0.025 );    // growth width
const vec3 h = vec3( 0.666, 0.666, 0.666 );     // growth strength
/**/

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

vec3 bell(vec3 x, vec3 m, vec3 s) {
    return exp(-(x-m)*(x-m)/s/s/2.);  // bell-shaped curve
}

vec3 getWeight(in float r) {
    if (r > 1.) return vec3(0.);
    vec3 Br = B * r;
    ivec3 iBr = ivec3(Br);
    vec3 height = 
      b[0] * vec3(equal(iBr, ivec3(0))) + 
      b[1] * vec3(equal(iBr, ivec3(1))) + 
      b[2] * vec3(equal(iBr, ivec3(2)));
    return height * bell(mod(Br, 1.), vec3(0.5), vec3(0.15));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  vec2 uv = fragCoord / iResolution;
  ivec2 px = ivec2(fragCoord);
  float newCell = 0.;
  float age = 0.;

  ivec3 c = ivec3(1, -1, 0);

  if (px.y > int(iResolution.y) / 2) {
    // Lenia (multi-kernel)
    vec3 sum = vec3(0.), total = vec3(0.);
    for (int x=-int(R); x<=int(R); x++)
    for (int y=-int(R); y<=int(R); y++)
    {
        float r = sqrt(float(x*x + y*y)) / R;
        vec3 weight = getWeight(r); //bell(r, 0.5, 0.15);
        float cell = getCellAge(px + ivec2(x,y)) / ageRatio;
        sum += cell * weight;
        total += weight;
    }
    vec3 avg = sum / total;
    vec3 growth = bell(avg, m, s) * 2. - 1.;
    float weightedGrowth = dot(h, growth);
    newCell = clamp(getCellAge(px) / ageRatio + dt * weightedGrowth, 0., 1.);
    age = newCell * ageRatio;

  } else if (px.y > int(iResolution.y) / 4) {
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
    //newCell = fract(103.2 * sin(5102.2 * uv.x + 983.87 * uv.y + 23.1) * (7. * uv.x + 923.2 * uv.y));
    age = 0.;
    if (px.y == 0 && px.x == int(iResolution.x) / 2) {
      newCell = 1.;
    }
  }
  
    
	float blue = 0.5 - 0.5 * cos(2.0 * newCell + 3.0 * age);
	fragColor = vec4(vec3(newCell, age, blue), 1.);
    
}
