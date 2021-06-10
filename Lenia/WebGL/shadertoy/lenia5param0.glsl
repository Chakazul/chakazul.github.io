// from https://www.shadertoy.com/view/7lsGDr
// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

// multiple species
// maximum 16 kernels by using 4x4 matrix
// when matrix operation not available (e.g. exp, mod, equal, /), split into four vec4 operations

#define EPSILON 0.000001
#define mult matrixCompMult
#define speciesNum 2

// change to other numbers (with nearest or linear filter in Buffer A) for funny effects :)
const float samplingDist = 1.;
// 1:normal, int>1:heavy phantom, 
// 0.1-0.2:dots, 0.3-0.9:smooth zoom out, 1.1-1.8,2.2-2.8:smooth zoom in, 
// 1.9,2.1,2.9,3.1,3.9(near int):partial phantom, >=3.2:minor glitch, increase as larger
// linear filter: smoother, nearest filter: more glitch/phantom

const ivec4 iv0 = ivec4(0), iv1 = ivec4(1), iv2 = ivec4(2), iv3 = ivec4(3);
const vec4 v0 = vec4(0.), v1 = vec4(1.);
const mat4 m0 = mat4(v0, v0, v0, v0), m1 = mat4(v1, v1, v1, v1);

const float R = 12.;  // space resolution = kernel radius
struct genome {
	float T;  // time resolution = number of divisions per unit time
	float baseNoise;
	mat4 betaLen;  // kernel ring number
	mat4 beta0;  // kernel ring heights
	mat4 beta1;
	mat4 beta2;
	mat4 mu;  // growth center
	mat4 sigma;  // growth width
	mat4 eta;  // growth strength
	mat4 relR;  // relative kernel radius
};
uniform genome genomes[speciesNum];

const mat4 src = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4 dst = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels
//const mat4 src = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., v0, v0 );  // source channels
//const mat4 dst = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., v0, v0 );  // destination channels

// precalculate
const vec4 kmv = vec4(0.5);    // kernel ring center
const mat4 kmu = mat4(kmv, kmv, kmv, kmv);
const vec4 ksv = vec4(0.15);    // kernel ring width
const mat4 ksigma = mat4(ksv, ksv, ksv, ksv);

const ivec4 src0 = ivec4(src[0]), src1 = ivec4(src[1]), src2 = ivec4(src[2]), src3 = ivec4(src[3]);
const ivec4 dst0 = ivec4(dst[0]), dst1 = ivec4(dst[1]), dst2 = ivec4(dst[2]), dst3 = ivec4(dst[3]);


// Noise simplex 2D by iq - https://www.shadertoy.com/view/Msf3WH

vec2 hash( vec2 p ) {
	p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p ) {
    const float K1 = 0.366025404; // (sqrt(3)-1)/2;
    const float K2 = 0.211324865; // (3-sqrt(3))/6;

	vec2  i = floor( p + (p.x+p.y)*K1 );
    vec2  a = p - i + (i.x+i.y)*K2;
    float m = step(a.y,a.x); 
    vec2  o = vec2(m,1.0-m);
    vec2  b = a - o + K2;
	vec2  c = a - 1.0 + 2.0*K2;
    vec3  h = max( 0.5-vec3(dot(a,a), dot(b,b), dot(c,c) ), 0.0 );
	vec3  n = h*h*h*h*vec3( dot(a,hash(i+0.0)), dot(b,hash(i+o)), dot(c,hash(i+1.0)));
    return dot( n, vec3(70.0) );
}


/*/
// high precision = value, low precision = species
const float highSize = 64.;  // 6 bits value
const float lowSize = 4.;  // 2 bits species = none + max 3 species
const float valMargin = 0.01;  // value 0..1 pack into 0.01..0.99
const float valBand = 1. - 2. * valMargin;
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(fract(texel * highSize) * lowSize + 0.5);
}
vec3 unpackValue(in vec3 texel) {
    return (floor(texel * highSize) / highSize -valMargin)/valBand;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (floor((val*valBand+valMargin) * highSize) + vec3(species) / lowSize) / highSize;
}
/*/
// high precision = species, low precision = value
const float highSize = 4.;  // 2 bits species = none + max 3 species
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(floor(texel * highSize));
}
vec3 unpackValue(in vec3 texel) {
    return (fract(texel * highSize) - 0.1) / 0.8;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (vec3(species) + val * 0.8 + 0.1) / highSize;
}
/**/

// bell-shaped curve (Gaussian bump)
mat4 bell(in mat4 x, in mat4 m, in mat4 s) {
    mat4 v = -mult(x-m, x-m) / s / s / 2.;
    return mat4( exp(v[0]), exp(v[1]), exp(v[2]), exp(v[3]) );
}

// get neighbor weights for given radius
mat4 getWeight(in float r, in genome g) {
    if (r > 1.) return m0;
    mat4 Br = g.betaLen / g.relR * r;  // scale radius by number of rings and relative radius
    ivec4 Br0 = ivec4(Br[0]), Br1 = ivec4(Br[1]), Br2 = ivec4(Br[2]), Br3 = ivec4(Br[3]);

    // get heights of kernel rings
    // for each Br: Br==0 ? beta0 : Br==1 ? beta1 : Br==2 ? beta2 : 0
    mat4 height = mat4(
        g.beta0[0] * vec4(equal(Br0, iv0)) + g.beta1[0] * vec4(equal(Br0, iv1)),// + g.beta2[0] * vec4(equal(Br0, iv2)),
        g.beta0[1] * vec4(equal(Br1, iv0)) + g.beta1[1] * vec4(equal(Br1, iv1)),// + g.beta2[1] * vec4(equal(Br1, iv2)),
        g.beta0[2] * vec4(equal(Br2, iv0)) + g.beta1[2] * vec4(equal(Br2, iv1)),// + g.beta2[2] * vec4(equal(Br2, iv2)),
        g.beta0[3] * vec4(equal(Br3, iv0)) + g.beta1[3] * vec4(equal(Br3, iv1)) );// + g.beta2[3] * vec4(equal(Br3, iv2)) );
    mat4 mod1 = mat4( mod(Br[0], 1.), mod(Br[1], 1.), mod(Br[2], 1.), mod(Br[3], 1.) );

    return mult(height, bell(mod1, kmu, ksigma));
}

// draw the shape of kernels
vec3 drawKernel(in vec2 uv, in genome g) {
    ivec2 ij = ivec2(uv / 0.25);  // divide screen into 4x4 grid: 0..3 x 0..3
    vec2 xy = mod(uv, 0.25) * 8. - 1.;  // coordinates in each grid cell: -1..1 x -1..1
    if (ij.x > 3 || ij.y > 3) return vec3(0.);
    float r = length(xy);
    mat4 weight = getWeight(r, g);
    vec3 rgb = vec3(weight[3-ij.y][ij.x]);
    return rgb;
}

// map values from source channels to kernels
vec4 mapKernels(in vec3 v, in ivec4 srcv) {
    // for each src: src==0 ? r : src==1 ? g : src==2 ? b : 0
    return
        v.r * vec4(equal(srcv, iv0)) + 
        v.g * vec4(equal(srcv, iv1)) +
        v.b * vec4(equal(srcv, iv2));
}

// reduce values from kernels to destination channels
float reduceKernels(in mat4 m, in ivec4 ch) {
    // (r,g,b) = sum of m where dst==(0,1,2)
    return 
        dot(m[0], vec4(equal(dst0, ch))) + 
        dot(m[1], vec4(equal(dst1, ch))) + 
        dot(m[2], vec4(equal(dst2, ch))) + 
        dot(m[3], vec4(equal(dst3, ch)));
}

// add to weighted sum and total for given cell, return cell
vec3 addSum(in vec2 p, in vec2 d, in mat4[speciesNum] weight, inout mat4[speciesNum] sum, inout mat4[speciesNum] total) {
    vec2 uv = mod((p + d * samplingDist) / iResolution.xy, 1.);
    vec3 texel = texture(iChannel0, uv).rgb;
	vec3 value = unpackValue(texel);
    mat4 valueK = mat4( mapKernels(value, src0), mapKernels(value, src1), mapKernels(value, src2), mapKernels(value, src3) );
	sum[0] += mult(valueK, weight[0]);
	total[0] += weight[0];
	return texel;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = fragCoord / iResolution.xy;

    // loop through the neighborhood, optimized: same weights for all quadrants/octants
    // calculate the weighted average of neighborhood from source channel
    mat4[speciesNum] sum, total, weight, avg;
	sum[0] = mat4(0.);
	total[0] = mat4(0.);
    // self
    float r = 0.;
	weight[0] = getWeight(r, genomes[0]);
    vec3 texel = addSum(fragCoord, vec2(0, 0), weight, sum, total);
    // orthogonal
	const int intR = int(ceil(R));
    for (int x=1; x<=intR; x++) {
        r = float(x) / R;
		weight[0] = getWeight(r, genomes[0]);
        addSum(fragCoord, vec2(+x, 0), weight, sum, total);
        addSum(fragCoord, vec2(-x, 0), weight, sum, total);
        addSum(fragCoord, vec2(0, +x), weight, sum, total);
        addSum(fragCoord, vec2(0, -x), weight, sum, total);
    }
    // diagonal
    const int diagR = int(ceil(float(intR) / sqrt(2.)));
    for (int x=1; x<=diagR; x++) {
        r = sqrt(2.) * float(x) / R;
		weight[0] = getWeight(r, genomes[0]);
        addSum(fragCoord, vec2(+x, +x), weight, sum, total);
        addSum(fragCoord, vec2(+x, -x), weight, sum, total);
        addSum(fragCoord, vec2(-x, +x), weight, sum, total);
        addSum(fragCoord, vec2(-x, -x), weight, sum, total);
    }
    // others
    for (int y=1; y<=intR-1; y++)
    for (int x=y+1; x<=intR; x++) {
        r = sqrt(float(x*x + y*y)) / R;
        if (r <= 1.) {
			weight[0] = getWeight(r, genomes[0]);
            addSum(fragCoord, vec2(+x, +y), weight, sum, total);
            addSum(fragCoord, vec2(+x, -y), weight, sum, total);
            addSum(fragCoord, vec2(-x, +y), weight, sum, total);
            addSum(fragCoord, vec2(-x, -y), weight, sum, total);
            addSum(fragCoord, vec2(+y, +x), weight, sum, total);
            addSum(fragCoord, vec2(+y, -x), weight, sum, total);
            addSum(fragCoord, vec2(-y, +x), weight, sum, total);
            addSum(fragCoord, vec2(-y, -x), weight, sum, total);
        }
    }
    avg[0] = sum[0] / (total[0] + EPSILON);    // avoid divided by zero

	ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    // calculate growth, add a small portion to destination channel
    mat4 growthK = mult(genomes[0].eta, bell(avg[0], genomes[0].mu, genomes[0].sigma) * 2. - 1.) / genomes[0].T;
    vec3 growth = vec3( reduceKernels(growthK, iv0), reduceKernels(growthK, iv1), reduceKernels(growthK, iv2) );
    value = clamp(growth + value, 0., 1.);

    // debug: uncomment to show list of kernels
    //rgb = drawKernel(fragCoord / iResolution.y, genomes[0]);

    // randomize at start, or add patch on mouse click
    if (iFrame == 0 || iMouse.z > 0.) {
		float speciesNoise = clamp((noise(fragCoord/R/samplingDist/8. + fract(iDate.w)*10.) + 1.) / 1.5, 0., 0.99);
        species = ivec3(speciesNoise * float(speciesNum+1) );
        vec3 valueNoise = vec3( 
            noise(fragCoord/R/samplingDist/2. + fract(iDate.w)*100.), 
            noise(fragCoord/R/samplingDist/2. + sin(iDate.w)*100.), 
            noise(fragCoord/R/samplingDist/2. + cos(iDate.w)*100.) );
        value = clamp(genomes[0].baseNoise + valueNoise, 0., 1.);
    }

	texel = packTexel(species, value);
    fragColor = vec4(texel, 1.);
}
