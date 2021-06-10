// from https://www.shadertoy.com/view/7lsGDr
// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

// multiple species, -1=none, 0=specie 1, 1=species 2...
// maximum 16 kernels by using 4x4 matrix
// when matrix operation not available (e.g. exp, mod, equal, /), split into four vec4 operations

#define EPSILON 0.000001
#define mult matrixCompMult
#define useHigherDigitsForSpecies 1
#define testRandomTexel 0

const int speciesNum = 2;  // number of coexist species
const float maxR = 12.;  // space resolution = kernel radius
const float samplingDist = 1.;
// change samplingDist to other numbers (with nearest or linear filter in Buffer A) for funny effects :)
// 1:normal, int>1:heavy phantom, 
// 0.1-0.2:dots, 0.3-0.9:smooth zoom out, 1.1-1.8,2.2-2.8:smooth zoom in, 
// 1.9,2.1,2.9,3.1,3.9(near int):partial phantom, >=3.2:minor glitch, increase as larger
// linear filter: smoother, nearest filter: more glitch/phantom
const vec3 aliveThreshold = vec3(-0.7);

struct genome {
    float R;  // space resolution = kernel radius
    float T;  // time resolution = number of divisions per unit time
    mat4 betaLen;  // kernel ring number
    mat4 beta0;  // kernel ring heights
    mat4 beta1;
    mat4 beta2;
    mat4 mu;  // growth center
    mat4 sigma;  // growth width
    mat4 eta;  // growth strength
    mat4 relR;  // relative kernel radius
    float baseNoise;
    float randomScale;
};
uniform genome genomes[speciesNum];

const ivec3 ic1 = ivec3(1);
const ivec4 iv0 = ivec4(0), iv1 = ivec4(1), iv2 = ivec4(2), iv3 = ivec4(3);
const vec4 v0 = vec4(0.), v1 = vec4(1.);
const mat4 m0 = mat4(v0, v0, v0, v0), m1 = mat4(v1, v1, v1, v1);

// kernel params
const vec4 kmv = vec4(0.5);    // kernel ring center
const mat4 kmu = mat4(kmv, kmv, kmv, kmv);
const vec4 ksv = vec4(0.15);    // kernel ring width
const mat4 ksigma = mat4(ksv, ksv, ksv, ksv);

// source/destination channels
const mat4 src = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4 dst = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels
//const mat4 src = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 1., 1., 2., 2., 0., v0 );  // source channels
//const mat4 dst = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 0., 2., 1., 0., 2., v0 );  // destination channels

const ivec4 src0 = ivec4(src[0]), src1 = ivec4(src[1]), src2 = ivec4(src[2]), src3 = ivec4(src[3]);
const ivec4 dst0 = ivec4(dst[0]), dst1 = ivec4(dst[1]), dst2 = ivec4(dst[2]), dst3 = ivec4(dst[3]);

const mat4[3] srcK = mat4[] (
    mat4( vec4(equal(src0, iv0)), vec4(equal(src1, iv0)), vec4(equal(src2, iv0)), vec4(equal(src3, iv0)) ),
    mat4( vec4(equal(src0, iv1)), vec4(equal(src1, iv1)), vec4(equal(src2, iv1)), vec4(equal(src3, iv1)) ),
    mat4( vec4(equal(src0, iv2)), vec4(equal(src1, iv2)), vec4(equal(src2, iv2)), vec4(equal(src3, iv2)) ) );
const mat4[3] dstK = mat4[] (
    mat4( vec4(equal(dst0, iv0)), vec4(equal(dst1, iv0)), vec4(equal(dst2, iv0)), vec4(equal(dst3, iv0)) ),
    mat4( vec4(equal(dst0, iv1)), vec4(equal(dst1, iv1)), vec4(equal(dst2, iv1)), vec4(equal(dst3, iv1)) ),
    mat4( vec4(equal(dst0, iv2)), vec4(equal(dst1, iv2)), vec4(equal(dst2, iv2)), vec4(equal(dst3, iv2)) ) );

#if useHigherDigitsForSpecies == 1
// higher digits = species, lower digits = value
const float highSize = 8.;  // 2 bits species = none + max 3 species
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(floor(texel * highSize)) - 1;
}
vec3 unpackValue(in vec3 texel) {
    return (fract(texel * highSize) - 0.1) / 0.8;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (vec3(species + 1) + val * 0.8 + 0.1) / highSize;
}
#else
// higher digits = value, lower digits = species
const float highSize = 64.;  // 6 bits value
const float lowSize = 4.;  // 2 bits species = none + max 3 species
const float valMargin = 0.01;  // value 0..1 pack into 0.01..0.99
const float valBand = 1. - 2. * valMargin;
ivec3 unpackSpecies(in vec3 texel) {
    return ivec3(fract(texel * highSize) * lowSize + 0.5) - 1;
}
vec3 unpackValue(in vec3 texel) {
    return (floor(texel * highSize) / highSize -valMargin)/valBand;
}
vec3 packTexel(in ivec3 species, in vec3 val) {
    return (floor((val*valBand+valMargin) * highSize) + vec3(species + 1) / lowSize) / highSize;
}
#endif

// bell-shaped curve (Gaussian bump)
mat4 bell(in mat4 x, in mat4 m, in mat4 s) {
    mat4 v = -mult(x-m, x-m) / s / s / 2.;
    return mat4( exp(v[0]), exp(v[1]), exp(v[2]), exp(v[3]) );
}

// Noise simplex 2D by iq - https://www.shadertoy.com/view/Msf3WH

vec2 hash( vec2 p ) {
    p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
    return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float simplexNoise( in vec2 p ) {
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

// Voronoi noise from https://www.shadertoy.com/view/XlB3zW

float voronoiNoise(in vec2 n, in float s) {
    float i = 0.0, dis = 2.0;
    for(int x = -1;x<=1;x++)
    for(int y = -1;y<=1;y++) {
        vec2 p = floor(n/s)+vec2(x,y);
        vec2 h2 = fract(cos(p*mat2(89.4,-75.7,-81.9,79.6))*343.42);
        float d = length(h2+vec2(x,y)-fract(n/s));
        if (dis > d) { dis = d; i = fract(cos(p.x*89.42-p.y*75.7)*343.42); }
    }
    return i;
}

vec3 randomTexel(in vec2 p) {
    const float pNoSpecies = 0.1;
    vec2 p1 = p / maxR / samplingDist;
    float noise = voronoiNoise(p1 + iDate.w * vec2(10., 10.), 8.);
    int species0 = int( noise * (float(speciesNum) + pNoSpecies) + 0.99 - pNoSpecies ) - 1;
    ivec3 species = ivec3(species0);
    vec3 value;
    if (species0 > -1) {
        genome g = genomes[species0];
        vec2 p2 = p / g.R / samplingDist / g.randomScale;
        vec3 valueNoise = vec3( 
            simplexNoise(p2 + sin(iDate.w) * vec2(10.1, 10.2)), 
            simplexNoise(p2 + sin(iDate.w) * vec2(10.3, -10.4)), 
            simplexNoise(p2 + sin(iDate.w) * vec2(-10.5, 0.6)) );
        value = clamp(g.baseNoise + 0.01 + valueNoise, 0., 1.);
    } else {
        value = vec3(0., 0., 0.);
    }
    return packTexel(species, value);
}

// get neighbor weights for given radius
mat4 getWeight(in float r, in genome g) {
    if (r > g.R) return m0;
    mat4 Br = g.betaLen / g.relR * r / g.R;  // scale radius by number of rings and relative radius
    ivec4 Br0 = ivec4(Br[0]), Br1 = ivec4(Br[1]), Br2 = ivec4(Br[2]), Br3 = ivec4(Br[3]);

    // get heights of kernel rings
    // for each Br: Br==0 ? beta0 : Br==1 ? beta1 : Br==2 ? beta2 : 0
    mat4 height = mat4(
        g.beta0[0] * vec4(equal(Br0, iv0)) + g.beta1[0] * vec4(equal(Br0, iv1)),// + g.beta2[0] * vec4(equal(Br0, iv2)),
        g.beta0[1] * vec4(equal(Br1, iv0)) + g.beta1[1] * vec4(equal(Br1, iv1)),// + g.beta2[1] * vec4(equal(Br1, iv2)),
        g.beta0[2] * vec4(equal(Br2, iv0)) + g.beta1[2] * vec4(equal(Br2, iv1)),// + g.beta2[2] * vec4(equal(Br2, iv2)),
        g.beta0[3] * vec4(equal(Br3, iv0)) + g.beta1[3] * vec4(equal(Br3, iv1)) );// + g.beta2[3] * vec4(equal(Br3, iv2)) );
    mat4 mod1 = mat4( fract(Br[0]), fract(Br[1]), fract(Br[2]), fract(Br[3]) );

    return mult(height, bell(mod1, kmu, ksigma));
}

// draw the shape of kernels
vec3 drawKernels(in vec2 uv, in genome g) {
    ivec2 ij = ivec2(uv / 0.25);  // divide screen into 4x4 grid: 0..3 x 0..3
    vec2 xy = mod(uv, 0.25) * 8. - 1.;  // coordinates in each grid cell: -1..1 x -1..1
    if (ij.x > 3 || ij.y > 3) return vec3(0.);
    float r = length(xy) * maxR;
    mat4 weight = getWeight(r, g);
    vec3 rgb = vec3(weight[3-ij.y][ij.x]);
    return rgb;
}

// map values from source channels to kernels
mat4 mapKernels(in vec3 v) {
    // for each src: src==0 ? r : src==1 ? g : src==2 ? b
    return 
        v.r * srcK[0] + 
        v.g * srcK[1] + 
        v.b * srcK[2];
}

float matrixDot(in mat4 m1, in mat4 m2) {
    return 
        dot(m1[0], m2[0]) + 
        dot(m1[1], m2[1]) + 
        dot(m1[2], m2[2]) + 
        dot(m1[3], m2[3]);
}

// reduce values from kernels to destination channels
vec3 reduceKernels(in mat4 m) {
    // sum of: dst==0 ? r : dst==1 ? g : dst==2 ? b
    return vec3( 
        matrixDot(m, dstK[0]),
        matrixDot(m, dstK[1]),
        matrixDot(m, dstK[2]) );
}

vec3 addSum(in vec2 p, in vec2 d, in mat4[speciesNum] weight, inout mat4[speciesNum] sum, inout mat4[speciesNum] total) {
    vec2 uv = (p + d * samplingDist) / iResolution.xy;  // either set texture repeat or use fract(...)
    vec3 texel = texture(iChannel0, uv).rgb;
    ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    for (int s=0; s<speciesNum; s++) {
        vec3 valueS = value * vec3(equal(species, ivec3(s)));
        mat4 valueK = mapKernels(valueS);
        sum[s] += mult(valueK, weight[s]);
        total[s] += weight[s];
    }
    /*for (int c=0; c<3; c++) {
        int s = species[c];
        if (s > -1) {
            sum[s] += value[c] * mult(srcK[c], weight[s]);
            total[s] += weight[s];
        }
    }*/    
    return texel;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = fragCoord / iResolution.xy;

    #if testRandomTexel == 0

    // loop through the neighborhood, optimized: same weights for all quadrants/octants
    // calculate the weighted average of neighborhood from source channel
    mat4[speciesNum] sum, total, weight;
    for (int s=0; s<speciesNum; s++) {
        sum[s] = mat4(0.);
        total[s] = mat4(0.);
    }
    // self
    float r = 0.;
    for (int s=0; s<speciesNum; s++) 
        weight[s] = getWeight(r, genomes[s]);
    vec3 texel = addSum(fragCoord, vec2(0, 0), weight, sum, total);
    // orthogonal
    const int intR = int(ceil(maxR));
    for (int x=1; x<=intR; x++) {
        r = float(x);
        for (int s=0; s<speciesNum; s++) 
            weight[s] = getWeight(r, genomes[s]);
        addSum(fragCoord, vec2(+x, 0), weight, sum, total);
        addSum(fragCoord, vec2(-x, 0), weight, sum, total);
        addSum(fragCoord, vec2(0, +x), weight, sum, total);
        addSum(fragCoord, vec2(0, -x), weight, sum, total);
    }
    // diagonal
    const int diagR = int(ceil(float(intR) / sqrt(2.)));
    for (int x=1; x<=diagR; x++) {
        r = sqrt(2.) * float(x);
        for (int s=0; s<speciesNum; s++) 
            weight[s] = getWeight(r, genomes[s]);
        addSum(fragCoord, vec2(+x, +x), weight, sum, total);
        addSum(fragCoord, vec2(+x, -x), weight, sum, total);
        addSum(fragCoord, vec2(-x, +x), weight, sum, total);
        addSum(fragCoord, vec2(-x, -x), weight, sum, total);
    }
    // others
    for (int y=1; y<=intR-1; y++)
    for (int x=y+1; x<=intR; x++) {
        r = sqrt(float(x*x + y*y));
        if (r <= maxR) {
            for (int s=0; s<speciesNum; s++) 
                weight[s] = getWeight(r, genomes[s]);
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

    // unpack current cell
    //ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    vec3[speciesNum] growth;
    vec3 maxGrowth = vec3(-99.);
    ivec3 argmax = ivec3(-99);
    for (int s=0; s<speciesNum; s++) {
        genome g = genomes[s];

        // weighted average = weighted sum / total weights, avoid divided by zero
        mat4 avg = sum[s] / (total[s] + EPSILON);    // avoid divided by zero

        // calculate growth (scaled by time step), reduce from kernels to destination channels
        mat4 growthK = mult(g.eta, bell(avg, g.mu, g.sigma) * 2. - 1.) / g.T;
        vec3 growthS = reduceKernels(growthK);
        growth[s] = growthS;

        // calculate max and argmax of growth
        ivec3 isMore = ivec3(greaterThan(growthS, maxGrowth));
        maxGrowth = maxGrowth * float(ic1-isMore) + growthS * float(isMore);
        argmax = argmax * (ic1-isMore) + s * isMore;
    }

    vec3 finalGrowth = vec3(-0.1);
    ivec3 finalSpecies = ivec3(-1);
    for (int s=0; s<speciesNum; s++) {
        // decide which species to be next, by maximum growth (or other criteria?)
        vec3 growthS = growth[s];
        ivec3 isSelect = ivec3(greaterThan(growthS, aliveThreshold)) * ivec3(equal(ivec3(s), argmax));

        finalGrowth = finalGrowth * float(ic1-isSelect) + growthS * float(isSelect);
        finalSpecies = finalSpecies * (ic1-isSelect) + s * isSelect;
    }

    value = clamp(finalGrowth + value, 0., 1.);
    ivec3 species = finalSpecies;

    /*
    mat4 avg = sum[0] / (total[0] + EPSILON);    // avoid divided by zero

	ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    genome g = genomes[0];
    // calculate growth, add a small portion to destination channel
    mat4 growthK = mult(g.eta, bell(avg, g.mu, g.sigma) * 2. - 1.) / g.T;
    vec3 growth = reduceKernels(growthK);
    value = clamp(growth + value, 0., 1.);
    */


    // debug: uncomment to show list of kernels
    //rgb = drawKernels(fragCoord / iResolution.y, genomes[0]);

    // randomize at start
    if (iFrame == 0)
        texel = randomTexel(fragCoord);
    else
        texel = packTexel(species, value);

    #else
    // test randomize
    vec3 texel = randomTexel(fragCoord);
    #endif

    fragColor = vec4(texel, 1.);
}
