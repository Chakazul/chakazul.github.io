// from https://www.shadertoy.com/view/7lsGDr
// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

// multiple species
// maximum 16 kernels by using 4x4 matrix (mat4)

#define EPSILON 0.000001
#define mult matrixCompMult

#define usePackSpecies 1
#define usePackSpeciesHigher 0
#define useDecideSpecies 1
#define useDecideByMaxGrowth 1
#define use3RingKernels 0

const int speciesNum = 2;
const float maxR = 15.;

const float biasNoise = speciesNum == 1 ? 0. : -0.06;
const float pNoSpecies = -0.1;
const float samplingDist = 1.;
// change samplingDist to other numbers (with nearest or linear filter in Buffer A) for funny effects :)
// 1:normal, int>1:heavy phantom, 
// 0.1-0.2:dots, 0.3-0.9:smooth zoom out, 1.1-1.8,2.2-2.8:smooth zoom in, 
// 1.9,2.1,2.9,3.1,3.9(near int):partial phantom, >=3.2:minor glitch, increase as larger
// linear filter: smoother, nearest filter: more glitch/phantom
const vec3 aliveThreshold = vec3(-0.7);

const ivec3 ic1 = ivec3(1);
const ivec4 iv0 = ivec4(0), iv1 = ivec4(1), iv2 = ivec4(2), iv3 = ivec4(3);
const vec4 v0 = vec4(0.), v1 = vec4(1.);
const mat4 m0 = mat4(0.), m1 = mat4(v1, v1, v1, v1);

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

const genome genomes[speciesNum] = genome[] (
    genome(
        15., 2.,
        mat4( 1., 1., 1., 2., 1., 2., 1., 1., 1., 1., 1., 2., 1., 1., 2., v0 ),  // kernel ring number
        mat4( 1., 1., 1., 1./12., 1., 5./6., 1., 1., 1., 1., 1., 1., 1., 1., 1., v0 ),  // kernel ring heights
        mat4( 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 11./12., 1., 0., 0., v0 ),
        mat4( v0, v0, v0, v0 ),
        mat4( 0.118, 0.174, 0.244, 0.114, 0.374, 0.222, 0.306, 0.449, 0.498, 0.295, 0.43, 0.353, 0.238, 0.39, 0.1, v0 ),  // growth center
        mat4( 0.0639, 0.159, 0.0287, 0.0469, 0.0822, 0.0294, 0.0775, 0.124, 0.1836, 0.1373, 0.0999, 0.0954, 0.0995, 0.1114, 0.0601, v1 ),  // growth width
        mat4( 0.082, 0.462, 0.496, 0.27, 0.518, 0.576, 0.324, 0.306, 0.544, 0.374, 0.33, 0.528, 0.498, 0.43, 0.26, v0 ),  // growth strength
        mat4( 0.85, 0.61, 0.5, 0.81, 0.85, 0.93, 0.88, 0.74, 0.97, 0.92, 0.56, 0.56, 0.95, 0.59, 0.58, v1 ),  // relative kernel radius
        0.175, 1.
    ),
    genome(
        10., 2.,
        mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 ),  // kernel ring number
        mat4( 1., 1., 1., 0., 1., 3./4., 1., 1., 1., 11./12., 3./4., 1., 1., 1./4., 1., v0 ),  // kernel ring heights
        mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 11./12., 0., 1., 0., v0 ),
        mat4( v0, v0, v0, v0 ),
        mat4( 0.175, 0.382, 0.231, 0.123, 0.398, 0.224, 0.193, 0.512, 0.427, 0.286, 0.508, 0.372, 0.196, 0.371, 0.246, v0 ),  // growth center
        mat4( 0.0682, 0.1568, 0.034, 0.0484, 0.0816, 0.0376, 0.063, 0.1189, 0.1827, 0.1422, 0.1079, 0.0724, 0.0934, 0.1107, 0.0672, v1 ),  // growth width
        mat4( 0.138, 0.544, 0.326, 0.256, 0.544, 0.544, 0.442, 0.198, 0.58, 0.282, 0.396, 0.618, 0.382, 0.374, 0.376, v0 ),  // growth strength
        mat4( 0.78, 0.56, 0.6, 0.84, 0.76, 0.82, 1.0, 0.68, 0.99, 0.72, 0.56, 0.65, 0.85, 0.54, 0.82, v1 ),  // relative kernel radius
        0.155, 1.
    )
);

// when matrix operation not available (e.g. exp, mod, equal, /), break down into four vec4 operations
mat4 intEqual4(in mat4 m, in ivec4 v) {
    return mat4( equal(ivec4(m[0]), v), equal(ivec4(m[1]), v), equal(ivec4(m[2]), v), equal(ivec4(m[3]), v) );
}
mat4 fract4(in mat4 m) {
    return mat4( fract(m[0]), fract(m[1]), fract(m[2]), fract(m[3]) );
}
mat4 exp4(in mat4 m) {
    return mat4( exp(m[0]), exp(m[1]), exp(m[2]), exp(m[3]) );
}

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

//const mat4[3] srcK = mat4[] ( intEqual4(src, iv0), intEqual4(src, iv1), intEqual4(src, iv2) ); 
//const mat4[3] dstK = mat4[] ( intEqual4(dst, iv0), intEqual4(dst, iv1), intEqual4(dst, iv2) ); 
const ivec4 src0 = ivec4(src[0]), src1 = ivec4(src[1]), src2 = ivec4(src[2]), src3 = ivec4(src[3]);
const ivec4 dst0 = ivec4(dst[0]), dst1 = ivec4(dst[1]), dst2 = ivec4(dst[2]), dst3 = ivec4(dst[3]);
const mat4[3] srcK = mat4[] (
    mat4( equal(src0, iv0), equal(src1, iv0), equal(src2, iv0), equal(src3, iv0) ),
    mat4( equal(src0, iv1), equal(src1, iv1), equal(src2, iv1), equal(src3, iv1) ),
    mat4( equal(src0, iv2), equal(src1, iv2), equal(src2, iv2), equal(src3, iv2) ) );
const mat4[3] dstK = mat4[] (
    mat4( equal(dst0, iv0), equal(dst1, iv0), equal(dst2, iv0), equal(dst3, iv0) ),
    mat4( equal(dst0, iv1), equal(dst1, iv1), equal(dst2, iv1), equal(dst3, iv1) ),
    mat4( equal(dst0, iv2), equal(dst1, iv2), equal(dst2, iv2), equal(dst3, iv2) ) );

const float valueSize = 64.;  // 6 bits value
const float speciesSize = 4.;  // 2 bits species: -1=none, 0=specie 1, 1=species 2...
const float valueMargin = 0.01;  // value 0..1 pack into 0.01..0.99
const float valueRange = 1. - 2. * valueMargin;
// higher digits = value, lower digits = species
ivec3 unpackSpecies(in vec3 texel) { return ivec3(fract(texel * valueSize) * speciesSize + 0.5) - 1; }
vec3  unpackValue  (in vec3 texel) { return texel; }  // return (floor(texel * valueSize) / valueSize - valueMargin) / valueRange; }
vec3 packTexel(in ivec3 species, in vec3 value) { return ( vec3(species+1) / speciesSize + floor((value*valueRange+valueMargin) * valueSize) ) / valueSize; }

// bell-shaped curve (Gaussian bump)
mat4 bell(in mat4 x, in mat4 m, in mat4 s) {
    return exp4( -mult(x-m, x-m) / s / s / 2. );
}

// Simplex 2D noise from https://www.shadertoy.com/view/Msf3WH
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
float voronoiNoise(in vec2 n, in float scale) {
    float i = 0.0, dis = 2.0;
    for(int x = -1;x<=1;x++)
    for(int y = -1;y<=1;y++) {
        vec2 p = floor(n/scale)+vec2(x,y);
        vec2 h2 = fract(cos(p*mat2(89.4,-75.7,-81.9,79.6))*343.42);
        float d = length(h2+vec2(x,y)-fract(n/scale));
        if (dis > d) { dis = d; i = fract(cos(p.x*89.42-p.y*75.7)*343.42); }
    }
    return i;
}

int randomSpecies(in vec2 p) {
    vec2 p1 = p / maxR / samplingDist;
    float noise = voronoiNoise(p1 + iDate.w * vec2(10., 10.), 8.);
    return int( noise * (float(speciesNum) + pNoSpecies) + 0.99 - pNoSpecies ) - 1;
}

vec3 randomValue(in vec2 p, in genome g) {
    vec2 p2 = p / g.R / samplingDist / g.randomScale;
    vec3 valueNoise = vec3( 
        simplexNoise(p2 + sin(iDate.w) * vec2(10.1, 10.2)), 
        simplexNoise(p2 + sin(iDate.w) * vec2(10.3, -10.4)), 
        simplexNoise(p2 + sin(iDate.w) * vec2(-10.5, 0.6)) );
    return clamp(g.baseNoise + biasNoise + valueNoise, 0., 1.);
}

vec3 randomTexel(in vec2 p) {
    int species = randomSpecies(p);
    vec3 value = species > -1 ? randomValue(p, genomes[species]) : vec3(0.);
    return packTexel(ivec3(species), value);
}

// get neighbor weights for given radius
mat4 getWeight(in float r, in genome g) {
    if (r > g.R) return m0;
    mat4 Br = g.betaLen * r / g.R / g.relR;  // scale radius by number of rings and relative radius
    #if use3RingKernels == 1
    mat4 height = mult(g.beta0, intEqual4(Br, iv0)) + mult(g.beta1, intEqual4(Br, iv1)) + mult(g.beta2, intEqual4(Br, iv2));
    #else
    mat4 height = mult(g.beta0, intEqual4(Br, iv0)) + mult(g.beta1, intEqual4(Br, iv1));
    #endif
    return mult(height, bell(fract4(Br), kmu, ksigma));
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
    return texel;
}

// simulation
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = fragCoord / iResolution.xy;

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


    mat4 avg_a = sum[0] / (total[0] + EPSILON);
    mat4 avg_b = sum[1] / (total[1] + EPSILON);

    // calculate growth (scaled by time step), reduce from kernels to destination channels
    mat4 growthK_a = mult(genomes[0].eta, bell(avg_a, genomes[0].mu, genomes[0].sigma) * 2. - 1.) / genomes[0].T;
    mat4 growthK_b = mult(genomes[1].eta, bell(avg_b, genomes[1].mu, genomes[1].sigma) * 2. - 1.) / genomes[1].T;
    vec3 growth_a = reduceKernels(growthK_a);
    vec3 growth_b = reduceKernels(growthK_b);

    // unpack current cell
    vec3 value = unpackValue(texel);

    // decide which species to be next, by maximum growth (or other criteria?)
    const vec3 aliveThreshold = vec3(-0.7);
    ivec3 select_a = ivec3(greaterThan(growth_a, aliveThreshold)) * ivec3(greaterThanEqual(growth_a, growth_b));
    ivec3 select_b = ivec3(greaterThan(growth_b, aliveThreshold)) * ivec3(greaterThan     (growth_b, growth_a));
    ivec3 species = 0 * select_a + 1 * select_b;

    // choose growth according to species, add to original value
    vec3 is_none = vec3(equal(species, ivec3(-1)));
    vec3 growth = growth_a * float(select_a) + growth_b * float(select_b) + vec3(-0.1) * is_none;
    value = clamp(growth + value, 0., 1.);

/*
    // unpack current cell
    //ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    #if useDecideByMaxGrowth == 1
    vec3[speciesNum] growthList;
    vec3 maxGrowth = vec3(-99.);
    ivec3 argmax = ivec3(-99);
    for (int s=0; s<speciesNum; s++) {
        genome g = genomes[s];

        // weighted average = weighted sum / total weights, avoid divided by zero
        mat4 avg = sum[s] / (total[s] + EPSILON);    // avoid divided by zero

        // calculate growth (scaled by time step), reduce from kernels to destination channels
        mat4 growthK = mult(g.eta, bell(avg, g.mu, g.sigma) * 2. - 1.) / g.T;
        vec3 growthS = reduceKernels(growthK);
        growthList[s] = growthS;

        // calculate max and argmax of growth
        ivec3 isMore = ivec3(greaterThan(growthS, maxGrowth));
        maxGrowth = mix(maxGrowth, growthS, float(isMore));
        //maxGrowth = maxGrowth * float(ic1-isMore) + growthS * float(isMore);
        argmax = argmax * (ic1-isMore) + s * isMore;
    }

    vec3 growth = vec3(-0.1);  // (-0.1)
    ivec3 species = ivec3(-1);
    for (int s=0; s<speciesNum; s++) {
        // decide which species to be next, by maximum growth (or other criteria?)
        vec3 growthS = growthList[s];
        ivec3 isSelect = ivec3(greaterThan(growthS, aliveThreshold)) * ivec3(equal(ivec3(s), argmax));

        growth = mix(growth, growthS, vec3(isSelect));
        //finalGrowth = finalGrowth * vec3(ic1-isSelect) + growthS * vec3(isSelect);
        species = species * (ic1-isSelect) + s * isSelect;
    }

    value = clamp(growth + value, 0., 1.);

    #elif useDecideByMaxGrowth == 0

    mat4[speciesNum] avgList;
    vec3 maxAvg = vec3(-99.);
    ivec3 argmax = ivec3(-99);
    for (int s=0; s<speciesNum; s++) {
        // weighted average = weighted sum / total weights, avoid divided by zero
        mat4 avgK = sum[s] / (total[s] + EPSILON);    // avoid divided by zero
        avgList[s] = avgK;
        vec3 avgS = reduceKernels(avgK);

        // calculate max and argmax of avg
        ivec3 isMore = ivec3(greaterThan(avgS, maxAvg));
        maxAvg = mix(maxAvg, avgS, float(isMore));
        argmax = argmax * (ic1-isMore) + s * isMore;
    }

    vec3 growth = vec3(-0.1);  // (-0.1)
    ivec3 species = ivec3(-1);
    for (int s=0; s<speciesNum; s++) {
        genome g = genomes[s];

        // decide which species to be next, by maximum neighbor avg (or other criteria?)
        mat4 avgK = avgList[s];
        mat4 growthK = mult(g.eta, bell(avgK, g.mu, g.sigma) * 2. - 1.) / g.T;
        vec3 growthS = reduceKernels(growthK);

        ivec3 isSelect = ivec3(greaterThan(growthS, aliveThreshold)) * ivec3(equal(ivec3(s), argmax));
        growth = mix(growth, growthS, vec3(isSelect));
        species = species * (ic1-isSelect) + s * isSelect;
    }

    value = clamp(growth + value, 0., 1.);

    #endif
*/
    // randomize at start
    if (iFrame == 0)
        texel = randomTexel(fragCoord);
    else
        texel = packTexel(species, value);

    fragColor = vec4(texel, 1.);
}
