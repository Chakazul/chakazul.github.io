// from https://smoothstep.io/anim/41ff4ba126d9

// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

// multiple species (A, B...), can have different spatial scales (R) and time scales (T) 
// maximum 16 kernels by using 4x4 matrix
// when matrix operation not available (e.g. exp, mod, equal, /), split into four vec4 operations

// common
#define EPSILON 0.000001
#define mult matrixCompMult
#define speciesNum 3
#define species_0 0
#define species_a 1
#define species_b 2

const ivec4 iv0 = ivec4(0);
const ivec4 iv1 = ivec4(1);
const ivec4 iv2 = ivec4(2);
const ivec4 iv3 = ivec4(3);
const vec4 v0 = vec4(0.);
const vec4 v1 = vec4(1.);
const mat4 m0 = mat4(v0, v0, v0, v0);
const mat4 m1 = mat4(v1, v1, v1, v1);

/// choose pixel size
const float pixelSize = 1.;

// choose species A, change numbers at end of xxx_?
const float R_a = 15.;  // space resolution = kernel radius
const float T_a = 2.;  // time resolution = number of divisions per unit time
#define baseNoise_a baseNoise_2
#define   betaLen_a   betaLen_2
#define     beta0_a     beta0_2
#define     beta1_a     beta1_2
#define     beta2_a     beta2_2
#define        mu_a        mu_2
#define     sigma_a     sigma_2
#define       eta_a       eta_2
#define      relR_a      relR_2
#define       src_a       src_2
#define       dst_a       dst_2

// choose species B, change numbers at end of xxx_?
const float R_b = 10.;  // space resolution = kernel radius
const float T_b = 2.;  // time resolution = number of divisions per unit time
#define baseNoise_b baseNoise_1
#define   betaLen_b   betaLen_1
#define     beta0_b     beta0_1
#define     beta1_b     beta1_1
#define     beta2_b     beta2_1
#define        mu_b        mu_1
#define     sigma_b     sigma_1
#define       eta_b       eta_1
#define      relR_b      relR_1
#define       src_b       src_1
#define       dst_b       dst_1


// species list

// species: VT049W Tessellatium (sometimes reproductive)
const float baseNoise_0 = 0.185;
const mat4    betaLen_0 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_0 = mat4( 1., 1., 1., 0., 1., 5./6., 1., 1., 1., 11./12., 3./4., 11./12., 1., 1./6., 1., v0 );  // kernel ring heights
const mat4      beta1_0 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 1., 0., 1., 0., v0 );
const mat4      beta2_0 = mat4( v0, v0, v0, v0 );
const mat4         mu_0 = mat4( 0.272, 0.349, 0.2, 0.114, 0.447, 0.247, 0.21, 0.462, 0.446, 0.327, 0.476, 0.379, 0.262, 0.412, 0.201, v0 );  // growth center
const mat4      sigma_0 = mat4( 0.0595, 0.1585, 0.0332, 0.0528, 0.0777, 0.0342, 0.0617, 0.1192, 0.1793, 0.1408, 0.0995, 0.0697, 0.0877, 0.1101, 0.0786, v1 );  // growth width
const mat4        eta_0 = mat4( 0.19, 0.66, 0.39, 0.38, 0.74, 0.92, 0.59, 0.37, 0.94, 0.51, 0.77, 0.92, 0.71, 0.59, 0.41, v0 );  // growth strength
const mat4       relR_0 = mat4( 0.91, 0.62, 0.5, 0.97, 0.72, 0.8, 0.96, 0.56, 0.78, 0.79, 0.5, 0.72, 0.68, 0.55, 0.82, v1 );  // relative kernel radius
const mat4        src_0 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_0 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: Z18A9R Tessellatium (highly reproductive) (modified for lower reproduction)
const float baseNoise_1 = 0.155;
const mat4    betaLen_1 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_1 = mat4( 1., 1., 1., 0., 1., 3./4., 1., 1., 1., 11./12., 3./4., 1., 1., 1./4., 1., v0 );  // kernel ring heights
const mat4      beta1_1 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 11./12., 0., 1., 0., v0 );
const mat4      beta2_1 = mat4( v0, v0, v0, v0 );
const mat4         mu_1 = mat4( 0.175, 0.382, 0.231, 0.123, 0.398, 0.224, 0.193, 0.512, 0.427, 0.286, 0.508, 0.372, 0.196, 0.371, 0.246, v0 );  // growth center
//const mat4      sigma_1 = mat4( 0.0682, 0.1568, 0.034, 0.0484, 0.0816, 0.0376, 0.063, 0.1189, 0.1827, 0.1422, 0.1079, 0.0724, 0.0934, 0.1107, 0.0712, v1 );  // growth width
const mat4      sigma_1 = mat4( 0.0682, 0.1568, 0.034, 0.0484, 0.0816, 0.0376, 0.063, 0.1189, 0.1827, 0.1422, 0.1079, 0.0724, 0.0934, 0.1107, 0.0672, v1 );  // growth width
const mat4        eta_1 = mat4( 0.138, 0.544, 0.326, 0.256, 0.544, 0.544, 0.442, 0.198, 0.58, 0.282, 0.396, 0.618, 0.382, 0.374, 0.376, v0 );  // growth strength
const mat4       relR_1 = mat4( 0.78, 0.56, 0.6, 0.84, 0.76, 0.82, 1.0, 0.68, 0.99, 0.72, 0.56, 0.65, 0.85, 0.54, 0.82, v1 );  // relative kernel radius
const mat4        src_1 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_1 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: G6G6CR Ciliatium (immune system) (modified for higher cilia production)
const float baseNoise_2 = 0.175;
const mat4    betaLen_2 = mat4( 1., 1., 1., 2., 1., 2., 1., 1., 1., 1., 1., 2., 1., 1., 2., v0 );  // kernel ring number
const mat4      beta0_2 = mat4( 1., 1., 1., 1./12., 1., 5./6., 1., 1., 1., 1., 1., 1., 1., 1., 1., v0 );  // kernel ring heights
const mat4      beta1_2 = mat4( 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 11./12., 1., 0., 0., v0 );
const mat4      beta2_2 = mat4( v0, v0, v0, v0 );
const mat4         mu_2 = mat4( 0.118, 0.174, 0.244, 0.114, 0.374, 0.222, 0.306, 0.449, 0.498, 0.295, 0.43, 0.353, 0.238, 0.39, 0.1, v0 );  // growth center
const mat4      sigma_2 = mat4( 0.0639, 0.159, 0.0287, 0.0469, 0.0822, 0.0294, 0.0775, 0.124, 0.1836, 0.1373, 0.0999, 0.0954, 0.0995, 0.1114, 0.0601, v1 );  // growth width
//const mat4      sigma_2 = mat4( 0.0639, 0.159, 0.0287, 0.0469, 0.0822, 0.0294, 0.0775, 0.124, 0.1836, 0.1373, 0.0999, 0.0754, 0.0995, 0.1144, 0.0601, v1 );  // growth width
const mat4        eta_2 = mat4( 0.082, 0.462, 0.496, 0.27, 0.518, 0.576, 0.324, 0.306, 0.544, 0.374, 0.33, 0.528, 0.498, 0.43, 0.26, v0 );  // growth strength
const mat4       relR_2 = mat4( 0.85, 0.61, 0.5, 0.81, 0.85, 0.93, 0.88, 0.74, 0.97, 0.92, 0.56, 0.56, 0.95, 0.59, 0.58, v1 );  // relative kernel radius
const mat4        src_2 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_2 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: tri-color ghosts
// set R = 10., T = 10.
const float baseNoise_3 = 0.185;
const mat4    betaLen_3 = mat4( 2., 3., 1., 2., 3., 1., 2., 3., 1., v0, v0 );  // kernel ring number
const mat4      beta0_3 = mat4( 1./4., 1., 1., 1./4., 1., 1., 1./4., 1., 1., v0, v0 );  // kernel ring heights
const mat4      beta1_3 = mat4( 1., 3./4., 0., 1., 3./4., 0., 1., 3./4., 0., v0, v0 );
const mat4      beta2_3 = mat4( 0., 3./4., 0., 0., 3./4., 0., 0., 3./4., 0., v0, v0 );
const mat4         mu_3 = mat4( 0.16, 0.22, 0.28, 0.16, 0.22, 0.28, 0.16, 0.22, 0.28, v0, v0 );  // growth center
const mat4      sigma_3 = mat4( 0.025, 0.042, 0.025, 0.025, 0.042, 0.025, 0.025, 0.042, 0.025, v1, v1 );  // growth width
const mat4        eta_3 = mat4( 0.666, 0.666, 0.666, 0.666, 0.666, 0.666, 0.666, 0.666, 0.666, v0, v0 );     // growth strength
const mat4       relR_3 = mat4( 1., 1., 1., 1., 1., 1., 1., 1., 1., v0, v0 );  // relative kernel radius
const mat4        src_3 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., v0, v0 );  // source channels
const mat4        dst_3 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., v0, v0 );  // destination channels

// more species

// species: KH97WU Tessellatium (courting, slightly reproductive)
const float baseNoise_4 = 0.185;
const mat4    betaLen_4 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 1., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_4 = mat4( 1., 1., 1., 0., 1., 5./6., 1., 1., 1., 11./12., 3./4., 1., 1., 1./6., 1., v0 );  // kernel ring heights
const mat4      beta1_4 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 0., 0., 1., 0., v0 );
const mat4      beta2_4 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_4 = mat4( 0.204, 0.359, 0.176, 0.128, 0.386, 0.229, 0.181, 0.466, 0.466, 0.37, 0.447, 0.391, 0.299, 0.398, 0.183, v0 );  // growth center
const mat4      sigma_4 = mat4( 0.0574, 0.152, 0.0314, 0.0545, 0.0825, 0.0348, 0.0657, 0.1224, 0.1789, 0.1372, 0.1064, 0.0644, 0.0891, 0.1065, 0.0773, v1 );  // growth width
const mat4        eta_4 = mat4( 0.116, 0.448, 0.332, 0.392, 0.398, 0.614, 0.448, 0.224, 0.624, 0.352, 0.342, 0.634, 0.362, 0.472, 0.242, v0 );  // growth strength
const mat4       relR_4 = mat4( 0.93, 0.59, 0.58, 0.97, 0.79, 0.87, 1.0, 0.64, 0.67, 0.68, 0.5, 0.85, 0.69, 0.87, 0.66, v1 );  // relative kernel radius
const mat4        src_4 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_4 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: XEH4YR Tessellatium (explosive)
const float baseNoise_5 = 0.175;
const mat4    betaLen_5 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 3., 1., v0 );  // kernel ring number
const mat4      beta0_5 = mat4( 1., 1., 1., 0., 1., 5./6., 1., 1., 1., 11./12., 3./4., 11./12., 1., 1./6., 1., v0 );  // kernel ring heights
const mat4      beta1_5 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 1., 0., 1., 0., v0 );
const mat4      beta2_5 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_5 = mat4( 0.282, 0.354, 0.197, 0.164, 0.406, 0.251, 0.259, 0.517, 0.455, 0.264, 0.472, 0.417, 0.208, 0.395, 0.184, v0 );  // growth center
const mat4      sigma_5 = mat4( 0.0646, 0.1584, 0.0359, 0.056, 0.0738, 0.0383, 0.0665, 0.1164, 0.1806, 0.1437, 0.0939, 0.0666, 0.0815, 0.1049, 0.0748, v1 );  // growth width
const mat4        eta_5 = mat4( 0.082, 0.544, 0.26, 0.294, 0.508, 0.56, 0.326, 0.21, 0.638, 0.346, 0.384, 0.748, 0.44, 0.366, 0.294, v0 );  // growth strength
const mat4       relR_5 = mat4( 0.85, 0.62, 0.69, 0.84, 0.82, 0.86, 1.0, 0.5, 0.78, 0.6, 0.5, 0.7, 0.67, 0.6, 0.8, v1 );  // relative kernel radius
const mat4        src_5 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_5 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: HAESRE Tessellatium (zigzaging)
const float baseNoise_6 = 0.185;
const mat4    betaLen_6 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_6 = mat4( 1., 1., 1., 0., 1., 3./4., 1., 1., 1., 11./12., 5./6., 1., 1., 1./4., 1., v0 );  // kernel ring heights
const mat4      beta1_6 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 11./12., 0., 1., 0., v0 );
const mat4      beta2_6 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_6 = mat4( 0.272, 0.337, 0.129, 0.132, 0.429, 0.239, 0.25, 0.497, 0.486, 0.276, 0.425, 0.352, 0.21, 0.381, 0.244, v0 );  // growth center
const mat4      sigma_6 = mat4( 0.0674, 0.1576, 0.0382, 0.0514, 0.0813, 0.0409, 0.0691, 0.1166, 0.1751, 0.1344, 0.1026, 0.0797, 0.0921, 0.1056, 0.0813, v1 );  // growth width
const mat4        eta_6 = mat4( 0.15, 0.474, 0.342, 0.192, 0.524, 0.598, 0.426, 0.348, 0.62, 0.338, 0.314, 0.608, 0.292, 0.426, 0.346, v0 );  // growth strength
const mat4       relR_6 = mat4( 0.87, 0.65, 0.67, 0.98, 0.77, 0.83, 1.0, 0.7, 0.99, 0.69, 0.7, 0.57, 0.89, 0.84, 0.76, v1 );  // relative kernel radius
const mat4        src_6 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_6 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: GDNQYX Tessellatium (stable)
const float baseNoise_7 = 0.205;
const mat4    betaLen_7 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_7 = mat4( 1., 1., 1., 0., 1., 5./6., 1., 1., 1., 11./12., 3./4., 1., 1., 1./6., 1., v0 );  // kernel ring heights
const mat4      beta1_7 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 11./12., 0., 1., 0., v0 );
const mat4      beta2_7 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_7 = mat4( 0.242, 0.375, 0.194, 0.122, 0.413, 0.221, 0.192, 0.492, 0.426, 0.361, 0.464, 0.361, 0.235, 0.381, 0.216, v0 );  // growth center
const mat4      sigma_7 = mat4( 0.061, 0.1553, 0.0361, 0.0531, 0.0774, 0.0365, 0.0649, 0.1219, 0.1759, 0.1381, 0.1044, 0.0686, 0.0924, 0.1118, 0.0748, v1 );  // growth width
const mat4        eta_7 = mat4( 0.144, 0.506, 0.332, 0.3, 0.502, 0.58, 0.344, 0.268, 0.582, 0.326, 0.418, 0.642, 0.39, 0.378, 0.294, v0 );  // growth strength
const mat4       relR_7 = mat4( 0.98, 0.59, 0.5, 0.93, 0.73, 0.88, 0.93, 0.61, 0.84, 0.7, 0.57, 0.73, 0.74, 0.87, 0.72, v1 );  // relative kernel radius
const mat4        src_7 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_7 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: 5N7KKM Tessellatium (stable)
const float baseNoise_8 = 0.175;
const mat4    betaLen_8 = mat4( 1., 1., 2., 2., 1., 2., 1., 1., 1., 2., 2., 2., 1., 2., 1., v0 );  // kernel ring number
const mat4      beta0_8 = mat4( 1., 1., 1., 0., 1., 3./4., 1., 1., 1., 11./12., 3./4., 1., 1., 1./6., 1., v0 );  // kernel ring heights
const mat4      beta1_8 = mat4( 0., 0., 1./4., 1., 0., 1., 0., 0., 0., 1., 1., 11./12., 0., 1., 0., v0 );
const mat4      beta2_8 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_8 = mat4( 0.22, 0.351, 0.177, 0.126, 0.437, 0.234, 0.179, 0.489, 0.419, 0.341, 0.469, 0.369, 0.219, 0.385, 0.208, v0 );  // growth center
const mat4      sigma_8 = mat4( 0.0628, 0.1539, 0.0333, 0.0525, 0.0797, 0.0369, 0.0653, 0.1213, 0.1775, 0.1388, 0.1054, 0.0721, 0.0898, 0.1102, 0.0749, v1 );  // growth width
const mat4        eta_8 = mat4( 0.174, 0.46, 0.31, 0.242, 0.508, 0.566, 0.406, 0.27, 0.588, 0.294, 0.388, 0.62, 0.348, 0.436, 0.39, v0 );  // growth strength
const mat4       relR_8 = mat4( 0.87, 0.52, 0.58, 0.89, 0.78, 0.79, 1.0, 0.64, 0.96, 0.66, 0.69, 0.61, 0.81, 0.81, 0.71, v1 );  // relative kernel radius
const mat4        src_8 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_8 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels

// species: Y3CS55 fast emitter
const float baseNoise_9 = 0.165;
const mat4    betaLen_9 = mat4( 1., 1., 1., 2., 1., 2., 1., 1., 1., 1., 1., 3., 1., 1., 2., v0 );  // kernel ring number
const mat4      beta0_9 = mat4( 1., 1., 1., 1./12., 1., 5./6., 1., 1., 1., 1., 1., 1., 1., 1., 1., v0 );  // kernel ring heights
const mat4      beta1_9 = mat4( 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 11./12., 0., 0., 1./12., v0 );
const mat4      beta2_9 = mat4( 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., v0 );
const mat4         mu_9 = mat4( 0.168, 0.1, 0.265, 0.111, 0.327, 0.223, 0.293, 0.465, 0.606, 0.404, 0.377, 0.297, 0.319, 0.483, 0.1, v0 );  // growth center
const mat4      sigma_9 = mat4( 0.062, 0.1495, 0.0488, 0.0555, 0.0763, 0.0333, 0.0724, 0.1345, 0.1807, 0.1413, 0.1136, 0.0701, 0.1038, 0.1185, 0.0571, v1 );  // growth width
const mat4        eta_9 = mat4( 0.076, 0.562, 0.548, 0.306, 0.568, 0.598, 0.396, 0.298, 0.59, 0.396, 0.156, 0.426, 0.558, 0.388, 0.132, v0 );  // growth strength
const mat4       relR_9 = mat4( 0.58, 0.68, 0.5, 0.87, 1.0, 1.0, 0.88, 0.88, 0.86, 0.98, 0.63, 0.53, 1.0, 0.89, 0.59, v1 );  // relative kernel radius
const mat4        src_9 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 0., 0., 1., 1., 2., 2., v0 );  // source channels
const mat4        dst_9 = mat4( 0., 0., 0., 1., 1., 1., 2., 2., 2., 1., 2., 0., 2., 0., 1., v0 );  // destination channels


// common

#define vecMax(v) max(max(max(v.x, v.y), v.z), v.w)
#define matMax(m) max(max(max(vecMax(m[0]), vecMax(m[1])), vecMax(m[2])), vecMax(m[3]))
// const float maxRelR_a = matMax(relR_a);
// const float maxRelR_b = matMax(relR_b);

const vec4 kmv = vec4(0.5);    // kernel ring center
const mat4 kmu = mat4(kmv, kmv, kmv, kmv);
const vec4 ksv = vec4(0.15);    // kernel ring width
const mat4 ksigma = mat4(ksv, ksv, ksv, ksv);

const ivec4 src0_a = ivec4(src_a[0]), src1_a = ivec4(src_a[1]), src2_a = ivec4(src_a[2]), src3_a = ivec4(src_a[3]);
const ivec4 src0_b = ivec4(src_b[0]), src1_b = ivec4(src_b[1]), src2_b = ivec4(src_b[2]), src3_b = ivec4(src_b[3]);


// bell-shaped curve (Gaussian bump)
mat4 bell(in mat4 x, in mat4 m, in mat4 s)
{
    mat4 v = -mult(x-m, x-m) / s / s / 2.;
    return mat4( exp(v[0]), exp(v[1]), exp(v[2]), exp(v[3]) );
}

// pack / unpack texels

/**/
// high precision = value, low precision = species
const float highSize = 64.;  // 6 bits value
const float lowSize = 4.;  // 2 bits species = none + max 3 species
const float valMargin = 0.01;  // value 0..1 pack into 0.01..0.99
const float valBand = 1. - 2. * valMargin;
ivec3 unpackSpecies(in vec3 texel)
{
    return ivec3(fract(texel * highSize) * lowSize + 0.5);
}
vec3 unpackValue(in vec3 texel)
{
    return (floor(texel * highSize) / highSize -valMargin)/valBand;
}
vec3 packTexel(in ivec3 species, in vec3 val)
{
    return (floor((val*valBand+valMargin) * highSize) + vec3(species) / lowSize) / highSize;
}
/*/
// high precision = species, low precision = value
ivec3 unpackSpecies(in vec3 texel)
{
    return ivec3(floor(texel * float(speciesNum)));
}
vec3 unpackValue(in vec3 texel)
{
    return (fract(texel * float(speciesNum)) - 0.1) / 0.8;
}
vec3 packTexel(in ivec3 species, in vec3 val)
{
    return (vec3(species) + val * 0.8 + 0.1) / float(speciesNum);
}
/**/

// ----------------------------------------

// Noise simplex 2D by iq - https://www.shadertoy.com/view/Msf3WH

vec2 hash( vec2 p )
{
	p = vec2( dot(p,vec2(127.1,311.7)), dot(p,vec2(269.5,183.3)) );
	return -1.0 + 2.0*fract(sin(p)*43758.5453123);
}

float noise( in vec2 p )
{
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


// get neighbor weights for given radius
mat4 getWeight(in float r, in mat4 relR, in mat4 betaLen, in mat4 beta0, in mat4 beta1, in mat4 beta2)
{
    if (r > 1.) return m0;
    mat4 Br = betaLen / relR * r;  // scale radius by number of rings and relative radius
    ivec4 Br0 = ivec4(Br[0]), Br1 = ivec4(Br[1]), Br2 = ivec4(Br[2]), Br3 = ivec4(Br[3]);

    // get heights of kernel rings
    // for each Br: Br==0 ? beta0 : Br==1 ? beta1 : Br==2 ? beta2 : 0
    mat4 height = mat4(
        beta0[0] * vec4(equal(Br0, iv0)) + beta1[0] * vec4(equal(Br0, iv1)),// + beta2[0] * vec4(equal(Br0, iv2)),
        beta0[1] * vec4(equal(Br1, iv0)) + beta1[1] * vec4(equal(Br1, iv1)),// + beta2[1] * vec4(equal(Br1, iv2)),
        beta0[2] * vec4(equal(Br2, iv0)) + beta1[2] * vec4(equal(Br2, iv1)),// + beta2[2] * vec4(equal(Br2, iv2)),
        beta0[3] * vec4(equal(Br3, iv0)) + beta1[3] * vec4(equal(Br3, iv1)) );// + beta2[3] * vec4(equal(Br3, iv2)) );
    mat4 mod1 = mat4( mod(Br[0], 1.), mod(Br[1], 1.), mod(Br[2], 1.), mod(Br[3], 1.) );

    return mult(height, bell(mod1, kmu, ksigma));
}

// draw the shape of kernels
vec3 drawKernel(in vec2 uv, in mat4 relR, in mat4 betaLen, in mat4 beta0, in mat4 beta1, in mat4 beta2)
{
    ivec2 ij = ivec2(uv / 0.25);  // divide screen into 4x4 grid: 0..3 x 0..3
    vec2 xy = mod(uv, 0.25) * 8. - 1.;  // coordinates in each grid cell: -1..1 x -1..1
    if (ij.x > 3 || ij.y > 3) return vec3(0.);
    float r = length(xy);
    mat4 weight = getWeight(r, relR, betaLen, beta0, beta1, beta2);
    vec3 rgb = vec3(weight[3-ij.y][ij.x]);
    return rgb;
}

// map values from source channels to kernels
vec4 mapK(in vec3 v, in ivec4 srcv)
{
    // for each src: src==0 ? r : src==1 ? g : src==2 ? b : 0
    return
        v.r * vec4(equal(srcv, iv0)) + 
        v.g * vec4(equal(srcv, iv1)) +
        v.b * vec4(equal(srcv, iv2));
}

// reduce values from kernels to destination channels
float reduceK(in mat4 m, in ivec4 ch, in mat4 dst)
{
    // (r,g,b) = sum of m where dst==(0,1,2)
    return 
        dot(m[0], vec4(equal(ivec4(dst[0]), ch))) + 
        dot(m[1], vec4(equal(ivec4(dst[1]), ch))) + 
        dot(m[2], vec4(equal(ivec4(dst[2]), ch))) + 
        dot(m[3], vec4(equal(ivec4(dst[3]), ch)));
}

// add to weighted sum and total for given cell, return cell
vec3 addSum(in vec2 xy, 
    in mat4 weight_a, inout mat4 sum_a, inout mat4 total_a, 
    in mat4 weight_b, inout mat4 sum_b, inout mat4 total_b)
{
    // get neighbor cell, unpack species and value
    vec2 uv = xy / iResolution.xy;
    vec3 texel = texture(iPrevFrame, uv).rgb;
    ivec3 species = unpackSpecies(texel);
    vec3 value = unpackValue(texel);

    // get value according to species
    vec3 is_a = vec3(equal(species, ivec3(species_a)));
    vec3 is_b = vec3(equal(species, ivec3(species_b)));
    vec3 value_a = value * is_a;
    vec3 value_b = value * is_b;

    // map values to kernels
    mat4 valueK_a = mat4( mapK(value_a, src0_a), mapK(value_a, src1_a), mapK(value_a, src2_a), mapK(value_a, src3_a) );
    mat4 valueK_b = mat4( mapK(value_b, src0_b), mapK(value_b, src1_b), mapK(value_b, src2_b), mapK(value_b, src3_b) );

    // add to weighted sum and total according to species
    sum_a += mult(valueK_a, weight_a);
    sum_b += mult(valueK_b, weight_b);
    total_a += weight_a;
    total_b += weight_b;
    return texel;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    if (any(greaterThan(mod(fragCoord, vec2(pixelSize)), vec2(0.6))))
        discard;

    vec2 uv = fragCoord / iResolution.xy;

    // loop through the neighborhood, optimized: same weights for all quadrants/octants
    // calculate the weighted average of neighborhood from source channel
    mat4 sum_a = mat4(0.), total_a = mat4(0.);
    mat4 sum_b = mat4(0.), total_b = mat4(0.);


    // current cell
    float r = 0.;
    mat4 weight_a = getWeight(r / R_a, relR_a, betaLen_a, beta0_a, beta1_a, beta2_a);
    mat4 weight_b = getWeight(r / R_b, relR_b, betaLen_b, beta0_b, beta1_b, beta2_b);
    vec3 texel = addSum(fragCoord + vec2(0, 0)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);

    // orthogonal neighbors, repeat for 4 quadrants
    const int maxR = int(ceil(max(R_a, R_b)));
    for (int x=1; x<=maxR; x++)
    {
        r = float(x);
        weight_a = getWeight(r / R_a, relR_a, betaLen_a, beta0_a, beta1_a, beta2_a);
        weight_b = getWeight(r / R_b, relR_b, betaLen_b, beta0_b, beta1_b, beta2_b);
        addSum(fragCoord + vec2(+x, 0)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(-x, 0)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(0, +x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(0, -x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
    }

    // diagonal neighbors, repeat for 4 quadrants
    const int diagR = int(ceil(float(maxR) / sqrt(2.)));
    for (int x=1; x<=diagR; x++)
    {
        r = sqrt(2.) * float(x);
        weight_a = getWeight(r / R_a, relR_a, betaLen_a, beta0_a, beta1_a, beta2_a);
        weight_b = getWeight(r / R_b, relR_b, betaLen_b, beta0_b, beta1_b, beta2_b);
        addSum(fragCoord + vec2(+x, +x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(+x, -x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(-x, +x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        addSum(fragCoord + vec2(-x, -x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
    }

    // others neighbors, repeat for 8 octants
    for (int y=1; y<=maxR-1; y++)
    for (int x=y+1; x<=maxR; x++)
    {
        r = sqrt(float(x*x + y*y));
        if (r <= float(maxR)) 
        {
            weight_a = getWeight(r / R_a, relR_a, betaLen_a, beta0_a, beta1_a, beta2_a);
            weight_b = getWeight(r / R_b, relR_b, betaLen_b, beta0_b, beta1_b, beta2_b);
            addSum(fragCoord + vec2(+x, +y)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(+x, -y)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(-x, +y)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(-x, -y)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(+y, +x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(+y, -x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(-y, +x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
            addSum(fragCoord + vec2(-y, -x)*pixelSize, weight_a, sum_a, total_a, weight_b, sum_b, total_b);
        }
    }
    // weighted average = weighted sum / total weights, avoid divided by zero
    mat4 avg_a = sum_a / (total_a + EPSILON);
    mat4 avg_b = sum_b / (total_b + EPSILON);

    // calculate growth (scaled by time step), reduce from kernels to destination channels
    mat4 growthK_a = mult(eta_a, bell(avg_a, mu_a, sigma_a) * 2. - 1.) / T_a;
    mat4 growthK_b = mult(eta_b, bell(avg_b, mu_b, sigma_b) * 2. - 1.) / T_b;
    vec3 growth_a = vec3( reduceK(growthK_a, iv0, dst_a), reduceK(growthK_a, iv1, dst_a), reduceK(growthK_a, iv2, dst_a) );
    vec3 growth_b = vec3( reduceK(growthK_b, iv0, dst_b), reduceK(growthK_b, iv1, dst_b), reduceK(growthK_b, iv2, dst_b) );

    // unpack current cell
    vec3 value = unpackValue(texel);

    // decide which species to be next, by maximum growth (or other criteria?)
    const vec3 aliveThreshold = vec3(-0.7);
    ivec3 select_a = ivec3(greaterThan(growth_a, aliveThreshold)) * ivec3(greaterThanEqual(growth_a, growth_b));
    ivec3 select_b = ivec3(greaterThan(growth_b, aliveThreshold)) * ivec3(greaterThan     (growth_b, growth_a));
    ivec3 species = species_a * select_a + species_b * select_b;
    //ivec3 species = unpackSpecies(texel);
    //ivec3 select_0 = max(ivec3(1) - select_a - select_b, ivec3(0));
    //species = species_a * select_a + species_b * select_b + species * select_0;

    // choose growth according to species, add to original value
    vec3 is_none = vec3(equal(species, ivec3(species_0)));
    vec3 growth = growth_a * float(select_a) + growth_b * float(select_b) + vec3(-0.1) * is_none;
    value = clamp(growth + value, 0., 1.);

    // debug: uncomment to show list of kernels
    //val = drawKernel(fragCoord / iResolution.y, relR_a, betaLen_a, beta0_a, beta1_a, beta2_a);
    //val = drawKernel(fragCoord / iResolution.y, relR_b, betaLen_b, beta0_a, beta1_b, beta2_b);
    //species = ivec3(0);

    // initialize with random species and values
    if (iFrame == 0)
    {
        float speciesNoise_a = noise(fragCoord/R_a/pixelSize/8. + 17.);
        float speciesNoise_b = noise(fragCoord/R_b/pixelSize/8. + 23.);
        bool is_a = (speciesNoise_a > 0.);
        bool is_b = (speciesNoise_b > -0.1) && !is_a;
        species = ivec3(species_a * int(is_a) + species_b * int(is_b));
        float R = R_a * float(is_a) + R_b * float(is_b);
        float baseNoise = baseNoise_a * float(is_a) + baseNoise_b * float(is_b);
        vec3 valueNoise = vec3( 
            noise(fragCoord/R/pixelSize), 
            noise(fragCoord/R/pixelSize + 13.), 
            noise(fragCoord/R/pixelSize + 27.) );
        value = clamp(baseNoise + valueNoise, 0., 1.);
    }
    else if (iFrame==1) {}

    texel = packTexel(species, value);
    fragColor = vec4(texel, 1.);
}
