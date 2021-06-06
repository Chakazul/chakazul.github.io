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


// modified from SmoothLife by davidar - https://www.shadertoy.com/view/Msy3RD

const float R = 15.;       // space resolution = kernel radius
const float T = 10.;       // time resolution = number of divisions per unit time
const float dt = 1./T;     // time step
const float mu = 0.14;     // growth center
const float sigma = 0.014; // growth width
const float rho = 0.5;     // kernel center
const float omega = 0.15;  // kernel width

float bell(float x, float m, float s)
{
    return exp(-(x-m)*(x-m)/s/s/2.);  // bell-shaped curve
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord / iResolution.xy;

    float sum = 0.;
    float total = 0.;
    for (int x=-int(R); x<=int(R); x++)
    for (int y=-int(R); y<=int(R); y++)
    {
        float r = sqrt(float(x*x + y*y)) / R;
        vec2 txy = mod((fragCoord + vec2(x,y)) / iResolution.xy, 1.);
        float val = texture(iChannel0, txy).r;
        float weight = bell(r, rho, omega);
        sum += val * weight;
        total += weight;
    }
    float avg = sum / total;

    float val = texture(iChannel0, uv).r;
    float growth = bell(avg, mu, sigma) * 2. - 1.;
    float c = clamp(val + dt * growth, 0., 1.);

    if (iFrame < 1) // || iMouse.z > 0.
        c = 0.013 + noise(fragCoord/R + mod(iDate.w,1.)*100.);
    if (iMouse.z > 0.) {
        float d = length((fragCoord.xy - iMouse.xy) / iResolution.xx);
        if (d <= R/iResolution.x)
        	c = 0.02 + noise(fragCoord/R + mod(iDate.w,1.)*100.);
    }

    fragColor = vec4(c,c,c,1.);
}
