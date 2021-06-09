// from https://www.shadertoy.com/view/XlB3zW

#define speciesNum 2

float voronoi(vec2 n, float s) {
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
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 coord = fragCoord.xy+iTime*32.0;  // speed
    float species = floor(voronoi(coord, 12.*8.) * 0.99 * float(speciesNum+1) ) / float(speciesNum);
	fragColor = vec4(vec3(species),1.0);  // scale
}
