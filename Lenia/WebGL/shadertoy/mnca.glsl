//	----    ----    ----    ----    ----    ----    ----    ----
//
//	Shader developed by Slackermanz:
//
//		https://slackermanz.com
//
//		Discord:	Slackermanz#3405
//		Github:		https://github.com/Slackermanz
//		Twitter:	https://twitter.com/slackermanz
//		YouTube:	https://www.youtube.com/c/slackermanz
//		Shadertoy: 	https://www.shadertoy.com/user/SlackermanzCA
//		Reddit:		https://old.reddit.com/user/slackermanz
//
//		Communities:
//			Reddit:		https://old.reddit.com/r/cellular_automata
//			Discord:	https://discord.com/invite/J3phjtD
//			Discord:	https://discord.gg/BCuYCEn
//
//	----    ----    ----    ----    ----    ----    ----    ----

#define txdata (iChannel0)
#define PI 3.14159265359

const uint MAX_RADIUS = 10u;

uint u32_upk(uint u32, uint bts, uint off) { return (u32 >> off) & ((1u << bts)-1u); }

float  tp(uint n, float s) 			{ return (float(n+1u)/256.0) * ((s*0.5)/128.0); }
float bsn(uint v, uint  o) 			{ return float(u32_upk(v,1u,o)*2u)-1.0; }
float vwm()							{ return 255.000336 * 1.0; }
float utp(uint v, uint  w, uint o) 	{ return tp(u32_upk(v,w,w*o), vwm()); }
    
vec4  gdv( ivec2 of, sampler2D tx ) {
	of 		= ivec2(gl_FragCoord) + of;
	of[0] 	= (of[0] + textureSize(tx,0)[0]) % (textureSize(tx,0)[0]);
	of[1] 	= (of[1] + textureSize(tx,0)[1]) % (textureSize(tx,0)[1]);
	return 	texelFetch( tx, of, 0); }
    
vec4 nbhd( vec2 r, sampler2D tx ) {
//	Precision limit of signed float32 for [n] neighbors in a 16 bit texture (symmetry preservation)
	uint	chk = 2147483648u /
			(	( 	uint( r[0]*r[0]*PI + r[0]*PI + PI	)
				- 	uint( r[1]*r[1]*PI + r[1]*PI		) ) * 128u );
	float	psn = (chk >= 65536u) ? 65536.0 : float(chk);
	vec4	a = vec4(0.0,0.0,0.0,0.0);
	for(float i = 0.0; i <= r[0]; i++) {
		for(float j = 1.0; j <= r[0]; j++) {
			float	d = round(sqrt(i*i+j*j));
			float	w = 1.0;
			if( d <= r[0] && d > r[1] ) {
				vec4 t0  = gdv( ivec2( i, j), tx ) * w * psn; a += t0 - fract(t0);
				vec4 t1  = gdv( ivec2( j,-i), tx ) * w * psn; a += t1 - fract(t1);
				vec4 t2  = gdv( ivec2(-i,-j), tx ) * w * psn; a += t2 - fract(t2);
				vec4 t3  = gdv( ivec2(-j, i), tx ) * w * psn; a += t3 - fract(t3); } } }
	return a; }
    
vec4 totl( vec2 r, sampler2D tx ) {
//	Precision limit of signed float32 for [n] neighbors in a 16 bit texture (symmetry preservation)
	uint	chk = 2147483648u /
			(	( 	uint( r[0]*r[0]*PI + r[0]*PI + PI	)
				- 	uint( r[1]*r[1]*PI + r[1]*PI		) ) * 128u );
	float	psn = (chk >= 65536u) ? 65536.0 : float(chk);
	vec4 	b = vec4(0.0,0.0,0.0,0.0);
	for(float i = 0.0; i <= r[0]; i++) {
		for(float j = 1.0; j <= r[0]; j++) {
			float	d = round(sqrt(i*i+j*j));
			float	w = 1.0;
			if( d <= r[0] && d > r[1] ) { b	+= w * psn * 4.0; } } }
	return b; }
                
vec4 bitring(vec4[MAX_RADIUS] rings_a, vec4[MAX_RADIUS] rings_b, uint bits, uint of) {
	vec4 sum = vec4(0.0,0.0,0.0,0.0);
	vec4 tot = vec4(0.0,0.0,0.0,0.0);
	for(uint i = 0u; i < MAX_RADIUS; i++) {
		if(u32_upk(bits, 1u, i+of) == 1u) { sum += rings_a[i]; tot += rings_b[i]; } }
	return sum / tot; }
    
//	Used to reseed the surface with lumpy noise
float get_xc(float x, float y, float xmod) {
	float sq = sqrt(mod(x*y+y, xmod)) / sqrt(xmod);
	float xc = mod((x*x)+(y*y), xmod) / xmod;
	return clamp((sq+xc)*0.5, 0.0, 1.0); }
float shuffle(float x, float y, float xmod, float val) {
	val = val * mod( x*y + x, xmod );
	return (val-floor(val)); }
float get_xcn(float x, float y, float xm0, float xm1, float ox, float oy) {
	float  xc = get_xc(x+ox, y+oy, xm0);
	return shuffle(x+ox, y+oy, xm1, xc); }
float get_lump(float x, float y, float nhsz, float xm0, float xm1) {
	float 	nhsz_c 	= 0.0;
	float 	xcn 	= 0.0;
	float 	nh_val 	= 0.0;
	for(float i = -nhsz; i <= nhsz; i += 1.0) {
		for(float j = -nhsz; j <= nhsz; j += 1.0) {
			nh_val = round(sqrt(i*i+j*j));
			if(nh_val <= nhsz) {
				xcn = xcn + get_xcn(x, y, xm0, xm1, i, j);
				nhsz_c = nhsz_c + 1.0; } } }
	float 	xcnf 	= ( xcn / nhsz_c );
	float 	xcaf	= xcnf;
	for(float i = 0.0; i <= nhsz; i += 1.0) {
			xcaf 	= clamp((xcnf*xcaf + xcnf*xcaf) * (xcnf+xcnf), 0.0, 1.0); }
	return xcaf; }
float reseed(int seed) {
	vec4	fc = gl_FragCoord;
	float 	r0 = get_lump(fc[0], fc[1],  2.0, 19.0 + mod(iDate[3]+float(seed),17.0), 23.0 + mod(iDate[3]+float(seed),43.0));
	float 	r1 = get_lump(fc[0], fc[1], 14.0, 13.0 + mod(iDate[3]+float(seed),29.0), 17.0 + mod(iDate[3]+float(seed),31.0));
	float 	r2 = get_lump(fc[0], fc[1],  6.0, 13.0 + mod(iDate[3]+float(seed),11.0), 51.0 + mod(iDate[3]+float(seed),37.0));
	return clamp((r0+r1)-r2,0.0,1.0); }

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {

    //	Parameters
	const	float 	mnp 	= 1.0 / 65536.0;			//	Minimum value of a precise step for 16-bit channel
	const	float 	s  		= mnp *  96.0 *  96.0;
	const	float 	n  		= mnp *  96.0 *  16.0;
    
	vec4 res_c = gdv( ivec2(0, 0), txdata );
    
    //	NH Rings
	vec4[MAX_RADIUS] nh_rings_c_a;
	vec4[MAX_RADIUS] nh_rings_c_b;
	for(uint i = 0u; i < MAX_RADIUS; i++) {
		nh_rings_c_a[i] = nbhd( vec2(i+1u,i), txdata );
		nh_rings_c_b[i] = totl( vec2(i+1u,i), txdata ); }
        
    uint[12] NB = uint[12] (
		1857775778u, 373976652u,  2050229723u, 3957459881u, 
		2971170932u, 179461069u,  3501247727u, 1512432443u, 
		956981535u,  2382362945u, 168996476u,  2033912016u  );
        
    uint[24] UD = uint[24] (
		69206018u,   887201330u,  529408915u,  4175580383u, 
		394674196u,  1588503026u, 90701929u,   2284734815u, 
		1556401948u, 2673158198u, 4228048446u, 716239805u, 
		1175339528u, 4070697602u, 3108524504u, 4271215171u, 
		1366665979u, 2104054760u, 3123817950u, 2010068087u, 
		2867880264u, 2061177276u, 3179780908u, 683022191u   );
        
    uint[ 2] SN = uint[ 2] (
		1631696705u, 2329356951u                            );

	uint[ 3] IO = uint[ 3] (
		2286157824u, 295261525u, 1713547946u                );


    for(uint i = 0u; i < 24u; i++) {
		float nhv = bitring( nh_rings_c_a, nh_rings_c_b, NB[i/2u], (i & 1u) * 16u )[u32_upk( IO[i/8u], 2u, (i*4u+0u) & 31u )];
		if( nhv >= utp( UD[i], 8u, 0u) && nhv <= utp( UD[i], 8u, 1u)) { res_c[u32_upk( IO[i/8u], 2u, (i*4u+2u) & 31u )] += bsn(SN[i/16u], ((i*2u+0u) & 31u))*s; }
		if( nhv >= utp( UD[i], 8u, 2u) && nhv <= utp( UD[i], 8u, 3u)) { res_c[u32_upk( IO[i/8u], 2u, (i*4u+2u) & 31u )] += bsn(SN[i/16u], ((i*2u+1u) & 31u))*s; } }
    
	res_c -= n;
    res_c  = clamp(res_c,0.0,1.0);
    
//	----    ----    ----    ----    ----    ----    ----    ----
//	Shader Output
//	----    ----    ----    ----    ----    ----    ----    ----

    if (iMouse.z > 0. && length(iMouse.xy - fragCoord) < 14.0) {
        res_c[0] = round(mod(float(iFrame),2.0));
        res_c[1] = round(mod(float(iFrame),3.0));
        res_c[2] = round(mod(float(iFrame),5.0)); }
    if (iFrame == 0) { res_c[0] = reseed(0); res_c[1] = reseed(1); res_c[2] = reseed(2); }
    fragColor=res_c;
}

