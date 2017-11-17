#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gf.h"

/*
 ~~~~~~~~ARITHMETIC FIELD ELEMENT CONSTRUCTION ~~~~~~~~~~~~~~~~

 We define arithmetic field in  F[2][x]/(f), where f is an m-irreducible polynomial.
 In our case, m=5 or m=6.
 For multiplication ,inversion,square,exponentiation field elements, we have adopted the "bitsliced-operation" technic.
 Addition field element, we used XOR between integers.
 */

static void gf_init_antilog_sf() {

	int i = 1;
	gf_antilog_sf = (gf_t *) malloc((1 << gf_extd_sf * sizeof(gf_t)));
	gf_antilog_sf[0] = 1;
	for (i = 1; i < gf_ord_sf; ++i) {
		gf_antilog_sf[i] = gf_antilog_sf[i - 1] << 1;
		if (((gf_antilog_sf[i - 1]) & (1 << (gf_extd_sf - 1)))) {
			gf_antilog_sf[i] ^= poly_primitif_subfield;
		}
	}
	gf_antilog_sf[gf_ord_sf] = 1;
}

static void gf_init_log_sf() {

	int i = 1;
	gf_log_sf = (gf_t *) malloc((1 << gf_extd_sf * sizeof(gf_t)));

	gf_log_sf[0] = -1;
	gf_log_sf[1] = 0;
	for (i = 1; i < gf_ord_sf; ++i) {
		gf_log_sf[gf_antilog_sf[i]] = i;
	}
}

//for F_2^12
/*static void gf_init_antilog() {
 int i=1;
 gf_antilog = (gf *) malloc((1<<gf_extd * sizeof(gf)));
 gf_antilog[0]=1;
 gf_antilog[1]= u_val;
 gf_antilog[2]= 100;
 gf_t x,y,temp;
 x=1;
 y=36;

 for (i=3;i<=gf_ord;++i){
 temp = x;
 x ^= y;
 y = gf_Mult_subfield(temp ,36);
 gf_antilog[i]= (x<<6)^y;

 }
 gf_antilog[gf_ord]=1;

 }*/

static void gf_init_antilog() {
	int i;
	gf_antilog = (gf *) calloc(gf_card, sizeof(gf));
	gf p = 1;
	gf_antilog[0] = 1;
	gf_antilog[1] = 0;
	for (i = 1; i < gf_card; i++) {

		p = gf_Mult1(p, 64);
		gf_antilog[i] = p;

	}

}

//for F_2^12
static void gf_init_log() {

	int i;
	gf_log = (gf *) malloc((1 << gf_extd * sizeof(gf)));
	gf_log[0] = -1;
	gf_log[1] = 0;
	for (i = 1; i < gf_ord; ++i) {
		gf_log[gf_antilog[i]] = i;
	}
}

gf gf_diff1(gf a, gf b) {
	uint32_t t = (uint32_t) (a ^ b);
	t = ((t - 1) >> 20) ^ 0xFFF;
	return (gf) t;
}

gf gf_Div1(gf a, gf b) {
	if (b == 0) {
		fprintf(stderr, "ERROR %d is not invertible", b);
		exit(-1);
	} else {
		gf res = gf_Mult(a, gf_Inv(b));
		return res;
	}
}

gf gf_Pow1(gf in, int n) {

	gf h, t;
	h = 1;
	t = in;
	while (n != 0) {
		if (n % 2 == 1) {
			h = gf_Mult(h, t);
		}
		n = (n / 2);
		t = gf_Mult(t, t);
	}
	return h;
}

gf gf_Mult1(gf x, gf y) {
	gf a1, b1, a2, b2, a3, b3;
	a1 = x >> 6;
	b1 = x & 63;
	a2 = y >> 6;
	b2 = y & 63;
	a3 = gf_Mult_subfield(gf_Mult_subfield(a1, a2),
			36)^gf_Mult_subfield(a1, b2)^gf_Mult_subfield(b1, a2);
	b3 =
			gf_Mult_subfield(gf_Mult_subfield(a1, a2),
					2)^gf_Mult_subfield(b1, b2);
	return (a3 << 6) ^ b3;
}

gf gf_sq1(gf x) {
	gf a1, b1, a3, b3;
	a1 = x >> 6;
	b1 = x & 63;

	a3 = gf_Mult_subfield(gf_Mult_subfield(a1, a1), 36);
	b3 =
			gf_Mult_subfield(gf_Mult_subfield(a1, a1),
					2)^gf_Mult_subfield(b1, b1);
	return (a3 << 6) ^ b3;
}

gf gf_Inv1(gf in) {
	gf tmp_11;
	gf tmp_1111;

	gf out = in;

	out = gf_sq1(out);
	tmp_11 = gf_Mult1(out, in);

	out = gf_sq1(tmp_11);
	out = gf_sq1(out);
	tmp_1111 = gf_Mult1(out, tmp_11);

	out = gf_sq1(tmp_1111);
	out = gf_sq1(out);
	out = gf_sq1(out);
	out = gf_sq1(out);
	out = gf_Mult1(out, tmp_1111);

	out = gf_sq1(out);
	out = gf_sq1(out);
	out = gf_Mult1(out, tmp_11);

	out = gf_sq1(out);
	out = gf_Mult1(out, in);

	return gf_sq1(out);
}

int init_done = 0;

int gf_init(int extdeg) {
	if (extdeg > gf_extd_sf) {

		exit(0);
	}

	if (init_done != extdeg) {
		if (init_done) {
			free(gf_antilog_sf);
			free(gf_log_sf);
			free(gf_antilog);
			free(gf_log);
		}
		init_done = extdeg;
		gf_init_antilog_sf();
		gf_init_log_sf();
		gf_init_antilog();
		gf_init_log();
	}

	return 1;
}

gf_t gf_rand(int (*u8rnd)()) {
	return (u8rnd() ^ (u8rnd() << 8)) & gf_ord_sf;
}

