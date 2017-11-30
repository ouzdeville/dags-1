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

gf_t gf_log_sf[64] = { 63, 0, 1, 6, 2, 12, 7, 26, 3, 32, 13, 35, 8, 48, 27, 18,
		4, 24, 33, 16, 14, 52, 36, 54, 9, 45, 49, 38, 28, 41, 19, 56, 5, 62, 25,
		11, 34, 31, 17, 47, 15, 23, 53, 51, 37, 44, 55, 40, 10, 61, 46, 30, 50,
		22, 39, 43, 29, 60, 42, 21, 20, 59, 57, 58 };
gf_t gf_antilog_sf[64] = { 1, 2, 4, 8, 16, 32, 3, 6, 12, 24, 48, 35, 5, 10, 20,
		40, 19, 38, 15, 30, 60, 59, 53, 41, 17, 34, 7, 14, 28, 56, 51, 37, 9,
		18, 36, 11, 22, 44, 27, 54, 47, 29, 58, 55, 45, 25, 50, 39, 13, 26, 52,
		43, 21, 42, 23, 46, 31, 62, 63, 61, 57, 49, 33 };

static void gf_init_antilog() {
	/*
	 MAINFIELD table.
	 Build table for faster calculation.
	 In memory access is faster than calculating.
	 */
	int i = 1;
	gf_antilog = (gf *) malloc(gf_card * sizeof(gf));
	gf p = 1;
	gf_antilog[0] = 1;
	// gf_card = 4096
	for (i = 1; i < gf_card; i++) {
		p = gf_mult(p, 64);
		gf_antilog[i] = p;
	}
}

static void gf_init_log() {
	/*
	 MAINFIELD table.
	 Build table for faster calculation.
	 In memory access is faster than calculating.
	 */
	int i = 1;
	gf_log = (gf *) malloc((gf_card * sizeof(gf)));
	gf_log[0] = -1;
	gf_log[1] = 0;
	for (i = 1; i < gf_ord; ++i) {
		gf_log[gf_antilog[i]] = i;
	}
}

gf_t gf_mul_fast_subfield(gf_t x, gf_t y) {
	if (y) {
		return gf_antilog_sf[gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])];
	}
	return 0;
}
gf_t gf_Mult_subfield(gf_t x, gf_t y) {
	if (x) {
		return gf_mul_fast_subfield(x, y);
	}
	return 0;
}
gf_t gf_Pow_subfield(gf_t x, int i) {

	return gf_antilog_sf[(gf_modq_1_sf(gf_log_sf[x] * i))];
}
gf_t gf_Inv_subfield(gf_t x) {
	return (gf_antilog_sf[gf_ord_sf - gf_log_sf[x]]);
}

gf gf_mul_fast(gf x, gf y) {
	if (y) {
		return gf_antilog[_gf_modq_1(gf_log[x] + gf_log[y])];
	}
	return 0;
}
// Correct gf_Div
// Use in poly.c
gf gf_div(gf a, gf b) {
	if (b == 0) {
		fprintf(stderr, "ERROR %d is not invertible", b);
		exit(-1);
	} else {
		gf res = gf_mult(a, gf_inv(b));
		return res;
	}
}

// Correct gf_Pow
gf gf_pow(gf in, int n) {

	gf h, t;
	h = 1;
	t = in;
	while (n != 0) {
		if ((n & 1) == 1) {
			h = gf_mult(h, t);
		}
		n = n >> 1;
		t = gf_mult(t, t);
	}
	return h;
}

// Correct gf_mult
gf gf_mult(gf x, gf y) {
	gf a1, b1, a2, b2, a3, b3;

	a1 = x >> gf_extd_sf;
	b1 = x & (u_val-1);
	a2 = y >> gf_extd_sf;
	b2 = y & (u_val-1);

	a3 = gf_Mult_subfield(gf_Mult_subfield(a1, a2), 36)
			^ gf_Mult_subfield(a1, b2) ^ gf_Mult_subfield(b1, a2);

	b3 = gf_Mult_subfield(gf_Mult_subfield(a1, a2), 2)
			^ gf_Mult_subfield(b1, b2);

	return (a3 << gf_extd_sf) ^ b3;
}

// Correct gf_sq
gf gf_sq(gf x) {
	gf a1, b1, a3, b3;

	a1 = x >> gf_extd_sf;
	b1 = x & (u_val-1);

	a3 = gf_Mult_subfield(gf_Mult_subfield(a1, a1), 36);

	b3 = gf_Mult_subfield(gf_Mult_subfield(a1, a1), 2)
			^ gf_Mult_subfield(b1, b1);

	return (a3 << gf_extd_sf) ^ b3;
}

// Correct gf_Inv
gf gf_inv(gf in) {
	gf tmp_11;
	gf tmp_1111;

	gf out = in;

	out = gf_sq(out);
	tmp_11 = gf_mult(out, in);

	out = gf_sq(tmp_11);
	out = gf_sq(out);
	tmp_1111 = gf_mult(out, tmp_11);

	out = gf_sq(tmp_1111);
	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_mult(out, tmp_1111);

	out = gf_sq(out);
	out = gf_sq(out);
	out = gf_mult(out, tmp_11);

	out = gf_sq(out);
	out = gf_mult(out, in);

	return gf_sq(out);
}

int init_done = 0;

int gf_init(int extdeg) {
	if (extdeg > gf_extd_sf) {

		exit(0);
	}

	if (init_done != extdeg) {
		if (init_done) {
			free(gf_antilog);
			free(gf_log);
		}
		gf_init_antilog();

		gf_init_log();
	}

	return 1;
}

