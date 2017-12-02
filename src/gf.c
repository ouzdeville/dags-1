#include "gf.h"

/*
 ~~~~~~~~ARITHMETIC FIELD ELEMENT CONSTRUCTION ~~~~~~~~~~~~~~~~
 We define arithmetic field in  F[2][x]/(f), where f is an m-irreducible polynomial.
 In our case, m=5 or m=6.
 For multiplication ,inversion,square,exponentiation field elements, we have adopted the "bitsliced-operation" technic.
 Addition field element, we used XOR between integers.
 */

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
	b1 = x & (u_val - 1);
	a2 = y >> gf_extd_sf;
	b2 = y & (u_val - 1);

	a3 = gf_mult_fast(gf_mult_fast(a1, a2),
			36) ^ gf_mult_fast(a1, b2) ^ gf_mult_fast(b1, a2);

	b3 = gf_mult_fast(gf_mult_fast(a1, a2), 2) ^ gf_mult_fast(b1, b2);

	return (a3 << gf_extd_sf) ^ b3;
}

// Correct gf_sq
gf gf_sq(gf x) {
	gf a1, b1, a3, b3;

	a1 = x >> gf_extd_sf;
	b1 = x & (u_val - 1);

	a3 = gf_mult_fast(gf_mult_fast(a1, a1), 36);

	b3 = gf_mult_fast(gf_mult_fast(a1, a1), 2) ^ gf_mult_fast(b1, b1);

	return (a3 << gf_extd_sf) ^ b3;
}

// Correct gf_Inv
gf gf_inv(gf in) {
	gf tmp_11;
	gf tmp_1111;

	gf out = in;
	out = gf_sq(out); //a^2
	tmp_11 = gf_mult(out, in); //a^2*a = a^3

	out = gf_sq(tmp_11); //(a^3)^2 = a^6
	out = gf_sq(out); // (a^6)^2 = a^12
	tmp_1111 = gf_mult(out, tmp_11); //a^12*a^3 = a^15

	out = gf_sq(tmp_1111); //(a^15)^2 = a^30
	out = gf_sq(out); //(a^30)^2 = a^60
	out = gf_sq(out); //(a^60)^2 = a^120
	out = gf_sq(out); //(a^120)^2 = a^240
	out = gf_mult(out, tmp_1111); //a^240*a^15 = a^255

	out = gf_sq(out); // (a^255)^2 = 510
	out = gf_sq(out); //(a^510)^2 =  1020
	out = gf_mult(out, tmp_11); //a^1020*a^3 = 1023

	out = gf_sq(out); //(a^1023)^2 = 2046
	out = gf_mult(out, in); //a^2046*a = 2047
	//gf t = gf_sq(out); //(a^2047)^2 = 4094
	/*
	 gf tmp = gf_pow(in, 4094);
	 //gf tmp = gf_pow(in, 4094);*/
	return gf_sq(out);
}

