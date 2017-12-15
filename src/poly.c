#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"

poly_t poly_alloc(int d) {
	poly_t p;
	p = (poly_t) malloc(sizeof(struct polynomial));
	p->deg = -1;
	p->size = d + 1;
	p->coeff = (gf *) calloc(p->size, sizeof(gf));
	return p;
}

poly_t poly_copy(poly_t p) {
	poly_t q;
	q = (poly_t) malloc(sizeof(struct polynomial));
	q->deg = p->deg;
	q->size = p->size;
	q->coeff = (gf *) calloc(q->size, sizeof(gf));
	memcpy(q->coeff, p->coeff, p->size * sizeof(gf));
	return q;
}

void poly_free(poly_t p) {
	free(p->coeff);
	free(p);
}

void poly_set_to_zero(poly_t p) {
	memset(p->coeff, 0, p->size * sizeof(gf));
	p->deg = -1;
}

poly_t poly_set_to_null() {
	poly_t p;
	p = poly_alloc(0);
	p->coeff[0] = 0;
	p->deg = -1;
	return p;
}
void poly_set_to_unit(poly_t p) {
	memset(p->coeff, 0, p->size * sizeof(gf));
	p->coeff[0] = 1;
	p->deg = 0;
}

//Compute the maximal degree of a given polynomial
int poly_compute_deg(poly_t p) {
	int d = p->size - 1;
	while ((d >= 0) && (p->coeff[d] == 0))
		--d;
	p->deg = d;
	return d;
}
//Set the polynomial p which is equal to q
void poly_set(poly_t p, poly_t q) {
	int d = p->size - q->size;
	if (d < 0) {
		memcpy(p->coeff, q->coeff, p->size * sizeof(gf));
		poly_compute_deg(p);
	} else {
		memcpy(p->coeff, q->coeff, q->size * sizeof(gf));
		memset(p->coeff + q->size, 0, d * sizeof(gf));
		poly_compute_deg(q);
		p->deg = q->deg;
	}
}
//Used by poly_eval
//TODO: Generate constant time for this function
gf poly_eval_aux(gf * coeff, gf a, int d) {
	gf b;
	b = coeff[d--];
	for (; d >= 0; --d) {
		if (b != 0) {
			b = (gf_mult(b, a) ^ coeff[d]);
		} else {
			b = coeff[d];
		}
	}
	return b;
}

//Returns the multiplication of p by q
poly_t poly_mul(poly_t p, poly_t q) {
	int i, j, dp, dq;
	poly_t r;

	poly_compute_deg(p);
	poly_compute_deg(q);
	dp = poly_deg(p);
	dq = poly_deg(q);
	r = poly_alloc(dp + dq);
	for (i = 0; i <= dp; ++i) {
		for (j = 0; j <= dq; ++j) {
			poly_addto_coeff(r, i + j, gf_mult(poly_coeff(p,i),poly_coeff(q,j)));
		}
	}
	poly_compute_deg(r);

	return r;
}


//TODO: implement in constant time : remove the if/else
//Returns the addition of p by q in r (used to avoid memory leak)
void poly_add_free(poly_t r, poly_t a, poly_t b) {
	int i;
	if (a->deg == -1) {
		r->deg = b->deg;
		memcpy(r->coeff, b->coeff, b->deg * sizeof(gf));
	} else if (b->deg == -1) {
		r->deg = a->deg;
		memcpy(r->coeff, a->coeff, a->deg * sizeof(gf));
	} else {
		if (a->deg == b->deg) {
			r->deg = a->deg;
			for (i = 0; i < a->deg + 1; i++) {
				r->coeff[i] = (a->coeff[i]) ^ (b->coeff[i]);
			}
		}
		if (a->deg > b->deg) {
			r->deg = a->deg;
			for (i = 0; i < b->deg + 1; i++)
				r->coeff[i] = (a->coeff[i]) ^ (b->coeff[i]);
			memcpy(r->coeff + (b->deg + 1), a->coeff + (b->deg + 1),
					(a->deg - b->deg) * sizeof(gf));
		}
		if (b->deg > a->deg) {

			r->deg = b->deg;
			for (i = 0; i < a->deg + 1; i++)
				r->coeff[i] = (a->coeff[i]) ^ (b->coeff[i]);
			memcpy(r->coeff + (a->deg + 1), b->coeff + (a->deg + 1),
					(b->deg - a->deg) * sizeof(gf));
		}
	}
}

gf poly_eval(poly_t p, gf a) {
	poly_compute_deg(p);
	return poly_eval_aux(p->coeff, a, p->deg);
}

void poly_rem(poly_t p, poly_t g) {
	int i, j, d;
	gf a, b;
	poly_compute_deg(p);
	poly_compute_deg(g);
	d = p->deg - g->deg;
	if (d >= 0) {
		a = gf_inv(poly_tete(g));
		for (i = p->deg; d >= 0; --i, --d) {
			if (poly_coeff(p, i) != 0) {
				b = gf_mult(a, poly_coeff(p, i));
				for (j = 0; j < g->deg; ++j) {
					poly_addto_coeff(p, j + d, gf_mult(b, poly_coeff(g, j))); //In F2^m, addition=soustraction
				}
				poly_set_coeff(p, i, 0);
			}
		}
		poly_set_deg(p, g->deg - 1);
		while ((p->deg >= 0) && (poly_coeff(p, p->deg) == 0)) {
			poly_set_deg(p, p->deg - 1);
		}
	}
}

poly_t poly_quo(poly_t p, poly_t d) {
	int i, j, dd, dp;
	gf a, b;
	poly_t quo, rem;

	dd = poly_compute_deg(d);
	dp = poly_compute_deg(p);
	rem = poly_copy(p);
	quo = poly_alloc(dp - dd);
	quo->deg = dp - dd;
	a = gf_inv(poly_coeff(d, dd));
	for (i = dp; i >= dd; --i) {
		b = gf_mult(a, poly_coeff(rem, i));
		quo->coeff[i - dd] = b;
		if (b != 0) {
			rem->coeff[i] = 0;
			for (j = i - 1; j >= i - dd; --j) {
				poly_addto_coeff(rem, j, gf_mult(b, poly_coeff(d, dd - i + j)));
			}
		}
	}
	poly_free(rem);
	return quo;
}


poly_t poly_srivastava(gf * W, int s, int t) {
	poly_t poly_de_goppa_srivasta, poly1, temp;
	int i, j;

	poly_de_goppa_srivasta = poly_alloc(s * t);
	poly_set_to_zero(poly_de_goppa_srivasta);
	poly_de_goppa_srivasta->coeff[0] = 1;

	poly1 = poly_alloc(1);
	poly_set_to_zero(poly1);

	for (i = 0; i < s; i++) {
		poly1->coeff[1] = 1;
		poly1->coeff[0] = W[i];
		for (j = 0; j < t; j++) {
			//poly_de_goppa_srivasta = poly_mul (poly_de_goppa_srivasta, poly1);
			temp = poly_mul(poly_de_goppa_srivasta, poly1);
			poly_free(poly_de_goppa_srivasta);
			poly_de_goppa_srivasta = temp;
		}
	}

	poly_free(poly1);
	return poly_de_goppa_srivasta;
}
