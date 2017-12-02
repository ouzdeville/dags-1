/*
 * gf.h
 *
 *  Created on: Dec 1, 2017
 *      Author: vader
 */

#ifndef GF_H_
#define GF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#include "gf_tables.h"
#include "types_def.h"
#include "param.h"

//int gf_extension_degree, gf_cardinality, gf_multiplicative_order;

//gf *gf_log;
//gf *gf_antilog;

#define gf_modq_1_sf(d) ((d) % (u_val-1))

#define _gf_modq_1(d) ((d) % (gf_ord))

#define gf_mult_fast_subfield(x, y) ((y) ? gf_antilog_sf[gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])] : 0)

// Multiplication in the field : apha^i*alpha^j=alpha^(i+j)
// Check x is zero, if zero, return 0, else compute
#define gf_mult_fast(x, y) ((x) ? gf_mult_fast_subfield(x, y) : 0)

// gf_Pow_subfield is always calculate 2^k
#define gf_pow_subfield(x, i) (gf_antilog_sf[(gf_modq_1_sf(gf_log_sf[x] * i))])

// Inverse in the subfield
#define gf_inv_subfield(x) gf_antilog_sf[gf_ord_sf - gf_log_sf[x]]

#define gf_mul_fast(x, y) ((y) ? gf_antilog[_gf_modq_1(gf_log[x] + gf_log[y])] : 0)

////////////////////////////////////////////////////////////////////
///////////////////////// MAIN FIELD OPERATION /////////////////////
// Check y is zero, if zero, return 0, else calculate

// Correct gf_Mult1 =>> will rename to gf_Mult
gf gf_mult(gf in0, gf in1);

// Correct gf_Inv1 =>> will rename to gf_Inv
// Use in poly, matrix, keygen, decoding
gf gf_inv(gf in);

// Correct gf_sq1 =>> will rename to gf_sq
gf gf_sq(gf in);

// Incorrect gf_Div1
//gf gf_Div1(gf a, gf b);

// Propose gf_Div2
gf gf_div(gf a, gf b);

// Propose gf_Ppw
gf gf_pow(gf f, int n);

#endif /* GF_H_ */

