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

#include "gf_tables.h"
#include "types_def.h"
#include "param.h"

//int gf_extension_degree, gf_cardinality, gf_multiplicative_order;

gf *gf_log;
gf *gf_antilog;

#define gf_unit() 1
#define gf_zero() 0
#define gf_add(x, y) ((x) ^ (y)) // Addition in the field

////////////////////////////////////////////////////////////////////
///////////////////////// SUBFIELD OPERATION ///////////////////////

/* we obtain a value between 0 and (q-1) included, the class of 0 is
 represented by 0 or q-1 (this is why we write _K->exp[q-1]=_K->exp[0]=1)*/



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

// In correct gf_Pow1
gf gf_Pow1(gf f, int n);

// Propose gf_Ppw
gf gf_pow(gf f, int n);

// Un-nessary use antilog and log table
int gf_init(int extdeg);

#endif /* GF_H_ */
