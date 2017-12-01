/*
 * gf_tables.h
 *
 *  Created on: Dec 1, 2017
 *      Author: vader
 */

#ifndef SRC_GF_TABLES_H_
#define SRC_GF_TABLES_H_

#include <stdint.h>
#include "types_def.h"


gf_t gf_log_sf[64];
gf_t gf_antilog_sf[64];

#define gf_modq_1_sf(d) ((d) % 63)

#define gf_mul_fast_subfield(x, y) ((y) ? gf_antilog_sf[gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])] : 0)

// Multiplication in the field : apha^i*alpha^j=alpha^(i+j)
// Check x is zero, if zero, return 0, else calculate
#define gf_Mult_subfield(x, y) ((x) ? gf_mul_fast_subfield(x, y) : 0)

// In direct way to calculate power in range 2^6.
// Only use in line
// 404:decoding.c: 				valeur_erreurs->coeff[i] = gf_Pow_subfield(2, k);
// gf_Pow_subfield is always calculate 2^k
#define gf_Pow_subfield(x, i) (gf_antilog_sf[(gf_modq_1_sf(gf_log_sf[x] * i))])

// Inverse in the subfield
#define gf_Inv_subfield(x) gf_antilog_sf[gf_ord_sf - gf_log_sf[x]]

#define gf_mul_fast(x, y) ((y) ? gf_antilog[_gf_modq_1(gf_log[x] + gf_log[y])] : 0)
#endif /* SRC_GF_TABLES_H_ */
