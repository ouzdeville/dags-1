#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <math.h>
#include "matrix.h"
#include "poly.h"
#include "fichier.h"
#include "param.h"


void polynome_syndrome_1(  binmat_t H, gf * mot, poly_t S);
binmat_t alternant_matrix(binmat_t H, gf * u);
int decoding_H(binmat_t H_alt, gf* c,gf* e,gf* mot );
