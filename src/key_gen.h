
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "decoding.h"
#include "decoding.h"
#include "matrix.h"
#include "fichier.h"
#include "poly.h"
#include "param.h"
#include "util.h"



int disjoint_test(gf * u, gf * v );
int Test_disjoint(gf * L,int n);
void Random_Vect(int m, gf *vect);

void Init_Random_U(gf *U);
void binary_quasi_dyadic_sig(int m, int n, int t, int * b, gf * h_sig, gf * w );
void Cauchy_Support(gf * Support, gf * W,gf * w);
int key_pair(unsigned char *pk, unsigned char *sk);
