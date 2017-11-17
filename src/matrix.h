#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gf.h"
#include "param.h"

typedef struct matrix{
 	unsigned int rown;
	unsigned int coln;
	gf **coeff;
} binmat_t;


#define mat_coeff(A,i,j) ((coeffe[i][j]))


int* test_mat(binmat_t A);
binmat_t mat_ini(int rown, int coln);
void mat_free(binmat_t A);
binmat_t mat_copy(binmat_t A);
void mat_rowxor(binmat_t A, int i, int j);
void mat_swaprow(binmat_t A, int i, int j);
void mat_swapcol(binmat_t A, int i, int j);
void mat_random_swapcol(binmat_t A);



void G_mat(binmat_t G,binmat_t H_syst);
binmat_t mat_Into_base(binmat_t H );

int syst_mat(binmat_t H);
void mat_line_mult_by_gf(binmat_t A, gf a, int i);

void mat_rowxor_with_another(binmat_t A,int i,gf * Line);


//int syst(binmat_t H, binmat_t P);
int syst(binmat_t H);
binmat_t mat_ini_Id(int rown);
binmat_t  Inversion_mat_pivot(binmat_t A);




void affiche_vecteur(gf * P, int taille);


binmat_t transpose(binmat_t A);
binmat_t produit_matrix(binmat_t A, binmat_t B);
binmat_t produit_matrix_subfield(binmat_t A, binmat_t B);
gf* produit_matrix_vector_subfield(binmat_t A, gf* v);
void aff_mat(binmat_t mat);
gf* produit_vector_matrix_subfield(gf* v,binmat_t A);
binmat_t permut_mat( int * P);
gf* produit_matrix_vector(binmat_t A, gf* v);
gf* produit_vector_matrix(gf* v,binmat_t A);
gf* produit_vector_matrix_Sf(gf* v,binmat_t A);
void Permut_vecteur(int * P, gf* v, int taille);
gf* produit_matrix_vector(binmat_t A, gf* v);

int inverse_matrice(binmat_t A,binmat_t S);

gf* Eltseq(gf a);
void secret_matrix(binmat_t H,gf * u, gf * v,gf * z);
void quasi_dyadic_bloc_mat(int s,binmat_t M,gf * sig,int ind_col,int ind_rown);

#endif

