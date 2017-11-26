/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                             *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ask permission.                             *
 **********************************************************************************************
 */
#include "decoding.h"

//Bulding of decoding fuction

/*
 * The polynome_syndrome_1 function compute the syndrome S in polynomial form
 *  with inputs a short IV and a Parity matrix H
 */
void polynome_syndrome_1(binmat_t H, gf * mot, poly_t S) {
	int i, j;

	gf tmp;
	for (j = 0; j < H.rown; j++) {
		tmp = 0;
		for (i = 0; i < H.coln; i++) {
			tmp ^= gf_mult(H.coeff[j][i], mot[i]);
		}
		S->coeff[j] = tmp;
	}
	poly_calcule_deg(S);
}

/*
 * The alternant_matrix function transform a Generalized Srivastava cheek
 * matrix H in alternant form
 */
binmat_t alternant_matrix(binmat_t H, gf * u) {
	gf_init(6);
	int i, j, k;
	int st = order * pol_deg;
	poly_t Srivastava;
	poly_t *g, pol1, pol2, temp;
	binmat_t H_alter, H_alt, B, C;;

	//Construction of the First intermediate matrix
//printf("\n   Construction of the First intermediate matrix \n");
	H_alter = mat_ini(st, code_length);
	for (i = 0; i < order; i++) {
		for (k = 0; k < pol_deg; k++) {
			for (j = 0; j < code_length; j++) {
				H_alter.coeff[i * pol_deg + k][j] = H.coeff[k * order + i][j];
			}
		}
	}

	//Construction of Srivastava's Polynome
//printf("\n  Construction of Srivastava's Polynome  \n");
	Srivastava = poly_srivastava(u, order, pol_deg);

	//Construction of the Second intermediate matrix
//printf("\n   Construction of the Second intermediate matrix \n");

	C = mat_ini(st, st);

	g = (poly_t *) malloc(st * sizeof(poly_t));
	for (i = 0; i < st; i++) {
		g[i] = poly_alloc(st - 1);
	}

	pol1 = poly_alloc(st - 1);
	pol2 = poly_alloc(1);

	for (i = 0; i < order; i++) {
		poly_set_to_unit(pol1);
		pol2->coeff[0] = u[i];
		pol2->coeff[1] = 1;
		for (k = 0; k < pol_deg; k++) {
			//Store in temp and free to prevent memory leaks
//			pol1 = poly_mul(pol1, pol2);
			temp = poly_mul(pol1, pol2);
			poly_free(pol1);
			pol1 = temp;
			//g[i * t + k] = poly_quo(Srivastava, pol1);
			 temp = poly_quo(Srivastava, pol1);
			 poly_free(g[i * pol_deg + k]);
			 g[i * pol_deg + k] = temp;
			poly_calcule_deg(g[i * pol_deg + k]);
		}
	}
	poly_free(Srivastava);
	for (i = 0; i < st; i++) {
		for (j = 0; j < g[i]->deg + 1; j++) {
			C.coeff[i][j] = g[i]->coeff[j];
		}
	}
	poly_free(pol1);
	poly_free(pol2);
	for(i = 0; i < st; i++){
		poly_free(g[i]);
	}
	free(g);

	//Construction of the Third intermediate matrix
//printf("\n  Construction of the Third intermediate matrix \n");
	B = mat_ini_Id(st);
	inverse_matrice(C, B);
	mat_free(C);

	//And finally construction of the matrix of Srivastion
//printf("\n And finally construction of the matrix of Srivastion \n");
	H_alt = produit_matrix(B, H_alter);
	//aff_mat(H_alt);
	mat_free(H_alter);
	mat_free(B);
	return H_alt;
}

/*
 * ALTERNANT DECODING: Take as input the matrix in alternate form just created
 * (it was Halt) and the received word c.
 */

int decoding_H(binmat_t H_alt, gf* c, gf* error, gf* code_word) {
	gf_init(6);
	int i, k, j, dr;
	int * LOG_12;
	int st = order * pol_deg;
	poly_t Syndrome;
	poly_t omega, sigma, re, copy_synd, uu, u, quotient, resto, app, temp;
	poly_t pol, pos, error_values;
	gf delta, * ver, pol_gf, tmp, tmp1, o;
	gf alpha;

	//Compute Syndrome normally
	Syndrome = poly_alloc(st - 1);
	polynome_syndrome_1(H_alt, c, Syndrome);
	//aff_poly(Syndrome);

	if (Syndrome->deg == -1) {
		return -1;
	}



	//Resolution of the key equation
	re = poly_alloc(st);
	re->coeff[st] = 1;
	poly_calcule_deg(re);

	uu = poly_alloc(st);
	u = poly_alloc(st);
	poly_set_to_zero(uu);
	poly_set_to_unit(u);
	app = poly_alloc(st);

	copy_synd = poly_copy(Syndrome);//TODO use Syndrome in place of copy_synd
	poly_calcule_deg(copy_synd);
	poly_free(Syndrome);

	resto = poly_copy(re);
	dr = copy_synd->deg;

	while (dr >= (st / 2)) {
		quotient = poly_quo(re, copy_synd);
		poly_rem(resto, copy_synd);

		poly_set(re, copy_synd);
		poly_set(copy_synd, resto);
		poly_set(resto, re);

		poly_set(app, uu);
		poly_set(uu, u);

		temp = poly_mul(u, quotient);
		poly_free(u);
		poly_free(quotient);
		u = temp;
		poly_calcule_deg(u);
		poly_add_free(u, u, app);
		poly_calcule_deg(copy_synd);
		dr = copy_synd->deg;
	}
	poly_free(re);
	poly_free(uu);
	poly_free(app);
	poly_free(resto);

	//Then we find error locator poly (sigma) and error evaluator poly (omega)
	delta = gf_inv(poly_eval(u, 0));
	pol = poly_alloc(0);
	pol->coeff[0] = delta;
	poly_calcule_deg(pol);
	omega = poly_mul(copy_synd, pol);
	sigma = poly_mul(u, pol);
	poly_free(pol);
	poly_free(u);
	poly_free(copy_synd);

	//Support
	ver = (gf*) calloc(code_length, sizeof(gf));
	for (i = 0; i < code_length; i++) {
		ver[i] = gf_mult(H_alt.coeff[1][i], gf_inv(H_alt.coeff[0][i]));
	}

	//Polynome POS gives the position of the errors
	pos = poly_alloc(st / 2);
	j = 0;
	for (i = 0; i < code_length; i++) {
		if ((ver[i] != 0) && (poly_eval(sigma, gf_inv(ver[i])) == 0)) {
			pos->coeff[j] = i;
			j += 1;
		}
	}
	poly_calcule_deg(pos);
	poly_free(sigma);

	//Element for determining the value of errors
	if (pos->deg == -1) {
		return -1;
	}

	// H_alt = produit_matrix(H_alt,PP);
	app = poly_alloc(pos->deg);
	for (j = 0; j <= pos->deg; j++) {
		pol_gf = 1;
		for (i = 0; i <= pos->deg; i++) {
			if (i != j) {
				tmp = gf_mult(ver[pos->coeff[i]],
					      gf_inv(ver[pos->coeff[j]]));
				tmp = gf_add(1, tmp);
				pol_gf = gf_mult(pol_gf, tmp);
			}
		}
		o = poly_eval(omega, gf_inv(ver[pos->coeff[j]]));
		tmp1 = gf_mult(H_alt.coeff[0][pos->coeff[j]], pol_gf);
		app->coeff[j] = gf_mult(o, gf_inv(tmp1));
	}
	poly_calcule_deg(app);
	free(ver);
	poly_free(omega);

	//Determining the value of the errors
	alpha = gf_pow(64, 65);

	LOG_12 = (int*) calloc(gf_card, sizeof(int));
	tmp = 1;
	LOG_12[0] = -1;
	LOG_12[1] = 0;
	for (i = 1; i < gf_card; i++) {
		tmp = gf_mult(tmp, 64);
		LOG_12[tmp] = i;
	}

	k = 0;
	//Reconstruction of the error vector
	for (i = 0; i <= app->deg; i++) {
		j = LOG_12[app->coeff[i]];
		k = j / LOG_12[alpha];
		error[pos->coeff[i]] = gf_Pow_subfield(2, k);
		//printf(" %d " ,valeur_erreurs->coeff[i]);
	}
	poly_free(app);
	free(LOG_12);
	poly_free(pos);

	//Reconstruction of code_word
	for (i = 0; i < code_length; i++) {
		code_word[i] = c[i] ^ error[i];
		//printf(" %d ",code_word[i]);
	}

	return 1;

}

