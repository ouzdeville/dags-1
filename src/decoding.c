/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                             *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ask permission.                             *
 **********************************************************************************************
 */
#include "decoding.h"

//Bulding of decoding fuction

/*************************************************************************************************************
 * The polynome_syndrome_1 function compute the syndrome S in polymial form with inputs a mot and a Parity matrix H                           *
 **************************************************************************************************************
 */
void polynome_syndrome_1(binmat_t H, gf * mot, poly_t S) {
	int i, j;

	gf tmp;
	for (j = 0; j < H.rown; j++) {
		tmp = 0;
		for (i = 0; i < H.coln; i++) {
			tmp ^= gf_Mult(H.coeff[j][i], mot[i]);
		}
		S->coeff[j] = tmp;
	}
	poly_calcule_deg(S);
}

/***************************************************************************************************************
              The alternant_matrix function transform a Generalized Srivastava cheek matrix H in alternant form
***************************************************************************************************************/

binmat_t alternant_matrix(binmat_t H, gf * u) {
	gf_init(6);
	int s = (order);
	int t = pol_deg, i, j, k;
	int st = s * t;
	poly_t Srivastava;
	Srivastava = poly_alloc(st);
	binmat_t H_alter, H_alt;

	/*                                      Construction of the First intermediate matrix
	 ***************************************************************************************************************/
//printf("\n   Construction of the First intermediate matrix \n");
	H_alter = mat_ini(st, code_length);
	H_alt = mat_ini(st, code_length);
	for (i = 0; i < s; i++) {
		for (k = 0; k < t; k++) {
			for (j = 0; j < code_length; j++) {
				H_alter.coeff[i * t + k][j] = H.coeff[k * s + i][j];
			}
		}
	}

	/*                                      Construction of Srivastava's Polynome
	 ***************************************************************************************************************/
//printf("\n  Construction of Srivastava's Polynome  \n");
	Srivastava = poly_srivastava(u, s, t);

	/*                                      Construction of the Second intermediate matrix
	 ***************************************************************************************************************/
//printf("\n   Construction of the Second intermediate matrix \n");
	binmat_t C;
	C = mat_ini(st, st);
	poly_t * g;
	g = (poly_t *) malloc(st * sizeof(poly_t));
	for (i = 0; i < st; i++) {
		g[i] = poly_alloc(st - 1);
	}
	poly_t pol1, pol2;
	pol1 = poly_alloc(st - 1);
	pol2 = poly_alloc(1);

	for (i = 0; i < order; i++) {
		poly_set_to_unit(pol1);
		pol2->coeff[0] = u[i];
		pol2->coeff[1] = 1;
		for (k = 0; k < pol_deg; k++) {
			pol1 = poly_mul(pol1, pol2);
			g[i * t + k] = poly_quo(Srivastava, pol1);
			poly_calcule_deg(g[i * t + k]);
		}
	}
	for (i = 0; i < st; i++) {
		for (j = 0; j < g[i]->deg + 1; j++) {
			C.coeff[i][j] = g[i]->coeff[j];
		}
	}
	poly_free(pol1);
	poly_free(pol2);

	/*                                      Construction of the Third intermediate matrix
	 ***************************************************************************************************************/
//printf("\n  Construction of the Third intermediate matrix \n");
	binmat_t B;
	B = mat_ini_Id(st);
	inverse_matrice(C, B);

	/*                                  And finally construction of the matrix of Srivastion
	 ***************************************************************************************************************/
//printf("\n And finally construction of the matrix of Srivastion \n");
	H_alt = produit_matrix(B, H_alter);
	//aff_mat(H_alt);
	mat_free(H_alter);
	return H_alt;
}

/***************************************************************************************************************
 * ALTERNANT DECODING:                                                                                          *
 * Take as input the matrix in alternant form just created (it was Halt) and the received word c.               *
 ****************************************************************************************************************/

int decoding_H(binmat_t H_alt, gf* c, gf* error, gf* code_word) {
		gf_init(6);
	int i, k, j;
	int s = (order), t = pol_deg;
	int st = s * t;

	/*                                      compute Syndrome normally
	 ***************************************************************************************************************/
	poly_t Syndrome = poly_alloc(st - 1);
	polynome_syndrome_1(H_alt, c, Syndrome);
	//aff_poly(Syndrome);

	if (Syndrome->deg == -1) {
		return -1;
	} else {

		/*                                    Resolution of the key equation
		 ***************************************************************************************************************/

		poly_t omega, sigma, x_st, re, copy_synd, uu, u, quotient, resto, app;
		omega = poly_alloc(st / 2);
		sigma = poly_alloc(st / 2);
		x_st = poly_alloc(st);
		x_st->coeff[st] = 1;
		poly_calcule_deg(x_st);

		re = poly_alloc(st);
		re = poly_copy(x_st);
		poly_calcule_deg(re);

		copy_synd = poly_alloc(st);
		uu = poly_alloc(st);
		u = poly_alloc(st);
		quotient = poly_alloc(st);
		resto = poly_alloc(st);
		app = poly_alloc(st);

		copy_synd = poly_copy(Syndrome);
		poly_calcule_deg(copy_synd);

		poly_set_to_zero(uu);
		poly_set_to_unit(u);

		resto = poly_copy(re);
		int dr = copy_synd->deg;

		while (dr >= (st / 2)) {
			quotient = poly_quo(re, copy_synd);
			poly_rem(resto, copy_synd);

			poly_set(re, copy_synd);
			poly_set(copy_synd, resto);
			poly_set(resto, re);

			poly_set(app, uu);
			poly_set(uu, u);

			u = poly_mul(u, quotient);
			poly_calcule_deg(u);
			poly_add_free(u, u, app);
			poly_calcule_deg(copy_synd);
			dr = copy_synd->deg;

		}

		/*                  Then we find error locator poly (sigma) and error evaluator poly (omega)
		 ***************************************************************************************************************/
		gf delta;
		delta = gf_Inv(poly_eval(u, 0));
		poly_t pol;
		omega = poly_alloc(st / 2);
		sigma = poly_alloc(st / 2);
		pol = poly_alloc(0);
		pol->coeff[0] = delta;
		poly_calcule_deg(pol);
		omega = poly_mul(copy_synd, pol);
		sigma = poly_mul(u, pol);
		poly_calcule_deg(omega);
		poly_calcule_deg(sigma);

		/*                                                 Support
		 ***************************************************************************************************************/
		gf * ver;
		ver = (gf*) calloc(code_length, sizeof(gf));
		for (i = 0; i < code_length; i++) {
			ver[i] = gf_Mult(H_alt.coeff[1][i], gf_Inv(H_alt.coeff[0][i]));
		}

		/*                                   Polynome POS gives the position of the errors
		 ***************************************************************************************************************/
		poly_t pos;
		pos = poly_alloc(st / 2);
		j = 0;
		for (i = 0; i < code_length; i++) {
			if ((ver[i] != 0) && (poly_eval(sigma, gf_Inv(ver[i])) == 0)) {
				pos->coeff[j] = i;
				j += 1;
			}
		}
		poly_calcule_deg(pos);

		/*                                    Element for determining the value of errors
		 ***************************************************************************************************************/
		if (pos->deg == -1) {
			return -1;
		} else {
			// H_alt = produit_matrix(H_alt,PP);
			poly_t app;
			app = poly_alloc(pos->deg);
			gf pol;
			gf tmp, tmp1;
			for (k = 0; k <= pos->deg; k++) {
				j = k;
				pol = 1;
				for (i = 0; i <= pos->deg; i++) {
					if (i != j) {
						tmp = gf_Mult(ver[pos->coeff[i]],
								gf_Inv(ver[pos->coeff[j]]));
						tmp = gf_add(1, tmp);
						pol = gf_Mult(pol, tmp);
					}
				}
				gf o = poly_eval(omega, gf_Inv(ver[pos->coeff[j]]));
				tmp1 = gf_Mult(H_alt.coeff[0][pos->coeff[j]], pol);
				app->coeff[k] = gf_Mult(o, gf_Inv(tmp1));
			}
			poly_calcule_deg(app);

			/*                                   Determining the value of the errors
			 ***************************************************************************************************************/
			gf alpha;
			alpha = gf_Pow(64, 65);

			int * LOG_12;
			LOG_12 = (int*) calloc(gf_card, sizeof(int));
			gf p = 1;
			LOG_12[0] = -1;
			LOG_12[1] = 0;
			for (i = 1; i < gf_card; i++) {

				p = gf_mult(p, 64);
				LOG_12[p] = i;

			}

			poly_t valeur_erreurs = poly_alloc(pos->deg);
			int k = 0;

			for (i = 0; i <= app->deg; i++) {
				j = LOG_12[app->coeff[i]];
				k = j / LOG_12[alpha];
				valeur_erreurs->coeff[i] = gf_Pow_subfield(2, k);
				//printf(" %d " ,valeur_erreurs->coeff[i]);
			}

			/*                                       Reconstruction of the error vector
			 ***************************************************************************************************************/
			for (i = 0; i <= pos->deg; i++) {
				error[pos->coeff[i]] = valeur_erreurs->coeff[i];
				//printf(" %d ",error[i]);
			}

			/*                                        Reconstruction of code_word
			 ***************************************************************************************************************/
			for (i = 0; i < code_length; i++) {
				code_word[i] = c[i] ^ error[i];
				//printf(" %d ",code_word[i]);
			}

			return 1;

		}

	}

}

