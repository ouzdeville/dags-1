/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                            *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.                             *
 **********************************************************************************************
 */

#include "encapsulation.h"

int encapsulation(const unsigned char *pk, unsigned char *ct, unsigned char *ss) {
	gf_init(6);

	unsigned char *m, *d, *rho, *sigma, *u, *e, *hash_sigma, *r, *m_extend,
			*sigma_extend;
	gf *c, *e1, *u1, *dd;

	int i;

	/* Memory's allocation.............................................................................
	 *************************************************************************************************/
	m = (unsigned char*) calloc(k_prime, sizeof(unsigned char));
	rho = (unsigned char*) calloc(k_sec, sizeof(unsigned char));
	sigma = (unsigned char*) calloc(code_dimension - k_sec,
			sizeof(unsigned char));

	hash_sigma = (unsigned char*) calloc(code_length, sizeof(unsigned char));

	//x=( unsigned char*)calloc(400,sizeof( unsigned char));
	u = (unsigned char*) calloc(code_dimension, sizeof(unsigned char));
	u1 = (gf*) calloc(code_dimension, sizeof(gf));
	r = (unsigned char*) calloc(code_dimension, sizeof(unsigned char));
	m_extend = (unsigned char*) calloc(code_dimension, sizeof(unsigned char));
	sigma_extend = (unsigned char*) calloc(code_dimension,
			sizeof(unsigned char));
	e = (unsigned char*) calloc(code_length, sizeof(unsigned char));
	e1 = (gf*) calloc(code_length, sizeof(gf));

	/*ETAPE_1:  Choose randomly  m ←  F_q^k, m is seen as a sequence of k_prime integer modulo 2^6.......
	 *****************************************************************************************************/

	m = random_m(k_prime, gf_card_sf);
	m_extend = extend(m, k_prime, code_dimension);

	/*ETAPE_2:  Compute r = G(m) and d = H(m) with  G(x) = sponge(x,k) and H(x) = sponge(x,k_prime).........
	 *******************************************************************************************************/
	dd = (gf*) calloc(k_prime, sizeof(gf));

	r = sponge(m_extend, code_dimension);
	for (i = 0; i < code_dimension; i++)
		r[i] = r[i] % gf_card_sf;

	d = sponge(m, k_prime);
	for (i = 0; i < k_prime; i++)
		d[i] = d[i] % gf_card_sf;

	for (i = 0; i < k_prime; i++) { //Conversion
		dd[i] = d[i];
	}

	//cfile_vec_F6("d_file.txt",k_prime, dd);

	/*ETAPE_3:  Parse r as (ρ||σ) then set u = (ρ||m).......................................................
	 *******************************************************************************************************/

	for (i = 0; i < code_dimension; i++) {
		if (i < k_sec)
			rho[i] = r[i] % gf_card_sf;  //rho recovery
		else
			sigma[i - k_sec] = r[i] % gf_card_sf; // sigma recovery
	}

	for (i = 0; i < code_dimension; i++) {
		if (i < k_sec)
			u[i] = rho[i];
		else
			u[i] = m[i - k_sec];
	}

	/*ETAPE_4: Generate error vector e of length n and weight w from sigma.........................................
	 *********************************************************************************************************/

	sigma_extend = extend(sigma, k_prime, code_dimension);
	hash_sigma = sponge(sigma_extend, code_length);

	e = random_e(code_length, gf_card_sf, n0_w, hash_sigma);

//           //cfile_vec_char("erreur.txt", code_length, e);

	/*ETAPE_5: Recovery of G and Compute c = uG + e................................................................
	 *************************************************************************************************************/
	gf* c2 = (gf*) calloc(code_length, sizeof(gf));

	c = (gf*) calloc(code_length, sizeof(gf));

	//binmat_t G = mat_ini(code_dimension,code_length-code_dimension);
	binmat_t G = mat_ini(code_dimension, code_length);
	//set_Public_matrix(pk, code_dimension, code_length-code_dimension, G);
	recup_pk(pk, G);

	for (i = 0; i < code_dimension; i++)
		u1[i] = (gf) u[i];

	//c1=produit_vector_matrix_Sf(u1,G);

	//for(i = 0;i<code_length;i++){
	//if(i<code_dimension) c[i] = u1[i] ;
	//else                 c[i] = c1[i-code_dimension];
	//}

	c = produit_vector_matrix_Sf(u1, G);

	for (i = 0; i < code_length; i++) { //Conversion
		e1[i] = e[i];
	}

	for (i = 0; i < code_length; i++)
		c2[i] = c[i] ^ e1[i];

	for (i = 0; i < code_length + k_prime; i++) {
		if (i < code_length)
			ct[i] = c2[i];
		else
			ct[i] = dd[i - code_length];

	}

	//cfile_vec_F6("chiffrer.txt",code_length,c2);

	/*ETAPE_6: Compute K = K(m)...................................................................................
	 ************************************************************************************************************/

	unsigned char* K = sponge(m_extend, ss_lenght);
	for (i = 0; i < ss_lenght; i++)
		ss[i] = K[i];

	return 0;
	/*END********************************************************************************************************/

}

