/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                            *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.                             *
 **********************************************************************************************
 */

#include "encapsulation.h"

int encapsulation(const unsigned char *pk, unsigned char *ct, unsigned char *ss) {
	gf_init(6);

	unsigned char *m, *d, *rho, *sigma, *error_array, *hash_sigma, *r, *m_extend,
			*sigma_extend;
	gf *c, *u, *dd;
	const unsigned char *custom = "DAGs"; // customization = "DAGs";
	int i;
	int test; // Catch error 
	/*
	 * Memory's allocation
	 */
	m = (unsigned char*) calloc(k_prime, sizeof(unsigned char));
	rho = (unsigned char*) calloc(k_sec, sizeof(unsigned char));
	sigma = (unsigned char*) calloc(code_dimension - k_sec,
			sizeof(unsigned char));

	hash_sigma = (unsigned char*) calloc(code_length, sizeof(unsigned char));//Should remove
	u = (gf*) calloc(code_dimension, sizeof(gf));
	r = (unsigned char*) calloc(code_dimension, sizeof(unsigned char));//Should remove
	sigma_extend = (unsigned char*) calloc(code_dimension,//Should remove
			sizeof(unsigned char));
	dd = (gf*) calloc(k_prime, sizeof(gf));

	/*
	 * Step_1:  Choose randomly  m ←  F_q^k, m is seen as a sequence of k_prime integer
	 * modulo 2^6
	 */

	m = random_m(k_prime, gf_card_sf);
	m_extend = extend(m, k_prime, code_dimension); // TODO: REMOVE

	/*
	 * Step_2:  Compute r = G(m) and d = H(m) with  G(x) = sponge(x,k) and
	 * H(x) = sponge(x,k_prime)
	 */

	// Replace by KangarooTwelve
	// r = sponge(m_extend, code_dimension);
	// m: input type unsigned char len k_prime | r: output type unsigned char len code_dimesion
	test = KangarooTwelve(m, k_prime, r, code_dimension, custom, cus_len);
	assert(test == 0); // Catch Error

	for (i = 0; i < code_dimension; i++)
		// Optimize modulo
		//r[i] = r[i] % gf_card_sf;
		r[i] = r[i] & (gf_card_sf - 1);

	// Replace by KangarooTwelve
	// d = sponge(m, k_prime);
	// m: input type unsigned char len k_prime | d: output type unsigned char len k_prime
	test = KangarooTwelve(m, k_prime, d, k_prime, custom, cus_len);
	assert(test == 0); // Catch Error

	// Type conversion 
	for (i = 0; i < k_prime; i++)
		// Optimize modulo 
		//dd[i]= (unsigned char)(d[i] % gf_card_sf);
		dd[i] = (unsigned char)(d[i] & (gf_card_sf -1));

	free(d);


	//cfile_vec_F6("d_file.txt",k_prime, dd);

	/*
	 * Step_3:  Parse r as (ρ||σ) then set u = (ρ||m)
	 */
	for (i = 0; i < code_dimension; i++) {
		if (i < k_sec)
			// Optimize modulo
			// rho[i] = (unsigned char)(r[i] % gf_card_sf);  //rho recovery
			rho[i] = (unsigned char)(r[i] & (gf_card_sf-1));	  //rho recovery
		else
			// Optimize modulo
			// sigma[i - k_sec] = (unsigned char)(r[i] % gf_card_sf); // sigma recovery
			sigma[i - k_sec] = (unsigned char)(r[i] & (gf_card_sf-1)); // sigma recovery
	}

	for (i = 0; i < code_dimension; i++) {
		if (i < k_sec){
			u[i] = ((unsigned char)rho[i]);
		}
		else{
			u[i] = ((unsigned char)m[i - k_sec]);
		}
	}
	free(rho);


	/*
	 * Step_4: Generate error vector e of length n and weight w from sigma
	 * TODO verify that extend is supposed to be size code_dimension or code_length
	 * originally code_dimension
	 */
	sigma_extend = extend(sigma, k_prime, code_dimension); // TODO: REMOVE
	// Replace by KangarooTwelve
	// hash_sigma = sponge(sigma_extend, code_length);
	// sigma: input type unsigned char len k_prime | hash_sigma: output type unsigned char len code_length
	test = KangarooTwelve(sigma, k_prime, hash_sigma, code_length, custom, cus_len);
	assert(test == 0); // Catch Error

	error_array = random_e(code_length, gf_card_sf, n0_w, hash_sigma);

//           //cfile_vec_char("erreur.txt", code_length, e);
	free(sigma);
	free(sigma_extend);

	/*
	 * Step_5: Recovery of G and Compute c = uG + e
	 */
	binmat_t G = mat_ini(code_dimension, code_length);

	//set_Public_matrix(pk, code_dimension, code_length-code_dimension, G);
	recup_pk(pk, G);

	c = produit_vector_matrix_Sf(u, G);
	mat_free(G);
	free(u);

	for (i = 0; i < code_length + k_prime; i++) {
		if (i < code_length){
			ct[i] = (unsigned char)((unsigned char)c[i] ^ 
				(unsigned char)error_array[i]);
		}
		else{
			ct[i] = dd[i - code_length];
		}

	}
	free(c);
	free(dd);
	free(error_array);
	//cfile_vec_F6("chiffrer.txt",code_length,c2);

	/*
	 * Step_6: Compute K = K(m)
	 */
	unsigned char* K;
	// Replace by KangarooTwelve
	// K = sponge(m_extend, ss_lenght);
	// m: input type unsigned char len k_prime | K: output type unsigned char len ss_length
	test = KangarooTwelve(m, k_prime, K, ss_lenght, custom, cus_len);
	assert(test == 0); // Catch Error

	for (i = 0; i < ss_lenght; i++){
		ss[i] = K[i];
	}
	free(m_extend);
	free(K);
	free(m);
	return 0;
	/*END*/
}

