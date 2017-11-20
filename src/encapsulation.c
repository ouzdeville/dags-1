/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                            *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ascode_dimmension permission.                             *
 **********************************************************************************************
 */
#include <assert.h>
#include "encapsulation.h"

void modulo_array(unsigned char *array, unsigned int array_size, unsigned int modulo)
{
	// TODO: Move this function to util.c, declare in util.h
	int i;
	// because gf_card_sf is power of two, so use the trick to perform fast modulo
	// https://stackoverflow.com/a/14997264
	for (i = 0; i < array_size; i++)
		array[i] = ( array[i] & (modulo - 1) );
}

int encapsulation(const unsigned char *pk, unsigned char *ct, unsigned char *ss)
{
	// Init subfield
	gf_init(6);

	/*
	 * Memory's allocation
	 */

	// TODO: Should move to header, easier to manage this fix variable
	const unsigned char *customization = "DAGs";
	/*
	 * Step_1:  Choose randomly  m ←  F_q^k, m is seen as a sequence of k_prime integer
	 * modulo 2^6
	 */

	// length of m is 32
	unsigned char *m = (unsigned char *)calloc(k_prime, sizeof(unsigned char));
	m = random_m(k_prime, gf_card_sf);

	/*
	 * Step_2:  Compute r = G(m) and d = H(m) with  G(x) = sponge(x,k) and
	 * H(x) = sponge(x,k_prime)
	 */

	// No need to use extend
	// m_extend = extend(m, k_prime, code_dimension);

	// KangarooTwelve( *input, size of input, *output, size of output, *customize, length of customize)

	// KangarooTwelve return 0 is success, 1 is fail. Need to catch it.
	unsigned char *r = (unsigned char *)calloc(code_dimension, sizeof(unsigned char)); // TODO: Should remove

	int test = 0;
	// m : input, k_prime : length of input, r: output, code_dimension: length of output, customization: domain-string, customization_length: length of customization
	test = KangarooTwelve(m, k_prime, r, code_dimension, customization, customization_length);
	assert(test == 0);

	modulo_array(r, code_dimension, gf_card_sf);

	// KangarooTwelve( *input, size of input, *output, size of output, *customize, length of customize)
	// d = sponge(m, k_prime)
	unsigned char *d = (unsigned char *)calloc(k_prime, sizeof(gf));

	// m: input, k_prime: length of input, d: output, k_prime: output length
	test = KangarooTwelve(m, k_prime, d, k_prime, customization, customization_length);
	assert(test == 0);

	modulo_array(d, k_prime, gf_card_sf);
	/*
	for (i = 0; i < k_prime; i++)
	{
		dd[i] = (unsigned char)(d[i] % gf_card_sf);
	}
	free(d);
	*/
	//cfile_vec_F6("d_file.txt",k_prime, dd);

	/*
	 * Step_3:  Parse r as (ρ||σ) then set u = (ρ||m)
	 */
	unsigned char *rho = (unsigned char *)calloc(k_sec, sizeof(unsigned char));
	unsigned char *sigma = (unsigned char *)calloc(code_dimension - k_sec, sizeof(unsigned char));

	int i;
	for (i = 0; i < code_dimension; i++)
	{
		if (i < k_sec)
			// optimized modulo
			// rho[i] = (unsigned char)(r[i] % gf_card_sf); //rho recovery
			rho[i] = (unsigned char)(r[i] & (gf_card_sf - 1)); //rho recovery
		else
			// optimized modulo
			// sigma[i - k_sec] = (unsigned char)(r[i] % gf_card_sf); // sigma recovery
			sigma[i - k_sec] = (unsigned char)(r[i] & (gf_card_sf - 1)); // sigma recovery
	}

	gf *u = (gf *)calloc(code_dimension, sizeof(gf));
	for (i = 0; i < code_dimension; i++)
	{
		if (i < k_sec)
			u[i] = ((unsigned char)rho[i]);
		else
			u[i] = ((unsigned char)m[i - k_sec]);
	}
	free(rho);

	/*
	 * Step_4: Generate error vector e of length n and weight w from sigma
	 * TODO verify that extend is supposed to be size code_dimension or code_length
	 * originally code_dimension
	 */

	// No need to extend
	// sigma_extend = extend(sigma, k_prime, code_dimension);

	// KangarooTwelve( *input, size of input, *output, size of output, *customize, length of customize)
	// hash_sigma = sponge(sigma_extend, code_length);
	unsigned char *hash_sigma = (unsigned char *)calloc(code_length, sizeof(unsigned char)); // TODO: Should remove

	// sigma : input, code_dimension: length input, hash_sigma : output, code_length: length output
	test = KangarooTwelve(sigma, code_dimension, hash_sigma, code_length, customization, customization_length);
	assert(test == 0);

	unsigned char *error_array = random_e(code_length, gf_card_sf, n0_w, hash_sigma);

	free(sigma);
	free(hash_sigma);

	/*
	 * Step_5: Recovery of G and Compute c = uG + e
	 */

	binmat_t G = mat_ini(code_dimension, code_length);

	//set_Public_matrix(pk, code_dimension, code_length-code_dimension, G);
	recup_pk(pk, G);

	gf *c = produit_vector_matrix_Sf(u, G);
	mat_free(G);
	free(u);

	for (i = 0; i < code_length + k_prime; i++)
	{
		if (i < code_length)
			ct[i] = (unsigned char)((unsigned char)c[i] ^ (unsigned char)error_array[i]);
		else
			ct[i] = d[i - code_length];
	}
	free(c);
	free(d);
	free(error_array);

	/*
	 * Step_6: Compute K = K(m)
	 */
	//unsigned char* K = sponge(m_extend, ss_lenght);
	

	// KangarooTwelve( *input, size of input, *output, size of output, *customize, length of customize)
	unsigned char *K = (unsigned char *)calloc(ss_lenght, sizeof(unsigned char));
	KangarooTwelve(m, k_prime, K, ss_lenght, customization, customization_length);

	// Now we have K

	free(K);
	free(m);
	return 0;
	/*END*/
}
