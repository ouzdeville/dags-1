/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                             *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ask permission.                             *
 **********************************************************************************************
 */

#include "decapsulation.h"

/*
 * decapsulation() fuction compute the shared secret (ss) of type
 * unsigned char* and the ciphertext *(ct) of type unsigned char* by using the
 * secret key (sk)                             *
 */
int
decapsulation (unsigned char *ss, const unsigned char *ct,
	       const unsigned char *sk)
{

  int i, test, decode_value;
  gf_init (6); //initialization of Log Antilog table
  const unsigned char *custom = (unsigned char*)"DAGs"; // customization = "DAGs";
  gf * u, *v, *z;
  gf *e, *mot, *c;
  unsigned char *m1, *rho1;
  unsigned char *r1, *d1, *rho2, *sigma, *e2, *hash_sigma1;
  unsigned char *d, *e_prime;
  unsigned char *K;
  binmat_t H, H_alt;

  /*
   * Recovery of the secret elements u, v and z from the secret key (sk)
   */


  u = (gf*) calloc (order, sizeof(gf));
  v = (gf*) calloc (code_length, sizeof(gf));
  z = (gf*) calloc (n0_val, sizeof(gf));
  set_uvz (u, v, z, sk);

  /*
   * Construction of secret matrix H from u, v and z
   */

  H = mat_ini (pol_deg * (order), code_length);
  secret_matrix (H, u, v, z);
  free(v);
  free(z);

  /*
   * Construction of alternant matrix H_alt from H and u
   */

  H_alt = mat_ini (pol_deg * (order), code_length);
  H_alt = alternant_matrix (H, u);
  mat_free(H);
  free(u);
  //aff_mat(H_alt);

  /*
   * Step_1 of the decapasulation :  Decode the noisy codeword C received as
   * part of the ciphertext ct = (c||d) with d is “ a plaintext confirmation”.
   * We obtain codeword mot = u1G and error e
   */
  e = (gf*) calloc (code_length, sizeof(gf));
  mot = (gf*) calloc (code_length, sizeof(gf));
  c = (gf*) calloc (code_length, sizeof(gf));

  //Could be a memcpy
  for (i = 0; i < code_length; i++){
      c[i] = ct[i];
  }

  decode_value = decoding_H (H_alt, c, e, mot);
  mat_free(H_alt);
  free(c);

  /*
   * Step_2 of the decapasulation :  Output ⊥ if decoding fails or wt(e) != n0_w
   */

  if (decode_value == -1 || weight (e, code_length) != n0_w){
      return -1;
  }

  e_prime = (unsigned char*) calloc (code_length, sizeof(unsigned char));

  //casting
  for (i = 0; i < code_length; i++){
    e_prime[i] = (unsigned char) e[i];
  }
  free(e);

  /*
   * Step_3 of the decapasulation :  Recover u_prime = mot and parse it as (rho1||m1)
   */
  m1 = (unsigned char*) calloc (k_prime, sizeof(unsigned char));
  rho1 = (unsigned char*) calloc (k_sec, sizeof(unsigned char));

  // Optimize modulo and removed copy to u1
  for (i = 0; i < code_dimension; i++){
      if (i < k_sec){
      	rho1[i] = ((unsigned char)mot[i]) & gf_ord_sf;  //rho1 recovery
      }
      else{
      	m1[i - k_sec] = ((unsigned char)mot[i]) & gf_ord_sf;     // m1 recovery
      }
  }
  free(mot);
  /*
   * Step_4 of the decapasulation :  Compute r1 = G(m1) and d1 = H(m1)
   */
  r1 = (unsigned char*) calloc (code_dimension,
  		sizeof(unsigned char));
  d1 = (unsigned char*) calloc (k_prime, sizeof(unsigned char));

  //Compute r1 = G(m1) where G is composed of sponge SHA-512 function and extend function.
  //m_extend is no longer required because we are using kangrooTwelve which handles sizing
  //r1 = sponge (m_extend, code_dimension);
  // Replace by KangarooTwelve
  // r = sponge(m_extend, code_dimension);
  // m: input type unsigned char len k_prime | r: output type unsigned char len code_dimesion
  test = KangarooTwelve (m1, k_prime, r1, code_dimension, custom,
  cus_len);
  assert (test == 0); // Catch Error

  //Already performed in copy to rho2 and sigma
//  for (i = 0; i < code_dimension; i++){
//  	// Optimize modulo
//  	//r1[i] = r1[i] % gf_card_sf;
//    r1[i] = r1[i] & gf_ord_sf;
//  }

  //Compute d1 = H(m1) where H is  sponge SHA-512 function
  //d1 = sponge (m1, k_prime);
  test = KangarooTwelve (m1, k_prime, d1, k_prime, custom, cus_len);
  assert (test == 0); // Catch Error
  //TODO remove. Fairly certain that this is not doing anything because
  //d1 is already a unsigned char*
  for (i = 0; i < k_prime; i++){
    d1[i] = d1[i] % gf_card_sf;
  }


  /*
   * Step_5 of the decapasulation: Parse r1 as (rho2||sigma1)
   */
  rho2 = (unsigned char*) calloc (k_sec, sizeof(unsigned char));
  sigma = (unsigned char*) calloc (code_dimension, sizeof(unsigned char));

  for (i = 0; i < code_dimension; i++){
      if (i < k_sec){
      	// Optimize modulo
      	//rho2[i] = r1[i] % gf_card_sf;  //rho2 recovery
      	rho2[i] = r1[i] & gf_ord_sf;  //rho2 recovery
      }
      else{
      	// Optimize modulo
      	//sigma[i - k_sec] = r1[i] % gf_card_sf; // sigma1 recovery
      	sigma[i - k_sec] = r1[i] & gf_ord_sf; // sigma1 recovery
      }
   }
  free(r1);
  /*
   * Step_6 of the decapasulation: Generate error vector e2 of length n and
   * weight n0_w from sigma1
   */
  e2 = (unsigned char*) calloc (code_length, sizeof(unsigned char));
  hash_sigma1 = (unsigned char*) calloc (code_length,	sizeof(unsigned char));

  //Hashing sigma_extend by using sponge SHA-512 function.
  // hash_sigma1 = sponge (sigma_extend, code_length);
  test = KangarooTwelve (sigma, k_prime, hash_sigma1, code_length,
  		custom, cus_len);
  assert (test == 0); // Catch Error
  free (sigma);

  //Generate error vector e2 of length code_length and weight n0_w from
  //hash_sigma1 by using random_e function.
  e2 = random_e (code_length, gf_card_sf, n0_w, hash_sigma1);


  /*
   * Step_7 of the decapasulation: Return ⊥ if e_prime distinct e2 or rho1
   * distinct rho2 or d distinct d1
   */
  d = (unsigned char*) calloc (k_prime, sizeof(unsigned char));


  //Recovery of “plaintext confirmation” d from ct.
  for (i = 0; i < k_prime; i++){
    d[i] = ct[i + code_length];
  }


  /*
   * Step_7 of the decapasulation: If the previous condition is not satisfied,
   * compute the shared secret ss by using sponge function and extend function
   */
  if (compare (e_prime, e2, code_length) == 0
      || compare (rho1, rho2, k_sec) == 0 || compare (d, d1, k_prime) == 0){
      return -1;
  }
  else{
      K = (unsigned char*) calloc (ss_lenght, sizeof(unsigned char));
      //unsigned char* K = sponge (m_extend, ss_lenght);
      test = KangarooTwelve (m1, k_prime, K, ss_lenght, custom, cus_len);
      assert (test == 0); // Catch Error
      for (i = 0; i < ss_lenght; i++){
      	ss[i] = K[i];
      }
  }
  free(rho1);
  free(rho2);
  free(d);
  free(d1);
  free(e_prime);
  free(e2);
  free(m1);
  free (K);


  return 0;

}
/*END**/

