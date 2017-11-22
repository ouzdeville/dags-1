/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                             *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ask permission.                             *
 **********************************************************************************************
 */

#include "decapsulation.h"

/**************************************************************************************************
 * decapsulation() fuction compute the shared secret (ss) of type  unsigned char* and the ciphertext
 (ct) of type unsigned char* by using the secret key (sk)                             *
 ***************************************************************************************************
 */
int
decapsulation (unsigned char *ss, const unsigned char *ct,
	       const unsigned char *sk)
{

  int i;
  gf_init (6); //initialization of Log Antilog table
  const unsigned char *custom = "DAGs"; // customization = "DAGs";
  int test; // Catch error

  /* Recovery of the secret elements u, v and z from the secret key (sk) ...........................
   *************************************************************************************************/

  gf * u, *v, *z;
  u = (gf*) calloc (order, sizeof(gf));
  v = (gf*) calloc (code_length, sizeof(gf));
  z = (gf*) calloc (n0_val, sizeof(gf));
  set_uvz (u, v, z, sk);

  /* Construction of secret matrix H from u, v and z  .........................................
   *************************************************************************************************/

  binmat_t H = mat_ini (pol_deg * (order), code_length);
  secret_matrix (H, u, v, z);

  /* Construction of alternant matrix H_alt from H and u  .........................................
   *************************************************************************************************/

  binmat_t H_alt = mat_ini (pol_deg * (order), code_length);
  H_alt = alternant_matrix (H, u);
  //aff_mat(H_alt);

  //gf * motcode=( gf*)calloc(code_length,sizeof( gf));
  //unsigned char* e1=( unsigned char*)calloc(code_length,sizeof( unsigned char));

  /*ETAPE_1 of the decapasulation :  Decode the noisy codeword C received as part of the ciphertext....
   ct = (c||d) with d is “ a plaintext confirmation”. We obtain codeword mot = u1G and error e.........
   *****************************************************************************************************/
  //variable declaration
  gf* e = (gf*) calloc (code_length, sizeof(gf));
  gf* mot = (gf*) calloc (code_length, sizeof(gf));
  gf* c = (gf*) calloc (code_length, sizeof(gf));

  for (i = 0; i < code_length; i++)
    {
      c[i] = ct[i];
    }
  int decode_value = decoding_H (H_alt, c, e, mot);

  /*ETAPE_2 of the decapasulation :  Output ⊥ if decoding fails or wt(e) != n0_w...........................
   *******************************************************************************************************/

  if (decode_value == -1 || weight (e, code_length) != n0_w)
    {
      return -1;
    }

  /*ETAPE_3 of the decapasulation :  Recover u1 such as mot = u1G and parse it as (rho1||m1)..............
   *******************************************************************************************************/
  //variable declaration
  unsigned char* u1 = (unsigned char*) calloc (code_dimension,
					       sizeof(unsigned char));
  unsigned char* m1 = (unsigned char*) calloc (k_prime, sizeof(unsigned char));
  unsigned char* rho1 = (unsigned char*) calloc (k_sec, sizeof(unsigned char));

  //Recover u1
  for (i = 0; i < code_dimension; i++)
    {
      u1[i] = (unsigned char) mot[i];
    }

  // parse u1 as (rho1||m1)
  unsigned char* m_extend = (unsigned char*) calloc (code_dimension,
						     sizeof(unsigned char));

  for (i = 0; i < code_dimension; i++)
    {
      if (i < k_sec)
	rho1[i] = u1[i] % gf_card_sf;  //rho1 recovery
      else
	m1[i - k_sec] = u1[i] % gf_card_sf;     // m1 recovery
    }

  /*ETAPE_4 of the decapasulation :  Compute r1 = G(m1) and d1 = H(m1)....................................
   *******************************************************************************************************/
  //variable declaration
  unsigned char* r1 = (unsigned char*) calloc (code_dimension,
					       sizeof(unsigned char));
  unsigned char* d1 = (unsigned char*) calloc (k_prime, sizeof(unsigned char));

  // extend m1 to m_extend of lenght code_dimension by using extend function.
  m_extend = (unsigned char*) extend (m1, k_prime, code_dimension);

  //Compute r1 = G(m1) where G is composed of sponge SHA-512 function and extend function.
  //r1 = sponge (m_extend, code_dimension);
  // Replace by KangarooTwelve
  // r = sponge(m_extend, code_dimension);
  // m: input type unsigned char len k_prime | r: output type unsigned char len code_dimesion
  test = KangarooTwelve (m_extend, k_prime, r1, code_dimension, custom,
  cus_len);
  assert (test == 0); // Catch Error
  for (i = 0; i < code_dimension; i++)
    r1[i] = r1[i] % gf_card_sf;

  //Compute d1 = H(m1) where H is  sponge SHA-512 function
  //d1 = sponge (m1, k_prime);
  test = KangarooTwelve (m1, k_prime, d1, k_prime, custom, cus_len);
  assert (test == 0); // Catch Error
  for (i = 0; i < k_prime; i++)
    d1[i] = d1[i] % gf_card_sf;

  /*ETAPE_5 of the decapasulation: Parse r1 as (rho2||sigma1).....................................................
   *********************************************************************************************************/
  //variable declaration
  unsigned char*rho2 = (unsigned char*) calloc (k_sec, sizeof(unsigned char));
  unsigned char* sigma1 = (unsigned char*) calloc (k_prime,
						   sizeof(unsigned char));

  for (i = 0; i < code_dimension; i++)
    {
      if (i < k_sec)
	rho2[i] = r1[i] % gf_card_sf;  //rho2 recovery
      else
	sigma1[i - k_sec] = r1[i] % gf_card_sf; // sigma1 recovery
    }

  /*ETAPE_6 of the decapasulation: Generate error vector e2 of length n and weight n0_w from sigma1............
   *************************************************************************************************************/
  //variable declaration
  unsigned char* e2 = (unsigned char*) calloc (code_length,
					       sizeof(unsigned char));
  unsigned char* sigma_extend = (unsigned char*) calloc (code_dimension,
							 sizeof(unsigned char));
  unsigned char* hash_sigma1 = (unsigned char*) calloc (code_length,
							sizeof(unsigned char));

  // extend sigma to sigma_extend of lenght code_dimension by using extend function.
  sigma_extend = extend (sigma1, k_prime, code_dimension);

  //Hashing sigma_extend by using sponge SHA-512 function.
  // hash_sigma1 = sponge (sigma_extend, code_length);

  test = KangarooTwelve (sigma_extend, k_prime, hash_sigma1, code_length,
			 custom,
			 cus_len);

  assert (test == 0); // Catch Error

  //Generate error vector e2 of length code_length and weight n0_w from hash_sigma1 by using random_e function.
  e2 = random_e (code_length, gf_card_sf, n0_w, hash_sigma1);

  /*ETAPE_7 of the decapasulation: Return ⊥ if e_prime distinct e2 or rho1 distinct rho2 or d distinct d1 .....
   ************************************************************************************************************/
  //variable declaration
  unsigned char* d = (unsigned char*) calloc (k_prime, sizeof(unsigned char));
  unsigned char* e_prime = (unsigned char*) calloc (code_length,
						    sizeof(unsigned char));

  //casting
  for (i = 0; i < code_length; i++)
    e_prime[i] = (unsigned char) e[i];

  //Recovery of “plaintext confirmation” d from ct.
  for (i = 0; i < k_prime; i++)
    d[i] = ct[i + code_length];


  /*ETAPE_7 of the decapasulation: If the previous condition is not satisfied, compute the shared secret ss by
   using sponge function and extend function
   ************************************************************************************************************/
  if (compare (e_prime, e2, code_length) == 0
      || compare (rho1, rho2, k_sec) == 0 || compare (d, d1, k_prime) == 0)
    {
      return -1;
    }
  else
    {
      unsigned char* K = (unsigned char*) calloc (ss_lenght,
						  sizeof(unsigned char));
      //unsigned char* K = sponge (m_extend, ss_lenght);
      test = KangarooTwelve (m_extend, k_prime, K, ss_lenght, custom, cus_len);
      assert (test == 0); // Catch Error
      for (i = 0; i < ss_lenght; i++)
	ss[i] = K[i];

      free (K);
    }
  free (sigma1);
  free (rho2);
  free (sigma_extend);

  return 0;

}
/*END********************************************************************************************************/

