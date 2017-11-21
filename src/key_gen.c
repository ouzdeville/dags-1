/*********************************************************************************************
 * DAGS: Key Encapsulation using Dyadic GS Codes.                       *
 * This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
 * For any other usage , contact the author(s) to ask permission.                             *
 **********************************************************************************************
 */

#include "key_gen.h"

////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////       Key Generation Elements   ////////////////////////
//////////////////QUASI-DYADIC-SIGNATURE 	AND  CAUCHY SUPPORT  GENERATION.////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

//***********************************************************************************
//                               Test_disjoint
////////////////////////////////////////////////////////////////////////////////////
int
disjoint_test (gf * u, gf * v)
{
  int i, j;
  for (i = 0; i < (order); i++)
    {
      for (j = 0; j < code_length; j++)
	{
	  if (u[i] == v[j])
	    {
	      return -1;
	    }
	}
    }
  return 0;
}

int Test_disjoint(gf * L,int n){
	int i,j;
	for (i=0;i<n;i++){
		for (j=i+1;j<n;j++){
			if(L[i]==L[j]){
				return -1;
			}
		}
	}
	return 0;
}


void
generate_random_vector (int m, gf *vect)
{
  int i, j;
  register int v;
  gf tmp;
  gf *U;
  U = (gf *) calloc (gf_card, sizeof(gf));
  unsigned char *random_bytes = malloc (gf_card);
  randombytes (random_bytes, gf_ord);
  U[0] = 1;
  for (i = 1; i < gf_card; i++)
    {
      U[i] = i;
    }
  for (j = 1; j < gf_card; j++)
    {
//		srand(time(NULL));
//		v=(rand()%(j+1));
      v = random_bytes[j] % (j + 1);
      tmp = U[j];
      U[j] = U[v + 1];
      U[v + 1] = tmp;
    }

  memcpy (vect, U + 1, (m) * sizeof(gf));

  free (U);
  free (random_bytes);
}

void
init_random_element (gf *U)
{
  int i, j;
  register int v;
  gf tmp;
  unsigned char *random_bytes = malloc (gf_ord);
  randombytes (random_bytes, gf_ord);
  for (i = 0; i <= gf_ord; i++)
    {
      U[i] = i;
    }
  //TODO verify that this is not supposed to be j < gf_ord
  //Also do we need two versions openssl rand() and randombytes()?
  for (j = 1; j <= gf_ord; j++)
    {
//		srand(time(NULL));
//		v=(rand()%(j+1));
      v = random_bytes[j] % (j + 1);
      tmp = U[j];
      U[j] = U[v + 1];
      U[v + 1] = tmp;
    }
}

void
Remove_From_U (gf elt, gf * U)
{
  int k;
  for (k = 0; k <= gf_ord; k++)
    {
      if (U[k] == elt)
	{
	  U[k] = 0;
	  break;
	}
    }
}

void
binary_quasi_dyadic_sig (int m, int n, int t, int * b, gf * h_sig, gf * w)
{
  int i, j, k, s, p, l, l1, c, r, consistent_root, consistent_support_block;
  int const C = ((gf_card) / t);
  gf * U, *V, *h;
  gf sum_inv_h_i_j_0, sum_inv_h_i_0;
  h = (gf *) calloc (gf_card, sizeof(gf));
  U = (gf *) calloc (gf_card, sizeof(gf));
  V = (gf *) calloc (gf_card, sizeof(gf));

  do
    {
      init_random_element (U);
      h[0] = U[1];
      U[1] = 0;

      for (s = 0; s < m; s++)
	{
	  i = 1 << s;
	  h[i] = U[i + 1];
	  Remove_From_U (h[i], U);
	  for (j = 1; j < i; j++)
	    {
	      h[i + j] = 0;
	      if ((h[i] != 0) && (h[j] != 0))
		{
		  sum_inv_h_i_j_0 = (gf_Inv(h[i])^gf_Inv(h[j])) ^(gf_Inv(h[0]));
		  if(sum_inv_h_i_j_0!=0)
		    {
		      h[i+j]=gf_Inv(sum_inv_h_i_j_0);
		      Remove_From_U(h[i+j],U);
		    }
		  else
		    {
		      h[i+j]=0;
		    }
		}
	      else
		{
		  h[i+j]=0;
		}
	    }
	}

      c = 0;
      init_random_element (V);
      consistent_root = 1;
      for (p = 0; p < t; p++)
	{
	  consistent_root = consistent_root & (h[p] != 0);
	}
      if (consistent_root)
	{
	  b[0] = 0;
	  c = 1;
	  for (r = 0; r < t; r++)
	    {
	      sum_inv_h_i_0 = (gf_Inv(h[r])) ^ (gf_Inv(h[0]));
	      Remove_From_U (gf_Inv(h[r]), V);
	      Remove_From_U (sum_inv_h_i_0, V);
	    }
	  for (j = 1; j < C; j++)
	    {
	      consistent_support_block = 1;
	      for (p = j * t; p < (j + 1) * t; p++)
		{
		  consistent_support_block = consistent_support_block
		      & (h[p] != 0);
		}
	      if (consistent_support_block)
		{
		  b[c] = j;
		  c = c + 1;
		  for (l = j * t; l < (j + 1) * t; l++)
		    {
		      sum_inv_h_i_0 = (gf_Inv(h[l])) ^ (gf_Inv(h[0]));
		      Remove_From_U (sum_inv_h_i_0, V);
		    }
		}
	    }

	}
    }
  while (c * t < n);

  // Computing w: We just one value of omega. So we stop at the first non-zero element of V.
  for (j = 0; j < gf_card; j++)
    {
      if (V[j])
	{
	  *w = V[j];
	  break;
	}
    }
  /******************************************
   We choose n0=33 consistent blocks from all the consistent blocks given by the vector  b;
   We then obtain
   ******************************************/
  for (j = 0; j < n0_val; j++)
    {
      for (k = 0; k < order; k++)
	{
	  l = (order) * j + k;
	  l1 = (order) * b[j] + k;
	  h_sig[l] = h[l1];
	}
    }

  free (U);
  free (V);
  free (h);
}

void
cauchy_support (gf * Support, gf * W, gf * w)
{
  int i;
  gf sum_inv_h_i_0;
  gf * h;
  int * b, test_u = 0, test_v = 0, test_u_inter_v = 0;
  do
    {
      b = (int*) calloc (gf_card, sizeof(int));
      h = (gf*) calloc (code_length, sizeof(gf));
      binary_quasi_dyadic_sig (gf_extd, code_length, order, b, h, w);
      for (i = 0; i < code_length; i++)
	{
	  sum_inv_h_i_0 = (gf_Inv(h[i])) ^ (gf_Inv(h[0]));
	  Support[i] = (sum_inv_h_i_0) ^ (w[0]);
	}
      for (i = 0; i < order; i++)
	{
	  W[i] = (gf_Inv(h[i])) ^ (w[0]);
	}
      test_u = Test_disjoint (Support, code_length);
      test_v = Test_disjoint (W, order);
      test_u_inter_v = disjoint_test (W, Support);
      //printf ("\n calcul\n");

    }
  while ((test_u != 0) || (test_v != 0) || (test_u_inter_v != 0));

  free (h);
  free (b);
}

/*
 * The function key_pair generates the public key and the
 * secret key which will be stored in files
 */
int
key_pair (unsigned char *pk, unsigned char *sk)
{
  gf * u, *v, *w, *z;
  int return_value;
  binmat_t H, H_syst;
  u = (gf*) calloc (order, sizeof(gf));
  v = (gf*) calloc (code_length, sizeof(gf));
  w = (gf*) calloc (code_length, sizeof(gf));
  z = (gf*) calloc (n0_val, sizeof(gf));

  gf_init (6);
  cauchy_support (v, u, w);
  free (w);
  cfile_vec_F12 ("omega.txt", order, u);
  generate_random_vector (n0_val, z);
  H = mat_ini (pol_deg * (order), code_length);
  // construction matrix H
  secret_matrix (H, u, v, z);

  //cfile_matrix_F12("secret_matrix.txt", H.rown, H.coln, H);

  /*
   * The matrix H_base is obtained by the projection of the matrix
   *  H into the base field through the function  'mat_Into_base'
   */
  H_syst = mat_Into_base (H);

// Transform H_syst into its systematic by the function "syst" .
  return_value = syst (H_syst);
  if (return_value != 0)
    {
      return return_value;
    }
  /*
   * H_syst is in the form (G | I) we determine G and store
   * it  in the file "pubkey.text" as  the public key
   */

//G=mat_ini(code_dimension,code_length-code_dimension);
//G_mat(G,H_syst);
//set_public_key( G,pk);
  store_pk (H_syst, pk);
  store_sk (u, v, z, sk);
  free (u);
  free (v);
  free (z);
  mat_free (H);
  mat_free (H_syst);
  return 0;

}

