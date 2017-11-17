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
int disjoint_test(gf * u, gf * v ){
	int i,j;
	for (i=0;i<(order);i++){
		for (j=0;j<code_length;j++){
			if( u[i]==v[j]){
				return -1;
			}
		}
	}
	return 0;
}


int Test_disjoint(gf * L,int n){
	int i,j;
	int val=1;
	for (i=0;i<n;i++){
		for (j=i+1;j<n;j++){
			if(L[i]==L[j]){
				return -1;
			}
		}
	}
	return 0;
}


gf* Random_Vect(int m){
	int i,j,k;
	register int v;
	gf tmp;
	gf * vect,* U;
	U=(gf *)calloc(gf_card,sizeof(gf));
	vect=(gf *)calloc(m,sizeof(gf));
	U[0]=1;
	for(i=1;i<gf_card;i++){
		U[i]=i;
	}
	for (j=1;j<gf_card;j++){
		srand(time(NULL));
		v=(rand()%(j+1));
		tmp=U[j];
		U[j]=U[v+1];
		U[v+1]=tmp;
	}
	for (j=1;j<=m;j++){
		vect[j-1]=U[j];
	}

	return vect;
}





gf* Init_Random_U(){
	int i,j,k;
	register int v;
	gf tmp;
	gf * U;
	U=(gf *)calloc(gf_card,sizeof(gf));
	for(i=0;i<=gf_ord;i++){
		U[i]=i;
	}
	for (j=1;j<=gf_ord;j++){
		srand(time(NULL));
		v=(rand()%(j+1));
		tmp=U[j];
		U[j]=U[v+1];
		U[v+1]=tmp;
	}
	return U;
}

void Remove_From_U(gf elt,gf * U){
	int k;
	for (k=0;k<=gf_ord;k++){
		if (U[k]==elt){
			U[k]=0;
			break;
		}
	}
}

void binary_quasi_dyadic_sig(int m, int n, int t, int * b, gf * h_sig, gf * w ){
	int i,k,j,s,p,l,c,r,ct,consistent_root,consistent_support_block;
	int const C=((gf_card)/t);
	gf * U, * V,* h;
	gf h_i_inv, h_j_inv,sum_inv_h_i_j_0,sum_inv_h_i_0;
	h=(gf *)calloc(gf_card,sizeof(gf));
	U=(gf *)calloc(gf_card,sizeof(gf));
	V=(gf *)calloc(gf_card,sizeof(gf));
	do{
		U=Init_Random_U();
		h[0]=U[1];
		U[1]=0;
		k=1;
		for (s=0;s<m;s++){
			i=1<<s;
			h[i]=U[i+1];
			Remove_From_U(h[i],U);
			for (j=1;j<i;j++){
				h[i+j]=0;
				if ( (h[i]!=0 ) && (h[j]!=0 ) ){
					sum_inv_h_i_j_0=0;
					sum_inv_h_i_j_0=gf_Inv(h[i])^gf_Inv(h[j]);
					sum_inv_h_i_j_0=(sum_inv_h_i_j_0)^ (gf_Inv(h[0]));
					if(sum_inv_h_i_j_0!=0){
						h[i+j]=gf_Inv(sum_inv_h_i_j_0);
						Remove_From_U(h[i+j],U);
					}
					else {
					h[i+j]=0;
		    		}
				}
				else {
					h[i+j]=0;
		    	}
			}
		}

		c=0;
		V=Init_Random_U();
		consistent_root=1;
		for (p=0;p<t;p++){
					consistent_root=consistent_root & (h[p]!=0);
				}
		if(consistent_root){
			b[0]=0;
			c=1;
			for (r=0;r<t;r++){
				sum_inv_h_i_0= (gf_Inv(h[r]))^(gf_Inv(h[0]));
				Remove_From_U(gf_Inv(h[r]),V);
				Remove_From_U(sum_inv_h_i_0,V);
			}
			for (j=1;j<C;j++){
				consistent_support_block=1;
				for (p=j*t;p<(j+1)*t;p++){
					consistent_support_block=consistent_support_block & (h[p]!=0);
				}
				if(consistent_support_block){
					b[c]=j;
					c=c+1;
					for (l=j*t;l<(j+1)*t;l++){
						sum_inv_h_i_0;
						sum_inv_h_i_0= (gf_Inv(h[l]))^(gf_Inv(h[0]));
						Remove_From_U(sum_inv_h_i_0,V);
					}
				}
			}

		}

		ct=c*t;
	}while(ct<n);
	 // Computing w: We just one value of omega. So we stop at the first non-zero element of V.
	for (j=0;j<gf_card;j++){
		if(V[j]){
			*w=V[j];
			break;
		}
	}
	/******************************************
	 We choose n0=33 consistents blocks from all the consistents blocks given by the vector  b;
	We then obtain
	******************************************/
	l=0;
	int l1=0;
	for (j=0;j<n0_val;j++){
		for (k=0;k<order;k++){
			l=(order)*j+k;
			l1=(order)*b[j] + k;
			h_sig[l]=h[l1];
		}
	}

	free(U);
	free(V);
	free(h);
}



void Cauchy_Support(gf * Support, gf * W,gf * w){
	int i;
	gf sum_inv_h_i_0;
	gf * h;
	int * b,test_u=0,test_v=0,test_u_inter_v=0;
do{
	b=(int*)calloc(gf_card,sizeof(int));
	h=(gf*)calloc(code_length,sizeof(gf));
	binary_quasi_dyadic_sig(gf_extd,code_length,order,b,h,w);
	for (i=0;i<code_length;i++){
		sum_inv_h_i_0=0;
		sum_inv_h_i_0=(gf_Inv(h[i]))^(gf_Inv(h[0]));
		Support[i]=(sum_inv_h_i_0)^(w[0]);
	}
	for (i=0;i<order;i++){
		W[i]=(gf_Inv(h[i]))^(w[0]);
	}
        test_u=Test_disjoint(Support,code_length);
        test_v=Test_disjoint(W,order);
        test_u_inter_v=disjoint_test(W,Support);
        printf("\n calcul\n");
        

}while((test_u!=0) || (test_v!=0) || (test_u_inter_v!=0));

	free(h);
	free(b);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//***********************************************************************************
//      The fonction key_pair generates the public key and the secret key which will be stored in files
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

int key_pair(unsigned char *pk, unsigned char *sk) {
	gf * u, *v, *w, *z;
	binmat_t H;
	u = (gf*) calloc(order, sizeof(gf));
	v = (gf*) calloc(code_length, sizeof(gf));
	w = (gf*) calloc(code_length, sizeof(gf));
	z = (gf*) calloc(n0_val, sizeof(gf));

	gf_init(6);
	Cauchy_Support(v, u, w);
	cfile_vec_F12("omega.txt", order, u);
	z = Random_Vect(n0_val);
	H = mat_ini(pol_deg * (order), code_length);

	// construction de la matrice H
	secret_matrix(H, u, v, z);

	//cfile_matrix_F12("secret_matrix.txt", H.rown, H.coln, H);

	binmat_t H_base, H_syst;

//*********************************************************************************************************
//   The matrix H_base is obtained by the projection of the matrix H_fin into the base field through the
//                                 	 function  'mat_Into_base'
///////////////////////////////////////////////////////////////////////////////////////////////////////////
	H_base = mat_Into_base(H);
	H_syst = mat_copy(H_base);

	binmat_t S;
	S = mat_ini_Id(H_syst.rown);

//**************************************************************************************************************
//  Transform H_syst into its systematic by the function "syst" .
//  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int valueretur = syst(H_syst);
	if (valueretur != 0) {
		return valueretur;
	}

//*************************************************************************************************************
//       H_syst is in the form (G | I) we determine G and store it  in the file "pubkey.text" as  the publci key
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//G=mat_ini(code_dimension,code_length-code_dimension);
//G_mat(G,H_syst);
//set_public_key( G,pk);
	store_pk(H_syst, pk);
	store_sk(u, v, z, sk);

	return 0;

}

