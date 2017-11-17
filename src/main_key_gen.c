/*


********************************************************************************************
* QD_sRIVASTAVA ALGORITHM.                              *
* This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ask permission.                             *
**********************************************************************************************



#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gf.h"
#include "param.h"
#include "matrix.h"
#include"poly.h"
#include "fichier.h"
#include "decodage.h"
#include "key_gen.h"
#include "encapsulation.h"
#include "decapsulation.h"
#include "api.h"


int main(){
	float dureeKeyGen;
    clock_t tpsFinalKeyGen;
    clock_t tpsInitialKeyGen=clock()/CLOCKS_PER_SEC;

                 KEY_GENERATION

	unsigned char *pk, *sk, *ct, *ss,*s;
	pk=(unsigned char *)calloc(CRYPTO_PUBLICKEYBYTES,sizeof(unsigned char));
	sk=(unsigned char *)calloc(CRYPTO_SECRETKEYBYTES,sizeof(unsigned char));
        ct=(unsigned char *)calloc(CRYPTO_CIPHERTEXTBYTES,sizeof(unsigned char));
        ss=(unsigned char *)calloc(CRYPTO_BYTES,sizeof(unsigned char));
	//crypto_kem_keypair
	crypto_kem_keypair(pk,sk);    //75.018 s
	


                ENCAPSULATION

    crypto_kem_enc(ct,ss,pk);


               DECAPSULATION
   decapsulation();  //171.771 s s


	tpsFinalKeyGen=clock()/CLOCKS_PER_SEC;
	dureeKeyGen= tpsFinalKeyGen-tpsInitialKeyGen ;
	printf("Temps_d'execution: %f \n",dureeKeyGen);


}
*/
