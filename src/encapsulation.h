
/*********************************************************************************************
* DAGS: Key Encapsulation using Dyadic GS Codes.                             *
* This code is exclusively intended for submission to the NIST Post=Quantum Cryptography.    *
* For any other usage , contact the author(s) to ascode_dimmension permission.                             *
**********************************************************************************************
*/


#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<stdint.h>


#include "round.h"

// replace sponge.h by KangarooTwelve --fastest branch of SHA3 
// #include "sponge.h"


#include "fichier.h"
#include "util.h"

#include <keccak/KangarooTwelve.h>

int encapsulation(const unsigned char *pk,unsigned char *ct,  unsigned char *ss);

#define customization_length 4 