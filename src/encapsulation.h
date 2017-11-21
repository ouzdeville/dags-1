
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
#include <assert.h>

#include "round.h"
#include "sponge.h"
#include "fichier.h"
#include "util.h"



int encapsulation(const unsigned char *pk,unsigned char *ct,  unsigned char *ss);

#include <keccak/KangarooTwelve.h>

#define cus_len 4 // I random pick this number to fullfill parameter 