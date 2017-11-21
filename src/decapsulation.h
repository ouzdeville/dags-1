

#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<stdint.h>

#include "decoding.h"
#include "round.h"
#include "sponge.h"
#include "fichier.h"
#include "util.h"


/**
 * @brief Method to perform the decapsulation using secret key
 * @param ss
 * @param ct
 * @param sk secret key to extract the key
 */
int decapsulation(unsigned char *ss,  const unsigned char *ct, const unsigned char *sk);





