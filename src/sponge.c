#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<math.h>
#include<stdint.h>

#include "round.h"
#include "sponge.h"
#include "keccak_f.h"





/* FUNCTION HACHAGE SHA3**KECCACK**ROUND SHA3******************************************************
  *************************************************************************************************/

const uint64_t r[5][5]={ {0,36,3,41,18},
                         {1,44,10,45,2},
                         {62,6,43,15,61},
                         {28,55,25,21,56},
                         {27,20,39,8,14}
                       };

int modulo (int a, int b){

     if(b < 0) //you can check for b == 0 separately and do what you want
     return modulo(-a, -b);
     int ret = a % b;
     if(ret < 0)
     ret+=b;
     return ret;
}


uint64_t **sha3_round(uint64_t **A, uint64_t RC){

  uint64_t C[5];
  uint64_t D[5];
  uint64_t B[5][5];

  /* Theta step */
  uint8_t x,y;
  for(x=0;x<5;x++){
    C[x]=A[x][0] ^ A[x][1] ^ A[x][2]^ A[x][3] ^ A[x][4];
   }
  for(x=0;x<5;x++){
    D[x]=C[(x + 4) % 5] ^ ((C[(x + 1) % 5] << 1) | (C[(x + 1) % 5] >> 63));
   }
  for(x=0;x<5;x++){
    for(y=0;y<5;y++){
      A[x][y]=A[x][y]^D[x];
    }
  }

  /*Rho and pi steps*/
  for(x=0;x<5;x++){
    for(y=0;y<5;y++){
      B[y][modulo((2*x+3*y),5)]=((A[x][y] << r[x][y]) | (A[x][y] >> (64-r[x][y])));
    }
  }

  /*Xi state*/
  for(x=0;x<5;x++){
    for(y=0;y<5;y++){
      A[x][y]=B[x][y]^((~B[modulo((x+1),5)][y]) & B[modulo((x+2),5)][y]);
     }
  }

  /*Last step*/
  A[0][0]=A[0][0]^RC;

  return A;
}


// Keccak ==SHA3 -
const uint64_t RC[24]={ 0x0000000000000001,
                       0x0000000000008082,
                       0x800000000000808A,
                       0x8000000080008000,
                       0x000000000000808B,
                       0x0000000080000001,
                       0x8000000080008081,
                       0x8000000000008009,
                       0x000000000000008A,
                       0x0000000000000088,
                       0x0000000080008009,
                       0x000000008000000A,
                       0x000000008000808B,
                       0x800000000000008B,
                       0x8000000000008089,
                       0x8000000000008003,
                       0x8000000000008002,
                       0x8000000000000080,
                       0x000000000000800A,
                       0x800000008000000A,
                       0x8000000080008081,
                       0x8000000000008080,
                       0x0000000080000001,
                       0x8000000080008008
};

uint64_t **keccak_f(uint64_t **A){
  int32_t i;
  for(i=0;i<24;i++){
    A=sha3_round(A,RC[i]);
  }
  return A;
}

///Sponge-----------------------------------------------------------------------------

uint8_t *sponge(uint8_t* M,int32_t size){
  int32_t r=72;
  int32_t w=8;
  /*Padding*/
  if((size%r)!=0){//r=72 bytes
    M=padding(M,&size);
  }
  uint64_t *nM;
  nM=(uint64_t *)M;
  /*Initialization*/
  uint64_t **S=(uint64_t **)calloc(5,sizeof(uint64_t*));
  uint64_t i;
  for(i = 0; i < 5; i++) S[i] = (uint64_t *)calloc(5,sizeof(uint64_t));
  /*Absorbing Phase*/
  int32_t y,x;
  for(i=0;i<size/r;i++){//Each block has 72 bytes
    for(y=0;y<5;y++){
      for(x=0;x<5;x++){
        if((x+5*y)<(r/w)){
          S[x][y]=S[x][y] ^ *(nM+i*9+x+5*y);
          //      S=keccak_f(S);
        }
      }
    }
    S=keccak_f(S);

  }
  /*Squeezing phase*/
  int32_t b=0;
  uint64_t *Z=(uint64_t *)calloc(size,sizeof(uint64_t));
  while(b<size){
  	int32_t y;
  for(y=0;y<5;y++){
    for(x=0;x<5;x++){
      if((x+5*y)<(r/w)){
        *(Z+b)^=S[x][y];
        b++;
      }
    }
  }
 }
  return (uint8_t *)Z;
}

uint8_t *padding(uint8_t* M, int32_t* S){
  int32_t i=*S;
  int32_t newS=(*S+72-(*S%72));;
  uint8_t *nM;
  nM=malloc(*S+(72-(*S%72)));
   /*Copy string*/
  int32_t j;
  for(j=0;j<*S;j++){
    *(nM+j)=*(M+j);
  }
  *(nM+i)=0x01;
  i++;
  while(i<(newS-1)){
    *(nM+i)=0x00;
    i++;
  }
  *(nM+i)=0x80;
  i++;
  *S=i;
  return nM;
}
//*END**********************************************************************************************/

