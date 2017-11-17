
#include<stdlib.h>
#include<stdio.h>
#include<sys/types.h>
#include<stdint.h>
#include"matrix.h"
#include"time.h"



unsigned char* extend(unsigned char* m, int size1 , int size2);
int  weight( gf* r,int size);

unsigned char* random_m(int size, int q);

int indice_in_vec(unsigned int * v, int j, int size);

unsigned char* random_e(int size, int q, int w,  unsigned char* sigma);


int  compare(unsigned char* tab1,unsigned char* tab2, int size);

unsigned char* gf_to_char(gf* a, int lenght);
void recup_pk(const unsigned char * pk,binmat_t G);
void store_pk(binmat_t M,unsigned char * pk);

void store_sk(gf * u, gf * v, gf * z,unsigned char * sk);
void set_uvz( gf * u, gf* v, gf * z,const unsigned char * sk );




