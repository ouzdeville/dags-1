#ifndef SPONGE_H
#define SPONGE_H

int modulo (int a, int b);
uint8_t *sponge(uint8_t*,int32_t);
uint8_t *padding(uint8_t*,int32_t*);
void swap(char *, int);
#endif
