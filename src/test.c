#include <stdlib.h>
#include <stdio.h>
#include "gf.h"
#include "util.h"



int main(){
    char *m = random_m(32, gf_card_sf);
    unsigned char * m_extend, *r;
    //printf('%d\n', sizeof(*m_extend));
    m_extend = extend(m, 32, 704);
    r = sponge(m_extend, 704);
    for (int i = 0; i < 704; i++)
    {
        r[i] = r[i] % gf_card_sf;
        printf("%x", r[i]);
    }
    return 0;
}
