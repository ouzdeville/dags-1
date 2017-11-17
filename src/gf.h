#include <stdint.h>

#ifndef GF_H
#define GF_H


typedef uint16_t gf;
typedef uint16_t gf_t;



#define gf_extd 12 // extension degree
#define gf_card (1<<gf_extd) // Field size
#define gf_ord ((gf_card) - 1) // Group order
#define m_val 2
#define gf_extd_sf gf_extd/m_val // subfield extension degree
#define gf_card_sf (1<<gf_extd_sf) // Field sizeml
#define gf_ord_sf ((gf_card_sf) - 1) // Group order


//gf_antilog_alpha[i] =gf_Mult_subfield(gf_antilog_alpha[i-1],alpha);

#define u_val  64
#define poly_primitif_field 4196


//#define alpha_val  524
//#define poly_primitif 4359
#define primitif_elt_field    gf_antilog_sf[34]
#define poly_primitif_subfield 67





//int gf_extension_degree, gf_cardinality, gf_multiplicative_order;
gf_t * gf_log_sf;
gf_t * gf_antilog_sf;

gf * gf_log;
gf * gf_antilog;









#define gf_unit() 1
#define gf_zero() 0
#define gf_add(x, y) ((x) ^ (y)) // Addition in the field

/*#define gf_antilog_sf(i) gf_antilog_sf[i] //alpha^i
#define gf_log_sf(x) gf_log_sf[x]


#define gf_antilog(i) gf_antilog[i] //alpha^i
#define gf_log(x) gf_log[x]*/



/* we obtain a value between 0 and (q-1) included, the class of 0 is
represented by 0 or q-1 (this is why we write _K->exp[q-1]=_K->exp[0]=1)*/


//#define _gf_modq_1_sf(d) (((d) & gf_ord_sf + ((d)>> gf_extd_sf)))
#define  _gf_modq_1_sf(d) ((d) % 63)
//#define gf_mul_fast_sf(x, y) ((y) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])] : 0)
#define gf_mul_fast_sf(x, y) ((y) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])] : 0)
#define gf_Mult_subfield(x, y) ((x) ? gf_mul_fast_sf(x, y) : 0) // Multiplication in the field : apha^i*alpha^j=alpha^(i+j)
#define gf_Pow_subfield(x,i) ( gf_antilog_sf[ ( _gf_modq_1_sf( gf_log_sf[x]*i))])
#define gf_sq_subfield(x) ((x) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] << 1)] : 0) // Square of an element of the field
#define gf_sqrt_subfield(x) ((x) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] << (gf_extd_sf-1))] : 0) // Squareroot of an element of the field

#define gf_div_subfield(x, y) ((x) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] - gf_log_sf[y])] : 0) // Division in the field : apha^i/alpha^j=alpha^(i-j)
#define gf_Inv_subfield(x) gf_antilog_sf[gf_ord_sf - gf_log_sf[x]] // Inverse in the field




//#define _gf_modq_1_sf(d) (((d) & gf_ord_sf + ((d)>> gf_extd_sf)))
#define  _gf_modq_1(d) ((d) % 4095)
//#define gf_mul_fast_sf(x, y) ((y) ? gf_antilog_sf[_gf_modq_1_sf(gf_log_sf[x] + gf_log_sf[y])] : 0)
#define gf_mul_fast(x, y) ((y) ? gf_antilog[_gf_modq_1(gf_log[x] + gf_log[y])] : 0)
#define gf_Mult(x, y) ((x) ? gf_mul_fast(x, y) : 0) // Multiplication in the field : apha^i*alpha^j=alpha^(i+j)
#define gf_Pow(x,i) ( gf_antilog[ ( _gf_modq_1( gf_log[x]*i))])
#define gf_sq(x) ((x) ? gf_antilog[_gf_modq_1(gf_log[x] << 1)] : 0) // Square of an element of the field
#define gf_sqrt(x) ((x) ? gf_antilog[_gf_modq_1(gf_log[x] << (gf_extd-1))] : 0) // Squareroot of an element of the field

#define gf_Div(x, y) ((x) ? gf_antilog[_gf_modq_1(gf_log[x] - gf_log[y])] : 0) // Division in the field : apha^i/alpha^j=alpha^(i-j)
#define gf_Inv(x) gf_antilog[gf_ord - gf_log[x]] // Inverse in the field







gf gf_Mult1(gf in0, gf in1);
gf gf_diff1(gf a, gf b);
gf gf_Inv1(gf in);
gf gf_sq1(gf in);
//gf gf_add(gf a, gf b);
gf gf_Div1(gf a, gf b );
gf gf_Pow1(gf  f, int n);







int gf_init(int extdeg);
gf_t gf_rand(int (*u8rnd)());
#endif
