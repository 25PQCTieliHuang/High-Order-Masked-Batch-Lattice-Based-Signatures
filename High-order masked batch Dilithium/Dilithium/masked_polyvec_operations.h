#ifndef POLYVEC_OP_H
#define POLYVEC_OP_H


#include "polyvec.h"
#include "masking_interface.h"


//A*y
uint64_t fast_matrix_mult(polyveck C[N_SHARES], const polyvecl mat[K], const polyvecl s1hats[N_SHARES]);
uint64_t fast_matrix_mult_precom(polyveck C[N_SHARES],const polyvecl mat[K],const polyvecl s1hats[N_SHARES],poly T3[3][K-1]);
uint64_t masked_polyvec_matrix_pointwise_montgomery_opt(masked_polyveck *t, const polyvecl mat[K], const masked_polyvecl *v);
uint64_t masked_polyvec_matrix_pointwise_montgomery_precom(masked_polyveck *t, const polyvecl mat[K], const masked_polyvecl *v, poly T3[3][K-1]);
uint64_t masked_polyvec_matrix_pointwise_montgomery(masked_polyveck *t, const polyvecl mat[K], const masked_polyvecl *v);
void masked_polyveck_reduce(masked_polyveck *v);
void masked_polyveck_invntt_tomont(masked_polyveck *v);
void masked_polyvecl_ntt(masked_polyvecl *v);
void masked_polyveck_ntt(masked_polyveck *v);

//decompose
void masked_polyveck_caddq(masked_polyveck *v);

//sc+y
void masked_polyvecl_pointwise_poly_montgomery(masked_polyvecl *r, const poly *a, const masked_polyvecl *v);
void masked_polyvecl_invntt_tomont(masked_polyvecl *v);
void masked_polyvecl_add(masked_polyvecl *w, const masked_polyvecl *u, const masked_polyvecl *v);
void masked_polyvecl_reduce(masked_polyvecl *v);

//w-cs_2
void masked_polyveck_pointwise_poly_montgomery(masked_polyveck *r, const poly *a, const masked_polyveck *v);
void masked_polyveck_sub(masked_polyveck *w, const masked_polyveck *u, const masked_polyveck *v);





#endif