#include "../Masking/random.h"
#include "./test/cpucycles.h"
#include "masked_polyvec_operations.h"
#include "masking_interface.h"
#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
#pragma GCC diagnostic ignored "-Wint-conversion"

#if DILITHIUM_MODE == 2 && N_SHARES > 1 && N_SHARES < 10

/**
 * @brief Fast matrix multiplication for optimized masked signing
 * @param C Output share matrix
 * @param mat Precomputed matrix A
 * @param s1hats Input masked shares
 * @return Cycle count used
 */
uint64_t fast_matrix_mult(
    polyveck C[N_SHARES],
    const polyvecl mat[K],
    const polyvecl s1hats[N_SHARES])
{
    poly add, add2, sub;
    poly t1, t2, t3;

    // Intermediate buffers
    poly T1[N_SHARES][2];
    poly T2[N_SHARES][2];
    poly T3[3][K-1];
    poly T4[N_SHARES][2][K-1];

    uint64_t cycles = cpucycles();

    // Compute T1 = B[2k] * (A0[2k] + B[2k+1])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_add(&add, &mat[0].vec[2*k], &s1hats[j].vec[2*k+1]);
            poly_pointwise_montgomery(&T1[j][k], &s1hats[j].vec[2*k], &add);
        }
    }

    // Compute T2 = B[2k+1] * (A0[2k+1] - B[2k])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_sub(&sub, &mat[0].vec[2*k+1], &s1hats[j].vec[2*k]);
            poly_pointwise_montgomery(&T2[j][k], &s1hats[j].vec[2*k+1], &sub);
        }
    }

    // Compute T3 = Ai[2k+1] * (A0[2k] + Ai[2k]) for i=1..K-1
    for (int k = 0; k < 2; k++) {
        for (int i = 1; i < K; i++) {
            poly_add(&add, &mat[0].vec[2*k], &mat[i].vec[2*k]);
            poly_pointwise_montgomery(&T3[k][i-1], &mat[i].vec[2*k+1], &add);
        }
    }

    // Compute T4 = (B[2k] + Ai[2k+1]) * (B[2k+1] + A0[2k] + Ai[2k])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_add(&t2, &s1hats[j].vec[2*k+1], &mat[0].vec[2*k]);
            for (int i = 1; i < K; i++) {
                poly_add(&t1, &s1hats[j].vec[2*k], &mat[i].vec[2*k+1]);
                poly_add(&t3, &t2, &mat[i].vec[2*k]);
                poly_pointwise_montgomery(&T4[j][k][i-1], &t1, &t3);
            }
        }
    }

    // Sum T3 for final combination
    for (int i = 1; i < K; i++) {
        poly_add(&T3[2][i-1], &T3[1][i-1], &T3[0][i-1]);
    }

    // Assemble output rows
    for (int j = 0; j < N_SHARES; j++) {
        // Row 0
        poly_add(&C[j].vec[0], &T1[j][1], &T1[j][0]);

        // Rows 1..K-1
        for (int i = 1; i < K; i++) {
            poly_add(&t3, &T4[j][0][i-1], &T4[j][1][i-1]);
            poly_sub(&C[j].vec[i], &t3, &C[j].vec[0]);
            poly_sub(&C[j].vec[i], &C[j].vec[i], &T3[2][i-1]);
        }

        // Finalize row 0
        poly_add(&C[j].vec[0], &C[j].vec[0], &T2[j][0]);
        poly_add(&C[j].vec[0], &C[j].vec[0], &T2[j][1]);
    }

    return cpucycles() - cycles;
}

/**
 * @brief Further optimized fast matrix multiplication (reuses precomputed T3)
 * @param C Output share matrix
 * @param mat Precomputed matrix A
 * @param s1hats Input masked shares
 * @param T3 Precomputed T3 terms
 * @return Cycle count used
 */
uint64_t fast_matrix_mult_precom(
    polyveck C[N_SHARES],
    const polyvecl mat[K],
    const polyvecl s1hats[N_SHARES],
    poly T3[3][K-1])
{
    poly add, add2, sub;
    poly t1, t2, t3;

    // Intermediate buffers
    poly T1[N_SHARES][2];
    poly T2[N_SHARES][2];
    poly T4[N_SHARES][2][K-1];

    uint64_t cycles = cpucycles();

    // Compute T1 = B[2k] * (A0[2k] + B[2k+1])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_add(&add, &mat[0].vec[2*k], &s1hats[j].vec[2*k+1]);
            poly_pointwise_montgomery(&T1[j][k], &s1hats[j].vec[2*k], &add);
        }
    }

    // Compute T2 = B[2k+1] * (A0[2k+1] - B[2k])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_sub(&sub, &mat[0].vec[2*k+1], &s1hats[j].vec[2*k]);
            poly_pointwise_montgomery(&T2[j][k], &s1hats[j].vec[2*k+1], &sub);
        }
    }

    // Compute T4 = (B[2k] + Ai[2k+1]) * (B[2k+1] + A0[2k] + Ai[2k])
    for (int j = 0; j < N_SHARES; j++) {
        for (int k = 0; k < 2; k++) {
            poly_add(&t2, &s1hats[j].vec[2*k+1], &mat[0].vec[2*k]);
            for (int i = 1; i < K; i++) {
                poly_add(&t1, &s1hats[j].vec[2*k], &mat[i].vec[2*k+1]);
                poly_add(&t3, &t2, &mat[i].vec[2*k]);
                poly_pointwise_montgomery(&T4[j][k][i-1], &t1, &t3);
            }
        }
    }

    // Assemble output rows
    for (int j = 0; j < N_SHARES; j++) {
        // Row 0
        poly_add(&C[j].vec[0], &T1[j][1], &T1[j][0]);

        // Rows 1..K-1
        for (int i = 1; i < K; i++) {
            poly_add(&t3, &T4[j][0][i-1], &T4[j][1][i-1]);
            poly_sub(&C[j].vec[i], &t3, &C[j].vec[0]);
            poly_sub(&C[j].vec[i], &C[j].vec[i], &T3[2][i-1]);
        }

        // Finalize row 0
        poly_add(&C[j].vec[0], &C[j].vec[0], &T2[j][0]);
        poly_add(&C[j].vec[0], &C[j].vec[0], &T2[j][1]);
    }

    return cpucycles() - cycles;
}

#endif

/**
 * @brief Wrapper: masked matrix multiplication (with fast path)
 * @param t Output masked polynomial vector
 * @param mat Precomputed matrix
 * @param v Input masked polynomial vector
 * @return Cycle saving from optimization
 */
uint64_t masked_polyvec_matrix_pointwise_montgomery(
    masked_polyveck *t,
    const polyvecl mat[K],
    const masked_polyvecl *v)
{
    uint64_t cyc = 0;
    cyc = cpucycles();
    for (int i = 0; i < N_SHARES; ++i) {
        polyvec_matrix_pointwise_montgomery(&t->shares[i], mat, &v->shares[i]);
    }
    cyc = cpucycles() - cyc;
    return cyc;
}
uint64_t masked_polyvec_matrix_pointwise_montgomery_opt(
    masked_polyveck *t,
    const polyvecl mat[K],
    const masked_polyvecl *v)
{
    uint64_t cyc = 0;

#if DILITHIUM_MODE == 2 && N_SHARES > 1 && N_SHARES < 10
   cyc = fast_matrix_mult(t, mat, v);
   
#else
    // Default implementation
    for (int i = 0; i < N_SHARES; ++i) {
        polyvec_matrix_pointwise_montgomery(&t->shares[i], mat, &v->shares[i]);
    }
#endif

    return cyc;
}
/**
 * @brief Wrapper: fully optimized masked matrix multiplication
 * @param t Output masked polynomial vector
 * @param mat Precomputed matrix
 * @param v Input masked polynomial vector
 * @param T3 Precomputed T3 terms
 * @return Cycle count used
 */
uint64_t masked_polyvec_matrix_pointwise_montgomery_precom(
    masked_polyveck *t,
    const polyvecl mat[K],
    const masked_polyvecl *v,
    poly T3[3][K-1])
{
    uint64_t cyc = 0;

#if DILITHIUM_MODE == 2 && N_SHARES > 1 && N_SHARES < 10
    cyc = fast_matrix_mult_precom(t, mat, v, T3);
#else
    for (int i = 0; i < N_SHARES; ++i) {
        polyvec_matrix_pointwise_montgomery(&t->shares[i], mat, &v->shares[i]);
    }
#endif
    return cyc;
}

void	masked_polyveck_reduce(masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_reduce(&(v->shares[i]));
}
void	masked_polyveck_invntt_tomont(masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_invntt_tomont(&(v->shares[i]));
}

void	masked_polyvecl_ntt(masked_polyvecl *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyvecl_ntt(&(v->shares[i]));
}

void	masked_polyveck_ntt(masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_ntt(&(v->shares[i]));
}

void	masked_polyvecl_invntt_tomont(masked_polyvecl *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyvecl_invntt_tomont(&(v->shares[i]));
}

void	masked_polyvecl_pointwise_poly_montgomery(masked_polyvecl *r,
		const poly *a, const masked_polyvecl *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyvecl_pointwise_poly_montgomery(&(r->shares[i]), a,
			&(v->shares[i]));
}

void	masked_polyvecl_add(masked_polyvecl *w, const masked_polyvecl *u,
		const masked_polyvecl *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyvecl_add(&(w->shares[i]), &(u->shares[i]), &(v->shares[i]));
}

void	masked_polyvecl_reduce(masked_polyvecl *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyvecl_reduce(&(v->shares[i]));
}

void	masked_polyveck_pointwise_poly_montgomery(masked_polyveck *r,
		const poly *a, const masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_pointwise_poly_montgomery(&(r->shares[i]), a,
			&(v->shares[i]));
}
void	masked_polyveck_sub(masked_polyveck *w, const masked_polyveck *u,
		const masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_sub(&(w->shares[i]), &(u->shares[i]), &(v->shares[i]));
}

void	masked_polyveck_caddq(masked_polyveck *v)
{
	for (int i = 0; i < N_SHARES; ++i)
		polyveck_caddq(&(v->shares[i]));
}
