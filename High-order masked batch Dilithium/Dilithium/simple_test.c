#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/resource.h>

#include "bench_masked_sign.h"
#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "symmetric.h"
#include "fips202.h"
#include "randombytes.h" 
#include "sign.h"
#include "masked_sign.h"
#include "masking_interface.h"
#include "masked_polyvec_operations.h"
#include "./test/cpucycles.h"
#include "gen_zr.h"
#pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
#pragma GCC diagnostic ignored "-Wint-conversion"



uint8_t g_msgs[BATCH_SIZE][MLEN] = {0};
uint8_t g_sigs[BATCH_SIZE][MLEN + CRYPTO_BYTES];

// -----------------------------------------------------------------------------
// Benchmark: masked batch signature generation (internal, used by batch benchmark)
// -----------------------------------------------------------------------------
int benchmark_masked_batch_sign_internal(
    size_t *sig_lens,
    size_t msg_len,
    uint8_t *seed_buf,
    polyvecl *s1,
    polyveck *s2,
    polyveck *t0,
    uint8_t *masked_key
)
{
    // BENCH INDEX: [0:NTT, 1:SAMPLE_Y, 2:AY, 3:DECOMPOSE, 4:Z=y+cs1, 5:REJECT, 6:w-cs2]
    uint64_t start, stop;
    start = cpucycles();

    unsigned int n;
    uint8_t *rho, *tr, *key;
    uint16_t nonce = 0;
    polyvecl mat[K], y, z;
    polyveck w1, w0, h;
    poly cp;
    keccak_state state;

    rho = seed_buf;
    tr = rho + SEEDBYTES;
    key = tr + CRHBYTES;

    uint8_t **mu = malloc(sizeof(uint8_t*) * BATCH_SIZE);
    uint8_t **rho_prime = malloc(sizeof(uint8_t*) * BATCH_SIZE);

    masked_polyvecl masked_y;
    masked_polyvecl masked_z;
    masked_polyveck masked_w1;
    masked_polyveck masked_w0;
    masked_polyveck masked_h;
    masked_polyvecl masked_s1;
    masked_polyveck masked_s2;
    poly T3[3][K-1], tmp_add;
    // Precompute T3 optimization terms
    for (int k = 0; k < 2; k++) {
        for (int i = 1; i < K; i++) {
            poly_add(&tmp_add, &mat[0].vec[2*k], &mat[i].vec[2*k]);
            poly_pointwise_montgomery(&T3[k][i-1], &mat[i].vec[2*k+1], &tmp_add);
        }
    }
    for (int i = 1; i < K; i++) {
        poly_add(&T3[2][i-1], &T3[1][i-1], &T3[0][i-1]);
    }
    // Precompute message hashes
    for (int i = 0; i < BATCH_SIZE; i++) {
        mu[i] = malloc(CRHBYTES);
        rho_prime[i] = malloc(CRHBYTES);

        shake256_init(&state);
        shake256_absorb(&state, tr, CRHBYTES);
        shake256_absorb(&state, g_msgs[i], msg_len);
        shake256_finalize(&state);
        shake256_squeeze(mu[i], CRHBYTES, &state);

        randombytes(rho_prime[i], CRHBYTES);
    }

    // Expand matrix
    polyvec_matrix_expand(mat, rho);
    polyveck_ntt(t0);
    int iter = 1;
    
    // Batch signing loop
    for (int i = 0; i < BATCH_SIZE; i++) {
        mask_polyvecl(&masked_s1, s1);
        mask_polyveck(&masked_s2, s2);
        masked_polyvecl_ntt(&masked_s1);
        masked_polyveck_ntt(&masked_s2);

        int flag = 0;
        while (!flag) {
            masked_sample_y(&masked_y);
            masked_z = masked_y;
            masked_polyvecl_ntt(&masked_z);
            
        
            // Matrix multiplication
            start = cpucycles();
            masked_polyvec_matrix_pointwise_montgomery_precom(&masked_w1, mat, &masked_z, T3);

            masked_polyveck_reduce(&masked_w1);
            masked_polyveck_invntt_tomont(&masked_w1);
            masked_polyveck_caddq(&masked_w1);

            masked_decompose(&w1, &masked_w0, &masked_w1);
            polyveck_pack_w1(g_sigs[i], &w1);

            // Challenge generation
            shake256_init(&state);
            shake256_absorb(&state, mu[i], CRHBYTES);
            shake256_absorb(&state, g_sigs[i], K * POLYW1_PACKEDBYTES);
            shake256_finalize(&state);
            shake256_squeeze(g_sigs[i], SEEDBYTES, &state);

            poly_challenge(&cp, g_sigs[i]);
            poly_ntt(&cp);

            // Compute z = y + c*s1
            masked_polyvecl_pointwise_poly_montgomery(&masked_z, &cp, &masked_s1);
            masked_polyvecl_invntt_tomont(&masked_z);
            masked_polyvecl_add(&masked_z, &masked_z, &masked_y);
            masked_polyvecl_reduce(&masked_z);

            // Reject z
            if (masked_rejection_sampling_z(&masked_z)) continue;

            // Compute w0 - c*s2
            masked_polyveck_pointwise_poly_montgomery(&masked_h, &cp, &masked_s2);
            masked_polyveck_invntt_tomont(&masked_h);
            masked_polyveck_sub(&masked_w0, &masked_w0, &masked_h);
            masked_polyveck_reduce(&masked_w0);

            // Reject w0
            if (masked_rejection_sampling_r(&masked_w0)) continue;

            // Unmask
            unmask_polyvecl(&masked_z, &z);
            unmask_polyveck(&masked_w0, &w0);

            // Process hint
            polyveck_pointwise_poly_montgomery(&h, &cp, t0);
            polyveck_invntt_tomont(&h);
            polyveck_reduce(&h);

            if (polyveck_chknorm(&h, GAMMA2)) continue;

            polyveck_add(&w0, &w0, &h);
            polyveck_caddq(&w0);
            n = polyveck_make_hint(&h, &w0, &w1);

            if (n > OMEGA) continue;

            // Pack signature
            pack_sig(g_sigs[i], g_sigs[i], &z, &h);
            sig_lens[i] = CRYPTO_BYTES + msg_len;
            flag = 1;
        }
    }

    // Cleanup
    for (int i = 0; i < BATCH_SIZE; i++) {
        free(rho_prime[i]);
        free(mu[i]);
    }
    free(rho_prime);
    free(mu);

    return iter;
}

// -----------------------------------------------------------------------------
// Benchmark: batch signature 
// -----------------------------------------------------------------------------
void benchmark_batch_sign(
    size_t *sig_lens,
    uint8_t *seed_buf,
    uint8_t *pk,
    uint8_t *sk)
{
    int count = 0;
    int ret;

    uint8_t *rho = seed_buf;
    uint8_t *tr = rho + SEEDBYTES;
    uint8_t *key = tr + CRHBYTES;

    polyvecl s1;
    polyveck s2, t0, t0_backup;
    const int ITERS = 50;

    // Key generation
    if ((ret = crypto_sign_keypair(pk, sk)) != 0) {
        fprintf(stderr, "crypto_sign_keypair failed: %d\n", ret);
        return;
    }

    memset(seed_buf, 0, 2 * SEEDBYTES + 3 * CRHBYTES);
    unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
    t0_backup = t0;

    uint8_t masked_key[SEEDBYTES * N_SHARES];
    mask_bitstring(masked_key, key, SEEDBYTES);
    // Benchmark loop
    uint64_t cycle_start = cpucycles();
    for (int i = 0; i < ITERS; i++) {
        t0 = t0_backup;
        benchmark_masked_batch_sign_internal(sig_lens, MLEN, seed_buf, &s1, &s2, &t0, masked_key);
    }
    uint64_t cycle_end = cpucycles();

    printf("\n=============================================\n");
    printf("        Benchmark: Batch Sign Result\n");
    printf("=============================================\n");
    printf("Avg repetitions: %f\n", (double)count / ITERS);
    printf("\n");
    printf("Batch serial total cycles: %.0f kcycles\n",
           (double)(cycle_end - cycle_start) / (ITERS * 1000));
    printf("=============================================\n\n");
}

// -----------------------------------------------------------------------------
// Benchmark: single masked signature (many iterations to average)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Benchmark: Serial masked signing (repeating single sign many times)
// -----------------------------------------------------------------------------
void benchmark_serial_masked_sign(void)
{
    int ret;
    uint8_t msg[MLEN] = {0};
    uint8_t sig[MLEN + CRYPTO_BYTES];
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    const int ITERS = 50;

    uint8_t seed_buf[2 * SEEDBYTES + 3 * CRHBYTES];
    uint8_t *rho = seed_buf;
    uint8_t *tr = rho + SEEDBYTES;
    uint8_t *key = tr + CRHBYTES;

    masked_polyvecl masked_s1;
    masked_polyveck masked_s2;
    polyvecl s1;
    polyveck s2, t0;

    crypto_sign_keypair(pk, sk);
    unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);

    uint8_t masked_key[SEEDBYTES * N_SHARES];
    mask_bitstring(masked_key, key, SEEDBYTES);

    uint64_t bench_vec[7] = {0};
    int count = 0;

    uint64_t cycle_start = cpucycles();
    for (int i = 0; i < ITERS; i++) {
        for (int j = 0; j < BATCH_SIZE; j++) {
            unpack_sk(rho, tr, key, &t0, &s1, &s2, sk);
            mask_polyvecl(&masked_s1, &s1);
            mask_polyveck(&masked_s2, &s2);
            randombytes(msg, MLEN);

            count += bench_masked_crypto_sign(
                sig, &ret, msg, MLEN,
                seed_buf, &masked_s1, &masked_s2, &t0, masked_key, bench_vec);
        }
    }
    uint64_t cycle_end = cpucycles();

    printf("\n=============================================\n");
    printf("     Benchmark: Serial Sign Performance\n");
    printf("=============================================\n");
    printf("Security order: %d\n", MASKING_ORDER);
    printf("Avg repetitions: %f\n", (double)count / (ITERS * BATCH_SIZE));
    printf("\n");
    printf("Total cycles: %.0f kcycles\n",
           (double)(cycle_end - cycle_start) / (ITERS * 1000));
    printf("=============================================\n\n");
}

// -----------------------------------------------------------------------------
// Matrix multiplication optimization test (without precomputation)
// -----------------------------------------------------------------------------
void benchmark_matrix_mult_optimized(void)
{
    polyvecl mat[K];
    polyvecl mat_shares[N_SHARES][K];
    masked_polyveck mw1, mw2;
    masked_polyvecl my;
    uint8_t rho[SEEDBYTES];

    randombytes(rho, SEEDBYTES);
    polyvec_matrix_expand(mat, rho);
    masked_sample_y(&my);
    masked_polyvecl_ntt(&my);
    const int ITERS = 100000;
    uint64_t cyc, orig_cycles = 0, opt_cycles = 0;

    // Benchmark
    for (int i = 0; i < ITERS; i++) {
        // Original
        cyc = cpucycles();
        for (int j = 0; j < N_SHARES; j++) {
            polyvec_matrix_pointwise_montgomery(&mw1.shares[j], mat, &my.shares[j]);
        }
        orig_cycles += cpucycles() - cyc;
        // Optimized
        opt_cycles += masked_polyvec_matrix_pointwise_montgomery_opt(&mw1, mat, &my);
    }

    printf("\n=============================================\n");
    printf("    Benchmark: Matrix Multiplication (Optimized)\n");
    printf("=============================================\n");
    printf("Optimized avg cycles: %f\n", (double)opt_cycles / ITERS);
    printf("Original avg cycles: %f\n", (double)orig_cycles / ITERS);
    printf("Speedup ratio: %.2f%%\n",
           (1.0 - (double)opt_cycles / orig_cycles) * 100);
    printf("=============================================\n\n");
}

// -----------------------------------------------------------------------------
// Matrix multiplication with precomputed T3 terms
// -----------------------------------------------------------------------------
void benchmark_matrix_mult_precomputed(void)
{
    polyvecl mat[K];
    polyvecl mat_shares[N_SHARES][K];
    masked_polyveck mw1, mw2;
    masked_polyvecl my;
    uint8_t rho[SEEDBYTES];

    randombytes(rho, SEEDBYTES);
    polyvec_matrix_expand(mat, rho);

    poly T3[3][K-1], tmp_add;
    masked_sample_y(&my);
    masked_polyvecl_ntt(&my);

    const int ITERS = 100000;

    // Precompute T3 optimization terms
    for (int k = 0; k < 2; k++) {
        for (int i = 1; i < K; i++) {
            poly_add(&tmp_add, &mat[0].vec[2*k], &mat[i].vec[2*k]);
            poly_pointwise_montgomery(&T3[k][i-1], &mat[i].vec[2*k+1], &tmp_add);
        }
    }
    for (int i = 1; i < K; i++) {
        poly_add(&T3[2][i-1], &T3[1][i-1], &T3[0][i-1]);
    }

    uint64_t cyc, orig_cycles = 0, opt_cycles = 0;

    // Benchmark
    for (int i = 0; i < ITERS; i++) {
        // Original
        cyc = cpucycles();
        for (int j = 0; j < N_SHARES; j++) {
            polyvec_matrix_pointwise_montgomery(&mw1.shares[j], mat, &my.shares[j]);
        }
        orig_cycles += cpucycles() - cyc;
        // Optimized
        opt_cycles += masked_polyvec_matrix_pointwise_montgomery_precom(&mw1, mat, &my, T3);
    }

    printf("\n=============================================\n");
    printf("    Benchmark: Matrix Multiplication (Precomputed)\n");
    printf("=============================================\n");
    printf("Optimized avg cycles: %f\n", (double)opt_cycles / ITERS);
    printf("Original avg cycles: %f\n", (double)orig_cycles / ITERS);
    printf("Speedup ratio: %.2f%%\n",
           (1.0 - (double)opt_cycles / orig_cycles) * 100);
    printf("=============================================\n\n");
}

// -----------------------------------------------------------------------------
// Main function
// -----------------------------------------------------------------------------
int main(void)
{
    srand(time(NULL));
    size_t sig_lens[BATCH_SIZE] = {0};
    uint8_t pk[CRYPTO_PUBLICKEYBYTES] = {0};
    uint8_t sk[CRYPTO_SECRETKEYBYTES] = {0};
    uint8_t seed_buf[2 * SEEDBYTES + 20 * CRHBYTES] = {0};

    // Initialize messages
    for (int i = 0; i < BATCH_SIZE; i++) {
        randombytes(g_msgs[i], MLEN);
        memcpy(g_sigs[i] + CRYPTO_BYTES, g_msgs[i], MLEN);
    }

    printf("=====================================================\n");
    printf("          Dilithium mode: %d\n", DILITHIUM_MODE);
    printf("=====================================================\n\n");

    // Run all benchmarks
    benchmark_batch_sign(sig_lens, seed_buf, pk, sk);
    benchmark_serial_masked_sign();
    benchmark_matrix_mult_optimized();
    benchmark_matrix_mult_precomputed();

    return 0;
}