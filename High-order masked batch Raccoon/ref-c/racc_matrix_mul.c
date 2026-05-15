#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "plat_local.h"
#include "racc_core.h"
#include "polyr.h"
#include "mont32.h"
#include "mont64.h"
#include "racc_matrix_mul.h"

#if defined(RACCOON_128_2)||defined(RACCOON_128_4)||defined(RACCOON_128_8)||defined(RACCOON_128_16)||defined(RACCOON_128_32)

#define HALF_N (RACC_ELL / 2)

/**
 * @brief Matrix multiplication for RACC_ELL=4 (even dimension)
 * @param C Output matrix [K×D×N]
 * @param A Input matrix A [K×ELL×N]
 * @param B Input matrix B [ELL×D×N]
 * @return Execution cycles
 */
int64_t racc_matrix_mul(
    int64_t C[RACC_K][RACC_D][RACC_N],
    int64_t A[RACC_K][RACC_ELL][RACC_N],
    int64_t B[RACC_ELL][RACC_D][RACC_N])
{
    int i, j, k;
    int64_t tmp_0[RACC_N], tmp_1[RACC_N], tmp_2[RACC_N];
    int64_t cycles_start = plat_get_cycle();

    // Initialize output matrix (required for accumulation mode)
    for (i = 0; i < RACC_K; i++) {
        for (j = 0; j < RACC_D; j++) {
            polyr_zero(C[i][j]);
        }
    }

    // Precompute shared terms from matrix A (computed once)
    int64_t a_share[HALF_N][RACC_K-1][RACC_N];
    int64_t a_share_sum[RACC_K-1][RACC_N];

    for (k = 0; k < HALF_N; k++) {
        for (i = 1; i < RACC_K; i++) {
            polyr_add(tmp_0, A[0][2*k], A[i][2*k]);
            polyr_ntt_cmul(a_share[k][i-1], A[i][2*k+1], tmp_0);
        }
    }

    // Sum of shared terms (maps to original T3[2])
    for (i = 1; i < RACC_K; i++) {
        polyr_add(a_share_sum[i-1], a_share[0][i-1], a_share[1][i-1]);
    }

    // Compute column by column, accumulate directly to output
    for (j = 0; j < RACC_D; j++) {
        int64_t even_share_sum[RACC_N];
        int64_t odd_share_sum[RACC_N];
        polyr_zero(even_share_sum);
        polyr_zero(odd_share_sum);

        // Compute even shared terms and cross terms
        for (k = 0; k < HALF_N; k++) {
            // Even block shared term
            polyr_add(tmp_0, A[0][2*k], B[2*k+1][j]);
            polyr_ntt_cmul(tmp_1, B[2*k][j], tmp_0);
            polyr_add(even_share_sum, even_share_sum, tmp_1);

            // Cross term (reuse common value for all i)
            polyr_add(tmp_2, B[2*k+1][j], A[0][2*k]);
            for (i = 1; i < RACC_K; i++) {
                polyr_add(tmp_0, B[2*k][j], A[i][2*k+1]);
                polyr_add(tmp_1, tmp_2, A[i][2*k]);
                polyr_ntt_cmul(tmp_1, tmp_0, tmp_1);
                polyr_add(C[i][j], C[i][j], tmp_1);
            }
        }

        // Combine results in original algorithm order
        // First row initialization
        polyr_add(C[0][j], C[0][j], even_share_sum);

        // Non-first rows
        for (i = 1; i < RACC_K; i++) {
            polyr_sub(C[i][j], C[i][j], C[0][j]);
            polyr_sub(C[i][j], C[i][j], a_share_sum[i-1]);
        }

        // Compute and add odd shared terms
        for (k = 0; k < HALF_N; k++) {
            polyr_sub(tmp_0, A[0][2*k+1], B[2*k][j]);
            polyr_ntt_cmul(tmp_1, B[2*k+1][j], tmp_0);
            polyr_add(odd_share_sum, odd_share_sum, tmp_1);
        }

        // Final first row combination
        polyr_add(C[0][j], C[0][j], odd_share_sum);
    }

    int64_t cycles_elapsed = plat_get_cycle() - cycles_start;
    return cycles_elapsed;
}

#else

// General matrix multiplication macros
#define MAX_PAIRS (((RACC_K) + 1) / 2)
#define EVEN_PART (RACC_ELL - 3)
#define EVEN_HALF (EVEN_PART / 2)

/**
 * @brief Optimized matrix multiplication workspace
 * @note Uses union for memory sharing between A1B1 and A2B2 blocks
 */
typedef struct {
    union {
        // A1B1 block: global shared terms
        struct {
            int64_t share_01_10[RACC_N];
            int64_t share_02_20[RACC_N];
            int64_t share_12_21[RACC_N];

            int64_t pair_share_0[MAX_PAIRS][RACC_N];
            int64_t pair_share_1[MAX_PAIRS][RACC_N];
            int64_t pair_share_2[MAX_PAIRS][RACC_N];
        } a1b1;

        // A2B2 block: even part terms
        struct {
            int64_t even_share[RACC_K][EVEN_HALF][RACC_N];
            int64_t even_pre[RACC_D][RACC_N];
            int64_t even_extra[RACC_D][RACC_N];
        } a2b2;
    } u;
} MatMulWorkspace;

/**
 * @brief General matrix multiplication (odd dimension)
 * @param C Output matrix [K×D×N]
 * @param A Input matrix A [K×ELL×N]
 * @param B Input matrix B [ELL×D×N]
 * @return Execution cycles, -1 on allocation failure
 */
int64_t racc_matrix_mul(
    int64_t C[RACC_K][RACC_D][RACC_N],
    int64_t A[RACC_K][RACC_ELL][RACC_N],
    int64_t B[RACC_ELL][RACC_D][RACC_N])
{
    MatMulWorkspace *ws = (MatMulWorkspace *)calloc(1, sizeof(MatMulWorkspace));
    if (!ws) return -1;

    int i, j, p;
    int64_t tmp_0[RACC_N], tmp_1[RACC_N], tmp_2[RACC_N];
    int64_t cycles_start = plat_get_cycle();

    // Initialize output matrix
    for (i = 0; i < RACC_K; i++) {
        for (j = 0; j < RACC_D; j++) {
            polyr_zero(C[i][j]);
        }
    }

    // Stage 1: A1B1 block (compute and combine on-the-fly)
    // Precompute global shared terms (once for all j/i)
    polyr_ntt_cmul(ws->u.a1b1.share_01_10, A[1][0], A[0][1]);
    polyr_ntt_cmul(ws->u.a1b1.share_02_20, A[2][0], A[0][2]);
    polyr_ntt_cmul(ws->u.a1b1.share_12_21, A[2][1], A[1][2]);

    // Precompute pair shared terms (once for all j)
    for (i = 3, p = 0; i < RACC_K - 1; i += 2, p++) {
        int r0 = i, r1 = i + 1;

        polyr_sub(tmp_0, A[0][1], A[r0][1]);
        polyr_sub(tmp_1, A[r0][0], A[1][0]);
        polyr_sub(tmp_1, tmp_1, A[r1][0]);
        polyr_ntt_cmul(ws->u.a1b1.pair_share_0[p], tmp_0, tmp_1);

        polyr_sub(tmp_0, A[0][2], A[r0][2]);
        polyr_sub(tmp_1, A[r1][0], A[2][0]);
        polyr_ntt_cmul(ws->u.a1b1.pair_share_1[p], tmp_0, tmp_1);

        polyr_add(tmp_0, A[1][2], A[r0][2]);
        polyr_sub(tmp_0, tmp_0, A[r1][2]);
        polyr_sub(tmp_1, A[r1][1], A[2][1]);
        polyr_ntt_cmul(ws->u.a1b1.pair_share_2[p], tmp_0, tmp_1);
    }

    // Column-wise computation: no temporary storage
    for (j = 0; j < RACC_D; j++) {
        int64_t pre_0[RACC_N], pre_1[RACC_N], pre_2[RACC_N];
        int64_t cross_01[RACC_N], cross_02[RACC_N], cross_12[RACC_N];

        // Precompute terms for current column j
        polyr_sub(tmp_0, A[0][0], A[1][0]);
        polyr_sub(tmp_0, tmp_0, A[2][0]);
        polyr_sub(tmp_0, tmp_0, B[1][j]);
        polyr_sub(tmp_0, tmp_0, B[2][j]);
        polyr_ntt_cmul(pre_0, B[0][j], tmp_0);

        polyr_sub(tmp_0, A[1][1], A[0][1]);
        polyr_sub(tmp_0, tmp_0, A[2][1]);
        polyr_sub(tmp_0, tmp_0, B[0][j]);
        polyr_sub(tmp_0, tmp_0, B[2][j]);
        polyr_ntt_cmul(pre_1, B[1][j], tmp_0);

        polyr_sub(tmp_0, A[2][2], A[0][2]);
        polyr_sub(tmp_0, tmp_0, A[1][2]);
        polyr_sub(tmp_0, tmp_0, B[0][j]);
        polyr_sub(tmp_0, tmp_0, B[1][j]);
        polyr_ntt_cmul(pre_2, B[2][j], tmp_0);

        // Cross terms for current column j
        polyr_add(tmp_0, B[0][j], A[0][1]);
        polyr_add(tmp_1, B[1][j], A[1][0]);
        polyr_ntt_cmul(cross_01, tmp_0, tmp_1);

        polyr_add(tmp_0, B[0][j], A[0][2]);
        polyr_add(tmp_1, B[2][j], A[2][0]);
        polyr_ntt_cmul(cross_02, tmp_0, tmp_1);

        polyr_add(tmp_0, B[1][j], A[1][2]);
        polyr_add(tmp_1, B[2][j], A[2][1]);
        polyr_ntt_cmul(cross_12, tmp_0, tmp_1);

        // Combine first 3 rows
        polyr_add(C[0][j], C[0][j], pre_0);
        polyr_add(C[0][j], C[0][j], cross_01);
        polyr_add(C[0][j], C[0][j], cross_02);
        polyr_sub(C[0][j], C[0][j], ws->u.a1b1.share_01_10);
        polyr_sub(C[0][j], C[0][j], ws->u.a1b1.share_02_20);

        polyr_add(C[1][j], C[1][j], pre_1);
        polyr_add(C[1][j], C[1][j], cross_01);
        polyr_add(C[1][j], C[1][j], cross_12);
        polyr_sub(C[1][j], C[1][j], ws->u.a1b1.share_01_10);
        polyr_sub(C[1][j], C[1][j], ws->u.a1b1.share_12_21);

        polyr_add(C[2][j], C[2][j], pre_2);
        polyr_add(C[2][j], C[2][j], cross_02);
        polyr_add(C[2][j], C[2][j], cross_12);
        polyr_sub(C[2][j], C[2][j], ws->u.a1b1.share_02_20);
        polyr_sub(C[2][j], C[2][j], ws->u.a1b1.share_12_21);

        // Pair-wise cross terms and combination
        for (i = 3, p = 0; i < RACC_K - 1; i += 2, p++) {
            int r0 = i, r1 = i + 1;
            int64_t pair_cross_0[RACC_N], pair_cross_1[RACC_N], pair_cross_2[RACC_N];

            // Compute pair cross terms
            polyr_add(tmp_0, B[0][j], A[0][1]);
            polyr_sub(tmp_0, tmp_0, A[r0][1]);
            polyr_sub(tmp_1, A[r0][0], B[1][j]);
            polyr_sub(tmp_1, tmp_1, A[1][0]);
            polyr_sub(tmp_1, tmp_1, A[r1][0]);
            polyr_ntt_cmul(pair_cross_0, tmp_0, tmp_1);

            polyr_add(tmp_0, B[0][j], A[0][2]);
            polyr_sub(tmp_0, tmp_0, A[r0][2]);
            polyr_sub(tmp_1, A[r1][0], B[2][j]);
            polyr_sub(tmp_1, tmp_1, A[2][0]);
            polyr_ntt_cmul(pair_cross_1, tmp_0, tmp_1);

            polyr_add(tmp_0, B[1][j], A[1][2]);
            polyr_add(tmp_0, tmp_0, A[r0][2]);
            polyr_sub(tmp_0, tmp_0, A[r1][2]);
            polyr_sub(tmp_1, A[r1][1], B[2][j]);
            polyr_sub(tmp_1, tmp_1, A[2][1]);
            polyr_ntt_cmul(pair_cross_2, tmp_0, tmp_1);

            // Combine to output
            polyr_add(C[r0][j], C[r0][j], cross_01);
            polyr_add(C[r0][j], C[r0][j], cross_02);
            polyr_add(C[r0][j], C[r0][j], pair_cross_0);
            polyr_add(C[r0][j], C[r0][j], pair_cross_1);
            polyr_sub(C[r0][j], C[r0][j], ws->u.a1b1.share_01_10);
            polyr_sub(C[r0][j], C[r0][j], ws->u.a1b1.share_02_20);
            polyr_sub(C[r0][j], C[r0][j], ws->u.a1b1.pair_share_0[p]);
            polyr_sub(C[r0][j], C[r0][j], ws->u.a1b1.pair_share_1[p]);

            polyr_add(C[r1][j], C[r1][j], cross_02);
            polyr_add(C[r1][j], C[r1][j], cross_12);
            polyr_add(C[r1][j], C[r1][j], pair_cross_1);
            polyr_add(C[r1][j], C[r1][j], pair_cross_2);
            polyr_sub(C[r1][j], C[r1][j], ws->u.a1b1.share_02_20);
            polyr_sub(C[r1][j], C[r1][j], ws->u.a1b1.share_12_21);
            polyr_sub(C[r1][j], C[r1][j], ws->u.a1b1.pair_share_1[p]);
            polyr_sub(C[r1][j], C[r1][j], ws->u.a1b1.pair_share_2[p]);
        }
    }

    // Stage 2: A2B2 even block computation
#define SRC_A(i, k) B[3 + (k)][i]
#define SRC_B(k, j) A[j][3 + (k)]

    // Even block shared terms
    for (i = 1; i < RACC_K; i++) {
        for (int k = 0; k < EVEN_HALF; k++) {
            int idx = 2 * k;
            polyr_add(tmp_0, SRC_B(idx, 0), SRC_B(idx, i));
            polyr_ntt_cmul(ws->u.a2b2.even_share[i][k], SRC_B(idx + 1, i), tmp_0);
        }
    }

    // Precompute even block terms
    for (j = 0; j < RACC_D; j++) {
        polyr_zero(ws->u.a2b2.even_pre[j]);
        polyr_zero(ws->u.a2b2.even_extra[j]);

        for (int k = 0; k < EVEN_HALF; k++) {
            int idx = 2 * k;
            polyr_add(tmp_0, SRC_B(idx, 0), SRC_A(j, idx + 1));
            polyr_ntt_cmul(tmp_2, SRC_A(j, idx), tmp_0);
            polyr_add(ws->u.a2b2.even_pre[j], ws->u.a2b2.even_pre[j], tmp_2);

            polyr_sub(tmp_0, SRC_B(idx + 1, 0), SRC_A(j, idx));
            polyr_ntt_cmul(tmp_2, SRC_A(j, idx + 1), tmp_0);
            polyr_add(ws->u.a2b2.even_extra[j], ws->u.a2b2.even_extra[j], tmp_2);
        }
    }

    // Combine even block results
    for (j = 0; j < RACC_D; j++) {
        // First row
        polyr_add(C[0][j], C[0][j], ws->u.a2b2.even_pre[j]);
        polyr_add(C[0][j], C[0][j], ws->u.a2b2.even_extra[j]);

        // Non-first rows
        for (i = 1; i < RACC_K; i++) {
            polyr_sub(C[i][j], C[i][j], ws->u.a2b2.even_pre[j]);

            for (int k = 0; k < EVEN_HALF; k++) {
                int idx = 2 * k;
                polyr_add(tmp_0, SRC_A(j, idx), SRC_B(idx + 1, i));
                polyr_add(tmp_1, SRC_B(idx, 0), SRC_B(idx, i));
                polyr_add(tmp_1, tmp_1, SRC_A(j, idx + 1));
                polyr_ntt_cmul(tmp_2, tmp_0, tmp_1);

                polyr_add(C[i][j], C[i][j], tmp_2);
                polyr_sub(C[i][j], C[i][j], ws->u.a2b2.even_share[i][k]);
            }
        }
    }

#undef SRC_A
#undef SRC_B

    int64_t cycles_elapsed = plat_get_cycle() - cycles_start;
    free(ws);
    return cycles_elapsed;
}

#endif
