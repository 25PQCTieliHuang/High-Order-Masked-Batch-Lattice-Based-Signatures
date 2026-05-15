#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "plat_local.h"
#include "racc_core.h"
#include "polyr.h"
#include "mont32.h"
#include "mont64.h"

int64_t racc_matrix_mul(int64_t C[RACC_K][RACC_D][RACC_N], 
    int64_t A[RACC_K][RACC_ELL][RACC_N],  // 输入矩阵A（7×5）
    int64_t B[RACC_ELL][RACC_D][RACC_N]    // 输入矩阵B（5×d）
    );
