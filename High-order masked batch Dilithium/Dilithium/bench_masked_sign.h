#include <stdint.h>
#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"
#include "symmetric.h"
#include "fips202.h"
#include "./test/cpucycles.h"

#include "masked_sign.h"
#include "masking_interface.h"
#include "masked_polyvec_operations.h"
#include <execinfo.h>
#include <stdio.h>


int bench_masked_crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen, uint8_t* seedbuf,
                          masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0, uint8_t* masked_key, uint64_t* bench_vector);



int bench_masked_crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                uint8_t* seedbuf,
                masked_polyvecl* ms1, masked_polyveck* ms2, polyveck* t0, uint8_t* masked_key, uint64_t* bench_vector);

                          
                          


