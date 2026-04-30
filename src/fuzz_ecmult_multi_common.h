/***********************************************************************
 * Distributed under the MIT software license, see the accompanying
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.
 ***********************************************************************/

/* Shared input parsing and helpers for the ecmult_multi differential
 * fuzz harnesses (Strauss, Pippenger, dispatcher).
 *
 * Each harness includes "secp256k1.c" first, then this header.
 */

#ifndef SECP256K1_FUZZ_ECMULT_MULTI_COMMON_H
#define SECP256K1_FUZZ_ECMULT_MULTI_COMMON_H

#include <stddef.h>
#include <stdint.h>

#include "scalar.h"
#include "field.h"
#include "group.h"

#define FUZZ_ECMULT_MULTI_MAX_POINTS 128
#define FUZZ_ECMULT_MULTI_PER_POINT_BYTES 65 /* 32 scalar + 32 x + 1 flag */
#define FUZZ_ECMULT_MULTI_HEADER_BYTES 35    /* 1 ctrl + 32 g_sc + 1 n + 1 pad */
#define FUZZ_ECMULT_MULTI_SCRATCH_SIZE (8 * 1024 * 1024)

typedef struct {
    secp256k1_scalar scalars[FUZZ_ECMULT_MULTI_MAX_POINTS];
    secp256k1_ge points[FUZZ_ECMULT_MULTI_MAX_POINTS];
    size_t n;
    int has_g;
    secp256k1_scalar g_sc;
} fuzz_ecmult_multi_input;

/* Parses a fuzzer input into the structured form. Returns 1 on success
 * (caller may run the targets), 0 if the input is too short. */
static int fuzz_ecmult_multi_parse(fuzz_ecmult_multi_input *out,
                                   const uint8_t *data, size_t size) {
    size_t i;
    size_t n_max_from_size;
    uint8_t ctrl;

    if (size < FUZZ_ECMULT_MULTI_HEADER_BYTES) return 0;

    ctrl = data[0];
    out->has_g = (ctrl & 1) ? 1 : 0;

    secp256k1_scalar_set_b32(&out->g_sc, data + 1, NULL);

    n_max_from_size = (size - FUZZ_ECMULT_MULTI_HEADER_BYTES)
                      / FUZZ_ECMULT_MULTI_PER_POINT_BYTES;
    if (n_max_from_size > FUZZ_ECMULT_MULTI_MAX_POINTS) {
        n_max_from_size = FUZZ_ECMULT_MULTI_MAX_POINTS;
    }
    out->n = data[33];
    if (out->n > n_max_from_size) out->n = n_max_from_size;

    for (i = 0; i < out->n; i++) {
        const uint8_t *p = data + FUZZ_ECMULT_MULTI_HEADER_BYTES
                         + i * FUZZ_ECMULT_MULTI_PER_POINT_BYTES;
        secp256k1_fe x;
        secp256k1_scalar_set_b32(&out->scalars[i], p, NULL);
        if (!secp256k1_fe_set_b32_limit(&x, p + 32) ||
            !secp256k1_ge_set_xo_var(&out->points[i], &x, p[64] & 1)) {
            /* Fall back to G when x is out of range or has no y. */
            out->points[i] = secp256k1_ge_const_g;
        }
    }
    return 1;
}

/* Callback that pulls (scalar, point) pairs out of a parsed input. */
static int fuzz_ecmult_multi_cb(secp256k1_scalar *sc, secp256k1_ge *pt,
                                size_t idx, void *cbdata) {
    const fuzz_ecmult_multi_input *in = (const fuzz_ecmult_multi_input *)cbdata;
    if (idx >= in->n) return 0;
    *sc = in->scalars[idx];
    *pt = in->points[idx];
    return 1;
}

#endif /* SECP256K1_FUZZ_ECMULT_MULTI_COMMON_H */
