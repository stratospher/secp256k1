/***********************************************************************
 * Distributed under the MIT software license, see the accompanying
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.
 ***********************************************************************/

/* Fuzz harness for secp256k1_ecmult (R = na*A + ng*G).
 *
 * Build with libFuzzer (requires Homebrew LLVM on macOS; Apple's CLT clang
 * does not ship the fuzzer runtime):
 *   brew install llvm   # macOS only; skip on Linux
 *   cmake -B build-fuzz -DSECP256K1_BUILD_FUZZ_TESTS=ON \
 *         -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang ..
 *   cmake --build build-fuzz --target fuzz_ecmult
 *
 * Input layout (INPUT_SIZE bytes; inputs shorter than this are skipped):
 *
 *   bytes   0– 31  na       scalar multiplier for point A
 *   bytes  32– 63  ng       scalar multiplier for generator G
 *   bytes  64– 95  x        x-coordinate of input point A
 *   byte        96  flags   bit 0 = y-parity of A; if x not on curve -> A = G
 *   bytes  97–128  z_bytes  Jacobian Z-rescale factor (skipped if 0 or >= p)
 *   bytes 129–160  na1      scalar split: na1 + na2 = na (na2 derived)
 */

#include <stdint.h>
#include <stdlib.h>

/* Pull in the full secp256k1 implementation and internal headers.
 * Follows the same pattern as bench_internal.c. */
#include "secp256k1.c"
#include "assumptions.h"
#include "util.h"
#include "field_impl.h"
#include "group_impl.h"
#include "scalar_impl.h"
#include "ecmult_impl.h"

/* INPUT_SIZE: minimum number of fuzz bytes required; inputs shorter than
 * this are skipped so the fuzzer focuses on fully-structured inputs. */
#define INPUT_SIZE 161

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    secp256k1_scalar na, ng, na1, na2;
    secp256k1_fe x, z;
    secp256k1_ge ge_a;
    secp256k1_gej a, r1, r2, r3, tmp;

    if (size < INPUT_SIZE) return 0;

    /* ------------------------------------------------------------------ *
     * Parse scalars.
     * secp256k1_scalar_set_b32 reduces mod n, so any 32 bytes are valid. *
     * ------------------------------------------------------------------ */
    secp256k1_scalar_set_b32(&na,  data + 0,   NULL);
    secp256k1_scalar_set_b32(&ng,  data + 32,  NULL);
    secp256k1_scalar_set_b32(&na1, data + 129, NULL);

    /* ------------------------------------------------------------------ *
     * Build input point A (Option C: x-lift to affine + fuzz-derived Z). *
     *                                                                     *
     * Step 1: lift x-coordinate to a curve point.                        *
     *   secp256k1_fe_set_b32_limit fails when x >= field prime p (~50%). *
     *   secp256k1_ge_set_xo_var fails when x has no square root for y.   *
     *   In either case fall back to G, which is always a valid point.    *
     * ------------------------------------------------------------------ */
    if (!secp256k1_fe_set_b32_limit(&x, data + 64) ||
            !secp256k1_ge_set_xo_var(&ge_a, &x, data[96] & 1)) {
        ge_a = secp256k1_ge_const_g;
    }

    /* Step 2: convert affine -> Jacobian with Z = 1. */
    secp256k1_gej_set_ge(&a, &ge_a);

    /* Step 3: rescale with a fuzz-derived Z to vary the Jacobian          *
     * representation.  This exercises magnitude/normalisation handling    *
     * inside ecmult without changing the affine point being multiplied.   *
     * secp256k1_gej_rescale requires Z != 0; also skip if Z >= p.        */
    if (secp256k1_fe_set_b32_limit(&z, data + 97) && !secp256k1_fe_is_zero(&z)) {
        secp256k1_gej_rescale(&a, &z);
    }

    /* ================================================================== *
     * Invariant 1: NULL ng equivalence                                   *
     *   ecmult(A, na, NULL) == ecmult(A, na, &zero)                      *
     *                                                                     *
     * NULL ng is a documented fast-path alias for the zero scalar.       *
     * Both calls must produce identical results.                          *
     * ================================================================== */
    secp256k1_ecmult(&r1, &a, &na, NULL);
    secp256k1_ecmult(&r2, &a, &na, &secp256k1_scalar_zero);
    CHECK(secp256k1_gej_eq_var(&r1, &r2));

    /* ================================================================== *
     * Invariant 2: Negation cancellation                                 *
     *   ecmult(A, na, ng) + ecmult(A, -na, -ng) == infinity              *
     *                                                                     *
     * For any point P, P + (-P) == O.  Tests the group inverse.          *
     * ================================================================== */
    {
        secp256k1_scalar neg_na, neg_ng;
        secp256k1_scalar_negate(&neg_na, &na);
        secp256k1_scalar_negate(&neg_ng, &ng);
        secp256k1_ecmult(&r1, &a, &na,     &ng);
        secp256k1_ecmult(&r2, &a, &neg_na, &neg_ng);
        secp256k1_gej_add_var(&tmp, &r1, &r2, NULL);
        CHECK(secp256k1_gej_is_infinity(&tmp));
    }

    /* ================================================================== *
     * Invariant 3: Scalar additivity for point multiplication            *
     *   na1*A + na2*A == na*A   (where na2 = na - na1)                   *
     *                                                                     *
     * Tests linearity in the point scalar.  na1+na2 == na holds exactly  *
     * because na2 is derived as na - na1 (no approximation).             *
     * ================================================================== */
    {
        secp256k1_scalar_negate(&na2, &na1);
        secp256k1_scalar_add(&na2, &na2, &na); /* na2 = na - na1 */

        secp256k1_ecmult(&r1,  &a, &na1, &secp256k1_scalar_zero);
        secp256k1_ecmult(&r2,  &a, &na2, &secp256k1_scalar_zero);
        secp256k1_gej_add_var(&tmp, &r1, &r2, NULL);
        secp256k1_ecmult(&r3,  &a, &na,  &secp256k1_scalar_zero);
        CHECK(secp256k1_gej_eq_var(&tmp, &r3));
    }

    /* ================================================================== *
     * Invariant 4: Scalar additivity for generator multiplication        *
     *   ng1*G + ng2*G == ng*G   (where ng1 = na1, ng2 = ng - ng1)       *
     *                                                                     *
     * Tests linearity in the generator scalar.  Passes NULL for the      *
     * point (with na=0) to isolate the ng*G path.                        *
     * ================================================================== */
    {
        secp256k1_scalar ng1, ng2;
        ng1 = na1; /* reuse fuzz bytes 129–160 */
        secp256k1_scalar_negate(&ng2, &ng1);
        secp256k1_scalar_add(&ng2, &ng2, &ng); /* ng2 = ng - ng1 */

        secp256k1_ecmult(&r1, NULL, &secp256k1_scalar_zero, &ng1);
        secp256k1_ecmult(&r2, NULL, &secp256k1_scalar_zero, &ng2);
        secp256k1_gej_add_var(&tmp, &r1, &r2, NULL);
        secp256k1_ecmult(&r3, NULL, &secp256k1_scalar_zero, &ng);
        CHECK(secp256k1_gej_eq_var(&tmp, &r3));
    }

    /* ================================================================== *
     * Invariant 5: Additive decomposition of the double multiply         *
     *   ecmult(A, na, 0) + ecmult(NULL, 0, ng) == ecmult(A, na, ng)     *
     *                                                                     *
     * The two independent single-scalar multiplications must sum to the  *
     * combined double-multiply result.                                    *
     * ================================================================== */
    secp256k1_ecmult(&r1, &a,   &na,                   &secp256k1_scalar_zero);
    secp256k1_ecmult(&r2, NULL, &secp256k1_scalar_zero, &ng);
    secp256k1_gej_add_var(&tmp, &r1, &r2, NULL);
    secp256k1_ecmult(&r3, &a,   &na,                   &ng);
    CHECK(secp256k1_gej_eq_var(&tmp, &r3));

    return 0;
}
