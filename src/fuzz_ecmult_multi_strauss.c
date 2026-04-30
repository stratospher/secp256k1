/***********************************************************************
 * Distributed under the MIT software license, see the accompanying
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.
 ***********************************************************************/

/* Differential fuzz harness:
 *   secp256k1_ecmult_strauss_batch  vs  secp256k1_ecmult_multi_simple_var
 *
 * Both compute r = inp_g_sc*G + sum(scalars[i]*points[i]) and must agree.
 *
 * Build / run instructions match fuzz_ecmult.c — see CMakeLists.txt
 * (target: fuzz_ecmult_multi_strauss).
 */

#include <stdint.h>
#include <stdlib.h>

#include "secp256k1.c"
#include "assumptions.h"
#include "util.h"
#include "field_impl.h"
#include "group_impl.h"
#include "scalar_impl.h"
#include "ecmult_impl.h"
#include "scratch_impl.h"

#include "fuzz_ecmult_multi_common.h"

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    static fuzz_ecmult_multi_input in;
    secp256k1_callback err_cb = {NULL, NULL};
    secp256k1_scratch *scratch;
    secp256k1_gej r_test, r_ref;
    const secp256k1_scalar *g_sc;
    int ret_test, ret_ref;

    if (!fuzz_ecmult_multi_parse(&in, data, size)) return 0;
    g_sc = in.has_g ? &in.g_sc : NULL;

    scratch = secp256k1_scratch_create(&err_cb, FUZZ_ECMULT_MULTI_SCRATCH_SIZE);
    if (scratch == NULL) abort();

    /* Reference: simple loop. Cannot fail under valid inputs. */
    ret_ref = secp256k1_ecmult_multi_simple_var(&r_ref, g_sc,
                                                fuzz_ecmult_multi_cb, &in,
                                                in.n);

    /* Target: Strauss batch directly (cb_offset=0). */
    ret_test = secp256k1_ecmult_strauss_batch(&err_cb, scratch, &r_test,
                                              g_sc, fuzz_ecmult_multi_cb, &in,
                                              in.n, 0);

    CHECK(ret_ref == 1);
    CHECK(ret_test == 1);
    CHECK(secp256k1_gej_eq_var(&r_test, &r_ref));

    secp256k1_scratch_destroy(&err_cb, scratch);
    return 0;
}
