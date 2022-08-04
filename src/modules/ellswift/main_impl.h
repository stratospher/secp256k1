/***********************************************************************
 * Copyright (c) 2022 Pieter Wuille                                    *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_MODULE_ELLSWIFT_MAIN_H
#define SECP256K1_MODULE_ELLSWIFT_MAIN_H

#include "../../../include/secp256k1.h"
#include "../../../include/secp256k1_ellswift.h"
#include "../../hash.h"

/** c1 = the square root of -3 ((-3)**((p+1)/4)). */
static const secp256k1_fe secp256k1_ellswift_c1 = SECP256K1_FE_CONST(0x0a2d2ba9, 0x3507f1df, 0x233770c2, 0xa797962c, 0xc61f6d15, 0xda14ecd4, 0x7d8d27ae, 0x1cd5f852);
/** c2 = -1/2 * (c1 - 1). */
static const secp256k1_fe secp256k1_ellswift_c2 = SECP256K1_FE_CONST(0x7ae96a2b, 0x657c0710, 0x6e64479e, 0xac3434e9, 0x9cf04975, 0x12f58995, 0xc1396c28, 0x719501ef);

/** Decode ElligatorSwift encoding (u, t) to a fraction xn/xd representing a curve X coordinate. */
static void secp256k1_ellswift_fe2_to_gexfrac_var(secp256k1_fe* xn, secp256k1_fe* xd, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe v1 = *u, v2 = *t;
    secp256k1_fe v3, v4, v5, v6, v7, v8;
    secp256k1_fe_normalize_var(&v1);
    secp256k1_fe_normalize_var(&v2);
    if (secp256k1_fe_is_zero(&v1)) v1 = secp256k1_fe_one;
    if (secp256k1_fe_is_zero(&v2)) v2 = secp256k1_fe_one;
    secp256k1_fe_sqr(&v3, &v1);
    secp256k1_fe_mul(&v3, &v3, &v1);
    secp256k1_fe_add(&v3, &secp256k1_fe_const_b);
    secp256k1_fe_sqr(&v4, &v2);
    v5 = v3;
    secp256k1_fe_add(&v5, &v4);
    if (secp256k1_fe_normalizes_to_zero_var(&v5)) {
        secp256k1_fe_add(&v2, &v2);
        secp256k1_fe_sqr(&v4, &v2);
        v5 = v3;
        secp256k1_fe_add(&v5, &v4);
    }
    secp256k1_fe_mul(&v6, &v1, &secp256k1_ellswift_c1);
    secp256k1_fe_negate(&v4, &v4, 1);
    secp256k1_fe_add(&v4, &v3);
    secp256k1_fe_mul(&v4, &v4, &v6);
    secp256k1_fe_mul(&v2, &v2, &v6);
    secp256k1_fe_sqr(&v2, &v2);
    secp256k1_fe_sqr(&v8, &v5);
    secp256k1_fe_mul(&v3, &v1, &v2);
    secp256k1_fe_add(&v3, &v8);
    secp256k1_fe_sqr(&v6, &v2);
    secp256k1_fe_sqr(&v6, &v6);
    secp256k1_fe_mul_int(&v6, 7);
    secp256k1_fe_sqr(&v7, &v3);
    secp256k1_fe_mul(&v7, &v7, &v3);
    secp256k1_fe_mul(&v7, &v7, &v2);
    secp256k1_fe_add(&v7, &v6);
    if (secp256k1_fe_jacobi_var(&v7) >= 0) {
        *xn = v3;
        *xd = v2;
        return;
    }
    secp256k1_fe_mul(&v1, &v1, &v5);
    secp256k1_fe_add(&v1, &v4);
    secp256k1_fe_half(&v1);
    secp256k1_fe_negate(&v1, &v1, 3);
    secp256k1_fe_sqr(&v6, &v8);
    secp256k1_fe_mul_int(&v6, 7);
    secp256k1_fe_sqr(&v7, &v1);
    secp256k1_fe_mul(&v7, &v7, &v1);
    secp256k1_fe_mul(&v7, &v7, &v5);
    secp256k1_fe_add(&v7, &v6);
    *xd = v5;
    secp256k1_fe_inv_var(&v5, &v5);
    if (secp256k1_fe_jacobi_var(&v7) >= 0) {
        *xn = v1;
        return;
    }
    secp256k1_fe_add(&v1, &v4);
    *xn = v1;
}

/** Decode ElligatorSwift encoding (u, t) to X coordinate. */
static void secp256k1_ellswift_fe2_to_gex_var(secp256k1_fe* x, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe xn, xd;
    secp256k1_ellswift_fe2_to_gexfrac_var(&xn, &xd, u, t);
    secp256k1_fe_inv_var(&xd, &xd);
    secp256k1_fe_mul(x, &xn, &xd);
}

/** Decode ElligatorSwift encoding (u, t) to point P. */
static void secp256k1_ellswift_fe2_to_ge_var(secp256k1_ge* p, const secp256k1_fe* u, const secp256k1_fe* t) {
    secp256k1_fe x;
    secp256k1_ellswift_fe2_to_gex_var(&x, u, t);
    secp256k1_ge_set_xo_var(p, &x, secp256k1_fe_is_odd(t));//todo: hmm?
}

/* Try to complete an ElligatorSwift encoding (u, t) for X coordinate x, given u and x.
 *
 * There may be up to 8 distinct t values such that (u, t) decodes back to x, but also
 * fewer, or none at all. Each such partial inverse can be accessed individually using a
 * distinct input argument i (in range 0-7), and some or all of these may return failure.
 * The following guarantees exist:
 * - Given (x, u), no two distinct i values give the same successful result t.
 * - Every successful result maps back to x through secp256k1_ellswift_fe2_to_gex_var.
 * - Given (x, u), all t values that map back to x can be reached by combining the
 *   successful results from this function over all i values, with the exception of:
 *   - this function cannot be called with u=0
 *   - no result with t=0 will be returned
 *   - no result for which u^3 + t^2 + 7 = 0 will be returned.
 */
static void print_buf_1(const unsigned char *buf, size_t n) {
    size_t i;
    for (i = 0; i < n; i++) {
        printf("%02x", buf[i]);
        if(i % 4 == 3) {
            printf(" ");
        }
    }
    printf("\n");
}
static void print_fe_1(const secp256k1_fe *inp) {
    unsigned char value[32];
    secp256k1_fe_get_b32(value, inp);
    print_buf_1(value, 32);

#ifdef VERIFY
    printf(", %d (mag), %d (normal)\n", inp->magnitude, inp->normalized);
#endif
}
static int secp256k1_ellswift_fegex_to_fe_var(secp256k1_fe* t, const secp256k1_fe* x, const secp256k1_fe* u, int i) {
    secp256k1_fe xm = *x, um = *u;
    secp256k1_fe g, s, w2, w;
    secp256k1_fe_normalize_weak(&xm);
    secp256k1_fe_normalize_weak(&um);
    secp256k1_fe_sqr(&g, u);
    secp256k1_fe_mul(&g, &g, u);
    secp256k1_fe_add(&g, &secp256k1_fe_const_b);
    if ((i & 2) == 0) {
        secp256k1_fe o;
        s = xm;
        secp256k1_fe_add(&s, &um);
        secp256k1_fe_sqr(&o, &s);
        secp256k1_fe_mul(&o, &o, &s);
        secp256k1_fe_negate(&o, &o, 1);
        secp256k1_fe_add(&o, &secp256k1_fe_const_b);
        if (secp256k1_fe_jacobi_var(&o) >= 0){
//            printf("hits case 1");
            return 0;
        }
        if (i & 1) {
            secp256k1_fe_add(&xm, &um);
            secp256k1_fe_negate(&xm, &xm, 2);
//            printf("hits case 2");
        }
        o = um;
        secp256k1_fe_add(&o, &xm);
        secp256k1_fe_sqr(&o, &o);
        secp256k1_fe_negate(&o, &o, 1);
        secp256k1_fe_mul(&w2, &um, &xm);
        secp256k1_fe_add(&w2, &o);
        secp256k1_fe_inv_var(&w2, &w2);
        secp256k1_fe_mul(&w2, &w2, &g);
//        printf("hits case 3");
    } else {
        secp256k1_fe r2, r;
        secp256k1_fe_negate(&w2, &um, 1);
        secp256k1_fe_add(&w2, &xm);
        if (secp256k1_fe_normalizes_to_zero_var(&w2)){
//            printf("hits case 4");
            return 0;
        }
        secp256k1_fe_normalize_weak(&g);
        secp256k1_fe_mul_int(&g, 4);
        secp256k1_fe_sqr(&r2, &um);
        secp256k1_fe_mul_int(&r2, 3);
        secp256k1_fe_mul(&r2, &r2, &w2);
        secp256k1_fe_add(&r2, &g);
        secp256k1_fe_mul(&r2, &r2, &w2);
        secp256k1_fe_negate(&r2, &r2, 1);
        if (!secp256k1_fe_sqrt(&r, &r2)){
//            printf("hits case 5");
            return 0;
        }
        if (i & 1) {
            if (secp256k1_fe_normalizes_to_zero_var(&r)){
//                printf("hits case 6");
                return 0;
            }
            secp256k1_fe_negate(&r, &r, 1);
//            printf("hits case 7");
        }
        secp256k1_fe_inv_var(&xm, &w2);
        secp256k1_fe_mul(&xm, &xm, &r);
        secp256k1_fe_add(&xm, &um);
        secp256k1_fe_half(&xm);
        secp256k1_fe_negate(&xm, &xm, 2);
//        printf("hits case 8");
    }
    if (!secp256k1_fe_sqrt(&w, &w2)){
//        printf("hits case 9");
        return 0;
    }
    if (i & 4){
//        printf("hits case 10");
        secp256k1_fe_negate(&w, &w, 1);
    }
    secp256k1_fe_normalize_var(&xm);
    secp256k1_fe_mul(&um, &um, &secp256k1_ellswift_c2);
    secp256k1_fe_add(&um, &xm);
    secp256k1_fe_mul(t, &w, &um);
    return 1;
}

/** Find an ElligatorSwift encoding (u, t) for X coordinate x.
 *
 * hasher is a SHA256 object which a incrementing 4-byte counter is added to to
 * generate randomness for the rejection sampling in this function. Its size plus
 * 4 (for the counter) plus 9 (for the SHA256 padding) must be a multiple of 64
 * for efficiency reasons.
 */
static void secp256k1_ellswift_gex_to_fe2_var(secp256k1_fe* u, secp256k1_fe* t, const secp256k1_fe* x, const secp256k1_sha256* hasher) {
    /* Pool of 3-bit branch values. */
    unsigned char branch_hash[32];
    /* Number of 3-bit values in branch_hash left. */
    int branches_left = 0;
    /* Field elements u and branch values are extracted from
     * SHA256(hasher || cnt) for consecutive values of cnt. cnt==0
     * is first used to populate a pool of 64 4-bit branch values. The 64 cnt
     * values that follow are used to generate field elements u. cnt==65 (and
     * multiples thereof) are used to repopulate the pool and start over, if
     * that were ever necessary. */
    uint32_t cnt = 0;
    VERIFY_CHECK((hasher->bytes + 4 + 9) % 64 == 0);
    while (1) {
        int branch;
        /* If the pool of branch values is empty, populate it. */
        if (branches_left == 0) {
            printf("%d:%d\n",cnt, branches_left);
            secp256k1_sha256 hash = *hasher;
            unsigned char buf4[4];
            buf4[0] = cnt;
            buf4[1] = cnt >> 8;
            buf4[2] = cnt >> 16;
            buf4[3] = cnt >> 24;
            ++cnt;
//            printf("bytesbefore=%d\n",hash.bytes);
            secp256k1_sha256_write(&hash, buf4, 4);
            secp256k1_sha256_finalize(&hash, branch_hash);
//            print_buf_1(&branch_hash, 32);
//            printf("bytes=%d\n",hash.bytes);
            branches_left = 64;
            printf("%d:%d\n",cnt, branches_left);
        }
        /* Take a 3-bit branch value from the branch pool (top bit is discarded). */
        --branches_left;
        branch = (branch_hash[branches_left >> 1] >> ((branches_left & 1) << 2)) & 7;
        /* Compute a new u value by hashing. */
        {
            secp256k1_sha256 hash = *hasher;
            unsigned char buf4[4];
            unsigned char u32[32];
            buf4[0] = cnt;
            buf4[1] = cnt >> 8;
            buf4[2] = cnt >> 16;
            buf4[3] = cnt >> 24;
            ++cnt;
            printf("%d:%d\n",cnt, branches_left);
            secp256k1_sha256_write(&hash, buf4, 4);
            secp256k1_sha256_finalize(&hash, u32);
//            print_buf_1(&u32, 32);
            if (!secp256k1_fe_set_b32(u, u32)) continue;
            if (secp256k1_fe_is_zero(u)) continue;
        }
        /* Find a remainder t, and return it if found. */
        if (secp256k1_ellswift_fegex_to_fe_var(t, x, u, branch)) {
            secp256k1_fe_normalize_var(t);
            break;
        }
    }
}

/** Find an ElligatorSwift encoding (u, t) for point P. */
static void secp256k1_ellswift_ge_to_fe2_var(secp256k1_fe* u, secp256k1_fe* t, const secp256k1_ge* p, const secp256k1_sha256* hasher) {
    secp256k1_ellswift_gex_to_fe2_var(u, t, &p->x, hasher);
    if (secp256k1_fe_is_odd(t) != secp256k1_fe_is_odd(&p->y)) {
        secp256k1_fe_negate(t, t, 1);
        secp256k1_fe_normalize_var(t);
    }
}

int secp256k1_ellswift_encode(const secp256k1_context* ctx, unsigned char *ell64, const secp256k1_pubkey *pubkey, const unsigned char *rnd321) {
//    const unsigned char rnd32[32] = {
//            0xf8, 0xe7, 0x50, 0x37, 0xaa, 0xfa, 0x82, 0x3c,
//            0x76, 0xaa, 0x22, 0xf0, 0xa3, 0xa4, 0x90, 0xda,
//            0x31, 0x67, 0x93, 0xac, 0x33, 0xb0, 0xb3, 0xe8,
//            0x3a, 0x58, 0xbf, 0xe0, 0xc0, 0x27, 0x49, 0x59,
//    };
    unsigned char rnd32[32] = {
            0xbf, 0x72, 0x23, 0x9a, 0x41, 0xa4, 0xb8, 0x46,
            0x20, 0x46, 0xe9, 0x99, 0x9d, 0x7a, 0x7f, 0xf5,
            0x33, 0x10, 0x66, 0xcf, 0x81, 0x67, 0x37, 0x35,
            0xa6, 0xf0, 0x15, 0xb2, 0x48, 0xd2, 0x0a, 0xbf,
    };
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(ell64 != NULL);
    ARG_CHECK(pubkey != NULL);
    ARG_CHECK(rnd321 != NULL);
    if (secp256k1_pubkey_load(ctx, &p, pubkey)) {
        static const unsigned char PREFIX[128 - 9 - 4 - 32 - 33] = "secp256k1_ellswift_encode";
        secp256k1_fe u, t;
        unsigned char p33[33];
//        unsigned char branch_hash[32];
        secp256k1_sha256 hash, hash2;

        /* Set up hasher state */
        secp256k1_sha256_initialize(&hash);
        secp256k1_sha256_write(&hash, PREFIX, sizeof(PREFIX));
//        printf("printing PREFIX\n");
//        print_buf_1(&PREFIX, 128 - 9 - 4 - 32 - 33);
//        printf("printing 32 bytes\n");
//        print_buf_1(&rnd32, 32);
        secp256k1_sha256_write(&hash, rnd32, 32);
//        printf("printing 33 bytes\n");
        secp256k1_fe_get_b32(p33, &p.x);
        p33[32] = secp256k1_fe_is_odd(&p.y);
//        print_buf_1(&p33, 33);
        secp256k1_sha256_write(&hash, p33, sizeof(p33));
//        printf("final hash\n");
//        hash2 = hash;
//        secp256k1_sha256_finalize(&hash2, branch_hash);
//        print_buf_1(&branch_hash, 32);
        printf("it's gonna be lit tonight\n");
        VERIFY_CHECK(hash.bytes == 128 - 9 - 4);
//        print_buf_1(&hash.buf, hash.bytes);

        /* Compute ElligatorSwift encoding and construct output. */
        secp256k1_ellswift_ge_to_fe2_var(&u, &t, &p, &hash);
        secp256k1_fe_get_b32(ell64, &u);
        secp256k1_fe_get_b32(ell64 + 32, &t);
//        print_buf_1(&ell64, 64);
        return 1;
    }
    /* Only returned in case the provided pubkey is invalid. */
    return 0;
}

int secp256k1_ellswift_create(const secp256k1_context* ctx, unsigned char *ell64, const unsigned char *seckey321, const unsigned char *rnd321) {

    unsigned char seckey32[32] = {
            0x07, 0xcf, 0x4f, 0x05, 0x15, 0xac, 0x8a, 0x95,
            0xb3, 0x3f, 0x75, 0x5f, 0xe0, 0xa8, 0xa1, 0x76,
            0xe2, 0x2a, 0xb0, 0xf1, 0x96, 0x89, 0x96, 0x59,
            0x2a, 0x8e, 0x95, 0x1b, 0xda, 0x72, 0x4a, 0x1d,
    };
    unsigned char rnd32[32] = {
            0x1d, 0x5d, 0xa9, 0xfa, 0x8b, 0xd2, 0xd7, 0x30,
            0x4b, 0x18, 0xe1, 0x46, 0xec, 0x74, 0xb1, 0x0c,
            0xe7, 0x58, 0xf9, 0x0d, 0x24, 0x3d, 0x8f, 0xde,
            0x07, 0x65, 0x56, 0xce, 0xf6, 0xd4, 0x3a, 0xe4,
    };
    secp256k1_ge p;
    secp256k1_fe u, t;
    secp256k1_sha256 hash;
    secp256k1_scalar seckey_scalar;
    static const unsigned char PREFIX[32] = "secp256k1_ellswift_create";
    static const unsigned char ZERO[32] = {0};
    int ret = 0;

    /* Sanity check inputs. */
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(ell64 != NULL);
    memset(ell64, 0, 64);
    ARG_CHECK(secp256k1_ecmult_gen_context_is_built(&ctx->ecmult_gen_ctx));
    ARG_CHECK(seckey321 != NULL);

    /* Compute (affine) public key */
    ret = secp256k1_ec_pubkey_create_helper(&ctx->ecmult_gen_ctx, &seckey_scalar, &p, seckey32);
    secp256k1_fe_normalize_var(&p.x);
    secp256k1_fe_normalize_var(&p.y);

    /* Set up hasher state */
    secp256k1_sha256_initialize(&hash);
    secp256k1_sha256_write(&hash, PREFIX, sizeof(PREFIX));
    secp256k1_sha256_write(&hash, seckey32, 32);
    secp256k1_sha256_write(&hash, rnd32, 32);
    secp256k1_sha256_write(&hash, ZERO, 32 - 9 - 4);

    /* Compute ElligatorSwift encoding and construct output. */
    secp256k1_ellswift_ge_to_fe2_var(&u, &t, &p, &hash);
    secp256k1_fe_get_b32(ell64, &u);
    secp256k1_fe_get_b32(ell64 + 32, &t);

    secp256k1_memczero(ell64, 64, !ret);
    secp256k1_scalar_clear(&seckey_scalar);

    return ret;
}

int secp256k1_ellswift_decode(const secp256k1_context* ctx, secp256k1_pubkey *pubkey, const unsigned char *ell64) {
    secp256k1_fe u, t;
    secp256k1_ge p;
    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(pubkey != NULL);
    ARG_CHECK(ell64 != NULL);

    secp256k1_fe_set_b32(&u, ell64);
    secp256k1_fe_normalize_var(&u);
    secp256k1_fe_set_b32(&t, ell64 + 32);
    secp256k1_fe_normalize_var(&t);
    secp256k1_ellswift_fe2_to_ge_var(&p, &u, &t);
    secp256k1_pubkey_save(pubkey, &p);
    return 1;
}

static int ellswift_xdh_hash_function_sha256(unsigned char *output, const unsigned char *x32, const unsigned char *ours64, const unsigned char *theirs64, void *data) {
    secp256k1_sha256 sha;

    (void)data;

    secp256k1_sha256_initialize(&sha);
    if (secp256k1_memcmp_var(ours64, theirs64, 64) <= 0) {
        secp256k1_sha256_write(&sha, ours64, 64);
        secp256k1_sha256_write(&sha, theirs64, 64);
    } else {
        secp256k1_sha256_write(&sha, theirs64, 64);
        secp256k1_sha256_write(&sha, ours64, 64);
    }
    secp256k1_sha256_write(&sha, x32, 32);
    secp256k1_sha256_finalize(&sha, output);

    return 1;
}

const secp256k1_ellswift_xdh_hash_function secp256k1_ellswift_xdh_hash_function_sha256 = ellswift_xdh_hash_function_sha256;
const secp256k1_ellswift_xdh_hash_function secp256k1_ellswift_xdh_hash_function_default = ellswift_xdh_hash_function_sha256;

int secp256k1_ellswift_xdh(const secp256k1_context* ctx, unsigned char *output, const unsigned char* theirs64, const unsigned char* ours64, const unsigned char* seckey32, secp256k1_ellswift_xdh_hash_function hashfp, void *data) {
    int ret = 0;
    int overflow;
    secp256k1_scalar s;
    secp256k1_fe xn, xd, px, u, t;
    unsigned char sx[32];

    VERIFY_CHECK(ctx != NULL);
    ARG_CHECK(output != NULL);
    ARG_CHECK(theirs64 != NULL);
    ARG_CHECK(ours64 != NULL);
    ARG_CHECK(seckey32 != NULL);

    if (hashfp == NULL) {
        hashfp = secp256k1_ellswift_xdh_hash_function_default;
    }

    /* Load remote public key (as fraction). */
    secp256k1_fe_set_b32(&u, theirs64);
    secp256k1_fe_normalize_var(&u);
    secp256k1_fe_set_b32(&t, theirs64 + 32);
    secp256k1_fe_normalize_var(&t);
    secp256k1_ellswift_fe2_to_gexfrac_var(&xn, &xd, &u, &t);

    /* Load private key (using one if invalid). */
    secp256k1_scalar_set_b32(&s, seckey32, &overflow);
    overflow = secp256k1_scalar_is_zero(&s);
    secp256k1_scalar_cmov(&s, &secp256k1_scalar_one, overflow);

    /* Compute shared X coordinate. */
    secp256k1_ecmult_const_xonly(&px, &xn, &xd, &s, 256, 1);
    secp256k1_fe_normalize(&px);
    secp256k1_fe_get_b32(sx, &px);

    /* Invoke hasher */
    ret = hashfp(output, sx, ours64, theirs64, data);

    memset(sx, 0, 32);
    secp256k1_fe_clear(&px);
    secp256k1_scalar_clear(&s);

    return !!ret & !overflow;
}

#endif
