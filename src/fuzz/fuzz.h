/***********************************************************************
 * Minimal Bitcoin Core-style fuzz framework for libsecp256k1.         *
 *                                                                     *
 * A single fuzz binary bundles many targets. Each target is defined   *
 * with FUZZ_TARGET() and auto-registers itself; the one to run is     *
 * chosen at runtime via the FUZZ=<name> environment variable, e.g.    *
 *                                                                     *
 *     FUZZ=field ./fuzz corpus/field                                  *
 *                                                                     *
 * Targets access libsecp256k1 internals, so they are compiled into    *
 * the same translation unit as secp256k1.c (see fuzz.c).              *
 ***********************************************************************/
#ifndef SECP256K1_FUZZ_H
#define SECP256K1_FUZZ_H

#include <stddef.h>
#include <stdlib.h>

typedef void (*secp256k1_fuzz_fn)(const unsigned char *data, size_t size);

/* Register a fuzz target by name. Returns 1 (see FUZZ_TARGET). */
int secp256k1_fuzz_register(const char *name, secp256k1_fuzz_fn fn);

/* Define and auto-register a fuzz target:
 *
 *     FUZZ_TARGET(field) {
 *         ... use data, size ...
 *     }
 *
 * Uses a constructor so registration happens before main(). This is a
 * clang/gcc extension, which is fine: fuzzing already requires clang. */
#define FUZZ_TARGET(tname)                                                   \
    static void fuzz_target_##tname(const unsigned char *data, size_t size); \
    __attribute__((constructor)) static void fuzz_reg_##tname(void) {        \
        secp256k1_fuzz_register(#tname, fuzz_target_##tname);                \
    }                                                                        \
    static void fuzz_target_##tname(const unsigned char *data, size_t size)

/* An invariant that must hold; a violation is a crash the fuzzer keeps. */
#define FUZZ_CHECK(cond) do { if (!(cond)) { abort(); } } while(0)

#endif /* SECP256K1_FUZZ_H */
