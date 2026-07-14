/***********************************************************************
 * libsecp256k1 fuzz harness entry point.                              *
 *                                                                     *
 * This is a unity build: it includes secp256k1.c (like tests.c) so    *
 * that fuzz targets can exercise the library's static internals, then *
 * includes each target file. A single binary is produced; the target  *
 * to run is selected with the FUZZ=<name> environment variable.       *
 ***********************************************************************/
#include <stdio.h>
#include <string.h>

#include "secp256k1.c"
#include "../include/secp256k1.h"
#include "util.h"
#include "modinv32_impl.h"
#include "modinv64_impl.h"
#include "int128_impl.h"

#include "fuzz.h"

/* ---- Target registry -------------------------------------------------- */

#define FUZZ_MAX_TARGETS 64
static struct { const char *name; secp256k1_fuzz_fn fn; } fuzz_targets[FUZZ_MAX_TARGETS];
static int fuzz_num_targets = 0;

int secp256k1_fuzz_register(const char *name, secp256k1_fuzz_fn fn) {
    if (fuzz_num_targets < FUZZ_MAX_TARGETS) {
        fuzz_targets[fuzz_num_targets].name = name;
        fuzz_targets[fuzz_num_targets].fn = fn;
        fuzz_num_targets++;
    }
    return 1;
}

/* ---- Targets (share secp256k1.c's internals via the unity build) ------ */

#include "silentpayments_scan.c"

/* ---- Target selection ------------------------------------------------- */

static secp256k1_fuzz_fn fuzz_selected = NULL;

static void fuzz_select_from_env(void) {
    const char *name = getenv("FUZZ");
    int i;
    if (name != NULL) {
        for (i = 0; i < fuzz_num_targets; i++) {
            if (strcmp(fuzz_targets[i].name, name) == 0) {
                fuzz_selected = fuzz_targets[i].fn;
                return;
            }
        }
        fprintf(stderr, "Unknown FUZZ target: %s\n", name);
    } else {
        fprintf(stderr, "Set the FUZZ environment variable to one of:\n");
    }
    for (i = 0; i < fuzz_num_targets; i++) {
        fprintf(stderr, "  %s\n", fuzz_targets[i].name);
    }
    exit(1);
}

/* ---- libFuzzer / OSS-Fuzz interface ----------------------------------- */

int LLVMFuzzerInitialize(int *argc, char ***argv);
int LLVMFuzzerInitialize(int *argc, char ***argv) {
    (void)argc;
    (void)argv;
    fuzz_select_from_env();
    return 0;
}

int LLVMFuzzerTestOneInput(const unsigned char *data, size_t size);
int LLVMFuzzerTestOneInput(const unsigned char *data, size_t size) {
    fuzz_selected(data, size);
    return 0;
}

#ifdef SECP256K1_FUZZ_STANDALONE
/* Built without a fuzzing engine: replay the given corpus files once each.
 * Lets CI run a corpus for regression without libFuzzer/AFL. */
static void fuzz_run_file(const char *path) {
    unsigned char buf[1 << 20];
    size_t len;
    FILE *f = fopen(path, "rb");
    if (f == NULL) {
        fprintf(stderr, "Cannot open %s\n", path);
        exit(1);
    }
    len = fread(buf, 1, sizeof(buf), f);
    fclose(f);
    fuzz_selected(buf, len);
}

int main(int argc, char **argv) {
    int i;
    fuzz_select_from_env();
    if (argc < 2) {
        fprintf(stderr, "usage: FUZZ=<target> %s <inputfile>...\n"
                        "(this is a replay build; for coverage-guided fuzzing "
                        "rebuild with -DSECP256K1_FUZZ_LIBFUZZER=ON)\n", argv[0]);
        return 0;
    }
    for (i = 1; i < argc; i++) {
        fuzz_run_file(argv[i]);
    }
    fprintf(stderr, "ran %d input(s), no invariant violated\n", argc - 1);
    return 0;
}
#endif
