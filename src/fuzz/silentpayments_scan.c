/***********************************************************************
 * Fuzz target for the silentpayments module's receiver scanning.      *
 *                                                                     *
 * Performs a real send -> scan roundtrip with fuzz-chosen keys, output *
 * counts (sweeping the label-candidate batch boundaries), label usage  *
 * and label_lookup behavior (absent / honest cache / hostile           *
 * always-match), then checks recipient-side invariants: every sent     *
 * output is found exactly once and is spendable with the returned      *
 * tweak, and nothing outside the transaction is ever reported.         *
 ***********************************************************************/
#ifdef ENABLE_MODULE_SILENTPAYMENTS

#include "../../include/secp256k1_silentpayments.h"

#define SP_FUZZ_MAX_SEND 4
#define SP_FUZZ_MAX_DECOYS 23
#define SP_FUZZ_MIN_SIZE 168

struct sp_fuzz_lookup_ctx {
    int always_match;                 /* hostile mode: claim every label is in the cache */
    unsigned char label33[33];        /* honest cache: the recipient's single label */
    unsigned char label_tweak[32];    /* tweak for the honest entry */
    unsigned char hostile_tweak[32];  /* arbitrary fuzz bytes returned in hostile mode */
};

static const unsigned char *sp_fuzz_label_lookup(const unsigned char *label33, const void *ctx_ptr) {
    const struct sp_fuzz_lookup_ctx *lctx = (const struct sp_fuzz_lookup_ctx *)ctx_ptr;
    if (lctx->always_match) {
        return lctx->hostile_tweak;
    }
    if (secp256k1_memcmp_var(label33, lctx->label33, 33) == 0) {
        return lctx->label_tweak;
    }
    return NULL;
}

FUZZ_TARGET(silentpayments_scan) {
    static secp256k1_context *ctx = NULL;
    unsigned char scan_seckey[32], spend_seckey[32], input_seckey[32];
    unsigned char outpoint_smallest[36];
    unsigned char ser_a[32], ser_b[32];
    secp256k1_pubkey scan_pubkey, spend_pubkey, address_spend_pubkey, input_pubkey;
    secp256k1_keypair input_keypair;
    secp256k1_silentpayments_label label;
    secp256k1_silentpayments_recipient recipients[SP_FUZZ_MAX_SEND];
    const secp256k1_silentpayments_recipient *recipient_ptrs[SP_FUZZ_MAX_SEND];
    secp256k1_xonly_pubkey generated_outputs[SP_FUZZ_MAX_SEND];
    secp256k1_xonly_pubkey *generated_output_ptrs[SP_FUZZ_MAX_SEND];
    secp256k1_xonly_pubkey tx_output_objs[SP_FUZZ_MAX_DECOYS + SP_FUZZ_MAX_SEND];
    const secp256k1_xonly_pubkey *tx_outputs[SP_FUZZ_MAX_DECOYS + SP_FUZZ_MAX_SEND];
    secp256k1_silentpayments_found_output found_output_objs[SP_FUZZ_MAX_DECOYS + SP_FUZZ_MAX_SEND];
    secp256k1_silentpayments_found_output *found_outputs[SP_FUZZ_MAX_DECOYS + SP_FUZZ_MAX_SEND];
    secp256k1_silentpayments_prevouts_summary prevouts_summary;
    struct sp_fuzz_lookup_ctx lookup_ctx;
    secp256k1_silentpayments_label_lookup lookup_fn;
    const secp256k1_keypair *keypair_ptrs[1];
    const unsigned char *seckey_ptrs[1];
    const secp256k1_pubkey *input_pubkey_ptrs[1];
    const secp256k1_xonly_pubkey *input_xonly_ptrs[1];
    secp256k1_xonly_pubkey input_xonly;
    uint32_t n_found = 0;
    size_t n_send, n_decoys, n_tx_outputs;
    int use_labels, lookup_null, hostile, taproot_input;
    unsigned int label_m;
    size_t i;
    int ret;

    if (size < SP_FUZZ_MIN_SIZE) {
        return;
    }
    if (ctx == NULL) {
        ctx = secp256k1_context_create(SECP256K1_CONTEXT_NONE);
        FUZZ_CHECK(ctx != NULL);
    }

    use_labels = data[0] & 1;
    lookup_null = (data[0] >> 1) & 1;
    hostile = !lookup_null && ((data[0] >> 2) & 1);
    taproot_input = (data[0] >> 3) & 1;
    n_send = 1 + ((data[0] >> 4) & 3);          /* 1..4 */
    n_decoys = data[1] % (SP_FUZZ_MAX_DECOYS + 1); /* 0..23: sweeps batch sizes 7/8/9, 15/16/17 */
    label_m = data[2];
    memcpy(scan_seckey, data + 3, 32);
    memcpy(spend_seckey, data + 35, 32);
    memcpy(input_seckey, data + 67, 32);
    memcpy(lookup_ctx.hostile_tweak, data + 99, 32);
    memcpy(outpoint_smallest, data + 131, 36);

    if (!secp256k1_ec_seckey_verify(ctx, scan_seckey)) return;
    if (!secp256k1_ec_seckey_verify(ctx, spend_seckey)) return;
    if (!secp256k1_ec_seckey_verify(ctx, input_seckey)) return;
    FUZZ_CHECK(secp256k1_ec_pubkey_create(ctx, &scan_pubkey, scan_seckey));
    FUZZ_CHECK(secp256k1_ec_pubkey_create(ctx, &spend_pubkey, spend_seckey));

    /* Recipient address: unlabeled, or labeled with label m (which is also the
     * single honest cache entry). */
    FUZZ_CHECK(secp256k1_silentpayments_recipient_label_create(ctx, &label, lookup_ctx.label_tweak, scan_seckey, label_m));
    FUZZ_CHECK(secp256k1_silentpayments_recipient_label_serialize(ctx, lookup_ctx.label33, &label));
    if (use_labels) {
        FUZZ_CHECK(secp256k1_silentpayments_recipient_create_labeled_spend_pubkey(ctx, &address_spend_pubkey, &spend_pubkey, &label));
    } else {
        address_spend_pubkey = spend_pubkey;
    }
    lookup_ctx.always_match = hostile;
    lookup_fn = lookup_null ? NULL : sp_fuzz_label_lookup;

    /* Sender: n_send outputs to our own address, from a single plain or taproot input. */
    for (i = 0; i < n_send; i++) {
        recipients[i].scan_pubkey = scan_pubkey;
        recipients[i].spend_pubkey = address_spend_pubkey;
        recipients[i].index = i;
        recipient_ptrs[i] = &recipients[i];
        generated_output_ptrs[i] = &generated_outputs[i];
    }
    if (taproot_input) {
        FUZZ_CHECK(secp256k1_keypair_create(ctx, &input_keypair, input_seckey));
        keypair_ptrs[0] = &input_keypair;
        ret = secp256k1_silentpayments_sender_create_outputs(ctx, generated_output_ptrs, recipient_ptrs, n_send,
            outpoint_smallest, keypair_ptrs, 1, NULL, 0);
    } else {
        seckey_ptrs[0] = input_seckey;
        ret = secp256k1_silentpayments_sender_create_outputs(ctx, generated_output_ptrs, recipient_ptrs, n_send,
            outpoint_smallest, NULL, 0, seckey_ptrs, 1);
    }
    /* Can only fail on (negligible-probability) hash edge cases; treat as bug. */
    FUZZ_CHECK(ret == 1);

    /* Transaction outputs: decoys first, real outputs last, so that with labels the
     * matches sit in the later label-candidate batch windows. */
    for (i = 0; i < n_decoys; i++) {
        unsigned char decoy_seckey[32];
        secp256k1_pubkey decoy_pubkey;
        size_t j;
        for (j = 0; j < 32; j++) {
            decoy_seckey[j] = data[167] ^ (unsigned char)(i * 17 + j);
        }
        if (!secp256k1_ec_seckey_verify(ctx, decoy_seckey)) return;
        FUZZ_CHECK(secp256k1_ec_pubkey_create(ctx, &decoy_pubkey, decoy_seckey));
        FUZZ_CHECK(secp256k1_xonly_pubkey_from_pubkey(ctx, &tx_output_objs[i], NULL, &decoy_pubkey));
    }
    for (i = 0; i < n_send; i++) {
        tx_output_objs[n_decoys + i] = generated_outputs[i];
    }
    n_tx_outputs = n_decoys + n_send;
    for (i = 0; i < n_tx_outputs; i++) {
        tx_outputs[i] = &tx_output_objs[i];
        found_outputs[i] = &found_output_objs[i];
    }

    /* Receiver: prevouts summary from the public form of the sender's input. */
    if (taproot_input) {
        FUZZ_CHECK(secp256k1_keypair_xonly_pub(ctx, &input_xonly, NULL, &input_keypair));
        input_xonly_ptrs[0] = &input_xonly;
        FUZZ_CHECK(secp256k1_silentpayments_recipient_prevouts_summary_create(ctx, &prevouts_summary,
            outpoint_smallest, input_xonly_ptrs, 1, NULL, 0));
    } else {
        FUZZ_CHECK(secp256k1_ec_pubkey_create(ctx, &input_pubkey, input_seckey));
        input_pubkey_ptrs[0] = &input_pubkey;
        FUZZ_CHECK(secp256k1_silentpayments_recipient_prevouts_summary_create(ctx, &prevouts_summary,
            outpoint_smallest, NULL, 0, input_pubkey_ptrs, 1));
    }

    FUZZ_CHECK(secp256k1_silentpayments_recipient_scan_outputs(ctx, found_outputs, &n_found,
        tx_outputs, n_tx_outputs, scan_seckey, &prevouts_summary, &spend_pubkey,
        lookup_fn, lookup_fn ? &lookup_ctx : NULL));

    /* Structural invariants (hold in every mode). */
    FUZZ_CHECK(n_found <= n_tx_outputs);
    for (i = 0; i < n_found; i++) {
        size_t j;
        int in_tx = 0;
        FUZZ_CHECK(secp256k1_xonly_pubkey_serialize(ctx, ser_a, &found_output_objs[i].output));
        for (j = 0; j < n_tx_outputs; j++) {
            FUZZ_CHECK(secp256k1_xonly_pubkey_serialize(ctx, ser_b, tx_outputs[j]));
            if (secp256k1_memcmp_var(ser_a, ser_b, 32) == 0) {
                in_tx = 1;
                break;
            }
        }
        FUZZ_CHECK(in_tx); /* scanning must never report an output not in the tx */
        FUZZ_CHECK(found_output_objs[i].found_with_label == 0 || found_output_objs[i].found_with_label == 1);
    }
    if (hostile) {
        /* A lying cache yields unpredictable (but crash-free and in-tx) results. */
        return;
    }

    if (use_labels && lookup_fn == NULL) {
        /* Labeled outputs are unfindable without a label cache. */
        FUZZ_CHECK(n_found == 0);
        return;
    }

    /* Honest modes: exactly the sent outputs are found, each spendable with
     * spend_seckey + tweak. */
    FUZZ_CHECK(n_found == n_send);
    for (i = 0; i < n_found; i++) {
        unsigned char full_seckey[32];
        secp256k1_pubkey full_pubkey;
        secp256k1_xonly_pubkey full_xonly;
        FUZZ_CHECK(found_output_objs[i].found_with_label == use_labels);
        memcpy(full_seckey, spend_seckey, 32);
        FUZZ_CHECK(secp256k1_ec_seckey_tweak_add(ctx, full_seckey, found_output_objs[i].tweak));
        FUZZ_CHECK(secp256k1_ec_pubkey_create(ctx, &full_pubkey, full_seckey));
        FUZZ_CHECK(secp256k1_xonly_pubkey_from_pubkey(ctx, &full_xonly, NULL, &full_pubkey));
        FUZZ_CHECK(secp256k1_xonly_pubkey_serialize(ctx, ser_a, &full_xonly));
        FUZZ_CHECK(secp256k1_xonly_pubkey_serialize(ctx, ser_b, &found_output_objs[i].output));
        FUZZ_CHECK(secp256k1_memcmp_var(ser_a, ser_b, 32) == 0);
        if (use_labels) {
            unsigned char found_label33[33];
            FUZZ_CHECK(secp256k1_silentpayments_recipient_label_serialize(ctx, found_label33, &found_output_objs[i].label));
            FUZZ_CHECK(secp256k1_memcmp_var(found_label33, lookup_ctx.label33, 33) == 0);
        }
    }
}

#endif /* ENABLE_MODULE_SILENTPAYMENTS */
