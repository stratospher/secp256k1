#!/usr/bin/env python3
"""Differential test for the silentpayments module against secp256k1lab.

Replays inputs in the same 168-byte format as the silentpayments_scan fuzz
target (share the corpus!): parses keys/config identically, runs the real C
library via ctypes AND theStack's secp256k1lab reference implementation, and
compares sender outputs and scan results byte-for-byte.

The hostile always-match mode is skipped (a lying label cache has no
reference-defined semantics); the C-side invariants in the fuzz target
already cover it.

Usage:
    python3 src/fuzz/diff_silentpayments_scan.py corpus/silentpayments_scan/*
    SECP256K1_LIB=path/to/libsecp256k1.dylib python3 src/fuzz/diff_silentpayments_scan.py <files>
"""
import ctypes
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "secp256k1lab", "src"))

from secp256k1lab.secp256k1 import G, GE, Scalar
from secp256k1lab import bip352

MIN_SIZE = 168

# ---- ctypes bindings ------------------------------------------------------

def _default_lib():
    ext = "dylib" if sys.platform == "darwin" else "so"
    for build_dir in ("build-ai", "build"):
        path = os.path.join(REPO, build_dir, "lib", f"libsecp256k1.{ext}")
        if os.path.exists(path):
            return path
    sys.exit("libsecp256k1 shared library not found; build with\n"
             "  cmake -B build -DSECP256K1_EXPERIMENTAL=ON -DSECP256K1_ENABLE_MODULE_SILENTPAYMENTS=ON\n"
             "  cmake --build build\n"
             "or point SECP256K1_LIB at it.")

LIB = os.environ.get("SECP256K1_LIB") or _default_lib()
C = ctypes.CDLL(LIB)

Buf64 = ctypes.c_ubyte * 64
Buf68 = ctypes.c_ubyte * 68
Buf96 = ctypes.c_ubyte * 96
Buf101 = ctypes.c_ubyte * 101

class Recipient(ctypes.Structure):
    _fields_ = [("scan_pubkey", Buf64), ("spend_pubkey", Buf64), ("index", ctypes.c_size_t)]

class FoundOutput(ctypes.Structure):
    _fields_ = [("output", Buf64), ("tweak", ctypes.c_ubyte * 32),
                ("found_with_label", ctypes.c_int), ("label", Buf68)]

LOOKUP_FN = ctypes.CFUNCTYPE(ctypes.c_void_p, ctypes.POINTER(ctypes.c_ubyte), ctypes.c_void_p)

C.secp256k1_context_create.restype = ctypes.c_void_p
CTX = C.secp256k1_context_create(1)  # SECP256K1_CONTEXT_NONE
assert CTX

def c_seckey_verify(sk):
    return C.secp256k1_ec_seckey_verify(ctypes.c_void_p(CTX), bytes(sk)) == 1

def c_pubkey_create(sk):
    pk = Buf64()
    assert C.secp256k1_ec_pubkey_create(ctypes.c_void_p(CTX), pk, bytes(sk)) == 1
    return pk

def c_xonly_serialize(xonly64):
    out = bytes(32)
    assert C.secp256k1_xonly_pubkey_serialize(ctypes.c_void_p(CTX), out, xonly64) == 1
    return out

def ptr_array(objs, ty):
    return (ctypes.POINTER(ty) * len(objs))(*[ctypes.pointer(o) for o in objs])

# ---- one differential run -------------------------------------------------

def run_one(data):
    if len(data) < MIN_SIZE:
        return "short"
    cfg = data[0]
    use_labels = cfg & 1
    lookup_null = (cfg >> 1) & 1
    hostile = (not lookup_null) and ((cfg >> 2) & 1)
    taproot_input = (cfg >> 3) & 1
    n_send = 1 + ((cfg >> 4) & 3)
    n_decoys = data[1] % 24
    label_m = data[2]
    scan_sk, spend_sk, input_sk = data[3:35], data[35:67], data[67:99]
    outpoint = bytes(data[131:167])
    if hostile:
        return "hostile-skipped"
    for sk in (scan_sk, spend_sk, input_sk):
        if not c_seckey_verify(sk):
            return "invalid-seckey"
    decoy_sks = []
    for i in range(n_decoys):
        dsk = bytes((data[167] ^ ((i * 17 + j) & 0xFF)) for j in range(32))
        if not c_seckey_verify(dsk):
            return "invalid-decoy"
        decoy_sks.append(dsk)

    # ---- C side ----
    scan_pk, spend_pk = c_pubkey_create(scan_sk), c_pubkey_create(spend_sk)
    label = Buf68()
    label_tweak = bytes(32)
    assert C.secp256k1_silentpayments_recipient_label_create(
        ctypes.c_void_p(CTX), label, label_tweak, bytes(scan_sk), ctypes.c_uint32(label_m)) == 1
    label33 = bytes(33)
    assert C.secp256k1_silentpayments_recipient_label_serialize(ctypes.c_void_p(CTX), label33, label) == 1
    if use_labels:
        addr_spend_pk = Buf64()
        assert C.secp256k1_silentpayments_recipient_create_labeled_spend_pubkey(
            ctypes.c_void_p(CTX), addr_spend_pk, spend_pk, label) == 1
    else:
        addr_spend_pk = spend_pk

    recipients = [Recipient(Buf64(*scan_pk), Buf64(*addr_spend_pk), i) for i in range(n_send)]
    gen_outputs = [Buf64() for _ in range(n_send)]
    if taproot_input:
        kp = Buf96()
        assert C.secp256k1_keypair_create(ctypes.c_void_p(CTX), kp, bytes(input_sk)) == 1
        ret = C.secp256k1_silentpayments_sender_create_outputs(
            ctypes.c_void_p(CTX), ptr_array(gen_outputs, Buf64), ptr_array(recipients, Recipient),
            ctypes.c_size_t(n_send), outpoint, ptr_array([kp], Buf96), ctypes.c_size_t(1), None, ctypes.c_size_t(0))
    else:
        sk_buf = ctypes.create_string_buffer(bytes(input_sk), 32)
        sk_ptrs = (ctypes.c_char_p * 1)(ctypes.cast(sk_buf, ctypes.c_char_p))
        ret = C.secp256k1_silentpayments_sender_create_outputs(
            ctypes.c_void_p(CTX), ptr_array(gen_outputs, Buf64), ptr_array(recipients, Recipient),
            ctypes.c_size_t(n_send), outpoint, None, ctypes.c_size_t(0), sk_ptrs, ctypes.c_size_t(1))
    assert ret == 1
    c_gen = [c_xonly_serialize(o) for o in gen_outputs]

    # ---- reference sender ----
    scan_ge = Scalar.from_bytes_checked(bytes(scan_sk)) * G
    spend_ge = Scalar.from_bytes_checked(bytes(spend_sk)) * G
    ref_label_ge, ref_label_tweak = bip352.silentpayments_recipient_label_create(bytes(scan_sk), label_m)
    assert ref_label_ge.to_bytes_compressed() == label33, "label point mismatch"
    assert ref_label_tweak.to_bytes() == label_tweak, "label tweak mismatch"
    addr_spend_ge = bip352.silentpayments_recipient_create_labeled_spend_pubkey(spend_ge, ref_label_ge) \
        if use_labels else spend_ge
    ref_recipients = [bip352.silentpayments_recipient(scan_ge, addr_spend_ge, i) for i in range(n_send)]
    ref_gen = bip352.silentpayments_sender_create_outputs(
        ref_recipients, outpoint,
        [bytes(input_sk)] if taproot_input else [],
        [] if taproot_input else [bytes(input_sk)])
    assert c_gen == ref_gen, f"sender outputs differ:\nC:   {[o.hex() for o in c_gen]}\nref: {[o.hex() for o in ref_gen]}"

    # ---- transaction outputs (decoys first, real last; same as fuzz target) ----
    tx_outputs32 = [(Scalar.from_bytes_checked(d) * G).to_bytes_xonly() for d in decoy_sks] + c_gen

    # ---- C scan ----
    tx_xonly = []
    for o in tx_outputs32:
        x = Buf64()
        assert C.secp256k1_xonly_pubkey_parse(ctypes.c_void_p(CTX), x, o) == 1
        tx_xonly.append(x)
    ps = Buf101()
    if taproot_input:
        input_xonly = Buf64()
        assert C.secp256k1_keypair_xonly_pub(ctypes.c_void_p(CTX), input_xonly, None, kp) == 1
        assert C.secp256k1_silentpayments_recipient_prevouts_summary_create(
            ctypes.c_void_p(CTX), ps, outpoint, ptr_array([input_xonly], Buf64), ctypes.c_size_t(1),
            None, ctypes.c_size_t(0)) == 1
        input_ges = ([GE.from_bytes_xonly(c_xonly_serialize(input_xonly))], [])
    else:
        input_pk = c_pubkey_create(input_sk)
        assert C.secp256k1_silentpayments_recipient_prevouts_summary_create(
            ctypes.c_void_p(CTX), ps, outpoint, None, ctypes.c_size_t(0),
            ptr_array([input_pk], Buf64), ctypes.c_size_t(1)) == 1
        input_ges = ([], [Scalar.from_bytes_checked(bytes(input_sk)) * G])

    tweak_keepalive = ctypes.create_string_buffer(label_tweak, 32)
    def lookup(label_ptr, _lctx):
        if bytes(label_ptr[0:33]) == label33:
            return ctypes.cast(tweak_keepalive, ctypes.c_void_p).value
        return None
    lookup_cb = LOOKUP_FN(lookup)

    found = [FoundOutput() for _ in range(len(tx_xonly))]
    n_found = ctypes.c_uint32(0)
    assert C.secp256k1_silentpayments_recipient_scan_outputs(
        ctypes.c_void_p(CTX), ptr_array(found, FoundOutput), ctypes.byref(n_found),
        ptr_array(tx_xonly, Buf64), ctypes.c_size_t(len(tx_xonly)), bytes(scan_sk), ps, spend_pk,
        None if lookup_null else lookup_cb, None) == 1

    # ---- reference scan ----
    ref_ps = bip352.silentpayments_recipient_prevouts_summary_create(outpoint, *input_ges)
    ref_cache = None if lookup_null else {label33: label_tweak}
    ref_found = bip352.silentpayments_recipient_scan_outputs(
        list(tx_outputs32), bytes(scan_sk), ref_ps, spend_ge, ref_cache)

    # ---- compare ----
    assert n_found.value == len(ref_found), f"n_found: C={n_found.value} ref={len(ref_found)}"
    for i, rf in enumerate(ref_found):
        cf = found[i]
        c_out = c_xonly_serialize(cf.output)
        assert c_out == rf.output, f"found[{i}].output: C={c_out.hex()} ref={rf.output.hex()}"
        assert bytes(cf.tweak) == rf.tweak.to_bytes(), \
            f"found[{i}].tweak: C={bytes(cf.tweak).hex()} ref={rf.tweak.to_bytes().hex()}"
        assert bool(cf.found_with_label) == rf.found_with_label, f"found[{i}].found_with_label differs"
        if rf.found_with_label:
            c_label33 = bytes(33)
            assert C.secp256k1_silentpayments_recipient_label_serialize(
                ctypes.c_void_p(CTX), c_label33, cf.label) == 1
            assert c_label33 == rf.label.to_bytes_compressed(), f"found[{i}].label differs"
    return "ok"

# ---- driver ----------------------------------------------------------------

def main():
    files = sys.argv[1:]
    if len(files) == 1 and os.path.isdir(files[0]):
        d = files[0]
        files = [os.path.join(d, f) for f in sorted(os.listdir(d))]
    if not files:
        print(__doc__)
        return 1
    stats = {}
    for path in files:
        with open(path, "rb") as f:
            data = f.read()
        try:
            result = run_one(data)
        except AssertionError as e:
            print(f"MISMATCH in {path}:\n{e}")
            return 1
        stats[result] = stats.get(result, 0) + 1
    print(f"differential OK over {len(files)} input(s): {stats}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
