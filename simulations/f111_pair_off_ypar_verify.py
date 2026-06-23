"""
F111 verify: classify all pairs in the off-y_par(D) sector of the diagonal
Klein cell at k=4 N=4. Expected: ALL SOFT (no hard).

For each dephase D, the off-y_par(D) sector has 32 templates. Pairs (with
self): 32 * 33 / 2 = 528 pairs.

Also classify pairs in the on-y_par(D) sector to confirm hard count = 228.

Also test the tensor-sum lemma: for a pair (P, Q) in off-y_par(D), build
H = chain(P) + chain(Q) and decompose by joint non-I support.
"""

import os
from datetime import datetime
from itertools import product as iprod, combinations_with_replacement
from collections import Counter, defaultdict

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
os.makedirs(RESULTS_DIR, exist_ok=True)
OUT_LOG = os.path.join(RESULTS_DIR, "f111_pair_off_ypar_verify.txt")
f_log = open(OUT_LOG, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


def bit_a(letter): return 1 if letter in "XY" else 0
def bit_b(letter): return 1 if letter in "YZ" else 0


def term_klein(letters):
    return (sum(bit_a(c) for c in letters) % 2,
            sum(bit_b(c) for c in letters) % 2)


def term_y_par(letters):
    return sum(1 for c in letters if c == "Y") % 2


def diagonal_klein(dephase):
    return (bit_a(dephase), bit_b(dephase))


def y_par_dephase(dephase):
    return bit_a(dephase) & bit_b(dephase)


def build_pauli_op_full(letters):
    mat = PAULI[letters[0]]
    for k in range(1, len(letters)):
        mat = np.kron(mat, PAULI[letters[k]])
    return mat


def build_chain_k_body(N, template):
    k = len(template)
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    if k > N:
        return H
    for shift in range(N - k + 1):
        full = ["I"] * N
        for i in range(k):
            full[shift + i] = template[i]
        H = H + build_pauli_op_full(full)
    return H


def build_L(H, gamma, N, dephase):
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    D_op = PAULI[dephase]
    for k in range(N):
        ops = [I2] * N
        ops[k] = D_op
        Dk = ops[0]
        for o in ops[1:]:
            Dk = np.kron(Dk, o)
        L += gamma * (np.kron(Dk, Dk.conj()) - np.eye(d2))
    return L


def canonical_pi_per_site(dephase):
    if dephase == "Z":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase == "X":
        return {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase == "Y":
        return {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    raise ValueError(dephase)


def build_Q_Nsite(per_site, N):
    d = 2 ** N
    d2 = d * d
    basis_N = list(iprod(LABELS, repeat=N))
    label_to_idx = {bl: i for i, bl in enumerate(basis_N)}
    mats_N = []
    for bl in basis_N:
        mat = PAULI[bl[0]]
        for k in range(1, N):
            mat = np.kron(mat, PAULI[bl[k]])
        mats_N.append(mat)
    vecs_N = [m.flatten() for m in mats_N]
    Q = np.zeros((d2, d2), dtype=complex)
    for idx, bl in enumerate(basis_N):
        phase = 1.0
        tgt_labels = []
        for sl in bl:
            ph, t = per_site[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


def spec_pairs(evals, sigma, tol):
    used = [False] * len(evals)
    max_err = 0.0
    for i in range(len(evals)):
        if used[i]:
            continue
        target = -evals[i] - 2 * sigma
        best_j, best_d = -1, np.inf
        for j in range(len(evals)):
            if used[j]:
                continue
            d = abs(evals[j] - target)
            if d < best_d:
                best_d = d
                best_j = j
        if best_j < 0:
            return False, np.inf
        used[i] = True
        if best_j != i:
            used[best_j] = True
        if best_d > max_err:
            max_err = best_d
    return max_err < tol, max_err


def classify_pair(N, templates, dephase, gamma=0.05,
                  op_tol=1e-10, spec_tol=1e-6):
    H = sum(build_chain_k_body(N, t) for t in templates)
    L = build_L(H, gamma, N, dephase)
    sigma = N * gamma
    Pi = build_Q_Nsite(canonical_pi_per_site(dephase), N)
    Pi_inv = np.linalg.inv(Pi)
    M = Pi @ L @ Pi_inv + L + 2 * sigma * np.eye(L.shape[0])
    op_norm = np.linalg.norm(M)
    if op_norm < op_tol:
        return "truly", op_norm, 0.0
    evals = np.linalg.eigvals(L)
    paired, max_err = spec_pairs(evals, sigma, spec_tol)
    cls = "soft" if paired else "hard"
    return cls, op_norm, max_err


def enumerate_k_terms_in_cell(k, klein, y_par):
    terms = []
    for s in iprod(LABELS, repeat=k):
        if all(c == "I" for c in s):
            continue
        seq = "".join(s)
        if term_klein(seq) == klein and term_y_par(seq) == y_par:
            terms.append(seq)
    return terms


def is_pure_D_template(template, dephase):
    """A template is 'pure-D' if it has only D and I letters (no other non-I)."""
    for c in template:
        if c != "I" and c != dephase:
            return False
    return True


def main():
    log("=" * 88)
    log("F111 pair verification: classify all pairs in OFF-y_par(D) at k=4 N=4")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}")
    log("Expected per F111 claim: off-y_par(D) sector has 0 hard pairs.")

    overall = {}
    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase} (diagonal Klein = {diagonal_klein(dephase)}, "
            f"y_par(D) = {y_par_dephase(dephase)})")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        on_y_par = y_par_dephase(dephase)

        off_templates = enumerate_k_terms_in_cell(k, diag, off_y_par)
        on_templates = enumerate_k_terms_in_cell(k, diag, on_y_par)

        # Classify on-y_par pairs
        log(f"\n-- ON-y_par(D) pairs (y_par={on_y_par}, {len(on_templates)} templates,"
            f" {len(on_templates) * (len(on_templates) + 1) // 2} pairs) --")
        on_counts = Counter()
        hard_on = []
        soft_on = []
        for i, t1 in enumerate(on_templates):
            for j in range(i, len(on_templates)):
                t2 = on_templates[j]
                cls, on, se = classify_pair(N, [t1, t2], dephase, gamma)
                on_counts[cls] += 1
                if cls == "hard":
                    hard_on.append((t1, t2, on, se))
                elif cls == "soft":
                    soft_on.append((t1, t2, on, se))
        log(f"  truly={on_counts['truly']}, soft={on_counts['soft']}, "
            f"hard={on_counts['hard']}")
        log(f"  EXPECTED per F106: hard={228}")

        # Classify off-y_par pairs
        log(f"\n-- OFF-y_par(D) pairs (y_par={off_y_par}, {len(off_templates)} templates,"
            f" {len(off_templates) * (len(off_templates) + 1) // 2} pairs) --")
        off_counts = Counter()
        hard_off = []
        for i, t1 in enumerate(off_templates):
            for j in range(i, len(off_templates)):
                t2 = off_templates[j]
                cls, on, se = classify_pair(N, [t1, t2], dephase, gamma)
                off_counts[cls] += 1
                if cls == "hard":
                    hard_off.append((t1, t2, on, se))
        log(f"  truly={off_counts['truly']}, soft={off_counts['soft']}, "
            f"hard={off_counts['hard']}")
        log(f"  EXPECTED per F111: hard=0")
        log(f"  F111 verification for dephase {dephase}: "
            f"{'PASS' if off_counts['hard'] == 0 else 'FAIL'}")

        # Analyze hard ON pairs: pure-D vs mixed?
        if hard_on:
            log(f"\n-- Hard ON pair analysis (first 10) --")
            pure_pure = 0
            mixed_pure = 0
            pure_mixed = 0
            mixed_mixed = 0
            for t1, t2, on, se in hard_on:
                p1 = is_pure_D_template(t1, dephase)
                p2 = is_pure_D_template(t2, dephase)
                if p1 and p2:
                    pure_pure += 1
                elif p1 and not p2:
                    pure_mixed += 1
                elif not p1 and p2:
                    mixed_pure += 1
                else:
                    mixed_mixed += 1
            log(f"  Pure-Pure pairs: {pure_pure}")
            log(f"  Pure-Mixed pairs (sum, both orderings): {pure_mixed + mixed_pure}")
            log(f"  Mixed-Mixed pairs: {mixed_mixed}")
            log(f"  Total: {pure_pure + pure_mixed + mixed_pure + mixed_mixed}")

        overall[dephase] = (on_counts, off_counts)

    log("\n" + "=" * 88)
    log("Summary across dephase letters:")
    log("=" * 88)
    log(f"{'D':<3} | {'On hard':>8} | {'On total':>9} | {'Off hard':>9} | {'Off total':>10}")
    for dephase in ["Z", "X", "Y"]:
        on_c, off_c = overall[dephase]
        on_tot = sum(on_c.values())
        off_tot = sum(off_c.values())
        log(f"{dephase:<3} | {on_c['hard']:>8d} | {on_tot:>9d} | "
            f"{off_c['hard']:>9d} | {off_tot:>10d}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
