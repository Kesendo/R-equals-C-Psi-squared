"""
F111 empirical exploration: extract the actual F87-hard pairs at k=4 N=4
in each diagonal Klein cell per dephase, and characterize their letter-count
profiles by y_par.

Statement under test (F111):
  At k=4 N=4, every F87-hard Pauli pair (P, Q) in the diagonal Klein cell
  for dephase letter D in {Z, X, Y} has y_par(pair) = y_par(D) = D.BitA() AND D.BitB().

Recall:
  D = Z: diagonal Klein (0, 1), y_par(D) = 0
  D = X: diagonal Klein (1, 0), y_par(D) = 0
  D = Y: diagonal Klein (1, 1), y_par(D) = 1

Empirical (F106): 228:0 / 228:0 / 0:228 splits across the three diagonals.

This script:
  1. Enumerates all 4248 k=4 Klein/y_par-homogeneous pairs.
  2. For each dephase, classifies each pair via the explicit Liouvillian
     palindrome residual + spectrum pairing (matches PauliPairTrichotomy).
  3. For HARD pairs in the diagonal cell, dumps letter-count profiles
     binned by y_par. Also dumps complete pair tuples to JSON.
  4. Counts soft/truly pairs in the off-y_par(D) sector of the diagonal cell
     (the sector that must be HARD-FREE per F111).
"""

import json
import os
from collections import Counter, defaultdict
from itertools import product as iprod
from datetime import datetime

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
OUT_JSON = os.path.join(RESULTS_DIR, "f111_hard_pairs_extraction.json")
OUT_LOG = os.path.join(RESULTS_DIR, "f111_hard_pairs_extraction.txt")
os.makedirs(RESULTS_DIR, exist_ok=True)
f_log = open(OUT_LOG, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


# ---- Pauli matrices ----
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


def bit_a(letter):
    return 1 if letter in "XY" else 0


def bit_b(letter):
    return 1 if letter in "YZ" else 0


def term_klein(letters):
    return (sum(bit_a(c) for c in letters) % 2,
            sum(bit_b(c) for c in letters) % 2)


def term_y_par(letters):
    return sum(1 for c in letters if c == "Y") % 2


def term_letter_counts(letters):
    return (sum(1 for c in letters if c == "I"),
            sum(1 for c in letters if c == "X"),
            sum(1 for c in letters if c == "Y"),
            sum(1 for c in letters if c == "Z"))


def diagonal_klein(dephase):
    return (bit_a(dephase), bit_b(dephase))


def y_par_dephase(dephase):
    """y_par(D) = (#Y in D) mod 2 = bit_a(D) AND bit_b(D)."""
    return bit_a(dephase) & bit_b(dephase)


# ---- Build N-qubit Pauli operator from letter sequence ----
def build_pauli_op(letters):
    mat = PAULI[letters[0]]
    for k in range(1, len(letters)):
        mat = np.kron(mat, PAULI[letters[k]])
    return mat


# ---- Chain k-body Hamiltonian: place k-letter template at each shift on chain ----
def build_chain_k_body(N, template):
    """Place k-letter template at sites (0..k-1), (1..k), ..., (N-k..N-1).
    Returns the d x d Hamiltonian matrix (no coefficient; coefficient is implicit 1.0
    matching PauliPairTrichotomy.ChainKBody convention)."""
    k = len(template)
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    if k > N:
        return H
    for shift in range(N - k + 1):
        # template at sites [shift, shift+1, ..., shift+k-1], identities elsewhere
        full = ["I"] * N
        for i in range(k):
            full[shift + i] = template[i]
        H = H + build_pauli_op(full)
    return H


def build_chain_pair_H(N, term1, term2):
    """H for a Pauli pair: chain k-body sum of term1 PLUS chain k-body sum of term2."""
    return build_chain_k_body(N, term1) + build_chain_k_body(N, term2)


# ---- Build L = -i [H, .] + sum_l gamma * D[D_l] in vec basis ----
def build_L(H, gamma, N, dephase_letter):
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    D_op = PAULI[dephase_letter]
    for k in range(N):
        ops = [I2] * N
        ops[k] = D_op
        Dk = ops[0]
        for o in ops[1:]:
            Dk = np.kron(Dk, o)
        L += gamma * (np.kron(Dk, Dk.conj()) - np.eye(d2))
    return L


# ---- Build palindrome residual M = Pi L Pi^-1 + L + 2 sigma I ----
def build_pi_canonical(N, dephase_letter):
    """Canonical PiOperator (NOT Pi_5bilinear). Per dephase letter, per-site action:
       Z: I<->X (sign +1), Y<->Z (sign +i)
       X: I<->Z (sign +1), X<->Y (sign -i)
       Y: I<->X (sign +1), Y<->Z (sign -i)
    """
    if dephase_letter == "Z":
        per_site = {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}
    elif dephase_letter == "X":
        per_site = {"I": (1, "Z"), "Z": (1, "I"), "X": (-1j, "Y"), "Y": (-1j, "X")}
    elif dephase_letter == "Y":
        per_site = {"I": (1, "X"), "X": (1, "I"), "Y": (-1j, "Z"), "Z": (-1j, "Y")}
    else:
        raise ValueError(dephase_letter)
    return build_Q_Nsite(per_site, N)


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


# ---- Spectrum pairing greedy match (matches PauliPairTrichotomy.SpectrumPairs) ----
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


def classify_pair(N, term1, term2, dephase_letter,
                  op_tol=1e-10, spec_tol=1e-6, gamma=0.05):
    H = build_chain_pair_H(N, term1, term2)
    L = build_L(H, gamma, N, dephase_letter)
    sigma = N * gamma
    Pi = build_pi_canonical(N, dephase_letter)
    Pi_inv = np.linalg.inv(Pi)
    M = Pi @ L @ Pi_inv + L + 2 * sigma * np.eye(L.shape[0])
    op_norm = np.linalg.norm(M)
    if op_norm < op_tol:
        return "truly", op_norm
    evals = np.linalg.eigvals(L)
    paired, err = spec_pairs(evals, sigma, spec_tol)
    if paired:
        return "soft", op_norm
    return "hard", op_norm


# ---- Enumerate all k-letter Pauli sequences ----
def enumerate_k_terms(k):
    """All non-trivial k-letter strings (excluding all-I)."""
    return ["".join(s) for s in iprod(LABELS, repeat=k)
            if any(c != "I" for c in s)]


def enumerate_homogeneous_pairs(k):
    """Unordered Klein-homogeneous + y_par-homogeneous pairs (with self-pairs)."""
    terms = enumerate_k_terms(k)
    # Index terms by (klein, y_par)
    by_grp = defaultdict(list)
    for t in terms:
        by_grp[(term_klein(t), term_y_par(t))].append(t)
    pairs = []
    for (klein, y_par), grp in by_grp.items():
        for i, t1 in enumerate(grp):
            for j in range(i, len(grp)):
                t2 = grp[j]
                pairs.append((t1, t2, klein, y_par))
    return pairs


def main():
    log("=" * 88)
    log("F111 empirical exploration: extract F87-hard pairs at k=4 N=4")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05

    pairs = enumerate_homogeneous_pairs(k)
    log(f"\nTotal k={k} Klein/y_par-homogeneous pairs: {len(pairs)}")

    # Classify all pairs under each dephase, only collect HARD pairs in diagonal cell.
    hard_records = defaultdict(list)  # dephase -> list of (term1, term2, y_par, op_norm)
    sector_counts = defaultdict(lambda: defaultdict(Counter))  # dephase -> y_par -> class -> count

    for dephase in ["Z", "X", "Y"]:
        log(f"\n--- Dephase {dephase} (diagonal cell = {diagonal_klein(dephase)},"
            f" y_par(D) = {y_par_dephase(dephase)}) ---")
        diag = diagonal_klein(dephase)
        diag_pairs = [p for p in pairs if p[2] == diag]
        log(f"  Klein-diagonal pairs at k={k}: {len(diag_pairs)}")
        # Classify
        # Bucket by y_par for clean reporting
        per_ypar = defaultdict(list)
        for (t1, t2, klein, y_par) in diag_pairs:
            per_ypar[y_par].append((t1, t2))
        for y_par_val in [0, 1]:
            log(f"  y_par={y_par_val}: {len(per_ypar[y_par_val])} pairs")

        # Run classification (this is the slow part: ~3000 pairs, each ~0.5s)
        for y_par_val in [0, 1]:
            for (t1, t2) in per_ypar[y_par_val]:
                cls, op_norm = classify_pair(N, t1, t2, dephase, gamma=gamma)
                sector_counts[dephase][y_par_val][cls] += 1
                if cls == "hard":
                    hard_records[dephase].append({
                        "term1": t1, "term2": t2,
                        "y_par": y_par_val,
                        "op_norm_M": float(op_norm),
                        "lc_t1": list(term_letter_counts(t1)),
                        "lc_t2": list(term_letter_counts(t2)),
                    })

        log(f"\n  Per-y_par classification counts:")
        for y_par_val in [0, 1]:
            cnt = sector_counts[dephase][y_par_val]
            log(f"    y_par={y_par_val}: truly={cnt['truly']:>5d}, "
                f"soft={cnt['soft']:>5d}, hard={cnt['hard']:>5d}")

    # ============================================================
    # Letter-count profile dump for HARD pairs per dephase
    # ============================================================
    log("\n" + "=" * 88)
    log("HARD PAIR LETTER-COUNT PROFILES (binned by y_par)")
    log("=" * 88)
    for dephase in ["Z", "X", "Y"]:
        records = hard_records[dephase]
        log(f"\nDephase {dephase}: {len(records)} hard pairs total")
        if not records:
            log("  (no hard pairs)")
            continue
        by_ypar = defaultdict(list)
        for r in records:
            by_ypar[r["y_par"]].append(r)
        for y_par_val in [0, 1]:
            sub = by_ypar[y_par_val]
            if not sub:
                continue
            log(f"  y_par={y_par_val} ({len(sub)} hard pairs):")
            # Aggregate letter-count profiles
            profile_counter = Counter()
            for r in sub:
                # Profile = sorted (lc1, lc2) (#I, #X, #Y, #Z)
                key = (tuple(r["lc_t1"]), tuple(r["lc_t2"]))
                profile_counter[key] += 1
            for prof, n in sorted(profile_counter.items(), key=lambda kv: -kv[1]):
                lc1, lc2 = prof
                log(f"    {n:>4d} pairs: t1=(I={lc1[0]},X={lc1[1]},Y={lc1[2]},Z={lc1[3]}) "
                    f"+ t2=(I={lc2[0]},X={lc2[1]},Y={lc2[2]},Z={lc2[3]})")

    # ============================================================
    # Sample list of off-y_par(D) candidates (the sector that is empirically HARD-FREE)
    # ============================================================
    log("\n" + "=" * 88)
    log("OFF-y_par(D) CANDIDATE PAIRS (the empirically HARD-FREE sector)")
    log("=" * 88)
    for dephase in ["Z", "X", "Y"]:
        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        diag_pairs = [(t1, t2) for (t1, t2, klein, y_par) in pairs
                      if klein == diag and y_par == off_y_par]
        log(f"\nDephase {dephase}, off-y_par sector (Klein={diag}, y_par={off_y_par}):")
        log(f"  {len(diag_pairs)} candidate pairs, classification:")
        cnt = sector_counts[dephase][off_y_par]
        log(f"    truly={cnt['truly']}, soft={cnt['soft']}, hard={cnt['hard']}")
        log(f"  Sample first 5 pairs:")
        for t1, t2 in diag_pairs[:5]:
            log(f"    ({t1}, {t2})")

    # ============================================================
    # Dump JSON
    # ============================================================
    out = {
        "N": N,
        "k": k,
        "gamma": gamma,
        "sector_counts": {
            dephase: {
                str(y_par): dict(sector_counts[dephase][y_par])
                for y_par in [0, 1]
            }
            for dephase in ["Z", "X", "Y"]
        },
        "hard_records": {dephase: hard_records[dephase] for dephase in ["Z", "X", "Y"]},
    }
    with open(OUT_JSON, "w") as fp:
        json.dump(out, fp, indent=2)
    log(f"\nWrote: {OUT_JSON}")
    log(f"Wrote log: {OUT_LOG}")


if __name__ == "__main__":
    main()
