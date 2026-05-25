"""
F111 structural decomposition: key insight is that at k=N=4, the chain k-body H
is a single Pauli string at fixed positions. If the template has internal I's,
those sites are free (no H, only dephasing). The Liouvillian factorizes.

We verify:
1. At k=N=4, the effective H of an off-y_par(D) template = a Pauli-string
   at a subset of sites (the non-I positions), surrounded by free dephasing
   sites.
2. The spec of L = L_free (on I-sites) + L_interacting (on non-I-sites).
3. For off-y_par templates, the interacting H is in the diagonal Klein cell
   of a SMALLER N (the non-I-site count) and is SOFT under D-dephasing on
   those sites.

If we can show: for every off-y_par(D) k=4 template, the SMALLER-N interacting
Hamiltonian is SOFT (or truly), then the full L spectrum is palindromic by
the tensor-sum lemma.

We then need to verify:
  (a) Tensor-sum spec palindrome lemma: if L1 has palindromic spec around
      -sigma1 and L2 has palindromic spec around -sigma2, does L = L1 (x) I2 +
      I1 (x) L2 have palindromic spec around -(sigma1+sigma2)?

  (b) That off-y_par(D) templates' reduced (non-I-position) Hamiltonians
      indeed have palindromic spec.

Strategy:
1. Enumerate all off-y_par(D) k=4 templates; classify each as soft/hard
   (already done in _f111_spec_palindrome_single_term.py).
2. For each, extract its non-I positions and the reduced template.
3. Check whether the reduced template at the reduced N matches an F108/F110
   known soft sub-case.

This may give us the closed-form: every off-y_par(D) k=4 template reduces
to a smaller-N sub-Hamiltonian that is provably soft via F108 + sub-N analysis.
"""

import os
from datetime import datetime
from itertools import product as iprod
from collections import Counter, defaultdict

import numpy as np

RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
os.makedirs(RESULTS_DIR, exist_ok=True)
OUT_LOG = os.path.join(RESULTS_DIR, "f111_structural_decomposition.txt")
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


def reduce_template(template):
    """Extract the non-I sub-template and the I-positions. Returns
    (reduced_template, non_I_positions, num_I_sites)."""
    non_I = [i for i, c in enumerate(template) if c != "I"]
    reduced = "".join(template[i] for i in non_I)
    return reduced, non_I, len(template) - len(non_I)


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


def build_L_reduced(reduced_template, gamma, dephase, full_N):
    """Build L on the REDUCED system: only the non-I sites of the template.
    H is the reduced Pauli string at sites 0..k_body-1; dephasing is on all
    full_N sites (free dephasing on the implicit I sites contributes
    -gamma*num_I_sites * I to the spec)."""
    k_body = len(reduced_template)
    if k_body == 0:
        return None, 0
    # Reduced H acts on k_body qubits
    d_red = 2 ** k_body
    H_red = build_pauli_op_full(list(reduced_template))
    L_red = -1j * (np.kron(H_red, np.eye(d_red)) - np.kron(np.eye(d_red), H_red.T))
    # Dephasing on the k_body sites of the reduced system
    D_op = PAULI[dephase]
    for s in range(k_body):
        ops = [I2] * k_body
        ops[s] = D_op
        Ds = ops[0]
        for o in ops[1:]:
            Ds = np.kron(Ds, o)
        L_red += gamma * (np.kron(Ds, Ds.conj()) - np.eye(d_red * d_red))
    return L_red, k_body


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


def enumerate_k_terms_in_cell(k, klein, y_par):
    terms = []
    for s in iprod(LABELS, repeat=k):
        if all(c == "I" for c in s):
            continue
        seq = "".join(s)
        if term_klein(seq) == klein and term_y_par(seq) == y_par:
            terms.append(seq)
    return terms


def main():
    log("=" * 88)
    log("F111 structural decomposition: off-y_par k=4 templates reduce to smaller N")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    N = 4
    k = 4
    gamma = 0.05
    sigma_full = N * gamma

    log(f"\nSetup: N={N}, k={k}, gamma={gamma}, sigma_full={sigma_full}")

    for dephase in ["Z", "X", "Y"]:
        log("\n" + "=" * 88)
        log(f"DEPHASE {dephase} (diagonal Klein = {diagonal_klein(dephase)}, "
            f"y_par(D) = {y_par_dephase(dephase)})")
        log("=" * 88)

        diag = diagonal_klein(dephase)
        off_y_par = 1 - y_par_dephase(dephase)
        off_templates = enumerate_k_terms_in_cell(k, diag, off_y_par)
        log(f"\nOff-y_par templates: {len(off_templates)}")

        # Group by reduced template
        by_reduced = defaultdict(list)
        for t in off_templates:
            reduced, non_I, num_I = reduce_template(t)
            by_reduced[reduced].append((t, non_I, num_I))

        log(f"\nUnique reduced templates: {len(by_reduced)}")
        for reduced, instances in sorted(by_reduced.items(), key=lambda kv: (len(kv[0]), kv[0])):
            log(f"  {reduced} (k_body={len(reduced)}): {len(instances)} placements")

        # For each unique reduced template, compute L on the REDUCED system
        # and check spec palindrome.
        log(f"\nReduced-system spec palindrome check:")
        all_reduced_soft = True
        for reduced in sorted(by_reduced.keys(), key=lambda r: (len(r), r)):
            k_body = len(reduced)
            sigma_red = k_body * gamma  # dephasing on k_body sites
            L_red, _ = build_L_reduced(reduced, gamma, dephase, N)
            evals = np.linalg.eigvals(L_red)
            paired, max_err = spec_pairs(evals, sigma_red, tol=1e-6)
            status = "PAL" if paired else "HARD"
            log(f"  {reduced} (k_body={k_body}, sigma_red={sigma_red}): "
                f"{status}, max_err={max_err:.3e}")
            if not paired:
                all_reduced_soft = False

        # Now verify: full L spec = sum of reduced L spec + offset from free sites
        log(f"\nTensor-sum lemma verification: full L spec palindrome via reduction")
        for reduced in sorted(by_reduced.keys(), key=lambda r: (len(r), r))[:3]:
            instances = by_reduced[reduced]
            for full_template, non_I, num_I in instances[:1]:
                log(f"  Full template {full_template} = reduced {reduced} at non-I "
                    f"positions {non_I}, {num_I} I-sites")
                H_full = build_chain_k_body(N, full_template)
                L_full = build_L(H_full, gamma, N, dephase)
                evals_full = np.linalg.eigvals(L_full)
                paired_full, max_err_full = spec_pairs(evals_full, sigma_full, tol=1e-6)
                status_full = "PAL" if paired_full else "HARD"
                log(f"    Full L spec: {status_full}, max_err={max_err_full:.3e}")

        log(f"\nAll off-y_par reduced templates spec palindromic? {all_reduced_soft}")

    log("\n" + "=" * 88)
    log(f"Done: {datetime.now()}")
    log(f"Log: {OUT_LOG}")
    log("=" * 88)


if __name__ == "__main__":
    main()
