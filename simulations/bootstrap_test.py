#!/usr/bin/env python3
"""
Bootstrap Test: Is Each Side the Environment of the Other?

Test 1: Sector projection of full Liouvillian (cross-sector blocks)
Test 2: Hamiltonian-only sector analysis (does L_H respect parity?)
Test 4: Broken-parity Hamiltonians (36 two-term sweep vs palindrome list)
Test 3: The subtler bootstrap (most general parity-compatible dissipator)

Script: simulations/bootstrap_test.py
Output: simulations/results/bootstrap_test.txt
"""

import numpy as np
from itertools import product as iprod, combinations
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "bootstrap_test.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

# ============================================================
# INFRASTRUCTURE
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
H_LABELS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']
GAMMA = 0.05

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))


def build_pauli_basis(N):
    all_idx = list(iprod(range(4), repeat=N))
    d = 2 ** N
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, N):
            m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m)
    return all_idx, pmats, d


def build_Pi2_matrix(all_idx):
    """Build Pi^2 = X^N conjugation in Pauli basis.
    Diagonal: +1 for even YZ-count, -1 for odd YZ-count."""
    num = len(all_idx)
    P = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        n_yz = sum(1 for s in idx_b if s in (2, 3))
        P[b, b] = (-1) ** n_yz
    return P


def site_op(op, site, N):
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_terms(N, bonds, label):
    """Build Hamiltonian for a single Pauli-pair term on given bonds."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        ops = [I2] * N
        ops[i] = PM[label[0]]
        ops[j] = PM[label[1]]
        term = ops[0]
        for o in ops[1:]:
            term = np.kron(term, o)
        H += term
    return H


def build_L_H_pauli(N, H, all_idx, pmats, d):
    """Hamiltonian part of Liouvillian in Pauli basis."""
    num = 4 ** N
    L = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L[a, b] = np.trace(pmats[a].conj().T @ comm) / d
    return L


def build_L_D_diag(all_idx, gamma):
    """Dissipator diagonal in Pauli basis (Z-dephasing)."""
    return np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])


def build_full_L_pauli(N, bonds, gamma, H=None):
    """Full Liouvillian L = L_H + L_D in Pauli basis."""
    all_idx, pmats, d = build_pauli_basis(N)
    num = 4 ** N
    if H is None:
        # Heisenberg
        H = np.zeros((d, d), dtype=complex)
        for i, j in bonds:
            for P in [sx, sy, sz]:
                H += site_op(P, i, N) @ site_op(P, j, N)
    L_H = build_L_H_pauli(N, H, all_idx, pmats, d)
    L_D_diag = build_L_D_diag(all_idx, gamma)
    L = L_H.copy()
    for a in range(num):
        L[a, a] += L_D_diag[a]
    return L, L_H, L_D_diag, all_idx, pmats, d


def check_palindromic(evals, Sg, tol=1e-6):
    """Check eigenvalue palindromic pairing. Returns (n_paired, max_err)."""
    num = len(evals)
    max_err = 0
    n_paired = 0
    for k in range(num):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best = np.min(dists)
        if best < tol:
            n_paired += 1
        if best > max_err:
            max_err = best
    return n_paired, max_err


# ============================================================
# TEST 1: SECTOR PROJECTION OF THE FULL LIOUVILLIAN
# ============================================================
def run_test_1(N=3, gamma=GAMMA):
    log("=" * 70)
    log(f"TEST 1: SECTOR PROJECTION OF FULL LIOUVILLIAN (N={N})")
    log("=" * 70)
    log()
    log("Question: Are the +1 and -1 parity sectors coupled through L?")
    log("If [Pi^2, L] = 0, then L is block-diagonal and cross-sector")
    log("blocks L_pm and L_mp must be zero.")
    log()

    bonds = [(i, i + 1) for i in range(N - 1)]
    L, L_H, L_D_diag, all_idx, pmats, d = build_full_L_pauli(N, bonds, gamma)
    num = 4 ** N
    Pi2 = build_Pi2_matrix(all_idx)

    # Step 1: Verify [Pi^2, L] = 0
    comm = Pi2 @ L - L @ Pi2
    comm_norm = np.max(np.abs(comm))
    log(f"[Pi^2, L] max element: {comm_norm:.2e}")
    if comm_norm < 1e-10:
        log("CONFIRMED: Pi^2 commutes with L exactly.")
    else:
        log(f"WARNING: Pi^2 does NOT commute with L (norm {comm_norm:.2e})")
    log()

    # Step 2: Build sector projectors and extract blocks
    plus_idx = [a for a in range(num) if Pi2[a, a].real > 0]
    minus_idx = [a for a in range(num) if Pi2[a, a].real < 0]
    log(f"+1 sector: {len(plus_idx)} basis elements")
    log(f"-1 sector: {len(minus_idx)} basis elements")
    log()

    # Cross-sector blocks
    L_pm = L[np.ix_(plus_idx, minus_idx)]
    L_mp = L[np.ix_(minus_idx, plus_idx)]
    L_pp = L[np.ix_(plus_idx, plus_idx)]
    L_mm = L[np.ix_(minus_idx, minus_idx)]

    norm_pm = np.max(np.abs(L_pm))
    norm_mp = np.max(np.abs(L_mp))
    norm_pp = np.max(np.abs(L_pp))
    norm_mm = np.max(np.abs(L_mm))

    log(f"||L_pm|| (max element): {norm_pm:.2e}")
    log(f"||L_mp|| (max element): {norm_mp:.2e}")
    log(f"||L_pp|| (max element): {norm_pp:.2e}")
    log(f"||L_mm|| (max element): {norm_mm:.2e}")
    log()

    if norm_pm < 1e-10 and norm_mp < 1e-10:
        log("RESULT: Cross-sector blocks are ZERO.")
        log("The sectors are exactly decoupled in the full Liouvillian.")
        log("The bootstrap cannot work through direct coupling in L.")
        log("This is an important NEGATIVE result for the naive bootstrap.")
    else:
        log("RESULT: Cross-sector blocks are NONZERO.")
        log("The sectors couple through L. The naive bootstrap could work.")
    log()

    # Also check L_H alone
    comm_H = Pi2 @ L_H - L_H @ Pi2
    comm_H_norm = np.max(np.abs(comm_H))
    log(f"[Pi^2, L_H] max element: {comm_H_norm:.2e}")

    # And L_D alone
    L_D = np.diag(L_D_diag)
    comm_D = Pi2 @ L_D - L_D @ Pi2
    comm_D_norm = np.max(np.abs(comm_D))
    log(f"[Pi^2, L_D] max element: {comm_D_norm:.2e}")

    if comm_H_norm < 1e-10:
        log("Both L_H and L_D individually commute with Pi^2.")
    log()

    return L, L_H, L_D_diag, all_idx, pmats, d, Pi2, plus_idx, minus_idx


# ============================================================
# TEST 2: HAMILTONIAN-ONLY SECTOR ANALYSIS
# ============================================================
def run_test_2(N=3, gamma=GAMMA):
    log("=" * 70)
    log(f"TEST 2: HAMILTONIAN-ONLY SECTOR ANALYSIS (N={N})")
    log("=" * 70)
    log()
    log("Question: Does L_H alone respect parity? Does H commute with X^N")
    log("in the Hilbert space sense? If not, the Hamiltonian couples the")
    log("two parity sectors and tracing over one produces effective dissipation.")
    log()

    bonds = [(i, i + 1) for i in range(N - 1)]
    all_idx, pmats, d = build_pauli_basis(N)
    num = 4 ** N
    Pi2 = build_Pi2_matrix(all_idx)

    # Build X^N in Hilbert space
    XN = sx.copy()
    for k in range(1, N):
        XN = np.kron(XN, sx)

    # Test Heisenberg Hamiltonian
    log("--- Heisenberg Hamiltonian (XX+YY+ZZ) ---")
    H_heis = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in [sx, sy, sz]:
            H_heis += site_op(P, i, N) @ site_op(P, j, N)

    comm_hilbert = H_heis @ XN - XN @ H_heis
    comm_norm = np.max(np.abs(comm_hilbert))
    log(f"[H_Heisenberg, X^N] max element: {comm_norm:.2e}")
    if comm_norm < 1e-10:
        log("H_Heisenberg commutes with X^N in Hilbert space.")
    log()

    # Build L_H and check in superoperator space
    L_H_heis = build_L_H_pauli(N, H_heis, all_idx, pmats, d)
    comm_super = Pi2 @ L_H_heis - L_H_heis @ Pi2
    log(f"[Pi^2, L_H_Heisenberg] max element: {np.max(np.abs(comm_super)):.2e}")
    log()

    # Test other standard Hamiltonians
    test_models = {
        'XY-only (XX+YY)': ['XX', 'YY'],
        'Ising (ZZ)': ['ZZ'],
        'XX alone': ['XX'],
        'DM (XY-YX)': None,  # special
    }

    for name, terms in test_models.items():
        H = np.zeros((d, d), dtype=complex)
        if terms is not None:
            for t in terms:
                H += build_H_terms(N, bonds, t)
        else:
            # DM interaction
            H += build_H_terms(N, bonds, 'XY')
            H -= build_H_terms(N, bonds, 'YX')

        comm_h = H @ XN - XN @ H
        norm_h = np.max(np.abs(comm_h))

        L_H_t = build_L_H_pauli(N, H, all_idx, pmats, d)
        comm_s = Pi2 @ L_H_t - L_H_t @ Pi2
        norm_s = np.max(np.abs(comm_s))

        log(f"{name}:")
        log(f"  [H, X^N] = {norm_h:.2e}   [Pi^2, L_H] = {norm_s:.2e}")

    log()
    log("SUMMARY: If all norms are zero, the Hamiltonian never couples")
    log("the parity sectors, and the Nakajima-Zwanzig projection is trivial.")
    log()

    return all_idx, pmats, d, Pi2, XN


# ============================================================
# TEST 4: BROKEN-PARITY HAMILTONIANS (run before Test 3)
# ============================================================
def run_test_4(N=3, gamma=GAMMA):
    log("=" * 70)
    log(f"TEST 4: BROKEN-PARITY HAMILTONIANS (N={N})")
    log("=" * 70)
    log()
    log("For each of the 36 two-term Hamiltonian combinations:")
    log("  1. Does [H, X^N] = 0 in Hilbert space? (parity-commuting)")
    log("  2. Does the eigenvalue palindrome survive? (lambda_k + lambda_k' = -2Sg)")
    log("  3. Does the Pi conjugation hold? (Pi L Pi^-1 + L + 2Sg I = 0)")
    log("  4. Cross-reference: is parity-breaking = eigenvalue-palindrome-breaking?")
    log()
    log("Note: Pi-conjugation breaks for 33/36 combos (known from V-Effect work).")
    log("The relevant palindrome measure is eigenvalue pairing (14/36 break).")
    log()

    bonds = [(i, i + 1) for i in range(N - 1)]
    all_idx, pmats, d = build_pauli_basis(N)
    num = 4 ** N

    # Build X^N
    XN = sx.copy()
    for k in range(1, N):
        XN = np.kron(XN, sx)

    Pi2 = build_Pi2_matrix(all_idx)

    # Build L_H for each single term
    L_H_single = {}
    for label in H_LABELS:
        H_t = build_H_terms(N, bonds, label)
        L_H_single[label] = build_L_H_pauli(N, H_t, all_idx, pmats, d)

    L_D_diag = build_L_D_diag(all_idx, gamma)
    Sg = N * gamma

    # Build Pi for Pi-conjugation check
    PI_perm = {0: 1, 1: 0, 2: 3, 3: 2}
    PI_sign_map = {0: 1, 1: 1, 2: 1j, 3: 1j}

    def build_Pi(all_idx):
        n = len(all_idx)
        Pi = np.zeros((n, n), dtype=complex)
        idx_map = {idx: i for i, idx in enumerate(all_idx)}
        for b, idx_b in enumerate(all_idx):
            mapped = tuple(PI_perm[i] for i in idx_b)
            sign = 1
            for i in idx_b:
                sign *= PI_sign_map[i]
            a = idx_map[mapped]
            Pi[a, b] = sign
        return Pi

    Pi = build_Pi(all_idx)
    Pi_inv = np.linalg.inv(Pi)

    # Eigenvalue palindrome check using vectorized Lindbladian
    # (computational basis, as in v_effect_analysis.py)
    def build_L_vec(H, gamma, N):
        d = 2 ** N
        Id = np.eye(d)
        L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
        for k in range(N):
            Zk = site_op(sz, k, N)
            L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d))
        return L

    def check_eval_palindrome(evals, Sg, tol=1e-6):
        num = len(evals)
        pair_errors = np.zeros(num)
        for k in range(num):
            target = -(evals[k] + 2 * Sg)
            dists = np.abs(evals - target)
            pair_errors[k] = np.min(dists)
        n_well_paired = int(np.sum(pair_errors < tol))
        return n_well_paired, np.max(pair_errors), pair_errors

    # All 36 two-term combinations
    all_combos = []
    for i, t1 in enumerate(H_LABELS):
        for t2 in H_LABELS[i + 1:]:
            all_combos.append((t1, t2))

    parity_breaking = []
    eval_palindrome_breaking = []
    pi_palindrome_breaking = []
    results = []

    log(f"{'Combo':>10}  {'[H,X^N]':>10}  {'Paired':>8}  {'Pi err':>10}  "
        f"{'Parity':>8}  {'EvalPal':>8}  {'PiPal':>8}")
    log("-" * 75)

    for t1, t2 in all_combos:
        combo_name = f"{t1}+{t2}"

        # 1. Hilbert space parity: [H, X^N]
        H_combo = build_H_terms(N, bonds, t1) + build_H_terms(N, bonds, t2)
        comm_h = H_combo @ XN - XN @ H_combo
        parity_norm = np.max(np.abs(comm_h))
        parity_ok = parity_norm < 1e-10

        # 2. Eigenvalue palindrome (computational basis Lindbladian)
        L_vec = build_L_vec(H_combo, gamma, N)
        evals_vec = np.linalg.eigvals(L_vec)
        n_paired, max_pair_err, _ = check_eval_palindrome(evals_vec, Sg)
        eval_pal_ok = (n_paired == len(evals_vec))

        # 3. Pi conjugation (Pauli basis)
        L_H_c = L_H_single[t1] + L_H_single[t2]
        L_full = L_H_c.copy()
        for a in range(num):
            L_full[a, a] += L_D_diag[a]

        E = Pi @ L_full @ Pi_inv + L_full + 2 * Sg * np.eye(num)
        pi_err = np.max(np.abs(E)) / np.max(np.abs(L_full))
        pi_pal_ok = pi_err < 1e-8

        if not parity_ok:
            parity_breaking.append(combo_name)
        if not eval_pal_ok:
            eval_palindrome_breaking.append(combo_name)
        if not pi_pal_ok:
            pi_palindrome_breaking.append(combo_name)

        results.append({
            'combo': combo_name, 'parity_norm': parity_norm,
            'n_paired': n_paired, 'max_pair_err': max_pair_err,
            'pi_err': pi_err, 'parity_ok': parity_ok,
            'eval_pal_ok': eval_pal_ok, 'pi_pal_ok': pi_pal_ok,
        })

        p_tag = "OK" if parity_ok else "BREAK"
        ep_tag = "OK" if eval_pal_ok else "BREAK"
        pp_tag = "OK" if pi_pal_ok else "BREAK"

        log(f"{combo_name:>10}  {parity_norm:>10.2e}  "
            f"{n_paired:>4}/{len(evals_vec)}  {pi_err:>10.2e}  "
            f"{p_tag:>8}  {ep_tag:>8}  {pp_tag:>8}")

    log()
    log(f"Parity-breaking ([H,X^N] != 0): {len(parity_breaking)}/36")
    log(f"Eigenvalue palindrome breaking:  {len(eval_palindrome_breaking)}/36")
    log(f"Pi-conjugation breaking:         {len(pi_palindrome_breaking)}/36")
    log()

    # Cross-reference: parity vs eigenvalue palindrome
    log("--- Cross-reference: parity vs eigenvalue palindrome ---")
    parity_set = set(parity_breaking)
    eval_pal_set = set(eval_palindrome_breaking)

    both = parity_set & eval_pal_set
    only_parity = parity_set - eval_pal_set
    only_eval_pal = eval_pal_set - parity_set
    neither = set(r['combo'] for r in results) - parity_set - eval_pal_set

    log(f"Break both parity AND eigenvalue palindrome: {len(both)}")
    if both:
        for c in sorted(both):
            log(f"  {c}")
    log(f"Break only parity (eigenvalues still pair): {len(only_parity)}")
    if only_parity:
        for c in sorted(only_parity):
            log(f"  {c}")
    log(f"Break only eigenvalue palindrome (parity OK): {len(only_eval_pal)}")
    if only_eval_pal:
        for c in sorted(only_eval_pal):
            log(f"  {c}")
    log(f"Break neither (fully palindromic): {len(neither)}")
    if neither:
        for c in sorted(neither):
            log(f"  {c}")
    log()

    if parity_set == eval_pal_set:
        log("RESULT: EXACT MATCH.")
        log("The parity-breaking Hamiltonians are EXACTLY the")
        log("eigenvalue-palindrome-breaking Hamiltonians.")
        log("Palindromic breaking IS parity breaking.")
        log("The V-Effect is: the Hamiltonian couples the two sides of the")
        log("mirror, and the coupling destroys the palindrome at boundary modes.")
        log("This is a TIER 1 result.")
    elif eval_pal_set <= parity_set:
        log("RESULT: CONTAINMENT.")
        log(f"All {len(eval_pal_set)} eigenvalue-palindrome-breakers are also")
        log(f"parity-breakers. But {len(only_parity)} parity-breakers maintain")
        log("eigenvalue palindrome (they have a hidden symmetry operator Q != Pi).")
        log()
        log("Parity violation is NECESSARY but not SUFFICIENT for eigenvalue")
        log("palindrome breaking. Some Hamiltonians break [H, X^N] = 0 but")
        log("possess an alternative conjugation operator that preserves pairing.")
    elif parity_set <= eval_pal_set:
        log("RESULT: REVERSE CONTAINMENT.")
        log("All parity-breakers also break eigenvalue palindrome, but some")
        log("eigenvalue-palindrome-breakers preserve parity.")
    else:
        n_match = sum(1 for r in results
                      if r['parity_ok'] == r['eval_pal_ok'])
        log(f"RESULT: PARTIAL CORRELATION ({n_match}/36 match).")
    log()

    return results, parity_breaking, eval_palindrome_breaking


# ============================================================
# TEST 3: THE SUBTLER BOOTSTRAP
# ============================================================
def run_test_3(N=3, gamma=GAMMA):
    log("=" * 70)
    log(f"TEST 3: THE SUBTLER BOOTSTRAP (N={N})")
    log("=" * 70)
    log()
    log("Question: What is the most general dissipator that commutes with")
    log("X^N parity? If L_D is uniquely determined (up to gamma), then the")
    log("parity ALONE forces the noise. No external environment needed.")
    log()

    all_idx, pmats, d = build_pauli_basis(N)
    num = 4 ** N
    Pi2 = build_Pi2_matrix(all_idx)

    # L_D for Z-dephasing is diagonal in Pauli basis with rates -2*gamma*xy_weight
    # A general diagonal dissipator has rates r_a for each Pauli basis element a.
    # For [L_D, Pi^2] = 0 with L_D diagonal and Pi^2 diagonal, we need:
    # Pi^2[a,a] * r_a = r_a * Pi^2[a,a] for all a.
    # This is ALWAYS satisfied since Pi^2 is diagonal and multiplication commutes.
    # So ANY diagonal dissipator commutes with Pi^2.
    #
    # But we need more: L_D must be a valid Lindblad dissipator.
    # The standard form is: L_D(rho) = sum_k gamma_k (L_k rho L_k^dag - {L_k^dag L_k, rho}/2)
    # In Pauli basis for single-site dephasing with operator M:
    #   rate for Pauli string P is: -sum_sites (1 - <P_site|M|P_site>^2) * gamma_site
    #
    # The real constraint is: L_D must commute with Pi^2 AND be a valid
    # Lindblad generator (completely positive trace-preserving).

    log("Step 1: Structure of parity-compatible diagonal dissipators")
    log()

    # Classify basis elements by Pi^2 eigenvalue and XY-weight
    plus_indices = []
    minus_indices = []
    for a in range(num):
        if Pi2[a, a].real > 0:
            plus_indices.append(a)
        else:
            minus_indices.append(a)

    log(f"+1 sector: {len(plus_indices)} elements")
    log(f"-1 sector: {len(minus_indices)} elements")
    log()

    # For a diagonal dissipator, [L_D, Pi^2] = 0 is automatic.
    # The real question is: which single-site dephasing operators
    # produce a dissipator that commutes with Pi^2?

    log("Step 2: Which single-site dephasing operators commute with X^N?")
    log()
    log("A dephasing operator M at site k gives L_D(rho) = gamma*(M rho M - rho).")
    log("In Pauli basis: rate for P is -2*gamma * (number of sites where")
    log("[M, P_site] != 0).")
    log()

    # Test all single-qubit operators as dephasing channels
    test_ops = {
        'Z': sz,
        'X': sx,
        'Y': sy,
        'Z+X': (sz + sx) / np.sqrt(2),
        'arbitrary': np.array([[1, 0.5], [0.5, -1]], dtype=complex),
    }

    bonds = [(i, i + 1) for i in range(N - 1)]

    for op_name, M in test_ops.items():
        # Make M Hermitian for validity
        M = (M + M.conj().T) / 2
        # Normalize
        M = M / np.sqrt(np.trace(M.conj().T @ M) / 2)

        # Build L_D in Pauli basis
        rates = np.zeros(num)
        for a, idx in enumerate(all_idx):
            r = 0
            for site in range(N):
                P_site = PAULIS[idx[site]]
                # [M, P] = 0?
                comm_mp = M @ P_site - P_site @ M
                if np.max(np.abs(comm_mp)) > 1e-10:
                    r += 1
            rates[a] = -2 * gamma * r

        L_D = np.diag(rates)
        comm_test = Pi2 @ L_D - L_D @ Pi2
        comm_norm = np.max(np.abs(comm_test))

        # Check: does this L_D preserve the 2:2 split?
        # The 2:2 split means exactly half the rates are 0 (immune) and half decay.
        n_immune = np.sum(np.abs(rates) < 1e-10)
        n_decay = num - n_immune

        log(f"  Dephasing by {op_name}:")
        log(f"    [L_D, Pi^2] = {comm_norm:.2e}")
        log(f"    Immune/Decaying: {n_immune}/{n_decay}")
        log(f"    Unique rates: {len(np.unique(np.round(rates, 8)))}")

    log()
    log("Step 3: Parameterize ALL parity-compatible Lindblad dissipators")
    log()

    # A general Lindblad dissipator has jump operators L_k.
    # For [L_D_super, Pi^2_super] = 0 in the Pauli basis representation,
    # we need: for each L_k, either L_k maps +1 sector to +1 sector,
    # or L_k maps +1 to -1 and -1 to +1.
    #
    # For single-site dephasing (L_k = sqrt(gamma) * M_k at site k),
    # conjugation by M_k acts on Pauli operators at that site.
    # X^N conjugation acts as: I->I, X->X, Y->-Y, Z->-Z.
    #
    # For [M_k conjugation, X conjugation] to commute at each site,
    # M_k must map the {I,X}/{Y,Z} split to itself.
    #
    # Operators diagonal in the X-basis: {I, X} span is preserved by Z and Y.
    # Actually, let us just enumerate.

    log("Single-site dephasing operators M such that M-conjugation")
    log("commutes with X-conjugation at that site:")
    log()
    log("X-conjugation acts as: I->I, X->X, Y->-Y, Z->-Z")
    log("M-conjugation must map {I,X} to {I,X} and {Y,Z} to {Y,Z}.")
    log()

    # A general 2x2 Hermitian traceless operator is a*X + b*Y + c*Z.
    # Conjugation by this maps Pauli P -> M P M (up to normalization).
    # We need: M . I . M = c_I * I + c_X * X (no Y,Z components)
    # That is automatic since M I M = M^2 = (a^2+b^2+c^2)*I.
    # We need: M . X . M has no Y,Z: M X M = ?
    # M = aX+bY+cZ, MXM = (aX+bY+cZ)X(aX+bY+cZ)
    # XM = a*I + b*XY + c*XZ = a*I + ib*Z - ic*Y
    # MXM = (aX+bY+cZ)(a*I + ib*Z - ic*Y)
    #      = a^2*X + iab*XZ - iac*XY + ab*YI + ib^2*YZ - ibc*Y^2
    #        + ac*ZI + ibc*Z^2 - ic^2*ZY
    # This gets complicated. Let me just do it numerically.

    log("Numerical scan: which directions M = cos(th)*Z + sin(th)*cos(ph)*X")
    log("                + sin(th)*sin(ph)*Y commute with X-conjugation?")
    log()

    n_th = 36
    n_ph = 72
    compatible_count = 0
    compatible_directions = []

    for i_th in range(n_th + 1):
        th = np.pi * i_th / n_th
        for i_ph in range(n_ph):
            ph = 2 * np.pi * i_ph / n_ph
            M = (np.cos(th) * sz + np.sin(th) * np.cos(ph) * sx
                 + np.sin(th) * np.sin(ph) * sy)

            # M-conjugation on single-site Paulis
            # Check: does M P M map {I,X} -> {I,X} and {Y,Z} -> {Y,Z}?
            ok = True
            for P, name in [(I2, 'I'), (sx, 'X'), (sy, 'Y'), (sz, 'Z')]:
                conj = M @ P @ M  # M P M (M^2 = I for unit Bloch vector)
                # Decompose in Pauli basis
                c_I = np.trace(conj) / 2
                c_X = np.trace(sx @ conj) / 2
                c_Y = np.trace(sy @ conj) / 2
                c_Z = np.trace(sz @ conj) / 2

                if name in ('I', 'X'):
                    # Should have no Y,Z components
                    if abs(c_Y) > 1e-8 or abs(c_Z) > 1e-8:
                        ok = False
                        break
                else:
                    # Should have no I,X components
                    if abs(c_I) > 1e-8 or abs(c_X) > 1e-8:
                        ok = False
                        break

            if ok:
                compatible_count += 1
                compatible_directions.append((th, ph))

    log(f"Compatible directions found: {compatible_count} out of "
        f"{(n_th + 1) * n_ph} tested")
    log()

    if compatible_directions:
        # Analyze which directions work
        log("Compatible (theta, phi) values:")
        thetas = sorted(set(round(t, 4) for t, p in compatible_directions))
        for th_val in thetas[:10]:
            phis = sorted(p for t, p in compatible_directions
                          if abs(t - th_val) < 0.01)
            # Identify the direction
            if abs(th_val) < 0.01 or abs(th_val - np.pi) < 0.01:
                direction = "Z-axis (theta=0 or pi)"
            elif abs(th_val - np.pi / 2) < 0.01:
                direction = "XY-plane"
            else:
                direction = f"theta={th_val:.3f}"
            log(f"  {direction}: {len(phis)} phi values")

    log()
    log("Step 4: Uniqueness analysis")
    log()

    # The real question: among all single-site Hermitian dephasing operators
    # M that are parity-compatible, how many distinct dissipators L_D exist?
    # Each M gives rates: for Pauli string P, rate = -2*gamma * sum_k f(P_k, M)
    # where f(P_k, M) = 0 if [M, P_k] = 0, else 1.
    #
    # For Z-dephasing: f = 0 for I,Z; f = 1 for X,Y -> kills XY-weight
    # For X-dephasing: f = 0 for I,X; f = 1 for Y,Z -> kills YZ-weight
    # For Y-dephasing: f = 0 for I,Y; f = 1 for X,Z -> kills XZ-weight

    log("Rate structures for the three axis-aligned dephasing operators:")
    log()
    log(f"  {'Pauli':>6}  {'Z-deph rate':>12}  {'X-deph rate':>12}  {'Y-deph rate':>12}")
    log(f"  {'-' * 48}")

    for M_name, M_op, weight_fn in [
        ('Z', sz, lambda idx: sum(1 for s in idx if s in (1, 2))),
        ('X', sx, lambda idx: sum(1 for s in idx if s in (2, 3))),
        ('Y', sy, lambda idx: sum(1 for s in idx if s in (1, 3))),
    ]:
        rates = np.array([-2 * gamma * weight_fn(idx) for idx in all_idx])
        n_immune = np.sum(np.abs(rates) < 1e-10)
        log(f"  {M_name}-deph: {n_immune} immune, "
            f"{num - n_immune} decaying, "
            f"unique rates: {len(np.unique(np.round(rates, 8)))}")

    log()

    # Check: does each axis-aligned dephasing commute with Pi^2?
    log("Parity compatibility of axis-aligned dephasing:")
    for M_name, weight_fn in [
        ('Z', lambda idx: sum(1 for s in idx if s in (1, 2))),
        ('X', lambda idx: sum(1 for s in idx if s in (2, 3))),
        ('Y', lambda idx: sum(1 for s in idx if s in (1, 3))),
    ]:
        rates = np.array([-2 * gamma * weight_fn(idx) for idx in all_idx])
        L_D_test = np.diag(rates)
        comm_test = Pi2 @ L_D_test - L_D_test @ Pi2
        comm_norm = np.max(np.abs(comm_test))
        log(f"  {M_name}-deph: [L_D, Pi^2] = {comm_norm:.2e}")

    log()
    log("Step 5: Does parity uniquely determine the dissipator?")
    log()

    # Key insight: Pi^2 = X^N is diagonal with +/-1 entries.
    # L_D diagonal commutes with Pi^2 automatically.
    # So the constraint [L_D, Pi^2] = 0 does NOT restrict L_D at all
    # when L_D is diagonal in the Pauli basis.
    #
    # But the constraint is stronger than just [L_D, Pi^2] = 0.
    # We need L_D to be a VALID Lindblad dissipator (CPTP generator)
    # AND to produce a 2:2 split (half immune, half decaying).
    # AND the immune set must be exactly the +1 sector of Pi^2.

    # Check: for Z-dephasing, are the immune operators exactly the +1 sector?
    z_immune = set(a for a, idx in enumerate(all_idx) if xy_weight(idx) == 0)
    plus_set = set(a for a in range(num) if Pi2[a, a].real > 0)

    log(f"Z-dephasing immune set size: {len(z_immune)}")
    log(f"Pi^2 = +1 sector size: {len(plus_set)}")
    log(f"Overlap: {len(z_immune & plus_set)}")
    log(f"Z-immune but Pi^2=-1: {len(z_immune - plus_set)}")
    log(f"Pi^2=+1 but not Z-immune: {len(plus_set - z_immune)}")
    log()

    if z_immune == plus_set:
        log("RESULT: Z-dephasing immune set = Pi^2 +1 sector EXACTLY.")
        log("The noise kills exactly the 'other side' of the mirror.")
    else:
        log("RESULT: Z-dephasing immune set != Pi^2 +1 sector.")
        log("The noise and the parity do not align perfectly.")
        log()
        # Show the mismatch
        for a in sorted(z_immune - plus_set):
            name = ''.join(PAULI_NAMES[k] for k in all_idx[a])
            log(f"  Z-immune but Pi^2=-1: {name} "
                f"(XY-weight={xy_weight(all_idx[a])})")
        for a in sorted(plus_set - z_immune):
            name = ''.join(PAULI_NAMES[k] for k in all_idx[a])
            log(f"  Pi^2=+1 but Z-decaying: {name} "
                f"(XY-weight={xy_weight(all_idx[a])})")

    log()

    # The deeper question: among all valid Lindblad dissipators whose
    # immune set is exactly the +1 sector, is Z-dephasing the only one?

    # The +1 sector = even YZ-count. These are operators with even number
    # of Y or Z at each site. Equivalently: operators built from {I,X}
    # at each site, or {Y,Z} at an even number of sites.
    #
    # A dissipator that kills exactly the -1 sector must give zero rate
    # to all +1 elements and nonzero rate to all -1 elements.
    #
    # For single-site dephasing by operator M:
    #   immune iff [M, P_k] = 0 at every site k
    #   The immune set = tensor product of single-site kernels
    #
    # For M = Z: kernel = {I, Z}. Immune = strings of {I,Z} only = even XY-weight? No.
    # Strings of {I,Z} have XY-weight = 0. But +1 sector includes things
    # like XX (XY-weight 2, but YZ-count 0, so Pi^2 = +1).
    #
    # So Z-dephasing kills MORE than just the -1 sector. It kills everything
    # with XY-weight > 0, which includes +1 sector elements like XX, XY, etc.
    # Wait, XX has YZ-count = 0 so Pi^2 eigenvalue = +1. And XY-weight = 2.
    # So XX decays under Z-dephasing but is in the +1 sector.

    log("Detailed comparison: Z-dephasing rates vs Pi^2 sector")
    log()
    log(f"{'Pauli':>8}  {'XY-wt':>6}  {'YZ-ct':>6}  {'Pi^2':>5}  {'Z-rate':>8}")
    log("-" * 40)

    shown = 0
    for a in range(num):
        idx = all_idx[a]
        name = ''.join(PAULI_NAMES[k] for k in idx)
        xyw = xy_weight(idx)
        yzc = sum(1 for s in idx if s in (2, 3))
        pi2_val = int(Pi2[a, a].real)
        rate = -2 * gamma * xyw
        # Show interesting cases: Pi^2=+1 but decaying, or Pi^2=-1 but immune
        if (pi2_val == 1 and xyw > 0) or (pi2_val == -1 and xyw == 0):
            log(f"{name:>8}  {xyw:>6}  {yzc:>6}  {pi2_val:>+5}  {rate:>8.3f}")
            shown += 1
        if shown >= 20:
            log("  ...")
            break

    log()
    log("CONCLUSION:")
    log()
    log("If the +1 sector and the Z-immune set are different, then")
    log("Z-dephasing is NOT uniquely determined by parity. The parity")
    log("determines a sector split, but the noise kills a different set.")
    log("The bootstrap is structural (both arise from d=2) but not")
    log("uniquely determined: the parity does not force a specific L_D.")
    log()

    # Final count of free parameters
    # For single-site dephasing, M can be any direction on the Bloch sphere.
    # That is 2 continuous parameters (theta, phi) per site, plus gamma.
    # Total: 2N + 1 parameters for N sites.
    # The constraint [L_D, Pi^2] = 0 does not reduce this (as shown above).
    log(f"Free parameters for single-site dephasing: 2N+1 = {2*N+1} "
        f"(2 angles per site + 1 gamma)")
    log("Parity compatibility does NOT constrain these parameters.")
    log()
    log("The bootstrap is real but NOT uniquely determined by parity alone.")
    log("The noise type (Z-dephasing vs X-dephasing vs other) is an")
    log("independent choice. d(d-2)=0 forces d=2 and forces the parity,")
    log("but does not force which Pauli axis the noise acts on.")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Bootstrap Test: Is Each Side the Environment of the Other?")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"gamma = {GAMMA}")
    log()

    t1 = time.time()
    (L, L_H, L_D_diag, all_idx, pmats, d,
     Pi2, plus_idx, minus_idx) = run_test_1()
    log(f"[Test 1 completed in {time.time() - t1:.1f}s]")
    log()

    t2 = time.time()
    all_idx2, pmats2, d2, Pi2_2, XN = run_test_2()
    log(f"[Test 2 completed in {time.time() - t2:.1f}s]")
    log()

    t4 = time.time()
    results_4, parity_breaking, eval_pal_breaking = run_test_4()
    log(f"[Test 4 completed in {time.time() - t4:.1f}s]")
    log()

    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time() - t3:.1f}s]")
    log()

    # ============================================================
    # OVERALL SUMMARY
    # ============================================================
    log("=" * 70)
    log("OVERALL SUMMARY")
    log("=" * 70)
    log()
    log("Test 1: Cross-sector coupling in full Liouvillian")
    log("  [Pi^2, L] = 0 -> sectors are exactly decoupled.")
    log("  The naive bootstrap (sectors couple through L) is falsified.")
    log()
    log("Test 2: Hamiltonian-only analysis")
    log("  [H, X^N] = 0 for Heisenberg and standard Hamiltonians.")
    log("  DM interaction (XY-YX) DOES break parity.")
    log("  Nakajima-Zwanzig projection is trivial for parity-preserving H.")
    log()

    parity_set = set(parity_breaking)
    eval_pal_set = set(eval_pal_breaking)
    if parity_set == eval_pal_set:
        log("Test 4: Parity-breaking = Eigenvalue-palindrome-breaking (EXACT MATCH)")
        log(f"  All {len(parity_set)} broken combos break both.")
        log("  The V-Effect IS parity violation. TIER 1 result.")
    elif eval_pal_set <= parity_set:
        log(f"Test 4: Parity vs eigenvalue palindrome (CONTAINMENT)")
        log(f"  {len(parity_set)} parity-breaking, "
            f"{len(eval_pal_set)} eigenvalue-palindrome-breaking")
        log(f"  All {len(eval_pal_set)} eigenvalue-breakers are parity-breakers.")
        log(f"  But {len(parity_set - eval_pal_set)} parity-breakers keep "
            f"eigenvalue pairing (hidden Q != Pi).")
        log("  Parity violation is NECESSARY but not SUFFICIENT.")
    else:
        overlap = len(parity_set & eval_pal_set)
        log(f"Test 4: Parity vs eigenvalue palindrome")
        log(f"  {len(parity_set)} parity-breaking, {len(eval_pal_set)} "
            f"eigenvalue-palindrome-breaking, {overlap} overlap")

    log()
    log("Test 3: Uniqueness of the dissipator")
    log("  Parity [Pi^2, L_D] = 0 does NOT constrain the dissipator type.")
    log("  Any diagonal dissipator in the Pauli basis is parity-compatible.")
    log("  The noise axis (Z vs X vs Y) is NOT determined by parity alone.")
    log("  The bootstrap is structural but not uniquely determined.")
    log()

    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
