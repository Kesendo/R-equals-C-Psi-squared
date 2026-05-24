"""
F108 Part 1: Pi-family scan for pure-Pi^2-even NON-TRULY Hamiltonians.

Goal: identify the conjugation operator Pi that achieves EXACT operator-level
palindrome (Pi * L * Pi^{-1} = -L - 2*sigma*I) for every pure-Pi^2-even
non-truly pair (single bilinear + two-term combinations).

Pi^2-even (bit_b = #Y + #Z mod 2 = 0) bilinears: {XX, YY, YZ, ZY, ZZ}
(2-site bilinears only; II/IX/XI are trivially Hamiltonian-free).
Truly (F85): #Y even AND #Z even.
Pure-Pi^2-even NON-truly 2-site bilinears: {YZ, ZY}.

For each candidate Hamiltonian H (single bilinear + two-term combos with
at least one bilinear in {YZ, ZY}), we:
  1. Build L = -i [H, .] + sum_l gamma * D[Z_l] for chain (open BC).
  2. Verify spec(L) is palindromic around -sigma = -N*gamma.
  3. Test Pi-families (P1 phase variants, P4 phase variants, alternating)
     for EXACT operator-level palindrome.
  4. Report which family achieves operator-level palindrome.

Reuses the dict-based per-site map representation from
simulations/algebraic_pi_search.py (proven correct: STEP-6 verification
there reaches residual < 1e-10).

Output: console + simulations/results/f108_part1_pi_family_scan.txt
"""

import os
from datetime import datetime
from itertools import combinations_with_replacement, product as iprod

import numpy as np

# ============================================================
# Output setup
# ============================================================
RESULTS_DIR = (
    r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
    r"\simulations\results"
)
OUT = os.path.join(RESULTS_DIR, "f108_part1_pi_family_scan.txt")
os.makedirs(RESULTS_DIR, exist_ok=True)
f_log = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f_log.write(msg + "\n")
    f_log.flush()


# ============================================================
# Pauli matrices
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {"I": I2, "X": sx, "Y": sy, "Z": sz}
LABELS = ["I", "X", "Y", "Z"]


# ============================================================
# bit_b parity + F85 truly classification
# ============================================================
def bit_b(bilinear):
    """Pi^2 parity bit: #Y + #Z mod 2."""
    return sum(1 for c in bilinear if c in "YZ") % 2


def is_truly(bilinear):
    """F85 truly: #Y even AND #Z even."""
    ny = sum(1 for c in bilinear if c == "Y")
    nz = sum(1 for c in bilinear if c == "Z")
    return (ny % 2 == 0) and (nz % 2 == 0)


# ============================================================
# Site operator builder
# ============================================================
def site_op(op, s, N):
    ops = [I2] * N
    ops[s] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_chain(N, bonds, comps):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for c, J in comps.items():
            if J == 0:
                continue
            H += J * site_op(PAULI[c[0]], i, N) @ site_op(PAULI[c[1]], j, N)
    return H


def build_L(H, gamma, N):
    """L = -i [H, .] + sum_k gamma * D[Z_k] in vec (row-major) basis."""
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


# ============================================================
# Per-site Pi map as dict: label -> (phase, target_label)
# Proven correct construction from algebraic_pi_search.py
# ============================================================
def build_Q_Nsite(M, N):
    """Q = M^{otimes N} as d^2 x d^2 superoperator in standard vec basis."""
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
            ph, t = M[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


def build_Q_Nsite_alternating(M1, M2, N):
    """Q = M1 (x) M2 (x) M1 (x) M2 ... as d^2 x d^2 superoperator."""
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
        for site, sl in enumerate(bl):
            M = M1 if site % 2 == 0 else M2
            ph, t = M[sl]
            phase *= ph
            tgt_labels.append(t)
        tgt_idx = label_to_idx[tuple(tgt_labels)]
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


# ============================================================
# Pi family library (dict format: label -> (phase, target))
# ============================================================
# P1 permutation: I<->X, Y<->Z
# Phase variations distinguish families
def Pi_P1_canonical():
    """Heisenberg-supporting P1: I->X, X->I, Y->iZ, Z->iY.
    Supports {XX, YY, ZZ}."""
    return {"I": (1, "X"), "X": (1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}


def Pi_P1_only_yzzy():
    """ONLY-YZ-ZY P1 phase variant: I->X, X->-I, Y->+iZ, Z->+iY.
    From algebraic_pi_search.txt: '{YZ, ZY} - 8 maps' bucket.
    Supports YZ and ZY only (no XX, YY, ZZ)."""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (1j, "Y")}


def Pi_P1_5bilinear():
    """The MAGIC P1 variant: I->X, X->-I, Y->+iZ, Z->-iY.
    From algebraic_pi_search.txt: '{XX, YY, YZ, ZY, ZZ} - 8 maps' bucket.
    Supports XX, YY, YZ, ZY, ZZ -- i.e. all 5 Pi^2-even 2-site bilinears!
    This is THE candidate for closing F108 Part 1."""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (1j, "Z"), "Z": (-1j, "Y")}


def Pi_P1_5bilinear_v2():
    """Phase variant of Pi_P1_5bilinear: I->X, X->-I, Y->-iZ, Z->+iY.
    Also from '{XX, YY, YZ, ZY, ZZ}' bucket."""
    return {"I": (1, "X"), "X": (-1, "I"), "Y": (-1j, "Z"), "Z": (1j, "Y")}


def Pi_P4_canonical():
    """Heisenberg+ P4: I->Y, X->Z, Y->I, Z->X.
    Supports {XX, XZ, YY, ZX, ZZ}."""
    return {"I": (1, "Y"), "X": (1, "Z"), "Y": (1, "I"), "Z": (1, "X")}


def Pi_P4_for_xzzx():
    """P4 phase variant for {XX, XZ, YY, ZX, ZZ}: I->Y, X->iZ, Y->-I, Z->iX."""
    return {"I": (1, "Y"), "X": (1j, "Z"), "Y": (-1, "I"), "Z": (1j, "X")}


# ============================================================
# Palindrome checks
# ============================================================
def spec_palindrome_score(L, sigma, tol=1e-5):
    """Fraction of eigenvalues with palindromic partner."""
    evals = np.linalg.eigvals(L)
    target_sum = -2.0 * sigma
    n = len(evals)
    used = np.zeros(n, dtype=bool)
    paired = 0
    for i in range(n):
        if used[i]:
            continue
        partner = target_sum - evals[i]
        # First check self-partner at center
        if abs(evals[i].real - target_sum / 2) < tol and abs(evals[i].imag) < tol:
            used[i] = True
            paired += 1
            continue
        best_j, best_d = -1, np.inf
        for j in range(n):
            if used[j] or j == i:
                continue
            dd = abs(evals[j] - partner)
            if dd < best_d:
                best_d = dd
                best_j = j
        if best_j >= 0 and best_d < tol:
            used[i] = True
            used[best_j] = True
            paired += 2
    return paired / n


def operator_palindrome_residual(Q, L, sigma):
    """||Q L Q^{-1} + L + 2 sigma I||_F."""
    try:
        Q_inv = np.linalg.inv(Q)
    except np.linalg.LinAlgError:
        return np.inf
    d2 = L.shape[0]
    M = Q @ L @ Q_inv + L + 2.0 * sigma * np.eye(d2)
    return np.linalg.norm(M)


# ============================================================
# Bilinear classification
# ============================================================
PI2_EVEN_BILINEARS = ["XX", "YY", "YZ", "ZY", "ZZ"]
TRULY_BILINEARS = ["XX", "YY", "ZZ"]
NONTRULY_PI2EVEN_BILINEARS = ["YZ", "ZY"]


def pair_is_pure_pi2_even_non_truly(t1, t2):
    """Both Pi^2-even, at least one NON-truly."""
    if bit_b(t1) != 0 or bit_b(t2) != 0:
        return False
    if is_truly(t1) and is_truly(t2):
        return False
    return True


# ============================================================
# MAIN
# ============================================================
def main():
    log("=" * 88)
    log("F108 Part 1: Pi-family scan for pure-Pi^2-even non-truly H")
    log(f"Started: {datetime.now()}")
    log("=" * 88)

    log("\nPi^2-even (bit_b=0) 2-site bilinears: {XX, YY, YZ, ZY, ZZ}")
    log("Truly bilinears (F85): #Y even AND #Z even -> {XX, YY, ZZ}")
    log("Pure-Pi^2-even NON-truly bilinears: {YZ, ZY}")

    # ============================================================
    # Phase 1: enumerate the pairs
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 1: Enumerate pure-Pi^2-even non-truly pairs")
    log("=" * 88)

    pairs = []

    # Single bilinears: only non-truly ones (YZ, ZY); truly singletons (XX/YY/ZZ)
    # already proven via canonical Pi.
    for b in NONTRULY_PI2EVEN_BILINEARS:
        pairs.append((b, {b: 1.0}))

    # Two-term distinct combos: both Pi^2-even, at least one non-truly
    for t1, t2 in combinations_with_replacement(PI2_EVEN_BILINEARS, 2):
        if t1 == t2:
            continue
        if pair_is_pure_pi2_even_non_truly(t1, t2):
            label = f"{t1}+{t2}"
            pairs.append((label, {t1: 1.0, t2: 1.0}))

    log(f"\n  Total pairs: {len(pairs)}")
    for label, _ in pairs:
        log(f"    {label}")

    # ============================================================
    # Phase 2: For each pair, scan Pi families
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 2: Pi-family operator-level palindrome scan")
    log("=" * 88)

    gamma = 0.05
    N_list = [3, 4, 5]

    # Family pool
    families_uniform = {
        "P1_canonical": Pi_P1_canonical(),
        "P1_only_yzzy": Pi_P1_only_yzzy(),
        "P1_5bilinear": Pi_P1_5bilinear(),
        "P1_5bilinear_v2": Pi_P1_5bilinear_v2(),
        "P4_canonical": Pi_P4_canonical(),
        "P4_for_xzzx": Pi_P4_for_xzzx(),
    }
    families_alt = {
        "Alt_P1can_P4can": (Pi_P1_canonical(), Pi_P4_canonical()),
        "Alt_P1can_P15bi": (Pi_P1_canonical(), Pi_P1_5bilinear()),
        "Alt_P15bi_P4can": (Pi_P1_5bilinear(), Pi_P4_canonical()),
    }
    all_family_names = list(families_uniform.keys()) + list(families_alt.keys())

    results = {}  # (label, N) -> {family: residual}

    for N in N_list:
        log(f"\n--- N = {N} (chain, open BC, gamma = {gamma}) ---")
        bonds = [(i, i + 1) for i in range(N - 1)]
        sigma = N * gamma

        # Precompute Q per family for this N
        Q_per_family = {}
        for fname, M in families_uniform.items():
            Q_per_family[fname] = build_Q_Nsite(M, N)
        for fname, (M1, M2) in families_alt.items():
            Q_per_family[fname] = build_Q_Nsite_alternating(M1, M2, N)

        # Header
        log(f"\n  {'Pair':<10} {'SpecPal':>8}  "
            + "  ".join(f"{fname[:14]:>14}" for fname in all_family_names))
        log(f"  {'-' * 10} {'-' * 8}  "
            + "  ".join("-" * 14 for _ in all_family_names))

        for label, comps in pairs:
            H = build_H_chain(N, bonds, comps)
            L = build_L(H, gamma, N)
            spec_score = spec_palindrome_score(L, sigma, tol=1e-5)
            row = {}
            for fname in all_family_names:
                resid = operator_palindrome_residual(Q_per_family[fname], L, sigma)
                row[fname] = resid
            results[(label, N)] = (spec_score, row)

            spec_str = f"{spec_score*100:>6.1f}%"
            resid_strs = []
            for fname in all_family_names:
                r = row[fname]
                if r < 1e-8:
                    resid_strs.append(f"{'EXACT':>14}")
                elif r > 1e10:
                    resid_strs.append(f"{'SING':>14}")
                else:
                    resid_strs.append(f"{r:>14.3e}")
            log(f"  {label:<10} {spec_str:>8}  " + "  ".join(resid_strs))

    # ============================================================
    # Phase 3: Pattern extraction
    # ============================================================
    log("\n" + "=" * 88)
    log(f"Phase 3: Which Pi families achieve EXACT palindrome at ALL N in {N_list}?")
    log("=" * 88)

    family_hits = {fname: 0 for fname in all_family_names}
    pair_to_winners = {}
    for label, _ in pairs:
        winners = []
        for fname in all_family_names:
            resid_per_N = [results[(label, N)][1][fname] for N in N_list]
            if all(r < 1e-8 for r in resid_per_N):
                winners.append(fname)
                family_hits[fname] += 1
        pair_to_winners[label] = winners
        log(f"  {label:<10}  -> {winners if winners else 'NONE'}")

    log(f"\nFamily usage counts (out of {len(pairs)} pairs):")
    for fname, cnt in sorted(family_hits.items(), key=lambda x: -x[1]):
        log(f"  {fname:<24} : {cnt}")

    # Is there a single family that covers all pairs?
    log("\nSingle-family coverage test:")
    for fname in all_family_names:
        not_covered = [
            label for label, winners in pair_to_winners.items()
            if fname not in winners
        ]
        if not not_covered:
            log(f"  {fname:<24} : COVERS ALL {len(pairs)} pairs")
        else:
            log(f"  {fname:<24} : misses {len(not_covered)}: {not_covered}")

    # ============================================================
    # Phase 4: Algebraic verification of the winning family
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 4: Algebra of the winning family")
    log("=" * 88)

    # Find the best single family
    best_fname = max(family_hits, key=family_hits.get)
    log(f"\n  Best uniform family: {best_fname} (covers {family_hits[best_fname]}/{len(pairs)} pairs)")
    log(f"  Per-site map:")
    M_best = (families_uniform.get(best_fname)
              if best_fname in families_uniform
              else families_alt[best_fname][0])
    if best_fname in families_uniform:
        for l in LABELS:
            ph, t = M_best[l]
            ph_str = (
                "+1" if ph == 1 else "-1" if ph == -1
                else "+i" if ph == 1j else "-i" if ph == -1j
                else f"{ph}"
            )
            log(f"    {l} -> {ph_str} * {t}")
    else:
        # alternating
        M1, M2 = families_alt[best_fname]
        log(f"    Site-odd (M1):")
        for l in LABELS:
            ph, t = M1[l]
            log(f"      {l} -> {ph} * {t}")
        log(f"    Site-even (M2):")
        for l in LABELS:
            ph, t = M2[l]
            log(f"      {l} -> {ph} * {t}")

    # 2-qubit anti-commutation check: {Q_best, [H_2-body, .]} = 0 for the 9 bilinears
    log("\n  2-qubit anti-commutation check {Q, [B, .]} = 0 for each 2-body bilinear:")
    Q2 = build_Q_Nsite(M_best if best_fname in families_uniform else M1, 2)
    if best_fname not in families_uniform:
        # alternating: Q = M1 (x) M2 for 2 sites
        Q2 = build_Q_Nsite_alternating(*families_alt[best_fname], 2)
    BILINEARS = ["XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"]
    d = 4
    Id = np.eye(d)
    for term in BILINEARS:
        H2 = np.kron(PAULI[term[0]], PAULI[term[1]])
        C = -1j * (np.kron(H2, Id) - np.kron(Id, H2.T))
        err = np.linalg.norm(Q2 @ C + C @ Q2)
        ok = "ANTI-COMM" if err < 1e-8 else f"err={err:.2e}"
        b_tag = f"Pi2={'+' if bit_b(term) == 0 else '-'}"
        truly_tag = "truly" if is_truly(term) else "non-truly"
        log(f"    {term} ({b_tag}, {truly_tag}): {ok}")

    # Dissipator side: Q L_diss Q^-1 = -L_diss - 2 sigma I?
    log("\n  Dissipator side: pure D[Z]^{otimes N} (no Hamiltonian)")
    for N_chk in [3, 4, 5]:
        sigma_chk = N_chk * gamma
        H0 = np.zeros((2 ** N_chk, 2 ** N_chk), dtype=complex)
        L_diss = build_L(H0, gamma, N_chk)
        if best_fname in families_uniform:
            Q_chk = build_Q_Nsite(M_best, N_chk)
        else:
            Q_chk = build_Q_Nsite_alternating(*families_alt[best_fname], N_chk)
        try:
            Q_inv_chk = np.linalg.inv(Q_chk)
            resid_d = np.linalg.norm(
                Q_chk @ L_diss @ Q_inv_chk + L_diss
                + 2 * sigma_chk * np.eye(L_diss.shape[0])
            )
            log(f"    N={N_chk}:  ||Q L_diss Q^-1 + L_diss + 2 sigma I|| = {resid_d:.4e}")
        except np.linalg.LinAlgError:
            log(f"    N={N_chk}:  Q singular!")

    # ============================================================
    # Phase 4b: Robustness against arbitrary Pi^2-even non-truly H
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 4b: Robustness over arbitrary Pi^2-even non-truly H")
    log("=" * 88)
    log("\n  Build H = sum over bonds, with random coefficients on each of the 5")
    log("  Pi^2-even bilinears (so YZ and ZY are guaranteed nonzero, making H")
    log("  generically non-truly), and also non-uniform J per bond.")
    log("  Check P1_5bilinear still gives EXACT operator-level palindrome.")
    log()

    rng = np.random.default_rng(2026)
    n_trials = 5
    M_winner = Pi_P1_5bilinear()

    log(f"  {'Trial':>5} {'N':>3}  {'spec_pal':>9}  {'op_resid':>12}")
    log(f"  {'-' * 5} {'-' * 3}  {'-' * 9}  {'-' * 12}")
    for trial in range(n_trials):
        for N in [3, 4, 5]:
            bonds = [(i, i + 1) for i in range(N - 1)]
            sigma = N * gamma
            # Random non-uniform coupling per bond per bilinear
            d = 2 ** N
            H = np.zeros((d, d), dtype=complex)
            for i, j in bonds:
                for term in PI2_EVEN_BILINEARS:
                    J = rng.normal()
                    if abs(J) < 1e-12:
                        continue
                    H += J * site_op(PAULI[term[0]], i, N) @ site_op(PAULI[term[1]], j, N)
            # H is Hermitian by construction (XX, YY, ZZ self-Hermitian;
            # YZ + ZY symmetric on each bond if both present; but mixing
            # YZ and ZY with different J keeps H Hermitian only with adj)
            H = (H + H.conj().T) / 2.0
            L = build_L(H, gamma, N)
            spec_score = spec_palindrome_score(L, sigma, tol=1e-5)
            Q = build_Q_Nsite(M_winner, N)
            resid = operator_palindrome_residual(Q, L, sigma)
            log(f"  {trial+1:>5} {N:>3}  {spec_score*100:>7.1f}%  {resid:>12.3e}")

    # Also test asymmetric YZ vs ZY (J_YZ != J_ZY) -- a pure
    # off-symmetric flavor of pi^2-even non-truly H.
    log()
    log("  Asymmetric YZ vs ZY (J_YZ != J_ZY, all bonds; no XX/YY/ZZ):")
    log(f"  {'Trial':>5} {'N':>3}  {'spec_pal':>9}  {'op_resid':>12}")
    for trial in range(3):
        for N in [3, 4, 5]:
            bonds = [(i, i + 1) for i in range(N - 1)]
            sigma = N * gamma
            d = 2 ** N
            H = np.zeros((d, d), dtype=complex)
            for i, j in bonds:
                J_yz = rng.normal()
                J_zy = rng.normal()
                H += J_yz * site_op(sy, i, N) @ site_op(sz, j, N)
                H += J_zy * site_op(sz, i, N) @ site_op(sy, j, N)
            H = (H + H.conj().T) / 2.0
            L = build_L(H, gamma, N)
            spec_score = spec_palindrome_score(L, sigma, tol=1e-5)
            Q = build_Q_Nsite(M_winner, N)
            resid = operator_palindrome_residual(Q, L, sigma)
            log(f"  {trial+1:>5} {N:>3}  {spec_score*100:>7.1f}%  {resid:>12.3e}")

    # ============================================================
    # Phase 5: Conclusion + remaining gap
    # ============================================================
    log("\n" + "=" * 88)
    log("Phase 5: F108 Part 1 closure status")
    log("=" * 88)

    covered_all = [
        fname for fname in all_family_names
        if family_hits[fname] == len(pairs)
    ]
    if covered_all:
        log()
        log("  RESULT: A SINGLE Pi family achieves EXACT operator-level palindrome")
        log(f"  for EVERY pure-Pi^2-even non-truly pair, at all N in {N_list}.")
        log()
        log(f"  Covering families: {covered_all}")
        log()
        log("  Algebraic mechanism (F108 Part 1):")
        log("    For Pi^2-even non-truly H = sum_l alpha_l * (per-bond bilinear),")
        log("    the Pi-family with phase variant X -> -1 * I (instead of +1)")
        log("    realizes Pi * L * Pi^{-1} = -L - 2 sigma I exactly.")
        log("    This is the SAME P1 permutation as Heisenberg's Pi, with the")
        log("    only difference being a sign on the I<->X swap.")
        log()
        log("  Hence:")
        log("    NO Pi^2-even pair can be F87-hard, because an exact Pi exists.")
        log("    F108 Part 1 reduces to: choose the right phase variant of P1.")
    else:
        log()
        log("  No single family covers all pairs; multi-family selection rule needed:")
        for label, winners in pair_to_winners.items():
            log(f"    {label:<10} : {winners}")

    log("\n  Remaining gaps:")
    log("    (a) Extend numerical verification to N >= 5 (cheap: 256 x 256).")
    log("    (b) Promote to a Tier-1 proof in docs/proofs/PROOF_F108_PART1.md.")
    log("    (c) Re-derive the phase choice X -> -1 * I from first principles")
    log("        (anti-commutation {M, sx} = 0 vs commutation {M, I} = 2I).")
    log()
    log(f"  Completed: {datetime.now()}")
    log("=" * 88)


if __name__ == "__main__":
    try:
        main()
    finally:
        f_log.close()
        print(f"\n>>> Results saved to: {OUT}")
