#!/usr/bin/env python3
"""
ALGEBRAIC PI OPERATOR SEARCH
==============================
Enumerate all per-site Pauli transformations for Z-dephasing.
Find which transformations anti-commute with each Hamiltonian term.
Build compatibility table explaining the 22/14 palindrome split.

Approach:
  - For Z-dephasing, Paulis split: undephased {I,Z}, dephased {X,Y}
  - Any per-site bijection swapping these groups (with phases ±1, ±i)
    automatically satisfies the dissipator condition
  - We enumerate all such maps and check which Hamiltonian terms
    each map anti-commutes with: {Q, [H,.]} = 0
  - Two terms are "compatible" iff they share a common Q

Output: simulations/results/algebraic_pi_search.txt
"""

import numpy as np
from itertools import product as iprod, combinations
from datetime import datetime
import os

# ================================================================
# Output setup (dual: console + file, line-buffered)
# ================================================================
RESULTS_DIR = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results"
OUT = os.path.join(RESULTS_DIR, "algebraic_pi_search.txt")
f = open(OUT, "w", buffering=1)

def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()

# ================================================================
# Pauli matrices and constants
# ================================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
LABELS = ['I', 'X', 'Y', 'Z']
PHASES = [1, -1, 1j, -1j]
H_TERMS = [a + b for a in 'XYZ' for b in 'XYZ']  # XX,XY,...,ZZ (9 terms)

def phase_str(p):
    """Human-readable phase."""
    for ref, s in [(1, "+1"), (-1, "-1"), (1j, "+i"), (-1j, "-i")]:
        if abs(p - ref) < 1e-10:
            return s
    return f"{p:.4f}"

def map_str(M):
    """Concise string for a per-site map."""
    return ", ".join(f"{l}->{phase_str(M[l][0])}{M[l][1]}" for l in LABELS)

# ================================================================
# Precompute 2-qubit Pauli basis (for 16x16 superoperator checks)
# ================================================================
BASIS_2Q = list(iprod(LABELS, repeat=2))           # 16 label-pairs
BASIS_2Q_IDX = {bl: i for i, bl in enumerate(BASIS_2Q)}
MATS_2Q = [np.kron(PAULI[l1], PAULI[l2]) for l1, l2 in BASIS_2Q]
VECS_2Q = [m.flatten() for m in MATS_2Q]           # row-major vectorization

# Precompute commutator superoperators [Pa⊗Pb, .] for all 9 terms
COMM = {}
for term in H_TERMS:
    H = np.kron(PAULI[term[0]], PAULI[term[1]])     # 4x4 Hamiltonian
    COMM[term] = np.kron(H, np.eye(4)) - np.kron(np.eye(4), H.T)  # 16x16


# ================================================================
# Step 1: Enumerate all valid per-site maps for Z-dephasing
# ================================================================

def enumerate_maps():
    """
    For Z-dephasing:
      Undephased: {I, Z}  (commute with Z_k)
      Dephased:   {X, Y}  (anti-commute with Z_k)

    Valid map: bijection on {I,X,Y,Z} that swaps the two groups,
    with independent phase ∈ {+1,-1,+i,-i} on each label.

    Returns list of dicts: label -> (phase, target_label)
    """
    maps = []
    # 4 permutation structures:
    #   I -> {X,Y}  (2 choices), Z -> other
    #   X -> {I,Z}  (2 choices), Y -> other
    for i_tgt in ['X', 'Y']:
        z_tgt = 'Y' if i_tgt == 'X' else 'X'
        for x_tgt in ['I', 'Z']:
            y_tgt = 'Z' if x_tgt == 'I' else 'I'
            perm = {'I': i_tgt, 'X': x_tgt, 'Y': y_tgt, 'Z': z_tgt}
            for phs in iprod(PHASES, repeat=4):
                maps.append({l: (phs[i], perm[l]) for i, l in enumerate(LABELS)})
    return maps


# ================================================================
# Step 2: Build Q superoperator and test anti-commutation
# ================================================================

def build_Q(M):
    """
    Build 16x16 superoperator Q = M⊗M for a 2-qubit system.

    Q maps vec(P_a⊗P_b) to φ_a·φ_b · vec(P_{M(a)}⊗P_{M(b)}).

    In standard (row-major vectorized) basis:
      Q = Σ_{a,b} (φ_a φ_b / d) |vec(P_{M(a)}⊗P_{M(b)})><vec(P_a⊗P_b)|
    where d = 4 (2-qubit dimension) = ||vec(P_a⊗P_b)||².
    """
    Q = np.zeros((16, 16), dtype=complex)
    for idx, (l1, l2) in enumerate(BASIS_2Q):
        ph1, t1 = M[l1]
        ph2, t2 = M[l2]
        tgt_idx = BASIS_2Q_IDX[(t1, t2)]
        Q += (ph1 * ph2 / 4.0) * np.outer(VECS_2Q[tgt_idx], VECS_2Q[idx].conj())
    return Q

def anticomm_norm(Q, term):
    """||{Q, [H_term, .]}|| = ||Q·C + C·Q||"""
    C = COMM[term]
    return np.linalg.norm(Q @ C + C @ Q)


# ================================================================
# N-site verification helpers
# ================================================================

def site_op(op, s, N):
    """Operator `op` on site s, identity elsewhere, for N qubits."""
    ops = [I2] * N
    ops[s] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def build_H_chain(N, bonds, comps):
    """Build Hamiltonian on N-site chain."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for c, J in comps.items():
            if J == 0:
                continue
            H += J * site_op(PAULI[c[0]], i, N) @ site_op(PAULI[c[1]], j, N)
    return H

def build_L(H, gammas, N):
    """Build Lindbladian L = -i[H,.] + Σ_k γ_k D[Z_k]."""
    d = 2 ** N
    d2 = d * d
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L

def build_Q_Nsite(M, N):
    """Build d²×d² superoperator Q = M^{⊗N} for N qubits."""
    d = 2 ** N
    d2 = d * d
    basis_N = list(iprod(LABELS, repeat=N))
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
        tgt_idx = basis_N.index(tuple(tgt_labels))
        Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
    return Q


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    log("=" * 90)
    log("ALGEBRAIC PI OPERATOR SEARCH")
    log(f"Started: {datetime.now()}")
    log("=" * 90)

    # ==================================================================
    # STEP 1: Enumerate all valid per-site maps
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 1: Enumerate valid per-site maps for Z-dephasing")
    log(f"{'=' * 90}")

    all_maps = enumerate_maps()
    log(f"  Undephased set: {{I, Z}} (commute with Z)")
    log(f"  Dephased set:   {{X, Y}} (anti-commute with Z)")
    log(f"  Valid permutations: 4")
    log(f"  Phase choices per permutation: 4^4 = 256")
    log(f"  Total candidate maps: {len(all_maps)}")

    log(f"\n  The 4 permutation structures:")
    perms_seen = set()
    perm_list = []
    for M in all_maps:
        p = tuple(M[l][1] for l in LABELS)
        if p not in perms_seen:
            perms_seen.add(p)
            perm_list.append(p)
            log(f"    P{len(perm_list)}: I->{p[0]}, X->{p[1]}, Y->{p[2]}, Z->{p[3]}")

    # Identify known Pi
    known_idx = None
    for idx, M in enumerate(all_maps):
        if (M['I'] == (1, 'X') and M['X'] == (1, 'I')
                and M['Y'] == (1j, 'Z') and M['Z'] == (1j, 'Y')):
            known_idx = idx
            break
    log(f"\n  Known Pi (index {known_idx}): {map_str(all_maps[known_idx])}")

    # ==================================================================
    # STEP 2: Test anti-commutation with each Hamiltonian term
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 2: Test anti-commutation {Q, [H,.]} = 0 for each map and term")
    log(f"{'=' * 90}")

    TOL = 1e-10
    term_valid = {t: [] for t in H_TERMS}   # term -> list of map indices
    map_support = {}                          # map_idx -> frozenset of terms

    for m_idx, M in enumerate(all_maps):
        Q = build_Q(M)
        support = set()
        for term in H_TERMS:
            if anticomm_norm(Q, term) < TOL:
                support.add(term)
                term_valid[term].append(m_idx)
        map_support[m_idx] = frozenset(support)

    log(f"\n  Valid maps per Hamiltonian term:")
    log(f"  {'Term':6s} {'#Maps':>6s}")
    log(f"  {'-' * 14}")
    for t in H_TERMS:
        log(f"  {t:6s} {len(term_valid[t]):6d}")

    # Check: all 9 individual terms have at least one valid map?
    all_covered = all(len(term_valid[t]) > 0 for t in H_TERMS)
    log(f"\n  All 9 terms have valid maps: {'YES' if all_covered else 'NO — problem!'}")
    if not all_covered:
        missing = [t for t in H_TERMS if not term_valid[t]]
        log(f"  MISSING: {missing}")

    # Verify known Pi supports exactly {XX, YY, ZZ}
    log(f"\n  Known Pi supports: {sorted(map_support[known_idx])}")

    # Group maps by their support pattern
    pat_groups = {}
    for m_idx, pat in map_support.items():
        pat_groups.setdefault(pat, []).append(m_idx)

    log(f"\n  Distinct support patterns: {len(pat_groups)}")
    for pat in sorted(pat_groups, key=lambda p: (-len(p), sorted(p) if p else [])):
        indices = pat_groups[pat]
        terms_str = ', '.join(sorted(pat)) if pat else '(none)'
        log(f"\n  {{{terms_str}}} — {len(indices)} maps")
        # Show representative examples (up to 4)
        for idx in indices[:4]:
            log(f"    {map_str(all_maps[idx])}")
        if len(indices) > 4:
            log(f"    ... ({len(indices)} total)")

    # ==================================================================
    # STEP 2b: Site-swap symmetry check
    # ==================================================================
    log(f"\n  Site-swap symmetry check:")
    log(f"  (Q=MxM commutes with site swap, so XY and YX should have same maps)")
    swap_pairs = [('XY', 'YX'), ('XZ', 'ZX'), ('YZ', 'ZY')]
    for t1, t2 in swap_pairs:
        s1 = set(term_valid[t1])
        s2 = set(term_valid[t2])
        log(f"  {t1} maps == {t2} maps: {s1 == s2} (|{t1}|={len(s1)}, |{t2}|={len(s2)})")

    # ==================================================================
    # STEP 2c: Non-uniform maps for missing terms (M1 != M2)
    # ==================================================================
    missing_terms = [t for t in H_TERMS if not term_valid[t]]
    if missing_terms:
        log(f"\n{'=' * 90}")
        log("STEP 2c: Non-uniform maps Q = M1 x M2 (M1 != M2) for missing terms")
        log(f"{'=' * 90}")
        log(f"  Missing terms with uniform Q: {missing_terms}")
        log(f"  Testing all 1024 x 1024 = {len(all_maps)**2} pairs...")

        # For chain compatibility, we need BOTH M1xM2 AND M2xM1 to work
        # (alternating pattern on chain)
        nonunif_valid = {t: [] for t in missing_terms}  # term -> list of (m1,m2) pairs

        # Precompute all Q matrices
        all_Qs = [build_Q(M) for M in all_maps]  # 1024 x 16 x 16

        # For non-uniform Q = M1_at_site1 x M2_at_site2:
        # We need to build Q differently - M1 at site 1, M2 at site 2
        def build_Q_nonuniform(M1, M2):
            """Build 16x16 superoperator with M1 at site 1 and M2 at site 2."""
            Q = np.zeros((16, 16), dtype=complex)
            for idx, (l1, l2) in enumerate(BASIS_2Q):
                ph1, t1 = M1[l1]
                ph2, t2 = M2[l2]
                tgt_idx = BASIS_2Q_IDX[(t1, t2)]
                Q += (ph1 * ph2 / 4.0) * np.outer(VECS_2Q[tgt_idx], VECS_2Q[idx].conj())
            return Q

        for term in missing_terms:
            C = COMM[term]
            count = 0
            for i, M1 in enumerate(all_maps):
                for j, M2 in enumerate(all_maps):
                    # Build Q = M1 x M2
                    Q12 = build_Q_nonuniform(M1, M2)
                    err12 = np.linalg.norm(Q12 @ C + C @ Q12)
                    if err12 < TOL:
                        # Also check M2 x M1 (needed for chain compatibility)
                        Q21 = build_Q_nonuniform(M2, M1)
                        err21 = np.linalg.norm(Q21 @ C + C @ Q21)
                        chain_ok = err21 < TOL
                        nonunif_valid[term].append((i, j, chain_ok))
                        count += 1
            log(f"\n  {term}: {count} valid non-uniform pairs (M1 x M2 where M1!=M2 allowed)")
            if count > 0:
                # Count chain-compatible pairs
                chain_count = sum(1 for _, _, ok in nonunif_valid[term] if ok)
                log(f"    Chain-compatible (both M1xM2 and M2xM1 work): {chain_count}")
                # Show examples
                for i, j, ok in nonunif_valid[term][:5]:
                    tag = " [chain-OK]" if ok else " [2-site only]"
                    log(f"    M1: {map_str(all_maps[i])}")
                    log(f"    M2: {map_str(all_maps[j])}{tag}")
                    log(f"    ---")
            else:
                log(f"    NO valid non-uniform Pauli-to-Pauli maps exist either!")
                log(f"    => Q for {term} is NOT a per-site Pauli permutation.")
                log(f"    => Must be a more general linear map or non-tensor-product.")

    # ==================================================================
    # STEP 3: Compatibility table (2-term combinations)
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 3: Compatibility table for 2-term combinations")
    log(f"{'=' * 90}")

    combo_shared = {}
    for t1, t2 in combinations(H_TERMS, 2):
        key = f"{t1}+{t2}"
        s1 = set(term_valid[t1])
        s2 = set(term_valid[t2])
        combo_shared[key] = s1 & s2

    log(f"\n  {'Combo':<10s} {'Compat':>8s} {'#Shared':>8s}")
    log(f"  {'-' * 28}")
    n_yes = n_no = 0
    for t1, t2 in combinations(H_TERMS, 2):
        key = f"{t1}+{t2}"
        n_shared = len(combo_shared[key])
        compat = n_shared > 0
        if compat:
            n_yes += 1
        else:
            n_no += 1
        tag = 'YES' if compat else 'NO'
        log(f"  {key:<10s} {tag:>8s} {n_shared:>8d}")

    log(f"\n  Compatible: {n_yes}/36")
    log(f"  Incompatible: {n_no}/36")

    # ==================================================================
    # STEP 4: Compare with numerical results from non_heisenberg_deep.py
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 4: Compare with numerical results")
    log(f"{'=' * 90}")

    # Broken combos from Section 2 of non_heisenberg_deep.txt (N=3, Z-deph, gamma=0.05)
    num_broken = {
        'XX+XY', 'XX+YX',
        'XY+XZ', 'XY+YY', 'XY+YZ', 'XY+ZX', 'XY+ZY',
        'XZ+YX', 'XZ+ZY',
        'YX+YY', 'YX+YZ', 'YX+ZX', 'YX+ZY',
        'YZ+ZX'
    }

    alg_broken = {k for k, v in combo_shared.items() if len(v) == 0}

    log(f"  Numerical broken combos: {len(num_broken)}")
    log(f"  Algebraic broken combos: {len(alg_broken)}")

    match = num_broken == alg_broken
    log(f"\n  MATCH: {'YES — Perfect agreement!' if match else 'NO — Discrepancy detected!'}")

    if not match:
        only_num = num_broken - alg_broken
        only_alg = alg_broken - num_broken
        if only_num:
            log(f"  Broken numerically but NOT algebraically: {sorted(only_num)}")
            log(f"    => These have a shared Q algebraically but break numerically")
        if only_alg:
            log(f"  Broken algebraically but NOT numerically: {sorted(only_alg)}")
            log(f"    => No per-site Pauli-to-Pauli Q exists, yet palindrome holds")
            log(f"    => May need non-Pauli-permutation Q or non-per-site Q")

    # ==================================================================
    # STEP 4b: Non-uniform Q for discrepant combos (Pauli-basis, batched)
    # ==================================================================
    nonunif_combo_valid = {}  # combo -> list of (i,j) pairs

    if not match and only_alg:
        log(f"\n{'=' * 90}")
        log("STEP 4b: Non-uniform Q = M1 x M2 for discrepant combos")
        log(f"{'=' * 90}")

        discrepant = sorted(only_alg)
        log(f"  Combos palindromic numerically but no uniform Q: {discrepant}")

        # Precompute per-site 4x4 matrices in Pauli label space
        M_mats = np.zeros((len(all_maps), 4, 4), dtype=complex)
        for m_idx, M in enumerate(all_maps):
            for a_idx, label in enumerate(LABELS):
                phase, target = M[label]
                b_idx = LABELS.index(target)
                M_mats[m_idx, b_idx, a_idx] = phase

        # Precompute commutator superops in Pauli basis: C_p = V^-1 C_std V
        V = np.column_stack(VECS_2Q)   # 16x16, columns = vec(Pauli strings)
        V_inv = V.conj().T / 4.0       # V^-1 = V^dag / 4 since V^dag V = 4I
        C_pauli = {}
        for term in H_TERMS:
            C_pauli[term] = V_inv @ COMM[term] @ V

        N_maps = len(all_maps)
        for combo in discrepant:
            t1, t2 = combo.split('+')
            C1 = C_pauli[t1]
            C2 = C_pauli[t2]

            chain_valid = []
            for i in range(N_maps):
                M1 = M_mats[i]

                # Batch build Q12 = kron(M1, M2) for all M2
                # Q12[n, 4a+b, 4c+d] = M1[a,c] * M2_all[n,b,d]
                Q12_batch = np.einsum('ac,nbd->nabcd', M1, M_mats).reshape(N_maps, 16, 16)

                # Check {Q12, C1} = 0 for all M2 (batch)
                anti1 = np.einsum('nij,jk->nik', Q12_batch, C1) + np.einsum('ij,njk->nik', C1, Q12_batch)
                norms1 = np.sqrt(np.sum(np.abs(anti1.reshape(N_maps, -1))**2, axis=1))
                pass1 = norms1 < TOL

                if not np.any(pass1):
                    continue

                # Check {Q12, C2} = 0 for passing M2s
                idx1 = np.where(pass1)[0]
                Q_sub = Q12_batch[idx1]
                anti2 = np.einsum('nij,jk->nik', Q_sub, C2) + np.einsum('ij,njk->nik', C2, Q_sub)
                norms2 = np.sqrt(np.sum(np.abs(anti2.reshape(len(idx1), -1))**2, axis=1))
                pass2 = norms2 < TOL
                idx2 = idx1[pass2]

                if len(idx2) == 0:
                    continue

                # Check reverse ordering Q21 for each passing pair
                for j in idx2:
                    Q21 = np.kron(M_mats[j], M1)
                    a1 = Q21 @ C1 + C1 @ Q21
                    a2 = Q21 @ C2 + C2 @ Q21
                    if np.linalg.norm(a1) < TOL and np.linalg.norm(a2) < TOL:
                        chain_valid.append((i, j))

            nonunif_combo_valid[combo] = chain_valid
            log(f"\n  {combo}: {len(chain_valid)} chain-compatible non-uniform pairs")
            if chain_valid:
                for ci, cj in chain_valid[:3]:
                    pi = ''.join(all_maps[ci][l][1] for l in LABELS)
                    pj = ''.join(all_maps[cj][l][1] for l in LABELS)
                    log(f"    M1[{pi}]: {map_str(all_maps[ci])}")
                    log(f"    M2[{pj}]: {map_str(all_maps[cj])}")
                    log(f"    ---")
                if len(chain_valid) > 3:
                    log(f"    ... ({len(chain_valid)} total)")
            else:
                log(f"    NO Pauli-to-Pauli Q (uniform or non-uniform)!")
                log(f"    => Q must involve Pauli MIXING or non-tensor-product structure")

        # Update match assessment
        all_resolved = all(len(nonunif_combo_valid.get(c, [])) > 0 for c in discrepant)
        log(f"\n  All discrepant combos resolved by non-uniform Q: {'YES' if all_resolved else 'NO'}")
        if all_resolved:
            log(f"  => The palindromic symmetry uses BOTH uniform and non-uniform")
            log(f"     per-site Pauli maps. XY/YX terms require alternating maps.")
            match = True  # Update for summary

    # ==================================================================
    # STEP 5: Shared Q details for palindromic combos
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 5: Shared Q operators for each palindromic 2-term combo")
    log(f"{'=' * 90}")

    for t1, t2 in combinations(H_TERMS, 2):
        key = f"{t1}+{t2}"
        shared = combo_shared[key]
        if not shared:
            continue
        # Group shared maps by permutation
        perm_groups = {}
        for m_idx in shared:
            M = all_maps[m_idx]
            pk = tuple(M[l][1] for l in LABELS)
            perm_groups.setdefault(pk, []).append(m_idx)
        perms_summary = ", ".join(
            f"P({''.join(pk)})×{len(idxs)}"
            for pk, idxs in sorted(perm_groups.items())
        )
        log(f"\n  {key}: {len(shared)} shared maps [{perms_summary}]")
        for idx in list(shared)[:2]:
            log(f"    {map_str(all_maps[idx])}")
        if len(shared) > 2:
            log(f"    ... ({len(shared)} total)")

    # ==================================================================
    # STEP 6: Confirm no shared Q for broken combos
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("STEP 6: Confirm no shared Q for broken combos")
    log(f"{'=' * 90}")

    for t1, t2 in combinations(H_TERMS, 2):
        key = f"{t1}+{t2}"
        if combo_shared[key]:
            continue
        # Show WHY: the valid sets are disjoint
        n1 = len(term_valid[t1])
        n2 = len(term_valid[t2])
        log(f"  {key}: 0 shared  (|{t1}|={n1}, |{t2}|={n2}, intersection empty)")

    # ==================================================================
    # VERIFICATION: Spot-check on actual N=3 Lindbladian
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("VERIFICATION: Spot-check algebraic Q on actual Lindbladian (N=3)")
    log(f"{'=' * 90}")

    N = 3
    gamma = 0.05
    gammas = [gamma] * N
    sg = sum(gammas)
    bonds = [(0, 1), (1, 2)]

    tests = []

    # Known Pi on Heisenberg
    tests.append(("Known Pi + Heisenberg", known_idx, {'XX': 1, 'YY': 1, 'ZZ': 1}))

    # One map per individual term
    for term in H_TERMS:
        if term_valid[term]:
            tests.append((f"Map for {term} alone", term_valid[term][0], {term: 1}))

    # DM interaction
    dm_shared = combo_shared.get('XY+YX', set())
    if dm_shared:
        tests.append(("Map for DM (XY-YX)", list(dm_shared)[0], {'XY': 1, 'YX': -1}))

    # ZZ+XY (palindromic cross-combo)
    zzxy = combo_shared.get('XY+ZZ', set())
    if zzxy:
        tests.append(("Map for ZZ+XY", list(zzxy)[0], {'ZZ': 1, 'XY': 1}))

    # All9 uniform
    all9_shared = None
    for m_idx, pat in map_support.items():
        if pat == frozenset(H_TERMS):
            all9_shared = m_idx
            break
    if all9_shared is not None:
        tests.append(("Map for All9", all9_shared,
                       {t: 1 for t in H_TERMS}))

    log(f"\n  N={N}, gamma={gamma}, bonds={bonds}")
    log(f"\n  {'Test':35s} {'||QLQ^-1+L+2Sg||':>18s}  Map(s)")
    log(f"  {'-' * 105}")

    for name, m_idx, comps in tests:
        M = all_maps[m_idx]
        Q_full = build_Q_Nsite(M, N)
        H = build_H_chain(N, bonds, comps)
        L = build_L(H, gammas, N)
        d2 = L.shape[0]
        Q_inv = np.linalg.inv(Q_full)
        residual = np.linalg.norm(Q_full @ L @ Q_inv + L + 2 * sg * np.eye(d2))
        status = 'OK' if residual < 1e-10 else f'{residual:.4e}'
        log(f"  {name:35s} {status:>18s}  {map_str(M)}")

    # Non-uniform verification: build alternating Q = M1 x M2 x M1 for N=3
    def build_Q_Nsite_alternating(M1, M2, N):
        """Build Q = M1 x M2 x M1 x M2 x ... for N sites."""
        d = 2 ** N
        d2 = d * d
        basis_N = list(iprod(LABELS, repeat=N))
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
            tgt_idx = basis_N.index(tuple(tgt_labels))
            Q += (phase / d) * np.outer(vecs_N[tgt_idx], vecs_N[idx].conj())
        return Q

    # Test non-uniform maps for XY, YX, DM, and discrepant combos
    if nonunif_combo_valid:
        log(f"\n  --- Non-uniform Q verification (alternating M1-M2-M1 for N=3) ---")
        nonunif_tests = []

        # XY alone (from Step 2c non-uniform search)
        if missing_terms and 'XY' in missing_terms and nonunif_valid.get('XY'):
            chain_pairs = [(i, j) for i, j, ok in nonunif_valid['XY'] if ok]
            if chain_pairs:
                ci, cj = chain_pairs[0]
                nonunif_tests.append(("Non-unif Q for XY", ci, cj, {'XY': 1}))

        # Discrepant combos
        for combo, pairs in nonunif_combo_valid.items():
            if pairs:
                ci, cj = pairs[0]
                t1, t2 = combo.split('+')
                nonunif_tests.append((f"Non-unif Q for {combo}", ci, cj,
                                      {t1: 1, t2: 1}))

        for name, ci, cj, comps in nonunif_tests:
            M1 = all_maps[ci]
            M2 = all_maps[cj]
            Q_full = build_Q_Nsite_alternating(M1, M2, N)
            H = build_H_chain(N, bonds, comps)
            L = build_L(H, gammas, N)
            d2 = L.shape[0]
            Q_inv = np.linalg.inv(Q_full)
            residual = np.linalg.norm(Q_full @ L @ Q_inv + L + 2 * sg * np.eye(d2))
            status = 'OK' if residual < 1e-10 else f'{residual:.4e}'
            m1s = map_str(M1)
            m2s = map_str(M2)
            log(f"  {name:35s} {status:>18s}  M1={m1s}")
            log(f"  {'':35s} {'':>18s}  M2={m2s}")

    # ==================================================================
    # SUMMARY: The Pi operator family
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("SUMMARY: The Pi Operator Family for Z-Dephasing")
    log(f"{'=' * 90}")

    log(f"\n  Each Pauli-pair Hamiltonian term has its own set of valid")
    log(f"  conjugation operators Q. Two terms are 'compatible' if and")
    log(f"  only if their valid Q-sets intersect (share a common Q).")

    log(f"\n  Pi families (by support set):")
    log(f"  {'Support set':<40s} {'#Maps':>6s}  Example")
    log(f"  {'-' * 90}")
    for pat in sorted(pat_groups, key=lambda p: (-len(p), sorted(p) if p else [])):
        if not pat:
            count = len(pat_groups[pat])
            log(f"  {'(no terms)':40s} {count:>6d}  {map_str(all_maps[pat_groups[pat][0]])}")
            continue
        terms_str = ', '.join(sorted(pat))
        count = len(pat_groups[pat])
        example_M = all_maps[pat_groups[pat][0]]
        log(f"  {terms_str:40s} {count:>6d}  {map_str(example_M)}")

    log(f"\n  Compatibility result:")
    log(f"    Compatible 2-term combos:   {n_yes}/36")
    log(f"    Incompatible 2-term combos: {n_no}/36")
    log(f"    Matches numerical data:     {'YES' if match else 'NO'}")

    if match:
        log(f"\n  CONCLUSION:")
        log(f"    The palindromic symmetry for Z-dephasing is generated by a")
        log(f"    FAMILY of per-site Pauli conjugation operators, not just one.")
        log(f"    Each operator swaps undephased {{I,Z}} <-> dephased {{X,Y}}.")
        log(f"    A 2-term combo is palindromic iff its terms share a common Q.")
        log(f"    The 14 broken combos have DISJOINT Q-sets => no single Q works")
        log(f"    => breaking at O(gamma^2) from Q-operator conflict.")

    # ==================================================================
    # BONUS: Compatibility matrix visualization
    # ==================================================================
    log(f"\n{'=' * 90}")
    log("BONUS: Compatibility matrix (. = compatible, X = broken)")
    log(f"{'=' * 90}")

    log(f"\n       {'  '.join(f'{t:2s}' for t in H_TERMS)}")
    for i, t1 in enumerate(H_TERMS):
        row = f"  {t1:2s}  "
        for j, t2 in enumerate(H_TERMS):
            if i == j:
                row += " . "
            elif i < j:
                key = f"{t1}+{t2}"
                row += " . " if combo_shared[key] else " X "
            else:
                key = f"{t2}+{t1}"
                row += " . " if combo_shared[key] else " X "
        log(row)

    log(f"\n{'=' * 90}")
    log(f"Completed: {datetime.now()}")
    log(f"{'=' * 90}")
    f.close()
    print(f"\n>>> Results saved to: {OUT}")
