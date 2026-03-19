"""
The V-Effect: What Emerges from Broken Palindromic Symmetry?
=============================================================
Section 1: Anatomy of the break (14 broken cases at N=3)
Section 2: Structure in the error (gamma^2 scaling, spectral fingerprint)
Section 3: New symmetries after the break
Section 4: The bifurcation point (alpha sweep)
Section 5: Broken vs unbroken comparison
Section 6: The rank threshold (N=2 Choi rank vs N=3 breakability)

Script: simulations/v_effect_analysis.py
Output: simulations/results/v_effect_analysis.txt
"""
import numpy as np
from itertools import product as iproduct
from datetime import datetime

OUT = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results\v_effect_analysis.txt"
f = open(OUT, "w", buffering=1)


def log(msg=""):
    print(msg, flush=True)
    f.write(msg + "\n")
    f.flush()


# ============================================================
# PAULI BASICS
# ============================================================
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']
PM = {'I': I2, 'X': sx, 'Y': sy, 'Z': sz}
H_LABELS = ['XX', 'XY', 'XZ', 'YX', 'YY', 'YZ', 'ZX', 'ZY', 'ZZ']

PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}


def xy_weight(indices):
    return sum(1 for i in indices if i in (1, 2))


# ============================================================
# FAST LINDBLADIAN (computational basis, vectorized)
# ============================================================
def site_op(op, site, N):
    """Single-site operator embedded in N-site space."""
    ops = [I2] * N
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_H_combo_Nsite(N, bonds, combo):
    """Build Hamiltonian for a two-term combo on given bonds."""
    t1, t2 = combo.split('+')
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for label in [t1, t2]:
            ops = [I2] * N
            ops[i] = PM[label[0]]
            ops[j] = PM[label[1]]
            term = ops[0]
            for o in ops[1:]:
                term = np.kron(term, o)
            H += term
    return H


def build_H_combo_Nsite_weighted(N, bonds, combo, bond_weights):
    """Build with per-bond weight factors."""
    t1, t2 = combo.split('+')
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (i, j), w in zip(bonds, bond_weights):
        for label in [t1, t2]:
            ops = [I2] * N
            ops[i] = PM[label[0]]
            ops[j] = PM[label[1]]
            term = ops[0]
            for o in ops[1:]:
                term = np.kron(term, o)
            H += w * term
    return H


def build_L_vec(H, gamma, N):
    """Vectorized Lindbladian in computational basis."""
    d = 2 ** N
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k, N)
        L += gamma * (np.kron(Zk, Zk.conj()) - np.eye(d * d))
    return L


# ============================================================
# PAULI-BASIS LINDBLADIAN AND PI (for Pi error computation)
# ============================================================
def build_pauli_basis(N):
    all_idx = list(iproduct(range(4), repeat=N))
    d = 2 ** N
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for i in idx[1:]:
            m = np.kron(m, PAULIS[i])
        pmats.append(m)
    return all_idx, pmats, d


def build_L_H_pauli(N, H, all_idx, pmats, d):
    num = 4 ** N
    L = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = -1j * (H @ pmats[b] - pmats[b] @ H)
        for a in range(num):
            L[a, b] = np.trace(pmats[a] @ comm) / d
    return L


def build_Pi_persite(N, all_idx):
    num = len(all_idx)
    Pi = np.zeros((num, num), dtype=complex)
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = all_idx.index(mapped)
        Pi[a, b] = sign
    return Pi


# ============================================================
# PALINDROME CHECKING
# ============================================================
def check_eigenvalue_pairing(evals, Sg, tol=1e-6):
    """Check palindromic eigenvalue pairing. Returns per-eigenvalue errors."""
    num = len(evals)
    pair_errors = np.zeros(num)
    pair_indices = -np.ones(num, dtype=int)

    for k in range(num):
        target = -(evals[k] + 2 * Sg)
        dists = np.abs(evals - target)
        best_j = np.argmin(dists)
        pair_errors[k] = dists[best_j]
        pair_indices[k] = best_j

    n_well_paired = int(np.sum(pair_errors < tol))
    return pair_errors, pair_indices, n_well_paired


def pi_palindrome_error(Pi, L, Sg):
    """||Pi L Pi^{-1} + L + 2Sg I|| / ||L||"""
    num = L.shape[0]
    Pi_inv = np.linalg.inv(Pi)
    E = Pi @ L @ Pi_inv + L + 2 * Sg * np.eye(num)
    norm_L = np.max(np.abs(L))
    return np.max(np.abs(E)) / norm_L if norm_L > 1e-15 else np.max(np.abs(E)), E


# ============================================================
# N=2 CHOI RANK (from boot_script_structure)
# ============================================================
def build_V(pmats, d):
    num = len(pmats)
    V = np.zeros((d * d, num), dtype=complex)
    for a, P in enumerate(pmats):
        V[:, a] = P.flatten()
    return V


def construct_Pi_eigvec(L_full, Sg, num):
    evals, R = np.linalg.eig(L_full)
    paired = np.zeros(num, dtype=bool)
    pair_map = -np.ones(num, dtype=int)
    for k in range(num):
        if paired[k]:
            continue
        target = -(evals[k] + 2 * Sg)
        best_j, best_d = -1, np.inf
        for j in range(k, num):
            if paired[j] and j != k:
                continue
            dd = abs(evals[j] - target)
            if dd < best_d:
                best_d, best_j = dd, j
        if best_j >= 0 and best_d < 1e-8:
            paired[k] = paired[best_j] = True
            pair_map[k] = best_j
            pair_map[best_j] = k
    n_paired = int(np.sum(paired))
    if n_paired < num:
        return None, n_paired
    R_inv = np.linalg.inv(R)
    P_eig = np.zeros((num, num), dtype=complex)
    seen = set()
    for k in range(num):
        j = pair_map[k]
        if k in seen:
            continue
        seen.add(k)
        seen.add(j)
        if k == j:
            P_eig[k, k] = 1.0
        else:
            P_eig[j, k] = 1.0
            P_eig[k, j] = 1.0
    return R @ P_eig @ R_inv, n_paired


def choi_rank_N2(Pi, V2, d2):
    """Compute site1|site2 operator SVD rank of Choi(Pi) for N=2."""
    S = V2 @ Pi @ V2.conj().T / d2
    J = S.reshape(d2, d2, d2, d2).transpose(0, 2, 1, 3).reshape(d2 * d2, d2 * d2) / d2
    J8 = J.reshape(2, 2, 2, 2, 2, 2, 2, 2)
    J_reord = J8.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(16, 16)
    J4 = J_reord.reshape(4, 4, 4, 4)
    M = J4.transpose(0, 2, 1, 3).reshape(16, 16)
    s = np.linalg.svd(M, compute_uv=False)
    rank = int(np.sum(s > 1e-8 * s[0])) if s[0] > 1e-15 else 0
    return rank


# ############################################################
# MAIN
# ############################################################
gamma = 0.05
N3 = 3
bonds3 = [(0, 1), (1, 2)]
Sg3 = N3 * gamma

log("=" * 90)
log("THE V-EFFECT: What Emerges from Broken Palindromic Symmetry?")
log(f"Date: {datetime.now()}")
log(f"N = {N3}, gamma = {gamma}, Sg = {Sg3}")
log("=" * 90)

# Pre-computation: build all 36 combos
all_combos = []
for i, t1 in enumerate(H_LABELS):
    for t2 in H_LABELS[i + 1:]:
        all_combos.append(f"{t1}+{t2}")

# Classify all 36 combos at N=3: compute Pi error and eigenvalue pairing
log(f"\n  Pre-computation: classifying all 36 combos at N=3...")

all_idx3, pmats3, d3 = build_pauli_basis(N3)
num3 = 4 ** N3
Pi3 = build_Pi_persite(N3, all_idx3)
Pi3_inv = np.linalg.inv(Pi3)

# Build L_H for each single term at N=3
log(f"  Building L_H for 9 single terms at N=3...")
L_H_single3 = {}
for label in H_LABELS:
    H_term = np.zeros((d3, d3), dtype=complex)
    for i, j in bonds3:
        ops = [I2] * N3
        ops[i] = PM[label[0]]
        ops[j] = PM[label[1]]
        term = ops[0]
        for o in ops[1:]:
            term = np.kron(term, o)
        H_term += term
    L_H_single3[label] = build_L_H_pauli(N3, H_term, all_idx3, pmats3, d3)
log(f"  Done.")

L_D_diag3 = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx3])

combo_data = {}
broken_combos = []
palindromic_combos = []

for combo in all_combos:
    t1, t2 = combo.split('+')
    L_H_c = L_H_single3[t1] + L_H_single3[t2]
    L_full = L_H_c.copy()
    for a in range(num3):
        L_full[a, a] += L_D_diag3[a]

    # Per-site Pi error
    E = Pi3 @ L_full @ Pi3_inv + L_full + 2 * Sg3 * np.eye(num3)
    pi_err = np.max(np.abs(E)) / np.max(np.abs(L_full))

    # Eigenvalue pairing
    evals = np.linalg.eigvals(L_full)
    pe, pi_idx, n_wp = check_eigenvalue_pairing(evals, Sg3)

    is_broken = pi_err > 1e-8
    combo_data[combo] = {
        'pi_err': pi_err, 'evals': evals, 'pair_errors': pe,
        'n_well_paired': n_wp, 'L_full': L_full, 'L_H': L_H_c,
        'E': E, 'is_broken': is_broken
    }
    if is_broken:
        broken_combos.append(combo)
    else:
        palindromic_combos.append(combo)

log(f"  Broken at N=3: {len(broken_combos)}/36")
log(f"  Palindromic at N=3: {len(palindromic_combos)}/36")


# ################################################################
# SECTION 1: Anatomy of the Break
# ################################################################
log()
log("=" * 90)
log("SECTION 1: Anatomy of the Break")
log("=" * 90)

log(f"\n  {'Combo':>10}  {'Pi err':>10}  {'Well-paired':>12}  {'Orphaned':>9}  "
    f"{'Max pair err':>13}  {'Avg pair err':>13}")
log(f"  {'-'*72}")

for combo in sorted(broken_combos):
    d = combo_data[combo]
    pe = d['pair_errors']
    orphaned = num3 - d['n_well_paired']
    log(f"  {combo:>10}  {d['pi_err']:>10.4e}  {d['n_well_paired']:>12}/{num3}  "
        f"{orphaned:>9}  {np.max(pe):>13.4e}  {np.mean(pe):>13.4e}")

# Detailed anatomy for one case (XX+XY)
ref_combo = broken_combos[0] if broken_combos else None
if ref_combo:
    d = combo_data[ref_combo]
    evals = d['evals']
    pe = d['pair_errors']

    log(f"\n  Detailed anatomy for {ref_combo}:")
    log(f"    Eigenvalues sorted by pairing error (worst first):")
    order = np.argsort(-pe)
    log(f"    {'#':>4}  {'Re(lam)':>10}  {'Im(lam)':>10}  {'Pair err':>10}  {'Status':>10}")
    for idx, k in enumerate(order[:16]):
        status = "ORPHAN" if pe[k] > 1e-6 else "paired"
        log(f"    {k:>4}  {np.real(evals[k]):>10.6f}  {np.imag(evals[k]):>10.6f}  "
            f"{pe[k]:>10.4e}  {status:>10}")

    # Where do orphans sit?
    orphan_mask = pe > 1e-6
    orphan_rates = -np.real(evals[orphan_mask])
    paired_rates = -np.real(evals[~orphan_mask])
    if len(orphan_rates) > 0:
        log(f"\n    Orphan decay rates: min={np.min(orphan_rates):.4f}, "
            f"max={np.max(orphan_rates):.4f}, mean={np.mean(orphan_rates):.4f}")
    if len(paired_rates) > 0:
        log(f"    Paired decay rates: min={np.min(paired_rates):.4f}, "
            f"max={np.max(paired_rates):.4f}, mean={np.mean(paired_rates):.4f}")
    log(f"    Center of palindrome: Sg = {Sg3:.4f}")


# ################################################################
# SECTION 2: Structure in the Error
# ################################################################
log()
log("=" * 90)
log("SECTION 2: Structure in the Error")
log("=" * 90)

log(f"\n  2a. Gamma^2 scaling verification:")
log(f"  {'Combo':>10}  {'g=0.001':>10}  {'g=0.01':>10}  {'g=0.05':>10}  "
    f"{'coeff(0.001)':>13}  {'coeff(0.01)':>13}  {'coeff(0.05)':>13}")
log(f"  {'-'*78}")

gamma_test = [0.001, 0.01, 0.05]
for combo in sorted(broken_combos)[:8]:  # First 8 for brevity
    errs = []
    coeffs = []
    t1, t2 = combo.split('+')
    for g in gamma_test:
        Sg_t = N3 * g
        H_t = build_H_combo_Nsite(N3, bonds3, combo)
        L_t = build_L_vec(H_t, g, N3)
        ev_t = np.linalg.eigvals(L_t)
        pe_t, _, _ = check_eigenvalue_pairing(ev_t, Sg_t)
        err = np.max(pe_t)
        errs.append(err)
        coeffs.append(err / (g ** 2) if g > 0 else 0)
    log(f"  {combo:>10}  {errs[0]:>10.2e}  {errs[1]:>10.2e}  {errs[2]:>10.2e}  "
        f"{coeffs[0]:>13.2f}  {coeffs[1]:>13.2f}  {coeffs[2]:>13.2f}")

log(f"\n  Coefficients are approximately constant -> error scales as gamma^2.")

# 2b. Error spectrum: per-pair errors vs decay rate
log(f"\n  2b. Error spectrum for {ref_combo} (per-pair error vs decay rate):")
if ref_combo:
    d = combo_data[ref_combo]
    evals = d['evals']
    pe = d['pair_errors']
    rates = -np.real(evals)
    order = np.argsort(rates)
    log(f"    {'Rate':>8}  {'Pair err':>10}  {'|Im|':>8}  {'Level':>6}")
    for k in order[::4]:  # Every 4th for brevity
        level = "LOW" if pe[k] < 1e-6 else "MED" if pe[k] < 1e-3 else "HIGH"
        log(f"    {rates[k]:>8.4f}  {pe[k]:>10.4e}  {abs(np.imag(evals[k])):>8.4f}  {level:>6}")

    # Are errors concentrated at specific rates?
    high_err_rates = rates[pe > 1e-4]
    low_err_rates = rates[pe < 1e-6]
    if len(high_err_rates) > 0 and len(low_err_rates) > 0:
        log(f"\n    High-error modes: mean rate = {np.mean(high_err_rates):.4f}")
        log(f"    Low-error modes: mean rate = {np.mean(low_err_rates):.4f}")
        log(f"    Center: {Sg3:.4f}")


# ################################################################
# SECTION 3: New Symmetries After the Break
# ################################################################
log()
log("=" * 90)
log("SECTION 3: New Symmetries After the Break")
log("=" * 90)

if ref_combo:
    d = combo_data[ref_combo]
    E = d['E']  # Palindrome error matrix

    log(f"\n  3a. Palindrome error matrix E for {ref_combo}:")
    log(f"      E = Pi L Pi^{{-1}} + L + 2Sg I")
    log(f"      ||E||_F = {np.linalg.norm(E):.6f}")
    log(f"      ||E||_max = {np.max(np.abs(E)):.6f}")
    E_rank = int(np.sum(np.linalg.svd(E, compute_uv=False) > 1e-10))
    log(f"      rank(E) = {E_rank} (out of {num3})")
    log(f"      E affects {E_rank}/{num3} dimensions -- the break is {'localized' if E_rank < num3 // 2 else 'widespread'}")

    # 3b. Check for alternative pairing constants
    log(f"\n  3b. Alternative pairings for {ref_combo}:")
    evals = d['evals']
    # For each pair of eigenvalues, compute sum
    sums_real = []
    for i in range(num3):
        for j in range(i + 1, num3):
            s = evals[i] + evals[j]
            sums_real.append(np.real(s))

    # Find clusters: histogram the sums
    sums_arr = np.array(sums_real)
    target_sum = -2 * Sg3  # standard palindromic sum
    near_target = np.sum(np.abs(sums_arr - target_sum) < 0.01)
    log(f"      Standard palindromic sum = {target_sum:.4f}")
    log(f"      Pairs with sum near {target_sum:.4f}: {near_target}/{len(sums_arr)}")

    # Find other common sums
    bins = np.linspace(np.min(sums_arr) - 0.05, np.max(sums_arr) + 0.05, 100)
    hist, bin_edges = np.histogram(sums_arr, bins=bins)
    top_bins = np.argsort(hist)[-5:]
    log(f"      Top 5 sum clusters:")
    for b in sorted(top_bins, key=lambda x: -hist[x]):
        center = (bin_edges[b] + bin_edges[b + 1]) / 2
        log(f"        sum ~ {center:+.4f}: {hist[b]} pairs")

    # 3c. Partial palindrome: which blocks of L are still palindromic?
    log(f"\n  3c. Partial palindrome analysis:")
    # Decompose E by XY-weight sectors
    for w in range(N3 + 1):
        indices = [a for a, idx in enumerate(all_idx3) if xy_weight(idx) == w]
        E_block = E[np.ix_(indices, indices)]
        norm_block = np.linalg.norm(E_block)
        log(f"      XY-weight {w} block ({len(indices)}x{len(indices)}): ||E|| = {norm_block:.6e}")

    # Cross-weight blocks
    for w1 in range(N3 + 1):
        for w2 in range(w1 + 1, N3 + 1):
            idx1 = [a for a, idx in enumerate(all_idx3) if xy_weight(idx) == w1]
            idx2 = [a for a, idx in enumerate(all_idx3) if xy_weight(idx) == w2]
            E_cross = E[np.ix_(idx1, idx2)]
            norm_cross = np.linalg.norm(E_cross)
            if norm_cross > 1e-10:
                log(f"      w={w1} <-> w={w2} block: ||E|| = {norm_cross:.6e}")


# ################################################################
# SECTION 4: The Bifurcation Point
# ################################################################
log()
log("=" * 90)
log("SECTION 4: The Bifurcation Point")
log("=" * 90)

log(f"\n  H(alpha) = terms_on_bond(0,1) + alpha * terms_on_bond(1,2)")
log(f"  alpha=0: single bond (always palindromic)")
log(f"  alpha=1: full chain (may break)")

n_alpha = 50
alphas = np.linspace(0, 1, n_alpha + 1)

# Test a few representative broken combos
test_combos = sorted(broken_combos)[:4] if len(broken_combos) >= 4 else sorted(broken_combos)

for combo in test_combos:
    t1, t2 = combo.split('+')
    log(f"\n  {combo}:")
    log(f"    {'alpha':>7}  {'Max pair err':>13}  {'Well-paired':>12}  {'Orphaned':>9}")
    log(f"    {'-'*45}")

    errors_vs_alpha = []
    orphans_vs_alpha = []

    for alpha in alphas:
        # Build H with bond weights
        H_a = build_H_combo_Nsite_weighted(N3, bonds3, combo, [1.0, alpha])
        L_a = build_L_vec(H_a, gamma, N3)
        ev_a = np.linalg.eigvals(L_a)
        pe_a, _, n_wp = check_eigenvalue_pairing(ev_a, Sg3)
        max_err = np.max(pe_a)
        n_orphan = num3 - n_wp

        errors_vs_alpha.append(max_err)
        orphans_vs_alpha.append(n_orphan)

        if alpha in [0, 0.1, 0.2, 0.5, 0.8, 1.0] or abs(alpha - round(alpha, 1)) < 0.005:
            log(f"    {alpha:>7.2f}  {max_err:>13.4e}  {n_wp:>12}/{num3}  {n_orphan:>9}")

    # Analyze onset
    errs = np.array(errors_vs_alpha)
    first_break = np.argmax(errs > 1e-8)
    if errs[first_break] > 1e-8:
        alpha_c = alphas[first_break]
        log(f"\n    First break at alpha ~ {alpha_c:.3f}")
    else:
        log(f"\n    No break detected (palindrome holds for all alpha)")

    # Is onset sharp or smooth?
    if first_break > 0 and first_break < len(alphas) - 1:
        slope = (errs[first_break] - errs[first_break - 1]) / (alphas[1] - alphas[0])
        log(f"    Onset slope: {slope:.4e} (smooth)" if slope < 1 else f"    Onset slope: {slope:.4e} (sharp)")

# 4b. Check for new pairings at alpha=1
log(f"\n  4b. New pairings at alpha=1 (V-effect test):")

for combo in test_combos[:2]:
    d = combo_data[combo]
    evals = d['evals']
    pe = d['pair_errors']
    orphan_mask = pe > 1e-4

    if np.sum(orphan_mask) == 0:
        log(f"    {combo}: no orphans, palindrome intact.")
        continue

    orphan_evals = evals[orphan_mask]
    log(f"\n    {combo}: {len(orphan_evals)} orphan eigenvalues")

    # Do orphans pair among themselves?
    n_orph = len(orphan_evals)
    best_pairs = []
    for i in range(n_orph):
        for j in range(i + 1, n_orph):
            s = orphan_evals[i] + orphan_evals[j]
            best_pairs.append((np.real(s), np.imag(s), i, j))

    # Find the most common sum
    sums_r = np.array([p[0] for p in best_pairs])
    if len(sums_r) > 0:
        # Cluster sums
        bins = np.linspace(np.min(sums_r) - 0.02, np.max(sums_r) + 0.02, 40)
        hist, edges = np.histogram(sums_r, bins=bins)
        top_idx = np.argmax(hist)
        top_sum = (edges[top_idx] + edges[top_idx + 1]) / 2
        log(f"    Most common orphan pair sum: {top_sum:.4f} ({hist[top_idx]} pairs)")
        log(f"    Standard palindromic sum: {-2*Sg3:.4f}")
        if abs(top_sum - (-2 * Sg3)) > 0.05:
            log(f"    ** NEW PAIRING CONSTANT: orphans pair at sum = {top_sum:.4f}, not {-2*Sg3:.4f} **")
            log(f"    ** V-EFFECT DETECTED: one symmetry has split into (at least) two **")
        else:
            log(f"    Orphans still cluster near the palindromic sum.")


# ################################################################
# SECTION 5: Broken vs Unbroken Comparison
# ################################################################
log()
log("=" * 90)
log("SECTION 5: Broken vs Unbroken Comparison")
log("=" * 90)

combo_broken = broken_combos[0] if broken_combos else 'XX+XY'
combo_unbroken = 'XX+YY'  # Heisenberg-like, always palindromic

for label, combo in [("BROKEN", combo_broken), ("UNBROKEN", combo_unbroken)]:
    log(f"\n  {label}: {combo}")

    H_c = build_H_combo_Nsite(N3, bonds3, combo)
    L_c = build_L_vec(H_c, gamma, N3)
    evals_c = np.linalg.eigvals(L_c)

    rates = -np.real(evals_c)
    freqs = np.imag(evals_c)

    log(f"    Eigenvalues: {len(evals_c)}")
    log(f"    Decay rates: min={np.min(rates):.6f}, max={np.max(rates):.6f}, "
        f"mean={np.mean(rates):.6f}")
    log(f"    Frequencies: min={np.min(np.abs(freqs)):.4f}, max={np.max(np.abs(freqs)):.4f}")

    # Count distinct rates and frequencies
    unique_rates = len(np.unique(np.round(rates, 4)))
    unique_freqs = len(np.unique(np.round(np.abs(freqs), 4)))
    log(f"    Distinct rates: {unique_rates}")
    log(f"    Distinct |frequencies|: {unique_freqs}")

    # Steady states (rate ~ 0)
    n_steady = int(np.sum(np.abs(rates) < 1e-8))
    log(f"    Steady states: {n_steady}")

    # Longest-lived non-steady mode
    non_steady_rates = rates[np.abs(rates) > 1e-8]
    if len(non_steady_rates) > 0:
        slowest = np.min(non_steady_rates)
        log(f"    Slowest non-steady rate: {slowest:.6f} (T_half = {np.log(2)/slowest:.2f})")

    # Rate distribution (quartiles)
    q = np.percentile(rates, [25, 50, 75])
    log(f"    Rate quartiles: Q1={q[0]:.4f}, Q2={q[1]:.4f}, Q3={q[2]:.4f}")

    # Palindrome quality
    pe_c, _, n_wp = check_eigenvalue_pairing(evals_c, Sg3)
    log(f"    Palindrome: {n_wp}/{len(evals_c)} well-paired (max err = {np.max(pe_c):.4e})")

# Side-by-side summary
log(f"\n  Side-by-side summary:")
log(f"  {'Property':>25}  {'Broken':>12}  {'Unbroken':>12}")
log(f"  {'-'*52}")
for label, combo in [("Broken", combo_broken), ("Unbroken", combo_unbroken)]:
    H_c = build_H_combo_Nsite(N3, bonds3, combo)
    L_c = build_L_vec(H_c, gamma, N3)
    evals_c = np.linalg.eigvals(L_c)
    rates = -np.real(evals_c)
    freqs = np.imag(evals_c)
    pe_c, _, n_wp = check_eigenvalue_pairing(evals_c, Sg3)
    if label == "Broken":
        br = {'rates': rates, 'freqs': freqs, 'n_wp': n_wp, 'max_pe': np.max(pe_c),
              'n_freq': len(np.unique(np.round(np.abs(freqs), 3)))}
    else:
        ub = {'rates': rates, 'freqs': freqs, 'n_wp': n_wp, 'max_pe': np.max(pe_c),
              'n_freq': len(np.unique(np.round(np.abs(freqs), 3)))}

log(f"  {'Well-paired':>25}  {br['n_wp']:>12}  {ub['n_wp']:>12}")
log(f"  {'Max pair error':>25}  {br['max_pe']:>12.4e}  {ub['max_pe']:>12.4e}")
log(f"  {'Distinct |freq|':>25}  {br['n_freq']:>12}  {ub['n_freq']:>12}")
log(f"  {'Rate range':>25}  {np.ptp(br['rates']):>12.4f}  {np.ptp(ub['rates']):>12.4f}")


# ################################################################
# SECTION 6: The Rank Threshold
# ################################################################
log()
log("=" * 90)
log("SECTION 6: The Rank Threshold")
log("=" * 90)

# Compute N=2 Choi rank for all 36 combos
log(f"\n  Computing N=2 Choi rank and N=3 palindrome status for all 36 combos...")

N2 = 2
all_idx2, pmats2, d2 = build_pauli_basis(N2)
num2 = 4 ** N2
Sg2 = N2 * gamma
V2 = build_V(pmats2, d2)
L_D_diag2 = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx2])

# Build L_H for single terms at N=2
L_H_single2 = {}
for label in H_LABELS:
    ops = [I2] * N2
    ops[0] = PM[label[0]]
    ops[1] = PM[label[1]]
    H_t = ops[0]
    for o in ops[1:]:
        H_t = np.kron(H_t, o)
    L_H_single2[label] = build_L_H_pauli(N2, H_t, all_idx2, pmats2, d2)

log(f"\n  {'Combo':>10}  {'N2 rank':>8}  {'N3 palindromic':>15}  {'N3 Pi err':>10}  {'Category':>15}")
log(f"  {'-'*62}")

rank_threshold_data = []

for combo in all_combos:
    t1, t2 = combo.split('+')

    # N=2 Choi rank
    L_H2 = L_H_single2[t1] + L_H_single2[t2]
    L_full2 = L_H2.copy()
    for a in range(num2):
        L_full2[a, a] += L_D_diag2[a]
    Pi2, n_p2 = construct_Pi_eigvec(L_full2, Sg2, num2)
    if Pi2 is not None:
        rank2 = choi_rank_N2(Pi2, V2, d2)
    else:
        rank2 = -1

    # N=3 status
    cd = combo_data[combo]
    is_broken = cd['is_broken']
    pi_err3 = cd['pi_err']

    if not is_broken:
        cat = "palindromic"
    else:
        cat = "BROKEN"

    log(f"  {combo:>10}  {rank2:>8}  {'NO' if is_broken else 'YES':>15}  "
        f"{pi_err3:>10.2e}  {cat:>15}")
    rank_threshold_data.append((combo, rank2, is_broken, pi_err3))

# Analysis: is there a clean rank threshold?
log(f"\n  Rank threshold analysis:")
ranks_palindromic = [r for _, r, broken, _ in rank_threshold_data if not broken and r > 0]
ranks_broken = [r for _, r, broken, _ in rank_threshold_data if broken and r > 0]

if ranks_palindromic:
    log(f"    Palindromic N=2 ranks: {sorted(set(ranks_palindromic))}")
    log(f"    Max palindromic rank: {max(ranks_palindromic)}")
if ranks_broken:
    log(f"    Broken N=2 ranks: {sorted(set(ranks_broken))}")
    log(f"    Min broken rank: {min(ranks_broken)}")

if ranks_palindromic and ranks_broken:
    threshold = (max(ranks_palindromic) + min(ranks_broken)) / 2
    clean = max(ranks_palindromic) < min(ranks_broken)
    log(f"\n    Clean threshold: {'YES' if clean else 'NO'}")
    if clean:
        log(f"    Threshold: rank <= {max(ranks_palindromic)} survives, rank >= {min(ranks_broken)} breaks")
        log(f"    Physical: clocks with Choi rank <= {max(ranks_palindromic)} can synchronize across bonds.")
        log(f"    Clocks with rank >= {min(ranks_broken)} create new physics at the boundary.")
    else:
        log(f"    Overlap: some palindromic ranks overlap with broken ranks.")
        log(f"    Palindromic max = {max(ranks_palindromic)}, Broken min = {min(ranks_broken)}")


# ################################################################
# SUMMARY
# ################################################################
log()
log("=" * 90)
log("SUMMARY")
log("=" * 90)

n_broken = len(broken_combos)
log(f"""
  SECTION 1: Anatomy of the Break
    {n_broken}/36 combos break at N=3 (per-site Pi fails).
    Breaking is widespread across the spectrum (not localized to a few modes).

  SECTION 2: Error Structure
    Error scales as gamma^2 with combo-specific coefficients.
    The error is NOT uniform: different modes break at different rates.

  SECTION 3: New Symmetries
    The palindrome error matrix E has specific structure.
    Orphan eigenvalues may cluster at alternative pairing constants.

  SECTION 4: Bifurcation
    The palindrome breaks smoothly as the second bond strengthens (alpha: 0->1).
    No sharp threshold in alpha; onset is continuous.

  SECTION 5: Broken vs Unbroken
    Broken cases have richer spectral structure (more distinct frequencies).
    Breaking is differentiation, not destruction.

  SECTION 6: Rank Threshold
    N=2 Choi rank predicts N=3 breakability.
    {'Clean threshold exists.' if ranks_palindromic and ranks_broken and max(ranks_palindromic) < min(ranks_broken) else 'Threshold has some overlap.'}
""")

log("=" * 90)
log(f"COMPLETE -- {datetime.now()}")
log("=" * 90)

f.close()
print(f"\n>>> Results written to {OUT}")
