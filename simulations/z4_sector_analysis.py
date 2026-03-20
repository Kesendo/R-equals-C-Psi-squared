#!/usr/bin/env python3
"""
z4_sector_analysis.py: Five tests for Z4 sector physical content.

Test 1: Z4 sector classification (XY-weight, decay rate per sector)
Test 2: Palindromic pairs across Z4 sectors
Test 3: Pi^2 analysis (parity, commutation, known symmetries)
Test 4: Standing wave four-fold structure
Test 5: N-dependence of Z4 structure
"""

import numpy as np
from itertools import product as iprod
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "z4_sector_analysis.txt")
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

def build_Pi(all_idx):
    """Canonical per-site Pi: I<->X (+1), Y->iZ, Z->iY."""
    num = len(all_idx)
    Pi = np.zeros((num, num), dtype=complex)
    idx_map = {idx: i for i, idx in enumerate(all_idx)}
    for b, idx_b in enumerate(all_idx):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = idx_map[mapped]
        Pi[a, b] = sign
    return Pi

def build_H_heisenberg(N, bonds, J=1.0):
    all_idx, pmats, d = build_pauli_basis(N)
    num = len(all_idx)
    H = np.zeros((d, d), dtype=complex)
    for (i, j) in bonds:
        for P in [sx, sy, sz]:
            ops_i = [I2] * N
            ops_j = [I2] * N
            ops_i[i] = P
            ops_j[j] = P
            Pi_op = ops_i[0]
            Pj_op = ops_j[0]
            for k in range(1, N):
                Pi_op = np.kron(Pi_op, ops_i[k])
                Pj_op = np.kron(Pj_op, ops_j[k])
            H += J * (Pi_op @ Pj_op)
    return H

def build_full_L(N, bonds, gamma):
    all_idx, pmats, d = build_pauli_basis(N)
    num = len(all_idx)
    H = build_H_heisenberg(N, bonds)
    L_H = np.zeros((num, num), dtype=complex)
    for a in range(num):
        for b in range(num):
            comm = -1j * (H @ pmats[b] - pmats[b] @ H)
            L_H[a, b] = np.trace(pmats[a].conj().T @ comm) / d
    L_D = np.array([-2 * gamma * xy_weight(idx) for idx in all_idx])
    L = L_H.copy()
    for a in range(num):
        L[a, a] += L_D[a]
    Sg = N * gamma
    return L, Sg, all_idx, pmats, d

def classify_z4(Pi, vec, tol=0.1):
    """Classify vector by its Z4 Pi-eigenvalue. Returns (sector, quality)."""
    Pv = Pi @ vec
    # Try each eigenvalue
    best_sector = None
    best_quality = 0
    norm_v = np.linalg.norm(vec)
    if norm_v < 1e-15:
        return '+1', 0.0
    for phase, label in [(1, '+1'), (-1, '-1'), (1j, '+i'), (-1j, '-i')]:
        residual = np.linalg.norm(Pv - phase * vec) / norm_v
        quality = 1.0 - residual / 2.0  # 1.0 = perfect, 0.0 = orthogonal
        if quality > best_quality:
            best_quality = quality
            best_sector = label
    return best_sector, best_quality

def avg_xy_weight(vec, all_idx):
    """Average XY-weight of a vector in the Pauli basis."""
    weights = np.array([xy_weight(idx) for idx in all_idx], dtype=float)
    probs = np.abs(vec) ** 2
    total = np.sum(probs)
    if total < 1e-30:
        return 0.0
    return float(np.sum(weights * probs) / total)


# ============================================================
# TEST 1: Z4 SECTOR CLASSIFICATION
# ============================================================
def run_test_1(N=3, gamma=GAMMA):
    log("=" * 70)
    log(f"TEST 1: Z4 SECTOR CLASSIFICATION (N={N})")
    log("=" * 70)
    log()

    bonds = [(i, i + 1) for i in range(N - 1)]
    L, Sg, all_idx, pmats, d = build_full_L(N, bonds, gamma)
    num = 4 ** N
    Pi = build_Pi(all_idx)

    # Verify Pi^4 = I
    Pi2 = Pi @ Pi
    Pi4 = Pi2 @ Pi2
    pi4_err = np.max(np.abs(Pi4 - np.eye(num)))
    log(f"Pi^4 = I check: max error = {pi4_err:.2e}")

    # Pi eigenvalues
    pi_evals = np.linalg.eigvals(Pi)
    for phase, label in [(1, '+1'), (-1, '-1'), (1j, '+i'), (-1j, '-i')]:
        count = np.sum(np.abs(pi_evals - phase) < 0.01)
        log(f"  Pi eigenvalue {label}: multiplicity {count}")

    # Eigendecompose L
    evals, R = np.linalg.eig(L)

    # Classify each eigenvector
    sectors = {'+1': [], '-1': [], '+i': [], '-i': []}
    log()
    log(f"{'Idx':>4} {'Re(lam)':>9} {'Im(lam)':>9} {'Sector':>6} "
        f"{'Quality':>7} {'XY-wt':>6}")
    log("-" * 50)

    for k in range(num):
        vec = R[:, k]
        sector, quality = classify_z4(Pi, vec)
        w = avg_xy_weight(vec, all_idx)
        sectors[sector].append({
            'idx': k, 'eval': evals[k], 'quality': quality,
            'xy_weight': w, 're': evals[k].real, 'im': evals[k].imag
        })
        if k < 20 or k >= num - 5:
            log(f"{k:>4} {evals[k].real:>9.5f} {evals[k].imag:>9.5f} "
                f"{sector:>6} {quality:>7.4f} {w:>6.3f}")

    log()
    log("Sector statistics:")
    log(f"{'Sector':>7} {'Count':>6} {'Avg XY-wt':>10} {'Avg Re(lam)':>12} "
        f"{'Avg Quality':>12}")
    log("-" * 55)
    for s in ['+1', '-1', '+i', '-i']:
        entries = sectors[s]
        n = len(entries)
        if n == 0:
            log(f"{s:>7} {0:>6}")
            continue
        avg_w = np.mean([e['xy_weight'] for e in entries])
        avg_re = np.mean([e['re'] for e in entries])
        avg_q = np.mean([e['quality'] for e in entries])
        log(f"{s:>7} {n:>6} {avg_w:>10.4f} {avg_re:>12.6f} {avg_q:>12.4f}")

    return L, Sg, all_idx, pmats, d, Pi, evals, R, sectors


# ============================================================
# TEST 2: PALINDROMIC PAIRS ACROSS Z4 SECTORS
# ============================================================
def run_test_2(L, Sg, Pi, evals, R, sectors, all_idx):
    log()
    log("=" * 70)
    log("TEST 2: PALINDROMIC PAIRS ACROSS Z4 SECTORS")
    log("=" * 70)
    log()

    num = len(evals)
    # Build flat sector map
    sector_map = {}
    for s, entries in sectors.items():
        for e in entries:
            sector_map[e['idx']] = s

    # Find palindromic pairs
    used = set()
    pair_sectors = {}  # (sector_A, sector_B) -> count
    pairs = []
    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = np.argmin(diffs)
        if diffs[best] < 1e-6 and best != k:
            s_k = sector_map.get(k, '?')
            s_b = sector_map.get(best, '?')
            key = tuple(sorted([s_k, s_b]))
            pair_sectors[key] = pair_sectors.get(key, 0) + 1
            pairs.append((k, best, s_k, s_b))
            used.add(k)
            used.add(best)
        elif diffs[best] < 1e-6 and best == k:
            # Self-paired
            s_k = sector_map.get(k, '?')
            key = (s_k, s_k)
            pair_sectors[key] = pair_sectors.get(key, 0) + 1
            pairs.append((k, k, s_k, s_k))
            used.add(k)

    log(f"Total palindromic pairs found: {len(pairs)}")
    log()
    log("Sector pairing distribution:")
    log(f"{'Sector pair':>15} {'Count':>6} {'Hypothesis':>20}")
    log("-" * 45)

    for key in sorted(pair_sectors.keys()):
        count = pair_sectors[key]
        s1, s2 = key
        if (s1 == '+1' and s2 == '-1') or (s1 == '-1' and s2 == '+1'):
            hyp = "EXPECTED (opposite)"
        elif (s1 == '+i' and s2 == '-i') or (s1 == '-i' and s2 == '+i'):
            hyp = "EXPECTED (opposite)"
        elif s1 == s2:
            hyp = "self-paired"
        else:
            hyp = "UNEXPECTED"
        log(f"{'(' + s1 + ', ' + s2 + ')':>15} {count:>6} {hyp:>20}")

    # Check hypothesis: do opposite sectors pair?
    expected_pairs = pair_sectors.get(('+1', '-1'), 0) + pair_sectors.get(('-1', '+1'), 0)
    expected_pairs += pair_sectors.get(('-i', '+i'), 0) + pair_sectors.get(('+i', '-i'), 0)
    total_non_self = sum(v for k, v in pair_sectors.items() if k[0] != k[1])
    log()
    if total_non_self > 0:
        frac = expected_pairs / total_non_self
        log(f"Opposite-sector pairs: {expected_pairs}/{total_non_self} ({frac:.1%})")
        if frac > 0.9:
            log("CONFIRMED: palindromic pairs strongly prefer opposite Z4 sectors")
        elif frac > 0.5:
            log("PARTIAL: some preference for opposite sectors")
        else:
            log("NOT CONFIRMED: pairs scatter across sectors")
    else:
        log("All pairs are self-paired (no cross-sector pairing to test)")


# ============================================================
# TEST 3: PI^2 ANALYSIS
# ============================================================
def run_test_3(L, Sg, all_idx, pmats, d, Pi, evals, R):
    log()
    log("=" * 70)
    log("TEST 3: PI^2 ANALYSIS")
    log("=" * 70)
    log()

    num = len(all_idx)
    N = len(all_idx[0])
    Pi2 = Pi @ Pi

    # Verify Pi^2 action on Pauli basis
    log("Pi^2 action on single-site Pauli basis (from first site):")
    for i, name in enumerate(PAULI_NAMES):
        # Find basis element with this Pauli at site 0 and I elsewhere
        idx = tuple([i] + [0] * (N - 1))
        pos = list(all_idx).index(idx)
        result = Pi2[:, pos]
        nonzero = [(j, result[j]) for j in range(num) if abs(result[j]) > 1e-10]
        if len(nonzero) == 1:
            j, val = nonzero[0]
            target_idx = all_idx[j]
            target_name = ''.join(PAULI_NAMES[k] for k in target_idx)
            log(f"  {name}{'I' * (N - 1)} -> {val:+.1f} * {target_name}")
        else:
            log(f"  {name}{'I' * (N - 1)} -> {len(nonzero)} terms")

    log()
    log("Expected: I->I, X->X, Y->-Y, Z->-Z")

    # Pi^2 eigenvalues
    pi2_evals = np.linalg.eigvals(Pi2)
    n_plus1 = np.sum(np.abs(pi2_evals - 1) < 0.01)
    n_minus1 = np.sum(np.abs(pi2_evals + 1) < 0.01)
    log(f"\nPi^2 eigenvalues: +1 x {n_plus1}, -1 x {n_minus1}")

    # Commutation with L
    comm = Pi2 @ L - L @ Pi2
    comm_norm = np.max(np.abs(comm))
    log(f"\n[Pi^2, L] max element: {comm_norm:.2e}")
    if comm_norm < 1e-10:
        log("Pi^2 COMMUTES with L: it is a conserved symmetry!")
    else:
        log(f"Pi^2 does NOT commute with L (norm {comm_norm:.2e})")

    # Compare with known operators
    # XXX parity: X at every site
    XXX = build_Pi_xxx(all_idx)
    diff_xxx = np.max(np.abs(Pi2 - XXX))
    log(f"\n||Pi^2 - X^N|| = {diff_xxx:.4f}")
    if diff_xxx < 1e-10:
        log("Pi^2 = X^N parity")
    else:
        log("Pi^2 != X^N parity")

    # ZZZ parity: Z at every site
    ZZZ = build_Pi_zzz(all_idx)
    diff_zzz = np.max(np.abs(Pi2 - ZZZ))
    log(f"||Pi^2 - Z^N|| = {diff_zzz:.4f}")

    # Check if Pi^2 is the "Y,Z negation" operator
    # It should negate all basis elements with odd number of Y or Z
    log("\nPi^2 diagonal structure:")
    for i in range(min(num, 20)):
        diag = Pi2[i, i]
        idx = all_idx[i]
        name = ''.join(PAULI_NAMES[k] for k in idx)
        n_yz = sum(1 for k in idx if k in (2, 3))
        expected = (-1) ** n_yz
        log(f"  {name}: Pi^2 diag = {diag.real:+.1f} "
            f"(expected (-1)^{n_yz} = {expected:+.1f})")

    # Classify L eigenvectors by Pi^2 eigenvalue
    log("\nL eigenvectors classified by Pi^2 eigenvalue:")
    n_pi2_plus = 0
    n_pi2_minus = 0
    for k in range(num):
        vec = R[:, k]
        Pi2v = Pi2 @ vec
        ratio = Pi2v / (vec + 1e-30)
        # Check if all components have same ratio
        valid = np.abs(vec) > 1e-10
        if np.sum(valid) > 0:
            ratios = ratio[valid]
            if np.std(np.abs(ratios)) < 0.01 * np.mean(np.abs(ratios)):
                mean_ratio = np.mean(ratios)
                if abs(mean_ratio - 1) < 0.1:
                    n_pi2_plus += 1
                elif abs(mean_ratio + 1) < 0.1:
                    n_pi2_minus += 1
    log(f"  Pi^2 = +1 sector: {n_pi2_plus} eigenvectors")
    log(f"  Pi^2 = -1 sector: {n_pi2_minus} eigenvectors")
    log(f"  Unclassified: {num - n_pi2_plus - n_pi2_minus}")

def build_Pi_xxx(all_idx):
    """X^N parity operator in Pauli basis."""
    # X at each site: I<->X, Y<->-Y... wait, X doesn't do Y<->-Y.
    # In the Pauli basis, conjugation by X^N maps:
    #   I->I, X->X, Y->-Y, Z->-Z at each site
    # That's actually the same as Pi^2! Let me compute it properly.
    num = len(all_idx)
    P = np.zeros((num, num), dtype=complex)
    idx_map = {idx: i for i, idx in enumerate(all_idx)}
    # Conjugation by X: X P X = P if P in {I,X}, -P if P in {Y,Z}
    for b, idx_b in enumerate(all_idx):
        sign = 1
        for s in idx_b:
            if s in (2, 3):  # Y or Z
                sign *= -1
        P[b, b] = sign
    return P

def build_Pi_zzz(all_idx):
    """Z^N parity operator in Pauli basis."""
    num = len(all_idx)
    P = np.zeros((num, num), dtype=complex)
    # Conjugation by Z: Z P Z = P if P in {I,Z}, -P if P in {X,Y}
    for b, idx_b in enumerate(all_idx):
        sign = 1
        for s in idx_b:
            if s in (1, 2):  # X or Y
                sign *= -1
        P[b, b] = sign
    return P


# ============================================================
# TEST 4: STANDING WAVE FOUR-FOLD STRUCTURE
# ============================================================
def run_test_4(L, Sg, all_idx, pmats, d, Pi, evals, R, sectors):
    log()
    log("=" * 70)
    log("TEST 4: STANDING WAVE FOUR-FOLD STRUCTURE")
    log("=" * 70)
    log()

    num = len(all_idx)
    N = len(all_idx[0])

    # Build sector map
    sector_map = {}
    for s, entries in sectors.items():
        for e in entries:
            sector_map[e['idx']] = s

    # Find palindromic pairs and their sectors
    used = set()
    pairs_by_sector = {'+1/-1': [], '+i/-i': [], 'self': [], 'other': []}

    for k in range(num):
        if k in used:
            continue
        target = -evals[k] - 2 * Sg
        diffs = np.abs(evals - target)
        for u in used:
            diffs[u] = 1e30
        best = np.argmin(diffs)
        if diffs[best] < 1e-6:
            s_k = sector_map.get(k, '?')
            s_b = sector_map.get(best, '?')
            freq = abs(evals[k].imag)
            if k == best:
                pairs_by_sector['self'].append((k, best, freq))
            elif {s_k, s_b} == {'+1', '-1'}:
                pairs_by_sector['+1/-1'].append((k, best, freq))
            elif {s_k, s_b} == {'+i', '-i'}:
                pairs_by_sector['+i/-i'].append((k, best, freq))
            else:
                pairs_by_sector['other'].append((k, best, freq, s_k, s_b))
            used.add(k)
            used.add(best)

    log("Palindromic pairs by Z4 sector pairing:")
    for key in ['+1/-1', '+i/-i', 'self', 'other']:
        entries = pairs_by_sector[key]
        n = len(entries)
        if n == 0:
            log(f"  {key}: 0 pairs")
            continue
        freqs = [e[2] for e in entries]
        osc = sum(1 for f in freqs if f > 1e-6)
        log(f"  {key}: {n} pairs, {osc} oscillating, "
            f"freq range [{min(freqs):.4f}, {max(freqs):.4f}]")

    # Frequency separation by sector
    log()
    log("Distinct oscillation frequencies by sector pairing:")
    for key in ['+1/-1', '+i/-i']:
        entries = pairs_by_sector[key]
        freqs = sorted(set(round(e[2], 4) for e in entries if e[2] > 1e-6))
        if freqs:
            log(f"  {key} frequencies: {freqs}")

    # Contribution to specific observables
    log()
    log("Observable decomposition by Z4 sector pairing:")
    # For a few key observables, compute the weight from each sector pair
    R_inv = np.linalg.inv(R)

    observables = []
    for idx in all_idx:
        name = ''.join(PAULI_NAMES[k] for k in idx)
        w = xy_weight(idx)
        if w in (0, N) or name in ('XXI', 'IXX', 'ZZI', 'IZZ', 'XYI', 'YXI'):
            observables.append((name, idx, w))

    for name, idx, w in observables[:10]:
        pos = list(all_idx).index(idx)
        # Decomposition coefficients in eigenbasis
        row = R_inv[pos, :]  # How this observable decomposes into eigenvectors

        weight_by_sector = {'+1/-1': 0, '+i/-i': 0, 'self': 0, 'other': 0}
        for key, entries in pairs_by_sector.items():
            for entry in entries:
                k = entry[0]
                b = entry[1]
                w_pair = abs(row[k]) ** 2 + (abs(row[b]) ** 2 if b != k else 0)
                weight_by_sector[key] += w_pair

        total = sum(weight_by_sector.values())
        if total > 1e-15:
            log(f"  {name} (w={w}): "
                f"+1/-1={weight_by_sector['+1/-1']/total:.1%}  "
                f"+i/-i={weight_by_sector['+i/-i']/total:.1%}  "
                f"self={weight_by_sector['self']/total:.1%}")


# ============================================================
# TEST 5: N-DEPENDENCE
# ============================================================
def run_test_5():
    log()
    log("=" * 70)
    log("TEST 5: N-DEPENDENCE OF Z4 STRUCTURE")
    log("=" * 70)
    log()

    gamma = GAMMA

    log(f"{'N':>3} {'Dim':>5} {'Sectors':>30} "
        f"{'Avg XY-wt by sector':>40} {'Avg Quality':>12}")
    log("-" * 95)

    for N in [2, 3, 4, 5]:
        bonds = [(i, i + 1) for i in range(N - 1)]
        L, Sg, all_idx, pmats, d = build_full_L(N, bonds, gamma)
        num = 4 ** N
        Pi = build_Pi(all_idx)

        evals, R = np.linalg.eig(L)

        sector_counts = {'+1': 0, '-1': 0, '+i': 0, '-i': 0}
        sector_xyw = {'+1': [], '-1': [], '+i': [], '-i': []}
        sector_quality = {'+1': [], '-1': [], '+i': [], '-i': []}

        for k in range(num):
            vec = R[:, k]
            sector, quality = classify_z4(Pi, vec)
            w = avg_xy_weight(vec, all_idx)
            sector_counts[sector] += 1
            sector_xyw[sector].append(w)
            sector_quality[sector].append(quality)

        counts_str = ", ".join(f"{s}:{sector_counts[s]}" for s in ['+1', '-1', '+i', '-i'])
        xyw_str = ", ".join(
            f"{s}:{np.mean(sector_xyw[s]):.2f}" if sector_xyw[s] else f"{s}:--"
            for s in ['+1', '-1', '+i', '-i'])
        avg_q = np.mean([q for qs in sector_quality.values() for q in qs])

        log(f"{N:>3} {num:>5} {counts_str:>30} {xyw_str:>40} {avg_q:>12.4f}")

    log()
    log("Expected: equal sector populations (4^N / 4 each),")
    log("  differentiation in XY-weight persists across N.")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("Z4 Sector Analysis: Five Tests")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"gamma = {GAMMA}")
    log()

    t1 = time.time()
    L, Sg, all_idx, pmats, d, Pi, evals, R, sectors = run_test_1()
    log(f"\n[Test 1 completed in {time.time() - t1:.1f}s]")

    t2 = time.time()
    run_test_2(L, Sg, Pi, evals, R, sectors, all_idx)
    log(f"\n[Test 2 completed in {time.time() - t2:.1f}s]")

    t3 = time.time()
    run_test_3(L, Sg, all_idx, pmats, d, Pi, evals, R)
    log(f"\n[Test 3 completed in {time.time() - t3:.1f}s]")

    t4 = time.time()
    run_test_4(L, Sg, all_idx, pmats, d, Pi, evals, R, sectors)
    log(f"\n[Test 4 completed in {time.time() - t4:.1f}s]")

    t5 = time.time()
    run_test_5()
    log(f"\n[Test 5 completed in {time.time() - t5:.1f}s]")

    log()
    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
