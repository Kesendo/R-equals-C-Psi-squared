#!/usr/bin/env python3
"""
N=4 Across the Boundary: Two N=2 Systems Meeting
====================================================
Test 1: Baseline comparison (chain vs isolated vs bridge)
Test 3: Symmetry of the bridge (bidirectionality requirement)
Test 2: The 70/30 split in eigenvalue spectrum (kappa sweep)
Test 4: New modes created by the bridge
Test 5: Does the bridge have a Pi?

Script: simulations/n4_across_boundary.py
Output: simulations/results/n4_across_boundary.txt
"""

import numpy as np
from itertools import product as iprod
import os, sys, time

# ============================================================
# OUTPUT
# ============================================================
OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "n4_across_boundary.txt")
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
PNAMES = ['I', 'X', 'Y', 'Z']

N = 4
d = 2 ** N   # 16
d2 = d * d   # 256


def site_op(op, k, N=4):
    ops = [I2] * N; ops[k] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r


def build_H_heisenberg(bonds, J=1.0):
    H = np.zeros((d, d), dtype=complex)
    for i, j in bonds:
        for P in [sx, sy, sz]:
            H += J * site_op(P, i) @ site_op(P, j)
    return H


def build_L_H(H):
    """Hamiltonian superoperator in computational vectorized basis."""
    Id = np.eye(d)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def build_L_D_zdeph(gammas):
    """Z-dephasing dissipator for given per-site gammas."""
    L_D = np.zeros((d2, d2), dtype=complex)
    for k in range(N):
        Zk = site_op(sz, k)
        L_D += gammas[k] * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L_D


def add_lindblad_jump(L, F, kappa):
    """Add Lindblad term kappa*(F rho F^dag - 0.5{F^dag F, rho}) to L."""
    Id = np.eye(d)
    FdF = F.conj().T @ F
    L += kappa * (np.kron(F, F.conj())
                  - 0.5 * np.kron(FdF, Id)
                  - 0.5 * np.kron(Id, FdF.T))
    return L


def palindrome_err(evals, Sg, tol=1e-6):
    n = len(evals); max_err = 0; n_paired = 0
    for k in range(n):
        target = -(evals[k] + 2 * Sg)
        best = np.min(np.abs(evals - target))
        if best < tol: n_paired += 1
        if best > max_err: max_err = best
    return max_err, n_paired, n


def build_bridge_jumps(model, kappa_AB, kappa_BA):
    """Build cross-boundary jump operators and their kappas.
    Boundary qubits: 1 (from pair A) and 2 (from pair B).
    Returns list of (F, kappa) tuples.
    """
    jumps = []

    if model == 'xz_boundary':
        # Coherence on A boundary -> Population on B boundary
        # F = X_1 tensor Z_2  and  Y_1 tensor Z_2
        F1 = site_op(sx, 1) @ site_op(sz, 2)
        F2 = site_op(sy, 1) @ site_op(sz, 2)
        jumps.append((F1, kappa_AB))
        jumps.append((F2, kappa_AB))
        # Reverse: B -> A
        F3 = site_op(sz, 1) @ site_op(sx, 2)
        F4 = site_op(sz, 1) @ site_op(sy, 2)
        jumps.append((F3, kappa_BA))
        jumps.append((F4, kappa_BA))

    elif model == 'zz_collective':
        # Collective dephasing: Z_1 Z_2
        F = site_op(sz, 1) @ site_op(sz, 2)
        jumps.append((F, kappa_AB))  # symmetric by construction

    elif model == 'amplitude_damping':
        # |0><1| on A boundary, |1><0| on B boundary (and reverse)
        sp = (sx + 1j * sy) / 2  # sigma_plus = |1><0|
        sm = (sx - 1j * sy) / 2  # sigma_minus = |0><1|
        # A decays, B excites
        F1 = site_op(sm, 1) @ site_op(sp, 2)
        jumps.append((F1, kappa_AB))
        # B decays, A excites
        F2 = site_op(sp, 1) @ site_op(sm, 2)
        jumps.append((F2, kappa_BA))

    elif model == 'xx_yy_boundary':
        # Heisenberg-like dissipative: XX + YY across boundary
        F1 = site_op(sx, 1) @ site_op(sx, 2)
        F2 = site_op(sy, 1) @ site_op(sy, 2)
        jumps.append((F1, kappa_AB))
        jumps.append((F2, kappa_AB))

    elif model == 'pure_dephasing_cross':
        # Cross-boundary Z-dephasing: each side dephases the other
        # Z on A boundary (qubit 1)
        F1 = site_op(sz, 1)
        jumps.append((F1, kappa_AB))
        # Z on B boundary (qubit 2)
        F2 = site_op(sz, 2)
        jumps.append((F2, kappa_BA))

    return jumps


def build_system(system_type, gamma=0.05, J=1.0, kappa=0.05,
                 bridge_model='xz_boundary',
                 kappa_AB=None, kappa_BA=None):
    """Build Liouvillian for different system types."""
    if kappa_AB is None: kappa_AB = kappa
    if kappa_BA is None: kappa_BA = kappa

    gammas = [gamma] * N

    if system_type == 'chain':
        # Standard N=4 chain: 0-1-2-3
        bonds = [(0,1), (1,2), (2,3)]
        H = build_H_heisenberg(bonds, J)
        L = build_L_H(H) + build_L_D_zdeph(gammas)
        Sg = sum(gammas)
        return L, Sg

    elif system_type == 'isolated':
        # Two isolated pairs: (0,1) and (2,3)
        bonds_A = [(0,1)]
        bonds_B = [(2,3)]
        H = build_H_heisenberg(bonds_A, J) + build_H_heisenberg(bonds_B, J)
        L = build_L_H(H) + build_L_D_zdeph(gammas)
        Sg = sum(gammas)
        return L, Sg

    elif system_type == 'bridge':
        # Two pairs with dissipative bridge
        bonds_A = [(0,1)]
        bonds_B = [(2,3)]
        H = build_H_heisenberg(bonds_A, J) + build_H_heisenberg(bonds_B, J)
        L = build_L_H(H) + build_L_D_zdeph(gammas)
        # Add bridge coupling
        jumps = build_bridge_jumps(bridge_model, kappa_AB, kappa_BA)
        for F, k in jumps:
            if k > 0:
                L = add_lindblad_jump(L, F, k)
        # Sg includes bridge contributions
        # Each jump F with rate kappa contributes to Sg through the
        # anti-commutator term. The effective Sg depends on the model.
        # For eigenvalue pairing, compute empirically.
        Sg = sum(gammas)  # base contribution
        return L, Sg


# ============================================================
# TEST 1: BASELINE COMPARISON
# ============================================================
def run_test_1(gamma=0.05, J=1.0):
    log("=" * 70)
    log("TEST 1: BASELINE COMPARISON")
    log("=" * 70)
    log()
    log(f"N=4, gamma={gamma}, J={J}")
    log()

    systems = {
        'chain':    ('Standard N=4 chain (0-1-2-3)', 'chain'),
        'isolated': ('Two isolated N=2 pairs (0-1)(2-3)', 'isolated'),
    }

    # Add bridge with different coupling models
    bridge_models = ['xz_boundary', 'zz_collective', 'amplitude_damping',
                     'xx_yy_boundary', 'pure_dephasing_cross']

    results = {}

    for label, (desc, stype) in systems.items():
        L, Sg = build_system(stype, gamma, J)
        evals = np.linalg.eigvals(L)
        err, n_p, n_t = palindrome_err(evals, Sg)
        rates = -evals.real
        n_ss = int(np.sum(np.abs(rates) < 1e-8))
        unique_rates = len(np.unique(np.round(rates, 6)))

        results[label] = evals

        log(f"--- {desc} ---")
        log(f"Sg = {Sg:.4f}")
        log(f"Palindrome: {n_p}/{n_t}, err = {err:.2e}")
        log(f"Steady states: {n_ss}")
        log(f"Distinct rates: {unique_rates}")
        log(f"Rate range: [{np.min(rates):.6f}, {np.max(rates):.6f}]")
        log()

    # Bridge models
    log("--- Bridge coupling models (kappa = gamma = 0.05, symmetric) ---")
    log()

    log(f"{'Model':>22}  {'Paired':>8}  {'MaxErr':>10}  "
        f"{'Steady':>7}  {'Rates':>6}  {'Rate range':>20}")
    log("-" * 80)

    for model in bridge_models:
        L, Sg_base = build_system('bridge', gamma, J, kappa=gamma,
                                   bridge_model=model)
        evals = np.linalg.eigvals(L)

        # Try multiple Sg candidates
        best_err = 1e30; best_Sg = Sg_base; best_np = 0
        for Sg_try in np.linspace(Sg_base * 0.5, Sg_base * 3, 50):
            err_try, np_try, _ = palindrome_err(evals, Sg_try)
            if np_try > best_np or (np_try == best_np and err_try < best_err):
                best_err = err_try; best_Sg = Sg_try; best_np = np_try

        rates = -evals.real
        n_ss = int(np.sum(np.abs(rates) < 1e-8))
        unique_rates = len(np.unique(np.round(rates, 6)))

        results[f'bridge_{model}'] = evals

        log(f"{model:>22}  {best_np:>4}/{d2}  {best_err:>10.2e}  "
            f"{n_ss:>7}  {unique_rates:>6}  "
            f"[{np.min(rates):.4f}, {np.max(rates):.4f}]")

        if best_np == d2:
            log(f"  -> PALINDROMIC at Sg = {best_Sg:.4f}")
        elif best_np > d2 * 0.9:
            log(f"  -> Nearly palindromic ({best_np}/{d2}) at Sg = {best_Sg:.4f}")

    log()
    return results


# ============================================================
# TEST 3: SYMMETRY OF THE BRIDGE
# ============================================================
def run_test_3(gamma=0.05, J=1.0):
    log("=" * 70)
    log("TEST 3: SYMMETRY OF THE BRIDGE (BIDIRECTIONALITY)")
    log("=" * 70)
    log()

    configs = {
        'symmetric':    (0.05, 0.05),
        'asymmetric':   (0.07, 0.03),
        'one_way_AB':   (0.10, 0.00),
        'one_way_BA':   (0.00, 0.10),
        'strong_sym':   (0.20, 0.20),
        'weak_sym':     (0.01, 0.01),
    }

    models_to_test = ['zz_collective', 'pure_dephasing_cross',
                      'amplitude_damping', 'xz_boundary']

    for model in models_to_test:
        log(f"--- Bridge model: {model} ---")
        log()
        log(f"{'Config':>14}  {'kAB':>5}  {'kBA':>5}  "
            f"{'Paired':>8}  {'MaxErr':>10}  {'Best Sg':>8}  {'Steady':>7}")
        log("-" * 65)

        for cname, (kAB, kBA) in configs.items():
            L, Sg_base = build_system('bridge', gamma, J,
                                       bridge_model=model,
                                       kappa_AB=kAB, kappa_BA=kBA)
            evals = np.linalg.eigvals(L)

            # Search for best Sg
            best_err = 1e30; best_Sg = Sg_base; best_np = 0
            for Sg_try in np.linspace(0, Sg_base * 5, 100):
                err_try, np_try, _ = palindrome_err(evals, Sg_try)
                if np_try > best_np or (np_try == best_np and err_try < best_err):
                    best_err = err_try; best_Sg = Sg_try; best_np = np_try

            n_ss = int(np.sum(np.abs(evals.real) < 1e-8))

            pal_tag = "YES" if best_np == d2 else f"{best_np}/{d2}"
            log(f"{cname:>14}  {kAB:>5.2f}  {kBA:>5.2f}  "
                f"{pal_tag:>8}  {best_err:>10.2e}  {best_Sg:>8.4f}  {n_ss:>7}")

        log()

    # Summary
    log("--- Bidirectionality analysis ---")
    log()
    log("For each model, check: does asymmetry/one-way break the palindrome?")
    log()

    for model in models_to_test:
        sym_L, _ = build_system('bridge', gamma, J, bridge_model=model,
                                 kappa_AB=0.05, kappa_BA=0.05)
        ow_L, _ = build_system('bridge', gamma, J, bridge_model=model,
                                kappa_AB=0.10, kappa_BA=0.00)
        asym_L, _ = build_system('bridge', gamma, J, bridge_model=model,
                                  kappa_AB=0.07, kappa_BA=0.03)

        ev_sym = np.linalg.eigvals(sym_L)
        ev_ow = np.linalg.eigvals(ow_L)
        ev_asym = np.linalg.eigvals(asym_L)

        # Check palindrome with best Sg
        def best_palindrome(evals, Sg_base):
            best = (1e30, 0, Sg_base)
            for Sg in np.linspace(0, Sg_base * 5, 100):
                err, np_, _ = palindrome_err(evals, Sg)
                if np_ > best[1] or (np_ == best[1] and err < best[0]):
                    best = (err, np_, Sg)
            return best

        sym_pal = best_palindrome(ev_sym, N * gamma)
        ow_pal = best_palindrome(ev_ow, N * gamma)
        asym_pal = best_palindrome(ev_asym, N * gamma)

        sym_ok = sym_pal[1] == d2
        ow_ok = ow_pal[1] == d2
        asym_ok = asym_pal[1] == d2

        log(f"  {model}:")
        log(f"    Symmetric:   palindromic = {'YES' if sym_ok else 'NO'} "
            f"({sym_pal[1]}/{d2})")
        log(f"    One-way:     palindromic = {'YES' if ow_ok else 'NO'} "
            f"({ow_pal[1]}/{d2})")
        log(f"    Asymmetric:  palindromic = {'YES' if asym_ok else 'NO'} "
            f"({asym_pal[1]}/{d2})")

        if sym_ok and not ow_ok:
            log(f"    -> BIDIRECTIONALITY REQUIRED for {model}")
        elif sym_ok and ow_ok:
            log(f"    -> Palindrome survives one-way coupling")
        log()


# ============================================================
# TEST 2: THE 70/30 SPLIT (KAPPA SWEEP)
# ============================================================
def run_test_2(gamma=0.05, J=1.0):
    log("=" * 70)
    log("TEST 2: THE 70/30 SPLIT (KAPPA SWEEP)")
    log("=" * 70)
    log()

    # Baseline: isolated pairs
    L_iso, Sg_iso = build_system('isolated', gamma, J)
    evals_iso = np.sort(np.linalg.eigvals(L_iso).real)
    spectral_width = np.max(evals_iso) - np.min(evals_iso)

    kappas = np.concatenate([np.array([0]),
                             np.logspace(-4, np.log10(gamma), 15),
                             np.linspace(gamma, 5*gamma, 10)])
    kappas = np.unique(np.round(kappas, 6))

    model = 'zz_collective'  # use ZZ as representative
    log(f"Model: {model}")
    log(f"Spectral width at kappa=0: {spectral_width:.6f}")
    log()

    log(f"{'kappa':>10}  {'Shifted':>8}  {'Frac':>6}  "
        f"{'Paired':>8}  {'MaxErr':>10}  {'Best Sg':>8}")
    log("-" * 60)

    for kappa in kappas:
        L_br, Sg_base = build_system('bridge', gamma, J, kappa=kappa,
                                      bridge_model=model)
        evals_br = np.sort(np.linalg.eigvals(L_br).real)

        # Count shifted eigenvalues
        threshold = 0.01 * spectral_width
        n_shifted = int(np.sum(np.abs(evals_br - evals_iso) > threshold))
        frac = n_shifted / d2

        # Palindrome check
        evals_full = np.linalg.eigvals(L_br)
        best_err = 1e30; best_Sg = Sg_base; best_np = 0
        for Sg_try in np.linspace(0, Sg_base * 5, 80):
            err_try, np_try, _ = palindrome_err(evals_full, Sg_try)
            if np_try > best_np or (np_try == best_np and err_try < best_err):
                best_err = err_try; best_Sg = Sg_try; best_np = np_try

        log(f"{kappa:>10.5f}  {n_shifted:>4}/{d2}  {frac:>6.1%}  "
            f"{best_np:>4}/{d2}  {best_err:>10.2e}  {best_Sg:>8.4f}")

    log()


# ============================================================
# TEST 4: NEW MODES
# ============================================================
def run_test_4(gamma=0.05, J=1.0, kappa=0.05):
    log("=" * 70)
    log("TEST 4: NEW MODES CREATED BY THE BRIDGE")
    log("=" * 70)
    log()

    model = 'zz_collective'

    # Isolated
    L_iso, Sg_iso = build_system('isolated', gamma, J)
    evals_iso, R_iso = np.linalg.eig(L_iso)

    # Bridge
    L_br, Sg_br = build_system('bridge', gamma, J, kappa=kappa,
                                bridge_model=model)
    evals_br, R_br = np.linalg.eig(L_br)

    # Match modes: for each bridge eigenvalue, find nearest isolated eigenvalue
    matched = np.zeros(d2, dtype=bool)
    shifts = np.zeros(d2)
    for k in range(d2):
        dists = np.abs(evals_br[k] - evals_iso)
        j = np.argmin(dists)
        shifts[k] = dists[j]

    threshold = 0.001  # modes shifted by more than this are "new"
    new_modes = np.where(shifts > threshold)[0]
    old_modes = np.where(shifts <= threshold)[0]

    log(f"Bridge model: {model}, kappa = {kappa}")
    log(f"Modes shifted by > {threshold}: {len(new_modes)}/{d2}")
    log(f"Modes essentially unchanged: {len(old_modes)}/{d2}")
    log()

    if len(new_modes) > 0:
        log("--- New (shifted) modes ---")
        log(f"{'Mode':>5}  {'Rate':>10}  {'|Freq|':>8}  {'Shift':>10}")
        log("-" * 40)

        order = np.argsort(-shifts)
        for rank in range(min(20, len(new_modes))):
            k = order[rank]
            if shifts[k] <= threshold: break
            log(f"{k:>5}  {-evals_br[k].real:>10.6f}  "
                f"{abs(evals_br[k].imag):>8.4f}  {shifts[k]:>10.6f}")

        log()

        # Spatial structure: which qubits do new modes involve?
        log("--- Spatial structure of new modes ---")
        log("Weight on each qubit (from eigenvector structure):")
        log()

        log(f"{'Mode':>5}  {'Rate':>10}  "
            f"{'Q0':>6}  {'Q1':>6}  {'Q2':>6}  {'Q3':>6}  {'Cross':>6}")
        log("-" * 50)

        for rank in range(min(10, len(new_modes))):
            k = order[rank]
            if shifts[k] <= threshold: break
            vec = R_br[:, k]
            rho_mat = vec.reshape(d, d)

            # Site participation: off-diagonal elements per qubit
            site_wt = np.zeros(N)
            for bit_r in range(d):
                for bit_c in range(d):
                    val = abs(rho_mat[bit_r, bit_c]) ** 2
                    for s in range(N):
                        if ((bit_r >> (N-1-s)) & 1) != ((bit_c >> (N-1-s)) & 1):
                            site_wt[s] += val

            total = np.sum(site_wt)
            if total > 1e-30: site_wt /= total

            # Cross-boundary weight: elements where A and B qubits differ
            cross_wt = 0
            for bit_r in range(d):
                for bit_c in range(d):
                    a_diff = any(((bit_r >> (N-1-s)) & 1) != ((bit_c >> (N-1-s)) & 1)
                                 for s in [0, 1])
                    b_diff = any(((bit_r >> (N-1-s)) & 1) != ((bit_c >> (N-1-s)) & 1)
                                 for s in [2, 3])
                    if a_diff and b_diff:
                        cross_wt += abs(rho_mat[bit_r, bit_c]) ** 2
            norm = np.sum(np.abs(rho_mat) ** 2)
            cross_frac = cross_wt / norm if norm > 1e-30 else 0

            log(f"{k:>5}  {-evals_br[k].real:>10.6f}  "
                f"{site_wt[0]:>6.3f}  {site_wt[1]:>6.3f}  "
                f"{site_wt[2]:>6.3f}  {site_wt[3]:>6.3f}  "
                f"{cross_frac:>6.3f}")

    log()


# ============================================================
# TEST 5: DOES THE BRIDGE HAVE A PI?
# ============================================================
def run_test_5(gamma=0.05, J=1.0, kappa=0.05):
    log("=" * 70)
    log("TEST 5: DOES THE BRIDGE HAVE A Pi?")
    log("=" * 70)
    log()

    model = 'zz_collective'
    L_br, Sg_base = build_system('bridge', gamma, J, kappa=kappa,
                                  bridge_model=model)

    # Build Pi operators in Pauli basis (superoperator on vectorized d^2)
    all_idx = list(iprod(range(4), repeat=N))
    num = 4 ** N  # 256

    # Build Pauli basis matrices
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, N):
            m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m)

    # Pi_chain: per-site I<->X, Y->iZ, Z->iY
    PI_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
    PI_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

    def build_Pi_pauli(all_idx):
        n = len(all_idx)
        Pi = np.zeros((n, n), dtype=complex)
        idx_map = {idx: i for i, idx in enumerate(all_idx)}
        for b, idx_b in enumerate(all_idx):
            mapped = tuple(PI_PERM[i] for i in idx_b)
            sign = 1
            for i in idx_b:
                sign *= PI_SIGN[i]
            a = idx_map[mapped]
            Pi[a, b] = sign
        return Pi

    Pi_pauli = build_Pi_pauli(all_idx)

    # Convert Pi from Pauli basis to computational vectorized basis
    # V: change-of-basis matrix from Pauli to computational
    V = np.zeros((d2, num), dtype=complex)
    for a in range(num):
        V[:, a] = pmats[a].flatten()
    V_inv = np.linalg.inv(V)

    Pi_comp = V @ Pi_pauli @ V_inv

    # Check: Pi L_bridge Pi^-1 + L + 2*Sg*I = ?
    Pi_inv = np.linalg.inv(Pi_comp)

    # Try different Sg values
    log("--- Check Pi_chain conjugation ---")
    log()

    best_Sg = Sg_base; best_err = 1e30
    for Sg_try in np.linspace(0, Sg_base * 5, 200):
        E = Pi_comp @ L_br @ Pi_inv + L_br + 2 * Sg_try * np.eye(d2)
        err = np.max(np.abs(E)) / max(np.max(np.abs(L_br)), 1)
        if err < best_err:
            best_err = err; best_Sg = Sg_try

    log(f"Pi_chain conjugation error (best Sg = {best_Sg:.4f}): {best_err:.2e}")

    if best_err < 1e-8:
        log("Pi_chain WORKS for the bridge system!")
    else:
        log("Pi_chain does NOT work for the bridge system.")
    log()

    # Check Pi^2 = X^N parity
    Pi2_comp = Pi_comp @ Pi_comp
    # X^N parity in computational vectorized basis
    XN = sx.copy()
    for k in range(1, N):
        XN = np.kron(XN, sx)

    # X^N conjugation as superoperator: rho -> X^N rho X^N
    XN_super = np.kron(XN, XN.conj())
    comm_Pi2_XN = np.max(np.abs(Pi2_comp - XN_super))
    log(f"||Pi^2 - X^N conjugation|| = {comm_Pi2_XN:.2e}")

    # Check [X^N, L_bridge] = 0
    comm_XN_L = np.max(np.abs(XN_super @ L_br - L_br @ XN_super))
    log(f"[X^N, L_bridge] = {comm_XN_L:.2e}")

    if comm_XN_L < 1e-10:
        log("X^N parity is CONSERVED by the bridge.")
    else:
        log("X^N parity is BROKEN by the bridge.")
    log()

    # Check Pi_A tensor Pi_B
    # Pi_A acts on qubits 0,1; Pi_B on qubits 2,3
    all_idx_2 = list(iprod(range(4), repeat=2))
    Pi_2 = np.zeros((16, 16), dtype=complex)
    idx_map_2 = {idx: i for i, idx in enumerate(all_idx_2)}
    for b, idx_b in enumerate(all_idx_2):
        mapped = tuple(PI_PERM[i] for i in idx_b)
        sign = 1
        for i in idx_b:
            sign *= PI_SIGN[i]
        a = idx_map_2[mapped]
        Pi_2[a, b] = sign

    # Pi_A x Pi_B in the 4-qubit Pauli basis
    Pi_AB = np.kron(Pi_2, Pi_2)  # 256 x 256 in Pauli basis
    Pi_AB_comp = V @ Pi_AB @ V_inv
    Pi_AB_inv = np.linalg.inv(Pi_AB_comp)

    best_Sg_AB = Sg_base; best_err_AB = 1e30
    for Sg_try in np.linspace(0, Sg_base * 5, 200):
        E = Pi_AB_comp @ L_br @ Pi_AB_inv + L_br + 2 * Sg_try * np.eye(d2)
        err = np.max(np.abs(E)) / max(np.max(np.abs(L_br)), 1)
        if err < best_err_AB:
            best_err_AB = err; best_Sg_AB = Sg_try

    log(f"Pi_A x Pi_B conjugation error (best Sg = {best_Sg_AB:.4f}): "
        f"{best_err_AB:.2e}")

    if best_err_AB < 1e-8:
        log("Pi_A x Pi_B WORKS! The palindromic operator is the tensor product.")
    else:
        log("Pi_A x Pi_B does NOT work.")
    log()

    # Eigenvalue palindrome as backup check
    evals_br = np.linalg.eigvals(L_br)
    for Sg_try in [Sg_base, best_Sg, best_Sg_AB]:
        err_ev, np_ev, _ = palindrome_err(evals_br, Sg_try)
        log(f"Eigenvalue pairing at Sg={Sg_try:.4f}: "
            f"{np_ev}/{d2}, err = {err_ev:.2e}")

    log()


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    t0 = time.time()
    log("N=4 Across the Boundary: Two N=2 Systems Meeting")
    log(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"N={N}, d={d}, d^2={d2}")
    log()

    t1 = time.time()
    results_1 = run_test_1()
    log(f"[Test 1 completed in {time.time() - t1:.1f}s]")
    log()

    t3 = time.time()
    run_test_3()
    log(f"[Test 3 completed in {time.time() - t3:.1f}s]")
    log()

    t2 = time.time()
    run_test_2()
    log(f"[Test 2 completed in {time.time() - t2:.1f}s]")
    log()

    t4 = time.time()
    run_test_4()
    log(f"[Test 4 completed in {time.time() - t4:.1f}s]")
    log()

    t5 = time.time()
    run_test_5()
    log(f"[Test 5 completed in {time.time() - t5:.1f}s]")
    log()

    # ============================================================
    # OVERALL SUMMARY
    # ============================================================
    log("=" * 70)
    log("OVERALL SUMMARY")
    log("=" * 70)
    log()
    log("Test 1: Three spectral architectures compared")
    log("  Chain, isolated pairs, and bridge with multiple coupling models.")
    log()
    log("Test 3: Bidirectionality requirement")
    log("  Whether palindromic structure requires symmetric coupling.")
    log()
    log("Test 2: Kappa sweep")
    log("  How eigenvalue shifts track the opening of the bridge.")
    log()
    log("Test 4: New modes")
    log("  Modes created by the bridge that did not exist in isolation.")
    log()
    log("Test 5: Pi operator")
    log("  Whether Pi_chain or Pi_A x Pi_B conjugates the bridge Liouvillian.")
    log()

    log(f"Total runtime: {time.time() - t0:.1f}s")
    log(f"Results: {OUT_PATH}")
    _outf.close()
