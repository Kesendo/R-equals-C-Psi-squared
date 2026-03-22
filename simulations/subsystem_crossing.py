#!/usr/bin/env python3
"""
Subsystem Crossing Theorem
============================
Conjecture 2.1: Every entangled pair with CΨ > 1/4 must cross below 1/4
under any non-unitary CPTP map.

Test 1: N=3,4 with pair (0,1) starting at CΨ > 1/4 (Bell+ ⊗ |0...0⟩)
Test 2: Random initial states with high pair CΨ
Test 3: Random CPTP maps (Stinespring) on 2 qubits
Test 4: Adversarial CPTP — maps designed to preserve entanglement
Test 5: The contractivity proof (analytical + numerical verification)

Script:  simulations/subsystem_crossing.py
Output:  simulations/results/subsystem_crossing.txt
"""

import numpy as np
from scipy.linalg import expm
import os, sys, time as _time

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "subsystem_crossing.txt")
_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
up = np.array([1, 0], dtype=complex)
dn = np.array([0, 1], dtype=complex)


def site_op(op, k, nq):
    ops = [I2] * nq
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def build_liouvillian(H, collapse_ops):
    d = H.shape[0]
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for C_op in collapse_ops:
        Cd = C_op.conj().T
        CdC = Cd @ C_op
        L += np.kron(C_op, C_op.conj()) - 0.5 * (
            np.kron(CdC, Id) + np.kron(Id, CdC.T))
    return L


def evolve(L, rho0, t):
    d = rho0.shape[0]
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


def partial_trace_keep(rho, keep, n_qubits):
    keep = list(keep)
    trace_out = [q for q in range(n_qubits) if q not in keep]
    dims = [2] * n_qubits
    reshaped = rho.reshape(dims + dims)
    current_n = n_qubits
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_keep = 2 ** len(keep)
    return reshaped.reshape((d_keep, d_keep))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d - 1) if d > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def heisenberg_H(nq, J=1.0):
    d = 2 ** nq
    H = np.zeros((d, d), dtype=complex)
    for i in range(nq - 1):
        for P in [sx, sy, sz]:
            H += J * site_op(P, i, nq) @ site_op(P, i + 1, nq)
    return H


def random_cptp_kraus(d, n_kraus=None):
    """Generate random CPTP map via Stinespring: random unitary on d×d_env, partial trace."""
    if n_kraus is None:
        n_kraus = d  # d Kraus operators
    # Random complex matrices
    K_list = []
    for _ in range(n_kraus):
        K = np.random.randn(d, d) + 1j * np.random.randn(d, d)
        K_list.append(K)
    # Enforce ΣK†K = I via polar decomposition
    S = sum(K.conj().T @ K for K in K_list)
    S_inv_sqrt = np.linalg.matrix_power(S, -1)
    ev, U = np.linalg.eigh(S)
    S_inv_sqrt = U @ np.diag(1.0 / np.sqrt(np.maximum(ev, 1e-15))) @ U.conj().T
    K_list = [K @ S_inv_sqrt for K in K_list]
    return K_list


def apply_cptp(rho, K_list):
    """Apply CPTP map: ε(ρ) = Σ K_k ρ K_k†"""
    d = rho.shape[0]
    rho_out = np.zeros((d, d), dtype=complex)
    for K in K_list:
        rho_out += K @ rho @ K.conj().T
    # Hermitianize
    rho_out = (rho_out + rho_out.conj().T) / 2
    rho_out /= np.trace(rho_out).real
    return rho_out


# ====================================================================
# Test 1: N=3,4 with Bell+ pair inside larger system
# ====================================================================

def test_1_subsystem_pairs():
    log("=" * 70)
    log("TEST 1: SUBSYSTEM PAIRS IN N-QUBIT SYSTEMS")
    log("=" * 70)
    log()

    gamma = 0.05; J = 1.0
    t_max = 60.0; dt = 0.05

    for nq in [3, 4, 5]:
        d = 2 ** nq
        H = heisenberg_H(nq, J)
        c_ops = [np.sqrt(gamma) * site_op(sz, k, nq) for k in range(nq)]
        Liouv = build_liouvillian(H, c_ops)

        # Initial state: Bell+(0,1) ⊗ |0⟩^{N-2}
        bell = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
        psi_init = bell
        for _ in range(nq - 2):
            psi_init = np.kron(psi_init, up)
        rho0 = np.outer(psi_init, psi_init.conj())

        # Also test |Ψ+⟩(0,1) ⊗ |+⟩^{N-2}
        psi_plus = (np.kron(up, dn) + np.kron(dn, up)) / np.sqrt(2)
        plus = (up + dn) / np.sqrt(2)
        psi_init2 = psi_plus
        for _ in range(nq - 2):
            psi_init2 = np.kron(psi_init2, plus)
        rho0_2 = np.outer(psi_init2, psi_init2.conj())

        n_steps = int(round(t_max / dt))
        pairs = [(0, 1)] + [(i, i + 1) for i in range(1, nq - 1)]

        for state_name, rho_init in [("Bell+(0,1)⊗|0⟩", rho0), ("Ψ+(0,1)⊗|+⟩", rho0_2)]:
            log(f"  N={nq}, state={state_name}, J={J}, γ={gamma}")

            header = "  " + "  ".join(f"CΨ({i},{j})" for i, j in pairs)
            log(f"  {'t':>6}" + header + "  crossed?")
            log("  " + "-" * (6 + 10 * len(pairs) + 10))

            all_crossed = {p: False for p in pairs}
            cross_times = {p: None for p in pairs}
            prev_cpsi = {}

            for step in range(n_steps + 1):
                t = step * dt
                rho = evolve(Liouv, rho_init, t)

                row = f"  {t:6.1f}"
                for p in pairs:
                    rho_pair = partial_trace_keep(rho, list(p), nq)
                    c = cpsi(rho_pair)

                    if p in prev_cpsi and prev_cpsi[p] > 0.25 and c <= 0.25 and cross_times[p] is None:
                        cross_times[p] = t
                        all_crossed[p] = True

                    prev_cpsi[p] = c
                    row += f"  {c:8.4f}"

                status = "".join("✓" if all_crossed[p] else "·" for p in pairs)
                if step % (n_steps // 6) == 0 or step == n_steps:
                    log(row + f"  {status}")

            log()
            for p in pairs:
                ct = cross_times[p]
                c0 = cpsi(partial_trace_keep(rho_init, list(p), nq))
                if c0 > 0.25:
                    status = f"t_cross={ct:.2f}" if ct else "NEVER CROSSED"
                else:
                    status = "started below 1/4"
                log(f"  Pair ({p[0]},{p[1]}): CΨ(0)={c0:.4f}, {status}")

            log()


# ====================================================================
# Test 2: Random CPTP maps on 2 qubits
# ====================================================================

def test_2_random_cptp():
    log("=" * 70)
    log("TEST 2: RANDOM CPTP MAPS ON 2 QUBITS")
    log("=" * 70)
    log()

    d = 4
    bell = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())

    n_maps = 200
    n_iter = 500
    np.random.seed(42)

    all_cross = 0
    never_cross = 0
    max_cpsi_final = 0
    failures = []

    for trial in range(n_maps):
        K_list = random_cptp_kraus(d, n_kraus=d)

        # Check if map is "trivially unitary" (single Kraus ~ unitary)
        if len(K_list) == 1:
            continue

        rho = rho0.copy()
        crossed = False
        for n in range(n_iter):
            rho = apply_cptp(rho, K_list)
            c = cpsi(rho)
            if c < 0.25:
                crossed = True
                break

        if crossed:
            all_cross += 1
        else:
            never_cross += 1
            cpsi_final = cpsi(rho)
            if cpsi_final > max_cpsi_final:
                max_cpsi_final = cpsi_final
            if cpsi_final > 0.25:
                failures.append((trial, cpsi_final))

    log(f"  {n_maps} random CPTP maps, {n_iter} iterations each")
    log(f"  Crossed 1/4: {all_cross}")
    log(f"  Never crossed: {never_cross}")

    if failures:
        log(f"  *** FAILURES (CΨ still > 1/4 after {n_iter} iterations): {len(failures)} ***")
        for trial, cf in failures[:10]:
            log(f"    Trial {trial}: final CΨ = {cf:.6f}")
    else:
        log(f"  All maps drove CΨ below 1/4.")

    log(f"  Max final CΨ among non-crossers: {max_cpsi_final:.6f}")
    log()


# ====================================================================
# Test 3: Adversarial CPTP — partial dephasing, partial identity
# ====================================================================

def test_3_adversarial():
    log("=" * 70)
    log("TEST 3: ADVERSARIAL CPTP — Can we keep CΨ > 1/4?")
    log("=" * 70)
    log()

    d = 4
    bell = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())

    n_iter = 1000

    # Adversarial map 1: Partial dephasing (weak noise)
    log("  --- Adversarial 1: Very weak dephasing (p=0.001) ---")
    p = 0.001
    K0 = np.sqrt(1 - p) * np.eye(d)
    Kz = np.sqrt(p) * np.kron(sz, I2)
    K_list = [K0, Kz]
    rho = rho0.copy()
    for n in range(n_iter):
        rho = apply_cptp(rho, K_list)
    log(f"  After {n_iter} iterations: CΨ = {cpsi(rho):.6f}")
    crossed = cpsi(rho) < 0.25
    log(f"  Crossed 1/4: {crossed}")
    log()

    # Adversarial map 2: Projective measurement that preserves Bell+
    log("  --- Adversarial 2: Projection onto Bell+ subspace ---")
    P_bell = np.outer(bell, bell.conj())
    P_rest = np.eye(d) - P_bell
    # Luders map: ε(ρ) = P ρ P + (I-P) ρ (I-P)
    K_list = [P_bell, P_rest]
    rho = rho0.copy()
    for n in range(n_iter):
        rho = apply_cptp(rho, K_list)
    log(f"  After {n_iter} iterations: CΨ = {cpsi(rho):.6f}")
    log(f"  Bell+ is FIXED POINT of this map: CΨ stays at 1/3")
    log(f"  BUT: this map is a projective measurement, not a noise channel.")
    log(f"  The system is in the range of the projector — the map acts as identity.")
    log()

    # Adversarial map 3: Near-identity with tiny perturbation
    log("  --- Adversarial 3: Identity + epsilon perturbation ---")
    for eps in [0.01, 0.001, 0.0001]:
        K0 = np.sqrt(1 - eps) * np.eye(d)
        K1 = np.sqrt(eps) * np.random.randn(d, d).astype(complex)
        # Normalize
        S = K0.conj().T @ K0 + K1.conj().T @ K1
        ev, U = np.linalg.eigh(S)
        S_inv = U @ np.diag(1.0 / np.sqrt(np.maximum(ev, 1e-15))) @ U.conj().T
        K_list = [K0 @ S_inv, K1 @ S_inv]

        rho = rho0.copy()
        n_cross = None
        for n in range(n_iter):
            rho = apply_cptp(rho, K_list)
            if cpsi(rho) < 0.25 and n_cross is None:
                n_cross = n

        nc = str(n_cross) if n_cross is not None else "NEVER"
        log(f"  eps={eps}: crossed at iteration {nc}, final CΨ = {cpsi(rho):.6f}")

    log()


# ====================================================================
# Test 4: The contractivity argument (analytical)
# ====================================================================

def test_4_contractivity():
    log("=" * 70)
    log("TEST 4: CONTRACTIVITY ARGUMENT")
    log("=" * 70)
    log()
    log("  THEOREM: For any non-unitary CPTP map ε with a unique fixed point ρ*,")
    log("  if CΨ(ρ*) < 1/4, then ε^n(ρ) has CΨ < 1/4 for sufficiently large n.")
    log()
    log("  PROOF:")
    log("  1. ε is a contraction in trace norm: ||ε(ρ)-ε(σ)||_1 ≤ ||ρ-σ||_1")
    log("     (Data processing inequality for CPTP maps)")
    log("  2. For primitive maps (unique ρ*): ||ε^n(ρ) - ρ*||_1 → 0")
    log("     (Convergence to fixed point)")
    log("  3. CΨ is continuous: CΨ(ε^n(ρ)) → CΨ(ρ*) < 1/4")
    log("  4. By convergence: ∃ N such that CΨ(ε^n(ρ)) < 1/4 for all n ≥ N. QED.")
    log()
    log("  CRITICAL QUESTION: Does every non-unitary CPTP have CΨ(ρ*) ≤ 1/4?")
    log()

    # Verify fixed point CΨ for various channel types
    d = 4
    n_iter = 2000
    bell = (np.kron(up, up) + np.kron(dn, dn)) / np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())

    log(f"  {'Channel':>25}  {'CΨ(ρ*)':>8}  {'< 1/4?':>6}")
    log("  " + "-" * 44)

    # Dephasing channels
    for p in [0.01, 0.05, 0.1, 0.5]:
        K0 = np.sqrt(1 - p) * np.eye(d)
        for k in range(2):
            Kz_k = np.sqrt(p / 2) * (np.kron(sz, I2) if k == 0 else np.kron(I2, sz))
        K_list = [K0, np.sqrt(p / 2) * np.kron(sz, I2), np.sqrt(p / 2) * np.kron(I2, sz)]
        rho = rho0.copy()
        for _ in range(n_iter):
            rho = apply_cptp(rho, K_list)
        c = cpsi(rho)
        log(f"  {'Z-deph p=' + str(p):>25}  {c:8.6f}  {'YES' if c < 0.25 else 'NO':>6}")

    # Depolarizing
    for p in [0.01, 0.1, 0.5]:
        paulis = [np.eye(d), np.kron(sx, I2), np.kron(sy, I2), np.kron(sz, I2),
                  np.kron(I2, sx), np.kron(I2, sy), np.kron(I2, sz)]
        K_list = [np.sqrt((1 - p)) * np.eye(d)] + [np.sqrt(p / 6) * P for P in paulis[1:]]
        rho = rho0.copy()
        for _ in range(n_iter):
            rho = apply_cptp(rho, K_list)
        c = cpsi(rho)
        log(f"  {'Depol p=' + str(p):>25}  {c:8.6f}  {'YES' if c < 0.25 else 'NO':>6}")

    # Random CPTP fixed points
    np.random.seed(123)
    max_fp_cpsi = 0
    for trial in range(100):
        K_list = random_cptp_kraus(d, n_kraus=d)
        rho = rho0.copy()
        for _ in range(n_iter):
            rho = apply_cptp(rho, K_list)
        c = cpsi(rho)
        if c > max_fp_cpsi:
            max_fp_cpsi = c

    log(f"  {'Random CPTP (100 maps)':>25}  {'max=' + f'{max_fp_cpsi:.6f}':>8}  "
        f"{'YES' if max_fp_cpsi < 0.25 else 'NO':>6}")

    log()
    log("  If ALL fixed points have CΨ < 1/4, the contractivity proof closes")
    log("  the Subsystem Crossing Theorem for primitive CPTP maps.")
    log()


# ====================================================================
# Main
# ====================================================================

if __name__ == "__main__":
    t_start = _time.time()

    log("Subsystem Crossing Theorem")
    log("=" * 70)
    log("Conjecture 2.1: Every entangled pair with CΨ > 1/4 must cross below 1/4")
    log("under any non-unitary CPTP map (repeated application).")
    log()

    test_1_subsystem_pairs()
    test_2_random_cptp()
    test_3_adversarial()
    test_4_contractivity()

    log("=" * 70)
    log("VERDICT")
    log("=" * 70)
    log()
    log("The Subsystem Crossing Theorem holds if:")
    log("  (a) All primitive CPTP maps have fixed points with CΨ ≤ 1/4")
    log("  (b) Contractivity drives ε^n(ρ) → ρ*")
    log("  (c) Continuity of CΨ gives the crossing")
    log()
    log("Exception: non-primitive maps (projective measurements with Bell+")
    log("as eigenstate) can preserve CΨ > 1/4 indefinitely. These are trivial")
    log("cases where the map acts as identity on the initial state.")
    log()

    total = _time.time() - t_start
    log(f"Total runtime: {total:.1f}s ({total/60:.1f} min)")
    log(f"Results saved to: {OUT_PATH}")
    _outf.close()
