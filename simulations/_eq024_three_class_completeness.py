#!/usr/bin/env python3
"""EQ-024 sub-question 1: Three-class completeness via random F71-symmetric sampling.

J_BLIND_RECEIVER_CLASSES.md identifies three classes of states whose J-Jacobian
vanishes (J-blind):
  Class 1: DFS of L_D + per-bond eigenstate          (|0⟩⁵, |1⟩⁵ at N=5)
  Class 2: H-degenerate L_D-closed block             (GHZ via {|0⟩⁵, |1⟩⁵})
  Class 3: M_α-polynomial subspace under Heisenberg  (|+⟩⁵, Dicke |S_k⟩)

The doc flags as open: are these three classes exhaustive over the J-blind set,
or does another mechanism exist? Direction 3 swept 39+4 F71-symmetric receivers
optimized for MAX capacity — all came in at 11.56-12.07 bits, none blind. But
that was a max-search, not a random sample.

This script: random-sample F71-symmetric product states, classify each as in/out
of Classes 1-3, and compute J-channel capacity for the OUT-of-class ones. Any
out-of-class state with C ≈ 0 would identify a fourth blind mechanism. Strong
evidence for completeness if all out-of-class samples are J-sensitive.

Setup matches J_BLIND_RECEIVER_CLASSES.md anchor:
  N=5, Heisenberg H = Σ_b (XX+YY+ZZ)_b, J=1.0 uniform, γ₀=0.05 per site,
  t-grid [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], 25 features × 6 times.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))


# Paulis
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

# Setup (matches J_BLIND_RECEIVER_CLASSES.md)
N = 5
d = 2 ** N
d2 = d * d
N_BONDS = N - 1
GAMMA_0 = 0.05
J_REF = 1.0
T_POINTS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
DJ = 1e-4
SPREAD = 0.02
SIGMA_NOISE = 0.01
P_TOTAL = N_BONDS * SPREAD ** 2

CLASS_TOL = 1e-6   # tolerance for class-membership checks
CAPACITY_BLIND = 0.05  # bits — below this we call "J-blind"


def site_op(op, k):
    ops = [I2] * N
    ops[k] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


# Heisenberg per-bond operators (XX+YY+ZZ)
_BOND_BLOCKS = []
for b in range(N_BONDS):
    block = (site_op(sx, b) @ site_op(sx, b + 1)
             + site_op(sy, b) @ site_op(sy, b + 1)
             + site_op(sz, b) @ site_op(sz, b + 1))
    _BOND_BLOCKS.append(block)


def build_H_perbond(J_vec):
    H = np.zeros((d, d), dtype=complex)
    for b in range(N_BONDS):
        H += J_vec[b] * _BOND_BLOCKS[b]
    return H


def build_L(H):
    Id = np.eye(d)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        Zk = site_op(sz, k)
        L += GAMMA_0 * (np.kron(Zk, Zk.conj()) - np.eye(d2))
    return L


def evolve(L, rho0, t):
    v = expm(L * t) @ rho0.flatten()
    rho = v.reshape(d, d)
    return (rho + rho.conj().T) / 2


_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def ptrace_keep(rho, keep):
    keep = list(keep)
    trace_out = [q for q in range(N) if q not in keep]
    dims = [2] * N
    reshaped = rho.reshape(dims + dims)
    current_n = N
    for q in sorted(trace_out, reverse=True):
        reshaped = np.trace(reshaped, axis1=q, axis2=q + current_n)
        current_n -= 1
    d_k = 2 ** len(keep)
    return reshaped.reshape((d_k, d_k))


def purity(rho):
    return float(np.trace(rho @ rho).real)


def psi_norm(rho):
    d_r = rho.shape[0]
    l1 = float(np.sum(np.abs(rho)) - np.sum(np.abs(np.diag(rho))))
    return l1 / (d_r - 1) if d_r > 1 else 0.0


def cpsi(rho):
    return purity(rho) * psi_norm(rho)


def extract_features(rho):
    f = []
    for q in range(N):
        f.append(purity(ptrace_keep(rho, [q])))
    for i in range(N):
        for j in range(i + 1, N):
            f.append(cpsi(ptrace_keep(rho, [i, j])))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sz @ rho_q))))
    for q in range(N):
        rho_q = ptrace_keep(rho, [q])
        f.append(float(np.real(np.trace(sx @ rho_q))))
    return np.array(f)


def extract_timeseries(rho0, J_vec):
    H = build_H_perbond(J_vec)
    L = build_L(H)
    f = []
    for t in T_POINTS:
        rho = evolve(L, rho0, t)
        f.extend(extract_features(rho).tolist())
    return np.array(f)


def compute_jacobian_J(rho0):
    J_ref = [J_REF] * N_BONDS
    y_ref = extract_timeseries(rho0, J_ref)
    A = np.zeros((len(y_ref), N_BONDS))
    for b in range(N_BONDS):
        J_plus = list(J_ref); J_plus[b] += DJ
        y_plus = extract_timeseries(rho0, J_plus)
        J_minus = list(J_ref); J_minus[b] -= DJ
        y_minus = extract_timeseries(rho0, J_minus)
        A[:, b] = (y_plus - y_minus) / (2.0 * DJ)
    return A


def waterfilling_capacity(sv):
    n = len(sv)
    if np.max(sv) < 1e-8:
        return 0.0
    gains = sv ** 2 / SIGMA_NOISE ** 2
    idx = np.argsort(gains)[::-1]
    gains_sorted = gains[idx]
    inv_gains = 1.0 / gains_sorted
    k_used = 1
    powers = None
    for k in range(n, 0, -1):
        mu = (P_TOTAL + np.sum(inv_gains[:k])) / k
        powers_try = mu - inv_gains[:k]
        if np.all(powers_try >= 0):
            k_used = k
            powers = powers_try
            break
    C = 0.0
    for i in range(k_used):
        snr_i = powers[i] * gains_sorted[i]
        C += 0.5 * np.log2(1 + snr_i)
    return float(C)


# =============================================================================
# F71-symmetric product state parametrization
# =============================================================================
def bloch_state(theta, phi):
    return np.array([np.cos(theta / 2),
                     np.exp(1j * phi) * np.sin(theta / 2)], dtype=complex)


def f71_product(angles):
    """F71-symmetric product state at N=5 from 3 (θ, φ) pairs.

    State = |a⟩|b⟩|c⟩|b⟩|a⟩, mirror-symmetric across center.
    angles = (θ_a, φ_a, θ_b, φ_b, θ_c, φ_c).
    """
    th_a, ph_a, th_b, ph_b, th_c, ph_c = angles
    a = bloch_state(th_a, ph_a)
    b = bloch_state(th_b, ph_b)
    c = bloch_state(th_c, ph_c)
    psi = a
    for v in [b, c, b, a]:
        psi = np.kron(psi, v)
    return psi


# =============================================================================
# Class membership tests
# =============================================================================
def is_class_1_or_2(psi):
    """Class 1 or 2 (within product F71-symmetric): ψ in span{|0⟩⁵, |1⟩⁵}.

    For pure product F71-symmetric states, this is equivalent to ψ being |0⟩⁵
    or |1⟩⁵ exactly (since a product can only span a 1D subspace, and
    {|0⟩⁵, |1⟩⁵} are the only product states inside the 2-dim Class-2 block).
    """
    p = np.abs(psi) ** 2
    weight_on_block = p[0] + p[d - 1]
    return weight_on_block > 1.0 - CLASS_TOL


# Hadamard for x-basis rotation: H X H = Z, so U_x = H^⊗N takes x-basis to z-basis
H1 = (1.0 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
# S H for y-basis: actual rotation Y → Z is U_y = (H @ S†)^⊗N where S = diag(1, i)
# We want a unitary U_y such that U_y σ_y U_y† = σ_z. Choose U_y = (1/√2) [[1, -i],[1, i]]
U_y_single = (1.0 / np.sqrt(2)) * np.array([[1, -1j], [1, 1j]], dtype=complex)


def kron_n(op):
    out = op
    for _ in range(N - 1):
        out = np.kron(out, op)
    return out


U_X = kron_n(H1)
U_Y = kron_n(U_y_single)
U_Z = np.eye(d, dtype=complex)


def project_on_dicke_in_z_basis(psi):
    """Compute ‖projection of ψ onto symmetric subspace (Dicke states)‖²."""
    # Symmetric subspace at N=5 is span of Dicke |S_k⟩, k=0..N (dim N+1=6).
    # Equivalently: ψ symmetric under all site permutations.
    # Test: project onto each Dicke state |S_k⟩ = (1/sqrt(C(N,k))) Σ_{popcount=k} |s⟩
    weight = 0.0
    for k in range(N + 1):
        # |S_k⟩ amplitude on basis state i with popcount(i)=k is 1/sqrt(C(N,k))
        # ⟨S_k|ψ⟩ = (1/sqrt(C(N,k))) Σ_{popcount(i)=k} ψ[i]
        amp_sum = 0.0 + 0.0j
        for i in range(d):
            if bin(i).count("1") == k:
                amp_sum += psi[i]
        n_pop = sum(1 for i in range(d) if bin(i).count("1") == k)
        norm = 1.0 / np.sqrt(n_pop)
        weight += np.abs(norm * amp_sum) ** 2
    return float(weight)


def is_class_3(psi):
    """Class 3 (Heisenberg): ψ in M_α-polynomial subspace for some α ∈ {x, y, z}.

    Equivalently: ψ becomes permutation-symmetric (Dicke) after rotating from
    α-basis to z-basis via U_α^† (the single-qubit rotation that takes α-basis
    eigenstates to z-basis eigenstates).
    """
    for U in [U_X, U_Y, U_Z]:
        psi_rot = U.conj().T @ psi
        weight = project_on_dicke_in_z_basis(psi_rot)
        if weight > 1.0 - CLASS_TOL:
            return True
    return False


def classify(psi):
    """Return ('in-class', class_label) or ('outside', None)."""
    if is_class_1_or_2(psi):
        return 'in-class', '1or2'
    if is_class_3(psi):
        return 'in-class', '3'
    return 'outside', None


# =============================================================================
# Main: random sample, classify, compute capacity for outside-class states
# =============================================================================
def main():
    n_samples = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    rng = np.random.default_rng(seed)

    print(f"EQ-024 Three-Class Completeness test (random F71-symmetric product states)")
    print(f"  N={N}, Heisenberg, γ₀={GAMMA_0}, J={J_REF}, dJ={DJ}")
    print(f"  Random samples: {n_samples}, seed={seed}")
    print(f"  Tolerance: class membership {CLASS_TOL:.0e}, "
          f"J-blind capacity threshold {CAPACITY_BLIND} bits")
    print()

    in_class = {'1or2': 0, '3': 0}
    outside_capacities = []
    blind_outside = []  # candidate fourth-class states

    t0 = time.time()
    for i in range(n_samples):
        # 6 random angles: θ in [0, π], φ in [0, 2π)
        thetas = rng.uniform(0, np.pi, 3)
        phis = rng.uniform(0, 2 * np.pi, 3)
        angles = (thetas[0], phis[0], thetas[1], phis[1], thetas[2], phis[2])
        psi = f71_product(angles)

        verdict, label = classify(psi)
        if verdict == 'in-class':
            in_class[label] += 1
            continue

        # Outside all classes — compute capacity
        rho0 = np.outer(psi, psi.conj())
        rho0 = (rho0 + rho0.conj().T) / 2.0
        rho0 /= np.trace(rho0).real
        A = compute_jacobian_J(rho0)
        sv = np.linalg.svd(A, compute_uv=False)
        C = waterfilling_capacity(sv)
        outside_capacities.append(C)

        if C < CAPACITY_BLIND:
            blind_outside.append({
                'angles_deg': [float(np.degrees(a)) for a in angles],
                'capacity_bits': C,
                'sv_max': float(sv[0]),
                'sv': [float(s) for s in sv],
            })

        if (i + 1) % 10 == 0:
            elapsed = time.time() - t0
            print(f"  [{i + 1}/{n_samples}] elapsed {elapsed:.1f}s, "
                  f"in-class {sum(in_class.values())}, "
                  f"outside {len(outside_capacities)}, "
                  f"blind-outside {len(blind_outside)}")

    print()
    print("=== Summary ===")
    print(f"Total samples: {n_samples}")
    print(f"In-class breakdown: Class 1or2: {in_class['1or2']}, Class 3: {in_class['3']}")
    print(f"Outside-class samples: {len(outside_capacities)}")
    if outside_capacities:
        arr = np.array(outside_capacities)
        print(f"  capacity range: {arr.min():.4f} – {arr.max():.4f} bits")
        print(f"  capacity mean: {arr.mean():.4f}, std: {arr.std():.4f}")
        print(f"  capacity median: {np.median(arr):.4f}")
    print()
    print(f"Blind outside (C < {CAPACITY_BLIND} bits, fourth-class candidates): "
          f"{len(blind_outside)}")
    if blind_outside:
        print("  *** WARNING: fourth-class candidate(s) found ***")
        for entry in blind_outside[:5]:
            print(f"    angles_deg={entry['angles_deg']}  C={entry['capacity_bits']:.6f}  "
                  f"sv_max={entry['sv_max']:.3e}")
    else:
        print("  None found — strong empirical support for three-class completeness")
        print("  within F71-symmetric product states.")

    print()
    print(f"Total runtime: {time.time() - t0:.1f}s")

    # Save raw data
    out = SCRIPT_DIR / "results" / "eq024_three_class_completeness.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', encoding='utf-8') as f:
        json.dump({
            'config': {
                'N': N, 'gamma_0': GAMMA_0, 'J_ref': J_REF, 'dJ': DJ,
                't_points': T_POINTS, 'spread': SPREAD, 'sigma_noise': SIGMA_NOISE,
                'n_samples': n_samples, 'seed': seed,
                'class_tol': CLASS_TOL, 'capacity_blind_threshold': CAPACITY_BLIND,
            },
            'in_class_counts': in_class,
            'outside_capacities': outside_capacities,
            'blind_outside_candidates': blind_outside,
        }, f, indent=1)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
