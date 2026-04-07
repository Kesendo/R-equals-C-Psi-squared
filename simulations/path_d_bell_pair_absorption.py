#!/usr/bin/env python3
"""
Path D - Bell pair ESD time as a direct application of the Absorption Theorem.

Question: Does Re(lambda) = -2 gamma <n_XY> (the Absorption Theorem) actually
predict the ESD time of a Bell pair on the central qubits of an N=5 Heisenberg
chain under uniform Z-dephasing?

Setup mirrors cockpit_scaling N=5 chain center_bell run:
  - N = 5 qubits, Heisenberg XXX chain (J = 1)
  - Uniform Z-dephasing on all qubits, gamma = 0.05
  - Initial state: |Phi+><Phi+|_{2,3} (x) |+><+|^{otimes 3}
    (The spectator qubits are in the pure |+> state, NOT maximally mixed.
    This is a non-obvious detail of the cockpit_scaling setup that is not
    documented elsewhere; see the build_initial_state docstring below.)
  - Compare predicted concurrence trajectory against
    cockpit_scaling_N5_chain.csv (center_bell pair)

Two predictions are made:
  (A) EXACT: full spectral evolution rho(t) = sum_k c_k R_k exp(lambda_k t)
      then partial trace -> reduced 2-qubit DM -> Wootters concurrence.
      RESULT: reproduces the empirical CSV to ~1e-3 (CSV write precision).
      Confirms the spectral pipeline matches the C# cockpit engine exactly.

  (B) ANALYTICAL: identify the slowest Bell-pair-bearing eigenmodes with
      <n_XY> = 2 from the spectrum alone, predict t_ESD as -log(threshold)
      divided by the slowest dominant absorption rate.
      RESULT: this naive single-mode picture is WRONG by a factor ~100.
      ESD is set by the destructive interference of ~50 degenerate n_XY=2
      modes through the Wootters concurrence (a nonlinear functional of rho),
      not by single-rate exponential decay. Predicting t_ESD from the spectrum
      alone remains an open analytical question.

Created: April 7, 2026 (path_d session)
"""

import numpy as np
from itertools import product
from pathlib import Path

# ----- Paths (relative to this script's location in simulations/) -----
SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results"
EMP_CSV = RESULTS_DIR / "cockpit_scaling" / "cockpit_scaling_N5_chain.csv"
PLOT_OUT = RESULTS_DIR / "cockpit_scaling" / "path_d_comparison.png"

# ----- Parameters -----
N = 5
J = 1.0
gamma = 0.05
BELL_QUBITS = (2, 3)        # central pair in N=5 chain
THRESHOLD = 1e-6            # concurrence threshold for ESD
T_GRID = np.linspace(0.0, 3.0, 301)

# ----- Pauli matrices and labels -----
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': Xm, 'Y': Ym, 'Z': Zm}
LABELS = ['I', 'X', 'Y', 'Z']


def kron_chain(ops):
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def pauli_string_matrix(s):
    """Build the 2^N x 2^N matrix for a Pauli string like 'IIXIZ'."""
    return kron_chain([PAULIS[c] for c in s])


def n_xy(s):
    """Count X or Y characters in a Pauli string."""
    return sum(1 for c in s if c in ('X', 'Y'))


# ----- Build Hamiltonian (Heisenberg XXX chain) -----
def build_H(N, J):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in (Xm, Ym, Zm):
            ops = [I2] * N
            ops[i] = P
            ops[i + 1] = P
            H += J * kron_chain(ops)
    return H


# ----- Build Liouvillian L (column-stacking vec convention) -----
# vec(A B C) = (C^T (x) A) vec(B)
# So d/dt vec(rho) = L vec(rho), with
#   L = -i (H (x) I - I (x) H^T)             [Hamiltonian part]
#     + sum_k gamma_k [ Z_k (x) Z_k^T - I ]  [pure dephasing on each qubit]
def build_L(N, H, gamma):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    L_H = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    L_D = np.zeros((d * d, d * d), dtype=complex)
    for k in range(N):
        ops = [I2] * N
        ops[k] = Zm
        Zk = kron_chain(ops)
        # Pure dephasing Lindblad with single Hermitian L = sqrt(gamma) Z_k:
        # D[rho] = gamma (Z rho Z - rho)
        L_D += gamma * (np.kron(Zk, Zk.T) - np.kron(Id, Id))
    return L_H + L_D


# ----- Pauli decomposition of an operator (vectorised) -----
# Build the full Pauli basis once: 4^N strings, each as a vec'd matrix.
def build_pauli_basis(N):
    d = 2**N
    strings = [''.join(p) for p in product(LABELS, repeat=N)]
    basis = np.zeros((d * d, len(strings)), dtype=complex)
    for j, s in enumerate(strings):
        M = pauli_string_matrix(s)
        # Hilbert-Schmidt normalised: <P, P>_HS = 2^N, so divide by sqrt(2^N)
        basis[:, j] = M.flatten(order='F') / np.sqrt(d)
    return strings, basis


def pauli_decompose(vec_op, basis):
    """Return complex coefficients of vec_op in the (orthonormal) Pauli basis."""
    return basis.conj().T @ vec_op


# ----- Mode <n_XY> from Pauli decomposition -----
def mode_n_xy(coeffs, strings):
    """Weighted average <n_XY> for an eigenoperator: sum |c_i|^2 n_xy(s_i)."""
    p = np.abs(coeffs) ** 2
    p_sum = p.sum()
    if p_sum < 1e-15:
        return 0.0
    weights = p / p_sum
    return sum(w * n_xy(s) for w, s in zip(weights, strings))


# =====================================================================
# MAIN
# =====================================================================
print(f"N = {N}, J = {J}, gamma = {gamma}")
print(f"Bell pair on qubits {BELL_QUBITS}")
print()

d = 2**N
H = build_H(N, J)
L = build_L(N, H, gamma)
print(f"Liouvillian shape: {L.shape}")

# Eigendecomposition of L (non-Hermitian: get right eigvecs and left eigvecs)
print("Diagonalising L ...")
eigvals, R = np.linalg.eig(L)            # L = R diag(eigvals) R^{-1}
R_inv = np.linalg.inv(R)                 # rows of R_inv are left eigvectors

# Build Pauli basis once
print("Building Pauli basis ...")
strings, P_basis = build_pauli_basis(N)
print(f"  {len(strings)} Pauli strings ({d*d})")

# For each eigenmode, compute Pauli decomposition and <n_XY>
print("Computing <n_XY> per mode ...")
mode_n = np.zeros(len(eigvals))
for k in range(len(eigvals)):
    coeffs = pauli_decompose(R[:, k], P_basis)
    mode_n[k] = mode_n_xy(coeffs, strings)

# ----- (1) Verify Absorption Theorem on this setup -----
predicted_re = -2 * gamma * mode_n
actual_re = eigvals.real
deviation = np.abs(actual_re - predicted_re)
print()
print("=== Absorption Theorem check ===")
print(f"  Re(lambda) vs -2 gamma <n_XY>")
print(f"  max |deviation| = {deviation.max():.3e}")
print(f"  median |deviation| = {np.median(deviation):.3e}")

# ----- (2) Initial state: |Phi+><Phi+|_{2,3} (x) |+><+|^{otimes 3} -----
def build_initial_state(N, bell_qubits):
    # |Phi+> = (|00> + |11>)/sqrt(2)  on the Bell qubits.
    # Spectator qubits are in |+> = (|0> + |1>)/sqrt(2),  NOT maximally mixed.
    # (Confirmed by reading cockpit_scaling N=5 chain CSV at t=0:
    # the 0_1_far_edge reduced state has phi_plus=0.5, psi_plus=0.5, Pur=1,
    # which is exactly |+>|+> reduced from a fully pure spectator product.)
    Phi_plus = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    rho_bell = np.outer(Phi_plus, Phi_plus.conj())  # 4x4
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    rho_plus = np.outer(plus, plus.conj())          # 2x2
    b0, b1 = bell_qubits
    assert b1 == b0 + 1, "this builder assumes adjacent Bell qubits"
    factors = []
    i = 0
    while i < N:
        if i == b0:
            factors.append(rho_bell)   # 4x4 covers two qubits
            i += 2
        else:
            factors.append(rho_plus)
            i += 1
    return kron_chain(factors)

rho0 = build_initial_state(N, BELL_QUBITS)
assert rho0.shape == (d, d)
print()
print(f"Initial state shape: {rho0.shape}, Tr = {np.trace(rho0).real:.6f}")

vec_rho0 = rho0.flatten(order='F')

# Project initial state onto eigenoperators: c_k = (left_k . vec_rho0)
c = R_inv @ vec_rho0   # length d^2

# ----- (A) EXACT spectral evolution -----
print()
print("=== (A) Exact spectral evolution ===")
print("Computing rho(t), reducing to 2-qubit DM, computing concurrence ...")

# Indices to trace OUT (everything except bell qubits)
trace_out = [q for q in range(N) if q not in BELL_QUBITS]


def partial_trace_iter(rho_full, N, keep):
    """Partial trace of an N-qubit density matrix, keeping only `keep` qubits.
    Iterates trace operations one qubit at a time, recomputing axis indices
    after each contraction (because np.trace shifts axes down)."""
    keep = sorted(keep)
    trace = sorted([q for q in range(N) if q not in keep], reverse=True)
    shape = [2] * (2 * N)
    t = rho_full.reshape(shape)
    cur_N = N
    cur_keep = list(range(N))
    for q in trace:
        idx = cur_keep.index(q)
        t = np.trace(t, axis1=idx, axis2=idx + cur_N)
        cur_keep.pop(idx)
        cur_N -= 1
    dim = 2 ** cur_N
    return t.reshape(dim, dim)


def wootters_concurrence(rho2):
    """Wootters concurrence for a 2-qubit density matrix."""
    YY = np.kron(Ym, Ym)
    rho_tilde = YY @ rho2.conj() @ YY
    M = rho2 @ rho_tilde
    # Eigenvalues of M can have small imaginary parts; take real of sqrt of max(0,re).
    eigs = np.linalg.eigvals(M)
    sqrt_eigs = np.sqrt(np.maximum(eigs.real, 0))
    sqrt_eigs = np.sort(sqrt_eigs)[::-1]
    return max(0.0, sqrt_eigs[0] - sqrt_eigs[1] - sqrt_eigs[2] - sqrt_eigs[3])


# Time-evolve and measure concurrence on the grid
conc_pred = np.zeros_like(T_GRID)
purity_pred = np.zeros_like(T_GRID)
for ti, t in enumerate(T_GRID):
    vec_rho_t = R @ (np.exp(eigvals * t) * c)
    rho_t = vec_rho_t.reshape(d, d, order='F')
    rho2 = partial_trace_iter(rho_t, N, BELL_QUBITS)
    conc_pred[ti] = wootters_concurrence(rho2)
    purity_pred[ti] = np.trace(rho2 @ rho2).real

# Find predicted ESD time
esd_pred = None
for ti, c_t in enumerate(conc_pred):
    if c_t < THRESHOLD:
        esd_pred = T_GRID[ti]
        break

print(f"  Predicted t_ESD (first crossing of {THRESHOLD}) = {esd_pred}")
print()
print("  Predicted concurrence at selected times:")
for tt in [0.0, 0.5, 0.9, 1.0, 1.5, 1.9, 2.0, 2.1, 2.3, 2.4, 3.0]:
    idx = np.argmin(np.abs(T_GRID - tt))
    print(f"    t={T_GRID[idx]:.2f}  C_pred={conc_pred[idx]:.6f}  Pur_pred={purity_pred[idx]:.4f}")

# ----- (B) Analytical prediction from spectrum alone -----
print()
print("=== (B) Analytical single-mode prediction ===")
overlap = np.abs(c) ** 2
xy2_mask = (mode_n > 1.5) & (mode_n < 2.5)
print(f"  Modes with <n_XY> in (1.5, 2.5): {xy2_mask.sum()}")
print(f"  Top 5 n_XY=2 modes by initial-state overlap:")
xy2_indices = np.where(xy2_mask)[0]
xy2_sorted = xy2_indices[np.argsort(-overlap[xy2_indices])]
for k in xy2_sorted[:5]:
    print(f"    Re(lambda)={eigvals[k].real:+.6f}  Im={eigvals[k].imag:+.4f}  "
          f"<n_XY>={mode_n[k]:.4f}  |overlap|^2={overlap[k]:.4e}")

# The "slowest dominant" mode is the one with the smallest |Re(lambda)|
# among modes with non-negligible overlap and n_XY ~ 2.
# Spoiler: this gives a result wildly off from empirical (factor ~100).
# Kept here as a documented negative result, not as a working prediction.
sig_overlap = overlap > 1e-6
candidates = np.where(xy2_mask & sig_overlap)[0]
if len(candidates) > 0:
    slowest = candidates[np.argmin(np.abs(eigvals[candidates].real))]
    rate = -eigvals[slowest].real
    t_esd_analytical = -np.log(THRESHOLD) / rate
    print(f"\n  Slowest dominant Bell-bearing absorption mode:")
    print(f"    Re(lambda) = -{rate:.6f}, <n_XY> = {mode_n[slowest]:.4f}")
    print(f"    NAIVE t_ESD ~ -log({THRESHOLD}) / rate = {t_esd_analytical:.4f}")
    print(f"    (This is wrong by ~100x. ESD comes from interference of many")
    print(f"     degenerate n_XY=2 modes through the Wootters functional.)")
else:
    print("  No candidates found.")

print()
print("=== Empirical (cockpit_scaling CSV) ===")
print("  ESD at t = 0.9, revival peak C ~ 0.158 at t = 2.1, final death t = 2.4")

# ---------------------------------------------------------------------
# Comparison plot vs empirical CSV (only if CSV exists in expected place)
# ---------------------------------------------------------------------
if EMP_CSV.exists():
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    print()
    print(f"=== Comparison vs {EMP_CSV.name} ===")
    df = pd.read_csv(EMP_CSV)
    emp = df[df['pair'] == '2_3_center_bell'].sort_values('t').reset_index(drop=True)
    emp = emp[emp['t'] <= T_GRID[-1]]

    # Residual statistics
    residuals = []
    for _, row in emp.iterrows():
        idx = np.argmin(np.abs(T_GRID - row['t']))
        if abs(T_GRID[idx] - row['t']) < 1e-3:
            residuals.append(row['concurrence'] - conc_pred[idx])
    residuals = np.array(residuals)
    print(f"  N points compared:       {len(residuals)}")
    print(f"  max |C_emp - C_pred|:    {np.max(np.abs(residuals)):.6f}")
    print(f"  RMS:                     {np.sqrt(np.mean(residuals**2)):.6f}")
    print(f"  median |residual|:       {np.median(np.abs(residuals)):.6f}")

    # Plot
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax = axes[0]
    ax.plot(emp['t'], emp['concurrence'], 'o', markersize=8,
            label='C# matrix-free engine (cockpit_scaling)',
            markerfacecolor='none', markeredgewidth=1.8, color='#cc3333')
    ax.plot(T_GRID, conc_pred, '-', linewidth=1.8,
            label='Spectral prediction (Path D)', color='#1f6fb3')
    ax.set_ylabel('Concurrence C(t)', fontsize=12)
    ax.set_title(f'Bell pair (qubits {BELL_QUBITS}) on N={N} Heisenberg chain, '
                 f'gamma={gamma}', fontsize=12)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.annotate('ESD', xy=(0.9, 0.02), xytext=(1.2, 0.15),
                arrowprops=dict(arrowstyle='->', color='gray'),
                color='gray', fontsize=10)
    ax.annotate('Revival', xy=(2.1, 0.158), xytext=(2.45, 0.25),
                arrowprops=dict(arrowstyle='->', color='gray'),
                color='gray', fontsize=10)

    ax = axes[1]
    ax.plot(emp['t'], emp['purity'], 'o', markersize=8,
            label='C# matrix-free engine (cockpit_scaling)',
            markerfacecolor='none', markeredgewidth=1.8, color='#cc3333')
    ax.plot(T_GRID, purity_pred, '-', linewidth=1.8,
            label='Spectral prediction (Path D)', color='#1f6fb3')
    ax.set_ylabel('Purity Tr(rho^2)', fontsize=12)
    ax.set_xlabel('Time t', fontsize=12)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(PLOT_OUT, dpi=130, bbox_inches='tight')
    print(f"  Plot saved to: {PLOT_OUT}")
else:
    print()
    print(f"NOTE: Empirical CSV not found at {EMP_CSV}")
    print("      Skipping comparison plot.")
