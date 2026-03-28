#!/usr/bin/env python3
"""
Quantitative test: does the sacrifice-zone formula predict the
IBM hardware measurement?

IBM data: 5-qubit chain [Q85, Q86, Q87, Q88, Q94] on ibm_torino.
Q85 has T2* = 3.7 us (natural sacrifice zone, 26x more dephasing).

Three configurations measured:
  - Selective DD (DD on all except sacrifice Q85)
  - Uniform DD (DD on all)
  - No DD (free evolution)

We run the Lindblad simulation with IBM's ACTUAL dephasing rates
and compare SumMI.
"""
import numpy as np
from scipy.linalg import expm
import json
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_list(mats):
    result = mats[0]
    for m in mats[1:]:
        result = np.kron(result, m)
    return result


def build_heisenberg_chain(N, J=1.0):
    """Heisenberg coupling XX+YY+ZZ on a chain."""
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for pauli in [sx, sy, sz]:
            ops = [I2] * N
            ops[i] = pauli
            ops[i + 1] = pauli
            H += J * kron_list(ops)
    return H


def build_liouvillian(H, gamma_list):
    """Build Liouvillian with per-qubit Z-dephasing."""
    d = H.shape[0]
    N = int(np.log2(d))
    eye = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, eye) - np.kron(eye, H.T))
    for i, gamma in enumerate(gamma_list):
        if gamma > 0:
            ops = [I2] * N
            ops[i] = sz
            c = np.sqrt(gamma) * kron_list(ops)
            cd = c.conj().T
            cdc = cd @ c
            L += (np.kron(c, c.conj())
                  - 0.5 * np.kron(cdc, eye)
                  - 0.5 * np.kron(eye, cdc.T))
    return L


def compute_mi(rho, i, j, N):
    """Mutual information between qubits i and j."""
    d = 2**N
    # Partial traces
    rho_i = partial_trace_to_qubit(rho, i, N)
    rho_j = partial_trace_to_qubit(rho, j, N)
    rho_ij = partial_trace_to_pair(rho, i, j, N)

    S_i = von_neumann_entropy(rho_i)
    S_j = von_neumann_entropy(rho_j)
    S_ij = von_neumann_entropy(rho_ij)

    return max(0, S_i + S_j - S_ij)


def von_neumann_entropy(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-15]
    return -np.sum(ev * np.log2(ev))


def partial_trace_to_qubit(rho, qubit, N):
    """Trace out everything except one qubit."""
    d = 2**N
    rho_reshaped = rho.reshape([2]*2*N)
    # Sum over all indices except qubit
    axes_to_trace = []
    for k in range(N):
        if k != qubit:
            axes_to_trace.append(k)
    # Need to trace pairs (k, k+N)
    result = rho.copy()
    d_curr = d
    traced = 0
    rho_q = np.zeros((2, 2), dtype=complex)
    # Simple approach: build projectors
    for a in range(2):
        for b in range(2):
            # |a><b| on qubit, I on rest
            proj = np.zeros((d, d), dtype=complex)
            for basis in range(d):
                for basis2 in range(d):
                    bits1 = [(basis >> (N-1-k)) & 1 for k in range(N)]
                    bits2 = [(basis2 >> (N-1-k)) & 1 for k in range(N)]
                    if bits1[qubit] == a and bits2[qubit] == b:
                        # Check all other qubits match
                        match = all(bits1[k] == bits2[k] for k in range(N) if k != qubit)
                        if match:
                            proj[basis, basis2] = 1.0
            rho_q[a, b] = np.trace(proj @ rho)
    return rho_q


def partial_trace_to_pair(rho, i, j, N):
    """Trace out everything except qubits i and j."""
    d = 2**N
    rho_ij = np.zeros((4, 4), dtype=complex)
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for dd in range(2):
                    val = 0
                    for basis in range(d):
                        for basis2 in range(d):
                            bits1 = [(basis >> (N-1-k)) & 1 for k in range(N)]
                            bits2 = [(basis2 >> (N-1-k)) & 1 for k in range(N)]
                            if (bits1[i] == a and bits1[j] == b and
                                bits2[i] == c and bits2[j] == dd):
                                match = all(bits1[k] == bits2[k]
                                          for k in range(N) if k != i and k != j)
                                if match:
                                    val += rho[basis, basis2]
                    rho_ij[2*a+b, 2*c+dd] = val
    return rho_ij


def sum_mi(rho, N):
    """Sum of MI over adjacent pairs."""
    total = 0
    for i in range(N - 1):
        total += compute_mi(rho, i, i + 1, N)
    return total


# ================================================================
# IBM Hardware Parameters
# ================================================================
N = 5
chain = [85, 86, 87, 88, 94]

# From IBM calibration data (March 24, 2026)
T2star = [3.73, 61.35, 97.54, 67.99, 95.03]  # microseconds
T1 = [2.84, 211.22, 272.92, 155.53, 295.90]

# Dephasing rates: gamma = 1/(2*T2*) for pure dephasing
# (factor 2 because T2* = 1/(1/(2T1) + gamma_phi), and we want gamma_phi)
# Simplified: use gamma = 1/T2* as effective dephasing rate
gamma_hw = [1.0/t for t in T2star]  # in MHz (1/us)

print("=" * 65)
print("IBM SACRIFICE ZONE: Formula vs Hardware")
print("=" * 65)
print(f"\nChain: {chain}")
print(f"T2*:   {['%.1f' % t for t in T2star]} us")
print(f"gamma: {['%.4f' % g for g in gamma_hw]} MHz")
print(f"Ratio gamma[0]/gamma[2]: {gamma_hw[0]/gamma_hw[2]:.1f}x")

# ================================================================
# Load IBM measurements
# ================================================================
ibm_data = {}
for label, fn in [
    ('selective', 'sacrifice_zone_hw_selective_dd_20260324_191523.json'),
    ('uniform', 'sacrifice_zone_hw_uniform_dd_20260324_191614.json'),
    ('no_dd', 'sacrifice_zone_hw_no_dd_20260324_191713.json'),
]:
    with open(os.path.join(SCRIPT_DIR,
              f'../data/ibm_sacrifice_zone_march2026/{fn}')) as f:
        ibm_data[label] = json.load(f)['results']

print(f"\nIBM measured SumMI:")
print(f"  {'t(us)':>6s}  {'Selective':>10s}  {'Uniform':>10s}  {'No DD':>10s}  {'Sel/Uni':>8s}")
print(f"  {'-'*50}")
for k in range(5):
    t = ibm_data['selective']['time_points'][k]
    s = ibm_data['selective']['sum_mi'][k]
    u = ibm_data['uniform']['sum_mi'][k]
    n = ibm_data['no_dd']['sum_mi'][k]
    ratio = s/u if u > 0 else 0
    print(f"  {t:6.1f}  {s:10.4f}  {u:10.4f}  {n:10.4f}  {ratio:8.2f}x")

# ================================================================
# Corrected gamma profiles (Echo = palindromic mirror in time)
# ================================================================
print(f"\n{'=' * 65}")
print("CORRECTED: Echo-aware gamma profiles")
print("=" * 65)

T2_echo = [5.22, 122.70, 243.85, 169.97, 237.57]  # T2 with echo (us)
gamma_T2star = [1.0/t for t in T2star]   # no DD
gamma_T2echo = [1.0/t for t in T2_echo]  # with DD

# Three profiles matching IBM configurations:
gamma_selective = [gamma_T2star[0]] + [gamma_T2echo[i] for i in range(1, N)]
gamma_uniform = [gamma_T2echo[i] for i in range(N)]
gamma_nodd = [gamma_T2star[i] for i in range(N)]

print(f"\n  Profile          Q85      Q86      Q87      Q88      Q94")
print(f"  {'Selective DD':15s}  {gamma_selective[0]:.4f}   {gamma_selective[1]:.4f}   {gamma_selective[2]:.4f}   {gamma_selective[3]:.4f}   {gamma_selective[4]:.4f}")
print(f"  {'Uniform DD':15s}  {gamma_uniform[0]:.4f}   {gamma_uniform[1]:.4f}   {gamma_uniform[2]:.4f}   {gamma_uniform[3]:.4f}   {gamma_uniform[4]:.4f}")
print(f"  {'No DD':15s}  {gamma_nodd[0]:.4f}   {gamma_nodd[1]:.4f}   {gamma_nodd[2]:.4f}   {gamma_nodd[3]:.4f}   {gamma_nodd[4]:.4f}")

print(f"\n  Sacrifice qubit difference: {gamma_selective[0]/gamma_uniform[0]:.2f}x "
      f"(selective vs uniform)")

# ================================================================
# Simulation with corrected profiles
# ================================================================
print(f"\n{'=' * 65}")
print("SIMULATION: Selective vs Uniform vs No-DD")
print("=" * 65)

# Initial state: |+>^N
d = 2**N
plus = (np.array([1, 1], dtype=complex) / np.sqrt(2))
psi0 = plus.copy()
for _ in range(N - 1):
    psi0 = np.kron(psi0, plus)
rho0 = np.outer(psi0, psi0.conj())

# Sweep J to find the best match
print(f"\nSweeping J to find best match with IBM data:")
print(f"  {'J (MHz)':>8s}  {'J/g_ctr':>8s}  {'Sel(t=1)':>9s}  {'Uni(t=1)':>9s}  "
      f"{'Ratio(1)':>8s}  {'Ratio(5)':>8s}  {'IBM r(1)':>8s}  {'IBM r(5)':>8s}")
print(f"  {'-'*75}")

ibm_ratio_1 = ibm_data['selective']['sum_mi'][0] / ibm_data['uniform']['sum_mi'][0]
ibm_ratio_5 = ibm_data['selective']['sum_mi'][4] / ibm_data['uniform']['sum_mi'][4]

best_J = None
best_err = np.inf

for J in [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]:
    H = build_heisenberg_chain(N, J=J)

    results = {}
    for label, gammas in [('sel', gamma_selective),
                          ('uni', gamma_uniform),
                          ('ndd', gamma_nodd)]:
        L = build_liouvillian(H, gammas)
        rho_vec = rho0.flatten()
        results[label] = {}
        for t_us in [1.0, 5.0]:
            rho_t_vec = expm(L * t_us) @ rho_vec
            rho_t = rho_t_vec.reshape(d, d)
            results[label][t_us] = sum_mi(rho_t, N)

    r1 = results['sel'][1.0] / results['uni'][1.0] if results['uni'][1.0] > 1e-10 else 0
    r5 = results['sel'][5.0] / results['uni'][5.0] if results['uni'][5.0] > 1e-10 else 0

    err = abs(r1 - ibm_ratio_1) + abs(r5 - ibm_ratio_5)
    if err < best_err and r1 > 0 and r5 > 0:
        best_err = err
        best_J = J

    jg = J / np.mean(gamma_selective[1:])
    print(f"  {J:8.2f}  {jg:8.1f}  {results['sel'][1.0]:9.4f}  {results['uni'][1.0]:9.4f}  "
          f"{r1:8.2f}x  {r5:8.2f}x  {ibm_ratio_1:8.2f}x  {ibm_ratio_5:8.2f}x")

print(f"\n  Best J: {best_J} MHz")

# Full time sweep at best J
if best_J:
    print(f"\n{'=' * 65}")
    print(f"FULL TIME SWEEP at J = {best_J} MHz")
    print("=" * 65)

    H = build_heisenberg_chain(N, J=best_J)

    print(f"\n  {'t(us)':>6s}  {'Sel(sim)':>9s}  {'Sel(IBM)':>9s}  "
          f"{'Uni(sim)':>9s}  {'Uni(IBM)':>9s}  {'r(sim)':>7s}  {'r(IBM)':>7s}")
    print(f"  {'-'*65}")

    for k in range(5):
        t_us = ibm_data['selective']['time_points'][k]
        smi_ibm_sel = ibm_data['selective']['sum_mi'][k]
        smi_ibm_uni = ibm_data['uniform']['sum_mi'][k]

        results_t = {}
        for label, gammas in [('sel', gamma_selective), ('uni', gamma_uniform)]:
            L = build_liouvillian(H, gammas)
            rho_t_vec = expm(L * t_us) @ rho0.flatten()
            rho_t = rho_t_vec.reshape(d, d)
            results_t[label] = sum_mi(rho_t, N)

        r_sim = results_t['sel'] / results_t['uni'] if results_t['uni'] > 1e-10 else 0
        r_ibm = smi_ibm_sel / smi_ibm_uni if smi_ibm_uni > 0 else 0

        print(f"  {t_us:6.1f}  {results_t['sel']:9.4f}  {smi_ibm_sel:9.4f}  "
              f"{results_t['uni']:9.4f}  {smi_ibm_uni:9.4f}  {r_sim:7.2f}x  {r_ibm:7.2f}x")
