#!/usr/bin/env python3
"""Z⊗N-Mirror Hardware-Test: Vorhersagen + Run-Spezifikation.

Plan: |+−+⟩ und |−+−⟩ auf 3 Marrakesh-Qubits präparieren, beide unter
Heisenberg evolvieren, Pauli-Tomographie, dann chain.zn_mirror_diagnostic
auf die rekonstruierten ρ_a, ρ_b anwenden. Ergebnis: Quantifizierung der
effektiven transverse-field-Komponente auf Heron r2.

Setup-Spec:
  Backend:    ibm_marrakesh (oder Kingston/Fez)
  Qubits:     3 contiguous (z.B. [48, 49, 50] wie soft_break April 26)
  Initial:    |+−+⟩ und |−+−⟩ — X-basis Néel-Pair (Z⊗N-Partner)
  Hamiltonian: H = J·Σ (X_i X_{i+1} + Y_i Y_{i+1} + Z_i Z_{i+1}) (Heisenberg)
  Parameter:  J = 1, t = 0.8, n_trotter = 3
  Tomographie: 9 Pauli-Bases auf (q0, q2), 4096 shots/basis (gleich wie soft_break)
  Total:      2 states × 9 bases = 18 circuits, ~3-4 min QPU

Vorhersagen:
  Z⊗N-violation auf clean Hardware ≈ 0.0 (Heisenberg + Z-Dephasing + T1 ALL preserve Z⊗N)
  Effektives h_x ≠ 0 → Z⊗N-violation ≈ 0.085 · h_x (linear scaling)
  Effektives h_y ≠ 0 → Z⊗N-violation ≈ 3.5 · h_y (40× stärker)

Bei Marrakesh-Readout-Niveau (~0.01 für 4096 shots):
  Eine messbare Z⊗N-violation > 0.05 deutet auf signifikantes transverse-field hin.
  Eine violation < 0.02 ist konsistent mit "rein Heisenberg + Z-Dephasing + T1".
"""
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

import warnings
warnings.simplefilter('ignore')
import framework as fw

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

N = 3
J = 1.0
t_final = 0.8
gamma_Z_eff = 0.05  # Marrakesh-calibrated (Confirmation gamma_0_marrakesh_calibration)

# Two Z⊗N-partner states in X-basis Néel form
plus  = np.array([1, 1], dtype=complex) / np.sqrt(2)
minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
psi_a = np.kron(np.kron(plus, minus), plus)   # |+−+⟩
psi_b = fw.zn_mirror_state(psi_a, N=3)          # |−+−⟩
rho_a_0 = np.outer(psi_a, psi_a.conj())
rho_b_0 = np.outer(psi_b, psi_b.conj())

# Verify partners at t=0
chain = fw.ChainSystem(N=N, J=J, gamma_0=gamma_Z_eff)
diag_t0 = chain.zn_mirror_diagnostic(rho_a_0, rho_b_0)
assert diag_t0['verdict'] == 'preserved', \
    f"Initial states should be exact Z⊗N partners; got {diag_t0}"

# Heisenberg dynamics
heisenberg = [('X','X'), ('Y','Y'), ('Z','Z')]

print("=== Z⊗N-Mirror Hardware-Test: Predictions ===\n")
print(f"N={N}, J={J}, t_final={t_final}, γ_Z_eff={gamma_Z_eff} (Marrakesh-calibrated)")
print(f"States: |+−+⟩ (rho_a) and |−+−⟩ (rho_b)")
print(f"Dynamics: Heisenberg XX+YY+ZZ\n")

# ----------------------------------------------------------------------
# Predictions for various noise/field scenarios
# ----------------------------------------------------------------------

scenarios = [
    ('clean Heisenberg + γ_Z=0.05',           {}),
    ('+ T1 = 0.005 (typical Marrakesh)',      {'T1_l': [0.005]*N}),
    ('+ T1 + ZZ-crosstalk J_zz=0.1',          {'T1_l': [0.005]*N, 'J_zz': 0.1}),
    ('+ Realistic stack: T1 + Tphi + ZZ',     {'T1_l':[0.005]*N, 'Tphi_l':[0.05]*N, 'J_zz':0.05}),
    ('+ small h_x = 0.01 (transverse leak)',  {'h_x_l':[0.01]*N}),
    ('+ medium h_x = 0.05',                   {'h_x_l':[0.05]*N}),
    ('+ small h_y = 0.01',                    {'h_y_l':[0.01]*N}),
    ('+ medium h_y = 0.05',                   {'h_y_l':[0.05]*N}),
    ('+ Mini-Magnetfeld h_z=[0.05,0.1,0.05]', {'h_z_l':[0.05, 0.1, 0.05]}),
    ('+ Realistic + small h_x = 0.005',       {'T1_l':[0.005]*N, 'Tphi_l':[0.05]*N,
                                               'J_zz':0.05, 'h_x_l':[0.005]*N}),
]

print(f"{'Scenario':<48s} {'max_viol':>14s} {'verdict':<12s} {'worst_string'}")
print("-" * 90)
for label, noise in scenarios:
    rho_a = chain.propagate_with_hardware_noise(rho_a_0, t=t_final, terms=heisenberg, **noise)
    rho_b = chain.propagate_with_hardware_noise(rho_b_0, t=t_final, terms=heisenberg, **noise)
    diag = chain.zn_mirror_diagnostic(rho_a, rho_b, tol=1e-9)
    print(f"  {label:<46s}  {diag['max_violation']:>13.4e}  "
          f"{diag['verdict']:<10s}  {diag['worst_string']}")

# ----------------------------------------------------------------------
# Expected Pauli-expectations on (q0, q2) for both states under realistic
# Marrakesh stack — these are what hardware tomography should return
# ----------------------------------------------------------------------

print("\n=== Expected Pauli-expectations on (q0, q2) under realistic Marrakesh stack ===\n")
print("(reduce 3-qubit ρ to (q0, q2) via tracing out q1; same as soft_break April 26 measurement)\n")

realistic_noise = {'T1_l':[0.005]*N, 'Tphi_l':[0.05]*N, 'J_zz':0.05}
rho_a_real = chain.propagate_with_hardware_noise(rho_a_0, t=t_final, terms=heisenberg, **realistic_noise)
rho_b_real = chain.propagate_with_hardware_noise(rho_b_0, t=t_final, terms=heisenberg, **realistic_noise)

def reduce_q0q2(rho_3q):
    rho_t = rho_3q.reshape(2,2,2,2,2,2)
    return np.einsum('ikjlkm->ijlm', rho_t).reshape(4, 4)

rho_a_2q = reduce_q0q2(rho_a_real)
rho_b_2q = reduce_q0q2(rho_b_real)

PAULI = {'I': np.eye(2, dtype=complex),
         'X': np.array([[0,1],[1,0]], dtype=complex),
         'Y': np.array([[0,-1j],[1j,0]], dtype=complex),
         'Z': np.array([[1,0],[0,-1]], dtype=complex)}

print(f"{'Pauli':<6s} {'<P>_a (rho_+−+)':>17s} {'<P>_b (rho_−+−)':>17s} "
      f"{'expected ratio':>16s} {'(-1)^n_XY':>10s}")
print("-" * 80)
for a_l in 'IXYZ':
    for b_l in 'IXYZ':
        if (a_l, b_l) == ('I', 'I'): continue
        P = np.kron(PAULI[a_l], PAULI[b_l])
        e_a = float(np.real(np.trace(rho_a_2q @ P)))
        e_b = float(np.real(np.trace(rho_b_2q @ P)))
        n_xy = (1 if a_l in 'XY' else 0) + (1 if b_l in 'XY' else 0)
        sign = (-1)**n_xy
        ratio = e_b/e_a if abs(e_a) > 1e-10 else 0
        match = '✓' if (abs(e_a)<1e-10 and abs(e_b)<1e-10) or abs(ratio - sign) < 0.05 else 'check'
        print(f"  {a_l},{b_l}  {e_a:>17.4f} {e_b:>17.4f} {ratio:>16.4f} {sign:>+10d}  {match}")

print("\nIf hardware respects Z⊗N: ⟨P⟩_b ≈ (-1)^n_XY · ⟨P⟩_a for every Pauli P.")
print("Deviations quantify transverse-field effects + readout asymmetry.\n")

# ----------------------------------------------------------------------
# Hardware Run Spec (copy-paste into AIEvolution.UI)
# ----------------------------------------------------------------------

print("=== Hardware Run Spec ===\n")
print("Backend:        ibm_marrakesh (or kingston / fez)")
print("Path:           3 contiguous qubits (e.g. [48, 49, 50] from soft_break April 26)")
print("Two states:")
print("  State A: |+−+⟩  prepare via H_q0 H_q1 X_q1 H_q1 H_q2 (or just H + X + H per qubit)")
print("           Cleaner: H on each, then Z on q1 (= H X H = Z) yields |+−+⟩")
print("  State B: |−+−⟩  prepare via H Z on q0,q2 (gives |-⟩) and H on q1 (gives |+⟩)")
print(f"Hamiltonian:    Heisenberg J=1.0, Trotter steps n=3, t={t_final}")
print(f"Tomography:     9 Pauli-bases on (q0, q2), 4096 shots/basis")
print("Total circuits: 2 states × 9 bases = 18, est. QPU ~3-4 min")
print("Output JSON keys (mirror soft_break format):")
print("  expectations.state_a.{X,X}, .{X,Y}, ..., 16 expectations")
print("  expectations.state_b.{X,X}, ..., 16 expectations")
print("  parameters.{N, J, t, n_trotter, shots}, backend, path, job_id")
