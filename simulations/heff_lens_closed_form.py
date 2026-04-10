"""
Phase A: Closed-form check for psi_opt via effective Hamiltonian H_eff.
Convention: ibm_april_predictions (np.kron(H, Id) for commutator).

H_eff[i,j] = -J * adjacency(i,j) - 1j * delta_{ij} * gamma_i

Authors: Thomas Wicht, Claude (Opus 4.6)
Date: April 10, 2026
"""
import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')
import json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ibm_april_predictions import heisenberg_H

# Canonical parameters from ibm_april_predictions
N = 5
J = 1.0
T2_us = np.array([5.22, 122.70, 243.85, 169.97, 237.57])
gamma_phys = 1.0 / (2.0 * T2_us)
gamma_min = 0.05
gamma_sacrifice = gamma_phys / gamma_phys.min() * gamma_min

PSI_OPT_REF = np.array([0.099, 0.239, 0.428, 0.572, 0.651])
PSI_OPT_REF = PSI_OPT_REF / np.linalg.norm(PSI_OPT_REF)

out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "results", "heff_lens")
os.makedirs(out_dir, exist_ok=True)

print("=" * 70)
print("PHASE A: H_eff closed-form check for psi_opt")
print("=" * 70)
print(f"  J = {J}")
print(f"  gamma = {np.round(gamma_sacrifice, 5)}")

# Build H_eff: -J * adjacency - i * diag(gamma)
adj = np.zeros((N, N))
for k in range(N - 1):
    adj[k, k+1] = adj[k+1, k] = 1.0

H_eff = -J * adj - 1j * np.diag(gamma_sacrifice)
print(f"\n  H_eff =")
for row in H_eff:
    print(f"    [{', '.join(f'{x.real:+.3f}{x.imag:+.3f}j' for x in row)}]")

# Diagonalize
evals, evecs = np.linalg.eig(H_eff)

# Sort by |Im(eigenvalue)| ascending (smallest imaginary part = longest-lived)
order = np.argsort(np.abs(evals.imag))
evals = evals[order]
evecs = evecs[:, order]

print(f"\n  Eigenvalues (sorted by |Im|):")
for k in range(N):
    ev = evals[k]
    print(f"    {k}: Re={ev.real:+.6f}  Im={ev.imag:+.6f}  |Im|={abs(ev.imag):.6f}")

# The slowest mode: smallest |Im|
slow_idx = 0
a_heff = evecs[:, slow_idx]
# Normalize magnitudes and fix phase (largest component real positive)
a_heff_mag = np.abs(a_heff)
a_heff_mag = a_heff_mag / np.linalg.norm(a_heff_mag)

# Also extract with consistent phase
phase = np.exp(-1j * np.angle(a_heff[np.argmax(np.abs(a_heff))]))
a_heff_phased = a_heff * phase
a_heff_phased = a_heff_phased / np.linalg.norm(a_heff_phased)

print(f"\n  Slowest eigenvector (index {slow_idx}, |Im| = {abs(evals[slow_idx].imag):.6f}):")
print(f"  {'site':>4}  {'|a_heff|':>10}  {'a_opt_ref':>10}  {'diff':>10}")
print(f"  {'-' * 38}")
for k in range(N):
    diff = abs(a_heff_mag[k] - PSI_OPT_REF[k])
    print(f"  {k:>4}  {a_heff_mag[k]:>10.6f}  {PSI_OPT_REF[k]:>10.6f}  {diff:>10.6f}")

cosine_sim = abs(a_heff_mag @ PSI_OPT_REF)
print(f"\n  Cosine similarity: {cosine_sim:.8f}")
print(f"  Max componentwise diff: {np.max(np.abs(a_heff_mag - PSI_OPT_REF)):.6f}")

if cosine_sim > 0.999:
    print(f"\n  VERDICT: CONFIRMED (cosine > 0.999)")
elif cosine_sim > 0.99:
    print(f"\n  VERDICT: CLOSE MATCH (cosine > 0.99)")
else:
    print(f"\n  VERDICT: MISMATCH (cosine = {cosine_sim:.6f})")

# Save results
spec_data = {
    'eigenvalues': [{'re': float(e.real), 'im': float(e.imag)} for e in evals],
    'eigenvectors': [[float(abs(evecs[i, k])) for i in range(N)] for k in range(N)],
    'slowest_index': int(slow_idx),
}
with open(os.path.join(out_dir, "heff_spectrum.json"), 'w') as f:
    json.dump(spec_data, f, indent=2)

comp_data = {
    'psi_opt_reference': PSI_OPT_REF.tolist(),
    'heff_candidate': a_heff_mag.tolist(),
    'componentwise_diff': (np.abs(a_heff_mag - PSI_OPT_REF)).tolist(),
    'cosine_similarity': float(cosine_sim),
}
with open(os.path.join(out_dir, "comparison.json"), 'w') as f:
    json.dump(comp_data, f, indent=2)

print(f"\n  Saved to {out_dir}/")
