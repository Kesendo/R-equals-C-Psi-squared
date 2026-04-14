"""
Check 1: SWAP-Artifact Check for Nested Mirror Hypothesis
==========================================================

Scientific question: In the minimal two-qubit nest, eigenmodes of L
appear to be simultaneous SWAP eigenmodes with expectation +-1. Because
[L, SWAP] != 0 and L has degenerate eigenvalues (3/10/3), this could be
a genuine inter-layer mirror symmetry OR a basis artifact (numpy is free
to pick any basis within a degenerate eigenspace).

Method: Break the degeneracy by adding small direct dephasing on S.
Sweep gamma_S_rel = gamma_S / gamma_B from 0.0 to 0.1 in 21 steps.
Track SWAP expectation values per eigenmode as degeneracy lifts.

Interpretation:
  - Sharp survival (SWAP stays within 0.05 of +-1): genuine structure
  - Continuous smear (drifts away): numpy basis artifact
  - Mixed: document which modes survive and which smear

Authors: Tom and Claude (Code)
Date: 2026-04-14
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

np.set_printoptions(precision=6, suppress=True)

# Pauli basis
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def liouvillian(H, jumps):
    """Vectorized Liouvillian using column-stacking vec convention."""
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def partial_trace_B(M_4x4):
    """Trace out B (second qubit)."""
    r = M_4x4.reshape(2, 2, 2, 2)
    return np.einsum('ibjb->ij', r)


def partial_trace_S(M_4x4):
    """Trace out S (first qubit)."""
    r = M_4x4.reshape(2, 2, 2, 2)
    return np.einsum('aiaj->ij', r)


# Parameters
J = 1.0
gamma_B = 0.1
gamma_S_rel_values = np.linspace(0.0, 0.1, 21)

# SWAP operator on 2 qubits (Hilbert space)
SWAP_hilb = np.zeros((4, 4), dtype=complex)
for i in range(2):
    for j in range(2):
        SWAP_hilb[2 * j + i, 2 * i + j] = 1

# SWAP superoperator
SWAP_super = np.kron(SWAP_hilb.T, SWAP_hilb)

# Hamiltonian
H = J * 0.5 * (kron(X, X) + kron(Y, Y))

# Storage for tracking
results_dir = Path("simulations/results/nested_mirror_swap_check")
all_swap_data = []

log_lines = []
log_lines.append("Check 1: SWAP-Artifact Check")
log_lines.append("=" * 72)
log_lines.append(f"Parameters: J={J}, gamma_B={gamma_B}")
log_lines.append(f"Sweep: gamma_S_rel = gamma_S/gamma_B from 0.0 to 0.1, 21 points")
log_lines.append("")

print("Check 1: SWAP-Artifact Check")
print("=" * 72)

for idx_g, gamma_S_rel in enumerate(gamma_S_rel_values):
    gamma_S = gamma_S_rel * gamma_B

    # Jump operators
    jumps = [np.sqrt(gamma_B) * kron(I2, Z)]
    if gamma_S > 0:
        jumps.append(np.sqrt(gamma_S) * kron(Z, I2))

    L = liouvillian(H, jumps)

    # Diagonalize
    eigvals, eigvecs = np.linalg.eig(L)

    # Sort by (Re, Im)
    order = sorted(range(16), key=lambda i: (round(eigvals[i].real, 8),
                                              round(eigvals[i].imag, 8)))

    # Compute SWAP expectation for each mode
    swap_exps = []
    re_vals = []
    im_vals = []

    for idx in order:
        lam = eigvals[idx]
        v = eigvecs[:, idx]
        v = v / np.linalg.norm(v)
        swap_exp = (v.conj() @ SWAP_super @ v).real
        swap_exps.append(swap_exp)
        re_vals.append(lam.real)
        im_vals.append(lam.imag)

    all_swap_data.append({
        'gamma_S_rel': gamma_S_rel,
        'swap_exps': swap_exps,
        're_vals': re_vals,
        'im_vals': im_vals,
    })

    # Check degeneracy: count distinct Re(lambda) values
    re_rounded = np.round(re_vals, 6)
    unique_re = np.unique(re_rounded)

    header = f"gamma_S_rel = {gamma_S_rel:.4f} (gamma_S = {gamma_S:.5f})"
    log_lines.append(header)
    log_lines.append(f"  Distinct Re(lambda) classes: {len(unique_re)}")
    log_lines.append(f"  Re(lambda) values: {unique_re}")

    if idx_g == 0 or idx_g == len(gamma_S_rel_values) - 1 or idx_g == 10:
        print(f"\n{header}")
        print(f"  {'Re(lam)':>9} {'Im(lam)':>10} {'<SWAP>':>10}")
        print("  " + "-" * 35)

    for k in range(16):
        line = f"  {re_vals[k]:9.5f} {im_vals[k]:10.5f} {swap_exps[k]:+8.4f}"
        log_lines.append(line)
        if idx_g == 0 or idx_g == len(gamma_S_rel_values) - 1 or idx_g == 10:
            print(line)

    log_lines.append("")

# Analysis: track how SWAP expectations evolve
print("\n" + "=" * 72)
print("Summary: SWAP expectation value statistics across sweep")
print("=" * 72)

log_lines.append("\n" + "=" * 72)
log_lines.append("Summary: SWAP expectation value statistics across sweep")
log_lines.append("=" * 72)

# For each gamma_S_rel, compute: how many modes have |<SWAP>| > 0.95
# and what is the range of SWAP values
summary_header = f"{'gamma_S_rel':>12} {'|<SWAP>|>0.95':>15} {'min|<SWAP>|':>13} {'max|<SWAP>|':>13} {'mean|<SWAP>|':>13}"
print(summary_header)
log_lines.append(summary_header)

for data in all_swap_data:
    abs_swaps = np.abs(data['swap_exps'])
    n_sharp = np.sum(abs_swaps > 0.95)
    line = f"{data['gamma_S_rel']:12.4f} {n_sharp:15d} {np.min(abs_swaps):13.6f} {np.max(abs_swaps):13.6f} {np.mean(abs_swaps):13.6f}"
    print(line)
    log_lines.append(line)

# Detailed tracking: for each mode, track SWAP vs gamma_S_rel
# To correlate modes across the sweep, sort by (Re, Im) at each point
# and track the trajectory
print("\n" + "=" * 72)
print("Per-mode SWAP trajectory (sorted by Re then Im at each point)")
print("=" * 72)

log_lines.append("\n" + "=" * 72)
log_lines.append("Per-mode SWAP trajectory (sorted by Re then Im at each point)")
log_lines.append("=" * 72)

mode_header = f"{'gamma_S_rel':>12}" + "".join(f"{'mode'+str(i):>9}" for i in range(16))
print(mode_header)
log_lines.append(mode_header)

for data in all_swap_data:
    vals = [f"{s:+8.4f}" for s in data['swap_exps']]
    line = f"{data['gamma_S_rel']:12.4f} " + " ".join(vals)
    print(line)
    log_lines.append(line)

# Plot: SWAP expectation values vs gamma_S_rel
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left panel: all 16 modes
ax = axes[0]
gamma_vals = [d['gamma_S_rel'] for d in all_swap_data]
for mode_idx in range(16):
    swaps = [d['swap_exps'][mode_idx] for d in all_swap_data]
    re0 = all_swap_data[0]['re_vals'][mode_idx]
    if abs(re0) < 1e-6:
        color, label = 'blue', 'conserved (0)' if mode_idx == 13 else None
    elif abs(re0 + 2 * gamma_B) < 1e-6:
        color, label = 'red', 'correlation (-2g)' if mode_idx == 0 else None
    else:
        color, label = 'green', 'mirror (-g)' if mode_idx == 3 else None
    ax.plot(gamma_vals, swaps, 'o-', color=color, markersize=3,
            alpha=0.6, label=label)

ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.3)
ax.axhline(y=-1.0, color='gray', linestyle='--', alpha=0.3)
ax.axhline(y=0.0, color='gray', linestyle=':', alpha=0.3)
ax.set_xlabel('gamma_S / gamma_B')
ax.set_ylabel('<SWAP>')
ax.set_title('SWAP expectation per eigenmode vs asymmetric dephasing')
ax.legend(loc='center right')
ax.set_ylim(-1.3, 1.3)

# Right panel: |<SWAP>| histogram at gamma_S_rel = 0 vs 0.1
ax2 = axes[1]
swaps_0 = all_swap_data[0]['swap_exps']
swaps_end = all_swap_data[-1]['swap_exps']
bins = np.linspace(-1.1, 1.1, 45)
ax2.hist(swaps_0, bins=bins, alpha=0.5, label=f'gamma_S_rel=0.00', color='blue')
ax2.hist(swaps_end, bins=bins, alpha=0.5, label=f'gamma_S_rel=0.10', color='red')
ax2.axvline(x=1.0, color='gray', linestyle='--', alpha=0.3)
ax2.axvline(x=-1.0, color='gray', linestyle='--', alpha=0.3)
ax2.set_xlabel('<SWAP>')
ax2.set_ylabel('Count')
ax2.set_title('SWAP distribution: symmetric vs asymmetric')
ax2.legend()

plt.tight_layout()
plt.savefig(results_dir / 'swap_vs_gamma_S.png', dpi=150)
print(f"\nPlot saved to {results_dir / 'swap_vs_gamma_S.png'}")

# Verdict
print("\n" + "=" * 72)
print("VERDICT")
print("=" * 72)

# Check: at gamma_S_rel = 0.1, are SWAP values still sharp?
final_swaps = np.array(all_swap_data[-1]['swap_exps'])
n_sharp_final = np.sum(np.abs(np.abs(final_swaps) - 1.0) < 0.05)
n_total = len(final_swaps)
mean_abs_final = np.mean(np.abs(final_swaps))

# Check at gamma_S_rel = 0 (baseline)
baseline_swaps = np.array(all_swap_data[0]['swap_exps'])
n_sharp_baseline = np.sum(np.abs(np.abs(baseline_swaps) - 1.0) < 0.05)
mean_abs_baseline = np.mean(np.abs(baseline_swaps))

verdict_lines = []

# Already at baseline, note the 0.995 values
baseline_not_exact = np.sum(np.abs(np.abs(baseline_swaps) - 1.0) > 0.001)
verdict_lines.append(f"Baseline (gamma_S_rel=0): {n_sharp_baseline}/16 modes within 0.05 of +-1")
verdict_lines.append(f"  Modes NOT exactly +-1 (deviation > 0.001): {baseline_not_exact}")
verdict_lines.append(f"  Mean |<SWAP>|: {mean_abs_baseline:.6f}")
verdict_lines.append(f"Final (gamma_S_rel=0.1): {n_sharp_final}/16 modes within 0.05 of +-1")
verdict_lines.append(f"  Mean |<SWAP>|: {mean_abs_final:.6f}")

if n_sharp_final >= 14:
    verdict_lines.append("\nVERDICT: SHARP SURVIVAL. SWAP +-1 pattern is genuine structure.")
elif n_sharp_final <= 6:
    verdict_lines.append("\nVERDICT: CONTINUOUS SMEAR. SWAP +-1 pattern was basis artifact.")
else:
    verdict_lines.append("\nVERDICT: MIXED RESULT. Some modes survive, others smear.")
    # Identify which
    for k in range(16):
        re0 = all_swap_data[0]['re_vals'][k]
        swap0 = all_swap_data[0]['swap_exps'][k]
        swapf = all_swap_data[-1]['swap_exps'][k]
        if abs(re0) < 1e-6:
            cls = "conserved"
        elif abs(re0 + 2 * gamma_B) < 1e-6:
            cls = "correlation"
        else:
            cls = "mirror"
        sharp = "sharp" if abs(abs(swapf) - 1.0) < 0.05 else "SMEARED"
        verdict_lines.append(f"  mode {k:2d} ({cls:>12s}): <SWAP> {swap0:+.4f} -> {swapf:+.4f}  [{sharp}]")

for line in verdict_lines:
    print(line)
    log_lines.append(line)

# Save text log
with open(results_dir / 'swap_vs_gamma_S.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines))
print(f"\nLog saved to {results_dir / 'swap_vs_gamma_S.txt'}")
