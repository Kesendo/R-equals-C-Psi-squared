"""
Check 4: Rebound-Mechanism Check for Nested Mirror Hypothesis
==============================================================

Scientific question: The hypothesis claims the mirror modes (Re = -gamma)
drive the non-Markovian rebound. Test by:
  (1) Decompose rho(0) in the eigenmode basis of L
  (2) Identify which class carries the initial amplitude
  (3) Construct modified initial states that exclude mirror-class modes
  (4) Evolve and verify rebound vanishes without mirror modes

Caveat: arbitrary superpositions of eigenmodes need not be valid density
matrices (positive semidefinite, trace 1). If projections produce
non-physical states, use nearest physical state.

Authors: Tom and Claude (Code)
Date: 2026-04-14
"""

import numpy as np
from scipy.linalg import expm
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

np.set_printoptions(precision=6, suppress=True)

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
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jumps:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * (np.kron(Idd, LdL) + np.kron(LdL.T, Idd)))
    return L


def partial_trace_B(M):
    r = M.reshape(2, 2, 2, 2)
    return np.einsum('ibjb->ij', r)


def vec_to_rho(v):
    return v.reshape((4, 4), order='F')


def rho_to_vec(rho):
    return rho.reshape(-1, order='F')


def is_valid_density_matrix(rho, tol=1e-8):
    """Check if rho is a valid density matrix."""
    # Hermitian
    if np.max(np.abs(rho - rho.conj().T)) > tol:
        return False, "not Hermitian"
    # Trace 1
    if abs(np.trace(rho) - 1.0) > tol:
        return False, f"trace = {np.trace(rho).real:.6f}"
    # Positive semidefinite
    eigvals = np.linalg.eigvalsh(rho)
    if np.min(eigvals) < -tol:
        return False, f"min eigenvalue = {np.min(eigvals):.6f}"
    return True, "valid"


def nearest_density_matrix(rho):
    """Project onto nearest valid density matrix."""
    # Hermitianize
    rho = 0.5 * (rho + rho.conj().T)
    # Eigendecompose, clip negatives
    eigvals, eigvecs = np.linalg.eigh(rho)
    eigvals = np.maximum(eigvals, 0)
    # Renormalize
    if np.sum(eigvals) > 0:
        eigvals /= np.sum(eigvals)
    rho_proj = eigvecs @ np.diag(eigvals) @ eigvecs.conj().T
    return rho_proj


# Parameters
J = 1.0
gamma_B = 0.1

H = J * 0.5 * (kron(X, X) + kron(Y, Y))
L_jump = np.sqrt(gamma_B) * kron(I2, Z)
L_super = liouvillian(H, [L_jump])

# Diagonalize L
eigvals, R = np.linalg.eig(L_super)  # R: right eigenvectors
L_inv = np.linalg.inv(R)  # L_inv @ R = I; rows of L_inv are left eigenvectors

# Classify modes
conserved_idx = []
mirror_idx = []
correlation_idx = []

for i in range(16):
    re = eigvals[i].real
    if abs(re) < 1e-6:
        conserved_idx.append(i)
    elif abs(re + 2 * gamma_B) < 1e-6:
        correlation_idx.append(i)
    else:
        mirror_idx.append(i)

print("=" * 72)
print("Check 4: Rebound-Mechanism Check")
print("=" * 72)
print(f"Mode classification:")
print(f"  Conserved (Re=0):    {len(conserved_idx)} modes, indices {conserved_idx}")
print(f"  Mirror (Re=-g):      {len(mirror_idx)} modes, indices {mirror_idx}")
print(f"  Correlation (Re=-2g): {len(correlation_idx)} modes, indices {correlation_idx}")
print()

# Initial state
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_mix = 0.5 * np.eye(2, dtype=complex)
rho0 = np.kron(rho_plus, rho_mix)
v0 = rho_to_vec(rho0)

# Decompose rho0 in eigenmode basis
coeffs = L_inv @ v0  # coefficients c_k such that v0 = sum_k c_k * R[:,k]

print("=" * 72)
print("Step 1: Decompose rho(0) in eigenmode basis")
print("=" * 72)

total_weight = {}
for cls_name, indices in [("conserved", conserved_idx),
                           ("mirror", mirror_idx),
                           ("correlation", correlation_idx)]:
    class_weight = np.sum(np.abs(coeffs[indices]) ** 2)
    total_weight[cls_name] = class_weight
    print(f"  {cls_name:>12}: weight = {class_weight:.6f}  "
          f"({len(indices)} modes)")
    for i in indices:
        print(f"    mode {i:2d}: |c|={abs(coeffs[i]):.6f}, "
              f"Re(lam)={eigvals[i].real:.5f}, Im(lam)={eigvals[i].imag:.5f}")

total = sum(total_weight.values())
print(f"\n  Total weight: {total:.6f}")
print(f"  Mirror fraction: {total_weight['mirror']/total:.4f}")
print()

# =========================================================================
# Step 2: Construct class-projected initial states
# =========================================================================
print("=" * 72)
print("Step 2: Construct class-projected initial states")
print("=" * 72)

projections = {}
for cls_name, keep_indices in [("conserved_only", conserved_idx),
                                ("mirror_only", mirror_idx),
                                ("correlation_only", correlation_idx),
                                ("no_mirror", conserved_idx + correlation_idx),
                                ("no_conserved", mirror_idx + correlation_idx),
                                ("no_correlation", conserved_idx + mirror_idx)]:
    # Zero out coefficients for modes NOT in keep_indices
    c_proj = np.zeros_like(coeffs)
    for i in keep_indices:
        c_proj[i] = coeffs[i]

    # Reconstruct density matrix
    v_proj = R @ c_proj
    rho_proj = vec_to_rho(v_proj)

    valid, reason = is_valid_density_matrix(rho_proj)
    if not valid:
        rho_nearest = nearest_density_matrix(rho_proj)
        valid2, reason2 = is_valid_density_matrix(rho_nearest)
        projections[cls_name] = {
            'rho': rho_nearest,
            'physical': valid2,
            'note': f"projected: {reason}; nearest DM: {reason2}"
        }
    else:
        projections[cls_name] = {
            'rho': rho_proj,
            'physical': True,
            'note': "valid DM from projection"
        }

    print(f"\n  {cls_name}:")
    print(f"    Direct projection valid: {valid} ({reason})")
    if not valid:
        print(f"    Nearest DM valid: {projections[cls_name]['physical']} ({projections[cls_name]['note']})")
    print(f"    Tr(rho): {np.trace(projections[cls_name]['rho']).real:.6f}")
    rho_S = partial_trace_B(projections[cls_name]['rho'])
    print(f"    rho_S |01| = {abs(rho_S[0,1]):.6f}")

# =========================================================================
# Step 3: Evolve each projection and measure rebound
# =========================================================================
print()
print("=" * 72)
print("Step 3: Time evolution and rebound measurement")
print("=" * 72)

ts = np.linspace(0, 30, 601)
results = {}

# Also evolve original rho0 for reference
test_states = {'original': rho0}
for name, data in projections.items():
    if data['physical']:
        test_states[name] = data['rho']

fig, ax = plt.subplots(figsize=(12, 7))

for name, rho_init in test_states.items():
    v_init = rho_to_vec(rho_init)
    coh = []
    for t in ts:
        vt = expm(L_super * t) @ v_init
        rho_t = vec_to_rho(vt)
        rho_S = partial_trace_B(rho_t)
        coh.append(abs(rho_S[0, 1]))

    # Detect rebound
    min_coh = coh[0]
    rebound_amp = 0
    rebound_time = 0
    for i, c in enumerate(coh[1:], 1):
        if c < min_coh:
            min_coh = c
        elif c > min_coh + 0.002:
            if c - min_coh > rebound_amp:
                rebound_amp = c - min_coh
                rebound_time = ts[i]

    has_rebound = rebound_amp > 0.005
    results[name] = {
        'initial_coh': coh[0],
        'min_coh': min(coh),
        'rebound_amp': rebound_amp,
        'rebound_time': rebound_time,
        'has_rebound': has_rebound,
        'coh_trace': coh,
    }

    print(f"\n  {name}:")
    print(f"    Initial |S01|: {coh[0]:.6f}")
    print(f"    Minimum |S01|: {min(coh):.6f}")
    print(f"    Rebound amplitude: {rebound_amp:.6f}")
    print(f"    Has rebound: {has_rebound}")

    # Plot
    style = '-' if name == 'original' else '--'
    lw = 2.5 if name in ('original', 'mirror_only', 'no_mirror') else 1.0
    ax.plot(ts, coh, style, linewidth=lw, label=name, alpha=0.8)

ax.set_xlabel('t')
ax.set_ylabel('|rho_S_{01}|(t)')
ax.set_title('Coherence rebound by eigenmode class projection')
ax.legend(fontsize=8)
ax.set_ylim(-0.01, 0.55)
ax.grid(True, alpha=0.3)

results_dir = Path("simulations/results/qubit_in_qubit_mode_projection")
plt.tight_layout()
plt.savefig(results_dir / 'rebound_by_class.png', dpi=150)
print(f"\nPlot saved to {results_dir / 'rebound_by_class.png'}")

# =========================================================================
# Step 4: Verdict
# =========================================================================
print()
print("=" * 72)
print("VERDICT")
print("=" * 72)

verdict_lines = []

orig = results.get('original', {})
mirror_only = results.get('mirror_only', {})
no_mirror = results.get('no_mirror', {})

if 'original' in results:
    verdict_lines.append(f"Original: rebound={orig['has_rebound']}, amp={orig['rebound_amp']:.4f}")

if 'mirror_only' in results:
    verdict_lines.append(f"Mirror-only: rebound={mirror_only['has_rebound']}, amp={mirror_only['rebound_amp']:.4f}")

if 'no_mirror' in results:
    verdict_lines.append(f"No-mirror: rebound={no_mirror['has_rebound']}, amp={no_mirror['rebound_amp']:.4f}")

# Key test: does removing mirror modes eliminate the rebound?
if 'no_mirror' in results and 'original' in results:
    if not no_mirror['has_rebound'] and orig['has_rebound']:
        verdict_lines.append("\nMECHANISM CONFIRMED: removing mirror modes eliminates the rebound.")
        verdict_lines.append("The mirror-class eigenmodes are necessary for non-Markovian dynamics.")
    elif no_mirror['has_rebound']:
        verdict_lines.append(f"\nMECHANISM UNCERTAIN: rebound persists without mirror modes (amp={no_mirror['rebound_amp']:.4f}).")
        verdict_lines.append("Mirror modes may contribute but are not the sole source.")
    else:
        verdict_lines.append("\nCANNOT DETERMINE: original has no rebound or projection failed.")

# Does mirror-only state produce rebound?
if 'mirror_only' in results:
    if mirror_only['has_rebound']:
        verdict_lines.append(f"Mirror-only state produces rebound (amp={mirror_only['rebound_amp']:.4f}): mirror modes ARE sufficient.")
    else:
        verdict_lines.append("Mirror-only state has no rebound: mirror modes alone are NOT sufficient.")

# Check correlation between class weight and rebound
if len(results) > 2:
    verdict_lines.append("\nRebound amplitude vs initial mirror-class weight:")
    for name, res in sorted(results.items()):
        if name in projections:
            rho_init = projections[name]['rho']
            v_init = rho_to_vec(rho_init)
            c_init = L_inv @ v_init
            mirror_weight = np.sum(np.abs(c_init[mirror_idx]) ** 2)
            verdict_lines.append(f"  {name:>20s}: mirror_weight={mirror_weight:.4f}, rebound_amp={res['rebound_amp']:.4f}")

for line in verdict_lines:
    print(line)

# Save
with open(results_dir / 'mode_projection_results.txt', 'w', encoding='utf-8') as f:
    f.write("Check 4: Rebound-Mechanism Check Results\n")
    f.write("=" * 72 + "\n\n")
    for name, res in results.items():
        f.write(f"{name}: initial_coh={res['initial_coh']:.6f}, "
                f"min_coh={res['min_coh']:.6f}, "
                f"rebound_amp={res['rebound_amp']:.6f}, "
                f"has_rebound={res['has_rebound']}\n")
    f.write("\n" + "\n".join(verdict_lines))
print(f"\nResults saved to {results_dir / 'mode_projection_results.txt'}")
