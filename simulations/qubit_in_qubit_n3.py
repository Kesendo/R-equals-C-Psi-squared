"""
Check 2: N=3 Scaling Check for Nested Mirror Hypothesis
========================================================

Scientific question: The hypothesis predicts (N+1) eigenvalue classes
at Re(lambda_k) = -k * 2*gamma_outer / N for N-layer nests. At N=3
this means 4 classes at 0, -2g/3, -4g/3, -2g.

Setup: Three-qubit chain S - M - B with:
  H = J * 0.5 * [(X_S X_M + Y_S Y_M) + (X_M X_B + Y_M Y_B)]
  Single jump: sqrt(gamma_B) * I (x) I (x) Z
  Parameters: J = 1.0, gamma_B = 0.1
  Liouville dim: 4^3 = 64 (L is 64x64)

Tests:
  (1) How many distinct Re(lambda) classes exist?
  (2) Do they fall at 0, -2g/3, -4g/3, -2g?
  (3) Partial-trace weights per class

Authors: Tom and Claude (Code)
Date: 2026-04-14
"""

import numpy as np
from scipy.linalg import expm
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


def partial_trace_single(M_8x8, keep, n_qubits=3):
    """Partial trace of a 3-qubit operator, keeping one qubit.

    M_8x8: reshaped eigenmode as 8x8 matrix
    keep: which qubit to keep (0=S, 1=M, 2=B)
    Returns: 2x2 reduced operator
    """
    # Reshape to (2,2,2, 2,2,2) for bra and ket indices
    M = M_8x8.reshape([2] * (2 * n_qubits))
    # Trace over all qubits except 'keep'
    trace_over = [i for i in range(n_qubits) if i != keep]
    for i in sorted(trace_over, reverse=True):
        M = np.trace(M, axis1=i, axis2=i + n_qubits)
        # Adjust n_qubits for remaining traces
        n_qubits -= 1
        # Find new position of 'keep' index
    return M


def partial_trace_general(M_reshaped, keep_qubit, total_qubits):
    """General partial trace: trace out all qubits except keep_qubit."""
    d = 2 ** total_qubits
    M = M_reshaped.reshape([2] * (2 * total_qubits))

    # We need to trace pairs of indices (i, i+total_qubits) for all i != keep_qubit
    # Work with einsum for clarity
    trace_indices = [i for i in range(total_qubits) if i != keep_qubit]

    # Build einsum string
    # Input indices: row indices (0..N-1), col indices (N..2N-1)
    # For traced qubits: row index = col index
    # For kept qubit: row and col are free
    input_idx = list(range(2 * total_qubits))
    for t in trace_indices:
        input_idx[t + total_qubits] = input_idx[t]  # tie col to row

    # Output indices: just the kept qubit's row and col
    output_idx = [keep_qubit, keep_qubit + total_qubits]

    # Adjust for einsum: use letters
    letters = 'abcdefghijklmnopqrstuvwxyz'
    in_str = ''.join(letters[i] for i in input_idx)
    out_str = ''.join(letters[i] for i in output_idx)

    return np.einsum(f'{in_str}->{out_str}', M)


# Parameters
J = 1.0
gamma_B = 0.1

# Three-qubit chain Hamiltonian: S-M and M-B nearest-neighbor XX+YY
H = (J * 0.5 * (kron(X, X, I2) + kron(Y, Y, I2)    # S-M coupling
               + kron(I2, X, X) + kron(I2, Y, Y)))   # M-B coupling

# Single jump: dephasing on B only
L_jump = np.sqrt(gamma_B) * kron(I2, I2, Z)

# Build Liouvillian (64x64)
L_super = liouvillian(H, [L_jump])

print("=" * 72)
print("Check 2: N=3 Scaling Check")
print("=" * 72)
print(f"System: 3-qubit chain S-M-B, XX+YY coupling, J={J}, gamma_B={gamma_B}")
print(f"Liouvillian dimension: {L_super.shape[0]}x{L_super.shape[1]}")
print(f"Hypothesis prediction: 4 classes at Re = 0, -2g/3, -4g/3, -2g")
print(f"  = 0, {-2*gamma_B/3:.6f}, {-4*gamma_B/3:.6f}, {-2*gamma_B:.6f}")
print()

# Diagonalize
eigvals, eigvecs = np.linalg.eig(L_super)

# Sort by Re then Im
order = sorted(range(64), key=lambda i: (round(eigvals[i].real, 8),
                                          round(eigvals[i].imag, 8)))

# =========================================================================
# Test 1: How many distinct Re(lambda) classes?
# =========================================================================
print("=" * 72)
print("Test 1: Eigenvalue class count")
print("=" * 72)

re_vals = [eigvals[i].real for i in order]
re_sorted = sorted(re_vals)

# Cluster Re values
tol = 1e-5
classes = []
current_class = [re_sorted[0]]
for r in re_sorted[1:]:
    if abs(r - current_class[-1]) < tol:
        current_class.append(r)
    else:
        classes.append(current_class)
        current_class = [r]
classes.append(current_class)

print(f"Number of distinct Re(lambda) classes: {len(classes)}")
print()
print(f"{'class':>5} {'Re(lambda)':>12} {'count':>6} {'predicted':>12}")
print("-" * 45)

predicted = [0.0, -2*gamma_B/3, -4*gamma_B/3, -2*gamma_B]
predicted_sorted = sorted(predicted)

for i, cls in enumerate(classes):
    re_mean = np.mean(cls)
    count = len(cls)
    # Find closest prediction
    closest_pred = min(predicted_sorted, key=lambda p: abs(p - re_mean))
    match = "MATCH" if abs(re_mean - closest_pred) < 1e-4 else "NO MATCH"
    print(f"{i:5d} {re_mean:12.6f} {count:6d} {closest_pred:12.6f}  {match}")

print()

# =========================================================================
# Test 2: Do classes match the predicted spacing?
# =========================================================================
print("=" * 72)
print("Test 2: Spacing analysis")
print("=" * 72)

class_centers = [np.mean(c) for c in classes]
print(f"Class centers: {[f'{c:.6f}' for c in class_centers]}")
print(f"Predicted:     {[f'{p:.6f}' for p in sorted(predicted)]}")
print()

# Check against prediction
if len(classes) == 4:
    errors = []
    for c in class_centers:
        closest = min(predicted, key=lambda p: abs(p - c))
        errors.append(abs(c - closest))
    max_error = max(errors)
    print(f"Max deviation from prediction: {max_error:.2e}")
    if max_error < 1e-4:
        print("RESULT: N=3 prediction CONFIRMED. 4 classes at predicted positions.")
    else:
        print(f"RESULT: 4 classes found but positions deviate. Max error: {max_error:.6f}")
else:
    print(f"RESULT: {len(classes)} classes found, hypothesis predicted 4.")
    if len(classes) > 4:
        print("  More classes than predicted: hypothesis over-simplifies N=3.")
    else:
        print("  Fewer classes than predicted: some predicted classes may be absent.")

    # Check which predictions are matched
    for p in sorted(predicted):
        matched = False
        for c in class_centers:
            if abs(c - p) < 1e-4:
                matched = True
                break
        status = "found" if matched else "MISSING"
        print(f"  Predicted Re={p:.6f}: {status}")

print()

# =========================================================================
# Test 3: Degeneracy pattern and partial-trace weights
# =========================================================================
print("=" * 72)
print("Test 3: Partial-trace weights per eigenvalue class")
print("=" * 72)

# For each eigenmode, compute partial-trace weights on S, M, B
print(f"  {'Re(lam)':>9} {'Im(lam)':>10} {'|M_S|':>7} {'|M_M|':>7} {'|M_B|':>7}")
print("  " + "-" * 55)

class_weights = {}  # re_class -> list of (nS, nM, nB)

for idx in order:
    lam = eigvals[idx]
    vec = eigvecs[:, idx]
    M = vec.reshape((8, 8), order='F')
    nTotal = np.linalg.norm(M)
    if nTotal < 1e-15:
        continue

    # Partial traces
    M_S = partial_trace_general(M, 0, 3)
    M_M = partial_trace_general(M, 1, 3)
    M_B = partial_trace_general(M, 2, 3)

    nS = np.linalg.norm(M_S) / nTotal
    nM = np.linalg.norm(M_M) / nTotal
    nB = np.linalg.norm(M_B) / nTotal

    # Classify
    re_class = round(lam.real, 5)
    if re_class not in class_weights:
        class_weights[re_class] = []
    class_weights[re_class].append((nS, nM, nB))

    # Print first few and last few per class
    print(f"  {lam.real:9.5f} {lam.imag:10.5f} {nS:7.4f} {nM:7.4f} {nB:7.4f}")

print()
print("=" * 72)
print("Summary: average partial-trace weights per class")
print("=" * 72)
print(f"{'Re(lambda)':>12} {'count':>6} {'mean|M_S|':>10} {'mean|M_M|':>10} {'mean|M_B|':>10}")
print("-" * 55)

for re_class in sorted(class_weights.keys()):
    weights = class_weights[re_class]
    n = len(weights)
    mean_S = np.mean([w[0] for w in weights])
    mean_M = np.mean([w[1] for w in weights])
    mean_B = np.mean([w[2] for w in weights])
    print(f"{re_class:12.5f} {n:6d} {mean_S:10.4f} {mean_M:10.4f} {mean_B:10.4f}")

# =========================================================================
# Non-Markovian rebound check for N=3
# =========================================================================
print()
print("=" * 72)
print("Bonus: Non-Markovian rebound in rho_S for N=3")
print("=" * 72)

# Initial state: S in |+>, M and B maximally mixed
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_mix = 0.5 * np.eye(2, dtype=complex)
rho0 = np.kron(rho_plus, np.kron(rho_mix, rho_mix))

ts = np.linspace(0, 30, 301)
coh_S = []
for t in ts:
    v = rho0.reshape(-1, order='F')
    vt = expm(L_super * t) @ v
    rho_t = vt.reshape((8, 8), order='F')
    rho_S = partial_trace_general(rho_t, 0, 3)
    coh_S.append(abs(rho_S[0, 1]))

# Check for rebound
max_after_min = 0
min_val = coh_S[0]
for c in coh_S[1:]:
    if c < min_val:
        min_val = c
    elif c > min_val + 0.01:
        max_after_min = max(max_after_min, c)

print(f"Initial |S01|: {coh_S[0]:.4f}")
print(f"Minimum |S01|: {min(coh_S):.4f} at t={ts[np.argmin(coh_S)]:.1f}")
if max_after_min > 0:
    print(f"Max rebound:   {max_after_min:.4f} (non-Markovian rebound confirmed)")
else:
    print("No significant rebound detected.")

# =========================================================================
# Save results
# =========================================================================
results_dir = Path("simulations/results/qubit_in_qubit_n3")

log_lines = []
log_lines.append("Check 2: N=3 Scaling Check Results")
log_lines.append("=" * 72)
log_lines.append(f"System: 3-qubit chain S-M-B, J={J}, gamma_B={gamma_B}")
log_lines.append(f"Liouvillian dim: {L_super.shape[0]}")
log_lines.append(f"Number of distinct Re(lambda) classes: {len(classes)}")
log_lines.append(f"Class centers: {[f'{np.mean(c):.6f}' for c in classes]}")
log_lines.append(f"Class sizes:   {[len(c) for c in classes]}")
log_lines.append(f"Predicted:     {[f'{p:.6f}' for p in sorted(predicted)]}")
log_lines.append("")
log_lines.append("Per-class partial-trace weights:")
for re_class in sorted(class_weights.keys()):
    weights = class_weights[re_class]
    n = len(weights)
    mean_S = np.mean([w[0] for w in weights])
    mean_M = np.mean([w[1] for w in weights])
    mean_B = np.mean([w[2] for w in weights])
    log_lines.append(f"  Re={re_class:10.5f}: count={n:3d}, "
                     f"|M_S|={mean_S:.4f}, |M_M|={mean_M:.4f}, |M_B|={mean_B:.4f}")

with open(results_dir / 'eigenvalue_analysis.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(log_lines))

print(f"\nResults saved to {results_dir / 'eigenvalue_analysis.txt'}")
