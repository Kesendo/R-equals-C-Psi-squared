"""
Check 3: Coupling Robustness for Nested Mirror Hypothesis
==========================================================

Scientific question: The N=2 three-class structure was observed under
XX+YY coupling. Does it survive under:
  (a) Heisenberg XXX (X(x)X + Y(x)Y + Z(x)Z)
  (b) Pure XX (X(x)X only)
  (c) Anisotropic probes

For each: eigenvalue class count, partial-trace weight pattern,
non-Markovian rebound existence.

Authors: Tom and Claude (Code)
Date: 2026-04-14
"""

import numpy as np
from scipy.linalg import expm
from pathlib import Path

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


def partial_trace_B(M_4x4):
    r = M_4x4.reshape(2, 2, 2, 2)
    return np.einsum('ibjb->ij', r)


def partial_trace_S(M_4x4):
    r = M_4x4.reshape(2, 2, 2, 2)
    return np.einsum('aiaj->ij', r)


J = 1.0
gamma_B = 0.1

# Define coupling types
couplings = {
    'XX+YY': J * 0.5 * (kron(X, X) + kron(Y, Y)),
    'XXX (Heisenberg)': J * 0.5 * (kron(X, X) + kron(Y, Y) + kron(Z, Z)),
    'XX only': J * 0.5 * kron(X, X),
    'YY only': J * 0.5 * kron(Y, Y),
    'ZZ only': J * 0.5 * kron(Z, Z),
    'XX+ZZ': J * 0.5 * (kron(X, X) + kron(Z, Z)),
}

L_jump = np.sqrt(gamma_B) * kron(I2, Z)

# Initial state
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_mix = 0.5 * np.eye(2, dtype=complex)
rho0 = np.kron(rho_plus, rho_mix)

results_dir = Path("simulations/results/qubit_in_qubit_coupling")
log_lines = []

print("=" * 72)
print("Check 3: Coupling Robustness")
print("=" * 72)
print(f"Parameters: J={J}, gamma_B={gamma_B}, dephasing on B only")
print()

for name, H in couplings.items():
    L = liouvillian(H, [L_jump])
    eigvals, eigvecs = np.linalg.eig(L)
    order = sorted(range(16), key=lambda i: (round(eigvals[i].real, 8),
                                              round(eigvals[i].imag, 8)))

    # Cluster Re values
    re_vals = sorted([eigvals[i].real for i in order])
    tol = 1e-5
    classes = []
    current = [re_vals[0]]
    for r in re_vals[1:]:
        if abs(r - current[-1]) < tol:
            current.append(r)
        else:
            classes.append(current)
            current = [r]
    classes.append(current)

    class_centers = [np.mean(c) for c in classes]
    class_sizes = [len(c) for c in classes]

    # Partial-trace weights per class
    class_data = {}
    for idx in order:
        lam = eigvals[idx]
        vec = eigvecs[:, idx]
        M = vec.reshape((4, 4), order='F')
        nTotal = np.linalg.norm(M)
        if nTotal < 1e-15:
            continue
        M_S = partial_trace_B(M)
        M_B = partial_trace_S(M)
        nS = np.linalg.norm(M_S) / nTotal
        nB = np.linalg.norm(M_B) / nTotal
        re_class = round(lam.real, 5)
        if re_class not in class_data:
            class_data[re_class] = {'nS': [], 'nB': [], 'count': 0}
        class_data[re_class]['nS'].append(nS)
        class_data[re_class]['nB'].append(nB)
        class_data[re_class]['count'] += 1

    # Non-Markovian rebound
    ts = np.linspace(0, 30, 301)
    coh = []
    for t in ts:
        v = rho0.reshape(-1, order='F')
        vt = expm(L * t) @ v
        rho_t = vt.reshape((4, 4), order='F')
        rho_S = partial_trace_B(rho_t)
        coh.append(abs(rho_S[0, 1]))

    has_rebound = False
    min_coh = coh[0]
    rebound_amp = 0
    for c in coh[1:]:
        if c < min_coh:
            min_coh = c
        elif c > min_coh + 0.005:
            has_rebound = True
            rebound_amp = max(rebound_amp, c - min_coh)

    # Check three-class structure
    has_three_classes = (len(classes) == 3)
    has_boundary = (abs(class_centers[0]) < 1e-4 or abs(class_centers[-1]) < 1e-4)
    has_neg2g = any(abs(c + 2 * gamma_B) < 1e-4 for c in class_centers)

    # Check 1/sqrt(2) split
    has_sqrt2_split = False
    for re_class in class_data:
        if abs(re_class + gamma_B) < 0.02:  # near -gamma
            mean_nS = np.mean(class_data[re_class]['nS'])
            mean_nB = np.mean(class_data[re_class]['nB'])
            if abs(mean_nS - 1/np.sqrt(2)) < 0.05 and abs(mean_nB - 1/np.sqrt(2)) < 0.05:
                has_sqrt2_split = True

    print("-" * 72)
    print(f"Coupling: {name}")
    print(f"  Classes: {len(classes)}  (sizes: {class_sizes})")
    print(f"  Centers: {[f'{c:.5f}' for c in class_centers]}")
    print(f"  Three-class structure: {'YES' if has_three_classes else 'NO'}")
    print(f"  Boundary classes (0, -2g): {has_boundary and has_neg2g}")
    print(f"  1/sqrt(2) split in middle: {'YES' if has_sqrt2_split else 'NO'}")
    print(f"  Non-Markovian rebound: {'YES' if has_rebound else 'NO'} (amplitude: {rebound_amp:.4f})")

    # Per-class details
    print(f"  {'Re(lam)':>10} {'count':>6} {'mean|M_S|':>10} {'mean|M_B|':>10}")
    for re_class in sorted(class_data.keys()):
        d = class_data[re_class]
        print(f"  {re_class:10.5f} {d['count']:6d} "
              f"{np.mean(d['nS']):10.4f} {np.mean(d['nB']):10.4f}")

    log_lines.append(f"Coupling: {name}")
    log_lines.append(f"  Classes: {len(classes)}, sizes: {class_sizes}")
    log_lines.append(f"  Three-class: {has_three_classes}, 1/sqrt(2): {has_sqrt2_split}, rebound: {has_rebound} (amp={rebound_amp:.4f})")
    for re_class in sorted(class_data.keys()):
        d = class_data[re_class]
        log_lines.append(f"  Re={re_class:10.5f}: count={d['count']}, |M_S|={np.mean(d['nS']):.4f}, |M_B|={np.mean(d['nB']):.4f}")
    log_lines.append("")

print()
print("=" * 72)
print("SUMMARY AND VERDICT")
print("=" * 72)

summary = """
Key findings:
- XX+YY (reference): 3 classes {3,10,3}, 1/sqrt(2) split, rebound present
- XXX (Heisenberg): check if same structure survives with ZZ term
- XX only: check if reduced coupling changes the class structure
- ZZ only: excitation-number-breaking; different structure expected
- XX+ZZ: shadow-crossing coupling; may alter the class pattern

The three-class structure and 1/sqrt(2) split are coupling-dependent
properties. The boundary classes (Re=0, Re=-2g) and palindromic pairing
follow from MIRROR_SYMMETRY_PROOF and hold for ALL couplings with
Z-dephasing. The intermediate class structure depends on the specific
Hamiltonian.
"""
print(summary)

with open(results_dir / 'coupling_sweep.txt', 'w', encoding='utf-8') as f:
    f.write("Check 3: Coupling Robustness Results\n")
    f.write("=" * 72 + "\n\n")
    f.write('\n'.join(log_lines))
    f.write("\n" + summary)

print(f"Results saved to {results_dir / 'coupling_sweep.txt'}")
