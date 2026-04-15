"""
Inside-Outside Correspondence: Probe 2, Inside Observer BLP Index
==================================================================

Question: an observer with access only to rho_S(t) = Tr_B rho(t).
What can this observer learn about gamma_B and J separately?

The BLP (Breuer-Laine-Piilo 2009) non-Markovianity index measures
information backflow: N_BLP = integral of d|rho_S|/dt over intervals
where the trace distance INCREASES (non-Markovian signature).

Tests:
  (a) Fix J, sweep gamma_B: does BLP scale with gamma_B?
  (b) Fix gamma_B, sweep J: does BLP depend on J?
  (c) 2D grid: is BLP = f(gamma) or f(J) or f(J/gamma)?

If BLP depends primarily on gamma: the inside observer can extract
gamma from rho_S alone (the outside leaks through).
If BLP depends primarily on J: the inside observer reads J, not gamma.
If BLP depends on J/gamma: the parameters are entangled from inside too.

Date: 2026-04-15
Authors: Tom and Claude (Code)
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


def partial_trace_B(rho_AB):
    r = rho_AB.reshape(2, 2, 2, 2)
    return np.einsum('ibjb->ij', r)


def trace_distance(rho, sigma):
    """Trace distance D(rho, sigma) = 0.5 * ||rho - sigma||_1."""
    diff = rho - sigma
    eigvals = np.linalg.eigvalsh(diff)
    return 0.5 * np.sum(np.abs(eigvals))


def compute_blp(L_super, rho0_a, rho0_b, t_max=50, n_steps=1000):
    """Compute BLP non-Markovianity index.

    Use two initial states rho0_a, rho0_b. Evolve both under L.
    Track trace distance D(rho_S_a(t), rho_S_b(t)).
    BLP = integral of dD/dt over intervals where dD/dt > 0.
    """
    ts = np.linspace(0, t_max, n_steps + 1)
    dt = ts[1] - ts[0]

    distances = []
    for t in ts:
        eLt = expm(L_super * t)

        va = eLt @ rho0_a.reshape(-1, order='F')
        rho_a = va.reshape((4, 4), order='F')
        rho_Sa = partial_trace_B(rho_a)

        vb = eLt @ rho0_b.reshape(-1, order='F')
        rho_b = vb.reshape((4, 4), order='F')
        rho_Sb = partial_trace_B(rho_b)

        distances.append(trace_distance(rho_Sa, rho_Sb))

    # BLP index: sum of positive increments
    blp = 0.0
    max_rebound = 0.0
    min_d = distances[0]
    for i in range(1, len(distances)):
        dd = distances[i] - distances[i-1]
        if dd > 0:
            blp += dd
        if distances[i] < min_d:
            min_d = distances[i]
        elif distances[i] - min_d > max_rebound:
            max_rebound = distances[i] - min_d

    return blp, max_rebound, distances


# Initial states: |+><+| (x) I/2 and |0><0| (x) I/2
rho_plus = 0.5 * np.array([[1, 1], [1, 1]], dtype=complex)
rho_zero = np.array([[1, 0], [0, 0]], dtype=complex)
rho_mix = 0.5 * I2

rho0_a = np.kron(rho_plus, rho_mix)
rho0_b = np.kron(rho_zero, rho_mix)

results_dir = Path("simulations/results/inside_observer_blp")

# =========================================================================
# Sweep A: Fix J=1.0, vary gamma_B
# =========================================================================
print("=" * 72)
print("Sweep A: Fix J=1.0, vary gamma_B")
print("=" * 72)

J_fixed = 1.0
gamma_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

print(f"\n{'gamma_B':>8} {'BLP':>10} {'max_reb':>10} {'BLP*gamma':>10} {'BLP/gamma':>10}")
print("-" * 55)

sweep_a = []
for gamma_B in gamma_values:
    H = J_fixed * 0.5 * (kron(X, X) + kron(Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, Z)])
    t_max = min(200, 30 / gamma_B)
    blp, reb, _ = compute_blp(L, rho0_a, rho0_b, t_max=t_max, n_steps=800)
    sweep_a.append({'gamma': gamma_B, 'J': J_fixed, 'blp': blp, 'reb': reb})
    print(f"{gamma_B:8.3f} {blp:10.6f} {reb:10.6f} {blp*gamma_B:10.6f} {blp/gamma_B if gamma_B>0 else 0:10.6f}")

# =========================================================================
# Sweep B: Fix gamma_B=0.1, vary J
# =========================================================================
print("\n" + "=" * 72)
print("Sweep B: Fix gamma_B=0.1, vary J")
print("=" * 72)

gamma_fixed = 0.1
J_values = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]

print(f"\n{'J':>8} {'BLP':>10} {'max_reb':>10} {'J/gamma':>8}")
print("-" * 45)

sweep_b = []
for J in J_values:
    H = J * 0.5 * (kron(X, X) + kron(Y, Y))
    L = liouvillian(H, [np.sqrt(gamma_fixed) * kron(I2, Z)])
    blp, reb, _ = compute_blp(L, rho0_a, rho0_b, t_max=50, n_steps=800)
    sweep_b.append({'gamma': gamma_fixed, 'J': J, 'blp': blp, 'reb': reb})
    print(f"{J:8.3f} {blp:10.6f} {reb:10.6f} {J/gamma_fixed:8.1f}")

# =========================================================================
# 2D Grid: is BLP = f(gamma), f(J), or f(J/gamma)?
# =========================================================================
print("\n" + "=" * 72)
print("2D Grid: BLP as function of (gamma_B, J)")
print("=" * 72)

gamma_grid = [0.01, 0.05, 0.1, 0.5, 1.0]
J_grid = [0.1, 0.5, 1.0, 5.0]

grid_blp = {}
print(f"\n{'':>10}", end="")
for J in J_grid:
    print(f" {'J='+str(J):>10}", end="")
print()
print("-" * (12 + 11 * len(J_grid)))

for gamma_B in gamma_grid:
    print(f"g={gamma_B:<8.2f}", end="")
    for J in J_grid:
        H = J * 0.5 * (kron(X, X) + kron(Y, Y))
        L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, Z)])
        t_max = min(200, 30 / gamma_B)
        blp, _, _ = compute_blp(L, rho0_a, rho0_b, t_max=t_max, n_steps=500)
        grid_blp[(gamma_B, J)] = blp
        print(f" {blp:10.6f}", end="")
    print()

# Test: BLP at constant J/gamma
print("\n" + "-" * 72)
print("Test: BLP along constant J/gamma diagonals")
print("-" * 72)

for ratio in [1.0, 10.0]:
    print(f"\nJ/gamma = {ratio}:")
    for gamma_B in [0.01, 0.05, 0.1, 0.5, 1.0]:
        J = ratio * gamma_B
        if J <= 10:
            H = J * 0.5 * (kron(X, X) + kron(Y, Y))
            L = liouvillian(H, [np.sqrt(gamma_B) * kron(I2, Z)])
            t_max = min(200, 30 / gamma_B)
            blp, _, _ = compute_blp(L, rho0_a, rho0_b, t_max=t_max, n_steps=500)
            print(f"  gamma={gamma_B:.3f}, J={J:.3f}: BLP={blp:.6f}")

# =========================================================================
# VERDICT
# =========================================================================
print("\n" + "=" * 72)
print("VERDICT: What can the inside observer learn?")
print("=" * 72)

# Check column invariance (fixed J, varying gamma)
print("\nColumn variation (fixed J, varying gamma):")
for J in J_grid:
    vals = [grid_blp[(g, J)] for g in gamma_grid]
    spread = max(vals) - min(vals)
    print(f"  J={J:<4}: BLP range = {spread:.6f}  {'small' if spread < 0.01 else 'LARGE'}")

# Check row invariance (fixed gamma, varying J)
print("\nRow variation (fixed gamma, varying J):")
for g in gamma_grid:
    vals = [grid_blp[(g, J)] for J in J_grid]
    spread = max(vals) - min(vals)
    print(f"  gamma={g:<4}: BLP range = {spread:.6f}  {'small' if spread < 0.01 else 'LARGE'}")

# Save
with open(results_dir / 'blp_results.txt', 'w', encoding='utf-8') as f:
    f.write("Probe 2: Inside Observer BLP Index\n")
    f.write("=" * 72 + "\n\n")
    f.write("Sweep A (fix J=1, vary gamma):\n")
    for d in sweep_a:
        f.write(f"  gamma={d['gamma']:.3f}: BLP={d['blp']:.6f}\n")
    f.write("\nSweep B (fix gamma=0.1, vary J):\n")
    for d in sweep_b:
        f.write(f"  J={d['J']:.3f}: BLP={d['blp']:.6f}\n")
    f.write("\n2D Grid:\n")
    for g in gamma_grid:
        for J in J_grid:
            f.write(f"  gamma={g:.2f}, J={J:.1f}: BLP={grid_blp[(g,J)]:.6f}\n")

print(f"\nResults saved to {results_dir / 'blp_results.txt'}")
