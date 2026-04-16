"""
Structure Points on the Dissipation Scale
==========================================

Question B from chat 2026-04-16: when we normalize alpha by gamma_0
and look at small-gamma limit (so Hamiltonian corrections vanish),
do certain rational fractions {0, 1/2, 3/4, 1, 5/4, 3/2, 2, ...}
appear consistently across N? Does the [0, 2*gamma_0] interval have
inherent structure points, or is the spectrum topology-specific dust?

Strategy: take N=3, 4, 5 chains, run with very small gamma_0
(0.001 vs J=1). The leading order Liouvillian eigenvalues are then
gamma_0 * f(structure) where f are pure rational fractions. Hamiltonian
corrections are O(gamma_0^2 / J), invisible at this scale.

For each N: list distinct alpha/gamma_0 values, identify rational
fractions, check which fractions are common across N.

Date: 2026-04-16
Authors: Tom and Claude (chat)
"""

import numpy as np
import sys
from scipy.linalg import eig
from fractions import Fraction

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

np.set_printoptions(precision=8, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron(*factors)


def liouvillian(H, jump_ops):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


def build_xy_chain(N, J=1.0):
    H = np.zeros((2**N, 2**N), dtype=complex)
    for i in range(N - 1):
        H += J * 0.5 * (site_op(X, i, N) @ site_op(X, i+1, N)
                        + site_op(Y, i, N) @ site_op(Y, i+1, N))
    return H


def distinct_alphas(L_super, tol=1e-8):
    eigs, _ = eig(L_super)
    alphas = -eigs.real
    sorted_alphas = np.sort(alphas)
    distinct = []
    counts = []
    for a in sorted_alphas:
        if distinct and abs(a - distinct[-1]) < tol:
            counts[-1] += 1
        else:
            distinct.append(a)
            counts.append(1)
    return distinct, counts


def best_fraction(x, max_denom=16):
    """Find simplest rational p/q with q <= max_denom that approximates x."""
    if abs(x) < 1e-9:
        return Fraction(0)
    f = Fraction(x).limit_denominator(max_denom)
    return f


# ============================================================
# Scan: N = 3, 4, 5 chains, single-site dephasing, small gamma_0
# ============================================================
J = 1.0
gamma_0 = 1e-3

per_N_results = {}

print("=" * 78)
print(f"Single-site dephasing on B (last site), J={J}, gamma_0={gamma_0}")
print("=" * 78)
print()

for N in [3, 4, 5]:
    H = build_xy_chain(N, J=J)
    L_super = liouvillian(H, [np.sqrt(gamma_0) * site_op(Z, N - 1, N)])
    distinct, counts = distinct_alphas(L_super, tol=gamma_0 * 1e-3)
    normalized = [a / gamma_0 for a in distinct]
    per_N_results[N] = (distinct, counts, normalized)

    print(f"N={N}: Liouville dim={4**N}, distinct rates={len(distinct)}")
    print(f"  {'alpha':>14} {'mult':>4} {'alpha/gamma_0':>14} {'rational':>12} {'err':>10}")
    print("  " + "-" * 60)
    for a, c, x in zip(distinct, counts, normalized):
        f = best_fraction(x, max_denom=16)
        f_val = float(f) if f != 0 else 0.0
        err = abs(x - f_val)
        rational_str = str(f) if f != 0 else "0"
        print(f"  {a:>14.10f} {c:>4d} {x:>14.6f} {rational_str:>12} {err:>10.2e}")
    print()


# ============================================================
# Cross-N comparison: which fractions appear at multiple N?
# ============================================================
print("=" * 78)
print("Cross-N comparison: rational fractions present at each N")
print("=" * 78)
print()

# Collect rational sets per N
N_to_fractions = {}
for N, (distinct, counts, normalized) in per_N_results.items():
    fracs = set()
    for x in normalized:
        f = best_fraction(x, max_denom=16)
        if f == 0 or abs(x - float(f)) < 1e-4:
            fracs.add(f)
    N_to_fractions[N] = fracs

all_fractions = sorted(set().union(*N_to_fractions.values()))

print(f"  {'fraction':>10} {'value':>10} ", end="")
for N in sorted(N_to_fractions.keys()):
    print(f"{'N='+str(N):>6}", end="")
print(f" {'in all N':>10}")
print("  " + "-" * 60)

universal = set.intersection(*N_to_fractions.values())

for f in all_fractions:
    val = float(f)
    in_all = f in universal
    print(f"  {str(f):>10} {val:>10.4f} ", end="")
    for N in sorted(N_to_fractions.keys()):
        present = "yes" if f in N_to_fractions[N] else "-"
        print(f"{present:>6}", end="")
    print(f" {'YES' if in_all else '':>10}")

print()
print(f"  Universal fractions (present at N=3, 4, 5): {sorted(universal)}")
print(f"  Total distinct fractions across all N: {len(all_fractions)}")


# ============================================================
# Sanity check: does the gamma_0 -> 0 limit produce sharper rationals?
# ============================================================
print()
print("=" * 78)
print("Sanity check: scaling with gamma_0 (N=3)")
print("=" * 78)
print()
print("If structure points are universal, alpha/gamma_0 should converge")
print("to fixed rationals as gamma_0 -> 0. Drift away as gamma_0 -> J.")
print()

print(f"  {'gamma_0':>10} {'distinct':>10} {'sample alpha/gamma_0 (sorted)':>40}")
print("  " + "-" * 70)
for g in [1e-5, 1e-3, 1e-1, 5e-1, 1.0]:
    H = build_xy_chain(3, J=J)
    L_super = liouvillian(H, [np.sqrt(g) * site_op(Z, 2, 3)])
    distinct, _ = distinct_alphas(L_super, tol=g * 1e-3)
    normalized = [a / g for a in distinct]
    sample = " ".join(f"{x:.4f}" for x in normalized[:6])
    print(f"  {g:>10.5f} {len(distinct):>10} {sample:>40}")

print()
print("  As gamma_0 grows, eigenvalues of L pick up Hamiltonian-induced")
print("  corrections. The pure rational structure survives only in the")
print("  perturbative limit gamma_0 << J.")
