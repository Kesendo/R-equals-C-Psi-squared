"""
Verify the structural decomposition of {L_H, L_Dc}
====================================================

PROOF_CROSS_TERM_FORMULA gives ||{L_H, L_Dc}||^2 = 4*gamma^2*(N-2)*||L_H||^2
via a bond-sum rule (overlap vanishes) + spectator variance (non-overlap remains).

This script makes the decomposition explicit at N=3 by splitting the
anticommutator per (bond, site) pair and verifying:

  (a) for k in bond, {L_H_bond, L_Dc^(k)} = 0 (overlap, vanishes by bond-sum rule)
  (b) for k not in bond, {L_H_bond, L_Dc^(k)} != 0 (spectator, generically nonzero)

This is the structural reason why Pythagoras holds at N=2 (no spectators)
but breaks at N>=3, while [L, Pi^2] = 0 holds at all N (per-term symmetry,
no global cross-term cancellation needed).

Date: 2026-04-15
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)


def kron(*args):
    out = args[0]
    for a in args[1:]:
        out = np.kron(out, a)
    return out


def bond_H(N, i, j, J=1.0):
    """XX+YY bond Hamiltonian on sites (i,j) of an N-qubit system."""
    ops_x = [I2]*N; ops_x[i]=X; ops_x[j]=X
    ops_y = [I2]*N; ops_y[i]=Y; ops_y[j]=Y
    return 0.5 * J * (kron(*ops_x) + kron(*ops_y))


def L_H_super(H):
    """Hamiltonian superoperator L_H[rho] = -i[H, rho], in vec convention."""
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    return -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))


def L_Dc_site(N, k, gamma=0.1):
    """Centered single-site dissipator L_Dc^(k) = gamma * (Z_k otimes Z_k) as superop.

    Total L_Dc = sum_k L_Dc^(k) for uniform Z-dephasing.
    Diagonal in Pauli basis with entries +gamma (P_k in {I,Z}) or -gamma (P_k in {X,Y}).
    """
    ops = [I2]*N; ops[k] = Z
    Z_k = kron(*ops)
    return gamma * np.kron(Z_k.conj(), Z_k)


def anticommutator(A, B):
    return A @ B + B @ A


def fnorm(M):
    return np.linalg.norm(M)


# ============================================================
# Test the structural decomposition at N=3
# ============================================================
N = 3
gamma = 0.1
print("=" * 64)
print(f"Decomposition of {{L_H, L_Dc}} at N={N}")
print("=" * 64)
print()
print("For each (bond, site) pair, compute ||{L_H_bond, L_Dc^(site)}||")
print("Expectation: vanishes when site is on the bond (overlap),")
print("             nonzero when site is a spectator.")
print()

# Bonds in N=3 chain: <0,1> and <1,2>
bonds = [(0, 1), (1, 2)]
sites = list(range(N))

print(f"  {'bond':>6} {'site':>5} {'on bond?':>10} | {'||{LH_bond, LDc_site}||':>26}")
print("  " + "-" * 55)

for (i, j) in bonds:
    H_bond = bond_H(N, i, j)
    LH_bond = L_H_super(H_bond)
    for k in sites:
        LDc_k = L_Dc_site(N, k, gamma)
        anti = anticommutator(LH_bond, LDc_k)
        norm = fnorm(anti)
        on_bond = k in (i, j)
        print(f"  ({i},{j}) {k:>5} {'YES' if on_bond else 'no':>10} | {norm:>26.6e}")

print()

# Now sum over sites for each bond, then over all bonds
print("Per-bond sums (summing over all 3 sites):")
print(f"  {'bond':>6} {'sum over k in bond':>22} {'sum over k not in bond':>26}")
print("  " + "-" * 55)
for (i, j) in bonds:
    H_bond = bond_H(N, i, j)
    LH_bond = L_H_super(H_bond)
    on_bond_sum = np.zeros_like(LH_bond)
    off_bond_sum = np.zeros_like(LH_bond)
    for k in sites:
        LDc_k = L_Dc_site(N, k, gamma)
        if k in (i, j):
            on_bond_sum += anticommutator(LH_bond, LDc_k)
        else:
            off_bond_sum += anticommutator(LH_bond, LDc_k)
    print(f"  ({i},{j}) {fnorm(on_bond_sum):>22.6e} {fnorm(off_bond_sum):>26.6e}")



# ============================================================
# Compare to total ||{L_H, L_Dc}|| and verify formula
# ============================================================
print()
print("Total cross-term verification:")
H_total = sum(bond_H(N, i, j) for (i, j) in bonds)
LH_total = L_H_super(H_total)
LDc_total = sum(L_Dc_site(N, k, gamma) for k in sites)
cross = anticommutator(LH_total, LDc_total)
cross_norm_sq = fnorm(cross) ** 2
LH_norm_sq = fnorm(LH_total) ** 2

predicted = 4 * gamma**2 * (N - 2) * LH_norm_sq
print(f"  ||{{L_H, L_Dc}}||^2 measured  = {cross_norm_sq:.6e}")
print(f"  4*gamma^2*(N-2)*||L_H||^2  = {predicted:.6e}")
print(f"  ratio                     = {cross_norm_sq / predicted:.10f}")
print()

# ============================================================
# Repeat at N=4 to confirm spectator scaling
# ============================================================
print()
print("=" * 64)
print("Same decomposition at N=4 (more spectators per bond)")
print("=" * 64)

N = 4
bonds_4 = [(0, 1), (1, 2), (2, 3)]
sites_4 = list(range(N))

print(f"  {'bond':>6}  {'on-bond contribution':>22}  {'off-bond contribution':>23}")
print("  " + "-" * 55)
for (i, j) in bonds_4:
    H_bond = bond_H(N, i, j)
    LH_bond = L_H_super(H_bond)
    on_bond_sum = np.zeros_like(LH_bond)
    off_bond_sum = np.zeros_like(LH_bond)
    for k in sites_4:
        LDc_k = L_Dc_site(N, k, gamma)
        if k in (i, j):
            on_bond_sum += anticommutator(LH_bond, LDc_k)
        else:
            off_bond_sum += anticommutator(LH_bond, LDc_k)
    print(f"  ({i},{j})  {fnorm(on_bond_sum):>22.6e}  {fnorm(off_bond_sum):>23.6e}")

H_total = sum(bond_H(N, i, j) for (i, j) in bonds_4)
LH_total = L_H_super(H_total)
LDc_total = sum(L_Dc_site(N, k, gamma) for k in sites_4)
cross_norm_sq = fnorm(anticommutator(LH_total, LDc_total)) ** 2
LH_norm_sq = fnorm(LH_total) ** 2
predicted = 4 * gamma**2 * (N - 2) * LH_norm_sq

print()
print(f"  ||{{L_H, L_Dc}}||^2  = {cross_norm_sq:.6e}")
print(f"  predicted          = {predicted:.6e}")
print(f"  ratio              = {cross_norm_sq / predicted:.10f}")

# ============================================================
# Summary
# ============================================================
print()
print("=" * 64)
print("CONCLUSION")
print("=" * 64)
print()
print("On-bond contributions vanish (bond-sum rule).")
print("Off-bond (spectator) contributions are nonzero (disjoint factor).")
print("Pythagoras N=2: no spectators, only on-bond, vanishes exactly.")
print("Pythagoras N>=3: spectators exist, off-bond non-cancelling, breaks.")
print()
print("Contrast with [L, Pi^2] = 0:")
print("  Each bond H AND each site dissipator individually commutes with X^N.")
print("  Per-term symmetry. No global cancellation needed. Holds for all N.")
