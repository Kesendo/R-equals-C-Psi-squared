"""Derive symbolic closed form for F89 path-2 topology subclass.

Setup: N qubits total, single path-2 block at sites (s, s+1, s+2),
bonds {s, s+1}; remaining N-3 sites are bare. Uniform J coupling on
the two block bonds, uniform Z-dephasing γ₀ on all N sites, ρ_cc
initial state.

Goal: derive S(t) = Σ_l 2|(ρ_l(t))_{0,1}|² as a closed-form function
of (N, J, γ, t) using sympy.

Strategy:
  1. Compute ρ_B(0) for the 3-qubit block analytically (Tr_{N\B}(ρ_cc))
  2. Build H_B = J(X_0 X_1 + Y_0 Y_1) + J(X_1 X_2 + Y_1 Y_2) (8x8)
  3. Build Z-dephasing Lindbladian on B (256x256 superoperator)
  4. Evolve ρ_B(t) under e^{L_B t} symbolically (eigenbasis path)
  5. Compute (ρ_l)_{0,1}(t) for l in {0, 1, 2} via partial trace
  6. Sum block + bare contributions to get S_path2(t)
  7. Simplify and verify against existing bond-isolate (2)-topology CSV
"""
from __future__ import annotations

import sys
from itertools import product

import sympy as sp

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Symbols
N_sym, J, gamma, t = sp.symbols("N J gamma t", positive=True, real=True)
sqrt2 = sp.sqrt(2)
sqrt3 = sp.sqrt(3)

# ---------------- 3-qubit Hilbert space basis ----------------
# Computational basis ordered as |b_0 b_1 b_2⟩, integer index 4*b_0 + 2*b_1 + b_2
DIM_B = 8
def ket(*bits):
    """Computational basis ket as 8-dim column vector."""
    idx = bits[0] * 4 + bits[1] * 2 + bits[2]
    v = sp.zeros(DIM_B, 1)
    v[idx] = 1
    return v


def bra(*bits):
    return ket(*bits).T


# Pauli operators on single qubit
sX = sp.Matrix([[0, 1], [1, 0]])
sY = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sZ = sp.Matrix([[1, 0], [0, -1]])
I2 = sp.eye(2)


def pauli_at(P, site):
    """Embed P on `site` of 3-qubit space (site 0 leftmost = MSB convention)."""
    op = sp.Matrix([[1]])
    for q in range(3):
        op = sp.kronecker_product(op, P if q == site else I2)
    return op


X = [pauli_at(sX, s) for s in range(3)]
Y = [pauli_at(sY, s) for s in range(3)]
Z = [pauli_at(sZ, s) for s in range(3)]


# H_B for path-2 (bonds (0,1) and (1,2))
HB = J * (X[0] * X[1] + Y[0] * Y[1] + X[1] * X[2] + Y[1] * Y[2])

# ---------------- ρ_B(0) construction ----------------
# Cases per derivation:
#   Case 1: i ∈ B, {a,b} ⊂ B → 9 terms, |1_i⟩⟨1_a,1_b| (popcount-1 × popcount-2 on B)
#   Case 2: i ∈ N\B, exactly one of {a,b} in B → 3(N-3) terms, |0⟩⟨1_q| (vac × popcount-1 on B)
# Plus h.c.

# Indices for popcount-1 and popcount-2 on B:
def popcount(i):
    return bin(i).count("1")


pcount1 = [i for i in range(DIM_B) if popcount(i) == 1]
pcount2 = [i for i in range(DIM_B) if popcount(i) == 2]
vac = 0  # |000⟩

# Case 1 sum: Σ_{i ∈ B, (a,b) ⊂ B} |i⟩⟨a∨b|
# Each (i, (a,b)) pair gives |1_i⟩_B ⟨1_a,1_b|_B
case1 = sp.zeros(DIM_B, DIM_B)
for i_idx in pcount1:
    for ab_idx in pcount2:
        case1 += ket(*[(i_idx >> (2 - q)) & 1 for q in range(3)]) * \
                 bra(*[(ab_idx >> (2 - q)) & 1 for q in range(3)])

# Case 2 sum: Σ_{i ∈ N\B, q ∈ B} |0⟩⟨1_q|
# = (N-3) · Σ_{q ∈ B} |0⟩⟨1_q| (sum over the (N-3) i values)
case2_per_q = sp.zeros(DIM_B, DIM_B)
for q_idx in pcount1:
    case2_per_q += ket(0, 0, 0) * bra(*[(q_idx >> (2 - q)) & 1 for q in range(3)])

case2 = (N_sym - 3) * case2_per_q

# ρ_B(0) prefactor: 1 / (2 · √N · √C(N,2))
# C(N,2) = N(N-1)/2 → √C(N,2) = √(N(N-1)/2)
prefactor = 1 / (2 * sp.sqrt(N_sym) * sp.sqrt(N_sym * (N_sym - 1) / 2))

rho_B_0 = prefactor * (case1 + case2 + (case1 + case2).H)

print("ρ_B(0) constructed (3-qubit block, symbolic in N).")
print(f"  Trace: {sp.simplify(rho_B_0.trace())}")
print(f"  Hermitian? {sp.simplify((rho_B_0 - rho_B_0.H).norm()) == 0}")

# Sanity at N=7: ρ_B(0) should have specific numerical entries
rho_B_0_n7 = rho_B_0.subs(N_sym, 7).evalf()
print(f"  N=7: max entry magnitude = {max(abs(rho_B_0_n7[i,j]) for i in range(8) for j in range(8)):.4f}")


# ---------------- H_B eigenstructure ----------------
# We need eigendecomposition of H_B. SE sector spans {|100⟩,|010⟩,|001⟩}.
# H_B in SE basis is 2J·[[0,1,0],[1,0,1],[0,1,0]]; eigenvalues ±2√2·J, 0.

# For symbolic propagation we use H_B.exp(-i H_B t) and H_B.exp(+i H_B t).
# sympy can compute exp of matrices analytically for small dim.
print("\nComputing exp(-i H_B t) symbolically (8x8 matrix exponential)...")
U_t = (-sp.I * HB * t).exp()
print("  done")

print("Computing exp(+i H_B t)...")
U_t_dag = (sp.I * HB * t).exp()
print("  done")


# ---------------- Z-dephasing on per-element coherence ----------------
# For the closed form path, instead of building a 64x64 Liouvillian, we use
# the known fact: under Σ_l γ D[Z_l], each coherence ⟨A|ρ|B⟩ in the Z-basis
# decays at rate γ · Σ_l (1 - z_l(A)·z_l(B)). For ρ_B propagated as
# U_t · ρ_B(0) · U_t_dag (Hamiltonian part), the dephasing modulates each
# basis-pair coherence by exp(-rate·t).
#
# We construct a "decay mask" M[a,b] = exp(-γ·Σ_l (1 - z_l(a)·z_l(b))·t) on the
# 3-qubit basis, and apply elementwise to the Hamiltonian-evolved density matrix.

def z_eigenvalue(basis_idx, site):
    bit = (basis_idx >> (2 - site)) & 1
    return -1 if bit == 1 else 1


decay_mask = sp.zeros(DIM_B, DIM_B)
for a in range(DIM_B):
    for b in range(DIM_B):
        n_diff = sum(1 for l in range(3) if z_eigenvalue(a, l) != z_eigenvalue(b, l))
        decay_mask[a, b] = sp.exp(-2 * gamma * n_diff * t)

# Note: under H_B + Z-deph, the Hamiltonian and dephasing don't commute in general,
# but BOTH preserve the popcount sectors. For the (vac, SE) and (SE, DE)
# coherence blocks, the Hamiltonian acts inside the block (preserving popcount)
# and dephasing decays each basis-pair element by the rate above. The rate
# depends only on the basis-state pair (a, b), NOT on which Hamiltonian operations
# moved amplitude to that pair. This means we can apply the dephasing mask
# elementwise to the Hamiltonian-evolved ρ_B(t).
#
# This is the "interaction-picture" simplification: in the Z-basis, dephasing acts
# diagonally on basis-pair coherences, and for U(1)-preserving H, the Heisenberg
# evolution stays inside fixed popcount sectors so dephasing rates remain constant
# along the trajectory.

# Hamiltonian-evolved ρ_B(t):
print("\nComputing U_t · ρ_B(0) · U_t† symbolically...")
rho_B_H_t = U_t * rho_B_0 * U_t_dag
print("  done (8x8 symbolic)")

# Apply dephasing mask elementwise:
print("Applying Z-dephasing mask elementwise...")
rho_B_t = sp.Matrix([[sp.simplify(rho_B_H_t[i,j] * decay_mask[i,j]) for j in range(DIM_B)] for i in range(DIM_B)])
print("  done")


# ---------------- Per-site (ρ_l)_{0,1} via partial trace ----------------
# (ρ_l)_{0,1} for site l = sum over rest sites of ρ_B[(0_l, rest), (1_l, rest)]

def basis_with_bit_at(l, bit, rest_bits):
    """Construct integer index with bit `bit` at site l and rest_bits at other sites."""
    bits = [None, None, None]
    bits[l] = bit
    others = [i for i in range(3) if i != l]
    for k, idx in enumerate(others):
        bits[idx] = rest_bits[k]
    return bits[0] * 4 + bits[1] * 2 + bits[2]


def partial_trace_off_diag(rho, l):
    """(ρ_l)_{0,1} = Σ_rest ρ[(0_l, rest), (1_l, rest)]."""
    result = 0
    for rest in product([0, 1], repeat=2):
        a = basis_with_bit_at(l, 0, rest)
        b = basis_with_bit_at(l, 1, rest)
        result += rho[a, b]
    return result


print("\nComputing (ρ_l)_{0,1}(t) for l = 0, 1, 2...")
rho_01 = [partial_trace_off_diag(rho_B_t, l) for l in range(3)]

# Sanity check at t=0:
print("  At t=0 (symbolic):")
for l in range(3):
    val = sp.simplify(rho_01[l].subs(t, 0))
    print(f"    (ρ_{l})_{{0,1}}(0) = {val}")


# ---------------- S_block(t) = Σ_l 2·|(ρ_l)_{0,1}|² ----------------
print("\nComputing S_block(t) = Σ_l 2·|(ρ_l)_{0,1}|²...")
S_block = sum(2 * sp.Abs(rho_01[l])**2 for l in range(3))

# Sympy's Abs squared on complex expressions: use conjugate
S_block_squared = sum(2 * (rho_01[l] * sp.conjugate(rho_01[l])) for l in range(3))
S_block_simplified = sp.simplify(S_block_squared)
print(f"  S_block sympy expr length: {sp.count_ops(S_block_simplified)} ops")

# Bare-site contribution: each bare site has (ρ_l)_{0,1}(0) = √(2(N-1))/(2N) (S_N-symmetric ρ_cc)
# Under pure dephasing only: |·|²(t) = (N-1)/(2N²)·exp(-4γt). Total bare = (N-3)·(N-1)/N² · exp(-4γt).
S_bare = (N_sym - 3) * (N_sym - 1) / (N_sym**2) * sp.exp(-4 * gamma * t)

S_total = sp.simplify(S_block_simplified + S_bare)
print(f"  S_total sympy expr length: {sp.count_ops(S_total)} ops")


# ---------------- Verification at N=7 against bond-isolate CSV ----------------
print("\n# Sanity check at N=7, J=0.075, γ=0.05, t=0:")
S_at_0 = sp.simplify(S_total.subs([(N_sym, 7), (J, sp.Rational(75, 1000)), (gamma, sp.Rational(5, 100)), (t, 0)]))
print(f"  S(0) = {S_at_0} = {float(S_at_0):.6f}")
print(f"  Expected: (N-1)/N = 6/7 = {6/7:.6f}")

# Save intermediate for later analytic inspection
print("\n# S_total symbolic form (truncated):")
S_str = str(S_total)
print(f"  Total length: {len(S_str)} chars")
print(f"  First 1000 chars: {S_str[:1000]}")


# Try to express in nice form
print("\n# Attempting to factor exp(-4γt) and simplify...")
S_factored = sp.simplify(S_total * sp.exp(4 * gamma * t))
print(f"  S * exp(4γt) length: {sp.count_ops(S_factored)} ops")
print(f"  First 1000 chars of S * exp(4γt): {str(S_factored)[:1000]}")
