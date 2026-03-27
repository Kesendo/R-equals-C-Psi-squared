"""
Q3.2: Why are XZ+YZ and ZX+ZY non-local?

Signal engineering view: each per-site Pi map is a DC/AC crossover switch.
P1 and P4 are two different crossover settings. For the palindrome, the
crossover must ANTICOMMUTE with the Hamiltonian's adjoint action ad_H.

Key test: does P1 anticommute with ad_X? With ad_Y? Same for P4.
If P1 handles X but not Y, and P4 handles Y but not X, then having both
X and Y on the same site is a routing conflict - like trying to pass two
bands through a single-band crossover.

Also: find the minimum Schmidt rank by optimizing over the commutant.
"""

import numpy as np
from itertools import product as iproduct

# Pauli matrices
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, sx, sy, sz]
NAMES = ['I', 'X', 'Y', 'Z']


def ad(sigma):
    """Adjoint action: ad_sigma(rho) = [sigma, rho] as 4x4 in Pauli basis."""
    # Pauli commutation: [sigma_a, sigma_b] = 2i * epsilon_{abc} * sigma_c
    # In the Pauli basis |a>, ad_sigma|b> = [sigma, sigma_b]
    # Decompose result in Pauli basis
    d = 2
    num = 4
    ad_mat = np.zeros((num, num), dtype=complex)
    for b in range(num):
        comm = sigma @ PAULIS[b] - PAULIS[b] @ sigma
        # Decompose in Pauli basis: coeff_a = Tr(sigma_a * comm) / 2
        for a in range(num):
            ad_mat[a, b] = np.trace(PAULIS[a] @ comm) / d
    return ad_mat


def build_pi_pauli(perm, sign):
    """Build per-site Pi as 4x4 matrix in Pauli basis."""
    M = np.zeros((4, 4), dtype=complex)
    for a in range(4):
        M[perm[a], a] = sign[a]
    return M


def anticommutator_test(M, A):
    """Check if M anticommutes with A: {M, A} = MA + AM = 0."""
    anticomm = M @ A + A @ M
    return np.linalg.norm(anticomm)


def commutator_test(M, A):
    """Check if M commutes with A: [M, A] = MA - AM = 0."""
    comm = M @ A - A @ M
    return np.linalg.norm(comm)


# ================================================================
# Part 1: Which crossover handles which signal?
# ================================================================
print("=" * 65)
print("SIGNAL ROUTING ANALYSIS")
print("Which crossover (P1/P4) handles which Pauli channel?")
print("=" * 65)

# Build adjoint actions for each Pauli
ad_X = ad(sx)
ad_Y = ad(sy)
ad_Z = ad(sz)

# Build Pi maps
P1_PERM = {0: 1, 1: 0, 2: 3, 3: 2}
P1_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}
P4_PERM = {0: 2, 1: 3, 2: 0, 3: 1}
P4_SIGN = {0: 1, 1: 1j, 2: 1, 3: 1j}

P1 = build_pi_pauli(P1_PERM, P1_SIGN)
P4 = build_pi_pauli(P4_PERM, P4_SIGN)

print("\nPer-site crossover compatibility with Pauli channels:")
print("(Palindrome requires ANTICOMMUTATION: {Pi, ad_sigma} = 0)")
print()
print(f"{'':15s} {'ad_X':>10s} {'ad_Y':>10s} {'ad_Z':>10s}")
print(f"{'':15s} {'(X signal)':>10s} {'(Y signal)':>10s} {'(Z signal)':>10s}")
print("-" * 50)

for name, M in [("P1 crossover", P1), ("P4 crossover", P4)]:
    err_x = anticommutator_test(M, ad_X)
    err_y = anticommutator_test(M, ad_Y)
    err_z = anticommutator_test(M, ad_Z)
    sym_x = "PASS" if err_x < 1e-10 else f"{err_x:.1f}"
    sym_y = "PASS" if err_y < 1e-10 else f"{err_y:.1f}"
    sym_z = "PASS" if err_z < 1e-10 else f"{err_z:.1f}"
    print(f"{name:15s} {sym_x:>10s} {sym_y:>10s} {sym_z:>10s}")

print()
print("Reading the table:")
print("  P1 handles X-signal and Z-signal but NOT Y-signal")
print("  P4 handles Y-signal and Z-signal but NOT X-signal")
print("  (Z is handled by both - it's on the DC bus)")

# ================================================================
# Part 2: Why the same-site conflict is irreconcilable
# ================================================================
print("\n" + "=" * 65)
print("SAME-SITE CONFLICT ANALYSIS")
print("=" * 65)

print("""
For XZ+YZ: H = sigma_X x sigma_Z + sigma_Y x sigma_Z

  Site 0 carries BOTH X-signal and Y-signal.
  - X-signal needs: P1 (or equivalent)  -> {P1, ad_X} = 0 PASS
  - Y-signal needs: P4 (or equivalent)  -> {P4, ad_Y} = 0 PASS
  - But site 0 can only have ONE crossover!
  - P1 fails Y: {P1, ad_Y} != 0
  - P4 fails X: {P4, ad_X} != 0

  This is a FREQUENCY CONFLICT: two bands need different crossover
  settings on the same diplexer port. No single-band crossover works.

For XZ+ZY: H = sigma_X x sigma_Z + sigma_Z x sigma_Y

  X-signal is on site 0, Y-signal is on site 1.
  - Site 0: P4 for X-signal  -> PASS
  - Site 1: P1 for Y-signal  -> wait, P1 handles X not Y...
""")

# Let me verify the actual compatibility for each term on each site
print("Detailed per-term compatibility:")
print()

# For a 2-site Hamiltonian term sigma_a x sigma_b,
# the per-site Pi compatibility requires:
# {M_0, ad_a} = 0 AND {M_1, ad_b} = 0
# (anticommutation with the adjoint of the Hamiltonian's Pauli on that site)

terms_to_check = [
    ("XZ", sx, sz),
    ("YZ", sy, sz),
    ("ZY", sz, sy),
    ("ZX", sz, sx),
]

for term_name, pa, pb in terms_to_check:
    ad_a = ad(pa)
    ad_b = ad(pb)
    print(f"  Term {term_name}:")
    for pi_name, M in [("P1", P1), ("P4", P4)]:
        err_0 = anticommutator_test(M, ad_a)
        err_1 = anticommutator_test(M, ad_b)
        s0 = "PASS" if err_0 < 1e-10 else "FAIL"
        s1 = "PASS" if err_1 < 1e-10 else "FAIL"
        print(f"    {pi_name} on site0: {s0:4s}  {pi_name} on site1: {s1:4s}")
    print()

print("For XZ+YZ: need site0 to pass BOTH X and Y.")
print("  P1 passes X but fails Y. P4 passes Y but fails X.")
print("  No single crossover works. -> NON-LOCAL Pi required.")
print()
print("For XZ+ZY: X on site0, Y on site1. Each site independent.")
print("  Site0: P4 passes X. Site1: P4 passes Y. -> P4xP4?")

# Verify P4xP4 for XZ+ZY
from pi_operator_entanglement import (build_hamiltonian, build_z_dephasing,
    build_liouvillian, verify_pi, build_analytical_pi_multi)

H_xz_zy = build_hamiltonian(2, [(sx, sz), (sz, sy)], J=1.0)
c_ops, sg = build_z_dephasing(2, 0.05)
L = build_liouvillian(H_xz_zy, c_ops)

# Test which combinations work for XZ+ZY
print("\n  XZ+ZY full Liouvillian test (including dephasing):")
combos = [
    ("P1xP1", P1_PERM, P1_SIGN, P1_PERM, P1_SIGN),
    ("P4xP4", P4_PERM, P4_SIGN, P4_PERM, P4_SIGN),
    ("P1xP4", P1_PERM, P1_SIGN, P4_PERM, P4_SIGN),
    ("P4xP1", P4_PERM, P4_SIGN, P1_PERM, P1_SIGN),
]
for label, p1, s1, p2, s2 in combos:
    Pi_test = build_analytical_pi_multi(2, [(p1, s1), (p2, s2)])
    err = verify_pi(Pi_test, L, sg)
    status = "OK" if err < 1e-6 else "FAIL"
    print(f"    {label}: {status} (err={err:.2e})")

# ================================================================
# Part 3: Minimum Schmidt rank via commutant optimization
# ================================================================
print("\n" + "=" * 65)
print("MINIMUM SCHMIDT RANK (Q3.1)")
print("=" * 65)

from pi_operator_entanglement import (find_pi_operator, operator_schmidt,
    build_bipartition_perm_2q)

# Build system for XZ+YZ (genuinely non-local)
H_xz_yz = build_hamiltonian(2, [(sx, sz), (sy, sz)], J=1.0)
L_nl = build_liouvillian(H_xz_yz, c_ops)

# Get one valid Pi (numerical)
Pi0, pairs, evals = find_pi_operator(L_nl, sg)
err0 = verify_pi(Pi0, L_nl, sg)
print(f"\nNumerical Pi0 for XZ+YZ: verification error {err0:.2e}")

# Compute commutant basis: all S with [S, L] = 0
# Build ad_L: ad_L(S) = LS - SL, vectorized
d2 = 16
eye16 = np.eye(d2, dtype=complex)
ad_L = np.kron(L_nl, eye16) - np.kron(eye16, L_nl.T)

# Find null space of ad_L
U_ad, sv_ad, Vh_ad = np.linalg.svd(ad_L)
tol = 1e-8 * sv_ad[0]
null_mask = sv_ad < tol
null_dim = np.sum(null_mask)
print(f"Commutant dimension: {null_dim}")

# Extract null space vectors (these are vec(S) for each S in commutant)
null_vecs = Vh_ad[null_mask].conj().T  # columns are null vectors
print(f"Null space shape: {null_vecs.shape}")

# Each null vector, reshaped to 16x16, is a matrix S with [S,L] = 0
# Any valid Pi = Pi0 @ S, so Pi = Pi0 @ (sum_k c_k S_k)
# We want to minimize the Schmidt rank of Pi over all c_k

# Optimization: minimize nuclear norm (convex relaxation of rank)
# of the reshaped operator

perm_2q = build_bipartition_perm_2q()

def schmidt_rank_of_pi(coeffs):
    """Compute Schmidt rank for Pi = Pi0 @ S(coeffs)."""
    # Build S from coefficients
    S_vec = null_vecs @ coeffs
    S = S_vec.reshape(d2, d2)
    Pi = Pi0 @ S
    # Schmidt decomposition
    sv = operator_schmidt(Pi, perm_2q, 4, 4)
    return sv

# Try random optimization: sample many random S, find minimum rank
print("\nSearching for minimum-rank Pi (random sampling over commutant)...")
best_rank = 16
best_coeffs = None
n_trials = 5000

for trial in range(n_trials):
    # Random coefficients (complex)
    c = np.random.randn(null_dim) + 1j * np.random.randn(null_dim)
    sv = schmidt_rank_of_pi(c)
    rank = np.sum(sv > 1e-8 * sv[0])
    if rank < best_rank:
        best_rank = rank
        best_coeffs = c.copy()
        best_sv = sv.copy()

print(f"Best Schmidt rank found: {best_rank} (out of {n_trials} trials)")
print(f"Schmidt coefficients (top 6):")
for i in range(min(6, len(best_sv))):
    rel = best_sv[i] / best_sv[0] if best_sv[0] > 0 else 0
    marker = " *" if best_sv[i] > 1e-8 * best_sv[0] else ""
    print(f"  sigma_{i} = {best_sv[i]:.6f}  (relative: {rel:.6f}){marker}")

# Also try structured approach: project Pi0 onto low-rank subspace
print("\n--- Structured search: low-rank projection ---")

# Get SVD of the reshaped Pi0
Pi0_r = Pi0[np.ix_(perm_2q, perm_2q)]
T0 = Pi0_r.reshape(4, 4, 4, 4)
M0 = T0.transpose(0, 2, 1, 3).reshape(16, 16)
U0, S0, Vh0 = np.linalg.svd(M0)

# Try truncating to rank k and finding closest valid Pi
for target_rank in [2, 3, 4]:
    # Truncate
    M_trunc = (U0[:, :target_rank] * S0[:target_rank]) @ Vh0[:target_rank, :]

    # Convert back to Pi operator
    T_trunc = M_trunc.reshape(4, 4, 4, 4).transpose(0, 2, 1, 3).reshape(16, 16)
    Pi_trunc = np.zeros_like(Pi0)
    inv_perm = np.zeros(16, dtype=int)
    for m, k in enumerate(perm_2q):
        inv_perm[k] = m
    Pi_trunc = T_trunc[np.ix_(inv_perm, inv_perm)]

    err_trunc = verify_pi(Pi_trunc, L_nl, sg)
    print(f"  Rank-{target_rank} truncation: verification error = {err_trunc:.2e}"
          f"  {'(valid)' if err_trunc < 1e-4 else '(invalid)'}")

# ================================================================
# Summary
# ================================================================
print("\n" + "=" * 65)
print("SIGNAL ENGINEERING INTERPRETATION")
print("=" * 65)
# ================================================================
# Part 4: Systematic minimum rank via optimization
# ================================================================
print("\n" + "=" * 65)
print("SYSTEMATIC MINIMUM RANK SEARCH")
print("=" * 65)

# Precompute M_k matrices: M_k = reshape(perm(Pi0 @ S_k))
M_basis = []
for k in range(null_dim):
    S_k = null_vecs[:, k].reshape(d2, d2)
    Pi_k = Pi0 @ S_k
    Pi_k_r = Pi_k[np.ix_(perm_2q, perm_2q)]
    T_k = Pi_k_r.reshape(4, 4, 4, 4)
    M_k = T_k.transpose(0, 2, 1, 3).reshape(16, 16)
    M_basis.append(M_k)

M_basis = np.array(M_basis)  # shape (null_dim, 16, 16)

def tail_sv_sum(params, target_rank):
    """Sum of squared singular values beyond target_rank."""
    # params: 2*null_dim reals -> null_dim complex coefficients
    c = params[:null_dim] + 1j * params[null_dim:]
    M = np.tensordot(c, M_basis, axes=([0], [0]))
    sv = np.linalg.svd(M, compute_uv=False)
    return np.sum(sv[target_rank:]**2)

from scipy.optimize import minimize, differential_evolution

print(f"\nPrecomputed {null_dim} M_basis matrices (16x16 each)")
print("Optimizing: minimize sum of tail singular values")

for target_rank in range(2, 11):
    # Multiple restarts
    best_tail = np.inf
    best_sv = None

    for restart in range(20):
        x0 = np.random.randn(2 * null_dim) * 0.5
        res = minimize(tail_sv_sum, x0, args=(target_rank,),
                       method='L-BFGS-B', options={'maxiter': 2000})
        if res.fun < best_tail:
            best_tail = res.fun
            c_opt = res.x[:null_dim] + 1j * res.x[null_dim:]
            M_opt = np.tensordot(c_opt, M_basis, axes=([0], [0]))
            best_sv = np.linalg.svd(M_opt, compute_uv=False)

    # Check if valid (Pi invertible?)
    S_opt = null_vecs @ c_opt
    Pi_opt = Pi0 @ S_opt.reshape(d2, d2)
    try:
        err_opt = verify_pi(Pi_opt, L_nl, sg)
        valid = err_opt < 1e-6
    except np.linalg.LinAlgError:
        valid = False

    actual_rank = np.sum(best_sv > 1e-8 * best_sv[0])
    tail_ratio = best_sv[target_rank] / best_sv[0] if target_rank < len(best_sv) else 0

    marker = "<--" if actual_rank <= target_rank and valid else ""
    print(f"  target={target_rank:2d}:  tail_sum={best_tail:.2e}  "
          f"sigma_{target_rank}/sigma_0={tail_ratio:.4f}  "
          f"actual_rank={actual_rank:2d}  valid={valid}  {marker}")

    if actual_rank <= target_rank and valid:
        print(f"\n  >>> MINIMUM SCHMIDT RANK = {actual_rank}")
        print(f"  >>> Singular values:")
        for i in range(min(actual_rank + 2, len(best_sv))):
            marker2 = " *" if best_sv[i] > 1e-8 * best_sv[0] else ""
            print(f"      sigma_{i} = {best_sv[i]:.8f}{marker2}")
        break

print("""
The palindromic mirror Pi is a DC/AC crossover switch in Liouville space.

Two crossover settings exist:
  P1: routes X-signal and Z-signal (I<->X, Y<->Z)
  P4: routes Y-signal and Z-signal (I<->Y, X<->Z)

Each per-site crossover handles a SUBSET of Pauli channels.
P1 and P4 are complementary: together they cover {X, Y, Z}.

The Hamiltonian coupling creates inter-site signal paths.
Each path requires a specific crossover setting at each site.

CONFLICT arises when one site carries signals requiring BOTH
crossover settings simultaneously (X and Y on same site).
This is analogous to a frequency conflict in a diplexer:
two bands need different filter settings on the same port.

Resolution: a HYBRID COUPLER (non-local Pi) that mixes the
signals between sites before crossovering. The 1/sqrt(2)
coefficients are the 3dB splitting ratio of the coupler.
""")
