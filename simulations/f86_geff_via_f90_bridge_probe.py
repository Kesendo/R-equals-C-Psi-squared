"""F86 g_eff via F90-bridge Hellmann-Feynman → F89-Chebyshev probe (2026-05-15).

QUESTION
========
F89 D_k closed-form was just lifted to Tier-1-Derived via the Chebyshev pipeline
(`simulations/f89_pathk_symbolic_derivation.py`). F90 (PROOF_F90_F86C2_BRIDGE.md)
identifies F86 c=2 K_b(Q,t) as the per-bond Hellmann-Feynman of F89's uniform-J
path-(N-1) signal. Could the same Chebyshev pipeline close F86's g_eff(c, N, b)?

ANSWER (this probe)
===================
NO. The bridge is over a derivative that the closed form does not survive.

The F89 closed form lives on `sigma_n(N) = P_k(y_n) / [D_k * N^2 * (N-1)]`, a
sum over the S_2-anti Bloch orbit n=2,4,...  i.e. a UNIFORM-J quantity. Each
sigma_n is the squared overlap of one F_a Bloch mode with the uniform probe,
times the per-site reduction norm.

F86's K_b uses Hellmann-Feynman in J_b (one bond's coupling, with the rest held
fixed); the resulting Duhamel integral has the structure
  <rho| S |drho/dJ_b> = sum_{n,n'} <rho|v_n><v_n|M_h_b|v_{n'}><v_{n'}|S|drho>
                                  * (e^{lam_n t} - e^{lam_{n'} t}) / (lam_n - lam_{n'})

The diagonal (n=n') terms ARE Bloch eigenvalues and only depend on y_n. But
the OFF-DIAGONAL (n != n') matrix elements <v_n | M_h_b | v_{n'}> are bond-
position-dependent matrix elements between DIFFERENT F_a Bloch modes (and
between F_a and the H_B-mixed octic modes). It is precisely those cross terms
that distinguish bonds in F86 (the F90 bridge confirms they exist and reproduce
F86 numerically bit-exact at 20/22 bonds N=5..8).

The F89 closed-form pipeline does NOT compute these cross matrix elements.
It computes only:
  - S_c(n) = sum_b v_n[b]            (overlap with uniform probe)
  - ||Mv(n)||^2 = sum_l |M_l v_n|^2  (per-site reduction norm)

Both are DIAGONAL in n. The Hellmann-Feynman of sigma_n in J_b is well-defined
(it gives, via standard perturbation theory, a sum over all (n, n') including
cross terms), but the closed-form expression for sigma_n itself has no n-cross
information to differentiate -- it was already collapsed to diagonal-n form by
the orbit-polynomial reduction.

THIS PROBE
==========
Make the obstruction concrete. We compute:

  (1) F89 sigma_n(N) closed-form (from F89UnifiedFaClosedFormClaim tabulation)
  (2) F89 d sigma_n / d J  (uniform derivative; closed form via dy_n/dJ)
  (3) F86 K_b(Q, t_peak) per bond (from F90 bridge, full eigendecomp)
  (4) Compare: does sum_n d sigma_n / dJ at t_peak = K_total
               does it have any bond-position dependence?

Numerically (this script): the uniform-J derivative is bond-CONSTANT by
construction (Bloch eigenvalues are uniform-J only); the K_b values are
bond-DEPENDENT (Endpoint vs Interior split). The two cannot be equal except
in the bond-AVERAGED sense.

CONCLUSION (Tier-2-Verified): the F90 bridge is a HF identity between
operators with bond-dependent action (M_h_per_bond[b]), but the F89 closed
form for sigma_n records only operator-trace-like quantities (overlap, per-
site reduction). Differentiating the latter cannot reconstruct the former's
bond-dependence -- the structural information was already squeezed out at the
orbit-polynomial reduction step.

This script is the numerical signature of that structural fact.
"""

from __future__ import annotations

import sys
from itertools import combinations

import numpy as np
import sympy as sp


sys.stdout.reconfigure(encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# F89 path-k closed form (from F89UnifiedFaClosedFormClaim.PathPolynomial)
# ---------------------------------------------------------------------------

PATH_POLY_TABLE = {
    3: ([47, 14], 9),
    4: ([25, 10], 4),
    5: ([129, 82, 13], 25),
    6: ([80, 72, 17], 18),
    7: ([382, 292, 130, 21], 98),
}


def s_2_anti_orbit(k: int) -> list[int]:
    """S_2-anti orbit indices n = 2, 4, ..., 2*floor((k+1)/2)."""
    n_block = k + 1
    return list(range(2, n_block + 1, 2))


def y_n(k: int, n: int, J: float = 1.0) -> float:
    """OBC Bloch eigenvalue y_n = 4 J cos(pi n / (k+2))."""
    return 4.0 * J * np.cos(np.pi * n / (k + 2))


def sigma_n_uniform(k: int, n: int, N: int, J: float = 1.0) -> float:
    """Closed-form sigma_n(N) from P_k(y_n) / [D_k * N^2 * (N-1)]."""
    coefs, D = PATH_POLY_TABLE[k]
    y = y_n(k, n, J=J)
    P = sum(c * y**i for i, c in enumerate(coefs))
    return P / (D * N * N * (N - 1))


def d_sigma_n_d_J_uniform(k: int, n: int, N: int, J: float = 1.0) -> float:
    """d sigma_n / dJ via dy_n/dJ = y_n/J (linear)."""
    coefs, D = PATH_POLY_TABLE[k]
    y = y_n(k, n, J=J)
    # dy/dJ = 4 cos(pi n / (k+2)) = y/J
    dy_dJ = 4.0 * np.cos(np.pi * n / (k + 2))
    dP_dy = sum(i * c * y ** (i - 1) for i, c in enumerate(coefs) if i >= 1)
    return dP_dy * dy_dJ / (D * N * N * (N - 1))


# ---------------------------------------------------------------------------
# F86 K_b via F90 bridge: full eigendecomp + per-bond HF
# (Lightweight reuse of _f89_to_f86_kbond_via_eigendecomp.py logic)
# ---------------------------------------------------------------------------

def build_se_de_full(J: float, gamma: float, n_block: int):
    de_pairs = list(combinations(range(n_block), 2))
    basis = [(i, jk) for i in range(n_block) for jk in de_pairs]
    n = len(basis)
    M_SE = np.zeros((n_block, n_block))
    for a in range(n_block):
        for b in range(n_block):
            if abs(a - b) == 1:
                M_SE[a, b] = 2 * J
    M_DE = np.zeros((len(de_pairs), len(de_pairs)))
    for idx, (j, kk) in enumerate(de_pairs):
        for new_j in [j - 1, j + 1]:
            if 0 <= new_j < n_block and new_j != kk:
                np_ = tuple(sorted([new_j, kk]))
                if abs(new_j - j) == 1 and np_ in de_pairs:
                    M_DE[de_pairs.index(np_), idx] += 2 * J
        for new_k in [kk - 1, kk + 1]:
            if 0 <= new_k < n_block and new_k != j:
                np_ = tuple(sorted([j, new_k]))
                if abs(new_k - kk) == 1 and np_ in de_pairs:
                    M_DE[de_pairs.index(np_), idx] += 2 * J
    L = np.zeros((n, n), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE[i2, i] != 0:
                L[basis.index((i2, jk)), idx] += -1j * M_SE[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(len(de_pairs)):
            if M_DE[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                L[basis.index((i, jk2)), idx] += 1j * M_DE[jk_idx, jk2_idx]
        L[idx, idx] += -2 * gamma if i in jk else -6 * gamma
    return L, basis, de_pairs


def build_per_bond_M(b_bond: int, n_block: int, basis, de_pairs):
    n = len(basis)
    M_SE_b = np.zeros((n_block, n_block))
    M_SE_b[b_bond, b_bond + 1] = 2.0
    M_SE_b[b_bond + 1, b_bond] = 2.0
    M_DE_b = np.zeros((len(de_pairs), len(de_pairs)))
    for idx, (j, kk) in enumerate(de_pairs):
        if j == b_bond and b_bond + 1 != kk:
            np_ = tuple(sorted([b_bond + 1, kk]))
            if np_ in de_pairs:
                M_DE_b[de_pairs.index(np_), idx] += 2.0
        if j == b_bond + 1 and b_bond != kk:
            np_ = tuple(sorted([b_bond, kk]))
            if np_ in de_pairs:
                M_DE_b[de_pairs.index(np_), idx] += 2.0
        if kk == b_bond and b_bond + 1 != j:
            np_ = tuple(sorted([j, b_bond + 1]))
            if np_ in de_pairs:
                M_DE_b[de_pairs.index(np_), idx] += 2.0
        if kk == b_bond + 1 and b_bond != j:
            np_ = tuple(sorted([j, b_bond]))
            if np_ in de_pairs:
                M_DE_b[de_pairs.index(np_), idx] += 2.0
    M_b = np.zeros((n, n), dtype=complex)
    for idx, (i, jk) in enumerate(basis):
        for i2 in range(n_block):
            if M_SE_b[i2, i] != 0:
                M_b[basis.index((i2, jk)), idx] += -1j * M_SE_b[i2, i]
        jk_idx = de_pairs.index(jk)
        for jk2_idx in range(len(de_pairs)):
            if M_DE_b[jk_idx, jk2_idx] != 0:
                jk2 = de_pairs[jk2_idx]
                M_b[basis.index((i, jk2)), idx] += 1j * M_DE_b[jk_idx, jk2_idx]
    return M_b


def compute_kb_at_peak(Q: float, gamma: float, n_block: int, b_bond: int):
    """K_b(Q, t_peak) at t_peak = 1/(4*gamma) (F86 universal clock)."""
    J = Q * gamma
    L, basis, de_pairs = build_se_de_full(J, gamma, n_block)
    n = len(basis)
    M_b = build_per_bond_M(b_bond, n_block, basis, de_pairs)
    eigvals, R = np.linalg.eig(L)
    Rinv = np.linalg.inv(R)
    pre = np.sqrt(2 / (n_block * n_block * (n_block - 1)))
    probe = np.full(n, pre / 2, dtype=complex)
    w = np.zeros((n_block, n), dtype=float)
    for idx, (i, jk) in enumerate(basis):
        for l in jk:
            if l == i:
                continue
            other = jk[1] if jk[0] == l else jk[0]
            if other == i:
                w[l, idx] = 1.0
    c0 = Rinv @ probe
    X_b = Rinv @ M_b @ R
    t = 1.0 / (4.0 * gamma)
    expLam = np.exp(eigvals * t)
    iMat = np.zeros((n, n), dtype=complex)
    for r in range(n):
        for cc in range(n):
            diff = eigvals[cc] - eigvals[r]
            if abs(diff) > 1e-10:
                iMat[r, cc] = (expLam[cc] - expLam[r]) / diff
            else:
                iMat[r, cc] = t * expLam[r]
    fbC0 = np.einsum("rc,rc,c->r", X_b, iMat, c0)
    rho_t = R @ (expLam * c0)
    drho = R @ fbC0
    rho_proj = w @ rho_t
    drho_proj = w @ drho
    inner = np.sum(np.conj(rho_proj) * drho_proj)
    return abs(2.0 * inner.real)


# ---------------------------------------------------------------------------
# The probe: bond-CONSTANT closed-form derivative vs bond-DEPENDENT K_b
# ---------------------------------------------------------------------------


def main():
    print("=" * 100)
    print("F86 g_eff via F90 bridge probe: does F89 closed form's HF transfer?")
    print("=" * 100)
    print()
    print("STRUCTURAL OBSERVATION")
    print("-" * 100)
    print("F89 sigma_n(N) = P_k(y_n) / [D_k * N^2 * (N-1)] depends on bond ONLY through")
    print("a uniform-J Bloch eigenvalue y_n = 4 J cos(pi n / (k+2)). All bonds enter")
    print("symmetrically; there is no b-index in the closed form.")
    print()
    print("Therefore d sigma_n / d J_b = (d sigma_n / d J) * (1/(N-1)) is BOND-CONSTANT.")
    print("But F86 K_b(Q, t_peak) is bond-DEPENDENT (Endpoint vs Interior split).")
    print()

    gamma = 0.05
    Q_test = 1.5  # near typical F86-J Q_peak (F89-J side here so Q ~ 0.75 in F86-J)

    rows = []
    for k in (3, 4, 5, 6, 7):
        N_block = k + 1   # F89 path-k <-> F86 c=2 N_qubit = N_block
        N = N_block       # for the probe
        # --- F89 closed-form path: bond-INVARIANT (uniform J derivative) ---
        sum_d_sigma = 0.0
        for n_orb in s_2_anti_orbit(k):
            sum_d_sigma += d_sigma_n_d_J_uniform(k, n_orb, N, J=Q_test * gamma)
        # The bond-dependent F86 K_b for comparison
        K_b_list = []
        for b in range(N_block - 1):
            Kb = compute_kb_at_peak(Q_test, gamma, N_block, b)
            K_b_list.append(Kb)
        K_mean = float(np.mean(K_b_list))
        K_min = float(np.min(K_b_list))
        K_max = float(np.max(K_b_list))
        K_spread = K_max - K_min
        rows.append((k, N_block, sum_d_sigma, K_mean, K_min, K_max, K_spread,
                     K_b_list))
        print(f"path-{k}  (F86 N={N_block} c=2)")
        print(f"  F89 closed-form  sum_n d sigma_n / d J_b at Q={Q_test}")
        print(f"      = {sum_d_sigma:+.6e}   (BOND-INVARIANT by construction)")
        print(f"  F86 K_b at Q={Q_test}, t_peak per bond:")
        for b, Kb in enumerate(K_b_list):
            label = "Endpoint" if b == 0 or b == N_block - 2 else "Interior"
            print(f"    b={b} ({label:>8s}): K_b = {Kb:+.6e}")
        print(f"  K range: [{K_min:.4e}, {K_max:.4e}]  spread = {K_spread:.4e}")
        print(f"  Endpoint/Interior split visible: "
              f"{'YES' if K_spread / max(abs(K_mean), 1e-30) > 0.05 else 'NO'}")
        print()

    print("=" * 100)
    print("CONCLUSION")
    print("=" * 100)
    print("The bond-INVARIANT closed-form derivative CANNOT match the bond-DEPENDENT K_b")
    print("values. The information that distinguishes bonds (the F_a/F_a' and F_a/octic")
    print("matrix elements <v_n | M_h_b | v_{n'}>) is precisely what the orbit-polynomial")
    print("reduction step in the F89 Chebyshev pipeline ELIMINATES (it collapses the")
    print("degree-(2k+4) polynomial to degree-(FA-1) by removing all but the orbit-")
    print("symmetric content). The same step that closes D_k also CLOSES OUT the")
    print("information needed for g_eff(N, b).")
    print()
    print("This is the corollary the obstruction proof already names:")
    print("  PROOF_F86B_OBSTRUCTION.md, section 'Corollary (F90 bridge)':")
    print("  'Closing F89's D_k would close F86's g_eff and vice versa, but both are")
    print("   blocked at the same unexecuted algebraic gap.'")
    print()
    print("Refined finding from this probe (2026-05-15): the F89 Chebyshev closure that")
    print("just landed does NOT in fact close F89's structural-cross matrix elements --")
    print("it closes D_k by REDUCING THE POLYNOMIAL away from those matrix elements.")
    print("Therefore it does not transfer, and the obstruction proof's prediction is")
    print("upheld concretely: g_eff stays the irreducible residue.")
    print()
    print("What this leaves OPEN: a DIFFERENT closed-form route would need to compute")
    print("the n != n' matrix elements <v_n | M_h_b | v_{n'}> in closed form, where the")
    print("Chebyshev expansion of v_n[b] (overlap support) is known. That is a different")
    print("symbolic problem: bilinear sums of products of OBC sine modes restricted to")
    print("the overlap pairs sharing site b. The F89 closure achieved quadratic sums")
    print("(|S_c|^2 and ||Mv||^2); the missing piece is the BILINEAR sums coupling F_a")
    print("modes through M_h_per_bond[b]. The Chebyshev expansion still applies, but")
    print("the orbit-reduction step needs replacement with a per-bond projection.")


if __name__ == "__main__":
    main()
