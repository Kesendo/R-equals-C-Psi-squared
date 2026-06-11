#!/usr/bin/env python3
r"""The antilinear triangle: one Klein four-group behind five proofs, self-validating.

On N-qubit operator space three involutions meet the identity: the transpose
θ(A) = Aᵀ (ℂ-linear, antiautomorphism), entrywise conjugation conj(A) = Ā
(antilinear, automorphism), and the adjoint †(A) = A† (antilinear,
antiautomorphism), with † = θ∘conj = conj∘θ. Together with id they form a
Klein four-group graded by two characters,

    ℓ(μ) = ±1   (linear / antilinear),
    m(μ) = ±1   (automorphism / antiautomorphism),

and ONE transport law, the engine, governs how all four move a commutator
generator: for ANY H (Hermitian or not) and L_H = −i[H, ·],

    μ ∘ L_H ∘ μ  =  ℓ(μ)·m(μ) · L_{μ(H)},

so id → +L_H, θ → −L_{Hᵀ}, conj → −L_{H̄}, † → +L_{H†}. On a Pauli string σ
the maps collapse to signs, θ(σ) = conj(σ) = (−1)^{n_Y(σ)}·σ and †(σ) = σ;
in the Pauli coefficient basis θ acts as the D₄ mirror D = diag((−1)^{n_Y})
(ℂ-linear), † is pure coefficient conjugation 𝒦, and conj = D∘𝒦.

Five existing proofs in this repo are corollary legs of the engine. Each block
below re-derives its leg from the engine and cross-checks the published form:

  (a) F114, the θ-leg: D·L_σ·D = −L_{σᵀ} = (−1)^{n_Y(σ)+1}·L_σ for every Pauli
      string σ ≠ I. Home: compute/RCPsiSquared.Core/Symmetry/
      CommutatorDConjugationSign.cs and the F114 entry in
      docs/ANALYTICAL_FORMULAS.md.
  (b) The reversal kill, θ at word level: for a word of Pauli strings,
      Tr(A_j···A_1) = (−1)^{n_Y}·Tr(A_1···A_j) with n_Y the total Y count over
      the factors, via Tr(Wᵀ) = Tr(W) and the antiautomorphism reversing the
      product. Home: docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md §4.
  (c) F112 Lemmas A+B, the †-leg: in the row-stacking vec representation
      (vec(AρB) = (A⊗Bᵀ)vec(ρ)) the matrix of L_H is M = −i(H⊗I − I⊗Hᵀ), its
      Hilbert-Schmidt adjoint is the conjugate transpose, (L_H)* = −L_{H†},
      Hermitian H gives a skew matrix, and the dagger map is an antilinear
      isometry that conjugates Π-conjugation eigenvalues,
      (ΠAΠ⁻¹)† = ΠA†Π⁻¹ for unitary Π. Home:
      docs/proofs/PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md.
  (d) F113 Lemma C, the collapse leg: H = H† ⟺ θ(H) = conj(H), i.e. Hᵀ = H̄
      (both directions), and the moment face: for Hermitian H and real Z_l the
      ket-leg moment is the conjugate of the bra-leg moment at the same
      indices, Tr(Z_{l1}(Hᵀ)^{a1}Z_{l2}(Hᵀ)^{a2}) = conj(Tr(Z_{l1}H^{a1}Z_{l2}H^{a2})).
      Home: docs/proofs/PROOF_F113_COEFFICIENT_DERIVATION.md; the F117
      Hermitian conjugacy rides on the same collapse.
  (e) The K_b conj-leg: a generator whose Pauli-basis matrix is REAL has
      conj-transport as a symmetry, 𝒯L𝒯 = L for the antilinear coefficient
      conjugation 𝒯, hence eigenpairs come in (λ, v), (λ̄, 𝒯v) pairs; the
      antilinearity is load-bearing for the λ ↔ λ̄ pairing (a linear symmetry
      would fix λ, only an antilinear one conjugates it). Home:
      docs/proofs/PROOF_F86B_OBSTRUCTION.md and
      docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md §8; the dressed version
      Σ₁∘conj on a coherence block is the F86 K_b mode mirror
      (simulations/f86_kb_chiral_mirror.py).

Plus two structural blocks: the antilinear double ⟨R, D, 𝒦⟩ closes at order
16 = D₄ × Z₂ with exactly 8 antilinear members; 𝒦 commutes with D, while R
carries ±i phases (Y·X = −iZ, Z·X = iY per site) so 𝒦 twists R by the D
mirror, 𝒦R𝒦 = conj(R) = D·R·D, and the CENTRAL antilinear involution is
conj = D𝒦 (the double is ⟨R, D⟩ × ⟨D𝒦⟩); the triangle {id, D, 𝒦, D𝒦} is
the D-axis Klein subgroup, the Pauli-basis image of {id, θ, †, conj}. And the
dial trio, the bridge to the S₃⋉D₄ circle completion: θ- and conj-transport
invert the rotation dial, μ ∘ Ad_{R_z(φ)} ∘ μ = Ad_{R_z(−φ)}, while
†-transport fixes Ad_U for every unitary U.

Block ledger
------------
  Block 1  the grading + the engine : Klein table of {id, θ, conj, †} with both characters
                                      multiplicative and realized (ℓ on scalars EXACT, m on
                                      products); Pauli-string action θ(σ) = conj(σ) =
                                      (−1)^{n_Y}σ, †(σ) = σ EXACT at N = 1..3; Pauli-basis
                                      representations D, 𝒦, D𝒦; the transport law EXACT at
                                      superoperator-matrix level (N = 2, 3, non-Hermitian H)
                                      and < 1e-12 on random states
  Block 2  F114, the θ-leg          : D·L_σ·D = −L_{σᵀ} = (−1)^{n_Y+1}·L_σ EXACT for all
                                      4^N − 1 strings at N = 2, 3 (integer-valued Pauli-basis
                                      matrices, bit-for-bit)
  Block 3  the reversal kill        : Tr(reversed word) = (−1)^{n_Y}·Tr(word) EXACT (Gaussian-
                                      integer traces) over random words plus constructed
                                      nonzero-trace words of BOTH n_Y parities (the odd ones
                                      flip sign: the kill is non-vacuous)
  Block 4  F112 Lemmas A+B, †-leg   : M = −i(H⊗I − I⊗Hᵀ) acts as −i[H,·]; HS adjoint =
                                      conjugate transpose (pairing check); (L_H)* = −L_{H†}
                                      EXACT; Hermitian ⟹ skew EXACT; (ΠAΠ⁻¹)† = ΠA†Π⁻¹
                                      EXACT for the monomial Pauli-basis Π and < 1e-11 for a
                                      random unitary (with a non-unitary control that FAILS);
                                      eigen-projection corollary: C(A_s) = i·A_s ⟹
                                      C(A_s†) = −i·A_s† (< 1e-13, the projector
                                      reassociates the four-term sum)
  Block 5  F113 Lemma C + moments   : Hermitian ⟺ Hᵀ = H̄ both directions EXACT (with a
                                      generic non-Hermitian control); moment face at N = 3,
                                      five index choices, < 1e-10 relative
  Block 6  the K_b conj-leg         : framework Z-dephasing Lindbladian of an XXZ-type H with
                                      real Pauli coefficients, transformed to the Pauli basis,
                                      is real to machine precision; 𝒯L𝒯 = L; eigenpair
                                      pairing (λ, v) ↔ (λ̄, 𝒯v) with the conjugated residual
                                      BIT-EQUAL to the original (antilinearity in action)
  Block 7  the antilinear double    : BFS closure of ⟨R, D, 𝒦⟩ on the 4^N Pauli coefficient
                                      space (N = 2): order 16 = D₄ × Z₂, 8 antilinear
                                      members, 𝒦 commutes with D but 𝒦R𝒦 = conj(R) = DRD
                                      (the D-mirror twist), conj = D𝒦 central, D·R of
                                      order 4 (the D₄ witness), the triangle {id, D, 𝒦, D𝒦}
                                      a Klein subgroup, all EXACT
  Block 8  the dial trio            : θ and conj invert Ad_{R_z(φ)} (matrix level to one ulp
                                      via conj(R_z(φ)) = R_z(−φ); NumPy's FMA complex multiply
                                      is not bit-commutative, so not bit-for-bit); † fixes
                                      Ad_U for every unitary (R_z and a random unitary);
                                      apply level < 1e-13

Run: python simulations/antilinear_triangle.py (< 1 s).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import (  # noqa: E402
    _k_to_indices, _vec_to_pauli_basis_transform, pauli_basis_vector,
    pauli_string, site_op)
from framework.symmetry import build_pi_full  # noqa: E402

RNG = np.random.default_rng(2026)

# The four maps with their grading characters: name -> (callable, l, m).
# l = +1 linear / -1 antilinear; m = +1 automorphism / -1 antiautomorphism.
MAPS = {
    'id':     (lambda A: A,          +1, +1),
    'theta':  (lambda A: A.T,        +1, -1),
    'conj':   (lambda A: A.conj(),   -1, +1),
    'dagger': (lambda A: A.conj().T, -1, -1),
}


def rand_complex(shape):
    return RNG.standard_normal(shape) + 1j * RNG.standard_normal(shape)


def L_apply(H, rho):
    """L_H(rho) = -i[H, rho]."""
    return -1j * (H @ rho - rho @ H)


def L_vec_matrix(H):
    """Row-stacking vec matrix of L_H: vec(A rho B) = (A kron B^T) vec(rho)."""
    d = H.shape[0]
    Id = np.eye(d)
    return -1j * (np.kron(H, Id) - np.kron(Id, H.T))


def vec_transpose_perm(d):
    """Permutation p with vec(A^T) = vec(A)[p] in row-stacking (the swap S)."""
    return np.arange(d * d).reshape(d, d).T.reshape(-1)


def n_y(k, N):
    """Number of Y letters of Pauli string k; Y is the symplectic pair (1, 1)."""
    return sum(1 for ab in _k_to_indices(k, N) if ab == (1, 1))


def string_index_of(Qmat, N):
    """Index m and phase c with Qmat = c * sigma_m, for Qmat a phase times a string."""
    d = 2 ** N
    for m in range(4 ** N):
        Pm = pauli_string(list(_k_to_indices(m, N)))
        c = np.trace(Pm.conj().T @ Qmat) / d
        if abs(c) > 0.5:
            return m, c
    raise AssertionError("product is not proportional to a Pauli string")


# ======================================================================
# BLOCK 1 -- the grading + the engine.
# ======================================================================
def block1_grading_and_engine():
    print("-" * 92)
    print("BLOCK 1  the grading + the engine  [Klein table, Pauli action, D/K/DK "
          "representations, transport law]")
    print("-" * 92)

    # (i) Klein four-group with both characters multiplicative: mu o nu is the
    # unique member with characters (l_mu*l_nu, m_mu*m_nu), bit-exact on matrices.
    probes = [rand_complex((4, 4)) for _ in range(4)]
    for n1, (f1, l1, m1) in MAPS.items():
        for n2, (f2, l2, m2) in MAPS.items():
            matches = [n for n, t in MAPS.items() if t[1] == l1 * l2 and t[2] == m1 * m2]
            assert len(matches) == 1, "characters do not separate the four maps"
            f3 = MAPS[matches[0]][0]
            for A in probes:
                assert np.array_equal(f1(f2(A)), f3(A)), \
                    f"{n1} o {n2} != {matches[0]} (Klein table broken)"
    print("  Klein table: mu o nu = the member with characters (l*l', m*m'), bit-exact; "
          "dagger = theta o conj = conj o theta  OK")

    # (ii) the characters are realized, not just labels: l on scalars (EXACT),
    # m on products (float matmul, < 1e-12).
    c = 0.37 - 1.21j
    for name, (f, l_chr, m_chr) in MAPS.items():
        for A in probes:
            want = c * f(A) if l_chr == +1 else np.conj(c) * f(A)
            assert np.array_equal(f(c * A), want), f"{name}: l character not realized"
        A, B = probes[0], probes[1]
        want = f(A) @ f(B) if m_chr == +1 else f(B) @ f(A)
        dev = np.max(np.abs(f(A @ B) - want))
        assert dev < 1e-12, f"{name}: m character not realized (dev {dev:.2e})"
    print("  characters realized: l on scalars EXACT, m on operator products < 1e-12  OK")

    # (iii) Pauli-string action: theta(sigma) = conj(sigma) = (-1)^{n_Y} sigma,
    # dagger(sigma) = sigma, EXACT for all strings at N = 1, 2, 3.
    for N in (1, 2, 3):
        for k in range(4 ** N):
            sigma = pauli_string(list(_k_to_indices(k, N)))
            sgn = (-1.0) ** n_y(k, N)
            assert np.array_equal(sigma.T, sgn * sigma), f"theta action fails at N={N} k={k}"
            assert np.array_equal(sigma.conj(), sgn * sigma), \
                f"conj action fails at N={N} k={k}"
            assert np.array_equal(sigma.conj().T, sigma), f"dagger action fails at N={N} k={k}"
    print("  Pauli action: theta(sigma) = conj(sigma) = (-1)^{n_Y} sigma, "
          "dagger(sigma) = sigma, EXACT for all strings N = 1..3  OK")

    # (iv) Pauli-basis representations: theta = D (C-linear sign mirror), dagger = K
    # (pure coefficient conjugation), conj = D o K. Float traces, < 1e-13.
    N = 2
    Dvec = np.array([(-1.0) ** n_y(k, N) for k in range(4 ** N)])
    for _ in range(4):
        A = rand_complex((4, 4))
        coeff = pauli_basis_vector(A, N)
        dev_th = np.max(np.abs(pauli_basis_vector(A.T, N) - Dvec * coeff))
        dev_dg = np.max(np.abs(pauli_basis_vector(A.conj().T, N) - np.conj(coeff)))
        dev_cj = np.max(np.abs(pauli_basis_vector(A.conj(), N) - Dvec * np.conj(coeff)))
        assert dev_th < 1e-13 and dev_dg < 1e-13 and dev_cj < 1e-13, \
            f"Pauli-basis representation broken ({dev_th:.2e}, {dev_dg:.2e}, {dev_cj:.2e})"
    print("  Pauli-basis representations: theta = D = diag((-1)^{n_Y}), dagger = K, "
          "conj = D o K  (< 1e-13)  OK")

    # (v) THE ENGINE at superoperator-matrix level, EXACT, non-Hermitian H, N = 2, 3.
    # In row-stacking, theta is the swap S; the three nontrivial conjugations are
    # S M S, conj(M), and S conj(M) S, all exact entrywise rearrangements.
    for N in (2, 3):
        d = 2 ** N
        H = rand_complex((d, d))
        M = L_vec_matrix(H)
        p = vec_transpose_perm(d)
        assert np.array_equal(M[np.ix_(p, p)], -L_vec_matrix(H.T)), \
            f"engine theta-leg not exact at N={N}"
        assert np.array_equal(M.conj(), -L_vec_matrix(H.conj())), \
            f"engine conj-leg not exact at N={N}"
        assert np.array_equal(M.conj()[np.ix_(p, p)], L_vec_matrix(H.conj().T)), \
            f"engine dagger-leg not exact at N={N}"
        print(f"  ENGINE (matrix level, N={N}): theta -> -L_(H^T), conj -> -L_(H_bar), "
              f"dagger -> +L_(H^dag)  EXACT (bit-for-bit)")

    # (vi) the engine on random states (apply level), all four maps, < 1e-12.
    for N in (2, 3):
        d = 2 ** N
        H = rand_complex((d, d))
        for name, (mu, l_chr, m_chr) in MAPS.items():
            sign = l_chr * m_chr
            worst = 0.0
            for _ in range(6):
                rho = rand_complex((d, d))
                lhs = mu(L_apply(H, mu(rho)))
                rhs = sign * L_apply(mu(H), rho)
                worst = max(worst, np.max(np.abs(lhs - rhs)))
            assert worst < 1e-12, f"engine apply-level {name} N={N}: dev {worst:.2e}"
        print(f"  ENGINE (apply level, N={N}): mu o L_H o mu = l(mu)*m(mu) * L_(mu(H)) "
              f"for all four mu  (< 1e-12)  OK")
    print("BLOCK 1 PASS")


# ======================================================================
# BLOCK 2 -- F114, the theta-leg: D L_sigma D = -L_(sigma^T) = (-1)^{n_Y+1} L_sigma.
# The Pauli-basis convention is the framework's own (grepped from
# framework/symmetry.py y_parity_panel): L_pauli = (Mb^dag @ L @ Mb) / 2^N with
# Mb = _vec_to_pauli_basis_transform(N). All entries are Gaussian integers over
# 2^N, so every assertion is bit-for-bit.
# ======================================================================
def block2_f114_theta_leg():
    print("-" * 92)
    print("BLOCK 2  F114, the theta-leg  [D L_sigma D = -L_(sigma^T) = (-1)^{n_Y+1} L_sigma, "
          "EXACT, N = 2, 3]")
    print("-" * 92)
    for N in (2, 3):
        d2 = 4 ** N
        Mb = _vec_to_pauli_basis_transform(N)
        # orthogonality of the basis matrix, exact (integer matmul)
        assert np.array_equal(Mb.conj().T @ Mb, (2 ** N) * np.eye(d2)), \
            f"Mb^dag Mb != 2^N I at N={N}"
        Dvec = np.array([(-1.0) ** n_y(k, N) for k in range(d2)])
        n_minus = n_plus = 0
        for k in range(1, d2):
            sigma = pauli_string(list(_k_to_indices(k, N)))
            Lp = (Mb.conj().T @ L_vec_matrix(sigma) @ Mb) / 2 ** N
            LpT = (Mb.conj().T @ L_vec_matrix(sigma.T) @ Mb) / 2 ** N
            DLD = Dvec[:, None] * Lp * Dvec[None, :]
            eps = (-1.0) ** (n_y(k, N) + 1)
            # engine: D L D = theta-transport = l(theta)*m(theta) * L_(sigma^T) = -L_(sigma^T)
            assert np.array_equal(DLD, -LpT), f"D L D != -L_(sigma^T) at N={N} k={k}"
            # collapse: sigma^T = (-1)^{n_Y} sigma, so the sign is (-1)^{n_Y+1}
            assert np.array_equal(DLD, eps * Lp), f"F114 sign fails at N={N} k={k}"
            if eps < 0:
                n_minus += 1
            else:
                n_plus += 1
        print(f"  N={N}: all {d2 - 1} strings EXACT; eps = -1 for {n_minus} (n_Y even), "
              f"eps = +1 for {n_plus} (n_Y odd)")
    print("BLOCK 2 PASS")


# ======================================================================
# BLOCK 3 -- the reversal kill: Tr(A_j...A_1) = (-1)^{n_Y} Tr(A_1...A_j).
# Via the engine's theta-leg at word level: Tr(W) = Tr(W^T), and the
# antiautomorphism reverses the product while each factor pays (-1)^{n_Y(factor)}.
# Traces of Pauli words are Gaussian integers, so equality is EXACT.
# ======================================================================
def block3_reversal_kill():
    print("-" * 92)
    print("BLOCK 3  the reversal kill  [Tr(reversed) = (-1)^{n_Y} Tr(word), EXACT, "
          "both n_Y parities non-vacuous]")
    print("-" * 92)
    N, d = 2, 4
    words = []
    for _ in range(40):
        j = int(RNG.integers(2, 7))
        words.append([int(RNG.integers(0, 4 ** N)) for _ in range(j)])
    # constructed nonzero-trace words [P, Q, S] with P*Q proportional to S, so the
    # product is proportional to I; collect 6 of each total-n_Y parity
    odd_words, even_words = [], []
    attempts = 0
    while (len(odd_words) < 6 or len(even_words) < 6) and attempts < 1000:
        attempts += 1
        kp, kq = int(RNG.integers(0, 4 ** N)), int(RNG.integers(0, 4 ** N))
        P = pauli_string(list(_k_to_indices(kp, N)))
        Q = pauli_string(list(_k_to_indices(kq, N)))
        ks, _ = string_index_of(P @ Q, N)
        wd = [kp, kq, ks]
        tot = sum(n_y(k, N) for k in wd)
        if tot % 2 == 1 and len(odd_words) < 6:
            odd_words.append(wd)
        elif tot % 2 == 0 and len(even_words) < 6:
            even_words.append(wd)
    assert len(odd_words) == 6 and len(even_words) == 6, "constructed-word search failed"
    words += odd_words + even_words

    nz_odd = nz_even = 0
    for wd in words:
        mats = [pauli_string(list(_k_to_indices(k, N))) for k in wd]
        fwd = np.eye(d, dtype=complex)
        for Wm in mats:
            fwd = fwd @ Wm
        rev = np.eye(d, dtype=complex)
        for Wm in reversed(mats):
            rev = rev @ Wm
        nY = sum(n_y(k, N) for k in wd)
        sgn = (-1.0) ** nY
        tr_f, tr_r = np.trace(fwd), np.trace(rev)
        assert abs(tr_r - sgn * tr_f) == 0.0, \
            f"reversal kill broken on word {wd} (n_Y = {nY})"
        if abs(tr_f) > 0:
            if nY % 2 == 1:
                nz_odd += 1
            else:
                nz_even += 1
    assert nz_odd >= 6 and nz_even >= 6, "non-vacuity counts missed"
    print(f"  {len(words)} words checked EXACT (Gaussian-integer traces); nonzero-trace "
          f"witnesses: {nz_even} even-n_Y (sign +1), {nz_odd} odd-n_Y (sign -1, the kill)")
    print("BLOCK 3 PASS")


# ======================================================================
# BLOCK 4 -- F112 Lemmas A+B, the dagger-leg (the Hilbert-Schmidt adjoint level).
# ======================================================================
def block4_f112_dagger_leg():
    print("-" * 92)
    print("BLOCK 4  F112 Lemmas A+B, the dagger-leg  [(L_H)* = -L_(H^dag); Hermitian => skew; "
          "Pi-conjugation eigenvalues conjugate]")
    print("-" * 92)
    N, d = 2, 4
    H_gen = rand_complex((d, d))
    G = rand_complex((d, d))
    H_herm = G + G.conj().T

    # Lemma A, step 1: M = -i(H kron I - I kron H^T) IS the vec matrix of -i[H, .]
    M = L_vec_matrix(H_gen)
    worst = 0.0
    for _ in range(5):
        rho = rand_complex((d, d))
        lhs = (M @ rho.reshape(-1)).reshape(d, d)
        worst = max(worst, np.max(np.abs(lhs - L_apply(H_gen, rho))))
    assert worst < 1e-12, f"vec representation broken (dev {worst:.2e})"
    print(f"  vec representation: M = -i(H kron I - I kron H^T) acts as -i[H, .] "
          f"(dev {worst:.1e})  OK")

    # Lemma A, step 2: the HS adjoint is the conjugate transpose:
    # <A, L(B)>_HS = <M^dag(A), B>_HS for random A, B
    worst = 0.0
    for _ in range(5):
        A, B = rand_complex((d, d)), rand_complex((d, d))
        lhs = np.trace(A.conj().T @ L_apply(H_gen, B))
        Lstar_A = (M.conj().T @ A.reshape(-1)).reshape(d, d)
        rhs = np.trace(Lstar_A.conj().T @ B)
        worst = max(worst, abs(lhs - rhs) / max(1.0, abs(lhs)))
    assert worst < 1e-12, f"HS pairing broken (rel dev {worst:.2e})"
    print(f"  HS adjoint = conjugate transpose: <A, L(B)> = <L*(A), B> "
          f"(rel dev {worst:.1e})  OK")

    # Lemma A, step 3: (L_H)* = -L_(H^dag), EXACT (the engine's dagger-leg read at
    # the adjoint level: dagger has sign +1, the adjoint contributes the minus)
    assert np.array_equal(M.conj().T, -L_vec_matrix(H_gen.conj().T)), \
        "(L_H)* != -L_(H^dag)"
    print("  (L_H)* = -L_(H^dag) for generic non-Hermitian H  EXACT")

    # Lemma B, step 1: Hermitian H => skew matrix, EXACT
    Mh = L_vec_matrix(H_herm)
    assert np.array_equal(Mh.conj().T, -Mh), "Hermitian H: M not skew"
    print("  Hermitian H: M* = -M (skew)  EXACT")

    # Lemma B, step 2: the dagger map conjugates Pi-conjugation eigenvalues:
    # (Pi A Pi^{-1})^dag = Pi A^dag Pi^{-1} for unitary Pi. First with the
    # Pauli-basis Pi mirror (a monomial unitary, so everything is EXACT).
    d2 = 4 ** N
    Pi = build_pi_full(N, dephase_letter='Z')
    Pi_inv = Pi.conj().T
    assert np.array_equal(Pi @ Pi_inv, np.eye(d2)), "Pi not unitary (exact)"
    for _ in range(4):
        A = rand_complex((d2, d2))
        lhs = (Pi @ A @ Pi_inv).conj().T
        rhs = Pi @ A.conj().T @ Pi_inv
        assert np.max(np.abs(lhs - rhs)) == 0.0, "(Pi A Pi^{-1})^dag != Pi A^dag Pi^{-1}"
    print("  (Pi A Pi^{-1})^dag = Pi A^dag Pi^{-1} with Pi = the Pauli-basis Pi mirror  EXACT")

    # eigen-projection corollary: C = Pi . Pi^{-1} has order 4 (Pi^4 = I); project a
    # random A onto the s = i eigenspace; then C(A_s) = i A_s and the dagger map sends
    # it to the s = -i eigenspace: C(A_s^dag) = -i A_s^dag. The monomial conjugation is
    # exact, but the projector reassociates the four-term sum on the two sides, so the
    # comparison is to rounding (< 1e-13), not bit-for-bit.
    A = rand_complex((d2, d2))
    conj_pi = lambda X: Pi @ X @ Pi_inv  # noqa: E731
    Ck = [A]
    for _ in range(4):
        Ck.append(conj_pi(Ck[-1]))
    assert np.array_equal(Ck[4], A), "C^4 != id"
    As = (Ck[0] - 1j * Ck[1] - Ck[2] + 1j * Ck[3]) / 4
    assert np.linalg.norm(As) > 1e-3, "i-eigenspace projection vanished (vacuous)"
    dev_s = np.max(np.abs(conj_pi(As) - 1j * As))
    assert dev_s < 1e-13, f"C(A_s) != i A_s (dev {dev_s:.2e})"
    Asd = As.conj().T
    dev_sd = np.max(np.abs(conj_pi(Asd) + 1j * Asd))
    assert dev_sd < 1e-13, f"C(A_s^dag) != -i A_s^dag (dev {dev_sd:.2e})"
    print(f"  eigen corollary: C(A_s) = i*A_s  ==>  C(A_s^dag) = -i*A_s^dag "
          f"(dagger conjugates the eigenvalue; devs {dev_s:.1e}, {dev_sd:.1e})")

    # ... and with a random unitary conjugation superoperator (Pi^{-1} computed
    # numerically, so float tolerance), plus a non-unitary control that FAILS.
    Q, _ = np.linalg.qr(rand_complex((d2, d2)))
    Q_inv = np.linalg.inv(Q)
    A = rand_complex((d2, d2))
    dev_u = np.max(np.abs((Q @ A @ Q_inv).conj().T - Q @ A.conj().T @ Q_inv))
    assert dev_u < 1e-11, f"random-unitary conjugation dev {dev_u:.2e}"
    T = rand_complex((d2, d2))
    T_inv = np.linalg.inv(T)
    dev_t = np.max(np.abs((T @ A @ T_inv).conj().T - T @ A.conj().T @ T_inv))
    assert dev_t > 1e-2, f"non-unitary control unexpectedly passed (dev {dev_t:.2e})"
    print(f"  random unitary Pi: dev = {dev_u:.1e}  OK; non-unitary control FAILS as it must "
          f"(dev = {dev_t:.1f}): unitarity is load-bearing")
    print("BLOCK 4 PASS")


# ======================================================================
# BLOCK 5 -- F113 Lemma C: the Hermitian collapse theta(H) = conj(H) <=> H = H^dag,
# and the moment face (ket-leg = conjugate of bra-leg at the same indices).
# ======================================================================
def block5_hermitian_collapse():
    print("-" * 92)
    print("BLOCK 5  F113 Lemma C  [H = H^dag <=> H^T = H_bar, both directions; "
          "the moment face]")
    print("-" * 92)
    N, d = 3, 8
    G = rand_complex((d, d))
    H = G + G.conj().T

    # direction 1: Hermitian => theta(H) = conj(H), EXACT
    assert np.array_equal(H.T, H.conj()), "Hermitian H: H^T != H_bar"
    print("  Hermitian H:  H^T = H_bar  EXACT")

    # direction 2: H^T = H_bar by construction => Hermitian, EXACT
    # (H2 = R + iS with R symmetric, S antisymmetric, both real)
    Ar = RNG.standard_normal((d, d))
    Br = RNG.standard_normal((d, d))
    H2 = (Ar + Ar.T) + 1j * (Br - Br.T)
    assert np.array_equal(H2.T, H2.conj()), "premise H2^T = H2_bar broken"
    assert np.array_equal(H2, H2.conj().T), "H2^T = H2_bar did not force Hermitian"
    print("  H with H^T = H_bar by construction:  H = H^dag  EXACT")

    # control: generic non-Hermitian H differs in both readings
    H3 = rand_complex((d, d))
    gap_collapse = np.max(np.abs(H3.T - H3.conj()))
    gap_herm = np.max(np.abs(H3 - H3.conj().T))
    assert gap_collapse > 1e-1 and gap_herm > 1e-1, "non-Hermitian control degenerate"
    print(f"  generic non-Hermitian control: max|H^T - H_bar| = {gap_collapse:.2f} > 0, "
          f"max|H - H^dag| = {gap_herm:.2f} > 0  OK")

    # the moment face: for Hermitian H and real Z_l,
    # Tr(Z_l1 (H^T)^a1 Z_l2 (H^T)^a2) = conj(Tr(Z_l1 H^a1 Z_l2 H^a2))
    assert all(np.max(np.abs(site_op(N, l, 'Z').imag)) == 0.0 for l in range(N)), \
        "Z site operators not real"
    Zs = [site_op(N, l, 'Z').real for l in range(N)]
    HT = H.T
    mp = np.linalg.matrix_power
    worst = 0.0
    for (l1, a1, l2, a2) in [(0, 1, 1, 2), (0, 2, 2, 3), (1, 3, 2, 1),
                             (2, 4, 0, 2), (1, 1, 1, 1)]:
        ket = np.trace(Zs[l1] @ mp(HT, a1) @ Zs[l2] @ mp(HT, a2))
        bra = np.trace(Zs[l1] @ mp(H, a1) @ Zs[l2] @ mp(H, a2))
        rel = abs(ket - np.conj(bra)) / max(1.0, abs(bra))
        worst = max(worst, rel)
        assert rel < 1e-10, f"moment face fails at (l1,a1,l2,a2)=({l1},{a1},{l2},{a2})"
    print(f"  moment face: ket-leg = conj(bra-leg) at 5 index choices, random Hermitian H, "
          f"N=3 (worst rel dev {worst:.1e})  OK")
    print("BLOCK 5 PASS")


# ======================================================================
# BLOCK 6 -- the K_b conj-leg: a generator with a REAL Pauli-basis matrix has
# conj-transport as a symmetry; eigenpairs come in (lambda, v), (lambda_bar, T v).
# The dressed version Sigma_1 o conj on a coherence block is the F86 K_b mode
# mirror (simulations/f86_kb_chiral_mirror.py); here we exhibit the bare law.
# NOTE on the basis convention: Mb columns are COLUMN-stacked strings while the
# framework Lindbladian matrix is row-stacking, so (Mb^dag L Mb)/2^N is the
# theta-conjugated Pauli matrix D L_pauli D. D is a real sign matrix, so realness
# and the lambda <-> lambda_bar pairing are unaffected; this is the framework's
# own established reading (framework/symmetry.py y_parity_panel).
# ======================================================================
def block6_kb_conj_leg():
    print("-" * 92)
    print("BLOCK 6  the K_b conj-leg  [real Pauli-basis Lindbladian, T L T = L, "
          "eigenpair pairing (lambda, v) <-> (lambda_bar, T v)]")
    print("-" * 92)
    N = 3
    # XXZ-type chain with REAL Pauli coefficients (+ a real transverse field term)
    H = np.zeros((2 ** N, 2 ** N), dtype=complex)
    for l in range(N - 1):
        H += site_op(N, l, 'X') @ site_op(N, l + 1, 'X')
        H += site_op(N, l, 'Y') @ site_op(N, l + 1, 'Y')
        H += 0.7 * (site_op(N, l, 'Z') @ site_op(N, l + 1, 'Z'))
    H += 0.3 * site_op(N, 1, 'X')
    gammas = [0.35, 0.6, 0.2]
    L = lindbladian_pauli_dephasing(H, gammas, dephase_letter='Z')
    Mb = _vec_to_pauli_basis_transform(N)
    Lp = (Mb.conj().T @ L @ Mb) / 2 ** N

    im_max = float(np.max(np.abs(Lp.imag)))
    assert im_max < 1e-12, f"Pauli-basis Lindbladian not real (max |Im| = {im_max:.2e})"
    print(f"  Pauli-basis matrix real: max |Im| = {im_max:.1e} "
          f"(real H-coefficients => real structure constants)  OK")
    Lr = Lp.real  # the certified-real matrix; T L T = conj(L) = L is now manifest

    evals, evecs = np.linalg.eig(Lr)
    n_cplx = int(np.sum(np.abs(evals.imag) > 1e-8))
    assert n_cplx > 0, "no complex eigenvalues: pairing test vacuous"
    # the spectrum is closed under conjugation
    worst_pair = 0.0
    for lam in evals:
        worst_pair = max(worst_pair, float(np.min(np.abs(np.conj(lam) - evals))))
    assert worst_pair < 1e-9, f"spectrum not conj-closed (worst {worst_pair:.2e})"
    print(f"  spectrum conj-closed: {n_cplx}/{len(evals)} complex eigenvalues, worst "
          f"|conj(lambda) - spectrum| = {worst_pair:.1e}  OK")

    # exhibit the antilinear pairing: T = coefficient conjugation; for real L,
    # L(T v) = T(L v) BIT-EXACT, so (lambda, v) maps to (lambda_bar, T v) with the
    # conjugated residual BIT-EQUAL to the original. Antilinearity is load-bearing:
    # a linear symmetry would fix lambda, only an antilinear one conjugates it.
    j = int(np.argmax(np.abs(evals.imag)))
    lam, v = evals[j], evecs[:, j]
    assert np.array_equal(Lr @ v.conj(), (Lr @ v).conj()), "T L T = L not bit-exact"
    res_orig = float(np.linalg.norm(Lr @ v - lam * v))
    res_pair = float(np.linalg.norm(Lr @ v.conj() - np.conj(lam) * v.conj()))
    assert res_pair == res_orig, "paired residual not bit-equal to the original"
    assert res_orig < 1e-10, f"eig residual too large ({res_orig:.2e})"
    print(f"  pairing exhibited at lambda = {lam:.4f}: T v is an eigenvector of "
          f"lambda_bar, residual {res_pair:.1e} BIT-EQUAL to the original")
    print("BLOCK 6 PASS")


# ======================================================================
# BLOCK 7 -- the antilinear double <R, D, K> on the 4^N Pauli coefficient space.
# Elements are (matrix, antilinear flag) with composition
# (M1, a1) o (M2, a2) = (M1 . conj(M2) if a1 else M1 . M2, a1 xor a2);
# K is (I, True), D = diag((-1)^{n_Y}), R is rho -> rho . X^{kron N} (a signed
# permutation on Pauli coefficients). BFS closure, N = 2.
# ======================================================================
def block7_antilinear_double():
    print("-" * 92)
    print("BLOCK 7  the antilinear double <R, D, K>  [order 16 = D4 x Z2, 8 antilinear, "
          "conj = DK central, the D-axis Klein triangle]")
    print("-" * 92)
    N, d = 2, 4
    Mdim = 4 ** N
    Dmat = np.diag(np.array([(-1.0) ** n_y(k, N) for k in range(Mdim)]).astype(complex))
    # R: rho -> rho . X^{kron N} as a signed permutation on Pauli coefficients
    Xn = pauli_string(['X'] * N)
    Rmat = np.zeros((Mdim, Mdim), dtype=complex)
    for k in range(Mdim):
        P = pauli_string(list(_k_to_indices(k, N)))
        m, c = string_index_of(P @ Xn, N)
        Rmat[m, k] = c
    Imat = np.eye(Mdim, dtype=complex)
    assert np.array_equal(Dmat @ Dmat, Imat), "D^2 != I"
    assert np.array_equal(Rmat @ Rmat, Imat), "R^2 != I (X^N squares to I)"

    def compose(e1, e2):
        """Apply e2 then e1; the antilinear flag conjugates what follows: K M K = conj(M)."""
        M1, a1 = e1
        M2, a2 = e2
        return ((M1 @ (M2.conj() if a1 else M2)), a1 ^ a2)

    def key(e):
        Mm, a = e
        return ((np.round(Mm, 9) + 0.0).tobytes(), a)  # +0.0 normalizes -0.0

    gens = [(Dmat, False), (Rmat, False), (Imat, True)]
    ident = (Imat, False)
    elems = [ident]
    seen = {key(ident)}
    frontier = [ident]
    while frontier:
        nxt = []
        for e in frontier:
            for g in gens:
                ne = compose(g, e)
                k_ = key(ne)
                if k_ not in seen:
                    seen.add(k_)
                    elems.append(ne)
                    nxt.append(ne)
        frontier = nxt
    n_anti = sum(1 for _, a in elems if a)
    assert len(elems) == 16, f"|<R, D, K>| = {len(elems)} (expected 16)"
    assert n_anti == 8, f"antilinear members: {n_anti} (expected 8)"
    print(f"  closure: |<R, D, K>| = 16 with {n_anti} antilinear members  OK")

    # direct product structure: the 8 linear matrices reappear verbatim with the
    # antilinear flag, i.e. the 16 elements are (linear part) x {id, K} as a set
    lin_keys = {key((Mm, False)) for Mm, a in elems if not a}
    anti_as_lin = {key((Mm, False)) for Mm, a in elems if a}
    assert lin_keys == anti_as_lin, "the 16 elements are not (8 linear) x {id, K}"
    # K commutes with D (D is real); but K is NOT central: R carries +-i phases
    # (Y.X = -iZ, Z.X = iY per site), so K twists R by the D mirror,
    # K R K = conj(R) = D R D. The CENTRAL antilinear involution is conj = DK:
    # it commutes with both generators and squares to id, so the double is the
    # direct product <R, D> x <DK> = D4 x Z2.
    K_el, D_el, R_el, DK_el = (Imat, True), (Dmat, False), (Rmat, False), (Dmat, True)
    assert key(compose(K_el, D_el)) == key(compose(D_el, K_el)), "K does not commute with D"
    assert not np.array_equal(Rmat.conj(), Rmat), "R unexpectedly real (twist vacuous)"
    assert np.array_equal(Rmat.conj(), Dmat @ Rmat @ Dmat), "K R K = conj(R) != D R D"
    for g, gname in ((D_el, 'D'), (R_el, 'R')):
        assert key(compose(DK_el, g)) == key(compose(g, DK_el)), \
            f"DK does not commute with {gname}"
    assert key(compose(DK_el, DK_el)) == key(ident), "(DK)^2 != id"
    print("  16 = 8 linear x {id, K} as a set; K commutes with D, K R K = conj(R) = D R D "
          "(the D-mirror twist)")
    print("  the CENTRAL antilinear involution is conj = DK: the double is "
          "<R, D> x <DK> = D4 x Z2  OK")

    # the D4 witness: D and R are involutions and D R has order 4
    DR = Dmat @ Rmat
    DR2 = DR @ DR
    assert not np.array_equal(DR2, Imat), "(DR)^2 = I: linear part not D4"
    assert np.array_equal(DR2 @ DR2, Imat), "(DR)^4 != I"
    print("  linear part: D^2 = R^2 = I, (D R)^4 = I != (D R)^2, so <R, D> = D4 "
          "(8 elements)  OK")

    # the D-axis Klein triangle {id, D, K, DK}: closed, every non-identity member an
    # involution; it is the Pauli-basis image of {id, theta, dagger, conj} from Block 1
    triangle = [ident, (Dmat, False), (Imat, True), (Dmat, True)]
    tri_keys = {key(e) for e in triangle}
    for e1 in triangle:
        for e2 in triangle:
            assert key(compose(e1, e2)) in tri_keys, "triangle not closed"
    for e in triangle[1:]:
        assert key(compose(e, e)) == key(ident), "triangle member not an involution"
    print("  the triangle {id, D, K, DK} is a Klein four-group inside the double: the "
          "Pauli-basis image of {id, theta, dagger, conj}  OK")
    print("BLOCK 7 PASS")


# ======================================================================
# BLOCK 8 -- the dial trio (single qubit): theta- and conj-transport invert the
# rotation dial, mu o Ad_{R_z(phi)} o mu = Ad_{R_z(-phi)}; dagger-transport fixes
# Ad_U for EVERY unitary. The bridge to the S3 x| D4 circle completion.
# ======================================================================
def block8_dial_trio():
    print("-" * 92)
    print("BLOCK 8  the dial trio  [theta, conj invert Ad_{R_z(phi)}; dagger fixes Ad_U]")
    print("-" * 92)
    phi = 0.7321
    Rz = np.diag([np.exp(-1j * phi / 2), np.exp(1j * phi / 2)])
    Rz_m = Rz.conj()  # R_z(-phi) = conj(R_z(phi)) for the diagonal matrix
    Rz_m_fresh = np.diag([np.exp(1j * phi / 2), np.exp(-1j * phi / 2)])
    assert np.max(np.abs(Rz_m - Rz_m_fresh)) < 1e-15, "conj(R_z(phi)) != R_z(-phi)"

    # matrix level: Ad_U = U kron conj(U) in row-stacking; S the transpose swap.
    # theta-transport: S K S = conj(U) kron U = Ad_{R_z(-phi)}
    # conj-transport:  conj(K) = conj(U) kron U = Ad_{R_z(-phi)}
    # dagger-transport: S conj(K) S = K, i.e. Ad_U is FIXED
    # One ulp tolerance, not bit-for-bit: NumPy's FMA-based complex multiply is not
    # bitwise commutative (z * conj(z) carries an O(1e-17) imaginary rounding residue),
    # so the kron products on the two sides differ at the last bit.
    p = vec_transpose_perm(2)
    K = np.kron(Rz, Rz.conj())
    K_minus = np.kron(Rz_m, Rz_m.conj())
    dev_th = np.max(np.abs(K[np.ix_(p, p)] - K_minus))
    dev_cj = np.max(np.abs(K.conj() - K_minus))
    dev_dg = np.max(np.abs(K.conj()[np.ix_(p, p)] - K))
    assert dev_th < 1e-15, f"theta-transport != Ad_(R_z(-phi)) (dev {dev_th:.2e})"
    assert dev_cj < 1e-15, f"conj-transport != Ad_(R_z(-phi)) (dev {dev_cj:.2e})"
    assert dev_dg < 1e-15, f"dagger-transport moved the dial (dev {dev_dg:.2e})"
    U, _ = np.linalg.qr(rand_complex((2, 2)))
    KU = np.kron(U, U.conj())
    dev_du = np.max(np.abs(KU.conj()[np.ix_(p, p)] - KU))
    assert dev_du < 1e-15, f"dagger-transport moved a generic unitary dial ({dev_du:.2e})"
    print(f"  matrix level: theta ({dev_th:.1e}) and conj ({dev_cj:.1e}) send "
          f"Ad_(R_z(phi)) -> Ad_(R_z(-phi)); dagger fixes Ad_(R_z) ({dev_dg:.1e}) AND a "
          f"random Ad_U ({dev_du:.1e})  (one ulp)")

    # apply level on random states
    ad = lambda V, r: V @ r @ V.conj().T  # noqa: E731
    worst = {'theta': 0.0, 'conj': 0.0, 'dagger_rz': 0.0, 'dagger_u': 0.0}
    for _ in range(6):
        rho = rand_complex((2, 2))
        worst['theta'] = max(worst['theta'],
                             np.max(np.abs(ad(Rz, rho.T).T - ad(Rz_m, rho))))
        worst['conj'] = max(worst['conj'],
                            np.max(np.abs(ad(Rz, rho.conj()).conj() - ad(Rz_m, rho))))
        worst['dagger_rz'] = max(worst['dagger_rz'],
                                 np.max(np.abs(ad(Rz, rho.conj().T).conj().T - ad(Rz, rho))))
        worst['dagger_u'] = max(worst['dagger_u'],
                                np.max(np.abs(ad(U, rho.conj().T).conj().T - ad(U, rho))))
    for k, v in worst.items():
        assert v < 1e-13, f"dial trio apply level {k}: dev {v:.2e}"
    print(f"  apply level (phi = {phi}): theta {worst['theta']:.1e}, conj "
          f"{worst['conj']:.1e} invert; dagger fixes R_z {worst['dagger_rz']:.1e} and "
          f"generic U {worst['dagger_u']:.1e}  (< 1e-13)  OK")
    print("BLOCK 8 PASS")


def main():
    t_start = time.time()
    print("=" * 92)
    print("THE ANTILINEAR TRIANGLE -- self-validating verification "
          "(theta, conj, dagger and the one transport law)")
    print("=" * 92)

    # convention cross-check 1: the symplectic encoding (Y is the pair (1, 1))
    assert _k_to_indices(3, 1) == ((1, 1),), "k = 3 at N = 1 is not Y"
    assert np.array_equal(pauli_string([(1, 1)]),
                          np.array([[0, -1j], [1j, 0]])), "(1, 1) is not the Y matrix"
    assert np.array_equal(pauli_string([(1, 0)]),
                          np.array([[0, 1], [1, 0]], dtype=complex)), "(1, 0) is not X"
    assert np.array_equal(pauli_string([(0, 1)]),
                          np.array([[1, 0], [0, -1]], dtype=complex)), "(0, 1) is not Z"
    print("cross-check: symplectic pairs X = (1,0), Y = (1,1), Z = (0,1)  OK")

    # convention cross-check 2: the Pauli basis matrix columns are the column-stacked
    # strings, Mb[:, k] = vec_F(sigma_k), with Mb^dag Mb = 2^N I (the framework's
    # L_pauli = (Mb^dag L Mb)/2^N convention, framework/symmetry.py y_parity_panel)
    Mb = _vec_to_pauli_basis_transform(2)
    for k in (1, 7, 11):
        sigma = pauli_string(list(_k_to_indices(k, 2)))
        assert np.array_equal(Mb[:, k], sigma.flatten('F')), f"Mb column {k} convention"
    assert np.array_equal(Mb.conj().T @ Mb, 4 * np.eye(16)), "Mb^dag Mb != 2^N I"
    print("cross-check: Mb columns = column-stacked Pauli strings, Mb^dag Mb = 2^N I  OK\n")

    block1_grading_and_engine()
    block2_f114_theta_leg()
    block3_reversal_kill()
    block4_f112_dagger_leg()
    block5_hermitian_collapse()
    block6_kb_conj_leg()
    block7_antilinear_double()
    block8_dial_trio()

    print("=" * 92)
    print(f"ALL BLOCKS PASS   ({time.time() - t_start:.1f}s)")
    print("=" * 92)


if __name__ == "__main__":
    main()
