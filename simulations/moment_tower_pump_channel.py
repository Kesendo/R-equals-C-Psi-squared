#!/usr/bin/env python3
r"""The moment-tower pump channel: the deg-1 girth ladder read LINEARLY (committed verifier).

THE PUMP-SLOPE LAW. For L = −i[H,·] + Σ_l γ^deph_l D[Z_l] + Σ_l γ↓_l D[σ⁻_l] + Σ_l γ↑_l D[σ⁺_l]
(framework conventions: σ⁻ = (X+iY)/2 = [[0,1],[0,0]], D[c](ρ) = cρc† − ½{c†c, ρ}), the slope of
any measured polynomial A at the maximally mixed state reads the Z_l-weighted moment tower
linearly:

    d/dt Tr(A ρ)|_{ρ = I/d} = (1/d) · Σ_l Δγ_l · Tr(A Z_l),    Δγ_l = γ↓_l − γ↑_l,

because D[σ⁻_l](I) = +γ Z_l and D[σ⁺_l](I) = −γ Z_l per site while the dephasing (unital) and
the unitary part kill I. With A = H^j the slope is (1/d) Σ_l Δγ_l t_j(l) with t_j(l) = Tr(Z_l H^j):
the deg-1 girth-ladder tower of f87_girth_dichotomy (its t_moment), read LINEARLY by a physical
pump. Three exact blindnesses follow: dephasing-blind (unital), evolution-blind (the slope does
not depend on the H generating the dynamics, only on the measured polynomial A), and the
detailed-balance closure (Δγ ≡ 0 ⟹ slope ≡ 0; the F84 vacuum-only statement,
docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md).

THE F113 BRIDGE. Within F113's scope (H = Σ_l (ω_l/2) Z_l, cooling γ_T1, heating γ_pump) the
STATIC Frobenius polarity asymmetry equals the DYNAMIC pump slope:

    asymmetry = ‖M₊‖² − ‖M₋‖² = (4^N/2) · Σ_l ω_l (γ_pump,l − γ_T1,l) = −4^N · slope⟨H⟩.

THE CURVATURE FINGERPRINT. The second derivative d²/dt² Tr(A ρ)|_{I/d} =
(1/d) Σ_l Δγ_l Tr(A L(Z_l)) is EXACTLY affine in the evolution generator (the pump is
generator-independent; L acts once). For a chip Hamiltonian H_c = H_p + δV the curvature
difference is exactly (δ/d) Σ_l Δγ_l · (−i) Tr(V [Z_l, A]): linear in the parasite V against the
commutator probes [Z_l, H_p^j]. Z-flavored parasites are EXACTLY invisible ([Z, Z_l] = 0); that
is F113's complementary territory (its balance channel reads exactly the Z-drives); X/Y-flavored
parasites with overlap on the probe set are read linearly.

THE GIRTH CERTIFICATE (one-sided, deg-1 face only). Scan j = 1..j_max; the first j where the
slope fires is the girth ℓ (per the dichotomy: t_ℓ ≠ 0 somewhere ⟹ m* = 2ℓ+1, deg = 1, hard at
every γ > 0). HONEST LIMIT: silence of the deg-1 tower is NOT softness; the k = 4 witness
IIXY+ZXZY has t_j = 0 everywhere yet is hard at m* = 11 = 2ℓ+5 with p₁₁ = 86507520·γ⁵ (deg-5,
pinned exactly in f87_girth_dichotomy). Site-sum caveat: Σ_l Δγ_l t_ℓ(l) can cancel accidentally;
site-resolved weights (per-qubit calibration Δγ_l, or selective per-site damping) resolve t_ℓ(l)
individually (demonstrated in Block D).

Block ledger
------------
  Block A  pump directions            : D[σ⁻_l](I) = +γZ_l, D[σ⁺_l](I) = −γZ_l per site BIT-EXACT
                                        (exact-square rates), unitary + dephasing kill I bit-exact,
                                        mixed channel nets Δγ_l; builder cross-checked bit-for-bit
                                        against framework lindbladian_z_plus_t1
  Block B  the slope law vs dense L   : random Hermitian H and a structured chain H (N = 3),
                                        site-dependent rates, all three channels together,
                                        A = H^j (j = 1..4) and random Hermitian A;
                                        dev ≤ 1e-14 · max(1, |slope|); t_j(l) cross-checked
                                        against f87_girth_dichotomy.t_moment
  Block C  the three blindnesses      : dephasing-rate change and evolution-H change leave the
                                        slope BIT-IDENTICAL (the diagonal columns of L carry no
                                        unitary or dephasing residue); detailed balance γ↑ = γ↓
                                        gives slope = 0 exactly
  Block D  the girth certificate      : girth-1 (Z-drive) fires at j = 1 + the site-sum caveat
                                        (balanced ω cancels the uniform-weight sum; per-site
                                        weights resolve it); girth-2 witness X₀ + X₀Z₁ + 0.7·X₁X₂
                                        (t₁ ≡ 0, t₂ = [0, 16, 0]) silent at j = 1, fires at j = 2;
                                        HONEST negative IIXY+ZXZY (N = 5): t_j = 0 EXACTLY through
                                        j = 5, all slopes silent, yet hard at m* = 11 (deg-5)
  Block E  the F113 bridge            : N = 2 and 3, F113-scope generator; closed form ==
                                        −4^N·slope(dense L) == polarity_coordinates_from_hc
                                        asymmetry (three-way, ≤ 1e-12; closed vs slope ≤ 1e-15);
                                        bridge unchanged under added Z-dephasing
  Block F  the curvature fingerprint  : curvature difference exactly affine in δ (two δ values,
                                        bit-exact doubling); probe law dev ≤ 1e-15 for Y₀, Y₀Z₁,
                                        Y₁X₂ parasites (all read, nonzero); Z₀ parasite invisible
                                        BIT-EXACT; the F113-complementarity line
  Block G  the finite-time protocol   : ρ(t) = expm(Lt)(I/d) at small t; finite-difference slope
                                        recovers the law within the window error; the curvature is
                                        the first correction ((t/2)·curv ratio → 1, Richardson
                                        extrapolation cancels it)
  Block H  detailed balance + F84 tie : γ↓ = γ↑ per site ⟹ L(I) = 0 BIT-EXACT ⟹ every slope
                                        zero exactly (battery over H^j and random A); the pump
                                        weight Δγ_l IS F84's net vacuum rate (cross-tied against
                                        predict_amplitude_damping_violation)

Provenance: 2026-06-11, hunt #2 of the connection-hunt (originating scout _pump_slope_scout.py, a WIP not retained in the repo).
Framework diagnostic: simulations/framework/diagnostics/f120_moment_tower.py (F120).
Run: python simulations/moment_tower_pump_channel.py (< 1 min).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_general, lindbladian_z_plus_t1  # noqa: E402
from framework.pauli import _build_bilinear, _build_kbody_chain, site_op  # noqa: E402
from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_hc  # noqa: E402
from f87_girth_dichotomy import t_moment  # noqa: E402


# ======================================================================
# Builders and readers. vec convention: row-stacking flatten() (C order), matching
# lindbladian_general's −i(kron(H, I) − kron(I, Hᵀ)).
# ======================================================================
SM = np.array([[0, 1], [0, 0]], dtype=complex)   # σ⁻ = (X+iY)/2 (framework convention)
SP = SM.conj().T                                  # σ⁺


def site_2x2(N, l, op2):
    """N-qubit operator with the given 2×2 block on site l, identity elsewhere."""
    ops = [np.eye(2, dtype=complex)] * N
    ops[l] = op2
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def build_L(H, N, g_deph, g_down, g_up):
    """L = −i[H,·] + Σ_l γ^deph_l D[Z_l] + Σ_l γ↓_l D[σ⁻_l] + Σ_l γ↑_l D[σ⁺_l] (dense, vec form)."""
    c_ops = []
    for l in range(N):
        if g_deph[l]:
            c_ops.append(np.sqrt(g_deph[l]) * site_op(N, l, 'Z'))
    for l in range(N):
        if g_down[l]:
            c_ops.append(np.sqrt(g_down[l]) * site_2x2(N, l, SM))
    for l in range(N):
        if g_up[l]:
            c_ops.append(np.sqrt(g_up[l]) * site_2x2(N, l, SP))
    return lindbladian_general(H, c_ops)


def vec_max_mixed(d):
    return np.eye(d, dtype=complex).flatten() / d


def slope_dense(A, L, d):
    """d/dt Tr(A ρ)|_{ρ = I/d} from the dense Liouvillian."""
    return complex(np.vdot(A.conj().T.flatten(), L @ vec_max_mixed(d)))


def curvature_dense(A, L, d):
    """d²/dt² Tr(A ρ)|_{ρ = I/d} from the dense Liouvillian (L applied twice)."""
    return complex(np.vdot(A.conj().T.flatten(), L @ (L @ vec_max_mixed(d))))


def slope_law(H, N, j, dg):
    """The closed form (1/d) Σ_l Δγ_l t_j(l), t_j(l) = Tr(Z_l H^j)."""
    ts = t_moment(H, N, j)
    return sum(dg[l] * ts[l] for l in range(N)) / (2 ** N)


# ======================================================================
# BLOCK A -- pump directions (bit-exact; exact-square rates so √γ·√γ is float-exact).
# ======================================================================
def block_A_pump_directions():
    print("-" * 92)
    print("BLOCK A  pump directions  [D[σ⁻](I) = +γZ_l, D[σ⁺](I) = −γZ_l, BIT-EXACT]")
    print("-" * 92)
    N, d = 3, 8
    vecI = np.eye(d, dtype=complex).flatten()
    H0 = np.zeros((d, d), dtype=complex)
    zeros = [0.0] * N
    for l in range(N):
        for g in (0.25, 1.0, 2.25):                       # exact squares: sqrt(g)² == g in float
            gl = list(zeros)
            gl[l] = g
            out_dn = (build_L(H0, N, zeros, gl, zeros) @ vecI).reshape(d, d)
            out_up = (build_L(H0, N, zeros, zeros, gl) @ vecI).reshape(d, d)
            Zl = site_op(N, l, 'Z')
            assert np.max(np.abs(out_dn - g * Zl)) == 0.0, f"site {l} γ={g}: D[σ⁻](I) != +γZ_l"
            assert np.max(np.abs(out_up + g * Zl)) == 0.0, f"site {l} γ={g}: D[σ⁺](I) != −γZ_l"
        print(f"  site {l}: D[σ⁻](I) = +γZ_{l}, D[σ⁺](I) = −γZ_{l}  BIT-EXACT (γ ∈ {{¼, 1, 2¼}})")

    # unitary + dephasing kill I bit-exact (random Hermitian H, generic dephasing rates)
    rng = np.random.default_rng(7)
    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    Hr = (M + M.conj().T) / 2
    L_ud = build_L(Hr, N, [0.21, 0.34, 0.15], zeros, zeros)
    assert np.max(np.abs(L_ud @ vecI)) == 0.0, "unitary + dephasing do not kill I bit-exact"
    print("  unitary part + Z-dephasing kill I: max|L(I)| = 0.0 BIT-EXACT (random H)")

    # mixed channel: γ↓ and γ↑ together net Δγ_l = γ↓_l − γ↑_l
    L_mix = build_L(H0, N, [0.5] * N, [0.25, 0, 0], [1.0, 0, 0])
    out = (L_mix @ vecI).reshape(d, d)
    assert np.max(np.abs(out - (0.25 - 1.0) * site_op(N, 0, 'Z'))) == 0.0, \
        "mixed channel does not net Δγ·Z_l"
    print("  mixed channel γ↓ = ¼, γ↑ = 1 (+ dephasing): L(I) = (γ↓ − γ↑)·Z_0 = −¾·Z_0 BIT-EXACT")

    # builder cross-check: σ⁻-only path == framework lindbladian_z_plus_t1, bit-for-bit
    g_deph, g_dn = [0.21, 0.34, 0.15], [0.11, 0.27, 0.05]
    diff = np.max(np.abs(build_L(Hr, N, g_deph, g_dn, zeros)
                         - lindbladian_z_plus_t1(Hr, g_deph, g_dn)))
    assert diff == 0.0, f"build_L != lindbladian_z_plus_t1 (max diff {diff})"
    print("  cross-check: build_L(γ↑ = 0) == lindbladian_z_plus_t1 BIT-FOR-BIT")
    print("BLOCK A PASS")


# ======================================================================
# BLOCK B -- the slope law vs the dense Liouvillian.
# ======================================================================
def block_B_slope_law():
    print("-" * 92)
    print("BLOCK B  the slope law: d/dt⟨A⟩|_{I/d} = (1/d) Σ_l Δγ_l Tr(A Z_l)  [vs dense L, N=3]")
    print("-" * 92)
    N, d = 3, 8
    rng = np.random.default_rng(11)
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    dg = [a - b for a, b in zip(g_dn, g_up)]

    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    H_rand = (M + M.conj().T) / 2
    H_chain = _build_bilinear(N, [(0, 1), (1, 2)],
                              [('X', 'X', 1.0), ('Y', 'Y', 1.0), ('Z', 'Z', 1.0)]) \
        + 0.4 * site_op(N, 0, 'Z')
    for name, H in (("random Hermitian", H_rand), ("Heisenberg chain + Z-field", H_chain)):
        L = build_L(H, N, g_deph, g_dn, g_up)
        for j in range(1, 5):
            A = np.linalg.matrix_power(H, j)
            lhs = slope_dense(A, L, d)
            rhs = slope_law(H, N, j, dg)
            dev = abs(lhs - rhs)
            assert dev <= 1e-14 * max(1.0, abs(lhs)), f"{name} j={j}: dev {dev:.2e}"
            # cross-check t_j(l) against the f87 girth-dichotomy convention via an
            # independent evaluation path (einsum vs trace(matmul))
            Hp = np.linalg.matrix_power(H, j)
            for l in range(N):
                t_indep = complex(np.einsum('ij,ji->', site_op(N, l, 'Z'), Hp))
                assert abs(t_indep - t_moment(H, N, j)[l]) < 1e-10, \
                    f"{name} j={j} l={l}: t_moment cross-check"
            print(f"  {name:26s} j={j}: slope = {lhs.real:+.6f}  law = {rhs.real:+.6f}  "
                  f"dev = {dev:.2e}  (t_j == f87 t_moment)")
        # random Hermitian measurement A (not a power of H): the law is Tr(A Z_l)-linear
        Ma = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        A = (Ma + Ma.conj().T) / 2
        lhs = slope_dense(A, L, d)
        rhs = sum(dg[l] * np.trace(A @ site_op(N, l, 'Z')) for l in range(N)) / d
        dev = abs(lhs - rhs)
        assert dev <= 1e-14, f"{name} random A: dev {dev:.2e}"
        print(f"  {name:26s} random Hermitian A: dev = {dev:.2e}")
    print("BLOCK B PASS")


# ======================================================================
# BLOCK C -- the three blindnesses (bit-identical slopes; exact zero at detailed balance).
# ======================================================================
def block_C_blindnesses():
    print("-" * 92)
    print("BLOCK C  the three blindnesses  [deph-blind, evolution-blind, detailed-balance closure]")
    print("-" * 92)
    N, d = 3, 8
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    H = site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')
    L = build_L(H, N, g_deph, g_dn, g_up)

    # (1) dephasing-blind: change all dephasing rates, slope BIT-IDENTICAL
    L_deph = build_L(H, N, [0.9, 0.01, 0.5], g_dn, g_up)
    # (2) evolution-blind: change the evolution H entirely (the measured polynomial A stays
    #     the PROGRAMMED H^j; only the pump enters the slope)
    H2 = H + 2.3 * site_op(N, 0, 'Y') @ site_op(N, 1, 'Y') + 0.9 * site_op(N, 2, 'Z')
    L_evo = build_L(H2, N, g_deph, g_dn, g_up)
    for j in (1, 2, 3):
        A = np.linalg.matrix_power(H, j)
        s0 = slope_dense(A, L, d)
        assert slope_dense(A, L_deph, d) == s0, f"j={j}: dephasing-rate change moved the slope"
        assert slope_dense(A, L_evo, d) == s0, f"j={j}: evolution-H change moved the slope"
    print("  dephasing rates changed         → slope BIT-IDENTICAL (j = 1, 2, 3)")
    print("  evolution H changed entirely    → slope BIT-IDENTICAL (j = 1, 2, 3)")

    # (3) detailed balance γ↑ = γ↓ (exact-square rates): slope exactly zero
    g_bal = [0.25, 1.0, 2.25]
    L_bal = build_L(H, N, g_deph, g_bal, g_bal)
    for j in (1, 2, 3):
        A = np.linalg.matrix_power(H, j)
        assert slope_dense(A, L_bal, d) == 0.0, f"j={j}: detailed balance slope != 0"
    print("  detailed balance γ↑ = γ↓        → slope = 0.0 EXACT (j = 1, 2, 3)")
    print("BLOCK C PASS")


# ======================================================================
# BLOCK D -- the girth certificate (one-sided) + the honest negative control.
# ======================================================================
def block_D_girth_certificate():
    print("-" * 92)
    print("BLOCK D  the girth certificate  [first firing j = girth ℓ; silence is NOT softness]")
    print("-" * 92)
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    dg = [a - b for a, b in zip(g_dn, g_up)]

    # (i) girth-1 witness: Z-drive H fires at j = 1
    N, d = 3, 8
    H1 = 0.5 * site_op(N, 0, 'Z') + 0.3 * site_op(N, 1, 'Z')
    L1 = build_L(H1, N, g_deph, g_dn, g_up)
    s1 = slope_dense(H1, L1, d)
    law1 = slope_law(H1, N, 1, dg)
    assert abs(s1) > 1e-3 and abs(s1 - law1) <= 1e-14, "girth-1 witness does not fire at j=1"
    print(f"  (i)   girth-1 (Z-drive): slope(H^1) = {s1.real:+.6f} fires at j = 1 (== law)")

    # the site-sum caveat: balanced ω with UNIFORM Δγ cancels the site sum although t_1(l) != 0;
    # site-resolved weights (selective per-site damping) recover the certificate
    Hc = 0.5 * site_op(N, 0, 'Z') - 0.5 * site_op(N, 1, 'Z')
    L_unif = build_L(Hc, N, g_deph, [0.25] * N, [0.0] * N)
    L_site = build_L(Hc, N, g_deph, [0.25, 0.0, 0.0], [0.0] * N)
    s_unif = slope_dense(Hc, L_unif, d)
    s_site = slope_dense(Hc, L_site, d)
    ts = t_moment(Hc, N, 1)
    assert ts[0].real == 4.0 and ts[1].real == -4.0, "caveat witness t_1 unexpected"
    assert s_unif == 0.0, f"uniform-weight site sum should cancel exactly ({s_unif})"
    assert abs(s_site - 0.25 * 4.0 / d) <= 1e-15 and abs(s_site) > 1e-3, \
        "site-resolved weight does not recover t_1(0)"
    print(f"  (i)   site-sum caveat: t_1 = [+4, −4, 0] cancels under uniform Δγ "
          f"(slope = {s_unif.real:.1f} exactly); per-site Δγ resolves it "
          f"(slope = {s_site.real:+.4f} = Δγ_0·t_1(0)/d)")

    # (ii) girth-2 witness: H = X_0 + X_0 Z_1 + 0.7·X_1 X_2; t_1 ≡ 0, t_2 fires at site 1
    H2 = site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')
    L2 = build_L(H2, N, g_deph, g_dn, g_up)
    t1 = t_moment(H2, N, 1)
    t2 = t_moment(H2, N, 2)
    assert all(t == 0 for t in t1), f"girth-2 witness: t_1 = {t1} != 0"
    assert t2[0] == 0 and t2[2] == 0 and t2[1] == 16.0, f"girth-2 witness: t_2 = {t2}"
    s_j1 = slope_dense(H2, L2, d)
    s_j2 = slope_dense(H2 @ H2, L2, d)
    assert s_j1 == 0.0, f"girth-2 witness: slope fired below the girth ({s_j1})"
    assert abs(s_j2 - dg[1] * 16.0 / d) <= 1e-14 and abs(s_j2) > 1e-3, \
        "girth-2 witness: slope(H^2) does not fire with Δγ_1·16/d"
    print(f"  (ii)  girth-2 (X₀ + X₀Z₁ + 0.7·X₁X₂): t_1 ≡ 0, t_2 = [0, 16, 0]; "
          f"slope(H^1) = 0.0 EXACT, slope(H^2) = {s_j2.real:+.4f} = Δγ_1·16/d  → ℓ = 2, m* = 5")

    # (iii) HONEST negative control: the k=4 pair IIXY+ZXZY (N = 5). The deg-1 tower is silent
    # through j = 5, yet the pair is HARD at m* = 11 = 2ℓ+5 with p_11 = 86507520·γ⁵ (deg-5;
    # pinned exactly in f87_girth_dichotomy Block 3). Silence certifies NOTHING about softness.
    N5, d5 = 5, 32
    H4 = _build_kbody_chain(N5, [('I', 'I', 'X', 'Y', 1.0), ('Z', 'X', 'Z', 'Y', 1.0)])
    L4 = build_L(H4, N5, [0.2] * N5, [0.11, 0.27, 0.05, 0.18, 0.07], [0.02] * N5)
    for j in range(1, 6):
        ts = t_moment(H4, N5, j)
        assert max(abs(t) for t in ts) == 0.0, f"IIXY+ZXZY: t_{j} != 0"
        s = slope_dense(np.linalg.matrix_power(H4, j), L4, d5)
        assert abs(s) <= 1e-12, f"IIXY+ZXZY: slope fired at j={j} ({abs(s):.2e})"
    print("  (iii) HONEST negative IIXY+ZXZY (N=5): t_j = 0 EXACTLY and slope silent for "
          "j = 1..5,")
    print("        yet the pair is HARD at m* = 11 = 2ℓ+5, p₁₁ = 86507520·γ⁵ (deg-5).")
    print("        The deg-1 tower's silence is NOT a softness certificate; the certificate "
          "is one-sided.")
    print("BLOCK D PASS")


# ======================================================================
# BLOCK E -- the F113 bridge: static polarity asymmetry == the dynamic pump slope.
# ======================================================================
def block_E_f113_bridge():
    print("-" * 92)
    print("BLOCK E  the F113 bridge  [asymmetry = (4^N/2)Σω(γ_pump − γ_T1) = −4^N·slope⟨H⟩]")
    print("-" * 92)
    for N in (2, 3):
        d = 2 ** N
        omega = [0.13, 0.29, 0.07][:N]
        gt1 = [0.10, 0.04, 0.16][:N]
        gpu = [0.03, 0.11, 0.02][:N]
        H = sum((omega[l] / 2.0) * site_op(N, l, 'Z') for l in range(N))

        # (i) the F113 closed form
        closed = (4 ** N / 2.0) * sum(omega[l] * (gpu[l] - gt1[l]) for l in range(N))
        # (ii) −4^N · slope⟨H⟩ from the DENSE Liouvillian (the dynamic side)
        L = build_L(H, N, [0.0] * N, gt1, gpu)
        bridge = -(4 ** N) * slope_dense(H, L, d).real
        # (iii) the ACTUAL Frobenius asymmetry from the framework polarity decomposition
        c_ops = [site_2x2(N, l, SM) for l in range(N)] + [site_2x2(N, l, SP) for l in range(N)]
        actual = polarity_coordinates_from_hc(H, c_ops, list(gt1) + list(gpu), N)['asymmetry']

        assert abs(closed - bridge) <= 1e-15, f"N={N}: closed vs slope-bridge {closed} {bridge}"
        assert abs(closed - actual) <= 1e-12, f"N={N}: closed vs Frobenius {closed} {actual}"
        print(f"  N={N}: closed form = {closed:+.12f}")
        print(f"        −4^N·slope  = {bridge:+.12f}   (|Δ| = {abs(closed - bridge):.1e})")
        print(f"        ‖M₊‖²−‖M₋‖² = {actual:+.12f}   (|Δ| = {abs(closed - actual):.1e})")

        # the bridge is dephasing-blind on BOTH sides: add Z-dephasing, nothing moves
        c_ops_z = [site_op(N, l, 'Z') for l in range(N)] + c_ops
        actual_z = polarity_coordinates_from_hc(H, c_ops_z, [0.05] * N + list(gt1) + list(gpu),
                                                N)['asymmetry']
        L_z = build_L(H, N, [0.05] * N, gt1, gpu)
        bridge_z = -(4 ** N) * slope_dense(H, L_z, d).real
        assert bridge_z == bridge, f"N={N}: slope side moved under added dephasing"
        assert abs(actual_z - actual) <= 1e-12, f"N={N}: Frobenius side moved under dephasing"
        print(f"        + Z-dephasing γ = 0.05: both sides unchanged (slope side bit-identical)")
    print("BLOCK E PASS")


# ======================================================================
# BLOCK F -- the curvature fingerprint: parasite reading on the commutator probes.
# ======================================================================
def block_F_curvature_fingerprint():
    print("-" * 92)
    print("BLOCK F  the curvature fingerprint  [diff = (δ/d) Σ_l Δγ_l (−i)Tr(V[Z_l, A])]")
    print("-" * 92)
    N, d = 3, 8
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    dg = [a - b for a, b in zip(g_dn, g_up)]
    H_p = site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')
    L_p = build_L(H_p, N, g_deph, g_dn, g_up)

    def probe_pred(V, A, delta):
        """(δ/d) Σ_l Δγ_l (−i) Tr(V [Z_l, A]): the parasite read on the commutator probes."""
        out = 0.0 + 0.0j
        for l in range(N):
            Zl = site_op(N, l, 'Z')
            out += dg[l] * (-1j) * np.trace(V @ (Zl @ A - A @ Zl))
        return delta * out / d

    # exact affinity in the generator: the curvature difference doubles bit-exactly with δ
    V = site_op(N, 0, 'Y')
    A = H_p @ H_p
    diff_1 = curvature_dense(A, build_L(H_p + 0.05 * V, N, g_deph, g_dn, g_up), d) \
        - curvature_dense(A, L_p, d)
    diff_2 = curvature_dense(A, build_L(H_p + 0.10 * V, N, g_deph, g_dn, g_up), d) \
        - curvature_dense(A, L_p, d)
    assert diff_2 == 2 * diff_1, "curvature difference not exactly affine in δ"
    print(f"  affinity: diff(2δ) == 2·diff(δ) BIT-EXACT (δ = 0.05, V = Y₀, A = H_p²)")

    # X/Y-flavored parasites with probe overlap are read LINEARLY
    delta = 0.05
    for name, V in (('Y₀', site_op(N, 0, 'Y')),
                    ('Y₀Z₁', site_op(N, 0, 'Y') @ site_op(N, 1, 'Z')),
                    ('Y₁X₂', site_op(N, 1, 'Y') @ site_op(N, 2, 'X'))):
        L_c = build_L(H_p + delta * V, N, g_deph, g_dn, g_up)
        diff = curvature_dense(H_p, L_c, d) - curvature_dense(H_p, L_p, d)
        pred = probe_pred(V, H_p, delta)
        assert abs(diff - pred) <= 1e-15, f"{name}: probe law dev {abs(diff - pred):.2e}"
        assert abs(pred) > 1e-3, f"{name}: parasite not actually read"
        print(f"  parasite {name:4s}: curvature diff = {diff.real:+.6f} = probe prediction "
              f"(dev = {abs(diff - pred):.1e}, read linearly)")

    # Z-flavored parasite: EXACTLY invisible ([Z, Z_l] = 0)
    V_z = site_op(N, 0, 'Z')
    for j in (1, 2):
        A = np.linalg.matrix_power(H_p, j)
        L_c = build_L(H_p + delta * V_z, N, g_deph, g_dn, g_up)
        diff = curvature_dense(A, L_c, d) - curvature_dense(A, L_p, d)
        assert diff == 0.0, f"Z₀ parasite visible at j={j} ({diff})"
    print("  parasite Z₀  : curvature diff = 0.0 BIT-EXACT (j = 1, 2)")
    print("  → F113 complementarity: the pump curvature is blind exactly where F113's balance")
    print("    channel reads (the Z-drives); X/Y parasites land on the [Z_l, H_p^j] probes.")
    print("BLOCK F PASS")


# ======================================================================
# BLOCK G -- the finite-time protocol (the practical reading).
# ======================================================================
def block_G_finite_time():
    print("-" * 92)
    print("BLOCK G  the finite-time protocol  [expm propagation; curvature = first correction]")
    print("-" * 92)
    # Preparation note: ρ(0) = I/d is the average over the 2^N computational basis states,
    # so the protocol is: prepare each basis state, evolve for t, measure A, average.
    N, d = 3, 8
    g_deph = [0.21, 0.34, 0.15]
    g_dn = [0.11, 0.27, 0.05]
    g_up = [0.02, 0.09, 0.13]
    H = site_op(N, 0, 'X') + site_op(N, 0, 'X') @ site_op(N, 1, 'Z') \
        + 0.7 * site_op(N, 1, 'X') @ site_op(N, 2, 'X')
    L = build_L(H, N, g_deph, g_dn, g_up)
    A = H @ H
    vecI = vec_max_mixed(d)
    a_vec = A.conj().T.flatten()
    f0 = np.vdot(a_vec, vecI).real
    s_true = slope_dense(A, L, d).real
    c_true = curvature_dense(A, L, d).real

    fd = {}
    for t in (0.005, 0.01, 0.02, 0.04):
        f_t = np.vdot(a_vec, expm(L * t) @ vecI).real
        fd[t] = (f_t - f0) / t
        corr = fd[t] - s_true
        pred_corr = (t / 2.0) * c_true
        # the curvature is the first correction: the remainder after subtracting it is O(t²)
        assert abs(corr - pred_corr) < 0.05 * t * t * max(1.0, abs(c_true)), \
            f"t={t}: remainder after the curvature correction not O(t²)"
        print(f"  t = {t:5.3f}: FD slope = {fd[t]:+.6f}   FD − law = {corr:+.3e}   "
              f"(t/2)·curv = {pred_corr:+.3e}   ratio = {corr / pred_corr:.4f}")
    assert abs(fd[0.005] - s_true) < abs(fd[0.04] - s_true), "window error does not shrink with t"

    # Richardson extrapolation 2·FD(t/2) − FD(t) cancels the curvature term
    rich = 2 * fd[0.01] - fd[0.02]
    assert abs(rich - s_true) < abs(fd[0.01] - s_true) / 10, "Richardson does not cancel curvature"
    print(f"  Richardson 2·FD(0.01) − FD(0.02) = {rich:+.6f} vs law {s_true:+.6f} "
          f"(|Δ| = {abs(rich - s_true):.1e}: the curvature term cancels)")
    print("  (preparation: I/d = average over the 2^N computational basis states)")
    print("BLOCK G PASS")


# ======================================================================
# BLOCK H -- the detailed-balance closure + the F84 tie.
# ======================================================================
def block_H_detailed_balance_f84():
    print("-" * 92)
    print("BLOCK H  detailed-balance closure + F84 tie  [Δγ = 0 ⟹ slope ≡ 0; Δγ_l IS F84's rate]")
    print("-" * 92)
    N, d = 3, 8
    rng = np.random.default_rng(23)
    M = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    H = (M + M.conj().T) / 2
    vecI = np.eye(d, dtype=complex).flatten()

    # γ↓ = γ↑ per site (exact-square rates): L(I) = 0 bit-exact, so EVERY slope is zero exactly
    g_bal = [0.25, 1.0, 2.25]
    L_bal = build_L(H, N, [0.21, 0.34, 0.15], g_bal, g_bal)
    assert np.max(np.abs(L_bal @ vecI)) == 0.0, "detailed balance: L(I) != 0 bit-exact"
    for j in (1, 2, 3, 4):
        assert slope_dense(np.linalg.matrix_power(H, j), L_bal, d) == 0.0, \
            f"detailed balance: slope(H^{j}) != 0"
    Ma = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
    A = (Ma + Ma.conj().T) / 2
    assert slope_dense(A, L_bal, d) == 0.0, "detailed balance: slope(random A) != 0"
    print("  γ↓ = γ↑ = (¼, 1, 2¼) per site: L(I) = 0 BIT-EXACT → slope ≡ 0 EXACT "
          "(H^1..H^4 + random A)")
    # generic (non-exact-square) balanced rates: float-fuzz level only
    L_bal2 = build_L(H, N, [0.21, 0.34, 0.15], [0.3, 0.1, 0.2], [0.3, 0.1, 0.2])
    s = abs(slope_dense(H, L_bal2, d))
    assert s < 1e-15, f"generic balanced rates: slope {s:.2e}"
    print(f"  generic balanced rates (0.3, 0.1, 0.2): slope = {s:.1e} (√γ·√γ float fuzz only)")

    # the F84 tie: the pump weight Δγ_l IS F84's net vacuum rate (the same quantity drives the
    # F81 violation √(Σ_l Δγ_l²)·2^(N−1); see docs/proofs/PROOF_F84_AMPLITUDE_DAMPING.md):
    # detailed balance closes BOTH readings simultaneously, vacuum-only opens both.
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chain = fw.ChainSystem(N=N)
    v_bal = fw.predict_amplitude_damping_violation(chain, g_bal, g_bal)
    assert v_bal == 0.0, "F84 tie: balanced rates should zero the F84 violation"
    g_dn, g_up = [0.11, 0.27, 0.05], [0.02, 0.09, 0.13]
    dg = [a - b for a, b in zip(g_dn, g_up)]
    v_net = fw.predict_amplitude_damping_violation(chain, g_dn, g_up)
    v_expect = float(np.sqrt(sum(x * x for x in dg)) * 2 ** (N - 1))
    assert abs(v_net - v_expect) < 1e-12, "F84 tie: net violation != √(ΣΔγ²)·2^(N−1)"
    print(f"  F84 tie: Δγ_l = (γ↓ − γ↑)_l drives BOTH the pump slope and the F84 violation")
    print(f"           balanced → F84 violation = 0.0;  net Δγ = {tuple(round(x, 2) for x in dg)}"
          f" → violation = {v_net:.4f} = √(ΣΔγ²)·2^(N−1)")
    print("BLOCK H PASS")


def main():
    print("=" * 92)
    print("THE MOMENT-TOWER PUMP CHANNEL -- the deg-1 girth ladder read linearly (F120)")
    print("=" * 92)
    block_A_pump_directions()
    block_B_slope_law()
    block_C_blindnesses()
    block_D_girth_certificate()
    block_E_f113_bridge()
    block_F_curvature_fingerprint()
    block_G_finite_time()
    block_H_detailed_balance_f84()
    print("=" * 92)
    print("ALL BLOCKS PASS")
    print("=" * 92)


if __name__ == "__main__":
    main()
