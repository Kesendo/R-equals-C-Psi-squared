"""Gate for the zeta^2 anti-protection law (the derived form, 2026-07-16).

THE LAW (docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md): for a Theta-mirror
collision pair of 3-magnon branches (tau, nu = N+1-tau) of the brickwork XX
Floquet step U(theta) on the open N-site chain, under ANY occupation-diagonal
perturbation phase A (site-dependent ZZ of any range, single-site detunings),
the quasi-energies obey exactly

    theta_tau(zeta) = -theta_nu(-zeta)        (all orders, isolated branch),

so the branch DIFFERENCE is EVEN in a global sign flip of zeta: odd orders
transfer to both branches equally (first-order protection), even orders push
the branches exactly apart (anti-protection, exact factor 2). The owner of the
static symmetry is ChiralKClaim (class BDI, K = diag((-1)^i), PROOF_K_
PARTNERSHIP.md); Theta = T.K is its Floquet-step lift. The symmetry fixes the
evenness and the factor 2, NOT the sign; one-sidedness of the flown b_zz2
budget is a same-sign-regime property and stays in the instrument doc
(experiments/IBM_F129_RAMSEY_FRINGE.md section 5).

GATES (exit 0 iff all PASS):
  G0  the lift itself: C U(theta)* C^-1 = U(theta) machine-zero; the mode comb
      is exactly antisymmetric and Theta realizes k <-> N+1-k (n = N+1 = 9),
      sending Slater (1,5,7) to (2,4,8).
  G1  CITED, not owned: chiral_certificate() from f129_ramsey_fringe_design
      (committed 2026-07-15) already certifies the first-order protection
      (all-range Z-diagonal equality, site-dependent) and the Slater-Condon
      branch-mixing zero. Re-run here as an imported precondition.
  G2  the clean-limit coefficient: the exact quasi-energy branch-difference
      bias is quadratic in a global zeta scale; its zeta -> 0 limit
      C_phys ~ 0.002522 (seed profile, tau = 1.2 us reference) is reproduced
      by second-order Floquet PT (the cot formula) to < 1e-3 relative.
      C_phys is the physical constant the flown, estimator-dressed
      C2LAW = 0.00257 (gate mode zz2_scan of f129_ramsey_7a_gate.py) dresses.
  G3  exact PT2 antisymmetry: eps2_tau + eps2_nu = 0 at machine precision for
      uniform, mirror-symmetric AND asymmetric random zeta profiles.
  G4  evenness: bias(+zeta) = bias(-zeta) exactly (two profiles).
  G5  Trotter-split consistency: eps_tau - eps_nu = 2 eps_tau (eps_nu =
      -eps_tau exactly), and it equals the PREDICTED A0 drift center +0.0213
      (fringe doc section 5; the measured flown slope is +0.0326).

FENCES (load-bearing, from the design review): the law exists ONLY for
Theta-mirror pairs; non-mirror collision pairs (A1, A2, and the n = 12 pair
(1,2,10)~(3,5,6)) carry FIRST-order ZZ shifts. Open chain + NN hopping only
(the fermionic compound is exact there; a ring or beyond-NN hopping
reintroduces JW strings). The pt2 helper below is fenced to the
Slater-Condon-protected pair (globally the comb has exactly degenerate pairs
with REAL coupling where the formula's degeneracy skip would silently eat an
O(zeta) term); it asserts the S-C zero before summing.

Runtime ~20 s. Doc: docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md.
"""

import sys
from itertools import combinations

import numpy as np

sys.path.insert(0, "simulations")
from f129_ramsey_fringe_design import (chiral_certificate, floquet,  # noqa: E402
                                       one_particle_step, phase_sum)

N = 8                       # sites; the comb index is n = N + 1 = 9
THETA = 0.5                 # the flown hop angle per step
TAU_REF_US = 1.2            # the law's normalization point (zz2_scan ref)
ZETA_BASE_KHZ, ZETA_SPREAD = 3.8, 0.30
SEED = 20260715             # the 7a gate's seed zeta profile, reproduced
TAU_MODES, NU_MODES = (1, 5, 7), (2, 4, 8)
C_PHYS_PIN = 0.002522       # the clean-limit coefficient (this gate derives it)

FAILURES = []


def check(name, ok, detail):
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}: {detail}")
    if not ok:
        FAILURES.append(name)


CONFIGS = list(combinations(range(N), 3))
IDX = {c: i for i, c in enumerate(CONFIGS)}
D3 = len(CONFIGS)


def compound3(A):
    """Third compound matrix (exact 3-fermion sector, open chain + NN only)."""
    M = np.empty((D3, D3), dtype=complex)
    for a, I in enumerate(CONFIGS):
        AI = A[np.ix_(I, range(N))]
        for b, J in enumerate(CONFIGS):
            M[a, b] = np.linalg.det(AI[:, J])
    return M


def zz_diag_phase(zetas_khz, tau_us):
    """Occupation-diagonal per-step phase: sum_bonds pi zeta_i tau z_i z_{i+1}
    (the rzz(2 pi zeta tau) convention of the 7a gate), kHz and us -> rad."""
    out = np.empty(D3)
    for a, I in enumerate(CONFIGS):
        z = np.ones(N)
        for s in I:
            z[s] = -1.0
        out[a] = sum(np.pi * zetas_khz[i] * tau_us * 1e-3 * z[i] * z[i + 1]
                     for i in range(N - 1))
    return out


def main():
    ph, V = floquet(THETA, N)
    U1 = one_particle_step(THETA, N)
    U3 = compound3(U1)
    W3 = compound3(V)
    assert list(combinations(range(N), 3)) == CONFIGS
    eps0 = np.array([sum(ph[list(t)]) for t in CONFIGS])
    s_tau = IDX[tuple(k - 1 for k in TAU_MODES)]
    s_nu = IDX[tuple(k - 1 for k in NU_MODES)]

    rng = np.random.default_rng(SEED)
    zeta_seed = ZETA_BASE_KHZ * (1 + ZETA_SPREAD * (2 * rng.random(N - 1) - 1))

    # ---------------- G0: the Floquet lift of the chiral symmetry ----------
    print("G0  Theta = T.K commutes with the brickwork step (the lift)")
    C = np.diag([(-1.0) ** i for i in range(N)])
    lift_res = np.max(np.abs(C @ U1.conj() @ C - U1))
    check("C U* C^-1 = U", lift_res < 1e-14, f"residual {lift_res:.2e}")
    comb_res = np.max(np.abs(ph + ph[::-1]))
    check("mode comb antisymmetric (k <-> N+1-k)", comb_res < 1e-12,
          f"max |ph_k + ph_(9-k)| = {comb_res:.2e}")
    # Theta maps the tau Slater to the nu Slater (overlap modulus 1)
    theta_tau = compound3(C) @ W3[:, s_tau].conj()
    ov = abs(np.vdot(W3[:, s_nu], theta_tau))
    check("Theta(Slater 1,5,7) = Slater (2,4,8)", abs(ov - 1) < 1e-12,
          f"|overlap| = {ov:.12f}")

    # ---------------- G1 (cited): the committed chiral certificate ---------
    print("G1  cited precondition: chiral_certificate() "
          "(f129_ramsey_fringe_design, committed)")
    occ, bd, mix = chiral_certificate(V, TAU_MODES, NU_MODES, N)
    check("first-order protection + Slater-Condon zero",
          max(occ, bd, mix) < 1e-12,
          f"occ {occ:.2e}, bond {bd:.2e}, mixing {mix:.2e}")

    # ---------------- exact bias machinery ---------------------------------
    def exact_bias(zetas, tau_us):
        Dph = zz_diag_phase(zetas, tau_us)
        Wstep = np.diag(np.exp(-1j * Dph)) @ U3
        w, X = np.linalg.eig(Wstep)
        th = np.angle(w)
        t_i = np.argmax(np.abs(X.conj().T @ W3[:, s_tau]) ** 2)
        n_i = np.argmax(np.abs(X.conj().T @ W3[:, s_nu]) ** 2)
        d0 = eps0[s_tau] - eps0[s_nu]
        return np.angle(np.exp(1j * (th[t_i] - th[n_i] - d0)))

    def pt2_pair_difference(zetas, tau_us):
        """Second-order Floquet PT branch shifts (eps2_tau, eps2_nu).
        FENCED to the Slater-Condon-protected pair: asserts the degenerate
        partner coupling is zero; do NOT reuse on other pairs (the comb has
        exact degeneracies with real coupling where the skip eats O(zeta))."""
        A = zz_diag_phase(zetas, tau_us)
        Amat = (W3.conj().T * A) @ W3
        assert abs(Amat[s_nu, s_tau]) < 1e-12, "S-C fence violated"

        def pt2(s):
            tot = 0.0
            for m in range(D3):
                if m == s:
                    continue
                g = np.angle(np.exp(1j * (eps0[s] - eps0[m])))
                if abs(g) < 1e-12:
                    assert abs(Amat[m, s]) < 1e-12, "degenerate coupling"
                    continue
                tot += 0.5 * abs(Amat[m, s]) ** 2 / np.tan(g / 2)
            return tot
        return pt2(s_tau), pt2(s_nu)

    # ---------------- G2: clean-limit coefficient = Floquet PT2 ------------
    print("G2  the clean-limit coefficient and its PT2 derivation")
    s_small = 0.0625
    C_phys = exact_bias(zeta_seed * s_small, TAU_REF_US) / s_small ** 2
    e2t, e2n = pt2_pair_difference(zeta_seed, TAU_REF_US)
    rel = abs((e2t - e2n) - C_phys) / abs(C_phys)
    check("bias quadratic, clean limit near pin",
          abs(C_phys - C_PHYS_PIN) < 5e-5,
          f"C_phys = {C_phys:+.6f} (pin {C_PHYS_PIN:+.6f}; "
          f"flown estimator-dressed C2LAW = +0.00257)")
    check("PT2 reproduces the clean limit", rel < 1e-3,
          f"PT2 diff = {e2t - e2n:+.7f}, relative error {rel:.1e}")

    # ---------------- G3: exact PT2 antisymmetry ----------------------------
    print("G3  PT2 shifts exactly opposite, arbitrary profiles")
    profiles = [("uniform 3.8", np.full(N - 1, 3.8)),
                ("seed (asymmetric +-30%)", zeta_seed),
                ("random asymmetric", np.random.default_rng(7).uniform(1, 8, N - 1)),
                ("random asymmetric 2", np.random.default_rng(11).uniform(1, 8, N - 1))]
    worst = 0.0
    for label, z in profiles:
        a, b = pt2_pair_difference(z, TAU_REF_US)
        worst = max(worst, abs(a + b))
    check("eps2_tau + eps2_nu = 0 (4 profiles)", worst < 1e-14,
          f"worst |sum| = {worst:.2e}")

    # ---------------- G4: evenness under zeta -> -zeta ----------------------
    print("G4  the branch difference is even in a global sign flip")
    worst = 0.0
    for z in (zeta_seed, np.random.default_rng(7).uniform(1, 8, N - 1)):
        worst = max(worst, abs(exact_bias(z, TAU_REF_US)
                               - exact_bias(-z, TAU_REF_US)))
    check("bias(+z) = bias(-z)", worst < 1e-13, f"worst diff = {worst:.2e}")

    # ---------------- G5: Trotter-split consistency --------------------------
    print("G5  the unperturbed split is the predicted A0 drift center, "
          "and eps_nu = -eps_tau")
    et, en = phase_sum(TAU_MODES, ph), phase_sum(NU_MODES, ph)
    check("eps_nu = -eps_tau exactly", abs(et + en) < 1e-13,
          f"eps_tau {et:+.6f}, eps_nu {en:+.6f}")
    check("split = predicted A0 drift center", abs((et - en) - 0.0213) < 2e-4,
          f"split {et - en:+.5f} vs +0.0213 (doc section 5)")

    print()
    if FAILURES:
        print(f"GATE FAIL: {FAILURES}")
        return 1
    print("GATE PASS: the zeta^2 anti-protection law holds as derived "
          "(Theta-mirror pairs only; see the fences in the header).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
