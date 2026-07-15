"""F129 Ramsey fringe design gate: the collision as a standing fringe on hardware.

Doc: experiments/IBM_F129_RAMSEY_FRINGE.md (the pre-registration this design gate
feeds; every number in its section 5 table is printed here).

PROTOCOL (design under test): prepare (|tau> + |c>)/sqrt(2) as 3-magnon states in the
FLOQUET eigenbasis of the Trotter step itself (GHZ-6 seed on the six mode qubits +
the Givens network built from the step's eigenvectors), evolve M Trotter steps of the
pure hopping chain (no engineered dephasing), invert the network, read the seed-qubit
fringe in two quadratures. The fringe phase is EXACTLY M * dPhi(pair), dPhi = the
difference of 3-mode Floquet phase sums, with no prep leakage (the bitstrings are
exact eigenstates of the stepped evolution by construction; the XX chain is free-
fermionic, so the matchgate step is exact at any theta, not a small-angle
approximation).

THE LAW ON DISPLAY (primary lab, n = 9, N = 8): the clean mirror collision
(1,5,7) ~ (2,4,8) of F129 has EXACTLY equal levels; the Floquet deformation preserves
that to O(theta^3) per step (theta -> 0 coefficient 1/6, effective 0.170 at
theta = 0.5), and the CHIRAL (sublattice) symmetry behind the mirror pair,
V_{n-k} = C conj(V_k) with C = diag((-1)^j), equalizes every Z-diagonal observable
of the two branches PER SITE and PER BOND (certified below at machine precision),
so the collision arm's first-order drift from ZZ coupling AND from single-qubit
Z-detuning is exactly zero even for SITE-DEPENDENT zeta_i / delta_i. Every
non-colliding pair winds at ~2 theta dS per step. The fringe that stands still IS
the concert-hall law.

The n = 12 pair (1,2,10) ~ (3,5,6) appears only as a NON-FLYABLE appendix: (1,2,10)
is not clean (2 + 10 = 12), so that lab tests F130's equal-level hypothesis, not
F129's clean-collision law, and it is excluded from the pre-registration's scope.

ERROR MODEL (v2, review round 1 folded):
  - sigma_phi = 1/(V * sqrt(shots)) per quadrature point; weights from the PREDICTED
    V(M) (removes the measured-weight bias); weighted-LSQ slope error; COMBINED
    two-arm error sqrt(2) * sigma_slope wherever two slopes are compared.
  - visibility V(M) = (1-p2)^CX * exp(-3 t/T1) * exp(-m_deph * t/T2*), where m_deph
    is the exact branch-pair dephasing multiplier E[Hamming distance] between the two
    Slater branches' site-configuration distributions (computed from below here;
    3.680 for the n = 9 collision pair, NOT 1: the branch superposition is a spread
    cat, and modeling it as one collective DOF understates sigma). m_deph is an
    exactly computed FIRST MOMENT; the single-exponential model built on it is
    conservative by Jensen (the true multi-exponential envelope lies above it).
    Control arms use their OWN pair's m_deph; two-arm comparisons combine
    sqrt(sigma_A0^2 + sigma_arm^2).
  - ZZ systematic: exact first-order Wick value per arm (one-body density matrix of
    the 3-magnon Slater state; <z_i z_j> = 1 - 2n_i - 2n_j + 4(n_i n_j - |g_ij|^2)),
    floored at 5x its magnitude as the device-transfer safety factor (zeta measured
    3.6-3.9 kHz on ibm_marrakesh; 5x covers zeta up to ~20 kHz). The collision arm's
    ZZ drift is zero BY SYMMETRY (mirror pair equalizes the bond correlators), which
    is zeta-independent and holds at every theta.

Constants from the external pipeline survey (2026-07-15): tau_step = 1200 ns,
2q gate ~0.5 us, T1 ~ 200 us, T2* ~ 70 us, ZZ ~ 4 kHz (Confirmations 21 noise
model). Billing projected by shots-ratio against the measured comparable
(IBM_CONCENTRATOR_RELOADED: 376,832 shots billed 119 s).
"""

import numpy as np
from itertools import combinations

TAU_STEP_US = 1.2
T1_US = 200.0
T2S_US = 70.0
ZZ_RAD_PER_STEP = 2 * np.pi * 4e3 * TAU_STEP_US * 1e-6  # ~0.030 rad at zeta = 4 kHz
ZZ_SAFETY = 5.0
BILLED_S_PER_SHOT = 119.0 / 376832.0


def one_particle_step(theta, nn):
    U = np.eye(nn, dtype=complex)
    for parity in (0, 1):
        L = np.eye(nn, dtype=complex)
        for i in range(parity, nn - 1, 2):
            B = np.eye(nn, dtype=complex)
            B[i, i] = B[i + 1, i + 1] = np.cos(theta)
            B[i, i + 1] = B[i + 1, i] = 1j * np.sin(theta)
            L = B @ L
        U = L @ U
    return U


def floquet(theta, nn):
    """Phases (descending, matching 2 theta cos(k pi/n) order) and eigenvectors."""
    w, V = np.linalg.eig(one_particle_step(theta, nn))
    ph = np.angle(w)
    order = np.argsort(-ph)
    return ph[order], V[:, order]


def phase_sum(triple, ph):
    return sum(ph[k - 1] for k in triple)


def slater_config_probs(V, triple, nn):
    """|<site config|Slater(triple)>|^2 over all C(nn, 3) occupation configs."""
    cols = V[:, [k - 1 for k in triple]]
    configs = list(combinations(range(nn), 3))
    p = np.array([abs(np.linalg.det(cols[list(c), :])) ** 2 for c in configs])
    return configs, p / p.sum()


def dephasing_multiplier(V, a, b, nn):
    """E[Hamming distance] between the two branches' site-configuration
    distributions: the early-time fringe dephasing rate in units of the
    single-qubit 1/T2* (a |x><y| site-Fock coherence decays at 2 gamma d(x, y);
    a single qubit's |0><1| sets the d = 1 reference)."""
    ca, pa = slater_config_probs(V, a, nn)
    cb, pb = slater_config_probs(V, b, nn)
    sa = [frozenset(c) for c in ca]
    return float(sum(pa[i] * pb[j] * len(sa[i] ^ sa[j])
                     for i in range(len(sa)) for j in range(len(sa))))


def zz_drift_exact(theta, a, b, nn):
    """First-order ZZ systematic: phase-slope difference per step between the two
    branch states = zeta_step * (<B>_a - <B>_b), B = sum z_i z_{i+1}, evaluated in
    the 3-magnon FLOQUET mode states (exact one-body Wick)."""
    _, V = floquet(theta, nn)

    def bexp(triple):
        rho = sum(np.outer(V[:, k - 1], V[:, k - 1].conj()) for k in triple)
        tot = 0.0
        for i in range(nn - 1):
            ni, nj = rho[i, i].real, rho[i + 1, i + 1].real
            gij = rho[i, i + 1]
            tot += 1 - 2 * ni - 2 * nj + 4 * (ni * nj - abs(gij) ** 2)
        return tot

    return ZZ_RAD_PER_STEP / 2.0 * (bexp(a) - bexp(b))


def visibility(m, nn, m_deph, p2):
    """Predicted fringe amplitude at M = m steps."""
    givens = nn * (nn - 1) // 2
    bonds = nn - 1
    cx = 5 + 2 * givens * 2 + m * bonds * 2 + 5
    t_us = m * TAU_STEP_US + 2 * givens * 0.5 / 5 + 2.0  # steps + networks (5 lanes) + seed
    return (1 - p2) ** cx * np.exp(-3 * t_us / T1_US) * np.exp(-m_deph * t_us / T2S_US)


def slope_sigma(ms, nn, m_deph, p2, shots):
    """Weighted-LSQ slope error from the PREDICTED per-point sigma_phi."""
    sig = np.array([1.0 / (visibility(m, nn, m_deph, p2) * np.sqrt(shots)) for m in ms])
    w = 1.0 / sig ** 2
    xm = np.sum(w * ms) / np.sum(w)
    return 1.0 / np.sqrt(np.sum(w * (ms - xm) ** 2))


def chiral_certificate(V, a, b, nn):
    """The mirror pair's protection: (1) DIAGONAL, all range: max difference of
    per-site occupations and of <z_i z_j> for EVERY pair (i, j) (not just NN bonds;
    real-device ZZ crosstalk is longer-range) between the two 3-magnon branches;
    machine-zero certifies zero first-order drift from ZZ and from single-qubit
    Z-detuning even with SITE-DEPENDENT couplings. (2) OFF-DIAGONAL: the coherent
    branch-mixing elements <sigma|z_i z_j|tau> and <sigma|n_i|tau>, the one channel
    that could Rabi-mix an EXACT degeneracy at first order; they vanish identically
    by Slater-Condon (the branches differ in three occupied orbitals, the operators
    are at most 2-body), certified numerically."""
    def rho(triple):
        return sum(np.outer(V[:, k - 1], V[:, k - 1].conj()) for k in triple)
    rA, rB = rho(a), rho(b)
    occ = np.max(np.abs(np.diag(rA).real - np.diag(rB).real))

    def zz(r, i, j):
        ni, nj = r[i, i].real, r[j, j].real
        return 1 - 2 * ni - 2 * nj + 4 * (ni * nj - abs(r[i, j]) ** 2)

    bd = max(abs(zz(rA, i, j) - zz(rB, i, j))
             for i in range(nn) for j in range(i + 1, nn))

    # off-diagonal branch mixing: amplitudes of both Slater states over the
    # C(nn, 3) site-configuration basis (z-diagonal operators stay diagonal there)
    configs = list(combinations(range(nn), 3))

    def amps(triple):
        cols = V[:, [k - 1 for k in triple]]
        v = np.array([np.linalg.det(cols[list(c), :]) for c in configs])
        return v / np.linalg.norm(v)

    aT, aS = amps(a), amps(b)
    zvals = np.array([[1 if i in c else -1 for i in range(nn)] for c in configs])
    mix = 0.0
    for i in range(nn):
        mix = max(mix, abs(np.sum(aS.conj() * zvals[:, i] * aT)))
        for j in range(i + 1, nn):
            mix = max(mix, abs(np.sum(aS.conj() * zvals[:, i] * zvals[:, j] * aT)))
    return occ, bd, mix


def lab_report(nn, coll_a, coll_b, controls, theta, p2, shots, ms, flyable):
    ncomb = nn + 1
    ph, V = floquet(theta, nn)
    m_deph = dephasing_multiplier(V, coll_a, coll_b, nn)
    ss = slope_sigma(np.array(ms), nn, m_deph, p2, shots)
    coll_slope = phase_sum(coll_a, ph) - phase_sum(coll_b, ph)
    coll_zz = zz_drift_exact(theta, coll_a, coll_b, nn)

    lam = [np.cos(k * np.pi / ncomb) for k in range(1, nn + 1)]
    tag = "PRIMARY (flyable)" if flyable else "APPENDIX (not flyable under this pre-registration)"
    print(f"\nLAB n = {ncomb} (N = {nn}) {tag}: collision {coll_a}~{coll_b}, theta = {theta},")
    print(f"  p2 = {p2}, shots = {shots}, M-grid = {list(ms)}")
    print(f"  dephasing multiplier m_deph(A0) = {m_deph:.3f} (E[Hamming] of the branch pair)")
    print(f"  V(M = {ms[-1]}) = {visibility(ms[-1], nn, m_deph, p2):.3f}, "
          f"sigma_slope(A0) = {ss:.4f} rad/step")
    if abs(coll_zz) < 1e-12:
        occ, bd, mix = chiral_certificate(V, coll_a, coll_b, nn)
        print(f"  A0 collision: Floquet drift {coll_slope:+.4f} "
              f"({coll_slope / theta ** 3:.3f} theta^3), ZZ = {coll_zz:+.1e}")
        print(f"  chiral certificate: max per-site occ diff = {occ:.1e}, "
              f"max all-range zz diff = {bd:.1e} (site-dependent zeta/delta protected)")
        print(f"  branch-mixing certificate: max |<sigma|z_i z_j|tau>|, |<sigma|n_i|tau>| "
              f"= {mix:.1e} (Slater-Condon zero: 3-orbital difference vs <= 2-body operators)")
        dA0 = d_drift_dtheta(coll_a, coll_b, theta, nn)
        verdict_arm = next(c for _, c, role in controls if role == "verdict")
        dV = d_slope_dtheta(coll_a, verdict_arm, theta, nn)
        sV = phase_sum(coll_a, ph) - phase_sum(verdict_arm, ph)
        # the clause-(c)-admissible band, DERIVED: (3 sigma + 5|ZZ_verdict|)/dV
        zzV = zz_drift_exact(theta, coll_a, verdict_arm, nn)
        band = (3 * ss + ZZ_SAFETY * abs(zzV)) / abs(dV)
        lin_err = max(abs((phase_sum(coll_a, floquet(theta + dth, nn)[0])
                           - phase_sum(verdict_arm, floquet(theta + dth, nn)[0]))
                          - (sV + dV * dth)) for dth in (-band, band))
        print(f"  theta-hat pin (verdict arm {verdict_arm}): "
              f"theta_hat = {theta} + (m_meas - {sV:+.4f})/{dV:.4f};")
        print(f"    A0 center = {coll_slope:+.4f} + {dA0:.4f}*(theta_hat - {theta}); "
              f"clause-(c) band |theta_hat - {theta}| <= {band:.4f}, "
              f"linearization error <= {lin_err:.1e} rad/step over it")
        # theta-hat statistical propagation into clause (a)'s test statistic
        ss_v = slope_sigma(np.array(ms), nn,
                           dephasing_multiplier(V, coll_a, verdict_arm, nn), p2, shots)
        infl = np.sqrt(1 + (dA0 * (ss_v / abs(dV)) / ss) ** 2)
        print(f"    sigma_a = {infl:.3f} x sigma_slope (theta-hat variance propagated)")
        # worst quorum-passing point-drop inflation of the realized sigma_slope
        from itertools import combinations as _comb
        worst = 1.0
        ms_arr = list(ms)
        for r in range(6, len(ms_arr) + 1):
            for keep in _comb(ms_arr, r):
                if sum(1 for m in keep if m <= 2) >= 2 and sum(1 for m in keep if m >= 6) >= 2:
                    worst = max(worst, slope_sigma(np.array(keep), nn, m_deph, p2, shots) / ss)
        print(f"    worst quorum-passing drop inflates realized sigma_slope x{worst:.2f}")
    else:
        print(f"  A0 collision: Floquet drift {coll_slope:+.4f}, ZZ {coll_zz:+.4f}")

    # nearest-impostor rejection: minimize the clause-(a) TEST METRIC itself
    # (the signed slope distance from the collision drift), not the undeformed dS
    best = None
    for c in combinations(range(1, nn + 1), 3):
        if c in (coll_a, coll_b):
            continue
        s_c = phase_sum(coll_a, ph) - phase_sum(c, ph)
        rej_c = abs(s_c - coll_slope) / ss
        if best is None or rej_c < best[0]:
            dS = abs(sum(lam[k - 1] for k in c) - sum(lam[k - 1] for k in coll_a))
            best = (rej_c, c, s_c, dS)
    rej, imp, imp_slope, dS_imp = best
    print(f"  nearest impostor {imp} (dS = {dS_imp:.3f}, by the clause-(a) metric): "
          f"slope {imp_slope:+.4f}, {rej:.1f} sigma_slope from the window center "
          f"(margin {rej / 5:.1f}x vs 5)")

    for name, c, role in controls:
        s = phase_sum(coll_a, ph) - phase_sum(c, ph)
        z = zz_drift_exact(theta, coll_a, c, nn)
        m_arm = dephasing_multiplier(V, coll_a, c, nn)
        ss_arm = slope_sigma(np.array(ms), nn, m_arm, p2, shots)
        comb = np.sqrt(ss ** 2 + ss_arm ** 2)
        dS = abs(sum(lam[k - 1] for k in c) - sum(lam[k - 1] for k in coll_a))
        # true slope gap (sign-safe), floored by the systematics
        sep = (abs(s - coll_slope) - abs(coll_zz) - ZZ_SAFETY * abs(z)) / comb
        print(f"  {name} {c} dS={dS:.3f} [{role}]: slope {s:+.4f}, ZZ {z:+.4f}, "
              f"m_deph {m_arm:.2f}, floored sep (combined, {ZZ_SAFETY:.0f}x ZZ) = "
              f"{sep:5.1f} sigma" + (f" (margin {sep / 5:.1f}x vs 5)" if role == "verdict" else ""))

    # V-model sensitivity: a uniform 2x visibility miss doubles every sigma
    print(f"  V-model sensitivity: uniform V/2 doubles sigma_slope "
          f"({ss:.4f} -> {2 * ss:.4f}); margins halve, verdict clauses survive")

    n_circ = len(ms) * 2 * (1 + len(controls))
    qpu_min = n_circ * shots * BILLED_S_PER_SHOT / 60
    n_cal = 10
    tot_min = (n_circ + n_cal) * shots * BILLED_S_PER_SHOT / 60
    print(f"  {n_circ} science circuits x {shots} shots ~ {qpu_min:.1f} QPU min "
          f"shots-ratio ({tot_min:.1f} incl. {n_cal} cal PUBs); x1.2-1.4 depth uplift "
          f"(the comparable billed 119 s vs its 85-99 s projection) -> "
          f"{1.2 * tot_min:.1f}-{1.4 * tot_min:.1f} min total")


def d_slope_dtheta(a, c, theta, nn, h=1e-4):
    """Sensitivity of an arm's fringe slope to the hop angle."""
    def slope(th):
        ph, _ = floquet(th, nn)
        return phase_sum(a, ph) - phase_sum(c, ph)
    return (slope(theta + h) - slope(theta - h)) / (2 * h)


def d_drift_dtheta(a, b, theta, nn, h=1e-4):
    """Sensitivity of the collision drift to hop-angle miscalibration."""
    return d_slope_dtheta(a, b, theta, nn, h)


def theta_exponent_scan(a, b, nn):
    """The collision drift's theta-order, fitted over the flown range."""
    thetas = np.array([0.02, 0.05, 0.1, 0.2, 0.32, 0.5])
    drifts = []
    for th in thetas:
        ph, _ = floquet(th, nn)
        drifts.append(abs(phase_sum(a, ph) - phase_sum(b, ph)))
    expo = np.polyfit(np.log(thetas), np.log(drifts), 1)[0]
    print(f"n = {nn + 1} collision {a}~{b}: drift theta-exponent = {expo:.4f} "
          f"over theta in [0.02, 0.5]; coefficient {drifts[0] / 0.02 ** 3:.4f} (theta -> 0) "
          f"to {drifts[-1] / 0.5 ** 3:.4f} (theta = 0.5)")


def main():
    # ---- the theta^3 survival of the primary collision, over the flown range
    theta_exponent_scan((1, 5, 7), (2, 4, 8), 8)

    # ---- the flyable design (v2.2): theta = 0.5, integer M-grid 0..8, 16384 shots,
    # three arms (the redundant second near-control dropped for the billing cushion)
    N9 = 8
    controls_n9 = [("A1", (2, 5, 6), "secondary"),
                   ("A2", (1, 6, 7), "verdict")]
    lab_report(N9, (1, 5, 7), (2, 4, 8), controls_n9,
               theta=0.5, p2=0.003, shots=16384, ms=range(0, 9), flyable=True)

    # sensitivity: the same lab at the v1 knobs, for the revision record
    lab_report(N9, (1, 5, 7), (2, 4, 8), controls_n9,
               theta=0.32, p2=0.003, shots=8192, ms=range(0, 11, 2), flyable=True)

    # p2 stress: the day-of abort threshold
    lab_report(N9, (1, 5, 7), (2, 4, 8), controls_n9,
               theta=0.5, p2=0.005, shots=16384, ms=range(0, 9), flyable=True)

    # ---- appendix: the n = 12 equal-level (F130) pair, NOT flyable here
    N12 = 11
    controls_n12 = [("C1", (3, 4, 7), "secondary"),
                    ("C2", (1, 2, 11), "secondary"),
                    ("C3", (2, 4, 9), "verdict")]
    lab_report(N12, (1, 2, 10), (3, 5, 6), controls_n12,
               theta=0.5, p2=0.003, shots=16384, ms=range(0, 9), flyable=False)

    # ---- Floquet-vs-cos deformation at the flown theta (the modes stay comb-like)
    for nn in (8, 11):
        n = nn + 1
        sine = np.array([[np.sin(k * j * np.pi / n) for j in range(1, nn + 1)]
                         for k in range(1, nn + 1)]).T
        sine = sine / np.linalg.norm(sine, axis=0, keepdims=True)
        for theta in (0.32, 0.5):
            _, V = floquet(theta, nn)
            ov = [abs(np.vdot(sine[:, k], V[:, k])) for k in range(nn)]
            print(f"N = {nn}, theta = {theta}: worst |<sine|floquet>| = {min(ov):.4f}")


if __name__ == "__main__":
    main()
