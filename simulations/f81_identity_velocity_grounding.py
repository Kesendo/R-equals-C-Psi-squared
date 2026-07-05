"""Grounding check (from below): what IS f81_violation physically? (the bridge, part 2)

Companion of experiments/F81_VIOLATION_HARDWARE_BRIDGE.md.

Candidate object (to be CONFIRMED here, not asserted): the per-site velocity
d<Z_l>/dt of the MAXIMALLY MIXED state. I/2^N is the fixed point of every
unital channel (all Pauli watching); non-unital (sigma+-) content gives the
identity a velocity toward the cold pole. Claim to confirm:

  (1) v_l := d<Z_l>/dt |_{rho = I/2^N} = gamma_down,l - gamma_up,l
      (computed by STATE PROPAGATION, independent of the Pauli-basis entry route)
  (2) f81_violation = 2^(N-1) * sqrt(sum_l v_l^2)   on the fitted Kingston L
  (3) v is temperature-independent (F84 vacuum reading): gamma_down = g0 + c,
      gamma_up = c gives v = g0 for every c
  (4) the measurement recipe a = z_inf / T1 recovers v from the <Z>(t)
      trajectory of the EXACT generator (the R2 leg validation)

Run: python -X utf8 simulations/f81_identity_velocity_grounding.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import expm

REPO_ROOT = Path(__file__).resolve().parents[1]
os.chdir(REPO_ROOT)
sys.path.insert(0, str(REPO_ROOT / "simulations"))

import f113_t1_extraction_kingston as f113  # noqa: E402
from framework.lindblad import lindbladian_general, lindbladian_z_plus_t1  # noqa: E402
from framework.pauli import site_op  # noqa: E402

I2 = np.eye(2, dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)
SP = SM.conj().T

# The framework's vec convention, read off lindbladian_general's kron structure:
# L @ vec(rho) with rho flattened the same way palindrome_residual's transform
# flattens Pauli strings ('F'). We PIN the convention empirically below (gate),
# rather than trusting either docstring.


def vec_F(rho):
    return rho.flatten("F")


def devec_F(v, d):
    return v.reshape((d, d), order="F")


def vec_C(rho):
    return rho.flatten("C")


def devec_C(v, d):
    return v.reshape((d, d), order="C")


def zdot_at_state(L_vec, rho, N, vec, devec):
    """d<Z_l>/dt at state rho: propagate rho through L once (no entry inspection)."""
    d = 2 ** N
    rdot = devec(L_vec @ vec(rho), d)
    return [float(np.real(np.trace(site_op(N, l, "Z") @ rdot))) for l in range(N)]


def build_L(H, gz_list, gdown_list, gup_list, N):
    c_ops = []
    for l, g in enumerate(gz_list):
        if g != 0:
            c_ops.append(np.sqrt(g) * site_op(N, l, "Z"))
    for mat, glist in ((SM, gdown_list), (SP, gup_list)):
        for l, g in enumerate(glist):
            if g != 0:
                ops = [I2] * N
                ops[l] = mat
                full = ops[0]
                for op in ops[1:]:
                    full = np.kron(full, op)
                c_ops.append(np.sqrt(g) * full)
    return lindbladian_general(H, c_ops)


def main():
    print("=" * 88)
    print("S4 grounding: is f81_violation the identity-state polarization velocity?")
    print("=" * 88)

    # Gate: pin the vec convention empirically. Known truth: pure sigma- at
    # rate g on N=1 from |1><1| must give d<Z>/dt = +2g (z=-1, dz/dt=g(1-z)).
    g = 0.01
    L1 = build_L(np.zeros((2, 2), dtype=complex), [0.0], [g], [0.0], 1)
    rho_e = np.array([[0, 0], [0, 1]], dtype=complex)  # |1><1|
    zF = zdot_at_state(L1, rho_e, 1, vec_F, devec_F)[0]
    zC = zdot_at_state(L1, rho_e, 1, vec_C, devec_C)[0]
    print()
    print(f"  vec-convention gate (truth: d<Z>/dt = +2g = {2*g}):  'F': {zF:+.6f}   'C': {zC:+.6f}")
    if abs(zF - 2 * g) < 1e-12:
        vec, devec, conv = vec_F, devec_F, "F"
    elif abs(zC - 2 * g) < 1e-12:
        vec, devec, conv = vec_C, devec_C, "C"
    else:
        print("  ABORT: neither convention reproduces the known velocity.")
        sys.exit(1)
    print(f"  -> framework builders act on '{conv}'-flattened rho. GATE PASS")

    # (1) + (2): the fitted Kingston L, identity-state velocity by propagation
    print()
    print("  (1)+(2) fitted Kingston L (all four pair-runs), rho = I/4:")
    print(f"      {'pair-run':<28} {'v_0':>10} {'v_1':>10} {'gt1_fit':>10} "
          f"{'2^(N-1)*rms->':>14} {'f81_violation':>14}")
    for omega, t_us, rhos, label, T1_calib in f113.load_f95():
        (gz, gt1), _ = f113.fit_z_t1(omega, t_us, rhos)
        H = (omega / 2.0) * sum(site_op(2, l, "Z") for l in range(2))
        L_fit = lindbladian_z_plus_t1(H, [gz, gz], [gt1, gt1])
        rho_id = np.eye(4, dtype=complex) / 4.0
        v = zdot_at_state(L_fit, rho_id, 2, vec, devec)
        pack = (2 ** (2 - 1)) * np.sqrt(sum(x * x for x in v))
        # independent reference from the scout: F82 closed form of gt1
        ref = np.sqrt(2 * gt1 * gt1) * 2
        print(f"      {label:<28} {v[0]:>10.6f} {v[1]:>10.6f} {gt1:>10.6f} "
              f"{pack:>14.8e} {ref:>14.8e}")

    # (3) temperature-independence at the identity (F84 vacuum reading)
    print()
    print("  (3) N=1 thermal ladder, g_down = 0.007 + c, g_up = c (same vacuum part 0.007):")
    H1 = 0.065 * Z2  # the omega=0.13 drive, present to keep it honest
    for c in [0.0, 0.003, 0.01, 0.05]:
        L1t = build_L(H1, [0.1], [0.007 + c], [c], 1)
        v = zdot_at_state(L1t, np.eye(2, dtype=complex) / 2, 1, vec, devec)[0]
        print(f"      c = {c:<6}: v = {v:+.10f}   (predicts +0.007 for every c)")

    # (4) the R2 measurement recipe on the exact generator: a = z_inf / T1
    print()
    print("  (4) R2 recipe check, N=1 exact <Z>(t) trajectory (g_down=0.007, g_up=0.003, gz=0.1):")
    gdn, gup = 0.007, 0.003
    L1t = build_L(H1, [0.1], [gdn], [gup], 1)
    rho0 = np.array([[0, 0], [0, 1]], dtype=complex)  # |1><1|, the T1 leg
    ts = np.linspace(0.0, 1200.0, 400)
    zs = []
    for t in ts:
        rt = devec(expm(L1t * t) @ vec(rho0), 2)
        zs.append(float(np.real(np.trace(Z2 @ rt))))
    zs = np.array(zs)
    # fit z(t) = z_inf + dz * exp(-b t), the same nonlinear form the hardware
    # analysis uses (a crude log-linear estimator biases b by ~5%)
    from scipy.optimize import curve_fit

    def bloch(t, z_inf, dz, b):
        return z_inf + dz * np.exp(-b * t)

    (z_inf, dz_fit, b_fit), _ = curve_fit(bloch, ts, zs, p0=[0.5, -1.5, 0.008], maxfev=20000)
    a_recipe = z_inf * b_fit
    print(f"      z_inf (asymptote)   = {z_inf:+.8f}   (theory (gdn-gup)/(gdn+gup) = {(gdn-gup)/(gdn+gup):+.8f})")
    print(f"      b (T1-leg rate)     = {b_fit:.8f}   (theory gdn+gup = {gdn+gup:.8f})")
    print(f"      a = z_inf * b       = {a_recipe:+.8f}   (theory net cooling gdn-gup = {gdn-gup:+.8f})")
    print(f"      -> a recovers the identity-state velocity: |a - v| = "
          f"{abs(a_recipe - (gdn - gup)):.2e}")


if __name__ == "__main__":
    main()
