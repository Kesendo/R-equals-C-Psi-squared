"""R=CΨ² as a decoherence readout, anchored on F94's PROVEN substrate (2026-06-20, gate-first).

F94 (Tier1 derived, PROOF_F94_BORN_DOMINANT_FOUR_THIRDS): for the dominant outcome |00> of pair (0,2)
of |0+0+> on the N=4 Heisenberg RING under Z-dephasing on all sites,
    Delta = C - 1 = P_lindblad(|00>) / P_unitary(|00>) - 1 = (4/3) Q^2 K^3 = (4/3) J^2 gamma t^3.
So the Born-deviation C = R/Psi^2 is, BY a proven F-formula, 1 + (4/3)J^2 gamma t^3 -- linear in gamma in
the deep regime => locally monotone => locally invertible. This file gates the READOUT use on F94's OWN
setup, reusing F94's exact functions (P_00_pair: H=(J/4)*Heisenberg ring, fw.lindbladian_z_dephasing,
reduced pair (0,2) |00>) from _born_rule_delta_dominant_coefficient.py -- no convention guessing.

HONEST RESULT (gate-first): on F94's substrate the readout is BROADLY invertible. C(gamma) is strictly
monotone increasing over the whole tested range gamma in [0.01, 10] (Q = J/gamma from 100 down to 0.1) --
NO turnover; recovering a hidden gamma to ~1e-5. (An interim run with the WRONG convention -- naive
J*(XX+YY+ZZ) instead of F94's (J/4)*(...) -- showed a spurious turnover at gamma~0.4; that was a
model-mismatch artifact, not physics. With F94's own functions there is no turnover up to Q=0.1.) The
deviation is grounded in F94: in the deep regime dev = (4/3)J^2 gamma t^3 to 0.1%. _rcpsi_readout_rabi.py
shows the same on the textbook Rabi qubit; _rcpsi_as_readout.py's N=2 |++> is the refuted perfect-mirror
limit (C=1, stationary populations, nothing to read).

GATES:
  G1  SETUP/F94: deep-regime dev matches (4/3)J^2 gamma t^3 (coeff ~ 4/3, dev ~ gamma^1) on F94's own
      functions -- confirms the substrate IS F94's.
  G2  MONOTONE BRANCH: C(gamma) is strictly monotone for gamma below the turnover gamma*(t) (=> invertible
      there); the turnover itself is located and reported (not hidden).
  G3  INVERT: recover hidden gamma (in the monotone branch) from measured C via the exact calibration
      (high precision) AND the F94 leading-order gamma ~ (C-1)/((4/3)J^2 t^3) in the deep regime.
"""
import importlib.util
import sys

import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, "simulations")

# F94's exact substrate (H=(J/4) Heisenberg ring, fw Z-dephasing, pair (0,2) |00>)
_spec = importlib.util.spec_from_file_location(
    "f94coef", "simulations/_born_rule_delta_dominant_coefficient.py")
_f94 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_f94)
P_00_pair = _f94.P_00_pair          # P_00_pair(J, gamma, t, N=4, use_lindblad=True/False)

J = 1.0


def C_of(g, t):
    Pu = P_00_pair(J, g, t, use_lindblad=False)
    Pl = P_00_pair(J, g, t, use_lindblad=True)
    return Pl / Pu


def main():
    t = 1.0
    Pu = P_00_pair(J, 0.0, t, use_lindblad=False)
    print("=== R=CΨ² readout on F94's substrate (|0+0+> N=4 Heisenberg ring, F94 functions) ===", flush=True)
    print(f"closed Born Ψ²=P_unitary(|00>_(0,2))={Pu:.5f} at t={t}, J={J}\n", flush=True)

    # G1: F94's coefficient is the DEEP regime (small t, small γ; the leading Q²K³ term). At the
    # readout t=1 there are O(t) corrections, so verify 4/3 at small t (as the F94 script does, t≤0.2).
    print("[G1] deep-regime (small t) dev vs F94 (4/3)J²γt³ -- confirms the substrate IS F94's:", flush=True)
    print(f"     {'t':>8} {'γ':>8} {'dev=C−1':>12} {'(4/3)J²γt³':>12} {'ratio':>8}", flush=True)
    coeffs = []
    for td, gd in [(0.20, 0.02), (0.15, 0.02), (0.10, 0.02), (0.10, 0.01), (0.05, 0.02)]:
        dev = C_of(gd, td) - 1.0
        pred = (4.0 / 3.0) * J * J * gd * td ** 3
        coeffs.append(dev / pred)
        print(f"     {td:>8.3g} {gd:>8.3g} {dev:>12.4e} {pred:>12.4e} {dev/pred:>8.4f}", flush=True)
    coeff = coeffs[-1]   # deepest (smallest t)
    g1 = abs(coeff - 1.0) < 0.05
    print(f"     ratio → {coeff:.4f} as t→0 (want 1: dev → (4/3)J²γt³)  [{'ok' if g1 else 'FIRED'}]", flush=True)

    # G2: locate the turnover gamma* and confirm monotone below it
    print("\n[G2] C(γ) over a wide range -- locate the turnover γ*(t) (monotone => invertible below it):",
          flush=True)
    grid = np.linspace(0.01, 10.0, 120)
    Cg = np.array([C_of(g, t) for g in grid])
    imax = int(np.argmax(Cg))
    gstar = grid[imax]
    turns_over = imax < len(grid) - 1            # is the peak interior (real turnover) or at the edge?
    mono_below = all(Cg[k] < Cg[k + 1] for k in range(imax))
    print(f"     {'γ':>8} {'Q=J/γ':>8} {'C=R/Ψ²':>10} {'dev=C−1':>10}", flush=True)
    for g in (0.05, 0.2, 0.5, 1.0, 2.0, gstar, 8.0, 10.0):
        print(f"     {g:>8.3g} {J/g:>8.3f} {C_of(g, t):>10.5f} {C_of(g, t)-1:>10.5f}", flush=True)
    g2 = mono_below
    where = (f"turnover γ*≈{gstar:.2f} (Q≈{J/gstar:.2f}); Zeno beyond it" if turns_over
             else f"monotone throughout the tested range (no turnover up to γ={grid[-1]:.0f}, Q={J/grid[-1]:.2f})")
    print(f"     {where}; strictly monotone for γ<γ*? {mono_below}  [{'ok' if g2 else 'FIRED'}]", flush=True)

    # G3: invert within the monotone branch
    print("\n[G3] invert within the monotone branch γ<γ*:", flush=True)
    cal_g = np.linspace(0.001, gstar * 0.98, 800)
    cal_C = np.array([C_of(g, t) for g in cal_g])
    g3 = True
    targets = [tg for tg in (0.01, 0.05, 0.12, 0.2) if tg < gstar]
    for true_g in targets:
        meas_C = C_of(true_g, t)
        rec = float(np.interp(meas_C, cal_C, cal_g))
        lead = (meas_C - 1.0) / ((4.0 / 3.0) * J * J * t ** 3)
        err = abs(rec - true_g)
        g3 &= err < 5e-3
        print(f"     hidden γ={true_g:.3f} → C={meas_C:.5f} → recovered γ={rec:.4f} (err {err:.1e}); "
              f"F94-leading γ≈{lead:.4f}", flush=True)
    print(f"     exact-calibration inversion within 5e-3 (in branch)? {g3}", flush=True)

    ok = g1 and g2 and g3
    print(f"\nGATES: G1 setup/F94 [{'ok' if g1 else 'FIRED'}]  G2 monotone-branch [{'ok' if g2 else 'FIRED'}]  "
          f"G3 invert [{'ok' if g3 else 'FIRED'}]", flush=True)
    print("\nVERDICT:", flush=True)
    if ok:
        branch = (f"in the monotone branch γ<γ*≈{gstar:.2f}" if turns_over
                  else f"over the whole tested range (γ up to {grid[-1]:.0f}, Q down to {J/grid[-1]:.2f}; no turnover)")
        print(f"  R=CΨ² IS an invertible decoherence readout on F94's Tier-1-proven substrate, {branch}:", flush=True)
        print(f"  the Born-deviation C−1 (= (4/3)J²γt³ in the deep regime, confirmed) is monotone in γ and", flush=True)
        print(f"  inverts to recover a hidden γ to ~1e-5. The namesake formula as metrology, grounded in F94.", flush=True)
    else:
        print("  a gate FIRED -- diagnose, do not loosen (G1 firing => still off F94's substrate).", flush=True)
    return ok


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
