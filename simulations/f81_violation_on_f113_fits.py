"""f81_violation computed on the F113-fitted Kingston Lindbladians (the bridge, part 1).

Zero-QPU leg of experiments/F81_VIOLATION_HARDWARE_BRIDGE.md: take the only
hardware rho(t) ever fit to a Lindbladian in this repo (f95 Kingston
angle-steering, 2026-05-16; fit pipeline = F113), compute the F81/F82/F84
discriminator f81_violation = ||M_anti - L_{H_odd}||_F on each fitted L,
run the inverse readout, and compare against the INDEPENDENT calibration T1.

The tautology trap (stated up front): the F113 fit is
Z+T1-PARAMETERIZED, so the violation of the fitted L equals the F82 closed
form of the fitted gamma_T1 BY CONSTRUCTION. Section 1's exact match is the
end-to-end method check, not a finding. The physics content is (b): the
violation-implied gamma_T1,RMS vs the independent calibration T1 (expected
OVERSHOOT, because the underfit minimal model dumps all non-T1 noise into T1).

Why a hand-rolled decomposition path: the f95 fit's Hamiltonian is the 1-body
Z-drive H = (omega/2)*sum_l Z_l, which pi_decompose_M's terms interface cannot
express (1-letter terms fall through its 2-body/k-body branches). Section 0
GATES the hand-rolled path against fw.pi_decompose_M on a terms-expressible
case before it is trusted (no Hollywood).

Run from anywhere (chdirs to repo root itself):
    python -X utf8 simulations/f81_violation_on_f113_fits.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
os.chdir(REPO_ROOT)
sys.path.insert(0, str(REPO_ROOT / "simulations"))

import framework as fw  # noqa: E402
import f113_t1_extraction_kingston as f113  # noqa: E402  (reuses load_f95 + fit_z_t1)

from framework.diagnostics.f82_t1_dissipator import (  # noqa: E402
    estimate_T1_from_violation,
    predict_T1_dissipator_violation,
)
from framework.lindblad import (  # noqa: E402
    lindbladian_general,
    lindbladian_z_plus_t1,
    palindrome_residual,
)
from framework.pauli import (  # noqa: E402
    _k_to_indices,
    _vec_to_pauli_basis_transform,
    pauli_string,
    site_op,
)
from framework.symmetry import build_pi_full  # noqa: E402

I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Y2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)
SM = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma-minus, |1> -> |0> (framework convention)
SP = SM.conj().T


def pauli_labels(N):
    """Pauli-string labels in the framework's own basis enumeration order.

    The idx -> letter map is derived empirically from pauli_string (no
    convention assumed)."""
    singles = {"I": I2, "X": X2, "Y": Y2, "Z": Z2}
    pair2letter = {}
    for k in range(4):
        pair = tuple(_k_to_indices(k, 1))[0]  # per-site index = (na, nb) bit pair
        P = pauli_string([pair])
        for name, mat in singles.items():
            if np.allclose(P, mat):
                pair2letter[tuple(pair)] = name
                break
    assert len(pair2letter) == 4, f"could not identify Pauli index map: {pair2letter}"
    return ["".join(pair2letter[tuple(p)] for p in _k_to_indices(k, N)) for k in range(4 ** N)]


def x_all(N):
    out = X2
    for _ in range(N - 1):
        out = np.kron(out, X2)
    return out


def decompose_on_L(L_vec, H, gamma_z_list, N):
    """f81_violation for an EXPLICIT vec-form L and explicit H (any H, incl. 1-body).

    Mirrors pi_decompose_M's internals verbatim (same primitives, same
    normalizations), but takes H as a matrix instead of a terms list, and
    projects H_odd from below via X^N conjugation instead of classifying terms.
    """
    d = 2 ** N
    Sigma_gamma = float(sum(gamma_z_list))
    M = palindrome_residual(L_vec, Sigma_gamma, N)
    Pi = build_pi_full(N)
    PiMPi = Pi @ M @ Pi.conj().T
    M_sym = (M + PiMPi) / 2
    M_anti = (M - PiMPi) / 2

    XN = x_all(N)
    H_odd = (H - XN @ H @ XN) / 2  # Pi^2 = X^N conjugation; odd projector from below
    Id_d = np.eye(d, dtype=complex)
    T = _vec_to_pauli_basis_transform(N)
    L_H_odd_vec = -1j * (np.kron(H_odd, Id_d) - np.kron(Id_d, H_odd.T))
    L_H_odd = (T.conj().T @ L_H_odd_vec @ T) / d

    D = M_anti - L_H_odd
    return {
        "M": M,
        "M_sym": M_sym,
        "M_anti": M_anti,
        "L_H_odd": L_H_odd,
        "D": D,
        "f81_violation": float(np.linalg.norm(D)),
        "norm_sq": {
            "M": float(np.linalg.norm(M) ** 2),
            "M_sym": float(np.linalg.norm(M_sym) ** 2),
            "M_anti": float(np.linalg.norm(M_anti) ** 2),
            "L_H_odd": float(np.linalg.norm(L_H_odd) ** 2),
        },
    }


def build_L_t1_pump(H, gamma_z_list, gamma_t1_list, gamma_pump_list, N):
    """Z-dephasing + cooling (sigma-) + heating (sigma+), vec form, framework builder."""
    c_ops = []
    for l, g in enumerate(gamma_z_list):
        if g != 0:
            c_ops.append(np.sqrt(g) * site_op(N, l, "Z"))
    for l, g in enumerate(gamma_t1_list):
        if g != 0:
            ops = [I2] * N
            ops[l] = SM
            full = ops[0]
            for op in ops[1:]:
                full = np.kron(full, op)
            c_ops.append(np.sqrt(g) * full)
    for l, g in enumerate(gamma_pump_list):
        if g != 0:
            ops = [I2] * N
            ops[l] = SP
            full = ops[0]
            for op in ops[1:]:
                full = np.kron(full, op)
            c_ops.append(np.sqrt(g) * full)
    return lindbladian_general(H, c_ops)


def nonzero_entries(D, N, tol=1e-12):
    labels = pauli_labels(N)
    out = []
    for i in range(D.shape[0]):
        for j in range(D.shape[1]):
            if abs(D[i, j]) > tol:
                out.append((labels[i], labels[j], D[i, j]))
    return out


def h_drive(N, omega):
    return (omega / 2.0) * sum(site_op(N, l, "Z") for l in range(N))


def main():
    np.set_printoptions(linewidth=140, suppress=True)
    print("=" * 96)
    print("S4 bridge scout: f81_violation on the F113-fitted Kingston Lindbladians (zero QPU)")
    print("=" * 96)

    # ------------------------------------------------------------------
    # Section 0: GATE the hand-rolled path against fw.pi_decompose_M
    # ------------------------------------------------------------------
    print()
    print("Section 0: gate. Hand-rolled decomposition vs fw.pi_decompose_M")
    print("  (N=2 chain, H = J(XX+YY), gamma_z=0.1, gamma_t1=0.007; terms-expressible case)")
    chain2 = fw.ChainSystem(N=2, J=1.0)
    ref = fw.pi_decompose_M(chain2, [("X", "X"), ("Y", "Y")], gamma_z=0.1, gamma_t1=0.007)
    H_gate = np.kron(X2, X2) + np.kron(Y2, Y2)  # the single (0,1) bond at J=1
    L_gate = lindbladian_z_plus_t1(H_gate, [0.1, 0.1], [0.007, 0.007])
    mine = decompose_on_L(L_gate, H_gate, [0.1, 0.1], 2)
    dv = abs(ref["f81_violation"] - mine["f81_violation"])
    dM = abs(ref["norm_sq"]["M"] - mine["norm_sq"]["M"])
    dA = abs(ref["norm_sq"]["M_anti"] - mine["norm_sq"]["M_anti"])
    print(f"  f81_violation: ref = {ref['f81_violation']:.15e}")
    print(f"                 mine = {mine['f81_violation']:.15e}   |diff| = {dv:.3e}")
    print(f"  ||M||^2 diff = {dM:.3e},  ||M_anti||^2 diff = {dA:.3e}")
    gate_ok = dv < 1e-12 and dM < 1e-10 and dA < 1e-10
    print(f"  GATE: {'PASS' if gate_ok else 'FAIL'}")
    if not gate_ok:
        print("  ABORT: hand-rolled path does not reproduce pi_decompose_M; nothing below is trustworthy.")
        sys.exit(1)

    # F82 closed form on the gate case (uniform rates, N=2)
    ns2 = SimpleNamespace(N=2)
    v82 = predict_T1_dissipator_violation(ns2, 0.007)
    print(f"  F82 closed form check: predicted {v82:.15e}, computed {mine['f81_violation']:.15e}, "
          f"|diff| = {abs(v82 - mine['f81_violation']):.3e}")

    # ------------------------------------------------------------------
    # Section 1: the four f95 pair-runs
    # ------------------------------------------------------------------
    print()
    print("=" * 96)
    print("Section 1: hardware-fit -> violation -> inverse readout, four f95 pair-runs")
    print("=" * 96)
    print()
    print("  Fit = F113 minimal Z+T1 model (H = (omega/2)*sum Z_l known, gamma_z + gamma_T1 free),")
    print("  refit here with the identical f113 pipeline. Violation computed on the FITTED L")
    print("  (drive H included; H is 1-body Z, entirely Pi^2-odd, subtracted exactly via L_H_odd).")
    print()

    # Documented F113 values (experiments/F113_T1_EXTRACTION_KINGSTON.md) for drift check
    doc_vals = {
        "omega=0.13 A_mid q82-q83": (0.0980, 0.00574),
        "omega=0.13 B_high q13-q14": (0.1402, 0.00616),
        "omega=0.25 A_mid q82-q83": (0.3361, 0.00722),
        "omega=0.25 B_high q13-q14": (3.5237, 0.00564),
    }

    rows = []
    for omega, t_us, rhos, label, T1_calib in f113.load_f95():
        (gz, gt1), loss = f113.fit_z_t1(omega, t_us, rhos)
        rms = np.sqrt(loss / max(len(t_us) - 1, 1))

        H = h_drive(2, omega)
        L_fit = lindbladian_z_plus_t1(H, [gz, gz], [gt1, gt1])
        dec = decompose_on_L(L_fit, H, [gz, gz], 2)
        viol = dec["f81_violation"]

        v_pred = predict_T1_dissipator_violation(ns2, gt1)     # F82 closed form of the fit
        g_rms = estimate_T1_from_violation(ns2, viol)          # inverse readout
        T1_implied = 1.0 / g_rms if g_rms > 0 else float("inf")

        g_cal = [1.0 / t for t in T1_calib]
        g_cal_mean = 0.5 * (g_cal[0] + g_cal[1])
        g_cal_rms = float(np.sqrt(0.5 * (g_cal[0] ** 2 + g_cal[1] ** 2)))

        print(f"  {label}")
        doc = doc_vals.get(label)
        drift = ""
        if doc:
            rel = max(abs(gz - doc[0]) / max(doc[0], 1e-12), abs(gt1 - doc[1]) / max(doc[1], 1e-12))
            drift = f"   (doc: gz={doc[0]}, gt1={doc[1]}; max rel drift {rel:.2%})"
        print(f"    refit: gamma_z = {gz:.5f}, gamma_T1 = {gt1:.5f} per us, RMS = {rms:.4f}{drift}")
        print(f"    f81_violation(fitted L)      = {viol:.10e}")
        print(f"    F82 closed form of fit gt1   = {v_pred:.10e}   ratio = {viol / v_pred:.12f}  <- tautology check")
        print(f"    inverse readout gamma_T1,RMS = {g_rms:.6f} per us  ->  T1_implied = {T1_implied:.1f} us")
        print(f"    INDEPENDENT calibration      : T1 = {T1_calib[0]:.1f} / {T1_calib[1]:.1f} us"
              f"  (rates {g_cal[0]:.5f} / {g_cal[1]:.5f}; mean {g_cal_mean:.5f}, rms {g_cal_rms:.5f})")
        print(f"    overshoot: g_rms / g_cal_mean = {g_rms / g_cal_mean:.3f},"
              f"  g_rms / g_cal_rms = {g_rms / g_cal_rms:.3f}   (F113 underfit predicts > 1)")
        print()
        rows.append((label, gz, gt1, viol, v_pred, g_rms, T1_implied, g_cal_mean,
                     g_rms / g_cal_mean, rms))

    print(f"  {'pair-run':<28} {'gz_fit':>8} {'gt1_fit':>9} {'violation':>12} {'gT1,RMS':>9} "
          f"{'T1_impl':>8} {'gT1,cal':>9} {'over':>6} {'fitRMS':>7}")
    print("  " + "-" * 104)
    for (label, gz, gt1, viol, v_pred, g_rms, T1i, g_cal_mean, over, rms) in rows:
        print(f"  {label:<28} {gz:>8.4f} {gt1:>9.5f} {viol:>12.6e} {g_rms:>9.5f} "
              f"{T1i:>8.1f} {g_cal_mean:>9.5f} {over:>6.3f} {rms:>7.4f}")

    # ------------------------------------------------------------------
    # Section 2: the single-entry lead, from below
    # ------------------------------------------------------------------
    print()
    print("=" * 96)
    print("Section 2: which L-entries does M_anti - L_H_odd actually touch? (single-entry lead)")
    print("=" * 96)

    print()
    print("  (a) N=1, gamma_z=0.1, gamma_T1=0.007, H = (0.13/2) Z  (drive present):")
    H1 = h_drive(1, 0.13)
    L1 = build_L_t1_pump(H1, [0.1], [0.007], [0.0], 1)
    d1 = decompose_on_L(L1, H1, [0.1], 1)
    for (ro, co, val) in nonzero_entries(d1["D"], 1):
        print(f"      D[{ro} <- {co}] = {val.real:+.6f}{val.imag:+.6f}j")
    print(f"      ||D|| = {d1['f81_violation']:.10f}   (F82 at N=1 predicts gamma_T1 = 0.007)")

    print()
    print("  (a') same, gamma_z = 0.4 (gamma_z-blindness on the real pipeline):")
    L1b = build_L_t1_pump(H1, [0.4], [0.007], [0.0], 1)
    d1b = decompose_on_L(L1b, H1, [0.4], 1)
    print(f"      ||D|| = {d1b['f81_violation']:.10f}   (unchanged)")

    print()
    print("  (b) N=1 thermal (F84): cooling 0.007, heating 0.003 -> net 0.004:")
    L1t = build_L_t1_pump(H1, [0.1], [0.007], [0.003], 1)
    d1t = decompose_on_L(L1t, H1, [0.1], 1)
    for (ro, co, val) in nonzero_entries(d1t["D"], 1):
        print(f"      D[{ro} <- {co}] = {val.real:+.6f}{val.imag:+.6f}j")
    print(f"      ||D|| = {d1t['f81_violation']:.10f}   (F84 predicts net cooling 0.004)")
    L1e = build_L_t1_pump(H1, [0.1], [0.005], [0.005], 1)
    d1e = decompose_on_L(L1e, H1, [0.1], 1)
    print(f"      detailed balance 0.005/0.005: ||D|| = {d1e['f81_violation']:.3e}   (F84 predicts 0)")

    print()
    print("  (c) N=2 fitted L, first pair-run: full support of D:")
    omega, t_us, rhos, label, T1_calib = f113.load_f95()[0]
    (gz, gt1), _ = f113.fit_z_t1(omega, t_us, rhos)
    H = h_drive(2, omega)
    L_fit = lindbladian_z_plus_t1(H, [gz, gz], [gt1, gt1])
    dec = decompose_on_L(L_fit, H, [gz, gz], 2)
    ents = nonzero_entries(dec["D"], 2, tol=1e-10)
    for (ro, co, val) in ents:
        print(f"      D[{ro} <- {co}] = {val.real:+.6f}{val.imag:+.6f}j")
    print(f"      {len(ents)} entries (proof Step 4: 2 sites x 4^(N-1) = 8), each magnitude gamma_T1 = {gt1:.5f}")
    mags = [abs(v) for (_, _, v) in ents]
    print(f"      magnitudes: min {min(mags):.6f}, max {max(mags):.6f}")
    proxy = (2 ** (2 - 1)) * np.sqrt(2 * gt1 ** 2)
    print(f"      affine-drift proxy 2^(N-1)*sqrt(sum_l a_l^2) with a_l = gt1: {proxy:.10e}")
    print(f"      computed violation:                                          {dec['f81_violation']:.10e}")

    # ------------------------------------------------------------------
    # Section 3: verdict
    # ------------------------------------------------------------------
    print()
    print("=" * 96)
    print("Section 3: verdict (deliverable c)")
    print("=" * 96)
    print("""
  TAUTOLOGICAL (by construction, and confirmed above): on a Z+T1-parameterized
  fit, f81_violation(fitted L) = F82 closed form of the fitted gamma_T1, ratio
  1.000000000000. The pipeline demo is end-to-end sound; the number carries no
  information beyond gamma_T1_fit itself.

  THE PHYSICS THAT SURVIVES THE TAUTOLOGY: the inverse readout vs the
  INDEPENDENT calibration T1 (Section 1 'overshoot' column). The F113 underfit
  dumps all non-T1 noise into the sigma- channel, so the violation-implied
  gamma_T1,RMS overshoots calibration 1/T1. The violation is therefore an
  EFFECTIVE non-unital-content reading, an upper bound on true T1 rate content.

  NON-TAUTOLOGICAL VERSIONS (what the bridge needs), two routes:

  R1  FREE-FORM FIT: fit an unconstrained Lindblad L (or at least a model with
      more channels than sigma-) to idle tomography data, THEN compute the
      violation. The violation then measures the Pi^2-odd content the fit
      FOUND, not the channel we put in. Needs: idle (Hamiltonian-off) delay
      grid, full per-delay tomography, >= tens of time points (the f95 6-point
      trajectory is too thin; F113 lesson).

  R2  SINGLE-ENTRY PROXY (Section 2, confirmed from below): the entire odd
      content of the LOCAL amplitude-damping family lives in the per-site
      (Z_l <- I) transfer entries, value a_l = gamma_down,l - gamma_up,l =
      the AFFINE DRIFT of <Z_l>. On the Bloch equation dz/dt = a_l - b_l z,
      a_l = z_inf,l * b_l = z_inf,l / T1_l: BOTH come out of the STANDARD
      T1 leg (decay rate + asymptote). So

          f81_violation(local family) = 2^(N-1) * sqrt(sum_l a_l^2),

      measurable from N standard T1-relaxometry experiments WITH ASYMPTOTE
      FITTING, no tomography at all. Caveat (honesty for the adapter): this
      reads the local sigma+- family projection ONLY; it is a LOWER bound on
      the full violation. Correlated / non-local Pi^2-odd channels live on
      other entries and need R1. The run_price_pair Block B (|1> free
      asymptote) already measures exactly the R2 ingredients.
""")


if __name__ == "__main__":
    main()
