"""Numerical Painter sector classification and axis-probe tail fits.

The selected N = 4 XX+YY ring has a transverse Y field and all-site local-Z
dephasing. Local Z is an F1-instance channel; a selected bond jump would be
outside F1's jump premise. Neither channel is assigned here to a Holstein or
material bath.

The script diagonalizes the Liouvillian numerically, decomposes each reduced
single-site eigenmode into {I, X, Y, Z}, and classifies Y/non-Y weights using
a 1e-8 ratio tolerance. A Y-only classification refers to the Pauli Y matrix,
which is purely imaginary and antisymmetric; it does not say a general complex
reduced eigenmode is real antisymmetric.

It also propagates selected X/Y/Z axis probes and reports their specified
single-exponential tail fits. These are selected-model probe-decay observations,
not T1/T2, FID, TROSY, EXSY, or material predictions. Relation of the tail fits
to the slow-mode rates requires a convergence check.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))


PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_op(letters):
    op = PAULI[letters[0]]
    for L in letters[1:]:
        op = np.kron(op, PAULI[L])
    return op


def site_op(N, site, letter):
    letters = ["I"] * N
    letters[site] = letter
    return pauli_op(letters)


def two_site_op(N, a, b, la, lb):
    letters = ["I"] * N
    letters[a] = la
    letters[b] = lb
    return pauli_op(letters)


def hueckel_ring_H(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a in range(N):
        b = (a + 1) % N
        H = H + two_site_op(N, a, b, "X", "X") + two_site_op(N, a, b, "Y", "Y")
    return H


def zeeman_y_total(N):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for l in range(N):
        H = H + site_op(N, l, "Y")
    return H


def commutator_superop_vec(H):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T))


def lindbladian_vec(H, c_list, gammas):
    L = commutator_superop_vec(H)
    for c, g in zip(c_list, gammas):
        L = L + g * dissipator_superop_vec(c)
    return L


def unvec_column(v, d):
    return v.reshape(d, d, order="F")


def vec_column(M):
    return M.reshape(-1, order="F")


def partial_trace_simple(rho, N, keep_site):
    rho_tensor = rho.reshape([2] * N + [2] * N)
    result = np.zeros((2, 2), dtype=complex)
    for i_kept in range(2):
        for j_kept in range(2):
            total = complex(0.0)
            for indices in np.ndindex(*([2] * (N - 1))):
                full_i = list(indices[:keep_site]) + [i_kept] + list(indices[keep_site:])
                full_j = list(indices[:keep_site]) + [j_kept] + list(indices[keep_site:])
                total += rho_tensor[tuple(full_i + full_j)]
            result[i_kept, j_kept] = total
    return result


def per_site_pauli_decomp(rho_2x2):
    """Decompose a 2×2 complex matrix in the Pauli basis {I, X, Y, Z}.
    Returns dict letter -> complex coefficient α such that ρ = Σ α_l · σ_l."""
    coeffs = {}
    for L, sigma in PAULI.items():
        coeffs[L] = np.trace(rho_2x2 @ sigma) / 2.0
    return coeffs


def mode_flavor_signature(rho_full, N):
    """For a full d × d operator ρ, partial-trace to each site, decompose into
    {I, X, Y, Z} per site, and return Y vs non-Y total weight (sum |α|² over
    sites). A pure-Re-flavor mode has Y-weight = 0; pure-Im-flavor has
    XZ-weight = 0 (I is in both columns; we exclude it from non-Y to keep
    the diagnostic sharp)."""
    y_w = 0.0
    xz_w = 0.0
    i_w = 0.0
    for s in range(N):
        rho_s = partial_trace_simple(rho_full, N, s)
        c = per_site_pauli_decomp(rho_s)
        i_w += abs(c["I"]) ** 2
        xz_w += abs(c["X"]) ** 2 + abs(c["Z"]) ** 2
        y_w += abs(c["Y"]) ** 2
    return i_w, xz_w, y_w


def fit_exp_tail(t, signal, frac_tail=0.5):
    """Single-exponential fit to the late-time tail of |signal(t)|.
    Returns (rate, prefactor) so that |signal(t)| ≈ prefactor · exp(-rate · t)
    on the tail. Skips zero / near-zero values."""
    nt = len(t)
    t_tail = t[int(nt * (1 - frac_tail)):]
    s_tail = np.abs(signal[int(nt * (1 - frac_tail)):])
    mask = s_tail > 1e-12
    if mask.sum() < 5:
        return float("nan"), float("nan")
    log_s = np.log(s_tail[mask])
    A = np.vstack([t_tail[mask], np.ones(mask.sum())]).T
    slope, intercept = np.linalg.lstsq(A, log_s, rcond=None)[0]
    return -slope, np.exp(intercept)


def run(N=4, h_zeeman=0.5, gamma=1.0, t_max=20.0, n_steps=400, n_slowest=8):
    print()
    print("=" * 96)
    print(f"N = {N} selected XX+YY ring, h_y = {h_zeeman}, γ_local-Z = {gamma}")
    print("=" * 96)
    print()

    H = hueckel_ring_H(N) + h_zeeman * zeeman_y_total(N)
    c_local_z = [site_op(N, l, "Z") for l in range(N)]
    L = lindbladian_vec(H, c_local_z, [gamma] * N)
    d = 2**N

    print("Q1: numerical Y/non-Y classification of slow modes (1e-8 tolerance)")
    print("-" * 96)
    eigvals, eigvecs = np.linalg.eig(L)
    order = np.argsort(eigvals.real)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    print(f"{'k':>3}  {'Re(λ)':>10}  {'Im(λ)':>10}  {'‖I‖²':>10}  {'‖X,Z‖²':>10}  {'‖Y‖²':>10}  flavor")
    flavor_table = []
    for k in range(n_slowest):
        v = eigvecs[:, k]
        rho = unvec_column(v, d)
        i_w, xz_w, y_w = mode_flavor_signature(rho, N)
        total = xz_w + y_w
        if total < 1e-10:
            flavor = "(trivial)"
        elif y_w / total < 1e-8:
            flavor = "Re (non-Y)"
        elif xz_w / total < 1e-8:
            flavor = "Im (Y-only)"
        else:
            mix = y_w / total
            flavor = f"mixed (Y-frac={mix:.4f})"
        flavor_table.append((k, eigvals[k].real, eigvals[k].imag, flavor, i_w, xz_w, y_w))
        print(f"{k:>3}  {eigvals[k].real:>+10.6f}  {eigvals[k].imag:>+10.6f}  "
              f"{i_w:>10.6f}  {xz_w:>10.6f}  {y_w:>10.6f}  {flavor}")
    print()
    print("Re-flavor = numerically non-Y by the Pauli-weight diagnostic")
    print("Im-flavor = numerically Y-only; σY is purely imaginary and antisymmetric")
    print()

    re_modes = [(k, dec) for (k, dec, om, fl, _, _, _) in flavor_table if fl.startswith("Re")]
    im_modes = [(k, dec) for (k, dec, om, fl, _, _, _) in flavor_table if fl.startswith("Im")]
    print(f"Slowest Re-flavor mode: k={re_modes[0][0]}, Re(λ) = {re_modes[0][1]:+.6f}")
    print(f"Slowest Im-flavor mode: k={im_modes[0][0]}, Re(λ) = {im_modes[0][1]:+.6f}")
    ratio_pred = abs(re_modes[0][1]) / abs(im_modes[0][1])
    print(f"Slow-mode Im/Re rate ratio = |Re(λ_im)| / |Re(λ_re)| = {1.0/ratio_pred:.6f}")
    print()

    # ---- Q2: axis-probe tail fits from full propagation ----
    print("Q2: selected-model axis-probe tail fits")
    print("-" * 96)

    # selected total Pauli-coordinate readouts along each axis
    Mx = sum(site_op(N, l, "X") for l in range(N))
    My = sum(site_op(N, l, "Y") for l in range(N))
    Mz = sum(site_op(N, l, "Z") for l in range(N))

    # initial states: maximally mixed plus a tiny single-site Pauli probe
    # along the desired axis. Use site 0 (the selected ring is translation-symmetric so any
    # site works). Use the linearised-response coefficient (don't add full ‖σ‖
    # for trace-1 reasons; just a small probe).
    rho0_mixed = np.eye(d, dtype=complex) / d
    eps = 1.0 / d  # small probe amplitude
    rho0_x = rho0_mixed + eps * site_op(N, 0, "X") / d
    rho0_y = rho0_mixed + eps * site_op(N, 0, "Y") / d
    rho0_z = rho0_mixed + eps * site_op(N, 0, "Z") / d

    # propagator: ρ(t) = unvec(exp(L · t) · vec(ρ0))
    ts = np.linspace(0.0, t_max, n_steps + 1)

    def propagate_and_track(rho0):
        v0 = vec_column(rho0)
        tr_x, tr_y, tr_z = [], [], []
        for t in ts:
            vt = np.linalg.matrix_power if False else None  # placeholder
            # use direct exponential per timestep; cheap enough for d=16
            U = (np.linalg.matrix_power if False else None)  # ignore
            del U, vt
            from scipy.linalg import expm
            vt = expm(L * t) @ v0
            rhot = unvec_column(vt, d)
            tr_x.append(np.trace(Mx @ rhot).real)
            tr_y.append(np.trace(My @ rhot).real)
            tr_z.append(np.trace(Mz @ rhot).real)
        return np.array(tr_x), np.array(tr_y), np.array(tr_z)

    print("Initial state: ρ_mixed + ε · X_0 / d  (probe along x at site 0)")
    mx_x, my_x, mz_x = propagate_and_track(rho0_x)
    rate_x_x, _ = fit_exp_tail(ts, mx_x)
    rate_x_y, _ = fit_exp_tail(ts, my_x)
    print(f"  ⟨Mx⟩(t) specified tail-fit rate: {rate_x_x:+.6f}")
    print(f"  ⟨My⟩(t) specified tail-fit rate: {rate_x_y:+.6f}")
    print()

    print("Initial state: ρ_mixed + ε · Y_0 / d  (probe along y at site 0)")
    mx_y, my_y, mz_y = propagate_and_track(rho0_y)
    rate_y_x, _ = fit_exp_tail(ts, mx_y)
    rate_y_y, _ = fit_exp_tail(ts, my_y)
    print(f"  ⟨Mx⟩(t) specified tail-fit rate: {rate_y_x:+.6f}")
    print(f"  ⟨My⟩(t) specified tail-fit rate: {rate_y_y:+.6f}")
    print()

    print("Initial state: ρ_mixed + ε · Z_0 / d  (probe along z at site 0)")
    mx_z, my_z, mz_z = propagate_and_track(rho0_z)
    rate_z_z, _ = fit_exp_tail(ts, mz_z)
    print(f"  ⟨Mz⟩(t) specified tail-fit rate: {rate_z_z:+.6f}")
    print()

    print("-" * 96)
    print("Numerical comparison (tail-fit convergence not established):")
    tol_rel = 0.02  # 2% tail-fit tolerance against asymptotic slow-mode prediction
    def ok(meas, pred):
        return abs(meas - pred) / max(abs(pred), 1e-12) < tol_rel

    print(f"  Slowest Re-flavor mode  |Re(λ_k=1)| = {abs(re_modes[0][1]):.6f}")
    print(f"  ⟨Mx⟩ tail-fit rate from x-probe     = {rate_x_x:.6f}   "
          f"(within {tol_rel*100:.0f}% threshold: {ok(rate_x_x, abs(re_modes[0][1]))})")
    print(f"  ⟨Mz⟩ tail-fit rate from z-probe     = {rate_z_z:.6f}   "
          f"(within {tol_rel*100:.0f}% threshold: {ok(rate_z_z, abs(re_modes[0][1]))})")
    print()
    print(f"  Slowest Im-flavor mode  |Re(λ_k=2)| = {abs(im_modes[0][1]):.6f}")
    print(f"  ⟨My⟩ tail-fit rate from y-probe     = {rate_y_y:.6f}   "
          f"(within {tol_rel*100:.0f}% threshold: {ok(rate_y_y, abs(im_modes[0][1]))})")
    print()
    print(f"  Axis-probe inverse-rate ratio       = "
          f"{(1.0/rate_x_x)/(1.0/rate_y_y):.4f}")
    print(f"  Slow-mode Im/Re rate ratio          = "
          f"{abs(im_modes[0][1])/abs(re_modes[0][1]):.4f}")
    print("  Their relation is a convergence question, not a causal assignment.")
    print()


def main():
    print("=" * 96)
    print("Numerical per-Painter Y/non-Y classification and axis-probe tail fits")
    print("=" * 96)
    print()
    print("This script (a) numerically classifies the slow-mode Pauli weights at a")
    print("1e-8 tolerance; and (b) propagates selected single-site Pauli probes under")
    print("the full Liouvillian to report specified tail fits.")
    print()
    print("These outputs are selected-model observations. They do not define T1/T2,")
    print("FID, TROSY, EXSY, or a material observable without a separately specified")
    print("degree of freedom, Hamiltonian, bath, preparation, and measurement operator.")

    run(N=4, h_zeeman=0.5, gamma=1.0, t_max=20.0, n_steps=200, n_slowest=8)


if __name__ == "__main__":
    main()
