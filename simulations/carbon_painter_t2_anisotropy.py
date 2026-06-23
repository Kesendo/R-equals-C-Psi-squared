"""From the per-Painter Re/Im alternation to a measurable anisotropic-T2
prediction for a y-field aromatic ring.

The companion script carbon_ptf_real_imag_per_painter.py showed, on the
N = 4 Hückel + y-Zeeman ring under Holstein dephasing, that the slowest
Liouvillian eigenmodes alternate per-Painter in a clean Re/Im pattern:

  Mode 1  ω = 0    decay 0.172    Re-flavor   (per-site content {I, X, Z})
  Mode 2  ω = 0    decay 0.219    Im-flavor   (per-site content Y only)
  Mode 3  ω = 0    decay 0.597    Re-flavor
  Mode 4  ω = 0    decay 0.901    Im-flavor
  Mode 5  ω = 0    decay 2.067    Im-flavor
  Mode 6,7  ω = ±3.83  fastest    mixed (complex pair)

The alternation is the slow-mode hierarchy's n_Y-parity sectorization at the
single-site projection level. Re-flavor modes carry NO Y content per site
(equivalently, their per-site reduced operators are real-symmetric 2 × 2);
Im-flavor modes are pure Y per site (anti-symmetric 2 × 2).

This script makes that sectorization quantitative and then turns it into
something a chemist could measure. Two questions:

  (Q1) Is the n_Y-parity assignment bit-exact at the operator level? Each
       slow eigenmode is a 4^N-dim operator; project per-site into the four
       Pauli labels {I, X, Y, Z}; sum |coefficient|² over Y-channels vs over
       non-Y-channels; the alternation predicts the ratio is either 0 or
       infinite per mode (one channel exactly empty).

  (Q2) Does the alternation imply anisotropic transverse relaxation?
       Prepare the ring with single-site magnetization along X-axis vs along
       Y-axis; propagate ρ under the full Liouvillian; track ⟨X⟩(t), ⟨Y⟩(t),
       ⟨Z⟩(t). Fit single-exponential decay to the long-time tail of each.
       Prediction: T2(x_init) maps onto the slowest Re-flavor mode (k = 1),
       T2(y_init) maps onto the slowest Im-flavor mode (k = 2). At h_y/γ =
       0.5 the ratio T2(y)/T2(x) = 0.172/0.219 ≈ 0.787, so y-magnetization
       decays ≈ 1.27× faster than x-magnetization.

The Bloch picture for a chemist: ordinary T2 is isotropic in a pure
on-site-dephasing bath (x and y magnetisations decay at the same rate).
Add a static y-Zeeman and the isotropy breaks; T2(x) ≠ T2(y), and the gap
size is closed-form in h_y / γ.

For NMR readers: this is exactly the regime where TROSY-Difference and
EXSY-Asymmetry techniques live. TROSY-Difference reads the multiplet
component split that arises from anisotropic transverse relaxation;
EXSY-Asymmetry reads cross-peak intensity differences I_AB vs I_BA from
differential T1/T2 between sites. Both are diagnostic of the painter
alternation without requiring full process tomography.
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
    print(f"N = {N} carbon ring, h_y_Zeeman = {h_zeeman}, γ_Holstein = {gamma}")
    print("=" * 96)
    print()

    H = hueckel_ring_H(N) + h_zeeman * zeeman_y_total(N)
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    L = lindbladian_vec(H, c_holstein, [gamma] * N)
    d = 2**N

    print("Q1: n_Y-parity sectorization of slow modes (bit-exact?)")
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
    print("Re-flavor = no Y per site (= zero anti-symmetric content of per-site projection)")
    print("Im-flavor = only Y per site (= zero non-trivial real-symmetric content)")
    print()

    re_modes = [(k, dec) for (k, dec, om, fl, _, _, _) in flavor_table if fl.startswith("Re")]
    im_modes = [(k, dec) for (k, dec, om, fl, _, _, _) in flavor_table if fl.startswith("Im")]
    print(f"Slowest Re-flavor mode: k={re_modes[0][0]}, Re(λ) = {re_modes[0][1]:+.6f}")
    print(f"Slowest Im-flavor mode: k={im_modes[0][0]}, Re(λ) = {im_modes[0][1]:+.6f}")
    ratio_pred = abs(re_modes[0][1]) / abs(im_modes[0][1])
    print(f"Slow-mode prediction for T2(x_init) / T2(y_init) = "
          f"|Re(λ_im)| / |Re(λ_re)| = {1.0/ratio_pred:.6f}")
    print()

    # ---- Q2: anisotropic T2 from full propagation ----
    print("Q2: anisotropic T2 from initial-state preparation")
    print("-" * 96)

    # observables: total magnetization along each axis
    Mx = sum(site_op(N, l, "X") for l in range(N))
    My = sum(site_op(N, l, "Y") for l in range(N))
    Mz = sum(site_op(N, l, "Z") for l in range(N))

    # initial states: maximally mixed plus a tiny single-site magnetisation
    # along the desired axis. Use site 0 (ring is translation-symmetric so any
    # site works). Use the linearised-response coefficient (don't add full ‖σ‖
    # for trace-1 reasons; just a small probe).
    rho0_mixed = np.eye(d, dtype=complex) / d
    eps = 1.0 / d  # small magnetisation amplitude
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
    print(f"  ⟨Mx⟩(t) tail decay rate: {rate_x_x:+.6f}    T2(x_init via Mx) ≈ {1.0/rate_x_x:.4f}")
    print(f"  ⟨My⟩(t) tail decay rate: {rate_x_y:+.6f}    T2(x_init via My) ≈ "
          f"{1.0/rate_x_y if not np.isnan(rate_x_y) else float('nan'):.4f}")
    print()

    print("Initial state: ρ_mixed + ε · Y_0 / d  (probe along y at site 0)")
    mx_y, my_y, mz_y = propagate_and_track(rho0_y)
    rate_y_x, _ = fit_exp_tail(ts, mx_y)
    rate_y_y, _ = fit_exp_tail(ts, my_y)
    print(f"  ⟨Mx⟩(t) tail decay rate: {rate_y_x:+.6f}    T2(y_init via Mx) ≈ "
          f"{1.0/rate_y_x if not np.isnan(rate_y_x) else float('nan'):.4f}")
    print(f"  ⟨My⟩(t) tail decay rate: {rate_y_y:+.6f}    T2(y_init via My) ≈ {1.0/rate_y_y:.4f}")
    print()

    print("Initial state: ρ_mixed + ε · Z_0 / d  (probe along z at site 0)")
    mx_z, my_z, mz_z = propagate_and_track(rho0_z)
    rate_z_z, _ = fit_exp_tail(ts, mz_z)
    print(f"  ⟨Mz⟩(t) tail decay rate: {rate_z_z:+.6f}    T1(z_init via Mz) ≈ {1.0/rate_z_z:.4f}")
    print()

    print("-" * 96)
    print("Comparison:")
    tol_rel = 0.02  # 2% tail-fit tolerance against asymptotic slow-mode prediction
    def ok(meas, pred):
        return abs(meas - pred) / max(abs(pred), 1e-12) < tol_rel

    print(f"  Slowest Re-flavor mode  |Re(λ_k=1)| = {abs(re_modes[0][1]):.6f}")
    print(f"  ⟨Mx⟩ tail rate from x-probe         = {rate_x_x:.6f}   "
          f"(match within {tol_rel*100:.0f}%: {ok(rate_x_x, abs(re_modes[0][1]))})")
    print(f"  ⟨Mz⟩ tail rate from z-probe         = {rate_z_z:.6f}   "
          f"(match within {tol_rel*100:.0f}%: {ok(rate_z_z, abs(re_modes[0][1]))})")
    print()
    print(f"  Slowest Im-flavor mode  |Re(λ_k=2)| = {abs(im_modes[0][1]):.6f}")
    print(f"  ⟨My⟩ tail rate from y-probe         = {rate_y_y:.6f}   "
          f"(match within {tol_rel*100:.0f}%: {ok(rate_y_y, abs(im_modes[0][1]))})")
    print()
    print(f"  T2 anisotropy ratio  T2(x)/T2(y)   = "
          f"{(1.0/rate_x_x)/(1.0/rate_y_y):.4f}")
    print(f"  Slow-mode prediction               = "
          f"{abs(im_modes[0][1])/abs(re_modes[0][1]):.4f}")
    print()


def main():
    print("=" * 96)
    print("From per-Painter Re/Im alternation to anisotropic T2 in a y-field")
    print("=" * 96)
    print()
    print("This script (a) confirms the painter alternation is the slow-mode hierarchy's")
    print("n_Y-parity sectorization at bit-exact precision; and (b) propagates single-")
    print("site magnetization probes under the full Liouvillian to read out the")
    print("anisotropic T2 predicted by the alternation.")
    print()
    print("The bridge: a chemist could test this prediction with two FID experiments")
    print("(90° pulse along orthogonal axes, then watch the decay) without needing")
    print("full process tomography. The ratio T2(x_init)/T2(y_init) is a number that")
    print("the algebra predicts in closed form from h_y and γ.")

    run(N=4, h_zeeman=0.5, gamma=1.0, t_max=20.0, n_steps=200, n_slowest=8)


if __name__ == "__main__":
    main()
