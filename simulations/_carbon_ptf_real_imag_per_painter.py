"""PTF-style 360-degree view of the carbon Lindbladian: full eigendecomposition,
then per-Painter (per-site partial trace) real- and imaginary-part separation
of each complex eigenmode.

The complex eigenvectors of L on the 4^N-dimensional operator space carry
information about the relaxation structure that the eigenvalue spectrum alone
does not show. Each eigenmode is a complex operator on the d × d state space.
We unvec each eigenmode, partial-trace it down to each site (one "Painter"
per site, like the Perspectival Time Field idea), then separate the resulting
2 × 2 reduced matrices into real and imaginary parts.

What this exposes: the geometric structure of the slowest relaxation modes
as seen from each individual Painter's viewpoint, with the imaginary part
carried by the off-diagonal coherent component, the real part by the
populations and dissipative content.

We focus on the carbon-relevant configuration: cyclobutadiene C₄ ring with
Hückel hopping plus a y-direction Zeeman perturbation (which gives the
relaxation structure a substantive Π-anti component), under Holstein on-site
dephasing.

Output per slowest mode:
  - Eigenvalue λ (real-part = decay rate, imag-part = oscillation frequency)
  - For each of N Painters (one per site):
      Re-part = local population / dissipative structure  (real 2×2)
      Im-part = local coherent / oscillatory structure   (real 2×2)
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
    """Reshape a length-d² complex vector (column-major vec) into a d × d matrix."""
    return v.reshape(d, d, order="F")


def partial_trace_site(rho, N, keep_site):
    """Partial trace of a d × d operator over all sites except keep_site.
    Returns a 2 × 2 complex matrix on the kept site."""
    rho_t = rho.reshape([2] * N + [2] * N)
    keep_axes = [keep_site, N + keep_site]
    sum_axes = list(range(N)) + list(range(N, 2 * N))
    for a in keep_axes:
        sum_axes.remove(a)
    for i in range(N):
        if i == keep_site:
            continue
        rho_t = np.trace(rho_t, axis1=0, axis2=N - len([s for s in range(i) if s != keep_site]))
    return rho_t


def partial_trace_simple(rho, N, keep_site):
    """Simpler partial trace via einsum over tensor-product indices."""
    rho_tensor = rho.reshape([2] * N + [2] * N)
    keep_pair = (keep_site, N + keep_site)
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


def fmt_2x2(M, prec=3):
    return "[" + " ".join(f"{x:+.{prec}f}" for x in M.flatten()) + "]"


def run(N, h_zeeman=0.5, gamma=1.0, n_slowest=8):
    print()
    print("=" * 96)
    print(f"PTF 360-degree view: N = {N} carbon ring")
    print(f"H = Hückel + {h_zeeman:.2f}·Σ_l Y_l (y-Zeeman)")
    print(f"bath: Holstein c_l = Z_l per site, γ = {gamma}")
    print("=" * 96)

    H = hueckel_ring_H(N) + h_zeeman * zeeman_y_total(N)
    c_holstein = [site_op(N, l, "Z") for l in range(N)]
    L = lindbladian_vec(H, c_holstein, [gamma] * N)

    d = 2**N
    print(f"Operator-space dim 4^{N} = {4**N}, state-space dim 2^{N} = {d}")
    print()

    eigvals, eigvecs = np.linalg.eig(L)
    order = np.argsort(eigvals.real)[::-1]  # slowest first (least-negative real part)
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    print("Slowest eigenmodes (sorted by Re(λ), least-negative first):")
    print(f"  showing {n_slowest} of {len(eigvals)}")
    print()

    for k in range(n_slowest):
        lam = eigvals[k]
        v = eigvecs[:, k]
        rho = unvec_column(v, d)

        decay = lam.real
        omega = lam.imag

        # operator-level magnitudes
        norm_total = np.linalg.norm(rho)
        re_norm = np.linalg.norm(rho.real)
        im_norm = np.linalg.norm(rho.imag)

        print(f"  Mode {k}: λ = {decay:+.4f} {omega:+.4f}i  (decay {decay:+.4f}, ω = {omega:+.4f})")
        print(f"           ‖ρ‖ = {norm_total:.4f}, ‖Re(ρ)‖ = {re_norm:.4f}, ‖Im(ρ)‖ = {im_norm:.4f}")

        for site in range(N):
            rho_i = partial_trace_simple(rho, N, site)
            re_i = rho_i.real
            im_i = rho_i.imag
            re_n = np.linalg.norm(re_i)
            im_n = np.linalg.norm(im_i)
            print(f"    Painter {site}: ‖Re‖={re_n:.4f} {fmt_2x2(re_i, 3)}    "
                  f"‖Im‖={im_n:.4f} {fmt_2x2(im_i, 3)}")
        print()


def main():
    print("=" * 96)
    print("PTF 360-degree view: Re/Im decomposition of complex eigenmodes per Painter")
    print("=" * 96)
    print()
    print("Each complex eigenmode of L is an operator on d × d state space. We partial-")
    print("trace each eigenmode to one site (a 'Painter' perspective), then separate the")
    print("resulting 2 × 2 reduced matrix into real and imaginary parts.")
    print()
    print("Re-part: local populations and dissipative structure.")
    print("Im-part: local coherent / oscillatory structure.")
    print()
    print("Read across the N Painters: does the slowest-mode structure show a symmetric")
    print("rotation around the ring, or a localised pattern, or something else?")

    run(N=4, h_zeeman=0.5, gamma=1.0, n_slowest=8)


if __name__ == "__main__":
    main()
