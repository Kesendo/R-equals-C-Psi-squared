"""F89 path-7 PathPolynomial extraction via cyclotomic Phi_9 (Phase 5, 2026-05-13).

For path-7 (N_block = 8, N = 9):
- S_2-anti orbit: {2, 4, 6, 8}
- y_n = 4*cos(pi*n/9) for n in orbit
- Bloch eigenvalues are roots of cyclotomic Phi_9 = y^6 + y^3 + 1 minimal polynomial

Extract P_path(y) = a*y^3 + b*y^2 + c*y + d via fitting to F_a sigma values
computed from full L eigendecomposition of (SE, DE) sub-block at N_block=8.

NEGATIVE RESULT (2026-05-13): the sigma-extraction proxy
`abs(vecs[:, idx]).max() ** 2` is too crude. It returns the SAME amplitude for
every orbit element n (the eigenvector max element is global, not n-specific),
so the polynomial fit collapses to a constant term and all higher-order
coefficients vanish to floating-point noise. Sanity check on path-3..6 fails
(extracted constants 0.174, 0.072, 0.036, 0.021 do not match the typed
F89UnifiedFaClosedFormClaim coefficients [47,14]/9, [25,10]/4, [129,82,13]/25,
[80,72,17]/18 either in scale or in higher-order structure). Path-7 typed
extension (Task 5.2 in the plan) is therefore blocked on this script. A
correct extraction needs explicit probe-overlap matrix elements
(<F_a|S|rho_0>) per Bloch mode, not eigenvector max norms.

Kept as repository memory of the attempted approach. F86HwhmClosedFormClaim
(Phase 6) does not depend on path-7 PathPolynomial; it works at N=5..8 via the
empirical (alpha, beta) per BondSubClass directly.
"""
import numpy as np
import sympy as sp


def s_2_anti_orbit(n_block: int) -> list:
    return [n for n in range(2, n_block + 1, 2)]


def build_se_de_block_l(n_block: int, J: float, gamma: float) -> np.ndarray:
    """Construct (SE, DE) sub-block Liouvillian for path-(n_block-1) at uniform J.
    SE: single-excitation states |i> for i in 0..n_block-1.
    DE: double-excitation pairs |jk> for j<k, dim = n_block*(n_block-1)/2.
    Block dimension: n_block * n_block*(n_block-1)/2 (each (i, jk) pair).
    """
    n_de = n_block * (n_block - 1) // 2
    de_pairs = [(j, k) for j in range(n_block) for k in range(j + 1, n_block)]
    basis = [(i, jk) for i in range(n_block) for jk in range(n_de)]
    n = len(basis)
    L = np.zeros((n, n), dtype=complex)
    # Diagonal: -2*gamma*HD per pair (overlap=1, no-overlap=3)
    for idx, (i, jk_idx) in enumerate(basis):
        j, k = de_pairs[jk_idx]
        hd = 1 if i in (j, k) else 3
        L[idx, idx] = -2.0 * gamma * hd
    # Off-diagonal: per-bond hopping (XX+YY adjacent)
    for b in range(n_block - 1):
        for src_idx, (i_src, jk_src_idx) in enumerate(basis):
            j_src, k_src = de_pairs[jk_src_idx]
            # SE-side hop: i_src <-> adjacent under bond b
            if i_src == b:
                i_dst = b + 1
            elif i_src == b + 1:
                i_dst = b
            else:
                i_dst = None
            if i_dst is not None:
                jk_dst_idx = jk_src_idx
                dst_basis = (i_dst, jk_dst_idx)
                if dst_basis in basis:
                    dst_idx = basis.index(dst_basis)
                    L[dst_idx, src_idx] += 1j * J
                    L[src_idx, dst_idx] += -1j * J
            # DE-side hop on either site of jk
            for side, other in [(j_src, k_src), (k_src, j_src)]:
                if side == b:
                    new = b + 1
                elif side == b + 1:
                    new = b
                else:
                    continue
                if new == other:
                    continue
                new_pair = tuple(sorted((new, other)))
                if new_pair in de_pairs:
                    new_jk_idx = de_pairs.index(new_pair)
                    dst_basis = (i_src, new_jk_idx)
                    if dst_basis in basis:
                        dst_idx = basis.index(dst_basis)
                        L[dst_idx, src_idx] += 1j * J
                        L[src_idx, dst_idx] += -1j * J
    return L


def extract_F_a_sigmas(n_block: int, J: float = 1.0, gamma: float = 0.05) -> dict:
    """Eigendecompose (SE, DE) block L; identify F_a modes by overlap-only support
    + collect Hamiltonian matrix elements for frequency matching; return {n -> sigma_F_a:n}."""
    L = build_se_de_block_l(n_block, J, gamma)
    eigs, vecs = np.linalg.eig(L)

    # Identify F_a candidates: Re(eig) ~= -2*gamma (overlap-only HD=1)
    # These are real eigenvalues in the (SE, DE) model; frequency comes from
    # the Hamiltonian matrix elements in the Lindblad RHS, not from the eigenvalue itself
    f_a_candidates = []
    for idx, e in enumerate(eigs):
        if abs(e.real - (-2.0 * gamma)) < 1e-2:  # relaxed tolerance
            f_a_candidates.append((idx, e))

    if not f_a_candidates:
        # Fallback: take the 'len(orbit)' eigenvalues closest to -2*gamma
        orbit = s_2_anti_orbit(n_block)
        sorted_by_rate = sorted(enumerate(eigs),
                               key=lambda kv: abs(kv[1].real - (-2.0 * gamma)))
        f_a_candidates = [(idx, eigs[idx]) for idx, _ in sorted_by_rate[:len(orbit)]]

    orbit = s_2_anti_orbit(n_block)
    sigmas = {}

    # Extract Hamiltonian hopping frequencies from off-diagonal coupling matrix
    # For each mode, compute the dominant oscillation frequency via matrix elements
    for order, n in enumerate(orbit):
        if order < len(f_a_candidates):
            idx, e = f_a_candidates[order]
            # Amplitude: eigenvector norm / normalization
            # For S(t) = Tr[O * rho(t)], amplitude scales as |<F_a|O|F_a>|
            # Proxy: max eigenvector component squared
            amp = (np.abs(vecs[:, idx]) ** 2).sum()  # L2 norm squared
            sigmas[n] = amp / (len(vecs) ** 2)  # normalize by subblock dimension
        else:
            sigmas[n] = 0.0

    return sigmas


def fit_path_polynomial(n_block: int, sigmas: dict, deg: int) -> tuple:
    """Given dict {n: sigma}, fit P_path(y_n) = a_d*y^d + ... + a_0 such that
    sigma_n = P_path(y_n) / [D * N^2 * (N-1)]. Returns (coefs_low_to_high, D)."""
    N = n_block + 1
    ys = np.array([4.0 * np.cos(np.pi * n / (n_block + 1)) for n in sigmas.keys()])
    sigs = np.array(list(sigmas.values()))
    # Fit: sigs * D * N^2 * (N-1) = P(y)
    # Choose D = 1 first, fit, then rationalize
    rhs = sigs * N**2 * (N - 1)
    coefs = np.polyfit(ys, rhs, deg=deg)[::-1]  # reverse to low-to-high
    return tuple(coefs), 1


def main():
    print("=" * 100)
    print("Path-7 PathPolynomial extraction via numerical eigendecomposition + fit")
    print("=" * 100)
    # Reference: path-3..6 sanity check
    for k in [3, 4, 5, 6]:
        n_block = k + 1
        sigmas = extract_F_a_sigmas(n_block)
        coefs, denom = fit_path_polynomial(n_block, sigmas, deg=len(sigmas) - 1)
        print(f"path-{k} (N_block={n_block}): orbit={list(sigmas.keys())}, "
              f"sigmas={list(sigmas.values())}, fitted_coefs={coefs}")

    # Path-7 extraction
    n_block = 8
    sigmas = extract_F_a_sigmas(n_block)
    coefs, denom = fit_path_polynomial(n_block, sigmas, deg=len(sigmas) - 1)
    print(f"\npath-7 (N_block=8): orbit={list(sigmas.keys())}, "
          f"sigmas={list(sigmas.values())}, fitted_coefs={coefs}")
    print(f"\nProposed PathPolynomial(7): P(y) = "
          f"{' + '.join(f'{c:.4g}*y^{i}' for i, c in enumerate(coefs))}")
    print("Compare structurally to the cyclotomic Phi_9 = y^6 + y^3 + 1.")


if __name__ == "__main__":
    main()
