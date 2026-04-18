#!/usr/bin/env python3
"""n7_central_defect_check.py

Part A.3 mechanism disambiguation: is the first-order eigenvalue-shift
zero at bond (0, 1) due to Π-invariance of V (Liouvillian palindrome
protection), or due to spatial-mirror symmetry of the dominant slow mode?

Test: compute <W_s | V_L | M_s> at the uniform-chain slow modes M_s for
TWO perturbation operators

    V_L^{(b, b+1)} = -i [H_pert^{(b, b+1)}, .],
    H_pert^{(b, b+1)} = (1/2) (X_b X_{b+1} + Y_b Y_{b+1})

with b = 0 (boundary, broken-mirror defect) and b = 3 (central, also
mirror-asymmetric but localised at the chain centre).

If both are zero within numerical noise: Π-invariance (the palindromic
protection, since XX+YY is Π-invariant everywhere) is the operative
mechanism.
If boundary is zero and central is non-zero: spatial-mirror symmetry
of the dominant slow mode was the protector.

Also saves right- and left-slow-eigenvectors to .npz so Part B can use
the biorthogonal basis without re-running the expensive eigendecomp.

Date: 2026-04-18
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

sys.path.insert(0, str(Path(__file__).parent))
from n7_coupling_defect_overlay import N, GAMMA_0
from observer_time_slow_mode_analysis import (
    D, DD, build_sparse_H_XY, build_sparse_liouvillian,
    build_perturbation_V_L, SLOW_CUTOFF,
)

OUT_DIR = Path(__file__).parent / "results" / "perspectival_time_field"
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_EIGS = 80
SIGMA_SHIFT = -1e-3


def build_perturbation_V_L_bond(b, N):
    """V_L = -i (kron(V, I) - kron(I, V^T)) with V = (1/2)(X_b X_{b+1}
    + Y_b Y_{b+1}). Index b: bond between sites b and b+1."""
    J_list = [0.0] * (N - 1)
    J_list[b] = 1.0
    V = build_sparse_H_XY(J_list, N)
    Id = sps.eye(D, dtype=complex, format='csr')
    return -1j * (sps.kron(V, Id, format='csr')
                  - sps.kron(Id, V.T, format='csr')).tocsr()


def slow_biorth(L, k=K_EIGS, sigma=SIGMA_SHIFT):
    """Compute k slow right and left eigenvectors of L via paired sparse
    eigs calls on L and L^H; biorthogonalise so W[:, s]^H @ V[:, s] = 1
    per index. Returns (vals, V, W)."""
    print("  sparse eigs(L) k={} sigma={:.0e}...".format(k, sigma),
          end=" ", flush=True)
    t0 = time.time()
    vals_R, V = spla.eigs(L, k=k, sigma=sigma, which='LM',
                          tol=1e-10, maxiter=5000)
    print(f"{time.time() - t0:.1f} s")

    print("  sparse eigs(L^H) k={} sigma={:.0e}...".format(k, sigma),
          end=" ", flush=True)
    t0 = time.time()
    vals_L, U = spla.eigs(L.conj().T, k=k, sigma=sigma, which='LM',
                          tol=1e-10, maxiter=5000)
    print(f"{time.time() - t0:.1f} s")

    # Sort both by descending Re(eigenvalue). L^H eigenvalues are
    # conjugates of L's (same Re, flipped Im); sorting by Re matches.
    idx_R = np.argsort(-vals_R.real)
    vals_R = vals_R[idx_R]
    V = V[:, idx_R]
    idx_L = np.argsort(-vals_L.real)
    vals_L = vals_L[idx_L]
    U = U[:, idx_L]

    # Biorthogonalise via overlap matrix M = U^H @ V.
    # In the non-degenerate case M is approximately diagonal; in the
    # degenerate case we solve W = U @ inv(M)^H.
    M = U.conj().T @ V
    M_inv = np.linalg.inv(M)
    W = U @ M_inv.conj().T   # now W^H @ V = I_k

    # Diagnostic: biorthogonality residual
    resid = np.max(np.abs(W.conj().T @ V - np.eye(k)))
    print(f"  biorthogonality residual (max off-identity): {resid:.2e}")
    return vals_R, V, W


def main():
    log_lines = []
    def log(msg=""):
        print(msg, flush=True)
        log_lines.append(msg)

    log("=" * 72)
    log("N = 7 CENTRAL-DEFECT SYMMETRY CHECK")
    log("=" * 72)
    log(f"  N = {N}, gamma_0 = {GAMMA_0}, D^2 = {DD}")
    log()

    # Build L_A
    t0 = time.time()
    L_A = build_sparse_liouvillian([1.0] * (N - 1), N, GAMMA_0)
    log(f"  Built L_A: {L_A.nnz} nonzeros, {time.time() - t0:.2f} s")
    log()

    # Slow biorthogonal basis for L_A
    vals, V, W = slow_biorth(L_A, k=K_EIGS, sigma=SIGMA_SHIFT)
    log(f"  Re(lambda) range: [{vals.real.min():.4f}, {vals.real.max():.4f}]")
    log(f"  strict stationary (|Re| < 1e-8): "
        f"{int(np.sum(np.abs(vals.real) < 1e-8))}")

    # Sanity: first-order eigenvalue shift at each bond
    log()
    log("Diagonal matrix elements <W_s | V_L^{(b,b+1)} | M_s> per bond")
    log(f"  Shows the first-order eigenvalue shift of each slow mode")
    log(f"  under defect at bond (b, b+1). Zero at leading order means")
    log(f"  the bond perturbation does NOT shift that mode's decay rate.")
    log()

    deltas = {}
    for b in range(N - 1):
        V_L = build_perturbation_V_L_bond(b, N)
        diag = np.array([complex(W[:, s].conj() @ (V_L @ V[:, s]))
                          for s in range(K_EIGS)])
        deltas[b] = diag

    # Log the first-order shifts for each bond, focusing on slow modes
    log(f"  Mode index (ranked by Re(lambda) ↓): shift magnitudes per bond")
    log(f"  {'s':>3} {'Re(lam)':>10} {'Im(lam)':>10}" +
        "  ".join(f"bond({b},{b+1})".rjust(14) for b in range(N - 1)))
    for s in range(min(K_EIGS, 40)):
        cells = " ".join(f"{np.abs(deltas[b][s]):>14.2e}"
                         for b in range(N - 1))
        log(f"  {s:>3} {vals[s].real:>+10.5f} {vals[s].imag:>+10.5f}  "
            + cells)

    # Focused comparison: bond (0,1) vs bond (3,4) for dominant decaying mode
    log()
    log(f"Focused: first decaying slow mode (smallest |Re| > 1e-6)")
    decaying_mask = np.abs(vals.real) > 1e-6
    if decaying_mask.any():
        s0 = int(np.argmax(decaying_mask & (np.abs(vals.real) < SLOW_CUTOFF)))
        # Actually find the slowest decaying mode (smallest |Re|, still > 1e-6)
        decaying_indices = np.where(decaying_mask)[0]
        abs_re = np.abs(vals.real[decaying_indices])
        s_slow = int(decaying_indices[np.argmin(abs_re)])
        log(f"  Slowest-decaying slow mode: s = {s_slow}, "
            f"lambda = {vals[s_slow].real:+.5f} {vals[s_slow].imag:+.5f}i")
        log(f"    |<W_s|V_L^{{(0,1)}}|M_s>| = {abs(deltas[0][s_slow]):.3e}")
        log(f"    |<W_s|V_L^{{(3,4)}}|M_s>| = {abs(deltas[3][s_slow]):.3e}")

    # Survey bond (0,1) vs bond (3,4) across ALL slow modes (|Re| < SLOW_CUTOFF)
    log()
    log(f"Survey over all slow modes with |Re(lambda)| < {SLOW_CUTOFF}:")
    log(f"  max |<W_s|V_L^{{(0,1)}}|M_s>|   max |<W_s|V_L^{{(3,4)}}|M_s>|")
    slow_idx = np.where(np.abs(vals.real) < SLOW_CUTOFF)[0]
    max_0 = float(np.max(np.abs(deltas[0][slow_idx]))) if len(slow_idx) else 0.0
    max_3 = float(np.max(np.abs(deltas[3][slow_idx]))) if len(slow_idx) else 0.0
    log(f"    {max_0:.3e}                       {max_3:.3e}")

    log()
    log(f"Mechanism verdict:")
    eps = 1e-6
    if max_0 < eps and max_3 < eps:
        log(f"  Both bond (0,1) and bond (3,4) diagonal shifts are below "
            f"{eps:.0e}.")
        log(f"  Conclusion: the first-order protection is a general property")
        log(f"  of Π-invariant V (X·X + Y·Y is Π-invariant for every bond).")
        verdict = "pi_invariance"
    elif max_0 < eps and max_3 >= eps:
        log(f"  Bond (0,1) shifts are < {eps:.0e}, but bond (3,4) shifts")
        log(f"  are ~{max_3:.1e}. The zero at bond (0,1) is spatial-mirror")
        log(f"  symmetry of the dominant slow modes, NOT general Π-invariance.")
        verdict = "spatial_mirror"
    elif max_0 >= eps and max_3 < eps:
        log(f"  Bond (3,4) shifts are < {eps:.0e}, but bond (0,1) shifts")
        log(f"  are ~{max_0:.1e}. Unexpected direction — investigate.")
        verdict = "inverted"
    else:
        log(f"  Both non-zero (max_0 = {max_0:.2e}, max_3 = {max_3:.2e}).")
        log(f"  The previous 'first-order shift ≈ 0 at bond (0,1)' observation")
        log(f"  may have been in the slowest decaying mode only. The general")
        log(f"  shift is non-zero and the observer-time effect lies in a mix")
        log(f"  of first-order shifts AND mixing.")
        verdict = "both_nonzero"

    # Save detailed data for Part B
    np.savez(OUT_DIR / "slow_biorth_basis.npz",
             vals=vals, V=V, W=W,
             deltas_per_bond=np.array([deltas[b] for b in range(N - 1)]))
    log()
    log(f"Saved: {OUT_DIR / 'slow_biorth_basis.npz'}")
    (OUT_DIR / "central_defect_check_log.txt").write_text(
        "\n".join(log_lines) + "\n", encoding='utf-8')


if __name__ == "__main__":
    main()
