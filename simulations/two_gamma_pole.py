#!/usr/bin/env python3
"""
two_gamma_pole.py -- Which sector produces alpha = 2*gamma_0?

The single-excitation sector never reaches 2*gamma_0 (max is 4/(N+1)).
But the full Liouvillian always has eigenvalues at exactly 2*gamma_0.
Which modes are these?

The Absorption Theorem says alpha = 2*gamma_0 * <n_XY>_B.
For alpha = 2*gamma_0, we need <n_XY>_B = 1 exactly, meaning
the eigenmode has X or Y at site B with certainty.

Question: which eigenmodes have this property, and what is their
Pauli structure?

Date: 2026-04-16
"""

import numpy as np
from itertools import product as iprod
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
OUT_PATH = RESULTS_DIR / "two_gamma_pole.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = [I2, X, Y, Z]
PAULI_NAMES = ['I', 'X', 'Y', 'Z']


def kron(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, N):
    factors = [I2] * N
    factors[site] = op
    return kron(*factors)


def liouvillian(H, jump_ops):
    d = H.shape[0]
    Idd = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Idd, H) - np.kron(H.T, Idd))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Idd, LdL)
              - 0.5 * np.kron(LdL.T, Idd))
    return L


def build_xy_chain(N, J=1.0):
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        H += J * 0.5 * (site_op(X, i, N) @ site_op(X, i + 1, N)
                        + site_op(Y, i, N) @ site_op(Y, i + 1, N))
    return H


def build_pauli_basis(N):
    """Build all 4^N Pauli strings as matrices and index tuples."""
    d = 2 ** N
    all_idx = list(iprod(range(4), repeat=N))
    pmats = []
    for idx in all_idx:
        m = PAULIS[idx[0]]
        for k in range(1, N):
            m = np.kron(m, PAULIS[idx[k]])
        pmats.append(m / d)  # normalize: Tr(P_i P_j) = delta_ij * d
    return all_idx, pmats


def pauli_label(idx):
    return ''.join(PAULI_NAMES[i] for i in idx)


def xy_weight_at_site(idx, site):
    """1 if Pauli string has X or Y at site, 0 otherwise."""
    return 1 if idx[site] in (1, 2) else 0


def total_xy_weight(idx):
    return sum(1 for i in idx if i in (1, 2))


# =====================================================================
if __name__ == "__main__":
    log("WHICH SECTOR PRODUCES alpha = 2*gamma_0?")
    log()

    gamma_0 = 1e-4
    J = 1.0

    for N in range(3, 7):
        log("=" * 70)
        log(f"N={N}, dim(L)={4**N}")
        log("=" * 70)

        B = N - 1  # endpoint

        H = build_xy_chain(N, J)
        L = liouvillian(H, [np.sqrt(gamma_0) * site_op(Z, B, N)])

        eigenvalues, eigenvectors = np.linalg.eig(L)
        alphas = -eigenvalues.real
        normalized = alphas / gamma_0

        # find modes at alpha = 2*gamma_0
        at_two = np.where(np.abs(normalized - 2.0) < 1e-3)[0]
        # find modes at alpha = 0
        at_zero = np.where(np.abs(normalized) < 1e-3)[0]

        log(f"Modes at alpha = 0: {len(at_zero)}")
        log(f"Modes at alpha = 2*gamma_0: {len(at_two)}")

        if N <= 5:
            # decompose the alpha=2*gamma_0 eigenvectors in Pauli basis
            all_idx, pmats = build_pauli_basis(N)
            d = 2 ** N

            log(f"\nPauli decomposition of alpha = 2*gamma_0 modes:")
            for mode_i, ev_idx in enumerate(at_two[:4]):  # show first 4
                evec = eigenvectors[:, ev_idx]
                # reshape to density matrix
                rho = evec.reshape(d, d, order='F')

                # project onto Pauli basis
                log(f"\n  Mode {mode_i} (alpha/gamma_0 = {normalized[ev_idx]:.6f}):")
                dominant = []
                for p_i, (idx, pm) in enumerate(zip(all_idx, pmats)):
                    coeff = np.trace(pm.conj().T @ rho) * d  # undo normalization
                    if abs(coeff) > 0.01:
                        label = pauli_label(idx)
                        w_B = xy_weight_at_site(idx, B)
                        w_total = total_xy_weight(idx)
                        dominant.append((abs(coeff), label, w_B, w_total, coeff))

                dominant.sort(reverse=True)
                for amp, label, w_B, w_total, coeff in dominant[:8]:
                    log(f"    {label:>8}  |c|={amp:.4f}  "
                        f"n_XY(B)={w_B}  w_total={w_total}")

        # statistics: what is <n_XY>_B for the modes at alpha=2*gamma_0?
        if N <= 5:
            log(f"\n<n_XY>_B for alpha=2*gamma_0 modes (via Absorption Theorem):")
            for ev_idx in at_two:
                alpha_val = alphas[ev_idx]
                nxy_B = alpha_val / (2 * gamma_0)
                log(f"  mode: alpha/gamma_0 = {normalized[ev_idx]:.6f}, "
                    f"<n_XY>_B = {nxy_B:.6f}")

        # what is the total XY-weight of these modes?
        # this requires the Pauli decomposition
        if N <= 4:
            all_idx, pmats = build_pauli_basis(N)
            d = 2 ** N

            log(f"\nXY-weight distribution of alpha=0 vs alpha=2*gamma_0 modes:")

            for label_name, mode_indices in [("alpha=0", at_zero),
                                              ("alpha=2g0", at_two)]:
                w_counts = {}
                for ev_idx in mode_indices:
                    evec = eigenvectors[:, ev_idx]
                    rho = evec.reshape(d, d, order='F')
                    # find dominant Pauli string
                    best_w = -1
                    best_amp = 0
                    for p_i, (idx, pm) in enumerate(zip(all_idx, pmats)):
                        coeff = abs(np.trace(pm.conj().T @ rho) * d)
                        if coeff > best_amp:
                            best_amp = coeff
                            best_w = total_xy_weight(idx)
                    w_counts[best_w] = w_counts.get(best_w, 0) + 1

                log(f"  {label_name}: total XY-weight of dominant component: {w_counts}")

        log()

    # === Pattern check ===
    log("\n" + "=" * 70)
    log("MULTIPLICITY PATTERN")
    log("=" * 70)
    log(f"{'N':>4} {'at 0':>8} {'at 2g0':>8} {'pattern?':>12}")
    for N in range(3, 8):
        H = build_xy_chain(N, J)
        L = liouvillian(H, [np.sqrt(gamma_0) * site_op(Z, N - 1, N)])
        eigs = np.linalg.eigvals(L)
        norm = -eigs.real / gamma_0
        n_zero = np.sum(np.abs(norm) < 1e-3)
        n_two = np.sum(np.abs(norm - 2.0) < 1e-3)
        log(f"{N:4d} {n_zero:8d} {n_two:8d} {n_zero}={N+1}? {'yes' if n_zero==N+1 else 'no'}")

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
