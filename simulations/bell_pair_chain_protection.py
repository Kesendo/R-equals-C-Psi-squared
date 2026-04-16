#!/usr/bin/env python3
"""
Bell-pair protection by a dephasing-decoupled chain.

Setup: external reference qubit R (isolated) + N-site XY chain Q_0..Q_{N-1}
with single-site Z-dephasing gamma_0 on the endpoint Q_{N-1}. R is coupled
to the chain only through initial entanglement (no Hamiltonian coupling,
no dephasing on R).

Three encoding variants:
  A: Bell pair |Phi+>_{R,Q_0} (x) |0...0>_rest    (localized at inner end)
  B: (|0>_R|vac>_C + |1>_R|psi_1>_C)/sqrt(2)      (delocalized on k=1 mode)
  C: Bell pair |Phi+>_{R,Q_{N-1}} (x) |0...0>_rest (localized at outer end,
     directly next to dephasing - baseline for "no protection")

Prediction:
  B decays as exp(-alpha_1 t) with alpha_1 = (4 gamma_0/(N+1)) sin^2(pi/(N+1))
    because |vac><psi_1| is a Liouvillian right eigenvector.
  A decays multi-exponentially with long-time tail also at alpha_1.
  C decays fast (essentially 2 gamma_0, directly exposed).

Metric: concurrence C(R,Q_0)(t) from the 4x4 reduced density matrix
of R and Q_0 (trace out the rest of the chain).

Date: 2026-04-16
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import curve_fit
from pathlib import Path
import sys

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

RESULTS_DIR = Path(__file__).parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)
OUT_PATH = RESULTS_DIR / "bell_pair_chain_protection.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()


# ---- Paulis and basic ops ----
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, Ntot):
    """Single-site operator at position `site` in a chain of `Ntot` qubits.
    Site 0 is most significant (leftmost in the kron product)."""
    factors = [I2] * Ntot
    factors[site] = op
    return kron_chain(*factors)


def chain_hamiltonian_on_extended(N, J=1.0):
    """XY Hamiltonian acting only on chain qubits (sites 1..N in the
    extended space, where site 0 is R)."""
    Ntot = N + 1
    d = 2**Ntot
    H = np.zeros((d, d), dtype=complex)
    for i in range(1, N):  # chain internal bonds between site i and i+1
        H += J * 0.5 * (site_op(X, i, Ntot) @ site_op(X, i+1, Ntot) +
                        site_op(Y, i, Ntot) @ site_op(Y, i+1, Ntot))
    return H


def liouvillian_superop(H, jump_ops):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for Lk in jump_ops:
        LdL = Lk.conj().T @ Lk
        L += (np.kron(Lk.conj(), Lk)
              - 0.5 * np.kron(Id, LdL)
              - 0.5 * np.kron(LdL.T, Id))
    return L


def formula_alpha_1(N, gamma_0=1.0):
    """Slowest single-excitation decay rate at endpoint B."""
    return (4.0 * gamma_0 / (N + 1)) * np.sin(np.pi / (N + 1))**2


# ---- Concurrence via Wootters formula ----

# sigma_y tensor sigma_y (constant across N since this is 2-qubit)
SIGMA_YY = np.kron(Y, Y)


def concurrence(rho2q):
    """Wootters concurrence for a 4x4 two-qubit density matrix."""
    rho_tilde = SIGMA_YY @ rho2q.conj() @ SIGMA_YY
    M = rho2q @ rho_tilde
    # eigenvalues of M are non-negative real in exact arithmetic;
    # numerically take real part and clip
    eigs = np.linalg.eigvals(M)
    eigs = np.sort(np.real(eigs))[::-1]  # descending
    eigs = np.clip(eigs, 0.0, None)
    sqrts = np.sqrt(eigs)
    C = sqrts[0] - sqrts[1] - sqrts[2] - sqrts[3]
    return max(0.0, float(np.real(C)))


def reduce_to_2qubit(rho_full, keep_site_a, keep_site_b, Ntot):
    """Trace out all qubits except keep_site_a and keep_site_b from rho_full
    (dim 2^Ntot x 2^Ntot), returning a 4x4 reduced density matrix.
    Qubit ordering in rho_full: site 0 is most significant."""
    # reshape rho into 2x2x...x2 x 2x2x...x2 (Ntot copies each)
    shape = [2] * (2 * Ntot)
    rho = rho_full.reshape(shape)

    # trace out site k means contracting index k (row) with index k+Ntot (col)
    # do it in descending order so earlier indices don't shift
    traced = rho
    current_Ntot = Ntot
    site_map = list(range(Ntot))  # current positions of original sites

    for k in reversed(range(Ntot)):
        if k == keep_site_a or k == keep_site_b:
            continue
        # find k's current position in traced tensor
        pos_row = site_map.index(k)
        pos_col = pos_row + current_Ntot
        traced = np.trace(traced, axis1=pos_row, axis2=pos_col)
        site_map.pop(pos_row)
        current_Ntot -= 1

    # now traced has shape [2, 2] * 2 with the two kept sites; bring into
    # the order (a, b) if necessary
    remaining = site_map  # length 2, e.g. [keep_a, keep_b] or reversed
    # reshape to 4x4: we need (a, b) in the row axes and (a, b) in col axes
    if remaining == [keep_site_a, keep_site_b]:
        rho2 = traced.reshape(4, 4)
    elif remaining == [keep_site_b, keep_site_a]:
        # swap axes 0<->1 in rows and 2<->3 in cols
        rho2 = traced.transpose(1, 0, 3, 2).reshape(4, 4)
    else:
        raise RuntimeError(f"unexpected remaining: {remaining}")
    return rho2


# ---- Initial state preparation for three variants ----

def prepare_variant_A(N):
    """|Phi+>_{R, Q_0} tensor |0...0>_{Q_1..Q_{N-1}}
    Extended ordering: [R=0, Q_0=1, Q_1=2, ..., Q_{N-1}=N] (site 0 is most sig).
    """
    Ntot = N + 1
    d = 2**Ntot
    psi = np.zeros(d, dtype=complex)
    # |0>_R |0>_Q0 |0...0>: all bits zero -> index 0
    psi[0] = 1.0 / np.sqrt(2)
    # |1>_R |1>_Q0 |0...0>: R bit at position 0 of the word, Q0 at position 1.
    # Because site 0 is most significant, bit pattern 110...0 has value
    # 2^(Ntot-1) + 2^(Ntot-2).
    idx = 2**(Ntot - 1) + 2**(Ntot - 2)
    psi[idx] = 1.0 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def prepare_variant_C(N):
    """|Phi+>_{R, Q_{N-1}} tensor |0...0>_rest"""
    Ntot = N + 1
    d = 2**Ntot
    psi = np.zeros(d, dtype=complex)
    psi[0] = 1.0 / np.sqrt(2)
    # R=1, Q_{N-1}=1: R at site 0 -> 2^(Ntot-1). Q_{N-1} at site Ntot-1
    # in extended -> bit 0 (least significant) -> 2^0 = 1.
    idx = 2**(Ntot - 1) + 2**0
    psi[idx] = 1.0 / np.sqrt(2)
    return np.outer(psi, psi.conj())


def prepare_variant_B(N):
    """(|0>_R |vac>_C + |1>_R |psi_1>_C) / sqrt(2)
    where |psi_1>_C = sum_i sqrt(2/(N+1)) sin(pi (i+1)/(N+1)) |single_i>_C
    and |single_i> has Q_i excited."""
    Ntot = N + 1
    d = 2**Ntot
    psi = np.zeros(d, dtype=complex)
    # |0>_R |vac>: index 0
    psi[0] = 1.0 / np.sqrt(2)
    # |1>_R |single_i>_chain: R at site 0, Q_i at site i+1 in extended space
    # bit pattern: R bit set (pos 0 from left in Ntot) + Q_i bit set (pos i+1)
    # index = 2^(Ntot-1) + 2^(Ntot - 1 - (i+1)) = 2^(Ntot-1) + 2^(N-1-i)
    norm = 1.0 / np.sqrt(2) * np.sqrt(2.0 / (N + 1))
    for i in range(N):
        amp = norm * np.sin(np.pi * (i + 1) / (N + 1))
        idx = 2**(Ntot - 1) + 2**(N - 1 - i)
        psi[idx] = amp
    return np.outer(psi, psi.conj())


# ---- Propagation and tracking (with negativity) ----


def partial_transpose_R(rho_rc, N):
    """Partial transpose w.r.t. R qubit (site 0 in extended space).
    rho_rc is 2^(N+1) x 2^(N+1)."""
    dimR = 2
    dimC = 2**N
    rho_reshape = rho_rc.reshape(dimR, dimC, dimR, dimC)
    rho_PT = rho_reshape.transpose(2, 1, 0, 3)
    return rho_PT.reshape(dimR * dimC, dimR * dimC)


def negativity(rho_rc, N):
    """Logarithmic negativity is log2(1+2N), plain negativity sums |lambda|
    over negative eigenvalues of the partial transpose."""
    rho_PT = partial_transpose_R(rho_rc, N)
    rho_PT_H = (rho_PT + rho_PT.conj().T) / 2.0
    eigs = np.linalg.eigvalsh(rho_PT_H)
    neg = float(np.sum(np.abs(eigs[eigs < 0])))
    return neg


def propagate_and_track(rho_0, L_super, times, N, track_sites):
    """Return (times, concurrence_RQ0, negativity_RC) arrays."""
    Ntot = N + 1
    d = 2**Ntot
    rho_vec0 = rho_0.flatten(order='F')
    concs = []
    negs = []
    for t in times:
        rho_vec_t = expm(L_super * t) @ rho_vec0
        rho_t = rho_vec_t.reshape(d, d, order='F')
        rho2 = reduce_to_2qubit(rho_t, track_sites[0], track_sites[1], Ntot)
        concs.append(concurrence(rho2))
        negs.append(negativity(rho_t, N))
    return np.array(concs), np.array(negs)


# =========================================================================
if __name__ == "__main__":
    log("BELL-PAIR PROTECTION BY DEPHASING-DECOUPLED CHAIN")
    log("=" * 70)
    log("Setup: R + N-chain, dephasing gamma_0 at endpoint Q_{N-1}.")
    log("Metrics: C(R,Q_0) = concurrence between R and inner qubit,")
    log("         N(R:C)  = negativity between R and entire chain.")
    log()

    gamma_0 = 0.05
    J = 1.0

    for N in [3, 5]:
        log("-" * 70)
        log(f"N = {N}  (Hilbert dim = {2**(N+1)}, Liouvillian dim = {4**(N+1)})")
        log("-" * 70)

        Ntot = N + 1
        H_ext = chain_hamiltonian_on_extended(N, J)
        L_jump = np.sqrt(gamma_0) * site_op(Z, Ntot - 1, Ntot)
        L_super = liouvillian_superop(H_ext, [L_jump])

        alpha_1 = formula_alpha_1(N, gamma_0)
        log(f"  gamma_0 = {gamma_0}, J = {J}")
        log(f"  formula alpha_1 = (4 gamma_0/(N+1)) sin^2(pi/(N+1)) = {alpha_1:.6f}")
        log(f"  predicted T_2 (= 1/alpha_1) = {1/alpha_1:.3f}")
        log(f"  baseline: 2 gamma_0 (direct dephasing) = {2*gamma_0:.6f}")

        t_max = 5.0 / alpha_1
        times = np.linspace(0, t_max, 40)

        variants = [
            ("A_inner", "Bell R-Q_0 (localized inner)",
             prepare_variant_A(N), (0, 1)),
            ("B_mode",  "Bell R-psi_1 (delocalized k=1 mode)",
             prepare_variant_B(N), (0, 1)),
            ("C_outer", "Bell R-Q_{N-1} (localized outer, no buffer)",
             prepare_variant_C(N), (0, Ntot - 1)),
        ]

        for tag, label, rho_0, track_sites in variants:
            log()
            log(f"  [{tag}] {label}")
            log(f"    initial trace        = {float(np.real(np.trace(rho_0))):.6f}")
            rho2_0 = reduce_to_2qubit(rho_0, track_sites[0], track_sites[1], Ntot)
            log(f"    initial C(R,Q_{track_sites[1]-1 if track_sites[1]>0 else 'R'})  = {concurrence(rho2_0):.4f}")
            log(f"    initial N(R:C)        = {negativity(rho_0, N):.4f}")

            concs, negs = propagate_and_track(rho_0, L_super, times, N, track_sites)

            log(f"    {'t':>8}  {'C(R,Q0)':>10}  {'N(R:C)':>10}")
            for idx in [0, 5, 10, 20, 30, -1]:
                log(f"    {times[idx]:8.2f}  {concs[idx]:10.6f}  {negs[idx]:10.6f}")

            # Fit exponential to the negativity tail (skip first 20%).
            # Negativity is the cleaner metric because for variant B
            # C(R,Q_0) is small throughout (info delocalized).
            tail_start = len(times) // 5
            t_tail = times[tail_start:]
            n_tail = negs[tail_start:]
            mask = n_tail > 1e-6
            t_fit = t_tail[mask]
            n_fit = n_tail[mask]

            def expf(t, alpha, A):
                return A * np.exp(-alpha * t)

            if len(t_fit) >= 3:
                try:
                    popt, _ = curve_fit(expf, t_fit, n_fit,
                                        p0=[alpha_1, n_fit[0]])
                    alpha_fit = popt[0]
                    ratio = alpha_fit / alpha_1
                    ratio_2g = alpha_fit / (2 * gamma_0)
                    log(f"    negativity tail fit: alpha_fit = {alpha_fit:.6f}")
                    log(f"      ratio to alpha_1    : {ratio:.4f}  (1.00 = matches F65 prediction)")
                    log(f"      ratio to 2 gamma_0  : {ratio_2g:.4f}  (1.00 = fully exposed)")
                except Exception as e:
                    log(f"    fit failed: {e}")

        log()

    log("=" * 70)
    log("EXPECTED VERDICT")
    log("=" * 70)
    log("  A_inner: multi-exponential decay; long-time tail at alpha_1")
    log("           with reduced amplitude (|c_1|^2 ~ 2/(N+1) at Q_0)")
    log("  B_mode : pure exponential at alpha_1 (k=1 is Liouvillian eigvec)")
    log("  C_outer: fast decay near 2 gamma_0 (direct dephasing)")
    log()
    log("If results match: chain acts as a dephasing-decoupled buffer;")
    log("the k=1 bonding mode is the optimal encoding for R-chain entanglement")
    log("and achieves T_2 ~ (N+1)^3/(2 pi^2 gamma_0) scaling.")

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
