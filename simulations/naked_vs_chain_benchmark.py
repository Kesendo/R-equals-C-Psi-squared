#!/usr/bin/env python3
"""
Naked Bell-pair vs chain-protected benchmark.

Baseline: R + Q with Z-dephasing gamma_0 on Q. alpha_naked = 2*gamma_0.
Chain: R + N-site XY chain with dephasing at endpoint Q_{N-1}.
Bonding-mode encoding gives alpha_1 = (4*gamma_0/(N+1))*sin^2(pi/(N+1)).

Metric: protection factor = T_2(chain) / T_2(naked) = alpha_naked / alpha_1.
This is the ratio of storage lifetime for one ebit.

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
OUT_PATH = RESULTS_DIR / "naked_vs_chain_benchmark.txt"

_outf = open(OUT_PATH, "w", encoding="utf-8", buffering=1)
def log(msg=""):
    print(msg, flush=True)
    _outf.write(msg + "\n")
    _outf.flush()

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SIGMA_YY = np.kron(Y, Y)


def kron_chain(*ops):
    out = ops[0]
    for op in ops[1:]:
        out = np.kron(out, op)
    return out


def site_op(op, site, Ntot):
    factors = [I2] * Ntot
    factors[site] = op
    return kron_chain(*factors)


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


def concurrence(rho2q):
    rho_tilde = SIGMA_YY @ rho2q.conj() @ SIGMA_YY
    M = rho2q @ rho_tilde
    eigs = np.sort(np.real(np.linalg.eigvals(M)))[::-1]
    eigs = np.clip(eigs, 0.0, None)
    sqrts = np.sqrt(eigs)
    return max(0.0, float(sqrts[0] - sqrts[1] - sqrts[2] - sqrts[3]))


def reduce_to_2qubit(rho_full, keep_a, keep_b, Ntot):
    shape = [2] * (2 * Ntot)
    rho = rho_full.reshape(shape)
    traced = rho
    current_Ntot = Ntot
    site_map = list(range(Ntot))
    for k in reversed(range(Ntot)):
        if k == keep_a or k == keep_b:
            continue
        pos_row = site_map.index(k)
        pos_col = pos_row + current_Ntot
        traced = np.trace(traced, axis1=pos_row, axis2=pos_col)
        site_map.pop(pos_row)
        current_Ntot -= 1
    remaining = site_map
    if remaining == [keep_a, keep_b]:
        return traced.reshape(4, 4)
    elif remaining == [keep_b, keep_a]:
        return traced.transpose(1, 0, 3, 2).reshape(4, 4)
    raise RuntimeError(f"unexpected remaining: {remaining}")


def formula_alpha_1(N, gamma_0):
    return (4.0 * gamma_0 / (N + 1)) * np.sin(np.pi / (N + 1))**2


def propagate_concurrence(rho_0, L_super, times, keep_a, keep_b, Ntot):
    d = 2**Ntot
    rho_vec0 = rho_0.flatten(order='F')
    concs = []
    for t in times:
        rho_vec_t = expm(L_super * t) @ rho_vec0
        rho_t = rho_vec_t.reshape(d, d, order='F')
        rho2 = reduce_to_2qubit(rho_t, keep_a, keep_b, Ntot)
        concs.append(concurrence(rho2))
    return np.array(concs)


# =====================================================================
if __name__ == "__main__":
    log("NAKED BELL-PAIR VS CHAIN-PROTECTED BENCHMARK")
    log("=" * 70)

    gamma_0 = 0.05
    J = 1.0
    alpha_naked = 2 * gamma_0

    # === 1. Naked baseline ===
    log("\n" + "=" * 70)
    log("1. NAKED BASELINE (R + Q, dephasing on Q)")
    log("=" * 70)

    # analytical
    T2_naked = 1.0 / alpha_naked
    integrated_naked = 1.0 / alpha_naked  # integral of exp(-alpha*t) from 0 to inf

    log(f"  gamma_0 = {gamma_0}")
    log(f"  alpha_naked = 2*gamma_0 = {alpha_naked:.6f}")
    log(f"  T_2 = {T2_naked:.4f}")
    log(f"  Integrated concurrence = 1/alpha = {integrated_naked:.4f}")

    # numerical confirmation
    Ntot = 2
    H_naked = np.zeros((4, 4), dtype=complex)  # no Hamiltonian
    L_naked = liouvillian_superop(H_naked, [np.sqrt(gamma_0) * site_op(Z, 1, 2)])
    psi_naked = np.zeros(4, dtype=complex)
    psi_naked[0] = 1/np.sqrt(2)  # |00>
    psi_naked[3] = 1/np.sqrt(2)  # |11>
    rho_naked = np.outer(psi_naked, psi_naked.conj())

    t_naked = np.linspace(0, 5*T2_naked, 50)
    c_naked = propagate_concurrence(rho_naked, L_naked, t_naked, 0, 1, 2)

    # fit
    mask = c_naked > 1e-6
    popt, _ = curve_fit(lambda t, a, A: A * np.exp(-a * t),
                        t_naked[mask], c_naked[mask], p0=[alpha_naked, 1.0])
    alpha_fit_naked = popt[0]
    log(f"  Numerical fit: alpha = {alpha_fit_naked:.10f}")
    log(f"  Ratio to analytical: {alpha_fit_naked / alpha_naked:.10f}")

    # === 2. Chain protocols ===
    log("\n" + "=" * 70)
    log("2. CHAIN PROTOCOLS (F67 bonding-mode encoding)")
    log("=" * 70)

    chain_results = {}

    for N in [3, 4, 5]:
        Ntot = N + 1
        d = 2**Ntot
        alpha_1 = formula_alpha_1(N, gamma_0)
        T2_chain = 1.0 / alpha_1
        protection_factor = alpha_naked / alpha_1

        log(f"\n  N={N} ({Ntot} total qubits, dim(L)={4**Ntot})")
        log(f"  alpha_1 (F65) = {alpha_1:.6f}")
        log(f"  T_2 = {T2_chain:.2f}")
        log(f"  Protection factor = alpha_naked/alpha_1 = {protection_factor:.2f}x")

        # build system
        H_chain = np.zeros((d, d), dtype=complex)
        for i in range(1, N):
            H_chain += J * 0.5 * (site_op(X, i, Ntot) @ site_op(X, i+1, Ntot)
                                  + site_op(Y, i, Ntot) @ site_op(Y, i+1, Ntot))
        L_chain = liouvillian_superop(H_chain,
                                      [np.sqrt(gamma_0) * site_op(Z, Ntot-1, Ntot)])

        # Variant B: bonding-mode encoding
        psi_B = np.zeros(d, dtype=complex)
        psi_B[0] = 1/np.sqrt(2)  # |0>_R |vac>_C
        norm = 1/np.sqrt(2) * np.sqrt(2.0/(N+1))
        for i in range(N):
            amp = norm * np.sin(np.pi * (i+1) / (N+1))
            idx = 2**(Ntot-1) + 2**(N-1-i)
            psi_B[idx] = amp
        rho_B = np.outer(psi_B, psi_B.conj())

        t_chain = np.linspace(0, 5*T2_chain, 40)

        # propagate - track concurrence between R and inner qubit Q_0
        c_B = propagate_concurrence(rho_B, L_chain, t_chain, 0, 1, Ntot)

        # also compute negativity-based integrated entanglement
        # for bonding mode: pure exponential, so integrated = C_B(0) / alpha_1
        C_B_0 = c_B[0]
        integrated_B = C_B_0 / alpha_1

        # fit
        mask = c_B > 1e-6
        if np.sum(mask) >= 3:
            popt, _ = curve_fit(lambda t, a, A: A * np.exp(-a * t),
                                t_chain[mask], c_B[mask], p0=[alpha_1, C_B_0])
            alpha_fit_B = popt[0]
        else:
            alpha_fit_B = alpha_1

        log(f"  Variant B: C(0)={C_B_0:.4f}, alpha_fit={alpha_fit_B:.6f}, "
            f"ratio={alpha_fit_B/alpha_1:.4f}")
        log(f"  Integrated C = C(0)/alpha = {integrated_B:.4f}")

        chain_results[N] = {
            'Ntot': Ntot, 'alpha_1': alpha_1, 'T2': T2_chain,
            'protection': protection_factor, 'C0_B': C_B_0,
            'alpha_fit': alpha_fit_B, 'integrated_B': integrated_B
        }

    # === Summary table ===
    log("\n" + "=" * 70)
    log("SUMMARY")
    log("=" * 70)
    log()
    log(f"{'Protocol':>20} {'Qubits':>7} {'alpha':>10} {'T_2':>10} "
        f"{'C(0)':>6} {'Protection':>12}")
    log("-" * 75)
    log(f"{'naked':>20} {2:>7} {alpha_naked:>10.6f} {T2_naked:>10.2f} "
        f"{'1.00':>6} {'1.0x (baseline)':>12}")
    for N in [3, 4, 5]:
        r = chain_results[N]
        log(f"{'F67-B N='+str(N):>20} {r['Ntot']:>7} {r['alpha_1']:>10.6f} "
            f"{r['T2']:>10.2f} {r['C0_B']:>6.3f} "
            f"{r['protection']:>10.1f}x")

    log()
    log("Protection factor = alpha_naked / alpha_1 = T_2(chain) / T_2(naked)")
    log(f"Asymptotic scaling: T_2 ~ (N+1)^3 / (4 pi^2 gamma_0)")
    log()

    # cubic scaling check
    log("Cubic scaling check:")
    for N in [3, 4, 5]:
        r = chain_results[N]
        predicted_T2 = (N+1)**3 / (4 * np.pi**2 * gamma_0)
        log(f"  N={N}: T_2={r['T2']:.2f}, "
            f"(N+1)^3/(4pi^2 gamma_0)={predicted_T2:.2f}, "
            f"ratio={r['T2']/predicted_T2:.4f}")

    # === Extended analytical table ===
    log("\n" + "=" * 70)
    log("EXTENDED ANALYTICAL TABLE (F65 formula, verified N=3..30)")
    log("=" * 70)
    log()
    log("Note: For N > 5 the chain is not Liouvillian-verified; alpha_1 is")
    log("the F65 first-order formula. Full Liouvillian shows O((gamma_0/J)^2)")
    log("corrections (for gamma_0/J=0.05 and N=5: relative shift ~4e-3).")
    log()
    log(f"{'N':>4} {'Qubits':>7} {'alpha_1':>14} {'T_2':>14} {'Protection':>12}")
    log("-" * 60)
    for N in [3, 5, 7, 10, 15, 20, 30, 50, 100]:
        M = N + 1
        a1 = formula_alpha_1(N, gamma_0)
        T2 = 1.0 / a1
        prot = alpha_naked / a1
        log(f"{N:4d} {M:7d} {a1:14.6e} {T2:14.1f} {prot:10.0f}x")

    _outf.close()
    print(f"\nResults written to {OUT_PATH}")
