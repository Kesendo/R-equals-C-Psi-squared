"""tau_max vs the spectral gap: is the extracted formula tau_max = h / sqrt(lambda2(L) * J^2) real?

Spec:     docs/superpowers/specs/2026-05-28-tau-max-spectral-gap-design.md
Write-up: experiments/TAU_MAX_SPECTRAL_GAP.md

Finding: the relaxation timescale is the INVERSE spectral gap,

    tau = 1 / lambda2 = 1 / (2 gamma)

set by gamma alone -- gamma is the timekeeper. At gamma = 0 the gap is 0 and the
clock stops (tau -> infinity, no decay; only Hamiltonian oscillation remains).
The coupling J sets oscillation frequencies, never decay rates (the decay
spectrum is exact multiples of gamma, J-independent). The extracted formula is
wrong twice over: it uses 1/sqrt(lambda2) where the timescale is 1/lambda2
(gamma-power -1/2 instead of -1), and it injects a spurious 1/J (J-power -1
instead of 0).

Self-contained (no cross-drive imports). The Liouvillian construction is
identical to simulations/decay_derivation.py.
"""
import numpy as np

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def op_at(op, qubit, n_q):
    ops = [I2] * n_q
    ops[qubit] = op
    out = ops[0]
    for o in ops[1:]:
        out = np.kron(out, o)
    return out


def heisenberg_chain_H(n_q, J=1.0):
    d = 2 ** n_q
    H = np.zeros((d, d), dtype=complex)
    for i in range(n_q - 1):
        for P in (X, Y, Z):
            H += J * (op_at(P, i, n_q) @ op_at(P, i + 1, n_q))
    return H


def build_L(H, gamma, n_q):
    d = H.shape[0]
    I_d = np.eye(d, dtype=complex)
    L = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for i in range(n_q):
        Lk = np.sqrt(gamma) * op_at(Z, i, n_q)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    return L


def all_rates(L):
    """Full decay-rate multiset (-Re of every eigenvalue), sorted ascending."""
    ev = np.linalg.eigvals(L)
    return np.sort(np.round(-ev.real, 6))


def spectral_gap(L, tol=1e-6):
    """Smallest nonzero decay rate lambda2 (the 2*gamma floor; exactly 0 at gamma=0)."""
    r = all_rates(L)
    pos = r[r > tol]
    return float(pos[0]) if pos.size else 0.0


def clock_tau(lam2):
    """The relaxation timescale = inverse spectral gap, 1/lambda2 = 1/(2 gamma).
    Infinite at gamma=0 (gap 0): no dephasing, no clock."""
    return float("inf") if lam2 <= 0 else 1.0 / lam2


def formula_tau(lam2, J, hbar=1.0):
    """The extracted claim tau_max = hbar / sqrt(lambda2 * J^2)."""
    denom = np.sqrt(lam2 * J ** 2)
    return float("inf") if denom == 0 else hbar / denom


def loglog_slope(xs, ys):
    """Fitted power-law exponent d(log y)/d(log x)."""
    return float(np.polyfit(np.log(xs), np.log(ys), 1)[0])


def measure(n_q, J, gamma):
    lam2 = spectral_gap(build_L(heisenberg_chain_H(n_q, J), gamma, n_q))
    return {"N": n_q, "J": J, "gamma": gamma, "lambda2": lam2,
            "clock": clock_tau(lam2), "formula": formula_tau(lam2, J)}


def run_analysis():
    """Sweep N, J, gamma; fit the power-law exponents of the true clock (1/lambda2)
    and of the extracted formula; return the comparison and verdict.

    The formula matches the clock only if its exponents are (J: 0, gamma: -1).
    """
    g0, J0 = 0.05, 1.0
    N_vals, J_vals, g_vals = (2, 3, 4, 5), (0.5, 1.0, 2.0, 4.0), (0.02, 0.05, 0.10, 0.20)
    N_sweep = [measure(N, J0, g0) for N in N_vals]
    J_sweep = [measure(3, J, g0) for J in J_vals]
    g_sweep = [measure(3, J0, g) for g in g_vals]

    exponents = {
        "clock_vs_J":       loglog_slope(J_vals, [r["clock"] for r in J_sweep]),    # ~ 0
        "formula_vs_J":     loglog_slope(J_vals, [r["formula"] for r in J_sweep]),  # ~ -1
        "clock_vs_gamma":   loglog_slope(g_vals, [r["clock"] for r in g_sweep]),    # ~ -1
        "formula_vs_gamma": loglog_slope(g_vals, [r["formula"] for r in g_sweep]),  # ~ -0.5
    }
    gap_at_zero = spectral_gap(build_L(heisenberg_chain_H(3, J0), 0.0, 3))  # 0 -> clock stops

    tol = 0.05
    formula_correct = (abs(exponents["formula_vs_J"]) < tol
                       and abs(exponents["formula_vs_gamma"] + 1.0) < tol)
    verdict = "confirmed" if formula_correct else "rejected"
    return {"N_sweep": N_sweep, "J_sweep": J_sweep, "gamma_sweep": g_sweep,
            "exponents": exponents, "gap_at_zero": gap_at_zero, "verdict": verdict}


def main(out_path="simulations/results/tau_max_spectral_gap.txt"):
    import os
    a = run_analysis()
    e = a["exponents"]
    out = []
    out.append("tau_max vs the spectral gap")
    out.append("formula under test:  tau_max = h / sqrt(lambda2(L) * J^2)")
    out.append("true clock:          tau = 1 / lambda2 = 1 / (2 gamma)")
    out.append("=" * 72)
    out.append(f"{'sweep':>6} {'N':>2} {'J':>5} {'gamma':>6} {'lambda2':>8} {'clock=1/gap':>12} {'formula':>9}")
    for tag, rows in (("N", a["N_sweep"]), ("J", a["J_sweep"]), ("gamma", a["gamma_sweep"])):
        for r in rows:
            out.append(f"{tag:>6} {r['N']:>2} {r['J']:>5.2f} {r['gamma']:>6.3f} "
                       f"{r['lambda2']:>8.4f} {r['clock']:>12.3f} {r['formula']:>9.3f}")
    out.append("=" * 72)
    out.append("fitted power-law exponents (log-log slope):")
    out.append(f"  clock   ~ J^({e['clock_vs_J']:+.3f})       truth J^0     (J sets frequency, not the rate)")
    out.append(f"  formula ~ J^({e['formula_vs_J']:+.3f})       WRONG         (formula injects a spurious 1/J)")
    out.append(f"  clock   ~ gamma^({e['clock_vs_gamma']:+.3f})   truth gamma^-1 (tau = 1/(2 gamma))")
    out.append(f"  formula ~ gamma^({e['formula_vs_gamma']:+.3f})   WRONG         (formula gives gamma^-1/2)")
    out.append(f"gamma=0:  spectral gap = {a['gap_at_zero']:.4f}  ->  clock stops (tau = infinity)")
    out.append("=" * 72)
    out.append(f"VERDICT: {a['verdict'].upper()}")
    out.append("The timescale is 1/lambda2 = 1/(2 gamma), set by gamma alone (the timekeeper).")
    out.append("The formula's square root and its J are both wrong: 1/sqrt(lambda2*J^2) is")
    out.append("neither the right functional form (1/lambda2) nor J-independent.")
    text = "\n".join(out)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(text + "\n")
    print(text)
    return a


if __name__ == "__main__":
    main()
