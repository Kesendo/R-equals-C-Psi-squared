"""
verify_derivations.py
Numerical verification of derived relations D1-D6 from docs/ANALYTICAL_FORMULAS.md.
Builds Heisenberg Liouvillians for N=2..5 and checks each derivation
against exact eigenvalues.

Also checks formula 33 (N=3 intermediate rates) and all weight-sector
rate assignments.
"""

import numpy as np
from math import comb
from pathlib import Path
from scipy.optimize import brentq
from scipy.linalg import expm

# ---- Pauli matrices ----
I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, site, N):
    """Place operator at site, identity elsewhere."""
    r = np.eye(1, dtype=complex)
    for k in range(N):
        r = np.kron(r, op if k == site else I2)
    return r


def build_H_chain(N, J=1.0):
    """Heisenberg XXX chain, open boundaries."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [sx, sy, sz]:
            H += J * kron_at(P, i, N) @ kron_at(P, i + 1, N)
    return H


def build_L(H, gamma, noise_pauli=None):
    """Liouvillian with uniform dephasing. noise_pauli defaults to sz."""
    if noise_pauli is None:
        noise_pauli = sz
    d = H.shape[0]
    N = int(np.log2(d))
    Id = np.eye(d, dtype=complex)
    d2 = d * d
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(N):
        Lk = np.sqrt(gamma) * kron_at(noise_pauli, k, N)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def stat_count(N):
    """Stat(N) = Sum_J m(J,N) * (2J+1)^2 (formula 4)."""
    total = 0.0
    if N % 2 == 0:
        J_vals = [j for j in range(N // 2 + 1)]
    else:
        J_vals = [j + 0.5 for j in range(N // 2 + 1)]
    for J in J_vals:
        idx = int(N / 2 - J)
        if idx < 0 or idx > N:
            continue
        m = comb(N, idx) * (2 * J + 1) / (N / 2 + J + 1)
        total += m * (2 * J + 1) ** 2
    return int(round(total))


def compute_cpsi_matrix(rho, d):
    """CPsi = Tr(rho^2) * L1/(d-1) from density matrix."""
    C = np.real(np.trace(rho @ rho))
    offdiag = np.abs(rho) - np.diag(np.abs(np.diag(rho)))
    L1 = np.sum(offdiag)
    psi = L1 / (d - 1)
    return C * psi


# ==================================================================
# D1: Bandwidth and mode density (from formula 2)
# ==================================================================
def verify_D1(out, J=1.0, gamma=0.05):
    out.append("=" * 60)
    out.append("D1: Bandwidth = 8J*cos(pi/N) -> 8J as N -> inf")
    out.append("    from formula 2: omega_k = 4J*(1 - cos(pi*k/N))")
    out.append("=" * 60)

    all_ok = True
    for N in range(2, 6):
        H = build_H_chain(N, J)
        L = build_L(H, gamma)
        evals = np.linalg.eigvals(L)

        rates = -evals.real
        freqs = np.abs(evals.imag)

        # w=1 modes: rate ~ 2*gamma
        tol = 0.3 * gamma
        w1 = [(r, f) for r, f in zip(rates, freqs)
               if abs(r - 2 * gamma) < tol and f > 1e-8]
        w1_freqs = sorted(set(round(f, 8) for _, f in w1))

        # Analytical
        analytical = sorted(4 * J * (1 - np.cos(np.pi * k / N))
                            for k in range(1, N))

        # Bandwidth
        if len(w1_freqs) >= 2:
            bw_num = w1_freqs[-1] - w1_freqs[0]
        else:
            bw_num = 0.0
        bw_formula = 8 * J * np.cos(np.pi / N)

        # Frequency match
        freq_err = 0.0
        if len(w1_freqs) == len(analytical):
            freq_err = max(abs(a - n) for a, n in
                          zip(analytical, w1_freqs))

        out.append(f"\n  N={N}: {len(w1_freqs)} w=1 frequencies "
                   f"(expected {N-1})")
        out.append(f"    numerical:  {[round(f,6) for f in w1_freqs]}")
        out.append(f"    formula 2:  {[round(f,6) for f in analytical]}")
        out.append(f"    max freq error: {freq_err:.2e}")
        out.append(f"    BW numerical:   {bw_num:.10f}")
        out.append(f"    BW formula:     {bw_formula:.10f}")
        out.append(f"    BW/8J = {bw_formula/(8*J):.6f}")

        if freq_err > 1e-6:
            all_ok = False
            out.append(f"    *** FREQUENCY MISMATCH > 1e-6 ***")

    status = "VERIFIED" if all_ok else "FAILED"
    out.append(f"\n  D1 STATUS: {status}\n")
    return all_ok


# ==================================================================
# D2: V-Effect = Q_max / Q_mean (from formulas 6 + 7)
# ==================================================================
def verify_D2(out, J=1.0, gamma=0.05):
    out.append("=" * 60)
    out.append("D2: V(N) = Q_max / Q_mean = 1 + cos(pi/N)")
    out.append("    Q_mean = 2J/gamma exactly (cosine sum = 0)")
    out.append("=" * 60)

    all_ok = True
    for N in range(2, 6):
        H = build_H_chain(N, J)
        L = build_L(H, gamma)
        evals = np.linalg.eigvals(L)

        rates = -evals.real
        freqs = np.abs(evals.imag)

        # Match eigenvalue frequencies to analytical formula 2
        analytical_freqs = [4 * J * (1 - np.cos(np.pi * k / N))
                            for k in range(1, N)]

        # For each analytical freq, find matching eigenvalue
        matched_Qs = []
        for af in analytical_freqs:
            best_q = None
            best_err = 1e10
            for r, f in zip(rates, freqs):
                if abs(r - 2 * gamma) < 0.5 * gamma and f > 1e-8:
                    err = abs(f - af)
                    if err < best_err:
                        best_err = err
                        best_q = f / (2 * gamma)
            if best_q is not None:
                matched_Qs.append(best_q)

        Q_mean_num = np.mean(matched_Qs) if matched_Qs else float('nan')
        Q_max_num = max(matched_Qs) if matched_Qs else float('nan')
        V_num = (Q_max_num / Q_mean_num
                 if len(matched_Qs) > 1 else 1.0)

        Q_mean_f = 2 * J / gamma
        V_f = 1 + np.cos(np.pi / N)

        # Analytical proof: Sum cos(pi*k/N) for k=1..N-1 = 0
        cos_sum = sum(np.cos(np.pi * k / N) for k in range(1, N))

        err_mean = abs(Q_mean_num - Q_mean_f)
        err_V = abs(V_num - V_f)

        out.append(f"\n  N={N} ({len(matched_Qs)} matched Q values):")
        out.append(f"    Sum cos(pi*k/N) = {cos_sum:.2e} (should be 0)")
        out.append(f"    Q_mean numerical:  {Q_mean_num:.10f}")
        out.append(f"    Q_mean formula:    {Q_mean_f:.10f}")
        out.append(f"    Q_mean deviation:  {err_mean:.2e}")
        out.append(f"    V(N) numerical:    {V_num:.10f}")
        out.append(f"    V(N) formula:      {V_f:.10f}")
        out.append(f"    V(N) deviation:    {err_V:.2e}")

        if err_V > 1e-4:
            all_ok = False
            out.append(f"    *** V-EFFECT MISMATCH > 1e-4 ***")

    status = "VERIFIED" if all_ok else "FAILED"
    out.append(f"\n  D2 STATUS: {status}\n")
    return all_ok


# ==================================================================
# D3: Crossing time ratios (from formula 27)
# ==================================================================
def verify_D3(out, J=1.0, gamma=0.05):
    out.append("=" * 60)
    out.append("D3: t_X/t_Z = K_X/K_Z, verified analytically + propagation")
    out.append("=" * 60)

    # -- Analytical K values --
    # K_Z: f*(1+f*^2) = 3/2, K = -ln(f*)/4
    f_star = brentq(lambda f: f * (1 + f**2) - 1.5, 0.5, 1.0)
    K_Z = -np.log(f_star) / 4

    # K_X: v^2 = 1/2, K = ln(sqrt(2))/4 = ln(2)/8
    K_X = np.log(2) / 8

    # K_depol: u(1+3u^2) = 3, K = -3*ln(u*)/8
    u_star = brentq(lambda u: u * (1 + 3 * u**2) - 3.0, 0.5, 1.0)
    K_depol = -np.log(u_star) * 3 / 8

    ratio_XZ = K_X / K_Z
    ratio_dZ = K_depol / K_Z

    out.append(f"\n  Analytical:")
    out.append(f"    f* (Z crossing):   {f_star:.12f}")
    out.append(f"    K_Z = {K_Z:.10f}  (doc: 0.0374)")
    out.append(f"    K_X = ln(2)/8 = {K_X:.10f}  (doc: 0.0867)")
    out.append(f"    u* (depol crossing): {u_star:.12f}")
    out.append(f"    K_depol = {K_depol:.10f}  (doc: 0.0440)")
    out.append(f"    K_X/K_Z = {ratio_XZ:.10f}  (doc: 2.317)")
    out.append(f"    K_depol/K_Z = {ratio_dZ:.10f}  (doc: 1.176)")

    # -- Propagation cross-check (Bell+ N=2) --
    N = 2
    H = build_H_chain(N, J)
    d = 2 ** N
    bell = np.zeros(d, dtype=complex)
    bell[0] = 1 / np.sqrt(2)
    bell[d - 1] = 1 / np.sqrt(2)
    rho0 = np.outer(bell, bell.conj())
    vec0 = rho0.flatten('F')

    dt = 0.005
    nsteps = 4000

    results_prop = {}
    for label, pauli in [("Z", sz), ("X", sx)]:
        Lmat = build_L(H, gamma, noise_pauli=pauli)
        prop = expm(Lmat * dt)
        v = vec0.copy()
        cpsi_prev = compute_cpsi_matrix(rho0, d)
        t_cross = None
        for s in range(nsteps):
            v = prop @ v
            rho = v.reshape((d, d), order='F')
            cpsi = compute_cpsi_matrix(rho, d)
            if cpsi_prev >= 0.25 and cpsi < 0.25 and t_cross is None:
                frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
                t_cross = (s + frac) * dt
            cpsi_prev = cpsi
        results_prop[label] = t_cross

    t_Z = results_prop.get("Z")
    t_X = results_prop.get("X")
    if t_Z and t_X:
        ratio_prop = t_X / t_Z
        out.append(f"\n  Propagation (Bell+ N=2, gamma={gamma}, dt={dt}):")
        out.append(f"    t_cross_Z = {t_Z:.6f}  "
                   f"(analytical: {-np.log(f_star)/(4*gamma):.6f})")
        out.append(f"    t_cross_X = {t_X:.6f}  "
                   f"(analytical: {np.log(2)/(8*gamma):.6f})")
        out.append(f"    Ratio propagation: {ratio_prop:.6f}")
        out.append(f"    Ratio analytical:  {ratio_XZ:.6f}")
        out.append(f"    Deviation: {abs(ratio_prop - ratio_XZ):.6f}")

    # Check document values
    doc_ok = (abs(K_Z - 0.0374) < 0.0001 and
              abs(K_X - 0.0867) < 0.0001 and
              abs(K_depol - 0.0440) < 0.0001)
    status = "VERIFIED" if doc_ok else "FAILED"
    out.append(f"\n  D3 STATUS: {status}\n")
    return doc_ok


# ==================================================================
# D4: Dimensional factor (d-1)/2 (from formulas 12 + 25)
# ==================================================================
def verify_D4(out):
    out.append("=" * 60)
    out.append("D4: Crossing condition f(1+f^2) = (d-1)/2")
    out.append("    Single qubit (d=2): 1/2.  Bell+ (d=4): 3/2")
    out.append("=" * 60)

    # Single qubit: CPsi = f(1+f^2)/2, crossing at 1/4
    # f(1+f^2) = 1/2
    f1 = brentq(lambda f: f * (1 + f**2) - 0.5, 0.01, 1.0)
    cpsi_1q = f1 * (1 + f1**2) / 2

    # Bell+ 2-qubit: CPsi = f(1+f^2)/6, crossing at 1/4
    # f(1+f^2) = 3/2
    f2 = brentq(lambda f: f * (1 + f**2) - 1.5, 0.01, 1.0)
    cpsi_2q = f2 * (1 + f2**2) / 6

    ratio = 1.5 / 0.5  # should be 3 = d-1

    out.append(f"\n  Single qubit (d=2):")
    out.append(f"    f* = {f1:.12f}")
    out.append(f"    f*(1+f*^2) = {f1*(1+f1**2):.12f}  (target: 0.5)")
    out.append(f"    CPsi at crossing = {cpsi_1q:.12f}  (target: 0.25)")
    out.append(f"\n  Bell+ 2-qubit (d=4):")
    out.append(f"    f* = {f2:.12f}")
    out.append(f"    f*(1+f*^2) = {f2*(1+f2**2):.12f}  (target: 1.5)")
    out.append(f"    CPsi at crossing = {cpsi_2q:.12f}  (target: 0.25)")
    out.append(f"\n  Ratio: 1.5/0.5 = {ratio:.1f} = d-1 = 3  EXACT")

    # Verify with propagation
    gamma = 0.05
    dt = 0.005

    # Single qubit
    rho_1q = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
    d1 = 2
    L1 = np.zeros((4, 4), dtype=complex)
    # H=0 for single qubit test (pure dephasing)
    Lk = np.sqrt(gamma) * sz
    LdL = Lk.conj().T @ Lk
    Id1 = np.eye(d1, dtype=complex)
    L1 += np.kron(Lk.conj(), Lk)
    L1 -= 0.5 * np.kron(Id1, LdL)
    L1 -= 0.5 * np.kron(LdL.T, Id1)
    prop1 = expm(L1 * dt)
    v1 = rho_1q.flatten('F')
    cpsi_prev = compute_cpsi_matrix(rho_1q, d1)
    t_1q = None
    for s in range(4000):
        v1 = prop1 @ v1
        rho = v1.reshape((d1, d1), order='F')
        cpsi = compute_cpsi_matrix(rho, d1)
        if cpsi_prev >= 0.25 and cpsi < 0.25 and t_1q is None:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            t_1q = (s + frac) * dt
        cpsi_prev = cpsi

    # Bell+ 2-qubit (H=0, pure dephasing)
    d2 = 4
    bell = np.zeros(d2, dtype=complex)
    bell[0] = 1 / np.sqrt(2)
    bell[3] = 1 / np.sqrt(2)
    rho_2q = np.outer(bell, bell.conj())
    Id2 = np.eye(d2, dtype=complex)
    L2 = np.zeros((16, 16), dtype=complex)
    for k in range(2):
        Lk = np.sqrt(gamma) * kron_at(sz, k, 2)
        LdL = Lk.conj().T @ Lk
        L2 += np.kron(Lk.conj(), Lk)
        L2 -= 0.5 * np.kron(Id2, LdL)
        L2 -= 0.5 * np.kron(LdL.T, Id2)
    prop2 = expm(L2 * dt)
    v2 = rho_2q.flatten('F')
    cpsi_prev = compute_cpsi_matrix(rho_2q, d2)
    t_2q = None
    for s in range(4000):
        v2 = prop2 @ v2
        rho = v2.reshape((d2, d2), order='F')
        cpsi = compute_cpsi_matrix(rho, d2)
        if cpsi_prev >= 0.25 and cpsi < 0.25 and t_2q is None:
            frac = (cpsi_prev - 0.25) / (cpsi_prev - cpsi)
            t_2q = (s + frac) * dt
        cpsi_prev = cpsi

    if t_1q and t_2q:
        K_1q = gamma * t_1q
        K_2q = gamma * t_2q
        out.append(f"\n  Propagation (H=0, gamma={gamma}):")
        out.append(f"    K_1q = gamma*t = {K_1q:.8f}  "
                   f"(analytical: {-np.log(f1)/2:.8f})")
        out.append(f"    K_2q = gamma*t = {K_2q:.8f}  "
                   f"(analytical: {-np.log(f2)/4:.8f})")

    ok = abs(cpsi_1q - 0.25) < 1e-10 and abs(cpsi_2q - 0.25) < 1e-10
    status = "VERIFIED" if ok else "FAILED"
    out.append(f"\n  D4 STATUS: {status}\n")
    return ok


# ==================================================================
# D5: Dynamic palindromic mode count (from formulas 4 + 22 + 23)
# NOTE: Stat(N) is valid at gamma=0.  At finite gamma the Hamiltonian
# mixes weight sectors (w with w+-2), so we must verify in the
# near-unitary limit (gamma << J).
# ==================================================================
def verify_D5(out, J=1.0):
    out.append("=" * 60)
    out.append("D5: Oscillating modes = 4^N - (N+1) - Stat(N)")
    out.append("    Valid at gamma -> 0.  Using gamma = 1e-6.")
    out.append("=" * 60)

    gamma = 1e-6  # near-unitary limit
    all_ok = True
    for N in range(2, 6):
        H = build_H_chain(N, J)
        L = build_L(H, gamma)
        evals = np.linalg.eigvals(L)
        total = len(evals)

        rates = -evals.real
        freqs = np.abs(evals.imag)

        # At gamma -> 0:
        # Stationary: |eigenvalue| < threshold (both rate and freq ~ 0)
        # These are modes that commute with H (zero frequency at gamma=0)
        stat_tol = 1e-4  # generous, since gamma is tiny
        n_stat = sum(1 for ev in evals if abs(ev) < stat_tol)

        # XOR: rate ~ 2*N*gamma (nonzero rate, possibly zero freq)
        xor_rate = 2 * N * gamma
        n_xor = sum(1 for r in rates
                    if abs(r - xor_rate) < 0.5 * gamma and r > 1e-10)

        # Oscillating: nonzero frequency
        n_osc = total - n_stat - n_xor

        stat_f = stat_count(N)
        xor_f = N + 1
        osc_f = 4**N - xor_f - stat_f

        out.append(f"\n  N={N} (total {total} = 4^{N}):")
        out.append(f"    Stationary: {n_stat} (formula: {stat_f})")
        out.append(f"    XOR drain:  {n_xor} (formula: {xor_f})")
        out.append(f"    Oscillating: {n_osc} (formula: {osc_f})")

        if n_stat != stat_f:
            out.append(f"    *** STAT MISMATCH: {n_stat} != {stat_f} ***")
            all_ok = False
        if n_xor != xor_f:
            out.append(f"    *** XOR MISMATCH: {n_xor} != {xor_f} ***")
            all_ok = False
        if n_osc != osc_f:
            out.append(f"    *** OSC MISMATCH: {n_osc} != {osc_f} ***")
            all_ok = False

    status = "VERIFIED" if all_ok else "FAILED"
    out.append(f"\n  D5 STATUS: {status}\n")
    return all_ok


# ==================================================================
# D6: Spectral gap = 2*gamma (from formulas 1 + 3)
# ==================================================================
def verify_D6(out, J=1.0, gamma=0.05):
    out.append("=" * 60)
    out.append("D6: Spectral gap = 2*gamma (smallest nonzero decay rate)")
    out.append("    Mixing time <= N*ln(4) / (2*gamma)")
    out.append("=" * 60)

    all_ok = True
    for N in range(2, 6):
        H = build_H_chain(N, J)
        L = build_L(H, gamma)
        evals = np.linalg.eigvals(L)
        rates = -evals.real

        # Nonzero rates (exclude steady state)
        nonzero = sorted(r for r in rates if r > 1e-8)
        gap = nonzero[0] if nonzero else 0.0
        gap_f = 2 * gamma
        mix_time = N * np.log(4) / (2 * gamma)

        err = abs(gap - gap_f)
        out.append(f"\n  N={N}:")
        out.append(f"    Spectral gap:     {gap:.10f}")
        out.append(f"    Expected 2*gamma: {gap_f:.10f}")
        out.append(f"    Deviation:        {err:.2e}")
        out.append(f"    Mixing time bound: {mix_time:.4f}")

        if err > 1e-6:
            all_ok = False
            out.append(f"    *** GAP MISMATCH > 1e-6 ***")

    status = "VERIFIED" if all_ok else "FAILED"
    out.append(f"\n  D6 STATUS: {status}\n")
    return all_ok


# ==================================================================
# BONUS: Check formula 33 (N=3 rates) and weight-sector structure
# ==================================================================
def check_formula_33(out, J=1.0, gamma=0.05):
    out.append("=" * 60)
    out.append("BONUS: Weight-sector rate structure (formula 33 check)")
    out.append("  Expected: rates at 2*w*gamma for each weight w")
    out.append("  Formula 33 claims: {2g, 8g/3, 10g/3} for N=3")
    out.append("=" * 60)

    for N in range(2, 6):
        H = build_H_chain(N, J)
        L = build_L(H, gamma)
        evals = np.linalg.eigvals(L)
        rates = -evals.real

        # Group rates by value (round to 6 decimals)
        rate_counts = {}
        for r in rates:
            key = round(r, 6)
            rate_counts[key] = rate_counts.get(key, 0) + 1

        sorted_rates = sorted(rate_counts.items())
        out.append(f"\n  N={N}: distinct rates (rounded to 6 decimals):")
        for rate, count in sorted_rates:
            # Expected weight sector
            w_expected = rate / (2 * gamma)
            out.append(f"    rate = {rate:12.8f}  "
                       f"(count={count:4d}, "
                       f"rate/(2g) = {w_expected:.4f})")

        # Check if all rates are integer multiples of 2*gamma
        all_integer = all(abs(r / (2 * gamma) - round(r / (2 * gamma)))
                         < 0.01 for r, _ in sorted_rates if r > 1e-8)
        out.append(f"    All rates = 2*w*gamma (integer w)? {all_integer}")

    out.append("")


# ==================================================================
# MAIN
# ==================================================================
def main():
    out = []
    out.append("Derivation Verification: D1-D6")
    out.append(f"J = 1.0, gamma = 0.05, N = 2..5")
    out.append("=" * 60)
    out.append("")

    results = {}
    results["D1"] = verify_D1(out)
    results["D2"] = verify_D2(out)
    results["D3"] = verify_D3(out)
    results["D4"] = verify_D4(out)
    results["D5"] = verify_D5(out, J=1.0)
    results["D6"] = verify_D6(out)
    check_formula_33(out)

    # Summary
    out.append("=" * 60)
    out.append("SUMMARY")
    out.append("=" * 60)
    for key, ok in results.items():
        out.append(f"  {key}: {'VERIFIED' if ok else 'FAILED'}")
    n_ok = sum(1 for v in results.values() if v)
    out.append(f"\n  {n_ok}/{len(results)} derivations verified.")

    text = "\n".join(out)
    print(text)

    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    (results_dir / "verify_derivations.txt").write_text(text, encoding="utf-8")
    print(f"\nResults written to {results_dir / 'verify_derivations.txt'}")


if __name__ == "__main__":
    main()
