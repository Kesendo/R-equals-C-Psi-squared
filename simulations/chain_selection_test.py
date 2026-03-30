"""
Chain Selection Test: Sacrifice vs Mean-T2 (No DD).

Compares two IBM Torino chains using real calibration data:
  Chain A (Sacrifice-Top): [80, 8, 79, 53, 85] -- score 2.86x
  Chain B (Mean-T2-Top):   [18, 89, 19, 90, 60] -- score 1.06x

Three analyses:
  1. Lindblad time evolution (|+>^5 and |01010>) -- SumMI(t)
  2. Liouvillian spectral analysis -- eigenvalues, palindrome, protection
  3. Cavity mode localization -- eigenvector decomposition

Convention: gamma_k = 1/T2echo_k (same as sacrifice_zone_mapping.py).
Under free decoherence, effective gamma would be ~2x larger (T2* ~ T2echo/2).
The chain COMPARISON is convention-independent (same J=1.0 for both).
"""

import numpy as np
from scipy.linalg import expm, eig
from scipy.stats import pearsonr
import csv
import itertools
from pathlib import Path


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}


def kron_at(op, target, n_qubits):
    result = np.eye(1, dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == target else I2)
    return result


def build_heisenberg_chain(n, J=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for pauli in [X, Y, Z]:
            H += J * kron_at(pauli, i, n) @ kron_at(pauli, i + 1, n)
    return H


def build_liouvillian(H, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n_qubits):
        if gammas[k] <= 0:
            continue
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n_qubits)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def partial_trace(rho, keep, n_qubits):
    dims = [2] * n_qubits
    rho_tensor = rho.reshape(dims + dims)
    trace_over = sorted([q for q in range(n_qubits) if q not in keep],
                        reverse=True)
    nq = n_qubits
    for q in trace_over:
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + nq)
        nq -= 1
    k = len(keep)
    return rho_tensor.reshape(2**k, 2**k)


def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


def mutual_information(rho, qi, qj, n_qubits):
    rho_ij = partial_trace(rho, [qi, qj], n_qubits)
    rho_i = partial_trace(rho, [qi], n_qubits)
    rho_j = partial_trace(rho, [qj], n_qubits)
    return (von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j)
            - von_neumann_entropy(rho_ij))


def compute_observables(rho, n_qubits):
    purity = np.real(np.trace(rho @ rho))
    pairs = [(i, i+1) for i in range(n_qubits - 1)]
    mi_values = {}
    cpsi_values = {}
    for qi, qj in pairs:
        mi_values[f"MI_{qi}{qj}"] = mutual_information(rho, qi, qj, n_qubits)
        rho_pair = partial_trace(rho, [qi, qj], n_qubits)
        pair_purity = np.real(np.trace(rho_pair @ rho_pair))
        l1 = np.sum(np.abs(rho_pair)) - np.sum(np.abs(np.diag(rho_pair)))
        cpsi_values[f"CPsi_{qi}{qj}"] = pair_purity * l1 / 3.0
    return {'purity': purity, 'sum_mi': sum(mi_values.values()),
            **mi_values, **cpsi_values}


def make_plus_state(n_qubits):
    """Build |+>^N state."""
    d = 2**n_qubits
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus
    for _ in range(n_qubits - 1):
        psi = np.kron(psi, plus)
    return np.outer(psi, psi.conj())


def make_neel_state(n_qubits, pattern='01'):
    """Build Neel state |01010...>."""
    d = 2**n_qubits
    idx = 0
    for k in range(n_qubits):
        bit = int(pattern[k % len(pattern)])
        idx += bit * (2**(n_qubits - 1 - k))
    psi = np.zeros(d, dtype=complex)
    psi[idx] = 1.0
    return np.outer(psi, psi.conj())


def time_evolution(L, rho0, times, n_qubits):
    """Run Lindblad time evolution and compute observables."""
    d = 2**n_qubits
    rho0_vec = rho0.flatten()
    results = []
    for t in times:
        if t == 0:
            rho_vec = rho0_vec.copy()
        else:
            rho_vec = expm(L * t) @ rho0_vec
        rho = rho_vec.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        obs = compute_observables(rho, n_qubits)
        obs['t'] = t
        results.append(obs)
    return results


def spectral_analysis(gammas, J=1.0):
    """Eigenvalue analysis: palindrome, protection, mode count."""
    N = len(gammas)
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    sum_g = sum(gammas)
    center = sum_g

    # Palindrome check
    rates = sorted(-ev.real for ev in evals if abs(ev.imag) > 1e-10)
    paired = 0
    used = [False] * len(rates)
    for i in range(len(rates)):
        if used[i]:
            continue
        partner = 2 * center - rates[i]
        for j in range(len(rates)):
            if not used[j] and j != i and abs(rates[j] - partner) < 1e-4:
                paired += 2
                used[i] = True
                used[j] = True
                break
    pal_score = paired / max(len(rates), 1)

    # Mode classification
    osc = [(-ev.real, abs(ev.imag)) for ev in evals if abs(ev.imag) > 1e-10]
    stat = sum(1 for ev in evals
               if abs(ev.real) < 1e-10 and abs(ev.imag) < 1e-10)

    slowest_rate = min(r for r, _ in osc) if osc else float('inf')
    fastest_rate = max(r for r, _ in osc) if osc else 0

    freqs = sorted(set(round(f, 4) for _, f in osc if f > 1e-8))

    # Protection vs uniform
    uniform_gammas = np.array([sum_g / N] * N)
    L_uni = build_liouvillian(H, uniform_gammas)
    evals_uni = np.linalg.eigvals(L_uni)
    osc_uni = [(-ev.real, abs(ev.imag)) for ev in evals_uni
               if abs(ev.imag) > 1e-10]
    slowest_uni = min(r for r, _ in osc_uni) if osc_uni else float('inf')
    protection = slowest_uni / slowest_rate if slowest_rate > 1e-15 else 0

    return {
        'palindrome': pal_score,
        'slowest_rate': slowest_rate,
        'fastest_rate': fastest_rate,
        'slowest_uni': slowest_uni,
        'protection': protection,
        'n_stationary': stat,
        'n_oscillating': len(osc),
        'n_frequencies': len(freqs),
        'sum_gamma': sum_g,
        'max_rate_ratio': fastest_rate / (2 * center) if center > 0 else 0,
        'evals': evals,
    }


def eigenvector_localization(gammas, J=1.0):
    """Eigenvector decomposition for mode localization profile."""
    N = len(gammas)
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    d = 2**N

    evals, evecs = eig(L)

    pauli_labels = ['I', 'X', 'Y', 'Z']
    pauli_strings = []
    for combo in itertools.product(pauli_labels, repeat=N):
        name = ''.join(combo)
        active = [k for k in range(N) if combo[k] != 'I']
        mat = np.eye(1, dtype=complex)
        for c in combo:
            mat = np.kron(mat, PAULIS[c])
        pauli_strings.append((name, active, mat))

    # Compute weights for each oscillating mode
    mode_data = []
    for i in range(len(evals)):
        ev = evals[i]
        rate = -ev.real
        freq = abs(ev.imag)
        if abs(ev.real) < 1e-10 and freq < 1e-10:
            continue  # skip stationary

        V = evecs[:, i].reshape(d, d)
        weights = np.zeros(N)
        total_w = 0.0
        for name, active, P in pauli_strings:
            if not active:
                continue
            c = np.trace(P.conj().T @ V) / d
            w = abs(c)**2
            total_w += w
            for k in active:
                weights[k] += w
        if total_w > 1e-15:
            weights /= total_w

        mode_data.append((rate, freq, weights))

    mode_data.sort(key=lambda x: x[0])  # slowest first
    return mode_data


def load_t2_data(csv_path, chain_qubits):
    """Load T2echo from the latest date in the CSV."""
    qubit_data = {}
    dates = set()
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = int(row['qubit'])
            if q in chain_qubits:
                date = row['date']
                t2 = float(row['T2_us']) if row['T2_us'] else None
                t1 = float(row['T1_us']) if row['T1_us'] else None
                dates.add(date)
                if t2 and t2 > 0 and t1 and t1 > 0:
                    qubit_data[(date, q)] = {'T2': t2, 'T1': t1}

    latest = max(dates)
    result = {}
    for q in chain_qubits:
        d = qubit_data.get((latest, q))
        if d:
            result[q] = d
        else:
            # Fallback to most recent date with data
            for dd in sorted(dates, reverse=True):
                d2 = qubit_data.get((dd, q))
                if d2:
                    result[q] = d2
                    break
    return result, latest


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "chain_selection_test.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("CHAIN SELECTION TEST: Sacrifice vs Mean-T2 (No DD)")
    out("=" * 70)

    N = 5
    J = 1.0

    chain_A = [80, 8, 79, 53, 85]   # Sacrifice-Top, score 2.86x
    chain_B = [18, 89, 19, 90, 60]  # Mean-T2-Top, score 1.06x

    # ================================================================
    # Step 1: Load IBM calibration data
    # ================================================================
    out("\n--- STEP 1: Load IBM calibration data ---")
    csv_path = (Path(__file__).parent.parent / "data" / "ibm_history"
                / "ibm_torino_history.csv")

    all_qubits = set(chain_A + chain_B)
    t2_data, latest_date = load_t2_data(csv_path, all_qubits)
    out(f"Date: {latest_date}")

    for label, chain in [("A (Sacrifice)", chain_A),
                         ("B (Mean-T2)", chain_B)]:
        out(f"\nChain {label}: {chain}")
        out(f"  {'Qubit':>6} {'T1(us)':>8} {'T2(us)':>8} {'r=T2/2T1':>10}"
            f" {'gamma':>10}")
        for q in chain:
            d = t2_data[q]
            r = d['T2'] / (2 * d['T1'])
            gamma = 1.0 / d['T2']
            out(f"  Q{q:>4d} {d['T1']:8.1f} {d['T2']:8.1f} {r:10.4f}"
                f" {gamma:10.6f}")

    gammas_A = np.array([1.0 / t2_data[q]['T2'] for q in chain_A])
    gammas_B = np.array([1.0 / t2_data[q]['T2'] for q in chain_B])

    out(f"\ngamma_A: [{', '.join(f'{g:.6f}' for g in gammas_A)}]"
        f"  sum={sum(gammas_A):.6f}")
    out(f"gamma_B: [{', '.join(f'{g:.6f}' for g in gammas_B)}]"
        f"  sum={sum(gammas_B):.6f}")

    contrast_A = max(gammas_A) / min(gammas_A)
    contrast_B = max(gammas_B) / min(gammas_B)
    out(f"\nContrast (max/min gamma): A={contrast_A:.1f}x  B={contrast_B:.1f}x")

    # ================================================================
    # Step 2: Spectral analysis
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 2: Spectral analysis (Liouvillian eigenvalues)")
    out("=" * 70)

    spec_A = spectral_analysis(gammas_A, J)
    spec_B = spectral_analysis(gammas_B, J)

    out(f"\n{'Metric':>25} {'Chain A':>12} {'Chain B':>12} {'A/B':>8}")
    out(f"{'-'*25} {'-'*12} {'-'*12} {'-'*8}")

    metrics = [
        ('Palindrome score', f"{spec_A['palindrome']:.0%}",
         f"{spec_B['palindrome']:.0%}", ""),
        ('Sum(gamma)', f"{spec_A['sum_gamma']:.6f}",
         f"{spec_B['sum_gamma']:.6f}",
         f"{spec_A['sum_gamma']/spec_B['sum_gamma']:.2f}x"),
        ('Stationary modes', f"{spec_A['n_stationary']}",
         f"{spec_B['n_stationary']}", ""),
        ('Oscillating modes', f"{spec_A['n_oscillating']}",
         f"{spec_B['n_oscillating']}", ""),
        ('Distinct frequencies', f"{spec_A['n_frequencies']}",
         f"{spec_B['n_frequencies']}", ""),
        ('Slowest osc. rate', f"{spec_A['slowest_rate']:.6f}",
         f"{spec_B['slowest_rate']:.6f}",
         f"{spec_B['slowest_rate']/spec_A['slowest_rate']:.2f}x"),
        ('Fastest osc. rate', f"{spec_A['fastest_rate']:.6f}",
         f"{spec_B['fastest_rate']:.6f}", ""),
        ('Max/2*sum_gamma', f"{spec_A['max_rate_ratio']:.4f}",
         f"{spec_B['max_rate_ratio']:.4f}", ""),
        ('Slowest (uniform)', f"{spec_A['slowest_uni']:.6f}",
         f"{spec_B['slowest_uni']:.6f}", ""),
        ('Protection factor', f"{spec_A['protection']:.2f}x",
         f"{spec_B['protection']:.2f}x",
         f"{spec_A['protection']/spec_B['protection']:.2f}x"),
    ]

    for name, va, vb, ratio in metrics:
        out(f"{name:>25} {va:>12} {vb:>12} {ratio:>8}")

    # ================================================================
    # Step 3: Time evolution
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 3: Time evolution (SumMI)")
    out("=" * 70)

    times = np.linspace(0, 10, 51)  # 0 to 10, 51 points

    H = build_heisenberg_chain(N, J)
    L_A = build_liouvillian(H, gammas_A)
    L_B = build_liouvillian(H, gammas_B)

    rho_plus = make_plus_state(N)
    rho_neel = make_neel_state(N, '01')

    out("\nRunning 4 time evolutions (2 chains x 2 initial states)...")
    out("  This may take a few minutes (50 matrix exponentials each).")

    res_A_plus = time_evolution(L_A, rho_plus, times, N)
    out("  Chain A, |+>^5: done")
    res_B_plus = time_evolution(L_B, rho_plus, times, N)
    out("  Chain B, |+>^5: done")
    res_A_neel = time_evolution(L_A, rho_neel, times, N)
    out("  Chain A, |01010>: done")
    res_B_neel = time_evolution(L_B, rho_neel, times, N)
    out("  Chain B, |01010>: done")

    # --- SumMI table: |+>^5 ---
    out(f"\n--- SumMI(t): |+>^5 ---")
    out(f"{'t':>6} {'A_SumMI':>10} {'B_SumMI':>10} {'A/B':>8}"
        f" {'A_Pur':>8} {'B_Pur':>8}")
    for rA, rB in zip(res_A_plus, res_B_plus):
        t = rA['t']
        ratio = (rA['sum_mi'] / rB['sum_mi']
                 if rB['sum_mi'] > 1e-8 else float('inf'))
        if t % 1.0 < 0.01 or t < 0.5:
            out(f"{t:6.2f} {rA['sum_mi']:10.6f} {rB['sum_mi']:10.6f}"
                f" {ratio:8.2f}x {rA['purity']:8.4f} {rB['purity']:8.4f}")

    # --- SumMI table: |01010> ---
    out(f"\n--- SumMI(t): |01010> ---")
    out(f"{'t':>6} {'A_SumMI':>10} {'B_SumMI':>10} {'A/B':>8}"
        f" {'A_Pur':>8} {'B_Pur':>8}")
    for rA, rB in zip(res_A_neel, res_B_neel):
        t = rA['t']
        ratio = (rA['sum_mi'] / rB['sum_mi']
                 if rB['sum_mi'] > 1e-8 else float('inf'))
        if t % 1.0 < 0.01 or t < 0.5:
            out(f"{t:6.2f} {rA['sum_mi']:10.6f} {rB['sum_mi']:10.6f}"
                f" {ratio:8.2f}x {rA['purity']:8.4f} {rB['purity']:8.4f}")

    # --- Per-pair MI at peak time ---
    # Find peak SumMI time for each
    peak_A_plus = max(res_A_plus, key=lambda r: r['sum_mi'])
    peak_B_plus = max(res_B_plus, key=lambda r: r['sum_mi'])
    peak_A_neel = max(res_A_neel, key=lambda r: r['sum_mi'])
    peak_B_neel = max(res_B_neel, key=lambda r: r['sum_mi'])

    out(f"\n--- Peak SumMI ---")
    out(f"  Chain A, |+>^5:  SumMI={peak_A_plus['sum_mi']:.6f}"
        f" at t={peak_A_plus['t']:.2f}")
    out(f"  Chain B, |+>^5:  SumMI={peak_B_plus['sum_mi']:.6f}"
        f" at t={peak_B_plus['t']:.2f}")
    out(f"  Chain A, |01010>: SumMI={peak_A_neel['sum_mi']:.6f}"
        f" at t={peak_A_neel['t']:.2f}")
    out(f"  Chain B, |01010>: SumMI={peak_B_neel['sum_mi']:.6f}"
        f" at t={peak_B_neel['t']:.2f}")

    # Per-pair MI at peak for |01010>
    out(f"\n--- Per-pair MI at peak (|01010>) ---")
    out(f"{'Pair':>6} {'A (t={:.1f})'.format(peak_A_neel['t']):>12}"
        f" {'B (t={:.1f})'.format(peak_B_neel['t']):>12} {'A/B':>8}")
    for pair in ['MI_01', 'MI_12', 'MI_23', 'MI_34']:
        va = peak_A_neel[pair]
        vb = peak_B_neel[pair]
        ratio = va / vb if vb > 1e-8 else float('inf')
        out(f"{pair:>6} {va:12.6f} {vb:12.6f} {ratio:8.2f}x")

    # ================================================================
    # Step 4: Cavity mode localization
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 4: Cavity mode localization (eigenvector analysis)")
    out("=" * 70)

    out("\nComputing eigenvectors (this takes ~1 min)...")
    modes_A = eigenvector_localization(gammas_A, J)
    modes_B = eigenvector_localization(gammas_B, J)
    out(f"  Chain A: {len(modes_A)} oscillating modes")
    out(f"  Chain B: {len(modes_B)} oscillating modes")

    # Slowest 10 modes
    for label, modes, chain in [("A (Sacrifice)", modes_A, chain_A),
                                ("B (Mean-T2)", modes_B, chain_B)]:
        out(f"\nChain {label}: Top 10 slowest modes")
        out(f"{'#':>3} {'Rate':>10} {'Freq':>8}"
            + "".join(f" {'Q'+str(q):>7}" for q in chain))
        for i, (rate, freq, w) in enumerate(modes[:10]):
            ws = "".join(f" {w[k]:7.3f}" for k in range(5))
            out(f"{i+1:3d} {rate:10.6f} {freq:8.3f}{ws}")

    # Correlation: edge qubit weight vs decay rate
    for label, modes in [("A (Sacrifice)", modes_A),
                         ("B (Mean-T2)", modes_B)]:
        if len(modes) > 20:
            # Use the end with the highest gamma as "sacrifice" qubit
            w_edge = [max(w[0], w[4]) for _, _, w in modes]
            rate_vals = [r for r, _, _ in modes]
            r_p, p_p = pearsonr(w_edge, rate_vals)
            out(f"\nChain {label}: Correlation(edge_weight, rate)"
                f" = {r_p:.4f} (p={p_p:.2e})")

    # Slowest vs fastest quartile profile
    for label, modes in [("A (Sacrifice)", modes_A),
                         ("B (Mean-T2)", modes_B)]:
        nq = max(1, len(modes) // 4)
        slow_q = np.mean([w for _, _, w in modes[:nq]], axis=0)
        fast_q = np.mean([w for _, _, w in modes[-nq:]], axis=0)
        out(f"\nChain {label}:")
        out(f"  Slowest quartile profile: "
            f"[{', '.join(f'{v:.3f}' for v in slow_q)}]")
        out(f"  Fastest quartile profile: "
            f"[{', '.join(f'{v:.3f}' for v in fast_q)}]")
        out(f"  Center-to-edge ratio (slow): "
            f"{slow_q[2]/max(slow_q[0], slow_q[4]):.2f}")
        out(f"  Center-to-edge ratio (fast): "
            f"{fast_q[2]/max(fast_q[0], fast_q[4]):.2f}")

    # ================================================================
    # Step 5: Summary and interpretation
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 5: Summary")
    out("=" * 70)

    out(f"\n{'Metric':>30} {'Chain A':>12} {'Chain B':>12} {'Winner':>8}")
    out(f"{'-'*30} {'-'*12} {'-'*12} {'-'*8}")

    summary = [
        ('Gamma contrast (max/min)',
         f"{contrast_A:.1f}x", f"{contrast_B:.1f}x",
         'A' if contrast_A > contrast_B else 'B'),
        ('Sum(gamma)',
         f"{sum(gammas_A):.4f}", f"{sum(gammas_B):.4f}",
         'B' if sum(gammas_B) < sum(gammas_A) else 'A'),
        ('Protection factor',
         f"{spec_A['protection']:.2f}x", f"{spec_B['protection']:.2f}x",
         'A' if spec_A['protection'] > spec_B['protection'] else 'B'),
        ('Palindrome',
         f"{spec_A['palindrome']:.0%}", f"{spec_B['palindrome']:.0%}",
         'tie' if abs(spec_A['palindrome'] - spec_B['palindrome']) < 0.01
         else ('A' if spec_A['palindrome'] > spec_B['palindrome'] else 'B')),
        ('Peak SumMI (|+>^5)',
         f"{peak_A_plus['sum_mi']:.6f}", f"{peak_B_plus['sum_mi']:.6f}",
         'A' if peak_A_plus['sum_mi'] > peak_B_plus['sum_mi'] else 'B'),
        ('Peak SumMI (|01010>)',
         f"{peak_A_neel['sum_mi']:.6f}", f"{peak_B_neel['sum_mi']:.6f}",
         'A' if peak_A_neel['sum_mi'] > peak_B_neel['sum_mi'] else 'B'),
    ]

    for name, va, vb, winner in summary:
        out(f"{name:>30} {va:>12} {vb:>12} {winner:>8}")

    # Interpretation
    out("\n--- INTERPRETATION ---")

    plus_ratio = (peak_A_plus['sum_mi'] / peak_B_plus['sum_mi']
                  if peak_B_plus['sum_mi'] > 1e-8 else float('inf'))
    neel_ratio = (peak_A_neel['sum_mi'] / peak_B_neel['sum_mi']
                  if peak_B_neel['sum_mi'] > 1e-8 else float('inf'))

    out(f"\nPeak SumMI ratio A/B:")
    out(f"  |+>^5:   {plus_ratio:.2f}x")
    out(f"  |01010>: {neel_ratio:.2f}x")

    if spec_A['protection'] > spec_B['protection'] * 1.5:
        out("\nSacrifice chain has STRONGER mode protection.")
        out("The gamma contrast creates a localization advantage.")
    elif spec_B['protection'] > spec_A['protection'] * 1.5:
        out("\nMean-T2 chain has STRONGER mode protection.")
        out("Low uniform noise beats high-contrast noise.")
    else:
        out("\nBoth chains have SIMILAR mode protection.")

    if plus_ratio > 2:
        out(f"\nFor |+>^5: Sacrifice chain wins by {plus_ratio:.1f}x.")
        out("Noise drives dynamics -- more noise contrast = more signal.")
    elif plus_ratio < 0.5:
        out(f"\nFor |+>^5: Mean-T2 chain wins by {1/plus_ratio:.1f}x.")
    else:
        out(f"\nFor |+>^5: Similar performance (ratio {plus_ratio:.2f}x).")

    if neel_ratio > 1.5:
        out(f"\nFor |01010>: Sacrifice chain wins by {neel_ratio:.1f}x.")
        out("Mode protection advantage persists under Hamiltonian dynamics.")
    elif neel_ratio < 0.67:
        out(f"\nFor |01010>: Mean-T2 chain wins by {1/neel_ratio:.1f}x.")
        out("Quieter qubits preserve coherence better when Hamiltonian drives.")
    else:
        out(f"\nFor |01010>: Similar performance (ratio {neel_ratio:.2f}x).")

    out(f"\nNote: Chain A has {sum(gammas_A)/sum(gammas_B):.1f}x more total"
        f" noise than Chain B.")
    out("Any SumMI advantage for A comes DESPITE having more noise,")
    out("not because of having less noise. This would confirm mode protection.")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
