"""
Time evolution with Neel initial state |01010> for fair comparison.

The previous simulation used |+>^5, a Heisenberg eigenstate that
does not evolve without noise. The Neel state evolves unitarily,
providing a fair comparison between chain types.
"""

import numpy as np
from scipy.linalg import expm
import json
import csv
from pathlib import Path


# === Physics (reuse from previous scripts) ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


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
    trace_over = sorted([q for q in range(n_qubits) if q not in keep], reverse=True)
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
    return von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j) - von_neumann_entropy(rho_ij)


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
    return {'purity': purity, 'sum_mi': sum(mi_values.values()), **mi_values, **cpsi_values}


def make_neel_state(n_qubits, pattern='01'):
    """Build Neel state |01010...> or |10101...>."""
    d = 2**n_qubits
    idx = 0
    for k in range(n_qubits):
        bit = int(pattern[k % len(pattern)])
        idx += bit * (2**(n_qubits - 1 - k))
    psi = np.zeros(d, dtype=complex)
    psi[idx] = 1.0
    return psi, idx


def run_scenario(label, gammas, times, n_qubits, J, rho0_vec, d, out):
    out(f"\n--- {label} ---")
    out(f"gammas: [{', '.join(f'{g:.4f}' for g in gammas)}], sum={sum(gammas):.4f}")

    H = build_heisenberg_chain(n_qubits, J)
    L = build_liouvillian(H, gammas)

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

    out(f"{'t':>5} {'Purity':>8} {'SumMI':>8} {'MI_01':>7} {'MI_12':>7} {'MI_23':>7} {'MI_34':>7}")
    for r in results:
        out(f"{r['t']:5.1f} {r['purity']:8.4f} {r['sum_mi']:8.4f} "
            f"{r['MI_01']:7.4f} {r['MI_12']:7.4f} {r['MI_23']:7.4f} {r['MI_34']:7.4f}")
    return results


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "time_evolution_neel.txt"
    csv_path = results_dir / "time_evolution_neel_plotdata.csv"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== TIME EVOLUTION: NEEL STATE |01010> ===")
    out("N=5 Chain, J=1.0")
    out()

    N = 5
    J = 1.0
    d = 2**N
    times = np.arange(0, 5.5, 0.5)

    # Build Neel states
    psi_01, idx_01 = make_neel_state(N, '01')
    psi_10, idx_10 = make_neel_state(N, '10')
    out(f"|01010> index: {idx_01} (binary: {idx_01:05b})")
    out(f"|10101> index: {idx_10} (binary: {idx_10:05b})")

    rho0_01 = np.outer(psi_01, psi_01.conj())
    rho0_10 = np.outer(psi_10, psi_10.conj())
    rho0_01_vec = rho0_01.flatten()
    rho0_10_vec = rho0_10.flatten()

    # === VALIDATION: gamma=0 ===
    out(f"\n{'='*70}")
    out("VALIDATION: gamma=0, |01010>")
    out(f"{'='*70}")

    gammas_zero = [0.0] * N
    H = build_heisenberg_chain(N, J)
    L_zero = build_liouvillian(H, gammas_zero)

    rho_vec_05 = expm(L_zero * 0.5) @ rho0_01_vec
    rho_05 = rho_vec_05.reshape(d, d)
    rho_05 = (rho_05 + rho_05.conj().T) / 2
    purity_05 = np.real(np.trace(rho_05 @ rho_05))
    obs_05 = compute_observables(rho_05, N)

    out(f"  t=0.5, gamma=0:")
    out(f"  Purity: {purity_05:.6f} (must be 1.0 for unitary)")
    out(f"  SumMI:  {obs_05['sum_mi']:.6f} (must be > 0 for Neel)")

    if purity_05 < 0.999:
        out("  WARNING: Purity < 0.999. Non-unitary evolution at gamma=0!")
    if obs_05['sum_mi'] < 1e-6:
        out("  FATAL: SumMI ~ 0. Neel state not evolving. STOPPING.")
        return

    out(f"  VALIDATION PASSED: Purity={purity_05:.6f}, SumMI={obs_05['sum_mi']:.6f}")

    # === gamma=0 full evolution for reference ===
    out(f"\n--- gamma=0 reference (|01010>) ---")
    ref_results = run_scenario("gamma=0 (unitary)", gammas_zero, times, N, J, rho0_01_vec, d, out)

    # === Load IBM data ===
    json_path = Path(__file__).parent.parent / "data" / "ibm_sacrifice_zone_march2026" / "sacrifice_zone_hardware_20260324_191713.json"
    with open(json_path) as f:
        hw = json.load(f)
    ci = hw['chain_info']
    sac_T2star = ci['T2star_values']
    sac_T2echo = ci['T2_values']

    csv_data_path = Path(__file__).parent.parent / "data" / "ibm_history" / "ibm_torino_history.csv"
    t2_chain = [18, 89, 19, 90, 60]
    qubit_t2 = {}
    latest_date = None
    with open(csv_data_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            q = int(row['qubit'])
            t2 = float(row['T2_us']) if row['T2_us'] else None
            date = row['date']
            if t2 and t2 > 0:
                if latest_date is None or date > latest_date:
                    latest_date = date
                qubit_t2[(date, q)] = t2
    t2_T2echo = [qubit_t2.get((latest_date, q), 200.0) for q in t2_chain]
    t2_T2star = [t / 2.0 for t in t2_T2echo]

    # === 6 scenarios with |01010> ===
    scenarios = [
        ("1 Mean-T2, No DD", [1/t for t in t2_T2star]),
        ("2 Mean-T2, Uniform DD", [1/t for t in t2_T2echo]),
        ("3 Mean-T2, Selective DD", [1/t2_T2star[0], 1/t2_T2echo[1], 1/t2_T2echo[2], 1/t2_T2echo[3], 1/t2_T2star[4]]),
        ("4 Sacrifice, No DD", [1/t for t in sac_T2star]),
        ("5 Sacrifice, Uniform DD", [1/t for t in sac_T2echo]),
        ("6 Sacrifice, Selective DD", [1/sac_T2star[0], 1/sac_T2echo[1], 1/sac_T2echo[2], 1/sac_T2echo[3], 1/sac_T2echo[4]]),
    ]

    out(f"\n{'='*70}")
    out("6 SCENARIOS with |01010>")
    out(f"{'='*70}")

    all_results = {}
    csv_rows = []

    for label, gammas in scenarios:
        results = run_scenario(label, gammas, times, N, J, rho0_01_vec, d, out)
        all_results[label] = results
        for r in results:
            csv_rows.append({'scenario': label, 'initial': '01010', 't': r['t'],
                           'SumMI': r['sum_mi'], 'Purity': r['purity'],
                           'MI_01': r['MI_01'], 'MI_12': r['MI_12'],
                           'MI_23': r['MI_23'], 'MI_34': r['MI_34']})

    # === |10101> for asymmetry check ===
    out(f"\n{'='*70}")
    out("ASYMMETRY CHECK: |10101> vs |01010> (Scenario 4 and 6)")
    out(f"{'='*70}")

    for label, gammas in [scenarios[3], scenarios[5]]:
        results_10 = run_scenario(f"{label} [|10101>]", gammas, times, N, J, rho0_10_vec, d, out)
        results_01 = all_results[label]

        out(f"\n  |01010> vs |10101> SumMI difference:")
        for r01, r10 in zip(results_01, results_10):
            diff = r01['sum_mi'] - r10['sum_mi']
            out(f"  t={r01['t']:.1f}: {r01['sum_mi']:.4f} vs {r10['sum_mi']:.4f} (diff={diff:+.4f})")

    # === Comparison table ===
    out(f"\n{'='*70}")
    out("SumMI EVOLUTION (all scenarios, |01010>)")
    out(f"{'='*70}")

    headers = [f"{'t':>5}", f"{'g=0':>8}"]
    for label, _ in scenarios:
        short = label.split(",")[0].replace("Mean-T2", "T2").replace("Sacrifice", "Sac")
        headers.append(f"{short:>10}")
    out(" ".join(headers))

    for ti, t in enumerate(times):
        row = [f"{t:5.1f}", f"{ref_results[ti]['sum_mi']:8.4f}"]
        for label, _ in scenarios:
            r = all_results[label][ti]
            row.append(f"{r['sum_mi']:10.4f}")
        out(" ".join(row))

    # === Comparison at t=2.5 ===
    out(f"\n{'='*70}")
    out("COMPARISON AT t=2.5 (|01010>)")
    out(f"{'='*70}")
    out(f"\n{'Scenario':>30} {'SumMI':>8} {'Purity':>8} {'vs Scen1':>9}")

    baseline_mi = None
    for label, results in all_results.items():
        r25 = [r for r in results if abs(r['t'] - 2.5) < 0.01]
        if r25:
            r = r25[0]
            if baseline_mi is None:
                baseline_mi = r['sum_mi']
            ratio = r['sum_mi'] / baseline_mi if baseline_mi > 1e-10 else 0
            out(f"{label:>30} {r['sum_mi']:8.4f} {r['purity']:8.4f} {ratio:9.2f}x")

    # === |+>^5 vs |01010> comparison ===
    out(f"\n{'='*70}")
    out("COMPARISON: |01010> vs |+>^5 at t=1.5 (peak)")
    out(f"{'='*70}")
    out(f"\n{'Scenario':>30} {'Neel SumMI':>12} {'Plus SumMI':>12} {'Neel/Plus':>10}")

    # Read plus results from previous run
    plus_path = results_dir / "time_evolution_6scenarios.txt"
    if plus_path.exists():
        out("  (|+>^5 data from previous run)")
    else:
        out("  (|+>^5 data not available)")

    # Write CSV
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['scenario', 'initial', 't', 'SumMI', 'Purity',
                                               'MI_01', 'MI_12', 'MI_23', 'MI_34'])
        writer.writeheader()
        writer.writerows(csv_rows)

    out(f"\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")
    print(f">>> Plot data saved to: {csv_path}")


if __name__ == "__main__":
    main()
