"""
Time evolution for 6 optimization scenarios.

Simulates rho(t) = expm(L*t) @ rho(0) for all 6 combinations of
chain selection x DD strategy. Extracts SumMI, purity, CPsi at
each time point. Produces the curves that will be measured on IBM.
"""

import numpy as np
from scipy.linalg import eig, expm
import json
import csv
from pathlib import Path


# === Physics ===
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
    """Trace out all qubits except those in 'keep' list."""
    d = 2**n_qubits
    dims = [2] * n_qubits
    rho_tensor = rho.reshape(dims + dims)

    # Determine which axes to trace over
    trace_over = [q for q in range(n_qubits) if q not in keep]
    # Sort in reverse order to trace from highest index first
    for q in sorted(trace_over, reverse=True):
        # Trace over axis q (row) and axis q+n_qubits (col)
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + n_qubits)
        n_qubits -= 1  # one fewer qubit now
        # Adjust: after tracing q, remaining qubit indices shift

    # Result shape should be (2^k, 2^k) for k = len(keep)
    k = len(keep)
    return rho_tensor.reshape(2**k, 2**k)


def von_neumann_entropy(rho):
    """S(rho) = -Tr(rho log2 rho)."""
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return -np.sum(evals * np.log2(evals))


def mutual_information(rho, qi, qj, n_qubits):
    """MI(i,j) = S(i) + S(j) - S(i,j)."""
    rho_ij = partial_trace(rho, [qi, qj], n_qubits)
    rho_i = partial_trace(rho, [qi], n_qubits)
    rho_j = partial_trace(rho, [qj], n_qubits)
    return von_neumann_entropy(rho_i) + von_neumann_entropy(rho_j) - von_neumann_entropy(rho_ij)


def compute_observables(rho, n_qubits):
    """Compute all observables for a density matrix."""
    purity = np.real(np.trace(rho @ rho))

    # Pairwise MI for adjacent pairs
    pairs = [(i, i+1) for i in range(n_qubits - 1)]
    mi_values = {}
    for qi, qj in pairs:
        mi_values[f"MI_{qi}{qj}"] = mutual_information(rho, qi, qj, n_qubits)

    sum_mi = sum(mi_values.values())

    # Pairwise CPsi
    cpsi_values = {}
    for qi, qj in pairs:
        rho_pair = partial_trace(rho, [qi, qj], n_qubits)
        pair_purity = np.real(np.trace(rho_pair @ rho_pair))
        # L1 coherence of pair
        l1 = np.sum(np.abs(rho_pair)) - np.sum(np.abs(np.diag(rho_pair)))
        d_pair = 4
        cpsi = pair_purity * l1 / (d_pair - 1)
        cpsi_values[f"CPsi_{qi}{qj}"] = cpsi

    return {
        'purity': purity,
        'sum_mi': sum_mi,
        **mi_values,
        **cpsi_values,
    }


def run_scenario(label, gammas, times, n_qubits, J, out):
    """Run time evolution for one scenario."""
    out(f"\n--- {label} ---")
    out(f"gammas: [{', '.join(f'{g:.4f}' for g in gammas)}]")
    out(f"sum(gamma): {sum(gammas):.4f}")

    H = build_heisenberg_chain(n_qubits, J)
    L = build_liouvillian(H, gammas)

    # Initial state: |+>^N
    d = 2**n_qubits
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi = plus
    for _ in range(n_qubits - 1):
        psi = np.kron(psi, plus)
    rho0 = np.outer(psi, psi.conj())

    # Validate
    assert abs(np.trace(rho0) - 1.0) < 1e-10, "rho0 trace != 1"
    assert abs(np.real(np.trace(rho0 @ rho0)) - 1.0) < 1e-10, "rho0 not pure"

    rho0_vec = rho0.flatten()

    # Time evolution via matrix exponential
    results = []
    for t in times:
        if t == 0:
            rho_vec = rho0_vec.copy()
        else:
            rho_vec = expm(L * t) @ rho0_vec

        rho = rho_vec.reshape(d, d)
        # Enforce hermiticity
        rho = (rho + rho.conj().T) / 2

        obs = compute_observables(rho, n_qubits)
        obs['t'] = t
        results.append(obs)

    # Print table
    out(f"\n{'t':>5} {'Purity':>8} {'SumMI':>8} {'MI_01':>7} {'MI_12':>7} {'MI_23':>7} {'MI_34':>7} {'CPsi01':>7} {'CPsi12':>7}")
    for r in results:
        out(f"{r['t']:5.1f} {r['purity']:8.4f} {r['sum_mi']:8.4f} "
            f"{r['MI_01']:7.4f} {r['MI_12']:7.4f} {r['MI_23']:7.4f} {r['MI_34']:7.4f} "
            f"{r['CPsi_01']:7.4f} {r['CPsi_12']:7.4f}")

    return results


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "time_evolution_6scenarios.txt"
    csv_path = results_dir / "time_evolution_plotdata.csv"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== TIME EVOLUTION: 6 SCENARIOS ===")
    out("N=5 Chain, J=1.0, |+>^5 initial state")
    out()

    N = 5
    J = 1.0
    times = np.arange(0, 5.5, 0.5)  # 0.0 to 5.0 in steps of 0.5

    # Load data
    json_path = Path(__file__).parent.parent / "data" / "ibm_sacrifice_zone_march2026" / "sacrifice_zone_hardware_20260324_191713.json"
    with open(json_path) as f:
        hw = json.load(f)

    ci = hw['chain_info']
    sac_T2star = ci['T2star_values']
    sac_T2echo = ci['T2_values']

    # Mean-T2 chain T2echo from CSV
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

    # Build 6 scenarios
    scenarios = [
        ("1 Mean-T2, No DD", [1/t for t in t2_T2star]),
        ("2 Mean-T2, Uniform DD", [1/t for t in t2_T2echo]),
        ("3 Mean-T2, Selective DD", [1/t2_T2star[0], 1/t2_T2echo[1], 1/t2_T2echo[2], 1/t2_T2echo[3], 1/t2_T2star[4]]),
        ("4 Sacrifice, No DD", [1/t for t in sac_T2star]),
        ("5 Sacrifice, Uniform DD", [1/t for t in sac_T2echo]),
        ("6 Sacrifice, Selective DD", [1/sac_T2star[0], 1/sac_T2echo[1], 1/sac_T2echo[2], 1/sac_T2echo[3], 1/sac_T2echo[4]]),
    ]

    all_results = {}
    csv_rows = []

    for label, gammas in scenarios:
        results = run_scenario(label, gammas, times, N, J, out)
        all_results[label] = results

        for r in results:
            csv_rows.append({
                'scenario': label,
                't': r['t'],
                'SumMI': r['sum_mi'],
                'Purity': r['purity'],
                'MI_01': r['MI_01'],
                'MI_12': r['MI_12'],
                'MI_23': r['MI_23'],
                'MI_34': r['MI_34'],
                'CPsi_01': r['CPsi_01'],
                'CPsi_12': r['CPsi_12'],
            })

    # Comparison at t=2.5 (midpoint)
    out(f"\n{'='*70}")
    out("COMPARISON AT t=2.5 us (midpoint)")
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

    # SumMI at all measured time points
    out(f"\n{'='*70}")
    out("SumMI EVOLUTION (all scenarios)")
    out(f"{'='*70}")
    headers = [f"{'t':>5}"]
    for label, _ in scenarios:
        short = label.split(",")[0].replace("Mean-T2", "T2").replace("Sacrifice", "Sac")
        headers.append(f"{short:>10}")
    out(" ".join(headers))

    for ti, t in enumerate(times):
        row = [f"{t:5.1f}"]
        for label, _ in scenarios:
            r = all_results[label][ti]
            row.append(f"{r['sum_mi']:10.4f}")
        out(" ".join(row))

    # Write CSV
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['scenario', 't', 'SumMI', 'Purity',
                                               'MI_01', 'MI_12', 'MI_23', 'MI_34',
                                               'CPsi_01', 'CPsi_12'])
        writer.writeheader()
        writer.writerows(csv_rows)

    out(f"\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")
    print(f">>> Plot data saved to: {csv_path}")


if __name__ == "__main__":
    main()
