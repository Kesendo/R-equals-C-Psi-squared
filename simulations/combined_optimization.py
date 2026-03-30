"""
Combined Optimization: Mapping + Selective DD for IBM April Run.

Six scenarios comparing chain selection and DD strategies:
1. Mean-T2 chain, no DD (baseline)
2. Mean-T2 chain, uniform DD
3. Mean-T2 chain, selective DD
4. Sacrifice chain, no DD (mapping only)
5. Sacrifice chain, uniform DD
6. Sacrifice chain, selective DD (the combination)
"""

import numpy as np
import json
import csv
from pathlib import Path


# === Physics (from ibm_cavity_analysis.py) ===
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


def analyze(gammas, J=1.0):
    """Spectral analysis for a gamma profile."""
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
            if not used[j] and j != i and abs(rates[j] - partner) < max(1e-4, 1e-3 * center):
                paired += 2
                used[i] = True
                used[j] = True
                break
    pal_score = paired / max(len(rates), 1)

    osc = [(-ev.real, abs(ev.imag)) for ev in evals if abs(ev.imag) > 1e-10]
    slowest = min(r for r, _ in osc) if osc else float('inf')

    # vs own uniform
    uni_g = np.array([sum_g / N] * N)
    L_uni = build_liouvillian(H, uni_g)
    evals_uni = np.linalg.eigvals(L_uni)
    osc_uni = [(-ev.real, abs(ev.imag)) for ev in evals_uni if abs(ev.imag) > 1e-10]
    slowest_uni = min(r for r, _ in osc_uni) if osc_uni else float('inf')

    return {
        'sum_gamma': sum_g,
        'slowest': slowest,
        'slowest_uni': slowest_uni,
        'vs_uniform': slowest_uni / slowest if slowest > 1e-15 else 0,
        'palindrome': pal_score,
    }


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "combined_optimization.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== COMBINED OPTIMIZATION: MAPPING + SELECTIVE DD ===")
    out("Preparation for IBM Hardware Run, April 9, 2026")
    out()

    # Load sacrifice chain data
    json_path = Path(__file__).parent.parent / "data" / "ibm_sacrifice_zone_march2026" / "sacrifice_zone_hardware_20260324_191713.json"
    with open(json_path) as f:
        hw = json.load(f)

    ci = hw['chain_info']
    sac_chain = hw['chain']  # [85, 86, 87, 88, 94]
    sac_T2star = ci['T2star_values']
    sac_T2echo = ci['T2_values']
    sac_ramsey = ci['ramsey_ratios']

    out(f"Sacrifice chain: {sac_chain}")
    out(f"  T2* (us):    [{', '.join(f'{v:.1f}' for v in sac_T2star)}]")
    out(f"  T2echo (us): [{', '.join(f'{v:.1f}' for v in sac_T2echo)}]")
    out(f"  Ramsey ratio:[{', '.join(f'{v:.1f}' for v in sac_ramsey)}]")
    out()

    # Load mean-T2 chain data from CSV
    csv_path = Path(__file__).parent.parent / "data" / "ibm_history" / "ibm_torino_history.csv"
    t2_chain = [18, 89, 19, 90, 60]  # best mean-T2 chain from mapping

    # Get latest T2 values for this chain
    qubit_t2 = {}
    latest_date = None
    with open(csv_path) as f:
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
    # Estimate T2* from T2echo using typical Ramsey ratio of 2.0
    default_ramsey = 2.0
    t2_T2star = [t2 / default_ramsey for t2 in t2_T2echo]

    out(f"Mean-T2 chain: {t2_chain}")
    out(f"  T2echo (us): [{', '.join(f'{v:.1f}' for v in t2_T2echo)}]")
    out(f"  T2* est (us):[{', '.join(f'{v:.1f}' for v in t2_T2star)}] (Ramsey ratio = {default_ramsey})")
    out()

    # === Build 6 scenarios ===
    scenarios = []

    # Scenario 1: Mean-T2, no DD
    g1 = [1.0 / t for t in t2_T2star]
    scenarios.append(("1 Base    ", "Mean-T2", "None   ", g1, t2_chain))

    # Scenario 2: Mean-T2, uniform DD (all qubits get T2echo)
    g2 = [1.0 / t for t in t2_T2echo]
    scenarios.append(("2 T2+DD   ", "Mean-T2", "All    ", g2, t2_chain))

    # Scenario 3: Mean-T2, selective DD (DD on inner 3 only)
    g3 = [1.0 / t2_T2star[0],   # edge: no DD
          1.0 / t2_T2echo[1],   # inner: DD
          1.0 / t2_T2echo[2],   # inner: DD
          1.0 / t2_T2echo[3],   # inner: DD
          1.0 / t2_T2star[4]]   # edge: no DD
    scenarios.append(("3 T2+Sel  ", "Mean-T2", "Select ", g3, t2_chain))

    # Scenario 4: Sacrifice, no DD
    g4 = [1.0 / t for t in sac_T2star]
    scenarios.append(("4 Sac     ", "Sacrif ", "None   ", g4, sac_chain))

    # Scenario 5: Sacrifice, uniform DD
    g5 = [1.0 / t for t in sac_T2echo]
    scenarios.append(("5 Sac+DD  ", "Sacrif ", "All    ", g5, sac_chain))

    # Scenario 6: Sacrifice, selective DD (DD on inner, no DD on sacrifice Q85)
    g6 = [1.0 / sac_T2star[0],   # sacrifice: no DD
          1.0 / sac_T2echo[1],   # inner: DD
          1.0 / sac_T2echo[2],   # inner: DD
          1.0 / sac_T2echo[3],   # inner: DD
          1.0 / sac_T2echo[4]]   # quiet edge: DD
    scenarios.append(("6 Sac+Sel ", "Sacrif ", "Select ", g6, sac_chain))

    # === Run analysis ===
    out(f"{'='*90}")
    out("SCENARIO ANALYSIS")
    out(f"{'='*90}")
    out()
    out(f"{'Scenario':>10} {'Chain':>8} {'DD':>8} {'sum_g':>8} {'Slowest':>10} {'Uni_slow':>10} {'vs_uni':>8} {'Palindrome':>10}")

    baseline_slowest = None
    results = []

    for label, chain_type, dd, gammas, chain in scenarios:
        r = analyze(gammas)
        if baseline_slowest is None:
            baseline_slowest = r['slowest']
        vs_baseline = baseline_slowest / r['slowest'] if r['slowest'] > 1e-15 else 0

        out(f"{label:>10} {chain_type:>8} {dd:>8} {r['sum_gamma']:8.4f} {r['slowest']:10.6f} {r['slowest_uni']:10.6f} {r['vs_uniform']:8.2f}x {r['palindrome']:10.0%}")

        results.append({
            'label': label.strip(),
            'chain_type': chain_type.strip(),
            'dd': dd.strip(),
            'chain': chain,
            'gammas': gammas,
            **r,
            'vs_baseline': vs_baseline,
        })

    # Summary table with vs_baseline
    out()
    out(f"{'='*90}")
    out("COMPARISON vs BASELINE (Scenario 1)")
    out(f"{'='*90}")
    out()
    out(f"{'Scenario':>10} {'vs_baseline':>12} {'vs_uniform':>12} {'sum_g':>8} {'Slowest':>10}")

    for r in results:
        out(f"{r['label']:>10} {r['vs_baseline']:12.2f}x {r['vs_uniform']:12.2f}x {r['sum_gamma']:8.4f} {r['slowest']:10.6f}")

    # Key comparisons
    out()
    out(f"{'='*90}")
    out("KEY COMPARISONS")
    out(f"{'='*90}")

    r4 = results[3]  # Sac no DD
    r6 = results[5]  # Sac+Sel
    r1 = results[0]  # Baseline
    r2 = results[1]  # T2+DD

    out(f"\nMapping only (4) vs Baseline (1): {r4['vs_baseline']:.2f}x")
    out(f"Best T2 + DD (2) vs Baseline (1): {r2['vs_baseline']:.2f}x")
    out(f"Mapping + Selective DD (6) vs Baseline (1): {r6['vs_baseline']:.2f}x")
    out(f"Mapping + Selective DD (6) vs Mapping only (4): {r4['slowest']/r6['slowest']:.2f}x")

    if r6['vs_baseline'] > r4['vs_baseline'] * 1.1:
        out(f"\nSELECTIVE DD ADDS {r6['vs_baseline']/r4['vs_baseline']:.2f}x ON TOP OF MAPPING.")
        out("The combination is synergistic.")
    else:
        out(f"\nSelective DD adds only {r6['vs_baseline']/r4['vs_baseline']:.2f}x on top of mapping.")
        out("Marginal improvement; mapping alone captures most of the benefit.")

    # IBM Experiment Plan
    out()
    out(f"{'='*90}")
    out("IBM EXPERIMENT PLAN (April 9, 2026)")
    out(f"{'='*90}")

    out("""
OBJECTIVE: Verify that sacrifice-zone chain selection + selective DD
outperforms standard practice (best-T2 chain + uniform DD) on IBM Torino.

CHAINS TO MEASURE:
  A. Sacrifice chain: {sac} (Q85 = sacrifice, T2* = 3.7 us)
  B. Mean-T2 chain:   {t2c} (all quiet, mean T2 ~ 200 us)

DD CONFIGURATIONS PER CHAIN:
  1. No DD (raw T2*)
  2. Uniform DD (XZX on all 5 qubits)
  3. Selective DD (XZX on inner 3-4, no DD on sacrifice/edge)

TROTTER EVOLUTION:
  dt = 0.5 us, steps = [2, 4, 6, 8, 10] (t = 1, 2, 3, 4, 5 us)
  J_coupling = 1.0 (standard Heisenberg)

SHOTS: 4000 per configuration per time point

TOTAL CIRCUITS: 2 chains x 3 DD configs x 5 time points = 30 circuits
TOTAL SHOTS: 30 x 4000 = 120,000

OBSERVABLE: Mutual information between adjacent pairs (SumMI)
  Compute from bitstring statistics via reduced density matrix tomography.

EXPECTED RESULTS (from simulation):
  Sacrifice + Selective DD: {r6_vs:.2f}x vs baseline
  Sacrifice + No DD:        {r4_vs:.2f}x vs baseline
  Mean-T2 + Uniform DD:     {r2_vs:.2f}x vs baseline
  Mean-T2 + No DD:          1.00x (baseline)

SUCCESS CRITERION:
  Sacrifice chain + Selective DD measurably outperforms Mean-T2 + Uniform DD.
  Minimum: 1.5x difference. Expected: ~{ratio:.1f}x.

FALLBACK:
  If Selective DD does not help: Mapping alone ({r4_vs:.2f}x) vs Mean-T2 + DD ({r2_vs:.2f}x)
  is still a testable prediction.
""".format(
        sac=sac_chain, t2c=t2_chain,
        r6_vs=r6['vs_baseline'], r4_vs=r4['vs_baseline'], r2_vs=r2['vs_baseline'],
        ratio=r6['vs_baseline'] / r2['vs_baseline']
    ))

    out("=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
