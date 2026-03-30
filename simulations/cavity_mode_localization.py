"""
Cavity Mode Localization: Where do the protected modes live?

For the N=5 chain under IBM sacrifice-zone noise, the slowest
oscillating modes (freq ~ 7.234J) survive 2.81x longer than under
uniform noise. This script determines WHERE those modes are
spatially localized by decomposing Liouvillian eigenvectors into
the Pauli basis and computing per-qubit weights.

Hypothesis: Protected modes are localized on interior qubits
(away from the sacrifice qubit Q85 = index 0).
"""

import numpy as np
from scipy.linalg import eig
from scipy.stats import pearsonr, spearmanr
from pathlib import Path
import itertools


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
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n_qubits)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def generate_pauli_strings(n_qubits):
    """Generate all 4^N Pauli strings with their active qubit sets."""
    labels = ['I', 'X', 'Y', 'Z']
    strings = []
    for combo in itertools.product(labels, repeat=n_qubits):
        name = ''.join(combo)
        active = [k for k in range(n_qubits) if combo[k] != 'I']
        # Build the matrix
        mat = np.eye(1, dtype=complex)
        for c in combo:
            mat = np.kron(mat, PAULIS[c])
        strings.append((name, active, mat))
    return strings


def compute_qubit_weights(eigvec, pauli_strings, n_qubits, d):
    """Compute per-qubit weight of a Liouvillian eigenvector."""
    # Reshape eigenvector to operator
    V = eigvec.reshape(d, d)

    # Pauli decomposition: c_j = Tr(P_j^dag @ V) / d
    weights = np.zeros(n_qubits)
    total_weight = 0.0

    for name, active, P in pauli_strings:
        if not active:  # skip IIIII (identity)
            continue
        c = np.trace(P.conj().T @ V) / d
        w = abs(c)**2
        total_weight += w
        for k in active:
            weights[k] += w

    if total_weight > 1e-15:
        weights /= total_weight

    return weights, total_weight


def analyze_profile(name, gammas, H, n_qubits, pauli_strings, out):
    """Full eigenvector localization analysis."""
    d = 2**n_qubits
    d2 = d * d

    out(f"\n{'='*70}")
    out(f"PROFILE: {name}")
    out(f"{'='*70}")
    out(f"Gammas: [{', '.join(f'{g:.4f}' for g in gammas)}]")

    L = build_liouvillian(H, gammas)
    evals, evecs = eig(L)

    # Classify modes
    eps = 1e-10
    modes = []
    for i in range(len(evals)):
        ev = evals[i]
        rate = -ev.real
        freq = abs(ev.imag)
        if abs(ev.real) < eps and freq < eps:
            mode_type = 'stationary'
        else:
            mode_type = 'oscillating'
        modes.append((i, rate, freq, mode_type))

    osc_modes = [(i, r, f) for i, r, f, t in modes if t == 'oscillating']
    osc_modes.sort(key=lambda x: x[1])  # sort by rate (slowest first)

    out(f"Total modes: {len(modes)}")
    out(f"Oscillating: {len(osc_modes)}")

    # Compute qubit weights for oscillating modes
    out(f"\nComputing qubit weights for {len(osc_modes)} oscillating modes...")
    mode_weights = []
    for idx, rate, freq in osc_modes:
        v = evecs[:, idx]
        weights, total = compute_qubit_weights(v, pauli_strings, n_qubits, d)
        mode_weights.append((idx, rate, freq, weights, total))

    # Top 20 slowest
    out(f"\nTop 20 slowest oscillating modes:")
    out(f"{'#':>3} {'Rate':>10} {'Freq':>8} {'Q0(sac)':>8} {'Q1':>8} {'Q2':>8} {'Q3':>8} {'Q4':>8}")
    for rank, (idx, rate, freq, w, _) in enumerate(mode_weights[:20]):
        out(f"{rank+1:3d} {rate:10.6f} {freq:8.3f} {w[0]:8.3f} {w[1]:8.3f} {w[2]:8.3f} {w[3]:8.3f} {w[4]:8.3f}")

    # Top 20 fastest
    out(f"\nTop 20 fastest decaying modes:")
    fastest = sorted(mode_weights, key=lambda x: -x[1])
    out(f"{'#':>3} {'Rate':>10} {'Freq':>8} {'Q0(sac)':>8} {'Q1':>8} {'Q2':>8} {'Q3':>8} {'Q4':>8}")
    for rank, (idx, rate, freq, w, _) in enumerate(fastest[:20]):
        out(f"{rank+1:3d} {rate:10.6f} {freq:8.3f} {w[0]:8.3f} {w[1]:8.3f} {w[2]:8.3f} {w[3]:8.3f} {w[4]:8.3f}")

    # Correlation: weight on Q0 vs decay rate
    if len(mode_weights) > 10:
        w0_vals = [w[0] for _, _, _, w, _ in mode_weights]
        rate_vals = [r for _, r, _, _, _ in mode_weights]

        r_pearson, p_pearson = pearsonr(w0_vals, rate_vals)
        r_spearman, p_spearman = spearmanr(w0_vals, rate_vals)

        out(f"\nCorrelation: weight_on_Q0 vs decay_rate")
        out(f"  Pearson  r = {r_pearson:.4f} (p = {p_pearson:.2e})")
        out(f"  Spearman r = {r_spearman:.4f} (p = {p_spearman:.2e})")

    # The 7.234J modes specifically
    target_freq = 7.234
    target_modes = [(idx, r, f, w, t) for idx, r, f, w, t in mode_weights
                    if abs(f - target_freq) < 0.1 and r < 0.1]
    if target_modes:
        out(f"\nThe ~7.234J modes (freq ~ {target_freq}, rate < 0.1):")
        out(f"{'#':>3} {'Rate':>10} {'Freq':>8} {'Q0(sac)':>8} {'Q1':>8} {'Q2':>8} {'Q3':>8} {'Q4':>8}")
        for rank, (idx, rate, freq, w, _) in enumerate(target_modes):
            out(f"{rank+1:3d} {rate:10.6f} {freq:8.3f} {w[0]:8.3f} {w[1]:8.3f} {w[2]:8.3f} {w[3]:8.3f} {w[4]:8.3f}")
    else:
        out(f"\nNo modes found near freq={target_freq} with rate < 0.1")

    # Mean qubit profile for slowest vs fastest quartile
    n_quartile = max(1, len(mode_weights) // 4)
    slow_q = np.mean([w for _, _, _, w, _ in mode_weights[:n_quartile]], axis=0)
    fast_q = np.mean([w for _, _, _, w, _ in mode_weights[-n_quartile:]], axis=0)
    out(f"\nMean qubit profile (slowest quartile, {n_quartile} modes):")
    out(f"  [{', '.join(f'{v:.3f}' for v in slow_q)}]")
    out(f"Mean qubit profile (fastest quartile, {n_quartile} modes):")
    out(f"  [{', '.join(f'{v:.3f}' for v in fast_q)}]")

    return mode_weights


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "cavity_mode_localization.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== CAVITY MODE LOCALIZATION ANALYSIS ===")
    out(f"N=5 Chain, J=1.0")
    out()

    N = 5
    J = 1.0
    H = build_heisenberg_chain(N, J)

    out("Generating Pauli basis (1024 strings)...")
    pauli_strings = generate_pauli_strings(N)
    out(f"Done. {len(pauli_strings)} strings generated.")

    # Three profiles
    profiles = [
        ("IBM SACRIFICE ZONE", np.array([0.2681, 0.0163, 0.0103, 0.0147, 0.0105])),
        ("UNIFORM", np.array([0.0640, 0.0640, 0.0640, 0.0640, 0.0640])),
        ("ZERO NOISE (cavity modes)", np.array([0.0, 0.0, 0.0, 0.0, 0.0])),
    ]

    all_results = {}
    for name, gammas in profiles:
        all_results[name] = analyze_profile(name, gammas, H, N, pauli_strings, out)

    # Cross-profile comparison
    out(f"\n{'='*70}")
    out("CROSS-PROFILE COMPARISON: 7.234J modes")
    out(f"{'='*70}")

    for name, gammas in profiles:
        modes = all_results[name]
        target = [(r, f, w) for _, r, f, w, _ in modes if abs(f - 7.234) < 0.2]
        if target:
            target.sort(key=lambda x: x[0])
            w_mean = np.mean([w for _, _, w in target[:4]], axis=0)
            out(f"\n{name}: {len(target)} modes near 7.234J")
            out(f"  Mean qubit profile (4 slowest): [{', '.join(f'{v:.3f}' for v in w_mean)}]")
            out(f"  Weight on Q0 (sacrifice): {w_mean[0]:.3f}")
            out(f"  Weight on interior (Q1-Q4): {sum(w_mean[1:]):.3f}")

    # Summary
    out(f"\n{'='*70}")
    out("SUMMARY")
    out(f"{'='*70}")

    ibm_modes = all_results["IBM SACRIFICE ZONE"]
    if ibm_modes:
        w0_vals = [w[0] for _, _, _, w, _ in ibm_modes]
        rate_vals = [r for _, r, _, _, _ in ibm_modes]
        r_pearson, p_pearson = pearsonr(w0_vals, rate_vals)
        out(f"\nIBM Sacrifice: Correlation(weight_Q0, rate) = {r_pearson:.4f} (p={p_pearson:.2e})")
        if r_pearson > 0.3:
            out("CONFIRMED: Modes with more weight on Q0 (sacrifice) decay faster.")
            out("The sacrifice zone protects modes localized on interior qubits.")
        elif r_pearson < -0.3:
            out("UNEXPECTED: Modes with more weight on Q0 decay SLOWER.")
            out("This would contradict the localization hypothesis.")
        else:
            out("INCONCLUSIVE: No strong correlation between Q0 weight and decay rate.")
            out("The protection mechanism may not be simple spatial localization.")

    out()
    out("=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
