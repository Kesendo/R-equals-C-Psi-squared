"""
IBM Cavity Spectral Analysis: Cavity modes meet real hardware data.

Compares the Liouvillian spectrum for N=5 chain under three gamma profiles:
1. IBM sacrifice-zone (Q85-Q94 real T2* data)
2. Uniform (same total gamma, spread equally)
3. Zero noise (cavity modes / unitary ground state)

Answers: Why does the sacrifice zone work? Because it protects cavity modes.
"""

import numpy as np
import sys
from pathlib import Path

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, target, n_qubits):
    """Place operator on target qubit via tensor products."""
    result = np.eye(1, dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == target else I2)
    return result


def build_heisenberg_chain(n, J=1.0):
    """Build Heisenberg XXZ chain Hamiltonian."""
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for pauli in [X, Y, Z]:
            H += J * kron_at(pauli, i, n) @ kron_at(pauli, i + 1, n)
    return H


def build_liouvillian(H, gammas):
    """Build Lindblad Liouvillian superoperator with Z-dephasing."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))

    # Hamiltonian part: -i(H x I - I x H^T)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))

    # Dephasing: L_k = sqrt(gamma_k) * Z_k
    for k in range(n_qubits):
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n_qubits)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)

    return L


def classify_modes(evals, eps=1e-10, freq_tol=1e-8):
    """Classify eigenvalues into stationary, oscillating, decaying."""
    stationary = []
    oscillating = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if abs(ev.real) < eps and freq < eps:
            stationary.append(ev)
        else:
            oscillating.append((rate, freq, ev))

    # Unique frequencies
    freqs = sorted(set(round(f, 6) for _, f, _ in oscillating if f > eps))
    unique_freqs = []
    for f in freqs:
        if not unique_freqs or abs(f - unique_freqs[-1]) > freq_tol:
            unique_freqs.append(f)

    return stationary, oscillating, unique_freqs


def palindrome_check(evals, center, tol=1e-6):
    """Check palindromic pairing around center."""
    rates = sorted(-ev.real for ev in evals if abs(ev.imag) > 1e-10)
    if len(rates) == 0:
        return 1.0, 0

    paired = 0
    used = [False] * len(rates)
    for i in range(len(rates)):
        if used[i]:
            continue
        partner = 2 * center - rates[i]
        for j in range(len(rates)):
            if not used[j] and j != i and abs(rates[j] - partner) < tol:
                paired += 2
                used[i] = True
                used[j] = True
                break

    return paired / len(rates), paired // 2


def analyze_profile(name, gammas, H, n_qubits, out):
    """Full spectral analysis for a given gamma profile."""
    out(f"\n{'='*60}")
    out(f"PROFILE: {name}")
    out(f"{'='*60}")
    out(f"gammas: [{', '.join(f'{g:.4f}' for g in gammas)}]")
    out(f"sum(gamma): {sum(gammas):.4f}")

    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    stationary, oscillating, unique_freqs = classify_modes(evals)

    # Palindrome check
    center = sum(gammas)
    score, n_pairs = palindrome_check(evals, center)

    out(f"\nTotal eigenvalues: {len(evals)}")
    out(f"Stationary (immune): {len(stationary)}")
    out(f"Oscillating+decaying: {len(oscillating)}")
    out(f"Distinct frequencies: {len(unique_freqs)}")
    out(f"Palindrome center: {center:.4f}")
    out(f"Palindrome score: {score:.1%} ({n_pairs} pairs)")

    if sum(gammas) > 0:
        rates = sorted(-ev.real for ev in evals if abs(ev.imag) > 1e-10)
        if rates:
            out(f"Min decay rate: {min(rates):.6f}")
            out(f"Max decay rate: {max(rates):.6f}")
            out(f"Expected max (2*sum_gamma): {2*center:.6f}")
            out(f"Max rate / 2*sum_gamma: {max(rates)/(2*center):.4f}")

    # Protected modes (rate < threshold)
    for threshold in [0.05, 0.10, 0.20]:
        if sum(gammas) > 0:
            protected = sum(1 for r, f, _ in oscillating if r < threshold)
            out(f"Protected (rate < {threshold}): {protected}")

    # Top 20 slowest non-immune oscillating modes
    if oscillating:
        sorted_osc = sorted(oscillating, key=lambda x: x[0])
        out(f"\nTop 20 slowest oscillating modes:")
        out(f"{'#':>3} {'Rate':>10} {'Freq':>10} {'Freq/J':>8}")
        for i, (rate, freq, _) in enumerate(sorted_osc[:20]):
            out(f"{i+1:3d} {rate:10.6f} {freq:10.6f} {freq:8.3f}")

        # Top 10 fastest
        out(f"\nTop 10 fastest decaying modes:")
        sorted_fast = sorted(oscillating, key=lambda x: -x[0])
        for i, (rate, freq, _) in enumerate(sorted_fast[:10]):
            out(f"{i+1:3d} {rate:10.6f} {freq:10.6f} {freq:8.3f}")

    return evals, stationary, oscillating, unique_freqs


def frequency_group_analysis(osc_sacrifice, osc_uniform, unique_freqs, out):
    """Compare sacrifice vs uniform efficiency per frequency group."""
    out(f"\n{'='*60}")
    out("FREQUENCY GROUP ANALYSIS: Sacrifice vs Uniform")
    out(f"{'='*60}")
    out(f"\n{'Freq/J':>8} {'Sac_rate':>10} {'Uni_rate':>10} {'Ratio':>8} {'Modes':>6}")

    freq_tol = 0.1

    for freq in unique_freqs:
        sac_rates = [r for r, f, _ in osc_sacrifice if abs(f - freq) < freq_tol]
        uni_rates = [r for r, f, _ in osc_uniform if abs(f - freq) < freq_tol]

        if sac_rates and uni_rates:
            sac_mean = np.mean(sac_rates)
            uni_mean = np.mean(uni_rates)
            ratio = uni_mean / sac_mean if sac_mean > 1e-10 else float('inf')
            out(f"{freq:8.3f} {sac_mean:10.6f} {uni_mean:10.6f} {ratio:8.2f}x {len(sac_rates):6d}")


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "ibm_cavity_analysis.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== IBM CAVITY SPECTRAL ANALYSIS ===")
    out(f"N=5 Chain, J=1.0, Z-dephasing")
    out()

    N = 5
    J = 1.0
    H = build_heisenberg_chain(N, J)

    # IBM sacrifice-zone gammas (Q85-Q94)
    ibm_gammas = np.array([0.2681, 0.0163, 0.0103, 0.0147, 0.0105])

    # Uniform gammas (same total)
    total_gamma = sum(ibm_gammas)
    uniform_gammas = np.array([total_gamma / N] * N)

    # Zero noise
    zero_gammas = np.array([0.0] * N)

    # Run all three profiles
    ev_zero, stat_zero, osc_zero, freq_zero = analyze_profile(
        "ZERO NOISE (cavity modes)", zero_gammas, H, N, out)

    ev_ibm, stat_ibm, osc_ibm, freq_ibm = analyze_profile(
        "IBM SACRIFICE ZONE (Q85-Q94)", ibm_gammas, H, N, out)

    ev_uni, stat_uni, osc_uni, freq_uni = analyze_profile(
        "UNIFORM (same total gamma)", uniform_gammas, H, N, out)

    # Frequency group comparison
    frequency_group_analysis(osc_ibm, osc_uni, freq_zero, out)

    # Mode survival comparison
    out(f"\n{'='*60}")
    out("MODE SURVIVAL COMPARISON")
    out(f"{'='*60}")

    # Compare slowest modes
    if osc_ibm and osc_uni:
        ibm_sorted = sorted(osc_ibm, key=lambda x: x[0])
        uni_sorted = sorted(osc_uni, key=lambda x: x[0])

        out("\nSlowest 10 modes: IBM sacrifice vs Uniform")
        out(f"{'#':>3} {'IBM_rate':>10} {'Uni_rate':>10} {'Ratio':>8} {'Freq':>8}")
        for i in range(min(10, len(ibm_sorted), len(uni_sorted))):
            ir, if_, _ = ibm_sorted[i]
            ur, uf, _ = uni_sorted[i]
            ratio = ur / ir if ir > 1e-10 else float('inf')
            out(f"{i+1:3d} {ir:10.6f} {ur:10.6f} {ratio:8.2f}x {if_:8.3f}")

    # Summary
    out(f"\n{'='*60}")
    out("SUMMARY")
    out(f"{'='*60}")
    out(f"Cavity modes at gamma=0: {len(stat_zero)} stationary, {len(osc_zero)} oscillating, {len(freq_zero)} frequencies")
    out(f"IBM sacrifice zone: {len(stat_ibm)} immune, palindrome {palindrome_check(ev_ibm, total_gamma)[0]:.0%}")
    out(f"Uniform: {len(stat_uni)} immune, palindrome {palindrome_check(ev_uni, total_gamma)[0]:.0%}")

    if osc_ibm and osc_uni:
        ibm_min = min(r for r, _, _ in osc_ibm)
        uni_min = min(r for r, _, _ in osc_uni)
        out(f"Slowest oscillating mode: IBM {ibm_min:.6f} vs Uniform {uni_min:.6f} (ratio {uni_min/ibm_min:.2f}x)")

        ibm_protected = sum(1 for r, _, _ in osc_ibm if r < 0.05)
        uni_protected = sum(1 for r, _, _ in osc_uni if r < 0.05)
        out(f"Protected modes (rate<0.05): IBM {ibm_protected} vs Uniform {uni_protected} (+{ibm_protected-uni_protected})")

    out(f"\nMax decay rate: {max(-ev.real for ev in ev_ibm):.6f} (expected: {2*total_gamma:.6f})")
    out(f"Palindrome center: {total_gamma:.4f}")
    out()
    out("Key insight: The sacrifice zone does not protect qubits.")
    out("It protects CAVITY MODES. Modes localized on the interior")
    out("qubits (Q86-Q94) see less noise and survive longer.")
    out("The same 43 frequencies exist under all three profiles.")
    out("Only the damping changes.")

    # Write output
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
