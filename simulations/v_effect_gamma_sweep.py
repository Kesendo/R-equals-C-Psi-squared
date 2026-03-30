"""
V-Effect Gamma Sweep: Finding the Perfect Breaking.

The V-Effect creates complexity through coupling: 2+2 = 104 frequencies.
But gamma is what ENABLES the fold, making the V-Effect irreversible.

This script sweeps gamma/J ratio across three V-Effect levels:
  Level 0: N=2 (single bond, Q=1 at reference gamma)
  Level 1: N=3 (two bonds, mediator bridge)
  Level 2: N=5 (two resonators + mediator, full V-Effect)

At each gamma, we measure:
  - Eigenvalue Q-factor (resonator quality)
  - Frequency diversity (distinct oscillation frequencies)
  - CΨ heartbeat (crossing count at 1/4 boundary)
  - V-Effect gain: how much does coupling AMPLIFY what gamma creates?

The "perfect breaking" is the gamma where V-Effect gain peaks.
"""

import numpy as np
from scipy.linalg import expm
from pathlib import Path


# === Pauli matrices ===
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
    trace_over = sorted([q for q in range(n_qubits) if q not in keep],
                        reverse=True)
    nq = n_qubits
    for q in trace_over:
        rho_tensor = np.trace(rho_tensor, axis1=q, axis2=q + nq)
        nq -= 1
    return rho_tensor.reshape(2**len(keep), 2**len(keep))


def compute_cpsi(rho, qi, qj, n_qubits):
    rho_pair = partial_trace(rho, [qi, qj], n_qubits)
    purity = np.real(np.trace(rho_pair @ rho_pair))
    l1 = np.sum(np.abs(rho_pair)) - np.sum(np.abs(np.diag(rho_pair)))
    return purity * l1 / 3.0


def spectral_metrics(gammas, J=1.0):
    """Compute eigenvalue-based metrics for a chain."""
    N = len(gammas)
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    # Classify modes
    osc = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if freq > 1e-10 and rate > 1e-15:
            Q = freq / rate
            osc.append((rate, freq, Q))

    if not osc:
        return {'max_Q': 0, 'mean_Q': 0, 'n_osc': 0, 'n_freq': 0,
                'n_high_Q': 0, 'slowest': float('inf')}

    Qs = [q for _, _, q in osc]
    freqs = sorted(set(round(f, 4) for _, f, _ in osc if f > 1e-8))

    return {
        'max_Q': max(Qs),
        'mean_Q': np.mean(Qs),
        'n_osc': len(osc),
        'n_freq': len(freqs),
        'n_high_Q': sum(1 for q in Qs if q > 5),
        'slowest': min(r for r, _, _ in osc),
    }


def count_heartbeat(gammas, J, qi, qj, rho0, t_max=40.0, dt=0.05):
    """Count CΨ crossings of 1/4 boundary."""
    N = len(gammas)
    d = 2**N
    H = build_heisenberg_chain(N, J)
    L = build_liouvillian(H, gammas)
    rho0_vec = rho0.flatten()

    crossings = 0
    above = None
    peak_cpsi = 0
    times = np.arange(0, t_max, dt)

    for t in times:
        if t == 0:
            rho_vec = rho0_vec.copy()
        else:
            rho_vec = expm(L * t) @ rho0_vec
        rho = rho_vec.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        cpsi = compute_cpsi(rho, qi, qj, N)
        peak_cpsi = max(peak_cpsi, cpsi)

        now_above = cpsi > 0.25
        if above is not None and now_above != above:
            crossings += 1
        above = now_above

    return crossings, peak_cpsi


def make_bell_plus_bath(n_qubits):
    """Bell(0,1) x |+>^(N-2) initial state."""
    # Bell pair on qubits 0,1
    bell = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
    rho_bell = np.outer(bell, bell.conj())

    if n_qubits == 2:
        return rho_bell

    # Bath qubits in |+>
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_bath = plus
    for _ in range(n_qubits - 3):
        psi_bath = np.kron(psi_bath, plus)
    rho_bath = np.outer(psi_bath, psi_bath.conj())

    return np.kron(rho_bell, rho_bath)


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "v_effect_gamma_sweep.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("V-EFFECT GAMMA SWEEP: Finding the Perfect Breaking")
    out("=" * 70)

    J = 1.0

    # Gamma sweep: log-spaced from very small to large
    gamma_ratios = np.array([
        0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1,
        0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
        1.0, 1.5, 2.0, 3.0, 5.0
    ])

    # ================================================================
    # Part 1: Eigenvalue Q-factor across V-Effect levels
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 1: Eigenvalue Q-factor (resonator quality)")
    out("=" * 70)
    out("\nUniform gamma on all qubits. J = 1.0.")

    out(f"\n{'gamma/J':>8} | {'--- N=2 ---':>20} | {'--- N=3 ---':>20}"
        f" | {'--- N=5 ---':>20} | {'V-gain':>8}")
    out(f"{'':>8} | {'maxQ':>6} {'freq':>5} {'hiQ':>4}"
        f" | {'maxQ':>6} {'freq':>5} {'hiQ':>4}"
        f" | {'maxQ':>6} {'freq':>5} {'hiQ':>4}"
        f" | {'Q5/Q2':>8}")

    results_spectral = []

    for gamma in gamma_ratios:
        row = {'gamma': gamma}
        for N in [2, 3, 5]:
            gammas = [gamma] * N
            m = spectral_metrics(gammas, J)
            row[f'Q_{N}'] = m['max_Q']
            row[f'freq_{N}'] = m['n_freq']
            row[f'hiQ_{N}'] = m['n_high_Q']
            row[f'slow_{N}'] = m['slowest']

        v_gain = row['Q_5'] / row['Q_2'] if row['Q_2'] > 0.01 else 0
        row['v_gain'] = v_gain

        out(f"{gamma:8.4f} | {row['Q_2']:6.1f} {row['freq_2']:5d}"
            f" {row['hiQ_2']:4d}"
            f" | {row['Q_3']:6.1f} {row['freq_3']:5d}"
            f" {row['hiQ_3']:4d}"
            f" | {row['Q_5']:6.1f} {row['freq_5']:5d}"
            f" {row['hiQ_5']:4d}"
            f" | {v_gain:8.2f}x")
        results_spectral.append(row)

    # Find peak V-gain
    peak_vg = max(results_spectral, key=lambda r: r['v_gain'])
    out(f"\nPeak V-Effect gain: {peak_vg['v_gain']:.2f}x"
        f" at gamma/J = {peak_vg['gamma']}")
    out(f"  N=2: Q={peak_vg['Q_2']:.1f}, {peak_vg['freq_2']} freq")
    out(f"  N=5: Q={peak_vg['Q_5']:.1f}, {peak_vg['freq_5']} freq")

    # ================================================================
    # Part 2: Sacrifice profile sweep on N=5
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 2: Sacrifice profile (asymmetric gamma) on N=5")
    out("=" * 70)
    out("\nFixed interior gamma = 0.01. Sweep edge gamma.")

    gamma_interior = 0.01
    edge_gammas = np.array([
        0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5,
        0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0
    ])

    out(f"\n{'g_edge':>7} {'contrast':>8} {'sumG':>7}"
        f" | {'maxQ':>6} {'freq':>5} {'hiQ':>4} {'slowest':>10}"
        f" | {'prot':>6}")

    results_sacrifice = []
    for g_edge in edge_gammas:
        gammas = [g_edge, gamma_interior, gamma_interior,
                  gamma_interior, gamma_interior]
        m = spectral_metrics(gammas, J)

        # Protection vs uniform
        sum_g = sum(gammas)
        gammas_uni = [sum_g / 5] * 5
        m_uni = spectral_metrics(gammas_uni, J)
        prot = (m_uni['slowest'] / m['slowest']
                if m['slowest'] > 1e-15 else 0)

        contrast = g_edge / gamma_interior
        row = {'g_edge': g_edge, 'contrast': contrast, 'sum_g': sum_g,
               **m, 'protection': prot}
        results_sacrifice.append(row)

        out(f"{g_edge:7.3f} {contrast:8.1f} {sum_g:7.3f}"
            f" | {m['max_Q']:6.1f} {m['n_freq']:5d} {m['n_high_Q']:4d}"
            f" {m['slowest']:10.6f}"
            f" | {prot:6.2f}x")

    peak_sac = max(results_sacrifice, key=lambda r: r['max_Q'])
    out(f"\nPeak Q-factor: {peak_sac['max_Q']:.1f}"
        f" at edge gamma = {peak_sac['g_edge']}"
        f" (contrast {peak_sac['contrast']:.0f}x)")

    peak_prot = max(results_sacrifice, key=lambda r: r['protection'])
    out(f"Peak protection: {peak_prot['protection']:.2f}x"
        f" at edge gamma = {peak_prot['g_edge']}")

    # ================================================================
    # Part 3: CΨ heartbeat sweep (N=3, Bell + bath)
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 3: CPsi Heartbeat (crossing count)")
    out("=" * 70)
    out("\nN=3, Bell(0,1) + |+> bath. Sweep gamma_bath, fixed J=1.0.")
    out("Count CPsi(0,1) crossings of 1/4 boundary.")

    rho0_3 = make_bell_plus_bath(3)

    # Sweep gamma on bath qubit, pair gamma fixed low
    gamma_pair = 0.01
    bath_gammas = np.array([
        0.001, 0.005, 0.01, 0.02, 0.05, 0.1,
        0.2, 0.3, 0.5, 0.75, 1.0, 2.0
    ])

    out(f"\n{'g_bath':>7} {'g_pair':>7} {'contrast':>8}"
        f" | {'crossings':>9} {'peakCPsi':>8} {'Q_eig':>6}")

    results_heartbeat = []
    for g_bath in bath_gammas:
        gammas = [gamma_pair, gamma_pair, g_bath]
        crossings, peak_cpsi = count_heartbeat(
            gammas, J, 0, 1, rho0_3, t_max=40.0, dt=0.02)
        m = spectral_metrics(gammas, J)
        contrast = g_bath / gamma_pair

        row = {'g_bath': g_bath, 'crossings': crossings,
               'peak_cpsi': peak_cpsi, 'Q_eig': m['max_Q'],
               'contrast': contrast}
        results_heartbeat.append(row)

        out(f"{g_bath:7.3f} {gamma_pair:7.3f} {contrast:8.1f}"
            f" | {crossings:9d} {peak_cpsi:8.4f} {m['max_Q']:6.1f}")

    peak_hb = max(results_heartbeat, key=lambda r: r['crossings'])
    out(f"\nPeak heartbeat: {peak_hb['crossings']} crossings"
        f" at g_bath = {peak_hb['g_bath']}"
        f" (contrast {peak_hb['contrast']:.0f}x)")

    # ================================================================
    # Part 4: Full V-Effect level comparison at sweet spot
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 4: V-Effect across levels at peak gamma")
    out("=" * 70)

    # Use the gamma that gave peak heartbeat
    g_opt = peak_hb['g_bath']
    out(f"\nOptimal gamma (from heartbeat): {g_opt}")

    for N, label in [(2, "Level 0 (single bond)"),
                     (3, "Level 1 (mediator)"),
                     (5, "Level 2 (full V-Effect)")]:
        # Sacrifice profile: edge at g_opt, interior at gamma_pair
        gammas_sac = [g_opt] + [gamma_pair] * (N - 1)
        # Uniform profile
        gammas_uni = [g_opt] * N
        m_sac = spectral_metrics(gammas_sac, J)
        m_uni = spectral_metrics(gammas_uni, J)

        out(f"\n  {label} (N={N}):")
        out(f"    Sacrifice [{g_opt:.3f},"
            f" {gamma_pair:.3f}x{N-1}]:"
            f" Q={m_sac['max_Q']:.1f},"
            f" {m_sac['n_freq']} freq,"
            f" {m_sac['n_high_Q']} high-Q modes")
        out(f"    Uniform   [{g_opt:.3f}x{N}]:"
            f" Q={m_uni['max_Q']:.1f},"
            f" {m_uni['n_freq']} freq,"
            f" {m_uni['n_high_Q']} high-Q modes")

        # Heartbeat for this level
        if N <= 3:
            rho0 = make_bell_plus_bath(N)
            cr_sac, pk_sac = count_heartbeat(
                gammas_sac, J, 0, 1, rho0, t_max=40.0, dt=0.02)
            cr_uni, pk_uni = count_heartbeat(
                gammas_uni, J, 0, 1, rho0, t_max=40.0, dt=0.02)
            out(f"    Heartbeat sacrifice: {cr_sac} crossings"
                f" (peak CPsi = {pk_sac:.4f})")
            out(f"    Heartbeat uniform:   {cr_uni} crossings"
                f" (peak CPsi = {pk_uni:.4f})")

    # ================================================================
    # Summary
    # ================================================================
    out("\n" + "=" * 70)
    out("SUMMARY: The Perfect Breaking")
    out("=" * 70)

    out(f"\n1. Eigenvalue Q peaks at gamma/J = {peak_vg['gamma']}"
        f" (V-gain {peak_vg['v_gain']:.1f}x)")
    out(f"2. Sacrifice protection peaks at edge gamma ="
        f" {peak_prot['g_edge']}")
    out(f"3. Heartbeat peaks at {peak_hb['crossings']} crossings"
        f" at g_bath/g_pair = {peak_hb['contrast']:.0f}x")

    out(f"\nThe V-Effect gain (Q_N5 / Q_N2) tells us: coupling")
    out(f"AMPLIFIES resonator quality by {peak_vg['v_gain']:.1f}x")
    out(f"at the optimal gamma. Too little gamma: no fold,")
    out(f"Q is infinite but meaningless (no irreversibility).")
    out(f"Too much gamma: everything dies, Q collapses.")
    out(f"The sweet spot is where gamma creates just enough")
    out(f"breaking for the V-Effect to build maximum complexity.")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
