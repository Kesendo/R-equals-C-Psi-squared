"""
V-Effect with thermal breaking.

Three types of symmetry breaking:
  1. Bond breaking (V-Effect): coupling creates new palindromic pairs
  2. Dephasing (gamma_z): lifts degeneracies, preserves pairing (exact)
  3. Thermal (n_bar > 0): amplitude damping adds energy and explodes the
     frequency count. It does NOT break the pairing: the palindrome survives
     at a centre of -Sum(gamma_down + gamma_up)/2 (F137, extended 2026-08-05).
     What breaks the pairing is amplitude damping beside CO-AXIAL Z-dephasing

Question: Does heat change the V-Effect gain (1.81x)?
Does thermal breaking create frequencies that dephasing alone cannot?
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from framework import f1_distance_in_eps  # noqa: E402


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma_plus = |1><0|
SM = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma_minus = |0><1|


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


def build_liouvillian_thermal(H, gamma_z, gamma_amp, n_bar):
    """Build Liouvillian with Z-dephasing + thermal amplitude damping.

    Three noise channels per qubit:
      1. Z-dephasing: L_z = sqrt(gamma_z) * Z  (preserves palindrome)
      2. Decay: L_down = sqrt(gamma_amp * (1 + n_bar)) * sigma_minus
      3. Excitation: L_up = sqrt(gamma_amp * n_bar) * sigma_plus

    At n_bar=0: pure dephasing + zero-temperature amplitude damping.
    At n_bar>0: thermal excitation adds energy from the bath.
    """
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))

    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))

    for k in range(n_qubits):
        gz = gamma_z[k] if hasattr(gamma_z, '__len__') else gamma_z
        ga = gamma_amp[k] if hasattr(gamma_amp, '__len__') else gamma_amp
        nb = n_bar

        # Z-dephasing
        if gz > 0:
            Lk = np.sqrt(gz) * kron_at(Z, k, n_qubits)
            LdL = Lk.conj().T @ Lk
            L += np.kron(Lk.conj(), Lk)
            L -= 0.5 * np.kron(Id, LdL)
            L -= 0.5 * np.kron(LdL.T, Id)

        # Amplitude damping (decay)
        if ga > 0 and (1 + nb) > 0:
            Ld = np.sqrt(ga * (1 + nb)) * kron_at(SM, k, n_qubits)
            LdLd = Ld.conj().T @ Ld
            L += np.kron(Ld.conj(), Ld)
            L -= 0.5 * np.kron(Id, LdLd)
            L -= 0.5 * np.kron(LdLd.T, Id)

        # Thermal excitation
        if ga > 0 and nb > 0:
            Lu = np.sqrt(ga * nb) * kron_at(SP, k, n_qubits)
            LuLu = Lu.conj().T @ Lu
            L += np.kron(Lu.conj(), Lu)
            L -= 0.5 * np.kron(Id, LuLu)
            L -= 0.5 * np.kron(LuLu.T, Id)

    return L


def spectral_metrics(evals):
    """Compute Q-factor and frequency metrics from eigenvalues."""
    osc = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if freq > 1e-10 and rate > 1e-15:
            Q = freq / rate
            osc.append((rate, freq, Q))

    if not osc:
        return {'max_Q': 0, 'n_osc': 0, 'n_freq': 0, 'n_high_Q': 0}

    Qs = [q for _, _, q in osc]
    freqs = sorted(set(round(f, 4) for _, f, _ in osc if f > 1e-8))

    return {
        'max_Q': max(Qs),
        'n_osc': len(osc),
        'n_freq': len(freqs),
        'n_high_Q': sum(1 for q in Qs if q > 5),
    }


def palindrome_reading(evals):
    """Canonical F1 distance about the spectrum's OWN centre, in eps * rho.

    Two repairs in one, both from 2026-08-05.

    THE SCORER. A greedy first-fit percentage inside an absolute tolerance of
    1e-3 stood here. It measured its own matcher: the spectrum is far denser
    than the tolerance, so a rate rarely grabbed its true mirror, and tightening
    the tolerance SATURATES the printed percentage long before it fixes the
    pairing. A 100% therefore carried less information than a 91%. See
    experiments/CONCENTRATOR_MAPPING.md; the check used instead is the repo's
    canonical one, fw.max_f1_pairing_distance.

    THE CENTRE. The old call sites passed the Z-dephasing centre even for rows
    with an amplitude-damping channel running, which is the "centre is
    estimated" caveat experiments/THERMAL_BREAKING.md warned about. It need not
    be estimated at all: a multiset closed under lambda -> 2c - lambda has
    sum(lambda) = n*c, so c = mean(lambda) EXACTLY. There is one candidate
    centre and it is read off the trace, whatever the channels are. A large
    distance at that centre means no centre works.

    Returns (distance_in_eps_rho, centre).
    """
    centre = float(np.mean(evals).real)
    dist, _radius = f1_distance_in_eps(evals, -centre)
    return dist, centre


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "v_effect_thermal.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("V-EFFECT WITH THERMAL BREAKING")
    out("=" * 70)

    J = 1.0
    gamma_z = 0.1      # Z-dephasing (preserves palindrome)
    gamma_amp_base = 0.0  # start without amplitude damping

    n_bar_values = [0.0, 0.01, 0.05, 0.1, 0.2, 0.5,
                    1.0, 2.0, 3.0, 5.0, 10.0]

    # ================================================================
    # Part 1: Pure thermal sweep (no Z-dephasing, only thermal noise)
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 1: Pure thermal noise (gamma_z=0, gamma_amp=0.1)")
    out("Sweep n_bar (thermal occupation)")
    out("=" * 70)

    gamma_amp = 0.1

    out(f"\n{'n_bar':>6} | {'-- N=2 --':>16} | {'-- N=3 --':>16}"
        f" | {'-- N=5 --':>16} | {'V52':>5} {'F1 eps*rho':>10}")
    out(f"{'':>6} | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'':>5} {'':>5}")

    for n_bar in n_bar_values:
        row = {}
        for N in [2, 3, 5]:
            H = build_heisenberg_chain(N, J)
            L = build_liouvillian_thermal(H, 0.0, gamma_amp, n_bar)
            evals = np.linalg.eigvals(L)
            m = spectral_metrics(evals)
            row[f'Q_{N}'] = m['max_Q']
            row[f'f_{N}'] = m['n_freq']
            row[f'h_{N}'] = m['n_high_Q']

            if N == 5:
                row['pal'], row['centre'] = palindrome_reading(evals)

        vg = row['Q_5'] / row['Q_2'] if row['Q_2'] > 0.01 else 0

        out(f"{n_bar:6.2f} | {row['Q_2']:5.1f} {row['f_2']:4d}"
            f" {row['h_2']:4d}"
            f" | {row['Q_3']:5.1f} {row['f_3']:4d}"
            f" {row['h_3']:4d}"
            f" | {row['Q_5']:5.1f} {row['f_5']:4d}"
            f" {row['h_5']:4d}"
            f" | {vg:5.2f} {row.get('pal', 0):10.2e}")

    # ================================================================
    # Part 2: Z-dephasing + thermal sweep
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 2: Z-dephasing (gamma_z=0.1) + thermal (gamma_amp=0.05)")
    out("Sweep n_bar")
    out("=" * 70)

    gamma_z_val = 0.1
    gamma_amp_val = 0.05

    out(f"\n{'n_bar':>6} | {'-- N=2 --':>16} | {'-- N=3 --':>16}"
        f" | {'-- N=5 --':>16} | {'V52':>5} {'F1 eps*rho':>10}")
    out(f"{'':>6} | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'maxQ':>5} {'frq':>4} {'hiQ':>4}"
        f" | {'':>5} {'':>5}")

    for n_bar in n_bar_values:
        row = {}
        for N in [2, 3, 5]:
            H = build_heisenberg_chain(N, J)
            L = build_liouvillian_thermal(
                H, gamma_z_val, gamma_amp_val, n_bar)
            evals = np.linalg.eigvals(L)
            m = spectral_metrics(evals)
            row[f'Q_{N}'] = m['max_Q']
            row[f'f_{N}'] = m['n_freq']
            row[f'h_{N}'] = m['n_high_Q']

            if N == 5:
                row['pal'], row['centre'] = palindrome_reading(evals)

        vg = row['Q_5'] / row['Q_2'] if row['Q_2'] > 0.01 else 0

        out(f"{n_bar:6.2f} | {row['Q_2']:5.1f} {row['f_2']:4d}"
            f" {row['h_2']:4d}"
            f" | {row['Q_3']:5.1f} {row['f_3']:4d}"
            f" {row['h_3']:4d}"
            f" | {row['Q_5']:5.1f} {row['f_5']:4d}"
            f" {row['h_5']:4d}"
            f" | {vg:5.2f} {row.get('pal', 0):10.2e}")

    # ================================================================
    # Part 3: What breaks the 1.81x constant?
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 3: V-Effect gain across noise types")
    out("=" * 70)
    out("\nFixed J=1.0. Compare V-gain = Q(N=5)/Q(N=2)")

    configs = [
        ("Pure Z-dephasing (gz=0.1)",
         {'gamma_z': 0.1, 'gamma_amp': 0.0, 'n_bar': 0}),
        ("Z-deph + cold amp (gz=0.1, ga=0.05, n=0)",
         {'gamma_z': 0.1, 'gamma_amp': 0.05, 'n_bar': 0}),
        ("Z-deph + warm (gz=0.1, ga=0.05, n=0.5)",
         {'gamma_z': 0.1, 'gamma_amp': 0.05, 'n_bar': 0.5}),
        ("Z-deph + hot (gz=0.1, ga=0.05, n=2.0)",
         {'gamma_z': 0.1, 'gamma_amp': 0.05, 'n_bar': 2.0}),
        ("Z-deph + very hot (gz=0.1, ga=0.05, n=5.0)",
         {'gamma_z': 0.1, 'gamma_amp': 0.05, 'n_bar': 5.0}),
        ("Pure thermal cold (ga=0.1, n=0)",
         {'gamma_z': 0.0, 'gamma_amp': 0.1, 'n_bar': 0}),
        ("Pure thermal warm (ga=0.1, n=1.0)",
         {'gamma_z': 0.0, 'gamma_amp': 0.1, 'n_bar': 1.0}),
        ("Pure thermal hot (ga=0.1, n=5.0)",
         {'gamma_z': 0.0, 'gamma_amp': 0.1, 'n_bar': 5.0}),
        ("Only excitation (gz=0, ga=0.1, n=10)",
         {'gamma_z': 0.0, 'gamma_amp': 0.1, 'n_bar': 10.0}),
    ]

    out(f"\n{'Config':>45} | {'Q_N2':>6} {'Q_N5':>6}"
        f" {'V-gain':>7} | {'f_N2':>5} {'f_N5':>5}"
        f" {'f-gain':>7} | {'F1 eps*rho':>10}")

    for label, cfg in configs:
        row = {}
        for N in [2, 5]:
            H = build_heisenberg_chain(N, J)
            L = build_liouvillian_thermal(
                H, cfg['gamma_z'], cfg['gamma_amp'], cfg['n_bar'])
            evals = np.linalg.eigvals(L)
            m = spectral_metrics(evals)
            row[f'Q_{N}'] = m['max_Q']
            row[f'f_{N}'] = m['n_freq']

            if N == 5:
                row['pal'], row['centre'] = palindrome_reading(evals)

        vg_q = row['Q_5'] / row['Q_2'] if row['Q_2'] > 0.01 else 0
        vg_f = row['f_5'] / row['f_2'] if row['f_2'] > 0 else 0

        out(f"{label:>45} | {row['Q_2']:6.1f} {row['Q_5']:6.1f}"
            f" {vg_q:7.2f}x | {row['f_2']:5d} {row['f_5']:5d}"
            f" {vg_f:7.1f}x | {row.get('pal', 0):10.2e}")

    # ================================================================
    # Part 4: Sacrifice + thermal (the full picture)
    # ================================================================
    out("\n" + "=" * 70)
    out("PART 4: Sacrifice profile + thermal noise on N=5")
    out("=" * 70)
    out("\nEdge gamma_z = 0.5, interior gamma_z = 0.01")
    out("Sweep thermal occupation n_bar with gamma_amp = 0.05")

    gz_sac = [0.5, 0.01, 0.01, 0.01, 0.01]  # sacrifice Z-dephasing
    gz_uni = [0.108, 0.108, 0.108, 0.108, 0.108]  # uniform (same total)

    out(f"\n{'n_bar':>6} | {'-- Sacrifice --':>22} | {'-- Uniform --':>22}"
        f" | {'S/U':>5}")
    out(f"{'':>6} | {'maxQ':>6} {'frq':>4} {'hiQ':>4} {'F1eps':>9}"
        f" | {'maxQ':>6} {'frq':>4} {'hiQ':>4} {'F1eps':>9}"
        f" | {'Qrat':>5}")

    for n_bar in n_bar_values:
        H = build_heisenberg_chain(5, J)
        ga = 0.05

        L_sac = build_liouvillian_thermal(H, gz_sac, ga, n_bar)
        ev_sac = np.linalg.eigvals(L_sac)
        m_sac = spectral_metrics(ev_sac)
        pal_sac, _ = palindrome_reading(ev_sac)

        L_uni = build_liouvillian_thermal(H, gz_uni, ga, n_bar)
        ev_uni = np.linalg.eigvals(L_uni)
        m_uni = spectral_metrics(ev_uni)
        pal_uni, _ = palindrome_reading(ev_uni)

        q_rat = (m_sac['max_Q'] / m_uni['max_Q']
                 if m_uni['max_Q'] > 0.01 else 0)

        out(f"{n_bar:6.2f} | {m_sac['max_Q']:6.1f} {m_sac['n_freq']:4d}"
            f" {m_sac['n_high_Q']:4d} {pal_sac:10.2e}"
            f" | {m_uni['max_Q']:6.1f} {m_uni['n_freq']:4d}"
            f" {m_uni['n_high_Q']:4d} {pal_uni:10.2e}"
            f" | {q_rat:5.2f}")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
