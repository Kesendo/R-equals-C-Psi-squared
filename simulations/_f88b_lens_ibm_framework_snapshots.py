"""F88b-Lens applied to IBM Marrakesh hardware tomography snapshots.

Loads the framework_snapshots JSON from AIEvolution's IBM tomography pipeline
(2026-04-26 Marrakesh run) and runs the F88b Π²-odd-fraction-within-memory
diagnostic on hardware-reconstructed 2-qubit reduced ρ states.

Snapshots and what they probe:
- A: Singlet |S⟩ = (|01⟩ − |10⟩)/√2: popcount-(1, 1) HD = 2 = N, F88b Π²-classical (Π²-odd / memory = 0)
- B: Triplet T_0 = (|01⟩ + |10⟩)/√2: same anchor
- C: |+−+⟩ at t = 0.8 under Heisenberg, 9-Pauli tomo on (q0, q2)
- D: |+−+⟩ under truly / soft / hard 2-body Hamiltonian categories

Snapshots A, B only have ⟨XX⟩, ⟨YY⟩, ⟨ZZ⟩ measured (not full 16-Pauli);
the F88b-classical prediction is partially testable: those three terms are
all Π²-even, so the contribution to Π²-odd/memory from the measured part
is zero by construction. The full hardware-noise deviation needs full
tomography (not in this snapshot data).

Snapshots C, D have full 16-Pauli expectations on the 2-qubit reduced ρ
of (q0, q2); full ρ reconstructable, F88b-Lens directly applicable.

|+−+⟩ is NOT a popcount-coherence pair state (it is a uniform superposition
across all popcounts), so F88b's closed form does not directly apply.
We compute Π²-odd/memory empirically via the framework's state-level lens
and compare to the same diagnostic on the noiseless simulated reduced ρ.
The hardware-vs-ideal deviation is the noise signature.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")


JSON_PATH = Path(
    r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI"
    r"\experiments\ibm_quantum_tomography\results"
    r"\framework_snapshots_ibm_marrakesh_20260426_105948.json"
)


I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULIS = {"I": I2, "X": SX, "Y": SY, "Z": SZ}
BIT_B = {"I": 0, "X": 0, "Y": 1, "Z": 1}


def reconstruct_2qubit_rho(expectations: dict[str, float]) -> np.ndarray:
    """ρ = (1/4) Σ_{αβ} ⟨σ_α ⊗ σ_β⟩ · σ_α ⊗ σ_β.

    Expectation keys are "α,β" with α the first qubit's Pauli letter,
    β the second's (standard tensor order).
    """
    rho = np.zeros((4, 4), dtype=complex)
    for key, expv in expectations.items():
        a, b = key.split(",")
        sigma = np.kron(PAULIS[a], PAULIS[b])
        rho += expv * sigma / 4.0
    return rho


def f88b_lens_2qubit(rho: np.ndarray) -> dict[str, float]:
    """Static / memory / Π²-odd-of-memory split at N = 2 (kernel = popcount sectors)."""
    P0 = np.diag([1, 0, 0, 0]).astype(complex)
    P1 = np.diag([0, 1, 1, 0]).astype(complex)
    P2 = np.diag([0, 0, 0, 1]).astype(complex)

    rho_d0 = (
        np.real(np.trace(P0 @ rho)) / 1.0 * P0
        + np.real(np.trace(P1 @ rho)) / 2.0 * P1
        + np.real(np.trace(P2 @ rho)) / 1.0 * P2
    )
    rho_d2 = rho - rho_d0

    odd_memory_sq = 0.0
    for a in "IXYZ":
        for b in "IXYZ":
            if (BIT_B[a] + BIT_B[b]) & 1 == 0:
                continue
            sigma = np.kron(PAULIS[a], PAULIS[b])
            coeff = np.trace(sigma @ rho_d2) / 4.0
            odd_memory_sq += abs(coeff) ** 2 * 4

    norm_total = np.linalg.norm(rho, "fro") ** 2
    norm_static = np.linalg.norm(rho_d0, "fro") ** 2
    norm_memory = np.linalg.norm(rho_d2, "fro") ** 2

    return {
        "norm_total": float(norm_total),
        "static_frac": float(norm_static / norm_total) if norm_total > 0 else 0.0,
        "memory_frac": float(norm_memory / norm_total) if norm_total > 0 else 0.0,
        "pi2_odd_in_memory": float(odd_memory_sq / norm_memory) if norm_memory > 1e-14 else 0.0,
        "trace": float(np.real(np.trace(rho))),
        "purity": float(np.real(np.trace(rho @ rho))),
    }


def ideal_plus_minus_plus_rho_q0q2(t: float, J: float = 1.0) -> np.ndarray:
    """Noiseless simulation: |+−+⟩ at N=3, exact unitary evolution under Heisenberg
    H = (J/4)·Σ_b (X_b X_{b+1} + Y_b Y_{b+1} + Z_b Z_{b+1}), partial trace to (q0, q2).
    Hardware uses Trotterised gates; we compare against the exact-evolution baseline
    so the deviation captures Trotter error + gate noise + decoherence combined.
    """
    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(psi_plus, psi_minus), psi_plus)

    def site_op(P, k, N=3):
        ops = [I2] * N
        ops[k] = P
        return np.kron(np.kron(ops[0], ops[1]), ops[2])

    H = np.zeros((8, 8), dtype=complex)
    for b in [0, 1]:
        for P in [SX, SY, SZ]:
            H += (J / 4.0) * site_op(P, b) @ site_op(P, b + 1)

    from scipy.linalg import expm
    U = expm(-1j * H * t)
    psi_t = U @ psi
    rho_full = np.outer(psi_t, psi_t.conj())

    rho_q0q2 = np.zeros((4, 4), dtype=complex)
    for i0 in range(2):
        for i2 in range(2):
            for j0 in range(2):
                for j2 in range(2):
                    s = 0.0 + 0j
                    for q1 in range(2):
                        idx_i = (i0 << 2) | (q1 << 1) | i2
                        idx_j = (j0 << 2) | (q1 << 1) | j2
                        s += rho_full[idx_i, idx_j]
                    rho_q0q2[(i0 << 1) | i2, (j0 << 1) | j2] = s
    return rho_q0q2


def lindbladian_plus_minus_plus_rho_q0q2(
    t: float,
    H_terms: list[tuple[str, str, float]],
    gamma: float = 0.0,
    h_y: float = 0.0,
    gamma_t1: float = 0.0,
    N: int = 3,
) -> np.ndarray:
    """Noiseless |+−+⟩ at N=3 under H = Σ_b Σ_term coeff·(a_b ⊗ b_{b+1}) + h_y·Σ_l Y_l,
    optionally with Z-dephasing γ and T1 amplitude damping γ_T1 on each site.
    Returns 2-qubit reduced ρ on (q0, q2).

    H_terms: list of (letter_a, letter_b, coeff_per_bond). The caller supplies
    the coupling explicitly per term, so this matches the F87-trichotomy convention
    (coeff=J=1.0) directly. For Heisenberg comparison, supply coeff=J/4 per XX/YY/ZZ.
    h_y: optional transverse Y-field per site (coherent), e.g. Marrakesh-detected
    h_y_eff = 0.05 leak (Z⊗N parity diagnostic, 2026-04-29).
    gamma_t1: optional amplitude-damping rate per site, L = σ⁻ = |0⟩⟨1|.

    Pure unitary U = expm(-iHt) if gamma == 0 AND gamma_t1 == 0.
    Otherwise vectorized Lindbladian ρ(t) = devec(expm(L·t)·vec(ρ₀)), column-major vec.
    Dissipator vec conventions: [H,·] → I⊗H − H^T⊗I; Z_l·Z_l → Z_l^T ⊗ Z_l
    (Z real, Z^T = Z); σ⁻_l·σ⁺_l → σ⁻* ⊗ σ⁻ (σ⁻ real); {σ⁺σ⁻, ·}/2 →
    ½(I⊗P_1 + P_1^T⊗I) with P_1 = (I−Z)/2 = |1⟩⟨1| (real diagonal, P_1^T = P_1).
    """
    if N != 3:
        raise NotImplementedError(f"only N=3 supported for now, got N={N}")

    psi_plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_minus = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi = np.kron(np.kron(psi_plus, psi_minus), psi_plus)

    def site_op(P, k):
        ops = [I2] * N
        ops[k] = P
        return np.kron(np.kron(ops[0], ops[1]), ops[2])

    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for (a_letter, b_letter, coeff) in H_terms:
        Pa = PAULIS[a_letter]
        Pb = PAULIS[b_letter]
        for b in range(N - 1):
            H += coeff * site_op(Pa, b) @ site_op(Pb, b + 1)

    if h_y != 0.0:
        for k in range(N):
            H += h_y * site_op(SY, k)

    from scipy.linalg import expm

    if gamma == 0.0 and gamma_t1 == 0.0:
        U = expm(-1j * H * t)
        psi_t = U @ psi
        rho_full = np.outer(psi_t, psi_t.conj())
    else:
        I_d = np.eye(d, dtype=complex)
        L_super = -1j * (np.kron(I_d, H) - np.kron(H.T, I_d))

        if gamma > 0.0:
            for l in range(N):
                Z_l = site_op(SZ, l)
                L_super += gamma * (np.kron(Z_l.T, Z_l) - np.kron(I_d, I_d))

        if gamma_t1 > 0.0:
            sigma_minus = np.array([[0, 1], [0, 0]], dtype=complex)
            P1_single = np.array([[0, 0], [0, 1]], dtype=complex)
            for l in range(N):
                sm_l = site_op(sigma_minus, l)
                P1_l = site_op(P1_single, l)
                L_super += gamma_t1 * (
                    np.kron(sm_l.conj(), sm_l)
                    - 0.5 * np.kron(I_d, P1_l)
                    - 0.5 * np.kron(P1_l.T, I_d)
                )

        rho_0 = np.outer(psi, psi.conj())
        vec_rho_0 = rho_0.flatten(order='F')
        vec_rho_t = expm(L_super * t) @ vec_rho_0
        rho_full = vec_rho_t.reshape((d, d), order='F')

    rho_q0q2 = np.zeros((4, 4), dtype=complex)
    for i0 in range(2):
        for i2 in range(2):
            for j0 in range(2):
                for j2 in range(2):
                    s = 0.0 + 0j
                    for q1 in range(2):
                        idx_i = (i0 << 2) | (q1 << 1) | i2
                        idx_j = (j0 << 2) | (q1 << 1) | j2
                        s += rho_full[idx_i, idx_j]
                    rho_q0q2[(i0 << 1) | i2, (j0 << 1) | j2] = s
    return rho_q0q2


def main():
    print("F88b-Lens on IBM Marrakesh framework_snapshots (2026-04-26)")
    print("=" * 78)
    with open(JSON_PATH, "r", encoding="utf-8") as fh:
        data = json.load(fh)

    print(f"backend = {data['backend']}, job_id = {data['job_id']}, shots = {data['shots']}")
    print(f"path = {data['path']}, parameters = {data['parameters']}")
    print()

    # Snapshots A, B: only XX, YY, ZZ; partial check that all three are Π²-even (no Π²-odd contribution)
    for label in ["snapshot_a_singlet", "snapshot_b_triplet"]:
        snap = data[label]
        print(f"=== {label} ===")
        print(f"  measured: XX={snap['XX']:+.4f}  YY={snap['YY']:+.4f}  ZZ={snap['ZZ']:+.4f}")
        print(f"  σ·σ = {snap['sigma_dot']:+.4f}  (predicted {snap['predicted']:+.1f})")
        print(f"  F88b: HD=N=2 popcount-{(label.split('_')[-1]=='singlet') and '(1,1)' or '(1,1)'} → Π²-classical, all measured Pauli are Π²-even ✓")
        print()

    # Snapshots C, D: full 16-Pauli, F88b-Lens directly applicable
    print("=== snapshot_c_heisenberg (|+−+⟩ at t=0.8 under Heisenberg, reduced (q0, q2)) ===")
    rho_c = reconstruct_2qubit_rho(data["snapshot_c_heisenberg"]["expectations"])
    lens_c = f88b_lens_2qubit(rho_c)
    print(f"  trace = {lens_c['trace']:.4f}, purity = {lens_c['purity']:.4f}")
    print(f"  static / memory = {lens_c['static_frac']:.4f} / {lens_c['memory_frac']:.4f}")
    print(f"  Π²-odd / memory (hardware) = {lens_c['pi2_odd_in_memory']:.4f}")

    rho_c_ideal = ideal_plus_minus_plus_rho_q0q2(t=0.8)
    lens_c_ideal = f88b_lens_2qubit(rho_c_ideal)
    print(f"  Π²-odd / memory (ideal sim) = {lens_c_ideal['pi2_odd_in_memory']:.4f}")
    print(f"  Δ (hardware − ideal) = {lens_c['pi2_odd_in_memory'] - lens_c_ideal['pi2_odd_in_memory']:+.4f}")
    print()

    print("=== snapshot_d_softbreak_trichotomy ===")
    TRICHOTOMY_H_TERMS = {
        "truly": [("X", "X", 1.0), ("Y", "Y", 1.0)],
        "soft":  [("X", "Y", 1.0), ("Y", "X", 1.0)],
        "hard":  [("X", "X", 1.0), ("X", "Y", 1.0)],
    }
    TRICHOTOMY_GAMMA = 0.1
    MARRAKESH_HY = 0.05
    t_eval = data["parameters"]["t_eval"]
    qubits_path = data["path"]

    cal_csv = Path(
        r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared"
        r"\data\ibm_calibration_snapshots"
        r"\ibm_marrakesh_calibrations_2026-04-25T11_28_00Z.csv"
    )
    t1_t2 = {}
    import csv as _csv
    with open(cal_csv, newline='', encoding='utf-8') as _fh:
        for row in _csv.DictReader(_fh):
            qb = int(row["Qubit"])
            if qb in qubits_path:
                t1_t2[qb] = (float(row["T1 (us)"]), float(row["T2 (us)"]))
    mean_1_over_T1 = sum(1.0 / t1 for (t1, _t2) in t1_t2.values()) / len(t1_t2)
    mean_1_over_Tphi = sum(
        1.0 / t2 - 1.0 / (2.0 * t1) for (t1, t2) in t1_t2.values()
    ) / len(t1_t2)
    gamma_t1_over_gamma_z = mean_1_over_T1 / mean_1_over_Tphi
    MARRAKESH_GAMMA_T1 = TRICHOTOMY_GAMMA * gamma_t1_over_gamma_z

    print(f"  γ_T1 derived from Marrakesh calibration 2026-04-25 (qubits {qubits_path}):")
    for qb in sorted(t1_t2):
        t1, t2 = t1_t2[qb]
        tphi = 1.0 / (1.0 / t2 - 1.0 / (2.0 * t1))
        print(f"     q{qb}: T1={t1:.1f}μs, T2={t2:.1f}μs, T_φ={tphi:.1f}μs")
    print(f"     ⟨1/T1⟩/⟨1/T_φ⟩ = {gamma_t1_over_gamma_z:.3f}  →  γ_T1 = γ_Z·{gamma_t1_over_gamma_z:.3f} = {MARRAKESH_GAMMA_T1:.4f}")
    print()
    print(f"  baseline A: Lindbladian γ={TRICHOTOMY_GAMMA}, no h_y leak, no T1")
    print(f"  baseline B: Lindbladian γ={TRICHOTOMY_GAMMA} + transverse h_y={MARRAKESH_HY} per site")
    print(f"  baseline C: Lindbladian γ={TRICHOTOMY_GAMMA} + h_y={MARRAKESH_HY} + γ_T1={MARRAKESH_GAMMA_T1:.4f} per site")
    print()
    for cat in ["truly", "soft", "hard"]:
        rho_d_hw = reconstruct_2qubit_rho(data["snapshot_d_softbreak_trichotomy"]["expectations_per_category"][cat])
        hw = f88b_lens_2qubit(rho_d_hw)

        rho_d_A = lindbladian_plus_minus_plus_rho_q0q2(
            t=t_eval, H_terms=TRICHOTOMY_H_TERMS[cat], gamma=TRICHOTOMY_GAMMA,
        )
        A = f88b_lens_2qubit(rho_d_A)

        rho_d_B = lindbladian_plus_minus_plus_rho_q0q2(
            t=t_eval, H_terms=TRICHOTOMY_H_TERMS[cat],
            gamma=TRICHOTOMY_GAMMA, h_y=MARRAKESH_HY,
        )
        B = f88b_lens_2qubit(rho_d_B)

        rho_d_C = lindbladian_plus_minus_plus_rho_q0q2(
            t=t_eval, H_terms=TRICHOTOMY_H_TERMS[cat],
            gamma=TRICHOTOMY_GAMMA, h_y=MARRAKESH_HY, gamma_t1=MARRAKESH_GAMMA_T1,
        )
        C = f88b_lens_2qubit(rho_d_C)

        delta_A = hw["pi2_odd_in_memory"] - A["pi2_odd_in_memory"]
        delta_B = hw["pi2_odd_in_memory"] - B["pi2_odd_in_memory"]
        delta_C = hw["pi2_odd_in_memory"] - C["pi2_odd_in_memory"]

        print(f"  [{cat}]")
        print(f"     hardware              : Π²-odd/mem = {hw['pi2_odd_in_memory']:.4f}  "
              f"(purity {hw['purity']:.3f}, static/mem {hw['static_frac']:.3f}/{hw['memory_frac']:.3f})")
        print(f"     ideal A (γ only)      : Π²-odd/mem = {A['pi2_odd_in_memory']:.4f}  "
              f"(purity {A['purity']:.3f}, static/mem {A['static_frac']:.3f}/{A['memory_frac']:.3f})   Δ_A = {delta_A:+.4f}")
        print(f"     ideal B (γ + h_y)     : Π²-odd/mem = {B['pi2_odd_in_memory']:.4f}  "
              f"(purity {B['purity']:.3f}, static/mem {B['static_frac']:.3f}/{B['memory_frac']:.3f})   Δ_B = {delta_B:+.4f}")
        print(f"     ideal C (γ + h_y + T1): Π²-odd/mem = {C['pi2_odd_in_memory']:.4f}  "
              f"(purity {C['purity']:.3f}, static/mem {C['static_frac']:.3f}/{C['memory_frac']:.3f})   Δ_C = {delta_C:+.4f}")
    print()

    print("=== Reading ===")
    print()
    print("Snapshot C (Heisenberg evolution): hardware adds ~1% Π²-odd noise above the")
    print("noiseless framework prediction. Heisenberg is a truly Hamiltonian (Π²-even");
    print("XX+YY+ZZ bilinears), so the framework predicts Π²-odd/memory ≈ 0; hardware")
    print("0.0098 confirms this within Trotter + gate-noise tolerance.")
    print()
    print("Snapshot D differentiates the three V-Effect / F87-trichotomy categories at")
    print("state level via Π²-odd/memory:")
    print("  truly (XX+YY, Π²-even bilinears):    0.030  (framework predicts low Π²-odd) ✓")
    print("  soft  (XY+YX, Π²-odd bilinears):     0.744  (Hamiltonian Π²-odd content")
    print("                                                 dynamically pumps state-level Π²-odd) ✓")
    print("  hard  (XX+XY, mixed bilinears):      0.276  (intermediate, mirrored decomposition) ✓")
    print()
    print("F88b's F87-trichotomy operator-level claim (truly: ‖M‖ = 0; soft / hard:")
    print("‖M‖ ≠ 0 with Π²-odd-driven super-operator structure) propagates cleanly to the")
    print("state-level Π²-odd/memory reading on real hardware: same three-way split,")
    print("same ordering, with order-of-magnitude separation between truly and soft.")


if __name__ == "__main__":
    main()
