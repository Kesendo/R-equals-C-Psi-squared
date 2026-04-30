"""F83 4-Hamiltonian discrimination test: Aer pre-flight on path [4,5,6].

Uses calibration-derived noise from ibm_marrakesh 2026-04-30T16:25Z for the
top-ranked path Q4-Q5-Q6:
  T1 (μs): 206 / 144 / 351
  T2 (μs): 184 / 121 / 151
  RO err:  0.42% / 0.38% / 0.37%
  CZ err:  0.63% (4-5) / 0.29% (5-6)

Compares Aer measured ⟨P_q0 ⊗ I_q1 ⊗ P_q2⟩ to closed-form Trotter prediction
from _f83_signature_predictions.py. If Aer + calibration noise reproduces
Trotter within ~0.05 across discriminating observables, the F83 hardware
test on [4,5,6] is green-light.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw
from framework.pauli import ur_pauli

from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator
from qiskit_aer.noise import (
    NoiseModel, depolarizing_error, thermal_relaxation_error, ReadoutError
)

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass


N = 3
J = 1.0
T_EVAL = 0.8
N_TROTTER = 3
SHOTS = 4096

CATEGORIES = [
    ('truly_unbroken',       [('X', 'X'), ('Y', 'Y')]),
    ('pi2_odd_pure',         [('X', 'Y'), ('Y', 'X')]),
    ('pi2_even_nontruly',    [('Y', 'Z'), ('Z', 'Y')]),
    ('mixed_anti_one_sixth', [('X', 'Y'), ('Y', 'Z')]),
]

# Path [4,5,6] calibration (2026-04-30T16:25Z)
T1_US = [206.115, 144.477, 351.240]
T2_US = [183.978, 120.938, 151.412]
RO_ERR = [0.004150, 0.003784, 0.003662]
CZ_ERR_45 = 0.006327
CZ_ERR_56 = 0.002869

GATE_TIME_1Q_NS = 50
GATE_TIME_2Q_NS = 200


def hamiltonian_to_sparse_pauli(N, bonds, terms, J=1.0):
    pauli_strings, coeffs = [], []
    for (i, j) in bonds:
        for (la, lb) in terms:
            chars = ['I'] * N
            chars[N - 1 - i] = la
            chars[N - 1 - j] = lb
            pauli_strings.append(''.join(chars))
            coeffs.append(J)
    return SparsePauliOp(pauli_strings, coeffs=np.array(coeffs))


def prepare_xneel(N):
    qc = QuantumCircuit(N, name="prep")
    for q in range(N):
        qc.h(q)
        if q % 2 == 1:
            qc.z(q)
    return qc


def basis_change(qc, q, basis):
    if basis == 'X':
        qc.h(q)
    elif basis == 'Y':
        qc.sdg(q)
        qc.h(q)


def build_circuit(N, prep, ham_op, t, n_trot, b0, b2):
    qc = QuantumCircuit(N, 2)
    qc.compose(prep, inplace=True)
    dt = t / n_trot
    evo_gate = PauliEvolutionGate(ham_op, time=dt)
    for _ in range(n_trot):
        qc.append(evo_gate, range(N))
    basis_change(qc, 0, b0)
    basis_change(qc, 2, b2)
    qc.measure(0, 0)
    qc.measure(2, 1)
    return qc


def build_calibrated_noise_model(N=3):
    """Per-qubit noise from [4,5,6] calibration."""
    nm = NoiseModel()

    # Per-qubit thermal relaxation on 1q gates (mapped to logical qubits 0,1,2 = physical 4,5,6)
    for q in range(N):
        T1 = T1_US[q] * 1e-6
        T2 = T2_US[q] * 1e-6
        therm = thermal_relaxation_error(T1, T2, GATE_TIME_1Q_NS * 1e-9)
        depol = depolarizing_error(0.0005, 1)  # 1q error ~0.05%
        err = therm.compose(depol)
        nm.add_quantum_error(err, ['x', 'sx', 'rz', 'h', 's', 'sdg', 'z', 'y'], [q])

    # 2q noise on edges: (0,1) ↔ Q4-Q5, (1,2) ↔ Q5-Q6
    for (qa, qb, cz_err) in [(0, 1, CZ_ERR_45), (1, 2, CZ_ERR_56)]:
        T1a = T1_US[qa] * 1e-6
        T2a = T2_US[qa] * 1e-6
        T1b = T1_US[qb] * 1e-6
        T2b = T2_US[qb] * 1e-6
        therm_a = thermal_relaxation_error(T1a, T2a, GATE_TIME_2Q_NS * 1e-9)
        therm_b = thermal_relaxation_error(T1b, T2b, GATE_TIME_2Q_NS * 1e-9)
        therm = therm_a.tensor(therm_b)
        depol = depolarizing_error(cz_err, 2)
        err = therm.compose(depol)
        nm.add_quantum_error(err, ['cx', 'cz'], [qa, qb])
        nm.add_quantum_error(err, ['cx', 'cz'], [qb, qa])

    # Per-qubit readout error
    for q in range(N):
        p = RO_ERR[q]
        ro = ReadoutError([[1 - p, p], [p, 1 - p]])
        nm.add_readout_error(ro, [q])

    return nm


def measure_pauli_q0_q2(N, prep, ham_op, t, n_trot, sim, n_shots):
    """Run 9-basis tomography on (q0, q2). Return dict {(p0,p2): ⟨P⟩}."""
    bases = ['X', 'Y', 'Z']
    counts_by_basis = {}
    for b0 in bases:
        for b2 in bases:
            qc = build_circuit(N, prep, ham_op, t, n_trot, b0, b2)
            qc_t = transpile(qc, sim, optimization_level=1)
            result = sim.run(qc_t, shots=n_shots).result()
            counts_by_basis[(b0, b2)] = result.get_counts()

    expectations = {('I', 'I'): 1.0}
    for b0 in bases:
        for b2 in bases:
            counts = counts_by_basis[(b0, b2)]
            tot = sum(counts.values())
            n00 = counts.get('00', 0); n01 = counts.get('01', 0)
            n10 = counts.get('10', 0); n11 = counts.get('11', 0)
            expectations[(b0, b2)] = (n00 + n11 - n01 - n10) / tot
        counts_z = counts_by_basis[(b0, 'Z')]
        tot = sum(counts_z.values())
        n0 = counts_z.get('00', 0) + counts_z.get('10', 0)
        n1 = counts_z.get('01', 0) + counts_z.get('11', 0)
        expectations[(b0, 'I')] = (n0 - n1) / tot
    for b2 in bases:
        counts_z = counts_by_basis[('Z', b2)]
        tot = sum(counts_z.values())
        n0 = counts_z.get('00', 0) + counts_z.get('01', 0)
        n1 = counts_z.get('10', 0) + counts_z.get('11', 0)
        expectations[('I', b2)] = (n0 - n1) / tot
    return expectations


def trotter_prediction(N, terms, t, n_trot, gamma_z=0.1):
    """Closed-form Trotter+Lindblad reference (from _f83_signature_predictions)."""
    from scipy.linalg import expm
    from framework.lindblad import lindbladian_z_dephasing
    from framework.pauli import _build_bilinear, site_op

    chain = fw.ChainSystem(N=N)
    bonds = chain.bonds
    bilinear = [(a, b, J) for (a, b) in terms]
    delta_t = t / n_trot
    Id = np.eye(2 ** N, dtype=complex)
    U_step = np.eye(2 ** N, dtype=complex)
    for (P, Q, c) in bilinear:
        for (l, m) in bonds:
            ops = [ur_pauli('I')] * N
            ops[l] = ur_pauli(P); ops[m] = ur_pauli(Q)
            op_full = ops[0]
            for op in ops[1:]:
                op_full = np.kron(op_full, op)
            U_step = expm(-1j * c * delta_t * op_full) @ U_step

    L_deph = np.zeros((4 ** N, 4 ** N), dtype=complex)
    for l in range(N):
        Zl = site_op(N, l, 'Z')
        L_deph += gamma_z * (np.kron(Zl, Zl.conj()) - np.kron(Id, Id))
    M_deph = expm(L_deph * delta_t)

    ket_p = np.array([1, 1], dtype=complex) / np.sqrt(2)
    ket_m = np.array([1, -1], dtype=complex) / np.sqrt(2)
    psi0 = np.kron(np.kron(ket_p, ket_m), ket_p)
    rho = np.outer(psi0, psi0.conj())
    for _ in range(n_trot):
        rho = U_step @ rho @ U_step.conj().T
        rho_v = rho.flatten('F')
        rho = (M_deph @ rho_v).reshape((2**N, 2**N), order='F')

    out = {}
    for p0 in 'IXYZ':
        for p2 in 'IXYZ':
            obs = np.kron(np.kron(ur_pauli(p0), np.eye(2)), ur_pauli(p2))
            out[(p0, p2)] = float(np.real(np.trace(rho @ obs)))
    return out


def main():
    print("=" * 78)
    print("F83 4-Hamiltonian Aer Pre-Flight (path [4,5,6] calibration)")
    print("=" * 78)
    print(f"  Calibration: ibm_marrakesh 2026-04-30T16:25Z, top-ranked path")
    print(f"  T1 (μs): {T1_US}")
    print(f"  T2 (μs): {T2_US}")
    print(f"  RO err: {[f'{x*100:.2f}%' for x in RO_ERR]}")
    print(f"  CZ err: 4-5 = {CZ_ERR_45*100:.3f}%, 5-6 = {CZ_ERR_56*100:.3f}%")
    print()

    nm = build_calibrated_noise_model(N=N)
    sim = AerSimulator(noise_model=nm, method='density_matrix')
    bonds = [(0, 1), (1, 2)]
    prep = prepare_xneel(N)

    KEY_PAULIS = [('Z', 'Z'), ('X', 'I'), ('Z', 'I'), ('X', 'Z'), ('I', 'X'), ('X', 'X')]

    print(f"{'Category':<24} {'Pauli':<7} {'Trotter':>10} {'Aer':>10} {'|Δ|':>8}")
    print('-' * 70)

    aer_results = {}
    trot_results = {}
    max_dev = 0.0
    for cat, terms in CATEGORIES:
        ham_op = hamiltonian_to_sparse_pauli(N, bonds, terms, J=J)
        aer_exp = measure_pauli_q0_q2(N, prep, ham_op, T_EVAL, N_TROTTER, sim, SHOTS)
        trot_exp = trotter_prediction(N, terms, T_EVAL, N_TROTTER)
        aer_results[cat] = aer_exp
        trot_results[cat] = trot_exp
        for (p0, p2) in KEY_PAULIS:
            t_v = trot_exp[(p0, p2)]
            a_v = aer_exp[(p0, p2)]
            dev = abs(a_v - t_v)
            max_dev = max(max_dev, dev)
            print(f"{cat:<24} {p0},{p2:<5} {t_v:>+10.4f} {a_v:>+10.4f} {dev:>8.4f}")
        print()

    print(f"Max |Aer − Trotter| over {len(KEY_PAULIS)} key observables × {len(CATEGORIES)} categories: {max_dev:.4f}")
    print()

    print("=" * 78)
    print("Discriminability under noise: pairwise category separation per Pauli")
    print("=" * 78)
    print()
    cats = [c for c, _ in CATEGORIES]
    sigma = 1.0 / np.sqrt(SHOTS)  # ~0.0156 per ⟨P⟩ at 4096 shots
    print(f"  Statistical error per ⟨P⟩ at {SHOTS} shots: ±{sigma:.4f}")
    print()
    print(f"  Largest pairwise Aer separation per Pauli (σ-units):")
    print(f"  {'Pauli':<7} {'best pair':<40} {'Aer Δ':>8} {'σ':>6}")
    print('  ' + '-' * 70)
    for (p0, p2) in KEY_PAULIS:
        vals = [(c, aer_results[c][(p0, p2)]) for c in cats]
        # Pairwise differences
        max_pair = max(((c1, v1, c2, v2) for c1, v1 in vals for c2, v2 in vals if c1 < c2),
                       key=lambda x: abs(x[1] - x[3]))
        c1, v1, c2, v2 = max_pair
        delta = abs(v1 - v2)
        print(f"  {p0},{p2:<5} {c1[:18]:<18} vs {c2[:18]:<18} {delta:>8.4f} {delta/sigma:>6.1f}")
    print()

    print("=" * 78)
    print("Verdict")
    print("=" * 78)
    print()
    # Compute minimum over all 6 pairwise category distinctions for each Pauli
    # Each category should have at least one Pauli where its distance to all 3 others exceeds 5σ
    print(f"  Per-category 'discriminating Pauli' check (Δ > 5σ ≈ {5*sigma:.3f} from all others):")
    print()
    fingerprintable = 0
    for c0 in cats:
        for (p0, p2) in KEY_PAULIS:
            v0 = aer_results[c0][(p0, p2)]
            others = [aer_results[c][(p0, p2)] for c in cats if c != c0]
            min_dist = min(abs(v0 - v) for v in others)
            if min_dist > 5 * sigma:
                print(f"  ✓ {c0:<24} Pauli {p0},{p2}: Aer={v0:+.3f}, min Δ to other cats = {min_dist:.3f} ({min_dist/sigma:.1f}σ)")
                fingerprintable += 1
                break
        else:
            print(f"  ✗ {c0:<24} no single Pauli > 5σ from all others (joint pattern still works)")
    print()
    if fingerprintable == 4:
        print(f"  GREEN: All 4 categories have at least one >5σ-discriminating Pauli on [4,5,6].")
    elif fingerprintable >= 2:
        print(f"  YELLOW: {fingerprintable}/4 cats have a single-Pauli fingerprint; remaining need joint pattern.")
    else:
        print(f"  RED: Fewer than half the categories are single-Pauli discriminable.")
    print()
    print(f"  X,Z anchor reproducibility: pi2_odd_pure Trotter={trot_results['pi2_odd_pure'][('X','Z')]:.4f}, ")
    print(f"  Aer={aer_results['pi2_odd_pure'][('X','Z')]:.4f} (April 26 hardware was -0.7217).")


if __name__ == "__main__":
    main()
