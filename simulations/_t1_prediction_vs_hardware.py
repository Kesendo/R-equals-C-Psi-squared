#!/usr/bin/env python3
"""EQ-030 retrospective verification: framework T1-prediction vs Snapshot D hardware.

The T1-extended Lindbladian (`lindbladian_z_plus_t1`) predicts an asymmetry:
  - truly's Π-protection is fragile to T1 (almost everything leaks)
  - soft's Π-protection is robust (only specific observables leak)
  - hard collapses completely under T1

For Snapshot D's 9-Pauli tomography on (q0, q2), the relevant Pauli strings
are of the form A_q0 ⊗ I_q1 ⊗ B_q2 with A, B ∈ {X, Y, Z} — 9 strings per
category per backend.

This script:
  1. Loads the Snapshot D JSON results from today's hardware runs
  2. Extracts ⟨A_q0 B_q2⟩ for each of the 9 (A, B) per category per backend
  3. Computes framework's pi_protected_observables for truly/soft/hard
     under (a) pure Z-dephasing and (b) Z + small T1
  4. Matches each (A, B) to the corresponding N=3 Pauli string and reports
     predicted protected/active vs measured magnitude
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw

RESULTS_DIR = Path(r"D:\Entwicklung\Projekte\.NET Projekte\AIEvolution\AIEvolution.UI\experiments\ibm_quantum_tomography\results")


def get_protected_set(H, gamma_l, gamma_t1_l, rho_0, N):
    """Return set of Pauli labels (e.g., 'XIZ') that are Π-protected."""
    if any(g != 0 for g in gamma_t1_l):
        L = fw.lindbladian_z_plus_t1(H, gamma_l, gamma_t1_l)
    else:
        L = fw.lindbladian_z_dephasing(H, gamma_l)

    M_basis = fw._vec_to_pauli_basis_transform(N)
    L_pauli = (M_basis.conj().T @ L @ M_basis) / (2 ** N)
    evals, V = np.linalg.eig(L_pauli)
    Vinv = np.linalg.inv(V)
    rho_pauli = fw.pauli_basis_vector(rho_0, N)
    c = Vinv @ rho_pauli

    cluster_tol = 1e-8
    n_eig = len(evals)
    used = np.zeros(n_eig, dtype=bool)
    clusters = []
    for i in range(n_eig):
        if used[i]:
            continue
        cl = [i]
        used[i] = True
        for j in range(i + 1, n_eig):
            if not used[j] and abs(evals[j] - evals[i]) < cluster_tol:
                cl.append(j)
                used[j] = True
        clusters.append(cl)

    threshold = 1e-9
    protected_labels = set()
    for alpha in range(1, 4 ** N):
        max_S = 0.0
        for cl in clusters:
            S = sum(V[alpha, k] * c[k] for k in cl)
            max_S = max(max_S, abs(S))
        if max_S < threshold:
            label = ''.join(fw.PAULI_LABELS[idx] for idx in fw._k_to_indices(alpha, N))
            protected_labels.add(label)
    return protected_labels


def main():
    N = 3
    GAMMA_DEPH = 0.1
    J = 1.0
    bonds = [(i, i + 1) for i in range(N - 1)]

    # |+−+⟩ initial state
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = np.kron(plus, np.kron(minus, plus))
    rho_0 = np.outer(psi, psi.conj())

    cases = {
        'truly': [('X', 'X', J), ('Y', 'Y', J)],
        'soft':  [('X', 'Y', J), ('Y', 'X', J)],
        'hard':  [('X', 'X', J), ('X', 'Y', J)],
    }

    # T1 ratios to consider: 0 (pure Z) and small T1 (representative of hardware)
    t1_ratios = [0.0, 0.1]

    # Load hardware Snapshot D data from all framework_snapshots files
    hw_results = {}  # backend -> {category -> {(A,B) -> measured}}
    for json_path in sorted(RESULTS_DIR.glob("framework_snapshots_*.json")):
        if not json_path.name.startswith("framework_snapshots_ibm_"):
            continue
        with open(json_path, encoding='utf-8') as f:
            data = json.load(f)
        if 'snapshot_d_softbreak_trichotomy' not in data:
            continue
        backend = data.get('backend', json_path.stem)
        # Use timestamp to differentiate multiple runs on same backend
        ts = json_path.stem.split('_')[-1]
        key = f"{backend}_{ts}"
        d_data = data['snapshot_d_softbreak_trichotomy']['expectations_per_category']
        hw_results[key] = d_data

    print(f"Loaded {len(hw_results)} hardware snapshot-D runs")
    print()

    # Compute framework's protected set for each category at each T1 ratio
    framework_predictions = {}
    for cat_name, terms in cases.items():
        H = fw._build_bilinear(N, bonds, terms)
        framework_predictions[cat_name] = {}
        for ratio in t1_ratios:
            gamma_t1 = ratio * GAMMA_DEPH
            prot = get_protected_set(
                H,
                gamma_l=[GAMMA_DEPH] * N,
                gamma_t1_l=[gamma_t1] * N,
                rho_0=rho_0,
                N=N,
            )
            framework_predictions[cat_name][ratio] = prot

    # The 9 Pauli strings of form (A, I, B) with A, B in {X, Y, Z}
    # In framework's convention: site 0 first, site 1 middle (I), site 2 last
    pauli_3strings = []
    for A in 'XYZ':
        for B in 'XYZ':
            label = f"{A}I{B}"
            pauli_3strings.append((A, B, label))

    # Mapping from framework's (A, I, B) label to qiskit's (b0, b2)
    # In qiskit's little-endian convention: q0 is rightmost in tensor.
    # Snapshot D measures qubit 0 (cr[0]) and qubit 2 (cr[1]).
    # The expectation labelled (A_qiskit_q0, B_qiskit_q2) corresponds to
    # the operator B ⊗ I ⊗ A in framework's tensor (left-to-right = q2, q1, q0).
    # So framework's "AIB" label for site (0=q2, 1=q1, 2=q0) means
    # ⟨A_q2 I_q1 B_q0⟩ = framework's AIB → qiskit measurement of (B, A).
    # Equivalently: framework's "AIB" matches qiskit's b0=B, b2=A.

    print(f"Framework's Π-protected set for each category:")
    print(f"  (only 'A I B' Pauli strings, comparing pure-Z vs +T1)")
    print()
    for cat_name in ['truly', 'soft', 'hard']:
        prot_pure = framework_predictions[cat_name][0.0]
        prot_t1 = framework_predictions[cat_name][0.1]
        aib_pure = sorted([s for s in prot_pure if len(s) == 3 and s[1] == 'I'])
        aib_t1 = sorted([s for s in prot_t1 if len(s) == 3 and s[1] == 'I'])
        print(f"  {cat_name:>5s}: pure-Z protects A·I·B = {aib_pure}")
        print(f"        with T1=0.1·γ: {aib_t1}")
    print()

    # Now compare: for each backend run, show measured |⟨B_q0 A_q2⟩| for each (A,B)
    # alongside the framework's protected/active prediction.
    print(f"Hardware ⟨A_q0 B_q2⟩ values per backend, with framework verdict:")
    print()
    print(f"  Convention: framework label = AIB on sites (0,1,2);")
    print(f"  qiskit measurement of (b0=B, b2=A) corresponds to framework's AIB.")
    print()

    # Build expected protected sets in qiskit-(b0, b2) convention
    # framework's AIB → qiskit (b0=B, b2=A)
    cat_qiskit_protected = {'truly': {}, 'soft': {}, 'hard': {}}
    for cat in ['truly', 'soft', 'hard']:
        prot_t1 = framework_predictions[cat][0.1]  # use T1-extended (more realistic)
        prot_pure = framework_predictions[cat][0.0]
        for A in 'XYZ':
            for B in 'XYZ':
                label = f"{A}I{B}"
                # framework AIB → qiskit (b0=B, b2=A)
                key_qiskit = (B, A)
                cat_qiskit_protected[cat][key_qiskit] = (
                    label in prot_pure,  # protected under pure Z?
                    label in prot_t1,    # protected under Z+T1?
                )

    for backend_key in sorted(hw_results.keys()):
        d_data = hw_results[backend_key]
        if 'truly' not in d_data:
            continue
        print(f"  === {backend_key} ===")
        print(f"  {'(b0,b2)':>9s}  {'framework':>12s}  {'truly':>10s}   {'soft':>10s}   {'hard':>10s}")
        for A in 'XYZ':
            for B in 'XYZ':
                key_str = f"{A},{B}"
                truly_val = d_data.get('truly', {}).get(key_str, 'N/A')
                soft_val = d_data.get('soft', {}).get(key_str, 'N/A')
                hard_val = d_data.get('hard', {}).get(key_str, 'N/A')
                # Framework predicts protected for this (b0=A, b2=B)
                # which corresponds to framework's BIA
                pure_t, t1_t = cat_qiskit_protected['truly'].get((A, B), (None, None))
                truly_marker = "P→A" if pure_t and not t1_t else ("P" if pure_t else " ")
                pure_s, t1_s = cat_qiskit_protected['soft'].get((A, B), (None, None))
                soft_marker = "P→A" if pure_s and not t1_s else ("P" if pure_s else " ")
                pure_h, t1_h = cat_qiskit_protected['hard'].get((A, B), (None, None))
                hard_marker = "P→A" if pure_h and not t1_h else ("P" if pure_h else " ")

                fmt = lambda v: f"{v:>+8.3f}" if isinstance(v, (int, float)) else f"{v:>8s}"
                print(f"  {key_str:>9s}  {'P=protected':>12s}  "
                      f"{truly_marker:>3s} {fmt(truly_val):>8s}   "
                      f"{soft_marker:>3s} {fmt(soft_val):>8s}   "
                      f"{hard_marker:>3s} {fmt(hard_val):>8s}")
        print()

    print("Verdict columns:")
    print("  ' '   = framework predicts ACTIVE under both pure-Z and Z+T1")
    print("  'P'   = framework predicts PROTECTED under both pure-Z and Z+T1")
    print("  'P→A' = protected under pure-Z, but ACTIVE under Z+T1 (the 'fragile' case)")


if __name__ == "__main__":
    main()
