"""Welle 6: polarity_fingerprint applied to Tier-B Marrakesh/Kingston datasets.

Tier-B datasets are single-time-snapshot 2-qubit tomography (16 Pauli
expectations per Hamiltonian, no t-series for L_eff fitting). The
polarity_fingerprint workflow operates on (chain, terms) inputs, not on
trajectories, so it applies to the Hamiltonians used in each dataset
without needing the snapshot data itself for the framework reading.

This script exercises the typed `fw.polarity_fingerprint` workflow on
the 11 Hamiltonians used across:
  - data/ibm_soft_break_april2026/...marrakesh...json     (3 Hamiltonians)
  - data/ibm_f83_signature_april2026/...                  (4 Hamiltonians)
  - data/ibm_soft_break_april2026/...kingston...json      (4 Hamiltonians)

(zn_mirror is excluded: state_a / state_b are state-prep variants, not
Hamiltonian-classification categories, so polarity_fingerprint doesn't
directly apply.)

For each Hamiltonian:
  1. Run polarity_fingerprint(chain, terms) → F87 class + F112 verdict
  2. Report joint reading + cross-check that framework F87 matches the
     dataset's category label (which is itself an F87 classification name)

What this was expected to find, and what it actually finds (corrected 2026-08-07):
the expectation was "F112 asymmetry = 0 BALANCED across all 11", on the reading
that chain.L is standard bit_b-homogeneous Z-dephasing whatever the F87 class.
Six of the eleven come out that way. The other five have a Pi^2-EVEN H, which
against that same homogeneous c leaves no polarity content to balance: they are
SILENT, and a BALANCED read off them would have confirmed nothing. The
orthogonality anchor this extends (polarity_probe_f87_connection.py, 7 synthetic
cases of which 4 are substantive) therefore reaches 6 substantive
hardware-Hamiltonian instances, all tested on Marrakesh/Kingston in April-May 2026.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import framework as fw  # noqa: E402

# Term mappings sourced from simulations/framework/tests/workflows/test_diagnose_hardware.py
# (F83_TERMS_PER_CATEGORY) and simulations/f80_ibm_soft_break_check.py (soft_break_marrakesh).

F83_TERMS = {
    'truly_unbroken':       [('X', 'X'), ('Y', 'Y')],
    'pi2_odd_pure':         [('X', 'Y'), ('Y', 'X')],
    'pi2_even_nontruly':    [('Y', 'Z'), ('Z', 'Y')],
    'mixed_anti_one_sixth': [('X', 'Y'), ('Y', 'Z')],
}

SOFT_BREAK_MARRAKESH_TERMS = {
    'truly_unbroken': [('X', 'X'), ('Y', 'Y')],
    'soft_broken':    [('X', 'Y'), ('Y', 'X')],
    'hard_broken':    [('X', 'X'), ('X', 'Y')],
}


def fingerprint_table(dataset_name, terms_per_category, chain):
    print(f"\n{'=' * 100}")
    print(f"  {dataset_name}")
    print(f"{'=' * 100}\n")
    print(f"  {'Category':<28} {'Terms':<26} {'F87':<8} {'F112':<14} {'rel asym':>12} {'In typed scope':<14}")
    print(f"  {'-' * 110}")

    rows = []
    for cat, terms in terms_per_category.items():
        result = fw.polarity_fingerprint(chain, terms)
        rows.append((cat, terms, result))
        terms_str = ' '.join(''.join(p) for p in terms)
        print(f"  {cat:<28} {terms_str:<26} {result['f87_class']:<8} "
              f"{result['f112_verdict']:<14} {result['f112_rel_asymmetry']:>12.4e} "
              f"{str(result['in_f112_typed_scope']):<14}")

    return rows


def main():
    print("Welle 6: polarity_fingerprint on Tier-B Marrakesh/Kingston Hamiltonians")
    print("=" * 100)
    print()
    print("Apply fw.polarity_fingerprint to the 11 Hamiltonians measured across:")
    print("  - soft_break Marrakesh (2026-04-26): 3 Hamiltonians")
    print("  - f83_signature Marrakesh (2026-04-30): 4 Hamiltonians")
    print("  - soft_break Kingston (2026-05-05): 4 Hamiltonians (same set as f83)")
    print()
    print("All datasets: N=3 chain, J=1.0, t=0.8μs, n_trotter=3, shots=4096.")
    print("Framework chain: ChainSystem(N=3, gamma_0=0.05) with default Z-dephasing.")
    print()

    chain = fw.ChainSystem(N=3, gamma_0=0.05)

    all_rows = []
    all_rows.extend(('soft_break Marrakesh', cat, terms, r) for cat, terms, r in
                    fingerprint_table('soft_break Marrakesh (2026-04-26, path [48, 49, 50])',
                                      SOFT_BREAK_MARRAKESH_TERMS, chain))
    all_rows.extend(('f83_signature Marrakesh', cat, terms, r) for cat, terms, r in
                    fingerprint_table('f83_signature Marrakesh (2026-04-30, path [4, 5, 6])',
                                      F83_TERMS, chain))
    all_rows.extend(('soft_break Kingston', cat, terms, r) for cat, terms, r in
                    fingerprint_table('soft_break Kingston (2026-05-05, path [43, 56, 63])',
                                      F83_TERMS, chain))

    print(f"\n{'=' * 100}")
    print("  Aggregate orthogonality reading")
    print(f"{'=' * 100}\n")

    n_total = len(all_rows)
    n_balanced = sum(1 for _, _, _, r in all_rows if r['f112_verdict'] == 'BALANCED')
    n_in_scope = sum(1 for _, _, _, r in all_rows if r['in_f112_typed_scope'])

    # Split by silence BEFORE aggregating. A silent row's relative asymmetry is NaN (0/0,
    # there is no polarity content to form a ratio of), and NaN loses every comparison, so
    # feeding it to max() hides it whenever a finite value happens to come first: with the
    # soft class ordered as it is, two silent rows would vanish behind a 0.0 and the class
    # would look fully measured. Count them instead.
    by_f87 = {}
    silent_by_f87 = {}
    for _, _, _, r in all_rows:
        by_f87.setdefault(r['f87_class'], [])
        silent_by_f87.setdefault(r['f87_class'], 0)
        if r['f112_verdict'] == 'DEGENERATE':
            silent_by_f87[r['f87_class']] += 1
        else:
            by_f87[r['f87_class']].append(r['f112_rel_asymmetry'])

    n_silent = sum(1 for _, _, _, r in all_rows if r['f112_verdict'] == 'DEGENERATE')
    print(f"  Total Hamiltonian-instances: {n_total}")
    print(f"  F112 BALANCED with genuine polarity content: {n_balanced} / {n_total}")
    print(f"  Structurally SILENT (no content to balance): {n_silent} / {n_total}")
    print(f"  In F112 typed scope (Hermitian H + bit_b-homog c): {n_in_scope} / {n_total}")
    print()
    print("  F112 reading per F87 class (max relative asymmetry over the MEASURABLE rows):")
    for f87_class, rels in sorted(by_f87.items()):
        n_silent_here = silent_by_f87[f87_class]
        if rels:
            print(f"    F87 {f87_class:<12}: {len(rels)} measurable, max rel asym = "
                  f"{max(rels):.4e}, {n_silent_here} silent")
        else:
            print(f"    F87 {f87_class:<12}: 0 measurable, {n_silent_here} silent, so this "
                  "class carries NO F112 evidence here")

    print()
    print(f"{'=' * 100}")
    print("  Conclusion")
    print(f"{'=' * 100}\n")
    print(f"  {n_balanced} of the {n_total} hardware-Hamiltonian instances are F112 BALANCED")
    print("  with a genuine polarity content to be balanced, under the framework's standard")
    print("  chain.L (Heisenberg-style H + single-Pauli Z-dephasing, which is bit_b-homogeneous")
    print("  on the c-side and Hermitian on the H-side, so F112 typed Tier1Derived predicts")
    print("  asymmetry = 0).")
    print()
    print(f"  The other {n_silent} are SILENT, not balanced, and this script counted them as")
    print("  confirmations until 2026-08-07. Their H is Pi^2-EVEN, so together with the")
    print("  homogeneous c the polarity content vanishes as a theorem: the asymmetry is 0 - 0")
    print("  and no verdict may be read off them. They are not evidence against F112 either;")
    print("  they are simply not evidence.")
    print()
    print("  This extends the F87↔F112 orthogonality empirical anchor from synthetic")
    print("  (simulations/polarity_probe_f87_connection.py: 7 cases, 3 of them silent and 4")
    print("  substantive, per that script's own printed conclusion)")
    print("  to 11 real-hardware-tested Hamiltonian instances across 3 datasets / 2 IBM")
    print("  backends (Marrakesh twice, Kingston once). The extension is narrower than that")
    print("  sentence reads, in two ways that have to be said in the same breath.")
    print()
    print("  First, the anchor rests on the measurable rows alone, and there are 6 of them.")
    print("  F87 classification varies across them while the F112 verdict stays BALANCED,")
    print("  which is the orthogonality. But the F87 'truly' class contributes NO measurable")
    print("  row. The Pi^2 parity of a term is (#Y + #Z) mod 2, and a truly term has #Y even")
    print("  AND #Z even, so its sum is even too: every truly H is Pi^2-EVEN, and every")
    print("  Pi^2-even H is silent against this homogeneous c. Not a biconditional, though:")
    print("  silence is the WIDER condition. YZ+ZY is soft and silent too, and over the 15")
    print("  two-letter bond strings other than II, of their 105 unordered pairs 11 are silent")
    print("  without being truly while 0 are truly without being silent. Only the silence half")
    print("  of that is a letter count and N-independent; `truly` here comes from")
    print("  classify_pauli_pair, a thresholded spectral verdict, so the 11 is an observation")
    print("  on this classifier (checked at N = 2, 3, 4 at unit coupling), not a derivation.")
    print("  The independence of the two axes is therefore anchored on 'soft' and 'hard' only,")
    print("  and a truly-class instance needs a c that is NOT bit_b-homogeneous before it can")
    print("  say anything at all.")
    print()
    print("  Second, the 11 rows are not 11 Hamiltonians. Every fingerprint here is taken at")
    print("  default rates and no declared drive, so within THIS run the reading is a function")
    print("  of (chain, terms) alone; the workflow itself also takes gamma_z, gamma_T1,")
    print("  gamma_pump, a drive profile and a tolerance. The three datasets reuse the same term")
    print("  sets, so the 6 measurable rows are 3 DISTINCT Hamiltonians: XY+YX counted three")
    print("  times, XX+XY once, XY+YZ twice. The soft class in particular is one Hamiltonian")
    print("  counted three times. The dataset multiplicity is real hardware provenance but")
    print("  carries no extra information for THIS reading, which never touches the")
    print("  measurements (see the 'does NOT do' section of the experiment writeup).")
    print()
    print("  Hardware-effective L (with non-standard noise channels measured on these")
    print("  backends) would require trajectory data to fit, which Tier-B snapshot datasets")
    print("  don't provide. Future hardware proposal: time-sweep variants of these same")
    print("  Hamiltonians would let us compute F112 asymmetry directly on the hardware-")
    print("  effective L (per Welle 2's Tier-A pattern), revealing any backend-specific")
    print("  polarity-axis structure beyond the framework prediction.")


if __name__ == '__main__':
    main()
