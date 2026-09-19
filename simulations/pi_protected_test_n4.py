#!/usr/bin/env python3
"""EQ-026 finite audit of the separate historical 120-member catalog.

The catalog is the combinations-with-replacement set from the 15 non-II
two-site Pauli words, including all 15 self-pairs.  It is not the N=3 census
of 36 distinct pairs (14 hard / 19 soft / 3 truly).  For a selected finite N
(default 4), each historical catalog Hamiltonian is diagnosed by:
  1. Operator residual ‖M‖ via framework.palindrome_residual
  2. Eigenvalue spectrum pairing under λ ↔ −λ − 2Σγ
  3. Number of Π-protected Pauli-string observables on |+−+−...⟩

The finite historical comparison is 15 truly / 46 soft / 59 hard.  These
operator and spectral diagnostics do not establish a causal mechanism or an
all-N universality claim.

The init state |+−+−...⟩ generalises Snapshot D's N=3 |+−+⟩.

Usage:
    python pi_protected_test_n4.py            # N=4 (default)
    python pi_protected_test_n4.py 5          # N=5; runtime depends on the host
"""
import math
import sys
from itertools import combinations_with_replacement

import numpy as np

try:
    from simulations import framework as fw
    from simulations.framework.diagnostics.f77_trichotomy import spectrum_pairing_error
except ModuleNotFoundError as exc:  # Direct execution from the simulations directory.
    if exc.name != "simulations":
        raise
    import framework as fw
    from framework.diagnostics.f77_trichotomy import spectrum_pairing_error


# Historical EQ-026 catalog: all 15 non-II two-site Pauli words and unordered
# pairs with replacement, C(15+1, 2) = 120, including 15 self-pairs.
TWO_SITE_PAULI_WORDS = (
    "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI",
    "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ",
)
_BPE = {"II", "XX", "YY", "ZZ"}


def historical_two_site_pair_class(term1, term2):
    """Return the historical 15/46/59 finite-catalog class.

    This is not the nine-fully-lit F87 router: identity letters are allowed,
    while the identity-only word ``II`` is excluded from the alphabet.
    """
    n1, n2 = term1.count("I"), term2.count("I")
    if n1 == 1 and n2 == 1:
        letters = {c for term in (term1, term2) for c in term if c != "I"}
        if "Z" in letters:
            return "hard"
        return "truly" if letters == {"X"} else "soft"
    if (n1 == 1) != (n2 == 1):
        single, double = (term1, term2) if n1 == 1 else (term2, term1)
        letter = next(c for c in single if c != "I")
        if letter == "Z":
            return "hard"
        if letter == "X":
            if double in _BPE:
                return "truly"
            return "soft" if double in {"YZ", "ZY"} else "hard"
        return "soft" if double in _BPE or double in {"XZ", "ZX"} else "hard"
    if term1 in _BPE and term2 in _BPE:
        return "truly"
    if term1 == term2:
        return "soft"
    if (term1 in _BPE) != (term2 in _BPE):
        bpe, other = (term1, term2) if term1 in _BPE else (term2, term1)
        axis = bpe[0]
        if axis == "Z":
            return "soft"
        partner = "Y" if axis == "X" else "X"
        if ((other[0] == axis and other[1] == partner)
                or (other[1] == axis and other[0] == partner)):
            return "hard"
        return "soft"
    if term1[::-1] == term2:
        return "soft"
    if ((term1[0] == "Z" and term2[0] == "Z")
            or (term1[1] == "Z" and term2[1] == "Z")):
        return "soft"
    return "hard"


def historical_two_site_pair_catalog():
    """Return the exact lexicographic 120-member historical catalog."""
    return tuple(
        (left, right, historical_two_site_pair_class(left, right))
        for left, right in combinations_with_replacement(TWO_SITE_PAULI_WORDS, 2)
    )


def format_historical_catalog_comparison(N, results):
    """Format the selected-N row from the verdicts that were actually computed."""
    n_truly = sum(1 for row in results if row["verdict"] == "truly")
    n_soft = sum(1 for row in results if row["verdict"] == "soft")
    n_hard = sum(1 for row in results if row["verdict"] == "hard")
    return "\n".join((
        f"Historical 120-member catalog at selected N={N}:",
        "  Distinct-pair N=3 control (different 36-member sample): "
        "3 truly, 19 soft, 14 hard",
        f"  N={N}: {n_truly} truly, {n_soft} soft, {n_hard} hard",
    ))


def spectrum_pair_max_err(L, sigma_gamma):
    """Exact one-use bottleneck distance from the spectrum to its mirror."""
    evals = np.linalg.eigvals(L)
    return spectrum_pairing_error(evals, sigma_gamma)


def compute_historical_catalog_results(N, gamma=0.1, coupling=1.0):
    """Compute the selected-N verdict rows consumed by the CLI summary."""
    sigma_gamma = N * gamma
    # |+−+−...⟩ initial state, alternating starting with +
    plus = np.array([1, 1], dtype=complex) / math.sqrt(2)
    minus = np.array([1, -1], dtype=complex) / math.sqrt(2)
    psi = plus.copy()
    for k in range(1, N):
        psi = np.kron(psi, minus if k % 2 == 1 else plus)
    rho_0 = np.outer(psi, psi.conj())

    bonds = [(i, i + 1) for i in range(N - 1)]
    # All 120 unordered pairs with replacement from the 15 non-II words.
    pairs = [
        ((tuple(left), tuple(right)), (tuple(left), tuple(right)))
        for left, right, _ in historical_two_site_pair_catalog()
    ]

    results = []
    for sorted_terms, terms in pairs:
        bilinear = [(t[0], t[1], coupling) for t in terms]
        H = fw._build_bilinear(N, bonds, bilinear)
        L = fw.lindbladian_z_dephasing(H, [gamma] * N)

        # Operator residual
        M = fw.palindrome_residual(L, sigma_gamma, N)
        op_norm = float(np.linalg.norm(M))

        # Spectrum pairing
        spec_err = spectrum_pair_max_err(L, sigma_gamma)

        # Π-protected observables
        result = fw.pi_protected_observables(H, [gamma] * N, rho_0, N)
        n_prot = len(result['protected'])
        n_act = len(result['active'])

        # Classify
        spec_ok = spec_err < 1e-6
        op_ok = op_norm < 1e-10
        if op_ok:
            verdict = "truly"
        elif spec_ok:
            verdict = "soft"
        else:
            verdict = "hard"

        label = f"{terms[0][0]}{terms[0][1]}+{terms[1][0]}{terms[1][1]}"
        results.append({
            'pair': sorted_terms, 'label': label,
            'op_norm': op_norm, 'spec_err': spec_err,
            'n_protected': n_prot, 'n_active': n_act,
            'verdict': verdict,
        })

    order = {'truly': 0, 'soft': 1, 'hard': 2}
    results.sort(key=lambda r: (order[r['verdict']], -r['n_protected']))
    return results


def main(selected_n=None):
    if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    N = int(sys.argv[1]) if selected_n is None and len(sys.argv) > 1 else (
        4 if selected_n is None else int(selected_n)
    )
    GAMMA = 0.1
    J = 1.0
    pair_count = len(historical_two_site_pair_catalog())

    init_label = "|" + "".join("+" if k % 2 == 0 else "−" for k in range(N)) + "⟩"
    print(f"Testing pi_protected_observables at N={N}, γ={GAMMA}, init {init_label}")
    print(f"{pair_count} unordered two-term Pauli-pair Hamiltonians "
          "(15 non-II words, combinations with replacement).")
    print(f"Liouvillian dim = 4^{N} = {4**N}; expected runtime ~ {pair_count * (4**N)**3 / 1e10:.1f} min")
    print()
    print(f"  {'H = J(...)':>12s}  {'‖M‖_op':>10s}  {'spec_err':>10s}  "
          f"{'protected':>10s}  {'active':>8s}  {'verdict':>10s}")
    print(f"  {'-' * 12}  {'-' * 10}  {'-' * 10}  {'-' * 10}  {'-' * 8}  {'-' * 10}")

    results = compute_historical_catalog_results(
        N, gamma=GAMMA, coupling=J
    )
    for r in results:
        print(f"  {r['label']:>12s}  {r['op_norm']:>10.2e}  {r['spec_err']:>10.2e}  "
              f"{r['n_protected']:>10d}  {r['n_active']:>8d}  {r['verdict']:>10s}")

    # Summary
    n_truly = sum(1 for r in results if r['verdict'] == 'truly')
    n_soft = sum(1 for r in results if r['verdict'] == 'soft')
    n_hard = sum(1 for r in results if r['verdict'] == 'hard')

    print()
    print(f"Summary at N={N}:")
    print(f"  truly:  {n_truly:>3d}  (operator equation exact, spectrum paired)")
    print(f"  soft:   {n_soft:>3d}  (operator equation broken, spectrum still paired)")
    print(f"  hard:   {n_hard:>3d}  (both broken)")
    print(f"  total:  {n_truly + n_soft + n_hard:>3d}")
    print()
    print(f"Protected-observable count distribution:")
    n_prot_truly = [r['n_protected'] for r in results if r['verdict'] == 'truly']
    n_prot_soft = [r['n_protected'] for r in results if r['verdict'] == 'soft']
    n_prot_hard = [r['n_protected'] for r in results if r['verdict'] == 'hard']
    if n_prot_truly:
        print(f"  truly:  {min(n_prot_truly):>3d} - {max(n_prot_truly):>3d}  (range)")
    if n_prot_soft:
        print(f"  soft:   {min(n_prot_soft):>3d} - {max(n_prot_soft):>3d}  (range)")
    if n_prot_hard:
        print(f"  hard:   {min(n_prot_hard):>3d} - {max(n_prot_hard):>3d}  (range)")

    print()
    print(format_historical_catalog_comparison(N, results))
    return results


if __name__ == "__main__":
    main()
