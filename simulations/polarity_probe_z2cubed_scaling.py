"""Polarity probe #12: verify the Z2-cubed sub-cell pattern at N >= 3.

Probe 11 (N=2) found two classes of Pauli-string pairs:
  Class A (same Z_2^3 sub-cell):     BALANCE for ALL coefficient choices
  Class B (cross Z_2^3 sub-cell):    BALANCE only for phase-matched coefficients

Z_2^3 cell = (Klein_a, Klein_b, y_par) = ((#X+#Y) mod 2, (#Y+#Z) mod 2, #Y mod 2)

This probe extends the test to N=3 and N=4 to verify the pattern scales.

Sweeps:
  Sweep A: within-cell c (k_pauli = 2, 3, 4 Pauli strings, all in same cell)
           with random complex coefficients. Conjecture: BALANCE always.
  Sweep B: cross-cell c (k_pauli >= 2 Pauli strings spanning >= 2 cells)
           with random complex coefficients. Conjecture: BROKEN mostly.
  Sweep C: same as A but with random REAL coefficients (sanity check that
           ANY coefs work for same-cell).

If pattern scales: F112 candidate (StandardLindbladZ2CubedPiBalance) is real.

WHAT A "BALANCED" ROW HERE IS, AND WHY THE VERDICT IS THREE-WAY (2026-08-06).
The probe predates the F112 proof and was written to find the structural axis
empirically. It cannot MEASURE balance on a bit_b-homogeneous c against the
Π²-even Heisenberg H it uses, because for that PAIR of conditions the polarity
content it would divide by is zero as a theorem, both halves and not merely
their difference (proof Steps 3 and 4: the dissipator contributes no Π-conj ±i
content, so M_{+i} and M_{−i} come from L_H alone, which a Π²-even H empties).
Neither condition suffices by itself, and the cross-bit_b rows below are the
proof of that: same H, genuine ‖M_anti‖² from 6.9e-03 to 83.2. Measured over the
192 within-cell rows at N=3 and N=4 (Sweep A and Sweep C both build c inside one
cell), ‖M_anti‖² runs from exactly 0.0 (one row; the next three are 9.3e-68,
2.1e-36 and 7.5e-36) up to 1.5e-29, against a ‖M‖² of 2.6 to 2.2e+03. It is float noise from the vec-to-Pauli transform, and the asymmetry on
top of it is noise as well, so the ratio is noise over noise and reports nothing
in either direction. Those rows are marked DEGEN.

That is not a detail of presentation. Until 2026-08-06 the ratio divided by
max(‖M‖², 1e-15), which made every such row read BALANCED at ~1e-34 and hid the
emptiness; dividing by the polarity content instead, with no structural test,
flips three of the same rows to a spurious BROKEN at 1.8e-02, 4.3e-03 and
2.0e-03. Neither number is a measurement. Note that an exact `== 0.0` guard on
the denominator does not rescue this: it fires on exactly 1 of those 192 rows.
The verdict is therefore decided FIRST by `is_bit_b_homogeneous`, an exact
integer test on the Pauli letters, and only the rows that survive it are read
as BAL or BREAK.

READ THE VERDICT COUNTS WITH THAT IN MIND. Over N=2, 3, 4 the sweep now returns
132 BREAK, 321 DEGEN and zero BAL: the BAL branch is unreached, because every
bit_b-homogeneous row is taken by DEGEN first and no bit_b-MIXED row has ever
balanced here. So "the split partitions by bit_b parity" is now half true by
construction. What remains empirical, and is the whole of what this probe
establishes, is the other half: every one of the 132 cross-bit_b rows breaks,
none of them balances, and the row set is unchanged from before this rewrite.
"""

import sys
sys.path.insert(0, 'simulations')

import numpy as np
from itertools import product
from collections import defaultdict

from framework.diagnostics.polarity_coordinates import polarity_coordinates_from_L
from framework.pauli import _vec_to_pauli_basis_transform, total_bit_b_parity

PAULI_2X2 = {
    'I': np.eye(2, dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_string(letters):
    op = PAULI_2X2[letters[0]]
    for l in letters[1:]:
        op = np.kron(op, PAULI_2X2[l])
    return op


def z2cubed_cell(p):
    """(Klein_a, Klein_b, y_par) for Pauli string p."""
    nx = p.count('X')
    ny = p.count('Y')
    nz = p.count('Z')
    return ((nx + ny) % 2, (ny + nz) % 2, ny % 2)


def L_vec_to_pauli(L_vec, N):
    T = _vec_to_pauli_basis_transform(N)
    return (T.conj().T @ L_vec @ T) / (2 ** N)


def heisenberg_terms(N):
    """The Pauli strings H is built from, as letter lists. H is the sum of these.

    Kept separate from `heisenberg_H` so the Π²-parity of the ACTUAL summands can
    be checked, rather than a literal restating what this function is believed to
    build. `assert_pi2_even_H` reads this list.
    """
    terms = []
    for b in range(N - 1):
        for letter in ['X', 'Y', 'Z']:
            letters = ['I'] * N
            letters[b] = letter
            letters[b + 1] = letter
            terms.append(letters)
    return terms


def heisenberg_H(N):
    """H = sum_b (XX + YY + ZZ) on bonds."""
    d = 2 ** N
    H = np.zeros((d, d), dtype=complex)
    for letters in heisenberg_terms(N):
        H = H + pauli_string(letters)
    return H


def bit_b_parity(p):
    """Π²-parity of a Pauli string. Thin alias for the framework primitive.

    `framework.pauli.total_bit_b_parity` is the canonical definition (Σ bit_b
    mod 2, bit_b of a letter being n_YZ); for a string that is (#Y + #Z) mod 2,
    the SECOND component of the Z_2^3 cell above. Exact: no float enters it.
    Aliased rather than re-derived so the two cannot drift apart.
    """
    return total_bit_b_parity(list(p))


def is_bit_b_homogeneous(selected):
    """True iff every Pauli string in c shares one Π²-parity.

    This is F112's scope test on c, decided by integer counts alone, and it is
    the EXACT companion the numeric relative asymmetry needs. It is HALF of the
    condition for a vanishing polarity content, not all of it: by Steps 3 and 4
    of PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md, a bit_b-homogeneous c puts no
    Π-conjugation ±i content into M, so M_{+i} and M_{−i} come entirely from the
    Hamiltonian part L_H, and they vanish too only when H is itself Π²-even.

    Both halves vanish there, not merely their difference: that is why the rows
    carry no reading rather than a balanced one. F112 as stated gives the
    difference; Step 3's "zero Π-conj +i and zero Π-conj −i content" gives the sum.

    The caller must therefore supply a Π²-even H. `main` uses the Heisenberg
    XX+YY+ZZ, whose every bilinear has bit_b parity 0, and `assert_pi2_even_H`
    below pins that rather than assuming it.
    """
    return len({bit_b_parity(p) for p in selected}) == 1


def is_pi2_even_H(h_terms):
    """True iff every summand of H carries even Π²-parity. Exact integer test.

    Read by `test_c` on every row, so the DEGEN verdict depends on the ACTUAL H
    in play rather than on a precondition asserted once in `main`.
    """
    return all(total_bit_b_parity(t) == 0 for t in h_terms)


def assert_pi2_even_H(N):
    """Pin the second half of the degeneracy condition: H is Π²-even.

    Exact integer check over the ACTUAL summands of H, read from
    `heisenberg_terms`. If any carried odd Π²-parity, L_H would contribute ±i
    content of its own, the polarity content would NOT vanish, and every DEGEN
    verdict below would be unsound.

    Written this way on purpose. A first version tested `bit_b_parity([l, l])`
    over the letters X, Y, Z, which is 0 for EVERY letter because a doubled
    letter has an even #Y + #Z: a guard that could not fail, and one whose
    mutation test mutated the parity helper rather than the H it was guarding.
    The object an assertion protects has to be the object it reads.
    """
    odd = [''.join(t) for t in heisenberg_terms(N) if total_bit_b_parity(t) != 0]
    assert not odd, f"H carries Π²-odd terms {odd}; the DEGEN verdict would be unsound"


def test_c(N, c, H, sigma, selected, h_terms):
    """Return (‖M‖², ‖M_anti‖², rel_asym, degenerate) for the jump operator c.

    `rel_asym` divides by the POLARITY CONTENT ‖M_anti‖² = ‖M_+‖² + ‖M_−‖², the
    whole of what the asymmetry is a difference of, not by ‖M‖², which is a
    different scale (see simulations/framework/workflows/polarity_fingerprint.py).
    `degenerate` is the exact structural verdict from `selected` and outranks the
    ratio: see the note in this module's docstring.
    """
    d = 2 ** N
    Id = np.eye(d, dtype=complex)
    L_vec = -1j * (np.kron(H, Id) - np.kron(Id, H.T)) + 0.1 * np.kron(c, c.conj())
    L_pauli = L_vec_to_pauli(L_vec, N)
    result = polarity_coordinates_from_L(L_pauli, N, sigma)
    ns_M = float(result['norm_sq']['M'])
    ns_anti = float(result['norm_sq']['M_plus_half'] + result['norm_sq']['M_minus_half'])
    rel_asym = 0.0 if ns_anti == 0.0 else abs(result['asymmetry']) / ns_anti
    # BOTH halves of the condition, decided here rather than in the caller: a c-only
    # test would call a Pi^2-ODD H's substantive reading DEGEN (N=2, H = 0.7*ZI with
    # c = ZZ gives a genuine ||M_anti||^2 = 15.68 and an exactly-zero asymmetry).
    degenerate = is_bit_b_homogeneous(selected) and is_pi2_even_H(h_terms)
    return ns_M, ns_anti, rel_asym, degenerate


def verdict(rel_asym, degenerate):
    """DEGEN outranks the ratio; BAL/BREAK only where there is content to read."""
    if degenerate:
        return "DEGEN"
    return "BAL" if rel_asym < 1e-10 else "BREAK"


def main():
    for N in [2, 3, 4]:
        d = 2 ** N
        sigma = 0.1
        h_terms = heisenberg_terms(N)
        H = heisenberg_H(N)
        assert_pi2_even_H(N)

        # Group all 4^N Pauli strings by Z_2^3 cell
        all_paulis = list(product(['I', 'X', 'Y', 'Z'], repeat=N))
        cells = defaultdict(list)
        for p in all_paulis:
            cells[z2cubed_cell(p)].append(p)

        print("=" * 100)
        print(f"N = {N}, d = {d}, total Pauli strings = {len(all_paulis)}, Z_2^3 cells = {len(cells)}")
        print("=" * 100)
        for cell, members in sorted(cells.items()):
            sample = members[:3] + (['...'] if len(members) > 3 else [])
            print(f"  cell {cell}: {len(members):>3} members, e.g. "
                  f"{['' .join(p) if isinstance(p, tuple) else p for p in sample]}")
        print()

        rng = np.random.default_rng(seed=2026 + N * 1000)

        # ---- Sweep A: within-cell c, k_pauli in {2, 3, 4} ----
        print(f"--- Sweep A (within-cell c, random complex coefs, 3 trials) ---")
        print(f"{'cell':<12}  {'k':<3}  {'trial':<5}  {'||M||^2':>12}  "
              f"{'||M_anti||^2':>14}  {'rel_asym':>12}  status")
        print("-" * 75)
        for cell, members in sorted(cells.items()):
            if len(members) < 2:
                continue
            for k_pauli in [2, 3, 4]:
                if k_pauli > len(members):
                    continue
                for trial in range(3):
                    indices = rng.choice(len(members), size=k_pauli, replace=False)
                    selected = [members[i] for i in indices]
                    coeffs = rng.normal(size=k_pauli) + 1j * rng.normal(size=k_pauli)
                    c = sum(coef * pauli_string(p) for coef, p in zip(coeffs, selected))
                    ns_M, ns_anti, rel_asym, degen = test_c(N, c, H, sigma, selected, h_terms)
                    marker = verdict(rel_asym, degen)
                    print(f"{str(cell):<12}  {k_pauli:<3}  {trial:<5}  {ns_M:>12.4f}  "
                          f"{ns_anti:>14.3e}  {rel_asym:>12.4e}  [{marker}]")
        print()

        # ---- Sweep B: cross-cell c (pick 1 from each of 2 cells), random coefs ----
        print(f"--- Sweep B (cross-cell c, 1 Pauli from each of 2 cells, random complex coefs, 5 trials) ---")
        print(f"{'cell_a':<10}  {'cell_b':<10}  {'p_a':<8}  {'p_b':<8}  {'trial':<5}  "
              f"{'||M||^2':>12}  {'||M_anti||^2':>14}  {'rel_asym':>12}  status")
        print("-" * 90)
        cell_list = list(cells.keys())
        for i, cell_a in enumerate(cell_list):
            for cell_b in cell_list[i+1:]:
                if not cells[cell_a] or not cells[cell_b]:
                    continue
                # Pick one representative each
                p_a = cells[cell_a][0]
                p_b = cells[cell_b][0]
                for trial in range(3):
                    a, b = (complex(rng.normal(), rng.normal()),
                            complex(rng.normal(), rng.normal()))
                    c = a * pauli_string(p_a) + b * pauli_string(p_b)
                    ns_M, ns_anti, rel_asym, degen = test_c(N, c, H, sigma, [p_a, p_b], h_terms)
                    marker = verdict(rel_asym, degen)
                    print(f"{str(cell_a):<10}  {str(cell_b):<10}  "
                          f"{''.join(p_a):<8}  {''.join(p_b):<8}  {trial:<5}  "
                          f"{ns_M:>12.4f}  {ns_anti:>14.3e}  {rel_asym:>12.4e}  [{marker}]")
        print()

        # ---- Sweep C: within-cell with random REAL coefs (extra sanity) ----
        print(f"--- Sweep C (within-cell c, random REAL coefs, 3 trials per cell, k=4) ---")
        print(f"{'cell':<12}  {'trial':<5}  {'||M||^2':>12}  "
              f"{'||M_anti||^2':>14}  {'rel_asym':>12}  status")
        print("-" * 60)
        for cell, members in sorted(cells.items()):
            if len(members) < 4:
                continue
            for trial in range(3):
                indices = rng.choice(len(members), size=4, replace=False)
                selected = [members[i] for i in indices]
                coeffs = rng.normal(size=4)  # REAL coefficients
                c = sum(coef * pauli_string(p) for coef, p in zip(coeffs, selected))
                ns_M, ns_anti, rel_asym, degen = test_c(N, c, H, sigma, selected, h_terms)
                marker = verdict(rel_asym, degen)
                print(f"{str(cell):<12}  {trial:<5}  {ns_M:>12.4f}  "
                      f"{ns_anti:>14.3e}  {rel_asym:>12.4e}  [{marker}]")
        print()
        print()


if __name__ == '__main__':
    main()
