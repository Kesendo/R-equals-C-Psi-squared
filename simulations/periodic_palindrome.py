"""F1-style palindrome test on periodic-table properties.

The framework's F1 palindrome theorem says spec(L) is invariant under
λ → −λ − 2σ for some center σ. Atomic analogue: do per-element values
across a period form pairs whose sums sit at a constant?

Tested:
  (1) First ionization energies (eV): energy to remove one electron.
      Single-atom property, no coupling. Tests F1 directly on the
      element's electronic structure.

  (2) Pauling electronegativities: derived from bond polarity (atoms
      coupled in molecules). One layer up from atomic IE: EN exists
      only because atoms couple to other atoms. Same palindrome test.

  (3) Allen electronegativities (1989): configuration energy =
      weighted-average valence-shell ionization energy. Includes
      noble gases by definition. Strictest palindrome of the three.

For each property the test is the same:
  pair_sum(k) = value_k + value_{N−k+1}
  CoV of pair_sums = how palindromic the property is across the period.
  null = randomly permute the values and compute CoV. A small fraction
  of shuffles dropping below the actual CoV is the p-value.

Pauling EN doesn't include noble gases, so EN periods are 7 (period 2)
or 17 (period 4 with d-block) elements wide. Odd-length periods have
a center element which is excluded from pair sums.
"""
from __future__ import annotations

import numpy as np

# First ionization energies in eV. NIST standard values, rounded.
# Periods are even-length (groups 1-2 and 13-18 for periods 2-3;
# groups 1-2, 3-12, 13-18 for periods 4-5).

periods = {
    'Period 2 (Li-Ne)': {
        'Z': list(range(3, 11)),
        'symbols': ['Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne'],
        'IE':     [5.392, 9.323, 8.298, 11.260, 14.534, 13.618, 17.423, 21.565],
    },
    'Period 3 (Na-Ar)': {
        'Z': list(range(11, 19)),
        'symbols': ['Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar'],
        'IE':     [5.139, 7.646, 5.986, 8.152, 10.487, 10.360, 12.968, 15.760],
    },
    'Period 4 (K-Kr)': {
        'Z': list(range(19, 37)),
        'symbols': ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr'],
        'IE':     [4.341, 6.113, 6.561, 6.828, 6.746, 6.767, 7.434, 7.902, 7.881, 7.640,
                   7.726, 9.394, 5.999, 7.900, 9.815, 9.752, 11.814, 13.999],
    },
    'Period 5 (Rb-Xe)': {
        'Z': list(range(37, 55)),
        'symbols': ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                    'In', 'Sn', 'Sb', 'Te', 'I', 'Xe'],
        'IE':     [4.177, 5.695, 6.217, 6.634, 6.759, 7.092, 7.280, 7.361, 7.459, 8.337,
                   7.576, 8.994, 5.786, 7.344, 8.640, 9.010, 10.451, 12.130],
    },
    'Period 6 (Cs-Rn)': {
        # 32 elements including the lanthanide f-block (La through Lu).
        # The pairing test asks whether F1 palindrome structure survives the
        # f-block insertion, the most complex orbital-class addition in the
        # periodic table.
        'Z': list(range(55, 87)),
        'symbols': ['Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
                    'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os',
                    'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'],
        'IE':     [3.894, 5.212, 5.577, 5.539, 5.473, 5.525, 5.582, 5.644, 5.670, 6.150,
                   5.864, 5.939, 6.022, 6.108, 6.184, 6.254, 5.426, 6.825, 7.890, 7.980,
                   7.880, 8.700, 9.100, 9.000, 9.226, 10.438, 6.108, 7.417, 7.286, 8.417,
                   9.318, 10.748],
    },
}


# Pauling electronegativities. Noble gases excluded (not assigned in
# standard Pauling scale). Periods are odd-length; the center element
# is excluded from the pair-sum test.
en_periods = {
    'EN Period 2 (Li-F)': {
        'symbols': ['Li', 'Be', 'B',  'C',  'N',  'O',  'F'],
        'EN':     [0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98],
    },
    'EN Period 3 (Na-Cl)': {
        'symbols': ['Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl'],
        'EN':     [0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16],
    },
    'EN Period 4 (K-Br, full d-block)': {
        'symbols': ['K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
                    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br'],
        'EN':     [0.82, 1.00, 1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88,
                   1.91, 1.90, 1.65, 1.81, 2.01, 2.18, 2.55, 2.96],
    },
    'EN Period 5 (Rb-I, full d-block)': {
        'symbols': ['Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',
                    'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I'],
        'EN':     [0.82, 0.95, 1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28,
                   2.20, 1.93, 1.69, 1.78, 1.96, 2.05, 2.10, 2.66],
    },
}


# Allen electronegativity scale (Allen 1989: configuration energy = average
# valence-shell ionization energy). Includes noble gases by definition,
# so periods are even-length (matching IE structure). Direct apples-to-
# apples comparison with the IE_1 palindrome on the same period.
allen_en_periods = {
    'Allen Period 2 (Li-Ne)': {
        'symbols': ['Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'],
        'EN':     [0.912, 1.576, 2.051, 2.544, 3.066, 3.610, 4.193, 4.789],
    },
    'Allen Period 3 (Na-Ar)': {
        'symbols': ['Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar'],
        'EN':     [0.869, 1.293, 1.613, 1.916, 2.253, 2.589, 2.869, 3.242],
    },
}


def pair_sums(values):
    """Return list of pair sums v_k + v_{N−k+1} for k = 1, ..., N/2.
    For odd N, the center element is excluded."""
    n = len(values)
    return [values[k] + values[n - 1 - k] for k in range(n // 2)]


def palindrome_score(values):
    """Coefficient of variation of pair sums (lower = more palindromic)."""
    pairs = pair_sums(values)
    return float(np.std(pairs) / np.mean(pairs))


def shuffle_null(values, n_shuffles=10000, seed=0):
    """Distribution of palindrome scores under random permutation of values."""
    rng = np.random.RandomState(seed)
    scores = np.zeros(n_shuffles)
    arr = np.array(values)
    for i in range(n_shuffles):
        perm = rng.permutation(len(arr))
        scores[i] = palindrome_score(arr[perm])
    return scores


def report(name, symbols, values, value_label='value', value_fmt='{:6.3f}'):
    n = len(values)
    pairs = pair_sums(values)
    score = palindrome_score(values)
    null = shuffle_null(values)
    p_value = float(np.mean(null <= score))

    print(f'{name}')
    print(f'  N = {n}, pair sums:')
    for k, ps in enumerate(pairs):
        left, right = symbols[k], symbols[n - 1 - k]
        v_l, v_r = values[k], values[n - 1 - k]
        print(f'    ({left:<3s}, {right:<3s}): {value_fmt.format(v_l)} + {value_fmt.format(v_r)} = {ps:7.3f}')
    if n % 2 == 1:
        c_idx = n // 2
        print(f'    center: {symbols[c_idx]:<3s} = {value_fmt.format(values[c_idx])} (excluded)')
    mean_ps = np.mean(pairs)
    spread = np.max(pairs) - np.min(pairs)
    print(f'  pair-sum mean = {mean_ps:.3f}, spread = {spread:.3f} ({100*spread/mean_ps:.1f}%)')
    print(f'  CoV (palindrome score) = {score:.4f}')
    print(f'  null shuffle median CoV = {np.median(null):.4f}, p = {p_value:.4f}')
    print()


def main():
    print('=' * 78)
    print('Layer 1: F1-style palindrome on first ionization energies (eV)')
    print('=' * 78)
    print()
    for name, p in periods.items():
        report(name, p['symbols'], p['IE'])

    print('=' * 78)
    print('Layer 2: same test on Pauling electronegativities (coupling property)')
    print('=' * 78)
    print()
    for name, p in en_periods.items():
        report(name, p['symbols'], p['EN'])

    print('=' * 78)
    print('Layer 3: Allen EN scale (configuration-energy, includes noble gases)')
    print('=' * 78)
    print()
    for name, p in allen_en_periods.items():
        report(name, p['symbols'], p['EN'])


if __name__ == '__main__':
    main()
