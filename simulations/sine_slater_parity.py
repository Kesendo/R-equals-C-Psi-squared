"""Parity pin: the frozen Slater-lift transcriptions vs the framework mirror of the C# canon.

Step 3b of the F89 de-monolith (the pattern of `weight_coherence_parity.py`): the rank-3 sine
Slater lift D_tau(z1,z2,z3) = det[sin(pi*k*(z+1)/(N+1))] is transcribed in three committed
engines of record; they stay byte-untouched and are pinned here against
`framework/sine_slater.py` (the Python mirror of the C# canonicals `XyJordanWignerModes` +
`JwSlaterPairBasis`) through each copy's documented adapter.

Adapter table (module raw = sine_mode_matrix(N, normalized=False), interior columns j = 0..N-1):

| copy | storage | call convention | normalization |
|------|---------|-----------------|---------------|
| A `cross_triple_orthogonality.umat/slater/slater_norm_sq` | wall-padded x = 0..N+1, walls zeroed | slater takes RAW columns (interior 1..N) | raw det; norm-SQUARED separate |
| B `y_zero_and_level_law.modes/slater/slater_norm`          | wall-padded z = -1..N (col = z+1)     | slater shifts +1 internally (interior 0..N-1) | raw det; norm (sqrt) separate |
| C `resonant_n_twinning` (closure in check_Y_is_zero)       | interior only j = 0..N-1              | raw columns  | vector PRE-normalized |

Copy C is a local closure inside `check_Y_is_zero` and is NOT importable; it is deliberately
not pinned here. Its construction line is the same interior raw sine as the module, its Slater
vector is normalized before use (the raw magnitude never enters), and its own downstream
asserts (kernel membership, isospectral twinning, Y = 0) are its guard. If that closure is ever
promoted to module level, add it to this gate.

Float honesty (unlike the integer-valued pencil gate): the sine entries are transcendental and
the three copies use different but mathematically equal multiplication orders into sin
(`k*x*pi/n` vs `pi*k*(j+1)/n`), so parity here is pinned to TIGHT TOLERANCES, not entry-exact:
matrix entries to 5e-15 absolute, 3x3 determinants to 1e-12, norm laws to 1e-9 relative.
The ((N+1)/2)^3 raw norm law is pinned on both copies (copy B also asserts it in-file).

Run: python simulations/sine_slater_parity.py   (~seconds; all gates must print PASS)
"""
from __future__ import annotations

import os
import sys
from itertools import combinations

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import framework as fw
import cross_triple_orthogonality as ct
import y_zero_and_level_law as yz

ENTRY_ATOL = 5e-15
DET_ATOL = 1e-12
NORM_RTOL = 1e-9


def gate_cross_triple(n_range=range(4, 9)):
    for n_sites in n_range:
        n = n_sites + 1
        u_ct = ct.umat(n_sites)
        u_m = fw.sine_mode_matrix(n_sites, normalized=False)
        assert u_ct.shape == (n_sites, n + 1)
        assert np.all(u_ct[:, 0] == 0.0) and np.all(u_ct[:, n] == 0.0), "walls not zeroed"
        assert np.abs(u_ct[:, 1:n] - u_m).max() <= ENTRY_ATOL, f"umat mismatch at N={n_sites}"
        for tau in ((1, 2, 3), (1, 3, n_sites)):   # distinct modes at every N >= 4
            for z in ((1, 2, 4), (2, 3, n_sites)):
                d_ct = ct.slater(u_ct, tau, *z)                       # raw columns 1..N
                d_m = fw.slater_det(u_m, tau, tuple(x - 1 for x in z))  # interior 0..N-1
                assert abs(d_ct - d_m) <= DET_ATOL, f"slater mismatch N={n_sites} tau={tau} z={z}"
            got = ct.slater_norm_sq(u_ct, tau, n_sites)
            want = fw.slater_norm_sq_law(n_sites, 3)
            assert abs(got - want) <= NORM_RTOL * want, f"norm-sq law N={n_sites} tau={tau}"
    print(f"PASS cross_triple_orthogonality.umat/slater/slater_norm_sq == sine_slater raw "
          f"(wall-padded adapter, raw-column calls), N={n_range.start}..{n_range.stop - 1}")


def gate_y_zero(n_range=range(4, 9)):
    for n_sites in n_range:
        u_yz = yz.modes(n_sites)
        u_m = fw.sine_mode_matrix(n_sites, normalized=False)
        assert u_yz.shape == (n_sites, n_sites + 2)
        assert np.abs(u_yz[:, 0]).max() <= 1e-13 and np.abs(u_yz[:, -1]).max() <= 1e-13, \
            "walls not (numerically) zero"      # sin(k*pi) is ~1e-16*k in float, not exact 0
        assert np.abs(u_yz[:, 1:-1] - u_m).max() <= ENTRY_ATOL, f"modes mismatch at N={n_sites}"
        for tau in ((1, 2, 3), (1, 3, n_sites)):   # distinct modes at every N >= 4
            for z in ((0, 1, 3), (1, 2, n_sites - 1)):
                d_yz = yz.slater(u_yz, tau, *z)                      # +1 shift internal
                d_m = fw.slater_det(u_m, tau, z)                     # interior direct
                assert abs(d_yz - d_m) <= DET_ATOL, f"slater mismatch N={n_sites} tau={tau} z={z}"
            got = yz.slater_norm(u_yz, tau, n_sites) ** 2
            want = fw.slater_norm_sq_law(n_sites, 3)
            assert abs(got - want) <= NORM_RTOL * want, f"norm law N={n_sites} tau={tau}"
    print(f"PASS y_zero_and_level_law.modes/slater/slater_norm == sine_slater raw "
          f"(+1-shift adapter), N={n_range.start}..{n_range.stop - 1}; "
          f"norm-squared law ((N+1)/2)^3 pinned on both sides")


def gate_full_vector(n_sites=6):
    # the whole rank-3 vector agrees across copy B's enumeration and the module's
    tau = (1, 3, 5)
    u_yz = yz.modes(n_sites)
    v_yz = np.array([yz.slater(u_yz, tau, *z) for z in combinations(range(n_sites), 3)])
    v_m = fw.slater_vector(n_sites, tau, normalized=False)
    assert np.abs(v_yz - v_m).max() <= DET_ATOL
    print(f"PASS full rank-3 Slater vector agrees (N={n_sites}, tau={tau}, "
          f"{len(v_m)} triples, combinations order)")


if __name__ == "__main__":
    gate_cross_triple()
    gate_y_zero()
    gate_full_vector()
    print("ALL SLATER PARITY GATES PASS: the frozen transcriptions match the one lift "
          "(tight-tolerance; transcendental entries, see the float-honesty note).")
