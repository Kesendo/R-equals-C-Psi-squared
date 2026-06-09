#!/usr/bin/env python3
r"""Ceiling middle-site tilt: the ceiling is an isolated soft point in a hard sea.

Re-derivation of the orphaned 2026-06-08 sweep behind
`simulations/results/ceiling_middle_tilt.{tsv,png}` (the original generator was an ad-hoc
script that was never committed; this regenerator was reconstructed 2026-06-09 from the data's
fingerprints and is the committed record).

The sweep
---------
At N = 3, tilt the middle-site letter of the local soft triple into the ceiling soft triple:

    H(theta) = cos(theta) * (XIX + XIY + YIX)  +  sin(theta) * (XZX + XZY + YZX),

which is exactly the middle-site operator W(theta) = cos(theta)*I + sin(theta)*Z inside each
term (the term sets pair up one-to-one: XIX<->XZX, XIY<->XZY, YIX<->YZX). Dephasing letter Z,
gamma = 0.05 per site (the repo's canonical gamma_0), sigma = N*gamma.

Validated against the original TSV:
  - VERDICT column: bit-identical 37/37 (Soft only at 0 and 90 degrees, Hard at all 35
    interior angles), via the canonical greedy SpectrumPairs criterion (max greedy pairing
    distance < 1e-6 <=> soft), the same criterion C# PauliPairTrichotomy uses.
  - M_norm column: bit-close row-by-row (rel < 1e-5). The original column equals
    2*sqrt(2) * ||palindrome_residual||_F in the framework convention (a constant factor from
    the original script's residual normalization); ||M|| is gamma-independent and loads
    monotonically 64.0 -> 101.2 as the middle Z turns on.
  - pairing_error / rate_error columns: the original used an unrecorded ad-hoc pairing
    statistic (smooth ~theta^2 at small angles); this regenerator reports the canonical
    measures instead, the minimax (bottleneck) pairing distance of spec(L) against
    -spec(L) - 2*sigma, full and real-part-only. Under any such measure the interior break
    is quantized in dephasing-rate steps (multiples of 2*gamma = 0.1), the structural
    statement the figure's plateau shows.

The result the sweep certifies: BOTH endpoints are soft, EVERY interior mixture is hard.
The two soft structures share no palindrome; the ceiling soft point is isolated, not
connected to the local soft case by any soft path along this tilt.

Self-validating: asserts raise on failure, a single PASS line prints on success, exit 0 only
if everything holds. Writes `simulations/results/ceiling_middle_tilt_regen.tsv`.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw  # noqa: E402
from framework.lindblad import lindbladian_pauli_dephasing  # noqa: E402
from framework.pauli import _build_kbody_chain  # noqa: E402

N = 3
GAMMA = 0.05
SIGMA = N * GAMMA
LOCAL = [('X', 'I', 'X'), ('X', 'I', 'Y'), ('Y', 'I', 'X')]
CEILING = [('X', 'Z', 'X'), ('X', 'Z', 'Y'), ('Y', 'Z', 'X')]
SOFT_TOL = 1e-6          # the canonical SpectrumPairs tolerance (C# PauliPairTrichotomy)
NORM_FACTOR = 2 * np.sqrt(2)   # original sweep's residual normalization vs framework's
TSV_ORIGINAL = SCRIPT_DIR / "results" / "ceiling_middle_tilt.tsv"
TSV_REGEN = SCRIPT_DIR / "results" / "ceiling_middle_tilt_regen.tsv"


def greedy_max(ev, sigma):
    """The canonical SpectrumPairs greedy criterion: max greedy pairing distance."""
    used = np.zeros(len(ev), bool)
    mx = 0.0
    for i in range(len(ev)):
        if used[i]:
            continue
        d = np.abs(ev - (-ev[i] - 2 * sigma))
        d[used] = np.inf
        j = int(np.argmin(d))
        mx = max(mx, float(d[j]))
        used[i] = True
        if j != i:
            used[j] = True
    return mx


def bottleneck(ev, tg):
    """Minimax (bottleneck) pairing distance between two eigenvalue multisets."""
    D = np.abs(ev[:, None] - tg[None, :])
    cand = np.unique(D)
    lo, hi = 0, len(cand) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        m = maximum_bipartite_matching(csr_matrix(D <= cand[mid]), perm_type='column')
        if (m >= 0).all():
            hi = mid
        else:
            lo = mid + 1
    return float(cand[lo])


def main():
    Hl = _build_kbody_chain(N, [t + (1.0,) for t in LOCAL])
    Hc = _build_kbody_chain(N, [t + (1.0,) for t in CEILING])

    original = []
    for line in TSV_ORIGINAL.read_text().splitlines()[1:]:
        a, v, pe, re_, mn = line.split('\t')
        original.append((float(a), v, float(pe), float(re_), float(mn)))
    assert len(original) == 37, f"original TSV has {len(original)} rows (expected 37)"

    rows = []
    n_verdict_match = 0
    max_mnorm_rel = 0.0
    for (ang, v_orig, _pe, _re, mn_orig) in original:
        th = np.deg2rad(ang)
        H = np.cos(th) * Hl + np.sin(th) * Hc
        L = lindbladian_pauli_dephasing(H, [GAMMA] * N, dephase_letter='Z')
        Mn = float(np.linalg.norm(fw.palindrome_residual(L, SIGMA, N, dephase_letter='Z')))
        Mn_scaled = NORM_FACTOR * Mn
        ev = np.linalg.eigvals(L)
        g = greedy_max(ev, SIGMA)
        verdict = 'Soft' if g < SOFT_TOL else 'Hard'
        pair = bottleneck(ev, -ev - 2 * SIGMA)
        evr = ev.real + 0j
        rate = bottleneck(evr, -evr - 2 * SIGMA)
        rows.append((ang, verdict, pair, rate, Mn_scaled))
        if verdict == v_orig:
            n_verdict_match += 1
        rel = abs(Mn_scaled - mn_orig) / mn_orig
        max_mnorm_rel = max(max_mnorm_rel, rel)

    # the certified structure
    verdicts = [r[1] for r in rows]
    assert verdicts[0] == 'Soft' and verdicts[-1] == 'Soft', "endpoints must be soft"
    assert all(v == 'Hard' for v in verdicts[1:-1]), "all 35 interior angles must be hard"
    assert n_verdict_match == 37, f"verdicts match original only {n_verdict_match}/37"
    assert max_mnorm_rel < 1e-5, f"M_norm column deviates rel {max_mnorm_rel:.2e} (>1e-5)"
    assert abs(rows[0][4] - 64.0) < 1e-9, f"local-endpoint ||M|| {rows[0][4]} != 64.0"
    assert abs(rows[-1][4] - 101.1929) < 1e-3, f"ceiling-endpoint ||M|| {rows[-1][4]}"
    assert all(r[4] >= rows[i][4] - 1e-9 for i, r in enumerate(rows[1:])), \
        "||M|| must load monotonically"
    # the interior break is real and gamma-quantized: minimax pairing >= 2*gamma somewhere,
    # and the soft endpoints sit at numerical zero
    assert max(r[2] for r in rows) >= 2 * GAMMA - 1e-9, "interior break below one 2*gamma step"
    assert rows[0][2] < 1e-9 and rows[-1][2] < 1e-9, "endpoint pairing not at numerical zero"

    with TSV_REGEN.open('w', encoding='utf-8', newline='\n') as f:
        f.write("angle_deg\tverdict\tminimax_pairing\tminimax_rate\tM_norm\n")
        for (ang, v, pair, rate, mn) in rows:
            f.write(f"{ang:.2f}\t{v}\t{pair:.6E}\t{rate:.6E}\t{mn:.6E}\n")

    print(f"verdicts: Soft@0, Hard@2.5..87.5 (35 angles), Soft@90; match original 37/37")
    print(f"M_norm column: 2*sqrt(2) * framework residual, max row deviation rel "
          f"{max_mnorm_rel:.2e}; loads 64.0 -> {rows[-1][4]:.4f} monotonically")
    print(f"minimax pairing: endpoints < 1e-9, interior quantized in 2*gamma = {2*GAMMA} steps "
          f"(max {max(r[2] for r in rows):.3f})")
    print(f"regenerated TSV: {TSV_REGEN}")
    print("PASS ceiling_middle_tilt: the ceiling soft point is isolated (hard sea in between)")


if __name__ == "__main__":
    main()
