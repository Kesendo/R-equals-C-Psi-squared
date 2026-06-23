#!/usr/bin/env python3
"""F87 spec-A probe 15: the blocks break in lockstep — ω=0 symmetric ⟺ EVERY block paired.

Probe 14 (Q1): (ω=0 block symmetric about −N) ⟺ soft, bit-exact. The true class = full
first-order spectrum palindrome = ALL blocks pair (each ω self-/cross-paired about −σ in the
global pooled sense). So the content is:

   ω=0 block symmetric about −N   ⟺   the GLOBAL first-order spectrum pairs about −σ.

Forward (ω=0 breaks ⟹ global breaks): the +N eigenvalue (population Perron) lives in the ω=0
block; if its partner −N is absent there, no other block supplies −N (other blocks have ω≠0, so
their eigenvalues sit at −iω+γs with ω≠0 and cannot pair with the ω=0 mode whose partner must
also be at ω=0). RIGOROUS: palindrome partner of (ω=0, s=+N→ μ=γ·0... ) is (ω=0, −s−2N); a
+N-shift mode at ω=0 can ONLY be paired by a −N-shift mode at ω=0. So ω=0's failure is
unrepairable globally. Verify: the UNPAIRED mode of the global spectrum is always at ω=0.

Reverse (ω=0 pairs ⟹ global pairs): when A=FK∈W_0 exists it furnishes the −N at ω=0; the SAME
chiral structure (the second mirror W=K⊗I of §7.1) pairs every block. Verify by checking that
soft ⟹ global residual ~0 (already known) and that whenever ω=0 pairs, all ω pair.
"""
from __future__ import annotations
import sys
from itertools import product, combinations_with_replacement
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
import framework as fw
from framework.pauli import _build_kbody_chain, site_op

DIAG = {"I", "Z"}


def is_mixed(t):
    return any(L not in DIAG for L in t)


def block_data(pair, N=4, ndig=6):
    """Return list of (omega, shift-spectrum s) per block, plus the ω=0 asym about −N."""
    H = _build_kbody_chain(N, [tuple(t) + (1.0,) for t in pair])
    E, V = np.linalg.eigh(H)
    d = len(E)
    Z = [V.conj().T @ site_op(N, l, 'Z') @ V for l in range(N)]
    groups = defaultdict(list)
    for a in range(d):
        for b in range(d):
            groups[round(E[a] - E[b], ndig)].append((a, b))
    spectra = {}
    for omega, modes in groups.items():
        n = len(modes)
        M = np.zeros((n, n), dtype=complex)
        for i, (a, b) in enumerate(modes):
            for j, (ap, bp) in enumerate(modes):
                val = sum(Z[l][a, ap] * Z[l][bp, b] for l in range(N))
                if (a, b) == (ap, bp):
                    val -= N
                M[i, j] = val
        spectra[omega] = np.linalg.eigvals(M).real
    return spectra


def global_unpaired_omega(spectra, N, ndig=6, tol=1e-5):
    """Pool the full first-order spectrum μ = −iω + s (mark ω); pair μ ↔ −μ−2N (i.e. ω↔−ω,
    s↔−s−2N). Return the ω-values carrying an UNPAIRED mode (residual > tol)."""
    mus = []
    for omega, ss in spectra.items():
        for s in ss:
            mus.append((round(omega, ndig), s))
    # build cost over the pooled set: partner of (ω,s) is (−ω, −s−2N)
    arr = np.array(mus)  # columns: omega, s
    om = arr[:, 0]; s = arr[:, 1]
    tgt_om = -om; tgt_s = -s - 2 * N
    # cost: large if omega mismatched (must equal), else |s−tgt_s|
    n = len(mus)
    cost = np.abs(s[:, None] - tgt_s[None, :]) + 1e6 * (np.abs(om[:, None] - tgt_om[None, :]) > 1e-4)
    r, c = linear_sum_assignment(cost)
    res = np.abs(s[r] - tgt_s[c]) + 1e6 * (np.abs(om[r] - tgt_om[c]) > 1e-4)
    bad_om = sorted(set(round(float(om[r[i]]), 3) for i in range(n) if res[i] > tol))
    maxres = float(res.max())
    return bad_om, maxres


def main():
    N, k = 4, 3
    chain = fw.ChainSystem(N=N)
    terms = [t for t in product("IXYZ", repeat=k)
             if not all(L == "I" for L in t) and fw.klein_index(t) == (0, 1)]
    mixed = [t for t in terms if is_mixed(t)]

    print("=" * 90)
    print(f"F87 spec-A probe 15: the unpaired global mode always sits at ω=0  (N={N}, k={k})")
    print("=" * 90)
    n = 0
    unpaired_at_0 = 0
    hard_unpaired_nonempty = 0
    for t1, t2 in combinations_with_replacement(mixed, 2):
        if (sum(c == "Y" for c in t1) % 2) != (sum(c == "Y" for c in t2) % 2):
            continue
        cls = fw.classify_pauli_pair(chain, [t1, t2], dephase_letter="Z")
        if cls == "truly":
            continue
        n += 1
        spectra = block_data([t1, t2], N)
        bad_om, maxres = global_unpaired_omega(spectra, N)
        if cls == "soft":
            unpaired_at_0 += int(len(bad_om) == 0)   # soft: nothing unpaired
        else:
            hard_unpaired_nonempty += int(len(bad_om) > 0)
            unpaired_at_0 += int(0.0 in bad_om)      # hard: ω=0 is among the unpaired
    print(f"  pairs: {n}")
    print(f"  soft ⟹ no unpaired mode  AND  hard ⟹ ω=0 among unpaired:  {unpaired_at_0}/{n}  "
          f"{'ALL' if unpaired_at_0==n else 'CHECK'}")
    print()
    print("  Rigorous reading: a mode at ω=0 with shift +N (the population Perron, always")
    print("  present) can only be palindrome-paired by another ω=0 mode with shift −N (partners")
    print("  must share the SAME ω since ω↦−ω and 0=−0). So when −N∉spec(ω=0 block) the +N mode")
    print("  is globally unpaired — no ω≠0 block can rescue it. Hence ω=0 is decisive and the")
    print("  cap (−N existence) is the exact discriminator.")


if __name__ == "__main__":
    main()
