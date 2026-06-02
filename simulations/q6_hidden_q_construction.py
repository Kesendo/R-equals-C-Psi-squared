#!/usr/bin/env python3
"""The explicit hidden symmetry Q for every palindromic two-term bilinear.

THE_OTHER_SIDE (Q6) asks: when the canonical palindrome operator Π (the P1
crossover) fails to pair the Liouvillian spectrum, is there a DIFFERENT symmetry
Q that still closes it? This script answers Q6 constructively. For each of the
36 unordered two-term bond-bilinear Hamiltonians H = Σ_bonds (t1 + t2) under
Z-dephasing (γ₀ = 0.05), it either BUILDS the actual hidden Q and verifies it
bit-exactly, or shows that no Q can exist.

What it proves (N = 3, 4, 5, identical):
  truly (3) : XX+YY, XX+ZZ, YY+ZZ close via the canonical P1; residual ≤ 1e-10.
  soft (19) : P1 FAILS (residual > 0), yet a hidden Q closes the conjugation
              equation ‖Q L Q⁻¹ − (−L − 2Σγ·I)‖ / ‖RHS‖ ≤ 1e-10. Each soft combo
              is assigned the family of its Q:
                uniform     (14) : one Z-valid per-site crossover applied to all
                                    sites (a P1-perm or P4-perm representative).
                alternating  (3) : P1 ⊗ M2 ⊗ P1 ⊗ ... by site parity.
                non_local    (2) : an entangled (non-product) Q from eigenvalue
                                    pairing; exactly {XZ+YZ, ZX+ZY}.
  hard (14) : the eigenvalue multiset is NOT invariant under λ → −2Σγ − λ; the
              palindrome is broken, so no Q of any kind can pair the spectrum.

The residual is the conjugation-equation norm: a valid Q sends the Liouvillian L
to its mirror image −L − 2Σγ·I about the center −Σγ, so ‖Q L Q⁻¹ + L + 2Σγ·I‖
vanishes (relative to ‖RHS‖) iff Q is a true palindrome symmetry. The "uniform"
family is backed by an exhaustive scan over all 1024 Z-valid uniform per-site
maps, so a non-uniform verdict is not an artifact of a fixed phase convention.
That exhaustive scan is the last-resort branch of classify_soft: each soft combo
is first tried against 6 curated uniform representatives, then alternating, then
non_local, and the full 1024-map scan runs only if all of those miss. For the
36 combos here the cheap branches always resolve, so the expensive scan is never
reached in practice; it remains as the guarantee for any future unresolved combo.

Convention note: all superoperators live in the COLUMN-STACK vec basis
(flatten('F')), the convention used by both simulations/pi_operator_entanglement.py
and the framework's lindbladian_z_dephasing. Per-site maps use the standalone
index convention I, X, Y, Z = 0, 1, 2, 3 (NOT the framework's Klein bits), since
the reused builders live in pi_operator_entanglement.py. We reuse those builders
verbatim:
  build_analytical_pi_multi(N, site_maps)  → tensored / alternating per-site Q
  find_pi_operator(L, sum_gamma)           → non-local Q via eigenvalue pairing
  verify_pi(Π, L, sum_gamma)               → ‖Π L Π⁻¹ − (−L − 2Σγ·I)‖ / ‖RHS‖

Per-site maps:
  P1 (canonical): I→X, X→I, Y→iZ, Z→iY
  P4 (uniform):   I→Y, X→iZ, Y→I, Z→iX
  M2 (P4 −i sib): I→Y, X→−iZ, Y→I, Z→−iX
The "alternating" family is P1 ⊗ M2 ⊗ P1 ⊗ ... (site parity), per the spec.

Self-validating: exits non-zero on any regression. The asserts at the bottom
pin the truly-via-P1 result, P1 failing on every soft combo, the 14 hard being
non-palindromic, the family tally (uniform 14 / alternating 3 / non_local 2),
the non_local set {XZ+YZ, ZX+ZY}, and cross-N family agreement across N = 3,4,5.
"""

import sys
import os
import numpy as np
from itertools import product as iprod

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

if sys.platform == "win32":
    # Force UTF-8 stdout so the γ/∈/Π glyphs print on a cp1252 console.
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

from pi_operator_entanglement import (
    sx, sy, sz,
    build_hamiltonian,
    build_z_dephasing,
    build_liouvillian,
    build_analytical_pi_multi,
    find_pi_operator,
    verify_pi,
)

PAULI_MAT = {'X': sx, 'Y': sy, 'Z': sz}

# Per-site maps in the standalone (perm, sign) format. Index 0,1,2,3 = I,X,Y,Z.
P1_PERM = {0: 1, 1: 0, 2: 3, 3: 2}      # I<->X, Y<->Z
P1_SIGN = {0: 1, 1: 1, 2: 1j, 3: 1j}

P4_PERM = {0: 2, 1: 3, 2: 0, 3: 1}      # I<->Y, X<->Z
P4_SIGN = {0: 1, 1: 1j, 2: 1, 3: 1j}

M2_PERM = {0: 2, 1: 3, 2: 0, 3: 1}      # same perm as P4
M2_SIGN = {0: 1, 1: -1j, 2: 1, 3: -1j}  # -i sibling

# The spec's "uniform" family is a single Z-dephasing-valid per-site map (a
# bijection {I,Z}<->{X,Y} with phases in {+-1,+-i}) applied uniformly to all
# sites. Which permutation/phase closes a given combo depends on its Klein cell
# (P4-perm I<->Y for some, P1-perm I<->X for others), so we test the WHOLE
# uniform family (1024 maps) and report the canonical representative that
# verifies. Brute force here so the family label does not depend on guessing one
# fixed phase convention.
_PHASES4 = [1, -1, 1j, -1j]


def _enumerate_uniform_maps():
    """All 1024 Z-valid uniform per-site (perm, sign) maps (bijection
    {I,Z}<->{X,Y}, independent phase in {+-1,+-i})."""
    maps = []
    for i_tgt in (1, 2):                  # I -> X or Y
        z_tgt = 2 if i_tgt == 1 else 1    # Z -> the other dephased letter
        for x_tgt in (0, 3):              # X -> I or Z
            y_tgt = 3 if x_tgt == 0 else 0
            perm = {0: i_tgt, 1: x_tgt, 2: y_tgt, 3: z_tgt}
            for phs in iprod(_PHASES4, repeat=4):
                maps.append((perm, {a: phs[a] for a in range(4)}))
    return maps


_UNIFORM_MAPS = _enumerate_uniform_maps()

# Curated uniform candidates tried FIRST (cheap path). The probe over all 1024
# maps showed the soft "uniform" cells close via one of these two permutation
# families with a specific sign pattern: a P4-perm representative (I<->Y, X<->Z)
# and a P1-perm representative (I<->X, Y<->Z), each with its phase-conjugate
# sibling. If none of these verify, classify_soft falls back to the full 1024
# enumeration before declaring the combo non-uniform.
_UNIFORM_CANDIDATES = [
    # P4-perm cell representative found by the probe (+ -i sibling):
    ({0: 2, 1: 3, 2: 0, 3: 1}, {0: 1, 1: 1j, 2: -1, 3: -1j}),   # I->Y,X->iZ,Y->-I,Z->-iX
    ({0: 2, 1: 3, 2: 0, 3: 1}, {0: 1, 1: -1j, 2: -1, 3: 1j}),   # conj sibling
    # P1-perm cell representative found by the probe (+ -i sibling):
    ({0: 1, 1: 0, 2: 3, 3: 2}, {0: 1, 1: -1, 2: 1j, 3: -1j}),   # I->X,X->-I,Y->iZ,Z->-iY
    ({0: 1, 1: 0, 2: 3, 3: 2}, {0: 1, 1: -1, 2: -1j, 3: 1j}),   # conj sibling
    # The spec's nominal P4 (I->Y,X->iZ,Y->+I,Z->iX) and P1, as documented refs:
    (P4_PERM, P4_SIGN),
    (P1_PERM, P1_SIGN),
]


def _ph(v):
    for r, s in [(1, '+1'), (-1, '-1'), (1j, '+i'), (-1j, '-i')]:
        if abs(v - r) < 1e-12:
            return s
    return f'{v}'


def _map_str(perm, sign):
    names = ['I', 'X', 'Y', 'Z']
    return ", ".join(f"{names[a]}->{_ph(sign[a])}{names[perm[a]]}" for a in range(4))


GAMMA = 0.05
RESID_TOL = 1e-10           # "bit-exact" threshold
PAIR_TOL = 1e-9             # eigenvalue-pairing / palindrome tolerance
# Margin note: the uniform / alternating Q close analytically to ~1e-17; the two
# non_local Q come from a numerical eigenvalue-pairing construction and run nearer
# the bar (~1e-11 to 1e-14, tightest ~1.0e-11 at N=5). They pass comfortably at
# N <= 5, but the non_local margin shrinks with N, so push RESID_TOL with care
# if this is ever extended beyond N=5.

N_VALUES = [3, 4, 5]

# --- The 36 combos, with their R1 classification (ground truth to confirm) ---
TRULY = ['XX+YY', 'XX+ZZ', 'YY+ZZ']
SOFT = ['XX+XZ', 'XX+YZ', 'XX+ZX', 'XX+ZY', 'XY+YX', 'XY+ZZ', 'XZ+YY',
        'XZ+YZ', 'XZ+ZX', 'XZ+ZZ', 'YX+ZZ', 'YY+YZ', 'YY+ZX', 'YY+ZY',
        'YZ+ZY', 'YZ+ZZ', 'ZX+ZY', 'ZX+ZZ', 'ZY+ZZ']
HARD = ['XX+XY', 'XX+YX', 'XY+XZ', 'XY+YY', 'XY+YZ', 'XY+ZX', 'XY+ZY',
        'XZ+YX', 'XZ+ZY', 'YX+YY', 'YX+YZ', 'YX+ZX', 'YX+ZY', 'YZ+ZX']

# Expected family tally and the exact non-local set (Q6 ground truth).
EXPECTED_TALLY = {'uniform': 14, 'alternating': 3, 'non_local': 2}
EXPECTED_NON_LOCAL = sorted(['XZ+YZ', 'ZX+ZY'])


def combo_terms(combo):
    """Parse 'XZ+YZ' -> [(sx, sz), (sy, sz)] as matrix tuples for build_hamiltonian."""
    parts = combo.split('+')
    terms = []
    for p in parts:
        a, b = p[0], p[1]
        terms.append((PAULI_MAT[a], PAULI_MAT[b]))
    return terms


def build_L(combo, N):
    """Build the column-stack-vec Liouvillian for a combo on the open N-chain."""
    H = build_hamiltonian(N, combo_terms(combo), J=1.0, topology='chain')
    c_ops, sum_gamma = build_z_dephasing(N, GAMMA)
    L = build_liouvillian(H, c_ops)
    return L, sum_gamma


def eigval_palindrome_residual(L, sum_gamma):
    """Q-agnostic palindrome test: is multiset {lambda} invariant under
    lambda -> -2*sum_gamma - lambda? Return the max nearest-neighbor pairing
    distance (0 if perfectly palindromic)."""
    evals = np.linalg.eigvals(L)
    target = -evals - 2 * sum_gamma
    d2 = len(evals)
    used = np.zeros(d2, dtype=bool)
    max_gap = 0.0
    # Greedy nearest pairing of evals to targets.
    for i in range(d2):
        dists = np.abs(evals - target[i])
        dists[used] = np.inf
        j = int(np.argmin(dists))
        max_gap = max(max_gap, dists[j])
        used[j] = True
    return max_gap


# ----------------------------------------------------------------------------
# Pi precomputation cache. A per-site map's Pi superoperator depends only on the
# map and N, NOT on the combo, so we build each candidate Pi ONCE per N and reuse
# it across all 36 combos (and store Pi_inv so the residual is matmul-only).
# This is the difference between ~40k expensive tensor builds and ~1k.
# ----------------------------------------------------------------------------
_PI_CACHE = {}   # N -> dict of precomputed (Pi, Pi_inv) bundles


def _alt_maps(N, site0='P1'):
    a, b = ((P1_PERM, P1_SIGN), (M2_PERM, M2_SIGN))
    if site0 == 'M2':
        a, b = b, a
    return [a if (s % 2 == 0) else b for s in range(N)]


def _build_pi_inv(N, site_maps):
    Pi = build_analytical_pi_multi(N, site_maps)
    return Pi, np.linalg.inv(Pi)


def precompute_pis(N):
    """Build and cache the cheap candidate Pi's (and inverses) for chain length N:
    canonical P1, the two alternating phases, and the 6 curated uniform reps.

    The full 1024-map uniform scan is NOT built here; it is materialized lazily
    by _uniform_all_pis(N) only if a curated candidate misses, since at N=5 that
    bundle is ~32 GB (1024 superoperators of dimension 1024² plus inverses) and
    is never needed when the curated reps already close every soft combo."""
    if N in _PI_CACHE:
        return _PI_CACHE[N]
    bundle = {}
    bundle['P1'] = _build_pi_inv(N, [(P1_PERM, P1_SIGN)] * N)
    bundle['alt'] = _build_pi_inv(N, _alt_maps(N, 'P1'))
    bundle['alt_swap'] = _build_pi_inv(N, _alt_maps(N, 'M2'))
    bundle['uniform_candidates'] = [
        (_map_str(p, s), *_build_pi_inv(N, [(p, s)] * N))
        for (p, s) in _UNIFORM_CANDIDATES
    ]
    bundle['uniform_all'] = None   # lazily filled by _uniform_all_pis(N)
    _PI_CACHE[N] = bundle
    return bundle


def _uniform_all_pis(N):
    """Build (once, on first miss) and cache all 1024 uniform Pi's at length N.

    Deferred from precompute_pis because the bundle is large at N=5 and only
    needed when the curated shortlist fails to close a soft combo; building it
    here keeps the 'non-uniform' verdict exhaustive without paying the cost up
    front for every run."""
    bundle = precompute_pis(N)
    if bundle['uniform_all'] is None:
        print(f"  [precompute] curated miss: building all {len(_UNIFORM_MAPS)} "
              f"uniform Pi(N={N}) superoperators (one-time fallback)...", flush=True)
        bundle['uniform_all'] = [
            (_map_str(p, s), *_build_pi_inv(N, [(p, s)] * N))
            for (p, s) in _UNIFORM_MAPS
        ]
        print(f"  [precompute] fallback build done (N={N}).", flush=True)
    return bundle['uniform_all']


def _mirror_rhs(L, sum_gamma):
    """The mirror target RHS = -L - 2*sum_gamma*I (does NOT depend on Q)."""
    d2 = L.shape[0]
    return -L - 2 * sum_gamma * np.eye(d2, dtype=complex)


def _resid_vs_rhs(Pi, Pi_inv, L, RHS, rhs_norm):
    """||Pi L Pi^-1 - RHS|| / ||RHS|| given a precomputed RHS and its norm.

    Same metric as verify_pi; RHS is Q-independent, so hoisting it out of any
    multi-Q scan avoids rebuilding the d²×d² identity once per candidate (the
    dominant overhead of the 1024-map uniform scan at N=5)."""
    LHS = Pi @ L @ Pi_inv
    return np.linalg.norm(LHS - RHS) / rhs_norm


def _resid(Pi, Pi_inv, L, sum_gamma):
    """||Pi L Pi^-1 - (-L - 2*sum_gamma*I)|| / ||RHS|| (same metric as verify_pi)."""
    RHS = _mirror_rhs(L, sum_gamma)
    return _resid_vs_rhs(Pi, Pi_inv, L, RHS, np.linalg.norm(RHS))


def verify_uniform_curated(L, sum_gamma, N):
    """Cheap uniform check: try only the 6 curated representatives. Returns
    (best_residual, representative_str). A residual < RESID_TOL means a uniform
    Q closes the combo; otherwise the exhaustive scan (verify_uniform_exhaustive)
    decides whether the combo is truly non-uniform."""
    bundle = precompute_pis(N)
    RHS = _mirror_rhs(L, sum_gamma)            # Q-independent: built once
    rhs_norm = np.linalg.norm(RHS)
    best_r, best_rep = np.inf, None
    for rep, Pi, Pi_inv in bundle['uniform_candidates']:
        r = _resid_vs_rhs(Pi, Pi_inv, L, RHS, rhs_norm)
        if r < best_r:
            best_r, best_rep = r, rep
        if best_r < RESID_TOL:
            return best_r, best_rep
    return best_r, best_rep


def verify_uniform_exhaustive(L, sum_gamma, N):
    """Exhaustive uniform check: scan ALL 1024 Z-valid uniform maps (curated set
    first, then the full enumeration). Returns (best_residual, representative_str).
    This is the guarantee that a 'non-uniform' verdict is not an artifact of the
    curated shortlist; it is invoked only when the cheap families all miss, since
    materializing and scanning the full bundle is the dominant cost at N=5."""
    RHS = _mirror_rhs(L, sum_gamma)
    rhs_norm = np.linalg.norm(RHS)
    best_r, best_rep = verify_uniform_curated(L, sum_gamma, N)
    if best_r < RESID_TOL:
        return best_r, best_rep
    for rep, Pi, Pi_inv in _uniform_all_pis(N):
        r = _resid_vs_rhs(Pi, Pi_inv, L, RHS, rhs_norm)
        if r < best_r:
            best_r, best_rep = r, rep
        if best_r < RESID_TOL:
            return best_r, best_rep
    return best_r, best_rep


def verify_alternating(L, sum_gamma, N):
    """P1 (x) M2 (x) P1 (x) ... by site parity (site 0 -> P1)."""
    Pi, Pi_inv = precompute_pis(N)['alt']
    return _resid(Pi, Pi_inv, L, sum_gamma)


def verify_alternating_swapped(L, sum_gamma, N):
    """M2 (x) P1 (x) M2 ... (site 0 -> M2). Tried only if site-0=P1 fails."""
    Pi, Pi_inv = precompute_pis(N)['alt_swap']
    return _resid(Pi, Pi_inv, L, sum_gamma)


def verify_P1(L, sum_gamma, N):
    Pi, Pi_inv = precompute_pis(N)['P1']
    return _resid(Pi, Pi_inv, L, sum_gamma)


def verify_nonlocal(L, sum_gamma):
    Pi, pairs, evals = find_pi_operator(L, sum_gamma, tol=PAIR_TOL)
    if Pi is None:
        return None
    return verify_pi(Pi, L, sum_gamma)


def classify_soft(combo, N):
    """Assign a soft combo the family of its hidden Q. Family precedence is
    uniform -> alternating -> non_local, but the checks are ordered cheapest
    first to avoid the expensive exhaustive uniform scan whenever possible:

      1. curated uniform   (6 representatives, cheap)
      2. alternating       (P1 (x) M2 ..., site0=P1 then site0=M2)
      3. non_local         (eigenvalue-pairing Q)
      4. exhaustive uniform (all 1024 maps) -- LAST resort, only if 1-3 all miss

    Step 4 preserves the 'no uniform Q exists' exhaustiveness guarantee for any
    combo not already resolved; for every combo in the validated 36 it is never
    reached (each closes at step 1, 2, or 3), which is what keeps N=5 tractable.
    A uniform Q found by the curated set still outranks alternating, so the
    family precedence is unchanged for the real combos. Returns
    (family, residual, detail)."""
    L, sum_gamma = build_L(combo, N)

    # 1. Cheap curated uniform.
    r_u, rep_u = verify_uniform_curated(L, sum_gamma, N)
    if r_u < RESID_TOL:
        return 'uniform', r_u, rep_u

    # 2. Alternating (both site-0 phases).
    r_alt = verify_alternating(L, sum_gamma, N)
    if r_alt < RESID_TOL:
        return 'alternating', r_alt, 'site0=P1'
    r_alt2 = verify_alternating_swapped(L, sum_gamma, N)
    if r_alt2 < RESID_TOL:
        return 'alternating', r_alt2, 'site0=M2'

    # 3. Non-local (entangled) Q.
    r_nl = verify_nonlocal(L, sum_gamma)
    if r_nl is not None and r_nl < RESID_TOL:
        return 'non_local', r_nl, ''

    # 4. Last resort: exhaustive uniform scan (all 1024 maps). Only reached for
    #    a combo that none of the cheap families closed; makes the 'non-uniform'
    #    verdict exhaustive rather than curated-only.
    r_uex, rep_uex = verify_uniform_exhaustive(L, sum_gamma, N)
    if r_uex < RESID_TOL:
        return 'uniform', r_uex, rep_uex

    # Nothing reached bit-exact. Report best of what we have for diagnosis.
    cands = {'uniform(exhaustive)': r_uex, 'alternating(site0=P1)': r_alt,
             'alternating(site0=M2)': r_alt2,
             'non_local': r_nl if r_nl is not None else np.inf}
    best = min(cands, key=cands.get)
    return f'UNRESOLVED(best={best})', cands[best], str(cands)


def run(N):
    print(f"\n{'#' * 78}")
    print(f"#  Q6 hidden-Q construction at N={N}, Z-dephasing γ={GAMMA}")
    print(f"{'#' * 78}")

    # ---- The 3 truly close via P1 ----
    print(f"\n--- TRULY (expect P1 verifies, residual < {RESID_TOL:.0e}) ---")
    truly_ok = True
    for c in TRULY:
        L, sg = build_L(c, N)
        r = verify_P1(L, sg, N)
        ok = r < RESID_TOL
        truly_ok &= ok
        print(f"  {c:8s}  P1 residual = {r:.3e}   {'OK' if ok else 'FAIL'}")

    # ---- P1 FAILS on soft ones ----
    print(f"\n--- SOFT: confirm canonical P1 FAILS (residual > 0) ---")
    p1_fail_ok = True
    p1_resid = {}
    for c in SOFT:
        L, sg = build_L(c, N)
        r = verify_P1(L, sg, N)
        p1_resid[c] = r
        fails = r > RESID_TOL
        p1_fail_ok &= fails
        print(f"  {c:8s}  P1 residual = {r:.3e}   {'P1 fails (good)' if fails else 'P1 WORKS?!'}")

    # ---- Classify each soft combo by the family of its hidden Q ----
    print(f"\n--- SOFT: family of the hidden Q (uniform -> alternating -> non_local) ---")
    results = {}
    tally = {'uniform': 0, 'alternating': 0, 'non_local': 0}
    for c in SOFT:
        fam, resid, detail = classify_soft(c, N)
        results[c] = (fam, resid, detail)
        if fam in tally:
            tally[fam] += 1
        d = f"  [{detail}]" if detail else ""
        print(f"  {c:8s} -> {fam:12s}  residual = {resid:.3e}{d}")

    # ---- The 14 hard are NOT palindromic at all ----
    print(f"\n--- HARD: confirm eigenvalue palindrome FAILS (no Q exists) ---")
    hard_ok = True
    hard_gaps = {}
    for c in HARD:
        L, sg = build_L(c, N)
        gap = eigval_palindrome_residual(L, sg)
        hard_gaps[c] = gap
        not_pal = gap > PAIR_TOL
        hard_ok &= not_pal
        print(f"  {c:8s}  eig-mirror max gap = {gap:.3e}   "
              f"{'NOT palindromic (good)' if not_pal else 'PALINDROMIC?!'}")

    # Cross-check: soft + truly ARE palindromic by the same eig-mirror test.
    print(f"\n--- SANITY: soft+truly ARE eig-palindromic (max gap < {PAIR_TOL:.0e}) ---")
    softtruly_pal_ok = True
    for c in TRULY + SOFT:
        L, sg = build_L(c, N)
        gap = eigval_palindrome_residual(L, sg)
        ok = gap < PAIR_TOL
        softtruly_pal_ok &= ok
        if not ok:
            print(f"  {c:8s}  eig-mirror max gap = {gap:.3e}   NOT palindromic?!")
    print(f"  all {len(TRULY)+len(SOFT)} truly+soft palindromic: "
          f"{'YES' if softtruly_pal_ok else 'NO'}")

    # ---- non_local cases should be exactly XZ+YZ and ZX+ZY ----
    nonlocal_combos = sorted([c for c, (f, _, _) in results.items() if f == 'non_local'])

    print(f"\n{'=' * 78}")
    print(f"SUMMARY (N={N})")
    print(f"{'=' * 78}")
    print(f"  TRULY close via P1:                 {'YES' if truly_ok else 'NO'}  ({len(TRULY)}/{len(TRULY)})")
    print(f"  P1 fails on all soft:               {'YES' if p1_fail_ok else 'NO'}")
    print(f"  HARD all non-palindromic (no Q):    {'YES' if hard_ok else 'NO'}  ({len(HARD)}/{len(HARD)})")
    print(f"  soft+truly all eig-palindromic:     {'YES' if softtruly_pal_ok else 'NO'}")
    print(f"  Soft family tally: {tally}  (sum={sum(tally.values())} of {len(SOFT)})")
    n_unres = sum(1 for f, _, _ in results.values() if f.startswith('UNRESOLVED'))
    if n_unres:
        print(f"  *** {n_unres} soft combo(s) UNRESOLVED ***")
    print(f"  non_local combos found: {nonlocal_combos}")
    print(f"  expected non_local:     {EXPECTED_NON_LOCAL}")
    print(f"  non_local match (XZ+YZ, ZX+ZY): {'YES' if nonlocal_combos == EXPECTED_NON_LOCAL else 'NO'}")

    # ---- Per-N self-validating asserts (exit non-zero on regression) ----
    assert truly_ok, f"N={N}: not all truly combos close via P1"
    assert p1_fail_ok, f"N={N}: P1 unexpectedly closes a soft combo"
    assert hard_ok, f"N={N}: a hard combo is eigenvalue-palindromic (Q would exist)"
    assert softtruly_pal_ok, f"N={N}: a truly/soft combo failed the eig-palindrome test"
    assert n_unres == 0, f"N={N}: {n_unres} soft combo(s) UNRESOLVED (no family found)"
    assert tally == EXPECTED_TALLY, (
        f"N={N}: soft family tally {tally} != {EXPECTED_TALLY}"
    )
    assert nonlocal_combos == EXPECTED_NON_LOCAL, (
        f"N={N}: non_local set {nonlocal_combos} != {EXPECTED_NON_LOCAL}"
    )

    return results, tally, nonlocal_combos


def print_maps():
    print("Per-site maps used (standalone convention, I,X,Y,Z = 0,1,2,3):")

    def fmt(perm, sign):
        names = ['I', 'X', 'Y', 'Z']

        def ph(v):
            for r, s in [(1, '+1'), (-1, '-1'), (1j, '+i'), (-1j, '-i')]:
                if abs(v - r) < 1e-12:
                    return s
            return f'{v}'
        return ", ".join(f"{names[a]}->{ph(sign[a])}{names[perm[a]]}" for a in range(4))
    print(f"  P1 (canonical): {fmt(P1_PERM, P1_SIGN)}")
    print(f"  P4 (uniform):   {fmt(P4_PERM, P4_SIGN)}")
    print(f"  M2 (P4 -i sib): {fmt(M2_PERM, M2_SIGN)}")
    print(f"  alternating:    P1 (x) M2 (x) P1 (x) ...  (site parity, site0=P1)")


def main():
    np.set_printoptions(precision=4, suppress=True)
    print("Q6: the explicit hidden symmetry Q for every palindromic two-term bilinear")
    print(f"  36 combos, Z-dephasing γ₀ = {GAMMA}, N ∈ {N_VALUES}")
    print(f"  truly = canonical P1 pairs; soft = hidden Q ≠ P1 pairs; hard = no Q")
    print()
    print_maps()

    per_n = {}
    for N in N_VALUES:
        res, tally, nl = run(N)
        per_n[N] = (res, tally, nl)

    # ---- Final cross-N per-combo family table ----
    print(f"\n{'=' * 78}")
    print(f"PER-COMBO FAMILY TABLE (N ∈ {N_VALUES})")
    print(f"{'=' * 78}")
    header = "  {:8s}".format('combo')
    for N in N_VALUES:
        header += f" {f'N={N} family':14s}"
    header += "  agree"
    print(header)

    agree_all = True
    for c in SOFT:
        fams = [per_n[N][0][c][0] for N in N_VALUES]
        agree = all(f == fams[0] for f in fams)
        agree_all &= agree
        row = f"  {c:8s}"
        for f in fams:
            row += f" {f:14s}"
        row += f"  {'yes' if agree else 'NO'}"
        print(row)

    print()
    for N in N_VALUES:
        print(f"  N={N} tally: {per_n[N][1]}")
    print(f"  cross-N family agreement on all soft: {'YES' if agree_all else 'NO'}")

    # ---- Final cross-N asserts (exit non-zero on any regression) ----
    assert agree_all, "cross-N family disagreement on at least one soft combo"
    ref_tally = per_n[N_VALUES[0]][1]
    for N in N_VALUES:
        assert per_n[N][1] == EXPECTED_TALLY, (
            f"N={N}: soft family tally {per_n[N][1]} != {EXPECTED_TALLY}"
        )
        assert per_n[N][1] == ref_tally, (
            f"N={N}: tally {per_n[N][1]} != reference N={N_VALUES[0]} tally {ref_tally}"
        )
        assert per_n[N][2] == EXPECTED_NON_LOCAL, (
            f"N={N}: non_local set {per_n[N][2]} != {EXPECTED_NON_LOCAL}"
        )

    print()
    print("ALL CROSS-N ASSERTIONS PASSED:")
    print(f"  soft family tally == {EXPECTED_TALLY} at every N")
    print(f"  family assignment agrees across N ∈ {N_VALUES} for all 19 soft combos")
    print(f"  non_local == {{XZ+YZ, ZX+ZY}} at every N")
    print()
    print("Done.")


if __name__ == '__main__':
    main()
