# The zz face of the dead set: the set jumps, the magnitudes cross over

**Date:** 2026-07-16 (the evening's second face; the h face landed hours earlier)
**Status:** verified from below. Gate: [`simulations/lattice_dead_set_zz_face.py`](../simulations/lattice_dead_set_zz_face.py) (12 checks, all PASS, exit 0). Quantifies the zz boundary of [F132](../docs/ANALYTICAL_FORMULAS.md); no new F number.
**Where it came from:** [LATTICE_DEAD_SET_RULE](LATTICE_DEAD_SET_RULE.md) gated the zz boundary at a single coupling (zz = 0.7): V and F die, layer K persists, the lone fermionic survivor is the conserved parity Z^N. The handover carried the open question: how does the dead set LIFT as zz grows from zero, jump or crossover? With the explicit warning to measure rather than headwork (two guesses died on data in the parent arc). Measured: the answer is both, cleanly split by level.

## The question and the split answer

A readout that F132 declares dead is identically zero at every h in the free
world. Turn on a ZZ coupling and the two non-kinematic killers (the
doubly-mirrored V, the conserved degree F) break EXACTLY at any zz ≠ 0; the
necessity arguments predict a jump of the SET. But an exactly-broken symmetry
can still leave arbitrarily small magnitudes. Both faces are real:

**The set jumps.** At every sampled zz in {10⁻⁴, 10⁻³, 10⁻², 0.1, 0.7} the
alive set is already the full zz-large set, identical across four decades of
coupling: everything K-readable revives except the conserved parity Z^N at
odd N (which starts at 0 for the population pair and stays there at every
zz). Nothing newly dies. Counts: N=3 population 15 → 30 alive (ZZZ the sole
holdout), N=4 coherence 71 → 127 (no holdout: Z^N starts nonzero at even N),
N=5 s=3 coherence 507 → 1022 (ZZZZZ the holdout).

**The magnitudes cross over, with exact integer orders.** Every revived
readout obeys the revival-order law

    max_t |⟨O⟩(t)| ~ zz^m,   m = min |d(O) − d0| / 2
                             over populated prep degrees d0 of the same parity,

where d is the left-Jordan-Wigner Majorana degree of the readout and the
populated prep degrees are those of ρ(0)'s own expansion (d0 ∈ {0, 4, 8, ...}
from the population pair, plus d0 = N from the cross coherence). The
mechanism reading: H_zz is quartic in Majoranas, and a CONTRIBUTING zz vertex
moves the degree by exactly ±2 (the commutator survives only when the shared
Majorana count is odd; the degree-preserving and ±4 terms commute and drop),
so m vertices bridge a gap of 2m and parity is preserved throughout. The gate checks the law PER READOUT (local log-log slopes on
zz = 10⁻⁴ → 10⁻³ → 10⁻²): worst deviation |slope − m| = 0.0006 across all
586 revived readouts of the three configs.

## The 20 ride second order

At N = 5 with the popcount-2 seed and coherence, the revival orders split
three ways and account for the pool exactly (495 + 20 = 515):

| group | count | degree gap | order |
|---|---|---|---|
| d ≡ 2 mod 4 (V-killed, population channel) | 255 | 2 (to d0 ∈ {0, 4, 8}) | zz¹ |
| d ∈ {3, 7} (V-killed, coherence channel) | 240 | 2 (to d0 = 5) | zz¹ |
| d ∈ {1, 9} (the F-killed 20: the a_k and their Γ-duals) | 20 | 4 (to d0 = 5) | zz² |

The 20 strings F alone killed in the free world are exactly the ones that
need TWO interaction vertices to wake: F's kill was the deepest, and the
revival order remembers it. (On pure population preps the odd-d readouts have
no same-parity populated degree at all: they are w-odd, hence already K-dead,
and the parity fence never becomes visible as a separate layer there.)

## The jump is detection-relative, gated

"The set jumps at zz = 10⁻⁴" is a statement about a threshold (alive =
max|⟨O⟩| > 10⁻¹⁰ on the RK4 window). The gate pins the honest version: at
zz = 10⁻⁵ the m = 2 strings sit at ~7·10⁻¹¹, just BELOW the threshold (a thin
1.4× margin; the crossing point moves with the window and the floor, which is
exactly the point: it is a property of the detector, not of the physics),
while every m = 1 string is already alive, so the still-dead set at zz = 10⁻⁵
is exactly {Z^N} ∪ {the 20}. Any apparent growth of the set with zz is a
threshold artifact of the zz^m magnitudes. The two true statements are: the
identically-zero property is lost at every zz ≠ 0 (the set, in the exact
sense, jumps), and every magnitude grows continuously as zz^m (the
observable, at any fixed detection floor, crosses over in order of m).

## What this adds to F132

- The zz boundary was a single-point control (zz = 0.7); it is now a
  quantified lift: WHERE each silence goes when the interaction turns on,
  and how fast. The free world's dead set is not erased by the interaction,
  it is graded by it: each dead readout revives at the perturbative order
  set by its Majorana-degree distance from the prep.
- The sole survivor is kinematic on both sides: layer K persists (it is
  h- and zz-blind), and Z^N survives not as a symmetry kill but as a
  conserved quantity pinned to its initial value 0 at odd N. Z^N is also the
  honest scope marker for the order law: the formula assigns it a finite m
  (at N = 5, d = 10, gap 2 to d0 = 8), yet it never revives, so m = min
  |d − d0|/2 is a LOWER BOUND on the revival order whose prefactor can
  vanish, and Z^N is the one witnessed vanishing (by conservation, at every
  order). Every other pool member's prefactor was nonzero in the gated
  configs.
- Scope, honestly: three configs (N = 3 population, N = 4 coherence, N = 5
  s = 3 coherence), one watching profile, one generic h direction (the
  parent gate's), zz sampled over four decades plus the 10⁻⁵ threshold
  probe. The revival-order law is a gated observation with a perturbative
  mechanism reading; no all-orders derivation is claimed.

## Files

- Gate: [`simulations/lattice_dead_set_zz_face.py`](../simulations/lattice_dead_set_zz_face.py)
  (T1 the zz = 0 law reference, T2 the full jump at every sampled zz incl.
  the Z^N holdouts, T3 the per-readout revival-order law and the 20-string
  m = 2 pin, T4 the detection-relativity probe at zz = 10⁻⁵).
- The law whose boundary this quantifies: [LATTICE_DEAD_SET_RULE](LATTICE_DEAD_SET_RULE.md)
  (F132, the zz boundary section).
- The other face of the same axis walk: [LATTICE_DEAD_SET_H_ZERO](LATTICE_DEAD_SET_H_ZERO.md)
  (sufficiency generic in h; the h = 0 bi-degree refinement).
