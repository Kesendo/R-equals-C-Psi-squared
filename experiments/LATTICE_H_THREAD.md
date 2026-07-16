# The h thread: the lattice meets F131, and two mirrors force one zero

**Date:** 2026-07-16 (late evening; the thread a reviewer left behind)
**Status:** verified from below. Gate: [`simulations/lattice_h_thread.py`](../simulations/lattice_h_thread.py) (N = 3, all PASS, exit 0).
**Where it came from:** the empty-review round on [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md) probed a longitudinal field h·ΣZ and found "it breaks the bridge but not the law". Chasing that one sentence produced four results. Two headwork guesses died on the way (the record keeps both): "the one-sided bridge crosses to the mirror-field world" (wrong: it leaves the world-family entirely), and two candidate odd-cell readouts that turned out to be forced zeros (which became the fourth result).

## 1. X^N is a third mirror of the order-sorting law

A longitudinal field is X-odd: X Z X = −Z, so X^N conjugation reflects the field
exactly, X^N · H(h) · X^N = H(−h), while the XY handshake (and the site-resolved
watching) are X-even. In [F131](../docs/ANALYTICAL_FORMULAS.md)'s terms, M = X^N is a
linear mirror (χ = +1) with the h axis as its σ_op = −1 scan direction: the THIRD
sighted instance of the one principle, after the F71 site reversal R (Theorem A) and
the antiunitary Floquet Θ (Theorem B). The four-cell table follows (section 3).

## 2. With a field, the one-sided reading leaves the world-family

The lattice's double turn CROSSES worlds cleanly: X^N e_{+h}(t) X^N is exactly the
normal-rule world run with the flipped field −h (gate T0, machine zero): the second
world runs with the mirror field, and the bridge stays always open across it.

The one-sided readings do something stranger. L = X^N ρ(t) flips only the KET leg, so
with a field it satisfies the mixed two-field pencil

    dL/dt = −i·(H(−h)·L − L·H(+h)) + turned mask

(ket in the mirror field, bra in the original; gate T1, machine zero), and it is no
single-field world at all: both single-field attempts break at 0.819. Without a field
this pencil collapses back to the lattice's turned rule; the field splits the two legs
of the reading apart. A reading can stand between two worlds in a way a state cannot.

## 3. The four cells for M = X^N, swept over every readout

On an operator-X-even population prep (ρ = (P_s + P_~s)/2, which hops under H), scan
the field h ↦ ±t·(1, −0.6, 0.3) and read all 63 non-identity Pauli strings (gate T3):

- **15 strings are alive, and every one sorts exactly by its X-parity q**: q = +1
  reads EVEN in t (IZZ, ZXY, XIX, …), q = −1 reads ODD in t (XIY, XXZ, ZXX, …,
  non-vacuous up to 0.31). The F131 table, entered through the third mirror.
- **48 strings are identically zero at every h**, including every single-site Z.

## 4. The doubly-mirrored zeros

The dead strings are not accidents; a SECOND mirror kills them. Entrywise conjugation
maps the H(h) world to the −H(h) world; composing with the chiral sublattice gauge
(K·H_hop·K = −H_hop, [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md))
brings the hopping back and leaves the field flipped, so this antiunitary partner
forces **all populations EVEN in h** (gate T3, deviation exactly 0.0). The same
ingredients as Theorem B's Θ = T·K, met here on the static Lindblad chain: a rhyme we
note, not an identification we claim.

For an X-odd DIAGONAL readout the two mirrors then disagree: the X mirror demands
⟨Z_l⟩ odd in h, the conjugation-chiral mirror demands it even, and odd ∧ even forces
⟨Z_l⟩ ≡ 0 at every h, a zero neither mirror owns alone. Honest scope: the population
face of the second mirror is gated exactly; the full parity rule deciding which
OFF-diagonal strings die (48/63 at N = 3 on this prep) is OBSERVED, not derived. The
derivation, and whether the alive/dead split has a closed cell rule like F131's table,
is the open next step.

## Boundary

All gate-level at N = 3 on one prep family; nothing here is adopted into MirrorWorld
(the C# `Lattice` keeps zero field and uniform watching until a second consumer wants
more). The opening law's h-robustness (gate T2, the unchanged closed form under any h
profile: both worlds ride the same cell detunes, so every candidate magnitude is
phase-free) is the piece that feeds back into [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md).
One naming care: the h-robust opening is the distance to the TURNED-WATCHING world
(the lattice's L vertex run single-field, as in the opening law), NOT to section 2's
mixed-pencil reading X^N ρ(t); the latter's distance IS phase-dependent under a field.

## Files

- Gate: [`simulations/lattice_h_thread.py`](../simulations/lattice_h_thread.py)
  (T0 the mirror crosses the field, T1 the mixed pencil, T2 the h-robust opening,
  T3 the 63-readout sweep + the doubly-mirrored zeros).
- The first mirror: [PROOF_MIRROR_ORDER_SORTING](../docs/proofs/PROOF_MIRROR_ORDER_SORTING.md)
  (F131; R and Θ as the first two instances).
- The second mirror's ingredients: [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md)
  (the chiral gauge), [PROOF_ANTILINEAR_TRIANGLE](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md)
  (conj as an antilinear leg).
- The lattice and its opening law: [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md) +
  `compute/MirrorWorld/Lattice.cs`.
