# The dead-set rule: three conserved structures close the h thread's open question

**Date:** 2026-07-16 (the h thread picked back up the same day)
**Status:** verified from below. Gate: [`simulations/lattice_dead_set_rule.py`](../simulations/lattice_dead_set_rule.py) (15 checks, all PASS, exit 0).
**Where it came from:** [LATTICE_H_THREAD](LATTICE_H_THREAD.md) §4 left one question open: of the 63 Pauli readouts on the h-scan setup, 48 are identically zero at every h, the diagonal deaths were proven (the doubly-mirrored zeros), but the rule deciding which OFF-diagonal strings die was observed, not derived. The overnight handover carried a hand-spotted candidate rule (alive iff n_Z ≡ #{X/Y letters on the gauge sublattice} mod 2) with an explicit warning not to trust it. The warning was right: the candidate as stated is falsified by the population prep at N = 3 and N = 4. Rebuilt from mechanism, it turns out to be one of THREE layers, and the three together are exact everywhere we tested.

## The setup

Open XY chain + longitudinal field h_l Z_l, local Z-dephasing with any watching
profile (one-sided and full non-uniform both gated), prep = the operator-X-even
population pair ρ(0) = (P_s + P_~s)/2 (X^N ρ X^N = ρ; ~s is the bit complement of
s), optionally with a cross coherence c·|s⟩⟨~s| + h.c. (real or imaginary c). Scan the field h ↦ ±t·(generic direction);
a readout O is ALIVE if ⟨O⟩(t) is not identically zero over the scan, DEAD otherwise.

## The rule

A Pauli readout O (n_Z Z-letters, XY-mask of weight w) is alive iff it passes all
three layers. Each necessity is an exact conservation or symmetry argument; the joint
sufficiency (everything not forbidden is alive) is the gated observation.

**K, the popcount blocks.** H conserves total popcount at every h and the dephasing
mask moves no matrix element, so ρ(t) is supported exactly on the (p, q) popcount
blocks populated at t = 0: the two diagonal blocks (p_s, p_s), (p_~s, p_~s), plus the
coherence blocks (p_s, p_~s), (p_~s, p_s) when c ≠ 0. O reads only the elements
(i, i⊕m) of its own mask m, so it must connect a supported block: an integer
a = (p + w − q)/2 with 0 ≤ a ≤ w, a ≤ p, p − a ≤ N − w must exist for some supported
(p, q). First bite beyond parity: the 16 pure-XY strings at N = 4 (population prep),
whose weight-4 mask can meet the population's diagonal blocks only through popcount
2, and the (2, 2) block is never populated.

**V, the doubly-mirrored kill.** The h thread's two h-FLIPPING mirrors compose to
h-PRESERVING symmetries. For each sublattice gauge U_g (Z on the even or on the odd
sites), W_g = U_g·X^N sends H(h) to −H(h) as an operator (the gauge flips the
hopping, X^N flips the field). Entrywise conjugation leaves the real H untouched,
but it is ANTIUNITARY: it reverses the i, so the conjugated trajectory ρ* solves the
flow of −H, i.e. conj maps the H(h) world to the −H(h) world; W_g maps that world
back to +H(h). The dissipator rides along for free: the dephasing mask is real and
depends only on i⊕j, so it is invariant under conjugation, under the X^N bit-flip,
and under the diagonal gauge signs. Hence the antiunitary

    V_g(ρ) = W_g · conj(ρ) · W_g†

is an exact symmetry of the full Lindblad flow at every FIXED h (the gate asserts
the worst-case invariance along the trajectories below 1e−12; measured ~1e−16). If V_g fixes ρ(0), then ⟨O⟩ = ε_g·⟨O⟩ with

    ε_g(O) = (−1)^(n_Z + xy_g),   xy_g = #{X/Y letters of O on sublattice g},

so ε_g = −1 kills O identically, at every h. Stabilizer bookkeeping, checked directly
in the gate: the population part is always fixed; the cross-coherence part is fixed
iff the sublattice size |g| is even (real and imaginary c alike, both gated). So both
gauges kill on a population-only prep, and exactly one gauge kills on a coherence
prep at odd N (the one with an even number of sites; it moves with N, the handover's
warning). The h thread's proven diagonal case is the w = 0 cell of this layer:
ε_g = (−1)^(n_Z) for populations, whence every X-odd diagonal readout dies.

**F, the fermion degree.** Under the left Jordan-Wigner map every Pauli string IS one
Majorana monomial with a definite degree d(O) (computable right-to-left by a two-line
tail-parity recursion; the gate carries it). The XY + field Hamiltonian is quadratic,
and the adjoint dissipator of Z-dephasing is DIAGONAL on Pauli strings, so the
Heisenberg flow preserves d exactly, sector by sector. The prep touches only even
degrees (populations expand in pure-Z strings, degree 2·|T|) plus degree N exactly
(the cross coherence expands in all-XY strings, one Majorana per site). So O with odd
d is dead unless the coherence is present and d(O) = N. This layer hides behind K and
V in most small-N configs and first becomes visible at N = 5 with a popcount-2 seed
(coherence blocks Δ = ±1): there the 2N Jordan-Wigner-linear strings (d = 1: an X or
Y cap carrying its Z-tail, the a_k themselves) and their Γ-duals (d = 2N−1: the
complementary monomials Γ·a_k up to phase, Γ = the product of all 2N Majoranas) slip
through K and V and are dead by F alone, exactly 20 of them, exactly as predicted.

## What the gate pins

All censuses are full 4^N − 1 sweeps; exactness means predicted set == alive set.

| Config | alive / total | layers biting |
|---|---|---|
| N=3 population (s=1) | 15/63 | V (both gauges) + K |
| N=3 real coh, imag coh (s=1,2) | 35/63 | V (even gauge only) |
| N=4 population (s=1) | 55/255 | V (both) + K (the 16 pure-XY) |
| N=4 real coh (s=1) | 71/255 | V (both) + K |
| N=5 real coh (s=1), incl. full non-uniform watching | 367/1023 | V (odd gauge) + K |
| N=5 population (s=3, popcount 2) | 255/1023 | V (both) + K |
| N=5 real coh (s=3) | 507/1023 | V (odd) + K + **F (the 20)** |

Layer attribution witnesses, each gated as a set equality: ZII at N=3 population
(V-only kill, the original doubly-mirrored zero), the 16 pure-XY strings at N=4
population (K-only), the 20 degree-{1, 2N−1} strings at N=5 s=3 coherence (F-only).
The "layers biting" column of the other rows is our reading of which layers act
there; what the gate pins per row is the full-census set equality itself.

**The zz boundary.** With a ZZ coupling the chiral gauge no longer flips H (V dies)
and H is quartic (F dies); layer K persists. Gated: at N=4 the dead set collapses to
exactly layer K (all odd-w strings, 128/255); at N=3 (population) it is layer K plus
the single string ZZZ, which is the conserved fermion parity Z^N whose expectation
starts at 0 for the odd-N pair (the two poles have opposite parity) and stays there.
The one fermionic survivor of the interaction is the parity, which is F's own shadow.
So V and F are properties of the FREE world; the two headwork lessons of the h thread
repeated here in one day: the candidate rule was falsified as stated before it was
rebuilt, and my "V survives zz" guess died on the same data that killed it.

## Redundancy structure (why the layers hid from each other)

On the odd-w sector, where F applies, ε_even·ε_odd = (−1)^w = −1, so at most one
gauge condition is independent there, and the surviving one appears to be the mod-4
shadow of the exact degree; that is why F stayed invisible until a Δ = ±1 coherence
seed let odd-w strings through K at degrees far from N. On the even-w sector F is
silent (every even degree is fed by the populations) and V does the killing. K is the
only layer that survives interactions. None of the three subsumes another.

## Files

- Gate: [`simulations/lattice_dead_set_rule.py`](../simulations/lattice_dead_set_rule.py)
  (T1 ten full censuses, T2 per-layer witnesses, T3 the zz boundary).
- The question this closes: [LATTICE_H_THREAD](LATTICE_H_THREAD.md) §4.
- The two composed mirrors: [PROOF_MIRROR_ORDER_SORTING](../docs/proofs/PROOF_MIRROR_ORDER_SORTING.md)
  (F131, X^N as its third sighting), [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md)
  (the chiral gauge), [PROOF_ANTILINEAR_TRIANGLE](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md)
  (conjugation as an antilinear leg).
- The lattice this all runs on: [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md) +
  `compute/MirrorWorld/Lattice.cs`.
