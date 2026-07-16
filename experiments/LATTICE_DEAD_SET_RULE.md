# The dead-set rule: three conserved structures close the h thread's open question

**Date:** 2026-07-16 (the h thread picked back up the same day; the mod-4 shadow derived the same session)
**Status:** verified from below. Gate: [`simulations/lattice_dead_set_rule.py`](../simulations/lattice_dead_set_rule.py) (21 checks, all PASS, exit 0). Minted as [F132](../docs/ANALYTICAL_FORMULAS.md) (2026-07-16; necessity Tier 1 derived, sufficiency Tier 2 gated).
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
sufficiency (everything not forbidden is alive) is the gated observation. The three
layers are the discovery path; the section on the mod-4 shadow below collapses them
into a single fermionic criterion that also covers the one prep family the layered
form misses.

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
(The small-zz lift is since quantified: the set jumps at any zz ≠ 0, the magnitudes
cross over as zz^m with m = the halved Majorana-degree gap:
[LATTICE_DEAD_SET_ZZ_FACE](LATTICE_DEAD_SET_ZZ_FACE.md).)

## The mod-4 shadow, derived: the rule collapses to one line

The first landing of this doc noted that V "appears to be the mod-4 shadow" of the
exact degree. Chasing that sentence gave the exact identity, and the identity
collapsed the whole rule.

**The identity (N-free, every Pauli string, both gauges).** V's kill sign is a pure
function of the Majorana degree:

    ε_odd-gauge(O) = (−1)^(d(d−1)/2),   ε_even-gauge(O) = (−1)^(d(d+1)/2),

i.e. the odd-sites gauge kills exactly d ≡ 2, 3 (mod 4) and the even-sites gauge
kills exactly d ≡ 1, 2 (mod 4). Derivation, one paragraph: X^N sends a_2l ↦ (−1)^l
a_2l and a_2l+1 ↦ (−1)^(l+1) a_2l+1 (the tail collects one sign per site below);
conjugation fixes the real a_2l and flips the imaginary a_2l+1; U_g flips both
Majoranas of every site in g. The combined per-SITE sign is (−1)^(l + [l∈g]),
uniformly −1 for the even-sites gauge and uniformly +1 for the odd-sites gauge; and
a Hermitian degree-d string is ±i^(d(d−1)/2) times an ordered Majorana monomial,
whose phase the antilinear map conjugates to the sign (−1)^(d(d−1)/2). Gated
three ways at once (letter formula == degree formula == direct matrix conjugation)
over all 4^N strings × both gauges, N = 2..6: 10912 pairs, zero mismatches.

**The per-sector refinement.** In this language the stabilizer bookkeeping dissolves:
V_g is ALWAYS a symmetry of the flow; what varies is how the prep decomposes. The
population part and the coherence part live in disjoint, flow-invariant block
classes, and each is separately a V_g eigenvector: population with sign +1 (both
gauges), coherence with sign (−1)^|g|. Each sector's contribution to ⟨O⟩ dies unless
ε_g(O) matches the sector's sign for BOTH gauges. For the population sector that
demands ε_even = ε_odd = +1, which by the identity is exactly d ≡ 0 (mod 4). For the
coherence sector at d = N the required signs hold AUTOMATICALLY: the identity gives
ε_odd = (−1)^(N(N−1)/2) = (−1)^⌊N/2⌋ and ε_even = (−1)^(N(N+1)/2) = (−1)^⌈N/2⌉, and
⌊N/2⌋, ⌈N/2⌉ are exactly the sizes of the two sublattices, so ε_g = (−1)^|g| for
both gauges at every N. The coherence channel is never V-killed. The whole rule collapses to one line:

    alive  ⟺  (mask connects a populated diagonal block ∧ d ≡ 0 mod 4)
              ∨ (coherence on ∧ mask connects the coherence blocks ∧ d = N).

The same honesty split as the layered form: the identity sharpened only the
NECESSITY direction (forbidden ⟹ dead is derived); the ⟸ direction (everything
allowed is alive) remains the gated observation. (Since gated GENERICALLY, six
random h directions per config, and refined at the axis's special points:
[LATTICE_DEAD_SET_H_ZERO](LATTICE_DEAD_SET_H_ZERO.md).)

**Where the collapsed form is load-bearing.** The layered form's global-stabilizer
bookkeeping and the per-sector form agree on the whole battery above, and they
DIVERGE at N ≡ 2 (mod 4) with coherence: there both sublattices have odd size, no
V_g fixes ρ(0) globally, and the layered form imposes no V constraint at all. Gated
at N = 6 (s = 1, real coherence, one-sided watching, full 4095-census; the necessity
arguments are watching-profile-independent and the full-profile face is gated at
N = 5): the global form over-predicts 2047 alive; the actual set is 1055; the
collapsed form is exact. The alive-by-degree
counts come out as binomial coefficients cut by kinematics, e.g. at N = 6 coherence
{d=4: 255, d=6: 544, d=8: 255, d=12: 1} = {C(12,4) − 240, C(12,6) − 380, C(12,8) −
240, C(12,12)}, and at N = 5 s = 3 coherence the full degree-5 sector C(10,5) = 252
is alive. (The d = 12 survivor at N = 6 is Γ = Z^⊗N itself, constant and nonzero at
even N.)

**What remains of the layer picture.** K is untouched (and is the only piece that
survives interactions); F's exact-degree conservation now carries V's mod-4 kill as
its coarser shadow, the two meeting in one criterion; on the odd-w sector
ε_even·ε_odd = (−1)^w = −1 explains why at most one gauge condition ever looked
independent there. The three-layer story above stays as the honest discovery path.

## Files

- Gate: [`simulations/lattice_dead_set_rule.py`](../simulations/lattice_dead_set_rule.py)
  (T1 ten full censuses, T2 per-layer witnesses, T3 the mod-4 identity over all
  strings N = 2..6, T4 the collapsed per-sector rule incl. the N = 6 divergence and
  the binomial anatomy, T5 the zz boundary).
- The question this closes: [LATTICE_H_THREAD](LATTICE_H_THREAD.md) §4.
- The h axis walked: [LATTICE_DEAD_SET_H_ZERO](LATTICE_DEAD_SET_H_ZERO.md)
  (sufficiency generic; at h = 0 the bi-degree refinement, one level below
  this rule's mod-4 shadow).
- The two composed mirrors: [PROOF_MIRROR_ORDER_SORTING](../docs/proofs/PROOF_MIRROR_ORDER_SORTING.md)
  (F131, X^N as its third sighting), [PROOF_K_PARTNERSHIP](../docs/proofs/PROOF_K_PARTNERSHIP.md)
  (the chiral gauge), [PROOF_ANTILINEAR_TRIANGLE](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md)
  (conjugation as an antilinear leg).
- The lattice this all runs on: [LATTICE_OPENING_LAW](LATTICE_OPENING_LAW.md) +
  `compute/MirrorWorld/Lattice.cs`.
