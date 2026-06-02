# Two-Term Palindrome: the Klein Routing of the Hidden Symmetry Q

**Status:** Computed (Tier 2); the routing rule is bit-exact across N=3,4,5. 2026-06-01.
**Authors:** Thomas Wicht, Claude (Opus 4.8)
**Answers:** [The Other Side](../hypotheses/THE_OTHER_SIDE.md) Q6 (the hidden symmetry Q) and Q7
(the second condition).
**Builds on:** [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) (the three Q-families),
[PTF under palindrome-breaking](PTF_PALINDROME_BREAKING_PERTURBATIONS.md) (the bit_a/bit_b axes).
**Scripts:** [`simulations/q6_klein_routing_two_term.py`](../simulations/q6_klein_routing_two_term.py)
(classification + routing rule + Q7), [`simulations/q6_hidden_q_construction.py`](../simulations/q6_hidden_q_construction.py)
(the explicit hidden Q, bit-exact). Reuses `fw.classify_pauli_pair` (the F87 trichotomy); adds the
primitive `fw.classify_two_term_palindrome`.

---

## The question

[The Other Side](../hypotheses/THE_OTHER_SIDE.md) left a sharp loose end. Of the 36 two-term spin
Hamiltonians it surveyed, 26 break the X-parity, and 14 of those break the eigenvalue palindrome.
That leaves 12 that break the parity yet keep the mirror: they must carry a hidden symmetry Q that
is not the usual mirror operator Π, and the question was, what is it, and why do exactly these 12
survive while 14 do not? The three Q-families were already known from
[Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md): a uniform per-site map, an alternating
odd/even map, and a non-local entangled operator. What was never built was the explicit map from each
surviving Hamiltonian onto its family, and the structural reason for the 14 that die. This note builds
both, and finds that the answer is written in the same two-bit alphabet that runs through the rest of
the framework.

## The setup

The 36 combinations are exactly the C(9,2) unordered pairs of the nine traceless bond bilinears
{XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ}. Each pair {t1, t2} becomes a Hamiltonian H = Σ_bonds (t1 + t2)
on an open chain under uniform Z-dephasing. The fate of each is read by the existing trichotomy
classifier (`fw.classify_pauli_pair`, which is exactly the F87 reading): **truly** if the canonical Π
already pairs the spectrum, **soft** if a hidden Q ≠ Π does, **hard** if nothing does. Each bilinear
also carries its Klein-Vierergruppe index (`fw.klein_index`): with M = (0,0) = {XX,YY,ZZ},
F_a = (0,1) = {XY,YX}, C = (1,0) = {YZ,ZY}, F_b = (1,1) = {XZ,ZX}.

The classification is N-invariant: the entire fate table, the Klein indices, and the parity flags are
bit-identical at N = 3, 4, 5. N=3 is not a small-N accident; it is representative.

| fate | count | meaning |
|------|:---:|---|
| truly | 3 | canonical Π (P1) pairs the spectrum |
| soft | 19 | a hidden Q ≠ Π pairs it (the benign ones) |
| hard | 14 | no Q; the palindrome is broken |

(26 break X-parity; the 14 hard are a strict subset, so 12 of the parity-breakers are benign, matching
The Other Side. The other 7 soft are parity-preserving.)

## Q6 answered: the explicit hidden Q

Constructing the actual conjugation operator for every palindromic combination and verifying
Q·L·Q⁻¹ = −L − 2Σγ·I to ‖·‖ ≤ 10⁻¹¹ (bit-exact), the 19 soft split cleanly into the three families,
identically at N = 3, 4, 5:

| family | count | combinations |
|--------|:---:|---|
| uniform (a single per-site Q, not the canonical P1) | 14 | XX+XZ, XX+YZ, XX+ZX, XX+ZY, XZ+YY, XZ+ZX, XZ+ZZ, YY+YZ, YY+ZX, YY+ZY, YZ+ZY, YZ+ZZ, ZX+ZZ, ZY+ZZ |
| alternating (odd/even site map) | 3 | XY+YX, XY+ZZ, YX+ZZ |
| non-local (entangled Q) | 2 | XZ+YZ, ZX+ZY |
| (truly, via canonical P1) | 3 | XX+YY, XX+ZZ, YY+ZZ |

([Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md) counts 17 uniform; those 17 are these 14
plus the 3 truly, which also close via a uniform per-site Q, the canonical P1, so the counts agree.)

So the 12 benign parity-breakers, and the 7 parity-preserving soft, all have a named, explicitly
constructed Q. The two non-local cases are exactly XZ+YZ and ZX+ZY, the pairs that force an X and a Y
onto the same qubit over a shared Z, the place where no per-site map can work and only an entangled
mirror closes the spectrum.

## The Klein routing (the viewpoint)

The fate is not a property the Hamiltonian wears on the outside; it is decided by the two bilinears'
letters alone, with no eigensolve, through how each routes the dephasing-damped channels. Z-dephasing
splits every site into a settled half {I, Z} and a humming half {X, Y}, the framework's near and far
banks. The mirror operator routes the humming half, and per site it comes in two crossovers: one
routes the X-channel, the other the Y-channel. A single uniform mirror has to serve both bond terms at
once, so the chain is palindromic exactly when its two terms admit a common crossover. Reading off each
term's admissible set,

    XX, YY  →  {X-router, Y-router}      ZZ  →  every valid router (it is dephasing-aligned)
    YZ, ZY  →  {X-router} only           XZ, ZX  →  {Y-router} only
    XY, YX  →  { }  (no single uniform router serves a bare X·Y)

the whole table falls out of one decision: the two terms share a router (truly if both are diagonal
Mother terms, soft-uniform otherwise); failing that, an XY/YX term with an XY/YX or ZZ partner closes
via the alternating router; failing that, the two same-site collisions XZ+YZ and ZX+ZY close
non-locally; failing that, it is hard. This reproduces the trichotomy classifier with zero
exceptions over all 36 combinations at N = 3, 4, 5.

Read in the framework's two bits, the active coordinate is clear, though neither bit alone is the
bit-exact rule. **Truly is the all-Mother corner**, both terms in Klein cell (0,0) (the F38/F88a
palindrome leg, where the mirror is already exact); it needs both bits zero on both terms, not bit_b
by itself. Among the rest, **bit_a, the light / dephasing-axis content, is the axis that turns**,
but the exact soft-versus-hard split is the cell-pair rule below, not bit_a alone. A pure
SU(2)-symmetric reading of the Klein classes is blind to that axis, which is why the class pair alone
fails: XY+YY and XY+ZZ are both Y-Father + Mother, yet XY+YY is hard and XY+ZZ is soft. The
difference is whether the Mother term is lit (YY carries a Y) or dark (ZZ is pure Z): a lit Mother adds
a conflicting routing demand on the qubit already carrying the Father's channel; a dark Mother is
dephasing-aligned and adds none.

## Q7 answered: the second condition

Breaking the X-parity is necessary (a term must be Π²-odd) but not sufficient. The second condition is
an **irreducible, unroutable same-qubit X/Y demand**: the palindrome breaks exactly when the two terms
force both an X-channel and a Y-channel onto one qubit in a way no router can reconcile and no entangled
mirror can rescue. Through the Klein cells the discriminator is a cell-pair rule with one refinement per
split cell: the two-Father cells {F_a, F_b} and the Child-plus-Y-Father cell {C, F_a} are always hard;
{F_a, M} is hard only when the Mother is lit; {C, F_b} is hard only when the two lit letters sit on
different sites. That flags exactly the 14, bit-exact at N = 3, 4, 5.

The entangled-mirror escape (the 2 non-local cases) is the boundary of the second condition: when the
same-qubit X/Y collision sits over a shared dark Z port, the conflict is clean enough that a non-product
Q absorbs it; when the conflicting demands are spread across both sites, nothing closes it. And the
frustration is genuinely many-body: at N = 2, with a single bond, all 36 are palindromic; the hardness
appears only at N ≥ 3, when the second bond's mirror must agree with the first on the shared middle
qubit. This is the V-Effect in its bare algebraic form, the coupling of two locally palindromic bonds
breaking the shared mirror.

## The routing rotates with the dephasing axis

Everything above is Z-dephasing, where the lit (damped) axis is {X, Y} and the dark axis is {I, Z}.
The three single-axis dephasing letters are related by the framework's Klein four-group of dephase-swaps
([`Pi2KleinV4DephaseSwapGroup`](../compute/RCPsiSquared.Core/Symmetry/Pi2KleinV4DephaseSwapGroup.cs),
proved bit-exact at universal N in
[PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE](../docs/proofs/PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md)),
so the routing should rotate rather than being a Z-accident. It does, exactly. Classifying all 36
combinations under X- and Y-dephasing and relabelling the bilinears by the matching swap, the
X-dephasing fate table is the Z table relabelled X↔Z and the Y table is the Z table relabelled Y↔Z,
with zero mismatches at N=3 and N=4; the counts {3, 19, 14} are dephasing-invariant; and the Z-only
routing rule reproduces every X- and Y-dephasing fate once the combo is relabelled. So the hidden-Q
routing rides on whichever axis the dephasing makes lit, not on Z in particular. The same Hamiltonian
can change fate with the noise axis: XY+YY is hard under Z-dephasing but soft under X-dephasing, its
lit and dark channels having rotated.

The two legs are not on the same footing, and the framework already proved why. The X↔Z swap is the
Q_zx Hadamard duality ([PROOF_BIT_A_TWIN_VIA_HADAMARD](../docs/proofs/PROOF_BIT_A_TWIN_VIA_HADAMARD.md)):
Q_zx is the operator-space lift of the physical Hadamard U_H^⊗N, a genuine Hilbert-space unitary that
carries the Z-dephased Lindbladian to the X-dephased one, so the X table is a true transport of the Z
table. The Y↔Z swap has no such Hilbert-space lift, since no unitary sends Y → −Y while fixing X and Z;
it is the operator-space D-swap at the palindrome-operator level, and the Y fate table matches the
relabelled Z table for the Route-1 reason (re-running the derivation with the lit axis set to {Y, Z}),
not by a Y↔Z Hadamard. The empirical zero-mismatch Y table is consistent with this; its mechanism is
just different from the X leg's.

This is the explicit k=2 face of the dissipator-resonance law
([THE_POLARITY_LAYER](../hypotheses/THE_POLARITY_LAYER.md), typed as
[`DissipatorResonanceLaw`](../compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs)), which
states the same axis-equivalence at k ≥ 3: F87-hardness lives in the Klein cell matching the dephasing
letter's index (Z → (0,1), X → (1,0), Y → (1,1)). One disambiguation: the SU(2) here is the discrete
Clifford / Klein-V₄ acting on the dephasing axis, not the Heisenberg chain's spin-rotation SU(2) (the
total-S² Casimir), which Z-dephasing breaks and which the Star / Schur-Weyl results use only as a tool
on H alone (see the clarifying note in [STAR_CONFOCAL_LIMIT](STAR_CONFOCAL_LIMIT.md)). Script:
[`simulations/q6_dephase_axis_rotation.py`](../simulations/q6_dephase_axis_rotation.py).

## What this is, and the tie to the rest

The routing rule is a viewpoint rotation on structure the framework already holds: the trichotomy
([F87](../docs/ANALYTICAL_FORMULAS.md)), the Klein cells ([F88a](../docs/ANALYTICAL_FORMULAS.md)), and
the Π²-even-always-palindromic result ([F108](../docs/ANALYTICAL_FORMULAS.md)). The mathematics falls
out exactly, which is the signature of a viewpoint, not new content; the contribution is the explicit
Q per combination (closing Q6), the second condition (closing Q7), and the bit-level reading that ties
them together.

That reading rhymes with the [PTF result](PTF_PALINDROME_BREAKING_PERTURBATIONS.md) without collapsing
into it. There, the per-site closure law was found to ride on bit_a (the light / U(1) axis), not on
the spectral palindrome (bit_b). Here too the load-bearing discriminator for the hard cases is bit_a,
while the truly line sits on bit_b. What is genuinely shared is the axis: under Z-dephasing the active
structural coordinate is bit_a. What is specific, and must not be merged, is the theorem: PTF's bit_a
governs a dynamical trajectory closure, this bit_a governs whether a spectral symmetry operator exists
at all. The same axis carries two different statements, a parallel, not one unified law.

No new typed claim is warranted. The two axes are already typed (F87/F88a/F108 and the trichotomy
classifier), and the result is a reading of them plus an explicit construction; like the PTF finding it
sharpens a Tier-2/3 hypothesis rather than anchoring a Tier-1 statement. The reusable piece is the
cockpit primitive `fw.classify_two_term_palindrome`, which returns the fate, the Q-family, and the
structural reason for any two-term bilinear pair.

## Cross-references

- [The Other Side](../hypotheses/THE_OTHER_SIDE.md): Q6/Q7, the X-parity vs palindrome containment.
- [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md): the three Q-families and the all-36 resolution.
- [PTF under palindrome-breaking](PTF_PALINDROME_BREAKING_PERTURBATIONS.md): the bit_a/bit_b dissociation.
- [F87](../docs/ANALYTICAL_FORMULAS.md), [F88a](../docs/ANALYTICAL_FORMULAS.md),
  [F108](../docs/ANALYTICAL_FORMULAS.md), [F38](../docs/ANALYTICAL_FORMULAS.md): the typed structure the rule reads.
- `fw.classify_two_term_palindrome` (`simulations/framework/diagnostics/q_family_routing.py`), and the
  two probes above.
