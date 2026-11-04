# The sterile↔birth-canal boundary and the odd↔even junction are one object

**Status:** Tier 1 candidate (gate-verified, N = 5, 6 numerical; the (0,1) reduction it rests on is
Tier 1 derived). The sterile↔birth-canal boundary and the coherence-horizon junction are two readings of
one quantity, `rate_slow(Q)` decomposed by the change-number Δn. The open chain's per-block floor
(the section on the Δn = 1 blocks) rests on PROOF_CHAIN_GAP_DOMINANCE and is exact block by block at
N = 2 … 6 (`FreeFermionFloorWitness`).
**Date:** 2026-06-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Verifier:** [`birth_canal_junction_nature.py`](../simulations/birth_canal_junction_nature.py)
(self-validating, gate-first).
**Builds on:** [`VacuumBlockReductionClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/VacuumBlockReductionClaim.cs)
(the (0,1) reduction), `HandoverFloorClaim` (the junction = the handover Q_h, = Q*(N) at N=2,3 only),
[`PROOF_STRUCTURAL_CEILING.md`](proofs/PROOF_STRUCTURAL_CEILING.md) §7 (the star ceiling); the Absorption
floor is `AbsorptionTheoremClaim`.

<!-- STAGE 0 SWEEP RECORD (2026-11-03, for the section on the free-fermion floor of the whole
Δn = 1 family). Question: what does the repo hold on the floor of the (n, n+1) blocks for n ≥ 1, and
on its relation to F86. docs/proofs/: PROOF_CHAIN_GAP_DOMINANCE §2 and §3 own the floor (every mode
at Re λ = −2γ is c_k^(†)·f(N_tot) from N = 4, with Δn = 0 extras at N = 2, 3; completeness a
dimension match at N = 4, 5, 6, 2N² modes in all);
PROOF_WEIGHT1_DEGENERACY is the Heisenberg real-mode analogue; PROOF_F86A_EP_MECHANISM and
PROOF_F86B_OBSTRUCTION use the HD = 1 channel rate −2γ₀ but never the floor family.
docs/ANALYTICAL_FORMULAS.md: the F2b corollary cites the floor modes, nothing per block; F153 pins
whole blocks only (min(p, q) = 0 or max(p, q) = N); F50 is Heisenberg. Typed layer:
VacuumBlockReductionClaim covers the (0,1) block only; PinnedBlockFloorClaim (F153);
F89PathKAtLockMechanismClaim counts floor(N/2) F_a modes at rate 2γ in the (1,2) block;
CoherenceHorizonWitness says in prose that the band edge repeats once in each of the 2N Δn = ±1
sectors; ClockHandLadderClaim types the proof. experiments/: VEFFECT_CAVITY_MODES (Heisenberg),
THE_FLOW_BETWEEN_TWO_SINGULARITIES (the flow EP, kept apart from F86 on purpose), Q_SCALE_THREE_BANDS
nothing. docs/GLOSSARY.md: the 2γ rung exists at every coupling at uniform γ. OpenArcs:
birth_canal_horizon_junction, birth_canal_surface and site_resolved_vacuum_block, none on the interior
blocks. docs/CAUGHT_ERRORS.md: CHAIN_GAP_DOMINANCE audited clean. fw.Confirmations: nothing.
Adjacency: four results about this one object that did not cite each other, the proof's floor family,
F89's F_a set (the floor modes even under the chain reflection of the coherence, checked by
FreeFermionFloorWitness), CoherenceHorizonWitness's 2N
repetition, and this document's sterile reading; and no birth-canal document named F86, no F86
document named the birth canal; docs/NAVIGATING_THE_DIMENSIONS.md names F86a's proof and the flow's
birth canal a few paragraphs apart and keeps them apart. New here is only the per-block reading and the bridge. -->

## What this is about

The slowest-decaying (longest-lived) coherence of a dephased spin network has been studied through two
boundaries that looked separate: the *sterile↔birth-canal* boundary (is that survivor's decay rate frozen
as you change the dephasing, or does it move?) and the *coherence-horizon junction* (does a second,
number-conserving mode overtake the first as the survivor?). They are two readings of one quantity. Sort
the coherences by Δn, how much each changes the excitation number: the number-changing modes (Δn = 1) sit
pinned at the Absorption floor, the number-conserving interior (Δn = 0) drifts with the dephasing. The
junction is the moment the interior overtakes the band edge; the birth canal is any departure from the
frozen, "sterile" rate. The junction always lands you in the birth canal, but the birth canal has a
second entrance the junction does not use.

## Abstract

Under XY hopping with Z-dephasing, the slowest non-kernel Liouvillian rate as a function of `Q = J/γ` is a
minimum over the joint-popcount sectors (p_col, p_row), organised by the change-number `Δn = |p_col −
p_row|`. The Δn = 1 number-changing band edge (the (0,1) |1-exc⟩⟨vac| class) sits at the −2γ Absorption
floor: Q-flat at uniform γ (its dissipator is scalar there), Q-drifting at non-uniform γ. The Δn = 0
number-conserving interior (the diagonal (p,p), the {0,2}-coherence) sits at −2γ⟨n_XY⟩(Q), brightening
toward the floor as Q rises. Two boundaries are then facets of this one object: the sterile↔birth-canal
boundary (`BirthCanalSurfaceWitness`) reads whether `rate_slow(Q)` is Q-flat; the odd↔even junction (the
handover Q_h, = Q*(N) at N=2,3 only) reads whether the Δn = 0 interior overtakes the Δn = 1 band edge. The seam is
**junction ⟹ birth canal, but not conversely**: the canal is also entered by the Δn = 1 survivor itself
drifting under a non-uniform profile (odd-drift, no Δn switch), so the junction is a strict sub-mechanism.
When the Δn = 0 interior is the survivor, its nature is topology-set: ring (2,2) is a frozen level
crossing (robust to the γ-profile), chain (2,2) is the oscillating filling-degenerate SE-EP (rigidity → 0
at Q*), and star (1,1) is the frozen commutant survivor (the structural ceiling g2 = 4/(N−1) ≤ 1 read
dynamically). Gate-verified N = 5, 6; Tier 1 candidate.

## The object: rate_slow(Q) = min over joint-popcount sectors

Under XY hopping + Z-dephasing the slowest non-kernel rate, as a function of the dephasing knob
Q = J/γ, is a minimum over the joint-popcount sectors (p_col, p_row). Organise them by the
change-number **Δn = |p_col − p_row|** (the n_diff parity of `JointPopcountSectors`):

- **Δn = 1, the number-changing band edge** (the (0,1) |1-excitation⟩⟨vacuum| class). Its rate is the
  −2γ Absorption floor (`AbsorptionTheoremClaim`). At **uniform γ** the block −iQh − 2γ·I is scalar in
  its dissipator, so Re = −2γ exactly, **Q-flat**. At **non-uniform γ** the hopping h mixes sites of
  different γ, so the slowest Δn=1 mode's rate **drifts with Q**.
- **Δn = 0, the number-conserving interior** (the diagonal (p,p), the {0,2}-coherence). Its rate is
  −2γ⟨n_XY⟩(Q), **Q-dependent** (it brightens toward the −2γ floor as Q rises).

## On the open chain the floor is every Δn = 1 block's, and F86 works in those blocks

The (0,1) block is the smallest member of the Δn = 1 family, not the only one: every (n, n+1) block
is number-changing, and the interior ones, chromaticity c ≥ 2, are the blocks F86 works in. On the
open XY chain at uniform γ every one of them holds the same kind of floor. Inside these blocks the
modes at Re λ = −2γ are the Jordan-Wigner family c_k^(†)·f(N_tot) of
[PROOF_CHAIN_GAP_DOMINANCE](proofs/PROOF_CHAIN_GAP_DOMINANCE.md) §2 and §3: in the (n, n+1) block that
is P_n·c_k·P_(n+1) = c_k·P_(n+1), k = 1 … N, with its adjoint in the transposed block, oscillating at
the single-particle comb ±2J·cos(πk/(N+1)) whatever the chromaticity (with H = (J/2)·Σ(XX+YY),
so one excitation hops with amplitude J). Each such mode carries exactly
one X or Y letter, so the dephasing charges it the single price 2γ. A single letter in general does
not stay single under the hopping, which can turn a background Z beside an I into an X/Y pair, three
letters in all (the proof's §1); the Jordan-Wigner string, Z on every site to one side of the letter and
I on every site to the other, has no such bond, the hopping carries it to strings of the same shape, and the
number projector that dresses it in the block commutes with the hopping, so on this family the hopping
only moves the letter along, and Q turns the phase and never the rate. The family gives N
floor modes in every block at every N; that there are no others is exact, an integer lower bound
meeting a mod-p upper bound, at N = 2 … 6 (`FreeFermionFloorWitness`, `inspect --root fermionfloor`). F86's K at the (1,2) block splits exactly by
this family into floor², rest² and their interference (`FloorSplitWitness`, `inspect --root floorsplit`;
the plain-words reading is [ON_WHAT_CLOSES_ALONE](../reflections/ON_WHAT_CLOSES_ALONE.md)).

The fence. A bond profile leaves the family in place, the comb then being the spectrum of the
weighted hopping matrix (the proof's scope). So does a γ profile, in its form: P_n·c_j·P_(n+1)
disagrees at site j alone, so the dephasing charges it 2γ_j, and every block carries the vacuum
block's own operator, so every block holds the (0,1) block's eigenvalues: the (n, n+1) family is the spectator
map W of [PROOF_CODIM1_BY_ADDITIVITY](proofs/PROOF_CODIM1_BY_ADDITIVITY.md) Theorem B applied n times to the
(0,1) block's coherences, and on the open chain W commutes with L for any bond profile and any γ profile (not on the ring, whose wrap
bond carries (−1)^N_tot: the ring's H is not quadratic in the fermions, W flips the parity on both sides, and the
two sides of a Δn = 1 coherence do not even see one and the same wrap bond). Sides whose particle numbers share a parity do see one, and on those blocks of the ring W, run twice, still carries each block two up ([PROOF_CODIM1_BY_ADDITIVITY](proofs/PROOF_CODIM1_BY_ADDITIVITY.md) §5; for the fillings, [PROOF_WEIGHT1_DEGENERACY](proofs/PROOF_WEIGHT1_DEGENERACY.md#the-zeno-end)). So the wrap bond's sign, being set by a parity, keeps W² on the ring's blocks of even height and takes this family's construction away from its interior Δn = 1 blocks, where the floor is indeed empty (below). What a profile takes away is the single rate: the family's rates then move with Q,
which is this document's odd-drift entrance to the birth canal, taken in every block at once. A ZZ
term breaks the interior blocks' N-member family: at the Heisenberg point each interior block keeps a
single floor mode (read at N = 3 … 7), at frequency 0, the real modes of
[PROOF_WEIGHT1_DEGENERACY](proofs/PROOF_WEIGHT1_DEGENERACY.md). So does the ring in
its interior: there the Δn = 1 floor is the 4N modes of the four outer blocks alone at N = 4 … 8
([PROOF_RING_GAP_DOMINANCE](proofs/PROOF_RING_GAP_DOMINANCE.md); at N = 3 the (1,2) block holds more),
and the ring's interior blocks, read at N = 5 and 6 and Q = 1, 2 and 5, decay faster than 2γ at a rate
that moves with Q (at N = 5 their slowest rate goes 2.16, 2.61, 3.40 γ).

So on the uniform open chain, wherever the floor is the slowest layer, the sterile reading reads a
layer that is Q-flat in every block by construction. It does not follow that Q leaves the channel alone. The family's span is invariant under
L and under its adjoint alike, so L splits into the family's part, which is normal, and an orthogonal
rest: no exceptional point of a block involves a mode of the family, at any N and any Q. Where the count is exact (N = 2 … 6), every other mode of a block decays strictly
faster than 2γ at every Q > 0. The uniform Dicke probe is silent, exactly, on the floor modes that are odd under the chain
reflection of the coherence, site j ↔ N−1−j applied to both of its sides, and reads the even ones as
nonzero in every block checked. That parity is φ_k's, (−1)^(k+1), times the sign (−1)^n the
Jordan-Wigner string picks up, so the visible k alternate from block to block; in the (1,2) block
they are the even k, the k-labels of F89's F_a set, which F89 names by φ_k's own parity, the
S₂-antisymmetric orbit. A slowest-rate reading cannot tell any of this apart: wherever the floor is
the slowest layer, the reading stays flat at uniform γ, whatever Q does above it.

## The two facets

- **sterile ↔ birth-canal** (`PostEpFlowField.BirthCanalDeviation` / `IsInBirthCanal`,
  `BirthCanalSurfaceWitness`, `inspect --root surface`): reads whether rate_slow(Q) is Q-flat. Sterile
  = Q-independent closed form (frozen, "no creation"); birth canal = Q-modulated.
- **odd ↔ even junction** (the handover Q_h, = Q*(N) at N=2,3 and just below it from N=4; arc `birth_canal_horizon_junction`): reads whether the
  Δn=0 interior overtakes the Δn=1 band edge as the survivor.

## The seam

**Junction ⟹ birth canal, but not conversely.** The birth canal (Deviation ≠ 0) is entered by two
distinct mechanisms:

| mechanism | what changes | example (verifier) |
|---|---|---|
| **odd-drift** | the same Δn=1 survivor drifts (non-uniform γ; no Δn switch) | N=5 canal [.25,1.5,1.5,1.5,.25], Deviation 0.085 |
| **junction** | Δn flips 0→1: the interior overtakes the band edge (= the handover Q_h, just below Q*(N) from N=4) | N=6 deep-edge, Deviation 0.408, Δn 0→1 |

Sterile = the Δn=1 band edge reigns and is Q-flat (the −2γ floor; e.g. uniform N=5, where the (0,1) band
edge stays the global slowest across the whole probe window because at uniform γ its rate is Q-invariant,
the flat-γ blindness of `VacuumBlockReductionClaim`). So the junction is a **strict sub-mechanism** of the
birth canal: it is the
*survivor-identity-change* way of leaving the sterile zone, while odd-drift is the
*same-survivor-drifts* way. `PostEpFlowField` sees the union; the junction arc sees only the junction.

## Which interior survivor sits on the birth-canal side

When the Δn=0 interior is the survivor, its nature is topology-set, not γ-profile-set:

- **Ring (2,2): a frozen level crossing** (|Im| ≈ 1e-15, a non-rotating V-Effect seam). Tracked by
  continuation, it **stays frozen** as the γ-profile turns on (uniform → deep-edge): the freezing is
  intrinsic to the wrap-bond two-fermion structure, not a uniform-γ symmetry artifact.
- **Chain (2,2): oscillating.** At uniform γ it is filling-degenerate with (1,1), i.e. the
  single-excitation √-EP coalescence (rigidity → 0 at Q*(N), the coherence horizon); it does **not**
  inherit the ring's frozen character. The EP-vs-crossing distinction is read by the defectiveness
  test (`EpCharacterWitness`, `inspect --root epcharacter`: geometric multiplicity 1 below algebraic
  2), with the rigidity's square-root law (`CoherenceHorizonWitness`) beside it, not by |Im| alone (the SE-EP is
  overdamped-real below Q*, so a single-Q |Im| is blind).
- **Star (1,1): frozen, the commutant route.** Its boundary survivor never un-freezes (|Im| ≈ 0 at all Q,
  N ≥ 5), because it is the darkest [H,A]=0 commutant coherence (commutes with H ⟹ no oscillation) and it
  undercuts the −2γ floor exactly when the structural ceiling g2 = 4/(N−1) ≤ 1; N=4 (g2 = 4/3 > 1) is the
  outlier that un-freezes. So the star's frozen seam IS the ceiling read dynamically, a different route to
  |Im|=0 than the ring's level crossing. See `docs/THE_STAR_FROZEN_SEAM.md`.

## See also

- Verifier: `simulations/birth_canal_junction_nature.py` (control + Stage 0 nature + Stage 1
  freezing-robustness continuation + Stage 2 the seam).
- Sibling cluster docs: `docs/THE_STAR_FROZEN_SEAM.md` (the star's commutant route to |Im|=0),
  `experiments/THE_HUB_KILLS_THE_HORIZON.md` (a dominant hub removes the horizon however dispersive the
  leaves), `docs/proofs/PROOF_STRUCTURAL_CEILING.md` §7 (the g2 = 4/(N−1) ceiling) and
  `docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md` (the chain's horizon Q*(N)); and the unified view
  [`docs/THE_TRICHOTOMY_SEEN.md`](THE_TRICHOTOMY_SEEN.md): this Δn-seam alongside the chain/ring/star
  un-freeze, as one figure + the live `inspect --root trichotomy` tree.
- Claims: `HandoverFloorClaim`, `CoherenceHorizonClaim`, `VacuumBlockReductionClaim`,
  `SecondClockRegimeClaim`; the floor is `AbsorptionTheoremClaim`, its free-fermion family
  `ClockHandLadderClaim` ([PROOF_CHAIN_GAP_DOMINANCE](proofs/PROOF_CHAIN_GAP_DOMINANCE.md)), its
  Dicke-visible modes in the (1,2) block `F89PathKAtLockMechanismClaim`.
- Code: `PostEpFlowField`, `BirthCanalSurfaceWitness` (`--root surface`), `SectorReductionWitness`
  (`--root reduction`), `IncompletenessSurvivorWitness` (`--root survivor`), `FreeFermionFloorWitness`
  (`--root fermionfloor`, the floor of every (n, n+1) block).
- F86, the same Δn = 1 blocks read by their response rather than their slowest rate:
  [PROOF_F86_QPEAK](proofs/PROOF_F86_QPEAK.md), [F86_VALUES_INVENTORY](F86_VALUES_INVENTORY.md).
- Prior ground truth: `simulations/carbon/handover_q.py` (the chain handover Q_h, = Q*(N) at N=2,3 only / ring = N√3/(2π), slope DERIVED in `docs/proofs/PROOF_RING_HANDOVER_SLOPE.md`),
  `simulations/birth_canal_n6_mode_crossing.py` (the N=6 deep-edge crossing).
- `docs/ANALYTICAL_FORMULAS.md` F2b corollary (the handover entry); `docs/HIERARCHY_OF_INCOMPLETENESS.md`
  (the C=0.5 / sterile reading).

The full sweep of the seam over (Q, profile, N, topology) is now the live witness `inspect --root trichotomy`
([`docs/THE_TRICHOTOMY_SEEN.md`](THE_TRICHOTOMY_SEEN.md)): the carbon un-freeze read and this absolute Δn-seam
read assembled into one object. The star is characterized (`docs/THE_STAR_FROZEN_SEAM.md`): its survivor is
frozen at all Q for N ≥ 5 (the commutant route), exactly the structural ceiling g2 = 4/(N−1) ≤ 1 read
dynamically; N=4 is the outlier.
