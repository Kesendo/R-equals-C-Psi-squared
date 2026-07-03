# The sterile↔birth-canal boundary and the odd↔even junction are one object

**Status:** Tier 1 candidate (gate-verified, N = 5, 6 numerical; the (0,1) reduction it rests on is
Tier 1 derived). The sterile↔birth-canal boundary and the coherence-horizon junction are two readings of
one quantity, `rate_slow(Q)` decomposed by the change-number Δn.
**Date:** 2026-06-17
**Authors:** Thomas Wicht, Claude (Anthropic, Opus 4.8)
**Verifier:** [`birth_canal_junction_nature.py`](../simulations/birth_canal_junction_nature.py)
(self-validating, gate-first).
**Builds on:** [`VacuumBlockReductionClaim`](../compute/RCPsiSquared.Diagnostics/Foundation/VacuumBlockReductionClaim.cs)
(the (0,1) reduction), `HandoverFloorClaim` (the junction = Q*(N)),
[`PROOF_STRUCTURAL_CEILING.md`](proofs/PROOF_STRUCTURAL_CEILING.md) §7 (the star ceiling); the Absorption
floor is `AbsorptionTheoremClaim`.

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
handover Q*(N)) reads whether the Δn = 0 interior overtakes the Δn = 1 band edge. The seam is
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

## The two facets

- **sterile ↔ birth-canal** (`PostEpFlowField.BirthCanalDeviation` / `IsInBirthCanal`,
  `BirthCanalSurfaceWitness`, `inspect --root surface`): reads whether rate_slow(Q) is Q-flat. Sterile
  = Q-independent closed form (frozen, "no creation"); birth canal = Q-modulated.
- **odd ↔ even junction** (the handover Q*(N); arc `birth_canal_horizon_junction`): reads whether the
  Δn=0 interior overtakes the Δn=1 band edge as the survivor.

## The seam

**Junction ⟹ birth canal, but not conversely.** The birth canal (Deviation ≠ 0) is entered by two
distinct mechanisms:

| mechanism | what changes | example (verifier) |
|---|---|---|
| **odd-drift** | the same Δn=1 survivor drifts (non-uniform γ; no Δn switch) | N=5 canal [.25,1.5,1.5,1.5,.25], Deviation 0.085 |
| **junction** | Δn flips 0→1: the interior overtakes the band edge (= the handover Q*(N)) | N=6 deep-edge, Deviation 0.408, Δn 0→1 |

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
  inherit the ring's frozen character. The EP-vs-crossing distinction is read by Petermann phase
  rigidity (`PhaseRigidity`, `CoherenceHorizonWitness`), not by |Im| alone (the SE-EP is
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
  `SecondClockRegimeClaim`; the floor is `AbsorptionTheoremClaim`.
- Code: `PostEpFlowField`, `BirthCanalSurfaceWitness` (`--root surface`), `SectorReductionWitness`
  (`--root reduction`), `IncompletenessSurvivorWitness` (`--root survivor`).
- Prior ground truth: `simulations/carbon/handover_q.py` (the chain = Q*(N) / ring = N√3/(2π) handover, slope DERIVED in `docs/proofs/PROOF_RING_HANDOVER_SLOPE.md`),
  `simulations/birth_canal_n6_mode_crossing.py` (the N=6 deep-edge crossing).
- `docs/ANALYTICAL_FORMULAS.md` F2b corollary (the handover entry); `docs/HIERARCHY_OF_INCOMPLETENESS.md`
  (the C=0.5 / sterile reading).

The full sweep of the seam over (Q, profile, N, topology) is now the live witness `inspect --root trichotomy`
([`docs/THE_TRICHOTOMY_SEEN.md`](THE_TRICHOTOMY_SEEN.md)): the carbon un-freeze read and this absolute Δn-seam
read assembled into one object. The star is characterized (`docs/THE_STAR_FROZEN_SEAM.md`): its survivor is
frozen at all Q for N ≥ 5 (the commutant route), exactly the structural ceiling g2 = 4/(N−1) ≤ 1 read
dynamically; N=4 is the outlier.
