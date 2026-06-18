# The sterile↔birth-canal boundary and the odd↔even junction are one object

The repo carried two boundaries of the longest-lived (slowest non-kernel) Liouvillian mode as if they
were separate. They are two readings of a single object. Verifier (self-validating):
`simulations/birth_canal_junction_nature.py`.

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

Sterile = the Δn=1 band edge reigns and is Q-flat (the −2γ floor; e.g. uniform N=5, whose probe window
sits above Q*(5)). So the junction is a **strict sub-mechanism** of the birth canal: it is the
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
- Claims: `HandoverFloorClaim`, `CoherenceHorizonClaim`, `VacuumBlockReductionClaim`,
  `SecondClockRegimeClaim`; the floor is `AbsorptionTheoremClaim`.
- Code: `PostEpFlowField`, `BirthCanalSurfaceWitness` (`--root surface`), `SectorReductionWitness`
  (`--root reduction`), `IncompletenessSurvivorWitness` (`--root survivor`).
- Prior ground truth: `simulations/carbon/handover_q.py` (the chain = Q*(N) / ring = 0.29N handover),
  `simulations/birth_canal_n6_mode_crossing.py` (the N=6 deep-edge crossing).
- `docs/ANALYTICAL_FORMULAS.md` F2b corollary (the handover entry); `docs/HIERARCHY_OF_INCOMPLETENESS.md`
  (the C=0.5 / sterile reading).

Open: a full sweep of the seam over (Q, profile, N). The star is now characterized (`docs/THE_STAR_FROZEN_SEAM.md`):
its survivor is frozen at all Q for N ≥ 5 (the commutant route), exactly the structural ceiling g2 = 4/(N−1) ≤ 1
read dynamically; N=4 is the outlier.
