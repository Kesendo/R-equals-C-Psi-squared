# The corner beat: a hardware prediction (pre-registration, DRAFT v7.4)

**Status: DRAFT v7.4, refreeze stage (2026-08-18). The runner
`run_corner_beat.py` exists in the external pipeline
(ibm_quantum_tomography, its own git repository), went through FOUR
empty review rounds (physics, spec-compliance, measurement statistics;
round 4: no physics or statistics blocker), and its local proofs pass:
certify at machine precision (Strang-step sector parity 7e-16,
trajectory parity vs the committed gate simulator 5e-14, fit parity
exactly 0.0), and the analyze chain proven end-to-end on
hardware-shaped synthetic artifacts at the frozen M = 1024. Submission
is structurally blocked by the runner's constants manifest until the
remaining gate work freezes every §8 quantity (θ_D itself carries
refreeze_required, §8a). Amendment 2 (2026-08-18) MEASURES θ_D, θ_W and
the s²(C) floor through the analyze-side chain and records four
systematics beside them; it does NOT freeze the manifest.

**PARKED BY TOM, 2026-08-19, before submission. This flight is not
scheduled and no money has been spent.** The reason is §10's own
arithmetic: the signal grows with depth at exponent γ̄T·(10/3) ≈ 2.95
while the state is lost at 970·p2·f_leak, which the §9 layered inflation
carries above the signal exponent at BOTH ends of the f_leak bracket. The
design has no free lever left (shots are measured saturating, more depth
costs the state, less depth costs discrimination), the parameter that
decides it cannot be narrowed before spending, and §1's own scope
sentence says the registered arms cannot distinguish the class law from
a dead C′ while the quantitative width law sits outside the registered
conjuncts. Paying for the weak statement when the closed form is what we
hold is the trade that was declined. The arc `corner_beat` in
`compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs` carries the
wake-up condition: an effective two-qubit error about two to three times
below today's, or a re-registration around the width law. **Everything
below stands as the design it was; nothing in it is retracted by the
parking**, and the closed form of §1 never depended on a device. Design and gate history: seven design rounds,
the gate's two review rounds, v2.1 runs in §8a. Verdict: provisionally
flyable, CONDITIONAL ON FRACTIONAL-RZZ, at the frozen 21-point grid
with 16384 shots, ~23.7 QPU min (worst-end power 2.9, P(detect) 0.92,
pre-refreeze). Tom's two budget decisions (2026-08-16) froze the
21 × 3 grid at 16384 shots: worst-end P(detect) = 0.920
at H0 false rate 0.000 held-out (0/50; ~6% binomial 95% bound, §8a);
both numbers MEASURED AT p2 = 0.5%/2q, which since Amendment 1.8
(Tom's decision, 2026-08-17) IS the Class-1 bound: the registered
operating characteristics are the bound's OWN, and the linear
interpolation that stood here while the bound was 0.6% is gone with
it, as is the measured point the refreeze owed at that bound. The
guard is an ISOLATED-gate number and the circuit runs three rzz
simultaneously per half-layer, so the day-of LAYER-FIDELITY gate (§9)
now carries the whole load that interpolation used to flag; these
numbers
become VERDICT-GRADE only at the ≥ 500-rep refreeze (§8a, the 0/50
caveat); the far bracket end
re-evaluates at the single registered θ_D in the refreeze (§7,
§8a).
Amendment 1 (below) records the analysis chain as implemented. Any
flight remains Tom's separate go.**

## What this is about

Take six spins in a row, all coupled equally, with exactly one spin flipped
up. That one flip can ride along the row in six different wave patterns,
each with its own pitch. Prepare two patterns at once and the row hums: the
site occupations beat at the difference of the two pitches. Three such
pattern pairs share the same difference, so the row offers three hums at
one frequency. Shine engineered noise on the row and the hums fade. This
experiment moves WHERE the noise shines without changing how much of it
there is in total. Spread evenly over all six sites, all three hums fade at
one single rate. Pushed to a corner, all of it on one site of each mirror
pair and none on the partner, the fading splits, and WHICH corner decides
by how much: half the corners give a wide split whose fastest rate sits at
the edge of what the proven containment allows, the other half a split
exactly 1/√3 as wide that stays clear of the edge. The machine is asked to
fly one corner of each kind against the even spread and show all three
numbers: the wide split, the narrow split, and the unchanged average.

## What the repo already holds

Swept 2026-08-16 (two store-sweep agents + the session's own reads), by
store:

- `experiments/THE_TWO_SPIN_ZEROS.md` (45e2f29) and
  `experiments/THE_MIRROR_TRANSVERSAL_CERTIFICATE.md` (e3dbab0): the theory
  this flight rides. The (1,1) ω = −2 room at N = 6, Δ = 1 is dim 3,
  breaking (no C_l protection on (1,1)), with gated w = 1/3, δ = 0,
  TIGHT at the UPPER edge of the centre interval
  (`simulations/results/locus_containment_gate.txt`, the
  `(1,1) om=-2 ... TIGHT(upper)` line). The gate prints
  `bound=0.3333` on the interval [0, 2] (the fraction 1/3 is our
  reading of that decimal, round-15 precision); with s_min = 0 and
  s_max = 2 from §1's fence, where this s is the certificate's SIZE
  coordinate and NOT the split statistic s of §1/§2 and §6.4, two
  different quantities sharing a letter (round-17 warning; w, the
  half-width, is defined ten lines below), bound = min(s_max − m̄,
  m̄ − s_min)
  admits m̄ = 1/3 and m̄ = 5/3, and TIGHT(upper) selects m̄ = 5/3
  (recomputed from below; m̄ is the mean of the compressed rates in
  units of **2γ̄**, per §1's own mapping rate = 2·s·γ̄, so m̄ = 5/3
  is the mean rate (10/3)γ̄ that §2 states. Round-16 repair of a
  round-15 clause that wrote "units of γ̄" and put a factor 2 on the
  only bridge between the gate's printed number and the flown rates;
  the check is m̄ ± w = 5/3 ± 1/3 = {4/3, 2} ↦ {8/3, 4}γ̄, exactly
  §2's extreme maximizing rates, and TIGHT(upper) at s = 2 ↦ 4γ̄).
  w = max_σ ‖Σ_l σ_l C_l‖ (l ranging over MIRROR PAIRS in this sum) is
  an ADVERSARIAL MAXIMUM over the transversal sign choices; §2's two
  transversal classes, the maximizer list, and the rate triples are
  an arithmetic consequence COMPUTED HERE (the certificate doc holds
  the definition and the bound, not the class split).
- `docs/ANALYTICAL_FORMULAS.md`: F154 (the locus saturation law; its own
  text still says the breaking interiors are "measured with slack", stale
  against e3dbab0's certificate, so this doc cites the certificate, not
  F154, for the containment); F140 (the frozen divisor: the Liouvillian
  EIGENVALUE λ = −4γ̄ pinned on the single-excitation corner block on
  the R₉₀ locus, F140's own phrase; "the N²-dimensional (1,1)-type
  corner block" is F153's phrasing for the same object; one frozen
  mode per balanced pair, the VALUE unconditional in the coupling,
  the COUNT at γ̄ ≠ 0 and generic J, Heisenberg-fenced; the SAME block label as this experiment's room and the
  SAME number as its certified edge rate 4γ̄, but a different object:
  F140 is a full-Liouvillian eigenvalue with locus-wide multiplicity,
  this doc's 4γ̄ is the fastest eigenvalue of the COMPRESSED
  dissipator on one ad_H eigenspace and appears only for maximizing
  transversals; the F-registry's own F153 note already flags the
  conflation risk, do not conflate); F91 (the R₉₀ locus definition and its "not the
  decay rates of anything" fence, which the certificate answers); F122
  (the ceiling; the compression machinery's home, "high-Q degenerate PT");
  F67 (spectral encoding is the protecting mechanism; receiver territory).
- `experiments/` prior flights: `IBM_CONCENTRATOR_RELOADED.md` (the genre
  parent: randomized-RZ dephasing injection, M-binding channel, billing
  laws measured: 257 s for 908k shots = 0.283 ms/shot on the A-vs-B
  jobs, and 119 s on the concentrator job itself, which is 0.316
  ms/shot over its 376,832 shots, that division being THIS document's
  recompute and not a figure the record states (round-15 precision;
  §10 already said so and this bullet did not),
  sign/magnitude verdict split, the stale-band and
  pooled-binding traps, the transport-dressed trap, the shared-skeleton
  invariant (the trap/invariant labels are this doc's shorthand for
  that RECORD's findings), and the channel-realization ladder quoted
  verbatim in §5);
  `IBM_F129_RAMSEY_FRINGE.md` (compile the preparation from the ACTUAL
  Trotter step's modes; its design-gate certificates, C1/C2 in this
  doc's shorthand, are free-fermion- and
  Slater-specific and do NOT transfer here, §5);
  `STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md` (two B-BLOCK-INVALID
  flights: mid-batch T1 telegraphing invisible to calibration-reading
  guards is the recorded killer; frozen CLEAN conditions, family error
  rates, VOID precedence, the calibration-only day-of re-gate with its
  anti-circularity fence, the DD-disabled guard, the ζ convention pin
  ζ_c = ζ_shift/4, and per-shot billing 0.309-0.327 ms on delay-bearing
  jobs); `RECORD_PARITY_HARDWARE_PREDICTION.md` (Class 1 pre-submit vs
  Class 2 in-job aborts; the guard bank = pooled statistics PLUS
  individual floors, the retired "pooled bank" shorthand; the
  factor-collapse floors under ratio statistics; per-device band sets for
  twin backends; the no-runtime-mitigation requirement with the
  resilience flag secondary);
  `PRICE_PAIR_HARDWARE_PREDICTION.md` (the conditional-Ramsey ζ
  measurement this pipeline's ZZ numbers come from).
- `fw.Confirmations` / C# `ConfirmationsRegistry` (24 entries): entry 23
  `concentrator_site_contrast_kingston_july2026` (engineered per-site
  dephasing, rate observable, differential arms, ibm_kingston); entry 24
  `f129_standing_fringe_kingston_july2026` (Givens-network-prepared
  two-branch 3-magnon cats in the Floquet eigenbasis, so the two-MODE
  single-magnon dyad prep here is the smaller sibling, §11a); `ibm_ep_onset_may2026` (per-site ⟨n_l⟩ from
  Z-basis counts only, no tomography). The nearest flown kin is
  entry 23 itself: slope(MP) − slope(E), a decay-rate DIFFERENCE
  between arms differing in where one engineered sink sits at equal
  dose. What has not flown is a rate split WITHIN one room read as a
  spectrum (three channels per arm) under a multi-site moved
  profile. "New" throughout this section means: absent from THIS
  repository's registry and flight records; no external-literature
  sweep is claimed.
- OpenArcs `compressed_density_laws` (Open): NextStep (3) closed both
  halves; open residue named there and inherited here: M0 diagonality and
  the 2160/2304 denominators measured, underived. `sideways_spin_ladder`:
  arithmetic kin (su(2) Clebsch-Gordan per its registry entry; the
  "different algebra, no overlap" judgment is this doc's own);
  its strip-retention precondition binds the registry's named future
  mod-p disc-layer device run at (1,4)@N=8, not this flight
  (round-4 scope precision).
- `docs/GLOSSARY.md`, `docs/CAUGHT_ERRORS.md`: no prior "corner beat"; the
  recorded dose-factor and convention error shapes feed §5 and §10.
- Swept again at the G7 stage (2026-08-17, the interpreter round): the
  flown circuit's ZZ-on-sector reduction IS the repo's **degree
  diagonal**, and **the repo already carries the inventory of who owns
  it**, so this bullet cites that inventory instead of rebuilding one.
  F152's own credit paragraph (`docs/ANALYTICAL_FORMULAS.md`, "What
  was already owned, and what is new here") reads: "D10 Step 3 derives
  the Laplacian on the open chain and F2 carries the
  ZZ-supplies-the-degree sentence, both at uniform γ;
  `VacuumBlockReductionClaim` carries the γ PROFILE but with the bare
  XY hopping and no degree term; PROOF_R90_FROZEN_DIVISOR Lemma 5
  derives the same N×N matrix independently, citing neither."
  Two holders that inventory does not list, both relevant here:
  `docs/proofs/PROOF_K_PARTNERSHIP.md` has the flown BOUNDARY, "in the
  single-excitation sector the ZZ-term contributes V_eff(ℓ) =
  (#bonds) − 2·deg(ℓ). Open chains have deg(0) = deg(N−1) = 1 <
  deg(interior) = 2, generating an effective on-site potential at the
  boundary", which is the free rz layer on qubits 0 and 5 that §5
  flies, appearing there in its "does not cover" list as the term
  that breaks **L1**, that proof's Lemma 1 (KHK = −H for bipartite
  nearest-neighbour hopping; "the chiral-conjugation lemma" was this
  doc's own gloss, unmarked, round-16 fix), K-partnership
  failing downstream of that (round-15 precision: the earlier
  wording said it breaks K, which that sentence does not say; the
  K statement is a separate one later in the proof). Mind the
  polarity either way: what disqualifies a chain in that proof is
  our circuit element. And `experiments/CONCENTRATOR_GEOMETRY.md`,
  this flight's own parent-genre record, in the same words ("the ZZ
  coupling contributes a degree diagonal, which turns the adjacency
  matrix into the graph Laplacian"). Note that D10 Step 3 is titled
  "On the (0,1) block the generator is the graph Laplacian" and closes
  at L|₍₀,₁₎ = −2iJ·𝓛 − 2γ·Id, so it holds the Liouvillian face and
  the Hamiltonian one together, the H₁ algebra being a step inside it.
  **What is new here is therefore neither the diagonal nor the
  boundary field, but only their per-Trotter-layer application** (odd
  layer constant, even layer a boundary field) and the hardware
  realization, which is what §11a element 1 claims and all it claims.
  A NOTE ON THIS BULLET, because it is the most-repaired passage in
  the document and the repairs are the lesson. Rounds 11, 12, 13 and
  14 each corrected it, and the first four versions all failed the
  same way: they built an ownership map by hand, sorting sources into
  buckets, when F152's credit paragraph already WAS the map, inside
  the document being cited. One of the four also concluded from a
  phrase search that "F2 holds neither", because F2's Note words it
  "supplies the diagonal shift that turns the adjacency matrix into
  that Laplacian" rather than "degree diagonal". Searching for the
  phrase produced a false denial about the STATEMENT, and the phrase
  itself is not the discriminator it looked like: it occurs in seven
  tracked files, this one included (round-15 correction of a
  round-14 sentence that called it unique to F152, which holds only
  INSIDE the formula registry). Cite the inventory; do not rebuild
  it, and never let a string search stand in for a claim about
  content.
  Two adjacent-genre mentions are NOT holders
  (`PROOF_ZETA2_ANTI_PROTECTION.md` / `PROOF_MIRROR_ORDER_SORTING.md`
  use "occupation-diagonal" for a ZZ perturbation, a different
  object). Precedent for
  this rediscovery shape: the OpenArcs entry
  `site_resolved_vacuum_block` ("it was not new, it was D10 reached
  by the same route").
- Local session scouts (gitignored, re-runnable, deliberately not named:
  a tracked document may not cite a gitignored path; every load-bearing
  number moves into the committed gate §8 before any flight). Numbers
  corrected across the rounds by recompute: v1's split fraction 96% and
  secular error 4% at Q = 10 (Q ≡ J/γ̄, defined in §4;
  fit-window artifacts); v2's ω = −4 null
  (unpopulated dyads) and secular sign; v3's envelope-immunity claim,
  (1,2)-null form, M-ladder indexing, and C-only corridor (round 3).

## 1. The claim, one line

On a six-site uniform-J open XXZ chain at Δ = 1, sliding engineered
per-site dephasing from the uniform profile to a mirror-transversal corner
profile at the SAME total dose splits the one beat-decay rate of the
(1,1) ω = −2 room, (10/3)γ̄, into three, with a width set by the
transversal's CLASS: γ̄·{8/3, 10/3, 4} for the maximizing class (split
statistic s = (2/3)γ̄, top rate at the certified centre-interval edge) and
γ̄·{10/3 − 2/√27, 10/3, 10/3 + 2/√27} for the non-maximizing class
(s = (2/√27)γ̄, exactly 1/√3 of the maximizing width; the middle rate is
exactly 10/3 in both classes), the mean over the room exactly unchanged in
every case. What the REGISTERED arms test, stated honestly (round-3
scope repair): D-sign tests THAT the corner splits the room; W tests
that the maximizing corner's split EXCEEDS the non-maximizing one's
(one-sided, s²(C) − s²(C′) > θ_W; s², the verdict statistic, is the
SAMPLE VARIANCE (ddof = 1) of an arm's three fitted rates,
equivalently one SIXTH of the summed squared pairwise differences,
defined in §6.4, and κ, the power condition's frozen factor, in §7),
which is an ORDERING, not the
1/√3 ratio: W is MAXIMIZED by C′ not splitting at all and D-sign
never reads C′, so the registered pair cannot distinguish the class
law from a dead C′
(the quantitative width law, the triples and the factor 3, lives in
D-mag, whose bands are CONFIRMED-conjunct but outside the
Confirmations registration, whose ¬D-mag branch is exactly the
width law's falsification face and the RECORD names it so, and
whose C′ side carries §7's discrimination caveat), plus the mean
anchor; the full rate
triples and
the edge value are reported context (§7), since the second-moment
statistic is blind to the triple's shape; the shape (one rate exactly at
the mean, the decoupled dyad) is reported through the CENTRED third
invariant (§6.4) with a frozen band, not as a CONFIRMED conjunct.

Wording fences:
- "Edge" means the edge of the size-class CENTRE interval
  [s_min, s_max]/2 · (−4γ̄) = [0, −4γ̄], NOT the edge of
  PROOF_CODIM1_BY_ADDITIVITY §6's block rate window (that proof
  derives its window at uniform γ only, at Δ = 0, and its "corner
  block" is the (p_c+1, p_c+1) object, not this room; the
  wider-under-a-profile
  comparison is THIS doc's own, not the proof's).
- The certificate's TIGHT(upper) attainment is an exact statement about
  the COMPRESSION; the finite-Q Liouvillian eigenvalue sits near the
  edge on EITHER side (measured: +0.06% to +0.09% above at the
  continuous candidate points; under the Strang step ~0.1% below at
  the RETIRED point and 0.61% below at the FROZEN one, round-3
  refresh), so a
  dressed top rate slightly off 4γ̄ is expected finite-Q behaviour, not
  a violation. No verdict arm reads absolute attainment.

## 2. The objects

- Chain: N = 6, open, uniform J, Δ = 1 (Pauli book H = Σ_bonds XX+YY+ZZ).
- One-magnon modes k = 0..5 ascending, E = {−2.4641, −1, 1, 3, 4.4641, 5}.
- The room: ω = −2 dyads (1,2), (2,3), (3,5) spanning the (1,1) block's
  dim-3 mixed eigenspace (modes 1,2,3,5 at E = −1, 1, 3, 5).
- Notation: mirror SITE pairs are written {a,b} ({0,5}, {1,4}, {2,3});
  mode DYADS are written (i,j) ((1,2), (2,3), (3,5)); the collision on
  "(2,3)"-shaped strings is real, so the brackets carry the type.
- Profiles on the R₉₀ locus (γ_l + γ_{5−l} = 2γ̄, total dose 6γ̄ in every
  engineered arm):
  - UNIFORM: γ_l = γ̄ on every site. Compressed matrix EXACTLY
    −(10/3)γ̄·Id (scalar; rests on Σ_l ψ_i(l)²ψ_j(l)² = 1/6, which holds
    for 13 of the 15 mode pairs at N = 6 including all three room dyads
    (only (0,4) and (1,3) give 1/12), an
    exact route per the house no-rounding rule, compared exactly in G1).
  - CORNER: a mirror transversal, one site of each pair
    {0,5}, {1,4}, {2,3} at 2γ̄, partner dark. **The eight transversals
    split into two classes** (recomputed from below in all three rounds):
    - MAXIMIZING (σ_outer ≠ σ_middle: the lit sites of pairs {0,5} and
      {2,3} on OPPOSITE ends): {0,1,3}, {0,3,4}, {1,2,5}, {2,4,5} →
      rates γ̄·{8/3, 10/3, 4}, s = (2/3)γ̄, attains w = 1/3. The split
      lives in the (2,3)↔(3,5) coupling; dyad (1,2) is EXACTLY decoupled
      and keeps the uniform rate.
    - NON-MAXIMIZING (σ_outer = σ_middle): {0,1,2}, {0,2,4}, {1,3,5},
      {3,4,5} → the rate triple γ̄·{10/3 − 2/√27, 10/3, 10/3 + 2/√27}
      (the MIDDLE rate is exactly 10/3: here the DECOUPLED dyad is
      (2,3)), s = (2/√27)γ̄ ≈ 0.385γ̄ (rate
      spread = 4γ̄ × the eigenvalue SPREAD of the transversal's odd
      part, which is 2× the OPERATOR NORM under the stated recipe
      C_p = comp(N_low − N_high): the recipe's difference matrix has
      eigenvalues {±1/6, 0}, spread 1/3 = the certificate's
      w = 0.3333, whose own C_l convention carries the factor 2
      (round-4 factor audit, symbolic); per mirror pair p
      with comp the certificate's NUMBER-OPERATOR compression (the
      dissipator's dose factor lives entirely in the 4γ̄ prefactor;
      building C_p from the rate compression double-counts it),
      under which C_1 ≡ 0 on this room, the maximizing combination is
      the DIFFERENCE with SPREAD 1/3 = w and the non-maximizing the
      SUM with spread 1/√27 (spreads, not operator norms, per the
      factor-2 sentence above; design-round-5 repair of a swapped
      sign); the numeric
      kinship with the census's resC = 1/√27 for (1,1) is an observation,
      not a provenance: the census norm is a single C_l, this one a
      two-pair combination), top rate 7% inside the edge. Here
      the coupling sits on (1,2)↔(3,5): the (1,2) dyad mixes, and WHERE
      the coupling sits is itself class physics.
    The pair {1,4} is free in both classes (C_1 ≡ 0 on this room).
  - Pinned profiles: **C = A = {0,3,4}** (maximizing) and
    **C′ = A′ = {0,1,2}** (non-maximizing). Class membership of the
    REALIZED profiles on the PHYSICAL qubits is a post-transpile hard
    assertion (§9).
- The conserved mean is the compressed trace: Tr = −10γ̄ for the uniform
  profile and all eight transversals. Off the locus at equal dose the
  trace ranges over [−10.67, −8.67]γ̄ (single-site vertices; −6.7%/+13.3%
  around −10γ̄), so the mean anchors the locus but does not confirm it
  alone.
- The two rooms are anti-correlated across classes: Σσ_l C_l vanishes on
  the ω = −4 room for exactly the four ω = −2-maximizing transversals.
  The ω = −4 room is NOT read by this experiment (no pinned preparation
  populates its dyads (1,3)/(2,5); the v2 "free null" was refuted).

## 3. Arms

Four shot-bearing arms plus one in-circuit null. All four arms share ONE
transpiled skeleton, differing only in the bound phase tables (φ = 0 for
N0; the concentrator's shared-skeleton invariant), so
they are gate- and duration-identical BY CONSTRUCTION, in the strict
sense: all four arms index the same transpiled circuit object, so the
identity holds the way an assignment holds, not the way a measurement
does. The runner's check compares those object identities, which
means it cannot fail today and is there as a REGRESSION guard: it
fires the day someone builds circuits per arm. Round-11 wording
repair, because "asserted in code" invited the reading that the
sameness had been tested. Interleave order is
DETERMINISTIC and pinned: within each depth, N0 → U → C′ → C, depths
ascending. The honest rationale (round-4 repair of a leftover pre-s²
argument): under the verdict statistic s², UNIFORM monotone drift has
exactly zero effect (shift invariance), and SCATTERED drift is unsigned
(G4 budgets it symmetrically), so no order can be claimed
drift-conservative; the order is pinned for determinism (record-parity
pins its intra-variant arm order and leaves inter-PUB order unpinned;
the pin-for-determinism move is what transfers), and the W pair C′/C is placed adjacent
because W carries the two-class law, the design's sharpest comparison.

- **U**: uniform engineered profile γ̄ per site.
- **C**: corner profile A = {0,3,4} (maximizing class).
- **C′**: corner profile A′ = {0,1,2} (non-maximizing class), same total
  dose, same depth grid: the live alternative hypothesis that makes the
  W-arm a discrimination between two from-below predictions.
- **N0**: no engineered dephasing (φ = 0 bindings). N0's fitted FULL 3×3
  room generator B̂ (not a profile, not a scalar: the site-profile
  inverse is rank-deficient, 5 of 6, and the real background need not be
  of dephasing form) enters the dressed predictions (§7). Because s² is
  shift-invariant, in the COMPRESSED MODEL s²_dressed(U) = s²(B̂) up
  to the §4 secular term (exact if the U engineered part were exactly
  scalar). MEASURED CORRECTION (round 8; the identity was registered
  before the hop-noise floor was found and the two cannot both
  stand): under the committed gate's within-sector hop error the U
  arm carries its OWN floor, s²(U) ≈ 0.0017 in the H1 chain (the
  0.0014 → 0.0003 hop-bracket pair is the gate's H0 configuration,
  a SECOND measurement of the same symbol, per §8a's
  reconciliation note), generated by hop
  error INTERACTING with the engineered dephasing, while N0, with no
  engineered dephasing, reads s²(N0) ≈ 4e-5 under the SAME noise
  model, a ~40× gap (0.0017/4e-5 = 42.5, round-17 arithmetic fix):
  so the U-side baseline and band are produced by
  the ANALYZE-CHAIN H0 of the θ_D refreeze, never by centring on
  s²(B̂); N0's B̂ carries the background PART of the dressed centres,
  var(B̂) still propagates into every
  band (G1), and N0 carries its own CLEAN conditions (§7).
- **In-circuit null**: the decoupled (1,2) dyad, read as a
  **C-minus-U generator difference** (notation: Ĝ(arm) is that arm's
  fitted full 3×3 room generator, the same fit shape whose N0 instance
  is named B̂ above): [Ĝ(C) − Ĝ(U)] must have (1,2) row
  and column consistent with zero (exact in the model for ANY common
  background B̂, which cancels in the difference; stated as a difference
  because under generic B̂ neither arm's own (1,2) row is zero and the
  (1,2) direction is not an eigenvector of either single-arm generator;
  caught in round 3). Design asset: sites 1 and 4 read dyad (1,2) nearly free of the
  other two dyads (rows (±0.2357, 0, 0) EXACTLY in the continuous-H
  basis; in the flown Strang basis the zeros become ±0.0023, 1%
  contamination, which the G2 margin carries). The DIFFERENCE reading is a wiring/class-membership check and
  may VOID (threshold G2); the (1,2) RATE itself is class physics and is
  REPORTED, never a VOID trigger (a genuine falsification must not be
  recorded as instrument failure). In C′ the (1,2) dyad mixes by class
  physics, but C′ has its OWN exactly decoupled dyad, (2,3): the
  analogous C′-minus-U difference null on the (2,3) row/column applies
  there (round-5 addition: the only in-job wiring check behind C′'s
  Class-1 class assertion), with the caveat that no site reads (2,3)
  contamination-free (rows 1 and 4 are zero on it), so it needs the
  full three-dyad inversion, unlike the (1,2) null.

## 4. The working point and the corridor (two factors, opposite directions, per arm)

DEFINITION, used everywhere in this document and never before stated
(round-4 repair): **Q ≡ J/γ̄**, the dimensionless working-point knob
(hops per engineered dephasing; the flight's γ̄ = J/10 is Q = 10; the
canonical hardware regime γ₀ = 0.05, J = 0.075 is Q = 1.5, §11).

The corridor is the PRODUCT of two effects (measured estimator-free in
rounds 1-3, agreeing on mechanism and magnitude; the committed gate
re-establishes the numbers through the committed estimator):

- **Secular factor, per arm and opposite in sign:** the finite-Q
  compression error moves the two arms in OPPOSITE directions: it
  INCREASES s(C) (+0.60% at Q = 20, +1.72% at Q = 12, +2.51% at Q = 10;
  re-signed in round 2) and DECREASES s(C′) by roughly twice as much
  (−1.3% at Q = 20, −3.5% at Q = 12, −5.1% at Q = 10; caught in round 4,
  where the v3 bullet was still C-only). The per-rate error is smaller
  (≤ 1.2% down to Q = 10; 1.116% exact), sits mostly on the lowest rate, and is NOT
  the band-setting quantity; the split's error is.
- **Trotter factor, per arm:** the Strang step detunes dyad (3,5) by
  δ_T ≈ 2.0·(J·dt)² (exact Floquet 0.020050 / 0.045250, implied
  coefficient 2.005-2.011), an ANGULAR FREQUENCY in units of J, commensurate
  with the rates (renamed from
  Δ in round 6: Δ is the XXZ anisotropy alone, and a reader taking the
  detuning as a per-step phase concludes falsely that the candidate is
  past the EP). For C the
  detuning competes with the coupling c = (2/3)γ̄ in its 2×2 block:
  half-split RATE √(c² − (δ_T/2)²) (its ratio to c is the factor
  ≤ 1) that reproduces C's rows
  to 2.4-4.4% across the four table rows (2.42 / 4.41 / 3.61 / 4.13,
  writer-recomputed from √(c²−(δ_T/2)²)/c against the table on the
  EXACT Floquet δ_T this paragraph pins two lines above, round 15;
  round 14's 2.4-4.3% came from the rounded 2.0·(J·dt)² instead,
  which is the same paragraph reading its own inputs two ways, and
  both replace a "~4-5%" that not even the worst row reached). Two
  rows, Q = 10 and Q = 20, exceed 4%. They are NOT thereby outside
  the "2-4% method spread" quoted elsewhere in this document
  (round-16 repair of a comparison that was never commensurable):
  that spread measures the compressed map against exact eigenvalues,
  while these four numbers measure the 2×2 Trotter MODEL against the
  Floquet table, and nothing bounds the second by the first. The
  table below is the authority. For C′ (c′ = c/√3, and the coupling
  on a DIFFERENT dyad pair) the 2×2 model fails at EVERY tested Q
  (round-3 recompute: at Q = 10 it predicts 0.81 where the table
  measures 1.03, at Q = 15 it predicts 0.47 vs the table's 0.82, and
  past its own model-EP it is undefined where the table still reads
  0.27): the 2×2 is C's mechanism only, never a C′ law. The table is the authority
  for the ROOM SPECTRUM; the committed estimator adds its own
  measured dressing on top (pinned-gauge, noiseless, the committed
  gate's frozen mode: **s² = 1.21× bare ideal for C, 0.54× for C′,
  ~0 for U**; from the class-basis channel projection and a
  prep-dependent secular population leak, deterministic, carried by
  f_C/f_{C′}; §7, where the centre functions live),
  so the model is C's mechanism, not a per-arm law, and the TABLE below
  (direct Floquet-Lindblad) is the authority. **Each arm has its own
  corridor and its own exceptional point**, C′'s EP at markedly lower
  Q; beyond an EP that arm's split VANISHES. Higher Q is NOT safer
  here.

Full Floquet-Lindblad values, in s (NOT s², the verdict statistic
defined in §6.4 as the SAMPLE VARIANCE (ddof = 1) of the three
fitted rates, equal to one sixth of the summed squared pairwise
differences; s is its square root, and this table is the one place
the doc quotes the root of a MEASURED quantity, §1 and §2 quoting
the ideal s = (2/3)γ̄ and (2/√27)γ̄ (round-17 scope fix of a "one
place" that was false as written), hence the warning. Round-16 repair of both
prose definitions: they said "the mean square of the pairwise
differences", which is the sum over THREE and therefore twice s².
The check that catches it is §2's own numbers, {8/3, 10/3, 4}γ̄ →
s² = 4/9 and s = 2/3, the registered split statistic. Round-14 forward pointer: s² appears in
§1, §3 and here before §6.4 defines it, the same gap Q and f_leak got
pointers for in earlier rounds; the ratio column is
context only since v6 retired the ratio arm; the C′ cell at Q = 20
reads 0.273-0.274 across independent recomputes, steep near its EP,
and the committed gate mode settles the digit):

| J·dt | Q | s(C)/ideal | s(C′)/ideal | (s(C)/s(C′))/√3 |
|------|----|-----------|-------------|------------------|
| 0.10 | 12 (retired) | 1.008 | 1.027 | 0.982 |
| 0.15 | 10 (CANDIDATE) | 0.984 | 1.030 | 0.955 |
| 0.15 | 15 | 0.893 | 0.816 | 1.094 |
| 0.15 | 20 | 0.766 | 0.274 | 2.80 (C′ near EP) |

EP locations at J·dt = 0.15 (the candidate's row): C at Q ≈ 30 exact
(2×2 model 29.5, conservative by ~0.5 in Q, not more: Q = 32 is
already collapsed); **C′ at Q ≈ 21, the BINDING margin at the flown
Q = 10 is C′'s, a factor 2**. (At the retired J·dt = 0.10 point the
margins were C ≈ 67 and C′ ≈ 42.5 exact; the earlier quoted ≈ 38
was the 2×2 model's value for C′, round-3 provenance repair.) All bands come from G1's counts-level sim,
never from this table alone (§7: the post-selection envelope biases
absolute s² and no bare-corridor band survives it).

Method spread: compressed-map and exact-eigenvalue readings differ by
2-4%; the committed gate's counts-level estimator is the arbiter.

Candidate working point (chosen by G1 v2's two-stage selection on the
VERDICT metric P(d > θ_D), held-out seed, both f_leak ends, where
f_leak is the OUT-OF-SECTOR fraction of the per-gate error weight,
bracket [8/15, ~0.9], counted and derived in §10; the frozen gate run
records its near end as 0.53 rather than 8/15 = 0.5333, a rounding in
the run's own input that is numerically immaterial and worth naming
only because the registered number and the frozen run's number are
then not the same string (round 15) (forward
definition, round 8); §8a): **Q =
10, J·dt = 0.15, grid 21 points × 3 steps (deep end 60), 16384 total
shots per (arm, depth, prep) (Tom's second budget decision, on the
post-pin numbers), FRACTIONAL-RZZ gates (2 two-qubit gates
per XY block, §5/§10; the corner-beat runner itself implements
use_fractional_gates=True with the CZ = NO-FLIGHT assertion; the
record-parity pre-registration REQUIRES the same flag but records
that no in-house runner had set it before, so this runner is the
flag's first in-house carrier)**: at ~23.7 QPU min
(the frozen grid's science cost 22.5 plus Amendment 1's aux/CAL
plan, inside the accepted 21-25 band; verdict-grade only AFTER the
≥ 500-rep refreeze, header caveat) under the GAUGE-PINNED committed
estimator: **worst-end P(detect) = 0.920 (power 2.9, H0 false 0.000
held-out); the other bracket end's (1.000/0.020) pair was measured
at its own 3σ value and is superseded by the single registered θ_D,
§7** (the 8192-shot variant: P(detect) 0.815 at ~11 min; the
pre-pin 0.985 was gauge-favored and is withdrawn). **Grid AND shots
FROZEN by Tom's two budget decisions (2026-08-16), the second taken
on the honest post-pin numbers.** The v4 fear that the
deep-end exponent kills the flight was WRONG in mechanism: the weighted
estimator devalues the dying deep points instead of dying with them.
Standard gates (3 per block) reach worst-end power only ~1.7 at 11
min and ~1.4 at 22 min under the pinned estimator (two different
stage-1 configurations, not one config at two budgets; the pairing
is the gate record's and lands in the committed G1-stage record at
the refreeze, round-5 note): fractional-RZZ is
LOAD-BEARING and its availability on the
flown backend is a Class-1 guard, as is the calibrated 2q error on the
used edges (p2 ≤ 0.5% since Amendment 1.8, which puts the bound ON
the measured p2 ladder's first rung rather than between rungs; the
ladder is §8a's). The former candidate Q = 12, J·dt = 0.10
underperforms at every tested budget and is retired.

**J-calibration sensitivity** (G5 v2.1, estimator level, nominal-J
preparation basis, pinned gauge, the flight's situation): the response
is EVEN, ∂s²/s² = +4.5% for both ±5% J (a one-sided inflation bias;
§8a); bond-J scatter σ_J = 1% costs ≤ 10% of s²(C) (95%, frozen-draw
systematic), σ_J = 2% ≤ 25%. The band carries these terms (G5). J is VERIFIED IN-ANALYSIS from
the science arms' own fitted dyad frequencies (the Floquet frequencies
are J-sensitive and fitted anyway; estimator frozen by G5; no extra
PUBs), a precondition on the magnitude/anchor arms (D-mag both and
A; the one registered scope is §7's SCOPE list).

## 5. Protocol

**Backend and line.** ibm_kingston or its pinned twin ibm_marrakesh
(Heron r2). OPERATIONAL NOTE (round 7): the runner pins
BACKEND_NAME = ibm_kingston and hard-aborts on any other chain
file; switching to the twin is a SOURCE EDIT plus a re-commit
through the Class-1 commit gate, a pre-day decision, never a
day-of selection. Six-qubit line by the concentrator uniform-line rule (§9,
Class 1). **Per-device band sets:** G1 freezes a band set per twin;
the flown device's set governs (the record-parity twin rule). HOW
THAT IS CARRIED, since round 15 found it unimplementable as
registered: every band key in the constants file is a SCALAR with no
device dimension, so a second set has nowhere to live. Because
switching twins is already a source edit plus a re-commit through
the Class-1 gate, the twin's set is carried the same way, as the
edited constants file committed at that hash, and the flown file is
the flown device's set by construction. G1 freezes both sets and
records the unflown twin's beside the flight, but only ONE set is
ever live in the file the gate reads. BOTH twins
must expose fractional rzz for the twin rule to hold; the pre-gate
properties query checks it on both, and a twin without it drops out
of the pinned pair (no silent standard-gate fallback). **Runtime
mitigation OFF:** the pinned requirement is NO runtime mitigation and
no scheduler dynamical decoupling (DD would rewrite the very dephasing
profile under test; offline CAL inversion would double-correct);
per the record-parity wording the resilience flag is evidence, not the
assertion (it is not guaranteed a Sampler knob), so the runner asserts
the requirement on the transpiled circuits. The FLAG itself is
fail-closed, not best-effort (round-16 alignment with §9 and the
code, where "set where exposed" was the pre-G7 design: the runner
raises if the options object refuses the setting, and then asserts
that DD and both twirling channels read disabled).

**Preparations.** Three, one per room dyad: (ψ_i + ψ_j)/√2 is itself a
SINGLE-PARTICLE orbital, prepared by N − 1 = 5 Givens rotations from
|100000⟩ (no cat block; v1's cat was F129's three-particle structure
mis-imported). Orbitals from the eigenvectors of the ACTUAL Strang step
(the f129 lesson that does transfer). The f129 design-gate certificates
(here labelled C1/C2, our shorthand) do NOT
transfer (C2 is a Slater-Condon zero needing 3 differing orbitals; C1 is
the pure-hopping chiral relation, and its non-transfer at Δ = 1 is
this doc's inference, not the f129 doc's statement). Replacement
certificates, committed in the gate: prep fidelity of the compiled
circuit against the Strang-step eigenvector, and a leakage bound out of
the one-magnon sector under the interacting step.

**Time evolution.** Second-order (Strang) odd/even Trotter at the frozen
J·dt. First-order is excluded by measurement (spurious uniform-arm
half-split 13-25% of ideal across the grid vs Strang's ~0.6-2.3%,
the top of that range at the frozen point, round-3 refresh).
**The flown circuit realization (G7, verified at machine precision
against the committed gate's strang(dt) on the sector):** each Strang
substep flies its XY parts as two-qubit blocks (2 fractional rzz per
block) while the ZZ part rides the ONE-MAGNON SECTOR as the repo's
**degree diagonal** (F152, whose phrase this is; D10 step 3 carries
the derivation; PROOF_K_PARTNERSHIP holds V_eff(l) = #bonds − 2·deg(l)
(with parentheses, `V_eff(ℓ) = (#bonds) − 2·deg(ℓ)`), there in its
"does not cover" list as the term that breaks **L1**, the lemma
KHK = −H for bipartite nearest-neighbour hopping, K-partnership
failing downstream of that in a separate statement later in the
proof (round-16 repair: round 15 corrected this attribution in §2's
sweep bullet and left THIS copy, the one an executor reads while
building the circuit, saying "breaks K")), applied per Trotter layer: the odd-bond layer is a perfect matching, its ZZ is a CONSTANT
(a global phase); the even-bond layer leaves the chain ends uncovered,
its ZZ is 2n₀ + 2n₅, a free rz layer on qubits 0 and 5. That boundary
rz IS the boundary condition: dropping it would fly a different
chain's spectrum (PROOF_R90_FROZEN_DIVISOR's "boundary clocks", its
Lemma 5 carrying both, the
N-vs-N+1 modulus between the Heisenberg and XY single-excitation
Hamiltonians; the registry's F2-vs-F2b pair shares the N-vs-N+1
denominator and nothing else, since F2b IS an SE Hamiltonian
spectrum by its own title, E_k = 2J·cos(πk/(N+1)) over N modes,
while F2 is the Liouvillian (0,1) block, 4J(1 − cos πk/N) over
N−1 oscillating ones: the contrast is Liouvillian-vs-Hamiltonian
AND Neumann-vs-Dirichlet, not one modulus shape on one kind of
object. Round-15 repair of a round-4 sharpening that had put both
on Liouvillian blocks.)
§2's pinned mode list is PROOF_R90_FROZEN_DIVISOR Lemma 5,
λ_k = 4cos(kπ/N) + N − 5, equivalently (N − 1) − 2μ_k with μ_k the
path-Laplacian eigenvalues (5 − 2μ_k at N = 6 only; F152 states its
spectrum via the same μ, not via this shifted form); the SET is
identical but
the index runs opposite (this doc's k is energy-ascending, Lemma 5's
descending, k ↦ N − 1 − k), so map values, never indices. The
flown unitary equals the full XXZ Strang step ON the sector and
differs off-sector, where post-selection discards the state; §10's
"2 gates per block" presupposes exactly this reduction. The final
injected-RZ layer is a Z-basis no-op (realized dose at depth n is n−1
layers), absorbed by the estimator's free amplitude FOR n ≥ 1 ONLY:
round 3 measured that fitting on the nominal axis t = n·dt leaves
depth 0 (dose 0, the fit's highest-weight point) sitting e^{−r·dt}
≈ 5% below the model, biasing every rate ~−4 to −6% low and giving
away 16.6% of s² noiselessly; the −1.4%/−3.4% previously booked as
"estimator dressing" (§8a) is dominantly THIS off-by-one (round-10
honesty: the −4 to −6%, the 16.6%, and the −1.4%/−3.4% are THREE
measurements from different constructions, weights and channels;
they do not close under one common factor, and the fixed-axis
re-fit is what measures the true residual dressing). The
refreeze re-fits on the REALIZED-dose axis t_eff = max(n−1, 0)·dt in
gate and runner together (registered here, executed as part of the
θ_D refreeze, pre-data; every dressed centre re-measures on the
fixed axis). The axis is now a SWITCH the frozen constant throws,
not a repair someone must remember to make (round 11): the runner
derives its fit axis from its `time_axis` constant, "nominal" =
n·dt and "realized_dose" = max(n−1, 0)·dt, and analyze reports the
axis in force unconditionally and audits it against the frozen
entry. THE GATE SIDE NOW HAS ITS MACHINE CHECK TOO (2026-08-17, the
first piece of the refreeze package built). It did not: `--certify`
handed ONE axis array to both fit implementations, so it certified
that the two FITS agree on a given axis and never that the committed
gate builds the same axis, and the gate built the nominal axis inline
at FOUR fit sites, none of them switchable. The gate now derives its
axis from the same switch the runner uses
(`corner_beat_gate.fit_time_axis`), and certify's check (8) compares
the two constructions on BOTH settings and the two switch constants
against each other, so a refreeze that throws the switch on one side
only is caught before a flight rather than inside the fits. The
comparison is EXACT, since both sides multiply the same integers.
The refactor moves no committed number, and that is checkable rather
than asserted: the new function reproduces the replaced expression
bit for bit at every (steps, k, dt) the gate sweeps, and `ts` enters
the fits only through it. What the machinery entry still holds is the
EXECUTION, re-fitting at the ≥ 500-rep refreeze with the switch
thrown, which is also what θ_D's `refreeze_required` and the unfrozen
`time_axis` entry block on. Before that,
freezing the manifest entry to "realized_dose" moved a string that
the estimator never read, so the manifest could record this repair
as landed while every fit still ran on n·dt. The depth grid is
UNIFORM in Trotter steps: 21 points at spacing k = 3 steps, starting at
step 0 (the realized dose at depth n ≥ 1 is n−1 injection layers;
depth 0 carries no injection layer at all, Amendment 1.4) and
ending at step 60 (§4/§8a; frozen by Tom's budget decision). The
committed estimator (§6) fits on that grid; the per-grid-interval room
phase is |ω|·k·dt ≈ 0.9 rad < π/3 (frequency identification safe for
the DETUNED dyad (3,5); the other two room dyads are EXACTLY
frequency-degenerate under the Strang step and no phase criterion
separates them, nor needs to: the eigenchannel fit never asks it
to, §6.3's withdrawal note). No
grid point is dropped after data exists (the concentrator rule: failing
points are dropped before pre-registration, never after).

**Dephasing injection.** The concentrator recipe verbatim, convention
pinned: RZ(φ) = exp(−iφZ/2); per step, per site, φ ~ N(0, σ_l²),
σ_l = √(4·γ_l·Δt_step), giving per-step retention e^{−2γ_l·Δt} per site
and e^{−2(γ_l+γ_m)·Δt} on a two-site coherence (recomputed, rounds 1-3).
**Channel-realization ladder, quoted verbatim from the parent flight**
(`IBM_CONCENTRATOR_RELOADED.md`: "single-realization slope gap 0.0208
(28.3% of the channel slope) ... the single-realization gap shrinks with
M (0.0208 → 0.0032 → 0.0010 at M = 256/1024/4096)", recorded there as an
interpretive decomposition, not a gate number): so ~28% at M = 256,
~4.4% at M = 1024, ~1.4% at M = 4096, on THAT flight's slope observable
(note the quoted ladder falls FASTER than the 1/√M that §8a's
cumulative-retention item reasons from: two different observables,
the parent's slope gap vs this design's retention floor, quoted
side by side deliberately, round-3 note).
**M = 1024 candidate** (the parent's slope observable shrinks to
~4.4% there, but the quantity §8 freezes M from is THIS design's
own s²-realization spread, measured 15.3% at M = 1024 by G3 v2,
§8a: price M from 15.3, never from 4.4, round-5 repair); G3 freezes M from the realization
error measured on s² itself, and the §10 cost band re-freezes with it
(certificates computed at a provisional M never carry; note the
concentrator's flown FINAL was M = 256, and its billed 0.316 ms/shot
sits inside the delay-bearing band with no gross per-binding term
visible, though that job's own aux delays make the separation
inexact, §10; the per-binding overhead question stays priced to
G3 anyway). Shots are TOTAL per (arm, depth, prep), split across the M
bindings (16384/1024 = 16 per binding at the frozen numbers; the
concentrator flew 32/binding); per-binding counts persisted (bypassing
get_counts() pooling, the concentrator's recorded instrument deviation).
Exactly one injected RZ per site per step, asserted post-transpile.

**Readout.** Z-basis counts only, per-site ⟨n_l(t)⟩. Per-site resolution
is MANDATORY: the three dyads' site-amplitude columns sum to exactly
zero (total magnetization carries no beat). Sites 1 and 4 are exact
nodes of mode 3 and blind to dyads (2,3) and (3,5): four of six sites
carry the split (and exactly those two sites carry the clean (1,2)-null
readout, §3). CAL0/CAL1 confusion PUBs fly in the same job (per Amendment 1.5 a
CAL pair additionally flies in EVERY job, and the HEAD pair plus
the per-job pairs pool into the inversion; the station CAL pairs
stay out as drift diagnostics, round-7 alignment with Amendment
1.5's letter); the confusion inversion is linear, so it commutes with
pooling (SUPERSEDED IN FORM by Amendment 1.1: the flown chain pools
the per-binding raw counts, inverts once, clips once; the clip is
the point where the chain becomes nonlinear, answering the question
§8 G2 carried), and G2 must show the
post-selection identity survives asymmetric per-qubit confusion.
Post-selection on total excitation = 1 is pinned in: it removes the
LEAKAGE component of gate error exactly. The within-sector component
survives (7 of 15 two-qubit Pauli errors preserve the excitation count
ON a differing pair: IZ/ZI/ZZ outright, XX/XY/YX/YY on the one-magnon
pair subspace only) and is a G1 line item. The residual post-selection
envelope (1−p)/(1−0.90625p) (0.90625 = 1 − 6/64, the out-of-sector
fraction under full depolarization) is common to all dyads and arms but
NON-EXPONENTIAL: it biases ABSOLUTE s² by up to ~−20% at deep grids
(measured, round 3) while cancelling in the s²-ratio (to < 0.1%) and in
r̄ differences; which is why D-mag's band comes from G1's counts-level
sim and never from §4's envelope-free corridor (§7).

## 6. Committed estimator (frozen in structure; constants from the gate)

(Reference convention: the numbered items below are cited as
§6.1..§6.6 throughout the document.)

1. Counts → per-binding RAW counts persisted → pooled over the
   (resampled) bindings → ONE confusion inversion → ONE clip →
   post-selected one-magnon site occupations n_l(t) (the order as
   flown, Amendment 1.1; the pre-amendment per-binding-inversion
   wording is superseded here too, round-5 alignment).
2. Site→dyad inversion: the 6×3 map A[l, d] = ψ_i(l)ψ_j(l) has exactly
   orthogonal columns (Gram = (1/6)·Id, cond = 1.000); implemented as
   the bare orthogonal projection 6·n·A, which EQUALS the
   intercept-fitted least-squares solution exactly (the intercept
   column is orthogonal to all three dyad columns, Σ_l ψ_i(l)ψ_j(l) = 0
   for i ≠ j, and a common offset is annihilated; verified to 1e-12);
   the affine remainder is absorbed by §6.3's free baseline b.
3. THE COMMITTED FIT (amended in v7; the formerly pinned complex
   one-step propagator fit was withdrawn for measured cause, recorded
   in §14: blind pole retrieval cannot separate three near-degenerate
   frequencies over ~2 periods, and quadrature demodulation is biased
   by the decaying envelope at the 8%-per-sample level, both measured
   in the gate's estimator trials): the EIGENCHANNEL damped-cosine
   fit. Per arm, the class-predicted eigenchannels (C: y₁ decoupled,
   v± = y₂ ± y₃; C′: y₂ decoupled, u± = y₁ ± y₃; U: y₁, y₂, y₃;
   N0 currently runs U's scalar channels, and the FULL 3×3 room
   generator B̂ that §3 requires, and that both §7 difference-null
   triggers and every dressed centre consume, has NO committed
   estimator yet: it is G2-pending machinery, so "frozen in
   structure" covers steps 1-6 for the four scalar-channel arms and
   NOT the generator fit, and the two difference-null triggers become
   evaluable only when G2 lands, §8a) are
   each fit as y(t) = b + a·e^{−rt}·cos((ω+δω)t + φ), linear in
   (b, a·cosφ, a·sinφ) at fixed (r, δω), a bounded 2-D grid search
   with refinement over (r, δω), traces weighted by kept counts
   WITHIN each trace (√w in the design matrix, each trace normalized
   by its own peak kept count, so jointly fitted preparations enter
   at equal overall scale regardless of absolute counts; the
   committed behaviour, certify fit-parity 0.0 on single AND on
   JOINT two-trace cases at 40× unequal peak weight, round-11
   repair: every parity case had been a single trace, where the
   per-trace normalization is a no-op, so the certified 0.0 rested
   on evidence that never touched the branch this sentence claims;
   the joint case now runs and is also exactly 0.0), joint over the
   preparations that carry
   the channel. The estimator ENCODES the arm's class hypothesis
   (legitimate: the profile is chosen, not measured). The earlier
   safety rationale "a wrong basis shrinks the split, the
   conservative direction" was measured BACKWARDS in the G7 rounds:
   the forced channels INFLATE (s²(C) reads 1.21× bare ideal
   noiseless, because the ± combinations are not eigenchannels of the
   Trotter-detuned block and the single-cosine fit trades the
   mismatch into rate and δω). The actual safety statement, measured:
   the forced C-channels applied to UNSPLIT data return
   s²(C) − s²(U) = +0.000000 exactly, and ≤ +0.00035 under scattered
   backgrounds to 0.5γ̄ (against θ_D = 0.00253), so the basis cannot
   manufacture a split; the inflation is deterministic dressing
   carried by f_C/f_{C′} per the corridor-dressing rule, and the
   designed-in |δω| excursion (up to ~0.03 noiseless) is why the
   fit-health δω margin must be frozen from the ideal construction's
   own excursion plus H0 (G2), not set at the search grid's edge.
4. **s² ≡ ((r₁−r₂)² + (r₁−r₃)² + (r₂−r₃)²)/6** on the three fitted
   channel rates: sort-free by construction (symmetric), and
   NON-NEGATIVE by construction (the v5 negative-s² machinery is
   withdrawn with the estimator that produced it; under H0 the
   statistic has a POSITIVE noise floor, measured by G1's H0 runs and
   absorbed into θ_D, which is why θ_D is frozen on the H0
   distribution of the DIFFERENCE and not at zero). Mean rate
   r̄ = (r₁+r₂+r₃)/3. The shape statistic Π(r_i − r̄) (shift-
   invariant, predicted exactly 0 for both classes' symmetric
   triples) is REPORTED with a frozen band, not a CONFIRMED conjunct.
   The fitted δω per channel is the detuning readout feeding the G5
   J-estimator.
5. Rates are invariant under any fixed LINEAR mis-calibration of the
   site→dyad map (similarity transform); the affine part is handled by
   the intercept in step 2.
6. Resampling unit for all error bars: the BINDING, resampled jointly
   across depths and preparations within an arm. s² resamples cleanly
   (no folding); the sorted rates do not and are presentation only.

## 7. Verdict rule (structure frozen; bands, thresholds, margins, and family rates from the gate)

Arm types are heterogeneous and each carries its own decision rule
(round-3 repair: the staircase's HOLDS/MARGINAL/VIOLATED grading is a
NULL-statistic shape and applies only where a null is tested):

- **D-sign, one-sided detection.** s²(C) − s²(U) > θ_D (θ_D frozen
  by G1 as E[d | H0] + 3σ̂ over the H0 draws of this difference,
  ONE formula (round-4 unification): the freezing run's H0 mean was
  consistent with zero so the frozen 0.00253 printed as the bare
  3σ̂, and the refreeze states the mean explicitly; σ̂ measured, not
  a quantile, §8a states the estimator; both arms carry a
  positive noise floor and the s²(B̂) offset, which cancel in
  expectation but not in variance; the detection is on the difference,
  never against an assumed-zero baseline). The REGISTERED θ_D is ONE
  number: the MAXIMUM of the per-end thresholds over the f_leak
  bracket (f_leak is unknown at flight time), the worst-end value
  **θ_D = 0.00253** (§8a; the frozen record also holds the other
  end's own 3σ value 0.00063, which the flight does NOT use). Both
  ends' operating characteristics must be read AT the registered
  threshold: worst end P(detect) 0.920 at held-out H0 false 0.000
  (n = 50, §8a caveat); the far end's published pair (1.000 / 0.020)
  was evaluated at its own 0.00063 and is superseded; at 0.00253 the
  far end's margin is (0.00554 − 0.00253)/0.00133 ≈ 2.3σ, and the
  exact re-evaluation at the registered threshold is part of the θ_D
  refreeze (§8a). DECISION RULE AS IMPLEMENTED (registered here):
  D-sign passes on the POINT estimate d > θ_D; the interval beside
  it is a 95% bootstrap percentile CI (R_BOOT replicates, resampling
  BINDINGS per arm, Amendment 1.4's nested tables as the unit; the
  NaN-replicate EXCLUSION rule is registered as implemented, round
  5: a replicate whose fit dies drops from the percentiles, the
  dropped count prints as a health number, PER KEY (round 11: it was
  D-sign's count alone while the percentiles drop replicates from
  every key, so W's and A's C′ line's exclusions were invisible), and
  its bound freezes with the G2 fit-health margins as the manifest
  entry `boot_nan_replicate_bound` (round 12: "freezes with the G2
  margins" named none of the three fit-health keys and no evaluator,
  so the count could be computed and persisted forever without
  anything able to block on it). The SAME MECHANISM one level down is
  counted too (round 13): a resampled pool can itself dip under the
  kept-count floor at the deep end, where ~208 kept counts are expected
  at the p2 guard bound (Amendment 1.8), which drops that grid point from that
  replicate's fit; those replicate-level floor trips are now counted,
  printed and persisted beside the NaN counts rather than happening
  in silence. Both matter for the same reason: silent exclusion
  narrows
  the CI in exactly the direction that turns INCONCLUSIVE into
  FALSIFIED; the bootstrap seed rides in every artifact; and the
  DEPTH-0 UNDER-COVERAGE is registered, round 6: the unswept
  depth-0 cell is single-PUB and contributes NO binding-resampling
  variance while being the fit's highest-weight point, so the CI
  systematically narrows in the FALSIFIED direction; the refreeze
  either adds a parametric shot-noise resample for depth-0 cells,
  a machinery-ledger item, or measures the term and WIDENS THE CI
  by it. ROUND-15 SIGN REPAIR, and it is the first finding in four
  rounds that sits in the original material rather than in a repair:
  the fallback used to read "inflates θ_F/κ by it", and both
  inflations point the wrong way. FALSIFIED is containment of the CI
  in [−θ_F, θ_F], so a LARGER θ_F makes containment easier and
  raises exactly the false-FALSIFIED rate this bullet exists to
  bound; and the power condition is SE ≤ κ × projected, so a LARGER
  κ loosens a test whose input SE is already biased low. The
  conservative direction is the CI, or equivalently a SMALLER θ_F
  and a TIGHTER κ. An executor implementing the old wording would
  have increased the false-falsification rate while believing the
  correction conservative), and
  the CI, not the point, carries the FALSIFIED / Anti-D /
  INCONCLUSIVE partitions below. If the refreeze returns a different
  θ_D, it supersedes 0.00253 as a NUMBERED pre-data amendment (§13),
  never silently.
- **W, the class law as a LEVEL discrimination (restated in v6: the
  ratio form failed its own G1 reachability check; re-measured in
  v7.1 with the committed gauge-pinned estimator).** The gate measures
  s²(C′) BELOW the U arm at budget (s²(C′) − s²(U) = −0.0015: the
  hop-noise floor of the degenerate U arm exceeds the whole C′
  prediction), so neither a ratio nor a C′-detection arm is freezable.
  The class law lives in the exceedance: s²(C) − s²(C′) > θ_W
  (one-sided; PRODUCER: θ_W freezes INSIDE the θ_D refreeze package
  from the same ≥ 500-rep analyze-side H0 chain, §8a, with the
  both-arms-null H0 for the family rates, and BRACKET RULE
  registered like θ_D's (round 9): θ_W = the MAXIMUM of the
  per-end values over the f_leak bracket, W's worst end being the
  OPPOSITE one from D-sign's; the gate's G1 stage
  measured the exceedance but froze no θ_W; measured at the frozen config
  d_W = +0.00634 ± 0.00134 (unrounded so the division reproduces
  the power), power 4.7 at f_leak = 8/15, D-sign's worst
  end; W's OWN worst end is the other one, power 3.8). The C′
  containment itself lives ONLY in D-mag (v7 repair: a
  duplicated band in two arms broke the branch partition and the
  registration rationale); W's verdict content is the one-sided
  exceedance alone, which keeps W a non-magnitude arm and the
  registration scope (D-sign + W) coherent. (The centre functions
  f_C and f_{C′} belong to D-mag, where they are registered; a
  sentence placing them here survived the v6 retirement of the ratio
  arm and read as if W consumed them, which these same three lines
  deny. Removed in round 15.) The s²(C) floor,
  DEFINED here (round 7: it was registered, scoped, and owed
  without ever being stated; its ratio-arm ancestry died with the
  ratio in v6 and its function under the LEVEL statistic is this):
  a frozen ABSOLUTE lower bound on the measured s²(C), frozen by G1
  from the H0 distribution of s²(C) itself through the analyze
  chain (the level below which C's own reading is indistinguishable
  from the no-split noise floor); it TRIPS when the measured
  s²(C) < floor; its purpose: the W exceedance is uninformative
  when C itself did not split, and the trip decides between the
  data path and VOID per the adjudication.
  ITS PRODUCER, one formula, registered here (round 11: the floor
  was defined in round 7, given two forcing powers in rounds 8 and
  9, and never given a producer, while its two sibling thresholds
  both have one; since round 9 the operative detection rule is the
  CONJUNCTION d > θ_D **and** s²(C) ≥ floor, so a floor left to
  G1's choice of quantile is a researcher degree of freedom on the
  detection criterion itself, which is exactly what pinning θ_F's
  base closed in round 4):
  **floor = E[s²(C) | H0] + 3·σ̂(s²(C) | H0)**, the same 3σ
  separation-from-the-null shape θ_D carries, evaluated through the
  analyze chain, and, like θ_D and θ_W, taken as the MAXIMUM over
  the f_leak bracket and the hop-fraction bracket. Two consequences
  are registered with it. (i) Under a TRUE NULL the floor trips by
  design and not by accident: with no split, s²(C) IS drawn from
  that H0 distribution, so the trip is the null's own signature and
  the registered route (clean dose, clean N0 ⇒ DATA ⇒ ¬D-sign
  partition ⇒ FALSIFIED reachable) is the route a true null is
  SUPPOSED to take. The floor is not an anomaly detector and must
  not be read as one. (ii) G1's freeze list gains the matching
  REACHABILITY check, since the formula alone cannot guarantee the
  floor sits below a real detection. THE CRITERION IS A RATE, NOT A
  RATIO, registered here in round 14, which found the requirement
  standing with no number at all: **P(s²(C) < floor | H1) ≤ 0.01 at
  BOTH bracket ends**, measured in the same ≥ 500-rep freeze that
  produces the floor. A ratio bar was the obvious candidate and is
  the wrong object: the deep-end ≥ 3× margin guards COUNTS against
  Poisson noise, where a ratio is the natural scale, whereas what
  this floor can do is trip on a true detection, and the size of
  THAT hazard is a false-trip rate under H1, which depends on the H1
  spread and not on the gap alone. The rate also puts the floor in
  the same shape as θ_D, whose H0 false rate is registered as a
  measured rate on held-out replicates rather than as a distance.
  The gap is still worth carrying as orientation: the frozen floor
  should fall below the predicted s²(C) at BOTH bracket ends,
  else the floor would trip on the very signal the flight is built
  to see, and the freeze moves the design rather than the number.
  MIND THE TWO CHAINS, the caution θ_F's base already carries (round
  12): E[s²(C)|H0] and σ̂ are ANALYZE-CHAIN quantities and carry the
  hop-noise floor, while the 0.00539 this check used to be written
  against is the NOISELESS dressed reference, so the two sides were
  not measured on one chain. The comparand is the H1 chain's own
  expected s²(C), which this document already carries in two pieces:
  s²(U) ≈ 0.0017 in the H1 chain (§3, the hop floor) plus the
  predicted difference d = 0.00482 (§8a), i.e. **≈ 0.0065**. Order of
  magnitude from this page's own σ̂ = θ_D/3 = 0.00084: a floor near
  0.0035-0.0042 leaves **1.55-1.86×**, on a threshold that since
  round 9 forces ¬D-sign. TWO SUBSTITUTIONS RIDE IN THAT
  ILLUSTRATION and round 14 named both, neither being visible in the
  numbers: the σ̂ used is θ_D/3, the spread of the DIFFERENCE d and
  not of s²(C), which is the conservative direction (σ̂(d) > σ̂(s²C))
  and is why the illustrative floor sits high; and 0.0065 is the
  NEAR bracket end (f_leak 8/15), while the far end gives
  0.0017 + 0.00554 = 0.00724, a larger comparand and a wider margin,
  so the end shown is the binding one. The freeze measures both ends
  and gates on the rate above, not on either gap.
  (Round 13 repair, and the failure is worth
  naming because it is subtle: round 12 stated the correct comparand
  in one sentence and then did the only arithmetic on the page
  against the number it had just disqualified, reporting ~1.3×. A
  freeze copying that example would have set the floor from a
  comparand ~20% low, in the conservative-LOOKING direction, on the
  threshold that can veto a detection.) Both producer inputs,
  E[s²(C)|H0] and its σ̂, are owed in §8a's refreeze package, where
  they come out of the same ≥ 500-rep H0 chain that produces θ_D,
  θ_W and the floor itself; a reachability requirement whose inputs
  appear nowhere on the page cannot be checked by the executor who
  has to act on it. The trip FORCES ¬W
  regardless of the computed exceedance (round 8: s² ≥ 0 makes a
  floored s²(C) yield some small d_W, and whether that clears θ_W
  would be a race between two independently frozen numbers;
  forcing removes the race, and the runner's floor evaluator is a
  machinery-ledger item). A FORCING NEVER OVERWRITES A VOID: if W's
  line is void by scope while the floor forces ¬W, the line stays
  VOID, since the forcing is a statement about a value and a void
  line has none (round 18; the verdict is unaffected, the printed
  line is not). It also forces ¬D-sign the same way
  (round 9: the identical race was still open on D-sign, since
  s²(C) < floor and d > θ_D are simultaneously satisfiable when
  s²(U) comes in anomalously low; a detection claim whose C-level
  is H0-indistinguishable from no-split is no detection, so the
  floor subordinates the point-estimate rule, "D-sign fails on its
  own terms" now BY REGISTRATION, and the verdict falls to the
  ¬D-sign partition exactly as the adjudication routes it). The floor
  and its frozen adjudication stay, REROUTED in round 3 (the prior
  routing sent a floor trip to "¬W class anomalous", a label whose
  own definition requires D-sign, which a floored s²(C) fails: the
  larger the falsification, the more certainly it would have read as
  instrument failure): a floor trip WITH dose certificates passing
  AND N0-CLEAN holding is DATA, never VOID; W reads ¬W, D-sign fails
  on its own terms, and the verdict falls to the ¬D-sign partition
  (FALSIFIED / Anti-D / INCONCLUSIVE) evaluated exactly as
  registered. Only a floor trip with a FAILING dose certificate or
  N0-CLEAN is VOID (there the instrument, not the room, is the
  suspect). The floor never converts a candidate falsification into
  instrument failure (§3's own rule, applied to itself); the former
  s²(C′) floor is retired (a low s²(C′) is
  expected, not pathological).
- **D-mag, two-sided agreement, BOTH corner magnitudes.** s²(C)
  inside its frozen band around f_C(B̂) AND s²(C′) inside its band
  around f_{C′}(B̂). s² is NOT additive in the generator
  (s²(Ĝ + B̂) carries cross terms unless B̂ is scalar), so both
  CENTRES are frozen FUNCTIONS committed as code with the same status
  as §9's width scaling functions, evaluated on the measured B̂
  (round-4 repair: a frozen NUMBER cannot exist pre-flight). The
  functions include the post-selection envelope from G1's counts-level
  sim (which biases absolute s² up to ~−20% at deep grids); a
  bare-corridor band would fire VIOLATED on a perfect device. G1
  verifies the reachability of the "D-sign ∧ W ∧ ¬D-mag" branch (the
  W exceedance and the D-mag band both read s²(C); the branch exists
  only if the band widths admit it), AND the C′ band's own
  DISCRIMINATION (the gate measures s²(C′) at signal-to-noise ~1 at
  budget, so an honest band around f_{C′}(B̂) risks containing any
  outcome, confirmation by imprecision; if the frozen band cannot
  separate the class alternatives at budget, the C′ containment
  reverts to REPORTED and leaves the CONFIRMED conjunction, by
  pre-data amendment at freeze). THE DEAD-C′ FORK, registered
  (round 4; §1 concedes the registered pair cannot distinguish the
  class law from a dead C′, so this is the only place that
  distinction lives): the freeze DECIDES it: with a discriminating
  C′ band, s²(C′) ≈ 0 fails D-mag(C′) and the flight reads "split
  confirmed, width off-prediction", the width law's falsification
  face; if the discrimination check fails and C′ reverts to
  REPORTED, s²(C′) ≈ 0 leaves CONFIRMED reachable under the §1
  honest scope; EITHER WAY the RECORD prints the raw s²(C′) beside
  the branch taken and names which regime the freeze selected.
  HONESTY, round 5: TODAY's gate numbers admit only the second
  regime (the frozen-config s²(C′) prediction is itself
  0.0002 ± 0.0002, indistinguishable from zero at budget, and no
  honest band can contain f_{C′}(B̂) while excluding 0); the first
  regime becomes available only if the refreeze's fixed-axis
  re-measurement separates f_{C′} from zero (the §8a
  reconciliation item). An executor should expect regime two. HOW
  REGIME TWO IS RECORDED, registered in round 12 because it was not
  representable: the freeze writes the SENTINEL value "REPORTED" into
  `dmag_band_Cp`, meaning the band is deliberately not registered by
  pre-data amendment and the C′ magnitude prints as a reading. Before
  that, the constants gate blocked on a missing entry, blocked on a
  null value and blocked on an unfrozen one, so the only way to reach
  submission in the regime this document tells the executor to EXPECT
  was to invent a number for a band the same paragraph says will not
  exist. Any evaluator consuming a sentinel key prints the reading and
  never a verdict. TWO keys are sentinel-eligible and no others,
  enumerated in code and enforced there: `dmag_band_Cp` and, with it,
  `f_Cp_dressing`, the C′ centre function whose only registered
  consumer is the band that stops existing. Round 13 closed the hole
  round 12's own repair had opened: the sentinel was DOCUMENTED and
  never implemented, so for one round the freeze gate accepted ANY
  non-null value on ANY key, and a string left in θ_D would have
  passed the last guard before a paid submission and then raised on
  the paid data. Numbers must now be numbers. PRECONDITION: the G5
  in-analysis J
  estimator inside its frozen pass band; a J-band failure VOIDs the
  magnitude/anchor arms per the SCOPE list below (D-sign and W stay
  evaluated with a named flag), never a silent widening.
- **A, equivalence test (the anchor; renamed from M in round 7
  to end the collision with the binding count M, G3's object).** The
  CIs on r̄(C) − r̄(U) and
  r̄(C′) − r̄(U) are CONTAINED in frozen equivalence margins (a TOST
  shape: "consistent with zero" alone would be confirmation by
  imprecision). CENTRES: the equivalence intervals are centred on the
  DRESSED A-centres frozen as functions beside f_C (the noiseless
  flown r̄ misses the ideal by −1.4%/−3.4%, §8a: estimator dressing,
  not error; zero-centred holds only in the compressed model);
  the model's own secular error is §4's PER-RATE figure, ≤ 1.2% down
  to Q = 10 (round-14 citation repair: this said "§4's 1-2%", which
  is neither of §4's two stated numbers; the split error there is
  +2.51%/−5.1%, and the per-rate one is what an r̄ margin carries),
  and the margin carries it.
- **In-circuit nulls, equivalence tests**, one per corner arm: the
  (1,2) row/column of Ĝ(C) − Ĝ(U), and the (2,3) row/column of
  Ĝ(C′) − Ĝ(U) (§3; C′'s wiring check, needing the full three-dyad
  inversion, so G2 freezes its margin SEPARATELY from the
  contamination-free (1,2) case). Both VOID-capable (wiring), while
  the decoupled-dyad RATES are report-only (physics).

**N0-CLEAN (frozen by G2; failing it is VOID):** B̂ physical (spectrum
in the admissible half-plane, magnitude envelope), fit residual bound,
and first-half vs second-half batch consistency of B̂ (the staircase
bracket shape). N0 is verdict-input data consumed ONLY through the
frozen committed functions (the dressings and f_C); no analyst-side
adjustment ever; the day-of re-gate (§9) is calibration-only and N0 is
not among its inputs.

**VOID (Class 2, §9).** Triggers, each with its freezing owner: N0-CLEAN
(G2), T1 CLEAN conditions (G4), dose certificates (G3), the (1,2) and (2,3) difference-nulls (G2,
separate margins), fit health: channel-fit residual, frequency-offset
excursion, or RATE SATURATION against the fit's own ceiling (the code
enforces all three, each behind its own freeze; the third is the one
§7 had left to §8a's refreeze list, round 12), plus
any NaN verdict statistic on a COMPLETE grid (a dead CHANNEL,
Amendment 1.2; a kept-count-floored CELL is NOT this trigger, it is
grid-incompleteness, single route, round-4 repair) reading as
fit-health VOID (G2),
the s²(C) floor under W AND under D-sign (rounds 8-9 gave it both
forcings; "under W" alone was the round-8 state, round-12 repair)
(G1; VOID only when dose certificates or
N0-CLEAN fail beside it, else it is the ¬D-sign DATA path, see the
W bullet's rerouted adjudication), band-validity window left (G1). The bank
combines POOLED statistics with INDIVIDUAL floors (the record-parity
correction of its own earlier shorthand); with the G5 J pass band added and Amendment 1's grid-incomplete
trigger (any science cell missing or binding-broken VOIDS D-sign and
W: the no-drop rule made executable) the bank holds TEN triggers, and
the union-bound false-VOID accounting is measured JOINTLY by G1 (the
verdict arms are built from the same FOUR fits, N0's B̂ entering every
dressed prediction, and are strongly correlated: family rates are
measured, never multiplied). SCOPE, registered for ALL TEN triggers
(round-4 completion; the consumption map that decides each scope:
D-sign consumes {C, U}, W {C, C′}, D-mag(C) {C, N0}, D-mag(C′)
{C′, N0}, A {C, C′, U}, each null its two arms):
- GLOBAL (voids every verdict arm): N0-CLEAN (B̂ enters every
  dressed centre), T1-CLEAN (round-5 re-scope: it is a DEVICE
  condition on the coupled chain, and decay on ANY site corrupts
  every arm's dynamics regardless of where engineered dephasing was
  injected, so a lighting predicate keyed to the injected profiles
  was the wrong key; the RECORD names the failing site(s) and
  face(s)), and the band-validity window (the design's
  validity, pre-verdict).
- GRID-INCOMPLETE: voids D-sign and W outright (Amendment 1.2);
  D-mag and A print as REPORTED EXCLUDED-CELL FITS (round-5
  precision: s² and r̄ are whole-grid channel fits with no per-cell
  version; on an incomplete grid they run on the surviving depths
  with the excluded cells NAMED, a deviation from the §5 no-drop
  rule that is permitted only for REPORTED content, which is
  exactly why the registered verdicts void instead of degrading).
- PER CONSUMING ARM (a trigger voids exactly the verdict arms that
  consume the failing arm's data, per the map above): fit health
  (per arm-pair as implemented: C/U unhealthy voids D-sign, C/C′
  voids W, each D-mag band its own arm, A its triple), dose
  certificates (the failing arm's consumers), the (1,2) null (a C wiring
  check: voids D-sign, W, D-mag(C), and A's C line), the (2,3) null
  (a C′ wiring check: voids W, D-mag(C′), and A's C′ line; D-SIGN
  STANDS, it never reads C′, the executor question answered).
- The s²(C) floor: per the rerouted adjudication above (DATA unless
  dose or N0-CLEAN fail beside it; when it IS void, scope = C's
  consumers).

PRECEDENCE AMONG THE TRIGGERS, registered here (round 11: the §7
precedence line ordered the ARMS, naming which failure leads the
verdict when several arms fail, but nothing ordered the ten triggers
INSIDE a VOID, so a flight where two fired had no determined line to
print, e.g. brake-truncation beside a (1,2)-null trip, or any global
trigger beside grid-incompleteness). **The line names the FIRST
firing trigger in this order and appends the others in the same
order:** (1) band-validity window, (2) N0-CLEAN, (3) T1-CLEAN,
(4) grid-incomplete (brake truncation included, carrying its
Amendment 1.6 label), (5) dose certificates, (6a) the (1,2)
difference-null, (6b) the (2,3) difference-null, (7) fit health,
(8) the s²(C) floor, (9) the G5 J pass band. The two nulls are
SEPARATE ranks (round 12: they shared one and a single layout
permutation fires both, leaving the line's name undetermined), and
(1,2) precedes (2,3) for the same reason the whole order is causal:
the (1,2) null voids D-sign itself while the (2,3) null leaves it
standing, so the one with the wider reach names the line. The order is causal, not alphabetical, and reads from the
most upstream cause outward: whether the design was valid at all,
then whether the background and the device were sound, then which
data exists, then whether the injected physics was realized, then
whether the wiring was right, then whether the estimator was, then
what level was read, and last the arm-scoped band. Two properties
worth stating, both consequences rather than extra rules: the s²(C)
floor can never LEAD a VOID line, since it is VOID-capable only
beside a dose or N0-CLEAN failure and both outrank it, which is what
makes its own VOID branch a second name for an already named state
rather than a competing one; and the J pass band can never lead
either, since it never voids D-sign or W at all. Appending, never
replacing, keeps the RECORD's diagnosis complete: precedence decides
what the line is CALLED, never what is written down. (The SCOPE list
above resumes below with its fifth and last member, the G5 J pass
band; this paragraph sits between them because precedence is what
that last member's resolution depends on.)
- The G5 J pass band: RESOLVED round 4 (it had been both a D-mag
  precondition and an unscoped bank member): a J-band failure VOIDS
  the magnitude/anchor arms (D-mag both, A), whose centres assume
  the frozen J, while D-sign and W stay EVALUATED and carry a named
  J-band flag in the RECORD (the even +4.5% J-inflation is
  multiplicative on a real split and creates none under H0,
  measured on the §6.3 zero-spurious construction and quantified in
  §8a G5; §9's bank lists the trigger with this scope).

HOW A VOID ARM ENTERS THE VERDICT TABLE (round-5 rule, round-6
completion). A VOID in one scope leaves every other arm computed and
reported; the LINE is then composed by the algorithm below, which
places each void arm in its own ARM position, and the flight-level
VOID lines of step 3 are the one place where triggers are appended,
in the form step 3 states.
Round 19 cut this paragraph back to that. It had carried its own
ordering rule through four rewrites: round 12 repaired its pointer,
round 17 repaired the pointer again, round 18 scoped the rule to step
3 without noticing that step 3 has no surviving verdict to "lead
with", and each of those repairs left the rule itself standing beside
an algorithm that ordered the line differently. Two live orderings
for one line is what a reader had to reconcile; now there is one, and
it lives in the algorithm.

**THE COMPOSITION ALGORITHM, stated once and normative (round 17),
AND BUILT AS CODE WITH TESTS (2026-08-17, after round 19).** The
algorithm below is the specification; `corner_beat_verdict.py` beside
the runner is its executable form, and `test_corner_beat_verdict.py`
is the specification's test suite, each test naming the round whose
finding it encodes. It walks EXHAUSTIVELY: all 972 assignments of
TRUE / FALSE / VOID across the arms' lines, in both D-mag regimes,
asserting that every state yields exactly one line and that no line
leaves the five-form vocabulary. The reason for building it is a
measurement, not a preference: rounds 12 through 19 each read these
rules and each found its defects in the previous round's PROSE
repair, while finding none in the physics, so the state machine was
moved somewhere it can be walked instead of re-read. Where prose and
code disagree, THIS SECTION governs and the code is the defect; the
tests exist so that such a disagreement is visible rather than
discovered by an executor at three in the morning.
Rounds 13 through 16 each patched the verdict-line rules and each
patch produced the next round's findings, all of them landing on one
ordinary cell; the rules are therefore written here as a procedure
that generates every line mechanically, and the registered label list
below is its OUTPUT, not a second source of truth. The composer:
1. Evaluates every arm LINE to TRUE, FALSE or VOID. An arm's LINE
   SET is read from the frozen constants, not assumed: D-sign and W
   have one line each; A has two; and D-mag has two only when
   `dmag_band_Cp` carries a band, since under the reverted-band
   freeze this document tells the executor to EXPECT, that key holds
   the REPORTED sentinel and D-mag is the C line alone (round 18).
   D-sign, being single-line, can never be PARTIAL.
2. Classes each arm: VOID if all its lines are void; PARTIAL if some
   but not all are; otherwise VALUED, and a valued arm holds iff all
   of its non-void lines hold. PARTIAL is a statement about VOIDNESS
   only. An arm whose lines all valued and one of which FAILS is
   ¬arm, never PARTIAL.
3. If a GLOBAL trigger fired, or a trigger voided D-sign itself:
   flight-level VOID (<leading trigger by the precedence order
   above>; <every other firing trigger, in that same order>), every
   surviving arm REPORTED. Step 3 is the one place the appended-void
   form applies, and it appends ALL of them (round 18: the step named
   only the leading trigger while the precedence rule and Amendment
   1.6's worked case both append the rest).
4. Else if D-sign holds: the line enumerates the arms in the
   reporting order D-sign, W, D-mag, A, each as "D-sign detected" for
   D-sign and "<arm> holds" for the other three when valued true,
   "¬<arm>" when valued false, optionally carrying its lines'
   verdicts as "¬<arm> (<line> fails; <line> holds)" on a two-line
   arm, "<arm> PARTIAL (<surviving line> holds|fails; <line> VOID
   (<trigger>))" when partial, and "<arm> VOID (<trigger>)" when
   void. These five forms are the whole vocabulary; any other
   spelling in this document is a defect (round 19, which found the
   ¬-with-lines form used in a registered example and listed
   nowhere, and the PARTIAL form registered in three places of which
   round 18 updated one). The PARTIAL form names the surviving
   line's TRUTH VALUE (round 18: without it, an A that is partial with
   its C line holding and one with its C line failing printed the
   same string, and that difference is the entire content of the
   anchor arm). ONLY IF ALL FOUR ARMS ARE VALUED does that
   conjunction carry a registered NAME (CONFIRMED or one of the three
   qualified split lines); otherwise the enumeration IS the line and
   no named line is claimed. This is what "no partial pattern is ever
   promoted" means operationally.
5. Else (D-sign valued false): the ¬D-sign partition, under
   preconditions (i)-(iii), κ first and then the CI, as registered
   below. (iii) fails in TWO distinct ways and they route
   differently: if A's C line FAILS, the line is "anchor failed, no
   split observed", instrument-suspect INCONCLUSIVE; if A's C line is
   VOID, the line is "NOT DETECTED; partition unavailable (A VOID
   (<trigger>))", because a void line is neither true nor false and
   cannot satisfy ¬A at all. Precondition (ii) needs no branch of its
   own, and the reason is (iii): any trigger that voids A's C line
   fails (iii) and is routed there. Round-19 repair of a round-18
   sentence that gave the reason as "every trigger that touches A's C
   line also touches D-sign's consumers", which the G5 J pass band
   falsifies in one lookup, since it voids the magnitude and anchor
   arms while D-sign and W stay evaluated with a named flag.
   THE LINE FORMAT IS STEP 4's, with the partition label in front.
   The label carries D-sign, so the enumeration starts at W and runs
   W, D-mag, A in step 4's forms (round 18 gave step 5 a format,
   having had a label and none, so six of the ten registered forks
   had a verdict an executor could name and not write down). For a
   (2,3)-null trip under a true null, in the expected reverted-band
   freeze: "FALSIFIED; W VOID ((2,3) null); D-mag holds; A PARTIAL (C
   line holds; C′ line VOID ((2,3) null))"; should the freeze keep a
   C′ band, D-mag is PARTIAL beside A. Round-19 repair: round 18 wrote
   this example with D-mag PARTIAL in the very regime where D-mag has
   no C′ line, seven lines above the paragraph in which it diagnosed
   exactly that error in the neighbouring example.
A worked example, the cell that four rounds kept mislabelling, a
(2,3)-null trip under a strong positive, IN THE EXPECTED
REVERTED-BAND FREEZE where D-mag is the C line alone: the null voids
W and A's C′ line, so W is VOID, A is PARTIAL and D-mag is untouched,
and the line reads "D-sign detected; W VOID ((2,3) null); D-mag
holds; A PARTIAL (C line holds; C′ line VOID ((2,3) null))" with no
named line claimed, because not all four arms are valued. Should the
freeze instead keep a C′ band, D-mag has two lines and is PARTIAL
alongside A (round 18: the earlier version of this example printed
D-mag PARTIAL in the regime where D-mag has no C′ line at all, and
applied the sub-arm rule to A while leaving D-mag out of it). The
other cell, a true null where A's C′ line FAILS while its C line
holds: A is ¬A, not PARTIAL, precondition (iii) is satisfied, the
partition is entered, and the line is the partition label followed by
the arms, "FALSIFIED; W ...; D-mag ...; ¬A (C′ line fails; C line
holds)". CONFIRMED and every CONFIRMED-TYPE named
conjunction line (defined here, round 15, having been load-bearing
and inferred: a registered line whose statement requires at least one
verdict arm to hold as TRUE, i.e. CONFIRMED and the three qualified
split lines; the ¬D-sign lines are not of this type, since they
require arms to FAIL) require ALL their arms non-void, and no partial
pattern is ever promoted to them; FALSIFIED follows the ¬D-sign
partition's OWN registered precondition list, whose A-requirement
is A's C line (the SUB-ARM carve-out below), the ONE registered
exception to the all-non-void rule (round-9 reconciliation: the
two paragraphs had contradicted each other on exactly the
(2,3)-null-under-true-null walk). None of this is "Any unlisted
pattern" INCONCLUSIVE (that
catch-all is for non-void patterns only). THE STATES THAT
NEED THEIR OWN LABEL (rounds 6-8): (i) a GLOBAL trigger (N0-CLEAN,
T1-CLEAN, band-validity) leaves no non-void arm, and the verdict IS
the flight-level "VOID (<trigger>)", the precedence line's first
word; (ii) a ¬D-sign pattern whose partition PRECONDITION is void
(e.g. A under a J-band failure) prints the registered label
"NOT DETECTED; partition unavailable (<arm> VOID (<trigger>))",
which is the existing failing-guard fallback given its printable
form; (iii) grid-incomplete prints "VOID (grid incomplete)" with
the REPORTED excluded-cell fits beside it, no conjunction needed;
(iv) a trigger that voids D-SIGN ITSELF (any per-consuming-arm
trigger touching {C, U}: the (1,2) null, a C- or U-arm dose
certificate, C/U fit-health) prints the flight-level
"VOID (<trigger>)" with every surviving non-void or partial
sub-arm REPORTED beside it, round 8: with D-sign void no named
line exists, and the catch-all stays closed, so this class needs
its own label exactly as the global one does.
SUB-ARM RESOLUTION (round 6; A and D-mag are two-line arms): A
holds iff BOTH its lines hold; A is VOID only if BOTH lines are
void; an arm with ONE void line and one non-void line is PARTIAL
(round 7): a PARTIAL arm never satisfies a conjunction (any
CONFIRMED-type line requiring it is unreachable) and prints in step 4's PARTIAL form (which names the surviving
line's truth value; the form is registered ONCE, in the algorithm); a
PARTIAL arm whose SURVIVING line FAILS counts as ¬<arm> wherever a
¬-condition is being EVALUATED, which after round 17 means inside the
¬D-sign partition's preconditions and nowhere else. It does NOT carry
a partial pattern into a NAMED conjunction line. **ROUND-18 REPAIR OF
A BLOCKER, and of round 17's own structural fix:** round 10's clause
used to read "D-sign ∧ A-partial-with-failing-C-line reads 'Split
observed, anchor failed' with the partiality named beside it, closing
the last unlabeled pattern". That IS one of the three qualified split
lines, so on a perfectly reachable state (a C′-arm dose-certificate
failure voids A's C′ line without touching D-sign's consumers {C, U},
and r̄(C) − r̄(U) then falls outside its margin) step 4 of the
algorithm emitted an unnamed enumeration while this clause emitted a
registered name: round 10's patch surviving underneath round 17's
algorithm, which is the arc's whole pattern in one sentence. The
algorithm governs, by its own declaration of being the single source,
and the pattern is not left unlabeled: under step 4 the enumeration
IS its label, "D-sign detected; W holds; D-mag holds; A PARTIAL (C
line fails; C′ line VOID (dose certificate))". The
one registered exception is the ¬D-sign partition, where "A holds"
means specifically A's
C LINE (the partition concerns D-sign's own pair, and the C′
line's state is reported beside it), so a (2,3)-null trip, which
voids A's C′ line and D-mag(C′), leaves the ¬D-sign partition
REACHABLE on a true null while CONFIRMED is not; under a dead-C′
freeze (C′ containment reverted to REPORTED), D-mag means
D-mag(C) alone. RECONCILIATION with grid-incompleteness (round 6:
why a C′-only HOLE voids D-sign although the (2,3) null does not):
the null is a wiring check on data D-sign never reads, but
grid-incompleteness breaks the NO-DROP guarantee θ_D was frozen
under (the H0/H1 statistics were measured on complete
four-arm grids), so a hole anywhere invalidates the frozen
thresholds themselves, not the arm's data.

**Verdicts.** ALL arms are always computed (round-4 repair: a
short-circuit would make every A-dependent branch below unreachable);
the order VOID → D-sign → W → D-mag → A is REPORTING precedence only,
naming which failure leads the verdict when several co-occur:

- CONFIRMED = D-sign ∧ W ∧ D-mag ∧ A.
- "Split confirmed, width off-prediction" = D-sign ∧ W ∧ ¬D-mag ∧ A
  (this branch IS the quantitative width law failing while the split
  stands: the width law's falsification face, and the RECORD names
  it as such, round-3 repair; the registered pair alone cannot
  falsify the width numbers, §1).
- "Split confirmed, class anomalous" = D-sign ∧ ¬W ∧ A (follow-up: the
  realized-profile assertion records which classes flew).
- "Split observed, anchor failed" = D-sign ∧ ¬A (any W/D-mag);
  the symmetric ¬D-sign ∧ ¬A **where A's C LINE fails** (the scope
  amended here in round 17, in place: round 16 rescoped it seventy
  lines below and left this bullet reading whole-arm ¬A, which is the
  document's own named failure mode, a repair landing in one of two
  places; a failing C′ line alone leaves the partition reachable, per
  the composition algorithm above) is named "anchor failed, no split
  observed" and likewise reports the failing anchor:
  instrument-suspect INCONCLUSIVE; follow-up: the dose or the locus
  condition did not realize; the dose certificates and N0 record decide
  which.
- ¬D-sign partitions (disjoint by construction, all three under the
  same preconditions, SCOPED in round 7 so the sub-arm rule and this
  list read as one: (i) no GLOBAL void and grid complete, (ii) every
  Class-2 trigger whose SCOPE (the §7 list) touches D-sign's
  consumers {C, U} or A's C line has not VOIDED (round-10 word
  repair, was "untripped": a tripped trigger whose adjudication
  routes to the DATA path, the s²(C) floor with clean dose and
  N0-CLEAN, does NOT bar the partition, that routing being its
  whole point; "guards clean" means
  exactly this, never the whole bank: a (2,3)-null trip voids
  W/D-mag(C′)/A's C′ line and leaves this partition REACHABLE),
  (iii) A's C line holds and is non-void. The power condition
  SE[s²(C) − s²(U)] ≤ κ × G1-projected, κ frozen, on that NAMED
  statistic, is **NOT a fourth precondition but the FIRST
  DISCRIMINATOR INSIDE the partition**, evaluated before CI
  containment (round-14 repair: it had been listed as a precondition
  while round 13's label split simultaneously made it the
  discriminator of INCONCLUSIVE (underpowered), which is
  self-cancelling, since a precondition that fails means the
  partition is never entered and the label inside it is unreachable.
  Reading it as a discriminator is the only version in which every
  ¬D-sign outcome has exactly one name). So: power fails ⇒
  INCONCLUSIVE (underpowered), full stop, whatever the CI does; power
  holds ⇒ the CI decides among the three below):
  **FALSIFIED** = the CI on s²(C) − s²(U) is contained in
  [−θ_F, θ_F], where θ_F is a separate PHYSICS-SIZED equivalence
  margin, REGISTERED here as the fraction AND the base (round-3
  repair, base pinned in round 4: a
  "producer plus reachability check" was a researcher degree of
  freedom on the falsification criterion itself, and an unregistered
  base re-opens it): **θ_F = 0.75 × the DRESSED PREDICTED DIFFERENCE
  f-composed at the frozen config** (the same quantity the D-mag
  centre chain predicts for d; today's gate value 0.00538,
  re-frozen on the fixed axis at the refreeze, never the raw H1
  sample mean), G1's role reduced
  to REPORTING the resolving power, never tuning fraction or base.
  The arithmetic is on the page: the 95% CI half-width at the frozen
  budget is 1.96 × 0.00165 ≈ 0.0032 = 60% of the 0.00538 base
  (θ_F ≈ 0.0040), so FALSIFIED at this budget can only ever mean
  "the effect is smaller than ~3/4 of prediction", and the RECORD
  states that scope WITH the base's assumption named (round 9: the
  base is the NOISELESS dressed difference, dressed C 0.00539 −
  dressed U 0.00001; against the H1 chain's own expected
  d = 0.00482, which carries the measured U hop floor, the same
  θ_F is ~5/6, and the refreeze re-states the fraction against the
  re-frozen base as part of the §8a reconciliation). (Round-5 history: containment in [−θ_D, θ_D] is
  scale-invariantly impossible when θ_D is the detection critical
  value ~1.5 SE against a ~3.9-SE-wide CI); **Anti-D** = the CI lies
  entirely below −θ_F, recorded
  as falsified-with-inverted-sign; and TWO DISTINCT INCONCLUSIVE
  labels, split in round 13 because one slash-joined string served two
  disjoint causes and the composer had no rule for which half to print
  (scenario: a floor trip with d > θ_D and the power condition
  satisfied): **INCONCLUSIVE (underpowered)** = the κ power condition
  fails, and **INCONCLUSIVE (indeterminate)** = the power condition
  holds but the CI sits in no registered partition (straddling ±θ_F,
  or lying entirely above θ_F, reachable when θ_F < θ_D and also,
  round 10, via a floor trip with large d now that the floor forces
  ¬D-sign). Both name their cause in the RECORD (and any ¬D-sign case
  whose preconditions (i)-(iii) fail reports the failing guard/anchor
  instead). THE THIRD INCONCLUSIVE ON THIS PAGE IS A DIFFERENT
  OBJECT and is not a third partition member: "instrument-suspect
  INCONCLUSIVE" above belongs to the ¬D-sign ∧ ¬A conjunction and
  reports the failing anchor. ROUND-16 SCOPE REPAIR, because the
  sentence that used to reconcile the two read A as a WHOLE ARM and
  round 9 had already narrowed precondition (iii) to A's C LINE: ¬A
  is therefore NOT ¬(iii), and the state where A's C line holds while
  its C′ line fails satisfies (iii) and also satisfies ¬D-sign ∧ ¬A,
  so two registered strings collided on one perfectly ordinary cell,
  the one a true null is most likely to produce. The composition
  algorithm above settles it and this paragraph no longer carries its
  own rule (round 17, after the round-16 patch produced five findings
  of its own): A with a holding C line and a FAILING C′ line is ¬A
  and not PARTIAL, since PARTIAL is about voidness; precondition
  (iii) asks only about A's C line, so it is satisfied and the
  partition is entered; and "anchor failed, no split observed" is
  reserved for a FAILURE of A's C line, while a VOID A's C line is
  the partition-unavailable state instead. Registered labels, whole list, so the composer
  has a closed set: CONFIRMED; the three qualified split lines;
  "anchor failed, no split observed", the ¬D-sign ∧ ¬A line, whose
  label IS instrument-suspect INCONCLUSIVE (round 15: the two were
  listed as separate members, leaving the composer two strings for
  one state, which is the defect round 13 fixed for the slash-joined
  INCONCLUSIVE reappearing one level up; they are one entry, a named
  line and its label); FALSIFIED; Anti-D; INCONCLUSIVE
  (underpowered); INCONCLUSIVE (indeterminate); VOID (<trigger>);
  NOT DETECTED; partition unavailable (<arm> VOID (<trigger>)); and
  <arm> PARTIAL (per step 4's form); and the ENUMERATED line
  that step 4 of the composition algorithm produces whenever not all
  four arms are valued, which carries no registered name by
  construction. That enumerated form is not a label to be listed
  separately, it is what the algorithm emits, and round 16's attempt
  to register it as one more string is what produced this round's
  findings: the algorithm above is the single source, this list is
  its output. The catch-all "any unlisted pattern: INCONCLUSIVE,
  named" is not a label at all but the composer's fallback
  obligation: if it ever fires, the pattern that reached it is a gap
  in the algorithm and belongs in the RECORD as such.
- Any unlisted pattern: INCONCLUSIVE, named, with follow-up stated.

**On CONFIRMED:** Confirmations registration scoped to D-sign + W only
(the concentrator precedent: magnitude arms are not registered); the
entry records that D-mag and A held as CONFIRMED's other conjuncts and
that W's predicted centre was the B̂-dressed function, since the
registry stores predicted-vs-measured.

## 8. Gate items before bands freeze (the committed gate: `simulations/corner_beat_gate.py`)

- G1: counts-level end-to-end sim of the flown construction (shots,
  readout confusion, gate-error as depolarizing per 2q gate including
  its within-sector component, the post-selection envelope, M-finite
  channel realization, var(B̂) propagation) through the §6 estimator;
  joint choice of (Q, J·dt, depth SPACING, shots, gate type), sixteen
  configurations, two-stage on the verdict metric. STATED EXACTLY,
  because this sentence used to claim the "depth grid endpoint" as an
  axis and the committed code never varied it (2026-08-17, read at the
  source): the endpoint is 60 Trotter steps in EVERY searched
  configuration, and the only grid freedom exercised is the spacing k,
  which ties the point count to the endpoint through one variable
  (13 × 5 vs 21 × 3). The endpoint is therefore an UNSEARCHED axis, not
  a settled one, and §12 carries it as a follow-up rather than this
  bullet carrying it as done;
  freezes: all bands, θ_D, θ_F's EVALUATION at the §7-registered
  fraction × base (G1 reports the resolving power only, round-5
  ownership alignment: §7 abolished the producer-plus-reachability
  construction), the "W ∧ ¬D-mag" branch reachability AND the "¬W ∧ D-mag" one
  (both magnitudes in band while the exceedance fails; named verdict
  outcomes must exist before they can be read, round-4 addition), the
  s²(C) floor, the A-arm equivalence margins (both layers: the statistical
  part and the §4 secular physics term the margin carries; the null
  margins are G2's), the projected SE of s²(C) − s²(U) at the frozen
  budget (the κ condition's reference), the frozen
  centre-functions f_C(B̂) AND f_{C′}(B̂) and the dressing
  compositions, κ, the family error rates (P(≥1 VOID),
  P(false CONFIRMED) under H0, measured jointly), the frozen
  depth-point list, the band-validity window (dimensions and its
  relation to the Class-1 line rule), the day-of band-WIDTH scaling
  functions (§9's re-gate has no other source), and
  the centred-third-invariant report band (§6.4) (the billing CAP
  itself is NOT G1's: it is Tom-frozen at 25.0, the accepted band's
  top, Amendment 1.6 and the manifest); and the H0
  definition: H0 = the class physics absent (the committed gate
  implements the D-sign H0 by flying the C arm with the U profile
  while analyzing it with the C estimator, C′ untouched, which is
  sufficient for θ_D; the BOTH-corner-arms-null H0 that the W false
  rate and P(false CONFIRMED) require is part of the refreeze's
  family-rate measurement), all arms dressed by the same measured B̂
  (the s²(B̂) contribution then cancels in the D-sign difference in
  expectation, and G1 measures its variance contribution).
- G2: estimator well-posedness at counts level: the N0-dressing
  residual calibration and N0-CLEAN thresholds; post-selection under
  asymmetric readout confusion; the U-arm demodulation treatment
  (§6.3, pre-decided); the (1,2)- and (2,3)-difference-null margins
  (separately: the (2,3) one rides the full inversion); fit-health
  thresholds; invariants-vs-sorted bias demonstration; where the
  linear chain becomes nonlinear (ANSWERED by Amendment 1.1: the
  single post-pooling clip).
- G3: dose certificates: per-site retention table AND the two-site
  realized retention |E[e^{−i(φ_l−φ_m)}]| over the M bindings for a
  lit/lit, a lit/dark, and a dark/dark pair (the two-site face has
  never been flown; a correlated-phase bug is invisible per-site and
  reads as a rate error); pass criteria with provenance (concentrator:
  every step within 0.02, a measured certify OUTPUT there, adopted
  here as a provisional criterion, round-3 label; artifact < 10%); M frozen from the
  realization error measured on s² itself, with the per-binding billing
  overhead question BOUNDED by measurement 2026-08-17 (§10: no gross
  per-circuit surcharge at M = 256; not separable from that job's own
  delay cost, M = 1024 unmeasured until flight, brake is the
  backstop); the exactly-one-RZ post-transpile assertion.
- G4: the scattered-background budget, BOTH arms: a site-scattered
  profile (T1 as Γ_l/4-equivalent dephasing, gate-error scatter)
  splits the degenerate U arm from zero AND moves s²(C) at first order
  (the v1-gate per-arm figures, superseded by G4 v2's
  difference-quantity bands in §8a, kept only as provenance: U
  spurious ~12% of s(C), s(C) 1σ ±7%, both at the 0.25γ̄ scatter
  row). G4 additionally quantifies the same scatter on s²(C′) and on the
  W exceedance (v2 measures the D-sign difference; the C′-level
  quantities remain for the committed gate). G4 sets the T1 CLEAN conditions (staircase shape:
  named, individually stated, in-situ, sampling the arms' own time
  window) and the scattered-background band widening. ROUND-3
  SHARPENING, both directions: (a) the first-order background
  sensitivity is the PROJECTION of the background profile onto C's
  mirror-odd part, not its scale (measured: a mirror-symmetric
  background of ANY size leaves d exactly unchanged, a C-aligned
  0.25γ̄ background moves d ±25%, a C′-aligned one 0%), and the aux
  stations measure the per-site Γ_l, i.e. exactly that alignment, so
  G4 gates on the measured projection and the blanket ±28% band
  narrows to the aligned component; (b) honestly, the T1-CLEAN scale
  bound is close to non-binding as scatter (T1 scatter moves d by
  < 0.2% in the round-3 propagation), so its real load is the
  staircase TELEGRAPH shape (time structure), which four stations
  sample coarsely, N0's halves-consistency being the only continuous
  check; the RECORD says which of the two did the work.
- G5: J-disorder and J-calibration: detuning-aware treatment under σ_J
  (priced; one frozen draw = a systematic widening, the f129 b_qs
  shape); the ∂s²/∂(J·dt) band term (§4); the in-analysis J estimator
  from the fitted dyad frequencies (§4) with its frozen pass band as a
  precondition on the magnitude/anchor arms (D-mag both and A; §7's
  SCOPE list is the one registered scope). ROUND-3 ADDITIONS: (a) a Δ_eff / residual
  static-ZZ line: the room's threefold degeneracy is a Δ = 1
  structure, DD is deliberately off (§5), residual ZZ of 5-20 kHz
  against the implied J/2π ≈ 34-52 kHz gives Δ_eff errors of 2-15%,
  and the round-3 propagation measured dΔ = ±5% moving d by
  ∓10%/+9%, a SIGNED systematic (it does not average out like J
  scatter), so G5 budgets it with its sign convention stated; (b)
  the DYAD-FREQUENCY AGREEMENT band, promoted to a precondition on
  D-mag exactly like the J pass band: the three room dyads being
  degenerate is the compression's validity condition, the fitted δω
  already measure it (spread 0.0066 at dΔ = 0.01, well inside the
  fit's range), and the frozen band doubles as (a)'s in-data guard.
- G7: the runner gate (the concentrator 7b lesson): runner built,
  per-binding raw counts persisted before any reduction, `--analyze`
  proven end-to-end on a hardware-shaped synthetic artifact, the
  shared-skeleton and one-RZ-per-site-per-step and DD-off assertions in
  code, then an empty review of runner + gate records as its own stage.
  **STATUS (2026-08-17): DONE through four empty rounds.**
  `run_corner_beat.py` (external pipeline repo) with modes
  certify/calibrate/aer/hardware/analyze; certify passes at machine
  precision (sector parity 7e-16, trajectory parity vs the committed
  gate simulator 5e-14, fit-search parity vs the committed gate exactly
  0.0, dose certificates over all swept cells); analyze proven
  end-to-end at M = 64 and at the frozen M = 1024, the latter now on a
  FINAL-CODE full-M artifact (2026-08-17, all 254 PUBs, ~11.5 h of Aer;
  sha and demodulation-basis audits OK, grid complete, per-binding
  structure intact, and the synthetic realization lands on the gate's
  frozen expectations: d = +0.00484 against the predicted
  +0.00482 ± 0.00165 and d_W = +0.00532 against +0.00634 ± 0.00134.
  That is a PLUMBING result on synthetic data, not evidence about the
  room); every Class-1 guard
  a hard abort surviving python -O; submission blocked by a complete
  §8 constants manifest until the gate work freezes it. The rounds'
  load-bearing catches and the chain decisions are recorded in
  Amendment 1 and in the sharpened outstanding list (§8a).

(G-numbering note: G6 existed in v2 as the ω = −4 second demodulation
and is deleted, refuted in round 2; the number is retired, not reused.)

## 8a. Gate runs (2026-08-16, `simulations/corner_beat_gate.py` v2.1 after two promotion review rounds; outputs in `simulations/results/corner_beat/`)

The gate's own review round found and the v2 re-run fixed: a
winner's-curse selection (now two-stage on the VERDICT metric
P(d > θ_D), held-out seed, both f_leak ends), a gate-error bracket with
no bias channel (now a stochastic XX/YY magnon-hop channel inside the
evolution plus the Z-like half at the gate-count profile), a G3 that
never measured what §8 freezes M from (now the s²-realization spread
itself), a G4 on the wrong quantity (now the difference under a common
background), a G5 whose preparation basis followed the disorder (now
nominal-J basis, the flight's situation), missing readout confusion
(now symmetric 1.5% per qubit in the counts layer), and a
√w-weighting bug (~3% on the centre).

- **G1 v2.1 (two-stage, GAUGE-PINNED estimator; the pre-pin v2
  numbers rode a J·dt-unstable mode-sign gauge and are superseded;
  terminology: "power" throughout this section is the z-score d/SE,
  while the detection PROBABILITY is always written P(detect) or
  P(d > θ_D), round-3 disambiguation):**
  stage-1, worst f_leak end: fractional beats standard in every cell;
  no standard-gate config reaches worst-end power 2 at 11 min.
  Stage-2 (200 reps, held-out seed), the frozen 21 × 3 grid, 8192
  shots, fractional: **d = +0.0049 ± 0.0017, power 2.9, θ_D = 0.0035,
  P(d > θ_D) = 0.815, H0 false rate 0.000 in-sample AND held-out,
  ~11 QPU min** for the 8192-shot variant. The C′ lines land where
  the promotion round predicted: s²(C′) = +0.0002 ± 0.0002, BELOW the
  hop-noise U floor (s²(C′) − s²(U) = −0.0015), while the W exceedance
  is robust. **THE FROZEN CONFIG (16384 shots, Tom's second decision),
  at 200 H1 reps per end and 100 H0 reps split 50/50 (θ_D = 3σ̂ of
  the first 50, the false rate measured on the held-out 50: "0.000"
  is 0/50 with a ~6% binomial 95% bound, "0.020" is 1/50; the
  ≥ 500-rep refreeze below is what makes these verdict-grade), both
  f_leak ends:
  worst end (f_leak = 8/15) d = +0.00482 ± 0.00165 (unrounded, so
  the division reproduces the power, round-4 note), power 2.9,
  θ_D = 0.00253 (the unrounded transcription of the gate output, 3σ̂
  over 50 H0 replicates, σ̂ = np.nanstd with ddof = 1; the
  registered threshold, §7),
  P(d > θ_D) = 0.920, H0 false 0.000; other end (f_leak = 0.9)
  d = +0.00554 ± 0.00133 (the §7 far-end margin arithmetic's
  source, round-4 traceability),
  power 4.2, P(detect) = 1.000, H0 false 0.020, both AT ITS OWN 3σ
  value 0.00063, which the flight does not use (§7; the re-evaluation
  at the registered threshold is part of the refreeze);
  W exceedance power 4.7 / 3.8.** p2 sensitivity at the frozen grid:
  power 2.7 at the recorded 0.5%/2q, 1.1 at 0.75%, dead at 1% (the
  Class-1 guard input; measured at the frozen 16384 shots, 40 reps;
  the 2.7 here vs the stage-2 row's 2.9 at the same nominal p2 is
  the 40-vs-200-rep spread, reconciled at the ≥ 500-rep refreeze.
  SINCE AMENDMENT 1.8 THE CLASS-1 BOUND IS THIS LADDER'S FIRST RUNG:
  the guard sits at 0.5%, where the ladder measures rather than
  interpolates, which is why the interpolated power that used to
  stand in the header and the measured point the refreeze owed at
  0.6% are both gone. Read the rung honestly: it is the 40-rep
  number, and the registered header power 2.9 is the 200-rep row at
  the same nominal p2, one reconciliation, not two ladders).
- **G3 v2:** retention pass probability over 40 table draws: M = 256
  passes 20%, M = 512 75%, M = 1024 100% (criterion 0.02, transplanted
  from the concentrator at a different dose and noted as such); the
  quantity §8 freezes M from, the s²-realization spread across frozen
  channels (re-run under the pinned estimator): **31% of centre at
  M = 256, 15.3% at M = 1024, 8.1% at M = 4096**. M = 1024 holds, and
  the 15% enters the D-mag band as a systematic (M = 4096 would halve
  it at 4 shots/binding; a G3 freeze-time trade).
- **G4 v2 (the verdict quantity, common scattered background):** the
  95% band on Δ(s²C − s²U) is 14% of the base difference at scatter
  scale 0.125γ̄, 28% at 0.25γ̄, 58% at 0.5γ̄, 111% at 1.0γ̄. **Proposed
  T1-CLEAN bound: scattered Γ_l/4-equivalent profile scale ≤ 0.25γ̄**
  (band under half the base).
- **G5 v2.1 (estimator level, nominal basis, pinned gauge):**
  σ_J = 1% costs ≤ 10% of s²(C) (95%), σ_J = 2% ≤ 25%; and the
  J-calibration response is **EVEN under the pinned estimator:
  ∂s²/s² = +4.5% for BOTH ±5% J** (any J mis-calibration INFLATES s²,
  a one-sided bias toward false D-sign, which is exactly why the G5
  J pass band is a precondition and enters θ_D's H0 model). The v6
  draft's compressed-level −27% was an artifact of letting the basis
  follow the disorder; the pre-pin ∓4% antisymmetric form was a
  gauge artifact.

Still outstanding for the committed gate (SHARPENED by the four G7
review rounds, 2026-08-17): **the θ_D refreeze package**, one unit:
re-freeze θ_D through the ANALYZE-SIDE chain (the runner's fit_rate,
whose dead-cell NaN differs from the gate's silent 0; pooled
inversion + single clip; fixed R_MAX_FIT; NESTED tables; the
REALIZED-DOSE time axis t_eff = max(n−1, 0)·dt in gate and runner
together, the round-3 off-by-one repair in §5, whose effect this
paragraph used to state as "restores ~16.6% of s²" and which the
refreeze has now MEASURED through this very chain, with the opposite
sign for the verdict (2026-08-17, 500 reps per corner): s²(C) rises
8.7%, not 16.6%, at the worst corner; d rises only 2.4%, because
s²(U) rises with it and the arms' difference is what the verdict
reads; and the H0 threshold rises 14% while the H1 spread rises 11%,
so P(detect) at the worst corner falls 0.934 → 0.900 and the power
3.12 → 2.87. THE REGISTERED AXIS CORRECTION COSTS DETECTION, IT DOES
NOT BUY IT. That is not a reason to withdraw it, since the nominal
axis fits a model the flown circuit does not realize, but the
executor must not read the switch as a gain. Round 10's own caution
was right and this is its measurement: the −4 to −6%, the 16.6% and
the −1.4%/−3.4% are three numbers from three constructions, and only
the fixed-axis re-fit through the analyze chain measures what the
verdict sees; the H0 mean E[d|H0] stated and θ_D defined as
E[d|H0] + 3σ̂ if that mean is nonzero, the two H0 sides using
different estimator bases; and the hop-fraction BRACKET: the gate
pins the within-sector error split 50/50 while s²(U)'s floor moves
0.0014 → 0.0003 across hop fractions 0.5 → 0, a 4.7× lever on a
θ_D-sized quantity. MEASURED ON θ_D ITSELF (2026-08-17, the same
refreeze): the lever is 2.2×, 0.00291 at hop 0.5 against 0.00133 at
hop 0 on the realized-dose axis, and the two statements are about
different quantities, so both stand. The direction is what the
bracket rule protects: the gate's pinned 0.5 is the WORST end, so
θ_D registers the MAXIMUM over the hop
bracket exactly as it does over f_leak) at ≥ 500 H0
reps with held-out splits; inside the same unit, freeze the fit-health
margins (the δω excursion margin from the IDEAL construction's own
excursion, up to ~0.03 noiseless, plus H0, and consider widening the
δω grid or a two-frequency channel model so gate and analyze share an
interior optimum; the rate-saturation fraction, WITH ITS BASE, which
the freeze states as one of exactly two registered strings, `r_max` or
`1.1*r_max` (the fit's bound, or the search's true ceiling after
refinement; the runner hard-aborts on anything else, and round 12
brought the enumeration onto this page, where G2 can read it, instead
of leaving it in code where a plausible third spelling would have
stopped the flight at the gate); the runner's current
0.06/0.10 flags are informational only and would VOID every healthy
run if frozen as-is), and measure the dw-flag firing rate under
H0/H1. Then: the N0 arm and var(B̂) propagation (G2), the dressed
centre functions f_C, f_{C′} as code (and the A-arm and Π₃ centres,
Π₃ = the centred third invariant Π(r_i − r̄) of §6 item 4, as
DRESSED functions: the noiseless flown r̄ misses the ideal by
−1.4%/−3.4%, the gate record's two corner-arm lines, the arm
mapping re-stated at the refreeze re-measurement (round-4 note),
ON THE NOMINAL AXIS; round 3 identified the dominant
share as the §5 depth-0 off-by-one, so every dressing re-measures
on the realized-dose axis, and what remains after that fix is the
true dressing), the f_{C′} RECONCILIATION (round 3: the 0.54×
dressing factor times the bare ideal gives ~0.0008 while the frozen
config measures s²(C′) = 0.0002 ± 0.0002, a ~3σ gap: as specified
the C′ D-mag band would be frozen at a centre the gate's own run
misses, so the refreeze either closes the gap on the fixed axis or
the C′ containment reverts to REPORTED per §7's discrimination
clause), the difference-null generator fits
and margins (G2), family rates, the per-M retention criterion
replacing the transplanted 0.02 AND the CUMULATIVE realized-retention
face (the 1/√M floor reaches +3100% of target at the deep end,
arm-asymmetric; measured NOT to move the fitted rates, so a
certificate gap, not a verdict bias; G3), the DEEP-END KEPT-COUNT
MARGIN (round 5, and CLOSED ON THE ISOLATED AXIS BY AMENDMENT 1.8,
which is the lever the freeze took: at the registered 0.5% bound the
honest survival range puts the deepest cells at ~208-1233 kept counts
of 16384 against the kept-count floor of 50, a **4.17× margin at the
worst corner** BEFORE the unpriced |2⟩ leakage and the inversion clip,
where the retired 0.6% bound gave ~87-736 and 1.74×. Twelve deep
cells, 4 arms × 3 preps, each hold veto power over D-sign/W via
grid-incompleteness, so G1/G3 still freeze the EXPECTED deep-end
kept-count profile and the ≥ 3× requirement at the worst modeled
corner stands as the criterion; what changed is that today's numbers
now MEET it there instead of failing it. WHAT IS NOT CLOSED, and it is
the same isolated-versus-layered gap §9 carries: kept counts are
monotone in the per-gate error, so the requirement is EXACTLY
equivalent to an upper bound on the EFFECTIVE p2, namely
16384·exp(−970·p2·0.9) ≥ 150 ⟺ p2 ≤ 0.538%. The registered isolated
bound 0.5% clears it; layered error at the §9 bracket's 1.3× does not
(≈ 56 kept, 1.12×) and at 2× the deepest cells FLOOR outright
(≈ 3 kept), which VOIDs D-sign and W by grid-incompleteness. That
0.538% is therefore the natural candidate for the owed
`layer_fidelity_bound`, an UPPER limit on it rather than a safe value
since the same "before leakage and clip" caveat applies, and freezing
it is G1/G3's, not this paragraph's), G3(b)/G5 re-runs at the
frozen 21 × 3 grid (they ran at the retired 13 × 5), asymmetric
readout confusion (the analyze side already implements it; the gate
model is symmetric), the within-sector error split between the Z-like
and hop faces (the gate pins 50/50; s²(U)'s floor moves 0.0014 →
0.0003 across hop fractions 0.5 → 0, so the C′-below-U statement is
model-set in NUMBER, though C′ stays below U across the whole
bracket; note the H1 chain's implied s²(U) is 0.0017, same symbol
at a different config, reconciled when the refreeze re-measures
both on one chain), the G5 in-analysis J-estimator machinery (a VOID trigger, the G5
J pass band; the fitted δω are persisted as its raw input), and the
committed G1-stage record: the s²(C′) − s²(U) line, the stage-1
screen (fractional vs standard), the 8192-shot stage-2 row, the p2
ladder, §6.3's zero-spurious measurements, the r̄ dressing values,
the 1/√M cumulative-retention face, and the hop-split bracket exist
only in uncommitted output today and land in the repo with the
refreeze record. Optional, flagged by the rounds: common random
numbers across arms' phase tables would make the M-realization term
common-mode in d (free power; changes the H0 model, so it belongs to
the refreeze decision, not after it).

THE MARGIN'S LEVER INVENTORY, priced (2026-08-17, after Amendment 1.8,
because "three levers" was never a closed list and the freeze should
choose from the real one). The requirement is
16384·exp(−Σ_b w_b·p2_b·f_leak) ≥ 3·50 with w = 240 on each of the
three ODD bonds and 120 on each of the two EVEN ones (verified at the
source: a Strang step gives the odd bonds two XY blocks and the even
bonds one, at 2 rzz per block, so 3·240 + 2·120 + 10 preparation
gates = 970). Every free quantity in it is a lever:
- **The guard bound.** TAKEN, Amendment 1.8. Zero minutes, 4.17×.
- **Shots, targeted at the deep cells.** ALIVE and previously
  mispriced (Amendment 1.8): +1.16 QPU min at the old bound, against
  the +16.3 the uniform form costs. NEW HERE, and unprecedented in
  this repo along the depth axis, though `RECORD_PARITY` already flies
  unequal shots across ARMS (its r = 0 arms at 2×) for an estimator
  reason. Three structural conditions, all checkable and none
  prohibitive: M stays 1024 at every depth (Amendment 1.4's nested
  tables are what would break, and they key on M, not on shots per
  binding), the estimator already carries unequal counts through its
  √w weighting, and billing is shots-proportional, so this is
  budget-NEUTRAL rather than free.
- **The guard's FUNCTIONAL FORM. Investigated and PARTLY BUILT
  2026-08-17; building it is what settled how far it may be used.** The
  Class-1 rule is a MAXIMUM over the five used edges while the quantity
  it protects is the weighted sum above, so the guard cannot see WHICH
  edge is bad, and an odd bond costs exactly twice an even one. The
  weights are derived, not chosen, and the preparation distributes as
  two more rzz on every bond, so per edge they are
  (242, 122, 242, 122, 242) and they sum to the flown 970: no new
  number enters. Two lines both sitting at max p2 = 0.5%, one with its
  bad edges odd and one even, are admitted identically and differ by a
  factor 5.7 in deep-end kept counts (501 against 2843).
  **WHAT WAS NOT BUILT, and the reason is the finding:** a budget GUARD
  cannot fire. Any line passing the per-edge ceiling has budget
  ≤ 970 × 0.005 = 4.850, and the ≥ 3× requirement written as a budget
  is Σ_b w_b·p2_b ≤ 5.2149, so the ceiling already implies it with room
  to spare. Adding that guard would have been machinery that only looks
  like a gate, which is the shape this document's own machinery ledger
  exists to catch. A test pins the redundancy and FAILS if an amendment
  ever raises the ceiling past the ratio, at which point the guard
  becomes real and should be built.
  **WHAT WAS BUILT**, both of them reported and neither gating. (i) The
  line-selection SCORE is now the budget rather than `mean(p2s)`, which
  is not a cosmetic change: two admissible lines exist whose mean order
  and whose deep-end kept counts disagree (mean prefers 0.00260 over
  0.00290 while the deep end reads 1076 against 1483 kept), so the old
  score could pick the worse line, and it treated an odd edge and an
  even one as equal cost. (ii) `--calibrate` and the day-of path now
  PROJECT the deepest cell's kept counts and margin from the chosen
  line's OWN edges instead of leaving them at the guard's worst case,
  print them, and persist them in the chain record, with a printed note
  when a line's own projection falls under the registered ≥ 3×. That
  is the number the requirement is stated in, produced for the line
  actually flown. It does not discharge the margin's evaluator, which
  is G3's and stays on the machinery ledger: this is the line-side
  half, at isolated p2, and the counts-level expected profile is the
  other half.
- **The depth grid endpoint.** UNSEARCHED (§8 G1, §12): it changes
  n_gates directly, and it changes the signal lever arm with it.
- **The floor, 50.** Priced at 50 → 29 and blocked on G3. Note for
  whoever freezes it: it enters the committed gate at
  `corner_beat_gate.py` line 181 with no comment and no derivation.
- **f_leak's worst end, 0.9, and the ≥ 3× itself.** Both would widen
  the margin, and NEITHER may be used for that. 0.9 sits deliberately
  past the counted 12/15; retreating to 0.8 would "buy" 6.8× by
  relaxing a conservatism. The 3 is a round-12 judgment call with no
  measurement and no precedent behind it in this repo (swept
  2026-08-17: nothing in `fw.Confirmations`, nothing in the OpenArcs
  registry, nothing in `docs/CAUGHT_ERRORS.md`, no prior flight with
  an analogous count-margin rule; the concentrator's power ≥ 3 and
  this document's own 3σ̂ thresholds are the nearest 3s and neither is
  its source). A criterion may be moved pre-data for a stated reason;
  it may never be moved BECAUSE the design fails it, and both of these
  would be exactly that.
- **Robustness instead of survival**, the concentrator's own answer to
  a fragile deep end, recorded here because it is the repo's precedent
  and it is not on the list above: that flight did not buy counts at
  the deep points, it required its verdict to hold WITHOUT them
  (`run_concentrator_reloaded.py`: the A-sign recomputed without
  depths 6 and 8, same sign required). Here the no-drop rule and the
  grid-incomplete VOID make that a registration change rather than an
  analysis choice, and a pre-registered reduced-grid secondary
  analysis is not the same object as dropping points after data
  exists. Named, not proposed.

COMPLETENESS (round-2 repair; this list is not the whole freeze).
The freeze ledger is the runner's CONSTANTS_MANIFEST (its LENGTH is
not quoted here on purpose: the analysis printout reports it on every
line of its manifest audit, and round 12 grew the list by three while
leaving a "34" behind in this sentence, which is the third prose
count of a moving quantity this document has had to retire.
Round 10 added layer_fidelity_bound, the day-of gate's frozen
threshold, and leakage_ket2_price, so the two owed quantities that
had NO ledger hook can now block; round 12 added kept_count_floor,
deep_end_kept_profile and boot_nan_replicate_bound;
round 3 added fit_residual_bound, round 5 added time_axis, round 6
split T1-CLEAN into its three faces and gave the dressed A and Π₃
CENTRES their own keys beside the bands; a
MISSING entry blocks submission exactly like an unfrozen one.
Every machinery entry NAMES the manifest keys it discharges, in
brackets, so that "empty the ledger by its names" is an operation an
executor can actually perform: a test built on 2026-08-17 found seven
keys covered by an evaluator whose entry spelled them differently
(`thetaF` for `theta_F`, `difference_null_evaluators_12_23` for
`null_margin_12`), which is the same defect round 10 named when it
said an executor emptying the ledger by names "would have shipped the
W half only". THAT SENTENCE WAS TRUE OF FIVE ENTRIES OF TWENTY-FOUR
WHEN IT WAS WRITTEN, and the gap was found the same evening by reading
the ledger against the manifest rather than reading the sentence
(2026-08-17): fifteen further entries plainly discharged a key and
named none, so an executor emptying the ledger by name would have
walked past `s2C_floor`, `N0_clean_thresholds`, `J_pass_band`, the
D-mag bands, the dressed A and Π₃ centres, the deep-end profile and
the |2⟩ price. They name their keys now, and the promise is enforced
from BOTH sides instead of asserted. The invariant, and it is the
honest form of "empty the ledger by its names": every manifest key is
either NAMED by a ledger entry or DECLARED as having a live consumer
already in the runner, a list of twelve given explicitly in the code
with the place each is consumed. A key in neither set is a number
nobody has to build anything for, which is the precise mechanism by
which a freeze could release a flight with a dead arm. Four ledger
entries name no key, correctly, because they discharge CODE and not a
number: the B̂ 3×3 generator fit, the depth-0 CI term, the verdict-line
composer's wiring, and the day-of gate-reproduces-§8a check. Overlap
between the two sets is
legitimate and real: `kept_count_floor` has a live floor comparison
AND an owed ≥ 3× margin evaluator, and the time axis is read by the
runner while the gate-side refit is still owed. Beside the number ledger the runner carries that
MACHINERY ledger,
round 5: a hard-abort list of unbuilt verdict consumers that the
gate work empties as the code lands, so a frozen number without
its evaluator never releases submission; round 11 closed the level
ABOVE the consumers, where the same hole had opened one storey up:
the ledger held an entry for every owed NUMBER's evaluator but none
for the two pieces of machinery that turn those evaluators into the
RECORD's single line, namely the VERDICT-LINE COMPOSER with the VOID
routing, trigger precedence and PARTIAL resolution of §7, and the
¬D-SIGN PARTITION EVALUATOR that separates FALSIFIED from Anti-D
from INCONCLUSIVE using θ_F and κ. The runner prints per-arm lines
and the strings "κ power condition PENDING"; naming a pending NUMBER
is not the same as owing the CODE, and an executor emptying the
ledger by its entries would have released a flight whose
falsification arm had no evaluator at all. Both are ledger entries
now, as is the parking rule's third condition, "the committed gate
fails to reproduce §8a", whose two siblings hard-abort in the runner
while it was a human act: --certify checks PARITY against the gate,
never that the gate reproduces the governing numbers, so the G1
committed-table mode owes that check), and
beyond the items above it still owes: θ_W (producer registered in
§7: inside the refreeze package), θ_F's evaluation at the registered
0.75 fraction and κ, both
D-mag bands and the A margins, the s²(C) floor and its adjudication
TOGETHER WITH ITS TWO PRODUCER INPUTS, E[s²(C)|H0] and its σ̂, which
§7 registers as owed here and which round 12 named without adding
(round-13 repair; they come out of the same H0 chain as θ_D),
the band-validity window and day-of width-scaling functions, the
N0-CLEAN thresholds, the T1-CLEAN promotion WITH its µs→γ̄ bridge
(§9), the readout-confusion ladder to the 2%
guard bound (the G1 model used symmetric 1.5%), the day-of
LAYER-FIDELITY gate machinery + its frozen bound (§9, round 3: the
isolated-p2 guard does not see simultaneous-layer error), the
leakage-to-|2⟩ pricing (§11), and a committed gate
mode emitting §4's Floquet-Lindblad table and EP locations (rounds
2 and 3 both reproduced all twelve cells and all four EPs from the
gate's own primitives, last-digit rounding only, but the authority
table is uncommitted today).

## 9. Guards

**Class 1 (pre-submit, prevents the spend):** line rule (every qubit
T2echo ≥ 150 µs, max/min T2echo ≤ 2, readout ≤ 2%; the runner reads
the properties-file T2 as the echo-calibrated PROXY, the first aux
station's Hahn echo being the in-data cross-check, round-9 naming); calibrated 2q
error on the used edges p2 ≤ 0.5% (§4, Amendment 1.8; the guard now
sits ON the measured 0.5/0.75/1% ladder's first rung, where power is
2.7-2.9, against 1.1 one rung up and dead at 1%, §8a), WITH the
round-3 caveat that properties-file p2 is an ISOLATED-gate number
while the flown circuit runs three disjoint rzz simultaneously per
half-layer, and simultaneous-layer error on Heron-class devices runs
~1.3-2× the isolated value (the reason layer fidelity/EPLG exists as
a metric). TIGHTENING THE ISOLATED BOUND DID NOT CLOSE THAT GAP and
must not be read as having done so: a line passing at isolated 0.5%
can be running layered at 0.65-1.0%, which spans the ladder's
power-1.1 rung and reaches the measured dead point, and §8a's
kept-count arithmetic puts the same span at 1.12× the floor falling
to a floored deep end. So
the day-of re-gate ADDS a pre-job layer-fidelity measurement of the
actual 6-qubit line in the actual layer structure (a
calibration-class input, inside the §9 anti-circularity fence;
machinery + its frozen bound are owed BEFORE freeze, §8a, where the
≥ 3× deep-end requirement now hands that bound a computed candidate,
an effective p2 ≤ 0.538%, to be frozen by G1/G3 and not here); the realized
physical-qubit profiles of C and C′ in their asserted classes
(post-transpile assertion); FRACTIONAL-RZZ exposed and taken
(use_fractional_gates honored, rzz in the target's operation names, no
CZ decomposition in the transpiled circuits, CZ = NO-FLIGHT, the
record-parity assertion; load-bearing per §4/§8a); DD and both
twirling channels disabled fail-closed (SamplerV2 carries no further
resilience knob; the runner raises if the options object refuses);
backend-operational abort plus a pending-queue depth WARNING (by
design the queue never blocks, only non-operational does; and the
operational check FAILS CLOSED when the status cannot be read at
all, registered in round 15 after round 14 made the code do it: a
Class-1 guard that cannot see its input has not passed); billing
projection under the cap; the dose certificates on the BUILT tables
also run pre-submit as a hard abort (the Class-2 face below reads
the FLOWN tables: belt and suspenders, round-5 class note); the
pre-registration, gate, and runner committed at a real hash BEFORE the
Batch opens; every CONSTANTS_MANIFEST entry frozen (a missing entry
blocks, the runner's list is the §8 freeze manifest) AND the
runner's MACHINERY ledger empty (round 5: unbuilt verdict
consumers hard-abort submission independently of the numbers).

**Class 2 (in-job, protects the verdict, can only VOID):** the §7
bank, all TEN: N0-CLEAN; T1 CLEAN on interleaved in-situ T1/T2* PUBs
(BRIDGE registered here, the numbers frozen by G4 under THREE
manifest keys, one per face (T1_clean_projection_band,
T1_clean_scale, T1_clean_station_band; round 6: one scalar cannot
carry three bands): each station T1_l converts to engineered units via
the transpiled per-Strang-step wall duration τ_step,
γ_l^hw = τ_step/(4·T1_l·dt), the Γ_l/4-equivalent dephasing in the
simulation's J units, and γ_l^hw/γ̄ to read it IN γ̄ (round-14 unit
repair: the formula returns J units, ≈ 0.0033 at τ_step = 0.5 µs and
T1 = 250 µs, which is ≈ 0.033 γ̄ at the flown Q = 10; two analysts
sharing one conversion is the whole point of the bridge, so both
steps are written), τ_step recorded in the day-of addendum; CLEAN then has THREE registered faces, round-4
unification of what §8 G4 and this section had stated as two
different tests: (i) the PROJECTION face, the binding one, the
γ_l^hw profile's projection onto C's mirror-odd part inside its G4
band, since a mirror-symmetric background of any size moves d
exactly zero; (ii) the SCALE face, site scatter ≤ the G4 scale,
proposed 0.25γ̄, a coarse sanity bound honestly near-non-binding;
(iii) the TIME-STRUCTURE face, the recorded staircase killer:
station-to-station consistency of each site's fitted T1 (two
delays, fittable) and single-delay echo SURVIVAL FRACTION (one
delay per station, a fraction, not a fit; round-6 precision)
within a frozen band across the four stations, beside N0's
halves-consistency); dose certificates on the flown phase tables, the (1,2) and (2,3)
difference-nulls,
fit health, the s²(C) floor (per §7's rerouted adjudication: usually
DATA, VOID only beside a failing dose certificate or N0-CLEAN), the
G5 J pass band (scoped to the magnitude/anchor arms, §7), grid-incomplete
(Amendment 1: any science cell missing, binding-broken, or
flown-but-kept-count-floored), and the
band-validity window (if the device leaves the
window the design is invalid: abort, redesign, re-register).

**Day-of re-gate (pre-committed rule, anti-circular):** consumes
PRE-SUBMISSION calibration-class inputs ONLY (backend properties
snapshot; at most a separate pre-job calibration PUB, never the science
job's own CAL or T1/T2* PUBs, which are Class-2 guard inputs and would
be circular; the staircase's own fence); replaces band WIDTHS per the
G1-frozen scaling functions, never centers, never thresholds; its
output is committed as the Day-of addendum BEFORE the Batch opens. The
committed document governs, even against the runner's own printout
(the concentrator stale-band lesson).

**The order:** gate green → empty review of runner + gate records
(G7) → document + runner + constants committed at a hash → day-of
`--calibrate` (the calibration pull and line selection) → **the dated
chain record it writes, COMMITTED** → day-of `--certify` (likewise
required from the same day, and not strippable) → **the DRY RUN,
`--hardware` without `--yes`, which prices the flight and prints the
realized job plan** → re-gate → Day-of addendum committed, BOTH
FORMS, carrying that job plan → Class 1 guards → submission. The two
commit steps inside the day are not bookkeeping: the submit gate
requires the runner, the constants, the dated chain file AND the
addendum JSON all tracked and clean, and the chain file cannot be
part of the pre-day commit because `--calibrate` writes it on flight
day. An executor following an order that ended at "addendum
committed" reached submission with an untracked chain file and was
hard-aborted with the Batch window open (round-13 repair; round 12
added the two day-of runs to this order without carrying their commit
obligation with them, the same shape round 11 closed for the addendum
JSON itself) → in-job Class 2 → analyze from persisted
counts → RECORD (with its own empty round, the house rule).

## 10. Costs (at the candidate point; billing anchors by measurement)

- Billing anchors: 0.283 ms/shot (the A-vs-B jobs, shallow, recorded
  in the concentrator doc), 0.316 ms/shot (the concentrator job
  itself, 8 aux delay PUBs aboard; recomputed HERE as 119 s/376,832
  shots, 12% above the 0.283 the same record measures on its A-vs-B
  jobs (that record calls neither figure a "model" and does not
  characterize those jobs as shallow; both descriptions are ours,
  round-15 repair), a tension that RECORD itself flags against its
  85-99 s
  projection), 0.309-0.327 ms/shot (the two
  staircase flights, delay-bearing like this design; per-shot rates
  recomputed HERE from the RECORDS' billed seconds, 324 s/1.049M and
  380 s/1.163M shots, the records store seconds): the band uses
  the delay-bearing anchors.
- Per-PUB budget table (frozen, §4/§8a; aux plan per Amendment 1): 4
  shot-bearing arms × 21 depths × 3 preps = 252 science PUBs at 16384
  total shots each (M = 1024 bindings × 16 shots on the 20 swept
  depths; the twelve depth-0 PUBs fly unswept single circuits at
  the same 16384, Amendment 1.4) ≈ 4.13M shots, plus
  2 head CAL PUBs at 16384, 24 aux station PUBs and 20 per-job CAL
  PUBs at 4096 → ~4.34M shots ≈ **23.7 QPU min** at the delay-bearing
  anchors. Band quoted 21-25 QPU min, accepted by Tom 2026-08-16
  (second decision, post-pin numbers); M is a G3 freeze output and
  the band re-freezes with it. Per-binding overhead: no GROSS
  per-circuit surcharge is visible (2026-08-17, read-only
  job.usage() query, re-run same day after a first note misquoted
  69 s: the concentrator's 376,832 shots and 119 s are what that
  record STATES; the 6144 bound circuits and the 24-sink/12-no-sink
  PUB split are OUR reconstruction of its payload, checked against
  the recorded shot total (24·256·32 + 12·8192 + 10·8192 = 376,832)
  and labelled here in round 16 as derived, the way §10 already
  labels the 0.316 division. In that reconstruction the third group,
  10 PUBs at 8192 shots, is the 8 delay-bearing aux PUBs plus the 2
  CAL PUBs (round-17 reconciliation of an 8-versus-10 that stood
  unexplained in two places). That division, 119 s over 376,832
  shots, gives 0.316 ms/shot, INSIDE the delay-bearing band
  0.309-0.327; that job carried those 8 delay-bearing aux PUBs,
  so binding overhead and delay cost are not separable there, and
  the M = 1024 overhead stays unmeasured until flight). The
  0.327 anchor's margin over the measured 0.316 is ~3%, and the
  remaining anchor caveat is per-shot DURATION: these science
  circuits carry up to 970 2q gates, deeper than any anchor
  flight. The runner's mid-flight billing brake (Amendment 1) is
  the backstop, not the projection.
- Circuit depth at the candidate point: the RZ injection layer sits
  between Strang steps and does not commute with XX+YY, so half-layers
  do NOT merge across steps: 8 two-qubit blocks per step (XY blocks under the §5 reduction) × 60 steps
  ≈ 480 blocks; at the pinned FRACTIONAL gates (2 per block; the flown
  total is 960 + 10 preparation gates = **970 two-qubit gates**,
  asserted post-transpile), deep-end error exponent = 4.85 exactly at
  the registered p2 bound of 0.5% (970·p2; Amendment 1.8, and the
  5.82 that stood beside it was the retired 0.6% bound's, while the
  1.3-2× layered inflation of §9 carries the exponent to 6.3-9.7,
  which is where §8a's deep-end margin goes from 1.12× to floored;
  standard gates were modeled at 3 per block, ~1440, the retired
  configuration) against a signal exponent γ̄T_real·(10/3) ≈ 2.95 at
  the frozen point (59 realized layers, §5). WALL DURATION,
  pre-registered estimate (round 3; the transpiled schedule is the
  authority and the day-of addendum records the measured value): 6
  rzz sub-layers per Strang step at Heron rzz durations ≈ 60-100 ns
  gives τ_step ≈ 0.36-0.60 µs and a bare deep end of 22-37 µs (360
  step sub-layers plus the preparation's ten, the count the runner
  books; round 13 caught round 12 quoting 22-36 from the step
  sub-layers alone while quoting the BUFFERED end from the full 370,
  so the pair was not the same model after all). Against the 150 µs
  T2echo line rule that is 0.15-0.25 × bare and **0.22-0.37 ×
  buffered**, and the background two-site coherence exponent
  2·τ_deep/T2echo is 0.30-0.49 bare and 0.44-0.74 buffered, i.e.
  10-17% of the engineered exponent ≈ 2.95 bare and **15-25%
  buffered**. THE COMPARISON'S COMPARAND IS THE BUFFERED ONE (round
  14: "the guard's" named a guard that does not exist, since §9's
  Class-1 line rule gates the DEVICE's T2echo ≥ 150 µs and nothing
  gates the circuit-duration ratio; this is reported context, as the
  paragraph's own last sentence says), because that
  is the circuit the executor actually schedules; the bare pair is
  the physical floor and is quoted beside it, never instead of it
  (round-13 repair: this sentence had made the T2 comparison on the
  bare number while the same section calls the buffered one the
  conservative estimate thirty lines later). Either way the exponent
  is carried by N0's B̂, not by assumption.
  ROUND-12 CORRECTION, and it is a correction of round 11's own
  repair: this estimate used to add ~100 ns per step for the
  injection rz layer, giving τ_step 0.46-0.70 µs and deep end 28-42
  µs (and the coherence exponent 13-19%). RZ is a VIRTUAL frame
  change on Heron and spends no time, so that allowance was buying
  6 µs of imaginary circuit across 60 steps. Dropping it also
  explains the doc-vs-runner gap, which round 11 asserted as "≲1.5%"
  by adopting a reviewer's estimate instead of computing it: the
  actual bare gap was **14-21%** of the with-rz estimate (28.2-43.0
  vs 22.2-37.0 µs, BOTH at the 370 count, the denominator being the
  larger of the pair; round 15's repair of round 14's repair, which
  wrote 12-20% by comparing a 60-step-only 27.6-42.0 against the
  370-count 22.2-37.0, i.e. by reintroducing at the level of the GAP
  the very asymmetry it had just removed from the two ends, and
  without naming its denominator), and it was
  almost entirely this rz allowance, the preparation-booking
  CORRECTION (five sub-layers, the difference between the old count
  and the right one) being ~1.4% of t_deep and 5-8% of the gap
  itself. Round-17 subject fix on top of round 16's denominator fix:
  the booking itself is TEN of 370 sub-layers, 2.7%, and calling the
  correction "the booking" put a factor 2 on it in the same clause
  where the denominators had just been separated. The two models now agree by construction. The RUNNER computes the
  same layered SUB-LAYER COUNT but applies a 1.5× scheduling buffer
  (round-8 repair of its earlier serial gate-count estimate, ~145
  µs, 4× physical, which would have pushed the long T1 station
  point past Heron T1. With the rz allowance gone the two are ONE
  model at two buffers, in the strict sense the round-9 phrase
  claimed and the round-11 text could not support: bare 22-37 µs
  here, buffered t_deep ≈ 33-55 µs in the runner at its 1.5×
  scheduling factor, both from n_sublayers = 10 + 6·60 = 370, the ten
  being the five Givens on OVERLAPPING pairs (l, l+1) at two rzz
  each, which cannot schedule in parallel and which the runner had
  booked as five. The coherence exponent above, 10-17% bare, reaches
  ~25% at the buffered end. (Round-14 repair of round 13's own: the
  bare end was corrected in the paragraph above and left at the
  365-sub-layer numbers here, thirty lines below, in the sentence
  that NAMES 370. A half-replaced passage is the failure mode this
  document has now produced four rounds running, and it is why the
  two ends are stated from one count, out loud.) Which estimate is the SAFE one depends on what it
  feeds, and the two directions are opposite (round-13 repair of a
  sentence that had them the same way): for the T2 comparison the
  buffered number is conservative, since it assumes MORE decoherence;
  for the aux-station delays it is the riskier one, since the long
  station sits at 1.5 × t_deep and a delay past T1 reads noise, which
  is exactly what the retired ~145 µs serial estimate would have
  caused. The stations derive from the buffered number and are
  granularity-aligned; the transpiled schedule at submit is what
  settles both. Both remain pre-registered ESTIMATES and the
  transpiled schedule at submit is the authority), and the
  aux-station delays derive from the
  runner's buffered
  number; the transpiled schedule at submit is the authority and
  the day-of addendum records the measured τ_step. (Why 16384 shots buy P(detect) 0.815 → 0.920 while
  power stays 2.9: the gain runs through θ_D shrinking, H0 being
  shot-dominated, while H1's spread is dominated by the M = 1024
  realization term, which shots cannot buy down; common random
  numbers across arms, §8a's optional note, is the lever that
  would.) **Mechanism relabel (G7, measured):** the 2-per-block
  count is the SECTOR REDUCTION of §5 (pure XY blocks; the ZZ part is
  the degree diagonal, one-qubit), not a fractional-gate identity: an
  exact XXZ block is Weyl (t,t,t) and needs 3 entanglers. A CZ-basis
  transpile of the SAME reduced circuit reaches ~2 CZ per block at
  optimization level 2+ (measured 106 vs 106 2q on a 6-step skeleton;
  level 1 gives 4 per block), so the fractional advantage at equal
  optimization is per-gate error and 1q overhead, not the gate count;
  G1's standard-gate power numbers were computed at the 3-per-block
  model and the standard branch's retirement inherits that model
  (re-examining it is a named §12 follow-up, not a flight question:
  the flown configuration is fractional either way). The bridge from exponent to survival runs through
  the out-of-sector fraction f_leak of the error weight (post-selection
  removes only what leaves the sector): §5's 7-of-15 in-sector count on
  a differing pair gives f_leak = 8/15 → survival e^{−4.85·8/15} ≈ 7.5%,
  while X-type errors on non-differing pairs also leave the sector
  (the same count on a NON-differing pair gives 12/15 = 0.80;
  the flown bracket end ~0.9 sits deliberately beyond the counted
  value, conservative, round-3 note); at ~0.9 → ~1.3%; the honest
  range at the 60-step
  point is **~1.3-7.5% post-selected survival AT the registered
  0.5% guard bound** (exact 1.27/7.53; Amendment 1.8 moved this
  range from the retired bound's 0.53/4.49%, which is the whole of
  where the deep-end kept-count margin came from, §8a; fractional;
  the retired standard-gate range was 0.2-2%). The v4-v5 reading that
  this kills the flight was WRONG in mechanism (§4, §8a): the deepest
  points do not need to survive at full SNR, because the estimator's
  kept-counts weights devalue them smoothly; G1's counts-level MC
  carries the whole depth profile and still resolves the split at the
  candidate budget WITH fractional-RZZ (2 gates per block ≈ 960 2q at
  the deep end). Standard gates (~1440 2q) remain unflyable in power
  terms at any tested budget AT THE RETIRED 3-PER-BLOCK PRICING
  (round 6: after the §10 relabel the park's power comparison
  awaits the honest ~2-CZ re-pricing, the named §12 follow-up;
  conservative in direction, since the honest count can only
  improve the unflown standard branch). If fractional gates are unavailable on
  the flown backend (Class-1 guard, §9), the candidate is PARKED
  (recorded, not deleted); shallower-grid redesigns reopen it.

## 11. Honesty notes

- For THIS room the class split is DERIVED in closed form (round-4
  symbolic verification of the overlap route: the pair-odd matrices
  have entries 1/12 and √3/36, eigenvalues {±1/6, 0} and
  {±√3/18, 0}, giving w = 1/3 and both rate triples exactly), and
  m̄ = 5/3 = 2·(1 − 1/6) likewise. What stays gated-numeric at the
  eigensolver floor is the GENERAL R₉₀-locus statement (M0
  diagonality and the 2160/2304 denominators, e3dbab0's open
  residue); the flown room no longer rests on it. A derivation landing
  before the flight upgrades the claim, not the protocol.
- The room is a uniform-chain resonance; bond disorder detunes it.
  The device requirement, remeasured at the estimator level under the
  pinned gauge (G5 v2.1): σ_J = 1% costs ≤ 10% of s²(C), 2% ≤ 25%,
  and the J response is even (any error inflates); the requirement
  stays σ_J ≲ 1-2%, owned by G5, cost carried in the D-mag band.
- No repo document states the Q-validity fence for this room's
  compression; this doc carries both corridor factors itself (§4),
  measured across three rounds with a 2-4% method spread, gate to
  arbitrate.
- LEAKAGE out of the qubit space (|2⟩) is named and only bounded,
  not budgeted (round 3): over 970 2q gates at 1e-4 to 1e-3
  per-gate leakage the cumulative leaked population is percent-scale
  to tens of percent at the deep end; a leaked qubit often
  discriminates as |1⟩, so the one-magnon post-selection does NOT
  remove it, and its site- and depth-dependence mimics the signal's
  own shape. The §5 sector-leakage certificate covers the PREP, not
  this. The refreeze prices it (a leakage term in the G1 counts
  model or a measured bound from the CAL/aux record); until then it
  is the largest UNPRICED systematic on the page.
- The engineered γ̄ = J/10 = 0.1·J is a resolvability choice. The
  canonical hardware-anchored point is a RATIO, γ₀ = 0.05 WITH
  J = 0.075 (Q = 1.5; the pin's owner is docs/Q_REGIME_ANCHORS.md,
  CAUGHT_ERRORS cites it); this flight sits at
  Q = 10, a deliberately far slower dephasing per hop, and no "2×"
  comparison of the bare numbers survives the ratio reading.
- v1 cited a ζ factor-2 from the price-pair doc; the pinned house form
  is ζ_c = ζ_shift/4 in the staircase doc. No ζ observable flies here;
  kept as a convention-trap reminder.

## 11a. The first-flight stack (reported readings, never verdict arms)

The verdict arms (§7) are D-sign, W, D-mag and A, four, which is what
the composition algorithm enumerates and gates its named lines on;
the two difference-nulls have equivalence-test SHAPE but are VOID
TRIGGERS, not arms, which is how the scope map and the precedence
ranks 6a/6b already treat them (round-18 correction of a "plus the
nulls" that made the arm count ambiguous). The Confirmations
registration scope stays D-sign + W only, §7.
But the flight is not one verdict: it is a STACK of repo objects each
of which reaches hardware here for the first time, and each stack
element stays readable even under a not-D-sign outcome or a
fit-health VOID. Discipline first: every reading below is REPORTED,
never a registered verdict arm, and none of them can rescue or
overturn §7 (no verdict inflation). The post-flight RECORD reads
this stack element by element, whatever the verdict line says.

Swept 2026-08-17 (one store-sweep agent, per element, against
`fw.Confirmations` (24 entries), the `experiments/` flight records,
the OpenArcs registry, `docs/ANALYTICAL_FORMULAS.md`, and
`docs/CAUGHT_ERRORS.md`; findings folded into the wording):

1. **The degree diagonal as a circuit element** (F152; D10 step 3;
   PROOF_K_PARTNERSHIP's V_eff(l) = #bonds − 2·deg(l)): the ZZ part
   flies as its one-magnon-sector equivalent, a boundary rz layer on
   qubits 0/5 (§5). First flight of the object as hardware; the
   nearest flown kin is the F95 angle-steering flight's per-chunk
   RZ(Ω·Δt) injection (Kingston 2026-05-16), uniform and not
   degree-derived.
2. **The Lemma 5 spectrum read from the beat** (PROOF_R90_FROZEN_DIVISOR
   Lemma 5, λ_k = 4cos(kπ/N) + N − 5; index caution: this doc's k is
   energy-ASCENDING while Lemma 5's runs descending, so the map is
   k ↦ N − 1 − k. That direction is DERIVED from the cosine over
   k = 0..N−1, not stated in the proof, and the same file's Lemma 3
   indexes the same operator's eigenvalues ascending, which is why
   the standing rule here is to map VALUES and never indices,
   round-14 precision): the fitted dyad frequencies ARE the λ_k
   differences, a spectrum test independent of every rate claim.
   `WHAT_THE_R90_LOCUS_BUYS.md` records for its own R₉₀/frozen-divisor
   territory that `fw.Confirmations` returned nothing and no flight
   touches it; the one-magnon spectrum read has no flight either
   (this doc's 2026-08-17 sweep).
3. **A site-VARYING engineered dephasing profile, on the R₉₀ locus**
   (the e3dbab0 mirror-transversal classes as flown profiles,
   including each class's decoupled dyad). Engineered injected
   dephasing per se has flown three times (EP-onset's uniform twirl,
   the concentrator's and A-vs-B's single-site sink, itself a
   degenerate site-varying profile); a MULTI-SITE non-uniform
   profile (one lit site per mirror pair), and any profile ON the
   R₉₀ locus, has not.
4. **The uniform-compression scalar**: the U arm must read all three
   dyads at ONE rate, (10/3)γ̄, the §2 compressed-dissipator scalar
   (this doc's own exact 1/6 route, G1-compared). Kin of the
   45e2f29 two-spin-zeros arc but NOT covered by it: that law's C_l
   silences decide the protected rooms, and this room is breaking
   (no C_l protection on (1,1), per the sweep section). The scalar
   reading rides along as the U arm, whatever C and C′ do.
5. **The mean anchor** (the A arm): r̄ invariance under the moved
   profile, the compression machinery's trace face; appears in no
   prior flight or confirmation entry.
6. **Gauge-pinned Floquet modes + Givens preparation of a two-mode
   dyad** (ψ_i + ψ_j)/√2. QUALIFIED first: F129
   (`f129_standing_fringe_kingston_july2026`, Confirmations 24) flew
   a Givens network compiled from the actual Trotter step's modes;
   new here are the two-mode DYAD preparation and the mode-sign
   gauge pin (§8a, the pre-pin 0.985 → 0.815 lesson).
7. **The M = 1024 randomized-RZ channel with per-binding counts.**
   QUALIFIED first: the channel itself is the concentrator genre
   (M = 256, and A-vs-B at K = 16), both pooled at save time;
   per-binding persistence is exactly the concentrator RECORD's
   named future-flight item, flown here, and M = 1024 is this doc's
   own G3-frozen step past that flight's M.

Elements 1 and 2 are Tier-1 objects (proof-backed); 3 and 5 are the
e3dbab0-adjacent gated-numeric results the honesty notes (§11)
already scope; 4 is this doc's own exact identity (G1-compared);
6 and 7 are instrument capabilities. A not-D-sign flight still
returns: the U-arm scalar test (4), the spectrum test (2), the anchor
(5), and the instrument records (6, 7). A fit-health VOID still
returns the raw per-binding record for all of them.

## 12. Scope fence (named follow-ups, out of scope here)

The mirrored maximizing transversal {1,2,5} (the l → 5−l image of the
pinned C) as a repeat arm; a second
non-maximizing corner as a deliberate positive control; other-N mixed
rooms (round-3 correction: degenerate ω ≠ 0 spaces of dim 2 exist at
EVERY N ≥ 4 in the 3..8 range; what is unique to N = 6 is dim ≥ 3,
the only such room in the range, so the dim-2 rooms at other N are
themselves a named follow-up family, beside ω = 0, Δ ≠ 1, and
topology variants, all unexplored); the DEPTH GRID ENDPOINT as a search axis (§8 G1: the committed gate
fixed it at 60 steps in all sixteen configurations and searched the
spacing alone, so shortening the deep end has never been priced
against its own cost, which is the signal lever arm: both the
engineered exponent and the error exponent scale with depth, so this
is a trade and not a free margin); a second backend / second day repeat; the derivation
upgrading the GENERAL-locus w statement from gated to derived (for
the flown room this is DISCHARGED, §11 first bullet, round 4); the standard-gate variant if a target without fractional RZZ ever
must fly (retired by G1: worst-end power ≤ 1.7 at any tested budget
under the pinned estimator; NOTE per §10's mechanism relabel, that
retirement was priced at the 3-per-block model, while the sector-
reduced circuit reaches ~2 CZ per block at optimization level 2+, so
a re-pricing at the honest CZ count and CZ error rates is the first
step of any such revival).

## 13. Amendment protocol

After freeze (§9's commit), changes ONLY as numbered pre-data
amendments appended here, each recording what changed, why, its
committed hash, and that no science data existed when it landed (the
staircase Amendments 1-2 shape). Post-data, nothing changes; the
RECORD applies the committed rules and names deviations as instrument
deviations.

## Amendment 1 (2026-08-17, pre-data; the analysis chain and PUB plan as flown)

No science data exists. This amendment records where the built runner
(G7) deviates from or sharpens the §5/§6/§10 letter, each with its
measured reason; the runner's commit gate refuses submission unless
this section is present in the committed document.

1. **Pooling before inversion, one clip.** §5/§6.1 said "confusion
   inversion applied to per-binding counts before post-selection". The
   flown chain persists per-binding RAW counts, pools them over the
   (resampled) bindings, applies ONE confusion inversion to the pooled
   vector, and clips negative quasi-counts ONCE on the six one-magnon
   entries. The inversion itself is linear, so pooling first is exactly
   equal; the CLIP is the point where the chain becomes nonlinear
   (§8 G2's question, hereby answered), and per-binding
   invert-then-clip at 16 shots/binding is depth-biased (measured:
   moves d by ~6% at flight scale). Negative mass over the full 64-dim
   inverted vector is reported as a health number.
2. **Kept-count floor and dead channels.** A pooled cell with < 50
   kept quasi-counts reads NaN (the committed gate's own floor, on a
   different quantity); a channel whose every trace has < 5 live points
   fits NaN (the gate silently returned r = 0 there; the runner's NaN
   is the repair, and the θ_D refreeze runs the runner's fit_rate for
   exactly this reason). An INCOMPLETE grid (any science cell
   missing, binding-broken, OR flown but under the kept-count floor
   on the full data, persisted as floored_cells; one route, round-4
   repair) VOIDS D-sign and W outright: §5's no-drop rule,
   made executable.
3. **Registered estimator constants.** R_MAX_FIT = 0.9 (the fit bound
   AND the rate lattice, terminal spacing r_max/240) and R_BOOT
   (registered TARGET 2000; the code's provisional default is 500,
   the committed synthetic analyses used 500, and 2000 becomes
   operative when the entry freezes) join `corner_beat_constants.json`;
   R_MAX_FIT freezes with the θ_D refreeze (the gate froze θ_D on a
   config-dependent r_max), R_BOOT freezes with θ_F (its percentile
   endpoints feed the containment tests; one owner each, round-4
   alignment with the manifest notes). WHAT R_BOOT COSTS IN WALL
   CLOCK, **measured 2026-08-17 on the final `--aer --full` artifact at
   the flight configuration** (M = 1024, all 252 science cells), which
   is the number this Amendment owed: **54 s at R_BOOT = 500 and
   214.5 s at the registered 2000.** Both ends were TIMED and neither
   was scaled from the other, because scaling is exactly how the
   earlier attempt went wrong: round 11 anchored on 48 s at
   R_BOOT = 500 from a FOUR-CELL M = 64 synthetic and multiplied by the
   cell count, but the per-replicate cost is 4 arms × 3 channels of
   fitting and is largely INDEPENDENT of the cell count, and a
   four-cell artifact does not exercise the fit path at all, since
   `fit_rate` returns immediately on a trace with fewer than five live
   points. The measurement settles it: the full-grid number at 500
   replicates is 54 s against that four-cell 48 s, so the whole cell
   count bought about a tenth, not the factor 63 the extrapolation
   assumed.
   TWO CONSEQUENCES, and the first RETIRES a warning this paragraph
   used to carry. The analysis is minutes, not hours, so nothing about
   its duration tempts an executor to cut R_BOOT; the sentence that
   said otherwise was written against the guessed number. The second is
   a from-below argument FOR the registered 2000 that this document did
   not have: on the same data the d interval reads [+0.00173, +0.00765]
   at 500 replicates and [+0.00151, +0.00804] at 2000, so the cheap
   setting returns the NARROWER interval, and a narrower interval is
   precisely what turns INCONCLUSIVE into FALSIFIED under the θ_F
   containment rule. R_BOOT stays a registered constant only θ_F's
   owner may move, and that reason is now measured instead of asserted.
   The number lands in the day-of addendum beside τ_step.
4. **Nested frozen tables.** One channel realization per (arm, prep,
   binding), drawn at the full depth and sliced per grid point, so
   binding m is THE SAME realization at every depth (what the gate's
   trajectory model computes, and what makes §6.6's "resampled jointly
   across depths" meaningful). The flown tables are persisted as .npz;
   N0 flies the zero table SWEPT (keeps the binding axis; its bindings
   carry shot noise only); depth 0 is unswept (no injection layers
   exist there).
5. **Aux plan.** v7.2's §10 budgeted "~8 in-situ T1/T2* PUBs"; the
   flown plan (and §10 as revised) is 24 aux PUBs at
   4096 shots: four stations (after steps 9/27/45/60), each flying the
   IDENTICAL set {T1 at 0.5× and 1.5× the deep-circuit duration,
   fixed-delay Ramsey canary, Hahn echo at the same total delay,
   CAL0/CAL1}. One delay per station would confound drift with delay
   (the staircase killer would be invisible); the Ramsey absolute value
   is detuning-scrambled and serves ONLY as a repeated-circuit drift
   canary, the echo is the fittable dephasing face. Additionally one
   CAL pair flies IN EVERY JOB (the 372 MB parameter payload chunks
   on depth-block boundaries under the runner's 48 MB-per-job cap,
   both recorded in every hardware payload; the job count is
   payload-determined at submit, 8 by size alone and up to ~11 with
   the boundary granularity, so the per-job CAL PUBs number
   2·(jobs − 1) = 14-20 (7-10 PAIRS; two PUBs per pair, round-5
   unit repair) after the head pair's job 0; the §10 budget
   uses the 20 top end, and the shift across the 14-20 PUB range is
   0.13 min, both ends under cap (round-7 unit fix). FOUR CODE-PINNED
   CONSTANTS SIT UNDER THAT CHAIN and carry no manifest key, named
   here because they are flight-determining and an executor cannot
   see them from this document otherwise (round 11): the per-job
   payload cap (48 MB), the per-job bound-circuit cap (80,000), the
   billing anchor (0.327 ms/shot) and the station placement (after
   depths 9/27/45/60). They are code constants rather than
   registered thresholds because no verdict consumes them. A FIFTH
   constant used to sit in that sentence and did not belong there
   (round-12 blocker): the KEPT-COUNT FLOOR of 50 is consumed by a
   verdict, and by the harshest route in the document, since a cell
   under it is FLOORED, a floored cell makes the grid incomplete, and
   an incomplete grid VOIDS D-sign and W. It is a registered
   threshold now, manifest-keyed and sync-checked against the code
   like the billing cap, with its deep-end margin (§8a: ≥ 3× the
   floor at the worst modeled corner) as a machinery item, because
   the margin as registered was NOT met by the numbers this section
   was written against: at the then-current 0.6% guard bound the
   honest survival put the deepest cells at ~87 kept counts of 16384,
   which is 1.74×, not 3×. That was a real open decision and not a
   bookkeeping gap, it named three levers (the floor, the budget, the
   guard bound), and Amendment 1.8 took the third: at 0.5% the same
   arithmetic gives ~208 kept and 4.17×, so the requirement is met on
   the isolated axis and the layered one stays §9's. The margin's
   EVALUATOR is still a machinery item, and an unfrozen constant
   still blocks the flight. The
   first two set the job count, the job count sets the per-job CAL
   PUB count, and that sets the projection this section budgets, so
   a change to either cap re-shapes the flight without touching a
   frozen number. All four are recorded in every hardware payload,
   and the REALIZED job plan (job count, per-job PUB and
   bound-circuit counts, total payload) is printed by the DRY RUN,
   which runs before the addendum is written, so the executor
   carries that plan into the day-of addendum and checks it there
   against the 8-to-11 range this paragraph pre-registers; a head-of-batch-only CAL would invert the deep half
   of the grid with stale confusion), and the analysis POOLS the
   HEAD pair and every PER-JOB pair into the inversion, recording
   each pair separately as the
   drift record; the STATION CAL pairs stay OUT of the inversion by
   design, they are aux drift diagnostics riding the stations
   (round-4 alignment with the implemented routing)
   (per-job pairs at 4096 shots: at 16384 they alone
   would push the projection to the cap; pooled over head + per-job
   pairs the
   confusion precision stays ~0.05% per qubit). Budget: ~4.34M shots
   ≈ 23.7 QPU min at the 0.327 ms/shot anchor, cap 25.
6. **Billing brake.** Jobs are submitted SERIALLY inside the one
   Batch; before each submission the brake compares
   max(measured billed seconds, shots-based fallback) + the next job's
   projection against the cap and stops submitting when it would
   cross (runtime's job.usage() can return None; the fallback keeps
   the brake armed). BRAKE ⇒ FLIGHT VOID, registered (round 6): jobs
   run depth-ASCENDING, so a firing brake truncates the DEEP END,
   and any missing science cell already VOIDs D-sign and W outright
   (Amendment 1.2): there is no partial-grid verdict to salvage, so
   the registered rule is that a brake stop VOIDS the flight and
   PARKS the candidate, the RECORD documenting the spent budget and
   the anchor's failure (the alternative, depth-interleaved job
   chunking so a truncation loses whole arms instead, was
   considered and not taken: the payload chunks on depth-block
   boundaries and re-chunking would break that; the residual risk
   is priced by the pre-submit projection at 23.66 vs cap 25.0 and
   the ~3% anchor margin, §10, and the brake exists precisely so
   the failure mode is a parked flight, never an unbounded bill).
   The verdict LINE for a brake-truncated artifact is "VOID (grid
   incomplete; brake-truncated per Amendment 1.6)": one label, the
   hw record carrying the brake string and the analyze print naming
   the cause (round-8 unification). That label is the GRID-INCOMPLETE
   RANK'S name, not an exemption from §7's precedence (round-14
   precision): if a higher-ranked trigger fires beside the brake, the
   higher one leads and this one is appended INSIDE the parenthesis,
   in step 3's form: "VOID (N0-CLEAN; grid incomplete,
   brake-truncated per Amendment 1.6)" (round-19 bracketing fix; this
   Amendment and step 3 had the same line with the appended triggers
   on opposite sides of the parenthesis, and step 3 cited this
   Amendment as its reason for appending at all). The runner's
   analyze embeds the brake note inside the grid-incomplete line
   unconditionally, which is this Amendment's behaviour and not §7's;
   reconciling the two is the verdict-line composer's job and is on
   the machinery ledger.
7. **Demodulation record.** The C-arm ± channels demodulate at the
   mean of ALL THREE dyad frequencies (the committed gate's choice,
   kept verbatim; 0.0075 off the mixing pair's own mean), C′ at its
   own pair mean: an asymmetry inside the W statistic, carried into
   the θ_D/θ_W refreeze rather than silently symmetrized. The flown
   demodulation basis (Floquet modes, dyad frequencies, site→dyad
   map, preps) is persisted IN the artifact and governs the analysis.

8. **The p2 Class-1 guard bound, 0.6% → 0.5%** (Tom's decision,
   2026-08-17, pre-data; no science data exists). WHAT CHANGED: the
   Class-1 guard on the calibrated 2q error of the used edges, in §4,
   §9 and the parking rule of Open question 2, and with it the runner
   constant the guard reads. WHY: §8a registers a ≥ 3× deep-end margin
   over the kept-count floor at the worst modeled corner, and at 0.6%
   the design gave 1.74×, failing its own precondition. Round 12
   surfaced that and named three levers; all three were priced from
   committed numbers before the decision. Lowering the floor 50 → 29
   needs a G3 justification that 29 counts stay unbiased, which does
   not exist. More shots, priced as a UNIFORM increase of all 252
   science PUBs, needs 28,244 each (150/exp(−5.238)) for ~40 QPU min
   against the 25-min cap, so THAT form of it is dead. **THE PRICING
   OF THAT LEVER WAS WRONG BY A FACTOR ~14, AND THE CORRECTION IS
   RECORDED HERE RATHER THAN QUIETLY DROPPED** (2026-08-17, recomputed
   from below after the decision): the requirement binds on the DEEPEST
   cells only, which at the old bound is 36 of the 252 (depths 54, 57
   and 60, each × 4 arms × 3 preps), and raising only those, only by
   their shortfall, costs **+1.16 QPU min for a 24.86 total**, inside
   the cap. That lever was alive, not dead. It does not reopen the
   decision, and the two reasons are worth stating because they would
   decide a similar fork again: the guard costs no minutes at all
   against that lever's 1.16 and leaves 0.14 min of headroom under the
   cap instead of none, and it moves the flight's registered operating
   point ONTO a measured ladder rung, which no shot allocation can do,
   since at 0.6% the power at the bound stays an interpolation whatever
   the counts are. What the correction does change is the RESERVE: a
   targeted deep-end allocation is the cheapest answer if the day-of
   layer fidelity comes in near the bracket's bad end, and §8a carries
   it priced there. Tightening the
   guard costs nothing in budget and nothing in signal, since the
   signal exponent γ̄T·(10/3) ≈ 2.95 does not contain p2 (§10). WHAT
   IT BUYS, all of it arithmetic on §10's own chain
   kept = 16384·exp(−970·p2·f_leak): the worst corner moves 87 → 208
   kept and 1.74× → 4.17×; the guard lands ON the measured p2 ladder's
   first rung (power 2.7 at 40 reps, 2.9 at 200) instead of between
   rungs, which retires both the header's interpolated P(detect) ≈ 0.7
   and the measured point §8a owed AT the bound. WHAT IT COSTS: a
   stricter day-of parking condition, since fewer 6-qubit lines pass
   at 0.5% than at 0.6%. That is a device-availability price, not a
   physics one, and parking is cheap because it always happens BEFORE
   the spend. WHAT IT DOES NOT FIX, registered so no reader takes the
   4.17× for a clean bill: the bound is an ISOLATED-gate number, the
   flown circuit runs three rzz per half-layer, and the §9 bracket's
   layered 1.3-2× puts the effective error at 0.65-1.0%, where the
   margin falls to 1.12× and then to a floored deep end. The ≥ 3×
   requirement is exactly equivalent to an effective p2 ≤ 0.538%, so
   the owed day-of layer-fidelity gate is what actually holds this
   design's deep end, and its frozen bound (G1/G3) is where that
   number belongs.

## Amendment 2 (2026-08-18, pre-data; the refreeze, the floor, and four systematics measured beside them)

No science data exists. Nothing has been submitted, no QPU minute has
been spent, and the thresholds this amendment moves are the ones §8a
registered as awaiting the refreeze. It lands as a numbered pre-data
amendment because §7 pre-authorized exactly this route ("If the
refreeze returns a different θ_D, it supersedes 0.00253 as a NUMBERED
pre-data amendment (§13), never silently") and because the floor's
producer formula, registered in §7 round 11, does not survive contact
with its own reachability criterion. Committed at `36a8dc3` (this paragraph names the commit that landed
the amendment; the freeze commit that moves the constants manifest is a
separate later step and has not happened).

**READ §2.8 BEFORE ACTING ON ANY NUMBER HERE.** This amendment has had
one round of three empty reviews and carries their surviving findings
as OPEN ITEMS rather than as repairs. It is written up in that state
deliberately: the numbers below are what was measured, and the places
where the argument around them is not yet closed are named on the page.
The constants freeze, the OpenArcs entry, and the flight go are three
separate later steps, and nothing here performs any of them.

**WHAT THE REPO WAS ASKED, and what it returned.** The sweep recorded
under "What the repo already holds" (near the top of this document)
covered the corner beat as an object; this amendment's items are
different objects, so they were swept again on 2026-08-18 by two
store-sweep agents, by store.

- `docs/ANALYTICAL_FORMULAS.md`: F140, the R90 frozen divisor that pins
  −4γ̄ at every coupling, which is the object §2.4b is about; and F154,
  which carries the fence "exceptional couplings expected" twice,
  quoting F122's own. NOTHING on decision thresholds, on NaN or
  kept-count policy, and NOTHING on the T1 → γ conversion.
- `docs/proofs/`: `PROOF_R90_FROZEN_DIVISOR` §9 and §12, which own
  §2.4b and are quoted there. For §2.3, two SIBLINGS and neither is the
  object: `PROOF_ASYMPTOTIC_SECTOR_PROJECTION` proves the sector
  populations are CONSERVED under Z-dephasing (p_w(∞) = p_w(0)), the
  in-sector identity and the opposite of a magnon that leaves;
  `PROOF_DEPHASING_FRONT_RENEWAL` carries a renewal representation for
  the front, WITHIN the sector, caught and released. NOTHING on a
  magnon number that leaves and returns.
- `experiments/`, including the null results and the prior flights:
  `RECORD_PARITY_HARDWARE_PREDICTION.md` round 28 holds the
  lower-envelope rule §2.2 uses, verbatim, inside an enumerated
  protection-interaction class whose instances are numbered and not
  closed (the seventh is that round, the eighth two rounds later);
  `STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md` is this repository's
  template for a numbered pre-data amendment; `IBM_F129_RAMSEY_FRINGE.md`
  and `IBM_CONCENTRATOR_RELOADED.md` hold the governing prior rule for
  §2.5, that any failed fit, NaN or guard trip is an instrument failure
  and never a verdict. NOTHING on re-entry.
- The OpenArcs registry: NO `corner_beat` entry, which the freeze commit
  opens. `gamma_book_enforced_nowhere`, opened 2026-08-18, governs §2.6
  and is answered there.
- `fw.Confirmations` and its C# mirror, 24 entries in both: NOTHING. No
  entry concerns a decision threshold or the T1 conversion.
- `docs/GLOSSARY.md`: the γ conversion table, from which §2.6's
  arithmetic follows in one step (a `D[Z]` channel at rate γ decays
  coherences at 2γ, so γ = 1/(2T₂); the Γ/4 form is derivable from the
  table's γ_Z row and is not written there). NOTHING on the corner beat.
- `docs/CAUGHT_ERRORS.md` (machine-local, so a clone cannot follow this
  citation): the rate-book confusion, logged twice.
- One store the earlier sweep never named, named here because §2.4b
  would otherwise collide with it:
  `experiments/THE_EXCEPTIONAL_COUPLINGS.md`. It is a DIFFERENT object
  and says so itself, pointing out that `PROOF_R90_FROZEN_DIVISOR` §9
  "spends the word exceptional on a different object". §2.4b carries the
  return fence.

**TWO ADDRESSES THIS SWEEP CORRECTED**, recorded so the next sweep does
not walk them again: `STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md` has no
§39 and no numbered sections at all; and MirrorWorld has no `Rung`
object, the renewal object is `Renewal.cs` and the Krawtchouk identity
lives in Core's `PopcountCoherencePi2Odd.cs`.

**AND ONE "CORRECTION" WITHDRAWN, which is worth more than the
correction would have been.** A draft of this amendment asserted that
F154 carries no "exceptional couplings expected" fence and that the
phrase belongs to `PROOF_R90_FROZEN_DIVISOR` §9. Both halves are false:
F154 carries it twice and the proof carries it zero times. The store the
sweep reported as silent was the store that was speaking, and the error
entered while correcting a note that had been right.

**AND ONE DIGIT COLLISION, fenced.** The value 0.00253 that §2.1
supersedes is this document's θ_D. The same digits appear in
`experiments/PRICE_PAIR_HARDWARE_PREDICTION.md` as a correlator
c₀₂ = +0.00253 ± 0.00049, from the flight of 2026-07-04, and that
reading is itself recorded there as SUPERSEDED, a non-recurring outlier
that did not repeat on the same line in run 3. Coincidence of digits,
different object, and the neighbour's own value is not a standing
measurement.

### 2.1 θ_D, refrozen through the analyze-side chain

**θ_D = 0.00306, superseding 0.00253**, and θ_W = 0.00366 beside it.
Both are frozen by the run of 2026-08-18 22:13, whose record is
committed at
`simulations/results/corner_beat/corner_beat_refreeze_20260818_221318.json`
and summarised in §2.7. Produced as §8a registers the package: the runner's own `fit_rate`
on the ANALYZE side, 500 replicates per corner with a held-out 500 more,
the four corners of the 2 × 2 bracket (f_leak at 8/15 and at 0.9, the
hop fraction at 0 and ½), and the fit axis on `realized_dose` on the
runner and the gate side together. θ_D and θ_W are UPPER cuts and take
the UPPER envelope over those four; the floor of §2.2 is a LOWER cut and
takes the lower one. The H0 mean is stated explicitly rather than
absorbed into the 3σ̂, per §7's round-4 unification: E[d | H0] = +0.00002
with sd 0.00022 at the corner quoted in §2.7.

**WHAT "CONSERVATIVE" MEANS HERE, and it needs saying because this
document has been bitten by the other direction.** A larger θ_D is
conservative against a false DETECTION. It is not conservative
everywhere: ¬D-sign routes to FALSIFIED, Anti-D or INCONCLUSIVE, so a
larger θ_D also makes ¬D-sign easier and raises the false-FALSIFIED
rate. That is the same sign trap §7's round-15 repair caught for θ_F and
κ. The direction is named rather than assumed, and the two-sided
consequence is an open item (§2.8).

**Sample size, registered rather than assumed.** Below **n = 299** a
0.01 ceiling is not demonstrable **at 95% one-sided confidence** by a
campaign that observes zero events: the Clopper-Pearson upper bound
1 − α^(1/n) at α = 0.05 first falls to 0.01 at n = 299 (0.010002 at
n = 298, 0.009969 at n = 299). The confidence level is part of the
claim; at 90% the same bound is reached at n = 230. 500 was chosen above
299, not below it, and 299 is the best case, since a corner that
observes trips needs far more.

### 2.2 The s²(C) floor: the producer formula is superseded, and the floor is inert on D-sign

§7 round 11 registered **floor = E[s²(C) | H0] + 3·σ̂(s²(C) | H0)**, and
round 14 registered the reachability criterion
**P(s²(C) < floor | H1) ≤ 0.01 at BOTH bracket ends**. The refreeze ran
both. The floor freezes at **0.00067**, the LOWER envelope over the four
bracket corners, with the set's maximum 0.00385 kept beside it in the
record as `s2C_floor_max_superseded` so the discarded end stays readable.

**IT PASSES REACHABILITY AND IT CHANGES NOTHING.** Held out, at every
one of the four corners, 0 of 500 H1 replicates fall below it, a
Clopper-Pearson 95% upper bound of 0.0060 against §7's 0.01. And it
removes 0 of 2000 held-out H0 draws and 0 of 2000 held-out H1 draws: at
every corner the detection rate under the CONJUNCTION equals the
detection rate on d alone, to the last digit (§2.7).

**WHY THAT IS FORCED, and not a coincidence of this sample.** Since
round 9 the operative rule is the CONJUNCTION d > θ_D **and**
s²(C) ≥ floor, with d = s²(C) − s²(U). A replicate can therefore pass
the first conjunct and die at the floor only when

    s²(U) + θ_D < s²(C) < a + b·s²(U),

that is, only when `a − θ_D > (1 − b)·s²(U)`. Since s²(U) is a sample
variance and therefore non-negative, any floor with `b ≤ 1` needs
**a > θ_D** to bite at all. The frozen floor is 0.00067 against
θ_D = 0.00306, so **no replicate can pass d > θ_D and fail this floor**,
on any sample, at any corner. The harness carries the lemma as
`is_provably_inert`. The argument is pathwise: the event set is empty
for every realization, so nothing about the joint distribution of
s²(C) and s²(U) enters, and their correlation is irrelevant to it.

**THREE HYPOTHESES THE LEMMA RESTS ON, stated as hypotheses because two
of them could be broken later by a change that looks unrelated.**
(i) d is literally s²(C) − s²(U) with the SAME s²(U) the floor reads,
which is true in the harness and would break silently if a replicate's
d and its floor ever saw different exclusion sets. (ii) s²(U) ≥ 0
pathwise, which is true of a raw sample variance and would fail for any
future debiasing that subtracts an estimated noise floor. (iii) b ≤ 1.
The registered coefficient grid runs to b = 1.50, so for the two members
above 1 the inertness is MEASURED by the 2000-replicate sweep and not
proved.

**AND THE SCOPE IS THE LEMMA'S LIMIT, named here rather than left to be
found.** §7 round 8 gives the floor a SECOND forcing: a trip forces ¬W
regardless of the computed exceedance. W's rule is
s²(C) − s²(C′) > θ_W, in which s²(U) does not appear at all, so the
window this lemma empties has no bearing on it. The floor is inert on
D-sign and **unaudited on W**, and the harness cannot see it either:
`conjunction_rate` and `incremental_guard_counts` take only
(d, s²(C), s²(U)). This is an open item of the freeze (§2.8), not a
result of it.

**WHERE THE LOWER-ENVELOPE RULE COMES FROM, and why the borrowing is
legitimate.** `experiments/RECORD_PARITY_HARDWARE_PREDICTION.md` round
28 reached it the expensive way: a cut frozen at the single worst
admitted basis "terminally VOIDed the BEST admitted device", VOID in 269
of 300 runs, counted there as the seventh instance of an enumerated
protection-interaction class, and the rule extracted was that **for a
LOWER cut the worst basis IS the lower envelope**. The index set matters
and is the same here: the envelope is taken over the four bracket
corners, which are NUISANCE states the device could be in and whose
f_leak and hop the flight cannot know, not over a grid of coefficients
the analyst chose. Taking the minimum over unknown device states is
protection; taking a minimum over one's own candidate knobs would not be,
and is not what happens. Note the direction is opposite for the upper
cuts by the same logic: θ_D and θ_W take the MAXIMUM over the same four.

**WHAT THE NEIGHBOUR DOES AND DOES NOT SAY.** RECORD_PARITY verified
"no floor trips post-fix", an empirical absence of trips on its own
sample. That is strictly weaker than inertness, which is the statement
that no verdict CAN change. The proof above is this document's own, and
the neighbour supplies the freezing rule, not the inertness.

§7's two registered consequences survive unchanged: the floor is not an
anomaly detector, and a floor trip on a clean dose with a clean N0 routes
to DATA and never to VOID. Note that under a true null the floor trips
by design and does so often: §2.7 records 69 to 497 trips per 500 H0
replicates across the four corners, which is §7's consequence (i)
measured, not an anomaly.

**A PROCESS NOTE, because the shape is a repeat offender.** An earlier
pass scanned the (a, b) coefficient plane AFTER seeing freeze data and
picked a pair with margin. Nothing was frozen, so it cost only time, but
a replacement gate that is also a number chosen after seeing the data is
exactly the researcher degree of freedom round 11 closed. The
lower-envelope rule takes no such choice.

### 2.3 The re-entrant pedestal: a named systematic, measured, and NOT modelled in the frozen chain

§10 models sector loss as a per-gate survival, `exp(−970·p2·f_leak)` at
the deep end: a shot that leaves the one-magnon sector is gone. That
expression is exact for independent per-gate errors and is not a
one-error approximation; what it omits is RE-ENTRY. A magnon number that
random-walks under a depolarizing channel walks BACK, and post-selection
keeps it, because post-selection reads the FINAL magnon number and not
the history. At the registered p2 bound the deep end carries 4.85
expected errors (§10), so the walk has room. The words "re-enter",
"Markov" and "multi-error" appear nowhere else in this document; this
item is what the document did not know.

**THE SIZE, from an exact Markov chain on the magnon number, driven by
the same depolarizing model §10's f_leak count comes from. The chain
runs at the q-weighted f_leak ≈ 0.713**, the prepared states' own bond
occupancy, and NOT at the 8/15 end that sets θ_D; the two tables in this
item are at different f_leak and that is why the shares and the costs
cannot be read against each other line by line:

    depth  3:  kept 0.841   never left 0.839   re-entrant  0.3%
    depth 24:  kept 0.325   never left 0.245   re-entrant 24.5%
    depth 36:  kept 0.228   never left 0.121   re-entrant 46.6%
    depth 60:  kept 0.147   never left 0.030   re-entrant 79.8%

An independent cold reviewer reproduced the chain on an exact
64-dimensional density-matrix construction and agreed to 0.5% (share
0.776 against 0.798, the chain slightly high because it takes one error
per step rather than sixteen gates at 0.005).

**AND THE TABLE IS ABOUT THE MAGNON NUMBER, NOT ABOUT THE BEAT.** Of
the seven in-sector Paulis the §5 count leaves, four HOP and three are
Z-type, so none of the seven preserves the beat undamaged. "Never left"
is therefore an upper bound on "carries the beat", and the deepest
never-left fraction of 0.030 is larger than the no-error-at-all
fraction. The cost table below is not optimistic for this reason: the
committed gate models the in-sector residue explicitly (half Z-like
dephasing at the Strang gate-count profile, half depolarization toward
Id/6). What the caveat costs is the headline: "four fifths of the
deepest kept shots carry no beat" is the magnon-number reading, and the
beat-carrying fraction is smaller still.

**WHERE THEY REACH THE ESTIMATOR, and it is exactly one channel.** A
uniform release is ρ = I/6, which is the exact stationary state of a
unital map (Trotter unitary plus diagonal dephasing) and therefore has
identically zero component on every decaying mode; the earlier reason
given here, that the Floquet modes are orthonormal, is not available,
since the block is non-normal. So a uniform release projects to exactly
zero on all three registered dyads, and re-entrants are invisible in the
NUMERATOR and present in the DENOMINATOR: the dyad amplitude per kept
shot falls by (1 − share(t)). "Pedestal" names where they sit, not what
they do. TWO FENCES ON THAT EXACTNESS. It is exact for exactly I/6,
while the measured re-entrant vector is within 2% of I/6 at depth 24 and
further off at shallow depths where the share is small; the product of
deviation and share is the residual numerator contamination and nobody
has bounded it (§2.8). And the loss is not a constant added rate: it
grows with depth, so it distorts the SHAPE, which is why r̄ rises while
s² falls, where a constant added rate would leave the variance of the
three rates alone.

**WHAT IT COSTS, measured under the FULL chain, 400 reps, paired, at the
corner that sets θ_D (f_leak = 8/15), three ways (data clean / data
contaminated / data contaminated and fit corrected):**

    θ_D        0.00287   0.00287   0.00327
    d          0.00512   0.00434   0.00502
    r̄          0.50097   0.55866   0.49994
    s²(C)      0.00733   0.00647   0.00749
    P(detect)  0.907     0.807     0.810   (each at its OWN threshold)

The corrected column is ×1.142 on θ_D and ×0.981 on d against the CLEAN
column; the 0.003 the correction buys is against the CONTAMINATED one.
An H1 signal-to-noise ratio d/sd(d | H1) reads 2.93 / 2.55 / 2.53 across
the three and is NOT the ratio the thresholds are built from, which is
an H0 scale; it is quoted here because it was quoted before and because
it points the OTHER way on the correction, which is exactly why it must
not be read as the decision margin.

**The signal comes back and the discrimination does not.** Folding the
known envelope (1 − f(t)) into the fit model beside the fitted
exponential, with no free parameter, returns d to 98% of the clean value
and closes 87% of the contamination deficit, and returns r̄ and s²(C)
essentially exactly. But θ_D RISES 14%, because the correction acts
under H0 too: there the signal is near zero and dividing by an envelope
that falls to 0.202 at the deep end amplifies the null's own noise, so
sd(d | H0) grows and the threshold walks up with it. **Operationally the
correction buys 0.003 of P(detect).**

**THEREFORE, REGISTERED HERE AS A SYSTEMATIC RATHER THAN REPAIRED:** the
frozen chain does NOT model re-entry, and the flight's realistic
detection probability is about ten points below whatever the clean chain
returns at the registered thresholds. Against §2.7's worst-end 0.890
that is **approximately 0.79**. The design REGISTERED 0.920 (§8a and the
header); 0.907 is this amendment's own clean A/B column on a different
chain and is not the registered figure. `env_known` on
`fit_rate`/`arm_rates` and `REENTRY_FIT` on the harness stay in the code
as the way any successor RE-MEASURES this; with them off the output is
bit-identical to before they existed (verified, and an all-ones envelope
reproduces the committed call exactly).

**WHAT THIS DOES NOT MOVE.** §10's per-gate numbers are per-gate numbers
and stand as written: the 4.85 exponent, and the 1.3 to 7.5%
post-selected survival band at the registered p2 guard bound. What
changes is the READING of "kept": the kept fraction exceeds the
never-left fraction, and only part of the never-left fraction carries an
undamaged beat. The H0 side is immune, and the mechanism is not the one
an earlier pass asserted: nothing "cancels in d", the pedestal never
enters either estimator's projection at all. The null is immune because
the compression is proportional to the RATE SPREAD and under H0 the
three rates are nearly equal. Measured on the table above, the H0 mean
moves 1.9 × 10⁻⁵ against a signal loss of 7.8 × 10⁻⁴, about forty times
smaller. `dmag_band_C` being null in the manifest is what keeps the move
in s²(C) from mattering here, and that is luck rather than design.

**A METHOD NOTE THAT COST A DAY AND IS WORTH THE LINE.** The correction
was first reported as recovering 99.5%, measured in a NOISELESS
isolation where no variance exists and power could not be measured. The
reviewer who produced it fenced it himself; the fence was written down
and the headline was passed on anyway. **A recovery number measured
without noise says nothing about a decision that is a signal-to-noise
ratio.**

### 2.4 One negative, one located finding, and three dead levers

**(a) THE DEPTH GRID: it stays, and the reason is weaker than it was
written.** The axis §12 of this document lists as unsearched, priced at
250 replicates at the θ_D-setting corner, clean chain and with the
pedestal:

    end 36, 13 pts, 13.93 QPU min:  P(detect) 0.856 / 0.780
    end 48, 17 pts, 18.22 QPU min:  P(detect) 0.880 / 0.828
    end 60, 21 pts, 22.50 QPU min:  P(detect) 0.907 / 0.807

On the CLEAN column the lever arm is monotone: 60 → 36 costs about five
points to free 8.6 QPU minutes. **On the pedestal column it is not.**
End 48 reads 0.828 against end 60's 0.807 and saves 4.3 QPU minutes, and
§2.3 makes the pedestal column the realistic one. At 250 replicates that
0.021 is inside the sampling spread of these estimates, so the inversion
is not evidence for shortening; it is also not evidence against it, and
the earlier reading of this table as "real and monotone" was reading the
clean column only. The grid stays at 21 points to depth 60 because
nothing here justifies moving a frozen flight parameter, not because the
lever was shown to be absent. Two questions are registered untried: the
inversion at a replicate count that could resolve it, and the 2-D
question of the fit correction (which makes the deepest points noise
amplifiers) combined with a shorter grid.

**(b) THE CORNER BLOCK'S EXCEPTIONAL COUPLINGS EXIST, ARE LOCATED, AND
THE FLOWN POINT IS CLEAR OF THEM. This RETRACTS the "clean negative"
an earlier pass recorded.** In this item, §6, §9 and §12 are sections of
`docs/proofs/PROOF_R90_FROZEN_DIVISOR`, not of this document.

§9 records that on a NON-uniform locus profile the frozen divisor gains
an extra algebraic dimension at isolated NONZERO real couplings, "one
small coupling away from anywhere", and §12 leaves their COUNT open,
reporting that at N = 6 two generic profiles on the same locus give four
real couplings and eight. It is the same block: C(6,1)² = 36 = N² IS the
corner beat's room, and the flown lighting profiles are exactly the
non-uniform locus profiles the statement is about (both verified on the
locus here: 2·avg − reverse returns the profile).

An earlier pass reported "no exceptional couplings, clean negative"
after seeing multiplicity jumps only at Q ≤ 0.25 and reading them as the
J → 0 root §9 excludes. **§9 excludes J = 0 EXACTLY** ("Exceptional
coupling below always means a nonzero root"), not small J, so that
reading was wrong twice over: it dismissed nonzero couplings as the zero
one, and the couplings are not where it put them. Recomputed in exact
arithmetic (mod-p characteristic polynomials over eight verified primes,
lifted by CRT, then Sturm over ℤ; no floating point at any step, since a
rank test at a tuned degeneracy is the least trustworthy measurement
there is):

- The generic algebraic multiplicity of the frozen root −4γ̄ is exactly
  ⌊N/2⌋ = 3 and **semisimple**, for C, C′ and the uniform control alike.
  F140 reproduces, which validates the construction.
- The cofactor whose vanishing marks an exceptional coupling has degree
  **30 = N(N−1)** in the coupling and is EVEN, a polynomial in J², both
  exactly as §6 and §9 require.
- **C, the maximizing arm: exactly ONE real nonzero exceptional
  coupling, in closed form, Q\* = √2 = 1.414214.**
- **C′, the non-maximizing arm: exactly ONE, at Q\* = 4.647043**, the
  single real root of 27w³ + 594w² + 288w + 1120 in w = −Q².
- **The uniform control has NONE.** This is an independent finding and
  not a reading of §9, whose statement (non-uniform implies at least
  one) says nothing about the uniform case.
- **At each of the two, geometric 3 against algebraic 4: one Jordan
  block of size exactly two**, with 3, 3, 3, 3 at neighbouring real
  couplings and at the flown point. That is §9's law, and §9 measured it
  at N = 3 and N = 4 while §12 leaves the block size at larger N open on
  this stratum, so this is a new data point for that question at N = 6.
- Counting convention, because §12 counts differently: the count above
  is of POSITIVE Q, each standing for the pair ±J. Both flown profiles
  are far from generic (three of six sites dark), which is the honest
  reason they carry one where §12's generic profiles carry four and
  eight, and it is a data point for §12's question rather than a
  contradiction of it.
- Model sensitivity: with the ZZ term removed (XY rather than
  Heisenberg) the couplings MOVE, to Q\* = 1.872530 for C and
  Q\* = 3.104212 for C′, and stay well below the flown point. The BLOCK
  also grows there: the XY defect is one Jordan block of size THREE at
  both profiles, geometric 3 against algebraic 5. §9 had seen a size-3
  block only on the zero-mean stratum, at N = 3, so the taxed stratum
  reaches it too once the ZZ term is removed. The flown model is
  Heisenberg and its block is size two; this is a sensitivity reading,
  not a statement about the flight.

**WHAT THIS MEANS FOR THE FLIGHT.** The flown point is Q = 10. The
nearest exceptional coupling is C′'s at 4.647, a distance of 5.35 in Q,
and C's sits at 1.414, 8.59 away. The corrected statement is stronger
than the negative it replaces, because it is a LOCATION rather than an
absence: the defect exists, it is where the proof says such defects are,
and the working point is not there. The conditioning at the flown point
is benign and is quoted as the right quantity for a non-normal operator:
the eigenvalue condition numbers of the frozen modes read 1.10 to 1.69,
where a defect would drive them up. The separation of 0.67 between the
frozen root and the nearest other eigenvalue is also measured, but it is
essentially profile-independent (the uniform control gives the same
0.667), so it is a property of the block's generic structure and cannot
carry a statement about this design's lighting.

**THREE LIMITATIONS, named rather than hidden.** This is the CONTINUOUS
Liouvillian, the object `PROOF_R90_FROZEN_DIVISOR` and F140 are about,
while the flown dynamics is Strang-Trotterized and its Floquet
generator's exceptional set has not been computed. The result is about
the frozen root's own Jordan structure, and a review round measured that
the frozen modes carry weight 4 × 10⁻¹⁵ on the population cells, which
if it survives scrutiny would mean this defect has no channel at all to
an estimator read from populations, a stronger statement than distance;
that reading is NOT yet established here, because the "dyads" §2.3 names
are Floquet-mode dyads while the probe used site pairs, and the two have
not been shown to be the same statement (§2.8). And the exact-arithmetic
route was itself repaired mid-computation: the first attempt used
GUESSED moduli, two of which were composite, so every CRT lift was
silently wrong and returned a degree of 40 where the a-priori bound is
33. The contradiction with the bound is what exposed it, not the
plausibility of the output.

**AND THE RETURN FENCE, which the repo owed in one direction only.**
`experiments/THE_EXCEPTIONAL_COUPLINGS.md` is a DIFFERENT object and
already says so: that document is the ceiling arc's, on the UNIFORM
profile and the block (2,2). This item is §9's object, on the
non-uniform locus profiles and the corner block. Neither is the Trotter
EP §4 of THIS document carries, which is a third object again: the arms'
own EPs at Q ≈ 30 (C at J·dt = 0.15) and Q ≈ 21 (C′), the binding margin
at the flown Q = 10 being C′'s factor of two. Three objects, one word,
and the flown point is clear of all three.

The verification is `simulations/corner_beat_exceptional_verify.py`.

**THREE DEAD LEVERS, recorded so nobody re-derives them.**

- **CAL shots: FLAT.** 45056 → 90112 → 180224 moves θ_D by ×1.001 and
  ×1.002. The estimator is NOT calibration-limited.
- **Science shots: SATURATING at the operating point.** Measured
  exponents per doubling, where 1/N would read +1 and 1/√N +0.5:
  8192 → 16384 gives θ_D +0.56 and sd(s²(C) | H1) +0.25;
  16384 → 32768 gives +0.22 and +0.18. Reaching a materially different
  threshold by shots alone needs about twice the science budget, roughly
  47 QPU min against a 25 min hard abort. Not viable. `--shots` exists
  in the harness for exploration only, and the house rule "measure, do
  not scale" is why these were measured; the 47 is itself an
  extrapolation of a local exponent on a saturating curve, and it errs
  toward "not viable", which is the safe direction.
- **M_BIND: the one axis that is NOT saturated, and it is frozen flight
  configuration.** 16384 shots per cell are M_BIND = 1024 lightings × 16
  questions each. The engineered γ is injected as one parameterized RZ
  per site per step, a virtual frame change on Heron: no duration, no
  pulse. The light exists only in the AVERAGE over the 1024 draws, so
  the measured saturation says the estimator is limited by the number of
  LIGHTINGS and not by how often each is asked. Freed minutes would buy
  bindings, not shots. UNTESTED, and out of scope here because M_BIND is
  frozen.

**AND ONE NARROWING THAT IS REFUTED, so that it is not attempted
again.** The §5 and §10 Pauli count sorts the survivors as well as the
leavers, so the same count fixes the hop fraction: re-derived exactly on
the six-site sector and re-verified at operator level by a reviewer, an
excitation ON the acted bond has 80 of 150 leave (f_leak = 8/15) and 40
of the 70 left move (hop = 4/7), while an excitation OFF the bond has
240 of 300 leave (f_leak = 4/5, hop = 0), giving
`f_leak = (12 − 4q)/15` and `hop = 4q/(3 + 4q)` in the bond occupancy q,
which is the prepared state's own (min 0.2363, mean 0.3280 over every
depth and preparation). The closed form is correct and is kept.
**Narrowing the flown bracket with it is refuted three ways:** a GENERAL
Pauli channel has three free directions rather than the depolarizing
3:4:8, and setting the XX/YY weight to zero, the smallest component on
real CZ gates, makes (f_leak = 8/15, hop = 0) exactly admissible; a
coherent ZZ over-rotation does the same with no Pauli assumption at all;
and T1 and |2⟩ leakage remove a magnon with probability 1 while not
being bond errors, which drops the claimed f_leak cap. Structurally the
harness's `hop` is a SWAP on a uniformly random bond while the physical
post-selected move is a PROJECTION onto the bond (L1 distance 0.612 at
depth 30), so there is no q in the harness to narrow. **Narrowing a
model-free bracket with a noise model, on a paid one-shot, at the corner
that sets the threshold, is the trade not to make.**

### 2.5 Fit-health margins, and what stays unchanged

The §8a package freezes the fit-health margins beside θ_D, and they
freeze here: `HOLE: FIT_HEALTH_MARGINS` (the δω excursion margin from
the IDEAL construction, the channel-fit residual margin, and the rate
saturation margin, whose BASE must be registered as one of the two
strings §8a admits rather than left to the code, since the runner
hard-aborts on anything else). The two-route split §7 registers is
unchanged: a NaN verdict statistic on a COMPLETE grid is a dead channel
and reads fit-health VOID, while a cell under the kept-count floor of 50
is grid-incompleteness and takes the single round-4 route. The governing
prior rule from the two IBM flights stands: any failed fit, NaN or guard
trip is an instrument failure and never a verdict.

### 2.6 The T1-CLEAN scale face is computed, not transcribed, and "scale" means the MEAN

§9 gives the bridge γ_l^hw = τ_step/(4·T1_l·dt), the Γ_l/4-equivalent
dephasing in the simulation's J units, and the SCALE face compares the
profile against ≤ 0.25 γ̄. Until 2026-08-18 no line in either repository
computed that conversion: the runner enforced that `tau_step_us` was
PRESENT in the day-of addendum and never read it again, so on flight day
an executor would have done the arithmetic by hand, from prose, across
two repositories, once. It is now code (`t1_clean_profile`,
`t1_list_from_snapshot`, `t1_clean_scale_face`, built test-first with 14
tests inside the flight suite's 135).

**ONE STANDING SENTENCE OF §9 IS SUPERSEDED HERE, by name.** §9 states
the scale face as "site scatter ≤ the G4 scale", and §8a as "scattered
Γ_l/4-equivalent profile scale ≤ 0.25 γ̄". The quantity the committed
gate actually draws and calls scale is the profile MEAN in units of γ̄:
the draw is `uniform(0, 2·scale·γ̄)`, whose mean is `scale·γ̄`, and the
gate's own comment says "'scale' = Gam_mean/4 in units gbar". **The
face is a MEAN, not a spread**, and a standard deviation is the natural
guess and would have been wrong while looking right. Since T1-CLEAN is a
GLOBAL VOID trigger, two readings of it on one page is not a wording
problem, which is why this supersedes rather than clarifies.

Two further properties are registered so a successor does not silently
undo them. **It PRINTS, it does not abort**: T1-CLEAN is a registered
GLOBAL VOID trigger at rank 3 of §7's void order, and a hard exit there
would convert a registered void into instrument failure, which this
document forbids in its own words elsewhere. **The book is Lindblad**,
Γ/4 and not Γ/2, pinned by a test with the reason in the docstring, so a
later simplification of the 4 has to argue with a test.

**MEASURED, and it is the day-of expectation this document owed:** at
realistic Kingston T1 (255 to 327 µs) the scale face reads **0.029 γ̄
against the registered 0.25**, about an order of magnitude of margin,
which matches §9's own note that the bound is close to non-binding as
scatter. A degraded chain at T1 = 25 µs correctly reads 0.333 and fails.

**AND IT ANSWERS THE ARC OPENED THE SAME DAY.** The OpenArc
`gamma_book_enforced_nowhere` records that a dephasing rate is written in
two books a factor of two apart, that several flights have been compared
across the seam, and that its first parked step is that a NEW flight must
not land without naming its book. The corner beat is in the Lindblad book
throughout (`sigma = sqrt(4*gamma*dt)` in both the runner and the gate),
and the T1-CLEAN bridge was the one place it could still have acquired
the other one, which is why it names its book in the docstring and pins
it with a test. That is the compliant behaviour the arc asks for, done
once by hand for one field; the arc's own remedy, a gate that refuses an
unnamed book, is not this amendment's work.

### 2.7 The freeze record, as measured

`simulations/results/corner_beat/corner_beat_refreeze_20260818_221318.json`
(committed with this amendment), 500 replicates per
corner plus a held-out 500, four corners, fit axis `realized_dose` on
both sides, re-entry NOT modelled per §2.3:

    theta_D    0.00306      (upper envelope over the four corners)
    theta_W    0.00366      (upper envelope)
    s2C floor  0.00067      (LOWER envelope; max superseded 0.00385)

    P(detect) at the REGISTERED thresholds, held out, on d alone and
    under the conjunction. They agree at every corner, which is §2.2's
    inertness measured rather than argued:

      f_leak 8/15, hop 1/2   0.890    H0 false 0.002   floor trips 69/500 on H0
      f_leak 8/15, hop 0     0.958    H0 false 0.000   floor trips 455/500
      f_leak 0.9,  hop 1/2   0.994    H0 false 0.000   floor trips 488/500
      f_leak 0.9,  hop 0     0.986    H0 false 0.000   floor trips 497/500

    floor reachability, held out: 0 of 500 at EVERY corner, 95% upper
    bound 0.0060 against §7's 0.01; the floor removes 0 of 2000 H0
    draws and 0 of 2000 H1 draws.

**THE WORST END IS 0.890, and it is the number the flight should be read
against**: not the 0.920 of the earlier freeze, and not §2.3's clean A/B
column of 0.907. The re-entry systematic of §2.3 sits on top of it,
giving approximately 0.79 as the realistic expectation.

### 2.8 What this amendment leaves open

Written as a list because each item is a piece of work, not a caveat.
The first three bear on numbers that appear above.

1. **The floor's W forcing is unaudited.** §2.2's inertness is proved
   for D-sign only, and a floor trip also forces ¬W by §7 round 8. The
   harness's guard counters do not take s²(C′) and so cannot measure it.
2. **The reachability criterion on the page is not the one in the
   code.** §7 states P(s²(C) < floor | H1) ≤ 0.01 as a rate; the harness
   enforces a Clopper-Pearson 95% UPPER BOUND on that rate, which is the
   stricter and better object. The α belongs on the page, and §2.1's
   n = 299 already depends on it.
3. **Two freeze-grade runs disagree.** This freeze returns
   θ_D = 0.00306 where the run of 2026-08-18 05:14 returned 0.00291 on
   what was believed to be the same configuration. Small, and
   unexplained; an unexplained difference between two freeze-grade runs
   is a finding about the harness rather than a rounding, and it is not
   resolved here.
4. **The earlier existence-condition argument is superseded and not
   re-derived.** A previous pass argued from
   P(s²(C) < θ_D | H1) = 0.0000 / 0.0140 / 0.0040 / 0.0060, measured on
   the run of 04:02, that no content-bearing floor is reachable at all.
   That is a different quantity at a different cut from what this freeze
   reports, and the general claim is not carried by this amendment.
   What is carried is the specific one: THIS floor, at 0.00067, is inert
   on D-sign because 0.00067 ≤ θ_D.
5. **The numerator contamination is unbounded.** §2.3's exact-zero
   projection holds for exactly I/6; the measured re-entrant vector is
   within 2% of it at depth 24 and further off where the share is small.
   The product has not been bounded.
6. **The population-cell reading of §2.4b is not established.** If the
   frozen modes really carry no weight on the observable, the
   exceptional-coupling result is about an object the flight cannot see,
   which is stronger than distance. Site pairs are not Floquet-mode
   dyads and the two have not been reconciled.
7. **The depth-grid inversion is unresolved** at a replicate count that
   could resolve it (§2.4a).
8. **θ_D's two-sided consequence.** A larger θ_D lowers the
   false-detection rate and raises the false-FALSIFIED rate; only the
   first direction is measured here.
9. **Inside the same §8a package and not frozen by this amendment:**
   θ_F, κ, R_MAX_FIT's freeze, the depth-0 under-coverage binary that
   §7 registers as a choice the refreeze must make, and the δω-flag
   firing rate under H0 and H1. θ_W IS frozen here, at 0.00366.

## Day-of addendum (empty until flight day)

Filled by the §9 re-gate BEFORE the Batch opens: day-of governing band
widths, calibration snapshot, line selection with scores, Class-1 guard
outcomes, and τ_step (the transpiled per-Strang-step wall duration,
the T1-CLEAN µs→γ̄ bridge's shared conversion). The runner's schema
gate enforces the keys calibration_snapshot, band_width_factors,
tau_step_us, and job_plan (round-7 alignment: an addendum filled from
this section alone must not abort on flight day; job_plan added in
round 14, which found Amendment 1.5 instructing the executor to carry
the dry run's realized job plan into an addendum that had no key for
it, while the dry run printed only the job COUNT and kept the per-job
PUB and bound-circuit numbers to itself).

ONE ADDENDUM IN TWO FORMS, and the executor owes both (round 11: the
document registered an addendum committed HERE while the runner
hard-aborts on a JSON FILE it never named, so an executor reading
only this section would have filled the prose, passed the commit
gate, and been stopped at submission by a file whose name, location
and format appear nowhere in the pre-registration). The machine
form is `corner_beat_dayof_addendum_YYYYMMDD.json`, dated for the
flight day, beside the runner in the pipeline directory, carrying
the three enforced keys above; the runner refuses to submit without
it and copies it verbatim into the flight artifact under
`dayof_addendum`. The prose form is THIS section, committed in the
pre-registration before the Batch opens. They are the same
addendum: the JSON is what the guards read, this section is what
the RECORD reads, and the numbers in them must agree. Where they
disagree, this section governs, per §9. ONLY THE JSON IS
MACHINE-ENFORCED, said plainly because this document is careful
elsewhere about which conditions a guard can actually hold (round
12): the runner refuses to submit without the dated file, its three
keys, and a committed-and-clean tree, and it requires the string
"## Amendment 1" in the committed pre-registration, but nothing
checks that THIS section was filled. A flight whose prose addendum
still reads "empty until flight day" passes every gate. That half is
executor discipline, and it is the half the RECORD is written from.

## 14. Revision notes

- v7.4 (2026-08-18): Amendment 2, pre-data. The analyze-side refreeze
  (θ_D 0.00253 → 0.00306, θ_W 0.00366, the s²(C) floor at 0.00067 and
  inert on D-sign), the re-entrant pedestal registered as a named
  systematic worth about ten points of P(detect), the depth-grid
  negative restated as the weaker claim the pedestal column actually
  supports, the corner block's exceptional couplings LOCATED in exact
  arithmetic (Q\* = √2 for C, 4.647043 for C′, none for the uniform
  control) which RETRACTS the earlier "no exceptional couplings"
  reading, the T1-CLEAN scale face superseding §9's "site scatter"
  wording with the MEAN it actually computes, and nine open items
  carried on the page rather than smoothed away (§2.8). One round of
  three empty reviews; the surviving findings are the open items.
- v1 (2026-08-16): first draft after the two-agent store sweep and the
  depth scout.
- v2 (2026-08-16): round 1 (physics recompute, spec audit, measurement
  statistics), ~35 findings folded; the two transversal classes, the
  Trotter EP, sort-free invariants, the prep 3× correction, the U-arm
  rebuilt.
- v3 (2026-08-16): round 2 (physics-math, spec-consistency), ~30
  findings folded; the ω = −4 null refuted and replaced, the secular
  factor re-signed, the split statistic pinned, D split into sign/mag,
  the N0 dressing moved to the full 3×3 generator, genre organs added.
- v4 (2026-08-16): round 3 (physics-math, spec-consistency), ~35
  findings folded, the load-bearing ones re-verified by the writer at
  the sources: **s² replaces s as the verdict statistic** (signed,
  real, resampleable; negative = conjugate pair, a reading not a
  failure); the envelope-immunity claim corrected (absolute s² biased
  up to −20% by the non-exponential post-selection envelope; ratios
  and r̄-differences immune; D-mag banded by G1 only); the (1,2) null
  restated as the C-minus-U generator difference (B̂ cancels) and
  split into a VOID-capable wiring check vs a report-only physics
  reading; the C′ arm given its own corridor and EP (Q ≈ 21 at
  J·dt = 0.15) and the W centre corridor-corrected (0.982·√3 ratio →
  2.9 on s²); the M-ladder re-anchored at the source verbatim (28.3%
  IS the M = 256 value; v3's "was M = 1" withdrawn); the candidate
  point declared priced-as-failing with the live alternatives named;
  arm-typed verdict rules replacing the mis-imported ternary grading,
  with equivalence tests for M and the null, floors under the W ratio,
  a partitioned ¬D-sign (FALSIFIED / Anti-D / underpowered), and named
  D∧¬M outcomes; N0-CLEAN conditions; the day-of re-gate made
  anti-circular (pre-submission inputs only) and its addendum moved
  before the Batch opens; DD/resilience-off, shared-skeleton, and
  deterministic arm order pinned; per-device twin band sets; billing
  re-anchored on delay-bearing flights.
- v4.1 (2026-08-16): round 4 (convergence check; every prior-round
  number re-verified to the digit, the seams caught): the secular
  factor stated per arm with C′'s OPPOSITE sign; the 2×2 Trotter model
  scoped to C (C′'s factor is ABOVE 1 at both live points; the table
  is the authority); the verdict evaluation made total with the order
  demoted to reporting precedence (FALSIFIED was unreachable);
  the ¬D-sign partition made disjoint and executable; D-mag's centre
  made a frozen function f_C(B̂); the arm-order rationale rewritten
  for s² (shift invariance kills the drift argument; W pair adjacent);
  s²'s sign statement made one-directional; the shape statistic moved
  to the centred third invariant; the s-vs-s² units pinned at the §4
  table and the J-term; the survival bridge (f_leak ∈ [8/15, ~0.9],
  0.2-2%) put on the page; the resilience pin softened to the house
  wording.
- v5 (2026-08-16): round 5 (final convergence check; physics verified a
  fifth time): the W arm's centres made frozen functions f_C(B̂),
  f_{C′}(B̂) (the last one-corner-arm number); s² pinned as
  Re[(e₁²−3e₂)/3] of the complex fitted R (decidable under the U-arm
  complex fit; Im parts = the detuning readout); FALSIFIED given its
  own physics-sized equivalence margin θ_F (containment in the
  detection band was scale-invariantly impossible); C′'s EPs at both
  live points on the page (Q ≈ 38 candidate, Q ≈ 21 vs flown 10 at
  fallback); the pair-odd C_p convention pinned and the swapped
  SUM/DIFFERENCE aside fixed (maximizing = difference, norm 1/3 = w;
  non-maximizing = sum, norm 1/√27); the non-maximizing triple written
  as three rates with (2,3) named as its decoupled dyad, and the (2,3)
  C′-minus-U difference null added as C′'s wiring check; the s²-floor
  VOID marked physics-or-instrument ambiguous with the RECORD naming
  both; H0 defined; the exact 1/6 route scoped to the room dyads; G4
  extended to s²(C′) and the ratio; survival 0.8% replaced by the
  f_leak bracket; s-vs-s² factor-2 notes at G4.
- v5.1 (2026-08-16): round 6 (propagation check; physics verified a
  sixth time, all repairs confirmed correct where present): the (2,3)
  C′ null propagated into §7 (its own equivalence test), G2 (separate
  margin), and §9 Class 2, trigger count 8; f_{C′}(B̂) and θ_F and the
  branch-reachability checks added to G1's freeze list; the W arm
  narrowed to the RATIO alone (the registered class arm) with both
  corner magnitudes moved to D-mag, resolving the collision with the
  registration rule; the Trotter detuning renamed δ_T with units (Δ is
  the anisotropy alone); site-pair {a,b} vs mode-dyad (i,j) notation
  pinned; the §12 mirror image corrected to {1,2,5}; the middle-rate
  ordinal, the number-operator-compression clause, and the M
  single-owner (G3) fixed.
- v5.2 (2026-08-16): round 7 (verdict-logic diff check, ZERO blockers):
  the floor-trip adjudication frozen (dose+N0 clean → ¬W class
  anomalous, else VOID; §3's partition restored); the G5 J pass band
  registered as VOID trigger nine; the anchor arm renamed M → A (the
  binding count keeps M, G3's object); the ¬D-sign partition made
  exhaustive (the above-θ_F case named); the named-verdict table
  declared governing over the reporting order; the null bullet made
  plural with separate margins; G1's equivalence margins scoped to the
  A arm and the projected SE added as its deliverable; the mirror-pair
  index clause; the ideal tag on the §1 ratio; the registration
  sentence extended to record the other conjuncts. Design stage
  CONVERGED: the next work is the committed gate (G1-G5, G7),
  Tom-gated.

- v6 (2026-08-16, gate-build stage): the gate's first G1/G3/G4/G5 runs
  folded (§8a): working point moved to Q = 10, J·dt = 0.15, 13 × 5
  grid, fractional-RZZ REQUIRED (now a Class-1 guard) and the
  priced-as-failing status REVERSED (the deep points are weight-
  devalued, not fatal); the W arm restated from the ratio to the level
  discrimination after failing its own reachability check (s²(C′)
  detection power 0.6 at budget); the s²(C′) floor retired; M = 1024
  confirmed from below (M = 256 FAILS the two-site retention
  criterion); budget re-tabled at 156 science PUBs ≈ 7 QPU min.
- v7 (2026-08-16, gate v2): the gate's own empty review round (four
  blockers: non-committed estimator, bias-free error bracket,
  winner's-curse freeze, G3 measuring the wrong M quantity) folded and
  the gate re-run; §6 AMENDED to pin the eigenchannel damped-cosine
  fit as the committed estimator (the complex one-step fit withdrawn
  for measured cause) and the negative-s² machinery withdrawn with it;
  the freeze moved to the 21 × 3 grid on the verdict metric
  (P(detect) 0.985 vs 0.84, ~11 min vs 6.8); §8a rewritten from the
  v2 run; p2 ≤ 0.6% and both-twins fractional checks added to
  Class 1; the T1-CLEAN scale ≤ 0.25γ̄ proposed from G4 v2; σ_J
  re-measured at the estimator level (8%/21% at 1%/2%); the v6
  consistency round's fix list (labels, γ̄ = 2×γ₀, grid remnant,
  exponent-7 remnants, W(ii)/D-mag duplication, twin rule) applied.
- v7.1 (2026-08-16): Tom's budget decision recorded (21 × 3 at ~11
  QPU min); the promotion round ran and BLOCKED promotion: the §7 C′
  numbers were from the withdrawn estimator (sign reverses under the
  committed one: s²(C′) reads below the hop-noise U floor, while the
  W exceedance s²(C) − s²(C′) survives at power ~4.7), and the
  committed estimator's dressing is large (s²-level ~1.5×/0.3× of
  ideal for C/C′) and belongs to f_C/f_{C′}; the writer's diagnostic
  resolved the dressing's mechanism (a J·dt-dependent mode-sign gauge,
  now PINNED to the continuous-H overlap in the gate, plus a
  prep-dependent secular population leak, real and deterministic);
  gate v2.1 adds the gauge pin, stage-2 C′/W lines, a held-out H0
  split, and the p2 row at the frozen grid; §7's W bullet is restated
  from the v2.1 run; em dashes swept, M→A completed, textual fixes
  from the promotion round applied.

- v7.2 (2026-08-16): Tom's SECOND budget decision on the honest
  post-pin numbers (option b: 21 × 3 at 16384 shots, ~22 QPU min);
  promotion round 2 blocked again on the same shape one level deeper
  (the gauge pin carried through G1 only), so G3/G5/dressing were
  re-run under the pinned gate and re-quoted (M-realization 15.3% at
  M = 1024; σ_J 10%/25%; dressing 1.21×/0.54×; the J response found
  EVEN, +4.5% both directions, a one-sided inflation bias entering
  the G5 pass band's rationale); the frozen-configuration numbers
  moved to the committed `--frozen` mode (outputs archived under
  `simulations/results/corner_beat/`); the pre-pin leftovers swept
  (status, §5 16/binding, §7 θ_D and W lines, §10 exponent 3.0,
  standard-gate ≤ 1.7, G4 provenance row relabeled); the hop-split
  bracket added to the outstanding list; §8/§8a retitled to the
  committed state.

- v7.3 (2026-08-17, G7 runner stage): the runner built and taken
  through FOUR empty review rounds (physics / spec-compliance /
  measurement-statistics per round, plus one interpreter round);
  round 4 returned no physics and no statistics blocker. Folded into
  this document: Amendment 1 (the analysis chain and PUB plan as
  flown: pooled inversion + single clip, kept-count floor, registered
  R_MAX_FIT and R_BOOT, nested tables, N0 swept zeros, the 24-PUB
  station plan with per-station and per-job CALs pooled into the
  inversion, the predictive billing brake, the demodulation-basis
  record); the §5 flown-circuit realization named as the repo's
  degree diagonal (F152/D10/K-partnership; the boundary rz = the
  F2/F2b boundary clock) with the sweep recorded in "What the repo
  already holds"; the §6.3 safety rationale corrected from measured
  numbers (the forced channels INFLATE, and safety is the measured
  zero-spurious result, not a shrink); the §10 gate-count mechanism
  relabel (2 per block is the sector reduction; standard CZ reaches
  ~2 per block at optimization 2+, so the fractional advantage is
  per-gate error, not count) with the flown 970 pinned; the §10
  billing bullet updated (per-binding overhead bounded, not zero: the
  read-only re-query returned 119 s = 0.316 ms/shot for the
  concentrator job, refuting the same day's first 69 s note; the
  duration caveat and the brake named);
  §8 G7 marked DONE with its evidence; the §8a outstanding list
  restructured around the θ_D refreeze package (analyze-side chain,
  fit-health margins from the ideal excursion + H0, G3's cumulative
  retention face, G3(b)/G5 at the frozen grid). Reviews recorded:
  round 1 found six guard-layer blockers (submission without frozen
  constants, vacuous realized-profile assertion, unreachable commit
  gate, best-effort DD-off, endpoint-only T1 sampling, silent
  get_counts pooling); round 2 found the -O strippability, the
  delay-station confound, the nested-tables requirement and the
  invert-then-clip bias; round 3 demonstrated the δω-rail
  false-VOID and the missing Amendment; round 4 the constants
  manifest, the CAL-per-job need, the brake's None-usage no-op and
  the grid-incomplete PASS. All folded; each finding verified from
  below before application.

- v7.3 doc round 2 (2026-08-17, two fresh agents: record audit with
  recompute + cold coherence, plus the same morning's billing
  re-query and the §11a store sweep). Corrections folded: the
  billing note's 69 s refuted by API re-query (119 s = 0.316
  ms/shot; §10 rewritten, "overhead zero" withdrawn); the
  concentrator's bound-circuit count corrected 9216 → 6144 (only
  the 24 sink PUBs were bound); the F140 fence had pointed at the
  WRONG block ((0,1); F140's own block is the (1,1)-type corner,
  this room's, so the fence now states the real distinction);
  §11a element 4 re-attributed (the U-arm scalar is this doc's own
  exact 1/6 route, NOT 45e2f29's C_l law, which covers only the
  protected rooms); fractional-RZZ provenance moved to this runner
  (record-parity records the flag was never set in-house);
  registrations added that the runner already implemented but the
  doc never froze (D-sign on the point estimate, 95% bootstrap
  percentile CIs over bindings, per-arm-pair VOID scoping, θ_W's
  producer = the refreeze package, dressed A/Π₃ centres, the
  T1-CLEAN µs→γ̄ bridge, θ_D 3σ̂-of-H0 wording); §9 completed
  (p2 Class-1 item, grid-incomplete as the tenth Class-2 trigger,
  DD/twirling and queue-warn wording matched to the code); the §8a
  COMPLETENESS paragraph now names the manifest complement; smaller
  repairs: canonical-γ₀ ratio note (Q = 1.5 vs Q = 10, no bare-number
  2×), Lemma-5 index-direction caution (twice), 1/6 holds for 13 of
  15 pairs, 0.90625 = 1 − 6/64 attributed, dose rule at depth 0,
  CAL count 20, header cost 23.7, Tom-note that Amendment 1.5 moved
  the projection ~22 → 23.7 inside the accepted band, billing-cap
  ownership (Tom, not G1), F152/F2 shared sentence, staircase
  anchors recomputed-here label, coined-label markers, entry-24
  characterization (3-magnon cats). Runner side: certify JSON
  persists the sector/prep/trajectory deviations; W now VOIDs by
  branch on incomplete grids (was an annotation).

- v7.3 doc round 3 (2026-08-17, two fresh agents with new lenses:
  mathematician with full recompute + external referee,
  FIXES-REQUIRED). The three referee blockers, folded: the s²(C)
  floor adjudication REROUTED (the old routing sent a floor trip to
  a D-sign-conditioned label, so the larger the falsification the
  more certainly it read as instrument failure; now: floor + clean
  dose + clean N0 = DATA, falls to the ¬D-sign partition); the p2
  Class-1 guard named as ISOLATED-gate while the circuit runs
  simultaneous layers (day-of layer-fidelity gate added to §9, owed
  before freeze); flown-but-kept-count-floored cells now COUNTED as
  grid-incomplete, persisted as floored_cells, VOIDing D-sign/W
  (they had been silently dropped from the fit). The mathematician's
  construction find: the DEPTH-0 OFF-BY-ONE (realized dose n−1 vs
  nominal fit axis n·dt; ~16.6% of s² given away, the booked
  −1.4%/−3.4% "dressing" dominantly this; the refreeze re-fits on
  t_eff = max(n−1,0)·dt, registered in §5). Also folded: the honest
  §1 scope (the registered pair tests split + ordering, not the
  width numbers; ¬D-mag named as the width law's falsification
  face; W-scope question 4 to Tom); θ_F REGISTERED as 0.75× the
  dressed prediction with the CI arithmetic on the page (round 4
  re-based it: 60% of the 0.00538 base); θ_D's
  H0-mean and hop-fraction BRACKET joined the refreeze; the C′
  2×2-model sentence corrected (fails at every tested Q); the
  retired point's C′ EP 38 → 42.5 exact (model/exact provenance
  split); the f_{C′}-vs-measured 3σ reconciliation item; B̂'s 3×3
  estimator marked G2-pending in §6 (two triggers evaluable only
  then); G4 re-aimed at the background PROJECTION onto C's
  mirror-odd part (symmetric backgrounds move d exactly zero,
  measured) with the T1-scale bound honestly near-non-binding as
  scatter; G5 gains the Δ_eff/static-ZZ signed line and the
  dyad-frequency agreement band as a D-mag precondition; leakage
  to |2⟩ named as the largest unpriced systematic (§11); τ_step and
  deep-end wall duration pre-registered as estimates (§10); §12's
  mixed-space claim corrected (dim-2 rooms at every N ≥ 4, N = 6
  unique in dim 3); exponent/survival refreshed at 970 gates and
  both p2 corners; smaller repairs (Lemma-5 shift form N-general,
  F140/F153 phrase attribution + count fence, entry-23 named as
  nearest kin + registry-scope note on "new", m̄ = 5/3 re-classed
  as derived, stale edge/spurious numbers refreshed at the frozen
  point, s²(U) two-config note, job/CAL-count arithmetic, f_leak
  12/15 count, M-ladder scaling note, power-vs-P(detect)
  disambiguation, "every step within 0.02" provenance label).
  Runner side: floored_cells tracking + VOID, fit_residual_bound
  as a manifest key, saturation-base pin note, T2-proxy comment,
  N0 zero tables in the npz, seeds + payload caps in the payload
  config.

- v7.3 doc round 4 (2026-08-17, two fresh agents: coherence
  executor walking the verdict tree on four scenarios + record
  auditor with symbolic recompute; both FIXES-REQUIRED, most
  findings sitting in the round-3 repairs, the known shape). Folded:
  the floored-cell DOUBLE ROUTE closed (one route: grid-incomplete,
  Amendment 1.2's definition extended; §7's fit-health sentence no
  longer claims the cell); the VOID SCOPE registered for ALL TEN
  triggers via the consumption map (global / grid-incomplete /
  per-consuming-arm; the (2,3) null leaves D-sign standing; the G5
  J-band conflict RESOLVED: voids the magnitude/anchor arms, D-sign
  and W stay evaluated with a named flag); T1-CLEAN unified as
  THREE registered faces (projection, scale, station-to-station
  time structure, the telegraph criterion that had no test);
  "verdict-grade" removed from the header and Open question 1
  (§8a's own 0/50 caveat governs); the §11 bullet that round 3's
  leakage insertion had DECAPITATED restored; Q ≡ J/γ̄ defined
  (round-4 repair: used ~40 times, never defined); the dead-C′
  fork registered in D-mag (the freeze decides which regime, the
  RECORD prints raw s²(C′) either way); θ_D as ONE formula
  (E[d|H0] + 3σ̂); θ_F's BASE pinned (the dressed predicted
  difference, 0.00538 today, CI half-width 60% of base); the
  round-4 factor audit named the §2 "norm" a SPREAD (2× the
  operator norm; the certificate's C_l carries the 2); **the flown
  room's class split UPGRADED to derived-in-closed-form** (symbolic
  overlap route: entries 1/12 and √3/36, eigenvalues {±1/6, 0} and
  {±√3/18, 0}; §12's derivation follow-up discharged for this room,
  the general locus stays gated); far-end d = 0.00554 ± 0.00133 put
  on the §8a page; the header's 0.7 labelled an interpolation
  estimate; CAL pooling aligned with the code (head + per-job pairs
  pooled, station pairs stay out as drift diagnostics); R_BOOT's
  owner unified (freezes with θ_F); the wall-duration coherence
  exponent corrected to 13-19% (was "of order 30%"); guard-bound
  survival corrected to 0.53/4.49%; unrounded ±0.00165 quoted so
  the division reproduces the power; the ¬W ∧ D-mag branch added to
  G1's reachability list; §6.x reference convention stated; smaller
  repairs (per-rate error ≤ 1.2% with 1.116% exact, table-cell
  0.273/0.274 note, sideways-arc scope precision, F2-vs-F2b
  same-shape sharpening, Q_REGIME_ANCHORS as the γ₀-pin owner,
  Open-questions retitle, 0.13-min CAL arithmetic). Runner side:
  stale 4-PUB station comment corrected to 6.

- v7.3 doc round 5 (2026-08-17, one fresh cold executor-referee
  walking the verdict tree on six scenarios against the runner;
  FIXES-REQUIRED, ten majors, most sitting in the rounds-3/4 repair
  layers plus four design-layer finds). Folded: the VOID
  three-valuedness got its ROUTING RULE (a void arm makes every
  line requiring it unreachable; partial patterns print
  "<verdict>; <arm> VOID (<trigger>)", never CONFIRMED or FALSIFIED,
  and never the unlisted-pattern INCONCLUSIVE); T1-CLEAN re-scoped
  GLOBAL (the lighting predicate was keyed to injected profiles
  while decay on any chain site corrupts every arm; N0 has no lit
  sites at all, so the old rule could never reach D-mag); the
  J-band scope now stated ONCE (§7's list; §4/§8-G5 point at it);
  the dead-C′ fork got its honesty rider (today's numbers admit
  only the revert-to-REPORTED regime); grid-incomplete D-mag/A
  defined as REPORTED excluded-cell fits (no per-cell s² exists);
  the DEEP-END KEPT-COUNT MARGIN joined the refreeze (deepest
  cells ~87-736 kept of 16384 against the 50 floor, 1.7× at the
  worst corner, twelve cells with veto power: freeze the expected
  profile, require ≥ 3× margin); the runner gained a MACHINERY
  gate (a hard-abort list of unbuilt verdict consumers, deleted as
  gate work lands: frozen numbers alone no longer release
  submission) and a TIME_AXIS sync constant (a manifest key:
  freezing "time_axis" to realized_dose forces the code change in
  runner and gate, closing the axis-unbound gap);
  fit_residual_bound wired into fit-health (it existed as a key the
  code never read); the bootstrap NaN-exclusion rule and seed
  registered (silent exclusion narrows the CI toward FALSIFIED);
  θ_F's G1 line aligned with §7's ownership; M priced from the
  15.3% s²-realization spread, never the parent's 4.4%; dose
  certificates' dual class (pre-submit belt + in-data suspenders)
  stated; §6.1 aligned with Amendment 1.1; §6.3's weighting
  sentence now describes the per-trace normalization the certify
  parity pins; smaller repairs (PUBs-vs-pairs units, queue-warning
  wording, unrounded d_W = +0.00634 ± 0.00134, δ_T coefficient
  2.005-2.011 exact, the 2.7-vs-2.9 rep-spread note, the
  standard-gate two-config note, 60%-of-base in the round-3 note,
  Q forward pointer at first use, one orphaned capital).

- v7.3 doc round 6 (2026-08-17, one fresh cold executor-referee,
  seven scenarios; FIXES-REQUIRED: one real design blocker four
  rounds had missed, plus routing edge cases in the rounds-4/5
  layer). Folded: **BRAKE ⇒ FLIGHT VOID registered** (jobs fly
  depth-ascending, a firing brake truncates the deep end, and any
  missing cell voids D-sign/W outright: the reachable "full budget
  spent, no verdict" outcome now has its registered rule, a parked
  flight, with the interleaving alternative considered and declined
  on payload-chunking grounds, Amendment 1.6); the three unlabeled
  guard states got their labels (global trigger ⇒ flight-level
  "VOID (<trigger>)"; void partition-precondition ⇒ "NOT DETECTED;
  partition unavailable (<arm> VOID)"; grid-incomplete ⇒ "VOID
  (grid incomplete)" beside the reported fits); SUB-ARM RESOLUTION
  registered (A and D-mag are two-line arms; the ¬D-sign
  partition's "A holds" means A's C line, so a (2,3)-null trip
  leaves FALSIFIED reachable on a true null; dead-C′ freeze makes
  D-mag mean D-mag(C)); the C′-hole-voids-D-sign reconciliation
  stated (grid-incompleteness breaks the no-drop guarantee θ_D was
  frozen under, not the arm's data); the DEPTH-0 CI UNDER-COVERAGE
  registered (single unswept PUB, no binding-resampling variance,
  highest-weight point: narrows the CI toward FALSIFIED; refreeze
  adds a shot-noise resample or widens the CI by the measured term,
  machinery-ledger item; the "inflates θ_F/κ" this note carried until
  round 15 was anti-conservative, see §7); the ledger grew to 32 keys (T1-CLEAN's three faces split,
  dressed A/Π₃ CENTRES keyed beside their bands); τ_step joined
  the day-of addendum's enforced schema; the echo face named a
  SURVIVAL FRACTION (one delay is not a fit); the standard-gate
  "unflyable" sentence scoped to the retired pricing; header bold
  balance repaired; the routing rule un-folded from the J-band
  bullet into its own paragraph; §2's second "norm" sentence says
  SPREAD; stale ordinal and pointer fixes; seed_boot into the
  analysis JSON.

- v7.3 doc round 7 (2026-08-17, one fresh cold executor, eight
  scenarios incl. the brake firing; FIXES-REQUIRED: three majors,
  no blocker; convergence continues). Folded: the ¬D-sign
  partition's "guards clean" SCOPED (round 7: it means the triggers
  whose §7 scope touches D-sign's consumers {C, U} or A's C line,
  never the whole bank, so a (2,3)-null trip leaves FALSIFIED
  reachable, closing the falsification-vs-no-verdict fork the
  round-6 sub-arm rule had opened one level up); the s²(C) floor
  DEFINED at last (a frozen absolute lower bound on measured s²(C)
  from its own H0 distribution; trips below; its ratio-arm ancestry
  died in v6 and the level-statistic function is now stated: W is
  uninformative when C itself did not split); the runner no longer
  OVERWRITES the brake status with COUNTS RETRIEVED (the governing
  artifact now reads BRAKE-TRUNCATED with the brake string
  preserved: the fact distinguishing "parked by rule" from "the
  device dropped cells"); PARTIAL arms defined (one void line +
  one non-void: never satisfies a conjunction, prints
  "<arm> PARTIAL (<line> VOID)"; the ¬D-sign A-rule is the one
  registered exception); §5's "all pairs pool" aligned with
  Amendment 1.5 (head + per-job pairs pool, station pairs are
  drift diagnostics); §4's stray "verdict-grade" removed and the
  science cost stated as 22.5; the day-of addendum section now
  names its ENFORCED schema incl. tau_step_us; the twin rule got
  its operational note (a source edit + re-commit, not a day-of
  selection); the "difference-nulls moved" orphan and the
  pairs-vs-PUBs unit slip fixed; runner DEFAULT_CONSTANTS brought
  to parity with the SHIPPED constants file (the manifest keys
  without entries have none yet BY DESIGN and MISSING blocks
  submission; the complement's size is a moving number the analyze
  printout reports, not something to quote twice, round-12 repair
  after this note's 15 and §14's 17 stood 120 lines apart); docstring version pin
  v7.2 → v7.3.

- v7.3 doc round 8 (2026-08-17, one fresh cold executor, nine
  scenarios; FIXES-REQUIRED: five majors, no blocker). Folded: §3's
  U-arm identity CORRECTED by its own measurement (s²(U) ≈ 0.0017
  hop-noise floor vs s²(N0) ≈ 4e-5 under one noise model, a ~40×
  gap: the U-side baseline/band come from the analyze-chain H0,
  never from centring on s²(B̂); the identity survives only as the
  compressed-model statement); routing state (iv) added (a trigger
  voiding D-SIGN ITSELF prints the flight-level "VOID (<trigger>)"
  with surviving sub-arms REPORTED: with D-sign void no named line
  exists and the catch-all stays closed); the s²(C) floor trip now
  FORCES ¬W regardless of the computed exceedance (else a race
  between two independently frozen numbers); the WALL-DURATION
  model unified (the runner's serial gate-count estimate, ~145 µs,
  was ~4× the layered physical time and pushed the long T1 station
  point past Heron T1; the runner now computes the layered model
  with a 1.5× buffer, t_deep ≈ 33-55 µs, stations derive from
  that, the transpiled schedule stays the authority); the in-job
  dose-certificate face and the floor-forces-¬W evaluator joined
  the machinery ledger (15 entries); the per-pair CAL drift record
  is now COMPUTED and persisted (Amendment 1.5's promise; it was
  recoverable but never produced); the brake-truncated verdict
  label unified ("VOID (grid incomplete; brake-truncated per
  Amendment 1.6)"); f_leak defined at first §4 use; smaller
  (DEFAULT_CONSTANTS wording, stale trigger-nine and entry
  ordinals, 2×2-model gap 0.5 not ~1).

- v7.3 doc round 9 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: two majors, both collisions between earlier
  round-repairs, five minors). Folded: the all-non-void rule
  SCOPED to the CONFIRMED-type lines with the ¬D-sign partition's
  A's-C-line carve-out named as the one registered exception (the
  two paragraphs had contradicted each other fourteen lines apart
  on the (2,3)-null-under-true-null walk); the s²(C) floor trip
  now forces ¬D-SIGN too (the race round 8 closed for W was still
  open on D-sign: s²(C) < floor and d > θ_D are simultaneously
  satisfiable when s²(U) comes in anomalously low, and a detection
  whose C-level is H0-indistinguishable is no detection); θ_W got
  its BRACKET rule (maximum over the f_leak ends, W's worst end
  the opposite one from D-sign's); the two wall-duration estimates
  named as one model at two buffers (bare 28-42 µs, buffered 33-55
  µs, exponent 13-19% extending to ~25%); θ_F's base assumption
  named (noiseless dressed difference; ~5/6 against the H1 chain's
  own d, refreeze re-states); the bootstrap NaN-replicate count
  now PERSISTED in the analysis JSON (it printed but did not
  survive into the record); depth-0 PUBs described correctly in
  §10's budget line; the T2-proxy named in §9's line rule; the
  survival pair corrected to the 970-gate exponent (4.85 → 7.5%);
  "THE STATES" header count un-staled.

- v7.3 doc round 10 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: three majors, one family: OWED QUANTITIES WITH
  NO LEDGER HOOK). Folded: layer_fidelity_bound and
  leakage_ket2_price joined the manifest (34 entries; §11's
  "largest unpriced systematic" and §9's day-of threshold could
  never have BLOCKED before); the floor-evaluator machinery entry
  renamed to carry BOTH halves (¬W and ¬D-sign, rounds 8-9; an
  executor emptying the ledger by names would have shipped the W
  half only); precondition (ii)'s "untripped" → "has not VOIDED"
  (the one-word fork: the floor's DATA-path trip must not bar the
  ¬D-sign partition it routes to); the PARTIAL rule completed (a
  PARTIAL arm whose surviving line FAILS counts as ¬arm for the
  ¬-conditioned lines, closing the last unlabeled pattern); the
  INCONCLUSIVE reachability clause updated for the round-9 floor
  rule; §5's three off-by-one magnitudes named as three
  measurements that do not close under one factor (the fixed-axis
  re-fit measures the truth); §3's two s²(U) configurations
  labeled; the 0.6% interpolation's round attribution neutralized;
  ddof = 1 on the page. The round also verified: every §8a digit
  against the committed records, the full room physics from below,
  runner-vs-doc on all counts and schemas, every verbatim citation,
  and the gate's RNG seeding (no PYTHONHASHSEED dependence).

- v7.3 doc round 11 (2026-08-17, one fresh cold executor on the ten
  scenarios; FIXES-REQUIRED: three blockers, eight majors, and no
  arithmetic error on the page, every number it recomputed
  reproducing). The three blockers were one shape at three heights.
  (B1) The machinery ledger had an entry for every owed NUMBER's
  evaluator but none for the VERDICT-LINE COMPOSER or the ¬D-sign
  PARTITION EVALUATOR, so emptying it by its names would have
  released a flight whose falsification arm has no code at all; both
  are entries now, and so is the parking rule's third condition,
  which the writer found unhooked while confirming it (§13 question
  2: its two siblings hard-abort, "the gate reproduces §8a" was a
  human act). (B2) W carried fit health as a SUFFIX where D-sign
  carries it as a BRANCH, so a dead channel made d_W = NaN and
  `NaN > θ_W` printed **FAIL**: a falsification-shaped verdict
  manufactured by instrument failure, the thing §3 forbids and the
  thing round 3 already repaired once for the floor. W now takes
  D-sign's branch structure, and the writer closed the mirror hole
  the reviewer did not report: a non-finite d had read as NOT
  DETECTED. (B3) The s²(C) floor, defined in round 7 and given two
  forcing powers in rounds 8 and 9, had no PRODUCER while both its
  sibling thresholds have one; since the detection rule is now a
  conjunction, that was a free hand on the detection criterion
  itself. Registered as E[s²(C)|H0] + 3σ̂, bracket-maximum, with the
  reachability check in G1's freeze list and one clarification the
  reviewer had backwards: a floor trip under a TRUE NULL is the
  route the null is supposed to take, not a pathology. Majors
  folded: trigger PRECEDENCE registered (the scope map said which
  arms each trigger voids, nothing said which trigger NAMES the line
  when several fire; nine ranks then, ten after round 12 split the
  two difference-nulls, causal order, appended never
  replaced, and the floor and the J band provably cannot lead);
  `time_axis` made a real switch (it was sync-checked but never
  read, so freezing it to "realized_dose" would have recorded the
  §5 off-by-one repair as landed while every fit still ran on n·dt);
  a mid-batch job failure now saves its counts (the paid artifact
  had been written without them and was unreadable by --analyze,
  with the partial dumps unreadable too); the day-of addendum
  declared as ONE object in TWO forms, the JSON file the runner
  demands finally named on the page; `r_saturation_frac` given a
  BASE field, since its own note said the frozen entry states the
  base and the schema had nowhere to put it; the analyze PENDING
  report driven off the manifest instead of ten hardcoded keys,
  which immediately showed that half the manifest is MISSING from
  the constants file outright and most of the rest unfrozen (the
  printout carries the live counts; they are not quoted here,
  because a count in prose goes stale the next time the manifest
  grows, which is exactly how two of them ended up disagreeing); per-key
  bootstrap NaN counts (the count was D-sign's alone while
  nanpercentile drops replicates from every key). §6.3's
  certify-parity claim was resting on single-trace evidence that
  could not touch the joint two-trace weighting it describes; the
  joint case now runs and is also exactly 0.0. Minors: the degree-
  diagonal citation corrected a third time in this arc, at the
  source (the sentence is F152's alone, F2 is its spectrum in the
  commuting case, D10 step 3 is the same Laplacian on the
  Liouvillian block, a different operator); the shared-skeleton
  invariant described as what it is, a construction identity with a
  regression guard, not a tested one; the wall-duration pair no
  longer called "one model at two buffers" (the edge terms differ
  too, ≲1.5%); the four code-pinned flight-determining constants
  named; and the analysis wall clock stated at last, measured.

- v7.3 doc round 12 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: one blocker, six majors, fifteen minors, and again
  no arithmetic error on the page). THREE OF THE FINDINGS WERE IN
  ROUND 11's OWN REPAIRS, which is this project's oldest pattern and
  worth naming rather than quietly fixing. (a) The wall-duration
  sentence: round 11 adopted a reviewer's "≲1.5%" instead of
  computing it, and the true bare gap between the two estimates was
  13-21%. Computing it found the cause: §10 charged ~100 ns per step
  for the injection rz layer, and RZ is a VIRTUAL frame change on
  Heron that spends no time, so the estimate was buying 6 µs of
  imaginary circuit. τ_step is 0.36-0.60 µs and the bare deep end
  22-36 µs; with the allowance gone the document's model and the
  runner's are one model at two buffers, which is what the round-9
  phrase had claimed all along. The runner's preparation booking was
  wrong in the other direction and is corrected with it (five Givens
  on OVERLAPPING pairs are ten rzz sub-layers, not five). (b) The
  analysis wall clock: the 48 s anchor came from a four-cell
  artifact that never enters the fit path, and the per-replicate
  cost barely depends on the cell count, so the ×63 extrapolation
  was meaningless. The number is now owed from the full-M analyze
  rather than guessed. (c) The degree-diagonal citation, corrected
  for the third time in this arc and now, finally, at the level that
  matters: every registered holder owns the LIOUVILLIAN (0,1)
  block, a fourth holder was sitting in this flight's own parent
  record `CONCENTRATOR_GEOMETRY.md`, and what the flight actually
  flies is a subordinate remark inside F152, not a registered
  result. The BLOCKER was the round-10 family's third live instance:
  the kept-count floor of 50 is the sole gate on the
  grid-incompleteness trigger, hence on a VOID of D-sign and W, and
  it sat in the "no verdict consumes them" list of code constants
  with no manifest key and no ledger hook. Registering it surfaced
  the substantive half: §8a requires ≥ 3× margin over that floor at
  the worst modeled corner and today's own survival numbers give
  1.74× at the p2 guard bound, so the design fails its own
  precondition there and the freeze must move the floor, the budget
  or the guard bound, pre-data. Majors folded: the constants note
  claiming certify-parity forces the committed gate onto the frozen
  time axis was false (certify passes ONE axis array into both
  implementations, so the gate's own axis is never compared; the
  gate side is ledger-only, the same shape round 11 hooked for the
  parking rule); `dmag_band_Cp` given a "REPORTED" SENTINEL, because
  the regime this document tells the executor to EXPECT was not
  representable in the constants schema and the only way to reach
  submission was to invent a number for a band that will not exist;
  the bootstrap NaN-replicate bound registered as a manifest entry
  and an evaluator, having been a sentence with no key; the two
  admissible rate-saturation BASES brought onto the page, where G2
  can read them, instead of living in code where a third plausible
  spelling would hard-abort the flight; and the s²(C) floor's
  reachability check corrected for CHAIN, since its producer is an
  analyze-chain quantity carrying the hop floor while the 0.00539 it
  was written against is the noiseless dressed reference, which
  leaves roughly 1.3× on a threshold that forces ¬D-sign. Minors:
  the two difference-nulls split into separate precedence ranks (one
  layout permutation fires both); the "precedence line below"
  pointer corrected, since it pointed at the arm order and
  contradicted its own example; the floor's scope line brought up
  from its round-8 state; rate saturation added to §7's fit-health
  definition, where only the code had it; the axis line printed
  unconditionally; §9's order given the two same-day runs that
  hard-abort without it; the prose addendum named as the half no
  guard enforces; and two prose counts of one complement, 120 lines
  apart, replaced by the live printout.

- v7.3 doc round 13 (2026-08-17, one fresh cold executor, pointed
  first at round 12's own repairs; FIXES-REQUIRED: no blocker, eight
  majors, and no arithmetic error on the page for the third round
  running). FIVE OF THE EIGHT MAJORS SAT IN ROUND 12's REPAIRS, which
  makes the pattern measured rather than suspected: across rounds 11,
  12 and 13 the original material has produced no arithmetic error at
  all, while each round's own fixes produced the next round's
  findings. The heaviest: (a) the degree-diagonal attribution, wrong
  a FOURTH time and this time wrong in the claim-bearing direction.
  Round 12 wrote that the flown object is "owned nowhere". It is
  owned twice, precisely: D10 step 3 writes the single-excitation
  Hamiltonian with its degree diagonal as a numbered step, and
  PROOF_K_PARTNERSHIP holds V_eff(ℓ) = #bonds − 2·deg(ℓ) together
  with the open-chain boundary discontinuity that IS the boundary rz
  layer we fly. All three wrong versions failed the same way, by
  characterizing ownership in the abstract instead of quoting the
  sources; the fourth version quotes them. What is new here is only
  the per-Trotter-layer application and the hardware realization,
  which is all §11a element 1 ever claimed. (b) The "REPORTED"
  sentinel was DOCUMENTED and not implemented, so for one round the
  freeze gate accepted any non-null value on any manifest key: a
  string left in θ_D would have passed the last guard before a paid
  submission and raised on the paid data. Now enumerated, restricted
  to `dmag_band_Cp` and `f_Cp_dressing`, and every other key must
  hold a number. (c) The floor's reachability example did its
  arithmetic against the very comparand the preceding sentence had
  disqualified, understating the margin by ~20% in the
  conservative-LOOKING direction; the H1-chain comparand is ≈ 0.0065
  and the margin 1.55-1.86×. (d) The wall-duration pair still was not
  one model: the bare end was quoted from the step sub-layers alone
  while the buffered end used the full 370, and the T2 guard sentence
  compared against the bare number although the executor schedules
  the buffered circuit. Both ends now stand, the guard's comparand is
  the buffered one, and the two directions of "safe" are separated,
  since a longer estimate is conservative for decoherence and RISKIER
  for a station delay that must not sit past T1. (e) §9's order gained
  the two day-of runs in round 12 without their commit obligation, so
  an executor following it verbatim met a hard abort with the Batch
  window open. Also folded: E[s²(C)|H0] and its σ̂ actually added to
  §8a's package, having been declared owed there and not added; the
  manifest length retired as a prose count, the third such; the
  kept-count floor read from the constants in analyze rather than
  from the module, so a re-freeze cannot be honoured by the gate and
  ignored by every analysis; replicate-level floor trips counted, the
  same silent-exclusion mechanism as the NaN replicates one level
  down; the two INCONCLUSIVE causes split into two labels; the
  billing anchor, station placement and floor added to the flight
  artifact, which §10 had claimed already; and the gate side of the
  axis freeze declared as having no machine check at all.

- v7.3 doc round 14 (2026-08-17, one fresh cold executor pointed at
  round 13's repairs; FIXES-REQUIRED: no blocker, seven majors, and
  for the fourth round running no arithmetic error in the original
  material. ALL SEVEN MAJORS sat in round-12 or round-13 repair
  layers). The finding that matters most is about method, not about
  the flight. **The degree-diagonal attribution was wrong a fourth
  time, and the reason all four versions failed is that the repo had
  already written the map I kept rebuilding**: F152's own credit
  paragraph, inside the document being cited, reads "D10 Step 3
  derives the Laplacian on the open chain and F2 carries the
  ZZ-supplies-the-degree sentence, both at uniform γ", and lists
  VacuumBlockReductionClaim and PROOF_R90 Lemma 5 besides. Round 13
  had gone further and DENIED F2 any holding, on the strength of a
  phrase search: F2's Note says "supplies the diagonal shift that
  turns the adjacency matrix into that Laplacian" rather than "degree
  diagonal", so the phrase is unique to F152 and the statement is
  not. Searching for the phrase produced a false denial about the
  statement. Round 13 also mis-sorted D10, whose Step 3 is titled "On
  the (0,1) block the generator is the graph Laplacian" and holds
  both faces at once. The bullet now cites the existing inventory,
  adds only the two holders it does not list (PROOF_K_PARTNERSHIP for
  the flown boundary, CONCENTRATOR_GEOMETRY in parallel words), and
  says out loud that rebuilding an inventory the repo already keeps
  is what produced four wrong versions. Other majors, all in prior
  repairs: the wall-duration paragraph had been corrected in one half
  and left at the 365-sub-layer numbers in the other, thirty lines
  below, in the very sentence that names 370 (now one count, both
  ends, and the derived gap restated as 12-20%); `_STRING_KEYS`
  accepted ANY non-empty string, so a British "realised_dose" or a
  hyphenated spelling would have passed the last guard and raised
  inside the estimator on paid data, which is round 13's own stated
  failure mode one key over (now a closed enumeration); three
  manifest keys had no reader and no ledger entry, `dmag_band_C`,
  `dmag_band_Cp` and `leakage_ket2_price`, the last being what §11
  calls the largest unpriced systematic, so freezing it would have
  released submission silently; §9's order still lacked the DRY RUN,
  while Amendment 1.5 told the executor to carry its realized job
  plan into an addendum that had no key for it, and the dry run in
  fact printed only the job count (the per-job PUB and bound-circuit
  numbers are printed and persisted now, `job_plan` joins the enforced
  schema, and the dry run joins §9's order); the floor's reachability
  requirement had no number at all, now registered as a RATE,
  P(s²(C) < floor | H1) ≤ 0.01 at both bracket ends, because a ratio
  bar guards counts against Poisson noise while what this floor risks
  is tripping on a true detection, whose size is a rate; and round
  13's INCONCLUSIVE split had left κ serving as both a precondition of
  the ¬D-sign partition and the discriminator of a label inside it,
  which is self-cancelling, so κ is now the partition's FIRST
  discriminator and the registered label set is written out closed.
  Minors: the backend-operational guard failed OPEN on an unreadable
  status and now fails closed on a paid submission; replicate-level
  floor trips recorded per arm, since d is a between-arm difference;
  "the guard's comparand" renamed, no such guard existing; the §4
  secular-error citation corrected to the per-rate figure; the 2×2
  model's agreement recomputed by the writer at 2.4-4.3% against a
  claimed 4-5%; the T1-CLEAN bridge given its missing ÷γ̄; Lemma 5's
  index direction marked as derived rather than stated, the same
  file's Lemma 3 running the other way; and Amendment 1.6's "one
  label" scoped as the grid-incomplete rank's name rather than an
  exemption from precedence.

- v7.3 doc round 15 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: no blocker, ten majors, fourteen minors). **The one
  finding in the ORIGINAL material after four rounds of finding none
  is a SIGN, and it points the wrong way**: the depth-0 CI
  under-coverage bullet (registered round 6) offered as its fallback
  "measures the term and inflates θ_F/κ by it". The CI is biased
  NARROW, FALSIFIED is containment inside [−θ_F, θ_F], and the power
  condition is SE ≤ κ × projected, so inflating either loosens
  exactly the test the bias already favours: an executor
  implementing the registered correction would have raised the
  false-falsification rate while believing the move conservative. The
  conservative direction is to widen the CI, equivalently to shrink
  θ_F and tighten κ; both copies repaired. Everything else again sat
  in prior repairs. The degree-diagonal bullet, fifth version: round
  14 wrote that the phrase "degree diagonal" is unique to F152 while
  quoting that same phrase from CONCENTRATOR_GEOMETRY eighteen lines
  earlier; the phrase is in seven tracked files, this one included,
  and is unique to F152 only INSIDE the formula registry. The lesson
  survives the correction and is sharper for it: a string search
  cannot stand in for a claim about content, which is how round 13
  came to deny F2 a holding it has. Also at the sources:
  PROOF_K_PARTNERSHIP's boundary sentence breaks **L1**, the
  chiral-conjugation lemma, not K directly; and the registry's
  F2-vs-F2b pair shares only the N-vs-N+1 denominator, since F2b IS a
  single-excitation Hamiltonian spectrum by its own title while F2 is
  the Liouvillian block, so the contrast is Liouvillian-vs-Hamiltonian
  AND Neumann-vs-Dirichlet. The wall-duration gap was still computed
  across two sub-layer counts, now 14-21% like for like at 370 with
  its denominator named; the "closed" label set had two strings for
  one state and an eleventh label loose, now one entry per state with
  the catch-all named as the composer's obligation rather than a
  label; `family_error_rates` had neither reader nor ledger entry,
  the same silent-release shape as `leakage_ket2_price` one round
  earlier; the per-twin band rule was unimplementable in a
  scalar-valued schema and is now carried the way a twin switch
  already is, as the committed constants file; a sentence placing
  the centre functions under W survived v6's retirement of the ratio
  arm; and the runner still carried a comment describing §10's
  superseded rz allowance. Minors: the 2×2 model's agreement
  recomputed a second time, because round 14 used the rounded δ_T
  while the same paragraph pins the exact Floquet values (2.42 /
  4.41 / 3.61 / 4.13, so TWO rows leave the declared spread, not
  one); the gate's printed `bound=0.3333` distinguished from our
  reading 1/3; the concentrator record's 0.316 named as our
  recompute in the sweep bullet as §10 already named it; the frozen
  gate run's f_leak = 0.53 distinguished from the registered 8/15;
  the informational 0.95 saturation flag marked as unregistered;
  the job plan persisted in the flight payload and not only the dry
  one; §9's Class-1 list told that the operational check fails
  closed; forward pointers for s² and κ at their first use; and
  "CONFIRMED-type" defined where it is load-bearing.

- v7.3 doc round 16 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: no blocker, five majors, five minors, seven nits).
  **Two of the majors are factor-2 errors in PROSE DEFINITIONS added
  by rounds 14 and 15**, both caught by recomputing against §2's own
  numbers, and both worth stating because the formulas they gloss
  were correct all along. (a) The forward pointers called s² "the
  mean square of the three pairwise rate differences". That is the
  sum over three; §6.4 and the code divide by six. The honest plain
  name is the SAMPLE VARIANCE (ddof = 1) of the three fitted rates,
  which reproduces §2's s = 2/3 on {8/3, 10/3, 4}γ̄ while the prose
  version gives twice s². (b) Round 15 wrote that m̄ is the mean of
  the compressed rates "in units of γ̄"; §1's own mapping is
  rate = 2·s·γ̄, so the unit is 2γ̄ and m̄ = 5/3 IS the (10/3)γ̄ that
  §2 states. A factor 2 sat on the only bridge between the gate's
  printed bound and the flown rates, which is the same shape the
  round-4 audit caught in §2's norm-versus-spread. The remaining
  three: §5 still said PROOF_K_PARTNERSHIP's boundary term "breaks
  K" 400 lines after round 15 corrected exactly that in §2, and §5
  is the section read while building the circuit; the registered
  label set omitted the PARTIAL-CONJUNCTION lead form that the
  routing rule's own worked example prints, so the document's
  example was an unregistered string; and the reconciliation
  sentence for "instrument-suspect INCONCLUSIVE" read A as a whole
  arm although round 9 narrowed the partition's precondition to A's
  C LINE, so two registered strings collided on the ordinary cell
  where A's C line holds and its C′ line fails, which is exactly
  where a true null puts you. That cell now belongs to the
  partition, with the C′ line printing as A PARTIAL beside it.
  DOC-side minors: §5's DD/twirling clause aligned to fail-closed;
  the concentrator's 6144 bound circuits and 24/12 PUB split labelled
  as our reconstruction rather than as recorded; the "~1.4%" prep
  booking given its denominator; and the 2×2 model's two >4% rows
  released from a "2-4% method spread" they were never commensurable
  with. RUNNER-side: the stale rz-allowance comment stopped restating
  a retired figure at all, and `kappa_projected_SE` named inside the
  ledger entry that owns it, the last manifest key without a reader.
  (Round 15's own note mixed doc and runner items without saying
  which was which; this note separates them, as rounds 2-6 did.)

- v7.3 doc round 17 (2026-08-17, one fresh cold executor;
  FIXES-REQUIRED: no blocker, five majors, nine minors, five nits;
  no arithmetic error in the original material for the fifth round
  running). **All five majors sat in round 16's label-set repair and
  all five landed on ONE cell**, the true-null cell where A's C line
  holds and its C′ line fails. That concentration is the finding:
  rounds 13 to 16 each patched the verdict-line rules and each patch
  produced the next round's defects, which is what a design being
  patched instead of stated looks like. §7 therefore now carries
  **the composition algorithm**, five steps, normative, generating
  every line mechanically, with the registered label list demoted to
  its OUTPUT rather than a second source of truth. The five defects
  it dissolves: "A PARTIAL (C′ line fails)" misused PARTIAL, which
  §7 defines in terms of VOIDNESS, so a failing line made an arm ¬A
  and never PARTIAL; precondition (iii) fails in two ways, by FAILURE
  and by VOID, and round 16 equated them although a void line is
  neither true nor false and cannot satisfy ¬A at all; the bullet
  round 16 rescoped was never amended at its own site seventy lines
  above, this document's own named failure mode; the worked example
  applied the sub-arm rule to A but not to D-mag, which the same
  trigger makes PARTIAL too, and dropped the trigger the registered
  form requires; and the lead clause of the partial-conjunction form
  was undefined, with no named line covering D-sign ∧ D-mag ∧ A when
  W is void, so the example printed a string the closed set did not
  contain. Under the algorithm that state simply enumerates its arms
  and claims no name, which is what "no partial pattern is ever
  promoted" always meant. DOC-side minors: the 8-versus-10 aux PUBs
  in the concentrator reconstruction reconciled (the 10 is the 8 aux
  plus 2 CAL); the preparation-BOOKING distinguished from the
  five-sub-layer CORRECTION, a factor 2 that round 16's denominator
  fix left in the subject; "the one place the doc quotes the root"
  scoped to measured quantities, §1 and §2 quoting the ideal s; the
  U-versus-N0 gap corrected from ~50× to ~40× (42.5); and the
  certificate's size coordinate s marked as a different quantity from
  the split statistic s, with w's forward pointer. RUNNER-side: the
  ledger comment for `kappa_projected_SE` moved below the entry that
  actually owns it, where it had read as if the composer discharged
  it.

- v7.3 doc round 18 (2026-08-17, one fresh cold executor pointed at
  round 17's structural rewrite; FIXES-REQUIRED: **TWO BLOCKERS**,
  seven majors, six minors, three nits; no arithmetic error in the
  original material for the sixth round running). **Both blockers
  were in the two structural repairs of rounds 14-17, which is the
  finding**: replacing a patch with a structure is the right move and
  it does not exempt the structure from being wrong.
  **BLOCKER 1, the admissibility guard.** Round 13 taught it "numbers
  must be numbers" after round 12 documented a sentinel it never
  implemented. But nine manifest keys are registered by this document
  as objects that CANNOT be scalars: §7 says of the dressed centres,
  in as many words, "a frozen NUMBER cannot exist pre-flight"; the
  day-of width scalings are functions; the family error rates are a
  PAIR; N0-CLEAN is THREE thresholds; the band-validity window has
  dimensions; the deep-end kept profile is one number per deep cell.
  The guard would have rejected the correctly frozen objects at the
  last gate before a paid submission, leaving no route forward except
  inventing scalars, which is post-data freezing under another name.
  Each key now declares its FORM and the guard checks THAT form. THE
  FORMS, registered here rather than left in code, the way round 12
  brought the saturation bases onto this page: **scalar** (a real
  number, and the DEFAULT for every key not named below);
  **positive_integer** (R_BOOT); **vector** (a non-empty list of
  reals: deep_end_kept_profile); **map** (named reals:
  family_error_rates, N0_clean_thresholds); **map_of_intervals**
  (named axes each with endpoints: band_validity_window, the one key
  the document calls dimensioned); and **code** (f_C_dressing,
  f_Cp_dressing, A_centres_dressed, pi3_centre_dressed,
  dayof_width_scaling), whose frozen object is the committed function
  itself, recorded as source, symbol and sha256, with the guard
  hashing the file and comparing. Every band and margin is a
  HALF-WIDTH around its frozen centre and therefore a scalar,
  dmag_band_C/Cp, A_margin_C/Cp, pi3_report_band, J_pass_band and the
  three T1-CLEAN faces included; round 19 removed a "scalar or
  interval" form that had accepted both a half-width and a pair of
  absolute endpoints on a conjunct of CONFIRMED, which is not a
  guard, and it removed an interval-only J_pass_band that would have
  hard-aborted a perfectly good scalar freeze.
  **BLOCKER 2, the algorithm against a surviving patch.** Round 10's
  clause promoted a partial pattern to the named line "Split
  observed, anchor failed"; step 4 forbids exactly that. On a
  reachable state, a C′-arm dose-certificate failure with a failing
  A C line, the algorithm printed an unnamed enumeration and round
  10's clause printed a registered name. The algorithm governs, and
  the pattern is not left unlabeled: under step 4 its enumeration IS
  its label. Majors folded: the PARTIAL print form now names the
  SURVIVING line's truth value, without which a partial A with a
  holding C line and one with a failing C line printed identically,
  the whole content of the anchor arm; step 5 gained a line FORMAT,
  having had a label and none, so six of the ten forks had a verdict
  an executor could name and not write; step 1 now reads each arm's
  LINE SET from the frozen constants, since under the reverted-band
  freeze this document tells the executor to expect, D-mag has one
  line and cannot be PARTIAL at all, which the worked example got
  wrong in exactly that regime; the append-the-voids rule scoped to
  step 3's flight-level lines instead of "the whole table", two live
  orderings for one line; step 3 now appends every firing trigger and
  not only the leading one; `A_margin` split into one key per LINE,
  the same ruling round 6 made for T1-CLEAN's three faces; the arms
  counted as FOUR everywhere, the nulls being VOID triggers with
  equivalence-test shape and not arms; and a forcing registered as
  never overwriting a VOID, since a void line has no value to force.

- v7.3 doc round 19 (2026-08-17; FIXES-REQUIRED: **three blockers,
  nine majors**, all three blockers inside round 18's two structural
  repairs). Eight consecutive rounds have now found their defects in
  the previous round's fixes while finding none in the original
  material, which is a measurement and is recorded as one.
  **BLOCKER 1 and 2 were one act:** round 18 added `A_margin_Cp` to
  the sentinel-eligible keys, in code, and explained it in a JSON
  note, while §7 registers TWO sentinel keys "and no others". The
  note also asserted something this document nowhere decides, that
  the D-mag C′ BAND reverting to REPORTED makes the A arm's C′ LINE
  reported too; those are different statistics and A's C′
  discriminability is assessed nowhere. A pre-data registration
  cannot be made in a code comment, so the key came back out, and
  with it the state the algorithm could not class at all: a REPORTED
  line is neither true, false nor void, so an A pinned at two lines
  with one of them reported fell through step 2 entirely.
  **BLOCKER 3:** round 18's own step-5 example printed "D-mag
  PARTIAL" in the regime where D-mag has no C′ line, seven lines
  above the paragraph in which the same round diagnosed exactly that
  error in the neighbouring example. Majors: the void-arm paragraph
  was cut back to a pointer, having carried its own ordering rule
  through four rewrites, each repairing the pointer and leaving the
  rule; Amendment 1.6 and step 3 had the same VOID line with the
  appended triggers on opposite sides of the parenthesis, while step
  3 cited that Amendment as its reason for appending; step 5's
  format now says the partition LABEL carries D-sign, so enumeration
  starts at W; the print vocabulary is FIVE forms stated once, the
  ¬arm-with-lines form having appeared in a registered example and
  nowhere else, and the PARTIAL form having been registered in three
  places of which round 18 updated one; precondition (ii)'s
  subsumption argument replaced, since the G5 J pass band falsifies
  its premise in one lookup; the value-form table moved onto the page
  with its default stated; `scalar_or_interval` deleted, having
  accepted a half-width and a pair of absolute endpoints as the same
  frozen object on a conjunct of CONFIRMED; `J_pass_band`'s
  interval-only form deleted, which would have hard-aborted a
  perfectly good scalar freeze and made a TENTH non-scalar key one
  round after nine were found; `band_validity_window` given a form
  that can express a dimensioned window; and the `code` form made a
  real check, hashing the committed file and comparing, where it had
  accepted `{"source":"todo","symbol":"todo","sha256":"todo"}` at the
  last gate before a paid submission.

- v7.3, the guard-bound decision (2026-08-17, Tom): the deep-end
  margin fork round 12 opened is CLOSED on the isolated axis by
  **Amendment 1.8**, the p2 Class-1 guard tightened 0.6% → 0.5%, the
  third of the three levers round 14 priced. Folded through the
  document rather than at the guard alone, because the bound is an
  INPUT to a chain: §4 and §9's Class-1 line, Open question 2's
  parking rule, §10's exponent (4.85 exactly) and survival
  (1.27-7.53%), §8a's kept-count margin (87 → 208 kept, 1.74× →
  4.17×) and its now-discharged owed measurement AT the bound, §7's
  replicate-level floor note, and the header, whose interpolated
  P(detect) ≈ 0.7 is retired because the guard now sits ON the
  measured ladder's first rung. TWO THINGS THE FOLD FOUND that the
  decision did not contain. The ≥ 3× requirement is exactly
  equivalent to an effective p2 ≤ 0.538%, which hands the owed
  `layer_fidelity_bound` a computed candidate where it had none. And
  the tightening does NOT reach the layered error: at §9's 1.3-2×
  bracket the same arithmetic gives 1.12× and then a floored deep
  end, so the layer-fidelity gate, not this bound, is what holds the
  deep end, and §9 now says so where it used to imply the opposite.
  A FOLLOW-UP PASS the same evening asked whether "three levers" was
  the whole list. It was not, and §8a now carries the inventory priced.
  Two corrections came out of it, both in text written that day. The
  "more shots" lever was priced UNIFORMLY and declared dead at ~40 QPU
  min; targeted at the 36 cells where the requirement actually binds it
  costs +1.16 min, inside the cap, so it was alive and Amendment 1.8
  says so now. And §8 claimed G1 jointly chose the depth grid ENDPOINT,
  which the committed gate never varied: 60 steps in all sixteen
  configurations, spacing alone searched. Corrected at the source, with
  the endpoint moved to §12 as an unsearched axis. The pass also built
  the weighted line budget as a SCORE and as a projected margin, and
  established by building it that a budget GUARD is redundant under the
  0.5% ceiling; and it recorded that the ≥ 3× and the floor 50 are both
  set numbers with no measurement behind them, together with the two
  levers (f_leak's 0.9, the 3 itself) that must never be pulled to make
  the design pass.

## Open questions for Tom (spend and registration decisions)

1. ~~The budget fork~~ **DECIDED TWICE, finally (Tom, 2026-08-16
   ~15:53, on the honest post-pin numbers): option (b), the 21 × 3
   grid at 16384 shots, then ~22 QPU min.** The frozen-config
   confirmation (200 reps, both f_leak ends, held-out H0; verdict-
   grade AFTER the ≥ 500-rep refreeze, §8a) is IN §8a: worst-end
   P(detect) 0.920 at H0 false 0.000.
   ~~NOTE for the flight go~~ **RE-CONFIRMED (Tom, 2026-08-17
   ~13:20): the 23.7 min stands.** Amendment 1.5 replaced the aux
   plan and moved the projection from the ~22 he priced to 23.7,
   inside his accepted 21-25 band and under the 25-min hard abort;
   he re-confirmed it on those terms. The day-of go re-checks the
   projection against the realized job plan, nothing more.
2. ~~The C′ arm stays~~ **CONFIRMED (2026-08-17, Tom's delegation:
   "das musst Du entscheiden").** The C′ arm stays, its cost is
   inside both numbers above, and the parking rule stands in SHAPE:
   park if fractional-RZZ is unavailable on the flown backend, if
   day-of p2 > 0.5%/2q on the used edges (**0.6% until Amendment 1.8,
   Tom, 2026-08-17**: the one leg of this rule that has moved, and it
   moved TIGHTER, so the rule parks more often and never less), or if
   the committed gate fails to reproduce §8a. All three are pre-data
   and day-of checkable, which is why the rule's shape stands. What
   confirming it surfaced (round 11): only the first two are
   MACHINE-enforced (fractional exposure on backend and twin, the p2
   guard on the used edges, both hard aborts); the third was a human
   act, since --certify checks parity against the committed gate and
   never that the gate reproduces §8a's governing numbers. It is a
   machinery-ledger entry now, owned by the G1 committed-table mode,
   so the rule's third leg blocks submission until its check exists.
3. Committing the gate + this document (the promotion rename and the
   first commit of the arc) happens after the next empty round on the
   re-run gate; no QPU spend before a separate go from you.
4. ~~W's registration scope~~ **DECIDED: option (a), W stays
   REGISTERED with the honest §1 scope** (2026-08-17, Tom's
   delegation: "das musst Du entscheiden"). The round-3 referee's
   point stands and stays on the page: W is an ORDERING test
   maximized by a dead C′, so the registered pair carries no
   quantitative class-law content. It is still the only registered
   line that tests the ORDER of the two transversal classes, where
   D-sign alone tests only that a split exists, and it carries
   measured power (4.7 at D-sign's worst f_leak end, 3.8 at its
   own). A limit that is stated is not a defect; the alternative is
   worse in the one way a pre-registration can be worse, since an
   unregistered arm that reads well afterwards becomes a post-hoc
   promotion, and the cost is identical either way. The §1 scope
   sentence is the load-bearing half of this decision and must
   travel with any quotation of the W result.
