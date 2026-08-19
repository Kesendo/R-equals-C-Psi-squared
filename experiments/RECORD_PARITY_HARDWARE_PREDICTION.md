# Pre-registration: the record parity trichotomy on hardware, the SIGNED angle law (v33, sim-gate stage, round-31 four-lens fold on the re-opened loop; not flown)

**Status:** v33, SIM-GATE STAGE (2026-08-04): the flight was called by Tom
("IBM ist freigeschaltet, lass uns das Uhren Experiment fahren"; advance go
08:12, quoted in full at the HARDWARE RECORD when it exists) and the loop
re-opened exactly where the bank note says it must, at the sim-gate stage.
The sim-gate machinery is built and REHEARSED at the pinned scales: ONE
sim-gate module (circuits, estimator, verdict, band freezing:
`simulations/record_parity_gates.py`, PROMOTED at the round-31 fold, prior-work
M4: the rehearsal numbers below are its outputs and the pinned empty round on the
frozen numbers needs a repo artifact to run against; F129's committed gate scripts
are the precedent; it enters the repo with the next commit, BEFORE the binding
freeze) plus the Aer-parity/η_nom reader in the external pipeline
home (its VALUES and provenance enter the frozen-constants table; the
external home is named under Runner home below). Rehearsal record (local,
regenerated deterministically by the module's freeze-full mode; scales qualified
round 31, stats m1: the run predates the round-30/31 control-axis extensions, so
its 1.2×10⁶ is 12 banks of 10⁵ against a pinned axis set that now demands more; the
UCB 3.0×10⁻⁵ is PER BANK at n = 10⁵; the binding freeze runs the full pinned set):
joint power
0.9834 at the worst-admitted corner, negative controls 0 CONFIRMED across those
12 banks. Round 28 folded TWO
pinned-rule amendments, both found from below by the gate itself (the guard
band coverage envelope and the η_nom profile arithmetic; see the Revision
notes and the amended sites). Rehearsal numbers are NOT frozen: the binding
freeze runs under these amended rules from a fresh calibration snapshot with
the representative line re-pinned, and appends the frozen-constants table
here. Design rounds 28–31 of the re-opened loop are FOLDED (v30–v33, the
Revision notes hold the ledger; round 31 ran FOUR lenses; the prior-work/
repo-consistency lens is a standing member from round 31 on); the loop stays open
until an empty round; round 32 runs on this state. The v29 status below stands
as history.
**Status (v29):** DESIGN CAMPAIGN BANKED (Tom's call, 2026-07-18: the loop that truth
cannot end, ended by decision from outside it; recorded as a decision, not as a clean
round). Rounds 1–27 folded; every reviewer finding recomputed from below before
folding; the physics is closed: exact from-below recomputes CLEAN in rounds 7 through
27, twenty-first consecutive, round 27 verifying BEYOND the claim (the double-ratio
verdict exact under static write misses to 0.3 rad: the r-independent miss cancels).
Round 27 found the SIXTH sign-class instance: the clause-(b) fallback band,
verdict-bearing and unswept ONE round after the sign accounting was declared
"complete", so the completeness claim is RETIRED: sign-sweep membership is ENUMERATED
(b_blind, b_curve's SUP, T̂, the floors, LR, V_S, the (b)-fallback band, the negative
controls), never again declared complete, and the enumeration gained a SECOND sign
AXIS in round 31: beside the readout-bias sign, the coherent watcher-angle sign
(ε_w-negative grid members; |ρ̂(1.75)| is first-order in it, everything else checked
second-order or symmetric). Round-27 stats verdict: "essentially CLEAN …
nothing blocks the freeze" (joint power ≥ 0.9946 reproduced from below;
false-CONFIRMED 0/2×10⁵ at both signs). Banked, not flown; if the flight ever
proceeds, the loop re-opens at the sim-gate stage and stays behind Tom's explicit GO.
The tested face is the
SIGNED law: magnitude = Proposition 1's record radius dressed by the watcher factor η,
sign = the signed-coherence
corollary (the r = 2 flip). TO-FREEZE constants get
numbers only under the pinned rules below: fluctuation widths from the sim gate's SIGNAL
model, deterministic centers from the systematics ledger, the NULL model as negative
control only. Flight only
after: design rounds converge → sim/null/Aer gates recorded → pre-reg committed → fresh
calibration + Class-1 pre-flight aborts → Tom's explicit go → ONE Batch job → Class-2
in-job VOID gates → verdict (terminal rule at the gate order).
**Backend (pinned from Tom's fresh pulls, 2026-07-18 08:35; the quoted calibration
figures are that pull's HISTORICAL values; the binding freeze re-pins the line from a
fresh snapshot and the Class-1 aborts consume only day-of values, round 30):** primary
**ibm_kingston**
(2q median 1.90e-3, readout 0.745%, T1/T2 238/136 µs; maintenance 7/21 15:00–19:00), twin
**ibm_fez**; `rzz` is a basis gate on both pinned Herons (native fractional; VERIFY at runner
stage that θ = π/2 is INCLUSIVE in the fractional range: already read green at the
2026-08-04 rehearsal transpile, re-verified at runner stage). Twin-submit, fly whichever queue
empties, cancel the other (normally unbilled; if both complete, both bill, tie-break
pinned, rounds 11/13/15, and SCOPED here (round-27): the tie-break governs ONLY the
both-complete case, so first-to-complete and earliest-submit are never rival selectors;
the job with the EARLIEST SUBMIT-RECORD TIMESTAMP governs,
fixed at submit, fully data-independent; the SUBMIT ORDER itself pinned round 31,
spec m8, the tie-break is data-independent only if the order is: the PRIMARY,
ibm_kingston, is always submitted first, before any day-of information is
consulted; the other is discarded unanalyzed, recorded). Twin-submit RECONCILED
with F129's posture (round 31, prior-work: [IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md)
§9 pins "one device, one job, one lab; no cross-device claims"; that is a SCOPE
sentence about claims, and this design keeps it: ONE analyzed job, the flown
device's own band set governs, no cross-device claim is ever assembled; the twin is
a queue hedge, not a second lab, and STAIRCASE_NULLTEST already cites this file as
the pinned-twin precedent, which this sentence makes explicit rather than
accidental). **Band rule: each device gets its own band set
(Kingston + Fez), frozen at the sim gate and committed before flight; the flown device's
set governs.**
**Runner home (when built):** external pipeline `...\ibm_quantum_tomography\run_record_parity.py`,
modeled on `run_f129_ramsey_fringe.py` + `run_price_pair.py`.
**Law under test:** the F135 record-radius angle law: Proposition 1's distinguishability
face D_j(t*) = β_j = |cos(r·π/2)| (in the proof D_j = β_j·|sin θ_w|; the flown write
θ_w = π/2 makes trace distance and record radius coincide; this is the single-watcher,
γ = 0 face of the general radius;
[PROOF_RECORD_PARITY_LAW](../docs/proofs/PROOF_RECORD_PARITY_LAW.md)).
**The independent sibling law, cited (round 31, prior-work: this document had
re-derived its sign discriminator without naming the repo's second gated law for
exactly this geometry):** F136, the record LETTER law
([PROOF_RECORD_LETTER_LAW](../docs/proofs/PROOF_RECORD_LETTER_LAW.md), 87/87
pre-registered gates; registry entry F136 in
[ANALYTICAL_FORMULAS](../docs/ANALYTICAL_FORMULAS.md)): the flown line k–j–S–j′ is
F136's pointer branch (dressers D = ∅, private watchers P = {k}, controls Q = {j′}),
for which F136 independently predicts a pointer record in the ZY CHANNEL with sign
Π(−1)^{r/2}: the r = 2 rotation-by-π, "F135's signed-coherence corollary made
observable" in F136's own words, i.e. the entire flown discriminator as a second
gated theorem; the direct correlator ⟨Z_S Y_j⟩ that F136 states is what this
design's conditional-Bloch estimator reconstructs. GAUGE-BASELINE divergence named:
F136 pins signs relative to the r = 1 baseline, F135's corollary relative to r = 0;
this design uses r = 0 (correct per F135, and the r = 0 arm is the flown
normalizer). One F136 consequence is used directly at the η bundle (its clause
(iv), the S-idle note there).
The flown carrier ρ̂ is SIGNED (round-4): its magnitude is the trace-distance face, and its
sign carries additional law content, pinned to the proof's **signed-coherence corollary**
(PROOF_RECORD_PARITY_LAW, the Signed-coherence corollary inside the Law A section): the
watcher multiplier cos(r·π/2) is
negative for r ∈ (1, 2], i.e. the record returns at the even parity ROTATED BY π in the
equatorial plane (ideal multiplier ρ(2) = −1, not +1; the flown fixpoint is −η(2)); the
sign is gauge-free referenced to the r = 0 record. The flip is a prediction no
magnitude-only impostor WITH A SUPERPOSED WATCHER can fake (scoped round 29, physics
M3: a Z-POLARISED, never-watching k, H on k missing, k reset, or fully Z-relaxed,
reproduces the ENTIRE signed curve exactly, flip included, with η ≡ 1: j's coherence
ROTATES by rπ/2 instead of attenuating, and the r = 0-fitted axis projection cannot
tell the two apart; the polarisation leaks ONLY into |⟨Z_k⟩| (m ≡ ⟨Z_k⟩, the watcher's
Z-polarisation) and the transverse
channel (T̂_j = −m·sin(rπ/2)); the |⟨Z_k⟩| guard is the superposition
certificate FOR THE m-POLARISED BRANCH of this class ONLY (scope narrowed round 30,
physics BLOCKER: at m = 0, the maximally MIXED never-watching k, the leak vanishes
identically, |⟨Z_k⟩| = 0 and T̂_j = 0, and every flown science/guard PUB is Z-diagonal
on k, so it sees only diag(ρ_k), where |+⟩⟨+| and I/2 COINCIDE: the m = 0 member is
invisible to the entire Z-diagonal PUB set and reproduces the law's curve exactly),
its band capped by the pinned z_k DISCRIMINATION CEILING (round 30, stats BLOCKER,
the eighth protection-interaction instance: the coverage rule makes a TO-FREEZE guard
band as LARGE as possible, which here WIDENS the impostor's escape window m < band;
the ceiling is the opposite-direction pin, guard_z_k ≤ 0.05, an envelope above it is
unfreezable; the boundary control below quantifies the escape at m = the frozen band);
the m ≈ 0 discrimination is carried by the WATCHER-TOMOGRAPH PUB pair (its own
paragraph in the guard section: ⟨Y_k Z_j⟩, +1 at r = 1 for a superposed watcher,
identically 0 for ANY Z-diagonal k; DIAGNOSTIC-only, it gates the mechanism sentence,
never the record verdict, so the flown VERDICT tests the record-radius law on j and
the watching MECHANISM is claimed only under the tomograph certificate); a
decoding-table row names the class). NOT a
mutual-information measurement. The MI
reading "the record survives watching" (Law A's I = 1 − h₂((1+β)/2)) is a corollary that
additionally needs the classical-quantum premise V_S = 0; that premise is GUARDED in-flight
(the S-coherence guard below), and the MI sentence is claimed only under that guard (a
derived reading, never a second measured claim). Law C is not tested and, since v11,
not invoked anywhere in this design. Higher parities than r = 1, 2
are extrapolations of the fitted curve, not flown. Hardware cannot falsify the theorem; a
deviation reads "the device does not realize the angle law."
**The seeing over it (labeled, not the claim):** Tom's fixpoint frame (2026-07-18): stability
comes from pinning to agreed reference points; the fixpoint endpoints cos ∈ {0, ±1}
are the Grenzen (Niven's rational set is {0, ±½, ±1}; the flown endpoints are its
extremes, labeled precisely round 31, physics NIT); watching at the agreed (even) angle PRESERVES the watched record; the odd angle
erases; between them "das Leben je nach Blickwinkel."

## The claim in one line

The signed double-ratio record ρ̂(r) follows cos(r·π/2): it floors at the first odd parity
(r = 1, the watcher at the odd angle erases) and recovers WITH FLIPPED SIGN at the
watcher-dressed magnitude η(2) at the first even parity (r = 2, the watcher at the agreed
angle forgives, returning the record rotated by π; ideal η → 1 = full magnitude: the
TESTED discriminators are the flip plus substantial recovery, clause (b), since "full"
is untestable under a known device attenuation), with the generic signed curve between. The mechanism, grounded
from below (round-10): at r = 2 the total watcher gate RZZ(π) is a deterministic Pauli
(−i·Z⊗Z): k ends PURE and unentangled, the even watcher learns NOTHING, which is why
the record survives; at r = 1 the watcher becomes maximally CORRELATED with the witness
j, classically: ⟨Y_k Z_j⟩ = 1, concurrence C(k, j) = 0, a perfect 1-bit which-path
record, while k's entanglement ENTROPY against the rest is 1 bit (rounds 25–26), and
j's record of S is gone. Watching costs what the watcher learns (equal at the
parity endpoints; complementary, not identical, between them).

Not tested: Law C, J>0 transport, the plaquette, R_δ (the record-redundancy count),
parities beyond r = 2. MI sits in a third, labeled category: neither tested nor
unconditionally claimed; the V_S guard's 4 PUBs buy only the conditionally-appended
report sentence, never a verdict.

## Geometry and circuit

Four qubits on a hardware line: **k–j–S–j′** (the line satisfies the proof's Law A
hypotheses: deg(S) = 2 with witnesses j and j′, triangle-free, k not adjacent to S, j′
not adjacent to j). Per arm r:

1. H on all four.
2. RZZ(π/2) on (S, j): the write. 3. RZZ(π/2) on (S, j′): the control write.
4. Watcher on (j, k): **two half-angle gates RZZ(r·π/4)²** for every r > 0 (native range;
   r-independent watcher depth in GATE COUNT; native fractional-rzz pulse DURATION may
   scale with angle; the runner-stage VERIFY list covers this, and it is what makes η worst
   at r = 2; if the device pads fractional angles to fixed duration, η is ~r-flat and
   freezing b_forgive at r = 2 stays conservative). Per-gate coherent budget doubles:
   ε_watch = 2·ε_cal. The
   r = 0 arm omits the watcher: the structurally privileged reference (2 gates vs 4), the
   reference arm of the double ratio (it contributes one numerator and one denominator
   factor).
5. Terminal measurement of ALL qubits. S in Z; witnesses pre-rotated: variant A = H (X);
   variant B = **Sdg then H** (Y; order load-bearing; sim gate asserts Y-recovery). k in Z.
   (The record lives EXACTLY on Y in the ideal circuit: zero X component; RZ frame
   offsets can tilt it, hence the axis fit: variant A reads ~zero record COMPONENT by
   design, but its X-data feeds every conditional Bloch 2-vector, the axis fit, and the
   transverse guard; no variant is discarded from ρ̂.)

RZZ(θ) = exp(−i(θ/2) Z⊗Z). No INSERTED delays in science circuits; DD OFF; runtime
resilience_level = 0 (V1 idiom: the pinned requirement is NO runtime mitigation, flag
names per the VERIFY list's API note; the CAL0/CAL1 readout inversion below is our own
offline mitigation, not runtime resilience).
**Gate-count budget (hard abort):** r > 0: exactly 4 two-qubit gates; r = 0: exactly 2. A
CZ-fallback transpile is **NO-FLIGHT** (round-7: every band is frozen for the native
fractional path only; a CZ flight would be judged against bands never derived for it;
there is no CZ verdict basis, so there is no CZ flight). (The budget also catches the
single-vs-double watcher compilation bug, which would fake the half-period impostor.)

## Observable and estimator

Mitigation BEFORE conditioning (CAL0/CAL1 tensor-product inversion incl. finite-shot noise;
quasi-probabilities kept, no clamping; all models, SIGNAL and NULL, condition
identically). Stated residual:
tensor-product mitigation cannot correct S↔witness readout correlations; the SIGNAL
model's readout-correlation variant WIDENS the committed bands, and because conditioning
is ON the S readout, an S↔witness correlation is a coherent BIAS on the conditional Bloch
vectors, not only spread (round-7; round-9 correction of the round-8 over-generalization:
the double ratio cancels only the SYMMETRIC, multiplicative part of an r-independent
correlation; the ASYMMETRIC additive part, a fixed S-outcome-dependent assignment offset
on the witness, does NOT cancel even when r-independent, because the offset is constant
while the true record varies across arms: from below, ρ̂(2) bias ≈ 2.75 SE per 1%
S-correlated assignment differential (per-% basis PINNED round-27: the earlier
"1–3.5 SE per 1–2%" headline was per-basis ambiguous and read low; the round-27 stats
MC lands at ~5.5 SE per 2% at the pinned shots; model estimates ranged 0.7–9 SE per 2%
across review
models; the two protections split cleanly, round-14: false-CONFIRMED protection rests
on the empirical negative-control UCB and is coefficient-FREE; the power-side
signal-model injection is the ONLY place a coefficient enters, and it uses the worst
reviewed one; the residual bias itself is read SIM-GATE-measured; the SINGLE frozen
injected coefficient is recorded in the ledger at freeze; the range is review history,
never a frozen input, round-20), and one
sign EASES clause (b)): **the injected QUANTITY pinned in one equation (round 29,
physics M1: "a percent of WHAT" was undefined and two quoted figures sat on bases 2×
apart: the then-current 0.065 b_blind envelope on the narrower, the 0.043 ρ̂(1.75)
shift on the wider: historical figures quoted for the RATIO argument only, round
31 spec n4; the operative envelope is the coverage-grid value): the ceiling, the injection, and the in-job crosstalk pair's abort compare are
ALL in the S-branch assignment differential Δ ≡ the shift of
[p(w=1 | S=0) − p(w=1 | S=1)], which equals the induced offset on Ŝ; "2%" means
Δ = 0.02, the WIDER (conservative) of the two readings; the crosstalk REDUCTION is a
FOUR-CIRCUIT combination, pinned (round 31, spec BLOCKER + physics minor, converging:
the pair alone cannot identify Δ; it flips S and the witness prep TOGETHER, so two
of the four (S, w-prep) cells are missing, and the missing cells come from CAL0 and
CAL1, which the reduction therefore consumes as inputs): per witness w, with
e01(S) = p(w reads 1 | w prepared 0, S) and e10(S) = p(w reads 0 | w prepared 1, S),
the pair's two circuits give e01(S=1) and e10(S=0), CAL0 gives e01(S=0), CAL1 gives
e10(S=1), and **Δ̂ = ½·[(e01(0) − e01(1)) + (e10(1) − e10(0))]** (the equal-weight
form matching the ~50/50 conditional branch populations of the science arms) is what
the ≤ 2% Class-2 compare consumes; historical per-% figures quoted in this document are review
MC on mixed bases and are superseded by the gate's own frozen values**, handled by
PINNED SIGNS, not center surgery (round-12 rebuild,
the rounds-8–10 band-center-shift mandate is DROPPED as redundant AND power-fatal: an
easing-sign MC returns 0/20000 false CONFIRMED for null and every impostor on UNSHIFTED
bands, because a constant additive bias cannot jointly fake the flip, the r = 1 minimum,
and the 7-arm shape; while a center shift comparable to the whole ~2.5–3σ (b) margin
(the p0.13 quantile margin, shot-independent, distinct from the arm-power figures in
Arms and shots) collapses valid-device power toward 0, and no shot count repairs a
center offset). The
pinned rule (round-24 GENERALIZATION, the third instance of the sign-class bug forced
the general law): the SIGNAL model's band runs inject the bias at the 2% ceiling at the
PER-STATISTIC WORST sign-combination. Concretely: each magnitude CEILING, b_blind,
b_curve's SUP, AND the pooled transverse guards T̂ (round-25, the FOURTH sign-class
instance: T̂ is Ŝ's perpendicular twin, S-conditioned, so the S-antisymmetric assignment
offset enters the variant-A/X channel at the same magnitude as the Y channel; both are
the same physical Z-readout after a pre-rotation; a single-sign-frozen T̂ band terminally
VOIDed admitted devices at the other signs with up to ~96% probability in review MC; the
Z-leak and ⟨Z⟩ sanity statistics are EXEMPT, being S-SYMMETRIC; the offset averages
out), freezes from the UPPER ENVELOPE over the sign axis {−2%, 0, +2%}, i.e.
the ε_w-center-ALIGNED worst case; the LR threshold joins the sign sweep round-26 (its
own section); the SIGN envelope and the SE edge are INDEPENDENT axes taken in OPPOSITE
directions: the envelope buys coverage, the SE edge buys conservatism; never conflate
them into one direction (round-26; magnitudes stated honestly, round 30, stats n16:
the two are not comparable in size: the sign envelope moves a band by tens of percent,
~45% on b_blind, while the SE edge moves it ~1%) (round-24 BLOCKER: the round-23 floor fix did not
reach b_blind; the additive offset shifts ρ̂(1) nearly 1:1 and the (b)-HARDENING sign is
the (a)-EASING sign, so the uniform hardening freeze gave b_blind ≈ 0.042 instead of
≈ 0.065 and routed valid aligned-sign devices to INCONCLUSIVE up to ~40%, end-to-end
joint power 0.57 at the ±2% endpoint: UNFREEZABLE as written, since the sign-swept gate
would read 0.57 and the ladder cannot move a mis-signed center; envelope fix: joint
0.996, false-CONFIRMED untouched at 0/2×10⁵); b_forgive freezes from the LOWER
ENVELOPE over the sign axis (the gloss resolved round 30, spec M5 + stats M6: the
sign wording had inverted a THIRD time, and the resolution retires the label rather
than rewording it again: for a LOWER cut like b_forgive the sign axis is a COVERAGE
axis, so the frozen threshold is the minimum over {−2%, 0, +2%} of the per-sign
p0.13 quantiles; the "CONFIRM-HARDENING sign" LABEL is RETIRED as a mislabel, the
envelope needs no sign name; the never-ease-CONFIRMED duty lives in the +SE edge,
which shifts the frozen value UP, and, the load-bearing point, clause (b)'s
false-CONFIRMED protection rests on the SIGN of ρ̂(2) (a non-law device cannot fake
the flip) and is CERTIFIED by the sign-swept negative controls' measured ≤ 0.1% UCB,
never by the threshold's sign choice; verified sign-robust: (b) power ≥ 0.995 at both
signs); the VOID FLOORS keep their
round-23 per-floor worst sign. Asymmetry stated (round-24): floors are ladder-EXEMPT and
terminal, hence direct worst-sign construction; b_forgive is ladder-remediable and
gate-watched. What is FROZEN is the pinned Δ-mechanism at the Δ = 0.02 ceiling: the
r-independent S-antisymmetric additive channel, whose sim-gate-measured coefficient
(~5.5 SE per 2%) is an OUTPUT of that model (restated round 31, stats m5: the old
sentence here, "all injections sit at the upper end of the reviewed coefficient
range", described the round-13 state; the reviewed 0.7–9-SE range is history, its
9-SE high end is NOT what is injected, and no domination of the other reviewed
models at equal Δ is claimed; the protection against a WORSE-than-modeled
correlation is the coefficient-FREE negative-control UCB, per the round-14 split)
(round-13 context: a benign model + a high true coefficient wastes the flight to
INCONCLUSIVE);
the NEGATIVE
controls SWEEP the injection sign {−2%, 0, +2%} and take the bound at the WORST
(highest-rate) sign (round-27, "swept, not assumed" extended from the round-26 LR rule;
the old single CONFIRM-EASING sign, maximal help to the non-law device, is the swept
axis's empirically expected worst point, review MC 94.0 vs 103.5 SE, impact nil; their
measured ≤ 0.1% UCB bound is THE model-independent false-CONFIRMED protection for this
mode); the in-job crosstalk pair below is the hardware bound (a Class-2 in-job VOID gate:
it protects the verdict, never moves a band; freeze-before-data holds); and the gate
REPORTS the residual bias on ρ̂(1) and
ρ̂(2) as a diagnostic. Post-selection novelty owned (both branches
kept; ~half shots per branch).

- Conditional Bloch components b_s(w) = (⟨X_w⟩_s, ⟨Y_w⟩_s), w ∈ {j, j′}, s ∈ {0, 1}.
- **Per-witness record axes** (round-3: one shared axis spuriously VOIDs on a benign j′
  frame offset): û_w = unit of b₀(w) − b₁(w) at the r = 0 arm, per witness, fit on a
  SPLIT SAMPLE (half the r = 0 shots PER VARIANT fit the axis, the other half evaluate the
  reference Ŝ(w; 0): removes the in-sample over-alignment bias at source; partition rule
  PINNED (round-11): EVEN shot indices in delivered memory order fit the axis, ODD
  indices evaluate the normalizer, deterministic, no RNG (primary purpose:
  out-of-sample axis fit; the even/odd drift-interleaving benefit is SECONDARY and
  unguaranteed; delivered order need not equal temporal order, round-19); the models
  replicate the exact rule; EXECUTION REQUIREMENT (round-18):
  the split needs PER-SHOT MEMORY (memory=True, stable delivered order); aggregated
  counts have no shot indices, so a counts-only runner would make this verdict code
  UNEXECUTABLE; per-shot memory is a hard Class-1 requirement and a runner-VERIFY line;
  pipeline order PINNED: split the raw memory even/odd AT THE TWO r = 0 ARMS (the only
  split-bearing PUBs, round-20 scope: all other arms use their full shots) → aggregate
  each half to counts → CAL0/CAL1 inversion per half → conditioning per half
  (mitigate-then-split is impossible; the inversion operates on distributions); under the split Ŝ(w; 0) is the eval-half difference projected
  onto the independently fitted axis: Gaussian, positive with overwhelming probability at
  the r = 0 SNR, NOT the folded magnitude and not positive by construction; this is THE
  definition of the normalizer, there is no other). Transverse
  guards T̂(w; r) = ½·(b₀(w) − b₁(w))·û_w⊥ (the Ŝ construction projected on the
  perpendicular, each witness's own û_w⊥; written out round-18); the SIGNAL model
  injects per-qubit coherent
  RZ frame offsets so the guards' bands cover the drift DIFFERENTIAL (statics are
  absorbed by the r = 0 axis fit; the guard's operative duty is inter-arm drift,
  round-19).
- **SIGNED record projection** (round-3 constructive fix; zero floor, Gaussian at the crux):
  Ŝ(w; r) = ½·(b₀(w) − b₁(w))·û_w. The folded D̂ = |Ŝ| appears only where a magnitude is
  required; in the VERDICT CARRIER that means the denominators of the inner ratio;
  the informational-only R̂ = D̂(j; r)/D̂(j′; r) additionally uses D̂(j) as a numerator
  (scope stated round 30, spec n2). (Notation: D̂(w; r) with a hat is this
  folded ESTIMATOR; the proof's un-hatted D_j is the physics distinguishability.) The r = 0 normalizer is the split-sample
  projection defined above: the sole definition (the stale "½|b₀−b₁| by construction"
  reading was removed in round 6).
- **VERDICT CARRIER, the signed double ratio:**
  ρ̂(r) = [Ŝ(j; r)/D̂(j′; r)] / [Ŝ(j; 0)/D̂(j′; 0)].
  The inner ratio cancels slow common drift; dividing by the r = 0 arm cancels the static
  write-edge/qubit asymmetry. (Cost stated, round 30, physics NIT: the inner j′ ratio
  buys the drift cancellation at a ~√2 SE price per arm relative to the bare Ŝ(j)
  projection; accepted, the drift systematic is the larger enemy.) What remains one-sided in ρ̂(r > 0) is the watcher factor η,
  which BUNDLES the two watcher gates AND the r-arm depth differential (j′ and S idle
  ~140 ns extra during the watcher layer relative to the 2-gate r = 0 arm: ~10⁻³ per side
  at calibration T2 ≈ 136 µs (the day-of floor is T2* ≥ 70 µs; at a 140 ns idle the two
  give 1.0e-3 vs 2.0e-3: same order; the S side's DEPHASING part is exactly ZERO by
  F136(iv): "γ on dressers, watchers, or any other traced site is exactly
  invisible" to the pointer record (round 31, prior-work: a gated law zeroes what
  this bundle budgeted; only S's T1 part and the j′ side survive, both Aer-absorbed
 : conservative padding, kept)), three orders above the write-edge duration residual
  ~10⁻⁶ (an order-of-magnitude assertion, round-24: non-load-bearing, it sits three
  orders below the idle term and both are Aer-absorbed); both absorbed
  into η_nom(r) by the Aer parity on the actual transpiled circuits, not double-counted).
- Factor collapse floors, D̂(j′; r) (min over the EIGHT r > 0 arms, the scope pinned
  at the guard bank), D̂(j′; 0), and Ŝ(j; 0), each
  factor of the double ratio floored separately (D̂(j′; 0) is algebraically a NUMERATOR
  factor: these are noise-quality floors on every NORMALIZER factor, not only
  divide-by-zero guards; the signed numerator Ŝ(j; r) is deliberately unfloored, it is
  the signal and may be zero or negative; the Ŝ(j; 0) floor is on the SIGNED value,
  Ŝ(j; 0) > floor > 0, which simultaneously guarantees the double ratio's sign
  convention: never |Ŝ|, round-19), are VOID triggers derived from the signal
  model's lower p0.03 quantile (round-20 propagation of the round-19 flip: every one of
  the 8 guard-bank tests, floors INCLUDED, sits at the 0.03% tail; a terminal VOID on a
  valid device deserves the deep tail, and the null routing at ~10²σ is indifferent to
  it) minus that quantile's own fitted SE (no free margin knob; σ disambiguation,
  round-25: the ≥ 10σ sanity pin below is in SHOT-SE units of Ŝ(j; 0), while the p0.03
  cut lives on the full systematics-broadened distribution: two different scales, both
  satisfiable at the ~1.2×10²σ separation; bare σ in power margins elsewhere = the
  statistic's own SE). **The floors' band runs
  inject the readout bias at the sign WORST FOR EACH FLOOR (round-23 BLOCKER fix: the
  uniform hardening-sign injection RAISED the Ŝ(j; 0) floor, so a valid device with
  clean or easing-sign readout correlation, an admitted configuration, sat below the
  frozen floor and was terminally VOIDed with probability up to ~0.8–1.0, with every
  gate structurally blind because no gate ever ran a valid device at that sign; the fix
  is the fewer-false-void rule applied to the injected-bias SIGN: for Ŝ(j; 0) the
  easing/lowering sign; equivalently, each floor freezes from the lower envelope over
  the sign axis {−2%, 0, +2%}).** **Round 29 extends the floors' lower envelope to the
  remaining axes (stats M3 + m3): the j-vs-j′ calibration-spread DRAW is a mixture and
  a device has ONE fixed j′; a fixed bad-j′ valid device against a mixture-frozen
  floor read ~0.26% terminal false-VOID per floor at a ±1% spread, and TWO of the three
  floors ride j′, so FIXED-j′-CORNER banks (j′ parameters at the draw's admitted
  extremes, no draw) and the grid freeze halves join the envelope, since round 31
  the floors' lower envelope runs over THE COVERAGE GRID, the one definition at the
  band procedure (box, η-top, and ε_w-sign members included); the per-arm
  D̂(j′; r) min is thereby also slope-covered. The round-25 sanity pin on Ŝ(j; 0)
  (≥ 10 SHOT-SE below the worst valid mean, ≥ 10 above the null) becomes a
  CONSTRUCTION, not a check, in its below-valid leg (framing repaired round 30,
  stats m13: a p0.03 cut sits ~3.4 of its own SDs below its source mean BY
  DEFINITION; that is what a 0.03% quantile IS, not a defect; the real reason to
  combine is a UNIT mismatch: the pin is in SHOT-SE of Ŝ(j; 0), several times
  smaller than the systematics-broadened bank SD the envelope lives on):
  floor_Ŝ(j;0) = min(the envelope, worst valid mean − 10·shot-SE), with the BANKS
  NAMED (round 30, spec m2): the worst valid mean is taken over the held-out
  worst-basis sign banks AND every grid freeze half; the shot-SE is read from the
  fixed-j′-corner banks (fixed device, spread 0; the systematics-broadened bank SD
  is NOT the pin's unit); the ≥ 10-above-null leg REMAINS A CHECK, measured against
  the null bank in the same shot-SE unit;
  lowering is the fewer-false-void side, and the null routing keeps ~1.2×10² shot-SE
  (review MC: a floor at 0.80–0.90 routes the null identically, 0/2×10⁵ escapes);
  the loosened floor is additionally certified against a NEAR-DEGENERATE VALID
  device, not only the null (round 30, stats m13: the fixed-j′-corner held-out banks
  run through the GATES: power and false-VOID measured there, see the sim gate
  section).**

## Prediction and systematics

    ρ(r) = cos(r·π/2) · η(r)   (SIGNED),   η(0) ≡ 1,   η(r > 0) = the one-sided watcher factor
    fixpoints: ρ(0) = +1, ρ(1) = 0, ρ(2) = −η(2) (flown; −1 ideal); generic ρ(1/2) = +(√2/2)·η(1/2), ρ(3/2) = −(√2/2)·η(3/2)

Reading aid (pinned): η = the TRUE device factor; η_nom = its Aer estimate (targets);
η_min = its worst-admitted bound (feeds b_forgive; the bound is a physics number, the
band is its p0.13 quantile; η is INCOHERENT-ONLY throughout, η_nom, η_min, η_max, the
λ-map, per the round-23 no-double-count guard, propagated round-24: coherent residuals
are covered SOLELY by the ε_w ledger, the ±s_res axis, and the sign-swept gate, NEVER by
η headroom; the 0.92 lower edge is margin for INCOHERENT modeling error below the ≈ 0.96
abort-implied nominal, with the RULE pinned (round 31, spec M5: 0.92 was asserted
without a derivation while feeding b_forgive directly): **η_min ≡ 1 − 2·(1 −
η_nominal-worst)**, i.e. the abort-implied worst-case incoherent LOSS is DOUBLED as
the modeling-error allowance (1 − 0.96 = 4% → 8% → 0.92); if the binding freeze's
abort-implied nominal differs from 0.96, η_min moves with it by this rule, and a
device below η_min lands INCONCLUSIVE, never
false-CONFIRMED);
η_max ≡ 1, the IDEAL ceiling (round 29 note on the sentence below: "ANGLE-DEPENDENT and
worst at r = 2" was the PREDICTION; the round-28 representative-line measurement read the
fractional-rzz duration ANGLE-INDEPENDENT at 68 ns, the ~r-flat case the VERIFY list
anticipated; freezing b_forgive from the r = 2 arm is conservative in BOTH cases, which
is why the pin does not move; the binding freeze re-reads this on the freeze-time line)
(round-21 BLOCKER fix: the round-20 calibration-MEDIAN
computation contradicted "best-admitted"; the lexicographic line selection flies the
BEST edges, so the flown device sat systematically ABOVE the swept range, silently
reopening the round-15 power leak with the gate structurally blind to it; no abort
bounds a device from ABOVE, so the ideal edge is the only coverage-complete upper end;
the sweep runs [η_min, 1], the mid-range target = its midpoint ≈ 0.96, minimax over the
full COVERAGE family (round-25 label; the admitted family's own midpoint would be 0.95),
and the freeze-time estimation this replaces can no longer be wrong). Hat convention throughout: ρ̂/Ŝ/D̂/σ̂/R̂
are ESTIMATORS; un-hatted ρ/D_j are predictions/physics. η(r) is angle-dependent and
worst at r = 2 IN THE PREDICTED CASE (each half-gate at π/2, the top of the range;
the round-28 representative-line read was ~r-FLAT at 68 ns, hedged here too, round
30: either case leaves the r = 2 freeze conservative):
b_forgive is frozen from the r = 2 arm specifically; the abort's RZZ-error metric is
referenced to the π/2 gate. η_min from worst-admitted calibration + differential
witness-idle/T2, via the SIGNAL model's incoherent injection; η_nom(r) read off the Aer
parity (the η_nom analog of the τ two-reading pattern, round-14: frozen from the
representative line; the flown-line difference rides clause (c)'s fail-safe direction,
stated there).

**The r = 1 point is first-order fragile:** signed ρ̂(1) = −sin(ε_watch + ε_ZZ)/cos(ε_ZZ)
EXACTLY + Gaussian noise (round 29 exactification of the "≈ −sin ε_w" gloss: the static
ε_ZZ acts on the r = 0 reference arm too and pays the 1/cos(ε_ZZ) factor; a pure-static
miss reads −tan; at the budgets the gap to −sin ε_w is ~1.5×10⁻⁶, ~0.0002 SE
(figure corrected round 30, physics m1: the earlier ~1.4×10⁻⁵ was computed at
ε_ZZ = 0.03 = 3× the pinned budget; at ε_ZZ = 0.01 the gap is
sin(0.03)·(1/cos(0.01) − 1) ≈ 1.5×10⁻⁶),
anti-conservative on b_blind's center, and the signal model injects the exact channel
anyway; the sign consistent with clause (a)'s center), ε_w = ε_ZZ + ε_watch
(ε_watch = 2·ε_cal, the circuit section's two half-gates; ε_ZZ ≡ θ_static = π·Δf·τ, the
static-ZZ equivalent watcher angle of the ζ book, this passage). **The ζ factor-2 book, pinned (round-3 physics):** the price-pair
conditional-Ramsey observable is the fringe-frequency DIFFERENCE Δf (kHz). For H = ζ_P·Z_jZ_k
the conditional angular frequencies are ±2ζ_P, so 2π·Δf = 4ζ_P, and the equivalent watcher
angle over window τ is θ_static = 2ζ_P·τ = **π·Δf·τ** (at Δf = 3.9 kHz, τ = 400 ns:
0.0049 rad; the wrong book, 2π·Δf·τ, would double it, and the wrong book is not
hypothetical, round 31 prior-work BLOCKER, reconciled here BY NAME: the price-pair
run-4 record
([PRICE_PAIR_HARDWARE_PREDICTION](PRICE_PAIR_HARDWARE_PREDICTION.md), RUN 4) defines
ζ as the conditional |0⟩-vs-|1⟩ fringe-frequency DIFFERENCE, which for H = c·Z⊗Z is
2c/π, so c = πζ/2, consistent with that record's own cos(πζt) attenuation and with
THIS book; but [IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md) pins
"H_ZZ = πζ·Z_iZ_j, the convention of the price_pair conditional-Ramsey measurement"
and its gate injects rzz(2πζτ) per step; at the transferred 3.6–3.9 kHz numbers
that DOUBLES the physical coupling, so F129's b_zz2 bias budget carries an
unintended extra 4× beyond its deliberate 2×-ζ headroom; re-read performed at this
fold: F129's clause-(a) verdict SURVIVES the correction (A0 deviation +0.0114
against a corrected window upper edge 3σ_a + b_zz2/4 + b_qs = 0.0192, and 0.0182
even at b_zz2 = 0), while its v2.6 gloss "b_zz2 ≈ 38% of the excess" reads ~9%
corrected; the F129-side repairs ride this file's freeze commit; the 3.9 kHz
context value is a MARRAKESH measurement (the price-pair confirmation) quoted here
as a representative magnitude only, never a frozen input; the flown-line ζ_jk is
measured in-job, per this section). ζ_jk is MEASURED in-job
(conditional-Ramsey PUB pair on (j, k)); **the FROZEN ε_ZZ budget is pinned to the
Class-2 abort ceiling itself, ε_ZZ budget ≡ 0.01 rad (round-27: a frozen band cannot
take the in-job Δf as an input; the in-job read feeds ONLY the Class-2 abort and the
decoding table, never a frozen band; the abort then guarantees every non-VOID flight
sits at or below the budgeted ceiling, ~2× the representative 0.0049 rad)**
(single-ε_ZZ assumption noted, round 30, physics NIT: the r = 0 arm runs 2 gates
against the r > 0 arms' 4, so its always-on-ZZ window is shorter and the effective
ε_ZZ differs by ~1.7×10⁻³ rad between the reference arm and the watched arms, inside
the ledger's worst-case sweep, which budgets the full window on every arm); **τ has two pinned readings (round-9: the flown
line does not exist at freeze):** τ_freeze from the REPRESENTATIVE-line transpiled
schedule's full H-to-pre-rotation window (one number; the band/ledger input), and
τ_dayof RE-READ from the FLOWN line's actual schedule (the abort input: the product
abort with day-of values guarantees the ledger ceiling regardless of line-to-line
duration differences). ZZ during readout is harmless: the pre-rotation freezes the
record into Z populations and ZZ is Z-diagonal.
Class-2 in-job VOID gate (see Abort gates): π·Δf_jk·τ_dayof ≤ 0.01 rad. ε_cal ≤ 0.01 rad
per gate, anchored IN-JOB by a
watcher-angle amplification PUB pair (N = 8 repeated RZZ(π/2) on (j, k), Ramsey-read on j:
amplifies a coherent angle miss 8×; Class-2 in-job VOID gate |8·ε_cal| ≤ 0.08 rad): the coherent-angle
guard the incoherent EPC metric cannot provide; coverage assumption stated (round-9): the
amplification is π/2-referenced, so the fractional-angle gates (r < 2) carry no separate
coherent anchor; their deviations are covered by the fail-safe direction (they inflate
|ρ̂(1)| toward INCONCLUSIVE) plus the duration-vs-angle VERIFY; the day-of "RZZ error
≤ 0.5%" is pinned as gate INFIDELITY (EPC), a separate quantity. Residual ε_cal drift risk is fail-safe:
it inflates |ρ̂(1)| toward INCONCLUSIVE, never toward a false CONFIRMED.
ρ̂(2) is second-order flat in ε_w. **The amplification pair's READ pinned (round 29,
physics M2, the k preparation was unstated, and the natural |+⟩ preparation reads
cos(8ε): second-order BLIND at ε = 0, the same blindness class the round-4 V_S fix
removed, 12.5× less sensitive at the threshold): the PAIR is k prepared in |0⟩ and in
|1⟩ (Z eigenstates), the j read is the conditional-Ramsey phase (sin route, SIGNED,
first-order at ε = 0), the two preparations combined in ONE pinned equation (round
31, physics MAJOR: the branches are ANTISYMMETRIC, ⟨Y_j⟩ = ±sin(8ε_cal), so their
raw DIFFERENCE is 2·sin(8ε_cal); the old wording here, "taken like the ζ pair's
branches", would read the abort quantity at twice its value and abort at HALF the
budget; the factor-2 trap is the ζ book's sibling and gets the same
one-equation treatment): **½·[⟨Y_j⟩_{k=|0⟩} − ⟨Y_j⟩_{k=|1⟩}] = sin(8·ε_cal)**, and
the Class-2 compare is |8·ε_cal| from THAT half-difference against 0.08 rad;
a runner-VERIFY line asserts the preparation, the signed extraction, AND the ½
factor (beside the Δf-is-a-difference line).**

**The V_S guard (the CQ premise, guarded not assumed, round-4 order fix; FRAMING
sharpened round 31, prior-work: in PROOF_RECORD_PARITY_LAW V_S = 0 is a DERIVED
consequence of Law A's hypotheses, equal write couplings, deg(S) ≥ 2, exact t*,
not an independent premise; guarding it in-flight therefore guards the CALIBRATION
INPUTS to that derivation, which is exactly what a hardware guard is for):** the premise
for the (S, j) record is the per-pair V_S^{(j)} = cos(θ_Sj′), FIRST order in a write-angle
miss δ (θ_Sw = the (S, w) write angle, nominal π/2: the proof's θ_w per pair); the unconditional ⟨X_S⟩ = cos(θ_Sj)·cos(θ_Sj′) is O(δ²) and would be ~1/δ blind:
~100× at the 0.01-rad abort scale (20× already at δ = 0.05; the quoted ratios assume
the miss on BOTH writes, giving sin²δ; a single-write miss leaves the unconditional
read identically 0, formally infinitely blind, round-22), so the guard measures the **j-conditional equatorial S-coherence**: guard
PUBs with S pre-rotated to X (and Y) and the witnesses measured in Z, at r ∈ {0, 2}; the
statistic |b(S | z_j)| is first-order in δ (prediction 0; threshold from the SIGNAL model
at the budgeted δ = ε_cal = 0.01 rad, round-25 pin: the write gates share the watcher's
coherent-angle budget; the guard bounds false MI-denial on a valid device; the V_S
threshold freezes from the upper-guard COVERAGE ENVELOPE: the round-28 rule, which
SUPERSEDES the round-26 worst-sign-only wording this site previously carried (site
refreshed round 30, spec m3); S-adjacent and
witness-readout-conditioned, hence plausibly sign-sensitive; it gates only the MI
corollary, but the one-clause rule costs nothing; its tail family is assigned to the
GEV/empirical branch explicitly, round 30, stats m12: a folded 2-vector norm, never
Gaussian-fitted).
Exceedance → the MI corollary is NOT claimed (the record-radius verdict stands on its own;
the flight report separates the two sentences).

**The watcher-tomograph pair (the mechanism certificate, round 30, physics BLOCKER):**
every science and guard PUB above is Z-DIAGONAL on k (all gates commute with Z_k; k is
measured in Z or not at all), so the entire flown set sees only diag(ρ_k), and
|+⟩⟨+| and I/2 share the diagonal: a maximally mixed, never-watching k (m = 0) is
INVISIBLE to the whole design and reproduces the signed curve exactly (a classical
mixture of ±rπ/2 phase kicks multiplies j's coherence by cos(rπ/2), the law's own
factor). The separating correlator is ⟨Y_k Z_j⟩: +1 at r = 1 for a SUPERPOSED watcher
(the claim section's perfect classical bit), identically 0 for ANY Z-diagonal k, at
every m. Two PUBs measure it: the r ∈ {0, 1} science circuits with k pre-rotated
**Sdg then H** (reads Y_k; same order convention as variant B) and j, S, j′ measured
raw in Z; statistic = ⟨Y_k Z_j⟩ at the r = 1 PUB, prediction +1 (dressed), and ROBUST, not
sensitive (order of response corrected round 31, physics m1, before the vocabulary
repeated the round-29 amplification confusion in the opposite direction: from below
⟨Y_k Z_j⟩ = sin(rπ/2 + ε_w) exactly, so at r = 1 it sits at the sine's extremum,
SECOND-order flat in ε_w, cos(0.03) ≈ 0.99955 at budget, and it is EXACTLY
write-miss-independent; flatness is a virtue for a lower cut pinned at +1);
~10 s of budget. DIAGNOSTIC-ONLY, exactly like V_S: the frozen
tomograph threshold (a LOWER cut, p0.03 − SE lower envelope over the coverage sources)
gates ONLY the mechanism sentence, "the watcher demonstrably watched in superposition",
never the record-radius verdict; if the pair were ever dropped, the flight report
carries a Not-tested line: the watcher's superposition is then unverified and the
flight tests the record-radius law on j alone. The r = 0 tomograph PUB is the
reference/consistency read (prediction 0).

## Arms and shots

r ∈ {0, 0.25, ..., 2.0} (9 arms) × 2 variants = 18 science circuits + CAL0/CAL1 + Z-leak
pair + ζ_jk conditional-Ramsey pair + V_S guard PUBs (S in X and Y × r ∈ {0, 2} = 4) +
watcher-angle amplification pair (the |8·ε_cal| abort input) + S↔witness crosstalk pair
(round-9: S = |1⟩ rest |0⟩ and S = |0⟩ rest |1⟩: the S-correlated witness assignment
differential, the day-of input of the readout-correlation abort) + watcher-tomograph
pair (round 30: r ∈ {0, 1} with k pre-rotated Sdg·H, the mechanism certificate, its
own paragraph above) + in-situ T1/Ramsey
PUBs (diagnostic/provenance ONLY, round-21: they cross-check the calibration values the
aborts consumed and feed the flight report; they gate nothing); ~36–38 total, exact
count pinned at freeze (ranges are draft placeholders; the
committed pre-reg carries the frozen numbers). **Interleave order PINNED (round-18: the drift-slot anchor needs a deterministic
schedule, not a rule family; slots of r = 1.75 and r = 0.25 SWAPPED round 30, physics
m4: the round-18 order put r = 1.75 at the MAXIMUM drift distance d = 4 from the
r = 0 anchor, which EASES the clause-(b) raw ordering |ρ̂(2)| > |ρ̂(1.75)| by ~+0.012
at budget drift and suppresses the fallback center ~2%; drift must never ease a
CONFIRMED-bearing comparison):** per variant, the fixed arm order
[1.0, 0.5, 1.5, 1.75, 0, 0.75, 2, 1.25, 0.25]: r = 0 at slot 5, r = 2 at slot 7 (the
two-slot separation the drift-leak figure anchors to), r = 1.75 at slot 4 (d = 1 from
the anchor), the maximum-distance slots d = 4 now carrying r = 1.0 and r = 0.25,
neither ENDPOINT fixpoint (r = 0, r = 2) first (wording repaired round 31, physics
m3: r = 1.0 IS a listed fixpoint, ρ(1) = 0, and it sits first, and at d = 4 the
budget drift LOWERS |ρ̂(1)| by ~6×10⁻⁴, an (a)-EASING direction, stated and
accepted: 0.7% of the b_blind envelope, and clause (a)'s false-CONFIRMED duty is
~zero with the swept controls certifying the residue). The
inter-PUB order (CAL/guard/diagnostic PUBs among the science arms) is DELIBERATELY
unpinned (round-23: delivered order need not equal temporal order anyway, the drift
benefit is secondary, and the drift systematic is day-of-anchored through the
transverse guard, not through scheduling). Shots
**16384/circuit** (round-8 stats: doubled from 8192, mean-estimator SEs ÷ √2, exact;
tail quantiles approximately; round-12 relabel
of the absolute figures: the ~2.7σ → ~3.8σ pair describes a SUB-ABORT stress ratio
(η(2)/η(1.75) ≈ 0.96, i.e. η(2) ≈ 0.74–0.80, below the abort floor); at the
abort-implied η the raw (b) margin is ~8–12σ for a TYPICAL device (round-20 recompute) (~3–6σ at the exact
abort corner with worst-case hardening bias, still resolvable; the measured joint power
governs), so the doubling is retained as cheap margin, not necessity, qualified round-17: margin on (a)/(b)/contrast only, NOT on
(c)'s corner coverage, which is systematic-limited (the corner-conditional b_curve is
that fix; more shots are not); the stress-floor "joint power ~96%" figure was
self-consistent stress-frozen bands; against the OPERATIVE abort-corner bands a
true-0.74 device lands INCONCLUSIVE and would abort day-of anyway; the operative gate is
the MEASURED joint power at freeze; ~1 extra QPU minute on a flight that is ~2% of the
annual gift); the two r = 0 science arms carry
**2× shots (32768)** (round-9: the split-sample normalizer runs on half its shots and
feeds EVERY band as a common factor; doubling r = 0 restores full nominal precision on
both split halves at ~+9 s); shots and cost re-freeze together if the power check
demands more (recorded).

## Verdict rule (VOID precedence; ONE conjunction; pinned freezing rules)

**Precedence (full order, evaluated in sequence): VOID → CONFIRMED → DEVICE-DEVIATION →
INCONCLUSIVE.** Any trip of the GUARD BANK below (named precisely, round 30: the bank
holds 5 POOLED statistics plus 3 INDIVIDUAL floors; "pooled guard bank" was the old
shorthand and misnamed the floors) voids regardless of the physics
clauses (the V_S guard is a separate track: it gates only the MI corollary, never the
record-radius verdict); DEVICE-DEVIATION is
assessed only after CONFIRMED has failed; the catch-all takes the residue. **The guard bank
is POOLED (round-5 stats: per-(stream × arm × variant) counting gives ~90 thin-tail tests
and a ~11% union-bound false VOID, unacceptable on a paid one-shot):** one pooled
transverse statistic per witness (the SUP of T̂(w; r) over ALL NINE arms, r = 0
included; its transverse read comes from the eval half against the fit-half axis;
arm set pinned round 31, spec M6: leaving it open moves both the frozen band and a
terminal trip probability, and the 8-test accounting is untouched because the SUP
is ONE pooled test whatever its arm count, banded as a SUP), one
pooled ⟨Z⟩ sanity statistic each for k and S ONLY (round-7 fix of a j↔S slip: k and S are
the two qubits MEASURED in Z, both ideal 0; H then Z-diagonal gates leave their
populations at ½; the pre-rotated witnesses j and j′ have no Z stream, their measured
axis IS X/Y; the |⟨Z_S⟩| stream doubles as the branch-balance sanity on the conditioning
qubit), the Z-leak pair as one statistic (DEFINED: two extra circuits = the r = 0 and
r = 2 arms with NO pre-rotation on any qubit, all four measured raw in Z; statistic = max
over witnesses {j, j′} and both circuits of |⟨Z_w⟩|, ideal 0 for the same Z-diagonal
reason (k and S are measured too but excluded from THIS statistic; their Z streams
already feed the pooled sanity guards from the science arms, and double-counting them
would break the 8-test accounting, scope stated round 30); leak = rotation/T1 error
moving population), and THREE
individual floors: D̂(j′; r) (its min over the EIGHT r > 0 arms, scope pinned round
29: r = 0 carries its own dedicated floor at 2× shots, and double-counting it would
break the 8-test accounting), D̂(j′; 0), and Ŝ(j; 0) (each factor of
the double ratio floored separately; a jointly-small pair must not hide). 8 tests (5
pooled statistics + 3 individual floors);
the sim gate MEASURES the aggregate false-VOID empirically (correlations included), pinned
ceiling 0.5% MEASURED at the p0.03 default (round-21: the old 2% was measured at the
p0.13 default and would re-admit the terminal loss the flip eliminated; naive union
0.24%, measured expectation ~0.3–0.5%; this 0.5% is the false-VOID ceiling; the day-of
RZZ EPC abort's 0.5% is an unrelated constant, round-22); guard quantiles DEFAULT to p99.97/p0.03 across
ALL 8 tests, floors included
(rounds 19–20: statistical-guard VOIDs are TERMINAL on a paid one-shot; the earlier
default lost ~2% of valid flights before any remedy engaged; provenance of the two
figures, round-20: the ~2% was MEASURED through the frozen code, sitting ABOVE the naive
1.04% one-sided union bound precisely because fitted tails under-covered the folded
statistics: the failure the GEV/empirical branch addresses; the ~0.24% is the naive
union at the new default, 8 × 0.03%, and the MEASURED ceiling still governs; the guards
are sanity checks, never verdict content); remedy ladder: fixpoint-arms-only
guards → re-freeze
at raised seeds (each step recorded; PINNED EXCEPTION, rounds 13–14: no remedy step may
lower the Ŝ(j; 0) floor; it is the guard that ROUTES the degenerate null to VOID before
CONFIRMED is evaluated (the null's double ratio is genuinely Cauchy: ~27–29% of
zero-tilt null runs land below −b_forgive unprotected: a worst case, injected frame
tilt shrinks it; the shape clause independently rejects the null with probability ≈ 1,
so the floor is routing, not the sole protection; law-vs-null separation on Ŝ(j; 0) is
~1–2×10²σ, overwhelming; the earlier 388 figure was the law's own SNR, a different
quantity); the floor is EXEMPT from the ladder entirely and stays where the signal
model froze it: lowering it walks toward the Cauchy regime, and the ladder's loosening
steps apply to the OTHER guards only (round-15 rewording of a self-contradictory v13
direction rule)). VOID: NaN/failed fit/mitigation
abort; factor floors; per-witness transverse guards T̂(w); Z-leak pair; |⟨Z_k⟩| and
|⟨Z_S⟩| sanity streams (ideal 0); the V_S guard trips only the MI corollary (stated above),
not the record-radius verdict.

**Band-source architecture (round-6 rebuild; the earlier null-sourcing was DEGENERATE:
with the write zeroed, the double ratio's denominator Ŝ(j; 0) is pure noise, its bands
unestimable at real SNR and Cauchy-like where estimable):**
- **SIGNAL model** = the law-holds construction with every injected systematic at its
  budgeted value (the readout-bias injection at the CONFIRM-HARDENING sign here; the
  easing sign lives only in the negative controls, round-12). Sources ALL fluctuation
  widths: the b_blind tail (read at the r = 1 arm:
  write present, watcher erasing, denominator intact, clean Gaussian, verified),
  SE_contrast (its freeze takes the LARGER, conservative side across the sign banks,
  inflation direction pinned round 30, stats n15), b_forgive, b_curve, b_sym, η_min,
  the denominator floors, the guard
  bands (a guard's job is bounding false VOIDs on a VALID device; b_sym listed round-22
 : its omission was an oversight, it sources like every width and stays
  informational-only), the LR thresholds' reference distribution, the clause-(b)
  fallback band, and the watcher-tomograph threshold (the last three listed round 30,
  spec m4: the sources list must be exhaustive over frozen constants).
- **Deterministic centers** (the b_blind center max over the budgets of
  |−sin(ε_watch + ε_ZZ)/cos(ε_ZZ)|: the exact r = 1 form, round 30 spec n1; the
  |−sin ε_w| gloss differs by 1.5×10⁻⁶ at budget): the
  systematics ledger: a deterministic worst-case sweep, no model sampling.
- **NULL model** (the no-record device: the (S, j) write's conditional phase zeroed, all
  else identical) sources NO band. Its one job is the **mandated NEGATIVE control**: the
  gate measures, through the FROZEN verdict code, P(CONFIRMED | null) and
  P(CONFIRMED | each impostor), each pinned ≤ 0.1% empirical, acceptance rule (round-8):
  the 95% Poisson upper confidence bound of the measured rate must sit below 0.1% (a raw
  measured 0.1% at 10⁵ seeds is under-resolved: ~100 events carry ~10% count error); the
  null and each impostor run WITH the worst-case S-correlated readout-bias injection at
  its 2% ceiling, SWEPT over the sign axis {−2%, 0, +2%} AND over the CONFIG axis
  {worst-admitted, ideal, fixed-j′ corner, box corner, η-top-null}: the control
  AXIS SET pinned (round 30, spec M6, EXTENDED round 31, spec BLOCKER: the round-30
  fold widened the verdict bands over the corner sources while the controls still
  swept sign × level ends only: a control certifies the SAME band set that is
  committed, so every config the bands were widened over joins the control sweep;
  the corner configs run device-FIXED, spread 0; slope's effect on a NON-LAW device
  remains covered by the level ends, pinned; interior levels are an OWNED assumption,
  round 31 spec m4, stated rather than silently monotonized: η only attenuates a
  non-law curve, every listed alternative dies at a sign/flip/guard requirement that
  attenuation cannot manufacture, and the certification is the measured zero-event
  record across the full config × sign sweep, never a monotonicity theorem), with
  the pinned ≤ 0.1% UCB
  required PER BANK at EVERY swept point (per-bank, never pooled: pooling across
  banks would let a hot bank hide, the juxtaposition stated round 30, stats
  n14) (round-27: "swept, not assumed"; rounds 11–12 identified the
  CONFIRM-EASING sign as the one mode easing clause (b) and ran it alone; the round-27
  review MC confirms that sign is empirically the more adverse, 94.0 vs 103.5 SE from
  the (b) boundary, so the sweep costs cheap extra banks and changes nothing;
  this measured bound IS the false-CONFIRMED protection for that
  mode; round-12 MC: 0/20000 across null + all impostors, both signs, unshifted bands;
  round-13 MC: still 0 at a 0.10 bias, above the 9-SE coefficient high end). The
  Z-POLARISED never-watching class (impz, H on k omitted) is the FOURTH control kind
  (round 30, spec M2+M3: the class previously sat in NO control bank and NO detection
  gate while the scoping sentence leaned on its guard): impz runs in the same
  sign × level control sweep at m = 1 under the same ≤ 0.1% per-bank UCB; because the
  LR cannot see it, the gate ADDITIONALLY measures its DETECTION, P(VOID via the z_k
  guard | impz, m = 1) ≥ 0.999, asserted, and a BOUNDARY control runs impz at
  m = the frozen z_k band value, its P(CONFIRMED) REPORTED as the honest escape rate
  of the discrimination window (round 30, stats BLOCKER; the m ≈ 0 member is scoped
  to the watcher-tomograph diagnostic and the mechanism sentence, never to the record
  verdict); seed floor
  pinned (round-13, exactified round 31 stats n2: below 2996 seeds, 2.9957/n <
  0.001, exact binomial agreeing, the strict < 0.1% UCB is unsatisfiable even at
  zero events); the 10⁵ default governs, never below 10⁴;
  symmetric with the positive control (a null run correctly lands in VOID/INCONCLUSIVE,
  never CONFIRMED).

**Band procedure:** full-estimator models per the source architecture above, default
**10⁵ seeds** (sim seeds are cheap; a 0.13% tail at 2000 seeds has ~62% Poisson error on
the tail COUNT, ~6% on the quantile value; the raised default plus the parametric fit
removes the extrapolation reliance); every PHYSICS band (b_blind, b_forgive, b_curve,
the LR threshold) is the
fitted tail law's **p99.87 / p0.13 quantile, inflated by that quantile's own fitted SE**,
while every GUARD-BANK statistic, the 5 pooled guards AND the 3 floors, sits at the
**p99.97 / p0.03 default** (round-20 propagation of the round-19 flip; the same fit and
SE machinery at the deeper tail; the per-guard quantile labels are NOMINAL-only,
round-23: a measured aggregate above the naive union proves individual guards trip
above nominal, which is why the MEASURED aggregate is the sole governing quantity)
(bootstrap of the fit; the goodness-of-fit gate must bound the fitted-quantile SE, not only
accept the family; GoF PINNED numerically, rounds 14/17: pass requires Anderson–Darling
p ≥ 0.01 on the fitted tail AND the fitted p99.87 within 2 bootstrap SEs of the
DEFAULT-SEED empirical quantile AT THE STATISTIC'S OPERATING POINT: p99.87/p0.13 for
physics bands against the 10⁵-seed bank, p99.97/p0.03 for guards and floors against
their ≥ 10⁶-seed bank (rounds 21–22: the GoF pin predated the guard flip and the seed
raise; the comparison bank is always the freeze bank, reconciled with the seed-cost
scheme round 31, spec m2: "the freeze bank" means the bank the frozen value was READ
FROM, and under the argmax-deep-seed rule the operative band value always comes from
a deep bank: the worst-basis banks by default, or the argmax grid point re-run at
≥ 10⁶; non-argmax 10⁵ sources never carry the frozen value, so the ≥ 10⁶ GoF demand
attaches exactly where a value is read); either failure → next branch; roles labeled: the
quantile-agreement condition is PRIMARY, the AD leg is only an egregious-misfit screen;
expected operative path (round-18): at the default seeds AD rejects nearly any imperfect
parametric family, so the empirical branch with bootstrap SE is the NORM and the
parametric apparatus the exception; the guard/floor freeze at the p0.03 tail runs at
≥ 10⁶ seeds (round-21: 10⁵ gives ~30 tail events at ~18% count error; sim seeds are
cheap, 10⁶ gives ~300 at ~6%; the physics bands stay at the 10⁵ default; the SEED-COST
SCHEME pinned, round 30, stats m11: a uniform ≥ 10⁶ across every envelope source
would be ≥ 2×10⁸ evaluations, unaffordable as written: the DEEP ≥ 10⁶ tail runs on
the worst-basis sign banks, the grid/corner sources enter the envelope at their
10⁵-per-point scale, the gate REPORTS the envelope-ARGMAX source per guard band, and
an argmax landing on a shallow source is re-run at the deep seeds before freezing:
the remedy-ladder line for freeze-bank cost; deep-seeding off-argmax buys nothing); at a 0.13% tail and 2000 seeds the tail COUNT carries ~62% Poisson
error and the raw empirical quantile VALUE ~6% (round-7 relabel: the two must not be
conflated), which the parametric fit + SE inflation absorbs); the denominator-floor
margins are this same quantile-SE (no free knob). **SE-inflation direction per decision
role (pinned; SCOPED round 31, spec M8: the old absolute sentence and the floors'
rule below self-cancelled on exactly one statistic, and this is the round-30 "self-
cancelling sentence" located and discharged: lowering the Ŝ(j; 0) floor is
simultaneously the fewer-false-void side AND a CONFIRMED-easing move, since that
floor routes the Cauchy null to VOID before CONFIRMED is evaluated): for every
CONFIRMED-bearing BAND, the side that makes CONFIRMED easier is
never taken; the VOID FLOORS follow the fewer-false-void rule instead, and on
Ŝ(j; 0), where the two rules pull opposite ways, the resolution is PINNED: the
floor takes the fewer-false-void (lower) side but is guarded from below by the
≥ 10-shot-SE-above-null CHECK, whose failure HALTS the freeze for diagnosis
(round 31, spec m7: the floor is ladder-exempt, so a failed above-null check
indicts the null bank or the construction, never a threshold to loosen), and the
null routing is CERTIFIED
empirically at every freeze (~1.2×10² shot-SE, 0/2×10⁵ escapes), so the easing
direction is measured shut, never argued shut** ("inflated" throughout = SE-shifted toward the conservative edge, which for
a ceiling like b_blind means the LOWER edge; b_blind takes its lower edge, b_forgive its
upper edge; the projected power absorbs the cost; acknowledged round-25: clause (a)
carries ~zero false-CONFIRMED duty, the impostors sit two orders above the ceiling,
so b_blind's lower edge is pure convention cost, kept for rule-uniformity and certified
by the swept gate; figure re-based round 30, physics m2, superseding the round-29 7.7×/10.9× which used
the pre-envelope 0.065: at r = 1 the impostors sit at
0.5 and 0.707 against the coverage-envelope ceiling, ≈ 6.2× and 8.8× at the
rehearsal envelope 0.0802 (the frozen value governs; 0.5/0.0802 = 6.2,
0.707/0.0802 = 8.8), not "two orders", and the
sign-flat impostor (iii) PASSES clause (a) exactly (its |cos(π/2)| = 0 matches the law
there); the conclusion stands because (iii) dies at clause (b)'s sign and the others at
(a), but the protection is the CONJUNCTION, never clause (a) alone); VOID floors take the fewer-false-void
side. Tail family per statistic: Gaussian ONLY
for statistics whose SOURCE model shows them symmetric with Gaussian tails (the signed
ρ̂(1) qualifies under the signal model, verified excess kurtosis ≈ 0; symmetry alone does
not imply Gaussian tails); any
|·|-folded or near-zero statistic (denominator floors, transverse guards, Z-leak, the
pooled |⟨Z_k⟩|/|⟨Z_S⟩| sanity statistics, and the V_S guard, its folded 2-vector norm
assigned explicitly, round 30, stats m12) starts at the GEV/empirical branch (empirical
fallback: SE = bootstrap of the empirical quantile; at 10⁴ seeds the 0.13% tail COUNT
still carries ~28% Poisson error; the ~2.7% VALUE figure is Gaussian-at-3σ illustrative
only and does not transfer to folded/GEV tails, whose operative quantile-value SE is the
bootstrap itself; inflate accordingly; the guards' operative tail is 0.03%: ~300 events
at the pinned ≥ 10⁶ seeds (round-27 consistency fix: this site still framed the tail at
10⁵ seeds after the round-21 raise), ~6% count error, absorbed by the bootstrap-SE
inflation with a FURTHER raise as the remedy if the SE dominates, round-20);
goodness-of-fit gate pinned (fit fails → next
branch → 10⁵-seed empirical; round-11: the folded guards get the full default seed count
too, removing the thin-tail reliance the 10⁴ figure above motivates). **Band basis SPLIT
BY ROLE (round-15; the single worst-admitted basis protected clause (b)'s magnitude but
POISONED clause (c)'s shape: with angle-dependent η the profile η(r)/η(2) depends on the
η LEVEL, so a worst-edge target makes a BETTER-than-worst valid device deviate
systematically: ~2× b_curve at η(2) = 0.98 vs a 0.92 target; the better the hardware,
the likelier INCONCLUSIVE, and a single-basis power gate is blind to it, the
destructive protection interaction the rounds hunt for):** the LOWER magnitude band
b_forgive stays frozen at the worst-admitted basis, as the LOWER ENVELOPE over the
readout-sign axis + SE up (round 31, spec M2: this site still carried the retired
"hardening sign" label one version after its retirement; the estimator section's
sign-envelope construction is the one rule, restated here by reference; for a lower
cut the worst basis is the level envelope); the FLOORS freeze from their lower envelope (sign axis +
the round-29 extensions stated at the floors section); b_blind, an UPPER cut whose
center |sin(ε_w)|·(record level) rides the level axis exactly like T̂'s drift term,
joins the coverage-range UPPER envelope (round 29, stats m1: the worst-basis freeze was
coverage-correct only through an unpinned cancellation between the center's rise and
the width's fall with device quality; an envelope needs no such coincidence; clause
(a)'s false-CONFIRMED duty is ~zero, so the envelope's easing side is certified by the
swept negative controls); the UPPER GUARD BANDS (the pooled per-witness transverse T̂, the
|⟨Z_k⟩| and |⟨Z_S⟩| sanity statistics, the Z-leak statistic, the V_S threshold) freeze
from the p99.97 + SE UPPER ENVELOPE over the FULL coverage range, i.e. the worst-admitted
basis AND every member of THE COVERAGE GRID (the one definition at the band
procedure; the |⟨Z_k⟩| envelope is additionally CAPPED by the
pinned discrimination ceiling guard_z_k ≤ 0.05, round 30 stats BLOCKER; the ceiling
PINNED PROPERLY round 31, spec M3, a bare number is not a rule: (i) DERIVATION:
0.05 is the accepted-set semantics made numeric, "an escaping Z-diagonal watcher
must be ≥ 95% maximally mixed", sized against the shot-noise-only pooled envelope
≈ 0.034–0.037 at pinned shots (headroom ~1.35×, measured and reported at every
freeze); (ii) SEMANTICS, the ceiling is a FREEZE GATE, never a clip: envelope >
ceiling HALTS the freeze, the band is never silently set to 0.05; (iii) REMEDIES,
named, all envelope-LOWERING: narrow the pooled z_k scope to the fixpoint arms
(the remedy ladder's existing first step), raise the guard seeds (shrinks the SE
edge), or tighten the k readout/T1 Class-1 aborts (shrinks the systematic part);
raising the ceiling is never a remedy) (round 28, sim-gate stage
2026-08-04: the SEVENTH instance of the protection-interaction class, and the first found
by the gate itself rather than a reviewer; T̂'s drift term scales with the RECORD LEVEL,
so the earlier worst-basis-only pin terminally VOIDed the BEST admitted device: exact
distribution-level centers sup|T̂_j| = 0.190 at the worst basis vs 0.235 at the
(η(2) = 1, slope −1%, sign −2%) corner against a frozen 0.244, corner joint power 0.10,
VOID:t_sup_j in 269 of 300 runs (decomposition flagged round 31, stats m4: the
level + sign terms of the stated sup-model account for ~0.204 of the 0.235; the
residue rides the slope term and the folded-tail offset the distribution-level
center carries beyond the sup-model; the pair is RE-READ at the binding freeze and
the frozen artifact carries the decomposition); the fewer-false-void rule thereby extends to the LEVEL
axis: a valid device lives anywhere in [η_min, 1], so the guard's envelope must too; the
FLOORS are unaffected, since for a LOWER cut the worst basis IS the lower envelope,
verified: the post-fix corner breakdown shows no floor trips); the clause-(c) TARGET η_nom(r)/η_nom(2) is frozen at the MID-RANGE
basis (RE-PINNED round-22 in OUTPUT-η space: the round-17 per-parameter recipe broke at
η_max ≡ 1, where a per-parameter midpoint is undefined for time constants: the target
basis is η(2) = (η_min + 1)/2 ≈ 0.96, realized by scaling the representative-line Aer
noise config by ONE loss factor λ: all incoherent rates multiplied by λ, λ solved so
the Aer-read η(2) hits the target; the SAME λ-map realizes every level point of the
sweep, each of the 7 uniform η(2) points over [η_min, 1] solving its own λ, with λ = 0
the ideal end; the map from target η(2) to noise config is thereby deterministic;
η_nom is read at the λ-scaled INCOHERENT-only config with the budgeted coherent
injections OFF, round-23: the coherent budgets live in the ε_w ledger and the injection
list, never inside η_nom: the no-double-count guard); and
b_curve is the CORNER-CONDITIONAL quantile over a TWO-DIMENSIONAL sweep (round-17 fixed
the pooled-mixture under-coverage: p99.87(mixture) ≤ max_η p99.87(·|η), and the corner
offset is shot-INDEPENDENT, so more shots tighten a pooled band and LOWER corner power,
97.1% → 95.7% on doubling; round-18 fixed the remaining dimension: a LEVEL-only sweep
does not cover the SHAPE mode the joint-power gate stresses, review MC: the floor
crossed at a ~1.4% slope residual, on the GOOD device again, and the earlier cited ~93%
shape figure was itself below the floor): b_curve = the MAX over the FULL 3-D grid,
the
η(2) LEVEL over [η_min, 1] (the coverage range, round-22 label: the sweep intentionally
overshoots the admitted ≈[0.92, 0.98] to the ideal, since no abort bounds from above) ×
the SHAPE/slope residual over ±s_res × the readout-bias SIGN {−2%, 0, +2%} (the sign
axis stated at b_curve's own site round 30, spec M1: this passage still read 2-D/35
while the general per-statistic sign rule and the machinery already ran 3-D), of that
point's
p99.87 of the SUP, each SE-inflated (SE-inflation direction pinned, round-22: b_curve
takes its LOWER edge per the global never-ease-CONFIRMED rule; the power protection
lives in the max-over-grid construction and is certified by the swept ≥ 95% gate, never
by loosening the SE step; the max of the coverage grid's ~120 finite-seed quantile
estimates carries a MODEST
POWER-SAFE upward selection bias, E[max] ≈ 2.5–2.6 × the per-point quantile SE
(re-sized round 30, stats m10 at 105 points → 2.53, and re-counted round 31, spec
n2, after the corner members joined; the exact member count is the freeze report's;
at ~210 it would be ≈ 2.76),
do not "correct" it, and this upward step is ACKNOWLEDGED as the one deliberate
exception to the never-ease-CONFIRMED rule (round 31, stats n3: net ≈ +1.5–1.8 SE
upward on an upper cut, under 1% of b_curve, moving the accepted scale window by
~0.03 pp, landing on the continuous neighborhood no negative control bounds, which
is exactly why it is named here rather than left implicit), and the finite grid suffices
because η's level/slope dependence is physically smooth) (grid resolution PINNED, round-20: 7 uniform level
points × 5 uniform slope points × 3 signs = 105 grid points, JOINED (round 30, stats
M2+M3 + physics m3) by the device-fixed corners: the fixed-j′ corners, the
off-diagonal BOX corner (best coherence λ = 0 × worst readout FIXED, the BASE
stated exactly, round 31 stats B1: a reviewer computed this corner from the
backend-median readout and read a 3.2× coverage gap; the actual base is the
WORST-ADMITTED basis, whose readout sits AT the Class-1 abort ceilings, S (1.5%,
2%) and witnesses/k (2%, 3%), and the box corner multiplies THOSE by the draw
factor 1.25, witness up to (2.5%, 3.75%); coverage on the readout axis therefore
EXCEEDS the abort-admitted set, no gap × slope extremes × signs: the λ one-knob
family is a CURVE through
the admitted box, it ties all incoherent rates together, and the box corner is an
admitted device the curve-with-draw only visits transiently; the guard-envelope
arithmetic read up to ~9.6% terminal VOID there under an unpinned width
decomposition), and the η-TOP point (a j′-only differential config, η(2) ≈ 1.002:
η = a_j/a_j′ and no abort enforces a_j ≤ a_j′, so η sits STRUCTURALLY above 1 by up
to ~0.001–0.002; the λ map cannot realize it, the j′-idle config does); a
max-over-grid band needs its grid
fixed; **this list, the 105 level × slope × readout-sign points, the fixed-j′
corners, the box corner, the η-top point, and the ε_w-sign members, IS "the
COVERAGE GRID", the ONE definition, and every band/threshold site in this document
that says "coverage grid" or "coverage envelope" means exactly this set (round 31,
spec M1: six sites had drifted into six different axis enumerations; the per-site
lists are hereby superseded by this reference)**; the FIXED-vs-DRAWN law pinned
with it (round 30, stats M2): every
device-FIXED nuisance enters the freeze as a CORNER at its admitted extremes, only
SHOT noise (and the explicitly-declared calibration-spread mixture) is drawn, the
declaration: edge EPCs, T2s, and witness readout ride the ±25% DRAW as a mixture AND
their extremes enter as corners; the readout-bias sign, the box corner, and the j′
asymmetry extremes are FIXED corners; RZ statics and the drift rate are FIXED at
budget (day-of-anchored through the transverse guard, never drawn); the residual axis is LINEAR slope, assumption stated, round-21: curvature
residuals ride η_nom's per-arm Aer read, so only the device-vs-Aer residual needs
sweeping and a slope axis spans it at the VERIFY-bounded magnitude; a non-monotone
duration scaling found at the runner VERIFY re-reads η_nom PER-ARM, folding the
measured shape into the deterministic target, a linear s_res rescale cannot cover a
non-monotone residual (round-26 remedy fix), and re-derives s_res pre-commit; the SUP's
reference distribution is the (a)∧(b)-PASSING subset (label matched to the code
round 31, stats m6: under the short-circuit the SUP is only ever evaluated after
BOTH (a) and (b) pass, so the conditional quantile's conditioning set includes
clause (a) and its order check; the difference is negligible on a valid device and
the label now says what the code does); the bands are sourced through the
frozen verdict code, which short-circuits (a)→(b)→(c); stated round-25); s_res has the τ-pattern TWO READINGS (round-19: the
runner stage sits AFTER the band freeze, so a runner-read input cannot feed a frozen
band): s_res_freeze = the REPRESENTATIVE-line transpiled schedule's duration-slope
reading AT THE SIM GATE × 1.25 (the 25% allowance is multiplicative HEADROOM: ×1.25,
never ×0.25; round-24 notation pin), floored at 1% absolute (the
band/gate input; the REPRESENTATIVE line itself is pinned, round-20: the line selected
by the SAME lexicographic rule applied to the calibration snapshot at freeze time,
recorded with the frozen artifact; it feeds τ_freeze, η_nom, and s_res_freeze alike); day-of reconciliation = the FLOWN line's slope × 1.25 must lie within
the frozen ±s_res; exceedance is a Class-1 abort (the transpile is local, pre-submit;
intent stated, round-25: algebraically this grants the flown line ZERO admissibility
headroom beyond the representative slope, DELIBERATE: the whole 25% lives in band
coverage, the band covers 25% more than the abort admits, and the best-edge flown line
is expected at or below the representative slope)
(review MC: a level+slope-swept b_curve ≈ 0.070 vs 0.052 level-only restores the
+1%-slope η(2)=0.98 device from 0.87 to 0.9999, and is SAFE on false-CONFIRMED; the
impostors die at clause (b)'s sign, and the null, which being Cauchy passes (b)'s
sign+magnitude ~27–31% of the time (model-dependent: review recomputes span ~26–32%,
rounds 23–26), is routed to VOID
by the Ŝ(j; 0) floor at ~1.2×10²σ:
widening (c) cannot create a false CONFIRMED through either path (round-19 precision);
the whole leak family vanishes
identically if the runner VERIFY finds r-flat/padded durations). Re-freeze trigger
PINNED (round-12 made the f129 rule local by DELETING the reference, restored
round 31, prior-work M3: the trigger is F129's, verbatim,
[IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md) "re-freeze trigger (bootstrap
SE > 1.3× analytic)"): a frozen band is re-derived only if its verification bootstrap returns
SE > 1.3× the fitted-model SE, or the goodness-of-fit gate fails post-freeze, or the
joint-power floor fails, or a runner-VERIFY finding invalidates a frozen INPUT
(round-22: e.g. non-monotone duration scaling re-deriving s_res: the fourth trigger the
"only" clause previously omitted); re-freezing exists only BEFORE the commit hash;
after commit, bands never change. Bands frozen per the ROLE-SPLIT basis rule above (the earlier "ONE basis for all bands"
was the round-15 leak's cause), per device (Kingston + Fez sets).

CONFIRMED iff all of (all clauses written in the SIGNED world; round-4 blocker fix:
the previous magnitude-typed clauses would have returned INCONCLUSIVE on a perfect
device; evaluation order PINNED (round-9): (a) → (b) → (c) with short-circuit on first
failure, so (c)'s σ̂ = ρ̂/|ρ̂(2)| is computed only after (b) has passed, which guarantees
|ρ̂(2)| > b_forgive > 0, no NaN path exists):
(a) **Blindness:** |ρ̂(1)| < b_blind, where b_blind = max over the budgets of
    |−sin(ε_watch + ε_ZZ)/cos(ε_ZZ)| (the exact signed center: ONE formula with the
    band procedure's, round 31 spec m1/physics m6; the |−sin ε_w| gloss this site
    carried differs by 1.5×10⁻⁶ and is retired; note the sign of the
    linearization) + the SIGNAL
    model's upper tail HALF-WIDTH about the model's own center, the whole construction
    taken as the UPPER ENVELOPE over the FULL coverage grid, the sign axis (round-24:
    the r = 1 tail is
    read at the ε_w-ALIGNED worst sign, never the (b)-envelope sign; that one is (a)'s
    easing sign; envelope ≈ 0.065 vs 0.042 single-sign in review MC: the sign-only
    figure; the full coverage envelope reads ≈ 0.080 at rehearsal, round 30 stats m9:
    arguments quoting 0.065 forward are superseded) AND every other member of THE
    COVERAGE GRID (the one definition at the band procedure; round 31, spec M1
    replaced this site's own axis list, which had drifted; round 29, stats m1 brought
    b_blind here; the derivation-basis
    clause UNIFIED round 30, spec BLOCKER: this site still said "the band's derivation
    point is the worst basis" one round after the band joined the coverage envelope;
    the envelope site in the band procedure is the ONE derivation basis, and the old
    worst-basis sentence is retired), pinned arithmetic
    (round-7): half-width = (p99.87 − model center), quantile-SE inflated; the coherent
    center enters ONCE (the signal model already injects ε_w, so adding its raw p99.87
    would double-count the center and loosen the ceiling: the one anti-conservative
    direction in the band set; algebraically the construction collapses to ONE number
    PER SOURCE BANK: each coverage-grid source's SE-shifted p99.87 quantile of |ρ̂(1)|
    under the signal model, the frozen band = the MAXIMUM over sources; it is stated as
    center + half-width only to pin the no-double-count arithmetic; the two-sources-
    one-number coincidence argument holds PER SOURCE, each bank injecting its own ε_w
    center, and stops being a single global statement on an envelope, round 30; at
    every source carrying the budget ε_w the injected center sits ~3–4σ off zero and
    folding is negligible;
    the folded-near-zero regime of a very good device only passes MORE easily under a
    ceiling, round-19): a magnitude ceiling
    built from the signed statistic's two-sided band; AND the parameter-free order check
    min(|ρ̂(0.75)|, |ρ̂(1.25)|) > |ρ̂(1)| (the floor is a local minimum in magnitude;
    its power is REPORTED per coverage-grid member like clause (b)'s, round 31 spec
    m3; no fallback is needed because the margin is structural, ~0.7·η against
    ~0.03, and the report measures rather than asserts that);
(b) **Forgiveness with the flip:** ρ̂(2) < −b_forgive (the signed even fixpoint:
    watcher-dressed magnitude η(2), flipped sign, ideal: full; b_forgive from the signal
    model at η_min, frozen from the
    r = 2 arm), AND |ρ̂(2)| > |ρ̂(1.75)| (neighbor recovery check into the fixpoint, NOT
    a law-enforced ordering: below ratio 0.924 a valid device fails the raw ordering,
    which is exactly what the fallback covers;
    the thinnest raw margin in the rule: the law's own gap is only ~8%, |cos(1.75·π/2)| =
    cos(π/8) ≈ 0.924 vs 1, and η eats it: the effective gap scales with η(2)/η(1.75),
    predicted < 1 above (η worst at r = 2), so the sim gate REPORTS this check's power under
    the SIGNAL model at frozen shots, computed at the Aer-MEASURED η_nom(2)/η_nom(1.75)
    (never at the ideal gap; two distinct thresholds, both stated: the raw check's POWER
    falls below the 99.5% trigger near ratio ≈ 0.94 at η(1.75) ≈ 0.9; the trigger ratio
    RISES as the absolute η falls, ≈ 0.958 at a STRESS-TEST η ≈ 0.74 (illustrative floor,
    round-11: NOT the abort-implied worst case; the abort thresholds RZZ EPC ≤ 0.5% +
    T2* ≥ 70 µs imply η ≈ 0.92–0.98), because the
    absolute gap is η(1.75)·(ratio − 0.924); the raw GAP vanishes entirely at cos(π/8) ≈
    0.924; the branch is routed on the MEASURED power, never on these illustrative
    numbers, the SELECTOR pinned, round 31 spec m5: the routing power is the
    MINIMUM over the coverage grid, the module's order_check_min_power, never a
    single-basis read); if < 99.5%, the PINNED fallback replaces the pairwise sample comparison:
    |ρ̂(1.75)| must lie inside its own signal-model two-sided band (p99.87, SE-inflated;
    round 29, stats M4 supersedes the round-27 single-worst-sign freeze: |ρ̂(1.75)| is
    an UNNORMALIZED magnitude whose center rides the LEVEL axis directly; the level
    shift over the coverage range, 0.074, is 1.7× the 0.043 sign shift round 27 graded
    MAJOR, and "worst sign" is undefined for a two-sided band, whose single-sign edges
    barely overlap; the band is therefore the UNION-OF-EDGES envelope: each edge takes
    its own extreme of the per-grid-point SE-inflated quantiles over the FULL
    coverage grid: level × slope × readout-sign plus the corners AND the coherent
    ε_w SIGN members (round 31, physics MAJOR: |ρ̂(1.75)| = cos(π/8 ∓ ε_w) is FIRST
    order in the watcher-angle sign; the raw-order gap swings ±0.0115 at budget and
    −ε_w is the EASING sign, the same class as the drift finding and the SECOND sign
    axis beside the readout one; ρ̂(2) is second-order flat in ε_w, so b_forgive is
    untouched; the ε_w-negative banks join the grid, so the fallback union, the
    order-power report, b_curve, and the LR reference all see the sign); round-27,
    the SIXTH sign-class instance, had found: ρ̂(1.75) is the
    same signed double ratio the round-24 b_blind fix protects, and the additive
    S-antisymmetric offset survives the double ratio here too (the watched arm's
    magnitude 0.85 and the r = 0 normalizer 0.92 weight it differently): |ρ̂(1.75)|
    shifts by ∓0.043 at δ = ±2% (η = 0.92 linearization), ODD in the sign: a band
    frozen at one sign is mis-centered by multiples of the shot SE at the other);
    the neighbor is where the law puts it, certified by band instead of by the noisy raw
    ordering (overlap owned, round-8: (c)'s SUP already bands r = 1.75 in normalized form;
    the fallback tests the UNNORMALIZED magnitude, differing by the η(2) factor, and its
    ~0.26% two-sided false-exclusion cost enters the assembled joint power, the cost
    re-measured on the SWEPT band at freeze); the
    fallback's own power is reported beside the raw check's, the branch taken is
    recorded at freeze, and the sim gate ASSERTS whether the fallback ever engages on
    the coverage grid (round-27 from below: never in coverage; the raw ordering held
    at ~5.5σ in the round-27 read, ~12σ in the round-29 from-below recompute; the pair
    is RE-MEASURED at freeze and the frozen value governs (round-29 physics m3 debt,
    recorded round 30); the < 99.5% trigger is reached only near
    η(2) ≲ 0.42, far below η_min; if it ever fires, the union-of-edges envelope above
    governs));
(c) **Shape:** the SUP statistic over the SEVEN INTERIOR arms r ∈ {0.25, …, 1.75}
    (ρ̂(0) ≡ 1 by construction, σ̂(0) = 1/|ρ̂(2)| carries no shape information, and
    r = 2 is the normalizer; both excluded) of
    |σ̂(r) − cos(rπ/2)·η_nom(r)/η_nom(2)| < b_curve, with σ̂(r) = ρ̂(r)/|ρ̂(2)|
    (magnitude-normalized so the target keeps the true sign; normalizer choice stated:
    |ρ̂(2)| injects η(2)'s fluctuation as a common factor, self-consistent because b_curve
    comes from the same signal model, at a modest power cost the projection absorbs);
    b_curve = the SUP's p99.87 quantile **under the SIGNAL model** (simultaneous; η_nom
    from Aer parity, pinned; MODEL-DEPENDENCE owned, round-6: the target's
    η_nom(r)/η_nom(2) is a simulator read; an η_nom error moves the target AWAY from a
    law-realizing device's true curve, so it can only GROW the SUP deviation, failing (c)
    toward INCONCLUSIVE, never easing CONFIRMED; fail-safe for CORRECTNESS, round-15
    sharpening: NOT fail-safe for POWER, since with angle-dependent η a worst-edge
    target systematically mismatches a better-than-worst valid device; hence the
    role-split band basis and the swept-range power gate, see the band procedure; note
    also, round-9: the r = 1 arm's coherent ε_ZZ offset is covered by the BAND, the
    signal model injects it at budget, not by the ideal-cos target, same fail-safe
    direction),
    AND
    (|ρ̂(2)| − |ρ̂(1)|)/SE_contrast > 5 (false-positive threshold; projected power
    ~95–150σ over the admitted η range at 16384 shots (review recomputes ranged ~84–170σ
    across models, rounds 13–24; the gate records the measured value; the earlier ~50
    low end was a lower-shot draft's), recorded as power;
    SE_contrast from the SIGNAL model, per the source architecture; a robustness ECHO of
    clause (b), not an independent discriminator, kept INSIDE the conjunction as a
    cheap backstop, its ~100σ margin making it effectively free, rounds 18–19;
    SE_contrast's sign-sensitivity is dominated by the swept end-to-end gate; both
    controls run both signs and the margin is non-binding, round-26).
**The ACCEPTED SET, stated (round 30, stats M5; made honest round 31, stats M2 +
physics m4):** the conjunction (a) ∧ (b) ∧ (c) accepts (i) the signed-cosine family
with the watcher-angle SCALE inside a TWO-SIDED, mildly asymmetric window: from
below ≈ [−3.9%, +3.75%] at b_curve = 0.070 and ≈ [−4.5%, +4.3%] at the rehearsal
0.0802; the frozen b_curve sets it (clause
(c)'s SUP is what binds the scale; a scale miss walks the interior arms off the
target curve before it moves the fixpoints; clause (a) alone would tolerate far more)
AND the r = 2 sign flip; and (ii), stated plainly, because the previous version's
"guard-separated impz class" wording overstated the separation, the ENTIRE
Z-DIAGONAL never-watching class at polarisation m BELOW the frozen z_k band
(window m ≲ 0.037 from below at rehearsal): such a device reproduces the signed
curve exactly with η ≡ 1 and reaches CONFIRMED with probability ≈ 1; the z_k guard
separates only the m-polarised branch (detection ≥ 0.999 asserted at m = 1, the
boundary escape at m = band measured and reported), and the ONLY discrimination
against the m ≈ 0 branch is the watcher-tomograph pair, which is DIAGNOSTIC-ONLY:
it gates the mechanism sentence, never the verdict. The negative-control ≤ 0.1%
headline is therefore a statement about the four named remote point alternatives
(null + three impostors) and the m-polarised impz branch, never about the
continuous scale neighborhood and never about the m ≈ 0 branch; the three
protections are different in kind and all stated.

Calibration check, not law-content: |R̂(2) − R̂(0)| < b_sym, with R̂(r) = D̂(j; r)/D̂(j′; r)
the unnormalized single ratio (defined here; b_sym from the same gate procedure;
INFORMATIONAL ONLY: its exceedance is recorded in the flight report and triggers no
verdict action; RELABELED round 29, stats n1: algebraically R̂(2) − R̂(0) =
(η(2) − 1)·a/a′, so the statistic reads the η(2) LEVEL, not j-vs-j′ symmetry; it
duplicates the decoding table's η-anomaly entry and is kept under its honest name,
"the unnormalized η(2) level check").
**DEVICE-DEVIATION:** trigger = a frozen likelihood-ratio statistic: LR(impostor vs law)
> the pinned threshold for any named impostor. **The threshold is pinned like every band
(round-6):** reference distribution = the LR values under the SIGNAL model; threshold =
the fitted p99.87 quantile inflated by its own fitted-quantile SE, tightened FAMILY-WISE
across the three impostor tests (per-test quantile raised until the MEASURED family
exceedance on signal-model runs is ≤ 0.13%; a spurious DEVICE-DEVIATION on a
law-realizing device is the error this bounds; the SE-inflation side is the
fewer-false-deviations side, consistent with the pinned direction rule). Named
impostors: (i) monotone erasure, PINNED to a curve (rounds 8–10; an LR needs a
functional form): ρ_i(r) = cos²(r·π/4), a PER-HALF-GATE-dephased watcher: the same two half-angle
gates with EACH gate's sign independently randomized; the 4-way average
¼·[cos(rπ/2) + 1 + 1 + cos(rπ/2)] = cos²(r·π/4), exact. Two near-misses pinned as
implementer WARNINGS (rounds 9–10, both reproduce the law's own curve and would give 0σ
separation): a 50/50 mixture of the two FULL watcher directions (cos is even), and
measure-and-forget on k (tracing or measuring k leaves j's reduced state unchanged);
the incoherence must live PER HALF-GATE. Starts at 1, decays monotonically, 0 at r = 2,
no recovery, no flip, parameter-free. (ii) half-period cos(rπ/4), (iii) the SIGN-FLAT magnitude curve +|cos(rπ/2)|
(matches every magnitude, cannot produce the r > 1 sign flip; separated by the signed arms
alone), each with ≥ 10σ separation IN THE PINNED LR STATISTIC required before freeze
(round-21: the metric named); recorded as "the
device does not realize the angle law." The impostor set is a LABELING AID, not an
exhaustive alternative space (round-12): an unlisted non-law curve fails CONFIRMED
(clause (b) dies without the flip) and falls to INCONCLUSIVE if no LR trips, never to a
false CONFIRMED. The LR statistic is computed on ρ̂ DIRECTLY (no
σ̂ normalizer), so a monotone-erasure device with ρ̂(2) → 0 reaches DEVICE-DEVIATION
without a NaN routing to VOID. **LR construction PINNED (rounds 16–17):** the Gaussian
log-likelihood ratio on the 8-ARM ρ̂ vector r ∈ {0.25, …, 2.0} (round-17 fix of the
round-16 pin: ρ̂(0) ≡ 1 identically, zero variance, so a 9-arm covariance is SINGULAR
and Σ⁻¹ does not exist as previously written; the 8 non-degenerate arms carry the
statistic), both hypotheses as POINT curves DRESSED WITH THE SAME η_nom(r) profile
(μ_law(r) = cos(rπ/2)·η_nom(r), μ_i(r) = impostor_i(r)·η_nom(r), with η_nom at the
MID-RANGE basis: the same basis as clause (c)'s target, round-18; the LR separates
SHAPE and SIGN, never the common attenuation), sharing ONE covariance Σ: the 8×8 empirical
covariance of ρ̂ from the held-out signal-model seed bank (the arms correlate through
the shared r = 0 normalizer; a diagonal Σ would misweight):
LR_i = ½[(ρ̂−μ_law)ᵀΣ⁻¹(ρ̂−μ_law) − (ρ̂−μ_i)ᵀΣ⁻¹(ρ̂−μ_i)]. The LR threshold is
CORNER-CONDITIONAL over the FULL coverage grid: the level × slope × sign points plus
the round-30 fixed-j′/box/η-top corners (round 29, stats M1: the
slope axis was the one axis the LR reference never swept, on the statistic whose whole
job is shape: a law-realising device at the admitted −1% slope corner read family
DEVICE-DEVIATION at ~1% against the 0.13% pin, WORSENING with shots: the b_curve
shot-independent-offset signature; pooling the slope into the reference does NOT repair
it, ~0.56% residual at the corner; only b_curve's max-over-grid construction does):
per grid point the p99.87 + SE of the LR values under the signal model, the threshold =
the grid MAXIMUM per impostor (history of the axes: round-17 added the η(2) sweep;
round-26, the FIFTH sign-class instance, added the sign axis after the ≤ 0.13% family
exceedance had been certified at the hardening sign only: swept, not assumed; round 29
completes the grid with the slope axis and replaces pooling by the corner-conditional
maximum).
**Catch-all:** anything matching neither CONFIRMED nor DEVICE-DEVIATION nor VOID →
INCONCLUSIVE, never a physics verdict.

Decoding table (diagnostic only, never upgrades): watcher-angle residual → measured ζ_jk +
curve-zero fit; η anomaly → unnormalized D̂(j) r = 2 vs 0; generic points off → Z-leak +
transverse guards; zeros displaced → transpiled-angle record; Z-polarised watcher
(the impz class, round 30) → |⟨Z_k⟩| against its capped band (the m-polarised branch)
+ the watcher-tomograph ⟨Y_k Z_j⟩ (+1 superposed, 0 Z-diagonal: the m ≈ 0 branch).

## Sim gate specification (pinned)

This gate is the house **7a** stage and the counts-level gate below is the house
**7b**, the two-gate standard of the prior flights (round 31, prior-work M3:
thirty rounds had re-derived the architecture anonymously after round 12 deleted
the one cross-reference instead of resolving it; the citations restored:
[IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md) §"7a"/§"7b", "7b (counts
level): the actual runner's circuits ... counts through the actual estimator code;
hard pre-flight abort if 7b contradicts 7a", and
[IBM_CONCENTRATOR_RELOADED](IBM_CONCENTRATOR_RELOADED.md), whose
representative-background band pin is the precedent for this document's
representative-line / τ_freeze-vs-τ_dayof pattern; the price-pair reviewer list is
the origin of the counts-level-null / Aer-parity / instrument-failure-path
requirements re-derived here).

Full estimator end-to-end: asymmetric confusion incl. CAL finite-shot noise → mitigation
inversion → distribution-level quasi-probability conditioning (as the pipeline) →
stochastic branch statistics → per-witness axis fits → signed Ŝ, guards, ρ̂, σ̂.
Injections: coherent always-on ZZ (full pre-measurement window, every NN bond; NN-only
injection SUFFICES; non-NN/spectator ZZ phases are r-independent up to the
watcher-layer duration differences already budgeted in η_nom/s_res (with angle-scaled
durations the accumulated spectator phase varies slightly among r > 0 arms too; the
same budget covers it, round-21), so they divide out in the double ratio
ρ̂(r) = [·]_r/[·]_0; round-6, verified) +
watcher-angle offset 2·ε_cal + **per-qubit coherent RZ frame offsets, static AND
arm-time-dependent drift** (the axis is fit at r = 0 and applied to later arms; the sign
estimator must be certified against inter-arm frame drift, not only statics;
drift budget PINNED, rounds 16–17: ≤ 0.05 rad PER-ARM increment, the guard covering the
ACCUMULATED tilt at each arm's interleave slot, the power cliff sits near a 0.2 rad
per-arm rate, clause (c) breaking first, and the drift IS day-of-anchored through the
transverse guard: sin(φ_r)·|record_r| leaks into T̂(w; r), the operative leak the SUP's
worst arm, the worked example RESTATED round 30, stats M4 (the round-22 example
multiplied the wrong arm by the wrong rate, two errors compensating to the same 0.18):
sup over arms of sin(rate·d_r)·|dressed record_r| with d_r = the arm's slot distance
from the r = 0 anchor; under the pinned order the maximum distance d = 4 carries
r = 1.0 (record ≈ 0, no leak) and r = 0.25 (dressed record ≈ 0.9), so at the 0.05
budget the worst accumulated tilt is 4 × 0.05 = 0.2 rad at r = 0.25 and the leak is
≈ 0.9·sin(0.2) ≈ 0.18; the trip-point algebra implied by the band is AMBIGUOUS by ~2×
(0.069 vs 0.137 rad/arm depending on which arm trips first), so the pinned criterion
is the MEASURED probe, never the algebra: the T̂ bank must trip with P ≥ 0.99 at the
0.2 rad/arm cliff rate (round 30, spec m5 pins the number); an UNTRIPPED 0.10 rad/arm
device carries a clause-(c) systematic of ~0.056 against b_curve ≈ 0.070: inside the
band, stated, and the INTERVAL ownership stated with it (round 31, spec m9: between
the 0.05 budget the band covers and the 0.2 cliff the probe certifies, coverage is
the ±s_res axis plus the swept power gate, and the 0.10 spot value is one point of
that composition, not an envelope; an intermediate-rate device that neither trips
nor stays inside the (c) systematics lands INCONCLUSIVE, the fail-safe direction);
so the frozen T̂ band is REQUIRED to cover the 0.05 budget while tripping beyond it;
fail-safe both ways) +
incoherent T1/T2/depolarizing at the worst-admitted basis with j-vs-j′ asymmetry drawn from
the calibration spread (10⁵ seeds PER GRID POINT of the 2-D sweep, never one bank
partitioned across the grid, round-19); readout-correlation variant injected in BOTH dangerous modes
(round-8 + round-9: the double ratio cancels only the symmetric multiplicative part,
MC bias/SE ≤ 0.03 for that mode, so the injections are (1) the r-INDEPENDENT ASYMMETRIC
additive mode, a fixed S-outcome-dependent witness assignment offset, which does NOT
cancel and can bias ρ̂(2) by ~2.75 SE per 1% differential (round-27 basis pin; ~5.5 SE
per 2% at the pinned shots) with one sign easing clause
(b), and (2) the ARM-DEPENDENT mode: the witness readout populations differ per arm,
~50/50 at r = 1 vs ~96/4 at r = 2, making state-dependent/T1-during-readout error
r-dependent, symmetric with the RZ-frame arm-drift treatment; the injection enters
BOTH measured channels, variant-A/X and variant-B/Y; the offset is a property of the
physical Z-readout, present after either pre-rotation; round-25: a Y-only injection
would freeze the transverse band blind at every sign); the variant widens bands
AND its residual bias on ρ̂(1)/ρ̂(2) is reported under the estimator section's
PINNED-SIGN rule (round-12: hardening sign in the signal-model band runs, easing sign in
the negative controls, NO center surgery; the tensor-product CAL structurally cannot
measure this correlation; the in-job crosstalk pair measures the S-correlated assignment
differential on hardware day-of, a Class-2 in-job VOID gate at > 2%, recorded for
provenance, and the
injection covers the model side). Aer parity on the ACTUAL
transpiled circuits asserts: variant B recovers Y; both writes exactly RZZ(π/2); 2q counts
equal budget; η_nom(r) read here (the profile ARITHMETIC pinned, round 28: η_nom(r) =
|ρ̂_Aer(r)/cos(r·π/2)| with the r = 1 point interpolated from its neighbors; the
RULE pinned (round 31, spec M7: the "0/0 is harmless" justification covers only
clause (c), but η_nom(1) ALSO dresses the LR hypotheses, where impostor_i(1) is 0.5
and 0.707, far from zero): **η_nom(1) = ½·(η_nom(0.75) + η_nom(1.25))**, the linear
neighbor mean in η, exactly what the module computes, where
cos = 0 makes the division 0/0 while the (c)-target term cos·η_nom is 0 regardless;
NEVER the raw |ρ̂_Aer| profile: the first freeze rehearsal stored |cos·η| as η_nom and
bent the clause-(c) target into a cos·|cos| curve, corner power 0.10, caught by the
sign-and-level-swept power gate before anything froze). **POSITIVE CONTROL (round-4, mandatory; round-6 two-branch split:
the ideal signal has η ≡ 1 and may legitimately sit outside bands frozen at the noisy
budget, so the two branches check different things):** (branch 1, flight-realistic
FEASIBILITY SMOKE TEST) the Aer-noisy signal is run once through the ACTUAL frozen
verdict code against the FROZEN bands; a non-CONFIRMED outcome halts with a
NON-DISCRETIONARY response (round-12: re-run once with the next pinned seed; two
consecutive non-CONFIRMED → treated as a joint-power failure → the re-freeze path,
recorded; a single draw at ~96–98% joint power fails a few percent of the time by design;
the statistical gate is the JOINT POWER below, not this draw); branch 1 additionally
carries a BANK comparison (round 30, stats M7: a single draw cannot certify the
model's correlation structure, and an 8×8 ρ̂-correlation misspecification is worth up
to ~14% of b_curve through the LR covariance and the SUP): an Aer seed bank's
first two ρ̂ moments per arm and its 8×8 correlation matrix are compared to the dm
model's bank, each entry within 5 combined finite-bank SEs (a loud-only screen: the
tolerance is the two banks' own sampling law, never a hand-set number; exceedance →
the model's injection bank is re-examined before freeze); (branch 2,
code-consistency) the ideal signal against η ≡ 1 targets
(the same verdict code with η_nom set to 1; the OTHER constants pinned round 31,
spec m6: branch 2 runs the WIRING-level band set: the module's provisional bands,
never the frozen ones, since an ideal device against noise-frozen bands can
legitimately miss clause (c); a check of clause logic and signs, not of the
bands) must return CONFIRMED. Each of the three named impostors must return
DEVICE-DEVIATION at the pinned LR threshold, asserted by the gate, so any
sign/normalization regression in the verdict rules is self-catching before flight. The gate also reports the AGGREGATE false-VOID probability of
the whole guard bank on valid runs (the 8 guard tests, 5 pooled + 3 floors, at tail
quantiles compound); pinned
ceiling 0.5% at the p99.97/p0.03 default (round-21 re-derivation; the old 2% was the
p0.13-era figure), exceedance escalating directly to the remedy ladder with recorded
provenance; **the aggregate is MEASURED AT EVERY POINT of the coverage grid and the
corner maximum GATED at the same pinned 0.5% ceiling (round 29, spec m5 measured it;
round 30, spec M4 makes the corner maximum a GATE; a reported-only corner could
exceed the ceiling exactly where round 28's blindness lived); the held-out
FIXED-j′-CORNER twins run through the gates too (round 30, stats M2: power and
false-VOID measured on the fixed corners, not only the draw mixture; this is also
the near-degenerate-valid certification of the floor construction, stats m13); the
gate additionally measures guard DETECTION power: a probe bank at the
0.2 rad/arm drift cliff must trip the T̂ bank with P ≥ 0.99 (the numeric criterion
pinned round 30, spec m5; round 29, spec
m6 + stats m2: the T̂ band's pinned duty is two-sided, cover the budget AND trip
beyond it, and no gate previously measured the second half), and the impz DETECTION
and BOUNDARY measurements plus the watcher-tomograph separation (signal claim rate
vs the Z-diagonal class's claim rate ≤ 0.1%) per the negative-control section; the
SE(ρ̂(1))/SE(ρ̂(2)) ratio is REPORTED with its provenance isolated (round 30, physics
MAJOR: the gate reports the ratio under BOTH provenances, the full chain AND a
no-readout-layer isolation bank; measured reads have spanned 0.88 (full chain,
round-29) through 0.94/0.94 (both provenances, the v32 validate) to the reviewer's
own shot-noise model's 1.20–1.26; the ratio is configuration-dependent, not
structural in either direction, and no assert rides it).** The gate
further reports: the (b)-order-check power at frozen shots (fallback rule in clause (b))
and the family exceedance of the LR thresholds on signal-model runs (≤ 0.13%).
**JOINT POWER, assembled and pinned (round-7 stats):** the negative control is a measured
probability, so the positive side must be too; branch 1's single draw proves feasibility,
never power. The gate runs the FROZEN verdict code over the full signal-model seed bank at EVERY
member of THE COVERAGE GRID (round 31, stats M1: this site's old "[η_min, 1] × sign"
wording described a strictly smaller space than the band grid; the box corner, the
design's hottest VOID point, was power-blind, and the η-top member sits at
η(2) ≈ 1.002, literally OUTSIDE an interval whose top is 1; the swept set and the
band-freezing set are now the SAME set by construction, and the level range's true
top is the η-top member, above 1) (round-23: the sign axis was the second
gate-blindness; no gate ever ran a valid device at the easing sign, where the
hardening-frozen floor terminally VOIDed it; round-15: a single-basis
measurement is blind to the clause-(c) device-vs-target shape mismatch; the gate must
exercise it; the floor must hold at EVERY sweep point of every axis, reported as the
minimum) and
reports P(CONFIRMED | valid device) as ONE number over the coverage grid
(round-22 label fix: not "the admitted range"; the sweep overshoots to the ideal,
and since round 31 past it, to the η-top member): the
END-TO-END measured CONFIRMED
fraction through the frozen code, never a product formula (round-21: the gloss
"clause power × (1 − false-VOID)" is NOT exact; the Ŝ(j; 0) floor and clause (b) are
negatively correlated through the shared normalizer, corr ≈ −0.36, biasing the product
UPWARD; the gloss is illustration only); the gate records the corner
MINIMUM over the sweep explicitly, since that minimum is what the floor gates; pinned
floor ≥ 95%: below → re-freeze, recorded, with the LEVER NAMED PER FAILURE MODE
(round-18): interior-margin shortfalls respond to shots; corner-limited shortfalls (the
(b) 3σ quantile margin, the (c) systematics) respond only to DESIGN: narrow the
admitted range or widen the band per its pinned construction; shots cannot move a
shot-independent offset; a DESIGN remedy that narrows the admitted range means TIGHTER
day-of aborts (stricter than EPC 0.5% / T2* 70 µs), which RAISES the day-of abort
probability; that operational cost is accepted and recorded, never hidden inside
"re-freeze" (round-19); the gate
additionally perturbs the η(r) SHAPE (an independent duration-slope offset, not only the
level sweep, round-16: the level sweep alone cannot span a profile-slope error; a
non-exponential profile at the worst level measured ~93% in review MC while an
Aer-referenced gate would read ~99%) and reports the min joint power under the shape
perturbation; the perturbation range IS ±s_res, the same 2-D grid b_curve is swept over
(round-18: the gate never stresses a space the band does not cover), the slope itself
bounded by the runner-VERIFY duration-vs-angle reading; the joint power, the false-VOID aggregate, and
the negative controls are measured on a seed bank HELD OUT from band fitting (round-14:
in-sample scoring is optimistically biased; measured ~0.04 pp at 5·10⁴ seeds,
negligible at the 10⁵ default, pinned anyway; the floor and its per-mode re-freeze lever
are pinned above). Reported beside it: the joint power's η(2) sensitivity (an optimistic
η(2) freeze is what collapses it; freezing the magnitude bands at the worst-admitted
basis is the protection, so the sensitivity is the number that proves the protection
held).
**The full gate, band freezing, BOTH controls (positive two-branch + negative
null/impostor), the false-VOID measurement, runs PER DEVICE band-set: Kingston AND Fez,
both gate records committed with their bands; the flown device's record governs. If one
device fails its own gate, single-submit the passing device (recorded). The sim-gate
seeds are recorded in the committed pre-reg (the verdict is RNG-free, but the frozen
band VALUES must be reproducible).**

## Abort gates (day-of, hard, no override, TWO CLASSES, round-15: a Batch has no
mid-batch conditional logic, so the in-job thresholds cannot prevent the spend; the
class NUMBERS name the evaluation position: 1 = pre-submit, prevents the spend;
2 = in-job, protects the verdict, never a severity ranking, round 30)

**Class 1: PRE-FLIGHT aborts (evaluated before submit; they PREVENT the spend).** All
inputs exist pre-submit: the fresh calibration values, the local transpile, the runner
state. Fresh same-day calibration; line selection PINNED (round-8; "house score" was undefined):
enumerate candidate 4-qubit lines from the fresh calibration, drop any failing the abort
thresholds below, score the survivors lexicographically: (1) max write-edge RZZ EPC
(minimize), (2) watcher-edge EPC (minimize), (3) S readout assignment error (minimize),
(4) min T2* over the four (MAXIMIZE; direction pinned round 29: read uniformly as
minimize, item 4 would select the worst line), (5) max WITNESS readout assignment
error (minimize; added round 31, stats B1's true kernel: witness readout was the
one abort-bounded nuisance the score never ranked; coverage was already complete
without it, the freeze corner sits above the abort ceiling, but zero selection
pressure on a ranked axis was a free give-away), with BOTH role assignments of each physical line scored (k–j–S–j′ vs its mirror;
the writes must sit on the two best edges); the runner records the full score table and
the chosen assignment. Flown-line duration-slope × 1.25 within the frozen ±s_res
(round-19: the s_res day-of reconciliation, see the band procedure). Abort-to-η mapping
one-liner (round-19; NOMINAL/approximate: the EPC-to-attenuation map is
channel-dependent, and the Aer-measured sweep is what governs): two watcher gates at
EPC ≤ 0.5% ≈ 1–2% amplitude loss ⇒ η(2) ≳ 0.96 nominal (idle-term sign corrected
round 31, physics m2: the old line ADDED the 140 ns idle's 0.2–0.4% as an η loss,
but that idle sits on j′ and S; j′ is in η's DENOMINATOR, so its idle RAISES η,
the η-top member's own mechanism, and S carries no η weight; the 0.96 figure keeps
its padding and is NOMINAL anyway); the 0.92 lower edge adds headroom for
INCOHERENT modeling error only (round-25 propagation of the round-24 incoherent-only
pin; this was its last surviving stale site; coherent residuals live solely in the
ε_w ledger, the ±s_res axis, and the sign-swept gate; the sweep's corner minimum makes
the exact value non-load-bearing). S readout ≤ 2% both
directions, all ≤ 3%; the native fractional-rzz path taken (CZ fallback = NO-FLIGHT, no
doubled-budget branch); median RZZ error (at π/2 reference) ≤ 0.5%; T2* ≥ 70 µs all
four; transpiled 2q counts exact; τ_dayof read from the flown line's schedule (local
transpile, pre-submit); PER-SHOT MEMORY enabled and verified on every PUB (memory=True,
round-18: the split-sample verdict machinery is unexecutable on bare counts; every-PUB
is deliberately BROADER than the need: only the two r = 0 arms consume shot indices,
a uniform conservative setting, round-27); submit
record + RAW PER-SHOT MEMORY persistence ARMED (counts derivable from it; the split
needs the memory; the counts-level GATE itself ran at runner stage, pre-commit;
failures loud).

**Class 2: IN-JOB VOID GATES (computed from the Batch's own characterization PUBs,
evaluated BEFORE the verdict; the spend is already made, the VERDICT is protected; any
trip → the flight is VOID, no verdict claimed, recorded).** Measured π·Δf_jk·τ_dayof
≤ 0.01 rad (Δf from the in-job ζ conditional-Ramsey pair, per the ζ book); watcher-angle
amplification |8·ε_cal| ≤ 0.08 rad (the N = 8 repeated-RZZ(π/2) PUB pair, Ramsey-read on
j: the coherent-angle guard the incoherent EPC metric cannot provide); S↔witness
crosstalk differential ≤ 2% (the in-job crosstalk pair; recorded for provenance; the
model side is covered by the PINNED-SIGN injection at the 2% ceiling (round-12; no
center surgery), so the day-of value never moves a band; rounds 9–12: the one readout
mode the double ratio does not cancel has a hardware bound).

## Cost (by measurement)

Shots-ratio vs f129 (1.05M shots ≈ 297 billed s): ~36–38 PUBs × 16384 (with the two
r = 0 arms at 32768, i.e. +16384 EXTRA each, not additive circuits) ≈ 620–655k shots → linear anchor ~176–185 s; with per-PUB overhead margin
~180–225 s ≈ **~3.0–3.75 QPU minutes** (round-8: shots doubled for power margin; round-9:
crosstalk pair + 2× r = 0; round 30: + the watcher-tomograph pair, ~10 s; the linear
anchor is the tight figure, the range carries
overhead). Account state reconciled at freeze via the EXTERNAL pipeline's
`_query_billing.py` (ibm_quantum_tomography, outside this repo, not a repo-file
reference; measured, not the stale May figure); the flight is ~2% of the
remaining annual gift; worst case if the cancel loses the race (both twin jobs
complete and bill): ~2× ≈ 6.0–7.5 min ≈ 4% of the gift, accepted and recorded
(round-17).

## Runner-stage VERIFY list (pinned; each check loud, each recorded)

- **Δf is extracted as the conditional fringe-frequency DIFFERENCE**: both conditional
  branches fit separately, the difference taken; the entire ζ factor-2 book
  (θ_static = π·Δf·τ) hinges on Δf being a difference, not a single-branch frequency.
- **The amplification statistic carries the ½ factor**: the abort input is the
  HALF-difference ½·[⟨Y_j⟩_{k=|0⟩} − ⟨Y_j⟩_{k=|1⟩}] = sin(8·ε_cal), never the raw
  branch difference (which is 2× it; round 31).
- **Fractional-rzz duration-vs-angle scaling** read from the transpiled schedule (decides
  whether η is worst at r = 2 or ~r-flat; either way b_forgive frozen at r = 2 stays
  conservative, recorded, not assumed).
- **θ = π/2 is INCLUSIVE** in the device's fractional-rzz range (both writes sit exactly
  there).
- **The native fractional-rzz path is actually taken** (flag set + verified; no in-house
  runner has set it before; a CZ fallback is NO-FLIGHT; no CZ bands exist).
- **Variant B pre-rotation order** (Sdg THEN H) and transpiled 2q gate counts vs budget.
- **Watcher-tomograph pre-rotation order on k** (Sdg THEN H, the same convention) and
  the raw-Z measurement of j in that pair (round 30, the ⟨Y_k Z_j⟩ sign convention is
  load-bearing: prediction +1, not −1).
- **τ read as ONE number** from the transpiled schedule (full H-to-pre-rotation window).
- **Per-shot memory ON for every PUB** (stable delivered order): the even/odd split is
  undefined on aggregated counts (round-18); persistence BYPASSES `get_counts()`
  aggregation, the Concentrator's burned trap, named (round 31, prior-work:
  [IBM_CONCENTRATOR_RELOADED](IBM_CONCENTRATOR_RELOADED.md) persisted pooled
  per-cell histograms because SamplerV2's `get_counts()` aggregated at save time,
  silently no-op-ing a pre-registered analysis leg; the raw per-shot arrays are
  what this design persists). Flag names verified against the runtime
  version at runner stage (round-19: on SamplerV2 per-shot data is the BitArray default
  and resilience_level is not a Sampler knob; memory=True is the V1 idiom; the pinned
  REQUIREMENT is per-shot data + no runtime mitigation, whatever the current API calls
  it).
- **The two r = 2 watcher gates remain two DISTINCT fractional RZZ(π/2) pulses**, not
  merged into RZZ(π) or Clifford-substituted (a count-preserving Clifford pass could
  alter the r = 2 noise profile without tripping the gate-count abort; round-18).
- **The MINIMUM half-gate angle RZZ(π/16) (r = 0.25)** is inside the native fractional
  range (the maximum, θ = π/2, has its own line above).

## Gate order (pinned)

design rounds converge → sim gate per spec → bands frozen (per device, per the
ROLE-SPLIT basis rule of the band procedure: worst-admitted for magnitude bands,
coverage-range envelopes for the upper guards and ceilings, mid-range for the (c)
target; round 29 fixed this line, which still carried the one-basis wording the
round-15 note names as the leak's cause, provenance recorded) → EMPTY ROUND on the frozen numbers themselves
(numbers-against-rules: fresh reviewers check every frozen constant against this
document's pinned rules, round-12; the design rounds covered the rules, the committed
artifact's numbers need their own outside cut) → runner + counts-level gate (the house 7b,
[IBM_F129_RAMSEY_FRINGE](IBM_F129_RAMSEY_FRINGE.md) §7b, DEFINED for this flight
round 29: the runner's analysis path fed with synthetic counts/per-shot memory from
the sim-gate models, same seeds, THROUGH the imported verdict module, must reproduce
the sim gate's estimator outputs bit-for-bit before any hardware data exists) → empty
rounds on runner + records
→ pre-reg committed (hash in banner; the commit carries or follows the signed-coherence
corollary edits in PROOF_RECORD_PARITY_LAW + the F135 registry clause; the load-bearing
dependency must be in the repo before the hash points at it, DISCHARGED (round 30,
spec m7): the corollary and registry edits are committed, base df63957 + the v29 bank
4d55b14; the freeze commit additionally carries the F135 registry entry's inbound
link to this file (round 30, spec n5: the pre-reg had no inbound link from the
registry); **CODE-IDENTITY PIN, round
29, stats M2: the runner IMPORTS the promoted sim-gate module's verdict function and
estimator; the measured joint power and negative-control rates are properties of that
code object, and a re-implementation would void their transfer; the promoted module
is `simulations/record_parity_gates.py` (round 31) and the freeze commit hash binds
it; a runner-VERIFY line asserts the import**)
→ fresh calibration + hard aborts → **Tom's explicit
go** (the GO POSITION pinned, round 30, spec m8: Tom's 2026-08-04 08:12 advance go,
"Du hast das Go für alles", quoted in full at the HARDWARE RECORD, DISCHARGES this
gate-order go PROVIDED every pinned gate before the submit is green; a red gate
ABORTS and reports, and the advance go never overrides an abort) → ONE ANALYZED Batch job (twin-submit: two submitted, at most one analyzed,
normally whichever queue empties, the earliest-submit tie-break governing only IF BOTH
COMPLETE, round-23 scoping; cancel loser normally unbilled) → IN-JOB VOID GATES evaluated
(the three Class-2 thresholds from the Batch's own characterization PUBs; any trip →
VOID, no verdict claimed, round-15: these cannot precede the spend, a Batch has no
mid-batch conditionals) → HARDWARE RECORD → post-flight
empty rounds against the COMMITTED bands (the committed rule governs even against the
runner's own printout) → on a CONFIRMED (or DEVICE-DEVIATION) record, the
Confirmations entry lands in BOTH registries, `fw.Confirmations` (Python) and
`ConfirmationsRegistry` (C#, count-asserted), per cockpit rule 3 (round 31,
prior-work m9: the registry is where future sessions look it up instead of
re-deriving).

**Terminal rule (rounds 16–17; the campaign-level false-CONFIRMED door closes here):**
this is a SINGLE TERMINAL flight. A VOID from a re-fly-ELIGIBLE cause (enumerated,
round-17: NaN/failed-fit/mitigation-pipeline abort, persistence failure, or a Class-2
in-job gate trip: an out-of-spec device state; the STATISTICAL guard VOIDs, factor
floors, transverse guards, Z-leak, ⟨Z⟩ sanity, are NOT re-fly-eligible and are
terminal) may be re-flown ONCE, behind a fresh explicit go from Tom, recorded; a second
VOID is terminal. INCONCLUSIVE is TERMINAL; the record stands, and any new attempt is a
NEW pre-registration cycle with its own design rounds and its own commit. The re-fly
inherits the SAME machinery unchanged (round-23: same committed bands, same twin-submit
+ tie-break, same aborts: no post-data discretion over the re-fly's device). CONFIRMED
and DEVICE-DEVIATION are terminal by nature. The campaign-level false-CONFIRMED rate is
thereby bounded at ~0.2% campaign-wide (two near-mutually-exclusive flights at ≤ 0.1%
each; round-18 arithmetic; SCOPED round 31, spec M9, exactly as the accepted-set
paragraph scopes the per-flight bound: this is a bound over the named point
alternatives and the m-polarised impz branch, never over the continuous scale
neighborhood or the m ≈ 0 Z-diagonal branch, which the accepted set states), never
"fly until it confirms."

## Revision notes

- **v33 (2026-08-04, round 31, FOUR lenses for the first time: physics / spec /
  stats plus the NEW prior-work/repo-consistency lens (Tom's call: "das Experiment
  erfindet in Teilen das Rad neu", confirmed, the lens is a standing member from
  here); NOT clean: spec 2B/9M, stats 1B/2M, physics 0B/2M (core re-confirmed to
  machine precision), prior-work 1B/4M; every finding recomputed from below, TWO
  reviewer numbers refuted in the process):** prior-work BLOCKER = the ζ-book
  collision reconciled BY NAME (price-pair's measured ζ is the conditional
  difference ⟹ H = (πζ/2)·ZZ, this book; F129's H_ZZ = πζ·ZZ pin doubled it; the
  F129 v2.7 correction, the gate-script comment, and the ζ²-proof convention note
  ride this fold; F129's flown verdict re-read and it stands, +0.0114 vs the
  corrected 0.0192 window); spec B1 = the negative controls now sweep the CONFIG
  axis {worst, ideal, fixed-j′, box, η-top-null} so the ≤ 0.1% UCB certifies the
  band set actually committed; spec B2 + physics m5 = the crosstalk Δ reduction
  pinned as the FOUR-circuit combination (the pair + CAL0/CAL1; the pair alone
  cannot identify Δ); stats B1 = REFUTED in its number (the box corner's readout
  base is the abort ceilings, not the backend median; coverage exceeds the
  admitted set), its true kernel kept: witness readout joins the lexicographic
  line score as item 5; physics M1 = the coherent ε_w SIGN is the second sign
  axis (|ρ̂(1.75)| = cos(π/8 ∓ ε_w) first-order, −ε_w eases (b) by +0.0114;
  ε_w-negative banks join the coverage grid); physics M2 = the amplification
  statistic pinned in one equation, ½·[⟨Y_j⟩_{k=|0⟩} − ⟨Y_j⟩_{k=|1⟩}] =
  sin(8ε_cal) (the raw branch difference is 2× and would abort at half budget);
  stats M1 = the power sweep IS the coverage grid now (box was power-blind, η-top
  sat outside "[η_min, 1]"); stats M2 = the accepted-set sentence made honest
  (the m < band Z-diagonal branch CONFIRMS; the tomograph is the only
  discrimination and gates no verdict); spec M1 = THE COVERAGE GRID defined once,
  six drifted per-site axis lists superseded by reference; spec M2 = b_forgive's
  stale second construction retired; spec M3 = the z_k ceiling got derivation /
  freeze-gate semantics / named envelope-lowering remedies; spec M5 = η_min
  pinned as 1 − 2·(1 − η_nominal-worst); spec M6 = T̂'s SUP over ALL NINE arms;
  spec M7 = η_nom(1) = the linear neighbor mean, stated where the LR consumes it;
  spec M8 = the round-30 "self-cancelling sentence" LOCATED and discharged (the
  never-ease rule scoped to CONFIRMED-bearing bands; the s_j0 floor's opposite
  pull resolved by the above-null check as a freeze-halting guard); spec M9 = the
  campaign ~0.2% bound scoped like the accepted set; prior-work M2 = F136 cited
  (the flown line is its pointer branch: the sign law as a second gated theorem;
  gauge baselines named; F136(iv) zeroes the S-idle dephasing the η bundle
  budgeted); prior-work M3 = the 7a/7b house architecture citations restored
  (F129/Concentrator; the re-freeze trigger is F129's verbatim); prior-work M4 =
  the module PROMOTED (`simulations/record_parity_gates.py`, F129's
  committed-gate precedent; enters the repo before the binding freeze);
  prior-work M5 = twin-submit reconciled with F129 §9 (one ANALYZED job, scope
  sentence kept; Kingston always submitted first, spec m8); minors: branch-1
  bands = the wiring set (m6), fallback selector = the grid-minimum power (m5),
  (a)-order-check power reported (m3), GoF-vs-seed-cost reconciled (m2), interior
  control levels owned (m4), drift-interval composition owned (m9), rehearsal
  banner scales qualified (stats m1), the frozen-injection sentence restated
  (stats m5), (a)∧(b)-passing label (stats m6), the 0.235 decomposition flagged
  for the freeze re-read (stats m4), seed floor 2996 (stats n2), the max-bias
  never-ease exception acknowledged (stats n3), tomograph = second-order ROBUST
  (physics m1), the idle sign in the abort-to-η line corrected (physics m2),
  "neither fixpoint first" repaired + the (a)-easing at d = 4 stated (physics
  m3), the accepted window two-sided (physics m4), b_blind's center = ONE exact
  formula at both sites (m1/m6), the write-miss cancellation recognized as an
  algebraic identity (physics N1), Niven's set labeled precisely (N2), the
  get_counts() trap named on the VERIFY list (prior-work m6), the Marrakesh
  provenance labeled (m8), the Confirmations double-entry added to the gate
  order (m9), GLOSSARY anchors minted for band / VOID / double ratio / watcher
  factor / 7a-7b / record radius (m10). REFUTED from below, folds refused: stats
  B1's 3.2× readout-coverage gap (wrong base) and stats m1's claim that the PUB
  range was left unbumped (34–36 → 36–38 had already moved). Machinery: control
  config axis, ε_w-sign banks, (a)-order power report, corner banks into
  b_curve/LR/fallback via grid_and_corners. Debt carried: em dashes now 544 (the
  v32 count 460 was true at its fold; the round-31 reviewer's 461 counted the
  post-validate edits; this fold added more: one style campaign, one sweep, later)
  and prior-work n12 (CLAUDE.md names docs/engineering/, which does not exist:
  not this file's defect, noted for the repo). Post-fix states are new work;
  round 32 follows, WITH the prior-work lens.
- **v32 (2026-08-04, round 30, three fresh lenses on v31, NOT clean: spec 1 BLOCKER /
  6 MAJOR, stats 1 BLOCKER / 6 MAJOR, physics 1 BLOCKER / 1 MAJOR; every finding
  recomputed from below before folding; the physics core re-confirmed exactly, the
  write-miss cancellation verified to 1.5 rad, RECOGNIZED in round 31 as an
  ALGEBRAIC IDENTITY, not a probed range: the write angle enters Ŝ(j; r) only as
  the common factor sin θ_Sj, which divides out against the r = 0 normalizer, for
  ALL θ_Sj ∈ (0, π) and for single-write misses too; the verified-range framing is
  retired so no future round "extends" it):** physics BLOCKER, the deepest
  find of the campaign: a maximally MIXED k (m = 0) never watches and is INVISIBLE to
  every flown PUB (all Z-diagonal on k; |+⟩ and I/2 share diag(ρ_k)), so the flown
  data alone cannot distinguish "watcher watched and forgave" from "classical noise
  that never watched", folded as the WATCHER-TOMOGRAPH PUB pair (⟨Y_k Z_j⟩, +1 at
  r = 1 superposed, 0 for any Z-diagonal k; diagnostic-only, gates the mechanism
  sentence like V_S gates MI), the round-29 "Z_k = THE superposition certificate"
  sentence re-scoped to the m-polarised branch; stats BLOCKER: impz's escape window
  was uncomputed while the coverage rule made its sole discriminating band as large
  as possible (the EIGHTH protection-interaction instance, the first BETWEEN two
  consecutive amendments): the z_k DISCRIMINATION CEILING pinned (guard_z_k ≤ 0.05,
  envelope-above-ceiling unfreezable) + the impz BOUNDARY control at m = the frozen
  band value with P(CONFIRMED) reported; spec BLOCKER: b_blind carried TWO
  derivation bases one round after joining the coverage envelope (clause (a)'s old
  worst-basis sentence retired; the per-source collapse argument restated). MAJORs:
  the fixed-j′ corners enter the GATES (power/false-VOID on held-out corner twins,
  the fixed-vs-drawn nuisance law pinned); the off-diagonal BOX corner (λ = 0 ×
  worst fixed readout) and the η-TOP point (j′-differential, η(2) ≈ 1.002; η can
  structurally exceed 1) join the coverage grid; b_curve's own site now names the
  full 3-D grid (105 points; E[max] re-sized 2.53×, 210 → 2.76×); the negative
  controls sweep sign × level with impz as the FOURTH control kind (detection
  ≥ 0.999 asserted) and per-bank UCB stated; the false-VOID corner maximum now GATES
  at 0.5%; the interleave slots of r = 1.75 and r = 0.25 swapped (drift was EASING
  the (b) raw ordering by ~+0.012); the T̂ worked example restated (wrong-arm ×
  wrong-rate had compensated; the sup's worst slot distance is 4, the probe P ≥ 0.99
  criterion is the pinned trip rule); the ACCEPTED-SET sentence stated (+3.7%
  watcher-angle scale ∧ the flip; clause (c) binds); branch 1 gains the Aer-vs-model
  BANK comparison (moments + 8×8 correlation, 5-SE loud-only screen); the v31
  "REFUTED" SE-ratio wording softened to configuration-dependent (the round-29 full
  chain 0.88, the reviewer's shot-only model 1.20–1.26, the v32 gate's no-readout
  isolation 0.94: three configurations, three answers; the gate reports both of its
  own, no assert rides it); b_forgive's sign gloss resolved by RETIRING the "CONFIRM-HARDENING" label
  (lower envelope over signs; never-ease lives in the +SE edge; (b)'s
  false-CONFIRMED duty rests on the sign of ρ̂(2), certified by the swept controls);
  the seed-cost scheme pinned (deep ≥ 10⁶ on worst-basis banks, argmax-source
  reported, shallow argmax re-run deep). Numbers: the −sin gap 1.4×10⁻⁵ → 1.5×10⁻⁶
  (was computed at 3× budget); impostor margins re-based 6.2×/8.8× at the rehearsal
  envelope 0.0802; the raw-order margin pair (5.5σ vs 12σ) marked re-measure-at-
  freeze. Round-29 line-item debt paid: the stale line-776 reference named by role;
  the guard bank named precisely (5 pooled + 3 floors, the floors were never
  pooled); the abort class numbers glossed (position, not severity); supersession
  tags added (v26 hardening label, v29 fallback sign); the Z-leak all-four scope
  stated; V_S assigned to the GEV/empirical branch; the F135 registry inbound link
  added (rides this file's freeze commit); the GO position pinned at the gate order
  (the 08:12 advance go discharges it only with every gate green). Debt carried,
  recorded: the em-dash count (460 at this fold, a style campaign patient) and round-30 spec n3
  (a self-cancelling sentence the banked verdict names without quoting: to be
  located or discharged in round 31). Machinery folds validated by rerun; post-fix
  states are new work; round 31 follows this fold.
- **v31 (2026-08-04, round 29: the FIRST empty round on the amended v30; three
  fresh lenses, NOT clean: spec 5 MAJOR / physics 3 MAJOR / stats 4 MAJOR, every
  finding recomputed from below before folding, and TWO reviewer claims REFUTED
  from below in the process):** stats M1 = the LR threshold was the one band never
  swept over the SLOPE axis (valid device at the −1% slope corner: family
  DEVICE-DEVIATION ~1% vs the 0.13% pin, worsening with shots); LR now
  corner-conditional max-over-grid like b_curve; stats M3+m3 = the j-vs-j′ draw is
  a MIXTURE and a device has one fixed j′; fixed-j′-corner banks + grid halves
  join the floors' lower envelope; stats M4 = the (b)-fallback band rides the
  level axis and "worst sign" is undefined two-sided; union-of-edges envelope
  over the full grid; stats m1 = b_blind joins the coverage envelope (upper cut,
  level-riding center); physics M1 = the readout-bias basis pinned in one
  equation (Δ = the S-branch assignment differential; "2%" = Δ = 0.02, the wider
  reading; abort compares in Δ); physics M2 = the amplification pair's k
  preparation pinned (Z-eigenstate pair, signed conditional read; the |+⟩
  reading is second-order blind); physics M3 = the Z-POLARISED never-watching k
  reproduces the whole signed curve exactly (|⟨Z_k⟩| re-scoped as the watcher's
  superposition certificate; the "no impostor can fake the flip" sentence
  scoped); spec M1 = the gate-order line still carried the one-basis wording;
  spec M2 = the "η worst at r = 2" prediction reconciled with the round-28 flat
  measurement; spec M3 = the v30 edit had eaten the v29 revision-note header
  (restored); spec M4 = rehearsal readings are NOT frozen-constants entries;
  spec M5 = WIP paths removed from this tracked file; plus: the counts-level
  gate defined; D̂(j′; r) floor scope = the eight r > 0 arms; lexicographic item
  4 direction = MAXIMIZE; false-VOID measured at every coverage point; the T̂
  detection probe added; the round-25 floor sanity pin became a CONSTRUCTION
  (floor = min(envelope, worst valid mean − 10·shot-SE): a p0.03 envelope
  alone sits only ~3.4 own-SDs below its source mean, so pin and envelope must
  combine); ρ̂(1) center exactified to −sin(ε_watch+ε_ZZ)/cos(ε_ZZ); the
  clause-(a) "two orders" figure corrected (7.7×/10.9×; impostor (iii) passes
  (a) and dies at (b)); b_sym relabeled (it reads the η(2) level); b_forgive's
  self-contradictory sign warning reworded. REFUTED FROM BELOW, folds refused:
  the reviewer-proposed SE(ρ̂(1)) > SE(ρ̂(2)) assert (true for the NUMERATOR
  records; the double ratio weights denominator/normalizer noise by the value
  and inverts it, measured 0.88 fixed-device, demoted to a report; SOFTENED in
  v32, round-30 physics MAJOR: "refuted" over-claimed in the other direction,
  the 0.88 is the FULL chain incl. confusion + CAL + mitigation, the reviewer's
  shot-noise-only model reads 1.20–1.26, and the v32 gate's own no-readout
  isolation reads 0.94; the ratio is CONFIGURATION-DEPENDENT,
  not structural either way, and the gate now reports both provenances) and the
  −tan(ε_w) center (exact only for a pure-static miss; the exact center is the
  mixed form above). Physics core re-confirmed exactly by the physics lens
  (its whole table, incl. write-miss cancellation now verified to 1.2 rad,
  round-27's 0.3 understated). Machinery folds validated by rerun; post-fix
  states are new work; the next empty round follows this fold.
- **v30 (2026-08-04, round 28: the SIM-GATE stage; the finder was the GATE,
  not a reviewer: both findings surfaced as a joint-power collapse to 0.10 at
  a corner of the sign-and-level-swept end-to-end power measurement, exactly
  the blindness class the round-23/24 sweeps were built to catch):**
  (1) GUARD-BAND BASIS AMENDED: the v29 role-split sentence froze the guard
  bands at the worst-admitted basis alongside the magnitude bands, but T̂'s
  drift term scales with the RECORD LEVEL (exact centers: sup|T̂_j| = 0.190
  worst-basis vs 0.235 at the best-admitted corner, band 0.244 → the BEST
  valid device terminally VOIDed in 269/300 runs); the upper guard bands
  (T̂ per witness, |⟨Z_k⟩|, |⟨Z_S⟩|, Z-leak, V_S) now freeze from the
  p99.97 + SE upper envelope over the FULL coverage range (worst basis +
  the level × slope × sign grid); floors unaffected (lower cut: the worst
  basis IS the envelope; verified, no floor trips post-fix). The SEVENTH
  protection-interaction instance; the class stays ENUMERATED, never closed.
  (2) η_nom PROFILE ARITHMETIC PINNED: η_nom(r) = |ρ̂_Aer(r)/cos(r·π/2)|
  with r = 1 neighbor-interpolated, never the raw |ρ̂| profile (the rehearsal
  briefly stored |cos·η|, bending the (c) target into cos·|cos|; corner power
  0.10 before the gate caught it). Also recorded at this stage (REHEARSAL
  readings, not frozen values: the binding freeze re-derives every one of
  them from the freeze-time snapshot and its re-pinned line, round 29): the
  representative line by the pinned rule on the 2026-08-04T06:07Z snapshot is
  Kingston k–j–S–j′ = 92–93–94–95; transpile VERIFY all-green against the
  real target (rzz-native, budgets exact, writes exactly RZZ(π/2), two
  distinct watcher pulses, no routing; noiseless Aer parity 5.8×10⁻¹⁶); rzz
  duration 68 ns per edge, ANGLE-INDEPENDENT → η ~r-flat (the conservative
  case the VERIFY list anticipated) and the duration-slope reading is ~0 →
  s_res_freeze takes the pinned 1% absolute floor.
- **v29 (2026-07-18, round 27, the BANK: physics CLEAN twenty-first consecutive,
  verified BEYOND the claim (the double-ratio verdict exact under static write misses
  to 0.3 rad; the r-independent miss cancels); stats "essentially CLEAN … nothing
  blocks the freeze" (joint power ≥ 0.9946 reproduced from below, false-CONFIRMED
  0/2×10⁵ at both signs, impostors 94–266 SE from the (b) boundary); spec NOT CLEAN,
  one MAJOR, and the campaign ENDED BY DECISION, Tom: "das Experiment war super, lass
  uns das so banken", the loop truth cannot end, ended from outside it):** spec
  MAJOR = the clause-(b) fallback band on |ρ̂(1.75)| was verdict-bearing,
  sign-sensitive, and NOT in the sign sweep: the SIXTH sign-class instance, found ONE
  round after the accounting was declared "complete", independently converged on by
  the stats lens; verified from below by hand: ∓0.043 shift at δ = ±2%, odd in the
  sign, the round-24 mechanism at r = 1.75, folded: the fallback band joins the sweep
  at the per-statistic worst sign (SUPERSEDED in v31, stats M4: union-of-edges
  envelope over the full grid; "worst sign" is undefined for a two-sided band), the sim gate asserts whether the branch ever
  engages on the η×sign grid (from below: never in coverage, trigger near
  η(2) ≲ 0.42), and the completeness claim is RETIRED; membership is ENUMERATED,
  never again "complete" (the lesson's third strike: never declare a class closed,
  only a round clean); stats M1 = the negative controls certified at an ASSUMED
  confirm-easing sign; "swept, not assumed" extended to them, UCB required at every
  sign (review MC: the assumption was empirically correct, 94.0 vs 103.5 SE, impact
  nil, but an assumption is not a rule); stats M2 = the ρ̂(2)-bias headline "1–3.5 SE
  per 1–2%" was per-basis ambiguous and read low, pinned ≈ 2.75 SE per 1% (round-27
  MC ~5.5 SE per 2% at the pinned shots, inside the reviewed 0.7–9-per-2% span whose
  WORST model the injection already uses); spec MINOR 1 = the frozen ε_ZZ budget had
  no pinned Δf input, pinned to the Class-2 abort ceiling itself, ε_ZZ budget ≡ 0.01
  rad; in-job Δf feeds ONLY the abort and the decoding table, never a frozen band;
  spec MINOR 2 = stale 10⁵-seed framing at the guards' operative-tail site updated to
  the pinned ≥ 10⁶ (~300 events, ~6% count error); NITs folded: earliest-submit
  tie-break scoped in place (governs only the both-complete case); every-PUB per-shot
  memory tagged deliberately broader than the two consuming arms; physics NITs
  recorded no-change (the abort-to-η "2%" end = the conservative channel case;
  impostor-(i) = the maximal-dephasing idealization, scoped as a labeling aid). The
  round ledger stays IN the banked file deliberately: for this bank the history IS the
  record.
- **v28 (2026-07-18, round 26; physics CLEAN twentieth consecutive: the r = 1 record
  grounded as (k, j) = ¼(I + Y_k Z_j), a perfect classical bit; spec MAJOR and stats
  MINOR converged INDEPENDENTLY on the same single finding):** the LR threshold, the
  fifth and, per the round-26 stats sign-audit, LAST unswept sign-sensitive statistic,
  joins the sign sweep (reference distribution over η(2) × sign {−2%, 0, +2%},
  threshold at the worst sign; review MC suggests hardening was already LR-adverse,
  0.0013 vs 0.0000, but a default sign is never assumed adverse; b_blind is the
  standing counterexample); non-monotone-duration remedy fixed (per-arm η_nom re-read
  folds the shape into the target; a linear s_res rescale cannot); V_S threshold
  inherits the worst-sign rule; sign-envelope vs SE-edge declared INDEPENDENT opposite
  axes (never conflate); SE_contrast sign dominated by the swept gate (stated); Cauchy
  headline ~27–31%; the watcher-record prose made exact (maximally CORRELATED with the
  witness, classically: entanglement entropy 1 bit against the rest).
- **v27 (2026-07-18, round 25; physics CLEAN nineteenth consecutive, one NIT, folded:
  the watcher's r = 1 record is CLASSICAL which-path, concurrence C(k, j) = 0; stats
  MAJOR = the FOURTH sign-class instance; spec MAJOR = the last stale coherent-headroom
  site):** the pooled transverse guards T̂ join the upper-envelope-over-sign rule (T̂ is
  Ŝ's perpendicular twin, S-conditioned, so the S-antisymmetric offset enters the
  variant-A/X channel at full magnitude; a single-sign-frozen T̂ band terminally VOIDed
  admitted devices at up to ~96% in review MC); the readout-bias injection pinned
  TWO-CHANNEL (X and Y: a Y-only injection freezes the transverse band blind at every
  sign); Z-leak/⟨Z⟩ exempted with reason (S-symmetric, the offset averages out); the
  η_min "coherent headroom" wording at the abort-to-η one-liner: the sole survivor
  of the round-24 "end to end" sweep; named by role, round 30 spec debt: the old
  line-number reference went stale with every edit,
  reconciled to INCOHERENT-only; s_res day-of zero-headroom algebra stated as
  deliberate; the V_S budgeted δ pinned = ε_cal = 0.01 rad; σ-scale disambiguation
  (shot-SE sanity pin vs systematics-broadened p0.03 cut); b_blind's lower-edge
  convention cost acknowledged (clause (a) carries ~zero false-CONFIRMED duty); (c)'s
  SUP reference = the (b)-passing subset, stated; "mild" → MODEST upward max-bias
  (≈ 2.1 × per-point SE); "coverage family" label (the admitted midpoint would be
  0.95); contrast range annotated with the cross-model spread (~84–170σ).
- **v26 (2026-07-18, round 24; physics CLEAN eighteenth consecutive, incl. the
  sign-direction assignments of the round-23 repair verified from below; spec no-BLOCKER
  with the η_min propagation MAJOR; stats found the THIRD BLOCKER):** the round-23
  per-floor worst-sign fix had NOT been extended to b_blind; the additive readout
  offset shifts ρ̂(1) nearly 1:1, and the (b)-HARDENING sign is precisely the
  (a)-EASING sign, so the uniform hardening freeze set b_blind ≈ 0.042 instead of the
  envelope ≈ 0.065 and routed valid aligned-sign devices to INCONCLUSIVE up to ~40%
  (end-to-end joint power 0.57 at the ±2% endpoint; the sign-swept gate would have
  refused the freeze, but the remedy ladder could not repair a mis-signed center: the
  pre-reg was UNFREEZABLE as written); GENERAL LAW pinned: every band and floor freezes
  at its PER-STATISTIC worst sign-combination: magnitude ceilings (b_blind, b_curve's
  SUP) from the UPPER ENVELOPE over {−2%, 0, +2%} (envelope fix verified: joint 0.996,
  false-CONFIRMED untouched 0/2×10⁵), b_forgive at hardening (verified sign-robust,
  ≥ 0.995 both signs; the "hardening" LABEL retired in v32; the rule is the lower
  envelope over signs), floors at their per-floor worst; the floor-vs-band asymmetry
  stated (ladder-exempt ⟹ direct construction; remediable ⟹ gate-watched); spec MAJOR
  = η pinned INCOHERENT-ONLY end to end (η_min's "coherent headroom" wording
  contradicted the round-23 no-double-count guard; coherent residuals live solely in
  ε_w, ±s_res, and the sign sweep); s_res allowance disambiguated (= ×1.25, never
  ×0.25, both sites); Cauchy fraction tagged model-dependent (~26–32% across review
  models); the 10⁻⁶ write-edge residual tagged as an order-of-magnitude assertion;
  stale v20/v21 note figures annotated.
- **v25 (2026-07-18, round 23; physics CLEAN seventeenth consecutive; spec no-MAJOR
  with a two-word scoping fix; stats found the campaign's SECOND BLOCKER):** the
  uniform CONFIRM-HARDENING injection sign, pinned in round-12 for the bands, also fed
  the VOID FLOORS, raising the frozen Ŝ(j; 0) floor so that a valid device with clean
  (0%) or easing-sign (−2%, admitted) readout correlation sat BELOW it: terminal
  statistical-guard VOID, not re-fly-eligible, up to ~0.8–1.0 probability, and every
  gate structurally blind (the joint-power/false-VOID gates ran valid devices at the
  hardening sign only; the negative controls ran easing only on null/impostors; the
  day-of pair measures the sign but never moves a band); FIXED: floors inject at the
  per-floor WORST sign (the fewer-false-void rule extended to the injected-bias sign;
  equivalently the lower envelope over {−2%, 0, +2%}), the frozen Ŝ(j; 0) floor must
  sit ≥ 10σ below the worst-sign valid mean AND ≥ 10σ above the null (trivially
  satisfiable at ~1.2×10²σ; review MC: a 0.80–0.90 floor routes the null identically
  and false-VOIDs nothing), and the gates sweep the readout SIGN as a third axis on
  valid devices; per-guard quantile labels marked NOMINAL-only (the measured aggregate
  governs); analyzed-job rule scoped ("tie-break governs only IF BOTH COMPLETE");
  η_nom read at the λ-scaled INCOHERENT-only config (no-double-count guard); Cauchy
  fraction unified ~27–29%; separation ~1.2×10²σ; inter-PUB order deliberately
  unpinned (acknowledged); re-fly inherits the machinery unchanged; EPC one-liner
  tagged NOMINAL.
- **v24 (2026-07-18, round 22; physics CLEAN sixteenth consecutive, stats CLEAN at
  blocker/major: BOTH lenses now close "ready to freeze"; all findings in the
  η_max ≡ 1 propagation layer):** spec MAJOR = the mid-range basis was double-defined
  after the round-21 fix (the round-17 per-parameter midpoint recipe is undefined for
  time constants at the ideal ceiling, and the 7 level points had no pinned map to Aer
  configs), re-pinned in OUTPUT-η space: target η(2) = (η_min + 1)/2, realized by ONE
  loss factor λ scaling all incoherent rates of the representative config, λ solved per
  target, λ = 0 the ideal end, the same map realizing every sweep point; sweep ranges
  relabeled COVERAGE [η_min, 1] (admitted stays ≈[0.92, 0.98]; the overshoot is the
  round-21 point); stats MINOR = b_curve's SE-inflation direction pinned (LOWER edge
  per the global rule; power lives in max-over-grid + the swept gate, never in
  loosening the SE step; the max-of-35 upward selection bias is POWER-SAFE; do not
  "correct" it; finite grid suffices by smoothness); GoF comparison bank = each
  statistic's freeze bank (10⁵ physics / ≥ 10⁶ guards); runner-VERIFY
  input-invalidation = the FOURTH re-freeze trigger (the "only" clause omitted it);
  b_sym added to the SIGNAL source list; 0.5% overload disambiguated; "ONE ANALYZED
  Batch job"; V_S both-writes assumption and the 0.18 = 0.9·sin(0.2) arithmetic made
  explicit.
- **v23 (2026-07-18, round 21; physics CLEAN fifteenth consecutive: "ready to freeze",
  flip gauge-free to ANY static magnitude, record verdict independent of the CQ premise
  to 0.2 rad; stats no-MAJOR, "one of the most carefully-constructed pre-registrations
  I have refereed"; spec found the campaign's FIRST BLOCKER):** η_max was DEFINED
  best-admitted but COMPUTED as the calibration MEDIAN (round-20 fold); the
  lexicographic line selection flies the BEST edges, so the flown device sat
  systematically ABOVE the swept range, silently reopening the round-15 power leak with
  the gate blind to it; fixed by η_max ≡ 1 (no abort bounds a device from above; the
  ideal edge is the only coverage-complete upper end, and the fix removes the
  freeze-time estimation entirely; sweep [η_min, 1], target at its midpoint, minimax
  intact); GoF quantile-agreement now at each statistic's OPERATING point (p0.03 guards
  were validated only at p99.87); false-VOID ceiling re-derived 0.5% at the new default
  (the old 2% was the p0.13-era figure and would re-admit the eliminated loss); joint
  power = the END-TO-END measured CONFIRMED fraction (the product gloss is upward-biased
 ; floor and clause (b) anti-correlate through the shared normalizer, corr ≈ −0.36);
  p0.03 guard freeze at ≥ 10⁶ seeds (30 → 300 tail events); spectator-ZZ wording honest
  (r > 0 duration differences ride the same budget); linear-slope residual axis
  assumption stated (curvature rides the Aer per-arm read); "factor floors" in the VOID
  list; in-situ T1/Ramsey PUBs = diagnostic only; ≥ 10σ metric named (the pinned LR);
  σ̂(0) justification precise.
- **v22 (2026-07-18, round 20; physics CLEAN fourteenth consecutive: the commuting-ZZ
  circuit is an EXACT realization of the continuous proof, no Trotter budget; stats
  no-MAJOR, "ready to freeze" pending the provenance line):** spec MAJOR = the round-19
  p99.97 flip was stated at two sites but three others still pinned floors and folded
  guards at p0.13/p99.87, and the ~0.24% arithmetic only reconciles if all 8 move,
  propagated: ALL 8 guard-bank tests (5 pooled + 3 floors) at p99.97/p0.03; PHYSICS
  bands (b_blind, b_forgive, b_curve, LR) carved out explicitly at p99.87/p0.13; the
  guards' 0.03% empirical tail honestly noted (~30 events at 10⁵ seeds, ~18% count
  error, bootstrap-SE absorbed); stats MINOR = provenance of the ~2%/~0.24% figures
  (the 2% was MEASURED above the naive 1.04% union because fitted tails under-covered
  folded statistics; the 0.24% is the naive union at the new default; the measured
  ceiling governs); representative line PINNED (same lexicographic rule on the
  freeze-time snapshot, recorded, feeds τ_freeze/η_nom/s_res_freeze); 2-D grid
  resolution pinned (7 × 5 = 35 points); η_max computed at freeze (ledger η(2) at the
  calibration-median basis); split scope corrected to the two r = 0 arms only; single
  frozen bias coefficient recorded in the ledger (the 0.7–9 SE range is review
  history); "~8–12σ typical"; resilience_level annotated V1-idiom; η_min headroom
  cross-ref.
- **v21 (2026-07-18, round 19; physics CLEAN thirteenth consecutive, incl. a NOISY
  from-below model independently confirming the η(r)/η(2) basis dependence is a real
  systematic; stats no-MAJOR with the deepest verification list yet):** spec MAJOR =
  s_res was sourced from the runner stage, downstream of the band freeze: now the
  τ-pattern split (s_res_freeze from the representative-line schedule AT the sim gate;
  day-of reconciliation: flown slope × 1.25 within the frozen ±s_res = Class-1 abort);
  guard quantiles flipped to p99.97 DEFAULT (statistical-guard VOIDs are terminal on a
  paid one-shot; union ~0.24% vs ~2% at the old default); Ŝ(j; 0) floor pinned SIGNED
  (> floor > 0, also carries the ratio's sign convention); the 2-D b_curve safety
  sentence made precise (impostors die at (b); the Cauchy null passes (b) ~25–29%,
  unified in v25 to ~27–29%, model-dependent, and
  is routed by the FLOOR at ~10²σ); design-remedy operational cost stated (tighter
  day-of aborts, higher abort probability, never hidden in "re-freeze"); drift
  interleaving = secondary unguaranteed benefit (delivered ≠ temporal order); stale
  "round-robin" naming fixed; API flag names deferred to runtime verification
  (SamplerV2 BitArray vs memory=True); "no INSERTED delays"; contrast = in-conjunction
  backstop; b_blind derivation point stated; abort-to-η mapping one-liner; 10⁵ seeds
  PER GRID POINT; transverse guard duty = drift differential (statics absorbed by the
  axis fit).
- **v20 (2026-07-18, round 18; physics CLEAN twelfth consecutive exact recompute):**
  stats MAJOR = the round-17 corner-conditional b_curve swept the LEVEL only while the
  round-16 gate stressed the SHAPE; band and gate covered different spaces, the floor
  crossed at a ~1.4% slope residual on the GOOD device, and the earlier cited ~93%
  shape figure was itself below the floor; now the b_curve sweep is 2-D (level ×
  ±s_res slope residual; s_res TO-FREEZE = VERIFY duration-slope × 25% allowance
  (= ×1.25, pinned v26),
  ≥ 1% absolute; review MC: 0.87 → 0.9999 for the +1%-slope 0.98 device, safe on
  false-CONFIRMED) and the gate's shape perturbation runs on the SAME grid; spec
  MAJOR = the even/odd split needs PER-SHOT MEMORY, but the doc worded everything
  "counts"; a faithful counts runner would make the verdict code unexecutable; pinned:
  memory=True as Class-1 abort + VERIFY line, RAW-memory persistence, pipeline order
  split-then-aggregate-then-mitigate-per-half; re-freeze lever named PER FAILURE MODE
  (shots for interior margins; DESIGN only for corner-limited; shots cannot move a
  shot-independent offset); interleave order pinned as a fixed permutation (r = 0 slot
  5, r = 2 slot 7); LR dressing basis = mid-range (same as (c)); v18 note's 9-arm LR
  annotated superseded; T̂ written out; campaign bound ~0.2% honest; contrast = echo of
  (b); GoF expected path = the 10⁵-empirical branch; Clifford-substitution VERIFY line
  (r = 2 gates stay two distinct pulses); min-angle RZZ(π/16) VERIFY line.
- **v19 (2026-07-18, round 17; physics CLEAN eleventh consecutive exact recompute,
  incl. gauge-freedom of the flip to 1.2 rad statics; two construction defects caught
  in the NEWEST layers, both mine from rounds 15–16):** the round-16 LR pin was
  mathematically SINGULAR (spec MAJOR: ρ̂(0) ≡ 1 identically → the 9×9 covariance has
  rank ≤ 8, Σ⁻¹ undefined: an unpinned implementation choice in a verdict branch; now
  the 8 non-degenerate arms, hypotheses η_nom-DRESSED so the LR separates shape/sign
  never attenuation, threshold reference swept over η); the round-15 pooled b_curve
  UNDER-COVERED the admitted corner (stats MAJOR: p99.87(mixture) ≤ max_η p99.87(·|η),
  corner offset shot-independent → doubling shots LOWERED corner power 97.1% → 95.7%;
  now CORNER-CONDITIONAL: max over sweep points of the per-point p99.87, SE-inflated,
  corner power 99.9% in review MC); re-fly-eligible VOID causes ENUMERATED (mechanical
  class + Class-2 trips only; statistical guard VOIDs terminal); mid-range basis pinned
  (arithmetic midpoint per parameter); "cheap margin" qualified ((a)/(b)/contrast only);
  GoF roles labeled (quantile-agreement primary, AD = egregious-misfit screen);
  both-bill worst-case cost stated (~2× ≈ 4% of gift); guard-count wording (5 pooled +
  3 floors); drift phrasing fixed (per-arm rate vs accumulated tilt, interleave slot
  anchored); hat convention line.
- **v18 (2026-07-18, round 16; physics CLEAN tenth consecutive exact recompute, incl.
  the V_S guard's correct scoping and gauge-freedom of the flip; stats VERIFIED the
  round-15 role-split fix from below, leak absent):** TERMINAL RULE pinned (spec MAJOR:
  the re-fly decision was post-data judgment, the campaign-level "fly until it
  confirms" door; now: single terminal flight, ONE conditional instrument-fault re-fly
  behind a fresh go, INCONCLUSIVE terminal, any new attempt = a new pre-reg cycle);
  Class-1/Class-2 language propagated out of the older sections (ζ book, estimator,
  banner sequence; "day-of abort" no longer names in-job VOID gates); LR construction
  pinned (Gaussian log-LR on the full 9-arm ρ̂ vector, point-curve hypotheses, ONE
  shared 9×9 held-out covariance; arms correlate through the r = 0 normalizer;
  SUPERSEDED in v19: 8 arms, 8×8; the 9-arm covariance is singular);
  inter-arm drift budget pinned ≤ 0.05 rad/arm and day-of-anchored through the
  transverse guard (T̂ band must cover the budget and trip beyond it; cliff at ~0.2
  rad/arm); the joint-power gate perturbs the η(r) SHAPE beside the level sweep (a
  slope error is invisible to the level sweep; ~93% vs a false ~99% read in review MC)
  and reports the corner MINIMUM net of guard drag; smaller: η_max defined, cost
  phrasing un-double-counted, contrast ~95–150σ.
- **v17 (2026-07-18, round 15; physics CLEAN ninth consecutive exact recompute, incl.
  the carrier surviving removal of the control write; two REAL design defects caught
  late):** (1) the clause-(c) SHAPE POWER LEAK (stats, BLOCKER-adjacent): with
  angle-dependent η the profile η(r)/η(2) depends on the η level, so the worst-admitted
  target basis, correct for clause (b)'s magnitude, systematically mismatches a
  BETTER-than-worst valid device (~2× b_curve at η(2) = 0.98 vs a 0.92 target; the
  better the hardware, the likelier INCONCLUSIVE), and the single-basis power gate was
  measured at the one point that hides it; fixed by the ROLE-SPLIT band basis
  (magnitude bands at worst-admitted; the (c) target at mid-range; b_curve's seed bank
  and the ≥ 95% joint-power gate SWEPT across the admitted η range, floor = the sweep
  minimum; swept MC restores 98.9–99.8%); (2) the ABORT TWO-CLASS split (spec MAJOR): a
  Batch has no mid-batch conditionals: Class 1 pre-flight aborts prevent the spend
  (calibration + local transpile inputs), Class 2 in-job VOID gates (ζ·τ, |8·ε_cal|,
  crosstalk) protect the verdict after the spend; gate order gained the post-Batch
  evaluation stage; (3) the v13 floor direction rule was self-contradictory; the floor
  is EXEMPT from the ladder entirely; (4) three stale revision notes annotated (388σ /
  "it alone", first-completed tie-break, ~50–100σ); smaller: "(b) margin =
  shot-independent quantile margin" tag, "normally unbilled" honest, ~1–2×10²σ,
  mean-estimator-SE precision on the √2 claim.
- **v16 (2026-07-18, round 14; NO MAJOR on any lens: physics CLEAN eighth consecutive
  exact recompute (2 no-action NITs), spec and stats number/framing MINORs only):**
  the "388σ" corrected to ~10²σ (388 was the law's own Ŝ(j;0) SNR, not the law-vs-null
  separation; P(null escapes the floor) = 0.000000 either way); "single guard" softened
  to ROUTING (the shape clause independently rejects the null with probability ≈ 1; the
  floor routes the degenerate case to VOID instead of INCONCLUSIVE); joint power /
  false-VOID / negative controls measured on a seed bank HELD OUT from band fitting
  (in-sample optimism measured ~0.04 pp: negligible, pinned anyway); the
  coefficient-tension sentence resolved (false-CONFIRMED = coefficient-free UCB;
  power-side injection = the one place a coefficient enters, worst reviewed one);
  "sim-gate-measured" where "gate" was ambiguous; GoF pinned numerically (AD p ≥ 0.01 +
  fitted p99.87 within 2 bootstrap SEs of the 10⁵-seed empirical); η_nom named as the τ
  two-reading analog; single-submit fallback written; sim-gate seeds recorded in the
  committed artifact; Law A hypotheses restated at the geometry; D̂-vs-D_j notation
  line; zero-tilt caveat on the Cauchy fraction (~27–29%, worst case).
- **v15 (2026-07-18, round 13; physics CLEAN seventh consecutive exact recompute,
  incl. sign-vs-frame-drift −cos(φ) robustness and the RZZ(4π) = identity check of the
  amplification sequence; stats CLEAN at blocker/major with the frozen-code MC table:
  LAW 1.00 CONFIRMED, null VOIDed with probability 1, every impostor 0 even at a 0.10
  easing bias):** the v14-leftover abort-gates line still justifying the crosstalk abort
  by the DROPPED center shift reworded to the pinned-sign justification (spec MAJOR:
  my round-12 fold missed the third site); v12 revision note annotated with its v14
  supersession; tie-break disambiguated to EARLIEST SUBMIT timestamp (data-independent);
  signal-model bias injection pinned at the UPPER end of the reviewed coefficient range
  (round-13: hardening protects power only at the conservative end); Ŝ(j; 0) floor
  BARRED from every remedy step (it alone VOIDs the Cauchy null, ~388σ separation;
  CORRECTED in v16: routing, not sole protection, and the separation reads ~10²σ);
  UCB seed floor pinned (10⁵ default, never below 10⁴: strict < 0.1% unsatisfiable
  under ~3001 seeds); contrast power relabeled ~95–120σ (recomputed); (b) corner margin
  honest (~3–6σ at the exact abort corner); "8 tests (5 pooled + 3 floors)"; r = 0 arm =
  "reference arm"; commit rides with the proof-corollary edits; η(1/2)/η(3/2) explicit;
  "Watching costs what the watcher learns" (the "exactly" was endpoint-only).
- **v14 (2026-07-18, round 12; physics CLEAN sixth consecutive exact recompute; spec
  no-MAJOR; stats found a REAL design interaction):** the rounds-8–10 readout-bias
  band-CENTER-shift mandate DROPPED; round-12 MC showed it redundant (0/20000 false
  CONFIRMED for null + all impostors with easing-sign injection on UNSHIFTED bands: a
  constant additive bias cannot jointly fake the flip, the r = 1 minimum, and the 7-arm
  shape) and power-fatal (the shift is comparable to the whole ~2.5–3σ (b) margin; no
  shot count repairs a center offset); replaced by the PINNED-SIGN rule: hardening sign
  in signal-model band runs, easing sign in negative controls (their ≤ 0.1% UCB is the
  model-independent protection), crosstalk pair stays the hardware abort; bias
  coefficient deferred to the gate-measured value (model estimates ranged 0.7–9 SE);
  re-freeze trigger self-contained (bootstrap SE > 1.3× model SE / GoF fail /
  joint-power fail, only before the commit hash); EMPTY ROUND on the frozen numbers
  added to the gate order; power labels honest (2.7σ→3.8σ is a sub-abort stress pair:
  abort-implied raw (b) margin ~8–14σ; the ~96% was stress-frozen-bands
  self-consistent); `_query_billing.py` marked external; floors "on every NORMALIZER
  factor" (Ŝ(j;r) deliberately unfloored); branch-1 response non-discretionary; impostor
  set = labeling aid, not exhaustive.
- **v13 (2026-07-18, round 11; FIRST no-MAJOR round on all three lenses: physics CLEAN
  fifth consecutive exact recompute, spec and stats MINOR/NIT only):** one-line claim
  re-worded to the tested face (flip + recovery at η(2); "IN FULL MAGNITUDE" was
  untestable under known attenuation); split-sample partition PINNED (even shot indices
  fit the axis, odd evaluate; the verdict is now a deterministic function of the
  counts); η ≈ 0.74 relabeled STRESS floor (physics: the abort thresholds imply η ≈
  0.92–0.98; the 0.74 was a round-9 MC sweep point promoted in folding); shots-line
  power figures carry the stress-floor caveat (~2.5σ / ~96% joint there); negative
  controls run WITH the worst-case readout-bias injection (the ≤ 0.1% bound formally
  covers the clause-(b)-easing pathway); branch 1 = feasibility smoke test, not a hard
  gate (single stochastic draw); twin-submit tie-break pinned (first-completed governs;
  SUPERSEDED in v15: earliest SUBMIT timestamp, data-independent);
  Law C sentence honest (not invoked anywhere); folded guards raised to 10⁵-seed
  empirical; contrast power ~50–100σ (RECOMPUTED in v16: ~95–120σ); "monotone recovery"
  → "neighbor recovery check,
  NOT a law-enforced ordering"; v11 note annotated with its v12 supersession; "a
  per-half-gate-dephased watcher" (indefinite article: two other dephasings reproduce
  the law).
- **v12 (2026-07-18, round 10; stats lens fully CLEAN for the first time: every figure
  reproduced or conservative; physics essentially clean, fourth consecutive exact
  recompute):** the day-of crosstalk measurement reclassified ABORT-ONLY (spec MAJOR: a
  measured day-of value was worded as a ledger CENTER input, contradicting the
  deterministic-ledger definition and freeze-before-data; the center shift is now frozen
  at the 2% worst case at the gate, the τ_dayof pattern; SUPERSEDED in v14: the frozen
  center shift itself was dropped, pinned signs replaced it); impostor-(i) passage rewritten
  (physics MINOR: my round-9 "50/50 mixture of the two watcher directions" clause read
  literally IS the law's curve, cos being even; both near-misses now pinned as
  implementer warnings; the correct mechanism is per-half-gate independent signs); the
  r = 2 mechanism grounded in the claim section (RZZ(π) is a deterministic Pauli; k ends
  pure, the even watcher learns nothing; watching costs exactly what the watcher learns);
  R_δ relabeled the record-redundancy count; banner band sentence rephrased as rule (was
  present-tense-as-done); abort checklist names τ_dayof; v10 note carries its v11
  narrowing inline; b_blind two-source coincidence stated; ε_watch used inside ε_w;
  resilience_level = 0 disambiguated from the offline CAL inversion; floors renamed
  Factor collapse floors (D̂(j′;0) is a numerator factor); record-on-Y stated exact for
  the ideal circuit; "(model estimate)" tag on the bias coefficient.
- **v11 (2026-07-18, round 9; physics lens CLEAN for the THIRD consecutive round: no
  MINOR left on the physics; spec lens returned its first no-MAJOR round):** the round-8
  "double ratio cancels ANY r-independent readout correlation" claim NARROWED to the
  symmetric multiplicative part (stats MAJOR, verified by algebra: an additive
  S-antisymmetric assignment offset does not cancel; ρ̂(2) bias ≈ 1–3.5 SE per 1–2%
  differential, one sign easing clause (b)); sim gate now injects BOTH modes; NEW in-job
  S↔witness crosstalk pair (S=|1⟩/|0⟩ swaps) = the hardware bound the mode lacked, day-of
  abort ≤ 2%, measured value into the ledger (SUPERSEDED in v12: ABORT-ONLY, the center
  is frozen at the 2% worst case; the day-of value never moves a band); r = 0 arms at
  2× shots (the split-sample
  normalizer feeds every band; +~9 s); τ split into τ_freeze (representative line, bands)
  and τ_dayof (flown line, the abort input); the flown line does not exist at freeze;
  amplification coverage assumption stated (π/2-referenced; fractional angles ride the
  fail-safe direction); clause evaluation order (a)→(b)→(c) pinned with short-circuit (no
  NaN path); impostor-(i) gloss corrected (NOT measure-and-forget; tracing k reproduces
  the coherent curve; the incoherence lives in the gate sign); the 0.94 power-trigger
  qualified as η-magnitude-dependent (≈ 0.958 at worst-admitted η); b_blind collapse to
  one quantile stated; variant-A role precise (no variant discarded); MI as labeled third
  category in Not-tested; 96/4 (was 92/8); 8 pooled tests unified; η_min/b_forgive
  unconflated; θ_Sw defined in place; γ = 0 single-watcher facet named; cost updated
  ~2.8–3.6 QPU min (~2% of gift).
- **v10 (2026-07-18, round 8; physics lens CLEAN again: second independent exact
  recompute, all arms to 1e-6, spectator-ZZ cancellation and common-write-miss
  cancellation verified):** impostor (i) PINNED to ρ_i(r) = cos²(r·π/4), the incoherent
  (sign-randomized) twin of the flown watcher: parameter-free, an LR needs a functional
  form; readout-correlation injection made ARM-DEPENDENT and the two sections reconciled
  (the double ratio provably cancels r-independent readout correlations, MC bias/SE
  ≤ 0.02, so only the arm-dependent residual is dangerous, and the injection is its sole
  line of defense; NARROWED in v11: the cancellation holds for the SYMMETRIC part only); shots DOUBLED to 16384 (joint power well off the 95% floor at ~1
  extra QPU minute; cost now ~2.5–3.3 min ≈ 2% of the gift); negative-control acceptance
  = 95% Poisson UCB < 0.1%; line selection pinned lexicographically ("house score" was
  undefined); (b)-fallback overlap with (c) owned; the Gaussian-only caveat on
  quantile-VALUE error figures; banner TO-FREEZE sentence reconciled with the ledger;
  R̂'s v4-note identity marked historical; ε_ZZ ≡ π·Δf·τ written; V_S blindness corrected
  to ~100× at abort scale (was 20×, loose); T2-vs-T2* "same order" wording; power-trigger
  0.94 vs gap-vanish 0.924 both stated; η/η_nom/η_min reading aid.
- **v9 (2026-07-18, round 7; physics lens returned CLEAN: exact statevector recompute
  reproduced all 9 arms, the π-flip, and every linearization):** guard-bank j↔S slip
  fixed (the ⟨Z⟩ sanity streams are k and S, the Z-measured ideal-0 qubits; j is
  pre-rotated and has no Z stream; |⟨Z_S⟩| doubles as branch-balance sanity); CZ fallback
  = NO-FLIGHT (no CZ bands exist, so no CZ flight); joint P(CONFIRMED | valid device)
  assembled over the seed bank and pinned ≥ 95% (a single positive-control draw proves
  feasibility, never power); the two mislabeled quantile-MC-error figures relabeled
  (count vs value: 62%/6% at 2000 seeds, 28%/2.7% at 10⁴); (b)-order power computed at
  the measured η_nom(2)/η_nom(1.75); b_blind half-width pinned about the model center
  (the center enters once); readout-correlation treated as BIAS (center shift), not only
  width; Z-leak pair defined; V_S guard = 4 PUBs + amplification pair enumerated, cost
  re-estimated ~1.3–1.8 QPU min; smaller: "both pinned Herons", ρ(2) ideal-vs-flown
  disambiguated, corollary cross-ref precise, "inflated" defined, D_j = β_j·|sin θ_w|
  coincidence stated, T2-vs-T2* label, variant-A-carries-no-record note.
- **v8 (2026-07-18, round 6):** band-source degeneracy fixed; the NULL zeroes the
  denominator's record, so the SIGNAL model now sources ALL fluctuation widths (b_blind
  tail read at the r = 1 arm), deterministic centers from the systematics ledger, NULL
  demoted to negative control (P(CONFIRMED | null/impostor) ≤ 0.1%); 10⁵ seeds default;
  split-sample normalizer as THE definition; LR threshold pinned family-wise; two-branch
  positive control; consolidated day-of aborts incl. |8·ε_cal|; per-device gate runs;
  runner VERIFY list; (b)-order power fallback.
- **v7 (2026-07-18, round 5):** signed-coherence corollary MINTED into
  PROOF_RECORD_PARITY_LAW (the sign discriminator now rests on a written theorem); guard
  bank POOLED ~90 → ~8 tests (union-bound false VOID was ~11%); SE-direction per decision
  role; ε_cal in-job amplification anchor; NULL/SIGNAL generative models first pinned.
- **v6 (2026-07-18, round 4):** the sign blocker, all three reviewers independently:
  clauses were magnitude-typed while the carrier goes ρ(2) = −1; every clause re-typed
  SIGNED; third (sign-flat) impostor added; positive control mandated; V_S guard
  redefined j-conditional (the unconditional read is O(δ²)-blind).
- **v5 (2026-07-18, round 3):** SIGNED projection Ŝ with per-witness axes; fitted-tail
  bands (never μ±3σ); VOID precedence; the ζ factor-2 book θ = π·Δf·τ; V_S guard added;
  shape clause σ̂; LR deviation statistic; interleave rule.
- **v4 (2026-07-18, round 2):** ibm_kingston pinned primary (Tom's fresh pulls); verdict
  carrier promoted to the DOUBLE ratio (then written R̂(r)/R̂(0); the round-3 signed
  projection later replaced its numerator; today's R̂ is the informational-only single
  ratio of the b_sym check); axis-projected estimator; parametric tail fits; ζ_jk
  measured in-job; ε_watch = 2·ε_cal.
- **v3 (2026-07-18, round 1):** sim gate must INJECT coherent ZZ (standard Aer has none;
  b_blind responds linearly to ε_w); half-angle watcher pair pinned; η one-sided budget;
  freezing rules + impostor targets; variant-B pre-rotation order caught.
- **v2 (2026-07-18):** full pre-registration skeleton: 9 arms × 2 variants, 8192 shots,
  conditional-Bloch estimator, verdict clauses drafted.
- **v1 (2026-07-18):** first draft from Tom's fixpoint frame: line k–j–S–j′, ratio
  carrier, prediction R(r) = |cos(r·π/2)|.
