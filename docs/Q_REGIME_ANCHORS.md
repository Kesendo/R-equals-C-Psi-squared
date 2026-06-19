# Q-Regime Anchor Map (Q = J/γ₀)

The framework's three Q-bands (onset, peak, plateau) materialized as a table of (Q, J)-pairs at the code-convention γ₀ = 0.05. Use this when reading any script that hard-codes a J value: divide by γ₀ to recover the Q-regime it sits in.

## Anchor table

| Q anchor | J = Q·γ₀ (at γ₀=0.05) | Band / role          | Source / tier                                    |
|----------|------------------------|----------------------|--------------------------------------------------|
| **0.2**  | 0.010                  | onset start          | Q-band edge (Tier2Empirical)                     |
| **0.35** | 0.0175                 | onset end            | Q-band edge (Tier2Empirical)                     |
| **1.0**  | 0.050 = γ₀             | Balance              | J = γ₀ exactly (Tier1Derived)                    |
| **1.2**  | 0.060                  | peak start           | Q-band edge (Tier2Empirical)                     |
| **1.5**  | 0.075                  | F86 Q_peak (c=2)     | **PolarityPairQPeakDecompositionClaim** (Tier1Derived schema 2 − 1/2); PerBlockQPeakClaim wobble 1.4-1.6 empirical witness |
| **1.6**  | 0.080                  | F86 Q_peak (c=3)     | PerBlockQPeakClaim, saturated N=5..9 (Tier2Empirical) |
| **√3**   | ≈ 0.0866               | canonical θ=60° (Lindblad-Absorption-Match) | **LindbladAbsorptionMatchAtSixtyDegreesClaim** (Tier1Derived: \|λ_±\|=2γ₀=α at θ=60°); empirical F86 c=3 K-peak suggestively close at N=7,8 (within ~2%) but Tier-3/4 hypothesis NOT typed |
| **1.8**  | 0.090                  | F86 Q_peak (c=4, c=5) + peak end | PerBlockQPeakClaim, saturated N≥7 (Tier2Empirical); fine-grid: c=4 N=7=1.748, N=8=1.804 |
| **2.0**  | 0.100                  | Q_EP idealized       | Q_EP = 2/g_eff at g_eff=1 (Tier1Derived)         |
| **2.5**  | 0.125                  | Endpoint orbit Q     | **PolarityPairQPeakDecompositionClaim** (Tier1Derived schema 2 + 1/2); PerF71OrbitObservation empirical witness (stable ~2% across c=2..4 N=5..8) |

Standard resonance-scan grid (`ResonanceScan.DefaultQGrid()`): `LinearQGrid(0.20, 4.00, 153)`. These 10 named anchors are the structurally-meaningful points; the grid samples the continuum between and past them.

**The peak-band 1.2-1.8 is c-mediated**: Q_peak saturates with chromaticity (1.5 c=2, 1.6 c=3, 1.8 c=4,5). The earlier linear-extrapolation hypothesis Q_peak(c=5)=2.0 was refuted by the 40h N=9 Phase-2 run (2026-04-24). The peak-band width reflects c-variation, NOT finite-size scatter.

Per-orbit Q_EP values, HWHM ratios, fitted (α, β) sub-class parameters, and off-grid escape Q values (Q≈7-17 at high-N orbit plateaus) are empirical observations rather than anchors. See [`F86_VALUES_INVENTORY.md`](F86_VALUES_INVENTORY.md) for the full inventory with caveats — they should not be promoted to anchor status without further structural derivation.

## Per-anchor viewpoints

Each Q-anchor is reached from multiple framework lenses (parallel to the multi-viewpoint edges in `FractionReferenceGraph`). The viewpoints listed below are sources that independently land on the same Q-value — Painter-Principle pluralism applied to Q.

### Q = 0.2 (onset start)
- **Q-band edge**: lower bound of onset band (per `project_q_middle_structure`: bands {onset 0.2-0.35, peak 1.2-1.8, plateau block-specific})
- **F95 angle**: θ = arctan(0.2) ≈ 11.3°
- **H-clock**: 0.2 H-rotations per γ₀-tick (carrier-decay dominates 5:1)
- **Resonance grid**: lower bound of `ResonanceScan.DefaultQGrid()` = `LinearQGrid(0.20, 4.00, 153)`
- **J at γ₀=0.05**: J = 0.010

### Q = 0.35 (onset end)
- **Q-band edge**: upper bound of onset band
- **F95 angle**: θ = arctan(0.35) ≈ 19.3°
- **H-clock**: 0.35 H-rotations per γ₀-tick
- **Empirical onset W>0**: dressed-mode weight emerges at Q ≈ 0.3-0.4 (per `project_q_middle_structure`: "Q ~ 0.4: onset of dressed-mode emergence")
- **J at γ₀=0.05**: J = 0.0175

### Q = 1.0 (Balance, J = γ₀)
- **Identity**: J = γ₀ exactly — the only point where the two clocks have identical numerical rates
- **F95 angle**: θ = arctan(1.0) = 45° — diagonal in the Re/Im plane of Liouvillian eigenvalues
- **H-clock synchron**: 1 H-rotation per 1 γ₀-decay (perfectly synchronized)
- **Wave-breaking-scan label**: "Balance (J=γ₀)" (`wave_breaking_q_anchor_scan.py:14`)
- **F99-anker hits at Q=1**: N=2 → 3/8 (KIntermediate), N=3 → 1/2 (Generic), N=5 → 1/4 (Silver-Dicke) — three independent N's land on F99 ankers, structurally privileged
- **C# typed**: `FractionReferenceGraph.QBasisAnkers[0]`
- **J at γ₀=0.05**: J = 0.050

### Q = 1.2 (peak start)
- **Q-band edge**: lower bound of peak band (per `project_q_middle_structure`)
- **F95 angle**: θ = arctan(1.2) ≈ 50.2°
- **H-clock**: 1.2 H-rotations per γ₀-tick (H slightly ahead)
- **J at γ₀=0.05**: J = 0.060

### Q = 1.5 (F86 Q_peak c=2)
- **F86 identification**: `PerBlockQPeakClaim` c=2 saturated value (wobble 1.4-1.6 across N=4..9) — Tier1Candidate, finite-size sensitive
- **K observable**: K_CC_pr (J-derivative of S over t) maximizes here at c=2
- **F95 angle**: θ = arctan(1.5) ≈ 56.3°
- **H-clock**: 1.5 H-rotations per γ₀-tick — the slight H-favor that defines the resonance maximum
- **Wave-breaking-scan label**: "F86 Q_peak" (`wave_breaking_q_anchor_scan.py:15`)
- **C# typed**: `FractionReferenceGraph.QBasisAnkers[1]`, `PerBlockQPeakClaim.Standard[0]`
- **J at γ₀=0.05**: J = 0.075
- **Birth Canal mouth (flow-EP reading)**: [`THE_FLOW_BETWEEN_TWO_SINGULARITIES`](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md) places the **single-excitation flow exceptional point** at Q_EP = 1.5 (J_EP = Q_EP·γ₀ = 0.075, the "EP-family coupling" named in `ClockHandLadderWitness`), the **Birth Canal's birth singularity** where the slowly-decaying oscillating memory is born (the canal is bracketed by this EP and the λ=0 / equipartition kernel). So Q=1.5 being the project-wide default operating point is not incidental: J=0.075 was chosen to sit on the flow EP. **Three caveats, do not conflate**: (i) this single-excitation flow-EP at 1.5 is a *different mode* from the F86/idealized **Q_EP = 2.0** anchor below (= 2/g_eff at g_eff=1); in that convention Q=1.5 = Q_EP − ½ = the Q_peak (= the "marginal regime near the EP", `ClockHandLadderWitness`). Two EP readings near one Q. (ii) the canal's *boundary* (sterile↔birth-canal) is a separate, **Q-independent** γ-profile surface (coordinate s*≈0.709, a rejected 1/√2 lookalike), not this point. (iii) operating-point-on-the-flow-EP does **not** imply the hardware chip sits on the EP: the chip's natural Q ≈ 30 (deep memory), reached down to the EP only by injected noise; the toy γ₀=1 units make the chip's coupling falsely read 1.5.

### Q = 1.6 (F86 Q_peak c=3)
- **F86 identification**: `PerBlockQPeakClaim.Standard[1]` c=3 saturated value (N=5..9, N-invariant) — Tier2Empirical
- **F95 angle**: θ = arctan(1.6) ≈ 58.0°
- **H-clock**: 1.6 H-rotations per γ₀-tick
- **J at γ₀=0.05**: J = 0.080

### Q = 1.8 (F86 Q_peak c=4, c=5 + peak end)
- **F86 identification**: `PerBlockQPeakClaim.Standard[2,3]` c=4 and c=5 saturated value — Tier2Empirical plateau ("saturates at 1.8 for c ≥ 4" per `project_q_middle_structure`)
- **Q-band edge**: upper bound of peak band
- **F95 angle**: θ = arctan(1.8) ≈ 60.9°
- **H-clock**: 1.8 H-rotations per γ₀-tick
- **J at γ₀=0.05**: J = 0.090

### Q = 2.0 (Q_EP idealized)
- **F86 identification**: Q_EP = 2/g_eff at g_eff=1 (`QEpLaw.cs:14-27`) — Tier1Derived idealized
- **F95 angle**: θ = arctan(2.0) ≈ 63.4°
- **H-clock**: 2 H-rotations per γ₀-tick
- **Wave-breaking-scan label**: "Q_EP (g_eff=1)" (`wave_breaking_q_anchor_scan.py:16`)
- **C# typed**: `FractionReferenceGraph.QBasisAnkers[2]`, `QEpLaw`
- **J at γ₀=0.05**: J = 0.100

### Q = 2.5 (Endpoint orbit, candidate)
- **F86 identification**: `PerF71OrbitObservation` orbit 0 (Endpoint) stable at Q ≈ 2.5 across (c=2..4, N=5..8); 2.39-2.61 range, ~2% N-variation
- **F95 angle**: θ = arctan(2.5) ≈ 68.2°
- **H-clock**: 2.5 H-rotations per γ₀-tick
- **Tier-1 candidate**: stability across 9 (c, N) combinations is reasonable Tier-1 evidence
- **C# typed**: absent from `QBasisAnkers`; needs structural derivation for formal anchor promotion
- **J at γ₀=0.05**: J = 0.125

## γ₀ as time-tick: why "20" appears in the headers

γ₀ is the framework's time-tick: every dimensionful timescale is integer × 1/γ₀ (see [`reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md)). At γ₀ = 0.05:

- 1/γ₀ = **20 sim-time-units** = the γ₀-period (always, regardless of J)
- This is a unit-arithmetic fact: γ₀ = 1/20 in sim units → 1/γ₀ = 20 sim units

The **structurally-invariant** content lives in Q = J/γ₀:

- 1/γ₀ = **Q Hamiltonian-rotation-periods** (since 1 H-period = 1/J sim-units, and Q·(1/J) = 1/γ₀)
- At framework J values (table above, J ∈ [0.01, 0.1]), Q ∈ [0.2, 2.0] → 1/γ₀ encompasses 0.2 to 2.0 H-rotations
- At Balance (Q=1): 1/γ₀ = 1 H-rotation (γ₀-period and H-period perfectly synchronized)
- At F86 peak (Q=1.5): 1/γ₀ = 1.5 H-rotations
- At Q_EP/g_eff=1 (Q=2): 1/γ₀ = 2 H-rotations

**Caveat on the Propagate cockpit default**: `compute/RCPsiSquared.Propagate/Program.cs:1118` sets `J = 1.0`, which gives Q = J/γ₀ = 20 — this is the *original deep-quantum baseline from before the framework Q-band structure was identified*, NOT a framework anchor. The `wave_breaking_q_anchor_scan.py` script explicitly labels Q=20 as "Deep-quantum (J=1, original default)". When this default is used, the system rotates 20 times per γ₀-tick, far outside the [0.2, 2.0] framework-anchor range. Any analysis at J=1 should explicitly note it's outside the framework's Q-band structure.

The "20" in the C# headers (`t=20` in F73, `tMax=20` in cockpit, `tMax = 1/γ₀` in InspectCommand) refers to **one γ₀-period in absolute sim-units** — that's invariant under J choice. What is NOT invariant: the H-rotation count packed into that γ₀-period, which is Q.

The "20" shows up explicitly in three canonical C# sites:

1. **F73 verified value** ([`F73SpatialSumPurityClosurePi2Inheritance.cs:123-126`](../compute/RCPsiSquared.Core/Symmetry/F73SpatialSumPurityClosurePi2Inheritance.cs)):
   > N=5, γ₀=0.05, **t=20**: closure = (1/2)·exp(−4·0.05·20) = 9.157819·10⁻³
   The verification anchors at exactly 4 dyadic-ladder dimensionless units (a_{−1} = 4 × γ₀·t = 4×1 = 4).

2. **Propagate cockpit default** ([`compute/RCPsiSquared.Propagate/Program.cs:1118`](../compute/RCPsiSquared.Propagate/Program.cs)):
   ```csharp
   double gamma = 0.05, J = 1.0, tMax = 20.0, dt = 0.05, sample = 0.1;
   ```
   Default simulation runs for exactly one γ₀-tick.

3. **InspectCommand t-max** ([`compute/RCPsiSquared.Cli/Commands/InspectCommand.cs:134`](../compute/RCPsiSquared.Cli/Commands/InspectCommand.cs)):
   ```csharp
   // t_max = 4·t_peak = 1/γ₀, t-points = 41
   double tMax = p.OptionalDouble("t-max") ?? (1.0 / block.GammaZero);
   ```
   Inspect mode auto-scales tMax = 1/γ₀ (= 20 at γ₀=0.05).

The γ₀=0.05 convention is documented as substrate-invariant in [`UniversalCarrierClaim.cs:88`](../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs): "ValidateAgainstPythonStepFTests covers γ₀ ∈ {0.025, 0.05, 0.10} with identical Q-values; structural slot is substrate-invariant. IBM hardware (T2* ~ 100μs → γ ~ 10⁴ Hz) different physical value, identical role." 0.05 was chosen as a convenient round number, not a physical constant.

## Q-anchor hits on F99 ankers (from `wave_breaking_q_anchor_scan.py`)

At Q=1 (Balance), the chain Liouvillian Im=0 fraction lands directly on F99 ankers at three N-values:

| N | Im=0 count / total | Fraction | F99 anker          |
|---|--------------------|----------|--------------------|
| 2 | 6/16               | **3/8**  | KIntermediate      |
| 3 | 32/64              | **1/2**  | Generic            |
| 5 | 256/1024           | **1/4**  | Silver-Dicke       |

N=4 and N=6 at Q=1 land on **off-Anker** fractions (3/16 and 3/32). At higher Q (1.5, 2, 20), all small-N values collapse to the combinatorial baseline (3/8)^⌊N/2⌋. This is why Q=1 is structurally privileged in `FractionReferenceGraph` — three F99-anker-hitting edges from it.

## Anchors and cross-refs

- **Detailed F86 inventory**: [`F86_VALUES_INVENTORY.md`](F86_VALUES_INVENTORY.md) — per-c, per-orbit, per-bond-class typed values + open derivations.
- C# Q-anchor representation: `FractionReferenceGraph.QBasisAnkers = { 1.0, 1.5, 2.0 }` ([`compute/RCPsiSquared.Core/Symmetry/FractionReferenceGraph.cs`](../compute/RCPsiSquared.Core/Symmetry/FractionReferenceGraph.cs)). Currently captures 3 of the 10 anchors above; absent: onset-band (0.2, 0.35), peak-band start (1.2), F86 Q_peak c=3 (1.6), peak-band end / F86 Q_peak c=4,5 (1.8), Endpoint candidate (2.5).
- Q-band reading: memory `project_q_middle_structure`.
- Per-orbit Q_EP structure: memory `project_q_peak_ep_structure` + [`F86/C2BlockCpsiQScan.cs`](../compute/RCPsiSquared.Core/F86/C2BlockCpsiQScan.cs).
- γ₀ tick reading: memory `project_tick_and_angle` + [`reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md`](../reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md).
- Q as clock-exchange-rate (Tier 4 reading): [`hypotheses/Q_AS_THE_EXCHANGE_RATE.md`](../hypotheses/Q_AS_THE_EXCHANGE_RATE.md).
- Substrate invariance: [`UniversalCarrierClaim.cs`](../compute/RCPsiSquared.Core/Symmetry/UniversalCarrierClaim.cs).
- Wave-breaking Q-anker hits: [`simulations/wave_breaking_q_anchor_scan.py`](../simulations/wave_breaking_q_anchor_scan.py) + JSON results in [`simulations/results/wave_breaking_q_anchor_scan/`](../simulations/results/wave_breaking_q_anchor_scan/).
