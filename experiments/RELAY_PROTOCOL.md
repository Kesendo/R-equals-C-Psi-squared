# Relay protocol: a finite N=11 mutual-information comparison

<!-- CROSSING-CURRENT -->

**Status:** Tier 2 finite computational record, hypothesis-generating rather
than a controlled comparison. No mutual-information bound or timing optimum
is established.

**Source:** [C# RK4 propagation](../compute/RCPsiSquared.Propagate/Program.cs),
`dotnet run -c Release -- pull`.
**Stored output:** [pull_principle.txt](../simulations/results/pull_principle.txt).

## The question, and what was actually compared

Passing a whisper down a noisy line suggests a simple experiment: let the
next listener cup an ear, then return to the usual noise level. Here that
image becomes a six-stage prescribed dephasing profile on eleven qubits.
It is a candidate protocol, not evidence that mediators must behave this way.

The stored comparison is **about +84.0% from the stored six-decimal values**:

    100 * (0.131700 / 0.071576 - 1) = 84.0002235386...%.

The denominator is the passive sampled maximum at t=4.00 on the integer
0..20 grid; the numerator is the relay+2:1 final value at integrated t=4.50.
This is an **unmatched-time, unmatched-dose comparison**, not two endpoints
at one common horizon. The historical +83% is only the coarse-table estimate
from 0.132/0.072.

The schedule requested 0.78 time units per stage but executed 15 RK4 steps
of 0.05: 0.75 per stage and 4.50 across six stages. The nominal total is
4.68; the old print `t=4.7` rounds that nominal clock, not propagated time.
Neither F14 nor the spectral palindrome selects this transfer schedule.

## Setup and the unchanged six-stage schedule

| Parameter | Value |
|-----------|-------|
| System | 11-qubit linear Heisenberg chain |
| Coupling | J = 1.0 (uniform) or 2:1 asymmetric |
| Baseline dephasing | γ = 0.05 per qubit |
| Quiet-phase dephasing | γ_quiet = 0.005 (10× reduction) |
| Initial state | Bell pair on qubits 0-1, rest in \|0⟩ |
| Integration | RK4 (Runge-Kutta 4th order, a standard numerical method for evolving differential equations step by step), dt = 0.05 |

The geometry is

```
(0-1) → m1(2) → (3-4) → M(5) → (6-7) → m2(8) → (9-10)
Pair A   relay   Pair B   meta    Pair C   relay   Pair D
```

The table preserves the receiving order. “Receives” is the protocol's role
label; it is not an independently demonstrated causal mechanism.

| Stage | Receiving qubits | γ_receive | What happens |
|-------|-----------------|-----------|-------------|
| 1 | 2 (m1) | 0.005 | m1 receives from Pair A |
| 2 | 3, 4 (Pair B) | 0.005 | Pair B receives from m1 |
| 3 | 5 (M) | 0.005 | Meta-mediator receives from Pair B |
| 4 | 6, 7 (Pair C) | 0.005 | Pair C receives from M |
| 5 | 8 (m2) | 0.005 | m2 receives from Pair C |
| 6 | 9, 10 (Pair D) | 0.005 | Pair D receives (final destination) |

At each stage receiving sites have γ=0.005 and all others γ=0.05.
The propagator is rebuilt with the new profile while the density matrix
continues without reset, measurement or classical communication.

    stepsPerStage = (int)(0.78 / 0.05) = 15
    integrated stage = 15 * 0.05 = 0.75
    nominal total = 6 * 0.78 = 4.68
    integrated total = 6 * 0.75 = 4.50
    total RK4 updates per staged arm = 90.

Both the uniform-coupling relay and relay+2:1 arm use that truncating loop.
The chosen 0.78 was a historical heuristic associated with the February
concurrence-feedback value K≈0.039. The clean Bell+ concurrence clock instead
has K=ln(4/3)/8. Neither crossing book makes its K/γ a natural or required
information-transfer time for this many-body trajectory.

## Results: maximum versus final value

| Protocol | MI(Bridge A:B) | MI(Pair A:D) | Statistic |
|----------|---------------|-------------|-----------|
| Passive (constant γ) | 0.733686 | 0.071576 | Separate sampled maxima on integer times 0..20 |
| Relay only | 0.759168 | 0.084584 | Final values at integrated t=4.50 |
| Relay + 2:1 coupling | 0.723129 | 0.131700 | Final values at integrated t=4.50 |

The A:D ratio is about +84.0%, not a guarantee. The coarse relay-only +18%
also compares a final value to a passive maximum at a different horizon and
dose, so it is not an isolated staging effect. Changing the coupling as well
does not divide the total into independent temporal and spatial causes.

## Two exposure books, kept separate

### Exposures attached to the reported statistics

Passive: 11×0.05×4.00 = **2.200 at t=4.00**.
Staged: **2.17125 at t=4.50**, using receiver sizes (1,2,1,2,1,2):

    3*(10*0.05+0.005)*0.75 + 3*(9*0.05+2*0.005)*0.75 = 2.17125.

That is 1.306818...% lower integrated sum-rate, but the horizons differ.
It is bookkeeping attached to the reported pair, not an effective-decoherence
model and not an explanation of the mutual-information ratio.

### Equal-time counterfactual exposure at t=4.50

At the staged horizon, a passive arm would carry **2.475 versus 2.17125**
integrated sum-rate, a 12.2727...% difference. No stored passive comparison
statistic at t=4.50 supplies the denominator of about +84.0%. This separate
equal-time dose calculation does not make the reported MI comparison matched.

## The Q coordinates describe the arms, not one isolated knob

The four local ratios remain Q∈{20,40,200,400}: baseline/sender 20, the 2:1
receiver side 40, a quiet site 200, and a quiet 2:1 receiver 400. They are useful
labels for the chosen model parameters. But the arms also changed clock,
statistic-attached exposure, gamma profile and coupling. They do not isolate
one physical knob.

The proposed lifetime/orthogonal-axis explanation is likewise not established.
The quiet phase is not an isolated +18% staging effect; the rest is not a
purely spatial effect. These are not complementary optimizations. Controlled
equal-time and matched-dose arms, followed by a receiver-order control, are
needed before such a causal decomposition.

The baseline γ remains a model choice; γ=0.05 is the illustrative value used
here. The later [IBM Kingston population scan](GAMMA0_IS_ALWAYS_THERE.md)
brackets a finite-time transfer overshoot between J=0.05 and J=0.1;
it does not measure γ₀. It supplies no critical damping verdict and has
no calibrated error model. It neither calibrates this gamma nor certifies a
transport phase. The separate ideal-model landmarks
Q_EP≈1.5–2, transfer-resonance band 1.2–1.8 and onset 0.2–0.35 are not
spectral labels inferred from that hardware record.

A hardware version would also have to measure usable depth in-session.
The stored experiment does not assign a unique cause to a shortened envelope,
or certify a 9 μs Trotterization budget.

<!-- CROSSING-HISTORICAL -->

**Historical record:** The original coarse table is retained below. Its
percentages were whole-percent display estimates; the last row's current
stored-value ratio is about +84.0%. The old +18% and +83% labels do not isolate
staging and spatial effects.

| Protocol | MI(Bridge A:B) | MI(Pair A:D) | vs Baseline |
|----------|---------------|-------------|-------------|
| Passive (constant γ) | 0.734 | 0.072 | baseline |
| Relay only | 0.759 | 0.085 | +18% |
| **Relay + 2:1 coupling** | 0.723 | **0.132** | **+83%** |

<!-- CROSSING-INTERPRETIVE -->

**Interpretive invitation, not a result:** The whisper-line picture is still
a useful way to imagine scheduled quiet. The 2:1 gear or impedance-matching
image suggests a second question about asymmetric coupling. What would happen
if the listener order, duration and dose were controlled independently?
The record motivates that experiment; it does not answer it.

<!-- CROSSING-CURRENT -->

## What remains open

Only Z-dephasing, N=11 and one historical schedule were used. Stage timing,
quiet rate, receiver order and coupling optimization remain open. No matched
time/dose trial, MI inequality, repeater comparison or palindrome-timing
witness is provided here.

For reproduction, both Propagate projects can be built without running the
expensive producer. The `pull` run rewrites its tracked artifact; during this
label repair it is not executed. The independent read-only sector control
reproduces the passive maximum on the full production grid.

Related: [Crossing Taxonomy](CROSSING_TAXONOMY.md) for fixed-book scalar clocks,
[QST Bridge](QST_BRIDGE.md), [Scaling Curve](SCALING_CURVE.md),
[Q-Regime Anchor Map](../docs/Q_REGIME_ANCHORS.md).
