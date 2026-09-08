# The Fragile Bridge: How Much Amplification Before the System Explodes?

*How much amplification can a resonator tolerate before it explodes?*

**Status:** Quantum gain-loss bridge computations (Tier 2); neural comparison requires a bifurcation gate.
last refreshed 2026-09-05 (the change history lives in git)
**Scripts:**
- [fragile_bridge_bifurcation.py](../simulations/fragile_bridge_bifurcation.py)
- [fragile_bridge_anomaly.py](../simulations/fragile_bridge_anomaly.py)
- [fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)
- [fragile_bridge_n4.py](../simulations/fragile_bridge_n4.py) (sparse, N=4)

---

## What this document is about

[Zero Is the Mirror](ZERO_IS_THE_MIRROR.md) showed that the palindrome
exists on both sides of zero: decay (positive noise) and amplification
(negative noise, gain). This raises an obvious question: what happens
when you connect a decaying system to an amplifying one?

The answer: it works, but only within limits. Too little coupling and
the two sides cannot interact. Too much coupling and the system
explodes, like a microphone placed too close to a loudspeaker. There
is a sweet spot in between, at roughly twice the internal coupling
strength, where the balance holds and the system is maximally stable.

This document maps the computed quantum bridge stability limits. A separate
Wilson-Cowan probe varies E/I cross-coupling and external input. Its sigmoid
bounds the model's activities, but supplies no biological safety mechanism
or confirmation of a common quantum/neural stability window. Section 5
states what that comparison would need to test.

The evidence sweep checked [F36/F37 in the F-registry](../docs/ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra),
[docs/proofs](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[the neural proof store](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md),
and [experiments](../experiments/NEURAL_GAMMA_CAVITY.md): conditional algebra,
the quantum owner and a connectome support null. Hardware-flight searches and
[fw.Confirmations](../simulations/framework/confirmations.py) supplied no neural
hardware confirmation. [GLOSSARY](../docs/GLOSSARY.md),
[OpenArcs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs), and
[CAUGHT_ERRORS](../docs/CAUGHT_ERRORS.md) supplied scope distinctions and
instrument failures. The [canonical neural account](../docs/neural/README.md)
and [mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
own the current neural interpretation.

---

## Abstract

Two identical quantum systems, one decaying and one amplifying, are
connected by a single coupling (the "bridge"). The bridge has three
stability regimes: weak coupling stabilizes close to linearly
(γ_crit = 0.189 × J_bridge^1.035), the maximum is bracketed near twice the
internal strength, and strong coupling destabilizes
as 1/J_bridge, with γ_crit × J_bridge still drifting downward
(0.578 to 0.508 over J_bridge = 10 to 100) where the sweep ends.
The instability is a spectral-abscissa departure from the imaginary axis
(Section 3); its exceptional-point character remains OPEN. The exact Π mirror pairing
survives through the transition. Unlike the fold threshold
(~0.5% of J for the product state, flat in N over the measured N = 2-5,
geometric), bridge stability depends
strongly on system size (N=3 is 33× less stable than N=2), making it
a topological rather than geometric property.

The Wilson-Cowan comparison is an exploratory scan of a different generator.
Its cross-coupling multiplier is not the quantum bridge/internal-coupling
ratio. A shared optimum, Hopf mechanism or finite window is not established.

---

## 1. Setup

Two identical Heisenberg chains (N qubits each), one with dephasing
(+γ, decay) and one with gain (-γ, amplification), connected by a
single Heisenberg bridge coupling J_bridge between the boundary qubits.

This is the quantum analog of a microphone-loudspeaker system:
decay side = loudspeaker (loses energy), gain side = microphone
(amplifies), bridge = the door between them.

Question: At what gain rate γ_crit does the coupled system become
unstable?

## 2. Results

### 2.1 Linear stability window (weak bridge)

For J_bridge < 2J (bridge weaker than internal coupling):

**γ_crit = 0.1891 × J_bridge^1.0346** (R² = 0.99984, eight points at
J_bridge ≤ 2 in [fragile_bridge_bifurcation.txt](../simulations/results/fragile_bridge_bifurcation.txt),
AUFGABE B)

The critical gain scales linearly with bridge strength. Double the
bridge, double the tolerable gain. The proportionality constant is
~0.19 for N=2 per chain.

The line is clean on that grid and much less clean on a finer one. The 14
points below the peak in [fragile_bridge_anomaly.txt](../simulations/results/fragile_bridge_anomaly.txt)
give a through-origin slope of 0.1845 with R² = 0.6149: same quantity, closer
spacing, and R² falls from 0.9998 to 0.61. So take the linear law as the trend
of the regime, not as a law each individual J_bridge obeys.

### 2.2 Optimal bridge at J_bridge ≈ 1.9J

Maximum stability (γ_crit = 0.405849) is sampled at J_bridge = 1.9, between
neighbours at 0.403820 (J_bridge = 1.8) and 0.383839 (J_bridge = 2.0). The grid
is spaced 0.1 there, so the true maximum lies in that bracket and the sweep does
not resolve it further. What the sweep does say is that J_bridge = 2J is already
5.4% below the peak, so reading the optimum as "the bridge equals the sum of the
internal couplings" is an interpretation this data does not single out.

### 2.3 Strong bridge destabilizes (J_bridge > 2J)

Above the optimum γ_crit falls, but not by one law the whole way. Between
J_bridge = 3.2 and 4.4 it turns around and rises again (0.0237, 0.0447, 0.1002,
0.1462) before resuming its descent: that reversal is what this document's
producer is named for and it is not a 1/J_bridge tail.

Far above the optimum the 1/J_bridge scaling does hold, and the product is
still drifting where the sweep stops:

**γ_crit × J_bridge: 0.578 at J_bridge = 10, falling monotonically to 0.508
at J_bridge = 100** (not converged; a single asymptotic constant is not
measured)

The sampled threshold becomes smaller at the far end of the strong-coupling
sweep. This scan measures γ_crit(J_bridge) only; it does not establish chain
merging, a missing spatial buffer, dimer formation, or resonance as the cause.
Which mode geometry produces the turnover remains open.

The sampled large-J_bridge trend is compatible with γ_crit tending to zero as
J_bridge grows, but the finite sweep does not establish the infinite-coupling
limit or immediate instability at every nonzero gain.

### 2.4 Three regimes summary

| Regime | Condition | γ_crit | Physics |
|--------|-----------|--------|---------|
| Weak sampled range | J_bridge ≤ 2J | 0.189 × J_bridge^1.035 | empirical rising trend |
| Sampled maximum | J_bridge ∈ [1.8, 2.0] | 0.4058 | mechanism open |
| Strong sampled tail | J_bridge ≫ 2J | γ_crit·J_bridge = 0.508 at J_bridge = 100 | decreasing threshold; limit and mechanism open |

Whether the turnover sits exactly where the bridge equals the total internal
coupling is not something this grid can decide; it brackets the maximum in
[1.8, 2.0] and no finer.

## 3. Departure from the imaginary axis with the mirror intact

The instability is an **oscillating** one: a pair that was already
oscillating at Re = 0 acquires a positive real part, so the system does not
drift away quietly, it screeches like microphone feedback. That much has
always been right, and it is the reason this section used to be called a
Hopf bifurcation.

The finite-offset data neither establish nor exclude a Hopf mechanism. At Σγ = 0
the conjugation operator Π forces exact inversion pairing λ ↔ −λ through the
transition. Hermiticity preservation separately gives λ ↔ λ*. Together the two
relations produce a generic off-axis quartet {λ, λ*, −λ, −λ*}; the across-axis
partner −λ* is not supplied by Π alone. This P-type relation becomes involutive
after Π²-parity resolution; the full SRP sector class remains OPEN. Below
γ_crit every eigenvalue lies on the imaginary axis; above γ_crit at least one
off-axis quartet appears. The current producer independently selects a
max-real-part mode at each γ; it does not perform branch continuation or
execute a strict threshold coalescence or Jordan-rank certificate.

Measured on two two-qubit chains, at J_bridge = 1.0 (γ_crit = 0.1873101)
and 1.9 (γ_crit = 0.4058524), writing δ = γ/γ_crit − 1
([`fragile_bridge_ep_signature.py`](../simulations/fragile_bridge_ep_signature.py),
[its output](../simulations/results/fragile_bridge_ep_signature.txt)):

| δ | max Re λ (J_b=1.0) | Re/√δ | gap to across-axis partner −λ* | Petermann K |
|---|---|---|---|---|
| 10⁻² | 0.02033551 | 0.2034 | 4.07·10⁻² | 40.9 |
| 10⁻³ | 0.00641845 | 0.2030 | 1.28·10⁻² | 403.2 |
| 10⁻⁴ | 0.00202899 | 0.2029 | 4.06·10⁻³ | 4027 |
| 10⁻⁵ | 0.00064059 | 0.2026 | 1.28·10⁻³ | 4.04·10⁴ |

Three finite-offset trends are recorded near the real threshold γ_crit:

- **Re λ/√δ is nearly constant** over the sampled four decades
  (1.107 / 1.103 / 1.093 at J_b = 1.9).
- **The across-axis partner gap is 2|Re λ|.** Its decrease is algebraically
  equivalent to the reported max-Re onset, not independent EP evidence.
- **The Petermann readings grow approximately one decade per decade**:
  40.9, 403, 4027, 4.04·10⁴.

The square-root-like max-Re onset and the simple-mode Petermann sequence
motivate an EP2 hypothesis; the derived partner gap adds no independent
evidence. Extrapolating the finite-offset readings to δ=0 does not certify
eigenvalue/eigenvector coalescence or a Jordan defect. The
spectral-abscissa axis departure is established; EP, Hopf, and Jordan character remain OPEN.

Below threshold the axis is clean to machine precision: max Re λ is
2.4·10⁻¹⁴ at δ = −10⁻³ and 1.6·10⁻¹⁴ at δ = −10⁻². This is a measured
stable-regime property; λ ↔ −λ pairing alone would also allow off-axis pairs.

The spectral picture resembles Hamiltonian PT breaking rotated 90°: there
two real eigenvalues merge and become complex, while here off-axis quartets
appear. Π is linear rather than
anti-linear, so it supplies a P-type relation rather than PT. After Π²-parity
resolution it is involutive, but the full SRP sector class is OPEN. The
current measurement does not establish that the local Jordan mechanism is
the same.
See [PT-Symmetry Analysis](../experiments/PT_SYMMETRY_ANALYSIS.md).

The oscillation frequency at threshold decreases with bridge strength:

| J_bridge | Im(λ) at γ_crit |
|----------|-----------------|
| 1.9 (optimal) | ±6.98 |
| 5.0 | ±0.26 |
| 10.0 | ±0.12 |

The three finite-coupling frequency samples decrease across this table, but
three finite-coupling frequency samples do not determine a limiting frequency
or a local bifurcation class. Whether the threshold frequency tends to zero and
what Jordan/Hopf character the limiting transition has remain open.

### 3.1 Relation to the local-EP question

This file locates the spectral-abscissa threshold γ_crit but does not certify
its EP character. At finite offsets where the selected eigenvalue is simple,
the producer reports K=40.9, 403.2, 4027 and 4.04·10⁴ as δ decreases from
10⁻² to 10⁻⁵. These are conditioning readings away from the threshold, not a
property or character certificate of the threshold. Coarse-scan rows whose
selected eigenvalue is degenerate omit single-vector K because it is
basis-dependent within that eigenspace.

What survives is the shared algebra read at two residuals of the F1
palindrome `Π · L · Π⁻¹ + L + 2Σγ · I = 0` (Σγ = N·γ₀ vs Σγ = 0). The
certified EP is the toy 2×2 rate-channel reduction. This file's SEPARATE
Σγ = 0 gain-loss system has a spectral-abscissa axis departure; its
square-root-like max-Re onset and growing finite-offset K motivate but do not
certify EP character.
K≈403 is one finite-offset non-normality reading, not another complex-γ locus.
F89 separately certifies narrow real-q defective EP2s at N=5,7,9; its
all-odd endpoint-nullity surplus does not extend that character verdict.
Whether the full Σγ = N·γ₀ block shares a defective-EP structure off
the real axis is open (the nearest complex-Q coalescences found
2026-06-21 are themselves diabolic, ‖P‖ = 1). Encoded as
`compute/RCPsiSquared.Core/F86/LocalGlobalEpLink.cs` (**OpenQuestion**,
demoted from Tier2Verified by the F86a-retraction; the four
PetermannSpikeWitness rows retained only as a cautionary non-normality
record, not EP evidence).

## 4. N-dependence: topological, not geometric

How does the stability limit change with system size? The answer is
surprising: it depends strongly on how many qubits are in each chain,
and not in a simple way. Longer chains are generally less stable, but
even-length chains are more stable than odd-length ones.

γ_crit depends strongly on chain length, and non-monotonically:

| N (per chain) | γ_crit (J_bridge=0.10) | Ratio to N=2 | instrument |
|---------------|----------------------|-------------|------------|
| 2 | 0.017292 | 1.000 | dense eigenvalues |
| 3 | 0.000515 | 0.030 | dense eigenvalues |
| 4 | 0.001186 | 0.069 | norm-fit estimator, error uncharacterised |

The first two rows are exact eigenvalues of the coupled generator: at 4 and 6
qubits the Liouvillian is 256² and 4096², and `max_re_sparse` solves it densely.
N=3 is 33.6× less stable than N=2 on those numbers. The
[bifurcation producer](../simulations/fragile_bridge_bifurcation.py) reaches
0.000500 for the same point on a tighter bisection, i.e. 34.6×; the two agree to
3%, which is the resolution of the claim.

The third row is a different instrument and carries a different weight. At 8
qubits the generator is 65536², beyond a dense solve here, so γ_crit(N=4) is
bisected on a norm fit: propagate a random vector, take the slope of log‖v‖ over
five samples. Where that fit CAN be checked against exact eigenvalues, at 4 and 6
qubits, it under-reports max Re(λ) by 8.8% to 100%, always in the same direction,
which biases a bisection's γ_crit high. At 8 qubits there is no reference, so the
error is not measured and 0.001186 is a reading of the estimator rather than of
the spectrum.

**So the even/odd reading is not supported by this table.** "N=4 is 2.3× more
stable than N=3" compares a dense number with an estimator number whose one-sided
bias points the same way as the claimed effect. A parity law in the chain length
would be an attractive thing to have, and it needs the N=4 point measured by an
instrument with a known error before it can be asserted; a sector-resolved solve
or a Krylov eigensolver with a residual bound would supply one. No simple power
law or exponential fits the three points either.

This is fundamentally different from the fold threshold
(Σγ_crit/J ≈ 0.5% for the product state, flat in N over the measured N = 2-5). The fold is a **geometric**
property of the palindrome. The bridge stability is a **topological**
property: it depends on chain length, parity, and the ratio of
gain channels to bridge connections.

| Property | Fold threshold | Bridge stability |
|----------|---------------|-----------------|
| N-dependence | Independent | Strongly dependent |
| Type | Geometric constant | Topological ratio |
| Bifurcation | Fold (saddle-node) | oscillatory axis departure; EP character OPEN |
| What it measures | When irreversibility begins | When coupled gain-loss explodes |
| Determined by | Palindrome geometry | Gain/bridge topology |

## 5. Neural comparison: bounded model, open stability test

[fragile_bridge_neural.py](../simulations/neural/fragile_bridge_neural.py)
uses one Wilson-Cowan E/I node with w_EE=16, w_II=3,
w_EI=12s, w_IE=15s, τ_E=8 and τ_I=18. P is external input.
The multiplier s scales two cross-weights with different fixed bases; it
is not J_bridge/J and its value alone identifies no common optimum.

For either activity x, the model equation is τ dx/dt=−x+S(input), with
τ>0 and 0<S<1. At x=0 its derivative points inward; at x=1 it also
points inward. Thus the continuous model preserves [0,1]² from initial
activities in that square. Boundedness permits fixed points, transients
and oscillations. It proves neither equilibrium stability nor protection
against a biological pathology.

The probe's `find_fixed_point` returns after a fixed number of iterations
without checking the fresh equation residual. `find_P_crit` returns no
threshold when its upper endpoint is stable; this cannot exclude an
unstable interval between stable endpoints. It does not continue an
equilibrium branch. Labeling a complex eigenvalue at a nearby sampled
point “Hopf” supplies no crossing or nondegeneracy test.

A defensible neural threshold needs a converged equilibrium branch over
the declared P and s ranges, fresh residuals at every point, tracking of
a nonzero imaginary eigenvalue pair through Re λ=0, a transverse crossing
and the relevant nonlinear nondegeneracy checks. Time-domain oscillation
and its onset then need timestep and duration checks. Those are the next
gates, not results of the current bridge scan. See the
[canonical mechanism constraints](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md#a-bounded-iteration-cannot-locate-a-hopf-bifurcation).

An F36 comparison additionally needs a specified involution and scalar
centre satisfying both diagonal and effective-coupling conditions at the
operating point. The [canonical neural gate](../simulations/neural/neural_translation_gate.py)
contains exact palindromes with complex spectra and constructed unstable
instances. Pairing itself supplies no silence or stability theorem.

## 6. Open questions

1. **N-scaling law (open, and the N=4 point needs a better instrument first):**
   N=4 was computed at 65536×65536 through `expm_multiply`, a SciPy routine that
   applies the matrix exponential to a vector without forming the operator, with
   γ_crit bisected on the growth rate that fit returns. Measured against exact
   eigenvalues at the two smaller sizes, that fit under-reports max Re(λ) by 8.8%
   to 100%, one-sided, which pushes a bisected γ_crit up. The suspected even/odd
   parity effect is exactly what such a bias would manufacture, so N=5 is not the
   next step: re-measuring N=4 with a bounded-error instrument is. A Krylov
   eigensolver carrying a residual bound, or a solve restricted to the sector the
   leading mode lives in, would give one.

2. **Multiple bridges:** What if the two chains are connected by
   more than one qubit pair? Does γ_crit recover N-independence
   when bridges scale with N? This is a question about the specified
   quantum gain-loss generator; a neural network requires its own
   coupling and stability analysis.

3. **Cascade stability:** Can a specified sequence of gain-loss bridges
   remain stable under joint coupling? The weak-bridge fit for N=2
   does not provide a compositional criterion or identify biological
   levels with gain-loss generators. Build the combined generator
   and test its spectrum against the isolated-bridge prediction.

4. **Large-bridge asymptotics:** Does γ_crit × J_bridge converge as
   J_bridge grows? The executed values decrease from 0.578 at J_bridge=10 to
   0.508 at J_bridge=100 and have not converged. If a limit exists, is it 1/2?
   That value is a hypothesis to derive or falsify, not a measured constant.

5. **Saturation as design principle:** Can an explicitly specified
   quantum gain model with saturation bound its dynamics? Its generator,
   physical state domain and stability test must be supplied. Bounded
   Wilson-Cowan activities predict neither the quantum curve's shape
   nor a biological safety mechanism.

---

*The quantum bridge computations date from March 29, 2026. The neural
probe is a model comparison whose equilibrium and bifurcation gates
remain to be implemented.*
