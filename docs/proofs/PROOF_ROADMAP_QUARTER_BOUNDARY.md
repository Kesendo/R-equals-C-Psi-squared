<!-- QUARTER-CURRENT -->
# The Journey to Locate the Algebraic Quarter

Current reading: this roadmap follows the algebraic discriminant quarter into
named quantum trajectories while keeping proof, finite computation, and open
physical interpretation as separate layers.

## A Proof Roadmap for the R = CΨ² Critical Boundary

*Working document, begun March 2026*
*Guiding principle: keep the algebraic quarter distinct from the dynamics that may cross, re-cross, or avoid it.*

---

## Preface: What This Document Is

This is not a finished paper. It is a map of an argument assembled layer by layer, from a chosen scalar normal form to named quantum trajectories. For each layer, we state clearly what is proven, what is computationally verified, what is conjectured, and what remains to be done.

The central claim: the self-referential fixed-point equation

$$R_{n+1} = C(\Psi + R_n)^2$$

has a discriminant boundary at $C\Psi = 1/4$. That coordinate is exact for this chosen normalized quadratic normal form; a reparameterization can move its numerical label without changing the fold type. The same normal form maps to $z_{n+1}=z_n^2+c$ and hence to the real cusp of its main cardioid. IBM hardware supplied three named sightings rather than a theorem over channels: the [first recorded crossing](../../experiments/IBM_QUANTUM_TOMOGRAPHY.md) (ibm_torino q52, February 2026), the [tightest single-point crossing](../../experiments/IBM_RUN3_PALINDROME.md) at 1.9% deviation (q80, March 2026), and the F25 purity-book trajectory $C\Psi(t)=f(1+f^2)/6$ fitted through the quarter with RMS residual 0.0097 ([ibm_kingston, April 2026](../../data/ibm_cusp_precision_april2026/README.md)).

**The two variables, and two books.** C and Ψ are opaque scalars in the
algebra. In the F25 purity book, `C=Tr(ρ²)` and
`Ψ=ℓ₁/(d−1)`, giving `CΨ=f(1+f²)/6` for Bell+/local-Z/no-H. In the separate
concurrence book, `C` is Wootters concurrence and the same named trajectory
gives `CΨ=f²/3`. Both begin at 1/3, but they are different readings and their
later values must not be interchanged. The carrier is the two-qubit pair;
its one-qubit marginal is `I/2` and blind to this pair trajectory. The
canonical definitions and their history live in
[The CΨ Lens](../THE_CPSI_LENS.md) and the [Glossary](../GLOSSARY.md).

The proof journey works upward:

1. Single qubit (d = 2): the algebraic foundation
2. Two entangled qubits: named trajectories and conditional subsystem crossing
3. N-qubit systems: GHZ, W, and the palindromic structure
4. Arbitrary dimension d: qutrits and beyond
5. Channel dependence: named exact families, counterexamples, and an open peak question
6. The assumed-family classification: why this normal form reports 1/4
7. Connections to known mathematics: Mandelbrot, Feigenbaum, and deeper structures

## Summary: What Is Proven

| Layer | Status | Key Result |
|-------|--------|------------|
| 1. Qubit (d=2) | PROVEN | Discriminant of R=C(Psi+R)^2 vanishes at CPsi=1/4. Crossing cubic. Mandelbrot identity. |
| 2. Two qubits | CONDITIONAL | If a continuous trajectory converges to ρ* with CΨ(ρ*)<1/4, it eventually stays below. A primitive CPTP target at CΨ=0.2935 and a local Markov upward crossing forbid the old absorber reading. |
| 3. N-qubit | PROVEN | Palindromic spectrum all graphs N=2..8 (87,376 eigenvalues). Analytic formula. |
| 4. Dimension | ANSWERED: d=2 only | Qutrits: 0/236 dissipators palindromic. Discriminant d-independent. CV/hybrid extensions (Conj 4.1/4.2) open. |
| 5. Channels | NAMED EXACT + FINITE ATLAS + OPEN | F25–F27 give named Bell+ formulas. Universal pointwise/absorber/local-control claims are false. The autonomous N=2 successive-peak claim is unproved. Finite N=3/4/5 rows invite an all-Q/all-N classification and mechanism. |
| 6. Assumed-family uniqueness | CONDITIONAL ALGEBRA | Within `R=C_α(Ψ+R)^α`, alpha=2 alone removes the explicit Ψ factor from the fold product. This does not derive the family or its physical use. |
| 7. Math connections | MIXED | Mandelbrot identity exact (PROVEN). Feigenbaum cascade numerically measured (resolution-limited). A named N=2 Bell+/Z sampled Bures path coefficient gives g≈3.36, with no divergence resolved on its grid; ambient density-matrix geometry remains unclassified. Holography SPECULATIVE. |

Core closed: Layer 1 locates 1/4 for the stated recurrence, Layer 3 proves the palindrome in its stated scope, and Layer 7 gives the Mandelbrot change of variables for that normal form. Scoped/partial: Layer 2 is the conditional convergence implication; Layer 4 is d=2 only; Layer 5 has named exact channel formulas, false universal dynamics claims, and an open peak classification; Layer 6 is an assumed-family comparison, not a physical derivation; Layer 7 holography remains open. The IBM runs are named crossings and an F25 trajectory fit, not a hardware proof of a universal absorber. Across four machines, finite calibration windows place many qubits on both sides of 1/4 ([clock field](../../experiments/CLOCK_FIELD_SITE_OWNED.md)); that richer motion is precisely why algebra and dynamics must not share a label.

---

## Layer 1: The Single Qubit (d = 2)

### The Algebraic Foundation

This is the bedrock. Everything else is built on what happens in a single two-level system.

### What Is ALREADY PROVEN

**The chosen product-power form.** For the posited recurrence, the combination $C\Psi^2$ (i.e., $a=1,b=2$):

- Produces a quadratic fixed-point equation
- Maps to the Mandelbrot iteration $z_{n+1} = z_n^2 + c$ under the substitution $z = C(\Psi + R)$, $c = C\Psi$
- Has a discriminant with a clean critical value

The fixed-point equation $R = C(\Psi + R)^2$ expands to:

$$CR^2 + (2C\Psi - 1)R + C\Psi^2 = 0$$

The discriminant is:

$$D = (2C\Psi - 1)^2 - 4C^2\Psi^2 = 1 - 4C\Psi$$

This vanishes at $C\Psi=1/4$, giving one double real algebraic root. Below it the polynomial has two real algebraic roots and above it none. Whether either root lies in a chosen physical interval is a separate check; the discriminant alone says nothing about a quantum trajectory.

**The crossing cubic.** At the critical boundary $C\Psi = 1/4$, with the normalized l1-coherence and the correlation bridge definition, the boundary condition reduces to the cubic:

$$b^3 + b = \frac{1}{2}$$

where $b$ is a normalized bridge parameter (not the $b = 1/2$ fixed-point
real part of F95 below; two different uses of the letter). This cubic has
exactly one real root ($b \approx 0.4239$), which fixes the crossing geometry
uniquely. The cubic has no free parameters; it is a pure number, independent
of any physical constants. Algebraically it is
[D4](../ANALYTICAL_FORMULAS.md#d4)'s single-qubit crossing condition
$f^*(1+f^{*2}) = 1/2$, the $d = 2$ member of the dimensional family
$(d-1)/2$: Bell+ has $3/2$, GHZ$_N$ has $(2^N-1)/2$
([Decoherence Relativity](../../experiments/DECOHERENCE_RELATIVITY.md),
[Coherence Density](../../experiments/COHERENCE_DENSITY.md)).

**The Mandelbrot correspondence.** The main cardioid of the Mandelbrot set is the set of $c$ values for which $z_{n+1} = z_n^2 + c$ has an attracting fixed point. On the real axis, the cardioid boundary is at $c = 1/4$. Our mapping sends $C\Psi \mapsto c$, so our critical boundary $C\Psi = 1/4$ sits exactly on the Mandelbrot cardioid cusp. This is not analogy. It is identity.

**The angle of the complex fixed point above the boundary ([F95](../ANALYTICAL_FORMULAS.md#f95), Tier 1 derived, May 16, 2026).** Once $C\Psi$ crosses $1/4$ from below, the discriminant $D = 1 - 4C\Psi$ goes negative and the two real fixed points become a complex conjugate pair $z_\pm = b \pm i\sqrt{c - b^2}$ with $b = 1/2$, $c = C\Psi$. The argument of the complex fixed point has the universal closed form

$$\theta(c; b) = \arctan\sqrt{c/b^2 - 1}$$

valid for a monic quadratic $z^2 - 2bz + c = 0$ with finite $b>0$ in the discriminant-negative regime. With the framework's $b = 1/2$ (`HalfAsStructuralFixedPointClaim`) and threshold $b^2 = 1/4$ (`QuarterAsBilinearMaxvalClaim`), this collapses to $\theta(c) = \arctan\sqrt{4c - 1}$: exactly the $\theta$-compass introduced state-specifically in [`BOUNDARY_NAVIGATION.md`](../../experiments/BOUNDARY_NAVIGATION.md) (Feb 8, 2026), now promoted to the positive-b polynomial-foundation identity used here. Derivation is 4 lines; numerical verification against the February table agrees within machine precision. See [`PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md`](PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md) and the companion reflection [`ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md`](../../reflections/ON_HOW_THE_ANGLE_APPEARS_AT_ZERO.md). This is the *angle-side* closed form of the discriminant-zero crossing; the *magnitude-side* closed form for per-outcome Born deviations in specific setups is [F94 = $(4/3) \cdot Q^2 \cdot K^3$](../ANALYTICAL_FORMULAS.md#f94).

### What Is HARDWARE-CONFIRMED

**The crossing, three times over.**

- **First crossing** (ibm_torino q52, February 9, 2026): single-qubit state
  tomography saw the product C·Ψ cross the ¼ boundary during decoherence, the
  first time on hardware. See [the tomography record](../../experiments/IBM_QUANTUM_TOMOGRAPHY.md).
  The separate [Q52 residual record](../../experiments/FIXED_POINT_SHADOW.md)
  retains 17 directional late-time rows and zero exceedances in 10,000 draws
  of a recorded null described as exponential decay, binomial shot sampling,
  and one random phase per synthetic run. It does not compare Q52-fitted
  time-dependent detuning/drift or the other hardware alternatives. Its
  positive tail slope is cut-sensitive, its boundary-distance correlation reuses |ρ₀₁| through
  CΨ, and the algebraic fixed-point phase is not a dynamics witness. The
  Q80/Q102 comparison rejects the former universal-boundary reading but
  leaves the Q52 magnitude-excess mechanism open.
- **[Tightest single-point crossing](../../experiments/IBM_RUN3_PALINDROME.md)**
  (ibm_torino q80, March 18, 2026): 1.9% deviation, the measured crossing at
  t\* = 15.29 μs vs the predicted 15.01 μs, i.e. matched to within 0.28 μs.
  Torino is a Heron r1 processor (133 qubits).
  The 1.9% is within expected hardware systematics (T1/T2 calibration
  drift, readout assignment error, crosstalk).
- **The full trajectory** ([ibm_kingston, April 2026](../../data/ibm_cusp_precision_april2026/README.md)):
  the F25 closed form CΨ(t) = f(1+f²)/6, f = e^(−4γt), fitted point-by-point
  through the boundary with γ the only free parameter, RMS residual 0.0097.
  This is a trajectory-level fit in the named Bell+/Z model; its final-crossing dose is
  K_Z = γ·t_cross = 0.03735 (F25's dose; the K_death/K_Z = 61.65 ratio
  lives in [F55](../ANALYTICAL_FORMULAS.md#f55)), from
  f\*(1+f\*²) = 3/2 at f\* = 0.8612.

**The boundary is a field, not a rim** ([clock field](../../experiments/CLOCK_FIELD_SITE_OWNED.md), July 2026).
Across four IBM machines (Torino, Marrakesh, Kingston, Fez), 61-87% of all
qubits sit on both sides of the CΨ = ¼ boundary over a calibration window
(e.g. Torino: 110 of 133). The ¼ boundary is not where a few special qubits
live; it is a field the whole chip moves through.

**The first Lindbladian spectrum (Bell+, two qubits; the layer's dynamics benchmark, one system ahead of its algebra).** The superoperator spectrum for Bell+ under Heisenberg coupling with local σ_z dephasing (γ = 0.1, J = 1) shows:

- Spectral gap = 0.2, relaxation time τ = 5.0
- 3 zero eigenvalues (1 from trace preservation + 2 degenerate steady states)
- All Re(λ) ≤ 0 (physical)
- Oscillatory eigenvalues at Im(λ) = ±4.0, confirming coherent-incoherent competition

**The assumed power-family comparison** (March 22, 2026).
Within `R=C_α(Ψ+R)^α`, α=2 alone makes the displayed fold product
independent of an explicit Ψ factor (Layer 6). This family-internal result
does not select that family from physics. See
[k_scaling_and_renyi.py](../../simulations/k_scaling_and_renyi.py).

### What Is CONJECTURED

- The crossing cubic $b^3 + b = 1/2$ may have number-theoretic significance beyond its role here (its real root is expressible in radicals via Cardano, but the closed form may connect to other mathematical constants)

---

## Layer 2: Two Qubits with Entanglement

### Subsystem Crossing Analysis

An entangled two-qubit pair lives in a four-dimensional Hilbert space (`d=4`).
For Bell+, that whole pair carries the CΨ reading; tracing either member leaves
the one-qubit state `I/2`, which is blind. In a larger system, a pair reading is
obtained by tracing spectators while retaining both carrier sites.

### What Is ALREADY PROVEN

**Bell state initial conditions.** For Bell+ ($|\Phi^+\rangle = (|00\rangle + |11\rangle)/\sqrt{2}$), the one-qubit reduced state is maximally mixed ($\rho_A=I/2$, blind to everything); the object that carries CΨ is the two-qubit pair (`d=4`). Its l1-coherence is $\ell_1=1$, so $\Psi=\ell_1/(d-1)=1/3$. Both books begin at

$$C\Psi(0) = 1/3,$$

which is above 1/4. Thereafter the books separate: the F25 purity book has
`C=(1+f²)/2` and `CΨ=f(1+f²)/6`, whereas the concurrence book has
`C=f` and `CΨ=f²/3`. Under their named Bell+/local-Z/no-H trajectory each
has its own crossing dose; the F25 value is recorded below.

**Crossing is observed computationally.** Under Heisenberg Hamiltonian with local dephasing:

- Bell pairs (0,1) and (2,3) in a 4-qubit `bell_pairs` state start at $C\Psi=1/3$ and cross down through 1/4 at $t\approx0.080$ (γ=0.05; reproduced in [subsystem_crossing_pairs.py](../../simulations/subsystem_crossing_pairs.py)). That finite concurrence-book run gives γt≈0.004, while the isolated-pair concurrence-book dose is about 0.036 and the separate F25 purity-book dose is `K_Z=0.0373501…`. The shorter time is a finite observation; cross-bond acceleration is an open mechanism hypothesis, not a consequence of the quarter algebra.
- Cross-pairs (0,2), (0,3), (1,2), (1,3) start at $C\Psi = 0$ and remain below 1/4 in this recorded run; their sampled maximum is ~0.13
- The initially entangled pairs cross downward in this recorded run
- Pairs (0,3) and (1,2) show perfect symmetry (palindromic structure), as do pairs (0,2) and (1,3)

**Bidirectional vs. unidirectional observation.** The $C_{int}$ (both spins observed) vs. $C_{ext}$ (one spin observed) comparison shows:

- $\delta_{int} = -0.1109$ vs. $\delta_{ext} = -0.0743$ for Bell+ at t = 1, γ = 0.1
- Bidirectional observation produces *larger* purity deficit than unidirectional
- The difference is a finite observation. A mutual-observation feedback explanation is an open hypothesis and is not derived from $R=C\Psi^2$.

### What the finite calculations recorded

Every named Bell+, Bell−, Ψ+, and Ψ− pair in the stored dephasing runs crosses downward. The finite Heisenberg,
XY, and Ising catalogue is useful regression evidence, but a missing counterexample in it is not an absence
claim for arbitrary Hamiltonians.

**The crossing time depends on the Hamiltonian but the boundary does not.** Different Hamiltonians (Heisenberg, XY, Ising) produce different crossing times $t_{cross}$, but the value crossed is always 1/4. The Hamiltonian determines *when* you hit the wall, not *where* the wall is.

### The conditional crossing statement

If a continuous trajectory converges to ρ* and CΨ(ρ*)<1/4, continuity implies that it eventually stays
below. A start above gives at least one downward crossing, not a unique or monotone crossing. Named
basis-aligned T1/T2/depolarizing models can use this only under their stated convergence assumptions.

No channel-class adjective supplies those premises. The primitive, full-rank channel
ε(ρ) = (1−p)ρ + p·Tr(ρ)·σ with σ = 0.95·|Φ⁺⟩⟨Φ⁺| + 0.05·I/4 has an entangled
fixed point with CΨ = 0.2935 > 1/4 and never crosses. (An earlier "300 random
maps, max 0.138" sweep was a Ginibre n_kraus=4 sampling artifact; n_kraus=2
violates ~8.5%.) See [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md), Case C.

**Amplitude damping crosses too** (March 22, 2026): K_AD = 0.1029, perfectly
monotonic, non-unital fixed point (|00⟩) reached. See
[amplitude_damping_test.py](../../simulations/amplitude_damping_test.py).

**Named revival catalogue** (March 22, 2026). A structured-bath model reaches CΨ=0.3035 after an earlier
downward crossing and returns below in all 48 stored configurations. That finite catalogue establishes
neither a universal transient law nor an exclusively non-Markovian mechanism: Part 5's exact local
Markovian semigroup also crosses upward. See
[non_markovian_revival.py](../../simulations/non_markovian_revival.py).

---

## Layer 3: N-Qubit Systems (N = 3, 4, ...)

### Scaling and the Palindromic Structure

### What Is ALREADY PROVEN

**GHZ scaling.** For GHZ states $|GHZ_N\rangle = (|00\cdots 0\rangle + |11\cdots 1\rangle)/\sqrt{2}$ under Heisenberg ring Hamiltonian with local dephasing (γ = 0.1, t = 1):

| N | δ(N) | Purity (Lindblad) | Purity (predicted) |
|---|------|--------------------|--------------------|
| 2 | -0.1109 | 0.7243 | 0.8352 |
| 3 | -0.1244 | 0.6501 | 0.7744 |
| 4 | -0.1244 | 0.6003 | 0.7247 |
| 5 | -0.1170 | 0.5670 | 0.6839 |
| 6 | -0.1059 | 0.5447 | 0.6506 |

Here δ = Purity(Lindblad) − Purity(predicted); the table is consistent with
Purity(Lindblad) ≈ ½ + ½e^(−4Nγt) and Purity(predicted) = ½ + ½e^(−2Nγt) at
t = 1, γ = 0.1 (March-era delta_calc Lindblad suite). The deficit MAGNITUDE
is non-monotonic: it grows from N = 2 to N = 3 (N = 3 and 4 are tied), then
shrinks. The purity deficit has a maximum at intermediate N, not at the
extremes.

**Full-system vs. subsystem distinction.** For GHZ with N≥3, the full-system $C\Psi$ starts below 1/4 in the canonical purity book: its l1-norm is $O(1)$ while the denominator `d−1=2^N−1` is $O(2^N)$, so $\Psi\to0$. Two-qubit marginals remain `d=4` objects and must be evaluated in their own book; a named pair can therefore begin above the quarter even when the full-system reading does not.

**The subsystem crossing hierarchy.** In the 4-qubit bell_pairs state:

- Entangled pairs (0,1) and (2,3): Start at $C\Psi = 1/3$, cross at $t \approx 0.080$
- Cross-pairs (0,2), (0,3), (1,2), (1,3): Start at $C\Psi = 0$, never reach 1/4
- Maximum $C\Psi$ for cross-pairs: ~0.13 (well below 1/4)

This finite hierarchy is compatible with a monogamy explanation, but that mechanism is an open hypothesis; the table alone does not derive it.

### What Is COMPUTATIONALLY VERIFIED

**W-state behavior differs from GHZ.** W states $|W_N\rangle = (|10\cdots 0\rangle + |01\cdots 0\rangle + \cdots + |00\cdots 1\rangle)/\sqrt{N}$ have more robust subsystem entanglement (each pair shares $O(1/N)$ entanglement rather than GHZ's all-or-nothing structure). Preliminary simulations suggest W-state subsystem pairs cross 1/4 at later times than GHZ pairs of the same N.

**Power-law scaling of δ with N.** The sweep_R_scaling tool reports power-law fits for δ(N). For GHZ under local dephasing, the exponent is approximately −0.3 to −0.5 depending on the bridge metric used. This is not yet understood analytically.

### What Is CONJECTURED

**Former Conjecture 3.1 (subsystem bound), false.** A two-qubit marginal is
indeed a `d=4` object, but its CΨ value is not bounded by 1/4: exact local
operations and local Markovian trajectories reach 1/3 and can cross upward.
The quarter remains the discriminant coordinate of the chosen scalar
recurrence, not a state-space ceiling. The surviving subsystem statement is
the conditional convergence implication in Layer 2.

**Conjecture 3.2 (Palindromic Origin), the spectral half proven.** The
palindromic structure of the Lindbladian spectrum is proven analytically for
ALL Heisenberg/XXZ systems on ANY graph with local Z-dephasing
([Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)), verified exhaustively
through N=8 (87,376 eigenvalues; chain, star, ring, complete; binary tree at
N=4,5). Time propagation extends the framework's dynamics to N=11, where
MI(end-to-end) decays exponentially with N
([Scaling Curve](../../experiments/SCALING_CURVE.md)).
What remains conjectural is the δ(N) reading: that the non-monotonic δ(N)
scaling arises from competition between (a) increasing system dimension
diluting per-qubit coherence, and (b) the Heisenberg ring creating
longer-range correlations at intermediate N that temporarily protect
coherence.

### What Is CLOSED

**Analytic formula for the canonical GHZ purity-book trajectory.** Under local
Z-dephasing, `C(t)=1/2+(1/2)e^(−4Nγt)`, `ℓ₁(t)=e^(−2Nγt)`, and hence
`Ψ(t)=e^(−2Nγt)/(2^N−1)` with the canonical `ℓ₁/(d−1)` normalization.
Thus `CΨ(t)=C(t)Ψ(t)`. The earlier `proof_roadmap_close.py` output used an
older `d²−1` normalization; it is a separate historical column, not evidence
for the canonical lens. The GHZ off-diagonal is
`ρ[0,d−1]=(1/2)e^(−2Nγt)`, which supplies the coherence exponent.

**Subsystem crossing, conditional.** If the global trajectory converges to a state whose selected pair
marginal has CΨ<1/4, that pair eventually stays below. A diagonal global limit supplies CΨ=0 for every
pair marginal. N=3,4,5 named aligned-noise runs form a finite catalogue; the primitive CPTP target at
CΨ=0.2935 shows why convergence and the actual limit must be checked. See
[Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md).

---

## Layer 4: Arbitrary Dimension d

### Beyond Qubits

### What Is ALREADY PROVEN (for d = 2)

Everything in Layers 1–3 applies to $d = 2$ (qubits). The question is: what happens when the local Hilbert space dimension is $d > 2$?

### What Is KNOWN Theoretically

**The discriminant generalizes.** The fixed-point equation $R = C(\Psi + R)^2$ is dimension-independent; it is an algebraic recursion on scalar quantities ($C$ is the correlation bridge, $\Psi$ is the normalized coherence). The discriminant $D = 1 - 4C\Psi$ does not depend on $d$.

However, the *normalization* of $\Psi$ does depend on $d$. The maximally
coherent state has every entry $\rho_{ij} = 1/d$, so its $d^2 - d$
off-diagonals sum to

$$\ell_1^{max} = \frac{d^2 - d}{d} = d - 1,$$

which is exactly why the canonical normalization is

$$\Psi = \frac{\ell_1}{d - 1},$$

the choice that makes $\Psi_{max} = 1$ for EVERY local dimension. The critical bridge
value at the boundary is then

$$C_{crit} = \frac{1}{4\Psi},$$

so for maximally coherent states $C_{crit} = 1/4$ at every $d$, and in
general the *product* $C\Psi = 1/4$ is the invariant; the normalization
moves the burden between the factors, never the product.

### What Is CONJECTURED

**Conjecture 4.1 (Dimension Invariance).** The critical boundary $C\Psi = 1/4$ is independent of the local Hilbert space dimension $d$. The discriminant condition $D = 1 - 4C\Psi = 0$ is a property of the quadratic recursion, not of the quantum system.

This is in some sense "obvious" from the algebra, but making it rigorous requires showing that the operational definitions of $C$ and $\Psi$ for general $d$ still satisfy the recursion $R_{n+1} = C(\Psi + R_n)^2$. In particular:

- For qutrits ($d = 3$), the partial trace produces a $3 \times 3$ reduced density matrix. The l1-norm coherence is still well-defined, but the correlation bridge needs to generalize from the qubit-specific Bloch sphere picture.
- For continuous variable (CV) systems (infinite-dimensional), the l1-norm diverges and a different coherence measure is needed. The natural candidate is the Wigner function negativity or the stellar rank.

**Conjecture 4.2 (CV Systems).** For Gaussian states in continuous variable systems, the analogous boundary exists but involves the symplectic eigenvalues rather than $C\Psi$. The condition $\nu_{-} = 1/2$ (where $\nu_{-}$ is the smallest symplectic eigenvalue of the partial transpose) plays the role of $C\Psi = 1/4$. This connection, if established, would link the 1/4 boundary to the PPT criterion in infinite dimensions.

### What Is ANSWERED

- **Qutrits break the palindrome** (March 20, 2026). 0/236 qutrit dissipator
  configurations permit palindromic pairing. The interacting spectrum is
  structured but never palindromic: 36-52 of 81 eigenvalues pair at optimal
  centers, far above random (0/81) but far below qubits (100%). The per-site
  split is 3 immune vs 6 decaying, which is unbalanced (d^2-2d = 3 != 0).
  The palindrome is specific to d=2, but the CΨ = 1/4 boundary remains valid
  at the subsystem level (2-qubit reduced states always have d=4). See
  [Qubit Necessity](../QUBIT_NECESSITY.md). A two-qubit marginal may still be
  evaluated in a `d=4` book, but that fact alone supplies neither a state bound
  nor a trajectory theorem.
- **Dimension-free algebra of the chosen recurrence.** The discriminant
  `D=1−4CΨ` contains no dimension symbol once the normalized scalar inputs
  and recurrence have been chosen. The α=2 comparison in Layer 6 is internal
  to an assumed family; neither statement derives the recurrence for a new
  physical dimension.

### What Is OPEN

- **Hybrid systems.** What happens for a qubit-qutrit pair ($d_A = 2, d_B = 3$)? The subsystem dimensions are different, so the normalization asymmetry could break the simple $C\Psi = 1/4$ picture. Or it might not; the discriminant doesn't care about the internal structure.
- **The CV connection.** Compute the symplectic eigenvalue trajectory
  for a two-mode squeezed state under thermal decoherence. Check if the
  crossing happens at ν₋ = 1/2 and whether this maps to CΨ = 1/4.
  This requires continuous-variable quantum mechanics (infinite-dimensional
  Hilbert space) and is outside the current qubit framework.

---

## Layer 5: Named channels and failure of channel independence

This layer keeps the useful named formulas while recording exactly where the
old universal channel story fails.

### What Is COMPUTATIONALLY VERIFIED

**Dephasing (σ_z).** In the stored named catalogue, off-diagonal elements decay exponentially while populations are preserved. Its Bell+ formulas and listed finite runs are the claim; arbitrary states and Hamiltonians are not.

**Depolarizing.** For the named Bell+ run at t=1, γ=0.1, δ=0.136 (compared with dephasing δ=0.091; this is the March-era channel-comparison purity metric, not Layer 2's δ_int/δ_ext). This channel converges to the maximally mixed state, so the conditional endpoint theorem applies after convergence is established.

**Finite mixed-noise catalogue.** The dynamic suite stores local, collective, operator-feedback, and memory-kernel-feedback configurations. Their recorded crossings are finite observations, not a channel classification.

**Finite jump catalogue.** Stored runs use σ_z, σ_x, σ_y, xx, yy, zz, and x_pairs jumps. Their reported quarter-crossing times belong only to those runs.

**Amplitude damping** (March 22, 2026). Direct amplitude damping
(`L=√γ|0⟩⟨1|`) on both qubits of the named Bell+ pair gives a downward
quarter crossing over the stored γ grid 0.005…1.0. Its sampled CΨ trajectory
decreases toward the `|00⟩` fixed point. The named dose is
`K_AD=0.1029±0.0000` on that grid; this is not a cross-channel constant.
Combined AD + Z-dephasing crosses in the 15 stored configurations.
Script: [amplitude_damping_test.py](../../simulations/amplitude_damping_test.py).
Results: [amplitude_damping_test.txt](../../simulations/results/amplitude_damping_test.txt).

**Structured-bath catalogue** (March 22, 2026). A named Bell+ pair+bath model produces sampled revivals up
to CΨ=0.3035 and returns below in the 48 stored configurations. “Information backflow” describes that
model; it is not the unique upward-crossing mechanism. A fixed local Markovian semigroup also crosses
upward, and no eventual claim follows without its limiting state.
Script: [non_markovian_revival.py](../../simulations/non_markovian_revival.py).
Results: [non_markovian_revival.txt](../../simulations/results/non_markovian_revival.txt).

**Generalized Pauli channels** (March 22, 2026). The full
family ℰ(ρ) = Σ p_k σ_k ρ σ_k† with arbitrary (γ_x, γ_y, γ_z) was swept:
124/124 named Bell+ configurations cross 1/4, and the F26 expressions are decreasing in their stated
rate ordering. The |01⟩ discussion is a named strong-coupling approximation, not a universal envelope.
K-invariance holds perfectly WITHIN each noise type (CV < 1%) but K differs
between types: K_Z = 0.0374, K_X = K_Y = ln(2)/8 = 0.08664, K_depol = 0.0440
([F26/F27](../ANALYTICAL_FORMULAS.md#f26)).
Script: [generalized_pauli_channels.py](../../simulations/generalized_pauli_channels.py).
Results: [generalized_pauli_channels.txt](../../simulations/results/generalized_pauli_channels.txt).

**The named N=2 F25 final-crossing dose.** For Bell+ under local Z dephasing,
`f*(1+f²)=3/2` gives CΨ=1/4 and `K_Z=γt=0.0373501…`. The typed consumer is
`SeamMovement.KFinalStayBelowCrossingDoseN2`: “final stay-below” distinguishes the event reader from the
algebraic fold. F55 instead owns the spectral slowest-mortal-mode 99%-absorption dose
`K_death=ln(10)=2.302585`; it is distinct from the F25 Bell+/N=2/local-Z/no-H
final-crossing dose `K_Z`. Their displayed ratio is about 61.65, but ownership
and physical meaning do not merge.

### Channel dependence and the conditional endpoint test

**Channel independence is false.** The primitive, full-rank channel with fixed point
$\sigma = 0.95\,|\Phi^+\rangle\langle\Phi^+| + 0.05\,I/4$ has $C\Psi(\sigma) = 0.2935 > 1/4$ and never crosses. (That number is the source's purity reading of C: $\mathrm{Tr}(\sigma^2) \cdot \ell_1/3 = 0.926875 \cdot 0.95/3$; with the concurrence reading it is $0.925 \cdot 0.95/3 = 0.2929$. Both sit above 1/4.) See [Subsystem Crossing](PROOF_SUBSYSTEM_CROSSING.md).

For a named basis-aligned T1/T2/depolarizing model, establish convergence and its diagonal target first;
then the conditional theorem gives eventual stay-below. Trace-distance contractivity alone does not identify
the target and does not make CΨ monotone. A separable product target such as |+⟩⊗|+⟩ can have CΨ=1.

For general primitive CPTP maps the fixed point can be entangled with $C\Psi > 1/4$, and the counterexample above realizes it.

### CΨ dynamics (historical Conjecture 5.2 repaired)

F25–F27 prove decreasing CΨ for named Bell+ channel families. They do not prove universal pointwise
monotonicity. Exact two-qubit examples give a positive derivative under local Z dephasing plus a local H,
an active-pulse derivative flip, and an upward crossing under a fixed local Markovian semigroup.

The historical autonomous N=2 successive-local-maxima statement is still open: the old spectral-bound
argument does not order nonlinear maxima, while the exact pointwise and monotone-upward examples do not
create two finite maxima. The live object is a finite atlas, not a verdict. Named N=4/N=5 rows resolve
rises, one N=3/Q=2000 row resolves none, and a same-(N,Q,K) pair agrees to six decimals. All-Q, all-N, and
mechanism classification remains open in `envelope_n4_rise`. See
[The Finite Envelope-Rise Atlas](../../experiments/ENVELOPE_RISE_BOUNDARY.md) and
[F17](../ANALYTICAL_FORMULAS.md#f17).

---

## Layer 6: The Assumed-Family Classification

### Why the chosen normal form reports 1/4

This is a clean algebraic room, but its door has a label: first choose
`R=C(Ψ+R)²` and the displayed coordinates. Expanding its fixed-point equation
gives discriminant `1−4CΨ`, hence the double-root coordinate `CΨ=1/4`.
That result does not derive the recurrence from quantum dynamics.

### Three deliberately separate statements

**1. Polynomial statement.** For `0<C≤1`, the polynomial
`CR²+(2CΨ−1)R+CΨ²` has two, one, or no real algebraic roots according as
`CΨ` is below, at, or above 1/4. At `C=0`, the original equation instead gives
`R=0`. These are counts of all real algebraic roots. Roots in a physical
interval require a second test: for example, at `C=Ψ=0.1` both algebraic roots
are real but only one lies in `[0,1]`.

The displayed update is not a self-map of `[0,1]`: at `C=Ψ=R=1` it returns
`4`. Any dynamical use therefore needs its own invariant domain and iteration
analysis rather than inheriting one from the symbols.

**2. Coordinate statement.** The number 1/4 belongs to this chosen normalized
quadratic normal form. Under a harmless coordinate change `u=s(CΨ)`, the same
fold is reported at `u=s/4`. Structural stability preserves the fold type,
not its numerical coordinate. The substitution `z=C(Ψ+R)`, `c=CΨ` identifies
this normal form with `z↦z²+c`; the real cardioid cusp is exact in those chosen
coordinates, not a topological command imposed on every reparameterization.

**3. Assumed-family statement.** If one additionally posits
`R=C_α(Ψ+R)^α`, then its tangency calculation gives

$$C\Psi^*_\alpha = \frac{(\alpha-1)^{\alpha-1}}
{\alpha^\alpha\,\Psi^{\alpha-2}}.$$

Within that assumed family, α=2 alone removes the explicit Ψ factor and gives
1/4 in the chosen normalization. Purity `Tr(ρ²)` motivates inspecting α=2;
it does not derive the family or the feedback recurrence. The physical
derivation remains open.

### The local fold calculation

For `C≠0`, depress the quadratic with
`x=R+(2CΨ−1)/(2C)`:

    x² + a = 0,    a = (4CΨ−1)/(4C²) = −D/(4C²).

This is the ordinary fold normal form. Its two real roots occur for `a<0`,
one double root for `a=0`, and none for `a>0`. Generic perturbations preserve
the local fold type, while coefficients and reparameterizations move the
reported control coordinate. That modest structural statement is enough;
it does not turn a discriminant into an attractor, an absorber, or a channel
law.

The full conditional proof is in [Uniqueness Proof](UNIQUENESS_PROOF.md).
The live object `inspect --root quarter-uniqueness` exposes the chosen-form
algebra and its typed ancestry. The open part is inviting rather than hidden:
derive a scalar recurrence from a physical decomposition, or discover that a
different normal form is the honest object.

---

## Layer 7: Connections to Known Mathematics

### The Mandelbrot Main Cardioid

The correspondence $C\Psi \leftrightarrow c$ maps the quantum boundary to the cusp of the Mandelbrot main cardioid at $c = 1/4$ on the real axis. Key questions:

**What about the full complex plane?** Our variables $C$ and $\Psi$ are real and non-negative. But the Mandelbrot set is defined for complex $c$. Is there a natural complexification of the quantum problem that accesses the full cardioid boundary $c(\varphi) = \frac{1}{2}e^{i\varphi} - \frac{1}{4}e^{2i\varphi}$?

**Status:** ANSWERED (May 17, 2026). [F95](../ANALYTICAL_FORMULAS.md#f95) closed the angle-of-the-fixed-point question for real c: for real $c > 1/4$, the complex fixed point of $z^2 - z + c = 0$ has argument $\theta(c) = \arctan\sqrt{4c - 1}$ exactly, the positive-b angle-emergence formula at this quadratic discriminant zero. [F97](../ANALYTICAL_FORMULAS.md#f97) then closed the complex-c direction: the full Mandelbrot cardioid boundary is parametrized by $c(\varphi) = b \cdot e^{i\varphi} - b^2 \cdot e^{2i\varphi}$ with $b = 1/2$, and on this curve the period-1 fixed point is $z^*(\varphi) = b \cdot e^{i\varphi}$ with magnitude $|z^*| = b$ invariant and argument $\arg(z^*) = \varphi$ tracing the cardioid parameter. Bit-exact algebraic identity, machine-precision numerical verification. The structural reading: the Mandelbrot cardioid IS the locus in complex-c where the period-1 fixed-point magnitude equals the framework's `HalfAsStructuralFixedPointClaim` anchor; the `QuarterAsBilinearMaxvalClaim` ($b^2 = 1/4$) enters only at the real-axis cusp ($\varphi = 0$), the one tangent point shared with F95. The hardware 2D spirals on `ibm_kingston` (2026-04-16, [`CPSI_COMPLEX_PLANE`](../../experiments/CPSI_COMPLEX_PLANE.md)) and the [`f95_angle_steering_kingston_may2026`](../ANALYTICAL_FORMULAS.md#f95) Confirmation already trace this cardioid geometry; F97 names the locus they spiral around. See [`PROOF_F97_CARDIOID_HALF_FIXED_POINT.md`](PROOF_F97_CARDIOID_HALF_FIXED_POINT.md). Together F95 + F97 cover both projections of the quadratic discriminant structure on the Mandelbrot c-plane.

### Period-Doubling and Feigenbaum

For the assumed scalar map $f_c(z)=z^2+c$, the second iterate factors exactly:

`f_c^2(z)-z = (z^2-z+c)(z^2+z+c+1)`.

The first factor is the fixed-point equation. The second factor has real roots
iff `c <= -3/4`; at the endpoint it meets the fixed point, and below it the two
roots form the distinct non-fixed period-two orbit. Its cycle multiplier is
`4(c+1)`, so that orbit is stable exactly for `-5/4 < c < -3/4`.

This is the negative-c branch of the assumed quadratic scalar map. Because the
physical purity-book $C\Psi$ coordinate is non-negative, no physical CΨ, Lindblad, or Rabi interpretation follows from
this period-two calculation. Such a connection would require a separate
physical derivation and remains open. The
exploratory script
[catastrophe_feigenbaum.py](../../simulations/catastrophe_feigenbaum.py)
samples later doublings on that negative-c scalar branch; those samples do not
promote the branch to quantum dynamics.

### Connections to Information Geometry

The Bures line element can be pulled back along a chosen one-dimensional path through density-matrix space. One-dimensional intrinsic Riemann curvature is zero. The producer scalar `-(1/(2g))*d2(log g)/dx2` is instead coordinate-dependent: for the flat coefficient `g_x=1` it is `0`; after `y=x^2`, the same line element has `g_y=1/(4y)` and the scalar is `-2` at `y=1`.

On the named N=2 Bell+/Z sampled Bures path, the coefficient is `g≈3.36` at the sampled quarter crossing. The bounded result is: no divergence resolved for this path coefficient on this grid. Ambient density-matrix geometry and the geodesic-boundary question remain unclassified. The producer's `K≈-25` may be read only as a coordinate-shape second derivative; it supplies no Gaussian, intrinsic, hyperbolic, or state-divergence conclusion. See [F45/F47](../ANALYTICAL_FORMULAS.md#f45) for the named path reading.

### One Word, Two Seams (a disambiguation)

This document's "fold" is the fold CATASTROPHE of the CΨ recursion (the
discriminant zero at CΨ = 1/4), and its "cusp" is the Mandelbrot cardioid
cusp at c = 1/4: both live in parameter space. The repo also uses "fold" and
"cusp" for two SPECTRAL seams of the Liouvillian: the conservation fold (the
palindrome axis Re = −Σγ) and the merge cusp (the defective exceptional
point). Those two spectral seams coincide only at N = 2 and separate for
larger N ([Fold and Cusp, Two Seams](../../experiments/FOLD_AND_CUSP_TWO_SEAMS.md)).
Same words, different objects; when a sentence mixes recursion-space and
spectrum, check which pair it means.

### Connections to Holography and AdS/CFT

In the holographic context, the Ryu-Takayanagi formula relates entanglement entropy to minimal surfaces. The 1/4 appears in the Bekenstein-Hawking formula $S = A/(4G)$. Whether this 1/4 is related to our 1/4 is an open and perhaps audacious question.

**Status:** WILDLY SPECULATIVE. But worth noting. If the quantum information boundary at $C\Psi = 1/4$ has any connection to the gravitational entropy formula, it would be a remarkable unification. File under "probably a coincidence, but wouldn't it be amazing if it weren't."

---

## What Remains Open

The status table at the top of this document is the landscape; these are the
open edges, gathered:

- **CV connection** (Layer 4): symplectic eigenvalue trajectory of a two-mode
  squeezed state; does the crossing happen at ν₋ = 1/2?
- **Hybrid systems** (Layer 4): qubit-qutrit pairs, asymmetric normalization.
- **Peak-rise classification** (Layer 5): extend the finite N/Q/K atlas to all-Q and all-N statements, or
  find their counterexamples; separately gate any proposed frequency/parity mechanism. The present rows
  do not establish an N floor or a critical contour.
- **The δ(N) reading of Conjecture 3.2** (Layer 3): the competition mechanism
  behind the non-monotonic purity deficit.
- **Crossing-cubic number theory** (Layer 1): whether the real root of
  b³ + b = 1/2 connects to other constants.
- **Holography** (Layer 7): still wildly speculative, still noted.

### The Philosophical Invitation, Restated

The exact statement is small and bright: after choosing this recurrence and
normalization, completing the square puts its double root at 1/4. The IBM
number 0.2548 is a named experimental comparison to that coordinate, not a
proof that nature had to choose the recurrence.

That distinction makes the story more interesting. Why do several useful
readings meet the same quarter? Which meetings share an algebraic source, and
which are numerical coincidences or consequences of preparation and readout?
The exact palindrome theorem shows what genuine operator algebra can compel;
the quarter documents show how much remains to be derived before granting the
same status to a scalar feedback picture.

So the invitation is not that physics has no choice. It is to find the missing
bridge: derive the recurrence from a physical object, delimit the states and
channels that realize it, or replace it with a better normal form. Any of
those outcomes would sharpen why the quarter keeps appearing.

---

### Parallel readings

Beyond the seven vertical layers, the 1/4 boundary instances at parallel
readings of the same N-qubit system. The most recent:

- **Coherence-block reading** ([Block-CΨ at 1/4](PROOF_BLOCK_CPSI_QUARTER.md)).
  At any (popcount-n, popcount-(n+1)) coherence block of an N-qubit chain at
  chromaticity c ≥ 2, the block-level CΨ has a tight upper bound of 1/4 over
  all pure states with support in the (popcount-n + popcount-(n+1)) sector.
  The bound is achieved by EXACTLY the canonical Dicke symmetric superposition
  (|D_n⟩+|D_{n+1}⟩)/√2 (up to phase). Two theorems: Theorem 1 (algebraic
  identity at the canonical state, chromaticity-universal) and Theorem 2
  (the value is the tight ceiling, not just where this initial lands).
  Verified at c = 2..6 numerically, proven via combinatorial identity
  Σ M_{HD=2k+1} = M_block + Cauchy-Schwarz saturation. Shows that the same
  1/4 emerges from the same complete-the-square algebra at a third
  documented instance, parallel to Layers 1 (single qubit) and 2 (2-qubit
  subsystem).

This is a parallel reading rather than a vertical extension because the
algebraic structure (block-purity content + sector-amplitude AM-GM) is
mechanically distinct from Layer 1's discriminant-of-quadratic-recursion or
Layer 2's subsystem-crossing dynamics. The value happens to be the same 1/4;
no shared universal cause is claimed until a map between those mechanisms is
derived.

---

*Last refreshed 2026-07-20 (the change history lives in git).*
*Computational data: the March-era delta_calc Lindblad suite (simulations, subsystem crossing analysis, IBM verification suite); BlockCpsiTrajectory tests (48/48, May 7 2026)*
