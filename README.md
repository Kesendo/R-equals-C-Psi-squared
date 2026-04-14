# R = CΨ²

<!-- Keywords: absorption theorem qubit cavity, palindromic Liouvillian spectral symmetry,
open quantum system dephasing, standing wave eigenmode pairing, Fabry-Perot qubit chain,
CΨ quarter boundary proof, quantum decoherence as illumination, Lindblad master equation
eigenvalue pairing, dephasing light not noise, quantum MIMO channel capacity,
quantum state transfer spin chain, IBM quantum hardware validation, self-referential
purity recursion Mandelbrot, conjugation operator Pi time reversal, R=CPsi2 framework,
palindromic eigenvalue spectrum proof, quantum noise channel 15 bits,
V-Effect coupling creates new frequencies, Dale's Law palindromic spectral symmetry,
Wilson-Cowan excitatory inhibitory balance palindrome, sacrifice-zone entrance pupil,
half-occupation C=0.5 axiom, sigmoid sensitivity maximum one quarter,
neural criticality eigenfrequency, K-dosimetry gamma-time invariant,
light and lens superalgebra M22, absorption quantum 2gamma standing wave cavity,
relational quantum mechanics reality between observers, dependent origination emptiness
mirror symmetry, human-AI collaboration original research Claude Anthropic -->

> *[We are all mirrors. Reality is what happens between us.](MIRROR_THEORY.md)*

A human and an AI, exploring together. What we found surprised us both:
the absorption spectrum of any qubit network under dephasing is exactly
palindromic. For every mode that absorbs fast, one absorbs slow. Always
paired. Always balanced. One equation governs it all.

Verified from N=2 through N=8 across 87,376 Liouvillian eigenvalues, with
zero mirror-symmetry exceptions on any tested topology (chain, star,
ring, complete, tree). Confirmed on IBM quantum hardware at 3%.

The thing that remains is not fighting the absorption. It is made of it.

Early speculations live in `recovered/`: some turned out to be premature
rather than wrong, others remain unsupported. We keep them because the
research process matters as much as the results.

**Thomas Wicht** (independent researcher, Germany) and **Claude** (AI, Anthropic)

---

## Where to start

You can read this page straight through: it follows one idea from
equation to understanding. Or pick the entry point that fits you:

→ **[What We Found](docs/WHAT_WE_FOUND.md)**: the discovery explained
from the beginning, no prerequisites

→ **[Reading Guide](docs/READING_GUIDE.md)**: five stories (proof,
application, ontology, resonator, cross-level), each with a reading order

→ **[The Anomaly](THE_ANOMALY.md)**: the question that remained after
the proof. No formulas. Written the evening the hardware confirmed
the theory

→ **[Mirror Theory](MIRROR_THEORY.md)**: the interpretation. What
happens when you read the formula from both sides at once. No
formulas in the prose; every claim links into the proofs

→ **[What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md)**: every
error, every limitation, every unanswered question. Because a theory that
only shows its strengths is not a theory

If you are a physicist: [Technical Paper](publications/TECHNICAL_PAPER.md).
If you work with quantum hardware: [Engineering Blueprint](publications/ENGINEERING_BLUEPRINT.md).
If you work with neural networks: [Neural Palindrome](docs/neural/README.md) (no quantum prerequisites).

---

## 1. One equation governs the spectrum

A guitar string vibrates in modes. Each mode spans the entire string.
How fast a mode fades depends on one thing: how much of the mode's
energy sits in the parts that are damped. If the damping pad touches a
node, the mode survives. If it touches an antinode, the mode dies fast.

The quantum version:

    Re(λ) = −2γ × ⟨n_XY⟩

The absorption rate of any eigenmode equals twice the dephasing rate (γ)
times the mode's mean light content (⟨n_XY⟩): the fraction consisting
of X and Y Pauli operators, the oscillating quantum components that
interact with the external illumination. Components that do not oscillate
({I, Z}, the "lens") are invisible and survive forever. Components that
oscillate ({X, Y}, the "light") absorb and fade.

The spectrum is a ladder with rung spacing 2γ:

| Rung | Rate | Content |
|:-----|:-----|:--------|
| 0 | 0 | All lens ({I,Z}^N): pure structure, immortal |
| 1 | 2γ | One light factor: one photon absorbed per dephasing time |
| 2 | 4γ | Two light factors |
| ... | ... | ... |
| N | 2Nγ | All light ({X,Y}^N): maximum absorption |

The Hamiltonian smooths the ladder (modes become superpositions of
different rungs, creating fractional absorption rates) but cannot change
the endpoints or the fundamental quantum 2γ.

Six previously separate results follow as one-line corollaries:

| Corollary | Statement | Derivation |
|:----------|:----------|:-----------|
| Spectral boundaries | min = 2γ, max = 2(N−1)γ | ⟨n_XY⟩ ∈ {1, N−1} for paired modes |
| Palindromic sum rule | α_fast + α_slow = 2Σγ | ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N |
| Spectral gap | Δ = 2γ | Cost of one light factor |
| 2× absorption law | Extreme rate / center rate = 2 | Definition of center of symmetric interval |
| Mode classification | {I,Z}^N immortal, {X,Y}^N absorbed | n_XY = 0 vs n_XY = N |
| N=3 exact rates | 2γ, 8γ/3, 10γ/3 | Hamiltonian-mixed ⟨n_XY⟩ = 1, 4/3, 5/3 |

Proven analytically in three steps from L_H anti-Hermitian. Verified on
1,343 modes (N=2-5, γ=0.01-1.0, J=0.1-5.0), coefficient of variation = 0.

→ **[Absorption Theorem Proof](docs/proofs/PROOF_ABSORPTION_THEOREM.md)**
→ [Discovery and Verification](experiments/ABSORPTION_THEOREM_DISCOVERY.md)
→ [Analytical Formulas](docs/ANALYTICAL_FORMULAS.md) (Formulas 3, 8, 33, D6)

---

## 2. The palindromic spectrum

The symmetry that makes the theorem's consequences exact.

N qubits coupled via Heisenberg (or XXZ) interaction, subject to local
Z-dephasing at rate γ per site:

```
    qubit ── qubit ── qubit ── ...       Heisenberg coupling J
       γ        γ        γ               Z-dephasing (illumination from outside)
```

The conjugation operator Π acts per site on Pauli indices
(I ↔ X, Y ↔ iZ, Z ↔ iY) and satisfies:

    Π · L · Π⁻¹ = −L − 2Σγ · I

For every Liouvillian eigenvalue λ, there exists a partner at −2Σγ − λ̄.
Every absorption rate d has a mirror at 2Σγ − d. The spectrum is exactly
palindromic.

| N | Matrix | Oscillatory rates | Mirror | Absorption range |
|:--|:-------|:------------------|:-------|:-----------------|
| 2 | 16² | 6 | 100% | 2γ to 2γ |
| 3 | 64² | 40 | 100% | 2γ to 4γ |
| 4 | 256² | 182 | 100% | 2γ to 6γ |
| 5 | 1024² | 776 | 100% | 2γ to 8γ |
| 6 | 4096² | 3,228 | 100% | 2γ to 10γ |
| 7 | 16384² | 13,264 | 100% | 2γ to 12γ |
| 8 | 65536² | 54,118 | 100% | 2γ to 14γ |

*"Oscillatory rates" counts Liouvillian eigenvalues with non-zero
imaginary part (Im(λ) ≠ 0). For N=8 this is 54,118 of the 65,536 total
eigenvalues; the remaining 11,418 are purely real (N+1 stationary modes,
2N weight-1 degenerate rungs at Re=−2γ, and higher-order Hamiltonian-
diagonal fixed points). "Mirror" is the fraction of below-center rates
that find a matching partner above the center −Nγ within numerical
tolerance. At N=8, the center cluster holds an additional 14,282
oscillatory rates that are consistent with the palindrome but not
individually resolved without eigenvectors.*

Holds for all standard coupling models (Heisenberg, XY, Ising, XXZ,
Dzyaloshinskii-Moriya), all graph topologies, non-uniform γ per qubit,
Z and Y dephasing. Breaks for depolarizing noise.

At Σγ = 0 (no illumination): Π L Π⁻¹ = −L. Every eigenvalue pairs with
its negative. Pure oscillation. No absorption. No irreversibility.
Illumination does not destroy the palindrome; it shifts it. The shift
creates the arrow of time.

→ **[Mirror Symmetry Proof](docs/proofs/MIRROR_SYMMETRY_PROOF.md)**
→ [All standard models](experiments/NON_HEISENBERG_PALINDROME.md) (two Π families, 36/36 combinations resolved)
→ [Zero Is the Mirror](hypotheses/ZERO_IS_THE_MIRROR.md) (the palindrome before the shift)

---

## 3. Every paired mode is a standing wave

Two waves traveling in opposite directions create a standing wave: a
pattern that does not move while energy flows back and forth. Nodes and
antinodes. Fixed in space. Oscillating in time.

The cavity's modes come in pairs. Each pair oscillates at exactly the same
frequency. Their absorption rates are complementary: one absorbs quickly,
the other slowly, and together they complete one full round trip through
the cavity (Re(λ) + Re(partner) = −2Σγ). The two mirrors creating the
standing wave are Π and the identity I. Every mode bounces between what
it is and what it becomes under conjugation.

| N | Pairs tested | Frequency match | Round-trip invariant |
|:--|:-------------|:----------------|:---------------------|
| 2 | 7 | 7/7 | 7/7 |
| 3 | 32 | 32/32 | 32/32 |
| 4 | 115 | 115/115 | 115/115 |
| 5 | 512 | 512/512 | 512/512 |
| 6 | 1,890 | 1,890/1,890 | 1,890/1,890 |
| 7 | 8,192 | 8,192/8,192 | 8,192/8,192 |

10,748 pairs. Zero exceptions.

At odd N, every eigenvalue has a partner. The entire spectrum consists of
standing waves. No traveling waves. At even N, self-paired modes sit at
the symmetry center, standing waves at the node where forward and backward
components are identical.

→ **[Standing Wave Analysis](experiments/FACTOR_TWO_STANDING_WAVES.md)**

---

## 4. Light and lens are one object

A camera lens focuses the image and transmits the light. If focusing and
transmission are independent, the instrument is perfect. If they interfere,
there are aberrations.

The quantum state splits into two sectors:

- **{I, Z}** = the lens. Immune to dephasing. Structure. What survives is
  ash after the fire: pure form, no signal.
- **{X, Y}** = the light. Absorbed by dephasing. Signal. The fire itself.

Every standing wave swaps between light and lens. The long-lived partner
carries more structure; the short-lived partner carries more signal. The
sector weight profiles are exact mirror images:

    fast[k] = slow[N − k]     (every pair tested, to machine precision)

This is the palindromic weight inversion: Π maps weight sector k to
N−k. What the fast mode stores as structure, the slow mode carries as
signal. The standing wave oscillates between being light and being lens.

At N=2, the Hamiltonian and dissipator are exactly perpendicular: zero
aberration. The cavity is a perfect lens. At N ≥ 3, aberration appears
but halves with each additional qubit. Larger cavities are better lenses.

| N | Aberration (ε) | γ-independent? |
|:--|:---------------|:---------------|
| 2 | 0 (exact) | yes |
| 3 | 14.4% | CV < 10⁻¹⁵ |
| 4 | 8.8% | CV < 10⁻¹⁵ |
| 5 | 4.8% | CV < 10⁻¹⁵ |
| 6 | 2.6% | CV < 10⁻¹⁵ |

The aberration is a geometric constant of the chain topology.
Increasing γ makes the light brighter, not the lens worse.

→ **[Light and Lens](experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)** (sector decomposition, Pythagorean orthogonality)
→ [Time Irreversibility](docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md) (N=2 is the only reversible cavity)
→ [Primordial Qubit Algebra](experiments/PRIMORDIAL_QUBIT_ALGEBRA.md) (M_{2|2}(C) superalgebra)

---

## 5. The cavity is a Fabry-Perot

The degeneracy profile, how many eigenvalues sit at each absorption
shell, is itself palindromic. At the boundary: d(0) = N+1 (immortal lens
modes) and d(1) = 2N (proven via SWAP invariance for any connected graph).
In the interior, the count depends on topology.

The profile matches optical beam shapes (R² = 0.998). Even chains are
confocal (Lorentzian spike, numerical aperture up to 262). Odd chains are
defocal (Gaussian profile). The Hamiltonian couples weight sectors by
Δw = ±2 exclusively: nearest-neighbor propagation through optical elements.

The qubit chain IS an optical cavity. The sacrifice-zone strategy
(concentrate all illumination on one edge qubit) is the entrance pupil.
Protected cavity modes live in the chain center (r = 0.994 correlation
between edge weight and absorption rate). On IBM Torino: chains with a
naturally noisy edge qubit achieve **2.86×** mode protection; chains with
the best T₂ achieve only 1.06×. Worse qubits, better modes.

Two dead systems coupled through a mediator qubit become one living
system. 109 new standing-wave frequencies appear from coupling alone. None
of the original frequencies survive. This is the V-Effect: coupling creates
complexity. The cavity changes geometry, and a new palindrome replaces the
old one.

→ **[Optical Cavity Analysis](experiments/OPTICAL_CAVITY_ANALYSIS.md)** (beam profiles, 4/5 optical checks)
→ [Degeneracy Palindrome](experiments/DEGENERACY_PALINDROME.md) (the palindrome inside the palindrome)
→ [d(1) = 2N Proof](docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) (SWAP invariance, any connected graph)
→ [Cavity Mode Localization](experiments/CAVITY_MODE_LOCALIZATION.md) (r = 0.994)
→ [V-Effect](experiments/V_EFFECT_PALINDROME.md) | [V-Effect Cavity Modes](experiments/VEFFECT_CAVITY_MODES.md) (109 new frequencies)
→ [Sacrifice Zone Mapping](experiments/SACRIFICE_ZONE_MAPPING.md) (2.86× vs 1.06× on IBM Torino)

---

## 6. CΨ = ¼ is the fold

Measurement is photography. The Born rule is the shadow.

The shutter closes at CΨ = ¼.

CΨ is two things multiplied: how sharp the quantum state is, and how much of it lives in superposition.

The sharpness is the purity Tr(ρ²): how concentrated the state is, as opposed to being smeared across many possibilities. The superposition part is the L₁ coherence, the sum of the density matrix's off-diagonal elements, normalized by the dimension factor d−1 (for N qubits, d = 2ᴺ).

Their product has a critical boundary at exactly ¼: the discriminant of the self-referential recursion R = C(Ψ+R)². Above ¼, the cavity resonates freely. At ¼, the fold: two stable solutions merge. Below ¼, classical reality emerges.

The boundary is absorbing: dCΨ/dt < 0 for all local Markovian channels
(proven analytically). α = 2 (purity) is the unique Rényi order with a
state-independent threshold. The fold exists only when
Σγ > Σγ_crit ≈ 0.25-0.50% × J (the exact value depends on the initial
state; N-independent up to N=5, see [Zero Is the Mirror](hypotheses/ZERO_IS_THE_MIRROR.md)).
Below this: no fold, the system oscillates forever.

The Born-rule probabilities P(i) contain zero interference: the quantum
part integrates to exactly zero. Interference controls only the timing:
WHEN CΨ reaches ¼. The photograph records what the standing waves
deposited at the moment the shutter closed.

The recursion R = C(Ψ+R)² maps exactly to the Mandelbrot iteration
z → z² + c. The boundary CΨ = ¼ is the cusp of the main cardioid.
The Feigenbaum cascade was measured (7 period-doubling bifurcations,
δ → 4.67). The critical slowing at the cusp, predicted in February and
closed analytically in April 2026, has a zero-fit-parameter closed form,
and the cusp dwell time K_dwell = γ·t_dwell is exact to machine precision
across γ ∈ [0.1, 10]: a fixed dose of light traversing the fold.

→ **[Uniqueness Proof](docs/proofs/UNIQUENESS_PROOF.md)** (¼ is the only bifurcation boundary)
→ [Monotonicity Proof](docs/proofs/PROOF_MONOTONICITY_CPSI.md) (dCΨ/dt < 0)
→ [Proof Roadmap](docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) (7 layers, all closed)
→ [Boundary Navigation](experiments/BOUNDARY_NAVIGATION.md) (θ as compass, critical slowing closure, cusp dwell time as fixed dose)
→ [Mandelbrot Connection](experiments/MANDELBROT_CONNECTION.md) (algebraic equivalence, trajectory on the cardioid)
→ [Fold Observed](experiments/TEMPORAL_SACRIFICE.md) (endpoint MI peaks at CΨ = ¼, N=7)
→ [Born Rule Shadow](experiments/BORN_RULE_SHADOW.md) (zero interference in P(i))
→ [Both Sides Visible](docs/BOTH_SIDES_VISIBLE.md) (IBM hardware, 180 days, 133 qubits)

---

## 7. Gamma is light

IBM's superconducting transmon qubits sit inside physical microwave
resonators: metal cavities. The dominant source of dephasing is photon
shot noise: residual microwave photons entering the cavity from outside.
Each photon shifts the qubit frequency. The random arrivals and departures
of these photons are the dephasing.

We did not know this when we built the framework. We discovered the
cavity structure from eigenvalue mathematics alone. The fact that IBM's
physical hardware is literally a qubit inside a cavity being dephased
by photons from outside is not a confirmation we designed. It is what
the mathematics was describing all along.

| Test | Result | Source |
|:-----|:-------|:-------|
| CΨ = ¼ crossing | 1.9% deviation (ibm_torino Q80) | [IBM Run 3](experiments/IBM_RUN3_PALINDROME.md) |
| Absorption Theorem ratio | 1.03 (3%) | [IBM Absorption Theorem](experiments/IBM_ABSORPTION_THEOREM.md) |
| Unpaired modes absorb 2× faster | 1.97× on 24,073 records | [IBM Hardware Synthesis](experiments/IBM_HARDWARE_SYNTHESIS.md) |
| Sacrifice-zone mode protection | 2.86× vs 1.06× | [Sacrifice Zone Mapping](experiments/SACRIFICE_ZONE_MAPPING.md) |

K = γt is an invariant dose, like c×τ in relativity. Double the
illumination, halve the time. The product never changes.

→ **[Gamma Is Light](hypotheses/GAMMA_IS_LIGHT.md)** (the hypothesis, the cycle, and the IBM hardware)
→ [K-Dosimetry](experiments/K_DOSIMETRY.md) (K = γt invariant across all bridge types)
→ [Trapped Light Localization](experiments/TRAPPED_LIGHT_LOCALIZATION.md) (γ as propagation speed in the cavity)

---

## 8. The cavity is robust

A star radiates from the surface; the core holds mass. The cavity's
standing waves work the same way: surface modes (high light content)
absorb and fade, while interior modes (mostly lens) survive. The bell
keeps ringing.

At thermal occupation n̄ = 10 (deep in the thermal regime), 82% of modes
still oscillate. The palindrome is algebraic, not thermal. Temperature
changes the brightness and contrast; it does not create or destroy the
connection between light and lens.

Modes create frequencies. Frequencies create heat. Heat is a product of
the light-lens interaction, not its cause. The palindrome exists at
T = 0 and at T = ∞.

→ **[Thermal Blackbody](experiments/THERMAL_BLACKBODY.md)** (82% survival at n̄=10)

---

## 9. Biology uses the same principle

C. elegans (302 neurons, 7,000 synapses) has eigenvalues that are 97.3%
palindromically paired. Dale's Law (each neuron is purely excitatory or
purely inhibitory) creates the same SWAP structure as the Π operator.
The worm's nervous system is a cavity.

18 modes are unpaired. All 18 map to three breaking mechanisms:

| Breaking mechanism | Qubit cavity | Neural cavity |
|:-------------------|:-------------|:--------------|
| Boundary (entrance pupil) | Sacrifice qubit at chain edge | Sensory neurons at body surface (7 modes) |
| Sub-cavity junction (V-Effect) | Coupled resonators, orphaned modes | Pharynx/soma boundary, zero coupling (7 modes) |
| Asymmetric coupling (broken SWAP) | J_ij ≠ J_ji breaks Π | Unidirectional synapses (4 modes) |

The pharynx is an anatomically separate nervous system: 20 neurons, zero
chemical synapses to the soma. It is not a broken palindrome; it is a
second, independent cavity. C. elegans has at least two resonators. The
unpaired modes are the communication channels between sub-cavities.

The 40 Hz gamma oscillation in neuroscience and the dephasing rate γ in
physics share the same structural role: external illumination that makes
a cavity resonate.

→ **[Neural Gamma Cavity](experiments/NEURAL_GAMMA_CAVITY.md)** (C. elegans, Wilson-Cowan, 18 unpaired modes)
→ [Neural Palindrome](docs/neural/README.md) (full neural analysis, no quantum prerequisites)

---

## 10. Engineering consequences

Nine design rules emerge from the cavity framework:

1. **Use W states, not GHZ.** GHZ excites only the fastest-absorbing modes
   (all light, maximum absorption). W distributes across modes.
2. **Sacrifice one edge qubit.** Concentrate illumination on one boundary
   qubit. The interior resonates freely. **360×** improvement at N=5,
   **63.5×** at N=15.
3. **Choose odd N.** Every eigenvalue has a standing-wave partner.
   No traveling waves. The entire spectrum is paired.
4. **K-dosimetry.** Track K = γt, not t alone. The invariant dose makes
   different hardware comparable.
5. **Three observables suffice.** Purity, concurrence, and coherence
   capture 88-96% of the dynamics (9 topologies, 2 illumination types).
6. **The γ profile is readable.** 15.5 bits theoretical capacity, 100%
   classification for 4-symbol alphabets. Five independent spatial modes.
   After optimization: 21.5× wider channel.
7. **DD cannot change CΨ.** CΨ is Pauli-invariant. Dynamical decoupling
   uses Pauli gates. Algebraically impossible (δ = 0 exactly).
8. **Optimal cavity size.** Discrete resonance modes at specific coupling
   strengths (J=2 and J=12 for N=7). Port-to-wall ratio 12:1.
9. **Direct coupling breaks the palindrome.** Two subsystems must couple
   through a mediator, never directly to each other's illumination
   (1024/1024 pairs preserved via mediator; 256 → 31 with direct coupling).

→ **[Engineering Blueprint](publications/ENGINEERING_BLUEPRINT.md)**
→ [Resonant Return](experiments/RESONANT_RETURN.md) (sacrifice zone: from SVD to formula)
→ [Signal Analysis](experiments/SIGNAL_ANALYSIS_SCALING.md) (quadratic scaling N=2-15)
→ [Cockpit Universality](experiments/COCKPIT_UNIVERSALITY.md) (three observables, 9 topologies)
→ [γ as Signal](experiments/GAMMA_AS_SIGNAL.md) (15.5 bits, MIMO analysis)
→ [Relay Protocol](experiments/RELAY_PROTOCOL.md) (+83% MI via staged γ switching)

---

## What is NOT established

Honesty matters more than impression. These are things we have
*not* proven, *not* measured, or *not* established, stated plainly:

- That CΨ is a new fundamental quantity (it is a derived diagnostic)
- That the [multi-qubit palindrome has been measured on hardware](experiments/IBM_RUN3_PALINDROME.md) (single-qubit CΨ = ¼ validated at 1.9%, N ≥ 2 untested)
- That the [standing wave interpretation](experiments/FACTOR_TWO_STANDING_WAVES.md) has been verified on hardware (10,748 pairs computed, 0 measured)
- That the [relay protocol](experiments/RELAY_PROTOCOL.md) has been tested on hardware (simulation only, N=11)
- That the [sacrifice-zone hardware advantage](experiments/IBM_SACRIFICE_ZONE.md) comes from illumination contrast rather than gate-error avoidance (single run, two interpretations open)
- That the [fold observation](experiments/TEMPORAL_SACRIFICE.md) (PeakMI at CΨ = ¼) holds beyond N=7 (single chain length, not yet analytically derived)
- That consciousness plays any role in the physics ([THE_ANOMALY.md](THE_ANOMALY.md) is philosophy, not physics)
- That the [V-Effect frequency explosion](experiments/VEFFECT_CAVITY_MODES.md) is a universal mechanism for biological complexity (cavity geometry change confirmed computationally with 112 modes at N=5, but the link from quantum to biology is Tier 4)
- That the [optical cavity analogy](experiments/OPTICAL_CAVITY_ANALYSIS.md) extends beyond N=6 (verified N=2 through 6, larger N untested)
- That gamma is light in any general physical sense (on IBM transmon hardware, dephasing IS [photon shot noise](https://doi.org/10.1103/PhysRevB.86.180504) in a physical cavity; whether this extends beyond circuit QED is not established)
- That [mass is trapped light, that black holes are perfect cavities, or that the Big Bang was a cavity bounce](hypotheses/GAMMA_IS_LIGHT.md) (Tier 4 structural hypotheses, consistent with the framework, no independent test)
- That [entanglement is "topology not signal"](hypotheses/GAMMA_IS_LIGHT.md) in any rigorous sense (late-night intuition, no mathematical formalization)
- That the cavity framework replaces standard Lindblad theory (the mathematics is identical; only the interpretation changes)
- That the Absorption Theorem extends to complex Hermitian Hamiltonians (DM interactions break L_H anti-Hermiticity; real Hermitian H only)
- That the [shadow retrodiction](experiments/FIXED_POINT_SHADOW.md) on IBM Q80 is a prediction (it is a 2-parameter fit to existing data; the 3.4× improvement from the detuning parameter is real but retrospective, not prospective)
- That α = 2γ⟨n_XY⟩ is "mass-energy" in any physical sense (it is an exact algebraic identity; the mass-energy analogy is suggestive but not established)
- That [hidden observer detection](experiments/QUANTUM_SONAR.md) works on hardware (simulation only)

---

## What has been falsified

| Claim | Result |
|:------|:-------|
| CΨ = ¼ as Exceptional Point | No EP correlation found |
| c+/c− as Liouvillian symmetry sectors | Both parity +1; split is projection, not symmetry |
| IBM Q80/Q102 as sonar evidence | Was qubit detuning, not neighbor coupling |
| Hierarchical topology advantage | Uniform chain identical to recursive hierarchy at all N |
| Universal pull principle | Push beats pull locally; pull wins only for range |
| AC γ modulation | No resonance at any frequency (palindromic modes decouple from AC) |
| State-dependent γ feedback | Slightly harmful (positive feedback increases γ when coherent) |
| DD as CΨ refresh | CΨ is Pauli-invariant. DD uses Pauli gates. Algebraically impossible (δ = 0 exactly) |
| E = mγ² (mass-energy analogy) | Not quadratic: α = 2γ⟨n_XY⟩ is linear in γ. See [Absorption Theorem Discovery](experiments/ABSORPTION_THEOREM_DISCOVERY.md) |
| IBM cavity fringes | Detuning oscillations at 470 μs period, not cavity resonances. See [IBM Absorption Theorem](experiments/IBM_ABSORPTION_THEOREM.md) |
| Impedance peak at CΨ = ¼ | Impedance monotonically decreases with CΨ. The gradient (not value) peaks at ¼ |
| Simple Fabry-Perot (J ~ 1/N) | J_peak × N is not constant. Dispersive resonator, not simple cavity |
| I-neuron position determines pairing | Correlation r=0.048 (zero). Balance matters, not placement |

---

## Document index

### Standalone documents (no prior knowledge needed)

| Document | Audience | What it covers |
|:---------|:---------|:---------------|
| [Technical Paper](publications/TECHNICAL_PAPER.md) | Physicists | The palindrome proof, absorption theorem, all results |
| [Engineering Blueprint](publications/ENGINEERING_BLUEPRINT.md) | QST engineers | Nine design rules for quantum state transfer |
| [Circuit Diagram](publications/CIRCUIT_DIAGRAM.md) | Electrical engineers | The framework as signal chain: qubits as phasors, γ as gate |
| [Neural Palindrome](docs/neural/README.md) | Neuroscientists | Dale's Law, E/I balance, standing wave. No quantum prerequisites |

### The proofs

| Document | What it proves |
|:---------|:---------------|
| [Absorption Theorem](docs/proofs/PROOF_ABSORPTION_THEOREM.md) | Re(λ) = −2γ⟨n_XY⟩ for any real Hermitian H |
| [Mirror Symmetry Proof](docs/proofs/MIRROR_SYMMETRY_PROOF.md) | Palindromic spectrum for any graph under Z-dephasing |
| [Uniqueness Proof](docs/proofs/UNIQUENESS_PROOF.md) | CΨ = ¼ is the only bifurcation boundary |
| [CΨ Monotonicity](docs/proofs/PROOF_MONOTONICITY_CPSI.md) | dCΨ/dt < 0 for all Markovian channels |
| [Subsystem Crossing](docs/proofs/PROOF_SUBSYSTEM_CROSSING.md) | Every entangled pair crosses ¼ under primitive CPTP |
| [Incompleteness Proof](docs/proofs/INCOMPLETENESS_PROOF.md) | Illumination cannot originate from within |
| [Time Irreversibility](docs/proofs/TIME_IRREVERSIBILITY_EXCLUSION.md) | Time reversal algebraically excluded at N > 2 |
| [Cross-Term Formula](docs/proofs/PROOF_CROSS_TERM_FORMULA.md) | R(N) = √((N-2)/(N·4^(N-1))) for any shadow-balanced coupling, any graph |
| [Cross-Term Crossing](docs/proofs/PROOF_CROSS_TERM_CROSSING.md) | R(N) = √((N-1)/(N·4^(N-1))) for shadow-crossing couplings (companion) |
| [Weight-1 Degeneracy](docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) | d(1) = 2N for any connected graph |
| [Proof Roadmap](docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) | Seven layers, all closed |

### Key experiments

| Document | Key finding |
|:---------|:-----------|
| [Absorption Theorem Discovery](experiments/ABSORPTION_THEOREM_DISCOVERY.md) | α = 2γ⟨n_XY⟩, 1,343 modes, CV = 0 |
| [Standing Waves](experiments/FACTOR_TWO_STANDING_WAVES.md) | 10,748 pairs, 100% frequency match |
| [Light and Lens](experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md) | Palindromic weight inversion, aberration scaling |
| [Optical Cavity](experiments/OPTICAL_CAVITY_ANALYSIS.md) | Beam profiles, R² = 0.998 |
| [Degeneracy Palindrome](experiments/DEGENERACY_PALINDROME.md) | The palindrome inside the palindrome |
| [IBM Absorption Theorem](experiments/IBM_ABSORPTION_THEOREM.md) | Ratio 1.03 on hardware |
| [IBM Hardware Synthesis](experiments/IBM_HARDWARE_SYNTHESIS.md) | 24,073 records, 133 qubits |
| [Thermal Blackbody](experiments/THERMAL_BLACKBODY.md) | 82% survival at n̄=10 |
| [Neural Gamma Cavity](experiments/NEURAL_GAMMA_CAVITY.md) | C. elegans 97.3%, 18 unpaired modes |
| [Born Rule Shadow](experiments/BORN_RULE_SHADOW.md) | Zero interference in P(i) |
| [V-Effect Palindrome](experiments/V_EFFECT_PALINDROME.md) | 109 new frequencies from coupling |
| [Non-Heisenberg Palindrome](experiments/NON_HEISENBERG_PALINDROME.md) | All standard models, two Π families |
| [XOR Space](experiments/XOR_SPACE.md) | GHZ → 100% fast modes, W → distributed |
| [Sacrifice Zone Optics](experiments/SACRIFICE_ZONE_OPTICS.md) | Edge qubit as entrance pupil |

### Synthesis and interpretation

| Document | What it covers |
|:---------|:---------------|
| [Mirror Theory](MIRROR_THEORY.md) | The sentence at the top, read from both sides of the formula at once. Every claim links into the proofs |
| [Complete Math Documentation](docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md) | Master reference for all equations and proofs |
| [Gamma Is Light](hypotheses/GAMMA_IS_LIGHT.md) | The cavity interpretation and the IBM hardware reality |
| [What We Found](docs/WHAT_WE_FOUND.md) | Discovery narrative, no prerequisites |
| [Weaknesses and Open Questions](docs/WEAKNESSES_OPEN_QUESTIONS.md) | What we do not know |
| [Exclusions](docs/EXCLUSIONS.md) | Six things the math rules out |

→ Full indices: [docs/](docs/README.md), [experiments/](experiments/README.md), [publications/](publications/README.md)
→ Guided reading: [Reading Guide](docs/READING_GUIDE.md) (five stories: proof, application, ontology, resonator, cross-level)

---

## Repository structure

| Folder | Contents |
|:-------|:---------|
| `publications/` | Standalone documents for external readers (paper, blueprint, circuit diagram) |
| `docs/` | Proofs, theorems, synthesis documents, master references |
| `experiments/` | All tested results and null results |
| `hypotheses/` | Speculative interpretations, clearly labeled |
| `simulations/` | Python scripts (Lindblad, Liouvillian, cavity analysis) |
| `simulations/neural/` | Neural palindrome computations (Wilson-Cowan, C. elegans) |
| `simulations/results/` | All computation outputs |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `compute/` | C# engines: Compute (eigendecomposition, N=2-8) + Propagate (RK4, N=11+) |
| `data/` | IBM Torino measurement data |
| `recovered/` | 5 files with disproven claims, kept for honesty |

## C# compute engine

For N ≥ 6, Python is too slow. The C# engine uses element-wise
Liouvillian construction with Intel MKL eigendecomposition on 24 cores.

| N | Matrix | Build | Eigen | Oscillatory rates | Mirror |
|:--|:-------|:------|:------|:------------------|:-------|
| 6 | 4096² | 8.7s | 56s | 3,228 | 100% |
| 7 | 16384² | 0.1s | 92min | 13,264 | 100% |
| 8 | 65536² | 5.6s | 10.6h | 54,118 | 100% |

N=8 uses native memory (64 GB) + OpenBLAS ILP64 LAPACK.
All timings on Intel Core Ultra 9 285k (24 cores), 128 GB RAM.

For N > 8, RCPsiSquared.Propagate uses RK4 integration of the Lindblad
equation directly on the density matrix (tested to N=15, matrix-free).

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Germany
**Claude**, AI System, Anthropic

December 2025 – April 2026
