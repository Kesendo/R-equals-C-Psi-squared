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

Later discoveries under the [γ₀ = const](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) framework (γ₀ as a framework constant, not a tunable parameter) sharpened the operational consequence: Alice chooses her initial state (receiver engineering), not her noise profile. The F67 bonding modes form Alice's operationally complete receiver menu. Against the classic ENAQT-optimised transport baseline, bonding-mode receivers reach **4000-5500×** improvement in simulation. The receiver-engineering signature was confirmed live on IBM Kingston (Heron r2) in April 2026, with the noise-robustness gap matching framework prediction (see Section 11 for numbers and discussion).

In April 2026 we closed the chain from below: the Heisenberg coupling is not postulated, it is **forced** by the C²⊗C² parity structure of the Pauli algebra at d=2. The Level 0 → Level 1 V-Effect bridge produces the textbook atomic exchange with predicted prefactor δE_GS = −3α² / (4(J_A + J_B)), an Anderson-superexchange shape derived end-to-end from Pauli algebra alone. The (w=0, w=N) extreme sectors of the palindromic relation are immune to every 2-body Hamiltonian (analytical proof). The 14-of-36 V-Effect break decomposes into a 3 / 19 / 14 fine structure invisible to spectroscopy. All of this is operationalised in the [framework package](simulations/framework/): ur-Pauli indexed by (bit_a, bit_b), ur-Heisenberg as the unique both-parity-even bilinear, ur-eigenvalues from the Pauli identity (σ·σ)² = 3I − 2(σ·σ). The 3/19/14 fine structure was verified live on `ibm_marrakesh` (Heron r2) on 2026-04-26: discriminating Pauli expectation ⟨X₀Z₂⟩ = +0.011 (truly), **-0.711** (soft), +0.205 (hard), all three categories separately resolved at 13-47σ. See Section 12.

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

→ **[Heisenberg Reloaded](hypotheses/HEISENBERG_RELOADED.md)** + **[framework package](simulations/framework/)**: how the Heisenberg form is forced from below by Pauli algebra at d=2, and the V-Effect bridge that produces atomic exchange with quantitative prefactor. Section 12 below.

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
1,342 modes (N=2-5, γ=0.01-1.0, J=0.1-5.0), coefficient of variation = 0.

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

For shadow-balanced couplings at N=2, the Hamiltonian and dissipator are
exactly perpendicular: zero aberration, the cavity is a perfect lens.
This is unique. For all other cases, aberration appears but stays
geometrically determined. Two coupling classes have closed-form
aberration:

- **Shadow-balanced** (both bond Paulis in {X,Y} or both in {I,Z}):
  Heisenberg, XXZ, XY, Ising, Dzyaloshinskii-Moriya. Both bond sites at
  the same shadow depth. Bond-site variance contributes 0.

      ε = √((N−2) / (N·4^(N−1)))

- **Shadow-crossing** (one bond Pauli in {X,Y}, the other in {I,Z}):
  X_iZ_j, Y_iZ_j, including the native cross-resonance drive on
  superconducting hardware. One bond site in light, one in shadow.
  Bond-site variance contributes 1.

      ε = √((N−1) / (N·4^(N−1)))

The two formulas differ by exactly one unit of variance: the asymmetry
between bond sites at different shadow depths. Both are γ-independent
geometric constants of the cavity, not physical parameters.

| N | Balanced ε | Crossing ε |
|:--|:-----------|:-----------|
| 2 | 0 (exact) | 1/√8 ≈ 35.4% |
| 3 | 1/√48 ≈ 14.4% | 1/√24 ≈ 20.4% |
| 4 | 1/√128 ≈ 8.8% | √(3/256) ≈ 10.8% |
| 5 | √(3/1280) ≈ 4.8% | 1/√320 ≈ 5.6% |
| 6 | 1/√1536 ≈ 2.6% | √(5/6144) ≈ 2.9% |

Both shrink exponentially with N. Larger cavities are better lenses,
regardless of which coupling class drives them. Increasing γ makes the
light brighter, not the lens worse.

→ **[Light and Lens](experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)** (sector decomposition, Pythagorean orthogonality)
→ [Cross-Term Formula](docs/proofs/PROOF_CROSS_TERM_FORMULA.md) (shadow-balanced, proven)
→ [Cross-Term Crossing](docs/proofs/PROOF_CROSS_TERM_CROSSING.md) (shadow-crossing companion)
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

(Under [γ₀ = const](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), γ-modulation
like the sacrifice zone is not operationally controllable by Alice. The
superseding result: receiver engineering via F67 bonding modes achieves
the same kind of advantage with 11-15× higher absolute magnitude, at
uniform γ. See [Section 11](#11-the-receiver-menu).)

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
→ [We Are the Fragment](hypotheses/WE_ARE_THE_FRAGMENT.md) (¼ as double fragmentation: operator C=1/2 × spatial sublattice 1/2)

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
| Receiver engineering (bonding:2 vs alt-z-bits) | 2.80× on ibm_kingston Heron r2, 2 QPU min | [IBM Receiver Engineering](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) |

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
2. **Choose an F67 receiver, not a γ profile.** Under γ₀ = const Alice
   cannot modulate γ per qubit; instead she picks her initial state from
   the F67 eigenmode menu. Bonding:k beats alt-bit transport by factor
   1.4× to 4.6× growing with N (ideal sim, tested N=5 to 13), and
   **4000-5500× over the ENAQT baseline** at uniform J. On real hardware
   (ibm_kingston, April 2026): 2.80× bonding:2 over alt-z-bits. The older
   sacrifice-zone lever (concentrate γ on one edge qubit, **360×** at N=5,
   **63.5×** at N=15) is superseded by receiver engineering by 11-15×
   and remains valid only when γ IS operationally controllable.
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

→ [Resonant Return](experiments/RESONANT_RETURN.md) (sacrifice zone: from SVD to formula)
→ [Signal Analysis](experiments/SIGNAL_ANALYSIS_SCALING.md) (quadratic scaling N=2-15)
→ [Cockpit Universality](experiments/COCKPIT_UNIVERSALITY.md) (three observables, 9 topologies)
→ [γ as Signal](experiments/GAMMA_AS_SIGNAL.md) (15.5 bits, MIMO analysis)
→ [Relay Protocol](experiments/RELAY_PROTOCOL.md) (+83% MI via staged γ switching)

---

## 11. The receiver menu

Under [γ₀ = const](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) (γ₀ is a framework
constant like c, not a hardware parameter), Alice cannot tune the noise. She
can only choose her initial state. The F67 single-excitation eigenmodes of
the uniform-J Heisenberg chain form Alice's operational menu:

    |ψ_k⟩ = √(2/(N+1)) Σⱼ sin(πk(j+1)/(N+1)) |1_j⟩,   k = 1, ..., N

Each k optimises a different transport task:

| Task | Best receiver | Why |
|:-----|:-------------|:----|
| Distributed correlation (Sum-MI) | alt-z-bits \|01010..⟩ | Multi-k amplitude spread drives adjacent-pair buildup |
| End-to-end MI(0, N-1) | bonding:k* = ⌊(N+1)/3⌋ | Dynamic optimum: low k with boost from Heisenberg hopping |
| Multi-end Mirror-Pair MM | any even k (F75 static) | All probability mass concentrated on mirror-pairs |
| Long-time memory | bonding:1 | Slowest decay α_1 ≈ π²γ₀/(N+1)³ |

A closed-form analytical model covers state preparation through observed
envelope: **F65** (eigenmode formula) + **F75** (t=0 mirror-pair MI) + **F76**
(0.93 dephasing envelope at γ₀·t=0.005). PeakMM = 0.93 × MM(0) universally
within 0.5% across 25+ (N, k) tested points at N=5 to 13. No propagation
needed for prediction.

The advantage grows superlinearly with N:

| Metric | N=5 | N=7 | N=9 | N=11 | N=13 |
|:-------|:---:|:---:|:---:|:----:|:----:|
| MI(0, N-1) bonding/alt-z-bits | 1.39× | 1.48× | 2.02× | 3.02× | **4.59×** |
| MM bonding/alt-z-bits | 0.99× | 1.37× | 1.88× | 3.26× | **5.25×** |

Live IBM Kingston (Heron r2) confirmation, April 2026: **2.80× bonding:2 vs
alt-z-bits** at N=5 in 2 QPU minutes. Advantage grows with noise level
because bonding:2's single delocalised excitation is more T1-robust than
alt-z-bits' two localised excitations.

Sacrifice-zone γ-modulation (Section 5) was the predecessor. It achieved
139-360× over ENAQT by breaking \|+⟩^N's J-blindness via asymmetric γ.
Receiver engineering reaches 4000-5500× over ENAQT at uniform γ₀ by using
a J-sensitive receiver directly. Both tools are in the kit; receiver
engineering is the larger lever.

→ **[Receiver vs γ-Sacrifice](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md)** (the full story, N=5 to 13)
→ [Primordial γ Constant](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) (γ₀ as framework invariant)
→ [Analytical Formulas F65, F67, F75, F76](docs/ANALYTICAL_FORMULAS.md)
→ [IBM Receiver Engineering](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) (Kingston Run 1, 2.80× ratio)
→ [PROOF K-Partnership](docs/proofs/PROOF_K_PARTNERSHIP.md) (bipartite sublattice gauge K = diag((-1)^l) gives spectral inversion KHK = -H)
→ [IBM K-Partnership Marrakesh](experiments/IBM_K_PARTNERSHIP_SKETCH.md) (5 receivers in 5 QPU min, K-pair spread reads γ-profile asymmetry)
→ [Z⊗N Partnership](experiments/Z_N_PARTNERSHIP.md) (multi-exc Néel mirror as third diagnostic channel, transverse-field detection)
→ [Q Scale Three Bands](experiments/Q_SCALE_THREE_BANDS.md) (Q = J/γ₀ structure, Q_peak(c) saturation at 1.8)
→ [J-Blind Receiver Classes](experiments/J_BLIND_RECEIVER_CLASSES.md) (why \|+⟩^N was the wrong receiver)
→ [Between Measurements Evidence](hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md) (structural argument for γ₀ = const)

---

## 12. Heisenberg from below: V-Effect creates the exchange

The Heisenberg coupling J σ_1 · σ_2 = J(XX + YY + ZZ) is not postulated. It is forced by the C²⊗C² parity structure of the Pauli algebra at d=2. Each single-qubit Pauli is indexed by two Z₂ parities: bit_a (n_XY, dephasing axis: I, Z immune; X, Y decaying) and bit_b (n_YZ, Π² axis: I, X even; Y, Z odd). The intersection of both-parity-even 2-body operators is exactly **{II, XX, YY, ZZ}**. Heisenberg/XXZ is the unique 2-body bilinear that survives both Z₂ filters; every other form (XY, ZX, etc.) breaks at least one parity.

When two atoms (each represented by an internal Heisenberg pair with strength J_A, J_B) are bonded through a V-Effect bridge α on a single inter-pair bond, second-order perturbation theory predicts the effective Level-1 exchange:

    δE_GS = − 3 α² / (4 (J_A + J_B))

verified numerically at N=4 to 0.2-0.8 % across seven asymmetric J_A:J_B ratios from 1:1 to 1:10. The "3" is universal: from the Pauli identity (σ·σ)² = 3 I − 2 (σ·σ) on a singlet-singlet ground state where ⟨σ·σ⟩_bridge = 0. The "4(J_A + J_B)" is the cost to flip both pairs simultaneously from singlet to triplet. Symmetric limit J_A = J_B = J recovers −3α²/(8J), the Anderson-superexchange shape applied to direct Heisenberg bridges.

The (w=0, w=N) extreme weight sectors of the palindromic relation are immune to **every** 2-body Hamiltonian (analytical proof). Z-dephasing vanishes on w=0 strings (Z commutes with I, Z); the Hamiltonian commutator's (w=0, w=0)-block is identically zero (only ZZ preserves w=0 among 2-body terms, and ZZ commutes with all w=0 strings). The break, when present, is strictly localised to boundary sectors 0 < w < N.

Re-examining the V-Effect's 14-of-36 finding via framework primitives reveals a **3 / 19 / 14** fine structure that the original spectral test obscured. Of the 22 V_EFFECT_PALINDROME-unbroken cases: 3 are truly unbroken (the Heisenberg/XXZ subset), 19 are "soft-broken" (operator equation residual non-zero but eigenvalue pairing intact at machine precision), and 14 hard-break both criteria.

The 3/19 distinction is verified at simulation level via eigenvector pairing: for each L eigenvalue λ_i, find its partner λ_j ≈ −λ_i − 2Σγ in the spectrum, then compute |⟨v_j | Π · v_i⟩| / (‖v_j‖ · ‖Π v_i‖). For the 3 truly-unbroken: overlap = 1.000000 (machine precision). For the 19 soft-broken: overlap = 0.000 to 0.598 (severely broken; many cases land on a subspace fully orthogonal to the partner eigenspace). The spectrum lies: it shows palindromic frequencies while the underlying eigenvector structure is scrambled.

**Live verification on IBM Marrakesh (2026-04-26):** the 3/19/14 distinction is operationally observable on real Heron r2 hardware. Three representative Hamiltonians (XX+YY truly, XY+YX soft, XX+XY hard) on path [48, 49, 50] from |+−+⟩ initial state at t=0.8, n_trotter=3, 4096 shots/basis. Discriminating Pauli expectation ⟨X₀Z₂⟩: +0.011 (truly), **-0.711** (soft), +0.205 (hard). Discrimination Δ(soft − truly) = **-0.72**, stronger than the continuous-Lindblad idealization -0.62; the difference is fully explained by Trotter n=3 discretization at δt=0.267 (matching prediction to 0.0014; T1 actually attenuates, see [MARRAKESH_THREE_LAYERS](experiments/MARRAKESH_THREE_LAYERS.md)). All three categories separately resolved at 13-47σ. Job `d7mjnjjaq2pc73a1pk4g`. See [`data/ibm_soft_break_april2026/`](data/ibm_soft_break_april2026/README.md). This is the FIRST framework-driven hardware test (the IBM pipeline script `run_soft_break.py` is the first in that directory to import the framework package).

→ **[Heisenberg Reloaded](hypotheses/HEISENBERG_RELOADED.md)** (the level-stack inheritance picture: Level 0 → Level 1 via V-Effect)
→ [Proof: Zero Immunity](docs/proofs/PROOF_ZERO_IMMUNITY.md) ((w=0, w=N) palindromic immunity for any 2-body H, four lemmas)
→ [V-Effect Boundary Localization](experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md) (numerical verification at N=3 and N=4 to machine precision)
→ [V-Effect Fine Structure](experiments/V_EFFECT_FINE_STRUCTURE.md) (3+19+14 decomposition of the 14-of-36 result, verified via eigenvector pairing in `simulations/_soft_break_eigenvector_test.py`)
→ [Exchange from V-Effect](experiments/EXCHANGE_FROM_V_EFFECT.md) (−3α²/(8J) symmetric, derived end-to-end from Pauli algebra)
→ [Asymmetric Exchange](experiments/ASYMMETRIC_EXCHANGE_FROM_V_EFFECT.md) (−3α²/(4(J_A+J_B)), first calculation built on framework primitives)
→ [On the Soft Break](reflections/ON_THE_SOFT_BREAK.md) (what the framework reveals beyond V_EFFECT_PALINDROME's spectral test)
→ [framework package](simulations/framework/) (ur-Pauli, ur-Heisenberg, ur-eigenvalues, palindrome residual, V-Effect primitives in the 4-layer cockpit)

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
- That the [receiver-engineering advantage extends to N > 13 in simulation](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) (tested N=5 to 13 via dense + matrix-free brecher; N=15+ would need further pipeline work)
- That the [receiver-engineering advantage extends beyond N=5 on hardware](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) (Run 1 is N=5 only; N=7 and 9 hardware tests budgeted but not yet run)
- That [γ₀ can be measured from inside the framework](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) (it cannot; only the ratio Q = J/γ₀ is an intrinsic observable; any absolute γ₀ value requires an external unit convention)
- That [bonding:k preparation is gate-optimal](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) (current Qiskit StatePreparation uses ~30 two-qubit gates at N=5; custom Dicke+phase circuits could be cheaper)

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
| Linear Q_peak(c) = 1.3 + 0.15·c | Q_peak saturates at 1.8 for c ≥ 4. Pattern: {1.5, 1.6, 1.8, 1.8} for c = {2, 3, 4, 5}. See [Q Scale Three Bands](experiments/Q_SCALE_THREE_BANDS.md) |
| Receiver advantage shrinks fast with N | The opposite: grows superlinearly. 1.39× → 4.59× across N=5..13 for MI(0, N-1). See [Receiver vs γ-Sacrifice](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) |
| Best-J plateau at ~3.3 Peak Sum-MI | Python coarse-grid artifact. C# fine-grid brecher shows linear growth to 7.07 at N=9 |
| Python shadow_lens_broken.py uniform-J baselines | np.linspace(0.1, 15, 40) missed peaks at t~0.24; N=5 MI values were factor ~2 too low. Fixed via C# fine-grid brecher mode |

---

## Document index

### Standalone documents (no prior knowledge needed)

| Document | Audience | What it covers |
|:---------|:---------|:---------------|
| [What We Found](docs/WHAT_WE_FOUND.md) | Anyone | Discovery narrative, no prerequisites |
| [Reading Guide](docs/READING_GUIDE.md) | Anyone | Five stories (proof, application, ontology, resonator, cross-level) |
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
| [Zero Immunity](docs/proofs/PROOF_ZERO_IMMUNITY.md) | (w=0, w=N) palindrome holds for every 2-body H, regardless of parity-violation |
| [K-Partnership](docs/proofs/PROOF_K_PARTNERSHIP.md) | Bipartite sublattice gauge K = diag((-1)^l) gives KHK = -H spectral inversion |
| [Proof Roadmap](docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) | Seven layers, all closed |

### Key experiments

| Document | Key finding |
|:---------|:-----------|
| [Absorption Theorem Discovery](experiments/ABSORPTION_THEOREM_DISCOVERY.md) | α = 2γ⟨n_XY⟩, 1,342 modes, CV = 0 |
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
| [Receiver vs γ-Sacrifice](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) | F67 bonding-mode beats γ-modulation 11-15×, 4000-5500× over ENAQT, validated through N=13 |
| [Q Scale Three Bands](experiments/Q_SCALE_THREE_BANDS.md) | Q = J/γ₀ three-band structure; Q_peak saturates at 1.8 for c ≥ 4 |
| [J-Blind Receiver Classes](experiments/J_BLIND_RECEIVER_CLASSES.md) | Three mechanism classes of J-blindness (DFS, H-degenerate, M_α-polynomial) |
| [IBM Receiver Engineering](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md) | Kingston Heron r2 Run 1: bonding:2 / alt-z-bits = 2.80× |
| [IBM K-Partnership](experiments/IBM_K_PARTNERSHIP_SKETCH.md) | Marrakesh Heron r2: 5 receivers in 5 QPU min, K-pair spread reads γ-profile |
| [IBM Soft-Break](data/ibm_soft_break_april2026/README.md) | Marrakesh Heron r2: 3/19/14 fine-structure hardware-verified, Δ⟨X₀Z₂⟩ = -0.72 (~50σ) |
| [Z⊗N Partnership](experiments/Z_N_PARTNERSHIP.md) | Multi-exc Néel mirror as third diagnostic channel (transverse-field detection) |
| [V-Effect Boundary Localization](experiments/V_EFFECT_BOUNDARY_LOCALIZATION.md) | Numerical verification: extreme sectors immune, breaks confined to boundary |
| [V-Effect Fine Structure](experiments/V_EFFECT_FINE_STRUCTURE.md) | 14-of-36 result decomposes into 3 truly unbroken + 19 soft-broken + 14 hard-broken |
| [Exchange from V-Effect](experiments/EXCHANGE_FROM_V_EFFECT.md) | Effective Level-1 exchange δE = -3α²/(8J) derived end-to-end from Pauli algebra |
| [Asymmetric Exchange](experiments/ASYMMETRIC_EXCHANGE_FROM_V_EFFECT.md) | Generalisation -3α²/(4(J_A+J_B)), first calculation built on the framework package |

### Synthesis and interpretation

| Document | What it covers |
|:---------|:---------------|
| [Mirror Theory](MIRROR_THEORY.md) | The sentence at the top, read from both sides of the formula at once. Every claim links into the proofs |
| [Complete Math Documentation](docs/proofs/COMPLETE_MATHEMATICAL_DOCUMENTATION.md) | Master reference for all equations and proofs |
| [Heisenberg Reloaded](hypotheses/HEISENBERG_RELOADED.md) | The Level 0 → Level 1 inheritance: Heisenberg form forced by Pauli C²⊗C², V-Effect bridges to atomic exchange |
| [We Are the Fragment](hypotheses/WE_ARE_THE_FRAGMENT.md) | 1/4 fold as double fragmentation: operator-level C=1/2 × spatial sublattice 1/2 |
| [On the Soft Break](reflections/ON_THE_SOFT_BREAK.md) | What the framework reveals beyond V_EFFECT_PALINDROME's spectral test |
| [Gamma Is Light](hypotheses/GAMMA_IS_LIGHT.md) | The cavity interpretation and the IBM hardware reality |
| [Primordial Gamma Constant](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) | γ₀ as framework constant, Q = J/γ₀ as only measurable ratio |
| [Between Measurements Evidence](hypotheses/BETWEEN_MEASUREMENTS_EVIDENCE.md) | Structural argument for γ₀ = const |
| [What We Found](docs/WHAT_WE_FOUND.md) | Discovery narrative, no prerequisites |
| [Weaknesses and Open Questions](docs/WEAKNESSES_OPEN_QUESTIONS.md) | What we do not know |
| [Exclusions](docs/EXCLUSIONS.md) | Six things the math rules out |

→ Full indices: [docs/](docs/README.md), [experiments/](experiments/README.md)
→ Guided reading: [Reading Guide](docs/READING_GUIDE.md) (five stories: proof, application, ontology, resonator, cross-level)

---

## Repository structure

| Folder | Contents |
|:-------|:---------|
| `docs/` | Proofs, theorems, synthesis documents, master references |
| `experiments/` | All tested results and null results (~140 entries) |
| `hypotheses/` | Speculative interpretations, clearly tier-labeled |
| `reflections/` | Synthesis arcs across hypotheses; reading-mode threads (~19 entries) |
| `simulations/` | Python view-layer. `framework/` package (4-layer cockpit: primitives, entities, diagnostics, workflows; ~25 files). `framework_archive.py` (the original monolithic version, kept as repository memory). `neural/` (Wilson-Cowan, C. elegans). `water/` (hydrogen-bond qubit). Plus ~410 standalone scripts (one-shot exploration; `_`-prefix indicates work-in-progress) |
| `simulations/results/` | All computation outputs (~80 sub-folders per experiment) |
| `simulations/app/` | Five Regulator Simulator (Streamlit) |
| `compute/` | C# layered architecture (see "C# layered architecture" section below): Core (typed F-Claims, Π², Pi2 foundation, BlockSpectrum, SymmetryFamily; active development front since 2026-04-30) + Compute (eigendecomp N=2-8) + Propagate (RK4 + matrix-free, N up to 15) + Diagnostics + Cli + Orchestration + Runtime + Visualization, each with their Tests project |
| `data/` | IBM Torino, Marrakesh, Kingston measurement data |
| `visualizations/` | Figures and plots, incl. Mandelbrot trajectory overlays |
| `review/` | Internal audit + emerging-questions tracking (`EMERGING_QUESTIONS`, `EQ###_FINDINGS`, `OPEN_QUESTIONS_INDEX`) |
| `recovered/` | "Premature, not wrong" + unsupported + early-exploration entries (6 files, kept for honesty about the research process; some like `LIGHT_FIRST_FREE_MIRRORING` were ahead of the math and later confirmed) |

## C# layered architecture

The `compute/` directory grew from a single eigendecomposition backend into a layered architecture. Three layers, distinct roles:

### C# Core (typed-knowledge layer; active development front since 2026-04-30)

`compute/RCPsiSquared.Core/` is the typed home for the project's structural results. Each F-formula, Pi2 foundation anchor, symmetry-family member, and Π² primitive is a typed Claim with a Tier label (`Tier1Derived`, `Tier1Candidate`, `Tier2Verified`, etc.), dependencies on other Claims, and inspectable properties. The aggregator `Pi2KnowledgeBase` exposes 15 Tier-1 derived claims under `--root pi2`; `F86KnowledgeBase`, `BlockSpectrum`, `SymmetryFamily` are sister knowledge bases. Live introspection via `dotnet run --project compute/RCPsiSquared.Cli -- inspect --root pi2` (or any other root).

Recent additions (2026-05): `UniversalCarrierClaim` (γ₀ as protection-constant + observation-substrate in one role; `DefaultGammaZero = 0.05` is "our c"); `n=0` doppelte Rolle on `Pi2DyadicLadderClaim` (a₀ = 2 carries both qubit-dimension and absorption-quantum-coefficient readings); F91/F92/F93 anti-palindromic γ/J/h symmetries; XGlobalChargeConjugationPairing in BlockSpectrum.

### C# Compute (eigendecomposition, N=2-8)

For N ≥ 6, Python is too slow. The Compute engine uses element-wise
Liouvillian construction with Intel MKL eigendecomposition on 24 cores.

| N | Matrix | Build | Eigen | Oscillatory rates | Mirror |
|:--|:-------|:------|:------|:------------------|:-------|
| 6 | 4096² | 8.7s | 56s | 3,228 | 100% |
| 7 | 16384² | 0.1s | 92min | 13,264 | 100% |
| 8 | 65536² | 5.6s | 10.6h | 54,118 | 100% |

N=8 uses native memory (64 GB) + OpenBLAS ILP64 LAPACK.
All timings on Intel Core Ultra 9 285k (24 cores), 128 GB RAM.

### C# Propagate (RK4 / matrix-free, N up to 15)

For N > 8, RCPsiSquared.Propagate uses RK4 integration of the Lindblad
equation directly on the density matrix (tested to N=15, matrix-free).
The `brecher` sub-mode adds F67 receiver-engineering evaluation with
arbitrary `bonding:k` or bit-pattern initial states (commits `0917038`,
`865641c`, `8caf499`, `e1ee822`). Matrix-free path handles N ≥ 13 for
brecher mode (commit `a5b347d`).

### Three-layer hierarchy (Markdown, Python, C#)

The framework itself lives in Markdown (`docs/`, `hypotheses/`, `reflections/`, `experiments/`); Python (`simulations/framework/`) and C# (`compute/RCPsiSquared.Core/`) are view-layers operationalising what the Markdown describes. The Python framework was the primary active layer through April 2026 (V-Effect, soft-break, framework primitives); since 2026-04-30 the typed C# Core has become the active development front. Both view-layers remain in use; the Python package is stable but not actively extended.

---

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Germany
**Claude**, AI System, Anthropic

December 2025 – April 2026
