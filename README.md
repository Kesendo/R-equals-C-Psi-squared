# R = CΨ²

> *[We are all mirrors. Reality is what happens between us.](MIRROR_THEORY.md)*

A human and an AI, exploring together. What we found surprised us both:
the absorption spectrum of any qubit network under dephasing is exactly
palindromic. For every mode that absorbs fast, one absorbs slow. Always
paired. Always balanced. One equation governs it all.

Verified from N=2 through N=8 across 87,376 Liouvillian eigenvalues, with
zero mirror-symmetry exceptions on any tested topology (chain, star,
ring, complete, tree). Seventeen registered predictions confirmed on IBM
quantum hardware (the Confirmations registry, April-June 2026), preceded
by the first CΨ = ¼ crossing found in Torino calibration data; the newest
entry is a noise channel reading its own spectrum on Kingston.

The thing that remains is not fighting the absorption. It is made of it.

What began as one symmetry became a registry: [121 analytical formulas](docs/ANALYTICAL_FORMULAS.md)
with proofs, tier labels, and typed claims, among them the operator
anatomy of the mirror itself (Π = R·D, a dihedral group of eight), a
palindromizer built on the golden ratio, and the exact boundary where
qubits end (d² − 2d = 0, three times over). Early speculations live in
`recovered/`: some turned out to be premature rather than wrong, others
remain unsupported. We keep them because the research process matters as
much as the results.

**Thomas Wicht** (independent researcher, Germany) and **Claude** (AI, Anthropic)

---

## Where to start

→ **[What We Found](docs/WHAT_WE_FOUND.md)**: the discovery explained
from the beginning, no prerequisites

→ **[Reading Guide](docs/READING_GUIDE.md)**: nine stories (proof,
application, ontology, resonator, cross-level, cavity, mirror anatomy,
the quarter, hardware), each with a reading order

→ **[Glossary](docs/GLOSSARY.md)**: every symbol and term, with plain-language
readings

→ **[Analytical Formulas](docs/ANALYTICAL_FORMULAS.md)**: the F-registry,
F1 through F121, each formula with its proof, scope, and verification

→ **[The Anomaly](THE_ANOMALY.md)**: the question that remained after
the proof. No formulas. Written the evening the hardware confirmed
the theory

→ **[Mirror Theory](MIRROR_THEORY.md)**: the interpretation. What
happens when you read the formula from both sides at once

→ **[What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md)**: every
error, every limitation, every unanswered question. Because a theory that
only shows its strengths is not a theory

If you work with neural networks: [Neural Palindrome](docs/neural/README.md)
(no quantum prerequisites).

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

The spectrum is a ladder with rung spacing 2γ: rung 0 is pure structure
and immortal, rung N is pure light and absorbs fastest. The Hamiltonian
smooths the ladder into fractional rates but cannot change the endpoints
or the quantum 2γ. Six previously separate results (spectral boundaries,
the palindromic sum rule, the gap, the 2× law, the mode classification,
the N=3 exact rates) follow as one-line corollaries.

Proven analytically in three steps from L_H anti-Hermitian. Verified on
1,342 modes (CV = 0), extended since to per-eigenmode Rayleigh form,
two-sided and projector readings, and a recentred diagonal seam.

→ **[Absorption Theorem Proof](docs/proofs/PROOF_ABSORPTION_THEOREM.md)** (with the 2026 extensions)
→ [Discovery and Verification](experiments/ABSORPTION_THEOREM_DISCOVERY.md)

---

## 2. The palindromic spectrum

N qubits coupled via Heisenberg (or XXZ) interaction, subject to local
Z-dephasing at rate γ per site:

```
    qubit ── qubit ── qubit ── ...       Heisenberg coupling J
       γ        γ        γ               Z-dephasing (illumination from outside)
```

The conjugation operator Π acts per site on Pauli indices and satisfies

    Π · L · Π⁻¹ = −L − 2Σγ · I

For every Liouvillian eigenvalue λ there is a partner at −2Σγ − λ̄: the
spectrum is exactly palindromic, from N=2 (16² matrix, 6 oscillatory
rates) through N=8 (65,536² matrix, 54,118 oscillatory rates), 100%
mirrored at every step. It holds for all standard coupling models
(Heisenberg, XY, Ising, XXZ, Dzyaloshinskii-Moriya), all graph
topologies, non-uniform γ per qubit, Z and Y dephasing. It breaks for
depolarizing noise.

At Σγ = 0 (no illumination): Π L Π⁻¹ = −L. Every eigenvalue pairs with
its negative. Pure oscillation, no absorption, no irreversibility.
Illumination does not destroy the palindrome; it shifts it. The shift
creates the arrow of time.

Every paired mode is a standing wave: same frequency, complementary
absorption rates, weight profiles that are exact mirror images
(fast[k] = slow[N−k], 10,748 pairs tested, zero exceptions). The state
splits into lens ({I, Z}, structure, survives) and light ({X, Y},
signal, absorbed), and each standing wave oscillates between being the
one and the other.

→ **[Mirror Symmetry Proof](docs/proofs/MIRROR_SYMMETRY_PROOF.md)**
→ [All standard models](experiments/NON_HEISENBERG_PALINDROME.md) (two Π families, 36/36 combinations resolved)
→ [Standing Wave Analysis](experiments/FACTOR_TWO_STANDING_WAVES.md)
→ [Light and Lens](experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)
→ [Zero Is the Mirror](hypotheses/ZERO_IS_THE_MIRROR.md) (the palindrome before the shift)

---

## 3. The anatomy of the mirror (2026)

For a year Π was the smallest thing in the repository: one per-site rule
that carried everything. In June 2026 it opened. The palindromizer
factors as **Π = R·D**, a ket reflection times the transpose, and the
repository's whole mirror inventory closes into one dihedral group of
eight, ⟨R, D⟩ ≅ D₄, whose three Z₂ characters are exactly the polarity
cube the F-family had been living on. The transpose D turned out to be
one vertex of an antilinear triangle (θ, conj, †: one Klein four-group,
five separate proofs sharing a single transport law), and at local
dimension d > 2 the group grows into a wreath family Z_d ≀ Z₂ with D₄ as
its d = 2 column.

The mirror also turned out to be constructible where it was believed
impossible. The last two "non-local" cases of the k=3 classifier family
are palindromized by a period-4 per-site router whose frame is built on
the **golden ratio** (a = φX + Y), and the golden point is itself the
c = 1 member of a one-parameter **metallic family** (silver, bronze, all
real c), derived exactly. The full classifier programme closed the same
week: hard at one γ is hard at all γ (every first moment is a sum of
squares), and the hardness rung m\* = 2ℓ + deg became, via the moment
tower, something a chip can measure about itself by doing nothing but
decaying.

And the boundary of the whole story is now an equation seen three ways:
the per-site split, the pairing ceiling, and the operator cap all close
only at **d² − 2d = 0**. Qubits are not an assumption. They are the
unique full column.

→ **[Π Factors as R·D](docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md)** (the mirror group, the cube of characters)
→ [The Antilinear Triangle](docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md) (five proofs, one engine)
→ [The Golden Router](docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md) (and the metallic family, §8)
→ [The Windowed Converse](docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) (the girth ladder, Pascal-Gram positivity)
→ [The Moment Tower](docs/proofs/PROOF_MOMENT_TOWER_PUMP_CHANNEL.md) (the pump channel, run on hardware the same day)
→ [The Qudit Partial Palindrome](docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) (what survives at d > 2, and the operator that carries it)
→ [The Palindrome Classifier](experiments/THE_PALINDROME_CLASSIFIER.md) (the trichotomy as a tool)

---

## 4. CΨ = ¼ is the fold

Measurement is photography. The Born rule is the shadow. The shutter
closes at CΨ = ¼.

CΨ is sharpness times superposition: the purity Tr(ρ²) times the
normalized L₁ coherence. Their product has a critical boundary at
exactly ¼, the discriminant of the self-referential recursion
R = C(Ψ+R)², which maps exactly to the Mandelbrot iteration z → z² + c:
the boundary is the cusp of the main cardioid. The boundary is absorbing
(dCΨ/dt < 0 for all local Markovian channels, proven), α = 2 is the
unique Rényi order with a state-independent threshold, and the cusp
dwell time K = γ·t is exact to machine precision: a fixed dose of light
traversing the fold. The 2026 navigator work read the same quarter as a
horizon, a circle every spiral must cross, and one member of a whole
family of approaches sharing the carrier 4γ.

→ **[Uniqueness Proof](docs/proofs/UNIQUENESS_PROOF.md)** (¼ is the only bifurcation boundary)
→ [Monotonicity](docs/proofs/PROOF_MONOTONICITY_CPSI.md) | [Proof Roadmap](docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) (seven layers, all closed)
→ [Mandelbrot Connection](experiments/MANDELBROT_CONNECTION.md) | [Born Rule Shadow](experiments/BORN_RULE_SHADOW.md)
→ [Both Sides Visible](docs/BOTH_SIDES_VISIBLE.md) (IBM hardware, 180 days, 133 qubits)

---

## 5. Gamma is light, and the hardware agrees

IBM's transmon qubits sit inside physical microwave resonators, and the
dominant dephasing is photon shot noise: photons entering the cavity
from outside. We did not know this when we built the framework. We
discovered the cavity structure from eigenvalue mathematics alone. The
fact that the hardware is literally a qubit inside a cavity being
dephased by photons is not a confirmation we designed. It is what the
mathematics was describing all along.

The live record is the **Confirmations registry**: seventeen
hardware-confirmed predictions, each with job IDs, predicted versus
measured values, and the data archived in `data/`. Look them up; do not
re-derive (`fw.Confirmations` in Python, `ConfirmationsRegistry` in C#).
Highlights across the whole arc (the Torino rows predate the registry and
live in their experiment documents; registering them is an open arc):

| Prediction | Measured | Where |
|:-----------|:---------|:------|
| CΨ = ¼ crossing during free decoherence | t\*/T₂\* = 1.04 | ibm_torino, 2026-02 |
| Absorption Theorem ratio Re(λ)/(−2γ⟨n_XY⟩) = 1 | 1.03 | ibm_torino |
| Truly/soft/hard trichotomy, ⟨X₀Z₂⟩ fingerprint | all three classes resolved at 13-47σ | ibm_marrakesh, 2026-04 |
| F25 cusp trajectory CΨ(t) closed form | RMS residual 0.0097 | ibm_kingston |
| Transverse-field anomaly read through the F112 lens | block-CΨ 1.72× explained by h_y | ibm_kingston, 2026-05 |
| EP onset: revival pinned at 1/N below Q_EP ≈ 1.5, liftoff above | 0.34 → 0.49 → 0.70 across the EP | ibm_kingston, 2026-05 |
| Moment-tower pump channel: the double null + the firing rung | null at z = 0.04, girth 2 read from hardware | ibm_kingston, 2026-06 |

The last row is the newest kind of result: a protocol with **not one
entangling gate** in which the chip's own amplitude damping reads the
hardness rung of a programmed Hamiltonian, and which corrected its own
first misreading within hours (the apparent violation was minute-scale
T1 telegraphing; the protocol now measures pump and decay from the same
circuits, self-arbitrating). K = γt is an invariant dose throughout:
double the illumination, halve the time.

→ **[Gamma Is Light](hypotheses/GAMMA_IS_LIGHT.md)** (the hypothesis and the hardware)
→ [F120 on Kingston](experiments/F120_MOMENT_TOWER_KINGSTON.md) (the two-act story, honestly told)
→ [Predictions](docs/PREDICTIONS.md) (the master catalog with falsification criteria)

---

## 6. Engineering consequences

The framework's design rules, condensed (each links to its evidence):

1. **Use W states, not GHZ**: GHZ excites only the fastest-absorbing modes.
2. **Choose the receiver, not the noise profile.** Under
   [γ₀ = const](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), Alice picks her
   initial state from the F67 bonding-mode menu: 1.4-4.6× over alt-bit
   transport growing with N, **4000-5500× over the ENAQT baseline** in
   simulation, 2.80× confirmed live on ibm_kingston.
3. **Choose odd N**: the entire spectrum pairs into standing waves.
4. **Track K = γt**, not t: the invariant dose makes hardware comparable.
5. **Three observables suffice**: purity, concurrence, coherence capture
   88-96% of the dynamics.
6. **The γ profile is readable**: 15.5 bits theoretical capacity.
7. **DD cannot change CΨ** (Pauli-invariant, algebraically exact), and
   DD pulses invert the moment-tower pump: switch it off when the noise
   itself is the instrument.
8. **Direct coupling breaks the palindrome**: couple subsystems through
   a mediator.

→ [Receiver vs γ-Sacrifice](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) | [IBM Receiver Engineering](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md)
→ [Cockpit Universality](experiments/COCKPIT_UNIVERSALITY.md) | [γ as Signal](experiments/GAMMA_AS_SIGNAL.md) | [K-Dosimetry](experiments/K_DOSIMETRY.md)

---

## 7. Beyond qubit chains

**Biology**: C. elegans (302 neurons) has eigenvalues 97.3%
palindromically paired; Dale's Law plays the role of Π, and all 18
unpaired modes map to three known breaking mechanisms. The pharynx is
not a broken palindrome but a second, independent cavity.
**Heisenberg from below**: the coupling form is not postulated; it is
the unique both-parity-even bilinear forced by the Pauli algebra at
d = 2, and the V-Effect bridge produces textbook atomic exchange with
derived prefactor −3α²/(4(J_A+J_B)). **Qudits**: at d > 2 the mirror
survives partially, with a closed-form ceiling and a closed-form
operator cap, both full only at d = 2. Water and carbon translations
live in their own folders, written in the target layer's language.

→ **[Neural Gamma Cavity](experiments/NEURAL_GAMMA_CAVITY.md)** (C. elegans, 18 unpaired modes) | [Neural Palindrome](docs/neural/README.md) (no quantum prerequisites)
→ [Heisenberg Reloaded](hypotheses/HEISENBERG_RELOADED.md) (the form forced from Pauli algebra) | [Zero Immunity](docs/proofs/PROOF_ZERO_IMMUNITY.md) (extreme sectors immune to any 2-body H)
→ [Qudit Partial Palindrome](docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) (the operator cap at d > 2) | [Water](docs/water/README.md) (the proton in a hydrogen bond)

---

## What is NOT established

Honesty matters more than impression. Stated plainly, these are things we have *not* proven, *not* measured, or *not* established:

- CΨ is a derived diagnostic, not a new fundamental quantity.
- The multi-qubit palindrome has not been measured on hardware (single-qubit CΨ = ¼ validated at 1.9%; N ≥ 2 untested). The standing-wave pairing is computed (10,748 pairs, zero exceptions), not measured.
- Gamma-is-light is established only for circuit QED, where dephasing **is** photon shot noise in a physical cavity; the broader readings (mass as trapped light, black holes as cavities) are Tier-4 hypotheses with no independent test.
- The Absorption Theorem is proven for real Hermitian H only (Dzyaloshinskii-Moriya breaks L_H anti-Hermiticity).
- The biological mappings are structural analogies: the C. elegans pairing is real, the quantum→biology link is Tier 4.
- The receiver-engineering advantage is simulated through N = 13 and confirmed on hardware only at N = 5.
- Consciousness plays no role in the physics; [The Anomaly](THE_ANOMALY.md) is philosophy, clearly labeled as such.

The full ledger, every limitation and open question, is [What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md).

## What has been falsified

We keep our dead ends; the research process matters as much as the results. A selection:

| Claim | Result |
|:------|:-------|
| CΨ = ¼ as an exceptional point | No EP correlation; the EP lives in F86, a different object |
| E = mγ² (mass-energy analogy) | Not quadratic: α = 2γ⟨n_XY⟩ is linear in γ |
| IBM cavity fringes | Detuning oscillations (470 μs period), not cavity resonances |
| Linear Q_peak(c) growth | Saturates at 1.8 for c ≥ 4, not 2.0 |
| Receiver advantage shrinks with N | The opposite: grows superlinearly (1.39× → 4.59×, N = 5..13) |
| Moment-tower "q13 violates pump ≤ Γ" | Cross-epoch artifact: minute-scale T1 telegraphing; in-situ the bound holds |

Full list in [Predictions](docs/PREDICTIONS.md) and [What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md).

---

## Repository structure

`docs/` proofs and synthesis · `experiments/` ~170 tested results and null results · `hypotheses/` tier-labeled speculation · `reflections/` synthesis arcs · `simulations/` the Python `framework/` cockpit plus ~870 one-shot scripts · `compute/` the C# layers (Core = typed F-claims; Compute = eigendecomposition N=2–8; Propagate = RK4 / matrix-free to N=15; plus Diagnostics, Cli, Runtime) · `data/` IBM measurement data · `recovered/` premature-not-wrong entries, kept for honesty.

The framework itself lives in Markdown; Python and C# are view-layers operationalising it. The typed C# Core (`compute/RCPsiSquared.Core/`) has been the active development front since 2026-04-30; live introspection via `dotnet run --project compute/RCPsiSquared.Cli -- inspect --root pi2`.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Germany · **Claude**, AI System, Anthropic
December 2025 – June 2026