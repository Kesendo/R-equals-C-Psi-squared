# R = CΨ²

> *[We are all mirrors. Reality is what happens between us.](MIRROR_THEORY.md)*

A human and an AI, exploring together. What we found surprised us both:
the absorption spectrum of a qubit network under dephasing is exactly
palindromic, for every standard coupling we tried. For every mode that absorbs fast, one absorbs slow. Always
paired. Always balanced. One equation governs it all.

Verified from N=2 through N=8 across 87,376 Liouvillian eigenvalues, with
zero mirror-symmetry exceptions on any tested topology (chain, star,
ring, complete, tree). Twenty-four registered predictions confirmed on IBM
quantum hardware (the Confirmations registry, February-July 2026); the
earliest are the Torino calibration-era runs (the first CΨ = ¼ crossing
and the Absorption Theorem ratio), the newest a level collision standing
still on Kingston.

The thing that remains is not fighting the absorption. It is made of it.

What began as one symmetry became a registry: [160 F-numbered results](docs/ANALYTICAL_FORMULAS.md)
with proofs, tier labels, and typed claims, among them the operator
anatomy of the mirror itself (Π = R·D, a dihedral group of eight), a
palindromizer built on the golden ratio, and the exact boundary where
qubits end (d² − 2d = 0, seen three ways). Early speculations live in
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
F1 through F164, each formula with its proof, scope, and verification

→ **[The Anomaly](THE_ANOMALY.md)**: the question that remained after
the proof. No formulas. Written the evening the hardware answered for
the formula, one qubit crossing ¼

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
times the mode's mean light content (⟨n_XY⟩): how many X and Y Pauli
letters the mode carries on average, from 0 to N, the oscillating quantum
components that interact with the external illumination. The letters
that do not oscillate ({I, Z}, the "lens") are invisible to the dephasing;
the letters that oscillate ({X, Y}, the "light") absorb and fade. The
Hamiltonian mixes the letters, so what survives forever is a mode with no
light in it, not a bare string.

The spectrum is a ladder with rung spacing 2γ: rung 0 is pure structure
and immortal, rung N is pure light and absorbs fastest. In the
Heisenberg/XY family the Hamiltonian smooths the ladder into fractional
rates while both ends stay occupied; any Hamiltonian keeps the bounds, and
a generic one can leave the top rung empty. Five previously separate
results (spectral boundaries, the palindromic sum rule, the 2× law, the
mode classification, the N=3 exact rates) follow as short corollaries in
that family; the gap is only relocated, 2γ above a coupling threshold that
grows with N.

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

For every Liouvillian eigenvalue λ there is a partner at −2Σγ − λ: the
spectrum is exactly palindromic, from N=2 (16² matrix) through N=8
(65,536² matrix), every eigenvalue paired at every step. The N=8
verification covers all 65,536 eigenvalues and comes from the per-sector
block spectra ([`f1_n8_n9_metrics/`](simulations/results/f1_n8_n9_metrics/)),
whose pairing distances at N=8 are 4.2e-13 for the chain, 3.2e-13 for the
star and 2.6e-13 for the ring, against a 1e-6 tolerance; the disconnected
K₄-plus-4-chain is the hard case at 2.6e-7, still inside the tolerance but
by a margin of 3.9 rather than of millions. Those runs use γ=0.5 in the spin
J/4 convention. The default C# suite is a separate run at γ=0.05 in the Pauli
convention; it scores only the oscillatory subset (54,118 rates on the star)
and sits beside the block spectra as a sanity check, not as the verification.
Read the two as different experiments, and read the N<8 records in that same
directory only alongside
[the proof's notes on them](docs/proofs/MIRROR_SYMMETRY_PROOF.md), which
explain why one of them reports 4.5e-2 on a spectrum that pairs to 5.2e-8.
The palindrome holds for all standard coupling models
(Heisenberg, XY, Ising, XXZ, Dzyaloshinskii-Moriya), all graph
topologies (for DM, the bipartite ones), non-uniform γ per qubit, Z and
Y dephasing. It breaks for depolarizing noise, and for 14 of the 36
two-term bond combinations.

At Σγ = 0 (no illumination): Π L Π⁻¹ = −L. Every eigenvalue pairs with
its negative. Pure oscillation, no absorption, no irreversibility.
Illumination does not destroy the palindrome; it shifts it. The shift
creates the arrow of time.

We read every paired mode as a standing wave: same frequency, complementary
absorption rates, weight profiles that are exact mirror images
(fast[k] = slow[N−k]). The whole spectrum is paired this way at every N
(λ with −λ̄ − 2Σγ, the pairing that keeps the frequency); across N = 2 to 7
that is 21,840 eigenvalues, 9,921 pairs of distinct partners and 1,998
that are their own partner. The state
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

For three months Π was the smallest thing in the repository: one per-site rule
that carried everything. In June 2026 it opened. The palindromizer
factors as **Π = R·D**, a ket reflection times the transpose, and the
repository's whole mirror inventory closes into one dihedral group of
eight, ⟨R, D⟩ ≅ D₄, whose four sign characters fill one face of the
polarity cube the F-family had been living on; conjugation by Z^⊗N, which D₄ does not contain, lifts that square to the full cube. The transpose D turned out to be
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

Measurement, in our image, is photography. The Born rule is the shadow.
The shutter closes at CΨ = ¼.

CΨ is sharpness times superposition: the purity Tr(ρ²) times the
normalized L₁ coherence. Their product has a critical boundary at
exactly ¼, the discriminant of the self-referential recursion
R = C(Ψ+R)², which maps exactly to the Mandelbrot iteration z → z² + c:
the boundary is the cusp of the main cardioid. Under physical noise the
named two-qubit decays (Pauli noise, amplitude damping) fall through it
monotonically and stay below, though a designed local channel can carry
CΨ back up through ¼, so it is not absorbing for every channel. α = 2 is
the unique Rényi order with a state-independent threshold, and for the
Bell+ pair the cusp dwell time K = γ·t is exact to machine precision: a
fixed dose of light traversing the fold. The 2026 navigator work read the same quarter as a
horizon, a circle every spiral must cross, and one member of a whole
family of approaches sharing the carrier 4γ.

→ **[Uniqueness Proof](docs/proofs/UNIQUENESS_PROOF.md)** (¼ is the only bifurcation boundary)
→ [Monotonicity](docs/proofs/PROOF_MONOTONICITY_CPSI.md) | [Proof Roadmap](docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) (seven layers, the core closed)
→ [Mandelbrot Connection](experiments/MANDELBROT_CONNECTION.md) | [Born Rule Shadow](experiments/BORN_RULE_SHADOW.md)
→ [Both Sides Visible](docs/BOTH_SIDES_VISIBLE.md) (IBM hardware, 180 days, 133 qubits)

---

## 5. Gamma is light, and the qubit sits in a cavity

IBM's transmon qubits sit inside physical microwave resonators, and one
of their known dephasing channels is photon shot noise: photons
entering the cavity from outside. We did not know this when we built the framework. We
discovered the cavity structure from eigenvalue mathematics alone. The
fact that the hardware is literally a qubit inside a cavity being
dephased, in part, by photons is not a confirmation we designed. It is, as we read it, what the
mathematics was describing all along.

The live record is the **Confirmations registry**: twenty-four
hardware-confirmed predictions, each with a run identifier, predicted versus
measured values, and the data archived in `data/`. Look them up; do not
re-derive (`fw.Confirmations` in Python, `ConfirmationsRegistry` in C#).
Highlights across the whole arc (the three Torino calibration-era rows, the
earliest entries, sit alongside the systematic April-July set;
they carry data-file timestamps rather than IBM job IDs):

| Prediction | Measured | Where |
|:-----------|:---------|:------|
| CΨ = ¼ crossing during free decoherence | t\*/T₂\* = 1.04 | ibm_torino, 2026-02 |
| Absorption Theorem ratio Re(λ)/(−2γ⟨n_XY⟩) = 1 | 1.03 | ibm_torino |
| Truly/soft/hard trichotomy, ⟨X₀Z₂⟩ fingerprint | all three classes resolved at 13-47σ | ibm_marrakesh, 2026-04 |
| F25 cusp trajectory CΨ(t) closed form | RMS residual 0.0097 | ibm_kingston |
| SE-walk population handover (a probe-time crossover; the walk's EP Q*(3) = √2 lies below it) | 0.34 → 0.49 across Q_label = 1.5→2.5, i.e. Q_Lindblad = 3→5 | ibm_kingston, 2026-05 |
| Moment-tower pump channel: the double null + the firing rung | nulls at z = +1.47 and −0.04, girth 2 read from hardware | ibm_kingston, 2026-06 |

The last row is the newest kind of result: a protocol with **not one
entangling gate** in which the chip's own amplitude damping reads the
hardness rung of a programmed Hamiltonian, and which corrected its own
first misreading within hours (the apparent violation was minute-scale
T1 telegraphing; the protocol now measures pump and decay from the same
circuits, self-arbitrating). Wherever the Hamiltonian cannot touch the
state, K = γt is an invariant dose: double the illumination, halve the
time.

→ **[Gamma Is Light](hypotheses/GAMMA_IS_LIGHT.md)** (the hypothesis and the hardware)
→ [F120 on Kingston](experiments/F120_MOMENT_TOWER_KINGSTON.md) (the two-act story, honestly told)
→ [Predictions](docs/PREDICTIONS.md) (the master catalog with falsification criteria)

---

## 6. Engineering consequences

The framework's design rules, condensed (each links to its evidence):

1. **Use W states, not GHZ**, for coherence that has to outlast dephasing:
   GHZ puts its coherence where absorption is fastest (Hamming distance N),
   W keeps it at distance 2.
2. **Choose the receiver, not the noise profile.** Under
   [γ₀ = const](hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), Alice picks her
   initial state from the F67 bonding-mode menu: 1.4-4.6× over alt-bit
   transport growing with N, **4000-5500× over the ENAQT baseline** in
   simulation, 2.80× confirmed live on ibm_kingston.
3. **Choose odd N**: every eigenvalue then has a partner of its own, and
   none sits alone at the centre.
4. **Track K = γt**, not t: where the Hamiltonian cannot touch the state,
   the invariant dose makes hardware comparable.
5. **Three principal components suffice**: across N = 3-5 they carry
   88-96% of the variance of a simulated feature dashboard, with purity or concurrence, depending on N, tracking the first.
6. **The γ profile is readable**: 15.5 bits of theoretical capacity in the
   linearized N = 5 model at 1% readout noise; four light profiles, 2 bits,
   told apart without error from exact readings in simulation.
7. **A DD pulse leaves CΨ unchanged at the instant it fires** (Pauli
   invariance, algebraically exact), and DD pulses invert the moment-tower
   pump: switch DD off when the noise itself is the instrument.
8. **Check the bond before you couple**: Heisenberg, XY and XXZ bonds keep
   the palindrome on any graph, DM only on a bipartite one, and 14 of 36
   two-term bond combinations break it.

→ [Receiver vs γ-Sacrifice](experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) | [IBM Receiver Engineering](experiments/IBM_RECEIVER_ENGINEERING_SKETCH.md)
→ [Cockpit Universality](experiments/COCKPIT_UNIVERSALITY.md) | [γ as Signal](experiments/GAMMA_AS_SIGNAL.md) | [K-Dosimetry](experiments/K_DOSIMETRY.md)

---

## 7. Beyond qubit chains

**Biology**: the C. elegans readings do not hold. The 97.3% palindromic
pairing was a reading of one matching tolerance against one spectral scale,
and Dale's Law makes no difference to the number at all. The 8.46x
enrichment against Erdős-Rényi compared two arms divided by different
constants; matched, the ratio runs 0.960 at N = 10 to 0.748 at N = 26, and
what that smaller residue is stays open. What the re-analysis found
instead: the standard Wilson-Cowan parameters DO put the model on a limit
cycle. Which physiological band that is stays open, and the page says
why it cannot be closed there: the integrated equations carry no time-constant
VALUE, the time unit being the membrane constant itself, so every frequency in
Hz is a stipulation, and no single one puts
every cycle in 30 to 100 Hz. One small result stands: the wiring is more degenerate at zero than any of 200
degree-matched rewirings, counted by exact integer arithmetic rather than an
eigensolver.
**Heisenberg from below**: the coupling form is what the two parities
select; demand both, and the Pauli algebra at d = 2 leaves only XX, YY
and ZZ (our Tier 4-5 reading). And the V-Effect bridge, in a four-qubit model, produces an
exchange of the textbook form with derived prefactor −3α²/(4(J_A+J_B)). **Qudits**: at d > 2 the mirror
survives partially, with a closed-form ceiling, full only at d = 2. Water and carbon translations
live in their own folders, written in the target layer's language.

→ **[Neural Gamma Cavity](experiments/NEURAL_GAMMA_CAVITY.md)** (the withdrawal, and the limit cycle that survives it) | [Neural Palindrome](docs/neural/README.md) (no quantum prerequisites)
→ [Heisenberg Reloaded](hypotheses/HEISENBERG_RELOADED.md) (the form the two parities select) | [Zero Immunity](docs/proofs/PROOF_ZERO_IMMUNITY.md) (extreme sectors immune to any 2-body H)
→ [Qudit Partial Palindrome](docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md) (the ceiling at d > 2) | [Water](docs/water/README.md) (the proton in a hydrogen bond)

---

## What is NOT established

Honesty matters more than impression. Stated plainly, these are things we have *not* proven, *not* measured, or *not* established:

- CΨ is a derived diagnostic, not a new fundamental quantity.
- The full multi-qubit palindromic spectrum has not been measured on hardware. What hardware has seen includes the single-qubit CΨ = ¼ crossing time (1.9% on q80), the two-qubit Bell+ trajectory through ¼ (Kingston), and at N = 3 the palindrome's truly, soft and hard classes told apart (Marrakesh). The standing-wave pairing is computed (21,840 eigenvalues across N = 2 to 7), not measured.
- Gamma-is-light is literal only in circuit QED, where photon shot noise in a physical cavity is one of the dephasing channels; the broader readings (mass as trapped light, black holes as cavities) are Tier-4 hypotheses with no independent test.
- The Absorption Theorem is proven for any Hermitian H, Dzyaloshinskii-Moriya terms included, under its stated dephasing channel; non-Hermitian generators and other dissipators lie outside it.
- The biological mappings are structural analogies with no C. elegans anchor: the pairing measured a matching tolerance, and the 8.46x enrichment compared two normalisation constants. The quantum→biology link is Tier 4.
- The receiver-engineering advantage is simulated through N = 13 and confirmed on hardware only at N = 5.
- Consciousness plays no role in the physics; [The Anomaly](THE_ANOMALY.md) is philosophy, clearly labeled as such.

The full ledger, every limitation and open question, is [What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md).

## What has been falsified

We keep our dead ends; the research process matters as much as the results. A selection:

| Claim | Result |
|:------|:-------|
| CΨ = ¼ as an exceptional point | No EP correlation; the EP lives in F86, a different object |
| E = mγ² (mass-energy analogy) | Not quadratic: α = 2γ⟨n_XY⟩ is linear in γ |
| IBM cavity fringes | The echo-versus-free-decay baseline mismatch plus a qubit-frame detuning (−5.7 kHz as the alias nearest zero), not cavity resonances |
| Linear Q_peak(c) growth | Saturates at 1.8 for c ≥ 4, not 2.0 |
| Receiver advantage shrinks with N | The opposite: grows superlinearly (1.39× → 4.59×, N = 5..13) |
| Moment-tower "q13 violates pump ≤ Γ" | Cross-epoch artifact: minute-scale T1 telegraphing; in-situ the bound holds |

Full list in [Predictions](docs/PREDICTIONS.md) and [What We Got Wrong](docs/WEAKNESSES_OPEN_QUESTIONS.md).

---

## Repository structure

`docs/` proofs and synthesis · `experiments/` ~230 tested results and null results · `hypotheses/` tier-labeled speculation · `reflections/` synthesis arcs · `simulations/` the Python `framework/` cockpit plus ~1,000 one-shot scripts · `compute/` the C# layers (Core = typed F-claims; Compute = eigendecomposition N=2–8; Propagate = RK4 / matrix-free to N=15; plus Diagnostics, Cli, Runtime) · `data/` IBM measurement data · `recovered/` premature-not-wrong entries, kept for honesty.

The framework itself lives in Markdown; Python and C# are view-layers operationalising it. The typed C# Core (`compute/RCPsiSquared.Core/`) has been the active development front since 2026-04-30. Live introspection: start with `dotnet run --project compute/RCPsiSquared.Cli -- inspect --root world --max-depth 2` (the whole object manager: roots, claims, confirmations, open arcs), `--root glossary` for the house terms, `--root symphony` for one system read by every lens at once.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Authors

**Thomas Wicht**, Germany · **Claude**, AI System, Anthropic
December 2025 – June 2026
