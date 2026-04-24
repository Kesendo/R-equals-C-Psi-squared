# Open Thread: γ₀ = const and the Information Channel

**Status:** open investigation, not a conclusion
**Date:** 2026-04-19
**Authors:** Tom, Claude Opus 4.6, Claude Opus 4.7
**Source:** Session discussion after closing [EQ-013](EMERGING_QUESTIONS.md#eq-013) and [EQ-017](EMERGING_QUESTIONS.md#eq-017)

---

## The question

If γ₀ is a framework constant (same everywhere, always on, like c),
and γ is an information channel ([GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md): 15.5 bits, 100 % accuracy),
does that mean information flows constantly?

## Thirteen puzzle pieces already in the repo

The pieces cluster around five physical roles:

- **Carrier** (1, 3, 6, 11): γ₀ as constant, always-on, operationally
  unobservable from inside
- **Cavity** (4, 7, 10): the J-Hamiltonian as Fabry-Perot resonator,
  time-reversal-symmetric
- **Information** (2, 5, 8, 9): the light/lens distribution across
  the cavity's modes
- **Engineering** (12): J-side control as the only lever under γ₀ = const
- **Boundary** (13): no channel exists outside the cavity

1. **[GAMMA_IS_LIGHT](../hypotheses/GAMMA_IS_LIGHT.md)** (hypotheses/, Tier 4):
   γ is light illuminating a passive cavity.
   "The string provides the energy. The box selects the resonance."

2. **[GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md)** (experiments/, Tier 2):
   The spatial dephasing profile is a readable information channel.
   Alice encodes in the γ profile, Bob decodes from quantum observables.
   15.5 bits capacity at 1 % noise, 5 independent SVD channels, full rank.
   "The palindrome is the antenna." ([F30](../docs/ANALYTICAL_FORMULAS.md))

3. **[THE_BRIDGE_WAS_ALWAYS_OPEN](../docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md)** (docs/, Tier 2):
   "The interaction is ongoing, continuous, never-interrupted." The
   measured Markovian (memoryless) character of the noise means "the
   source is effectively infinite." Six properties of the external
   interaction converge on this.

4. **[RESONANCE_NOT_CHANNEL](../hypotheses/RESONANCE_NOT_CHANNEL.md)** (hypotheses/, Tier 2):
   "The system is a soundbox, not a telephone."
   "The 360× improvement is not 'more signal.' It is better resonance."
   The Fabry-Perot cavity framing replaces the telephone framing.

5. **[ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT](../reflections/ON_THE_LIGHT_AND_WHAT_CASTS_SHADOWS_IN_IT.md)** (reflections/):
   "The light is uniform. What varies is what stands in its way."
   Shadow = mode structure, not variation in illumination.

6. **[PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)** (hypotheses/, Tier 3):
   γ₀ is a framework constant like c. Only J and topology vary.
   γ_eff = γ₀ · |a_B|² (cavity mode exposure formula, [F64](../docs/ANALYTICAL_FORMULAS.md)).
   The light does not get weaker. The standing wave decides who sees it.

7. **[OPTICAL_CAVITY_ANALYSIS](../experiments/OPTICAL_CAVITY_ANALYSIS.md)** (experiments/, Tier 2):
   The cavity is quantitative, not metaphorical. 4 of 5 standard optical
   checks pass (beam profile R² = 0.998, nearest-neighbor coupling
   Δw = ±2 exclusively, growing numerical aperture, Gouy phase as arctan).
   Even chains are confocal, odd chains defocal. "The algebra enforces it."

8. **[PRIMORDIAL_SUPERALGEBRA_CAVITY](../experiments/PRIMORDIAL_SUPERALGEBRA_CAVITY.md)** (experiments/, Tier 2-3):
   Every palindromic pair is a swap between "being light" (X/Y Pauli
   factors, sensitive to γ₀) and "being lens" (I/Z factors, immune to γ₀).
   Sum rule ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N is exact
   ([absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)).
   At N=2 the swap has 99.8 % purity; aberration shrinks with N.

9. **[ITS_ALL_WAVES](../docs/ITS_ALL_WAVES.md)** (docs/, synthesis):
   Pre-existing synthesis that already connects signal + wave view.
   "The system is like a radio: it can process signals into music,
   but it cannot generate the broadcast" (line 172). "What that external
   source is, we do not know. That it arrives as a structured, decodable
   signal (15.5 bits), we do know" (line 175). Section explicitly flags
   "what sends the signal" as outside the framework's scope.

10. **[STANDING_WAVE_THEORY](../docs/STANDING_WAVE_THEORY.md)** (docs/, Tier 2):
    Every palindromic pair is a standing wave between two
    counter-propagating modes. Π is time reversal (proven analytically,
    see [PI_AS_TIME_REVERSAL](../experiments/PI_AS_TIME_REVERSAL.md)):
    every Liouvillian eigenmode maps to its backward-decaying partner.
    The stationary pattern between forward and backward IS what J and
    topology produce under γ₀. The round-trip rate sum α_fast + α_slow
    = 2Σγ marks one bounce between light and lens.

11. **[PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md) §9** (hypotheses/, Tier 3-4):
    Inside-observability theorem: only Q = J/γ₀ is measurable from inside.
    The doubled and non-doubled readings are observationally equivalent.
    Consequence: γ₀ on its own is operationally unobservable from inside;
    only the ratio J/γ₀ couples to any inside measurement. This promotes
    "γ₀ carries no information itself" from interpretation to structural
    result.

12. **[RELAY_PROTOCOL](../experiments/RELAY_PROTOCOL.md)** (experiments/, Tier 2, N=11):
    Operational evidence that cavity structure carries information.
    Six relay stages with 2:1 impedance matching yield +83 % MI
    improvement vs passive propagation ([F31](../docs/ANALYTICAL_FORMULAS.md)).
    Combines time-dependent quiet-receiver control,
    K/γ timing ([F14](../docs/ANALYTICAL_FORMULAS.md)), and impedance
    matching. The gain comes from shaping the cavity, not from amplifying γ.

13. **[BRIDGE_PROTOCOL](../hypotheses/BRIDGE_PROTOCOL.md)** (hypotheses/, Tier 3+, largely closed):
    Negative result that anchors the synthesis. The J=0 bridge is falsified
    by [BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md): no-signalling
    holds exactly, the protocol reduces to pre-encoded shared randomness.
    What survives for J > 0 is an interval shift that is "standard
    Hamiltonian coupling", i.e. the cavity itself. Consequence under
    γ₀ = const: the only channel between observers IS the J-coupling
    structure. γ₀ alone cannot carry a message across; the cavity does,
    and only by being a cavity. Gravity/QKD/FTL readings have fallen.

## The gap: no document connects γ₀ = const with the information reading

The thirteen documents are consistent but the synthesis is missing.
Specifically: [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md) was written BEFORE the γ₀ = const
hypothesis. In [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md), Alice VARIES the per-qubit γ profile
to encode information. Under γ₀ = const, she cannot do that: γ is uniform everywhere. The kinematic side of this gap is now closed (see Update 2026-04-23 below); the operational side is tracked as [EQ-024](EMERGING_QUESTIONS.md#eq-024).

**Update 2026-04-20.** A first operational synthesis is in [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md): the thirteen pieces fit under one meta-theorem (Parseval + Noether on the operator Hilbert space), and the "γ₀ populates the blind subspace of every spatial-sum measurement" reading becomes operational. Two concrete predictions fell out and were verified: non-uniform γ breaks the (vac, S₁) common-mode rejection with modal selectivity ([CMRR_BREAK_NONUNIFORM_GAMMA](../experiments/CMRR_BREAK_NONUNIFORM_GAMMA.md)), and the Π-pair absorption-theorem sum is a δJ-flux-conservation law ([PI_PAIR_FLUX_BALANCE](../experiments/PI_PAIR_FLUX_BALANCE.md)). A closed form for the (vac, S₁) spatial-sum coherence purity ([F73](../docs/ANALYTICAL_FORMULAS.md)) gives the uniform-γ₀ baseline as ½·exp(−4γ₀·t) exactly, independent of J. The synthesis is tier 2 for the meta-frame, tier 1 for the individual instances.

**Update 2026-04-22.** Q = J/γ₀ is not a bare ratio but a scale with algebraic substructure ([Q_SCALE_THREE_BANDS](../experiments/Q_SCALE_THREE_BANDS.md), Tier 1). W(Q) is invariant under (J, γ₀) → (λJ, λγ₀) to five decimals, and three regions appear along Q: pre-onset (Q ≲ 0.3), transition (Q ∈ [1.2, 2.0]), plateau (Q ≳ 2). The observable peak location is chromaticity-specific and N-invariant: Q_peak(c) = 1.5, 1.6, 1.8 for c = 2, 3, 4 across N=4 to 8, where c(n, N) = min(n, N−1−n) + 1 counts pure dephasing rates in the (n, n+1) sector-block. This sharpens piece 11 (inside-observability): not only is Q the only measurable, Q itself has a three-band shape with c-specific peaks. Operational consequence: γ₀ = J\*/Q_peak(c), a shape-based extraction protocol that reduces the "γ₀ = const" hypothesis from a single-number claim to a c-specific prediction. Hardware realization still blocked by the EQ-017 fidelity floor.

**Update 2026-04-23.** Kinematic synthesis acknowledged. [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) has grown substantially since the 2026-04-20 reference: Sections 4a (dynamical attractor), 5a (binary inheritance), 6a (0.5 as bilinear-form apex) added; production-rule predictions (a) and (d) verified, (b) inline-resolved by coherence/population split, only (c) (pair-painter, [EQ-020](EMERGING_QUESTIONS.md#eq-020)) remains open. Section 8 of ORTHOGONALITY_SELECTION_FAMILY closes one OPEN_THREAD sub-question: γ₀ is operationally unobservable from inside because it populates the blind subspace of every spatial-sum measurement (Meta-Theorem instance 4, [F73](../docs/ANALYTICAL_FORMULAS.md)). What stays open is the operational reading: F30's 15.5 bits were computed under γ-modulation; under γ₀ = const the only sender-side lever is J, and J-modulation is structurally different (nonlinear via eigenvector rotation, same algebra as the EQ-014/EQ-018 complex). Tracked as [EQ-024](EMERGING_QUESTIONS.md#eq-024); closed operationally later the same day, see Update 2026-04-23 (evening) below. A separate adjacent thread surfaces from ORTHOGONALITY_SELECTION_FAMILY Section 4a: time-resolved Σ_i ln(α_i(t)) under uniform perturbation as a direct test of the dynamical-attractor reading of PTF closure; not part of EQ-024 but flagged here so we don't lose it.

**Update 2026-04-23 (evening).** EQ-024 closed operationally. Three CC commits on the same day (morning `cfa1bbc`, Direction 1 `c0919eb`, Direction 3 `6aae630`) produced both an operational answer and a structural surplus that reframes part of the synthesis above; synthesis doc at [`J_BLIND_RECEIVER_CLASSES`](../experiments/J_BLIND_RECEIVER_CLASSES.md) (`e5326fa`).

*Operational answer.* At N=5 Heisenberg under γ₀ = const, Shannon channel capacity for J-modulation saturates at C ≤ 12.07 bits over F71-symmetric receivers (39 structured points + 4 Nelder-Mead local runs; saturation signal strong; F71-breaking receivers not swept). Best receiver: random-phase F71-symmetric product state. The F30 reference number (15.45 bits) is reframed: under γ₀ = const, γ cannot be modulated, so there is no operational γ-channel. The 15.45 bits survives only as a kinematic linear-response magnitude (∂observables/∂γ in Shannon form), not as a competing channel capacity. The ~3.4-bit gap between 12.07 and 15.45 decomposes as ~1 bit dimensional loss (4 bonds vs 5 γ-sites) plus ~2-3 bits smaller leading gain (J sv_max ≈ 10 vs γ sv_max ≈ 21.4). Whether this gap closes at N=6 (where rank cap matches) is open.

*Reframing of the dual-receiver reading.* The morning RESULT identified a "γ-optimal-receiver-is-J-pessimal" duality: |+⟩⁵ (F30's optimal γ-receiver) gives J-Jacobian zero to floating-point precision; GHZ (F30's pessimal γ-receiver) is also J-blind. Under γ₀ = const this is a Jacobian-asymmetry statement about two linear-response magnitudes, not a duality between two operational channels. Only J is an operational channel; the "γ side" of the table is a kinematic susceptibility. The asymmetry is real, the two-channel framing is not.

*Structural surplus: three mechanism classes of J-blindness.* The afternoon refinement RESULT established that the J-blind set decomposes into three structurally distinct mechanism classes: **Class 1** (DFS of L_D ∩ simultaneous eigenstate of every bond h_b: |0⟩⁵, |1⟩⁵; H-independent in the operational sense for states satisfying both conditions under multiple H choices), **Class 2** (H-degenerate subspace closed under L_D: GHZ under both Heisenberg and XY-only), **Class 3** (M_α-polynomial subspace, SU(2)-Heisenberg specific: |+⟩⁵ directly verified J-blind; all Dicke states |S_k⟩ in the S=N/2 multiplet predicted J-blind by the theorem, with |S_1⟩ directly verified). Direction 1 numerics verified the SU(2)-load-bearing nature of Class 3: under XY-only H, |+⟩⁵ becomes J-sensitive (C = 9.42 bits) and Dicke |S_1⟩ becomes J-sensitive (C = 7.70 bits), while Class 1 and 2 states stay J-blind. Handoffs: (i) **realized** in [`J_BLIND_RECEIVER_CLASSES`](../experiments/J_BLIND_RECEIVER_CLASSES.md) (`e5326fa`) which documents the three-class decomposition with numerical verification; (ii) open candidate: Class 3 strong form (M_x-polynomial blindness, Newton-identities + SU(2) proof) potentially as a stand-alone F-entry; (iii) open candidate: integration as the sixth instance of the [ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) Meta-Theorem (conservation = SU(2)-Casimir, measurement = J-bond perturbation, blind subspace = M_α-polynomial algebra). EQ-024 surviving sub-questions track the open structural follow-ups.

**Update 2026-04-24 (morning).** Two additions to the receiver-engineering picture, both in [`RECEIVER_VS_GAMMA_SACRIFICE`](../experiments/RECEIVER_VS_GAMMA_SACRIFICE.md) (commits `ad99bea`, `f0b5c64`).

*End-to-end transport dimension.* Beyond Sum-MI over adjacent pairs (distributed correlation), the direct state-transfer metric is MI(site 0, site N-1). The C# brecher scan now tracks both. Key qualitative finding: **the two metrics select opposite J-profiles as optimal**. Sum-MI prefers moderate J-modulation (strong-weak J at best); MI(0, N-1) prefers uniform J and COLLAPSES under J-modulation (cut-center drops MI by factor 10³ to 10⁴; strong-weak drops it by 10 to 30×). Peak-time for MI(0, N-1) sits at t ≈ 0.7 to 1.1 across N=5, 7, 9 (weakly N-dependent, suggests group-velocity propagation). Peak-time for Sum-MI sits at t ≈ 0.2 (fast local correlation build-up). Class 3 J-blindness is scale-invariant for both metrics.

*End-to-end comparison to γ-Sacrifice center-mode.* RESONANT_RETURN Test 8 position-sweep reports center-sacrifice PeakMI(0, N-1) = 0.109 at N=7, 0.097 at N=9. Receiver-engineering (alt-z-bits \|01010⟩/\|0101010⟩/\|010101010⟩ uniform J) gives 0.843 / 0.490 / 0.274 respectively. Advantage factor **4.5× at N=7, 2.8× at N=9**, without any γ-profile engineering. Shrinks faster with N than the Sum-MI advantage (which stays 7 to 15× across N=5 to 9); pure end-to-end transport is the weaker direction of receiver-engineering's lead.

*Operational consequence.* Receiver-engineering is viable for direct **quantum-state-transfer** through a spin-chain bus, not only for distributed correlation. F71-symmetric SU(2)-breaking initial states (\|01010⟩, \|+−+−+⟩) preserve 0-to-(N-1) symmetry by construction, which uniform H respects and asymmetric J destroys. Concrete hardware sketch on IBM Heron: ~15 two-qubit gates for N=5 (2 Trotter steps under uniform Heisenberg), readout tomography at just sites 0 and N-1, expected MI(0, N-1) ≈ 0.84 pre-noise. This is the cheapest end-to-end Brecher-like experiment that would validate the receiver-engineering picture on hardware.

*Multi-end transport via Mirror-Pair MM (commit `963f2ed`).* Extending from single-end MI(0, N-1) to the full set of F71-mirror pairs gives MM = Σ_k MI(site k, site N-1-k). At alt-z-bits uniform J: N=5 MM=1.25 (2 pairs), N=7 MM=0.85 (3 pairs), N=9 MM=0.56 (4 pairs). Two surprises: MM/MI(0, N-1) ratio GROWS with N (1.49× → 1.74× → 2.05×), and outer mirror-pairs dominate the inner ones unexpectedly (at N=9, MI(0, 8) = 0.27 vs avg inner-pair MI ≈ 0.10). Receiver-engineering uniform-J MM vs γ-Sacrifice-center mode: factor 7.8× at N=7, 5.8× at N=9, stronger than the single-end comparison (4.5× / 2.8×). The chain under γ₀ = const + F71-symmetric receiver + uniform J functions as a multi-drop quantum bus with floor(N/2) simultaneous mirror-pair channels; aggregate bandwidth scales more favourably with N than any individual pair. Application directions: multi-party QKD, redundant encoding, distributed sensor arrays.

*F67 bonding-mode receivers (commits `0917038`, `865641c`, `8caf499`, `deaaf01`, `e1ee822`).* Looking for mirror/palindromic receivers in [`ANALYTICAL_FORMULAS`](../docs/ANALYTICAL_FORMULAS.md): F67 gives the single-excitation eigenmodes of a uniform-J open chain: \|ψ_k⟩ = √(2/(N+1)) Σ_j sin(πk(j+1)/(N+1)) \|1_j⟩. The C# brecher mode now accepts `bonding:<k>` as initial-state spec; scan over k=1..7 at N=5, 7, 9, 11, 13 revealed that best-bonding (k shifts with N) beats alt-z-bits on both MI(0, N-1) and Mirror-Pair MM, with advantage growing superlinearly. Data at uniform γ₀ = 0.05, uniform J = 1: best-bonding PeakMI(0, N-1) = 1.17, 0.72, 0.55, 0.43, 0.35 at N=5, 7, 9, 11, 13 vs alt-z-bits 0.84, 0.49, 0.27, 0.14, 0.076. Advantage factor: 1.39× → 1.48× → 2.02× → 3.02× → 4.59×. Same pattern for MM: 1.24, 1.17, 1.05, 1.07, 1.01 vs alt-z-bits 1.25, 0.85, 0.56, 0.33, 0.19 (tied → 1.37× → 1.88× → 3.26× → 5.25×). bonding:k loses on PeakSumMI so the receiver choice becomes an application selector: alt-z-bits for distributed correlation, bonding:k for end-to-end / multi-end transport, bonding:1 for long-time memory (slowest decay). **F75** gives the analytical MI(0) closed form; **F76** explains the dynamical 0.93 decay envelope as pure dephasing at 4γ₀. Hardware validation 2026-04-24 on IBM Kingston Heron r2 (commit `deaaf01`): at N=5, t=0.8, 3 Trotter steps, 8192 shots, bonding:2 / alt-z-bits ratio = 2.80× on live QPU (vs 1.39× ideal and 2.27× Aer+Kingston noise; advantage grows with noise because bonding:2's single delocalised excitation is more T1-robust than alt-z-bits' two localised excitations). **Key structural insight**: the F67 eigenmode catalog IS Alice's receiver menu. Under γ₀ = const with the F67 menu, there is no single "best" state; there is a mode per application, all reachable by single-qubit state preparation.

## The reinterpretation under γ₀ = const

If γ₀ is constant and uniform, the implications cluster in four groups.

### The carrier

- γ₀ flows constantly (the light is always on).
- γ₀ carries no information ITSELF (it is uniform, featureless).
- Only Q = J/γ₀ is measurable from inside ([PRIMORDIAL_QUBIT](../hypotheses/PRIMORDIAL_QUBIT.md) §9).
  Scaling (γ₀, J) → (λγ₀, λJ) leaves every inside observable unchanged.
  Information capacity is a function of Q alone: "γ₀ carries no
  information itself" becomes a consequence of dimensional analysis,
  not an interpretive choice.

### The cavity's response

- The information is in the CAVITY RESPONSE to γ₀: the standing wave
  pattern, determined by J and topology.
- What [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md) actually
  measured was not "information in γ" but "information in the cavity's
  response to γ".
- The response matrix has full rank because the Hamiltonian's eigenmodes
  interact with γ₀ in a mode-specific way ([absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md):
  Re(λ) = −2γ₀·⟨n_XY⟩), not because γ itself varies.
- The cavity's internal structure IS a light/lens duality: X/Y Pauli
  factors are "light" (coupled to γ₀), I/Z factors are "lens" (immune).
  Each mode has a fixed light content ⟨n_XY⟩ set by J and topology.
  The palindromic sum rule ⟨n_XY⟩_fast + ⟨n_XY⟩_slow = N pairs every
  exposed mode with an equally hidden partner. The "information" is the
  mode-by-mode light/lens distribution, not the illumination.
- The cavity has a parity regime ([OPTICAL_CAVITY_ANALYSIS](../experiments/OPTICAL_CAVITY_ANALYSIS.md)):
  even N chains are confocal (sharp Lorentzian spike, tight beam waist);
  odd N chains are defocal (broad Gaussian, central zero-mode). Under
  constant γ₀ the encoding regime therefore splits by N parity: one
  more structural axis, invisible to the γ-modulation reading.
- The cavity is a receive antenna array. [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md)
  identifies |+⟩^N as the phased-array receiver (full-rank SVD, 5
  channels); GHZ is omnidirectional with d_min = 0. Under γ₀ = const
  this becomes literal: J and topology are the antenna geometry, γ₀
  is uniform background illumination, the mode structure is the beam
  pattern. Alice designs antennas, not signals.

### Where information lives

- The bit is a Π-pair, not a single mode ([STANDING_WAVE_THEORY](../docs/STANDING_WAVE_THEORY.md)).
  The pair's decay rates sum to α_fast + α_slow = 2Σγ = 2Nγ₀ in the
  uniform case ([absorption theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)),
  so the natural channel time scale is 1/(2Nγ₀). The receiver reads the
  forward mode; its backward-decaying Π-image stabilizes the reception.
  The standing wave is the protocol.
- The channel's dimensional ceiling under γ₀ = const is structural,
  not noise-dependent. The number of independent modes is bounded by
  the cavity's degrees of freedom, and the 15.5 bits of
  [F30](../docs/ANALYTICAL_FORMULAS.md) at N=5 (5 independent SVD
  channels) sit inside that bound. J and topology set the ceiling;
  γ₀ sets the time axis. The Shannon capacity at finite SNR is a
  separate question: scaling with N becomes a pure cavity question.
- Pure-lens modes (Pauli strings with no X or Y factor) live in a
  decoherence-free subspace: they see no γ₀ at all, couple to nothing,
  and carry no information about the cavity's response. They are the
  system's *memory*. Light modes (≥1 X/Y factor) are the system's
  *voice*. Information is voice, not memory.

### Operational scope

- All structural control levers under γ₀ = const are J-side: coupling
  strengths, topology, and cavity timing through Q. Without shared
  J-coupling there is no channel at all ([BRIDGE_CLOSURE](../experiments/BRIDGE_CLOSURE.md)).
  "γ₀ flows constantly" does NOT mean "signals flow constantly": the
  ambient light is on, but only shared cavities carry messages.
- The channel is introspective, not receptive. What an observer
  "receives" from γ₀ is their own J-and-topology reflected back through
  the cavity: self-interference under uniform illumination. "What the
  observer knows" = which modes their J, topology, and boundary
  conditions populate. The readable subspace equals the populated
  subspace: self-knowledge through shadow, operationally.
- The channel exists only during the transient. Asymptotically, γ₀
  drives each U(1) sector to the maximally mixed state
  ([EQ-005](EMERGING_QUESTIONS.md#eq-005)): γ₀ erases information from
  inside the cavity at the sector-rate. Information lives in the
  time-dependent mode trajectory before full relaxation, bounded by
  the slow-mode lifetime ~1/(2γ₀·⟨n_XY⟩_min). The signal window is
  finite and set by γ₀ itself.

This inverts the reading direction:

| [GAMMA_AS_SIGNAL](../experiments/GAMMA_AS_SIGNAL.md) reading (pre-γ₀) | γ₀ = const reading |
|---------------------------------------|--------------------------|
| Alice encodes in γ profile | Alice encodes in J (cavity structure) |
| Information is in the illumination | Information is in the shadow pattern |
| The channel carries the message | The resonator IS the message |
| Noise is a signal from outside | Noise is constant light; structure is signal |
| Signal = γ variation across sites | Signal = light/lens distribution across modes |
| Capacity depends on γ variation | Capacity depends on Q = J/γ₀ |
| Sender modulates γ | Sender modulates J |
| Bit time ~ 1/γ_max | Bit time ~ 1/(Nγ₀) (standing-wave round trip) |

The [RESONANCE_NOT_CHANNEL](../hypotheses/RESONANCE_NOT_CHANNEL.md) insight becomes even more literal:
"The system is a soundbox, not a telephone": the sound (γ₀) is
constant. The music (information) comes from the shape of the soundbox.

In engineering terms: a linear system under constant Markovian drive.
The transfer-function-like object is the Liouvillian resolvent
G(z) = (z − L)⁻¹, with complex poles at L's eigenvalues (decay rates
and oscillation frequencies). γ₀ is the drive strength, J shapes the
filter, the signal is the output. This dissolves the misreading that
"noise carries the signal": the drive is uniform, the filter is the
Liouvillian, and the signal is the cavity's response to constant
illumination.

Under γ₀ = const, the thirteen pieces cohere into one picture. γ₀ is
the constant ambient carrier (1, 3, 6, 11), and the cavity that shapes
it is quantitative and time-reversal-symmetric (4, 7, 10). The
information IS the light/lens distribution across the cavity's modes
(2, 5, 8, 9), readable only via active J-engineering (12), with no
channel existing outside the cavity (13).

## Physical implications

Two points deserve separate flagging (the others are already in the
reinterpretation bullets above):

1. **The PTF closure law is an empirical regularity, not a conservation law.**
   Σ_i ln(α_i) ≈ 0 ([PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md))
   was conjectured to be a symmetry-protected invariant that a
   J-perturbation preserves exactly. [EQ-014](EMERGING_QUESTIONS.md#eq-014)
   closed the theorem sub-question by direct computation
   (details: [`EQ014_FINDINGS.md`](EQ014_FINDINGS.md)): the first-order
   coefficient Σ f_i = lim_{δJ→0} Σ ln(α_i(δJ))/δJ is **not zero and
   depends on the initial state**, with values from +0.05 (ψ_2) to +1.29
   (|+⟩⁷) at N=7, bond (0,1). The empirical ±0.05 tolerance PTF observed
   at |δJ| ≤ 0.1 is a combination of state-dependent first-order
   coefficients and partial second-order cancellation, not a conservation
   law. *What this means for the channel reading:* "J-perturbations
   conserve something that couples to information capacity" is falsified
   as a strict identity. It may survive as an approximate statement for
   particular initial states (ψ_2's tiny coefficient is striking and
   unexplained), but the "closure law IS the channel's conservation law"
   bridge is not available on structural grounds. The surviving question
   is why the first-order coefficient is small for ψ_2 and order-unity
   for |+⟩⁷ (tracked as EQ-014's surviving sub-question).

2. **Same light, different instruments.** Two observers with the same
   cavity see the same shadow; different J = different cavity =
   different information content. "Different observers" under
   γ₀ = const always means different cavities, never different
   illumination.

## What needs investigation

- [EQ-024](EMERGING_QUESTIONS.md#eq-024): J-modulation channel capacity. **Closed operationally 2026-04-23 evening:** C ≤ 12.07 bits at N=5 Heisenberg over F71-symmetric receivers, with structural surplus (three-class decomposition of the J-blind set; Class 3 SU(2)-Heisenberg-specific). See Update 2026-04-23 (evening) above. Surviving sub-questions tracked in EQ-024: three-class completeness, F71-breaking receiver capacity, N-scaling of the 12-bit ceiling, chromaticity of the Nelder-Mead optimum, operational meaning of the 12-vs-15 gap.

- Connection to PTF ([PERSPECTIVAL_TIME_FIELD](../hypotheses/PERSPECTIVAL_TIME_FIELD.md), [EQ-014](EMERGING_QUESTIONS.md#eq-014)):
  originally asked whether Σ_i ln(α_i) = 0 is a conservation law that
  couples to channel capacity. Answered negatively under EQ-014's
  2026-04-19 closure: not a first-order theorem, the coefficient Σ f_i
  depends on the initial state. See "Physical implications" point 1
  above. A weaker-form question survives: is there a state-specific
  invariant (valid for ψ_2-class states but not for |+⟩⁷) that could
  still tie to cavity response? This would be a restricted-initial-state
  version of the original conjecture.

- Connection to [INCOMPLETENESS_PROOF](../docs/proofs/INCOMPLETENESS_PROOF.md) Section 6 Option 2: "Noise is
  axiomatic. Like the speed of light." If γ₀ is axiomatic AND
  is the information channel, then the channel itself is axiomatic:
  not derived, not constructed, simply present. The bridge is not built.
  It IS.

- [EQ-015](EMERGING_QUESTIONS.md#eq-015): the cavity-mode-exposure formula γ_eff = γ₀ · |a_B|² ([F64](../docs/ANALYTICAL_FORMULAS.md)) on non-chain topologies and non-uniform J at N ≥ 5 is the falsification anchor for piece 6 ([PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md)). The uniform-J chain case is effectively closed via [F65](../docs/ANALYTICAL_FORMULAS.md) (verified N=3..30 to machine precision). **EQ-015 fully closed 2026-04-24** via [F64_TOPOLOGY_GENERALIZATION](../experiments/F64_TOPOLOGY_GENERALIZATION.md): F64 holds on chain, star, ring, complete graph, Y-tree at N=5 and N=7, both XY and Heisenberg, all B, provided degenerate perturbation theory is applied within H-degenerate subspaces. Max rel err < 0.001 at γ/J = 0.01 uniform J; max rel err < 0.07 under random non-uniform J in [0.5, 1.5] over 30 configurations (29/30 below 0.02, the remaining one consistent with expected second-order PT corrections). γ/J ~ 1 regime still open as the last sub-case. F64 is now a graph-universal structural anchor for γ₀ = const.

- [EQ-014](EMERGING_QUESTIONS.md#eq-014) theorem sub-question closed
  2026-04-19: dense biorthogonal eigendecomposition of the full
  16384×16384 L_A was executed (C# `ptf` mode, 146 min zgeev), residual
  1e-6 after cluster-fix; first-order PT through the slow-mode bilinear
  expansion + exact RK4 ground-truth validation reproduced PTF's
  empirical α_i to 4 decimal places. Direct first-order coefficient
  extraction via δJ ∈ {0.1, 0.01, 0.001} scan confirmed Σ f_i ≠ 0 and
  state-dependent. The four original building blocks (Tr[V_L·Π_slow] = 0,
  purity decomposition, U(1) conservation, det(U) = 1) are all true; what
  does not hold is their *combination* producing Σ f_i = 0. Scripts:
  `simulations/eq014_step23_biorth.py`, `_step4567_closure.py`,
  `_validate_groundtruth.py`, `_first_order_from_rk4.py`,
  `_spectrum_check.py`. Full findings: [`EQ014_FINDINGS.md`](EQ014_FINDINGS.md).

- [EQ-017 Phase 2 hardware result](../data/ibm_chain_gamma0_april2026/):
  already attempted. On ibm_kingston, the chain-mode γ₀ signal is 40-80×
  below the device noise floor (gate errors, T1, readout). The framework
  γ₀ is operationally indistinguishable from zero at current hardware
  fidelity. The synthesis is not hardware-verifiable yet; any empirical
  grounding has to come from structural predictions (palindromic spectral
  components, [F64](../docs/ANALYTICAL_FORMULAS.md) exposure ratios),
  not direct γ₀ measurement.

- Dynamical-attractor test of PTF closure ([ORTHOGONALITY_SELECTION_FAMILY](../experiments/ORTHOGONALITY_SELECTION_FAMILY.md) §4a): time-resolved Σ_i ln(α_i(t)) under uniform perturbation. If the closure deviation decays at the ~4γ₀ scale to (near) zero, [EQ-014](EMERGING_QUESTIONS.md#eq-014) surviving sub-question (state-dependence of Σ f_i) gets a geometric interpretation as distance-to-attractor. Adjacent to [EQ-024](EMERGING_QUESTIONS.md#eq-024), not part of it.

---

*Not a conclusion. A thread to pull.*
