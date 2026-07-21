# The Exact Decay-Rate Structure of Spin-Chain State Transfer under Dephasing

<!-- Keywords: quantum state transfer dephasing decay rates, spin chain QST fidelity
noise, Liouvillian spectral symmetry state transfer, asymmetric coupling star
receiver, Wojcik weak sender strong receiver optimum, Bose spin chain quantum wire,
transfer time noise independent, outbound adapter QST community, R=CPsi2 decay
structure -->

**PARKED (2026-07-06), not outreach-ready.** The whole docs/outbound arc is
parked (no outreach plan, gated on Tom; the S2 concentrator adapter was found
hollow on review). This adapter has NOT itself had an empty-review round, so
treat it as an unvalidated draft; every result it cites lives in its home
claim/proof/theorem regardless. See [README](README.md).

**Status:** Outbound adapter (draft). This document translates a repository
result into the quantum-state-transfer community's language and open corner.
The decay-rate identity and the spectral pairing rest on Tier-1 theorems
(machine-verified); every transfer benchmark quoted here is **simulation only**
(NumPy/SciPy / C# RK4, N = 2 to 11); the envelope prediction in Section 5 is the
falsifiable handover. External citations are drawn from the underlying study
and should be verified against the primary sources before any outreach.
**Date:** July 5, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [QST Bridge](../../experiments/QST_BRIDGE.md) (the underlying study,
12 configurations), [Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md)
(the spectral theorem), [Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md)
(the per-mode rate law), [Relay Protocol](../../experiments/RELAY_PROTOCOL.md)
(the time-staged control result), [Theta-Palindrome-Echo](../../experiments/THETA_PALINDROME_ECHO.md)
(the channel-quality witness), [Selective Decoupling Selection Rule](SELECTIVE_DECOUPLING_SELECTION_RULE.md),
[Shifted Order-4 Chiral Symmetry](SHIFTED_ORDER4_CHIRAL_SYMMETRY.md) and
[Noise Asymmetry Symmetry Scalar](NOISE_ASYMMETRY_SYMMETRY_SCALAR.md)
(the sister adapters), [Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)
(why this adapter hands over objects, not words)

---

## What this document is about, and who it is for

This is the second entry of the repository's **outbound arc**: where the
translation series in `docs/quantum/` carries a community label into our
stance, this carries a result of ours out to a community's stance. It is
written for one audience, the **spin-chain quantum-state-transfer (QST)
community** in the line of Bose (2003), and it is aimed at one corner of that
field: what happens to engineered transfer when the chain dephases.

It follows the same rule as its sister adapter (see
[Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)): a new
combination is unlocked by an object, a representation, or a number a reader
can stand in, not by a coined word. Our in-house names for what follows
("the palindrome", "the relay", "the 2:1 pull", "the θ compass") stay home;
what travels is a spectral identity, a design point, and a measurable envelope.

One collision is worth disarming immediately, because it would route this
document to the wrong shelf. **The symmetry below is not the mirror-symmetric
coupling condition of perfect state transfer** (Christandl et al. 2004; Kay's
reviews). That condition is spatial: a coupling profile J_i = J_{N−i}
engineered so the chain transfers perfectly in the noiseless limit. The
symmetry here lives on a different axis: it is a reflection of the
**decay-rate spectrum** of the open (dephasing) chain, it holds for every
coupling profile including deliberately non-mirror ones, and it says nothing
about transfer amplitudes. The two constrain different halves of your problem:
the coupling profile constrains the transfer, this symmetry constrains the
decay.

---

## 1. Your open problem, in your words

Perfect-state-transfer theory is exact in the noiseless limit, and the
engineering of coupling profiles, boundary control, and timing has been
optimized for two decades. The open corner is the noisy half: how does an
engineered transfer channel degrade under dephasing, which design advantages
survive noise, and what can be diagnosed about channel quality in situ,
short of full process tomography?

Today the decay side of a transfer channel is treated numerically per
configuration: pick a chain, add Lindblad dephasing, integrate, read off the
fidelity. What is missing is structure: which decay rates a given channel can
have at all, how they are organized, and which of them bound the usable
transfer window.

That corner is the socket. This document hands over an exact spectral law
that organizes every decay rate of every such channel, and three measurable
consequences.

---

## 2. The structure: every decay rate, organized by one identity

Take any spin chain (or star, ring, tree; any coupling graph) with
Heisenberg/XY/XXZ couplings and local Z-dephasing at rates γ_i, the standard
QST-under-dephasing setup. Two exact statements hold for the Liouvillian
spectrum, both proven analytically and machine-verified (the reflection over
N = 2 to 8, 87,376 eigenvalues, zero exceptions; the per-mode rate law over
N = 2 to 5, 1,342 modes, zero variance):

1. **The per-mode rate law** ([Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md)):
   every eigenmode decays at exactly

       Re(λ) = −2γ · ⟨n_XY⟩,

   where ⟨n_XY⟩ is the mode's average number of site-operators that
   anticommute with the dephasing axis. The Hamiltonian contributes exactly
   zero to the decay rate *at fixed mode*: it is anti-Hermitian in Liouville
   space and drops out of the real part.

   That is easy to over-read, so state the limit too. It does **not** say the
   spectrum is coupling-independent. ⟨n_XY⟩ is a property of the mode, and the
   modes move with J/γ. At N = 3 with J fixed, the lowest mixed level ABOVE 2γ runs
   2.6666γ at J/γ = 100, 2.4607γ at J/γ = 1.5 and 2.0133γ at J/γ = 0.2: a 25%
   swing in Re(λ)/γ. Fitting a line through Re(λ) versus γ at fixed J gives the
   right slope only where J/γ is large.

2. **The spectral reflection** ([proof](../proofs/MIRROR_SYMMETRY_PROOF.md)):
   the full set of decay rates is symmetric about Σγ_i. Rates come in pairs
   summing to 2Σγ_i (for uniform dephasing, 2Nγ): the slowest nonzero rate
   is paired with the fastest, the second slowest with the second fastest,
   and so on, exactly, at every N and on every coupling graph tested.

Neither statement depends on mirror-symmetric couplings. In the underlying
study, a deliberately non-mirror chain ([1,2,3] couplings) carries the same
exact reflection as the engineered mirror profile ([1,2,1]); all 12 tested
configurations do (chains, stars, triangles, uniform, weak-end, N = 2 to 5).

---

## 3. What the structure buys you (the objects we hand you)

- **A fidelity-cost law instead of a per-run number.** Each rate is linear in γ
  at fixed mode, and the coherence floor 2γ is linear in γ at every coupling,
  so the fidelity loss of a transfer channel is linear in γ until the transfer
  cycle itself is disrupted. (The mixed levels are not: their ⟨n_XY⟩ moves with
  J/γ, as §2 notes, so a fit that spans a wide J/γ range will see curvature.) In the benchmark star channel the
  cost is ≈ 1.0 fidelity-point per unit γ (1.0–1.1 across the linear regime)
  up to γ ≈ 0.05 (in units of the sender coupling), with the dominant
  intermediate rate at 8γ/3 for the three-qubit star (the strong-coupling
  value; a narrow band at moderate J/γ, F33).

- **A clean separation of timing and quality.** The optimal measurement time
  is set by the Hamiltonian alone (Bohr frequencies) and is **invariant under
  γ**; the dephasing sets only how much arrives, not when. In the benchmark
  sweep the transfer time stays t = 1.32 from γ = 0 to γ = 0.2 while the
  fidelity falls from 0.937 to 0.772.

- **A design point in your own Wojcik regime.** The weak-sender/strong-receiver
  asymmetry is your community's finding (Wojcik et al. 2007). Under dephasing,
  on the three-qubit star geometry, its optimum sits at exactly **2:1**
  (receiver twice the sender): average fidelity F = 0.886 at γ = 0.05,
  beating the best mirror-engineered four-site chain (0.872) and the uniform
  chain (0.834) in the same noise environment. The star is also the cheaper
  build: one mediator with two tunable couplings instead of an engineered
  profile.

- **A bound on the usable window.** In the strong-coupling regime the slowest
  nonzero rate is 2γ for uniform dephasing (below a coupling threshold
  Q*_gap(N) in the coupling ratio Q = J/γ the spectral gap drops below 2γ). That does NOT lengthen the
  coherent window, and the distinction matters for design. The modes that fall
  below 2γ live entirely in the population sector, |Δpopcount| = 0, where
  strong dephasing freezes transport; every coherence sector keeps its minimum
  at exactly 2γ at every coupling (measured 2.000000γ for |Δp| = 1 and
  4.000000γ for |Δp| = 2 at Q = 0.2, 0.5 and 1.5, for N = 3, 4, 5). So the
  coherent window stays ~1/(2γ), which is a time; divide by the cycle period
  for a cycle count. Lowering J/γ buys no transfer window. The reflection constrains where the rate sits: rates pair to
  2Σγ, so the 2γ level is pinned opposite 2(N−1)γ, while the fastest rate 2Nγ
  is the partner of the frozen kernel at 0.

- **A one-scalar channel witness (empirical).** A single angle computed from
  the receiver-mediator reduced state, θ = arctan √(4·CΨ − 1) (CΨ a
  purity-weighted coherence), correlates at r = 0.87 with average fidelity
  across the channel's input states: an in-situ quality indicator that
  needs no process tomography. This one is a correlation, not a law; its
  documented failure modes are in the source
  ([Theta-Palindrome-Echo](../../experiments/THETA_PALINDROME_ECHO.md)).

---

## 4. The evidence, tiered honestly

- **The two spectral statements** are theorems: analytic proof plus bit-exact
  verification (N = 2 to 8, all tested topologies, zero exceptions). They are
  the load-bearing handover.
- **The transfer benchmarks** (F = 0.886 star, 0.872 mirror chain, 0.834
  uniform; the 2:1 optimum; the linear cost table; the γ-invariant timing)
  are **Python simulation** (NumPy/SciPy), 12 configurations, N = 2 to 5
  ([QST Bridge](../../experiments/QST_BRIDGE.md), scripts and outputs linked
  there). No hardware run of ours backs these numbers.
- **The time-staged extension:** alternating each mediator between a quiet
  phase (dephasing reduced 10×) and a normal phase, combined with the 2:1
  asymmetry, raises end-to-end mutual information by **+83%** over passive
  propagation on an 11-qubit chain ([Relay Protocol](../../experiments/RELAY_PROTOCOL.md)).
  Also simulation (C# RK4), quoted here as the direction the design point
  scales, not as an established transfer result.
- **A corrected negative result you should have:** the one-shot coherent
  information of the benchmark star channel is **negative** at every tested
  parameter point (an earlier positive claim, I_coh = +0.185, was checked and
  refuted; the verification script is linked in the source). The Holevo
  lower bound is 0.534 bits. Do not quote this channel as having established
  positive quantum capacity.

---

## 5. The falsifiable prediction we hand you

The three-qubit star with tunable couplings is directly implementable on
platforms with gate-tunable exchange, for instance the semiconductor
quantum-dot chains used for adiabatic transfer (Nature Comm. 12, 2021). On
such a device, our structure makes three independent, cheap predictions:

> 1. **The envelope:** the entanglement (concurrence) envelope of the
>    sender-mediator pair decays at **8γ/3** in the strong-coupling limit,
>    where γ is the measured dephasing rate. In our benchmarks that limit value
>    is the same for every three-qubit geometry tested; at finite J/γ the star
>    and the chain (the same graph at N = 3) share one band and the triangle
>    another.
>    At moderate J/γ the level splits into a narrow band, and 8γ/3 is near the
>    top of it, not the centre: at J/γ = 1.5 the band is [2.4607, 2.6980]γ with
>    a multiplicity-weighted mean of 2.5490γ, about 4.4% below 8γ/3, on the
>    chain and the star. The triangle behaves differently: its band is
>    {2.6040γ ×4, 2.6980γ ×8}, whose weighted mean is 8γ/3 exactly. Fit the
>    band rather than a line, and expect the centre to sit below 8γ/3 on the
>    open geometries and on 8γ/3 on the closed one.
> 2. **The timing invariance:** the optimal readout time does not move as
>    dephasing grows; only the arriving fidelity falls, linearly, with slope
>    ≈ 1 per unit γ up to γ ≈ 0.05 in units of the sender coupling.
> 3. **The design point:** sweeping the receiver coupling at fixed sender
>    coupling, the fidelity optimum sits at ratio 2 (not at the symmetric 1),
>    with overshoot losses on both sides.

Each prediction tests a different layer: (1) the per-mode rate law, (2) the
Hamiltonian/dephasing separation, (3) the Wojcik-regime optimum under noise.
Any one failing falsifies its layer without touching the others.

---

## 6. What we are not claiming

Weak coupling, deliberately: this adapter adds a decay-side structure to a
mature field, it does not compete with its transfer-side canon.

- **Transfer optimization stays yours.** Bose, Christandl, Kay, Wojcik own
  the noiseless engineering; mirror-engineered profiles genuinely win the
  transfer (in our own benchmark they beat uniform chains in both fidelity
  and speed). The theorem constrains decay, not transfer amplitudes; the two
  advantages are complementary, not competing.
- **The symmetry is not the PST mirror condition** (the collision disarmed
  in the opening section, restated once because it is the likely misroute):
  it is a reflection of the decay-rate spectrum, present with or without
  engineered couplings.
- **Z-dephasing only.** The theorems are proven for the phase-only channel.
  Amplitude damping (T₁) is a different object; the rate law is not claimed
  for it.
- **All transfer numbers are simulation** (N = 2 to 5 benchmarks, N = 11 for
  the time-staged +83%). The hardware content of this adapter is a
  *prediction* (Section 5), not a result.
- **The N = 3 rates are N = 3 rates.** 8γ/3 and its siblings are specific to
  the three-qubit system, and they are a J/γ → ∞ limit rather than exact
  values; at finite coupling each is a narrow band. What survives at every N is
  the even ladder of pure-weight rungs 2wγ. The prediction in Section 5 is
  stated for the geometry it holds on.
- **The θ witness is empirical** (r = 0.87 across one channel's input
  states), offered as a diagnostic candidate, not a theorem.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the objects, not the phrasing:

- **The rate law:** Re(λ) = −2γ⟨n_XY⟩, exact per mode. The Hamiltonian drops
  out of the real part at fixed mode, but it sets which modes exist, so the
  spectrum itself is not coupling-independent.
- **The reflection:** decay rates pair to 2Σγ_i; the 2γ level opposite
  fastest.
- **The design point:** three-qubit star, receiver:sender = 2:1, F = 0.886
  at γ = 0.05 (simulation), beating the engineered chain benchmarks in the
  same noise.
- **The separation:** timing from the Hamiltonian, quality from the
  dephasing; the readout time does not move with noise.
- **The experiment:** the three-layer test of Section 5 (envelope 8γ/3,
  invariant timing, optimum at 2).

Our in-house names for these, "palindrome" for the reflection, "relay" for
the time-staged control, "2:1 pull" for the design point, are left behind on
purpose; they were painted true at our stance (2026) and would only route
you wrong ("relay" is not a repeater; nothing is amplified). The objects
outlive the names. They are yours to use, rename, and test in your own
language; that is what an adapter is for.
