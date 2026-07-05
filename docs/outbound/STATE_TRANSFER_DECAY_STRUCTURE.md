# The Exact Decay-Rate Structure of Spin-Chain State Transfer under Dephasing

<!-- Keywords: quantum state transfer dephasing decay rates, spin chain QST fidelity
noise, Liouvillian spectral symmetry state transfer, asymmetric coupling star
receiver, Wojcik weak sender strong receiver optimum, Bose spin chain quantum wire,
transfer time noise independent, outbound adapter QST community, R=CPsi2 decay
structure -->

**Status:** Outbound adapter (draft). This document translates a repository
result into the quantum-state-transfer community's language and open corner.
The decay-rate identity and the spectral pairing rest on Tier-1 theorems
(machine-verified); every transfer benchmark quoted here is **simulation only**
(QuTiP / C# RK4, N = 2 to 11); the envelope prediction in Section 5 is the
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
(the channel-quality witness), [Selective Decoupling Selection Rule](SELECTIVE_DECOUPLING_SELECTION_RULE.md)
(the sister adapter), [Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)
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
   anticommute with the dephasing axis. The real part is **linear in γ and
   independent of the couplings**: the Hamiltonian moves the oscillation
   frequencies (Im λ) and the mode composition, but contributes exactly zero
   to the decay rate itself.

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

- **A fidelity-cost law instead of a per-run number.** Because every rate is
  linear in γ, the fidelity loss of a transfer channel is linear in γ until
  the transfer cycle itself is disrupted. In the benchmark star channel the
  cost is ≈ 1.0 fidelity-point per unit γ (1.0–1.1 across the linear regime)
  up to γ ≈ 0.05 (in units of the sender coupling), with the dominant
  intermediate rate at 8γ/3 for the three-qubit star.

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

- **A bound on the usable window.** The slowest nonzero rate (2γ for uniform
  dephasing) caps the number of coherent cycles available to any echo-based
  or repeated-readout scheme at ~1/(2γ); the reflection tells you this
  slowest rate is not free to sit anywhere, it is pinned opposite the fastest.

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
  are **QuTiP simulation**, 12 configurations, N = 2 to 5
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
>    sender-mediator pair decays at **8γ/3**, where γ is the measured
>    dephasing rate; in our benchmarks this value is the same for every
>    three-qubit geometry tested (star, chain, triangle).
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
  the three-qubit system (only the boundary rate 2γ stays universal at
  N ≥ 4). The prediction in Section 5 is stated for the geometry it holds on.
- **The θ witness is empirical** (r = 0.87 across one channel's input
  states), offered as a diagnostic candidate, not a theorem.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the objects, not the phrasing:

- **The rate law:** Re(λ) = −2γ⟨n_XY⟩ (coupling-independent decay, exact).
- **The reflection:** decay rates pair to 2Σγ_i; slowest pinned opposite
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
