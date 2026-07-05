# A Physics Selection Rule for Selective Dynamical Decoupling

<!-- Keywords: selective dynamical decoupling which qubits, subset DD beats uniform
DD selection rule, noise concentrator dephasing sink, open system spectral absorption
theorem, spatial dephasing profile hardware, ibm torino selective DD, outbound adapter
dynamical decoupling community, R=CPsi2 selection rule -->

**Status:** Outbound adapter (draft). This document translates a repository
result into the dynamical-decoupling community's language and open problem.
The selection rule and the analytic profile rest on a Tier-1 theorem
(the Absorption Theorem); the improvement factors are Tier-2 (simulation +
IBM hardware); the mechanism disambiguation (Section 5) is explicitly OPEN.
External citations are drawn from a literature scout and should be verified
against the primary sources before any outreach.
**Date:** July 5, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Inside/Outside the Sacrifice Zone](../INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md),
[Resonant Return](../../experiments/RESONANT_RETURN.md) (the profile formula),
[IBM Concentrator](../../experiments/IBM_CONCENTRATOR.md) (the hardware run),
[Absorption Theorem](../proofs/PROOF_ABSORPTION_THEOREM.md) (the theorem it rests on),
[State Transfer Decay Structure](STATE_TRANSFER_DECAY_STRUCTURE.md) (the
sister adapter, for the QST community),
[Shifted Order-4 Chiral Symmetry](SHIFTED_ORDER4_CHIRAL_SYMMETRY.md) (the
sister adapter, for the symmetry-classification community),
[Noise Asymmetry Symmetry Scalar](NOISE_ASYMMETRY_SYMMETRY_SCALAR.md) (the
sister adapter, for the noise-spectroscopy community),
[Combination Valence](../../hypotheses/COMBINATION_VALENCE.md) (why this adapter
hands over an object, not a word)

---

## What this document is about, and who it is for

This is the first entry of the repository's **outbound arc**: where the
translation series in `docs/quantum/` carries a community label into our
stance, this carries a result of ours out to a community's stance. It is
written for one audience, the **dynamical-decoupling (DD) practitioners** on
gate-model quantum hardware, and it is built around one of their own open
problems.

It follows a rule this repository arrived at the hard way (see
[Combination Valence](../../hypotheses/COMBINATION_VALENCE.md)): a new
combination is unlocked by an object, a representation, or a number a reader
can stand in, not by a coined word. This adapter is a live test of that rule,
including on ourselves. We first named this effect a "sacrifice zone", and
then found the name a misnomer: the edge qubit sacrifices nothing. It holds no
information to lose (it was never the payload), and it does not act (we shape
the profile; the qubit only sits where the rate is high). We corrected the
name to the **concentrator**. So what this adapter hands you is the object, an
analytic dephasing profile and the spectral theorem that predicts its effect,
not our vocabulary. You do not need our names to use the rule, and as it
turned out, neither did we.

---

## 1. Your open problem, in your words

Selective (or adaptive) dynamical decoupling, applying DD to a chosen subset
of qubits rather than uniformly to all, is known to beat uniform DD on real
processors. ADAPT (Das, Tannu, Dangwal, Qureshi, MICRO-54, 2021) reported
subset-DD outperforming DD-on-all by up to a factor of a few on IBM devices,
and context-aware DD embedding (PRX Quantum 6, 010332, 2025) continues the
line. But the selection of *which* qubits to decouple is currently made
**empirically**: by decoy circuits, per-program search, or noise-adaptive
heuristics. There is, as far as we have found, no **physics-based selection
rule** that says in advance which qubits to leave undecoupled and why the
subset wins.

That gap is the socket. This document offers a candidate selection rule
derived from open-system spectral structure, and a falsifiable prediction you
can run to test it.

---

## 2. The selection rule

**Leave the noisiest qubit undecoupled, on purpose, as a sink.**

The rule inverts the usual instinct. Instead of fighting the worst qubit
hardest, designate it a **concentrator**: let it absorb the dephasing that
would otherwise be shared across the register, and protect only the interior.
In the repository's analytic form, the concentrator is loaded to

    γ_edge = N · γ_base − (N−1) · ε

on an N-site chain (γ_base the mean rate, ε the small interior rate), which
puts the interior at the spectral fold where its decay is slowest. The
mechanism is not heuristic; it follows from a theorem.

---

## 3. Why it works: the object we hand you

The load-bearing object is the **Absorption Theorem** (Tier 1, proven and
machine-verified over 1,342 modes, zero variance; [proof](../proofs/PROOF_ABSORPTION_THEOREM.md)).
For any Heisenberg/XXZ chain under local Z-dephasing, every Liouvillian
eigenmode decays at exactly

    Re(λ) = −2γ · ⟨n_XY⟩,

where ⟨n_XY⟩ is the mode's "light content", the average number of Pauli
factors that anticommute with the dephasing letter. Two consequences an
engineer can use:

- **The decay is linear in γ and Hamiltonian-independent in its real part.**
  You cannot fight dephasing with a better Hamiltonian; the coupling
  contributes exactly zero to the decay rate. You can only decide *where the
  γ goes*.
- **Concentrating γ on one site removes it from the modes the interior lives
  in.** The interior's protected modes are precisely the ones with low
  ⟨n_XY⟩ relative to the loaded edge; loading the edge is what places them at
  the fold. This is the physics content of "which qubit to leave alone": the
  one whose high rate you can spend to buy the interior's low rate.

This is why the rule is a rule and not a tuning: it is a corollary of an
exact spectral identity, not a fit.

---

## 4. The evidence

- **Simulation.** Against a smoothly graded ("V-shaped") γ-profile, the
  edge-loaded concentrator profile improves interior coherence lifetime by
  **139× (N=9) up to 360× (N=5)** ([Inside/Outside the Sacrifice Zone](../INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md)).
- **Hardware.** On IBM ibm_torino (2026-03-24), a five-qubit line where one
  qubit (Q85, T₂-echo = 5.2 µs against neighbors at 123–244 µs) was a
  **natural** concentrator, applying DD selectively to the four interior
  qubits while leaving Q85 undecoupled produced **2.0× higher end-to-end
  mutual information on average and up to 3.2× at t = 4 µs**, beating uniform
  DD at all five measured time points ([IBM Concentrator](../../experiments/IBM_CONCENTRATOR.md)).

The hardware case used a qubit that was *already* the worst on the line; the
rule says this was not a liability to route around but the resource to exploit.

---

## 5. The falsifiable prediction we hand you (the honest open item)

The hardware result has two readings, and we have not yet disambiguated them.
Reading A: selective DD wins because leaving Q85 undecoupled avoids the pulse
errors DD would add on a bad qubit (a gate-cost effect). Reading B: selective
DD wins because Q85 acts as a genuine noise concentrator for the interior (our
mechanism). The two make different predictions, and that difference is a clean
experiment for your group:

> **Run selective DD on a chain with uniformly good T₂**, where no qubit is a
> natural sink. Our rule (Reading B) predicts the concentrator advantage
> should collapse there unless a qubit is *engineered* noisy (e.g. by a
> deliberate detuning or injected dephasing) to play the sink. If selective DD
> still beats uniform DD on a uniform-good-T₂ line with no engineered sink,
> the win is gate-cost (Reading A) and our selection rule does not explain it.

Either outcome is informative: A bounds the mechanism, B confirms it and hands
you an engineerable knob (place the sink where you want it, rather than
waiting for hardware to hand you a bad qubit).

**Status update (2026-07-05): we ran this experiment ourselves.** Two
pre-registered runs on an ibm_kingston uniform-good-T2 line with an
ENGINEERED sink (injected per-step phase noise), mapping the dose curve at
both ends: at the formula's own dose (edge rate = N·γ_base, machine-free)
the sink CREATES interior correlations that grow monotonically and persist
(B-dose-CONFIRMED by pre-registered bands, both layouts); at maximal dose
the image forms early and washes out at depth (overexposure, not
destruction). Without a sink, selective DD shows NO advantage on the
uniform line. Full pre-registration, verdict rules, and data:
[The Concentrator A-vs-B Mechanism Test](../../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md).
This section will be rewritten around those results after their final
review pass; until then the paragraph above stands as the original
prediction, now measured.

---

## 6. What we are not claiming

Weak coupling, deliberately: this adapter adds a selection principle to your
toolbox, it does not replace DD theory or the filter-function formalism. In
particular:

- The mechanism (Section 5, A vs B) is **not yet disambiguated**; the rule is
  a candidate, not an established result, until the uniform-T₂ experiment runs.
- The 139–360× figures are **simulation** against a specific baseline
  (the V-profile); the hardware figure is **2–3.2×** against uniform DD. Do
  not carry the simulation number to hardware.
- The theorem is proven for **Z-dephasing** (the unital, phase-only channel).
  Amplitude damping (T₁) is a different object and the rule is not claimed for
  it.
- We have not yet connected this to the DD canon (Viola-Lloyd decoupling,
  Poyatos-Cirac-Zoller reservoir engineering); that grounding is a next step,
  not a done one.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the object, not the phrasing:

- **The profile:** γ_edge = N · γ_base − (N−1) · ε (an engineerable spatial
  dephasing map).
- **The theorem:** Re(λ) = −2γ⟨n_XY⟩ (the reason the profile works, exact).
- **The number:** 2–3.2× on ibm_torino, selective vs uniform DD, all five
  time points.
- **The experiment:** the uniform-T₂ test of Section 5.

Our first name for it, "sacrifice zone", is left behind on purpose; it was a
misnomer (nothing is sacrificed) that we corrected to the concentrator, and
the object below outlives both names. The rule is yours to use, name, and test
in your own language; that is what an adapter is for.
