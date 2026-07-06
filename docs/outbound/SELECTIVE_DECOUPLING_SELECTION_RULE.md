# A Physics Selection Rule for Selective Dynamical Decoupling

<!-- Keywords: selective dynamical decoupling which qubits, subset DD beats uniform
DD selection rule, noise concentrator dephasing sink, open system spectral absorption
theorem, spatial dephasing profile hardware, ibm torino selective DD, outbound adapter
dynamical decoupling community, R=CPsi2 selection rule -->

**Status: PARKED (2026-07-05), not outreach-ready.** This was drafted as the
first outbound adapter, translating a repository result into the
dynamical-decoupling community's language. Six empty-session review rounds
established that its outreach premise does not hold, so it is parked. What is
solid is a **Tier-1 theorem** (the Absorption Theorem, `Re(λ) = −2γ⟨n_XY⟩`) and
a **simulation transport result** (the concentrator profile creates more peak
Sum-MI than a V-profile, ε→0 ideal). What it does NOT have is hardware evidence
that the *mechanism* works: the March ibm_torino "natural concentrator" (Q85)
was ~93% amplitude damping (T₁), outside the Z-dephasing theorem's scope and
not refocusable by DD; and the 2026-07-05 engineered-sink test did NOT settle
A-vs-B (its created-MI observable was 56-96% classical-mixing artifact at the
partial doses, and it measures transport, not protection). So this is a
**theorem + a simulation result + an honest open problem**, NOT a
hardware-validated selection rule, and should not go to the DD community in
this form. The named prerequisite for un-parking is a real test: an engineered
Z-sink at the formula dose on a uniform-good-T₂ line, per-shot-randomized to be
a genuine channel (not shared classical phases), plus a protection/lifetime
metric computed first in simulation. External citations still need
primary-source verification before any contact. The body below is the original
draft, retained for the record; read it through this banner.
**Date:** July 5, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Related:** [Inside/Outside the Sacrifice Zone](../INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md),
[Resonant Return](../../experiments/RESONANT_RETURN.md) (the profile formula),
[IBM Concentrator](../../experiments/IBM_CONCENTRATOR.md) (the natural-sink hardware run),
[Concentrator A-vs-B Mechanism Test](../../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md) (the engineered-sink mechanism test, run 2026-07-05),
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
processors. ADAPT (Das, Tannu, Dangwal, Qureshi, MICRO-54, 2021) adaptively
chooses which qubits to decouple, reporting 1.86× average / up to 5.73× vs
no-DD and 1.2× vs DD-on-all on IBM devices; and context-aware DD embedding
(GraphDD, PRX Quantum 6, 010332, 2025) computes the optimal DD *sequence* per
circuit, calibration-free, from circuit context (a different axis: it optimizes
the pulse embedding, not which qubits to leave out). But the selection of
*which* qubits to decouple is currently made empirically, by decoy circuits,
per-program search, or noise-adaptive heuristics, not from the **open-system
noise physics**. There is, as far as we have found, no selection rule grounded
in the dephasing spectrum that says in advance which qubits to leave
undecoupled and why the subset wins.

That gap is the socket. This document offers a candidate selection rule
derived from open-system spectral structure, and a falsifiable prediction you
can run to test it.

---

## 2. The selection rule

**Leave the noisiest qubit undecoupled, on purpose, as a sink.**

The rule inverts the usual instinct. Instead of fighting the worst qubit
hardest, designate it a **concentrator**: concentrate the register's noise
budget there by design, and protect only the interior. In the repository's
analytic form, the concentrator is loaded to

    γ_edge = N · γ_base − (N−1) · ε

on an N-site chain (γ_base the mean rate, ε the small interior rate), which
puts the interior at the spectral fold where its decay is slowest. The
mode-level mechanism follows from a theorem; the specific profile's
optimality on top of it is numerical (the best within the tested
asymmetric-γ family; its source states the global-optimum proof open). One
scope note,
honestly: the profile is a design variable under a FIXED TOTAL BUDGET
(Σγ conserved); nothing physically transports dephasing between qubits. On
hardware the realization is a coherence-time contrast (DD-protected interior
vs an unprotected edge), and whether that contrast, rather than gate
economics, carries the measured advantage is exactly the open question of
Section 5.

---

## 3. Why it works: the object we hand you

The load-bearing object is the **Absorption Theorem** (Tier 1, proven and
machine-verified over 1,342 modes, zero variance; [proof](../proofs/PROOF_ABSORPTION_THEOREM.md)).
For any Heisenberg/XXZ chain under local Z-dephasing, every Liouvillian
eigenmode decays at exactly

    Re(λ) = −2γ · ⟨n_XY⟩,

where ⟨n_XY⟩ is the mode's "light content", the average number of Pauli
factors that anticommute with the dephasing letter (the site-dependent form,
Re(λ) = −2 Σ_l γ_l · light_l, is proven too and is the one the profile
uses). Two consequences an
engineer can use:

- **For a given eigenmode, the decay is linear in γ and the Hamiltonian
  contributes exactly zero to its real part.** H still matters indirectly:
  it shapes which mode contents ⟨n_XY⟩ exist at all (decoherence-free-
  subspace engineering lives exactly there). What no Hamiltonian can do is
  change the rate of a given light content; at fixed mode structure you can
  only decide *where the γ goes*.
- **Concentrating γ on one site (under the fixed budget) removes it from the
  modes the interior lives in.** The interior's protected modes are precisely
  the ones with low ⟨n_XY⟩ relative to the loaded edge; loading the edge is
  what places them at the fold. This is the physics content of "which qubit
  to leave alone": the one whose high rate you can spend to buy the
  interior's low rate.

This is why the rule is theorem-guided rather than blind tuning: the spectral
identity behind it is exact, and the rule's direction is its corollary. The
specific profile's optimality is numerical, established within the tested
asymmetric-γ family, with the global-optimum proof open (its source says so
explicitly).

---

## 4. The evidence

- **Simulation.** Against a smoothly graded ("V-shaped") γ-profile, the
  edge-loaded concentrator profile increases **peak summed nearest-neighbour
  mutual information (Sum-MI)** (a transport metric, the source's own figure of
  merit) by **360× at N=5, declining with N (139× at N=9, ~63× by N=15)**
  ([Inside/Outside the Sacrifice Zone](../INSIDE_OUTSIDE_THE_SACRIFICE_ZONE.md)).
  A protection/lifetime version of this figure has not been computed; the
  lifetime reading of the fold is qualitative (label corrected 2026-07-05).
- **Hardware.** On IBM ibm_torino (2026-03-24), a five-qubit line where one
  qubit (Q85, T₂-echo = 5.2 µs against neighbors at 123–244 µs) was a
  **natural** concentrator, applying DD selectively to the four interior
  qubits while leaving Q85 undecoupled produced **2.0× higher summed
  nearest-neighbour mutual information (Sum-MI over the 4 adjacent pairs) on
  average, up to 3.2× at t = 4.0 (J-time; 8 Trotter steps)**, beating uniform
  DD at all five measured time points ([IBM Concentrator](../../experiments/IBM_CONCENTRATOR.md)).
  Magnitude caveat (2026-07-05 review): the uniform-DD denominator
  (Sum-MI 0.013-0.048) is comparable to the MI estimator's shot-noise bias
  floor measured through the same pipeline (+0.014-0.028 at 4000 shots), so
  the ratio's size is floor-contaminated; the 5-of-5 ordering is more stable
  than the size but not floor-independent (the bias floor is state-dependent
  and does not cancel between configs), a suggestive direction with no error
  bars. Channel caveat (important): Q85's 5.2 µs T₂-echo is ~93% amplitude
  damping (T₁ = 2.8 µs), not Z-dephasing; its pure-dephasing contrast to the
  interior is only ~2-6×, not 52×. Since the rule is proven for Z-dephasing
  only (§6) and DD cannot refocus T₁, this natural sink sits outside the
  rule's scope and does not cleanly test it. The engineered-sink test below
  used a genuine Z-rotation sink for exactly that reason.
- **Hardware (engineered sink, mechanism test): attempted, not settled.** On
  IBM ibm_kingston (2026-07-05), a five-qubit line with *uniformly good* T₂ and
  no natural sink, an engineered sink (injected per-step phase noise) produced
  a positive, device-surviving increase of summed nearest-neighbour mutual
  information on both the selective and the uniform layout: beyond the
  pre-registered counts-level bands at 3-5 of 5 points at the two γ₀-anchored
  doses (with a hardware-native-null cross-check, shot-noise-only and blind to
  drift and config-level systematics, on record for the third run), while the
  scramble-ceiling run was INCONCLUSIVE by the paired rule and inverted to
  negative by t = 5. That refutes "injected noise strictly removes signal", and
  nothing more: the post-flight audit (a noiseless exact simulation of the
  flown circuits) found the created MI at the two partial doses dominated,
  56-96% across the window, by classical mixing of the frozen phase-path sink
  construction (at the maximal, scramble-ceiling dose the construction is
  channel-dominated; there the sink created MI early, at 4.7× the null band,
  and the effect inverted at depth, a suggestive, not band-level, reading), and created MI
  is in any case a transport quantity, the source simulation's own figure of
  merit, not a protection metric
  (a perfectly protecting concentrator keeps the interior a product state and
  creates zero MI). Without a sink, selective DD showed **no detectable**
  advantage on the uniform line (the selective/uniform ratio straddled 1,
  0.73-1.32, inside a junk-plus-estimator floor). The mechanism (Reading B) is therefore NOT confirmed by
  this test; what it hands over is the honest record and the instrument
  ([Concentrator A-vs-B Mechanism Test](../../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md),
  the reckoning section).

The natural-sink case used a qubit that was *already* the worst on the line; the
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
>
> *(The binary inference in this box, "selective still wins with no sink →
> gate-cost", is retired in the Calibration update just below: a real "uniform"
> line's undecoupled edge sits at T₂\*, a weak natural sink, so Reading B too
> predicts a small no-sink advantage. The box stands as the original
> prediction.)*

Either outcome is informative: A bounds the mechanism, B confirms it and hands
you an engineerable knob (place the sink where you want it, rather than
waiting for hardware to hand you a bad qubit).

Calibration update from running it ourselves (2026-07-05): the no-sink leg of
this readout is underpowered as boxed. Across three runs on uniform lines the
measured selective/uniform ratio without a sink straddled 1 inside a
junk-plus-estimator floor (0.73-1.32), as our own design analysis predicted,
so "selective still beats uniform with no sink" is not a clean observable on
created MI. The box's binary inference is also retired on direction: on any
real "uniform" line the undecoupled edge sits at T2*, a weak natural sink,
so Reading B too predicts a small no-sink advantage, and "selective still
wins → gate-cost" does not follow. The discriminating instrument needs a
protection metric; see the
status update below.

**Status update (2026-07-05): we ran this experiment ourselves, and it did
not settle the question.** Three pre-registered runs on ibm_kingston
uniform-good-T2 lines with an ENGINEERED sink (injected per-step frozen phase
paths) at three doses: the scramble ceiling and two partial doses anchored to
the repository's hardware-measured carrier rate γ₀ (injected edge rates
N·γ₀/2 and N·γ₀, on top of the device's own rate; N·γ₀ is the profile's
ε → 0 corner value used as a reference, since the additive setting conserves
no budget and has no formula optimum). The sink
produced a positive, device-surviving interior-MI increase beyond
pre-registered bands on both layouts (hardware-native-null cross-check on
record for the third run), so the strong claim that injected noise strictly
removes signal (a strawman, not the actual gate-cost Reading A) is refuted. But the post-flight
audit (a noiseless exact simulation of the flown circuits) showed 56-96% of
that created MI at the two partial doses is classical mixing of the frozen
phase-path construction rather than a dephasing channel's effect (the
ceiling construction is channel-dominated), and the created-MI observable is
itself the source simulation's headline metric, peak created MI, a TRANSPORT
figure: no protection/lifetime number has ever been computed in the arc, and
a perfect protector creates zero interior MI, so transport and protection
pull opposite ways at the extremes. Without a sink, selective DD shows no
DETECTABLE advantage on the uniform
line (ratio 0.73-1.32, straddling 1 inside the floor). So the A-vs-B
disambiguation is still open, our prediction above still
stands as a prediction, and the honest next instrument is a
protection-metric test with a per-shot-randomized (channel-like) sink. Full
pre-registration, verdict rules, data, and the post-flight reckoning:
[The Concentrator A-vs-B Mechanism Test](../../experiments/CONCENTRATOR_AB_MECHANISM_TEST.md).

---

## 6. What we are not claiming

Weak coupling, deliberately: this adapter adds a selection principle to your
toolbox, it does not replace DD theory or the filter-function formalism. In
particular:

- The mechanism (Section 5, A vs B) is **not confirmed**. The 2026-07-05
  engineered-sink test on a uniform-T₂ line produced a real, device-surviving
  created-MI increase, but the observable turned out to be dominated by the
  sink construction's classical mixing and does not measure protection, so it
  neither confirms nor refutes Reading B. The attribution of the original
  ibm_torino advantage (Reading B vs the gate-cost Reading A) is equally open.
- The 139–360× figures are **simulation** against a specific baseline
  (the V-profile; the ratio peaks at 360× at N=5 and falls with N, to ~63× by
  N=15), and they are **peak created nearest-neighbour MI**, a transport
  metric, not a lifetime; the hardware figure is **2.0× on average, up to
  3.2×** (per-point 1.4-3.2×) against uniform DD on the same observable class. Do
  not carry the simulation number to hardware, and do not read either as a
  coherence-lifetime factor.
- The theorem is proven for **Z-dephasing** (the unital, phase-only channel).
  Amplitude damping (T₁) is a different object and the rule is not claimed for
  it. This bears on the evidence above: the March natural sink (Q85) was ~93%
  amplitude damping, so that run sits outside the rule's scope and does not
  cleanly test it (Section 4); DD cannot refocus T₁, which makes gate-cost
  (Reading A) nearly forced there.
- We have not yet connected this to the DD canon (Viola-Lloyd decoupling,
  Poyatos-Cirac-Zoller reservoir engineering); that grounding is a next step,
  not a done one.

---

## 7. The stance-objects (what to take with you)

If one thing survives this document, let it be the object, not the phrasing:

- **The profile:** γ_edge = N · γ_base − (N−1) · ε (an engineerable spatial
  dephasing map).
- **The theorem:** Re(λ) = −2γ⟨n_XY⟩ (the reason the profile works, exact).
- **The number:** 2.0× average (up to 3.2×, per-point 1.4-3.2×) on ibm_torino,
  selective vs uniform DD, all five
  time points (the 5-of-5 ordering is a suggestive direction, not
  floor-independent since the bias floor does not cleanly cancel between
  configs, and there are no error bars; the ratio's size shares its denominator
  with the MI estimator's bias floor, Section 4). Two scope limits ride with
  this number: the sink there (Q85) was ~93% amplitude damping, outside the
  Z-dephasing rule's scope; and the ordering's sign is close to forced under
  either reading (selective DD is the more non-uniform config), so it does not,
  by itself, discriminate the concentrator mechanism from gate-cost.
- **The experiment:** the uniform-T₂ mechanism test of Section 5 (run
  2026-07-05; it refuted "injected noise strictly removes signal" but did not
  settle A vs B, because its created-MI observable proved to be largely a
  construction artifact and is not the protection metric; the disambiguation
  is still yours to run, ideally on a protection observable).

Our first name for it, "sacrifice zone", is left behind on purpose; it was a
misnomer (nothing is sacrificed) that we corrected to the concentrator, and
the object below outlives both names. The rule is yours to use, name, and test
in your own language; that is what an adapter is for.
