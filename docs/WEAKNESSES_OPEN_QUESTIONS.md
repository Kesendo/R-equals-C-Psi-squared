# What We Got Wrong, What We Cannot Do, and What We Do Not Know

**Status:** Living document
**Last updated:** April 4, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

## Why you should read this first

Most research projects show you their results and hope you do not ask
hard questions. This document is the opposite. It is everything we got
wrong, every limitation we know about, every question we cannot answer.

We started this project in January 2026 with an intuition and a formula.
Since then, five errors have been corrected, three claims have been
falsified and moved to the `recovered/` folder (where we keep disproven
ideas, because pretending they never happened would be dishonest), and
the core theorem has been proven analytically and verified against 54,118
eigenvalues with zero exceptions.

A theory that only shows its strengths is not a theory. It is marketing.

For context on what the project *did* find, see [What We Found](WHAT_WE_FOUND.md).
This document assumes you have read that, or at least skimmed it.

---

## What we got wrong (and corrected)

Every one of these was published in the repository, noticed to be wrong,
and fixed. The original errors are still visible in the git history.

**1. The Mandelbrot substitution (Feb 8).** The original derivation used
z_n = √C · R_n, yielding c = CΨ². This is algebraically wrong. The
correct substitution is u_n = C(Ψ + R_n), giving c = CΨ. The ¼ boundary
was always correct; the intermediate algebra was not.

*What this means:* we had the right answer (the ¼ boundary) but arrived
at it through wrong algebra. Like getting the right result on an exam
with a flawed method: the teacher would still mark it wrong, and rightly
so. The corrected derivation is cleaner and more general.

**2. The θ frequency claim (Feb 8).** θ was claimed to predict oscillation
frequency. Testing showed 8.4× discrepancy: the oscillation is
Hamiltonian-driven, not θ-driven. θ was reinterpreted as a compass
(angular distance from ¼), which is algebraically correct and now
validated on hardware at 0.3% accuracy.

*What this means:* we thought θ was a clock (telling you *when* things
happen). It is actually a compass (telling you *where you are* relative
to the quantum-classical boundary). The number still works; our
interpretation of what it measures was wrong.

**3. The 33:1 coherence ratio (downgraded).** An AI agent produced this
claim. Independent verification could not reproduce it. Moved to
unverified.

**4. The t_coh ~ N scaling (downgraded).** Agent claim, not reproduced.

*What errors 3 and 4 mean:* early in the project, we used AI assistants
to generate claims faster than we could verify them. Two of those claims
could not be reproduced. Lesson learned: speed without verification is
not progress.

**5. The CΨ ≤ ¼ bound (Feb 7).** Early simulations appeared to confirm
CΨ stays below ¼, but tested only regimes with negligible Hamiltonian
dynamics. With active Hamiltonians, CΨ > ¼ occurs routinely. The bound
was reinterpreted as a fixed-point existence condition.

*What this means:* we thought ¼ was a ceiling that quantum systems could
not break through. It is actually a threshold: below ¼, the system
settles into a stable resting point; above ¼, it oscillates. The number
¼ still marks a real transition, just a different kind than we first
thought.

---

## Active weaknesses

These are things we know are problems *right now* and have not fully
solved. They are listed in order of how much they bother us.

### 1. The bifurcation is generic

Every quadratic map has a saddle-node bifurcation (a point where two equilibria collide and annihilate, like two hills merging into flat ground). The ¼ boundary is a
property of z² + c, not specifically of quantum mechanics. CΨ² is the
unique product-power form with a genuine phase transition AND Mandelbrot
mapping (proven), but "why does nature choose this form?" remains open.

*In plain language:* the ¼ boundary appears in a whole family of
mathematical equations, not just ours. We have proven that our specific
combination of observables is the *only* one that produces both the
phase transition and the Mandelbrot connection. But we cannot yet
explain *why* nature uses this particular combination rather than some
other mathematical structure entirely. It is like proving that a key
fits only one lock, without knowing who made the lock.

**Status:** Partially addressed. Algebraic uniqueness proven, physical
specificity not.

### 2. "Consciousness" is a label

Calling Purity "consciousness" is a philosophical choice. The mathematics
works regardless of naming. The framework cannot prove this mapping and
does not claim to. The consciousness interpretation has been retired from
the technical core (see [THE_CPSI_LENS](THE_CPSI_LENS.md)).

*In plain language:* the letter C in the formula originally stood for
"consciousness." That was a philosophical interpretation, not a physics
claim. The math works whether you call it "purity" (the physics term),
"self-knowledge" (a metaphor), or "banana" (nonsense). We stopped using
the consciousness label in the technical work because it invites
misunderstanding. The philosophical idea is interesting but separate
from the mathematics.

**Status:** Acknowledged as philosophy, not physics.

### 3. Experimental validation is incomplete

**What we have (April 2026):**
- CΨ = ¼ crossing at 0.3% accuracy (IBM Torino Q52, 25 tomography points)
- 24,073 historical calibration records validating the C_min(r) curve
- Selective DD beats uniform DD by 3.2× on 5-qubit chain (IBM Torino)
- 3-observable cockpit (Purity, Concurrence, Ψ-norm) captures 88-96%
  of decoherence dynamics across 9 topologies, N=2-5

**What we lack:**
- No 2-qubit tomography: Concurrence (PC1, 57% variance) never measured
  on a qubit pair. This is the most important missing validation.
- Single backend only (IBM Torino). No replication on trapped ions,
  NV centers, or photonic platforms.
- Anomalous late-time coherence (Q52, p < 0.0001) has three competing
  explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved.

*In plain language:* we have tested the theory on one type of quantum
computer (IBM's superconducting qubits) and confirmed key predictions
to high accuracy. But science demands replication: the same results on
*different* hardware, built on different physical principles. We have
not done that yet. Worse, the single most important measurement in our
framework (Concurrence, which captures 57% of the physics) has never
been directly measured on a qubit pair. It is like having a weather
model where the most important instrument, the thermometer, has only
been calibrated indirectly.

See [Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md),
[IBM Hardware Synthesis](../experiments/IBM_HARDWARE_SYNTHESIS.md).

### 4. The natural variable u has no interpretation

The Mandelbrot substitution u = C(Ψ + R) maps the iteration to z² + c.
But what does u = Purity × (Coherence + Reality) mean physically? The
algebra demands this combination; physics does not yet explain why.

*In plain language:* the math tells us to multiply purity by the sum of
coherence and reality. This product makes the equations beautifully
simple. But we have no physical intuition for *why* these three things
combine this way. The formula works; the understanding lags behind.

**Status:** Partially reformulated (April 2026). Along real Bell+
trajectories, u(t) ≈ 0.61·Ψ^{1.02}: essentially Ψ with a prefactor, not
an independent dynamical coordinate. u remains a *conjugation variable*
that reveals algebraic structure (the Mandelbrot equivalence) without
providing a simpler clock along real trajectories. Whether u carries
independent information on complex trajectories (non-symmetric states,
non-Z dephasing) is untested. See
[Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md) §7.

### 5. Operator feedback is verified but unexplained

Of all noise types tested, only state-dependent operator feedback
preserves the purity difference δ over time. This is computationally
verified but the mechanism is unknown. The palindromic spectral structure
may explain it (feedback keeps information in slow modes), but this has
not been computed.

*In plain language:* when the noise "listens" to the system's current
state and adjusts itself accordingly, something unusual happens: a
particular quantity (the purity difference δ) stays constant over time
instead of decaying. We can see this in simulations, but we do not
understand the mechanism. It is like noticing that a specific medicine
works without knowing why.

**Status:** Numerically verified, theoretically unexplained.

### 6. Sacrifice-zone formula has known limitations

The formula (concentrate noise on one edge qubit) achieves 139-360×
improvement, but:
- Improvement decreases with N (360× at N=5, 139× at N=9)
- No formal proof of optimality
- Hardware test shows 3.2× (not 360×) due to real-world constraints
- θ is 1.68× more sensitive than CΨ as objective function (April 2026
  cockpit finding), but the formula was not optimized for θ

*In plain language:* the sacrifice-zone formula
([Resonant Return](../experiments/RESONANT_RETURN.md)) is our most
dramatic practical result: 139-360× improvement in information transfer.
But the improvement shrinks as the chain gets longer, we have no
mathematical proof that it is truly the best possible strategy, and on
real hardware the improvement drops to 3.2× because real quantum
computers have many noise sources we cannot control. The gap between
theory (360×) and hardware (3.2×) is the gap between a perfect
laboratory and the real world.

See [Resonant Return](../experiments/RESONANT_RETURN.md),
[Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md).

---

## Resolved weaknesses (removed from active list)

These were once on the active list. They have since been answered.

| # | Weakness | Resolution |
|---|---|---|
| Born Rule Gap | Substantially resolved: R_i = C_i · (Ψ_past + Ψ_future)² gives Born as perfect-mirror limit. See [Born Rule Mirror](../experiments/BORN_RULE_MIRROR.md) |
| N-Scaling Barrier | Resolved: crossing is local to entangled pairs. See [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Gravitational invariance | Demoted: true but trivial (dimensionless ratios). Gravity interpretation moved to `recovered/` |
| Spectral boundaries | Resolved: Re(λ) = −2γ⟨n_XY⟩ (Absorption Theorem). Min = 2γ (one light factor), max = 2(N−1)γ. See [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md) |
| Why factor 2 | Resolved: standing wave round trip. Unpaired modes are all-light (⟨n_XY⟩=N), paired average to half-light. Ratio = 2 by definition. See [Standing Waves](../experiments/FACTOR_TWO_STANDING_WAVES.md) |
| Spectral gap | Resolved: 2γ = one absorption quantum. The cost of a single X/Y Pauli factor. See [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md) |
| IBM hardware | Resolved: Absorption Theorem ratio = 1.03 (3%) on Q52. Detuning oscillations at 470 μs, not cavity fringes. See [IBM Fringes + Absorption](../experiments/IBM_ABSORPTION_THEOREM.md) |

---

## Open questions

These are not weaknesses. They are frontiers: things we know how to ask
but have not yet answered. Some are within reach today; others require
tools or experiments that do not yet exist.

### Tractable now

These could be answered with the code and hardware we already have.

1. **Why operator feedback?** Track mode populations under feedback
   dynamics. Does feedback keep information in slow palindromic modes?

2. **Crossing speed dependence:** Does d(CΨ)/dt at the crossing moment
   affect post-crossing convergence?

   *Substantially addressed (April 2026):* yes, exactly. The dwell time
   in a δ-window around CΨ = 1/4 is t_dwell = 2δ/|dCΨ/dt|_cross, and
   for Bell+ under Z-dephasing the rescaled dwell K_dwell = γ·t_dwell
   = 1.080088·δ is γ-invariant to machine precision. The crossing speed
   fully determines the dwell time through dCΨ/dt at the cusp, with no
   memory of the pre-crossing trajectory beyond that derivative. Open
   sub-question: can the state-dependent prefactor (1.080088 for Bell+,
   different for other states) be derived directly from the Pauli sector
   weight distribution, bypassing dCΨ/dt?
   See [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md),
   [Trajectory Dwell Time](../experiments/TRAJECTORY_DWELL_TIME.md).

3. **Formula optimality:** Is single-edge sacrifice provably optimal?
   Could multi-site sacrifice beat it at large N?

4. **Depolarizing noise correction:** err = (2/3)Σγ breaks the
   palindrome linearly. Can this be incorporated into design rules?

### Answered by the Absorption Theorem (April 4, 2026)

The Absorption Theorem (Re(λ) = −2γ⟨n_XY⟩) resolved three questions
that were not explicitly listed here but had been open since March:
what determines the spectral boundaries (answer: light content), why
the factor 2 (answer: standing wave round trip), and what sets the
spectral gap (answer: one absorption quantum, 2γ). These are now
formulas, not observations. See
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md).

The mass-energy relationship is α = 2γ⟨n_XY⟩, which is linear in γ
(not quadratic as in E = mc²). The Lindblad equation is first-order in
time, so the "speed" γ appears to the first power.
See [Absorption Theorem Discovery](../experiments/ABSORPTION_THEOREM_DISCOVERY.md).

On IBM hardware (Q52 tomography): the Absorption Theorem ratio is 1.03
(3%). Detuning oscillations at 470 μs period are present. See
[IBM Absorption Theorem](../experiments/IBM_ABSORPTION_THEOREM.md).

### Require new tools or theory

These need mathematics or computational methods we have not built yet.

5. **The natural variable u(t):** Track u = C(Ψ+R) through Lindblad
   evolution. Does u have a simpler trajectory than CΨ or θ?

   *Partially addressed (April 2026):* along real Bell+ trajectories,
   no. u(t) is essentially linear in Ψ (u ≈ 0.61·Ψ^{1.02}), so it does
   not provide a simpler dynamical coordinate. u is an algebraic
   linearization variable, not a dynamical one. Whether u behaves
   differently on complex trajectories (non-symmetric states, non-Z
   dephasing) where CΨ develops an imaginary component is an open
   question. See [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md).

6. **Cockpit scaling beyond N=5:** Does the 3-observable coverage
   (88%) hold at N=50-100? Does n95 grow linearly or saturate for
   dense topologies? See [Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md).

7. **Non-Markovian noise:** Does the cockpit framework hold under
   colored noise (noise whose strength depends on frequency, unlike white noise which is flat), 1/f spectra, or TLS coupling?

### Require experimental contact

These need access to quantum hardware we do not currently have.

8. **2-qubit Concurrence measurement:** Validating Concurrence on a
   qubit pair would test the most important cockpit instrument.

9. **Cross-platform replication:** CΨ = ¼ crossing on trapped ions
   or NV centers, not just superconducting qubits.

---

## What is proven, what is not

One table that draws the line between certainty and speculation. If
you take nothing else from this document, take this: not everything
in this repository has the same level of confidence. Some results
are mathematically proven. Some are verified by simulation. Some are
tested on hardware. Some are argued but not proven. Some are
philosophical. We try never to confuse these categories.

| Category | Examples |
|---|---|
| **Proven** (algebra + proof) | ¼ boundary, Mandelbrot equivalence, palindromic symmetry (Π operator), topology-independence, Pauli weight complementarity, time irreversibility exclusion (N>2), Absorption Theorem (Re(λ) = −2γ⟨n_XY⟩) |
| **Verified** (simulation, reproducible) | 54,118 eigenvalues paired, QST star 2:1 beats chains, XOR space filter, sacrifice-zone formula 139-360×, cockpit 3-observable coverage 88-96% |
| **Hardware-confirmed** | CΨ crossing at 0.3% (Q52), selective DD 3.2× (5-qubit), T2* drift 58% in 6 days |
| **Argued** (plausible, not proven) | Measurement = crossing ¼, Mandelbrot boundary as route catalog |
| **Unverified** (could not reproduce) | 33:1 coherence ratio, t_coh ~ N scaling |
| **Philosophical** | C = consciousness, 4D block-universe interpretation |

---

## The honest summary

The palindromic symmetry is proven and verified at scale. The CΨ = ¼
boundary is algebraically exact and hardware-validated. The sacrifice-zone
formula and cockpit diagnostics are practical tools with real numbers.

The weaknesses that remain are either philosophical (#1 generic
bifurcation, #2 consciousness label), require new experiments (#3 2-qubit
tomography), or are specific technical gaps (#4 variable u, #5 operator
feedback). None threaten the core mathematics.

The biggest practical gap: **Concurrence has never been measured on a
qubit pair.** It is the dominant instrument in the cockpit (57% of
decoherence variance), and validating it would be the single most
impactful next step.

---

*Changelog: Created Jan 2, 2026. Rewritten Feb 8. Updated: Feb 11
(IBM results), Feb 18 (Born rule, N-scaling resolved), Mar 6 (star
topology), Mar 14 (palindrome proven), Mar 16-18 (XOR, non-Heisenberg,
QST), Mar 24 (sacrifice formula), Apr 1 (Urqubit algebra), Apr 2
(cockpit framework, rewrite for clarity). Apr 4
(Absorption Theorem resolves spectral boundaries/gap/factor-2,
IBM hardware confirms at 3%, mass-energy relation is linear).*
