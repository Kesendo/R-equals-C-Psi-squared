# What We Got Wrong, What We Cannot Do, and What We Do Not Know

**Status:** Living document
**Last refreshed:** September 6, 2026 (the change history lives in git)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

## Why you should read this first

Most research projects show you their results and hope you do not ask
hard questions. This document is the opposite: a front door to the
limitations and open questions that still affect the current result set.
The running forensic record is [Caught Errors](CAUGHT_ERRORS.md).

We started this project in January 2026 with an intuition and a formula.
Since then, claims have repeatedly been corrected, narrowed, or deleted;
git keeps their history and [Caught Errors](CAUGHT_ERRORS.md) records the
material failures. The core theorem has been proven analytically and
verified against 87,376 eigenvalues with zero exceptions.

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
(angular distance from ¼), which is algebraically correct. θ itself has
not been measured as an angle on hardware; what has is the argument of
the complex CΨ, which wears the same square-root form
(`f95_angle_steering_kingston_may2026`, predicted against measured to
within 6.8° to 15.7°).

*What this means:* we thought θ was a clock (telling you *when* things
happen). It is actually a compass (telling you *where you are* relative
to the quantum-classical boundary). The number still works; our
interpretation of what it measures was wrong.

**3. The 33:1 coherence ratio (refuted).** The claimed ratio was not a
property of the stated protocol. Proper Lindblad simulations over 21
distributions contradicted it.

**4. The t_coh ~ N scaling (downgraded).** Agent claim, not reproduced.

*What errors 3 and 4 mean:* early in the project, we used AI assistants
to generate claims faster than we could verify them. Two of those claims
failed independent reproduction, and one was directly refuted. Lesson
learned: speed without verification is not progress.

**5. The CΨ ≤ ¼ bound (Feb 7).** Early simulations appeared to confirm
CΨ stays below ¼, but tested only regimes with negligible Hamiltonian
dynamics. With active Hamiltonians, CΨ > ¼ occurs routinely. The bound
was reinterpreted as a fixed-point existence condition.

*What this means:* we thought ¼ was a ceiling that quantum systems could
not break through. It is actually the real fixed-point boundary of the
scalar recurrence: in the physical range 0 ≤ c = CΨ < ¼ the attracting
real fixed point exists; above ¼ the real roots disappear and real
iterates grow rather than settle.
That statement is about the recurrence, not a universal physical phase
change. Physical trajectories can cross CΨ = ¼ in either direction.

---

## Active weaknesses

These are things we know are problems *right now* and have not fully
solved. They are listed in order of how much they bother us.

### 1. The bifurcation is generic

The one-parameter family z -> z² + c has a fold at c = ¼, where its two
real fixed points merge. That is a generic mathematical mechanism, not
something specific to quantum mechanics. Within the repository's stated
product-power/Renyi family, alpha = 2 is singled out by the combination
of the quadratic Mandelbrot form and a state-independent ¼ threshold.
Why that construction should be physically privileged remains open.

*In plain language:* the fold itself is standard mathematics. The
uniqueness result applies only inside a specified family of candidate
recurrences; it does not prove that nature selected that family.

**Status:** Partially addressed. Algebraic uniqueness proven, physical
specificity not.

### 2. "Consciousness" is a label

Calling Purity "consciousness" is a philosophical choice. The mathematics
works regardless of naming. The framework cannot prove this mapping and
does not claim to. The consciousness interpretation has been retired from
the technical core (see [The CΨ Lens](THE_CPSI_LENS.md)).

*In plain language:* the letter C in the formula originally stood for
"consciousness." That was a philosophical interpretation, not a physics
claim. The math works whether you call it "purity" (the physics term),
"self-knowledge" (a metaphor), or "banana" (nonsense). We stopped using
the consciousness label in the technical work because it invites
misunderstanding. The philosophical idea is interesting but separate
from the mathematics.

**Status:** Acknowledged as philosophy, not physics.

### 3. Experimental validation is incomplete

**What we have:**
- On IBM Torino, Q52 crossed qualitatively but 10.7% later than its T2*
  prediction; the separate Q80 crossing agreed to 1.9%
- 24,073 historical calibration records provide a mathematical
  consistency check for C_min(r), not an independent validation
- Selective DD beats uniform DD on the tested 5-qubit chain
- In the entangled-observer cockpit class, the first three principal
  components of 8-9 active standardized features capture 91.9-99.0%
  across N=5-11 over eight `center_bell` configurations (an adjacent
  central pair on the chain or two symmetry-equivalent leaves on the
  star, |+>^(N-2), uniform local Z-dephasing); Purity is the
  active feature most strongly correlated with PC1

**What we lack:**
- No direct two-qubit Concurrence tomography for the cockpit claim
- Results span several IBM superconducting backends, but not an
  independent trapped-ion, NV-centre, or photonic platform

*In plain language:* several results have hardware contact on IBM
superconducting devices, with accuracies that depend strongly on the
specific test. Independent physical platforms and direct Concurrence
tomography are still missing.

See [Cockpit Scaling](../experiments/COCKPIT_SCALING.md),
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
providing a simpler clock along real trajectories. Whether u supplies
independent information beyond that real recurrence remains untested.
Any complex-plane extension must use the separately defined CΨ_com;
non-symmetric states or non-Z dephasing do not by themselves make the
original CΨ complex. See
[Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md)
§8, "Active Weakness #4 (the natural variable u)".

### 5. Sacrifice-zone formula has known limitations

The ε→0 single-edge-sacrifice formula produced the following
ideal-simulation gains in peak created nearest-neighbour Sum-MI relative
to the corresponding V-shape baseline profile for N=5,7,9,11,13,15:
360×, 180×, 139×, 97×, 105×, 68×. But:
- The trend is broadly downward with N, with a non-monotone N=13 point
- No formal proof of optimality
- One Torino run on one 5-qubit chain gives five time-point contrasts from
  1.4× to 3.2× (average 2.0×); the ratio sizes are floor-sensitive, while
  the registered result is the selective-over-uniform ordering
- θ improves by 1.68× against CΨ's 1.28× under edge sacrifice (April 2026
  cockpit finding), and the formula was not optimized for θ. Whether that
  makes θ the better objective is open on its own terms: θ is a function
  of CΨ, its zero sits at the ¼ boundary, and the two ratios are taken
  from different baselines, so the larger one is not yet a sharper
  separation

*In plain language:* the sacrifice-zone protocol strongly increases one
specified simulation metric: the peak nearest-neighbour mutual
information created above its initial value. That is not yet a proof of
longer lifetime, protected end-to-end transfer, or optimality. One Torino run
shows a smaller but positive contrast at its five tested time points; it is
not a replication claim.

See [Resonant Return](../experiments/RESONANT_RETURN.md),
[Cockpit Universality](../experiments/COCKPIT_UNIVERSALITY.md).

---

## Resolved weaknesses (removed from active list)

These were once on the active list. They have since been answered.

| Weakness | Resolution |
|---|---|
| N-Scaling Barrier | Resolved: crossing is local to entangled pairs. See [Subsystem Crossing](../experiments/SUBSYSTEM_CROSSING.md) |
| Gravitational invariance | Scoped result retained in [Gravitational Invariance](../experiments/GRAVITATIONAL_INVARIANCE.md): the Bell+ Lindblad runs exhibit the stated dimensionless rate/time scaling. The interpretation of γ as gravitational field strength or a GR metric coefficient is retired there. |
| Spectral boundaries | Scoped result: for the uniform-Z Heisenberg chain above Q*_gap(N), the generic band has min 2γ and max 2(N−1)γ; the kernel and XOR drain lie outside it. The identity Re(λ) = −2γ⟨n_XY⟩ is broader than those edge formulas. See [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md) |
| Why factor 2 | Resolved within the F8 scope: for local Z-dephasing instances with Σγ > 0 satisfying the F1 palindromizer hypotheses whose spectrum reaches both endpoints, it is the ratio of the full range (0 to 2Σγ) to the centre (Σγ), not a ratio between two sets of modes. See [Absorption Theorem](proofs/PROOF_ABSORPTION_THEOREM.md) §4.4 and [Standing Waves](../experiments/FACTOR_TWO_STANDING_WAVES.md) Result 3 |
| Spectral gap | Resolved above Q*_gap(N): 2γ = one absorption quantum, the cost of a single X/Y Pauli factor. Below the threshold the gap is Zeno-suppressed and the theorem supplies no lower bound; see [Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md) §4.3 |
| IBM hardware | Resolved: Absorption Theorem ratio = 1.03 (3%) on Q52. Detuning oscillations at 470 μs, not cavity fringes. See [IBM Fringes + Absorption](../experiments/IBM_ABSORPTION_THEOREM.md) |

---

## Open questions

These are not weaknesses. They are frontiers: things we know how to ask
but have not yet answered. Some are within reach today; others require
tools or experiments that do not yet exist.

### Tractable now

These could be answered with the code and hardware we already have.

1. **Crossing speed dependence:** Does d(CΨ)/dt at the crossing moment
   control post-crossing convergence beyond a local neighbourhood?

   *Locally addressed (April 2026):* in a symmetric δ-window around
   CΨ = 1/4, t_dwell = 2δ/|dCΨ/dt|_cross is the leading small-δ estimate,
   with higher-order δ corrections. For Bell+ under Z-dephasing its
   leading rescaled value is K_dwell = γ·t_dwell = 1.080088·δ and is
   γ-invariant to machine precision. This local derivative does not yet
   determine the later post-crossing convergence. For the documented
   two-sector states, F58/F59 already derive the state-dependent local
   prefactor from sector weights; the open extension is the general
   multi-sector case and the genuinely post-crossing timescale.
   See [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md),
   [Critical Slowing at the Cusp](../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md) (Section 6).

2. **Formula optimality:** Is single-edge sacrifice provably optimal?
   Could multi-site sacrifice beat it at large N?

3. **Depolarizing noise correction:** Let γ_l be the total depolarizing rate
   on site l, divided as γ_l/3 among X, Y and Z. F5's extreme pair-sum
   shortfall is (2/3)Σ_lγ_l, equivalently the spectral norm of the centered
   diagonal residual. The F1-centered Frobenius diagnostic instead obeys
   ‖M_F1‖²_F = 4^(N−1)(16/9)Σ_lγ_l² with σ=Σ_lγ_l. Can either precisely
   named diagnostic be incorporated into design rules?

### Answered by the Absorption Theorem (April 4, 2026)

For uniform local Z-dephasing, the Absorption Theorem
(Re(λ) = −2γ⟨n_XY⟩) resolved three questions
that were not explicitly listed here but had been open since March:
what determines the spectral boundaries in the theorem's stated scope
(answer: light content), why the factor 2 where Σγ > 0, the F1 palindromizer
hypotheses hold and both endpoints are reached (answer: the full decay
range 0...2Σγ divided by its centre Σγ), and what sets the
spectral gap (answer: one absorption quantum, 2γ, above an N-dependent
coupling threshold; below it the gap is Zeno-suppressed and the theorem
supplies no lower bound). These are now formulas, not observations. See
[Absorption Theorem Proof](proofs/PROOF_ABSORPTION_THEOREM.md).

Under that uniform-rate specialization, the decay-rate identity is
α = 2γ⟨n_XY⟩. For a nonuniform rate profile the theorem uses the
site-weighted sum instead. Its linearity in the rates follows from the
dissipator; it is not a mass-energy relation.
See [Absorption Theorem Discovery](../experiments/ABSORPTION_THEOREM_DISCOVERY.md).

On IBM hardware (Q52 tomography): the Absorption Theorem ratio is 1.03
(3%). Detuning oscillations at 470 μs period are present. See
[IBM Absorption Theorem](../experiments/IBM_ABSORPTION_THEOREM.md).

### Require new tools or theory

These need mathematics or computational methods we have not built yet.

4. **The natural variable u(t):** Track u = C(Ψ+R) through Lindblad
   evolution. Does u have a simpler trajectory than CΨ or θ?

   *Partially addressed (April 2026):* along real Bell+ trajectories,
   no. u(t) is essentially linear in Ψ (u ≈ 0.61·Ψ^{1.02}), so it does
   not provide a simpler dynamical coordinate. u is an algebraic
   linearization variable, not a dynamical one. Whether u behaves
   differently after an explicitly defined complex extension is an open
   question; that extension is CΨ_com, not the original real CΨ. See
   [Boundary Navigation](../experiments/BOUNDARY_NAVIGATION.md).

5. **Cockpit scaling beyond N=11:** Does the first-three-PC cumulative
   coverage remain high at N=50-100? Does n95 grow linearly or saturate for
   dense topologies? The N=5-11 baseline is in [Cockpit Scaling](../experiments/COCKPIT_SCALING.md).

6. **Non-Markovian noise:** Does the cockpit framework hold under
   colored noise (noise whose strength depends on frequency, unlike white noise which is flat), 1/f spectra, or TLS coupling?

### Require experimental contact

These need access to quantum hardware we do not currently have.

7. **2-qubit Concurrence measurement:** Directly measuring Concurrence on
   a qubit pair would test a cockpit observable not yet covered by the
   hardware record.

8. **Cross-platform replication:** CΨ = ¼ crossing on trapped ions
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
| **Proven** (algebra + proof) | ¼ boundary, Mandelbrot equivalence, palindromic symmetry (Π operator), topology-independence, Pauli weight complementarity, Frobenius orthogonality of oscillation and cooling blocks at N=2 and its exact loss for N>2 under uniform local Z-dephasing and a nonzero shadow-balanced bond Hamiltonian, uniform-rate Absorption Theorem Re(λ) = −2γ⟨n_XY⟩ (site-weighted form for profiles) |
| **Verified** (simulation, reproducible) | 87,376 eigenvalues paired, QST star 2:1 beats chains, XOR space filter, ε→0 sacrifice-zone gains 68-360× versus the corresponding V-shape baselines on the stated peak-created-Sum-MI metric, cockpit first-3-PC coverage 91.9-99.0% over 8-9 standardized features in eight `center_bell` configurations for N=5-11 |
| **Hardware-confirmed** | CΨ crossings: Q52 qualitative (10.7% timing offset), Q80 within 1.9%; selective-DD and sacrifice-zone contrasts on IBM devices; T2* drift 58% in 6 days |
| **Argued** (plausible, not proven) | Measurement = crossing ¼, Mandelbrot boundary as route catalog |
| **Unverified** (could not reproduce) | t_coh ~ N scaling |
| **Refuted** | 33:1 coherence ratio |
| **Philosophical** | C = consciousness, 4D block-universe interpretation |

---

## The honest summary

The palindromic symmetry is proven and verified at scale. The CΨ = ¼
boundary is algebraically exact; two IBM crossing runs provide hardware
contact of unequal precision. The sacrifice-zone and cockpit results are
reproducible within their stated observables and scopes.

The weaknesses that remain are either philosophical (#1 generic
bifurcation, #2 consciousness label), require new experiments (#3 2-qubit
tomography), or are specific technical gaps (#4 variable u, #5 sacrifice-zone
optimality). None threaten the core mathematics.

The biggest direct measurement gap is two-qubit Concurrence tomography.
Current large-N cockpit scaling instead identifies Purity as the dominant
PC1 proxy, so Concurrence should be tested as a missing observable rather
than advertised as the dominant instrument.
