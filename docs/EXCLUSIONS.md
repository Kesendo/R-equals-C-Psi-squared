# What the Mathematics Excludes

<!-- Keywords: d squared minus 2d equals zero qubit necessity, CΨ quarter boundary
structurally stable, noise external origin incompleteness proof, palindromic symmetry
algebraically exact, irreversible quantum classical transition, fold catastrophe
topological invariant, time origin external Lindblad, R=CPsi2 logical exclusions -->

**Status:** Derived from proven results; each exclusion cites its proof
and carries its own grade (Tier-1 algebra in Exclusions 1, 4, 5, 6;
computational and conditional in 2; state-class-scoped in 3;
structural in 7; the closing γ-chain is graded link by link where it
stands).
**Date:** March 25, 2026, last refreshed 2026-07-16 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## Purpose

This document lists what the proven mathematics of this project
**rules out**. Not what we believe. Not what we interpret. What is
logically excluded by the proofs and verified computations.

Each exclusion is numbered, states the mathematical basis, and names
what is ruled out. At the end: what remains after the exclusions,
and what we cannot exclude.

---

## Exclusion 1: The structure requires d = 2

**Basis:** Under single-axis dephasing (for general d: the
full-Cartan, diagonal case; single non-Cartan generators split
differently), each qudit has d² Pauli-like
operators: d are immune (diagonal), d²-d decay. The palindromic
pairing requires equal counts: d = d²-d, giving d²-2d = 0. The only
non-trivial solution is d = 2 (the qubit).

**Verification:** Qutrits (d=3, the three-level analogue of qubits) tested against 236 qutrit
dissipators (single Gell-Mann jump operators, all pairs, and 200 random linear
combinations, all dephasing-type; amplitude damping and other
non-dephasing channels lie outside the sweep, the source's named open
avenue; the sweep lives in
[Qubit Necessity](QUBIT_NECESSITY.md) §8a-§9): 0/236 palindromic.
The equation d² - 2d = 9 - 6 = 3, not zero.

**The law behind the null (F121):** for d > 2 the spectrum is
not random but partially palindromic, with the closed-form dissipator ceiling
paired(d, N) = Σₖ d^N · C(N, k) · (d−1)^min(k, N−k), equal to d^{2N} **iff d = 2**; the
operator realization caps any per-site mirror at (2d)^N paired coherences, full again iff
d² − 2d = 0, the trunk equation's third independent appearance. See
[F121 in Analytical Formulas](ANALYTICAL_FORMULAS.md) (the ceiling, and
the product cap in its operator-realization paragraph) and the live witness
(`inspect --root qudit` in the Object Manager CLI), which recomputes the ceiling at
inspect time.

**Operational handle (2026-05-01):** The d=0 substrate is now directly readable via `fw.stationary_modes(chain)` (the kernel of the Liouvillian L), `fw.d_zero_decomposition(rho, chain)` (substrate/decohering split of any state), and `fw.sector_populations(rho_or_psi, N)` (the natural d=0 observables p_n = Tr(P_n · ρ), measurable on hardware via Z-basis tomography alone). For uniform XY/Heisenberg + Z-dephasing the kernel has dimension N+1, spanned by the F4 sector projectors P_n in the {I, Z}^N Pauli sublattice (see [F4 in Analytical Formulas](ANALYTICAL_FORMULAS.md)).

**Proof:** [Qubit Necessity](QUBIT_NECESSITY.md),
[The Mirror That Looked Non-Local](../hypotheses/THE_BOOT_SCRIPT.md) Section 5
(the host document's non-local headline is superseded; the §5 counting
cited here is preserved).

**Ruled out:**
- d = 1 (a single-state system, nothing to superpose; a classical bit
  is the diagonal d = 2 case): d² - 2d = -1, no palindrome possible
- d = 3, 4, ... (higher dimensions): d² - 2d > 0, split is unbalanced
- d = infinity (continuous variables): no finite balanced split
- "Qubits emerged from something simpler": the palindromic structure
  does not degrade gracefully, it requires d = 2 exactly

---

## Exclusion 2: The time arrow cannot originate from within

**Basis:** Five internal candidates for the origin of dephasing noise (γ)
have been examined; four are eliminated, and the first is, in the
proof's own grading, "a structural constraint, not yet an elimination"
(the elimination of internal self-generation rests on Candidates 2-3):

| Candidate | Test | Result |
|-----------|------|--------|
| Internal bootstrap | [bootstrap_test.py](../simulations/bootstrap_test.py) | Cross-sector coupling = 0 exactly (structural constraint) |
| Qubit decay | [failed_third.py](../simulations/failed_third.py) | Non-Markovian; 0/16 palindromic (wrong structure) |
| Qubit bath | Requires its own γ | Circular (shifts the problem) |
| Nothing (d=0) | Has no properties | Cannot generate anything |
| Other dimensions (d>2) | d² - 2d > 0 (the balance test fails) | Only d=2 exists in the framework |

**Proof:** [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md),
Section 2, Candidates 1-5. The proof is conditional on the d ∈ {0, 2}
ontology (its own wording: IF d(d−2) = 0 is the correct description,
THEN the source cannot be internal; its Open Directions concede that a
more complete theory could contain its own source).

**Corollary:** γ provides the time arrow (without γ: unitary oscillation,
no irreversibility, no before/after); the formal parameter t exists
trivially without γ, which is why this exclusion is about the arrow,
not about t. If γ cannot originate internally, then the time arrow
cannot originate internally. At Σγ = 0: Π·L·Π⁻¹ = -L
(exact time reversal, no fold, no crossing, no irreversibility).
The fold at 1/4 emerges only above Σγ_crit/J ≈ 0.25-0.50%
(state-dependent, but N-independent; these numbers live in
[Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)).
See also [Π as Time Reversal](../experiments/PI_AS_TIME_REVERSAL.md)
and [γ-Time Distinction](GAMMA_TIME_DISTINCTION.md).

**Ruled out (within this framework):**
- "The system created its own γ": the internal candidates eliminated
  (four outright, the bootstrap as a structural constraint)
- "γ emerged from an internal event": no internal mechanism produces
  the Lindblad dephasing parameter (this is a statement about γ in
  the Lindblad equation, not about cosmological time in general)
- "The time arrow from dCΨ/dt < 0 is an approximation": it is strict
  under Markovian dynamics within Exclusion 3's proven scope (the
  Markov assumption itself is a separate question)

---

## Exclusion 3: The quantum past is irrecoverable

**Basis:** dCΨ/dt < 0 strictly for all t > 0, proven in closed form
for the Bell+ state channel by channel (the five channels in the
Proof line below), and as
the Envelope Theorem (the local maxima of CΨ form a strictly
non-increasing sequence) for ANY 2-qubit state under local
Z-dephasing. Within that proven scope, CΨ (for Bell+) and the CΨ
envelope (for general 2-qubit states; individual oscillations may
briefly cross) cross ¼ exactly once, downward, and never return under
Markovian dynamics. The N ≥ 3 full-state envelope is open, and it
genuinely rises at N ≥ 4 under strong coupling (see the proof's scope
section and
[Envelope Rise Boundary](../experiments/ENVELOPE_RISE_BOUNDARY.md)).
At Σγ = 0 (no noise): CΨ oscillates and never crosses 1/4. The
exclusion requires noise.
See [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md).

**Proof:** [CΨ Monotonicity](proofs/PROOF_MONOTONICITY_CPSI.md).
Proven for Z dephasing, X and Y noise (bit-flip, bit-phase-flip),
depolarizing, and amplitude damping.

**What survives a crossing:** Classical correlations (mutual information).
What does not survive: quantum coherence (superposition, the content CΨ
tracks) and entanglement (tracked by concurrence; CΨ = purity ×
normalized L₁-coherence is a coherence measure, not an entanglement
measure, see the definition in Exclusion 6).

**Ruled out (within the proven scope above):**
- "We can reconstruct the initial quantum state from present
  measurements": the quantum information has been irreversibly
  converted to classical correlations at each ¼ crossing
- "Quantum information is preserved in the reduced system": under
  Markovian dynamics of the reduced system, it is converted to
  classical correlation, not stored (the full system+environment
  evolves unitarily, but the reduced system does not)
- "The past is still quantum": everything that crossed ¼ is
  classically decided. The doors are closed.

**Caveat:** Non-Markovian dynamics (structured bath) can temporarily
push CΨ back above 1/4 (max revival: CΨ = 0.3035, the structured-bath
table in [CΨ Monotonicity](proofs/PROOF_MONOTONICITY_CPSI.md) at
J_SB = 5.0, γ_B = 0.50). But revivals are
always transient. CΨ goes to zero eventually. The past is delayed,
not reversed. See [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md)
for the oscillation data (81 crossings with Bell+bath at J=5).

**The general principle: information is not stored, it is converted.**

Exclusion 3 describes one level of this principle: quantum coherence
is irreversibly converted to classical correlation at the 1/4 crossing.
The same principle operates at every level of the framework:

| Level | What is converted | Into what | Evidence |
|-------|------------------|-----------|----------|
| 1/4 crossing | Quantum coherence | Classical correlation | CΨ monotonicity (this exclusion) |
| V-Effect (experiment) | Individual frequencies | New coupled frequencies | 100% NEW-NEW: no N=2 frequency survives in N=5 ([pairing_structure_n5.txt](../simulations/results/pairing_structure_n5.txt)) |
| Energy partition (hypothesis) | Unstructured (unpaired) modes | Entropy | 2x decay law: unpaired modes die twice as fast ([Energy Partition](../hypotheses/ENERGY_PARTITION.md)) |

At each level:
- The original information is irrecoverable (not stored, not preserved)
- The conversion is irreversible (no mechanism to undo it)
- The result has more structure than the input (109 frequencies from 4,
  palindromic structure survives while noise dies 2x faster)

Two N=2 resonators (each 2 frequencies, Q=1, no oscillation) couple
through a mediator and produce an N=5 system with 109 frequencies and
Q=19+ (at least 19, still rising with J). Not one of the original 4
frequencies survives. They are
converted into something qualitatively new. Like the quantum coherence
that becomes classical correlation at the fold: the old form is gone,
but what replaces it is richer, not poorer.
([V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md),
[Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md))

---

## Exclusion 4: The structure is not from randomness

**Basis:** The palindromic pairing is algebraically exact.
87,376 eigenvalues tested (N=2 through N=8). Zero exceptions.
Pairing error at machine precision (the identity itself being
algebraically exact). Holds across all topologies (chain, star,
ring, complete, binary tree, the last verified to N=5) and all standard coupling models
(Heisenberg, XY, Ising, XXZ, Dzyaloshinskii-Moriya, the last under the
site-alternating Π variant).
See [Non-Heisenberg Palindrome](../experiments/NON_HEISENBERG_PALINDROME.md)
for the coupling models; the topology sweep is in the
[Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md).

**Proof:** [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md).
Three-step algebraic identity: Π flips XY-weight, Π anti-commutes
with the Hamiltonian's action [H, ·] (that is, Π L_H Π⁻¹ = −L_H on
operator space; not {Π, H} = 0), combined: Π L Π⁻¹ = -L - 2Σγ I.

**The operator mechanism below the spectrum:** the palindromizer
factors as Π_Z = R·D with
⟨R, D⟩ ≅ D₄, eight signed permutations of the Pauli basis compared
exactly ([Π factors as R·D](proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md)),
and on the Δ = 0 block lattice the mirror is an exact entry-wise
rearrangement: the fold-lattice group of eight, every leg
Jordan-preserving and holomorphic at the price λ → −λ − 2Nγ (the
"fold lattice, derived" section of
[Codim-1 by Additivity](proofs/PROOF_CODIM1_BY_ADDITIVITY.md); its
separate antiunitary cross-fold carries −λ̄ instead, the conjugations
cancelling in the lattice composition), residuals exactly 0.0 with no
eigensolver involved.
A random process would now have to reproduce not merely a paired
spectrum but an entry-wise rearrangement identity.

**Hardware:** the single-qubit consequence of the palindrome (the
CΨ = 1/4 crossing equation) confirmed on IBM Torino at 1.9% deviation
([IBM Run 3](../experiments/IBM_RUN3_PALINDROME.md)); the N-qubit
pairing itself is algebraic + numerical, not yet a hardware readout.

**Cross-domain (tested 2026-05-01, hardened 2026-06-28):** the F1 pair-sum-constant signature was tested on the periodic table of the elements, across three property scales and five periods (first ionization energies, Pauling and Allen electronegativities; the sweep and its numbers live in [`simulations/periodic_palindrome.py`](../simulations/periodic_palindrome.py) and §3 of [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md)). The current verdict is the hardened one: the pair-sum-constant statistic is satisfied exactly by any linear ramp, the original shuffle-null significance was largely a monotonic-smoothness artifact, and the sign-flip-null re-analysis finds a residual mirror-respecting signal that cannot be pinned to F1 specifically, neither validating nor refuting it ([the hardened re-analysis](carbon/PERIODIC_PALINDROME_HARDENED.md)). The cross-domain transport from quantum F1 to atomic shell Hamiltonians remains empirical, not derived.

**Ruled out:**
- "The symmetry is approximate and breaks at some scale": it is
  algebraically exact and has been verified up to N=8 (65,536×65,536
  matrices)
- "The symmetry emerged from random fluctuations": random processes
  do not produce algebraically exact symmetries with zero exceptions
- "This symmetry is statistical": the palindromic pairing is
  algebraic, not the result of averaging or large numbers

---

## Exclusion 5: The ¼ boundary is not removable

**Basis:** The recursion R = C(Ψ+R)² has discriminant D = 1 - 4CΨ.
D = 0 at CΨ = ¼. This is the normal form of the fold catastrophe
(the simplest in Thom and Arnold's classification of qualitative transitions that cannot be removed by small perturbations), which is structurally stable: small
perturbations cannot remove it, split it, or move it qualitatively.

**Caveat:** The fold exists only when Σγ > Σγ_crit (≈ 0.25-0.50% of J,
state-dependent but N-independent).
At Σγ = 0: no fold, CΨ oscillates without crossing 1/4. The fold is
not removable by perturbation, but it does not EXIST without noise.
Noise creates the fold. See [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md).

The value ¼ = (½)² is the degree-2 signature of purity (Tr(ρ²)); the
formal site of ¼ is the discriminant zero above, and the uniqueness of
the threshold is the Uniqueness Proof's α = 2 result. It is not a
fitted parameter.
The same PRODUCT STRUCTURE appears in classical neural networks:
σ(1-σ) reaches its maximum 1/4 at the sigmoid midpoint (decided x
undecided = ½ x ½); wherever two complementary halves form a product,
the maximum is (½)² = ¼. The cited analysis is explicit that only the
structure transfers: CΨ = 1/4 does not appear as a boundary in the
neural case. See [Neural 1/4](neural/ALGEBRAIC_PALINDROME_NEURAL.md).

**Proof:** [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md),
[Proof Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md).

**Ruled out:**
- "The ¼ boundary is an artifact of our model": it is a topological
  invariant of any quadratic recursion with the purity structure
- "At higher energies or larger scales, the boundary shifts": the
  fold catastrophe cannot be deformed away, only the approach to
  it can change
- "There are other boundaries at other values": α = 2 (purity) is
  the unique Rényi order with a state-independent threshold
  ([Uniqueness Proof](proofs/UNIQUENESS_PROOF.md))
- "The fold point connects classical and quantum regimes": the fold
  point is the algebraic site where the recursion's solution-multiplicity
  changes (0 real fixed points → 1 degenerate → 2 real, going
  downward through CΨ = 1/4). All measured variables (CΨ trajectory,
  sector populations, off-diagonal coherence, per-site Bloch
  components) cross the fold point smoothly. The fold point does not
  connect two worlds; it is the internal site where stability becomes
  algebraically available. (The repo's informal name for this site,
  "the cusp", survives in the script filename below; in the
  Thom-Arnold classification the site is a fold, A2, not a cusp, A3.)
  See [`simulations/cusp_investigation.py`](../simulations/cusp_investigation.py)
  for the full reading on Bell+ under Z-dephasing

---

## Exclusion 6: DD (Dynamical Decoupling) cannot restore CΨ

**Basis:** CΨ = Tr(ρ^2) x L₁(ρ)/(d-1), where d = 2^N is the full
Hilbert-space dimension (not Exclusion 1's single-site d). Purity Tr(ρ^2) is invariant
under all unitaries. L₁ coherence is invariant under the Pauli group
(Paulis permute computational basis states with phases, absolute value
absorbs phases). DD uses Pauli gates. Therefore the DD pulses leave CΨ
exactly invariant, and the free evolution between pulses only
decreases it (Exclusion 3), so the sequence as a whole cannot restore
CΨ. Delta = 0.00e+00 for all 16 two-qubit Pauli strings
(I,X,Y,Z)⊗(I,X,Y,Z).

**Proof:** [CΨ Monotonicity, Part 7](proofs/PROOF_MONOTONICITY_CPSI.md).

**Additionally:** No local unitary (Pauli or non-Pauli) can push CΨ
back above ¼ once crossed. The only mechanism for CΨ revival: non-
Markovian backflow through J-coupling to a coherent reservoir.

**A derived sibling (2026-07-16, a different symmetry, the same
lesson):** the counter-intuition about symmetry protection generalizes
beyond Π. The chiral mirror Θ = T·K of the F129 fringe protects a
Θ-mirror mode pair at first order (zero first-order drift
from ZZ crosstalk and detuning) and anti-protects it at second order,
pushing the pair apart by exactly twice one branch's shift: one
symmetry read at odd and at even order. Derived in
[the ζ² anti-protection law](proofs/PROOF_ZETA2_ANTI_PROTECTION.md);
on hardware it entered the flown F129 experiment as the budgeted
second-order drift, sign-consistent within the blended budget
(the registry's F129 Kingston-fringe entry, Confirmation 24, confirms
[the standing fringe](../experiments/IBM_F129_RAMSEY_FRINGE.md),
not the anti-protection sign separately). Protection is a parity
statement, not a shield.

**Ruled out:**
- "DD can refocus coherence and restore CΨ": DD uses Pauli gates,
  CΨ is Pauli-invariant, delta is exactly zero
- "Phase refocusing restores the quantum state": for the CΨ metric
  specifically, no Pauli-gate sequence has any effect
- "Active error correction can push CΨ back above ¼": only energy
  exchange (J-coupling) works, not phase refocusing (DD)

---

## Exclusion 7: The classical/quantum distinction is not a world-separation

**Basis:** d=0 (sector populations, kernel of L) and d=2 (off-diagonal coherences, decohering content) are not two ontologies but two indices of one ρ on the operator space of dimension (2^N)² = 4^N (this 2^N is the full Hilbert dimension, not Exclusion 1's single-site d). R = CΨ² and d² − 2d = 0 are the same polynomial at two parameter regimes: setting Ψ=0 and C=1/2 in R = C(Ψ + R)² collapses to R(R−2) = 0, exactly the dimension equation. The qubit on IBM hardware is one piece of aluminum with simultaneous classical (geometry, temperature, position in a fridge) and quantum (superposition, entanglement when the operations ask) readings; there is no separated quantum sphere from which it descends.

**Verification:** The algebraic identity Ψ(ρ_d0) = 0 for the kernel projection is exact for uniform XY/Heisenberg + Z-dephasing, verified numerically by `test_psi_vanishes_on_d_zero_substrate` in [`simulations/framework/tests/diagnostics/test_d_zero.py`](../simulations/framework/tests/diagnostics/test_d_zero.py). Any state's coherence (Ψ) vanishes on its d=0 projection by construction; both d=0 and d=2 are diagonal-vs-off-diagonal entries of the same Hilbert-space density matrix.

**Proof:** Structural identity documented in [On What the Formula Already Knew](../reflections/ON_WHAT_THE_FORMULA_KNEW.md) (R=CΨ² and d²−2d=0 as one quadratic family). The motto "Reality is what happens between us / mirrors / sectors" reads at the operator level: sectors P_n are static, the d=2 off-diagonal content between them is dynamic. The canonical poetic form is in [reflections/TRANSMISSION](../reflections/TRANSMISSION.md) ("We are all mirrors. Reality is what happens between us").

**Ruled out (as readings of this framework's mathematics; this
exclusion is structural, about what the framework's own objects are,
not a metaphysical theorem):**
- "There is a quantum world somewhere separate from the classical world we inhabit": the framework's structures are operator-space distinctions, not world-separations
- "We are external observers of a quantum system from a classical world": we are inside the same operator space whose dynamics we describe; the framework is self-description from within
- "Classical physics emerges from quantum physics through some metaphysical transition": both readings exist on the same substrate at all times; the only "transition" is the algebraic phase change at the fold point (Exclusion 5)
- "d=0 is the classical part and d=2 is the quantum part": both are entries of the same density matrix on the same Hilbert space; the labels mark *reading mode*, not ontology

---

## What Remains After the Exclusions

The exclusions leave a narrow picture:

1. The fundamental unit is the qubit (d=2). Not by choice. By
   algebraic necessity.

2. The palindromic symmetry was not created. It exists whenever
   d = 2 and local dephasing is present. It is a mathematical
   consequence, not a physical event.

3. The time arrow (γ) comes from outside the framework. The system
   cannot explain its own irreversibility.

4. Everything that has happened since: CΨ crossings of ¼.
   Quantum possibility becoming classical fact. Irreversibly.
   Each crossing deposits a piece of decided reality.

5. What we measure classically is consistent with the irreversible
   result of quantum-to-classical conversion at the fold catastrophe.
   Each crossing deposits classical correlations that persist.

---

## What the exclusions imply together

Individually, each exclusion removes one possibility. Together, they
form a chain that constrains what γ can be, each link graded below:

1. γ must be external (Exclusion 2: four internal candidates
   eliminated, the first, the internal bootstrap, reduced to a
   structural constraint)
2. The qubit chain behaves as a passive optical cavity (four of five
   standard cavity tests satisfied quantitatively, R² = 0.998 for the
   beam profile, Gouy marginal at R² = 0.81; the failing fifth is the
   even = confocal identification, which fails at the small N tested;
   [Optical Cavity Analysis](../experiments/OPTICAL_CAVITY_ANALYSIS.md))
3. Therefore γ is whatever enters the cavity from outside
4. In optics, what enters a passive cavity from outside is light
5. On IBM transmon hardware, the physical carrier behind γ IS light:
   microwave photons in the qubit's readout resonator cause the
   dephasing through photon shot noise (Sears et al., Phys. Rev. B 86,
   180504, 2012); γ itself is the rate this light imprints

Step 1 is proven. Step 2 is verified quantitatively on four of five
checks (the fifth, the even = confocal identification, fails at the
small N tested). Step 3 is deductive within the cavity reading. Step 4 is the
physical identification, the analogy's content rather than a
derivation. Step 5 is independent experimental confirmation on a
different platform than our theoretical framework.

This does not prove "γ is light" as a universal statement. It proves
that the mathematics forces γ to be external; it shows that the system
γ enters matches an optical cavity quantitatively on four of five
checks; and on the one hardware platform where we can check, the
external input is literally photons. On other platforms (ion traps, neutral atoms), γ has a
different physical origin, but the algebraic role (external input
to a cavity) is the same.

See [What If Gamma Is Light?](../hypotheses/GAMMA_IS_LIGHT.md) for the
full discussion, including speculative extensions (Tier 4) that go
beyond what the exclusions support.

---

## What We Cannot Exclude

- What the external source of γ is *in general* (on IBM hardware it is
  photon shot noise; on other platforms it may be different)
- Whether the palindromic mirror (Π side) has any form of experience
- Whether our framework is embedded in a larger structure
- Whether non-Markovian revivals have deeper significance: they ARE
  the mechanism by which coupled resonators sustain oscillation.
  Each subsystem is the other's coherent reservoir. Backflow through
  a mediator is what produces Q > 1
  ([Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md))
- Whether the V-Effect (two 2-frequency systems producing 109
  frequencies through coupling) is the mechanism of biological
  evolution and complexity growth
  ([Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md))
- Whether the structured patterns in IBM hardware calibration data
  ([Q98 crossing pattern](BOTH_SIDES_VISIBLE.md)) are the "bridge" or
  mundane hardware drift
- The nature of what existed "before" γ began (if that question is
  even meaningful from inside, see [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md))

These are open questions at the boundary of the framework's reach.
They are not excluded. They are undecidable from within.

---

## References

- [Mirror Symmetry Proof](proofs/MIRROR_SYMMETRY_PROOF.md): the palindrome
- [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md): four eliminations plus a structural constraint
- [CΨ Monotonicity](proofs/PROOF_MONOTONICITY_CPSI.md): irreversibility
- [Qubit Necessity](QUBIT_NECESSITY.md): d² - 2d = 0
- [Uniqueness Proof](proofs/UNIQUENESS_PROOF.md): ¼ is unique
- [Mathematical Connections](MATHEMATICAL_CONNECTIONS.md): fold catastrophe
- [Proof Roadmap](proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md): seven layers
- [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md): γ as mediator
