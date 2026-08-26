# The Universal Palindrome Condition

<!-- Keywords: universal palindrome QXQ+X+2S antisymmetry, selective damping
two populations swap operator, Dale's Law Pauli algebra generalization,
hierarchy incompleteness palindrome mechanism, V-Effect coupled dead systems,
R=CPsi2 universal palindrome -->

**Status:** Hypothesis (Tier 4): proven in one domain, translated into a second
where no system satisfying it has yet been found
**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## What this document is about

The palindromic symmetry was discovered in quantum systems. Then it
appeared in neural networks. Two completely different domains, built
from different physics, studied by different communities, with
different mathematics. And yet both satisfy the same equation.

This document asks: is that a coincidence, or is there one underlying
rule? The answer proposed here is a single algebraic condition with
three ingredients: two types (fast and slow, quantum and classical,
excitatory and inhibitory), a swap that exchanges them, and coupling
that flips sign under the swap. Whenever all three are present, the
palindrome follows. It does not matter whether the system is made
of qubits, neurons, or protons in a hydrogen bond.

If this is right, the hierarchy of reality (atoms → molecules →
cells → brains) is not designed. It is algebraically forced: each
level builds palindromic units, couples them (breaking the local
symmetry), and the new frequencies that emerge ARE the next level.

---

## The Claim

Any dynamical system with three properties has palindromic spectral
symmetry:

1. **Two populations whose decay rates pair to one and the same sum**, that
   is, 1/τ_{Q(i)} + 1/τ_i identical at every seat, since the pairing constant
   is a single scalar. This is the condition; "different decay rates" is
   not. A Q that exchanges the two populations
   satisfies it at ANY rates, equal or not; at equal rates every permutation
   satisfies it and it imposes nothing. What different rates do is force Q to
   exchange the populations rather than preserve them.
2. **A swap operator Q (involution: applying it twice gives back the
   original, Q² = I) that exchanges them.**
3. **Coupling that is antisymmetric under Q.**

When all three hold, the evolution operator X satisfies:

```
Q · X · Q⁻¹ + X + 2S = 0
```

where Q⁻¹ is the inverse of Q. (In neural networks, Q is a permutation
with Q² = I, so Q⁻¹ = Q. In quantum, Π has complex phases and Π² ≠ I,
so Π⁻¹ ≠ Π. The Q⁻¹ form covers both.)

Every eigenvalue μ has a partner μ' with μ + μ' = const. The
decay rates mirror around a center. The eigenmodes swap character
between the two populations. Coupling two such systems changes the number
of oscillatory modes; that the change is caused by breaking the local
symmetry holds in the quantum case and is withdrawn in the neural one
([Proof: V-Effect Mechanism](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)).

This is one equation. It is proven at the smallest scale (qubits) and it can
be WRITTEN at the scale of neural networks, where it holds on networks built to
satisfy it and fails on the one connectome tested. Those are not two
confirmations; they are a theorem and a translation whose biological instance is
still missing. We conjecture the equation may hold at scales in between, but
this is untested: every intermediate candidate examined below
(antiferromagnets, atoms) so far fails to meet all three conditions.

---

## The Evidence

One domain proves the condition; the other supplies its translation and, so
far, no system that satisfies it.

| Component | Quantum | Neural |
|-----------|---------|--------|
| Evolution operator X | L (Liouvillian: the matrix that governs quantum decay) | J (Jacobian: the matrix that governs how neurons influence each other) |
| Two populations | Immune {I,Z} vs Decaying {X,Y} | Excitatory vs Inhibitory |
| Split | 2:2 (C = 0.5) | N/2 : N/2 (balanced) |
| Swap Q | Π (Pauli weight swap, Π² ≠ I) | Q (E-I permutation, Q² = I) |
| Selective damping | γ (dephasing) | 1/τ_E ≠ 1/τ_I, which forces Q to swap types rather than enabling the pairing (corrected 2026-08-26) |
| Antisymmetry source | Commutator [H, ρ] | Dale's Law |
| Eigenvalue pairing | λ + λ' = -2Σγ | μ + μ' = -(1/τ_E + 1/τ_I) |
| Character swap | Population ↔ Coherence (100%) | E-dominant ↔ I-dominant (96%), on SYNTHETIC balanced networks; never measured on a connectome |
| V-Effect | 2+2 = 109 frequencies | 0+0 = 48 frequencies, but the two zeros are a property of the seeds used, not of the symmetry (see below) |
| Status | **Algebraically proven** | **Algebra verified on constructed networks; the one connectome tested does not satisfy the condition** |

Sources: [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (quantum),
[Neural Palindrome Proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md) (neural),
[V-Effect Neural](../docs/neural/V_EFFECT_NEURAL.md) (0+0=48)

**Important difference:** In quantum, the commutator [H, ρ] = Hρ - ρH
provides antisymmetry of BOTH signs and magnitudes automatically. In
neural, Dale's Law provides only the signs. The magnitudes require an
additional condition (W[Q(i),Q(j)] = -(τ_{Q(i)}/τ_i) · W[i,j]).
The quantum palindrome is unconditional. The neural one is conditional.

In quantum: 87,376 eigenvalues verified (N=2..8), zero exceptions. That is
computation, not measurement: the multi-qubit palindrome has never been measured
on hardware. The 1.9% deviation on ibm_torino Q80 is a different quantity, the
single-qubit CΨ = 1/4 crossing (F24, IBM Run 3); N ≥ 2 is untested.

**Application to chemistry: the hydrogen bond as qubit.**
The proton in O-H...O tunnels between |L⟩ (donor) and |R⟩ (acceptor).
d = 2. This IS the quantum palindrome applied to a physical system.
Computed: palindrome exact, V-Effect creates 104 new frequencies when
two water molecules couple through a H-bond. Three regimes: classical
(J/γ << 1), fold (J/γ ~ 1, enzymes), quantum (J/γ >> 1). Bulk water was
assigned to the first from a lower bound; its supported band reaches at its top
the value quoted as quantum for the Zundel cation (see [Q Belongs to No Substance](../docs/Q_BELONGS_TO_NO_SUBSTANCE.md)).
See [Hydrogen Bond Qubit](../docs/water/HYDROGEN_BOND_QUBIT.md).

In neural: synthetic networks with exact condition give residual = 0.
There is no 8.46× against Erdős-Rényi: that measurement divided the two arms
by different constants,
and the worm arm's mean ‖W_eff‖ stands to the control's in a RATIO of 8.50 at
N = 10 (5.75 at N = 20, 4.82 at N = 26), so the ratio tracked coupling magnitude. Matched, it runs 0.960 at N = 10, 0.841 at N = 20 and
0.748 at N = 26, and what that smaller residue is has not been settled. The degree-preserving rewiring that scores an
identical 1.0 went with it: that null cannot move such a metric
([Algebraic Palindrome Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md)).

---

## Where to Look Next

### Confirmed: hydrogen bond (quantum application)

The proton in a hydrogen bond is a qubit (d=2). The palindrome is
proven. 104 new frequencies via V-Effect when two molecules couple.
This is not a new domain; it is the quantum palindrome applied to
a physical system. It bridges quantum (Level 0) and chemistry.
See [Hydrogen Bond Qubit](../docs/water/HYDROGEN_BOND_QUBIT.md).

In the Zundel configuration (proton centered between two oxygens,
H₅O₂⁺), J/γ = 4.8 at 300K. CΨ crosses 1/4 six times in 21
femtoseconds (~6 crossings per proton transfer event). Water is
not near the fold. Water IS the fold.

A classical model of the same system (donor/acceptor modes as coupled
oscillators) shows NO palindrome (residual 1.33). The palindrome lives
at the quantum level of the proton, not at the classical level of the
bond. This teaches: look for the QUBIT inside the system, not for
classical analogs.

### Further candidates

#### Antiferromagnets (difficulties identified)

An antiferromagnetic crystal has two sublattices (spin-up, spin-down).
A natural swap Q exists (sublattice exchange). But:

- **Condition 1 (selective damping): unclear.** Both sublattices are
  the same material with the same intrinsic decay rates. In quantum,
  selective damping comes from operator structure ({I,Z} immune,
  {X,Y} decaying). In neural, τ_E ≠ τ_I plays a different
  role than the quantum γ: it does not create the asymmetry, it constrains
  which swaps are admissible (corrected 2026-08-26). In antiferromagnets,
  the analog is not obvious. Magnon modes on different sublattices
  might couple differently to phonon baths due to staggered order,
  but this is not guaranteed.

- **Condition 3 (antisymmetric coupling): problematic.** The Heisenberg
  exchange H = J Σ S_i · S_j is SYMMETRIC under sublattice swap
  (Q·H·Q = H, not -H). The antisymmetric quantity is the staggered
  magnetization (order parameter), not the Hamiltonian. This is a
  different kind of antisymmetry than the commutator (quantum) or
  Dale's Law (neural).

**Status:** The three conditions may not be met as straightforwardly
as initially assumed. The test requires identifying the correct Q and
the correct sense of "antisymmetry" for magnetic systems.

**Answered 2026-06-01** (verified in [`simulations/palindrome_general.py`](../simulations/palindrome_general.py)
+ this session's antiferromagnet probe): the Condition-3 diagnosis above was at the wrong level, but
the intuition was right. The Heisenberg EXCHANGE under Z-dephasing IS palindromic about −Σγ (h = 0,
bit-exact N = 2..4) by F1, which does not use the sublattice swap at all; the F1 conjugation makes
−i[H,·] mirror-odd for the coupling regardless of any sublattice symmetry, so looking for QHQ = −H at
the Hamiltonian level was the wrong test. But the AFM ORDER, a longitudinal field h·Σ(±1)^k Z_k (the
staggered magnetisation this section already flagged as "the antisymmetric quantity"), BREAKS the
palindrome: both the full-eigenvalue palindrome and the rate-only (real-part) palindrome fail for any
h ≠ 0, uniform or staggered (verified). So the palindrome lives in the DISORDERED (paramagnetic)
phase and dies the moment the order sets in. **Order is the death of the mirror**: the magnetic
reading of the project's recurring theme (the full shell is dead, the V-Effect needs broken symmetry,
the incomplete carries the structure). Condition 1 (whether the physical AFM dissipation is
Z-dephasing-like) stays the genuine open part; Condition 3 is resolved (the exchange is fine, the
order is the breaker).

### Atoms in external fields (candidate, one difficulty)

An atom in a magnetic field has Kramers-degenerate pairs (states that
are related by time-reversal: if one spins clockwise, its partner spins
counterclockwise, and without a magnetic field they have identical
energy). The Zeeman term flips sign under time reversal (condition 3).
A natural Q exists (Kramers conjugation, condition 2).

- **Condition 1 (selective damping): unclear, and it does not say what it
  used to.** Kramers partners typically have identical decay rates (same
  selection rules). The sentence that used to stand here, that without
  selective damping the palindrome has no centre to mirror around, is wrong:
  the centre is trace over dimension and is defined whether or not the
  spectrum pairs (F137). What the damping half of the neural theorem actually
  asks, transposed to rates r_i, is that r_{Q(i)} + r_i be the SAME number at
  every seat, since the pairing constant is one scalar. Three regimes follow,
  and the middle one is this candidate's. Rates uniform throughout: every Q
  passes and the condition imposes nothing. Rates equal WITHIN a doublet but
  differing across doublets: the natural Kramers swap, which stays inside a
  doublet, FAILS, since the sums it produces differ from doublet to doublet.
  Whether ANY admissible Q exists is then a matching question on the rates, and
  with exactly two doublets a Q pairing the fast one with the slow one does it
  (checked on r = (0.7, 0.7, 0.2, 0.2): the within-doublet swap needs the sums
  1.4 and 0.4, two different numbers, while the cross-doublet ones need only 0.9
  everywhere). With three doublets at a, b, c it needs a + b = 2c or equal
  rates, so existence is not general.
  Rates all distinct: the condition becomes a matching problem on the rate
  multiset. So equal decay rates are not an obstruction in general, but for
  the conjugation this section proposes they can be exactly the obstruction.

**Status:** Conditions 2 and 3 are cleanly met. Condition 1 is unclear, and it
is the condition that would decide the natural Q here, since a within-doublet
swap needs every doublet to carry the same rate. What has not been checked for
atoms is the half that binds everywhere else: whether the full operator equation
holds, magnitudes and zero pattern included.

### The key question for all candidates

In quantum, the three conditions arise AUTOMATICALLY from the
mathematical structure (Pauli algebra, commutator). In neural, the damping
condition is satisfied by any type-swapping Q at any time constants, so it
costs nothing; Dale's Law comes free with neurotransmitter identity but
fixes only the SIGNS, and only where the weight is nonzero. What is left
conditional is the zero pattern and the magnitudes, and on the one
connectome tested it is the zero pattern that fails, on a count.

For other systems: which conditions are automatic and which require
fine-tuning? If all three require fine-tuning, the palindrome is a
coincidence. If at least two are structural, the hypothesis gains
support.

### The V-Effect as the sharpest test

If a candidate system shows the V-Effect (coupling two locally
palindromic subsystems creates new frequencies), that is stronger
evidence than the palindrome alone. The V-Effect requires exact
local symmetry AND its breaking through coupling. It is harder to
achieve by accident.

At both levels:
- Quantum: 2 frequencies each → 109 coupled
- Neural: 0 frequencies each → 48 coupled, for the two seeds used

The neural row no longer carries a mechanism. The palindrome does not
forbid oscillation, so the two zeros are not symmetry-protected silence:
at the same coupling, 24 of 200 exactly palindromic draws oscillate
([Proof: V-Effect Mechanism](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)).
Finding the V-Effect at a third level would still be evidence for
universality, but the neural level would first have to establish what its
frequency counts track.

---

## The Consequence

If the palindromic condition is universal, then the entire hierarchy
of reality follows from one algebraic identity applied recursively:

1. Build palindromic systems at level N (the "atoms" of that level)
2. Couple them (the "bonds")
3. The coupling breaks local palindromes and creates new frequencies
4. These new frequencies ARE level N+1
5. Repeat

This means the hierarchy is not designed. It is algebraically forced.
Given two populations, a swap, and antisymmetric coupling, palindromic
symmetry is inevitable. Given palindromic symmetry, the V-Effect is
inevitable. Given the V-Effect, the next level of complexity is
inevitable.

The quantum palindrome pairs 87,376 eigenvalues with zero exceptions.
Dale's Law provides the sign antisymmetry across 300 identified neurons
in C. elegans, and that is as far as it goes: the 8× is measured against
Erdős-Rényi and shrinks to a much smaller residue of open origin once both arms
are normalised the same way, and on the connectome itself the palindrome condition does not hold
at all, its condition (b) failing on a count of 253 non-empty excitatory rows
against 18 inhibitory ones
([Neural Gamma Cavity](../experiments/NEURAL_GAMMA_CAVITY.md), 2026-08-25).
The same equation, then, but the worm is not an instance of it. The gap between
a qubit and a neuron is vast, but the algebra that organizes both is
one line:

```
Q · X · Q⁻¹ + X + 2S = 0
```

---

## The Limitations

### What does NOT transfer between domains

| Feature | Quantum | Neural | Universal? |
|---------|---------|--------|------------|
| Palindromic pairing | Exact (proven) | Exact if magnitudes match | Yes |
| Character swap | 100% | 96% at moderate coupling, synthetic networks only | Yes |
| V-Effect | From oscillation (2+2=109) | 0+0=48 on the seeds used; the palindrome does NOT enforce the silence | No: the neural mechanism is withdrawn |
| 2× decay law | Exact (2.00, N=2..5) | Not exact (0.84-1.63) | No |
| CΨ = 1/4 threshold | Fold catastrophe of R = CΨ² | σ(1-σ) = 1/4 at sigmoid inflection | Yes: (0.5)² |
| V-Effect from heat | n_bar > 0 creates 2 new modes | Drive P creates no modes | No |
| Exponential state space | 2^N (tensor product) | N (linear) | No |

The 2× decay law and the heat-driven V-Effect require the d² Liouville
structure with exact Pauli algebra. Neural systems have no analog.

### The value 1/4

In quantum: CΨ = Purity x Coherence. At the fold: 1/2 x 1/2 = 1/4.

In neural: σ(θ)(1-σ(θ)) = Decided x Undecided. At the sigmoid
inflection point: 1/2 x 1/2 = 1/4. This holds for every sigmoid,
every parameter set, every network size. It is the MAXIMUM of the
neural sensitivity. Parameterindependent. Computed and verified.

Both are the product of two complementary halves. Both are (0.5)^2.
Both are upper bounds (quantum: maximum CΨ for real fixed points;
neural: maximum sensitivity of the sigmoid response). The structure
is identical. The value is identical. 1/4 = the axiom squared.

### How to falsify

Find a system with all three conditions (two populations, swap,
antisymmetric coupling) where the palindrome does NOT hold. This
would mean the conditions are necessary but not sufficient, and
additional structure specific to qubits or neurons is required.

---

## Open Questions

1. **Antiferromagnet test:** The most accessible intermediate level,
   but conditions 1 and 3 are problematic. The Heisenberg exchange
   is SYMMETRIC under sublattice swap (not antisymmetric), and both
   sublattices have the same intrinsic decay rate. The test requires
   identifying a different sense of antisymmetry for magnetic systems,
   or a dissipative environment that breaks sublattice symmetry.
   (See "Where to Look Next" section for details.)

2. **What is the correct Q for atoms?** Candidates: Kramers
   conjugation (time-reversal), parity, spin-flip. Condition 1
   (selective damping) is unclear for Kramers partners, which
   typically have identical decay rates, but it is not the deciding
   one; see the candidate section above.

3. **What does the condition cost a system dynamically?** Not stability:
   exact magnitude matching buys no such thing. The condition makes the spectrum symmetric under μ ↦ −μ − 2s, which
   forces only Re μ + Re μ′ = −2s: stability then holds exactly while the
   spectrum fits inside a fixed-width strip, and the strip does not grow
   with the coupling. Measured on 200 exactly palindromic draws at N = 10,
   24 oscillate at coupling 0.5 and 45 are outright unstable at coupling
   10 ([Proof: V-Effect Mechanism](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)).
   So the condition is not a stability guarantee and "exact symmetry is
   dead" was never measured, it was inferred from a single seed. The
   question of what a system pays for satisfying it is open again.

4. **The inheritance mechanism.** How does the qubit palindrome
   propagate through atoms, molecules, chemistry, biochemistry to
   neurons? Dale's Law is the inherited form of the commutator.
   Through what chain of physical mechanisms? Can the intermediate
   steps be identified?

5. **V-Effect scaling with N.** Neural: 0+0=6 (N=10, coupling=0.01),
   0+0=12 (N=10, coupling=0.05), 0+0=48 (N=20, coupling=0.01),
   0+0=62 (N=20, coupling=0.05). The 62 is at the resolution of its own
   frequency binning and should be recomputed before it is used.
   How does the number of V-Effect frequencies scale with N?

6. **Universal coupling window.** Quantum V-Effect peaks at
   J/γ ~ 2-5. Neural V-Effect peaks at coupling 0.01-0.05.
   Is there a dimensionless ratio that is the same in both?
   **Answered 2026-06-01 (the question asks for the wrong kind of thing):** there is likely no
   universal dimensionless ratio. The peak location is system-specific, for three structural
   reasons. The couplings have different denominators (the neural g is relative to the intra-network
   weights, the quantum J relative to the dissipation γ). The local palindrome's robustness differs
   (the neural exact palindrome dissolves at g ~ 0.5, the quantum one is robust). And the starting
   points differ (quantum 2+2, neural 0+0, already noted in the limitations table). What IS universal
   is not a number but the SHAPE: the marginal window, silent → peak → silent (verified on the neural
   V-Effect sweep this session: oscillating modes rise to a peak at g ~ 0.05-0.1, then fall to zero
   as the coupling dissolves the two-subsystem structure). It shares the "alive in the middle" shape
   with the neural thermal window (both are two-sided dynamical windows). Reading it as ONE form with
   everything else would over-unify, and on a closer look it does not hold: the formation map is a
   ONE-sided threshold (clustering rises 0 → 1 and does not fall), and the coherence ceiling is a
   STATIC bilinear, not a dynamical window. Two honest tiers separate here. The genuinely INHERITED,
   structural form is the bilinear p(1−p) (the currency, d=2-borne, the same across real-qubit
   substrates). The reading "life sits in the marginal, between two deaths", grouping the dynamical
   windows, is Tier 5, the work's motor, not a verified single structure. So: the universal coupling
   window is a shape, not a ratio, and that shape is a two-sided dynamical window, kin to the thermal
   window but not identical to the inherited bilinear.

7. **Can the condition Q·X·Q⁻¹ + X + 2S = 0 be derived from a
   single axiom set** rather than proven separately in each domain?

8. **What is (0.5)^2 at intermediate levels?** The sigmoid maximum
   σ(1-σ) = 1/4 is neural. The purity fold CΨ = 1/4 is quantum.
   Both give (0.5)^2. What is the (0.5)^2 of an atom? A crystal?

---

## Seen again 2026-06-01: the abstract leg, and the two faces of 1/4

Returning with sharper tools, two of the open questions above moved. Verification:
[`simulations/palindrome_general.py`](../simulations/palindrome_general.py).

**Open Question 7 (derive once, not per domain): answered for the implication.** The condition was
confirmed in two domains; the question was whether the palindrome follows from the ingredients alone,
once, rather than being re-proven per substrate. It does. The probe builds a FRESH generator that is
neither quantum nor neural: a coupling part A made mirror-odd (Q A Q⁻¹ = −A) and a bath diagonal B
paired to a constant (B + Q B Q⁻¹ = −2c·I). The two sub-conditions ALONE force Q X Q⁻¹ = −X − 2c·I,
hence the palindrome: if X v = λ v then X(Q⁻¹v) = (−2c − λ)(Q⁻¹v), so Q⁻¹v is the partner mode. Bit-exact,
no reference to qubits or neurons; the same script re-confirms the quantum Liouvillian (about −Σγ) and
the neural Jacobian (‖J + Q J Q + 2S·I‖ = 0). So the equation is not "the same in two domains by
analogy"; it is one substrate-free algebraic fact, the implication derived once. What stays
domain-specific is only WHY each system satisfies the ingredients (the commutator in quantum, Dale's
Law in neural, construction in the abstract case). The centre −c is bath-set because A is mirror-odd
and contributes nothing to it: the substrate-independent form of the Absorption Theorem (rate centre
dephasing-set) and the neural Takt identity (centre membrane-set), one fact, the trace.

**Open Question 8 (what is (0.5)² at other levels): there are TWO 1/4, degenerate at d=2.** The doc
reads 1/4 = (0.5)² = "the axiom squared", and that reading is correct and d-INDEPENDENT: the maximal
single coherence obeys |ρ_ab|² ≤ ρ_aa·ρ_bb ≤ 1/4 (AM-GM), so the equal superposition of any two states
gives exactly 1/4 at every dimension (verified d = 2..5: the value stays 0.2500). This is the
quadratic 1/4, the fold of x², the coherence ceiling, the currency. But a SECOND 1/4 hides under it:
the polarity reading 1/d² = (1/d)², which is 1/4 only at d=2 and diverges (1/9, 1/16, 1/25) at higher
d. At d=2 the two coincide because (1/2)² = 1/4 = the quadratic maxval, a degeneracy that makes them
indistinguishable from d=2 alone. So the answer to "what is (0.5)² at intermediate levels": the
quadratic 1/4 (fold, ceiling) stays 1/4 at any dimension because it is about the square, not the
dimension; the polarity 1/d² is the face that changes, but it diverges from the framework's 1/4 at
d ≠ 2, and the qubit necessity (d² − 2d = 0) keeps the framework itself at d=2. The currency is the
quadratic 1/4; the 1/d² is the polarity coincidence. See
[Quarter and Half in Carbon](../docs/carbon/QUARTER_HALF_IN_CARBON.md) for the currency anchors.

---

## References

- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): quantum palindrome
- [Neural Palindrome Proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md): neural derivation
- [Algebraic Palindrome Neural](../docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md): C. elegans, standing wave
- [V-Effect Neural](../docs/neural/V_EFFECT_NEURAL.md): 0+0=48, thermal window
- [V-Effect Palindrome](../experiments/V_EFFECT_PALINDROME.md): 2+2=109
- [Hierarchy of Incompleteness](../docs/HIERARCHY_OF_INCOMPLETENESS.md): the levels
- [Energy Partition](ENERGY_PARTITION.md): 2× law, thermal window
- [Exclusions](../docs/EXCLUSIONS.md): what the math rules out

---

*One equation, proven in one domain and written down in a second. Whether
anything in that second domain satisfies it is the open question.*
