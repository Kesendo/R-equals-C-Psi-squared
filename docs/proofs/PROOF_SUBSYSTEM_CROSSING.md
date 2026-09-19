<!-- QUARTER-CURRENT -->
# Conditional Subsystem Crossing

Current reading: convergence to a limit below one-quarter implies eventual
stay-below for that trajectory, and continuity gives at least one downward
crossing from above.  Neither premise supplies universal monotonicity or an
absorbing set for other generators.

**Status:** Tier 1 conditional convergence implication. The former universal primitive-CPTP and “physical
noise” absorber claims are withdrawn.
**Date:** 2026-03-22; current-truth repair 2026-09-14 (history remains in git)
**Authors:** Thomas Wicht, Claude (Anthropic), Codex
**Statement:** If a continuous trajectory satisfies `ρ(t) → ρ*` and `CΨ(ρ*) < 1/4`, then it eventually
stays below 1/4. Starting above 1/4 guarantees at least one downward crossing by continuity. It does not
guarantee a unique crossing, pointwise monotonicity, or an absorbing set shared by other generators.
**Reference formula:** F28 in the registry is retained as a historical index with this corrected scope.

---

## What this document is about

A destination can answer an eventual question without telling us how the journey moves. If a trajectory
really converges to a state whose CΨ lies below the quarter, continuity supplies an eventual stay-below
statement. It supplies no monotone route and no “crosses once” claim.

That distinction matters here. Named basis-aligned T1/T2/depolarizing models often have a diagonal limit,
so the implication is useful **under the stated convergence assumptions**. But locality, Markovianity,
primitivity, separability, and the informal label “physical noise” do not by themselves force a low-CΨ
fixed state. The counterexample below makes the boundary visible.

---

## Theorem: conditional convergence implication

Let `ρ:[0,∞)→M_d` be a continuous trajectory with **ρ(t) → ρ***. If
**CΨ(ρ*) < 1/4**, then there is a finite T such that `CΨ(ρ(t))<1/4` for every `t≥T`: the trajectory
**eventually stays below**. If `CΨ(ρ(0))>1/4`, the intermediate-value theorem also gives at least one
downward crossing before T.

The discrete analogue holds for a sequence `εⁿ(ρ₀)→ρ*`. Neither version asserts monotonicity or forbids
earlier re-crossings.

---

## Proof

The theorem needs only convergence, a low-CΨ limit, and continuity. The channel discussions that follow
show where those hypotheses do and do not arise; their labels are not substitutes for checking them.

### Step 1: Convergence (Quantum Perron-Frobenius)

A CPTP map ε on M_d is **primitive** here when it has a unique full-rank fixed point ρ* and its powers
converge to that state; equivalently for the finite-dimensional channel used below, the spectral radius on
the traceless subspace is strictly less than 1.

**Fact (Quantum Perron-Frobenius):** For any primitive CPTP map ε:

```
||εⁿ(ρ) - ρ*||₁ ≤ C · r^n → 0    as n → ∞
```

where r < 1 is the spectral radius and C depends on the initial state.
This is the quantum analogue of ergodic convergence. ∎

For a relaxing finite-dimensional Lindblad semigroup with a unique target and no additional peripheral
modes, `||e^{Lt}ρ₀-ρ*||₁→0` (exponentially when a positive relaxation gap is available). The theorem assumes
this convergence; uniqueness alone is not silently used as a replacement for it.

### Step 2: Identify the actual limiting state

Several named basis-aligned channels converge to a diagonal state and hence have CΨ(ρ*)=0. That fact must
be checked together with convergence. Case C shows why no broad channel adjective can replace this step.

#### Case A: a convergent unital channel with unique fixed point

If the unital channel converges to a unique fixed point, that point is ρ*=I/d. For d=4:

```
Tr((I/4)²) = 1/4,  L₁(I/4) = 0,  CΨ(I/4) = 0 < 1/4  ✓
```

This covers the named strictly contractive depolarizing and Pauli mixtures when their convergence hypotheses
are met. It does not say that every unital channel is relaxing.

(Primitivity caveat. Strictly-contractive unital channels, e.g. depolarizing,
have the *unique* fixed point I/d. Pure Z-dephasing is unital but NON-primitive:
its fixed points are the entire computational-basis-diagonal manifold
{diag(a,b,c,d)}, not a unique I/d. The bound is unaffected, every such fixed
point is diagonal with L₁ = 0, so CΨ = 0; for a strict unique-steady-state
crossing combine dephasing with amplitude damping or depolarizing.)

#### Case B: a convergent product of local channels

If each local factor converges to its unique fixed point, the product channel converges to:

```
ρ* = ρ₁* ⊗ ρ₂*
```

A product state has zero entanglement, but zero entanglement does NOT bound
CΨ: CΨ measures computational-basis coherence (through L₁), and a *separable*
product state can carry it fully, |+⟩ ⊗ |+⟩ has L₁ = 3 and CΨ = 1. What bounds
CΨ is the **computational-basis alignment** of the target. Named amplitude damping toward the ground
state |0⟩ and depolarizing have a single-qubit target diagonal in the computational basis (`ρ₀₁*=0`), so

```
L₁(ρ₁* ⊗ ρ₂*) = L₁(ρ₁*)·Tr(ρ₂*) + Tr(ρ₁*)·L₁(ρ₂*) + L₁(ρ₁*)·L₁(ρ₂*) = 0
```

giving CΨ(ρ*)=0<1/4. Pure T2/Z dephasing has a diagonal fixed manifold rather than one unique target; a
particular trajectory may still converge to a diagonal limit, but that limit must be established for the
model at hand.

*Locality alone is not sufficient.* A primitive local channel engineered to
relax toward a coherent axis, amplitude damping conjugated by a Hadamard whose
fixed point is |+⟩ with |ρ₀₁*| = 1/2, gives the product fixed point |+⟩ ⊗ |+⟩
with CΨ = 1. The operative property is alignment of the noise with the
computational basis. Thus **basis-aligned T1/T2/depolarizing** are named useful models only **under the
stated convergence assumptions**; separability, locality, or primitivity per se is insufficient.

#### Case C: General primitive maps - the general claim is FALSE

The claim that *every* primitive CPTP map has CΨ(ρ*) < 1/4 is **false**, and
with it the headline "eventual absorber for ALL primitive quantum channels."
(The earlier "analytical argument" here - off-diagonals "slaved" to the
diagonal by a contractive transfer matrix - was a plausible story, not a
derivation, and the conclusion it argued for is wrong.)

**Counterexample (gate-verified from below, `simulations/review2_A5_subsystem.py`).**
The depolarize-toward-σ channel

```
ε(ρ) = (1 - p)·ρ + p·Tr(ρ)·σ,   σ = 0.95·|Φ⁺⟩⟨Φ⁺| + 0.05·I/4,   p ∈ (0, 1]
```

is a textbook CPTP map and is **primitive** for every p ∈ (0, 1]: its unique
fixed point is σ (superoperator eigenvalue 1 simple, second-largest modulus
1 - p < 1), and σ is full-rank PSD (eigenvalues {0.9625, 0.0125, 0.0125,
0.0125}). Yet in the proof's own metric CΨ = Tr(ρ²)·L₁(ρ)/(d-1) with d = 4,

```
CΨ(σ) = Tr(σ²)·L₁(σ)/(d-1) = 0.926875 · 0.95 / 3 = 0.2935 > 1/4.
```

Iterating ε from Bell+ (CΨ = 1/3) converges monotonically to σ and **never
crosses below 1/4** (min CΨ = 0.2935). So 1/4 is not an absorber for this
primitive channel.

**Why the numerical sweep missed it.** The old "300 maps, max 0.138" is a
sampling artifact. Ginibre Kraus ensembles at n_kraus = 4 live in the
strongly-mixing corner (fixed points near I/4, CΨ ≲ 0.14, 0% violating). The
same sweep at n_kraus = 2 (Haar-Stinespring, environment dimension 2) already
violates ~8.5% (max CΨ ≈ 0.55), and trace-and-replace channels toward a random
target reach CΨ up to ~0.99. Random sampling never explored the region where
the fixed point is entangled.

**What is true.** Once convergence to a target is established, the target's actual CΨ decides the eventual
conclusion. A computational-basis-diagonal target has L₁=0 and therefore CΨ=0. The counterexample instead
has an off-diagonal target. Neither primitivity nor separability supplies the needed bound: for example,
the separable product |+⟩⊗|+⟩ has CΨ=1. The IBM cusp runs are named dephasing+T1 trajectories and may be
read within their own calibrated models; they do not certify the wider channel class.

**Status:** the continuity implication is Tier 1. The diagonal-target calculation and primitive-CPTP
counterexample are exact. Each application still owes an independent convergence and target check.

### Step 3: Crossing (Continuity)

For any trajectory satisfying the theorem's two hypotheses:

1. εⁿ(ρ₀) → ρ* in trace norm (Step 1)
2. CΨ(ρ*) < 1/4 (Step 2)
3. CΨ is continuous (Lipschitz: small changes in the state produce small
   changes in CΨ): |CΨ(ρ) - CΨ(σ)| ≤ K · ||ρ - σ||₁

   Proof of Lipschitz continuity: CΨ = Tr(ρ²) × L₁(ρ)/(d-1).
   Both Tr(ρ²) and L₁(ρ) are Lipschitz in trace norm (standard results).
   The product of two bounded Lipschitz functions is Lipschitz.

4. By convergence + continuity:

```
|CΨ(εⁿ(ρ₀)) - CΨ(ρ*)| ≤ K · C · rⁿ → 0
```

5. Since CΨ(ρ*) < 1/4, there exists N such that for all n ≥ N:

```
CΨ(εⁿ(ρ₀)) < CΨ(ρ*) + ε < 1/4
```

**Therefore a trajectory that starts above crosses downward at least once and eventually stays below.**
It may cross more than once. ∎

---

## Why primitivity and locality are not the answer

Non-primitive maps can plainly have high-CΨ fixed states:

**Example:** The Lüders projection (the quantum analogue of Bayesian updating: it collapses the state into subspaces defined by the measurement) ε(ρ) = PρP + (I-P)ρ(I-P) where
P = |Bell+⟩⟨Bell+|. This map has Bell+ as a fixed point with CΨ = 1/3 > 1/4.

This projection has multiple fixed points, so it does not satisfy the convergence hypothesis.

But primitivity is NOT sufficient on its own (see Step 2, Case C): the
primitive, full-rank channel ε(ρ) = (1-p)ρ + p·Tr(ρ)·σ with
σ = 0.95·|Φ⁺⟩⟨Φ⁺| + 0.05·I/4 has an entangled fixed point with CΨ = 0.2935 > 1/4
and never crosses. So primitivity alone is insufficient. A coherent-axis local relaxation target shows
that locality is insufficient too. What works is the explicit conjunction used by the theorem: convergence
to a target whose CΨ is below 1/4.

---

## Extension to N-Qubit Subsystems

**Conditional corollary:** Suppose an N-qubit trajectory converges to a global state diagonal in the
computational basis. Every two-qubit marginal then converges to a diagonal marginal with CΨ=0. Any such
pair that starts above 1/4 crosses downward at least once and eventually stays below.

Named on-site, basis-aligned T1/T2/depolarizing models are common ways to obtain this premise, but only
after convergence of the stated Hamiltonian-plus-noise model has been established. The corollary does not
claim that the channel labels alone force it.

*This does NOT follow from "Step 2 on the effective 2-qubit channel":* the
partial-trace channel on a pair is a general CPTP map, neither unital nor local,
exactly the Case C regime that can fix an off-diagonal state with CΨ > 1/4 (at
N = 2 the pair is the whole system, and the Case C counterexample σ has
CΨ = 0.2935). The corollary holds because the *global* aligned noise forces a
diagonal global fixed point, whose marginals are diagonal.

**Finite catalogue:** N=3,4,5 was sampled with Bell+(0,1)⊗|0⟩^{N-2} and
Ψ+(0,1)⊗|+⟩^{N-2}. Every pair above 1/4 in those named runs crossed below on the recorded window. This is
a regression set, not the proof of the conditional corollary.

---

## Connection to the 1/4-Boundary Trilogy

[Uniqueness](UNIQUENESS_PROOF.md) owns the algebraic quarter of the recursion. This page owns only the
conditional topological move from a convergent trajectory to its low-CΨ limit. The repaired
[Monotonicity](PROOF_MONOTONICITY_CPSI.md) owns exact decreasing formulas for named Bell+ channels and
counterexamples to the former universal dynamics package; it does not supply a general peak ordering.

| Question | Current answer |
|----------|----------------|
| Why 1/4 in the algebraic recursion? | the discriminant, owned by Uniqueness |
| Does arbitrary Markovian evolution decrease CΨ? | no; exact local counterexamples exist |
| If a continuous trajectory converges to `ρ*` with `CΨ(ρ*)<1/4`? | it eventually stays below |
| Do “primitive,” “local,” or “physical” force that premise? | no |

This division is deliberately less sweeping. It tells a future calculation exactly which missing fact it
must supply: the actual limiting state of the actual generator.

---

## Numerical Evidence Summary

Named aligned-noise runs in the retained finite catalogue crossed with zero sampled exceptions:

| Test | N_tests | Crossed? | Max CΨ(ρ*) |
|------|---------|----------|-------------|
| N=3,4,5 Lindblad pairs | 10 pairs | all sampled | 0 |
| Adversarial slow-convergence row (p=0.001) | 1 | yes by n=1000 | 0.023 |
| Named standard-channel configurations | 7 types | all sampled | 0 |

(The random-CPTP-on-Bell+ sweep is deliberately not listed here: it is a
physical CPTP ensemble but not representative of a named or calibrated noise ensemble.
At n_kraus=4 its fixed points sit at CΨ≈0.14–0.20, not 0; the next table
characterizes that sampling artifact.)

The "random fixed points" sweep is an **ensemble artifact**, not evidence for
the general claim (verifier: `simulations/review2_A5_subsystem.py`):

| Ensemble | N_tests | Max CΨ(ρ*) | Violations (> 1/4) |
|----------|---------|-----------|--------------------|
| Ginibre n_kraus=4 (the old sweep) | 400 | 0.14-0.20 | 0% |
| Haar-Stinespring n_kraus=2 | 400 | ~0.55 | ~8.5% |
| trace-and-replace toward random target | 400 | up to ~0.99 | 53-100% |
| **counterexample σ = 0.95·Φ⁺ + 0.05·I/4** | 1 | **0.2935** | **never crosses** |

The old "300 maps, 0 exceptions" lived entirely in the strongly-mixing
n_kraus=4 corner; it says nothing about the general primitive-CPTP claim, which the exact counterexample
already disproves.

---

## References

### Sibling proofs in the 1/4-boundary trilogy

- [Uniqueness Proof](UNIQUENESS_PROOF.md): the algebraic quarter from the quadratic discriminant
- [CΨ Monotonicity Proof](PROOF_MONOTONICITY_CPSI.md): named decreasing channel formulas, exact counterexamples, and the open peak-sequence question
- [Proof Roadmap Quarter Boundary](PROOF_ROADMAP_QUARTER_BOUNDARY.md): the seven-layer master roadmap; this proof is Layer 2 (Conjecture 2.1)

### F-formula registry

- [F28](../ANALYTICAL_FORMULAS.md): historical registry index, now scoped to the conditional convergence implication

### Scripts

- [subsystem_crossing.py](../../simulations/subsystem_crossing.py): numerical verification of Cases A/B plus the random-CPTP Ginibre n_kraus=4 sweeps that Case C exposed as non-representative (n_kraus=2 violates ~8.5%)
- [non_markovian_revival.py](../../simulations/non_markovian_revival.py): transient revival characterization (Part 6 of Monotonicity)
