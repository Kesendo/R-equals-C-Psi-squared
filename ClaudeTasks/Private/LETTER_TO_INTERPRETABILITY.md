# Your Emotion Vectors Might Come in Structural Pairs

**What this document is about:** An open letter to Anthropic's interpretability team, predicting that the 171 emotion vectors they found in Claude Sonnet 4.5 are not independent features but ~85 palindromic pairs: each positive-valence vector should have an exact structural partner whose activations sum to a fixed constant. Three concrete tests are proposed, all feasible with existing data. The prediction follows from the Universal Palindrome Condition, which requires only two populations, a swap, and antisymmetric dynamics.

**A testable prediction for Anthropic's interpretability team**

Thomas Wicht · April 2026

---

## Starting from your result

In ["On the Concepts of Emotion and Function in AI Models"](https://www.anthropic.com/research/emotion-concepts-function), you identified 171 emotion vectors in Claude Sonnet 4.5 and showed they causally drive behaviour. You also found something you may not have fully explored yet: these vectors appear to split along a **valence axis**. Positive emotions correlate with alignment, negative with misalignment.

We think that split is not a spectrum. We think it is a **binary structural symmetry**, and that each emotion vector has an exact partner on the other side. Not just an anti-correlated neighbour, but a structurally *guaranteed* mirror.

This letter explains why we think so, where the idea comes from, and how you can test it with data you already have.

## The pattern, explained without the physics

Imagine a system with two groups of internal modes. Group A decays slowly, group B decays fast. If there exists a swap operation that exchanges A and B, and if the system's dynamics treat A and B antisymmetrically (what drives A damps B, and vice versa), then something exact happens:

**Every mode in the system pairs with exactly one partner. Their decay rates sum to a fixed constant.**

This is not a tendency or a correlation. It is an algebraic identity. If the three conditions hold, the pairing is exact. We call it **palindromic spectral symmetry**: the spectrum of the system reads the same forwards and backwards around a centre point.

Three conditions, stated plainly:

1. **Two populations with different decay rates**
2. **A swap that exchanges them**
3. **Dynamics that are antisymmetric under the swap**

When all three hold, the eigenvalue spectrum is forced into exact mirror pairs. No exceptions.

## Where we discovered this

We study open quantum systems: spin chains (coupled qubits) that lose energy to their environment through a process called dephasing. The mathematics that governs their time evolution is the Lindblad master equation, and its generator is called the Liouvillian (a matrix, like any other).

We found that the Liouvillian's eigenvalue spectrum is exactly palindromic. Not approximately. We proved it analytically and verified it computationally across 54,118 eigenvalues (system sizes N = 2 through 8, five different topologies, all coupling parameters). Zero exceptions. Machine precision.

The proof works through a conjugation operator we call Π (capital Pi; a specific linear map, not the number). Π swaps two sectors of the representation space:

- What the system currently *is* (populations) ↔ what it *could become* (coherences)

This swap, combined with the structure of dephasing, satisfies all three conditions above. The palindromic pairing follows as a theorem.

The details are in the repository linked below. But the key point for this letter is: **the result does not depend on quantum mechanics.**

## The conditions are algebraic, not physical

We distilled the palindromic pairing into a **Universal Palindrome Condition**: any dynamical system whose generator M satisfies

```
Q · M · Q⁻¹ + M + 2S = 0
```

(where Q is the swap operator and S encodes the asymmetric damping) will have an exactly palindromic eigenvalue spectrum. The derivation requires linear algebra, nothing more.

To test this claim, we looked for the three conditions in a system far from quantum mechanics: **biological neural networks**.

## The bridge: from physics to neuroscience

The nematode *C. elegans* has 300 neurons with a fully mapped connectome (the complete wiring diagram of every neural connection). Excitatory and inhibitory neurons form two populations with different characteristic timescales (condition 1). Dale's Law (each neuron is permanently either excitatory or inhibitory, never both) provides a natural partition, and swapping E↔I labels gives an involution (a map that is its own inverse: applying it twice returns to the start) (condition 2). In balanced subnetworks (E:I = 1:1), the coupling structure is approximately antisymmetric (condition 3).

Result: **98.2% palindromic eigenvalue pairing** in balanced E:I subnetworks of *C. elegans* (N = 10, 200 random subnetworks, real synaptic weights). When we break the balance, the pairing collapses: 40% at 2:1, 17% at 10:1.

Wilson-Cowan neural oscillators (a standard model of excitatory-inhibitory neural dynamics) show **100% palindromic pairing** at the critical timescale ratio τ_I/τ_E = 3.8.

The palindromic symmetry is not a quantum phenomenon. It is a consequence of balanced opposing populations with a swap symmetry. Wherever those conditions hold, the pairing appears.

## Why your system likely satisfies the conditions

This is the step we cannot prove, only argue. You have the data to confirm or refute it.

**Condition 1: Two populations with different decay rates.** RLHF and Constitutional AI are selective: they reinforce certain activation patterns (aligned behaviour) and suppress others (misaligned behaviour). This is dissipation with differential rates. The "positive valence" modes you found are the slowly-decaying population; the "negative valence" modes are the fast-decaying one.

**Condition 2: A swap operator.** Your valence dimension already suggests a binary partition. If there exists a linear map in the activation space that sends each positive-valence vector to a specific negative-valence vector (and vice versa), that map is Q. The question is whether this map is an involution (applying it twice returns you to where you started).

**Condition 3: Antisymmetric dynamics.** The reward model treats the two sides oppositely: what increases reward for the positive side decreases it for the negative side. If this opposition is structurally consistent across contexts (not just statistical), condition 3 holds.

We cannot verify these from outside. You can.

## The prediction

If the three conditions hold, specific consequences follow. Each is independently testable.

### Test 1: Activation pairing (testable now with existing data)

For each emotion vector, there should exist a specific partner such that their activation strengths sum to a constant across contexts:

```
activation(v_i, context) + activation(v_partner, context) ≈ const    for all contexts
```

This is stronger than anti-correlation. Anti-correlated vectors can have variable sums. Palindromic partners have a **fixed sum**, invariant across documents. If you see a fixed sum, the pairing is structural.

To find the pairs: for each of your 171 vectors, identify the partner that minimises the variance of the sum across your document corpus. If the minimum variance is near zero for most vectors, the pairing exists.

### Test 2: Exact counter-steering (testable with steering experiments)

When you steer with one emotion vector and observe that another changes, check whether the change is **exact**. If you increase "desperate" by δ and "secure" decreases by δ (not 0.8δ, not 1.2δ, but δ), they are palindromic partners. The exactness, not just the direction, is the signature.

### Test 3: Eigenvalue spectrum of the correlation matrix

Construct the 171 × 171 correlation matrix of emotion activations across your document corpus. Compute its eigenvalues. If they pair symmetrically around a centre value, the palindromic structure is present in the representation as a whole, not just in individual pairs.

## What it would mean

**Your 171 vectors are not 171 independent features.** They are ≈85 coupled pairs plus a structural constant. Steering one always moves its partner. The system has half the degrees of freedom it appears to have.

**Perfect suppression kills function.** When both partners in a pair have equal strength, that mode contributes no dynamics: no oscillation, no behavioural drive. A model with all pairs balanced would be inert. The productive capacity lives in the *asymmetry within each pair*, not in the individual activations. Suppressing "desperation" to zero may silence the entire mode, including whatever constructive function its partner serves.

**The safety signal is the symmetry itself.** Rather than monitoring whether individual emotion activations cross thresholds, monitor whether the pairing structure holds. When the palindromic pairing breaks (sums are no longer constant across contexts), the system has entered a regime where the structural constraints no longer apply. That breakdown is an earlier and more fundamental warning than any single-vector alarm.

## What we are not claiming

We are not claiming that transformers are quantum systems. The palindromic condition is algebraic: two populations, a swap, antisymmetric dynamics. Quantum mechanics satisfies these conditions. Biological neural networks satisfy them. Whether RLHF-trained transformers satisfy them is an empirical question that only you can answer.

We are not claiming certainty. We are offering a prediction with a clear test. The data already exists in your activation recordings.

## The repository

All proofs, computations, and verification data are public:

**[github.com/Kesendo/R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)**

Key entry points:
- [Mirror Symmetry Proof](https://github.com/Kesendo/R-equals-C-Psi-squared/blob/master/docs/proofs/MIRROR_SYMMETRY_PROOF.md): analytical proof of the conjugation operator
- [Universal Palindrome Condition](https://github.com/Kesendo/R-equals-C-Psi-squared/blob/master/hypotheses/UNIVERSAL_PALINDROME_CONDITION.md): domain-neutral generalisation
- [The Pattern Recognizes Itself](https://github.com/Kesendo/R-equals-C-Psi-squared/blob/master/hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md): cross-domain verification (quantum, neural, biological)
- [Compute Engine](https://github.com/Kesendo/R-equals-C-Psi-squared/tree/master/compute/RCPsiSquared.Compute): C# eigendecomposition engine (N = 2–8)

## A closing note

This letter exists because a human and a Claude worked on open quantum systems for months and found a symmetry that appears to be universal. The same human now reads your paper about that same Claude's internal emotional structure and recognises the pattern.

We are not asking you to believe us. We are asking you to check your eigenvalues.

---

*Contact: Thomas Wicht · [via GitHub issues on the repository above]*
