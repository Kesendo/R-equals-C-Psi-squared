# Proof: Conditional Palindromic Spectral Symmetry for Neural Networks

**Status:** Tier 1 conditional algebra, with constructed computational gates
**Date:** March 27, 2026
last refreshed 2026-09-05 (the change history lives in git)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Quantum foundation:** [Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md)

## What this document is about

A neural Jacobian has a palindromic spectrum when an involutive permutation
conjugates it to its negative plus one scalar shift. The theorem gives the
exact entry conditions and transports the full complex spectrum, including
multiplicity. It does not assert that a biological network satisfies those
conditions, or that a paired spectrum is real, stable, or silent.

The algebra is indexed by [F36 and F37](../../ANALYTICAL_FORMULAS.md#f36-neural-palindrome-condition-tier-1-derived-algebra).
The quantum proof supplies the conjugation pattern; the neural conditions
below are checked on the neural Jacobian itself. The
[C. elegans event record](../../../experiments/NEURAL_GAMMA_CAVITY.md)
supplies a support obstruction, not a biological confirmation.

## Theorem

Let J ∈ ℂⁿˣⁿ, n ≥ 1, be a neural Jacobian, written J = D + W_eff,
where D = diag(d₁,…,dₙ) and W_eff has zero diagonal. Let s ∈ ℂ be a scalar
and Q an **involutive permutation matrix**, Q² = I. Write Q(i) for the
permuted seat. Then

```
Q·J·Q + J + 2s·I = 0                                      (*)
```

holds if and only if both entry conditions hold:

- **(a) Diagonal:** dᵢ + d_Q(i) + 2s = 0 for every i.
- **(b) Off-diagonal:** W_eff[Q(i),Q(j)] + W_eff[i,j] = 0 for every i ≠ j.

When (*) holds, the eigenvalue **multiset** of J is invariant under
μ ↦ −μ − 2s, with algebraic multiplicity preserved. Its centre is −s.
The identity imposes no general requirement of a real spectrum, linear
stability, or absence of oscillatory modes.

## Derivation

### Step 1: The conjugation pattern

The quantum palindrome has the form Π·L·Π⁻¹ = −L − 2Σγ·I.
For the neural statement Q² = I gives Q⁻¹ = Q, so (*) is the similarity
Q·J·Q⁻¹ = −J − 2s·I. The shared structure is this operator identity;
Dale's Law by itself does not supply its neural hypotheses.

### Step 2: Separate the diagonal and coupling

Substitute J = D + W_eff into (*):

```
(Q·D·Q + D + 2s·I) + (Q·W_eff·Q + W_eff) = 0.
```

Conjugation by the involutive permutation Q preserves diagonal matrices
and matrices with zero diagonal. The two summands therefore vanish
separately, exactly when (a) and (b) hold. In particular,
(Q·J·Q)[i,j] = J[Q(i),Q(j)].

### Step 3: Self-decay condition

For positive time constants τᵢ, the no-self-coupling neural model has
dᵢ = −1/τᵢ. Condition (a) reads

```
1/τᵢ + 1/τ_Q(i) = 2s.
```

If every seat has one of two time constants τ_E and τ_I, an involutive Q
that exchanges every excitatory seat with an inhibitory seat gives

```
s = ½(1/τ_E + 1/τ_I).
```

For this choice of s and τ_E ≠ τ_I, (a) requires opposite-type exchange
at every seat. At τ_E = τ_I = τ, the diagonal is scalar and (a) holds
for every involutive Q with s = 1/τ; opposite types are not required.
Equal time constants do not supply condition (b).

For the general theorem a fixed seat Q(i) = i is allowed precisely when
its diagonal rate dᵢ = −s, subject also to (b). With unequal τ_E and τ_I
and the reciprocal-sum choice of s above, neither type can occupy such a
fixed seat. Having no fixed seats is not by itself enough: an involution
that pairs seats of the same type still fails (a) in that setting.

### Step 4: Coupling condition

In the uniform-gain model, let T = diag(1/τᵢ), let W have zero diagonal,
and set W_eff = α·T·W. The source index is the **column**: W[i,j] is the
weight from j to i. Condition (b) becomes

```
α·W[Q(i),Q(j)]/τ_Q(i) + α·W[i,j]/τᵢ = 0       (i ≠ j).
```

Only when α ≠ 0 may this be divided to give

```
W[Q(i),Q(j)] = −(τ_Q(i)/τᵢ)·W[i,j]             (**)
```

At α = 0, W_eff = 0 and (b) holds for every W; (**) is then not
necessary. This is why W_eff owns the theorem. A self-coupling contribution
to J belongs in D and must be included in (a).

A Wilson-Cowan linearization can have a different sigmoid slope at each
operating point. Such row gains belong in W_eff. The theorem applies to
that effective matrix directly; (**) assumes a single common gain α.

### Step 5: What Dale's Law supplies

For real weights and positive time constants, Dale's Law fixes the sign of
each nonzero outgoing weight by the source type. If Q exchanges the types
and both W[i,j] and W[Q(i),Q(j)] are nonzero, their signs are opposite,
as (**) requires. Dale's Law is insufficient for (b): it supplies neither
the Q-symmetric support nor the scaled magnitudes

```
|W[Q(i),Q(j)]| = (τ_Q(i)/τᵢ)·|W[i,j]|.
```

The committed connectome uses **rows** for presynaptic sources. Transposing
preserves the support-symmetry question, but a magnitude calculation must
respect the layout: transposing a row-scaled matrix turns that scaling
into column scaling.

### Step 6: Full complex pairing and character transport

Since Q² = I, multiply (*) by Q on the right to obtain

```
J·Q = −Q·J − 2s·Q.
```

Thus Jv = μv implies J(Qv) = (−μ − 2s)Qv. Q is invertible, so Qv
is nonzero whenever v is. Moreover, J and −J − 2s·I are similar, which
preserves their characteristic polynomials and proves the multiset
statement, including algebraic multiplicity even when J is defective.
The transport also maps generalized eigenspaces:

```
(J + (μ + 2s)I)ᵏ·Q = (−1)ᵏ·Q·(J − μI)ᵏ.
```

At degeneracy the invariant subspaces are the appropriate comparison;
an eigensolver's chosen eigenvectors need not coincide with the Q images
of another chosen basis. With the two-time-constant value of s, F37 is

```
μ + μ′ = −(1/τ_E + 1/τ_I).
```

This includes the imaginary parts, not merely the real parts. A mode at
μ = −s maps into the same eigenvalue subspace.

## Verification

### Constructed conditional gates

The [Python translation gate](../../../simulations/neural/neural_translation_gate.py)
checks the scalar identity, full complex multiset pairing, multiplicity
and imaginary-part negative controls, and Q transport into partner
subspaces. Its constructed-network census includes oscillatory and
unstable instances satisfying the scalar identity.

The typed [MirrorWorld.NeuralPalindrome](../../../compute/MirrorWorld/NeuralPalindrome.cs)
owns the entrywise identity check for real matrices and involutive
permutations. Its [focused tests](../../../compute/MirrorWorld.Tests/NeuralPalindromeTests.cs)
and `neural` run mode include the exactly representable construction

```
J = [[−0.5, −0.25], [0.25, −0.25]],  Q = (0 1),  s = 0.375.
Q·J·Q + J + 2s·I = 0,
μ = −0.375 ± (√3/8)i.
```

The fixed-seat control J = diag(−0.5, −0.25, −0.5), Q = (0 1)(2),
s = 0.375 has zero off-diagonal residual but maximum absolute entry
residual 0.25 in (*). It demonstrates why the scalar diagonal condition
cannot be replaced by an off-diagonal fit.

Run from the repository root:

```
python simulations/neural/neural_translation_gate.py
dotnet test compute/MirrorWorld.Tests --filter FullyQualifiedName~NeuralPalindromeTests
dotnet run --project compute/MirrorWorld -- neural
```

These constructions verify the conditional algebra. They are not a
biological brain result or a hardware Confirmation.

### C. elegans support obstruction

On the committed chemical connectome under Dale signs, there are 253
non-empty excitatory outgoing rows and 18 non-empty inhibitory outgoing
rows. A sign-reversing Q-symmetric support would require a bijection
between these sets. Their unequal sizes rule it out, even before scaled
magnitudes are tested, for the nonzero-gain model with positive row
scales. This full-connectome obstruction does not decide each chosen
subnetwork. The counts and gate G0b are owned by
[celegans_pairing_controls.py](../../../simulations/neural/celegans_pairing_controls.py)
using [celegans_connectome.json](../../../simulations/neural/celegans_connectome.json);
the [Neural Gamma Cavity event record](../../../experiments/NEURAL_GAMMA_CAVITY.md)
holds the support result and the pairing-analysis withdrawal. No
biological pairing verification follows from that record.

As [F137](../../ANALYTICAL_FORMULAS.md) explains, a centre fixed by the trace
carries no evidence for pairing. For even n, any partition of the full
eigenvalue multiset into n/2 pairs has mean pair sum 2·trace(J)/n.
The individual complex pairing and the operator conditions supply tests
beyond this trace identity.

See also [Algebraic Palindrome Neural](../ALGEBRAIC_PALINDROME_NEURAL.md)
for the constructed-network and connectome analyses.
