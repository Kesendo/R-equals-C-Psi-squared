# Proof: Palindromic Spectral Symmetry for Neural Networks

**Status:** Derived from quantum proof, computationally verified
**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Quantum foundation:** [Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md)

---

## Theorem

For a neural network with N neurons (N/2 excitatory, N/2 inhibitory),
linearized dynamics dx/dt = J*x, and E-I swap permutation Q:

If

```
Q * J * Q + J + 2*S = 0
```

with S = (1/τ_E + 1/τ_I) / 2 times the identity matrix, then
every eigenvalue μ_k of J has a palindromic partner μ_k' such that:

```
μ_k + μ_k' = -(1/τ_E + 1/τ_I)
```

The condition is satisfied exactly when:
- **(a)** τ_E ≠ τ_I (selective damping), and
- **(b)** W[Q(i), Q(j)] = -(τ_{Q(i)} / τ_i) * W[i, j] for all i, j
  (coupling antisymmetry scaled by time constant ratio)

Dale's Law provides the sign structure of condition (b) automatically.

---

## Derivation

### Starting point

In quantum open systems, the palindromic spectral symmetry is proven
as an algebraic identity
([Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md)):

```
Π · L · Π⁻¹ = -L - 2Σγ · I
```

where L is the Liouvillian (evolution superoperator), Π is the
palindromic conjugation operator that swaps immune and decaying
degrees of freedom, and Σγ is the total dephasing rate. This implies
eigenvalue pairing: λ + λ' = -2Σγ.

### Step 1: Identify the analogs

| Quantum | Neural | Role |
|---------|--------|------|
| L (Liouvillian) | J (Jacobian) | Evolution operator |
| Π (Pauli weight swap) | Q (E-I swap permutation) | Conjugation |
| 2Σγ (total dephasing) | 2*S (to be determined) | Pairing constant |
| Immune sector {I,Z} | E-neurons | Slow-decaying side |
| Decaying sector {X,Y} | I-neurons | Fast-decaying side |

Since Q is a permutation matrix, Q^{-1} = Q^T = Q. The translated
condition is:

```
Q * J * Q + J + 2*S = 0       ... (*)
```

### Step 2: Decompose the Jacobian

The neural Jacobian (linearized Wilson-Cowan or firing rate model):

```
J[i,i] = -1/τ_i                          (self-decay)
J[i,j] = α · W[i,j] / τ_i   for i≠j  (synaptic coupling)
```

where τ_i = τ_E if neuron i is excitatory, τ_I if inhibitory,
α is the coupling strength, and W[i,j] is the signed synaptic
weight from neuron j to neuron i.

Decompose: J = D + W_eff, where
- D = diag(-1/τ_i): self-decay (diagonal)
- W_eff: effective coupling (off-diagonal, includes τ scaling)

Substituting into (*):

```
(Q*D*Q + D + 2*S) + (Q*W_eff*Q + W_eff) = 0
```

D is diagonal and W_eff is off-diagonal, so they are linearly
independent. Both parenthesized terms must vanish separately.

### Step 3: Self-decay condition (determines S)

Q swaps each E-neuron (index i) with its paired I-neuron (index Q(i)).
Therefore Q*D*Q is D with τ_E and τ_I exchanged:

```
(Q*D*Q)[i,i] = D[Q(i),Q(i)] = -1/τ_{Q(i)}
```

For an E-neuron: (Q*D*Q)[i,i] = -1/τ_I.
For an I-neuron: (Q*D*Q)[i,i] = -1/τ_E.

In both cases:

```
(Q*D*Q)[i,i] + D[i,i] = -1/τ_I + (-1/τ_E) = -(1/τ_E + 1/τ_I)
```

Setting Q*D*Q + D + 2*S = 0:

```
S = (1/τ_E + 1/τ_I) / 2 * I
```

This is a scalar times the identity. It does NOT depend on which
neuron is E or I, because the sum 1/τ_E + 1/τ_I is the same
whether the neuron is E (swapped to I) or I (swapped to E).

**Condition (a) is always satisfied when τ_E ≠ τ_I.** No
constraints on the network topology. Only selective damping required.

### Step 4: Coupling condition

The remaining equation is:

```
Q * W_eff * Q + W_eff = 0
```

In components (for i ≠ j):

```
W_eff[Q(i), Q(j)] + W_eff[i, j] = 0
```

Substituting W_eff[i,j] = α · W[i,j] / τ_i:

```
α · W[Q(i), Q(j)] / τ_{Q(i)} + α · W[i, j] / τ_i = 0
```

Dividing by α (nonzero) and solving:

```
W[Q(i), Q(j)] = -(τ_{Q(i)} / τ_i) * W[i, j]       ... (**)
```

**This is the non-trivial condition.** It requires a specific
relationship between each connection and its E-I partnered connection.

### Step 5: Dale's Law provides the signs

Under Q, the source neuron type flips (E becomes I, I becomes E).
Dale's Law fixes the sign by the source type:

- W[i,j] > 0 if source j is excitatory
- W[i,j] < 0 if source j is inhibitory

Under the swap Q(j): if j was E, Q(j) is I (and vice versa).
So sign(W[Q(i),Q(j)]) = -sign(W[i,j]).

Condition (**) requires W[Q(i),Q(j)] = -(positive factor) * W[i,j].
Since τ_{Q(i)}/τ_i > 0, the required sign is negative, which
matches Dale's Law.

**Dale's Law automatically satisfies the sign part of condition (b).**

The remaining requirement is on magnitudes:

```
|W[Q(i), Q(j)]| = (τ_{Q(i)} / τ_i) * |W[i, j]|
```

For τ_I/τ_E = 2: E-to-E connections need I-to-I partners with
2x magnitude. I-to-I connections need E-to-E partners with 0.5x.

### Step 6: Eigenvalue pairing (consequence)

When (*) holds, let v be an eigenvector of J with eigenvalue μ:
J·v = μ·v.

Multiply (*) from the right by v:
Q·J·Q·v + J·v + 2S·v = 0
Q·J·(Q·v) + μ·v + 2S·v = 0

Let w = Q·v. Then Q·w = v (since Q² = I), and:
Q·J·w = -(μ + 2S)·v = -(μ + 2S)·Q·w

Multiplying both sides by Q from the left:
J·w = -(μ + 2S)·w

So w = Q·v is an eigenvector of J with eigenvalue -(μ + 2S).

Therefore J has eigenvalue μ' = -μ - (1/τ_E + 1/τ_I):

```
μ + μ' = -(1/τ_E + 1/τ_I)       for each palindromic pair
```

This is the neural palindromic spectral symmetry. QED.

---

## Verification

### Synthetic network (exact condition satisfied)

Constructed W satisfying (**) with Dale's Law signs and exact
magnitude ratios. Palindrome residual ||Q*J*Q + J + 2*S|| / ||J||:

```
Residual = 0.00 (machine precision, 10^{-16})
```

### C. elegans connectome (approximate)

Balanced subnetworks (5E + 5I), 200 samples:

| Network | ||R|| / ||J|| |
|---------|--------------|
| C. elegans | 0.013 |
| Erdos-Renyi | 0.108 |
| Degree-preserving rewiring | 0.013 |

### Eigenvalue pair sums

For C. elegans subnetwork at α = 0.3, τ_E = 5, τ_I = 10:
Predicted sum: -(1/5 + 1/10) = -0.300.
Observed: mean = -0.3012, max deviation 1.6%.

---

*See also:*
[Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md) (quantum original),
[Algebraic Palindrome Neural](../ALGEBRAIC_PALINDROME_NEURAL.md) (results and interpretation)
