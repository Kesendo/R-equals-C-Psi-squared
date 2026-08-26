# Proof: Palindromic Spectral Symmetry for Neural Networks

**Status:** Derived from the quantum proof; the algebra is verified, the
C. elegans verification below is a NULL result and says so
**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Quantum foundation:** [Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md)

---

## What this document is about

This is the formal proof that the palindromic spectral symmetry found
in quantum spin chains also exists in neural networks. The key insight:
where an E-I swap Q turns the wiring into minus itself, every eigenvalue of the
network's dynamics matrix has a mirror partner. That is condition (b), and it is
the whole of the content. The two ingredients usually named beside it are weaker
than they look. Selective damping (different time constants) does not enable the
pairing: the diagonal half of the equation asks only that Q swap the two types,
which is automatic at uniform τ, and what τ_E ≠ τ_I adds is that only a genuine
type-swapping Q will then do (Step 3). Dale's Law (excitatory neurons always
excite, inhibitory always inhibit) fixes only the SIGNS of (b), and only where W
is nonzero, leaving the zero pattern and the magnitudes to be checked (Step 5).
On the C. elegans connectome no admissible Q exists at all, on a count. The proof translates
the quantum Π conjugation into a neuroscience E-I swap operator Q
and derives the exact conditions under which the pairing holds.

---

## Theorem

For a neural network with N neurons (N/2 excitatory, N/2 inhibitory),
linearized dynamics dx/dt = J*x (where J is the Jacobian, the matrix
of partial derivatives that governs small perturbations), and a permutation Q:

Q is written below as the E-I swap because that is the case of interest, but
the conditions are stated for an arbitrary permutation on purpose. Reading Q as
type-swapping by definition makes condition (a) true by definition and hides
what Step 3 derives, which is that (a) is exactly the demand that Q be
type-swapping when the two time constants differ, and no demand at all when
they do not.

If

```
Q * J * Q + J + 2*S = 0
```

with S = (1/τ_E + 1/τ_I) / 2 times the identity matrix, then
every eigenvalue μ_k of J has a palindromic partner μ_k' such that:

```
μ_k + μ_k' = -(1/τ_E + 1/τ_I)
```

The condition is satisfied exactly when both of the following hold, and this is an
if-and-only-if: Step 2 splits (*) into a diagonal and a hollow equation living in
complementary subspaces, each of which Q maps into itself, so (*) is equivalent to
the two of them together.
- **(a)** the diagonal condition: 1/τ_{Q(i)} + 1/τ_i = 1/τ_E + 1/τ_I at every
  seat i. This is exactly the diagonal half of (*), no more and no less
  (Step 3). A Q that swaps the two types satisfies it for any time constants;
  the converse holds only when τ_E ≠ τ_I, since at uniform τ both sides are 2/τ
  for EVERY permutation and (a) imposes nothing at all. So selective damping is
  not a hypothesis that buys the pairing, it is what makes (a) bite on Q, and
- **(b)** W[Q(i), Q(j)] = -(τ_{Q(i)} / τ_i) * W[i, j] for all i ≠ j
  (coupling antisymmetry scaled by time constant ratio)

The i ≠ j is load-bearing and is the quantifier Step 4 derives. The diagonal
entries W[i,i] never enter J, which carries -1/τ_i there by construction, so a
version of (b) quantified over all i, j would be strictly stronger than (*) and
the if-and-only-if would fail: a network whose only nonzero weight is a self-loop
satisfies (*) exactly while failing (b) at i = j.

Dale's Law provides the sign structure of condition (b) automatically, but
only where W is nonzero; the zero pattern and the magnitudes are separate
requirements and are the ones that bind in practice. Step 5 works this out.

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
| Immune sector {I,Z} | E-neurons | One side |
| Decaying sector {X,Y} | I-neurons | The other |

The table names the two sides, not which is slower: that is set by the τ values
and the theorem does not fix their order. Step 5 below works its example at
τ_I/τ_E = 2, where the inhibitory population is the slower one.

Since Q is a permutation matrix Q^{-1} = Q^T, and since an E-I swap is an
involution, Q^T = Q. The translated
condition is:

```
Q * J * Q + J + 2*S = 0       ... (*)
```

### Step 2: Decompose the Jacobian

The neural Jacobian (linearized Wilson-Cowan model, the standard
equations describing excitatory-inhibitory population dynamics):

```
J[i,i] = -1/τ_i                          (self-decay)
J[i,j] = α · W[i,j] / τ_i   for i≠j  (synaptic coupling)
```

where τ_i = τ_E if neuron i is excitatory, τ_I if inhibitory,
α is the coupling strength, and W[i,j] is the signed synaptic
weight from neuron j to neuron i.

Two conventions to keep straight, because both bite later.

**The gain is uniform here.** A Wilson-Cowan Jacobian carries the sigmoid slope
at each neuron's own operating point, S'(x_i*), which is a PER-ROW factor; the α
above absorbs it into one constant and so assumes every neuron sits at the same
slope. Where the slopes differ, condition (b) becomes
`W[Q(i),Q(j)] = −(τ_{Q(i)}·S'_i)/(τ_i·S'_{Q(i)})·W[i,j]`, strictly stronger.
The C. elegans block of `neural_gamma_cavity.py` uses a single f' = 0.3 and is
covered; its two-population block has f_E ≠ f_I and is not.

**The source index is the COLUMN here**, W[i,j] running from j to i. The
connectome file uses the opposite layout, rows presynaptic, so an argument
phrased in "outgoing rows" is speaking about the transpose of the W defined
here. The SUPPORT half of (b), which is all the count argument on the connectome
uses, does not care: a zero pattern is Q-symmetric exactly when its transpose is.
The MAGNITUDE half does care, and the unscaled equivalence Q W Q = −W ⟺
Q Wᵀ Q = −Wᵀ does not license it, because (b) carries the factor τ_{Q(i)}/τ_i,
which scales ROWS, and transposition turns that into a column scaling. At N = 2
with τ_E = 10 on seat 0, τ_I = 20 on seat 1 and Q the swap, the Dale-legal
W = [[0, −1], [2, 0]] (the excitatory source's column positive, the inhibitory
source's negative) satisfies (b) exactly while its transpose violates it at both
off-diagonal seats. So an
argument about magnitudes must name the layout.

Decompose: J = D + W_eff, where
- D = diag(-1/τ_i): self-decay (diagonal)
- W_eff: effective coupling (off-diagonal, includes τ scaling)

Substituting into (*):

```
(Q*D*Q + D + 2*S) + (Q*W_eff*Q + W_eff) = 0
```

D is diagonal and W_eff is hollow, and a permutation conjugation maps each of
those two sets into itself, so the parenthesized terms lie in complementary
subspaces. Both must vanish separately.

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

**The self-decay equation, seat by seat, reads
1/τ_{Q(i)} + 1/τ_i = 1/τ_E + 1/τ_I.** A Q that sends every seat to one of the
opposite type satisfies it for any τ_E and τ_I, equal or not, with
S = ((1/τ_E + 1/τ_I)/2)·I. No constraint on the topology enters. The converse,
that only such a Q satisfies it, is true when τ_E ≠ τ_I and FALSE when they are
equal: at uniform τ both sides read 2/τ for every permutation whatsoever,
type-preserving ones included.

Two readings of that are easy and both are wrong, in opposite directions.

It does NOT need τ_E ≠ τ_I. Nothing above divides by (1/τ_E − 1/τ_I) or breaks
at τ_E = τ_I; put τ_E = τ_I = τ and D = −(1/τ)·I is scalar, so Q·D·Q = D for
EVERY permutation and S = (1/τ)·I solves it whatever Q does. Selective damping
is not what makes the diagonal half true.

And fixed-point-freeness is NOT enough to replace it, once τ_E ≠ τ_I.
Type-PRESERVING involutions have no fixed points either, and they fail there: take τ_E = 5 on seats 0,1 and
τ_I = 10 on seats 2,3 with Q = (0 1)(2 3). Then Q·D·Q = D, and the residual
Q·D·Q + D + 2S is diag(−0.1, −0.1, +0.1, +0.1), not zero. One can build an
8-seat instance where (a)-as-fixed-point-freeness and (b) both hold exactly and
the spectrum does not pair, while the trace still gives the right mean pair sum,
which is why F137 says a centre check cannot catch it.

What τ_E ≠ τ_I does is make (a) BITE. At uniform τ every permutation satisfies
the diagonal half; at τ_E ≠ τ_I only a type-swapping one does, and fixed points
are forbidden as a special case, since a fixed seat would need 2/τ_i to equal the
reciprocal sum 1/τ_E + 1/τ_I. So a uniformly damped network is not outside this
theorem on account of its damping. What it still has to supply is condition (b).

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

**Dale's Law satisfies the sign part of condition (b) ON THE SUPPORT of W.**
The qualifier is not decoration. The argument above assumes both `W[i,j]` and
`W[Q(i),Q(j)]` are nonzero; condition (b) also demands
`W[i,j] = 0 ⟹ W[Q(i),Q(j)] = 0`, that is, that the zero PATTERN be
Q-symmetric. Dale's Law says nothing about the zero pattern, and in practice
that is the binding half: on the C. elegans connectome it fails on a count, with
253 neurons whose outgoing row is non-empty and excitatory against 18 non-empty
and inhibitory, so no Q can pair them
([Neural Gamma Cavity](../../../experiments/NEURAL_GAMMA_CAVITY.md), gate G0b).

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

| Network | \|\|R\|\| / \|\|J\|\| |
|---------|--------------|
| C. elegans | 0.013 |
| Erdos-Renyi | 0.108 |
| Degree-preserving rewiring | 0.013 |

Both comparison rows were withdrawn on 2026-08-26 and the section is a NULL
result under a heading called Verification. The Erdos-Renyi row compares two arms
normalised by different constants (the connectome block globally by max|W| = 37,
the control by its own maximum), and the ratio tracks the arms' coupling
magnitude to half a percent; matched, the ratio runs 0.960 at N = 10, 0.841 at
N = 20 and 0.748 at N = 26, so parity is the smallest block's answer and the
smaller residue at the larger sizes is an open question rather than a result in
either direction. The degree-preserving row
cannot move the number at all, because its rewire keeps each weight in its own
row and the metric weights by that row. What stands is that the 8.46 was a
difference of constants; the smaller matched residue at the larger sizes is open,
and the sibling page splits it by the closed form's coverage rather than
resolving it. Details, including the closed
form and its sparsity condition, in
[Algebraic Palindrome Neural](../ALGEBRAIC_PALINDROME_NEURAL.md).

### Eigenvalue pair sums

For C. elegans subnetwork at α = 0.3, τ_E = 5, τ_I = 10:
Predicted sum: -(1/5 + 1/10) = -0.300.
Observed: mean = -0.3012, max deviation 1.6%.

**The mean pair sum is an identity and confirms nothing.** For a 5E + 5I
network with no self-coupling, trace(J) = -(5/5 + 5/10) = -1.5, and ANY
partition of the ten eigenvalues into five pairs, palindromic or arbitrary, has
mean pair sum trace(J)/5 = -0.300 exactly. The quantity cannot distinguish a
paired spectrum from an unpaired one; this is F137 (`ANALYTICAL_FORMULAS.md`),
which states the centre is trace over dimension and is equally well defined for
a spectrum that does not pair. The only informative figure in the block is the
max deviation, 1.6%. The observed mean of -0.3012 is not a reading of how well
the spectrum pairs, and cannot be: the identity fixes the mean at -0.300 exactly
for EVERY full partition of the ten eigenvalues. A mean that differs therefore
says the reported figure was not computed over a full partition, was not
computed on this J, or was not computed at the time constants printed beside it
(the committed successor script runs at τ_E = 10, τ_I = 20, where the same
identity gives -0.150). Which of the three it was is not recoverable from the
March record.

---

*See also:*
[Mirror Symmetry Proof](../../proofs/MIRROR_SYMMETRY_PROOF.md) (quantum original),
[Algebraic Palindrome Neural](../ALGEBRAIC_PALINDROME_NEURAL.md) (results and interpretation)
