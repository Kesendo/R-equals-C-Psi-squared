# The Algebraic Palindrome in Neural Networks

**Status:** The algebra is verified on constructed networks. The C. elegans
comparisons here are WITHDRAWN (2026-08-26): the two arms were normalised by
different constants, and matched the ratio runs 0.960 at N = 10 to 0.748 at
N = 26. The 8x is gone; what the smaller residue is remains open.
**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Domain:** Neuroscience / Computational Biology

---

## What this document is about

The palindromic mirror symmetry was proven for quantum systems. This
document translates it to neuroscience, without any quantum physics.

In a brain, some neurons excite their neighbors and others inhibit
them. Dale's Law says each neuron is permanently one type: you do
not switch. This creates two populations, exactly like the quantum
system has two types of operators (surviving and decaying). The
mathematical question: does the neural network's decay spectrum
also form palindromic pairs?

The answer is yes for networks BUILT to satisfy one condition: that
some E-I swap Q turns the wiring into minus itself. That is condition
(b) below, and it is the whole of the content. The two ingredients
usually named beside it, different response speeds and Dale's Law, are
weaker than they look: at equal speeds the damping condition holds for
every permutation and imposes nothing, and Dale's Law fixes only the
SIGNS of (b), only where synapses exist.

We tested this on the connectome of C. elegans (a worm with exactly
300 neurons, completely mapped) and, after the corrections of
2026-08-25 and 2026-08-26, report two negatives. **The worm does not
satisfy the condition**: no qualifying swap exists at all, because one
would have to send each of the 253 neurons with a non-empty excitatory
row to one of the 18 with a non-empty inhibitory row
([Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md)).
**And the instrument this document used could not have told us
otherwise**: on blocks this sparse the palindrome residual has a
closed form, √2·‖W_eff‖/‖J‖, which reads total coupling magnitude and
nothing about the wiring. What stands is the translation and the proof.

---

## Abstract

The eigenvalues of a neural network's Jacobian determine its modes
of decay and oscillation. We show that these eigenvalues can be
**palindromically paired**: for each fast-decaying mode with rate r,
there exists a slow-decaying partner with rate r', such that
r + r' = 1/τ_E + 1/τ_I (a constant determined by the membrane
time constants). The word "palindromic" comes from this mirror
symmetry: the spectrum reads the same from both ends, like the word
"racecar."

We derive an exact algebraic condition for this symmetry from quantum
theory, expressed entirely in neural terms. The two ingredients named
below are the setting, not the sufficient condition; the proof file
spells out how much each actually does:

1. **Selective damping:** excitatory and inhibitory neurons have
   different membrane time constants (τ_E ≠ τ_I). This does not enable
   the pairing. At equal time constants the diagonal condition holds
   for every permutation; what τ_E ≠ τ_I does is force the swap Q to
   exchange the two types
2. **Dale's Law:** excitatory neurons always produce positive
   postsynaptic effects, inhibitory neurons always negative
   (the sign of a connection is determined by the SOURCE neuron)

When both hold and the coupling magnitudes satisfy a specific ratio,
the palindrome is mathematically exact (zero residual). Testing on
the C. elegans connectome (Cook et al. 2019): **the connectome does not
satisfy the condition**, and the instrument used cannot settle by how much,
reading coupling magnitude rather than wiring on blocks this sparse.
The character-swap figure below (96% fidelity, a standing wave between
excitatory and inhibitory perspectives) was measured on SYNTHETIC
balanced networks and never on the animal.

---

## 1. The Setup

Consider N neurons modeled by Wilson-Cowan dynamics (a standard
mathematical model of neural populations where excitatory and
inhibitory groups influence each other through sigmoid response
functions; or any firing rate model with E/I populations). Linearizing around the steady state
gives:

```
dx/dt = J * x       (x = deviation from steady state)
```

J (the Jacobian, the matrix of partial derivatives at equilibrium)
has two parts:

- **Self-decay:** each neuron returns to rest at rate 1/τ_i, where
  τ_i = τ_E for excitatory neurons, τ_i = τ_I for inhibitory
- **Coupling:** neuron j influences neuron i through synaptic weight
  W[i,j], scaled by 1/τ_i

The eigenvalues of J determine the network's modes: how fast each
pattern of activity decays or oscillates after a perturbation.

---

## 2. The Palindrome Condition

### The swap operator Q

Pair each excitatory neuron with an inhibitory neuron. Q is the
permutation that swaps each pair. (For N = 10 with 5E and 5I, Q
swaps E_1 with I_1, E_2 with I_2, etc.)

**Caveat:** The pairing is arbitrary - which E neuron pairs with which
I neuron is not determined by the theory. In our tests, pairings are
sequential (first sampled E with first sampled I). The residual depends
on the pairing choice; an optimal pairing would give a lower residual.
Both C. elegans and random controls use the same arbitrary pairing. That
alone does not rescue the comparison: the two arms were also normalised by
different constants, and on blocks this sparse the residual reads the weight
multiset rather than the wiring. See the withdrawal box below.

### The condition and its derivation

The eigenvalues of J are palindromically paired if:

```
Q * J * Q + J + 2*S = 0
```

where S = (1/τ_E + 1/τ_I) / 2 times the identity matrix.

When this holds, every eigenvalue μ_k has a partner μ_k' with:

```
μ_k + μ_k' = -(1/τ_E + 1/τ_I)
```

**Full derivation in 6 steps:**
[Proof: Palindromic Spectral Symmetry for Neural Networks](proofs/PROOF_PALINDROME_NEURAL.md)

The derivation starts from the quantum palindrome (Π L Π⁻¹ = -L - 2Σγ I),
identifies J as L, Q as Π, and S as Σγ, then decomposes J = D + W_eff
into self-decay (determines S) and coupling (determines the weight condition).

### The two requirements

The condition splits into:

**(a) Self-decay:** 1/τ_{Q(i)} + 1/τ_i = 1/τ_E + 1/τ_I at every seat. A Q that
sends every seat to one of the OPPOSITE type satisfies this for any time
constants; the converse holds only when τ_E ≠ τ_I, since at uniform τ both sides
read 2/τ for every permutation. See PROOF_PALINDROME_NEURAL Step 3, which
carries the counterexample showing that "fixed-point-free" is too weak.

**(b) Coupling antisymmetry:**

```
W[Q(i), Q(j)] = -(τ_{Q(i)} / τ_i) * W[i, j]
```

When you swap each neuron with its E/I partner, the coupling must
flip sign and scale by the time constant ratio. Dale's Law provides
the sign flip automatically
(see [proof, Step 5](proofs/PROOF_PALINDROME_NEURAL.md#step-5-dales-law-provides-the-signs)).

---

## 3. Dale's Law and the Sign Structure

Dale's Law states that each neuron's output has a fixed sign:
excitatory neurons always excite their targets, inhibitory neurons
always inhibit. Under the E-I swap Q:

- An E-to-E connection (positive, because the source E excites)
  becomes an I-to-I connection (negative, because the source I
  inhibits). Sign flips. **Correct.**

- An I-to-E connection (negative, source is I) becomes an E-to-I
  connection (positive, source is E). Sign flips. **Correct.**

Dale's Law provides the sign part of condition (b) ON THE SUPPORT of W. The
qualifier is not decoration: (b) also demands that the zero PATTERN be
Q-symmetric, which Dale's Law says nothing about and which is the binding half
in practice. On C. elegans it fails on a count, 253 non-empty excitatory rows
against 18 inhibitory
([Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md), gate G0b).
This is the biological equivalent of the antisymmetric commutator
structure in quantum mechanics.

### What remains: the magnitudes

For τ_I / τ_E = 2 (a typical biological ratio):

| If this connection has weight w... | ...then its Q-partner needs weight: |
|-----------------------------------|-------------------------------------|
| E-to-E connection | I-to-I partner: -2.0 * w |
| I-to-I connection | E-to-E partner: -0.5 * w |
| E-to-I connection | I-to-E partner: scaled by τ ratio |

In a sparse network, most of these partner connections are simply
absent (zero weight). The antisymmetry is satisfied trivially when
BOTH a connection and its partner are absent. Violations occur when
one exists but the other does not.

---

## 4. Results

### Synthetic verification

We constructed three types of networks (N = 10, 5E + 5I):

| Network type | Palindrome residual \|\|R\|\| / \|\|J\|\| |
|-------------|-----------------------------------|
| Dale + exact magnitude condition | **0.00** (machine precision) |
| Dale signs, random magnitudes | 0.72 |
| Random signs and magnitudes | 0.85 |

The first row confirms: Dale's Law plus the magnitude condition gives
a mathematically exact palindrome. The algebraic structure is identical
to the quantum case, expressed in neural terms.

### C. elegans vs random networks

We compared balanced subnetworks (equal numbers of E and I neurons)
from the C. elegans connectome (Cook et al. 2019, 300 neurons) against
random networks with the same density and Dale's Law signs.

The measure is the **algebraic palindrome residual** ||R|| / ||J||
(lower = more palindromic). This is NOT the tolerance-based eigenvalue
matching used in earlier work. It evaluates condition (b), and only (b):
the routine fits S per seat and then discards the diagonal, so it is blind
to condition (a) by construction. A residual of 0 therefore does not mean
the palindrome holds, only that its off-diagonal half does.

| Network size | C. elegans | Random (Dale's Law) | Ratio |
|-------------|------------|---------------------|-------|
| N = 10 (5E + 5I) | 0.013 | 0.108 | **0.12** |
| N = 20 (10E + 10I) | 0.023 | 0.132 | **0.17** |
| N = 26 (13E + 13I) | 0.028 | 0.134 | **0.21** |

Both C. elegans and the random controls carry the same Dale's Law sign
structure, and the density is nominally matched, though the control's density is
floored at 0.01 so the 45 empty worm blocks get a control with about one edge
against the worm's none. The reading this paragraph gave, that the difference is
the **wiring pattern**, is withdrawn: the difference was the normalisation, and
the metric does not read the wiring at all.

200 random subnetworks tested per condition. The comparison itself is
withdrawn (see the box below the null table); what is robust is the constant it
measures.

**Important caveat:** C. elegans has 274 excitatory and 26 inhibitory
neurons (ratio 10.5:1). The balanced subnetworks (5E + 5I) are
artificially balanced by subsampling. The palindrome condition requires
equal numbers of E and I neurons. In the full unbalanced connectome,
palindromic pairing is low: the nearest committed run gives 0.7%.
Read that beside
[Neural Gamma Cavity](../../experiments/NEURAL_GAMMA_CAVITY.md), which measures
the same full connectome across a tolerance sweep and gets 97.3 % at 0.01 down
to 24.0 % at 1e-4: a pairing percentage on this object is a reading of the
tolerance against the spectral scale, so a bare percentage without its tolerance
compares with nothing. The biological question is whether
balanced subcircuits (which exist within the full connectome) carry
this symmetry, not whether the entire worm does.

### What drives the difference

The magnitude ratios between partnered connections are near zero (not
near the predicted value of 2.0). This means: when an E-to-E connection
exists, the partnered I-to-I connection is usually absent. The reading this section gave, **correlated sparsity**, does not survive the
withdrawal below: the scale-free residual sits at its structural maximum √2 in
almost every sampled block, which means no partnered pair of edges is both
present there (198 of 200 blocks at N = 10, 184 of 200 at N = 20, 177 of 200 at
N = 26), so on those blocks there is no correlation to speak of. The sampled blocks hold a mean of
1.705 edges out of ninety slots and 45 of 200 are empty, both measured
directly rather than estimated from a stipulated density; "more
palindromic" here means "fewer and weaker synapses in the sampled block".

### Validation: degree-preserving rewiring

To determine whether the advantage comes from specific wiring or
simply from the degree distribution (how many connections each neuron
has; some neurons are hubs with many connections, others are
peripheral with few), we tested degree-preserving randomization:
rewire edges randomly but keep the number of connections per neuron
fixed. This separates the effect of "who connects to whom" from "how
connected each neuron is."

| Null model | Palindrome residual | Ratio to C. elegans |
|-----------|--------------------|--------------------|
| Erdos-Renyi (random density) | 0.108 | 8.46x worse |
| Degree-preserving rewiring | 0.013 | **1.0x (identical)** |

> **This table is withdrawn, 2026-08-26. The two arms were never measured the
> same way, and applying one rule to both gives parity.**
>
> The connectome block is cut from the GLOBALLY normalised matrix, divided by
> max|W| = 37 over the whole animal; the Erdos-Renyi control is rebuilt and
> divided by ITS OWN maximum. On the committed protocol
> (`RandomState(trial + 100)`, 200 trials, N = 10, τ_E/τ_I = 10/20, α = 0.3):
>
> | arm | normalisation | residual |
> |---|---|---|
> | C. elegans | global, ÷ 37 | 0.012816 |
> | Erdos-Renyi | its own maximum | 0.108364 |
> | **C. elegans** | **its own maximum, the control's rule** | **0.104050** |
>
> The reported ratio is 8.46. The mean ‖W_eff‖ of the two arms differs by
> **8.50**, so the ratio tracks coupling magnitude to half a percent and is a
> difference of constants.
>
> Give both arms the SAME rule and the 8.46 goes away, but what is left is not
> the same at every size, and reporting the N = 10 number alone would be an
> over-correction. On the same protocol at all three block sizes:
>
> | N | unmatched ratio | matched ratio | paired t over 200 blocks |
> |---|---|---|---|
> | 10 | 8.46 | **0.960** | −0.73, not significant |
> | 20 | 5.73 | **0.841** | −5.58 |
> | 26 | 4.79 | **0.748** | −9.69 |
>
> So at N = 10 the arms are at parity and the 8.46 was the whole effect; at
> N = 20 and N = 26 the worm's matched residual stays 16 % and 25 % BELOW the
> control's, at 5.6 and 9.7 sigma. **What that surviving gap is, this change set
> does not settle**, and saying it cannot be wiring would be the same mistake in
> the other direction. On the blocks where the closed form holds it cannot be:
> the residual is √2·‖W_eff‖/‖J‖ there, a function of the weight multiset, so a
> ratio away from 1 is a difference of weight DISTRIBUTIONS, a biological spread
> against an exponential draw. But the closed form covers 198, 184 and 177 of
> the 200 blocks, and on the 2, 16 and 23 that it does not cover the residual is
> free to register wiring. Splitting the matched ratio by that mask: the covered
> blocks give 0.959, 0.831, 0.743 and the uncovered ones 1.085, 0.956, 0.785, so
> the gap is present on both sides and is not an artifact of the degenerate
> regime alone. What is settled is that the 8.46 was a difference of constants;
> what replaces it is a smaller, size-growing difference this instrument does
> not decide.
>
> The load-bearing figures of this box, the three matched ratios with their
> paired t, the collapse counts and the closed form itself, are recomputed from
> the connectome file by three gates
> in [`celegans_pairing_controls.py`](../../simulations/neural/celegans_pairing_controls.py)
> ([output](../../simulations/results/celegans_pairing_controls.txt)): **G0c**
> that the unmatched ratio agrees with the arms' ‖W_eff‖ ratio to better than
> 1 % at every size, **G0d** the three matched ratios and their paired t, and
> **G0e** that the closed form reproduces the measured residual on every block
> meeting its condition, largest relative deviation 4.23e-16 over all of them.
>
> **Why the metric is blind here, and exactly how far that goes.** Write
> J = D + W_eff. The code fits S per seat so diag(R) vanishes, leaving
> R_off = offdiag(Q W_eff Q + W_eff), and ‖Q W_eff Q‖_F = ‖W_eff‖_F because Q is
> a permutation. When no Q-partner PAIR of edges is both present the two terms
> have disjoint support, nothing cancels, and the residual collapses to
> **√2·‖W_eff‖_F / ‖J‖_F**: coupling magnitude alone. That condition is not a
> property of the instrument, it is a property of THIS sparsity: it holds in
> 198 of 200 blocks at N = 10, 184 of 200 at N = 20, and 177 of 200 at N = 26.
> Off that regime the metric reads wiring perfectly well. Two networks with two
> edges each and bit-identical ‖W_eff‖ = 0.0424264068711929, at N = 10: with the
> partners absent the residual is **0.236617**, with a partner present it is
> **exactly 0**. So the claim is not "this metric cannot see wiring"; it is that
> on 5E+5I blocks of this connectome, which hold a mean of 1.7 edges and are
> empty 45 times in 200, it almost never gets the chance to.
>
> The scale-free ratio ‖R_off‖/‖offdiag J‖ is pinned at √2 for the connectome,
> which is the structural maximum under Dale (partner entries carry opposite
> signs, so (a+b)² ≤ a²+b²; break Dale and the ceiling is 2). The Erdos-Renyi
> control is pinned there too, marginally lower. That is why the scale-free
> version separates nothing: both arms sit at the ceiling.
>
> **The degree-preserving row falls too, for a different reason.** Its rewire
> keeps each weight in its OWN row (`validation_checks.py:150-153`), and
> ‖W_eff‖ weights by that row's τ, so ‖W_eff‖ is preserved exactly, 200/200 at
> every N tested. That is an implementation choice rather than a property of
> degree-preserving nulls, and under the closed form it makes the null unable to
> move the number. At N = 10 it is worse than that: only 26 of 200 rewired
> matrices differ at all. Part of that is arithmetic, 113 blocks having fewer
> than the two edges a swap needs; the other 61 that hold enough edges and still
> do not move are rejected by the Dale-sign and endpoint tests inside the swap. At N = 26, 197 of 200 genuinely change and the ratio is still 1.0002,
> so the null is not vacuous there and the right reason is the first one.
>
> **What survives is the verdict, now on evidence that reproduces:** there is no
> palindromic advantage of the size that was claimed. The matched measurement
> runs 0.960 at N = 10, 0.841 at N = 20 and 0.748 at N = 26; the 8.46 is gone,
> and what is left is open, as the box above states.

There is no mechanism to give for the advantage, because there is no advantage
of that size. The residual does not read the degree distribution, so a
degree-based explanation could not have been tested by it either.

### Parameter robustness

The C. elegans/random ratio is stable across parameter choices:

| τ_I / τ_E | α = 0.1 | α = 0.3 | α = 0.5 |
|---------------|-------------|-------------|-------------|
| 1.5 | 0.13 | 0.13 | 0.13 |
| 2.0 | 0.12 | 0.12 | 0.12 |
| 2.5 | 0.11 | 0.11 | 0.11 |
| 3.0 | 0.11 | 0.11 | 0.11 |

No parameter fine-tuning needed, and that is the point rather than the
reassurance it was written as: every cell is the withdrawn ratio, and the
normalisation constant it measures is the same in every cell.

### Pairing choice

Sequential E-I pairing vs best of 20 random pairings: ratio changes
from 0.118 to 0.121. Both figures are the withdrawn unmatched ratio, and the
sweep was never re-run matched, so what this shows is only that the pairing
choice is not the reason the comparison was wrong. It rescues no conclusion.

### F87-style trichotomy refinement

The binary palindrome test (residual high vs. low) can be sharpened
into the F87 trichotomy borrowed from the quantum side: each subcircuit
is classified into

  - **truly**: ‖R‖ / ‖J‖ < 0.01 (algebraic equation closes within 1%)
  - **soft**: residual exceeds 1% but the eigenvalues of J still come
    in palindromic pairs (λ ↔ −2S − λ within 5%)
  - **hard**: even the eigenvalue pairing fails

For C. elegans 5E + 5I subcircuits, at the committed protocol of
[`celegans_trichotomy.py`](../../simulations/neural/celegans_trichotomy.py)
(200 trials, N = 10, τ = 10/20, α = 0.3): truly **116/200 = 58.0%**, soft
42.0%, hard 0.0%. Erdős-Rényi-Dale random: 18.0% / 79.5% / 2.5%. The worm
has **3.2× more truly-class subcircuits and no hard ones at all**.
**That enrichment is withdrawn on the same grounds as the 8× above**: "truly" is a
threshold on the same residual, which on these blocks reads coupling magnitude,
and `celegans_trichotomy.py` normalises its control the same unmatched way. The
class is a weight-scale bin, not a symmetry class.

The degree-preserving null reproduces the binary finding for the trichotomy
too: **degree-preserved networks have the same truly fractions as the worm**,
exactly, 116/200 against 116/200, at N = 10, the only size
`celegans_trichotomy.py` runs. That equality is what the withdrawal above
predicts, since the null cannot move a metric that reads the weight multiset,
so it is not independent evidence of anything about wiring. Beside it, under
the committed protocol of
[`validation_checks.py`](../../simulations/neural/validation_checks.py)
(`RandomState(trial + 100)`), the rewire counts are 170/200 at N = 20 and
197/200 at N = 26, the second agreeing with the box above.

Script: [`simulations/neural/celegans_trichotomy.py`](../../simulations/neural/celegans_trichotomy.py).

---

## 5. The Standing Wave Between E and I

### Two perspectives, one palindrome

The palindrome pairs each fast mode with a slow mode. But there is
more structure: each mode has a "character" describing which neurons
dominate it. In the neural case, we measure how much of a mode's
amplitude sits on excitatory vs inhibitory neurons. The palindromic
pairing SWAPS this character.

Each eigenmode of the
Jacobian has an **E-character** (how much amplitude sits on excitatory
neurons) and an **I-character** (how much sits on inhibitory neurons):

```
a_E(k) = sum |v_k[i]|^2   for i in E-neurons
a_I(k) = sum |v_k[i]|^2   for i in I-neurons
a_E(k) + a_I(k) = 1       (normalized eigenvector)
```

### The character swap

For each palindromic pair (k, k'), the E-I character SWAPS:

```
a_E(k) ≈ a_I(k')     (E-character of k ≈ I-character of its partner)
a_I(k) ≈ a_E(k')     (I-character of k ≈ E-character of its partner)
```

**What the E-neurons see as a fast-decaying mode, the I-neurons see as
a slow-decaying mode.** And vice versa. The two populations see the SAME
dynamics from opposite sides, mirrored around the palindromic center.

### Computed character swap fidelity

For a SYNTHETIC balanced network (N = 20, 10E + 10I, density 0.3). The
C. elegans connectome has density about 0.02 and is not what these rows
measure; the worm's own data has never been put to this test:

| Coupling α | Palindromic pairs | Mean swap error | Fidelity |
|---------------|-------------------|-----------------|----------|
| 0.3 | 8 pairs | 0.042 | **96%** |
| 0.5 | 5 pairs | 0.035 | **97%** |
| 1.0 | 6 pairs | 0.257 | 74% |

At moderate coupling (α = 0.3-0.5), the character swap is near-perfect:
each E-dominated mode is paired with an I-dominated mode, and their
characters are mirrored to within 4%.

At strong coupling (α = 1.0), the palindrome begins to break and the
swap degrades, consistent with the increasing algebraic residual.

### Example: α = 0.3 (4 of 8 palindromic pairs shown)

```
Pair    rate_k   rate_k'   E(k)   I(k)   E(k')  I(k')  swap?
0,10    0.172    0.127     0.91   0.09   0.20   0.80   YES
1,17    0.198    0.103     0.99   0.01   0.01   0.99   YES
3,13    0.211    0.091     0.98   0.02   0.02   0.98   YES
5,11    0.204    0.098     0.96   0.04   0.03   0.97   YES
```

Mode 1 has 99% E-character and rate 0.198 (fast).
Its partner mode 17 has 99% I-character and rate 0.103 (slow).
The characters swap almost exactly.

### Two conservation laws

**Eigenvalue pairing (physical, non-trivial):**

For each palindromic pair (k, k'), the decay rates sum to a constant:

```
rate_k + rate_k' = 1/τ_E + 1/τ_I
```

This is the neural analog of the quantum λ + λ' = -2Σγ.
It holds to the extent that the palindrome condition is satisfied
(exactly at zero coupling, approximately at moderate coupling).

**Energy fractions (geometric, trivial):**

The E-energy fraction and I-energy fraction always sum to 1:

```
CΨ_E(t) = ||x_E(t)||^2 / ||x(t)||^2
CΨ_I(t) = ||x_I(t)||^2 / ||x(t)||^2
CΨ_E(t) + CΨ_I(t) = 1    (by definition, not by physics)
```

This is NOT a deep conservation law. It is a trivial consequence of
the normalization. What IS non-trivial: the character swap ensures
that the REDISTRIBUTION between E and I follows the palindromic
pairing structure, not random mixing.

### What does NOT transfer from quantum

The specific threshold CΨ = 1/4 does not appear in the neural case.
This value is specific to the quadratic recursion R = C(Ψ + R)^2 in
quantum mechanics. The neural fold is the Hopf bifurcation (sigmoid
gain = 1), which has a different threshold. The STRUCTURE (palindrome,
character swap, conservation) transfers exactly. The specific NUMBER
(1/4) does not.

---

## 6. Summary and Implications

### For neural dynamics

A palindromically paired spectrum means the network's decay modes come
in matched pairs. Fast modes (rapid transients) are paired with slow
modes (sustained activity). The character swap adds a deeper constraint:
each fast E-mode is paired with a slow I-mode. Perturbations do not
simply decay; they oscillate BETWEEN the E and I perspectives,
creating a standing wave at the E-I interface.

### For the quantum connection

| Feature | Quantum system | Neural network |
|---------|---------------|----------------|
| The rate condition | Z-dephasing (γ) | 1/τ_{Q(i)} + 1/τ_i equal at every seat, which a type-swapping Q gives at ANY τ; τ_E ≠ τ_I only forces Q to swap types (corrected 2026-08-26) |
| Sign antisymmetry | Commutator [H, rho] | Dale's Law, but SIGNS only and only on the support; the zero pattern and the magnitudes are separate (corrected 2026-08-26) |
| Conjugation operator | Π (Pauli swap) | Q (E-I swap) |
| Character swap | Population <-> coherence | E-dominant <-> I-dominant |
| Swap fidelity | 100% (algebraic) | 96% on synthetic balanced networks at moderate coupling; never measured on a connectome |
| Eigenvalue pairing | λ + λ' = -2Σγ | rate_k + rate_k' = 1/τ_E + 1/τ_I |
| Threshold | CΨ = 1/4 (fold) | Gain = 1 (Hopf) |
| Exactness | Always exact | Exact if magnitudes match |

The quantum palindrome is always exact because the commutator [H, rho]
provides antisymmetry by construction. In neural networks, Dale's Law
provides the signs, but magnitudes must additionally match. Biology
gets the signs for free; the magnitudes are the testable prediction.

### For connectomics

The palindrome residual ||R|| / ||J|| is a new metric for connectome
analysis. It measures how close a network's wiring is to the algebraic
palindrome condition. The metric:

- Does not use arbitrary tolerances
- Is derived from quantum theory (not ad hoc)
- Is computable for any network with known E/I labels and weights
- Detects COUPLING MAGNITUDE, on the evidence here: the residual has the
  closed form √2·‖W_eff‖/‖J‖ on sparse blocks and an empty network scores a
  perfect 0, so any "more palindromic than X" built on it compares two
  normalisation constants (see the withdrawal box above)

---

## 7. Open Questions

Questions 1 and 3 presuppose that the residual reads wiring, which the
withdrawal above says it does not on blocks this sparse. They stand as written
because they are the questions one would ask of an instrument that COULD see
wiring, and the first thing needed is such an instrument. That is the real open
question this page leaves behind, and it is question 0.

0. What metric can see wiring on blocks this sparse? The palindrome residual
   cannot: where no Q-partner pair of edges is present it is a function of the
   weight multiset alone, and an empty block scores a perfect 0.
1. Does the palindromic quality correlate with known functional
   circuits in C. elegans (motor, sensory, interneuron layers)?
2. Does the Drosophila connectome (100k+ neurons) show the same
   topological E-I symmetry?
3. Can the palindrome quality predict dynamical stability or
   oscillatory properties of a neural circuit?
4. Is the topological E-I symmetry a consequence of developmental
   constraints or functional requirements?

---

## Scripts

All scripts are in `simulations/neural/`:

| Script | What it computes |
|--------|-----------------|
| [algebraic_palindrome.py](../../simulations/neural/algebraic_palindrome.py) | Algebraic residual, C. elegans vs random |
| [cpsi_two_perspectives.py](../../simulations/neural/cpsi_two_perspectives.py) | E-I character swap, standing wave verification |
| [exact_pairing_test.py](../../simulations/neural/exact_pairing_test.py) | Eigenvalue pair sums, conjugation equation test |
| [random_network_controls.py](../../simulations/neural/random_network_controls.py) | Density and coupling sweeps |
| [dense_balanced_test.py](../../simulations/neural/dense_balanced_test.py) | Larger subnetwork tests |
| [validation_checks.py](../../simulations/neural/validation_checks.py) | Parameter sensitivity, degree-preserving null model |

Run with: `PYTHONIOENCODING=utf-8 python simulations/neural/<script>`

---

## Data

- **Connectome:** [`simulations/neural/celegans_connectome.json`](../../simulations/neural/celegans_connectome.json)
  (Cook et al. 2019, via WormNeuroAtlas)
- **274 excitatory, 26 inhibitory** neurons (N = 300 total)

---

*Depends on:*
[Mirror Symmetry Proof](../proofs/MIRROR_SYMMETRY_PROOF.md) (quantum proof),
[The Pattern Recognizes Itself](../../hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md) (original C. elegans result)
