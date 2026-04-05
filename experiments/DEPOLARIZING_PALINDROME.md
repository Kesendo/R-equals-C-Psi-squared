# Why Depolarizing Noise Breaks the Palindrome: The 2:2 vs 1:3 Split

<!-- Keywords: depolarizing noise palindrome breaking, dephasing spectral symmetry
condition, Pauli rate pairing conjugation operator, immune decaying sector balance,
two axis dephasing palindrome, per-site operator split qubit, Lindblad noise type
spectral effect, palindromic symmetry necessary condition, depolarizing error formula
linear N, open quantum system noise classification, R=CPsi2 depolarizing palindrome -->

**Status:** Proven (theorem + numerical verification)
**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Π as Time Reversal](PI_AS_TIME_REVERSAL.md)

---

## What this document is about

The palindromic mirror works perfectly when noise attacks along one
or two axes. But what if the noise comes from every direction at once?
This is called depolarizing noise, and it breaks the palindrome.

The reason is surprisingly simple. The mirror Π works by swapping two
halves of the quantum world: the part that survives noise and the part
that decays. Under normal dephasing, these two halves are equal in
size (two surviving operators, two decaying ones at each qubit). Under
depolarizing noise, only one operator survives and three decay: 1 versus 3.
You cannot build a mirror between one thing and three things, for the
same reason you cannot pair up five people for a dance when there are
only two chairs.

This is not just a technical limitation. In the time-reversal language
of [Π as Time Reversal](PI_AS_TIME_REVERSAL.md), it means the future
is exponentially larger than the past. There are overwhelmingly more
ways to be undecided than decided.

---

## Abstract

The palindromic Liouvillian symmetry holds for Z-dephasing, X-dephasing,
Y-dephasing, and all two-axis combinations. It breaks for depolarizing noise
(all three axes simultaneously). The reason is a per-site counting argument:
Z-dephasing splits the four Pauli indices {I, X, Y, Z} into 2 immune and
2 decaying (a balanced 2:2 split). Depolarizing noise splits them 1:3 (only I
is immune). A bijective mirror Π requires equal numbers of immune and decaying
indices at each site. With a 1:3 split, no such bijection exists. The palindrome
error under depolarizing noise is exactly (2/3)Nγ = (2/3)Σγ, Hamiltonian-independent
and linear in both γ and N. For typical hardware dephasing rates (γ ∼ 0.001),
this is < 0.1%. The general condition: the palindrome holds if and only if at
least one of γ_X, γ_Y, γ_Z is zero (at most two dephasing axes).

---

## The Question

The palindromic symmetry holds for Z-dephasing. It holds for X-dephasing
and Y-dephasing. It even holds for two-axis dephasing (Z+X, Z+Y, X+Y).
But it breaks for depolarizing noise, which is all three axes at once.

Why? What is special about three axes versus two?

---

## 1. The Answer in One Sentence

The palindrome requires a bijective mirror between "immune" and "decaying"
Pauli indices at each qubit site. Single-axis dephasing splits the four
indices {I, X, Y, Z} into two immune and two decaying (a 2:2 split).
Depolarizing noise splits them 1:3. You cannot build a bijection between
1 and 3. The mirror does not fit.

---

## 2. The Per-Site Rate Table

Each type of noise assigns a dephasing rate to each Pauli index at a
single qubit site. The rate depends on whether the Pauli operator commutes
or anti-commutes with the jump operators (the mathematical objects
that describe how the environment disturbs each qubit):

| Pauli | Z-deph | X-deph | Y-deph | Depol | Z+X | Z+Y | X+Y |
|---|---|---|---|---|---|---|---|
| I | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| X | 2γ | 0 | 2γ | 4γ/3 | γ | 2γ | γ |
| Y | 2γ | 2γ | 0 | 4γ/3 | 2γ | γ | γ |
| Z | 0 | 2γ | 2γ | 4γ/3 | γ | γ | 2γ |

The identity I always commutes with everything. It is always immune.
Under single-axis dephasing, one additional Pauli is immune (the one
that commutes with the jump operator). Under depolarizing noise, no
additional Pauli is immune. X, Y, and Z all anti-commute with at least
two of the three jump operators.

The split:

| Noise | Immune | Decaying | Split |
|---|---|---|---|
| Z-dephasing | {I, Z} | {X, Y} | 2:2 |
| X-dephasing | {I, X} | {Y, Z} | 2:2 |
| Y-dephasing | {I, Y} | {X, Z} | 2:2 |
| Depolarizing | {I} | {X, Y, Z} | 1:3 |
| Z+X | {I} | {X, Y, Z} | 1:3 |
| Z+Y | {I} | {X, Y, Z} | 1:3 |
| X+Y | {I} | {X, Y, Z} | 1:3 |

Two-axis dephasing also has a 1:3 split, yet the palindrome survives
for Z+X, Z+Y, and X+Y. The split alone is not the full story.

---

## 3. What the Mirror Actually Needs

The conjugation operator Π satisfies Π L_D Π^{-1} = -L_D - c I for
some constant c. This means: Π must permute the four Pauli indices such
that every rate r maps to c-r, where c is the same for all four.

This is a rate-pairing condition: the four per-site rates [r_I, r_X, r_Y, r_Z]
must split into two pairs that each sum to c.

For Z-dephasing the rates are [0, 2γ, 2γ, 0]. The pairs are (0, 2γ) and
(0, 2γ), both summing to c = 2γ. Four valid permutations exist:

```
I->X, X->I, Y->Z, Z->Y  (our known Pi)
I->X, X->Z, Y->I, Z->Y
I->Y, X->I, Y->Z, Z->X
I->Y, X->Z, Y->I, Z->X
```

For depolarizing noise the rates are [0, 4γ/3, 4γ/3, 4γ/3]. The only way
to pair the rate-0 index (I) is with some c, requiring the other three
indices to all have rate c - 4γ/3. But the only index with rate 0 is I
itself. If I maps to I, then c = 0 and all other rates would need to be
negative. If I maps to X, then c = 4γ/3 and X, Y, Z would all need rate 0.
Neither works. No valid permutation exists. Proof by exhaustion over all
24 permutations of four elements confirms: zero solutions.

For two-axis dephasing (e.g., Z+X), the rates are [0, γ, 2γ, γ]. The
pairs are (0, 2γ) and (γ, γ), both summing to c = 2γ. Two valid
permutations exist. The palindrome survives because the rates, while
having a 1:3 immune/decaying split, still admit a bijective pairing.

**The general condition:** a palindromic Π_D exists if and only if the
four per-site rates can be partitioned into two pairs with equal sums.

For noise with rates (γ_X, γ_Y, γ_Z) along three axes, the per-site
Pauli rates are [0, 2(γ_Y+γ_Z), 2(γ_X+γ_Z), 2(γ_X+γ_Y)]. The pairing
exists if and only if at least one of γ_X, γ_Y, γ_Z equals zero.
Equivalently: dephasing along at most two axes.

---

## 4. How Badly Does It Break?

Not all palindrome failures are equal. For N=3 Heisenberg under depolarizing
noise, the palindrome error at center Sγ is exactly 0.1000. For N=4, exactly
0.1333. The pattern:

```
error = (2/3) Sγ = (2/3) N γ
```

This is exact and Hamiltonian-independent. It comes from a simple gap: the
steady state at rate 0 needs a palindromic partner at rate 2Sγ. But under
depolarizing noise, the maximum achievable rate is (4/3)Sγ (where every
site carries a decaying Pauli). The gap between what the palindrome demands
(2Sγ) and what the noise can provide ((4/3)Sγ) is exactly (2/3)Sγ.

The mirror is too short. It cannot reach the far end.

---

## 5. The Interpolation: No Threshold

Define a mixed noise channel: (1-α) Z-dephasing + α depolarizing.
At α = 0, the palindrome is exact. At α = 1, the error is (2/3)Sγ.

The transition is perfectly linear: error = α (2/3) Sγ. There is no
critical threshold, no phase transition, no gradual softening. At any
α > 0, no matter how small, the palindrome breaks immediately. The Z
index acquires a nonzero rate (4γα/3), destroying its membership in the
immune set. The 2:2 split becomes 1:3 instantly.

This is a topological fact, not a quantitative one. The split is either
2:2 or it is not. There is no "almost 2:2."

---

## 6. The Counting Argument: Why the Future is Bigger Than the Past

The deepest way to see the obstruction is through counting. The idea:
at each qubit site, you choose either an "immune" or a "decaying"
Pauli operator. The total number of strings of each type depends on
how many immune vs decaying choices you have per site. If the numbers
are equal (2 and 2), the sectors balance perfectly. If not (1 and 3),
the imbalance grows exponentially with the number of qubits.

Under Z-dephasing, XY-weight w (how many sites carry X or Y) classifies
the Pauli strings. Weight w has C(N,w) 2^w 2^(N-w) = C(N,w) 2^N strings.
Its palindromic partner at weight N-w has C(N,N-w) 2^N = C(N,w) 2^N.
The sectors are always the same size. The factor 2^N comes from having
2 immune and 2 decaying choices per site: 2^w decaying choices times
2^(N-w) immune choices equals 2^N regardless of w.

Under depolarizing noise, the "weight" is the non-I count (sites with
X, Y, or Z). Weight w has C(N,w) 3^w 1^(N-w) = C(N,w) 3^w strings.
Its would-be partner at weight N-w has C(N,w) 3^(N-w). The ratio
between sectors is 3^(N-2w), which is 1 only when w = N/2.

For N=3:

| Weight | Z-deph count | Depol count | Partner weight | Partner count | Ratio |
|---|---|---|---|---|---|
| 0 | 8 | 1 | 3 | 8 / 27 | 1.0 / 27.0 |
| 1 | 24 | 9 | 2 | 24 / 27 | 1.0 / 3.0 |
| 2 | 24 | 27 | 1 | 24 / 9 | 3.0 / 1.0 |
| 3 | 8 | 27 | 0 | 8 / 1 | 27.0 / 1.0 |

Under Z-dephasing, every pair of partner sectors has equal count. Under
depolarizing noise, the asymmetry is exponential.

This is the root cause. The 2:2 per-site split gives 2 immune choices
and 2 decaying choices, so the immune and decaying sectors always
balance. The 1:3 split gives 1 immune choice and 3 decaying choices,
making the decaying sector exponentially larger.

---

## 7. What This Means for Time Reversal

[Π as Time Reversal](PI_AS_TIME_REVERSAL.md) established that Π is a
time-reversal operator: it swaps populations (past, classical, persistent)
with coherences (future, quantum, fragile). The palindrome exists because
past and future have equal weight at each site.

Under depolarizing noise, this balance is destroyed:

- **Past** (immune sector): only {I} per site, giving 1^N = 1 string
  for the entire system. For N=3, one string (III). For N=10, still
  one string.

- **Future** (decaying sector): {X, Y, Z} per site, giving 3^N strings.
  For N=3, twenty-seven strings. For N=10, fifty-nine thousand and forty-nine.

The ratio past/future = (1/3)^N. It is not just unequal. It is
exponentially unequal. The future is exponentially larger than the past.

No bijective mirror can exist between a set of size 1 and a set of size 3^N.
The mirror Π requires exactly as much past as future at each site.
Depolarizing noise breaks this by making almost everything quantum and
almost nothing classical.

In the time-reversal language: depolarizing noise is a universe where
there are exponentially more ways to be undecided than decided. The
arrow of time is not merely present. It is overwhelming.

---

## 8. The Theorem

Everything above condenses into one precise statement. The palindrome
lives or dies based on a single condition about the noise:

**Theorem (Palindromic noise condition).** For N qubits with Heisenberg
coupling on any graph, the Liouvillian palindrome holds under dephasing
noise if and only if the noise has at most two Pauli axes. Equivalently,
at least one of γ_X, γ_Y, γ_Z must be zero.

**Proof.** The palindrome requires a per-site bijection on {I, X, Y, Z}
that pairs dephasing rates symmetrically (Section 3). Such a bijection
exists if and only if the four rates can be partitioned into two pairs
with equal sums. For three-axis noise with rates (γ_X, γ_Y, γ_Z), the
per-site Pauli rates are [0, 2(γ_Y+γ_Z), 2(γ_X+γ_Z), 2(γ_X+γ_Y)].
Exhaustive enumeration shows pairing exists iff at least one rate is
zero. When all three are nonzero, the rate-0 index (I) cannot be paired:
mapping I to any decaying index forces c > 0, requiring the remaining
three indices to have matching partners, but three odd elements cannot
be partitioned into pairs. QED.

**Corollary.** The palindrome error under depolarizing noise (γ_X = γ_Y = γ_Z = γ/3)
is exactly (2/3) N γ, independent of the Hamiltonian.

---

## 9. What Is Proven

### Tier 1 (Theorem):
- Per-site rate table for all noise types (algebraic, from commutation relations)
- Rate-pairing condition: 4 rates must partition into 2 equal-sum pairs
- Depolarizing has 0 valid permutations (exhaustive enumeration)
- Two-axis dephasing has 2 valid permutations each (Z+X, Z+Y, X+Y)
- Error formula: (2/3) Sγ exactly, Hamiltonian-independent
- Interpolation: perfectly linear, no threshold
- Counting argument: weight sectors balance iff per-site split is 2:2

### Tier 1-2 (Standard physics + algebra):
- The palindrome condition reduces to a combinatorial property of the noise
- The obstruction is 1 ≠ 3: no bijection between immune and decaying sets
- Under depol, the future is exponentially larger than the past

---

## References

- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): Π operator and the palindrome theorem
- [Π as Time Reversal](PI_AS_TIME_REVERSAL.md): populations = past, coherences = future
- [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md): the standing wave requires a mirror
- [Non-Heisenberg Palindrome](NON_HEISENBERG_PALINDROME.md): palindrome across all standard models
- Script: [`simulations/depolarizing_analysis.py`](../simulations/depolarizing_analysis.py)
- Results: [`simulations/results/depolarizing_analysis.txt`](../simulations/results/depolarizing_analysis.txt)

---

*The palindrome is a mirror between past and future. Depolarizing noise
breaks the mirror because it makes the future three times larger than
the past at every site. You cannot reflect what does not fit.*
