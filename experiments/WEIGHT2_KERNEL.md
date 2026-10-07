# The Weight-2 Commutator Kernel: Where the Simple Story Ends

<!-- Keywords: weight-2 commutator kernel, Liouvillian degeneracy topology dependence,
SWAP representation theory, Heisenberg chain higher weight, d_real(2) formula,
multi-weight eigenvalue mixing, R=CPsi2 weight-2 kernel -->

**Status:** Numerically characterized (N = 3 through 6, three topologies); at N = 5..8 the chain's weight-2 count is the multiset bound plus the corner commutant, exactly ([DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md#a-lower-bound-at-every-k-and-the-chains-commutant-counts))
**Date:** April 3, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Proof: Weight-1 Degeneracy](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md),
[Degeneracy Palindrome](DEGENERACY_PALINDROME.md)
**Verification:** [`simulations/weight2_kern_analysis.py`](../simulations/weight2_kern_analysis.py)

---

## What this means

With one dancer, everything is simple: swap two positions, the dancer
moves, the sum does not change. With two dancers at once, the symmetry
breaks. When you swap the two dancers themselves, their types can mix
(XY becomes YX). Different networks mix differently. A chain has few
connections and few silent two-dancer modes. A complete network has many
connections and many more.

The universal simplicity ends at one dancer. Beyond that, the specific
wiring of the network shapes how many silent modes exist. The world
becomes complicated, and complex.

How this degeneracy structure shapes the geometry of the quantum state
flow is in [Bures Degeneracy](BURES_DEGENERACY.md). The optical
interpretation is in [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md).

---

## What this document is about

At weight 1, the Liouvillian degeneracy has a clean formula: d_real(1) ≥ 2N
on every graph (proven, by SWAP-invariant Z-count operators), with equality
on the chain (measured) and the K_3 exception at N = 3. This document asks: does the same structure extend
to weight 2?

The answer is no. The weight-2 commutator kernel is more complex: its dimension
depends on the topology, and beyond its 3(N − 1) SWAP-invariant vectors (the
multiset bound) it holds vectors without definite SWAP eigenvalues, so the
simple group-theoretic argument from weight 1 covers only part of it. This is not a failure of the method;
it is a structural discovery about where the symmetry protection ends.

---

## Results

### Weight-2 kernel dimension: topology-dependent

| N | Chain | Star | Complete | w=2 sector dim |
|---|-------|------|----------|----------------|
| 3 | 6 | n/a | **8 (= K_3, +2)** | 24 |
| 4 | 13* | 16 | 36 | 96 |
| 5 | 14 | 30 | 54 | 320 |
| 6 | 19 | n/a | n/a | 960 |

*Chain N=4: d_real(2) = 14, but ker([H,·]|_{w=2}) = 13. One eigenvalue
at Re = −0.2 comes from multi-weight mixing (not a pure w=2 mode).

**2026-05-17 added: K_3 (= ring = complete at N=3) gives ker = 8 at weight-2,
not 6 like chain. This is the palindromic partner (w → N − w = 3 − 2 = 1)
of the K_3 N=3 weight-1 anomaly recently identified as the 2-dim standard
irrep of S_3 = Aut(K_3). The F1 palindrome forces the +2 excess at w=1 and
the +2 excess at w=2 to pair up. See [`PROOF_WEIGHT1_DEGENERACY § Appendix
(2026-05-17)`](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) for the
matrix-commutator framework and per-weight centralizer breakdown that
exposes the pairing.**

All other (non-K_3-N=3) values match d_real(k=2), except the chain at N = 4, where one real mode on that line mixes weights (13 against 14; the ring at N = 4 has one too, DEGENERACY_PALINDROME).

### SWAP eigenvalue structure: the key difference from weight 1

**Weight 1:** On the chain every kernel vector satisfies SWAP_{i,i+1}(v) = v
for every bond (measured; the weight-1 proof's construction gives the 2N SWAP-invariant ones,
its triangle-inequality step, which would exclude others, has a gap, and K_3 at N = 3 has two more that are not).

**Weight 2 (N ≥ 4):** In the numerical kernel basis the solver returns, no basis vector has a
definite SWAP eigenvalue: SWAP_{i,i+1}(v) is neither +v nor −v. The kernel itself does contain
SWAP-invariant vectors, 3(N − 1) of them at every N (the symmetrised Pauli multisets of
[DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md#a-lower-bound-at-every-k-and-the-chains-commutant-counts)),
so the table below reads that basis, not the kernel's representation content. For the chain the kernel
is not S_N-stable (swapping sites 0 and 1 carries kernel vectors out of it), so no S_N representation is read
on it; on graphs whose automorphism group is S_N the irreps apply. The table's column "Mixed" counts basis vectors that
are not SWAP eigenvectors (the symmetric group of all N! permutations; "trivial" means every permutation acts as +1, "mixed" means permutations act as matrices that are neither all +1 nor all ±1).

| N | Trivial (+1) | Alternating (−1) | Mixed | Total |
|---|-------------|-----------------|-------|-------|
| 3 | 6 | 0 | 0 | 6 |
| 4 | 0 | 0 | 13 | 13 |
| 5 | 0 | 0 | 14 | 14 |
| 6 | 0 | 0 | 19 | 19 |

At N = 3 the whole kernel is SWAP-invariant (6 = 3(N − 1)), matching the weight-1 pattern. At N = 4..8
the kernel exceeds its SWAP-invariant part by the corner commutant (see Excess below).

### Type class decomposition: no antisymmetric modes

The weight-2 Pauli strings have four type classes based on the active
pair: XX, XY, YX, YY. The SWAP between the two active positions maps
XY ↔ YX, creating symmetric (XY + YX) and antisymmetric (XY − YX)
combinations.

**Hypothesis (from task):** The kernel contains both symmetric and
antisymmetric type-class vectors.

**Result:** All kernel vectors have purely symmetric type pairing
(c_{XY} = c_{YX} for every Z-dressing). Zero antisymmetric pairs
observed. The XY − YX sector contributes nothing to the kernel.

### Why the weight-1 proof breaks at weight 2

The weight-1 proof uses three steps:
1. [H, v] = 0 → Σ SWAP_i(v) = (N−1)·v
2. Triangle inequality → each SWAP_i(v) = v individually
3. S_N-invariance → unique invariant per orbit → dim = 2N

At weight 2, step 2 fails for part of the kernel. The 3(N − 1) multiset vectors satisfy
SWAP(v) = v for each bond; the corner commutant vectors beyond one orbit sum per corner do not, and for them
[H, v] = 0 is achieved through cancellation: different SWAPs pull v in different directions, but
the sum cancels to zero. This is possible because the weight-2 sector
has richer representation structure.

The triangle inequality still applies, but its conclusion is different.
At weight 1, ||SWAP(v)|| = ||v|| and there are N−1 terms, giving
||(N−1)v|| = (N−1)||v|| which forces parallelism. At weight 2, the
equation Σ SWAP_i(v) = (N−1)·v still holds, but the individual
SWAP_i(v) vectors are not all parallel to v; they lie in a subspace
that allows cancellation of the non-parallel components.

### Excess over naive count

The naive count from symmetric types (XX, YY, XY+YX) with Z-count
orbits gives 3 × (N − 1) kernel vectors (3 types × (N−1) Z-counts
for N−2 passive positions):

| N | Kernel dim | Naive 3(N−1) | Excess |
|---|-----------|-------------|--------|
| 3 | 6 | 6 | 0 |
| 4 | 13 | 9 | 4 |
| 5 | 14 | 12 | 2 |
| 6 | 19 | 15 | 4 |

The naive count 3(N − 1) is the lower bound (k + 1)(N − k + 1) at k = 2, which holds for every
graph. The excess is the corner commutant: the (1,1) and (N−1,N−1) blocks carry ⌊N/2⌋ weight-2
commutant vectors each, the frozen modes of the R90 frozen divisor at uniform γ, one of each
SWAP-invariant, so the excess is at least 2(⌊N/2⌋ − 1) for N ≥ 5 and equal to it at N = 5..8 (2 at N = 5, 4 at N = 6); at N = 4,
k = 2 = N − 2 puts both corner pairs on one line and the excess is 4
([DEGENERACY_PALINDROME](DEGENERACY_PALINDROME.md#a-lower-bound-at-every-k-and-the-chains-commutant-counts)).

---

## The multi-weight mixing phenomenon

For Chain N = 4, d_real(2) = 14 but dim(ker([H,·]|_{w=2})) = 13.
One purely-real eigenvalue at Re = −4γ = −0.2 is not a pure weight-2
mode. It arises from Hamiltonian mixing between weight sectors: the
eigenoperator has components at multiple weights (w=0, w=2, w=4, ...),
and the combination of different dephasing rates averages to exactly
Re = −4γ.

The same `d_real(k) > ker(w=k)` discrepancy was originally suspected
to also explain the **Ring/Complete N = 3 weight-1 anomaly** (d_real(1)
= 8 > 2N = 6). **2026-05-17 update / correction:** native C# verification
of F50 + explicit per-weight Pauli decomposition of the 2 K_3 extras
shows they are **pure weight-1** (`|c_α|² in weight-1 sector = 1.000000`,
all other weight sectors exactly zero). The K_3 N=3 anomaly is NOT
multi-weight mixing; the extras live entirely inside the weight-1
sector. The correct attribution is the **2-dim standard irreducible
representation of S_3 = Aut(K_3)** acting on the weight-1 c=1 Pauli
strings; see [`PROOF_WEIGHT1_DEGENERACY § Appendix (2026-05-17)`](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md)
for the proof's Step-5 gap analysis (matrix-commutator vs conjugation-
action) and the structural identification.

The Chain N=4 weight-2 `+1` excess at Re = −4γ remains an honest
multi-weight-mixing case (verified there explicitly via the eigenoperator's
multi-weight Pauli content). Multi-weight mixing and irrep-induced extras
are two distinct mechanisms for `d_real(k) > ker(w=k)`; the K_3 N=3 case
exemplifies the latter.

---

## Null results

- **No antisymmetric kernel vectors.** The XY − YX type class was
  hypothesized to contribute to the kernel via the alternating
  representation of S_N. This was not observed at any N.

- **No known combinatorial family for the raw sequence.** The sequence [6, 13, 14, 19]
  for Chain N = 3, ..., 6 (and 14 for d_real(2) including the mixed
  mode at N=4) matches none of the families tried; its structure is the
  multiset bound plus the corner commutant, c_2 = 3N − 5 + 2⌊N/2⌋ at N = 5..8
  (a lower bound at every N ≥ 5; DEGENERACY_PALINDROME).

- **Burnside counting insufficient.** The orbit-counting approach (Burnside's lemma: count distinct patterns by averaging fixed points over all group elements)
  (number of S_N-invariant vectors per orbit) accounts for only
  3(N−1) of the kernel vectors. The excess is the corner commutant of
  the (1,1) and (N−1,N−1) blocks (Excess above).

- **Triangle inequality does not force individual SWAP fixation.**
  At weight 1, the triangle inequality argument aims at SWAP(v) = v for
  each bond individually. At weight 2 (N ≥ 4), the vectors of the numerical kernel basis have
  SWAP ratios that are neither +1 nor −1 but continuously varying
  between bonds. The sum cancels to zero through destructive
  interference, not through individual fixation.

---

## What this means for the project

1. **Weight 1 is special, except K_3 N=3.** The clean `d_real(1) = 2N`
   formula with its topology-independent SWAP proof is the exception
   for k ≥ 2 (this document) but is itself violated at one specific
   small-graph case: the K_3 (= triangle, ring at N=3, complete on 3
   vertices) gives `d_real(1) = 8` instead of `2N = 6`. The 2 extras
   are weight-1 operators in the S_3 standard 2-dim irrep, not multi-
   weight mixing; see [`PROOF_WEIGHT1_DEGENERACY § Appendix (2026-05-17)`](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md).
   At F50's exceptional couplings, a finite set of γ/J, real weight-mixing
   modes of the diagonal blocks (p, p) do reach rung 1
   ([the count at exceptional couplings](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-count-at-exceptional-couplings)).
   So even weight 1 has irrep structure beyond the trivial; it just
   needs one specific graph (K_3) to surface. Weight 2 and beyond are
   governed by the corner commutant on the chain and by the automorphism group where it is large.

2. **The topology matters.** A universal formula `d_real(k, N)` for k ≥ 2
   does not exist. Any formula must incorporate the bond structure
   of the graph. At k = 1 the same holds in degenerate form: K_3 N=3
   is the unique counterexample to topology-independence.

3. **The palindrome is deeper than the kernel.** The palindrome
   `d_real(k) = d_real(N − k)` holds for every topology (proven via Π
   conjugation). But the individual values `d_real(k)` depend on the
   topology. The palindrome is a spectral symmetry; the degeneracy
   counts are dynamical.

4. **The Trivial/Alternating/Mixed table reads a basis.** For the chain the count is the multiset bound plus the corner
   commutant at N = 5..8 (DEGENERACY_PALINDROME); the table format stays useful where the automorphism group acts: the K_3 N=3 weight-1 anomaly slots into it
   with one new row and a "Standard (2-dim, S_3)" column. On graphs with a large automorphism group an
   irrep-by-irrep refinement at weight 2 is open.

---

## Reproduction

- Weight-2 kernel analysis: [`simulations/weight2_kern_analysis.py`](../simulations/weight2_kern_analysis.py)
- Output: [`simulations/results/weight2_kern_analysis.txt`](../simulations/results/weight2_kern_analysis.txt)
- Topology eigenvalues: `dotnet run -c Release -- rmt {chain|star|ring|complete}`
- Topology comparison: [`simulations/topology_degeneracy_compare.py`](../simulations/topology_degeneracy_compare.py)
