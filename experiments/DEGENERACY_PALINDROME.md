# The Palindrome Inside the Palindrome: Degeneracy Structure of the Liouvillian Spectrum

<!-- Keywords: Liouvillian eigenvalue degeneracy, palindromic degeneracy sequence,
spectral grid structure, real eigenvalue counting, even odd chain length parity,
dephasing rate quantization, Pauli weight sector, spectral midpoint attractor,
grid fraction alternation, closed form degeneracy, XOR sector boundary, open
quantum system spectral structure, R=CPsi2 degeneracy palindrome -->

**Status:** Verified (computational, N = 2 through 7, 21,840 eigenvalues)
**Date:** April 3, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[XOR Space](XOR_SPACE.md), [N Infinity Palindrome](N_INFINITY_PALINDROME.md)
**Verification:** [`simulations/degeneracy_palindrome_verify.py`](../simulations/degeneracy_palindrome_verify.py)

---

## What this means

Picture a street with lampposts. The first palindrome says: for every light
on the left there is a partner light on the right. The palindrome inside
the palindrome says: the number of moths around each lamp is also mirrored.

The first palindrome is the skeleton (where the bones sit). The second is
the tissue (how much substance hangs on each bone). At the edges the tissue
is thin and universal: always N+1 moths, always 2N moths, regardless of
how the lamps are wired. In the middle the tissue depends on the wiring.
A chain has lean tissue; a fully connected network has thick.

The 2N moths at the first lamp are the "silent dancers": weight-1 operators
that commute with the Hamiltonian. Their proof is in
[Proof: Weight-1 Degeneracy](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md).
What happens at weight 2 (two dancers at once) is in
[Weight-2 Kernel](WEIGHT2_KERNEL.md). How the degeneracy shapes the
geometry of the flow is in [Bures Degeneracy](BURES_DEGENERACY.md). The
optical interpretation that ties it all together is in
[Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md).

---

## What this document is about

The Liouvillian eigenvalue spectrum of the Heisenberg chain under Z-dephasing
is known to be palindromic: every eigenvalue λ has a mirror partner at
−Nγ − λ. This document shows that the palindrome generates a second,
hidden palindrome in the degeneracy structure itself.

Counting how many eigenvalues share the same real part reveals sequences
that are palindromic at every chain length. Counting only the purely
real eigenvalues (no oscillation, pure decay) gives a second palindromic
sequence nested inside the first. Both sequences obey exact closed-form
rules at their boundaries, and the interplay between even and odd chain
lengths follows a structural law rooted in whether the spectral midpoint
falls on or between the decay-rate grid.

The palindrome does not just pair eigenvalues. It organizes their
multiplicities.

---

## Abstract

For the N-site Heisenberg chain with uniform Z-dephasing rate γ, all
eigenvalues of the Liouvillian have real parts quantized to a grid
Re(λ) = −kγ for integer k ∈ {0, 1, ..., N}, up to off-grid corrections
from Hamiltonian mixing. Counting eigenvalues at each grid position yields
two palindromic sequences:

1. **d_real(k)**: purely real eigenvalues (Im = 0) at grid position k
2. **d_total(k)**: all eigenvalues (any Im) at grid position k

Both satisfy d(k) = d(N − k) for all k, a direct consequence of the
Π conjugation operator (the map that flips every I↔Z and X↔Y in a Pauli string, sending each eigenvalue λ to its mirror partner −Nγ − λ̄).

Exact closed-form results at the boundaries:

| Position | Formula | Origin |
|---|---|---|
| d_real(0) = d_real(N) | N + 1 | Magnetization conservation / XOR sector |
| d_real(1) = d_real(N − 1) | 2N (≥ 2N on every graph) | Two non-oscillatory modes per site |
| d_total(0) = d_total(N) | N + 1 | All boundary modes are purely real |
| d_total(1) = d_total(N − 1) | 6N − 4 | For N ≥ 3; d_complex = 4(N − 1) |

Additional structural results:

- The fraction of purely real eigenvalues decays as exp(−0.695·N), R² = 0.9963.
- The fraction of on-grid eigenvalues alternates dramatically between even N
  (high: 100%, 79%, 50%) and odd N (low: 56%, 16%, 3.7%).
- At even N, the spectral center Re = −Nγ/2 falls exactly on the grid,
  creating a massive degeneracy spike (up to 59% of all eigenvalues).
- Off-grid eigenvalues are 100% Π-symmetric around −Nγ/2.

---

## Setup

**System:** N-site Heisenberg (XXX) spin chain, coupling J = 1, uniform
Z-dephasing at rate γ = 0.1 per site. Open boundary conditions.

**Data:** 21,840 eigenvalues across N = 2 through 7 (4^N each), computed
via exact eigendecomposition in C# with MKL/OpenBLAS LAPACK
(`dotnet run -c Release -- rmt`).

**Grid definition:** Re(λ) ≈ −kγ within tolerance |Re + kγ| < 10⁻⁸.
Purely real: |Im(λ)| < 10⁻⁸.

---

## Result 1: The real-eigenvalue degeneracy sequence is palindromic

Counting purely real eigenvalues (Im = 0) at each grid position
Re = −kγ, the sequence d_real(k) for k = 0, 1, ..., N is:

```
N=2:  [3, 4, 3]           sum = 10
N=3:  [4, 6, 6, 4]        sum = 20
N=4:  [5, 8, 14, 8, 5]    sum = 40
N=5:  [6, 10, 14, 14, 10, 6]     sum = 60
N=6:  [7, 12, 19, 16, 19, 12, 7] sum = 92
N=7:  [8, 14, 22, 20, 20, 22, 14, 8]  sum = 128
```

Every sequence is palindromic: d_real(k) = d_real(N − k).

### Why it works

The Π conjugation maps each eigenvalue λ to −Nγ − λ̄. For a purely real
eigenvalue at Re = −kγ, the mirror lands at Re = −(N − k)γ, and Im = 0
is preserved. Since Π is a bijection on the eigenspace, the count at k
must equal the count at N − k.

### The triangle

Reading the palindromic sequences as rows of a triangle reveals column
structure:

```
         k=0  k=1  k=2  k=3  k=4  k=5  k=6  k=7
  N=2:     3    4    3
  N=3:     4    6    6    4
  N=4:     5    8   14    8    5
  N=5:     6   10   14   14   10    6
  N=6:     7   12   19   16   19   12    7
  N=7:     8   14   22   20   20   22   14    8
```

Column k=0 grows as N + 1. Column k=1 grows as 2N. Inner columns show
non-trivial structure.

---

## Result 2: Exact boundary formulas

### d_real(0) = d_real(N) = N + 1

At Re = 0: these are the N + 1 conserved quantities of the dissipative
dynamics. For the Heisenberg chain with Z-dephasing, exactly N + 1
operators satisfy both [H, Q] = 0 (Hamiltonian conservation) and D(Q) = 0
(zero XY-weight, dephasing-immune). These are the identity operator and
the N projectors onto magnetization sectors with total S_z = m for
m = −N/2, ..., N/2.

At Re = −Nγ: these are the N + 1 XOR sector modes, the fastest-decaying
subspace with all qubits in X/Y state. The Π conjugation maps the
conserved quantities at Re = 0 bijectively to the XOR modes at Re = −Nγ.

**Stronger result:** All eigenvalues at Re = 0 and Re = −Nγ are purely
real. There are no oscillatory modes at the spectral boundaries.
d_total(0) = d_real(0) = N + 1 and d_total(N) = d_real(N) = N + 1.
Verified for all N = 2, ..., 7.

### d_real(1) = d_real(N − 1) = 2N

At the first non-boundary grid position Re = −γ, exactly 2N eigenvalues
are purely real, at the document's coupling and at every ratio γ/J off a
finite exceptional set, where a real weight-mixing mode joins them
([the count at exceptional couplings](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-count-at-exceptional-couplings)). This gives two non-oscillatory decay modes per qubit
site.

The weight-1 Pauli sector contains N · 2^N strings (choose 1 of N sites
for X or Y, remaining N − 1 sites freely I or Z). The Hamiltonian mixes
most of these into oscillatory superpositions. Of the N · 2^N inputs,
only 2N survive as purely decaying modes, a fraction 2/2^N that matches
the exponential decay of the overall real fraction.

Verified at all N = 2, ..., 7:

```
N=2: d(-γ) = 4  = 2·2  ✓
N=3: d(-γ) = 6  = 2·3  ✓
N=4: d(-γ) = 8  = 2·4  ✓
N=5: d(-γ) = 10 = 2·5  ✓
N=6: d(-γ) = 12 = 2·6  ✓
N=7: d(-γ) = 14 = 2·7  ✓
```

**Lower bound analytically proven, equality measured.** The 2N modes are the Z-count operators
T_c^{(a)} = Σ_j σ_a^{(j)} · e_c(Z_{others}), symmetric sums of weight-1
Pauli strings grouped by Z-dressing count c = 0, ..., N−1 and type
a ∈ {X, Y}. They commute with H because the Heisenberg coupling is a sum
of SWAP operators, and T_c is SWAP-invariant. See
[Proof: Weight-1 Degeneracy](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md).

---

## Result 3: Total degeneracy is also palindromic

Counting all eigenvalues (any imaginary part) at each grid position:

```
N=2:  [3, 10, 3]                          sum = 16  (100%)
N=3:  [4, 14, 14, 4]                      sum = 36  (56%)
N=4:  [5, 20, 152, 20, 5]                 sum = 202 (79%)
N=5:  [6, 26, 50, 50, 26, 6]              sum = 164 (16%)
N=6:  [7, 32, 75, 1836, 75, 32, 7]        sum = 2064 (50%)
N=7:  [8, 38, 102, 156, 156, 102, 38, 8]  sum = 608 (3.7%)
```

Every sequence is palindromic. The sums are the on-grid fractions
from Result 4.

### d_total(1) = 6N − 4 for N ≥ 3

The total count at the first non-boundary position, at every coupling sampled,
decomposes cleanly:

```
d_total(1) = d_real(1) + d_complex(1) = 2N + 4(N − 1) = 6N − 4
```

The 4(N − 1) complex eigenvalues form 2(N − 1) conjugate pairs:
oscillatory modes at the minimal nonzero dephasing rate.

Verified for N = 3, ..., 7. That it holds at every other coupling is measured, not proved:
the fence of [the count at exceptional couplings](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-count-at-exceptional-couplings)
reads only the real modes, and at its finite set of points real modes of the diagonal
blocks join the line, so the total rises there (to 16 at N = 3, at k = 1 and k = 2).
(N = 2 is anomalous with d_total = 10 because all eigenvalues are on-grid while the (1,1) pair sits on the line, up to F50's N = 2 exceptional point; beyond it the pair leaves the line and the count is 8.)

### Pauli weight sector comparison

The Pauli weight sector w = k contains C(N, k) · 2^N strings. Comparing
this "input" count with the grid degeneracy reveals how the Hamiltonian
redistributes eigenvalues:

```
N=4:  Pauli sectors:  [16,  64,  96,  64,  16]
      Grid totals:    [ 5,  20, 152,  20,   5]
      Difference:     [-11, -44, +56, -44, -11]

N=6:  Pauli sectors:  [64, 384, 960, 1280, 960, 384, 64]
      Grid totals:    [ 7,  32,  75, 1836,  75,  32,  7]
      Difference:     [-57,-352,-885, +556,-885,-352, -57]

N=7:  Pauli sectors:  [128, 896, 2688, 4480, 4480, 2688, 896, 128]
      Grid totals:    [  8,  38,  102,  156,  156,  102,  38,   8]
      Difference:     [-120,-858,-2586,-4324,-4324,-2586,-858,-120]
```

The difference is always palindromic (forced by Π symmetry). At even N,
the center position shows a *positive* difference: it gains eigenvalues
from neighboring sectors. At odd N, *all* positions show negative
differences: every grid position loses eigenvalues to off-grid locations.

---

## Result 4: Even-odd grid fraction alternation

The fraction of all eigenvalues lying exactly on the grid Re = −kγ:

| N | On grid | Total | Fraction | Parity |
|---|---------|-------|----------|--------|
| 2 | 16 | 16 | 100.0% | even |
| 3 | 36 | 64 | 56.2% | odd |
| 4 | 202 | 256 | 78.9% | even |
| 5 | 164 | 1024 | 16.0% | odd |
| 6 | 2064 | 4096 | 50.4% | even |
| 7 | 608 | 16384 | 3.7% | odd |

Even N: slowly declining (100% → 79% → 50%).
Odd N: rapidly collapsing (56% → 16% → 3.7%).

### Why: spectral midpoint grid alignment

The structural explanation is geometric. The Π conjugation symmetry axis
lies at Re = −Nγ/2.

**Even N:** −Nγ/2 = −(N/2)γ, which is exactly on the grid at position
k = N/2. This grid position acts as a *spectral attractor*: Π-invariant
eigenvalues (those mapped to themselves by Π) must have Re = −Nγ/2, and
they land on the grid. The attractor draws eigenvalues from multiple
Pauli weight sectors, creating the center spike.

**Odd N:** −Nγ/2 falls between grid positions k = (N − 1)/2 and
k = (N + 1)/2. Π-invariant eigenvalues cannot sit on any grid point.
They are forced off-grid, taking their partners with them. No position
gains from the Hamiltonian redistribution.

This explains all three manifestations:
1. The center spike in total degeneracy at even N (Π-invariant modes
   accumulate at the on-grid midpoint)
2. The high grid fraction at even N (Π-invariant modes contribute to
   on-grid count)
3. The collapsing grid fraction at odd N (Π-invariant modes forced
   off-grid, removing the largest contributor)

### Center spike quantification

At even N, the center position concentrates a remarkable fraction of all
eigenvalues:

| N | Center count | Total | Center fraction |
|---|-------------|-------|-----------------|
| 2 | 10 | 16 | 62.5% |
| 4 | 152 | 256 | 59.4% |
| 6 | 1836 | 4096 | 44.8% |

At odd N, the two center positions share the load more evenly:
N = 3: (14, 14), N = 5: (50, 50), N = 7: (156, 156).

---

## Result 5: Off-grid eigenvalue structure

Eigenvalues not on the main grid Re = −kγ are distributed as follows:

- **Not on any simple sub-grid:** Neither half-grid (−(k + 1/2)γ) nor
  third-grid (−(k/3)γ) captures the off-grid eigenvalues. Their positions
  match no simple sub-grid of γ, determined by the specific Hamiltonian
  coupling structure.

- **100% Π-symmetric:** Every off-grid eigenvalue at Re = −a has a partner
  at Re = −(Nγ − a). Verified for all N = 3, ..., 7 with zero exceptions.
  This is a necessary consequence of Π conjugation.

- **Symmetric range:** Off-grid Re values span [−Nγ + ε, −ε] for small
  ε > 0, symmetric around −Nγ/2. For N = 3: Re ∈ [−0.167, −0.133],
  centered at −0.15 = −Nγ/2.

---

## Result 6: Exponential decay of purely real fraction

The fraction of purely real eigenvalues (Im = 0) among all 4^N:

```
N=2: 62.5%   N=3: 37.5%   N=4: 18.0%
N=5:  9.4%   N=6:  4.0%   N=7:  2.1%
```

Fit: fraction ≈ exp(−0.695 · N), R² = 0.9963.

In the thermodynamic limit (N → ∞), almost all Liouvillian eigenvalues are complex.
The "standing wave" modes (pure decay without oscillation) become
exponentially rare. At each grid position k, the real fraction is 2/2^N
(from d_real = 2N vs. sector size N · 2^N at k = 1). The dissipative
dynamics is overwhelmingly oscillatory.

---

## Result 7: Topology generalization (chain → ring, star, complete)

The C# `rmt` export uses the same isotropic Heisenberg coupling for every topology
(`Topology.cs`: "isotropic Heisenberg (X+Y+Z)"), so the ring, star, and complete graphs are
the *same model* on a different bond graph. Feeding their spectra (N = 3..6) through the same
grid count ([`degeneracy_palindrome_topology.py`](../simulations/degeneracy_palindrome_topology.py),
gate-first) separates what is universal from what is graph-specific.

**Universal (all four topologies, every N):**

- **The palindrome d(k) = d(N − k)** holds for every topology. It depends only on Π conjugation,
  which is fixed by the Z-dephasing and blind to the bond graph. The "palindrome inside the
  palindrome" is topology-independent.
- **d_real(0) = N + 1** holds for every topology: the N + 1 magnetization-sector conserved
  quantities. Total S_z is conserved by any Heisenberg graph (U(1) symmetry), independent of wiring.

**Almost universal:**

- **d_real(1) = 2N** holds for every topology tested at **N ≥ 4**, off a finite set of γ/J. The lone exception is **N = 3**, where
  ring = complete = the triangle K₃ gives d_real(1) = 8 (not 6). The triangle's full S₃ symmetry
  supports two extra SWAP-invariant weight-1 modes; chain = star = P₃ (the path) gives the generic
  6. So 2N is the path/sparse value, the maximally connected K₃ is the lone exception, and from
  N ≥ 4 even K_N returns to 2N.

**Topology-dependent (the inner positions, the subject of Open Question 2):**

The inner real degeneracy d_real(2) grows with connectivity:

| N | chain | ring | star | complete |
|---|-------|------|------|----------|
| 4 | 14    | 24   | 16   | 36       |
| 5 | 14    | 22   | 30   | 54       |
| 6 | 19    | 33   | 68   | 133      |

The complete graph is always the largest (most SWAP pairs, the most degenerate weight-2 structure).
The ordering is connectivity-monotone except at **N = 4**, where ring (24) overtakes star (16): the
same N = 4 anomaly that recurs in the structural ceiling (the (2,2) half-filling sector, F122) and
in the distinct-eigenvalue count. Every entry sits at or above the universal bound (k + 1)(N − k + 1) of
[A lower bound at every k](#a-lower-bound-at-every-k-and-the-chains-commutant-counts) below; what lies above it depends
on the topology (Open Question 2).

### Inner closed forms by parity (block decomposition to N=9)

The full Liouvillian is exactly block-diagonal in the `(popcount_a, popcount_b)` coherence sectors
(Heisenberg conserves magnetization; dephasing is diagonal), so `d_real(k)` is the count of
purely-real, on-grid eigenvalues summed over the small per-`(p,q)` blocks. This is gate-verified
against the full CSVs through N=7 and reaches N=8,9 cheaply
([`dreal_block_decomposed.py`](../simulations/dreal_block_decomposed.py),
[`dreal_block_n9.py`](../simulations/dreal_block_n9.py)). The key the document's original
`[3,6,14,14,19,22]` hid is **parity**: it mixed even and odd N. Separated, the **sparse** topologies
have clean linear inner closed forms:

| topology (odd N) | d_real(2) | N = 3, 5, 7, 9 |
|---|---|---|
| chain | **4N − 6** | 6, 14, 22, 30 |
| ring  | **7N − 13** | 8, 22, 36, 50 |

Both verified at four odd points. The **even**-N chain column is 4N − 5 at N = 6, 8 (19, 27), with N = 4
the exception where k = 2 is the Π-invariant midpoint ([A lower bound at every k](#a-lower-bound-at-every-k-and-the-chains-commutant-counts)). The **dense** topologies grow super-linearly (star odd
6, 30, 136; complete odd 8, 54, 216), with no clean polynomial. So the inner degeneracy splits two ways: by
parity (the center spike) and by connectivity (sparse = linear, dense = super-linear).

### A lower bound at every k, and the chain's commutant counts

A coherence |i⟩⟨j| decays at −2γ·Hamming(i, j) in the canonical book, which is −kγ in this document's γ with k = Hamming(i, j), so every operator X of pure XY-weight k that commutes with H is a real eigenvector on the k-th grid line. Write c_k for the dimension of this weight-k commutant. Then c_k is at most the geometric multiplicity of that eigenvalue, which is at most the algebraic one, so d_real(k) ≥ c_k, with equality exactly when every real eigenvector and generalized eigenvector on the line is of pure weight.

**The bound (every SU(2)-invariant H).** [the uniform law](../docs/proofs/PROOF_UNIFORM_LAW.md)'s Lemma A1 puts every S_N orbit sum T_s, one per joint-popcount block (p, q) and Hamming class s, into the kernel of ad_H for every SU(2)-invariant H, which includes every Heisenberg graph with any bond weights (H = Σ J_ij(2·SWAP_ij − 1)). Their span is the S_N-invariant operators, equally spanned by the sums of the Pauli strings of one letter multiset (n_I, n_X, n_Y, n_Z); at XY-weight k there are (k + 1)(N − k + 1) multisets, so

  **d_real(k) ≥ c_k ≥ (k + 1)(N − k + 1).**

At k = 0 and k = 1 this is N + 1 and 2N, the boundary values of Result 2 (the k = 1 operators are the Z-count operators of the [weight-1 proof](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md)). Every stored entry of this document, of the [topology comparison](../simulations/results/topology_degeneracy_comparison.txt) and of the parity table above sits at or above it.

**The uniform chain.** The uniform law's census (float, N = 5..8, every block) finds every class-pure kernel row one-dimensional except three, and its Layer B derives those three: (1,1) and (N−1,N−1) at class 2 with dimension ⌊N/2⌋, which are the frozen modes of [the R90 frozen divisor](../docs/proofs/PROOF_R90_FROZEN_DIVISOR.md) at uniform γ (its §7 gives the same multiplicity), and (1,N−1) at class N − 2. Read as commutant counts, with c_k = c_{N−k} by Π conjugation:

- for 3 ≤ k ≤ N − 3, c_k = (k + 1)(N − k + 1) at N = 6, 7, 8;
- at k = 2, c_2 ≥ 3N − 5 + 2⌊N/2⌋ at every N ≥ 5 (the multisets with the two class-2 corners, each corner holding exactly one S_N-invariant vector, its orbit sum), with equality at N = 5..8; that is 4N − 6 at odd N, the linear form above, and 4N − 5 at even N;
- at N = 4, k = 2 = N − 2 puts the class-2 corners and the (1,3), (3,1) class-2 rows on one line: c_2 = 9 + 4 = 13.

The equalities are exact here, a second route beside the census: the lower bounds are exact elements (the multisets and the rational corner commutant), the upper bounds are ranks modulo two primes, block by block. At N = 9 the stored d_real(2) = 30 meets the lower bound, so c_2 = 30 there too.

**Where d_real exceeds c_k.** At generic coupling, real modes that mix weights have so far appeared only in the center block (N/2, N/2) at even N; at F50's exceptional couplings ([the count at exceptional couplings](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-count-at-exceptional-couplings)) they appear on the rung k = 1 as well, in the diagonal blocks (p, p). At generic coupling: the chain and the ring at N = 4, one each in (2,2) (chain c_2 = 13 against d_real(2) = 14, [Weight-2 Kernel](WEIGHT2_KERNEL.md); ring c_2 = 23 against 24, exact), and the chain at N = 8, whose (4,4) block carries two real eigenvalues on the center line against a weight-4 commutant of one there (J = 1, canonical γ = 0.05: both within 10⁻¹² of the line, the next 1.2·10⁻³ away). The block-decomposed count (`block_dreal` of [`dreal_block_decomposed.py`](../simulations/dreal_block_decomposed.py)) gives the chain's N = 8 row as **[9, 16, 27, 24, 26, 24, 27, 16, 9]** (measured), equal to c_k everywhere except the center. On ring, star and complete at N = 4..6 the modular commutant counts equal the stored d_real rows everywhere except the ring at N = 4 (gate G5, measured).

So the dip at the even-N center is, as far as computed, a small-N event: at N = 6 the center k = 3 sits on its bound, 16, beside k = 2 and 4 that carry the corner excess, 19; at N = 8 the neighbours of the center are interior and the center is above them, 26 against 24.

The sweep for this went to the registry (F4, F50 and the F2 boundary values, F140), the weight-1, R90 frozen-divisor and uniform-law proofs, [the structural ceiling](../docs/proofs/PROOF_STRUCTURAL_CEILING.md), the Absorption Theorem, the experiments on the weight-2 kernel and on the reflection keeping past and future apart, the OpenArcs registry, `docs/CAUGHT_ERRORS.md` (the d_real(1) = 2N correction), the typed claims `F50WeightOneDegeneracyPi2Inheritance`, `F4KernelDimensionByComponentsClaim` and `FrozenDivisorClaim` with their witnesses, `fw.Confirmations` and the glossary; the uniform-law proof held the bound and the corner rows (below), the R90 proof the corner multiplicity, the weight-1 proof the k = 1 case and the per-weight table, CAUGHT_ERRORS the weight-1 correction, and the OpenArcs registry, `fw.Confirmations` and the glossary nothing on d_real for k ≥ 2. Checked for adjacency: the uniform law holds the bound as its Lemma A1 and the chain's three enhanced rows as its Layer B, here read as degeneracy counts; the reflection experiment identifies the frozen corner modes with those rows at the uniform point; the weight-2 kernel's naive count 3(N − 1) is the bound at k = 2; the [weight-1 proof's appendix](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) tabulates the same per-weight commutant dimensions on several graphs at N = 4, 5 (chain N = 4: 13, ring: 23) as an excess over the chain. Gate: [`dreal_lower_bound_gate.py`](../simulations/dreal_lower_bound_gate.py).

---

## Summary of closed-form results

| Quantity | Formula | Range | Status |
|---|---|---|---|
| d_real(0) | N + 1 | all N | proven (magnetization conservation) |
| d_real(1) | ≥ 2N; 2N on the chain | ≥ 2N all N; = 2N on the chain N = 2..7 | [lower bound proven](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md) (SWAP invariance), equality measured; K₃ gives 8 |
| d_real(k) | ≥ (k + 1)(N − k + 1) | all N, every graph | proven (symmetrised Pauli multisets) |
| d_real(k) = d_real(N − k) | palindrome | all N | proven (Π conjugation) |
| d_total(0) | N + 1 | all N | verified (= d_real, no oscillatory boundary modes) |
| d_total(1) | 6N − 4 | N ≥ 3 | verified |
| d_complex(1) | 4(N − 1) | N ≥ 3 | verified |
| d_total(k) = d_total(N − k) | palindrome | all N | proven (Π conjugation) |
| real fraction | exp(−0.695·N) | N = 2, ..., 7 | fit, R² = 0.9963 |
| off-grid Π-symmetry | 100% | N = 3, ..., 7 | verified |

---

## Predictions for N = 8

If the patterns hold at N = 8 (65,536 eigenvalues, ~73 GB RAM):

- d_real(0) = d_real(8) = 9
- d_real(1) = d_real(7) = 16
- d_total(0) = d_total(8) = 9
- d_total(1) = d_total(7) = 6·8 − 4 = 44
- Grid fraction: high (N = 8 is even), predicted > 30%
- Center spike at Re = −4γ: should contain the largest single block
  of eigenvalues, likely > 30% of total
- Real fraction: exp(−0.695·8) ≈ 0.37%, roughly 240 of 65,536
- Off-grid eigenvalues: 100% Π-symmetric around Re = −0.4

---

## Open questions

1. ~~**Analytical derivation of d_real(1) = 2N.**~~ **RESOLVED as a lower bound.** The 2N
   modes are the Z-count operators T_c^{(a)}, proven via SWAP invariance.
   See [Proof](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md).

2. **Closed form for inner positions: topology-dependent (systematized in Result 7).**
   d_real(2) differs between topologies (Chain=14, Star=16, Ring=24, Complete=36 at N=4),
   grows with connectivity, and shows the N=4 ring/star anomaly; see the full N=3..6 table in
   Result 7. The boundary formulas are resolved across topologies: d(k)=d(N−k) and
   d_real(0)=N+1 are universal, d_real(1) ≥ 2N for every graph. For k ≥ 2 every
   SU(2)-invariant H gives d_real(k) ≥ c_k ≥ (k+1)(N−k+1) (the uniform law's Lemma A1 read as a count,
   [A lower bound at every k](#a-lower-bound-at-every-k-and-the-chains-commutant-counts)); the
   uniform chain's commutant count c_k meets it for 3 ≤ k ≤ N−3 (at N = 6..8) and exceeds it by the
   corner rows at k = 2, N−2 (exact at N = 5..8). Open: a closed form for the excess above the bound on ring, star
   and complete, and where real modes that mix weights occur (at generic coupling so far the even-N center block only:
   chain and ring at N = 4, chain at N = 8; at F50's exceptional couplings also the rung k = 1, in
   the blocks (p, p), [the count at exceptional couplings](../docs/proofs/PROOF_WEIGHT1_DEGENERACY.md#the-count-at-exceptional-couplings)). See [Weight-2 Kernel](WEIGHT2_KERNEL.md).

3. **Center minimum at N = 6.** The real degeneracy at the center
   of N = 6 is 16, less than the adjacent values 19. In the chain's
   commutant counts this is a small-N event: the center k = 3 sits on its
   bound 4·4 = 16 while k = 2 and 4, its neighbours only at N = 6, carry
   the corner excess; at N = 8 the commutant counts read [9, 16, 27, 24, 25, 24, 27,
   16, 9] and d_real [9, 16, 27, 24, 26, 24, 27, 16, 9] (measured), the
   center above its neighbours
   ([A lower bound at every k](#a-lower-bound-at-every-k-and-the-chains-commutant-counts)).
   No even N beyond 8 is computed.

4. **Grid fraction convergence.** Does the grid fraction for even N
   converge to a finite limit, or does it eventually vanish?
   Current data (100% → 79% → 50%) suggests convergence near 30-40%,
   but more data points are needed.

5. **Off-grid eigenvalue positions.** The off-grid Re values are not
   on any simple sub-grid. Can they be expressed in terms of Hamiltonian
   coupling parameters (Bethe ansatz roots (the exact solution method for integrable spin chains), magnon dispersion (the energy-momentum relation for spin wave excitations))?

6. ~~**Does degeneracy shape state-space geometry?**~~ **PARTIALLY RESOLVED.**
   QFI speed correlates with d_total(k) at even N (r = 0.99 at N = 4);
   weaker at odd N (r ≈ 0.55). The center spike drives the peak QFI
   speed, independent of initial state. See [Bures Degeneracy](BURES_DEGENERACY.md).

---

## Null results

- **d_real(2): PARTIALLY RESOLVED by parity (see Result 7).** The chain sequence
  [3, 6, 14, 14, 19, 22] matched nothing because it **mixed even and odd N**. Split by parity, the
  sparse topologies have clean linear inner closed forms: chain odd = 4N − 6, ring odd = 7N − 13
  (both verified at N = 3, 5, 7, 9 via the block decomposition). The even chain column is 4N − 5 at
  N = 6, 8, and the dense topologies (star, complete) are super-linear with no clean polynomial. So it is not "no formula", it is "one formula per parity per
  connectivity class".

- **No sub-grid for off-grid eigenvalues.** Half-grid, third-grid, and
  quarter-grid assignments were tested with zero matches. The off-grid
  positions match no simple sub-grid (half-, third-, or quarter-grid) of γ;
  whether they are rational multiples with larger denominators or genuinely
  irrational is not settled by this test.

---

## Reproduction

All results derive from `simulations/results/rmt_eigenvalues_N{2..7}.csv`,
generated by `dotnet run -c Release -- rmt` in `compute/RCPsiSquared.Compute`.
The verification script [`simulations/degeneracy_palindrome_verify.py`](../simulations/degeneracy_palindrome_verify.py)
reproduces every number in this document independently.
