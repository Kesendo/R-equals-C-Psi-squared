# Why d(−γ) = 2N: The Commutator Kernel in the Weight-1 Sector

<!-- Keywords: commutator kernel, Heisenberg Hamiltonian weight-1 sector, Liouvillian
degeneracy proof, SWAP invariance Pauli strings, Z-count conserved operators,
palindromic degeneracy analytical derivation, d(-gamma) = 2N proof, open quantum
system spectral structure, SU(2) symmetry weight-1, R=CPsi2 commutator kernel -->

**Status:** Tier 1 derived (SWAP-invariance + triangle-inequality proof for all N on any connected graph) + Tier 2 verified (bit-exact numerical match N = 2 through 7)
**Date:** 2026-04-03
**Authors:** Thomas Wicht, Claude (Anthropic)
**Statement:** `d_real(Re = −2γ) = 2N` for the isotropic Heisenberg coupling on any connected graph under uniform Z-dephasing. The 2N kernel operators are the symmetric sums `T_c^{(a)} = Σⱼ σ_a^{(j)} ⊗ Z_S(c) ⊗ I_rest` grouped by active type a ∈ {X, Y} and Z-count c ∈ {0, ..., N−1}.
**Typed claim:** [`F50WeightOneDegeneracyPi2Inheritance.cs`](../../compute/RCPsiSquared.Core/Symmetry/F50WeightOneDegeneracyPi2Inheritance.cs) (Tier 1 derived; both 2s in 2N and 2γ inherit from the Pi2 dyadic ladder's a₀ anchor; F50 entry in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md)).
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Degeneracy Palindrome](../../experiments/DEGENERACY_PALINDROME.md),
[Mirror Symmetry Proof](MIRROR_SYMMETRY_PROOF.md)
**Verification:** [`simulations/kommutator_kern_analysis.py`](../../simulations/kommutator_kern_analysis.py),
[`simulations/verify_triangle_inequality.py`](../../simulations/verify_triangle_inequality.py),
`dotnet run -c Release -- eigvec`

---

## What this means

In a chain of N mirrors, exactly one mirror can "dance" (X or Y) while the
others stand still (I or Z). There are many such configurations. But only
2N of them die silently, without oscillation, without echo. Why? Because the
Heisenberg Hamiltonian is a swap operator. If you swap two mirrors, the
dancer moves, but the uniform sum over all positions does not change. What
cannot change does not oscillate. What does not oscillate dies silently.
Two costumes (X and Y) times N outfits (0 to N−1 silent companions dressed
in Z) equals 2N.

This works for any connected graph, not just a chain. It is pure algebra,
independent of the specific wiring.

What happens when two mirrors dance at once is fundamentally different:
see [Weight-2 Kernel](../../experiments/WEIGHT2_KERNEL.md).

---

## What this document is about

The Liouvillian spectrum of the Heisenberg chain with Z-dephasing has exactly
2N purely-real eigenvalues at the first non-zero grid position Re = −2γ.
This document proves that these 2N modes are weight-1 Pauli operators that
commute with the Hamiltonian, and constructs all 2N of them explicitly via a
SWAP-invariance argument.

The construction works for any isotropic Heisenberg model on any graph, not
just the chain. The key insight: the 2N conserved operators are symmetric sums
of weight-1 Pauli strings grouped by Z-count.

---

## Abstract

For the N-qubit isotropic (XXX) Heisenberg Hamiltonian
H = J Σ_{(i,j)} (X_i X_j + Y_i Y_j + Z_i Z_j) on any graph, define the
Z-count operators:

**T_c^{(a)} = Σ_{k} Σ_{S ⊂ complement(k), |S|=c} σ_a^{(k)} ⊗ Z_S ⊗ I_{rest}**

for a ∈ {X, Y} and c ∈ {0, 1, ..., N−1}, where σ_a is the Pauli-a matrix at
site k, Z_S places Z on every site in S, and I on all remaining sites.

**Theorem.** [H, T_c^{(a)}] = 0 for all c and a.

**Proof.** H = J Σ_{bonds} (2·SWAP_{ij} − I). Each SWAP_{ij} permutes the
sites of a weight-1 Pauli string but preserves both the active type (X or Y)
and the Z-count. Since T_c^{(a)} sums over all weight-1 strings of type a with
Z-count c, and SWAP is a bijection on this set, SWAP leaves T_c^{(a)} invariant.
Hence [SWAP_{ij}, T_c^{(a)}] = 0, and [H, T_c^{(a)}] = 0.

The 2N operators {T_c^{(X)}, T_c^{(Y)} : c = 0, ..., N−1} are linearly
independent and orthogonal in the Hilbert-Schmidt inner product (the operator analog of the dot product: Tr(A†B), zero when operators have no Pauli strings in common), since they
involve disjoint sets of Pauli strings. This gives dim(ker) ≥ 2N.

The reverse inequality dim(ker) ≤ 2N follows from a triangle inequality
argument: any kernel element must be fixed by every individual SWAP (not just
their sum), reducing the invariant space to exactly one vector per transitive
orbit. This establishes:

**dim(ker([H, ·]|_{w=1})) = 2N.**

**Corollary.** The Liouvillian of the Heisenberg chain with uniform Z-dephasing
γ has exactly 2N purely-real eigenvalues at Re = −2γ for all N.

---

## The chain of reasoning

### Step 1: Purely-real eigenvalues imply commutation

The Liouvillian acts on operators as L(ρ) = −i[H, ρ] + D(ρ), where the
Z-dephasing dissipator D is diagonal in the Pauli basis: D(P) = −2wγ · P for
a Pauli string P with XY-weight w.

For a weight-1 eigenoperator v with eigenvalue λ:

L(v) = −i[H, v] − 2γ · v = λ · v

If λ is purely real (Im = 0), then λ = −2γ and:

−i[H, v] = 0, hence [H, v] = 0.

The 2N modes at Re = −2γ with Im = 0 are exactly the weight-1 operators
that commute with the Hamiltonian.

### Step 2: The SWAP structure of the Heisenberg Hamiltonian

The isotropic Heisenberg coupling between sites i and j is:

X_i X_j + Y_i Y_j + Z_i Z_j = 2 · SWAP_{ij} − I

where SWAP_{ij} exchanges the quantum states of sites i and j. The Hamiltonian
is therefore H = J Σ_{bonds} (2 · SWAP_{ij} − I).

### Step 3: SWAP preserves Z-count

A weight-1 Pauli string has one "active" site (X or Y) and N−1 "passive" sites
(I or Z). The Z-count c is the number of passive sites carrying Z.

When SWAP_{ij} acts on such a string by conjugation (permuting site labels),
it maps the string to another weight-1 string with the same active type and
the same Z-count:

| Original site i | Original site j | After SWAP |
|---|---|---|
| I | I | I, I (unchanged) |
| Z | Z | Z, Z (unchanged) |
| I | Z | Z, I (Z moves, count preserved) |
| X/Y | I | I, X/Y (active site moves) |
| X/Y | Z | Z, X/Y (active moves, Z moves) |

In every case, the Z-count c is preserved.

### Step 4: T_c is the unique SWAP-invariant vector

Fix a type a ∈ {X, Y} and Z-count c. The weight-1 Pauli strings of type a
with Z-count c form a set of size N · C(N−1, c). The symmetric group S_N acts
on this set by permuting site labels, and this action is transitive (any two
such strings are related by a permutation).

For the chain, nearest-neighbor SWAPs generate the full symmetric group S_N.
The unique (up to scalar) S_N-invariant vector on a transitive orbit is the
uniform sum: T_c^{(a)}.

Since [H, T_c^{(a)}] = 0 follows from SWAP-invariance, and there is exactly
one invariant vector per (type, Z-count) pair, we get exactly
2 × N = 2N kernel vectors.

### Step 5: The upper bound -- no additional kernel elements (triangle inequality)

Steps 1–4 prove dim(ker) ≥ 2N by constructing 2N independent kernel vectors.
We now prove dim(ker) ≤ 2N, completing the proof.

Let v ∈ ker([H, ·]|_{w=1}). The Hamiltonian for the open chain is:

H = J Σ_{i=1}^{N−1} (2·SWAP_{i,i+1} − I)

so [H, v] = 0 becomes:

Σ_{i=1}^{N−1} (SWAP_{i,i+1}(v) − v) = 0

where SWAP_{i,i+1}(v) denotes the adjoint action (conjugation: SWAP · v · SWAP†, which permutes the sites inside the operator v). Equivalently:

Σ_{i=1}^{N−1} SWAP_{i,i+1}(v) = (N−1) · v

Now apply the triangle inequality in the Hilbert-Schmidt norm. Each SWAP is
unitary under adjoint action, so ||SWAP_{i,i+1}(v)|| = ||v|| for all i. Thus:

||(N−1) · v|| = ||Σ_{i=1}^{N−1} SWAP_{i,i+1}(v)|| ≤ Σ_{i=1}^{N−1} ||SWAP_{i,i+1}(v)|| = (N−1) · ||v||

This is equality in the triangle inequality. Equality holds if and only if all
vectors in the sum are parallel, i.e., proportional to the same direction.
Since each SWAP_{i,i+1}(v) has the same norm as v, this requires:

**SWAP_{i,i+1}(v) = v  for every bond (i, i+1).**

The sum condition forces each individual SWAP to fix v. This is strictly
stronger than the sum being zero.

For the open chain, the nearest-neighbor transpositions {(1,2), (2,3), ...,
(N−1,N)} generate the full symmetric group S_N (this is a standard result
in group theory: adjacent transpositions generate S_N). Therefore v must be
invariant under ALL permutations:

π(v) = v  for all π ∈ S_N.

On the weight-1 sector, S_N acts on Pauli strings by permuting site labels.
The orbits of this action are indexed by (type a ∈ {X,Y}, Z-count c ∈
{0,...,N−1}), and each orbit is transitive. On a transitive S_N-orbit, the
space of S_N-invariant vectors is exactly one-dimensional: spanned by the
uniform sum over the orbit.

There are 2 types × N Z-counts = 2N transitive orbits, each contributing
one invariant vector T_c^{(a)}. Therefore:

dim(ker([H, ·]|_{w=1})) ≤ 2N.

Combined with Step 4 (dim ≥ 2N), this gives **dim(ker) = 2N exactly**.

**Remark on graph generality.** The triangle inequality argument requires that
the SWAP operators on the bond set generate a transitive group action on each
orbit. For the chain, this gives the full S_N. For any connected graph, the
SWAP operators also generate S_N (since transpositions corresponding to edges
of a connected graph generate S_N). The proof therefore applies to all
connected graphs without modification.

### Step 6: Numerical verification

Independent numerical verification confirms every step of the proof for
N = 2, ..., 7. Five tests were run ([`simulations/verify_triangle_inequality.py`](../../simulations/verify_triangle_inequality.py),
results in [`simulations/results/verify_proof_weight1.txt`](../../simulations/results/verify_proof_weight1.txt)):

**V1 -- Each individual SWAP fixes each kernel vector (Step 5 core claim):**
2N × (N−1) checks per N. Max deviation: 3.84 × 10⁻¹⁵. All PASS.

**V2 -- Triangle inequality saturation:**
Kernel vectors saturate the inequality (lhs = rhs). Non-kernel vectors show
strict inequality (gap grows with N: 0.74 at N=3 to 2.98 at N=7). PASS.

**V3 -- No additional SWAP-invariant vectors (independent upper bound check):**
The stacked (I − SWAP) matrix over all bonds has invariant subspace dimension
exactly 2N at every N. This confirms the upper bound without the triangle
inequality argument. PASS.

| N | w=1 dim | Rank of [H,·] | Kernel dim | SWAP-inv dim | Expected | Status |
|---|---------|---------------|------------|--------------|----------|--------|
| 2 | 8 | 4 | 4 | 4 | 4 | ✓ |
| 3 | 24 | 18 | 6 | 6 | 6 | ✓ |
| 4 | 64 | 56 | 8 | 8 | 8 | ✓ |
| 5 | 160 | 150 | 10 | 10 | 10 | ✓ |
| 6 | 384 | 372 | 12 | 12 | 12 | ✓ |
| 7 | 896 | 882 | 14 | 14 | 14 | ✓ |

**V4 -- Orbit transitivity:** All (type, Z-count) orbits are transitive under
nearest-neighbor SWAPs for N = 2, ..., 7. PASS.

**V5 -- Analytical basis matches numerical eigenvectors:**
The analytically constructed T_c^{(a)} operators span the same subspace as the
numerically computed Liouvillian eigenvectors from `dotnet run -- eigvec`.
rank(T_c ∪ eigvecs) = 2N for all N = 2, ..., 6. PASS.

---

## The 2N conserved operators, explicitly

### c = 0: Global transverse spin

T_0^{(X)} = Σ_j X_j ⊗ I_{rest} = 2 S_x

T_0^{(Y)} = Σ_j Y_j ⊗ I_{rest} = 2 S_y

These are the well-known SU(2) conserved quantities. They commute with any
SU(2)-invariant Hamiltonian.

### c = 1: Single-Z-dressed transverse spin

T_1^{(X)} = Σ_j Σ_{k≠j} X_j Z_k ⊗ I_{rest}

This sums over all N(N−1) strings with X at one site and Z at exactly one
other site. Equivalently: T_1^{(X)} = Σ_j X_j · (Σ_{k≠j} Z_k) = Σ_j X_j · (M_z − Z_j)
where M_z = Σ_k Z_k is the total magnetization operator.

### c = N−1: Fully Z-dressed transverse spin

T_{N-1}^{(X)} = Σ_j X_j ⊗ Z_{all others}

This is the "Jordan-Wigner-type" operator (named after the transformation that maps spins to fermions by attaching a string of Z's) where X acts on one site and every other site carries a Z string.

### General c

T_c^{(a)} = Σ_j σ_a^{(j)} · e_c(Z_1, ..., Z_{j-1}, Z_{j+1}, ..., Z_N)

where e_c is the c-th elementary symmetric polynomial (the sum of all products of c variables chosen from the set: e_1 = sum, e_2 = sum of pairs, etc.) applied to the Z
operators at the passive sites (interpreting Z as a variable and I as its
absence).

---

## Why the proof applies beyond the chain

The SWAP-invariance argument uses only two properties:

1. H is a sum of SWAP operators (isotropic Heisenberg coupling)
2. The SWAP operators on the bond set generate a transitive action on each
   (type, Z-count) orbit

Property 1 holds for any isotropic Heisenberg model. Property 2 holds whenever
the bond graph is connected (since SWAP operators generate the symmetric group
if and only if the graph is connected).

**Corollary.** For any connected graph with N sites and isotropic Heisenberg
coupling, the Liouvillian with uniform Z-dephasing has exactly 2N purely-real
eigenvalues at Re = −2γ.

This extends the d(−γ) = 2N result from the chain to star, ring, complete,
binary tree, and any other connected topology.

---

## What breaks the result

The proof relies critically on isotropy (Δ = 1 in XXZ notation). For the
anisotropic Heisenberg model H = Σ (X_i X_j + Y_i Y_j + Δ Z_i Z_j) with
Δ ≠ 1, the ZZ term does not have the SWAP structure:

X_i X_j + Y_i Y_j + Δ Z_i Z_j = (2 − Δ) SWAP_{ij} + (Δ − 1) Z_i Z_j + const

The residual Z_i Z_j term does not preserve the Pauli type (X ↔ Y mixing via
[Z_i, X_i] = 2iY_i), so T_c^{(X)} and T_c^{(Y)} are no longer individually
conserved. The d(−γ) = 2N formula may not hold for Δ ≠ 1.

---

## Connection to the palindrome

The Π conjugation maps weight w ↔ N−w. The 2N modes at w = 1 (Re = −2γ) pair
with 2N modes at w = N−1 (Re = −2(N−1)γ), giving the palindromic edge:

d_real(1) = d_real(N−1) = 2N

This is [Result 2 of DEGENERACY_PALINDROME](../../experiments/DEGENERACY_PALINDROME.md).
The present document provides the analytical explanation for the left side;
the palindrome provides the right side for free.

---

## Open questions

1. **Anisotropic case.** For XXZ with Δ ≠ 1, how does d(−γ) depend on Δ?
   The SWAP argument breaks, but a modified argument might work for specific
   Δ values.

2. **Higher grid positions: topology-dependent (resolved).**
   The d_real(2) sequence differs between topologies (Chain=14, Star=16,
   Complete=36 at N=4). The weight-2 kernel vectors transform under mixed
   (non-trivial) S_N representations -- the triangle inequality argument
   from Step 5 does not apply. A universal formula for d_real(k≥2) does
   not exist; any formula must incorporate the graph's bond structure.
   See [Weight-2 Kernel](../../experiments/WEIGHT2_KERNEL.md).

3. **Topology dependence: confirmed (resolved).** d_real(k) depends on
   the graph structure for k ≥ 2. Tested on Chain, Star, Ring, Complete
   for N = 3, ..., 6. Only k = 0 (N+1) and k = 1 (2N) are universal.

---

## Reproduction

- Commutator kernel analysis: `python` [`simulations/kommutator_kern_analysis.py`](../../simulations/kommutator_kern_analysis.py)
- Triangle inequality verification: `python` [`simulations/verify_triangle_inequality.py`](../../simulations/verify_triangle_inequality.py)
  (results: [`simulations/results/verify_proof_weight1.txt`](../../simulations/results/verify_proof_weight1.txt))
- Eigenvector export + Pauli projection: `dotnet run -c Release -- eigvec`
  (results: `simulations/results/eigvec_at_minus_gamma_N{2..6}.csv`)
- Eigenvalue CSVs: `simulations/results/rmt_eigenvalues_N{2..7}.csv`
  (generated by `dotnet run -c Release -- rmt`)
