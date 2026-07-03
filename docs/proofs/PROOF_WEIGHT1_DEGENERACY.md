# Why d(−γ) = 2N: The Commutator Kernel in the Weight-1 Sector

<!-- Keywords: commutator kernel, Heisenberg Hamiltonian weight-1 sector, Liouvillian
degeneracy proof, SWAP invariance Pauli strings, Z-count conserved operators,
palindromic degeneracy analytical derivation, d(-gamma) = 2N proof, open quantum
system spectral structure, SU(2) symmetry weight-1, R=CPsi2 commutator kernel -->

**Status:** Lower bound `dim(ker) ≥ 2N` Tier 1 derived for any connected graph (SWAP-invariant construction, Steps 1-4). Upper bound `dim(ker) ≤ 2N` holds for the chain and every tested connected graph EXCEPT K_3 at N=3 (where it is 2N+2); the original triangle-inequality upper-bound argument (Step 5) has a matrix-commutator vs conjugation-action gap, fully characterized in the Appendix. Tier 2 verified bit-exact for chain N = 2 through 7.
**Date:** 2026-04-03
**Reviewed:** 2026-07-03 (readability + numbers + em-dash pass: Δ-decomposition formula corrected, header / Step 5 / beyond-chain claims reconciled with the K_3 N=3 Appendix finding)
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
orbit. This establishes, for the chain:

**dim(ker([H, ·]|_{w=1})) = 2N.**

(Caveat: the triangle-inequality step as originally written has a
matrix-commutator vs conjugation-action gap. It gives the correct answer for
the chain at every tested N but fails for K_3 at N=3, where the kernel is 2N+2.
See Step 5 and the Appendix.)

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

### Step 5: The upper bound, no additional kernel elements (triangle inequality)

Steps 1–4 prove dim(ker) ≥ 2N by constructing 2N independent kernel vectors.
We now prove dim(ker) ≤ 2N, completing the proof.

> **Caveat (added 2026-05-17).** The argument below identifies `[H, v] = 0` with
> the conjugation-action fixed-point equation `Σ_b (SWAP_b(v) − v) = 0`. This
> identification is not correct in general: `[H, v] = 0` reduces to the
> matrix-commutator sum `Σ_b [SWAP_b, v] = 0`, and the two coincide for the chain
> at every tested N but diverge for K_3 at N=3, where the true kernel is 2N+2.
> The conclusion `dim(ker) = 2N` is correct for the chain; the derivation is
> retained here as originally written, with the gap and its resolution documented
> in the Appendix (§ "Where Step 5 of the original proof breaks down").

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

**V1, each individual SWAP fixes each kernel vector (Step 5 core claim):**
2N × (N−1) checks per N. Max deviation: 3.84 × 10⁻¹⁵. All PASS.

**V2, triangle inequality saturation:**
Kernel vectors saturate the inequality (lhs = rhs). Non-kernel vectors show
strict inequality (gap grows with N: 0.74 at N=3 to 2.98 at N=7). PASS.

**V3, no additional SWAP-invariant vectors (independent upper bound check):**
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

**V4, orbit transitivity:** All (type, Z-count) orbits are transitive under
nearest-neighbor SWAPs for N = 2, ..., 7. PASS.

**V5, analytical basis matches numerical eigenvectors:**
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
eigenvalues at Re = −2γ, with the single tested exception of K_3 at N=3 (the
triangle = ring = complete graph on 3 vertices), which has 2N+2. The Appendix
shows K_3 N=3 to be the small-N face of a general central-weight excess in
high-symmetry topologies; the "exactly" here inherits the Step 5 gap, while the
lower bound "≥ 2N" is unconditional.

This extends the d(−γ) = 2N result from the chain to star, ring, complete,
binary tree, and any other connected topology (K_3 N=3 excepted).

---

## What breaks the result

The proof relies critically on isotropy (Δ = 1 in XXZ notation). For the
anisotropic Heisenberg model H = Σ (X_i X_j + Y_i Y_j + Δ Z_i Z_j) with
Δ ≠ 1, the ZZ term does not have the SWAP structure:

X_i X_j + Y_i Y_j + Δ Z_i Z_j = 2 · SWAP_{ij} + (Δ − 1) Z_i Z_j − I

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

The only original open question that remains is the **anisotropic case**: for
XXZ with Δ ≠ 1, how does d(−γ) depend on Δ? The SWAP argument breaks, but a
modified argument might work for specific Δ values.

The topology-dependence questions for higher grid positions k ≥ 2 are resolved:
d_real(k) depends on the graph's bond structure for k ≥ 2 (only k = 0 (N+1)
and k = 1 (2N) are universal). See [Weight-2 Kernel](../../experiments/WEIGHT2_KERNEL.md).

The K_3 N=3 anomaly + max-spin and sub-max closed-form results from the
2026-05 work, plus the residual open questions on sub-max non-sign-rep
contributions and graph extensions, are catalogued in the Appendix below
(§ Open questions, refined post-resolution).

---

## Reproduction

- Commutator kernel analysis: `python` [`simulations/kommutator_kern_analysis.py`](../../simulations/kommutator_kern_analysis.py)
- Triangle inequality verification: `python` [`simulations/verify_triangle_inequality.py`](../../simulations/verify_triangle_inequality.py)
  (results: [`simulations/results/verify_proof_weight1.txt`](../../simulations/results/verify_proof_weight1.txt))
- Eigenvector export + Pauli projection: `dotnet run -c Release -- eigvec`
  (results: `simulations/results/eigvec_at_minus_gamma_N{2..6}.csv`)
- Eigenvalue CSVs: `simulations/results/rmt_eigenvalues_N{2..7}.csv`
  (generated by `dotnet run -c Release -- rmt`)
- K_3 N=3 anomaly investigation: [`simulations/f50_topology_anomaly_sweep.py`](../../simulations/f50_topology_anomaly_sweep.py)
- Per-weight breakdown across topologies: [`simulations/f50_per_weight_breakdown.py`](../../simulations/f50_per_weight_breakdown.py)
- Spin-isotypic decomposition: [`simulations/f50_spin_isotypic_decomposition.py`](../../simulations/f50_spin_isotypic_decomposition.py)
- Max-spin closed-form via Dicke endpoint ladder rungs: [`simulations/f50_max_spin_closed_form.py`](../../simulations/f50_max_spin_closed_form.py)
- Sub-max sign-rep antisymmetric Pauli orbits: [`simulations/f50_sub_max_antisym.py`](../../simulations/f50_sub_max_antisym.py)
- Native C# verification + anomaly test: `compute/RCPsiSquared.Core.Tests/Symmetry/F50NativeEigenvalueCountTests.cs`

---

## Appendix: Step-5 derivation gap and the K_3 N=3 anomaly

### What we found

Native C# verification of F50 across topologies × N (sweep at J = γ = 1, counting pure-real Liouvillian eigenvalues at Re = -2γ, fully consistent across J ∈ [0.01, 5.0]):

| Graph at given N | bonds | count at -2γ | 2N | diff |
|---|---|---|---|---|
| Chain N=2..5 | 1, 2, 3, 4 | 4, 6, 8, 10 | 4, 6, 8, 10 | 0 |
| Ring C_n at n=4, 5 | 4, 5 | 8, 10 | 8, 10 | 0 |
| Star K_{1,n−1} at n=3, 4, 5 | 2, 3, 4 | 6, 8, 10 | 6, 8, 10 | 0 |
| Complete K_n at n=4, 5 | 6, 10 | 8, 10 | 8, 10 | 0 |
| paw, bowtie, book (N=4, 5 with triangles) | 4, 6, 7 | 8, 10, 10 | 8, 10, 10 | 0 |
| **N=3 K_3 (= ring = triangle = complete on 3 vertices)** | **3** | **8** | **6** | **+2** |

Every connected-graph + N combination tested gives `d_real = 2N` EXCEPT the single case **N=3 K_3**, which gives 8 = 2N + 2.

The two K_3-specific extras are weight-1 operators that commute with H_K_3 globally (verified to machine precision: `||[H_K_3, A]|| < 10^{-14}` for both extras) but do NOT commute with H_chain alone (where `||[H_chain, A]|| = 1.0` for the same operators). They survive only when all three K_3 bonds are present: the bond commutators cancel pairwise.

### Where Step 5 of the original proof breaks down

The proof writes (for chain at N):

> `H = J Σ_{i=1}^{N−1} (2·SWAP_{i,i+1} − I)`, so `[H, v] = 0` becomes `Σ_{i=1}^{N−1} (SWAP_{i,i+1}(v) − v) = 0`, where `SWAP_{i,i+1}(v)` denotes the adjoint action (conjugation: `SWAP · v · SWAP^dagger`).

This identification is **incorrect**. The matrix commutator and the conjugation action are different:

    H_b v − v H_b = (2·SWAP_b − I) v − v (2·SWAP_b − I)
                  = 2·(SWAP_b · v − v · SWAP_b)
                  = 2·[SWAP_b, v]    (matrix commutator)

But the conjugation action `SWAP_b(v) = SWAP_b · v · SWAP_b^dagger`, and

    SWAP_b · v · SWAP_b^dagger − v   ≠   SWAP_b · v − v · SWAP_b   in general.

So `[H, v] = 0` reduces to `Σ_b [SWAP_b, v] = 0` (matrix-commutator sum), not to `Σ_b (SWAP_b(v) − v) = 0` (conjugation-action sum). The triangle inequality argument as written then doesn't apply.

### Why the empirical 2N still holds for chain (and most graphs)

For the chain, the conclusion is empirically correct (V1, V3, V5 in Step 6 all pass at N=2..7), but the proof's derivation chain has a gap. The 6 chain N=3 SWAP-invariants (`T_c^{(a)}` for `(a, c) ∈ {X, Y} × {0, 1, 2}`) are honestly in `ker[H_chain, ·]`, but additional non-SWAP-invariant operators might satisfy `Σ_b [SWAP_b, v] = 0` without each `[SWAP_b, v] = 0` individually. For chain at N=2..5 such operators don't exist (the constraint is tight); for ring K_3 at N=3 there are 2 of them.

The correct upper-bound statement is:

    dim(ker[H_G, ·] | weight-1) = dim{v in weight-1 sector : Σ_{b in E(G)} [SWAP_b, v] = 0}

which is the dimension of the joint solution space of the matrix-commutator equation, NOT of the conjugation-fixed-point equation. The two coincide for chain (in tested cases) but diverge for K_3 N=3 by 2 dimensions.

### Structural identification of the K_3 extras

The 2 extras live in the 12-dim "weight-1 c=1 sector" (one X or Y, one Z, one I across N=3 sites). Their Pauli decomposition is a mixture of X- and Y-strings with complex phases of the form `α·X-string + β·Y-string` where `(α, β)` realize the 2-dim **standard irreducible representation of S_3** acting on the cyclic site-permutation orbits.

For chain N=3, `Aut(chain) = Z_2` (reflection through site 1), whose irreps are only trivial + sign (both 1-dim); the standard 2-dim rep of S_3 doesn't survive under the symmetry reduction.

For K_3 N=3, `Aut(K_3) = S_3` (full permutation symmetry), and the standard 2-dim irrep contributes 2 extra weight-1 invariants beyond the F50 SWAP-invariant count.

Adding any external bond to K_3 (paw, bowtie, book) breaks the S_3 symmetry; even paw at N=4 (triangle + leaf vertex 3) gives `count = 2N = 8` per the sweep above.

### Matrix-commutator framework: the right structural angle

The right view of `[H, A] = 0`: it is the centralizer condition `A ∈ Centr(H) = {A : HA = AH}`, and `dim(Centr(H)) = Σ_λ m_λ²` where m_λ is the multiplicity of H-eigenvalue λ (this is the standard centralizer-of-a-matrix formula).

For the Heisenberg ring on small N, K_n is fully `S_N`-permutation-symmetric, so H eigenvalues split by total-spin sector. The H spectrum is therefore strongly degenerate. For K_3 N=3 the spectrum is exactly 2 eigenvalues, each with multiplicity 4 (= S = 3/2 multiplet, dim 4; and S = 1/2 multiplet doubled, dim 4):

- K_3 N=3: H spec = (+3J/4 mult 4) ⊕ (−3J/4 mult 4). `dim(Centr) = 16 + 16 = 32`.
- Chain N=3: H spec = (−1 mult 2) ⊕ (0 mult 2) ⊕ (+1/2 mult 4). `dim(Centr) = 4 + 4 + 16 = 24`.

K_3's full centralizer is 8 dimensions larger than chain's. Of these 8 extra centralizer dimensions, how are they distributed across pure-weight sectors? Empirical per-weight breakdown at N=3:

| weight w | chain dim(ker[H, ·]\|_w) | K_3 dim(ker[H, ·]\|_w) | excess |
|---|---|---|---|
| 0 | 4 | 4 | 0 |
| 1 | 6 | 8 | **+2** |
| 2 | 6 | 8 | **+2** |
| 3 | 4 | 4 | 0 |
| Σ (pure-weight) | 20 | 24 | **+4** |
| Centralizer total | 24 | 32 | **+8** |

**Two key structural facts:**

1. **Palindrome pair.** The weight-1 anomaly is paired with a weight-2 anomaly of the same size (+2). This is the F1 palindrome symmetry `w ↔ N − w` (here `N = 3`, so `w = 1 ↔ w = 2`); the 2 K_3 N=3 extras at weight-1 have palindromic partners at weight-2. Together they contribute 4 extra centralizer dimensions in pure-weight sectors.

2. **Multi-weight tail.** The remaining 4 extra centralizer dimensions live in multi-weight operators (operators with components in multiple weight sectors simultaneously). Chain's centralizer is dominated by pure-weight operators (20 of 24 = 83%); K_3's centralizer has 4 multi-weight operators that chain lacks (24 of 32 = 75% pure-weight). The multi-weight mechanism observed in [`WEIGHT2_KERNEL.md`](../../experiments/WEIGHT2_KERNEL.md) for Chain N=4 weight-2 (`d_real(2) = 14 > ker(w=2) = 13`) is the same kind of phenomenon: H's spin symmetry allows operators whose dephasing decay rates average to a specific target via cross-weight cancellation.

The matrix-commutator picture also clarifies why `Aut(G)`-irrep arguments (Schur class-sum scalar = 0) don't directly predict ker contributions: those arguments use the group-algebra LEFT-multiplication structure, while `[H, A] = 0` is matrix commutation. The centralizer `Centr(H)` decomposes as `⊕_λ M(d_λ)` (matrix algebras on each H-eigenspace), and the weight-w intersection with this depends on how Pauli weight aligns with the eigenspace projectors, a graph-and-N specific algebraic question.

### Cross-link to WEIGHT2_KERNEL.md (April 3, 2026)

The empirical K_3 N=3 anomaly was first recorded in [`experiments/WEIGHT2_KERNEL.md`](../../experiments/WEIGHT2_KERNEL.md) (lines 151-155):

> "This phenomenon was also observed at weight 1 for Ring/Complete at N = 3 (d_real(1) = 8 > 2N = 6). It is rare but structurally important: d_real(k) counts ALL eigenvalues at a grid position, not just those from weight-k."

WEIGHT2_KERNEL tentatively attributed the K_3 N=3 weight-1 anomaly to **multi-weight mixing** (by analogy with the chain N=4 weight-2 case at Re = -4γ, where one eigenvalue at that grid position does come from a multi-weight operator). The 2026-05-17 native verification of F50 checks this attribution explicitly: the 2 K_3 N=3 extras are **pure-weight-1** (verified: `|c_α|² in weight-1 sector = 1.000000`, all other weight sectors exactly zero). They are not multi-weight operators. The correct attribution is the **S_3 standard irrep on the weight-1 c=1 sector**, as derived above.

WEIGHT2_KERNEL also introduced a **Trivial / Alternating / Mixed S_N-representation table** (lines 78-87) for weight-2 kernel vectors:

| N | Trivial (+1) | Alternating (−1) | Mixed | Total |
|---|-------------|-----------------|-------|-------|
| 3 | 6 | 0 | 0 | 6 |
| 4 | 0 | 0 | 13 | 13 |
| 5 | 0 | 0 | 14 | 14 |
| 6 | 0 | 0 | 19 | 19 |

Today's K_3 N=3 finding extends the same decomposition format **to weight-1**, with the "Mixed" column refined for S_3 into its precise irreps (Standard = the 2-dim non-trivial non-sign irrep):

| N | topology | Trivial (+1) | Sign (−1) | Standard (2-dim, S_3) | Total |
|---|----------|-------------|-----------|----------------------|-------|
| 3 | Chain | 6 | 0 | 0 | 6 |
| 3 | **K_3 (= ring = complete = triangle)** | **6** | **0** | **2** | **8** |
| 3 | Star (= chain by relabeling) | 6 | 0 | 0 | 6 |
| 4 | Chain, Ring, Star, Complete, paw, K_4 − e | 8 | 0 | 0 | 8 |
| 5 | Chain, Ring, Star, K_5, bowtie, book | 10 | 0 | 0 | 10 |

The Trivial column matches the F50 SWAP-invariant operators (T_c^{(a)} for `a ∈ {X, Y}, c ∈ {0, ..., N−1}`). At weight-1, the Sign and Standard columns are typically empty; **K_3 N=3 is the unique tested exception**, where the S_3 standard 2-dim irrep contributes 2 invariants.

Connecting to WEIGHT2_KERNEL's weight-2 table: there, "Mixed" at N ≥ 4 captures non-trivial reps of S_N (standard 3-dim of S_4, etc.). At weight-1, the analogous non-trivial reps almost never produce kernel vectors (the dispersion structure rules them out), with K_3 N=3 as the unique witness.

### Resolution of the open question (full topology sweep)

Extending the per-weight breakdown to chain/ring/star/K_4-e/K_4 at N=4 and chain/ring/K_5 at N=5 surfaces a clean universal pattern:

| (graph, N) | per-weight ker (w=0..N) | central-w excess vs chain | palindrome location |
|---|---|---|---|
| chain N=3 | (4, 6, 6, 4) | baseline | – |
| **K_3 N=3** | (4, 8, 8, 4) | **+2 at w=1, +2 at w=2** | pair around N/2 = 1.5 |
| chain N=4 | (5, 8, 13, 8, 5) | baseline | – |
| star N=4 | (5, 8, 16, 8, 5) | +3 at w=2 | self-palindromic at N/2 = 2 |
| K_4-e N=4 | (5, 8, 22, 8, 5) | +9 at w=2 | self-palindromic |
| ring N=4 | (5, 8, 23, 8, 5) | +10 at w=2 | self-palindromic |
| **K_4 N=4** | (5, 8, 36, 8, 5) | **+23 at w=2** | self-palindromic |
| chain N=5 | (6, 10, 14, 14, 10, 6) | baseline | – |
| ring N=5 | (6, 10, 22, 22, 10, 6) | +8 at w=2, +8 at w=3 | pair around N/2 = 2.5 |
| **K_5 N=5** | (6, 10, 54, 54, 10, 6) | **+40 at w=2, +40 at w=3** | pair around N/2 = 2.5 |

**The pattern:** every connected graph with non-trivial automorphism beyond chain shows centralizer excess at the **central weights** `w ∈ {floor(N/2), ceil(N/2)}`, palindromic pair when N is odd, self-palindromic single value when N is even. The excess magnitude is topology-dependent (K_N largest among K_N / ring / star / K_N − e). At edge weights w ∈ {0, 1, N-1, N} the count is topology-independent (matches chain = 2N at w=1 and 2 at w=0 etc.).

**Why K_3 N=3 surfaces as a "weight-1 anomaly":** at N=3, floor(N/2) = 1, so the central palindromic pair is (w=1, w=2). F50 specifically tracks weight-1, so the K_3 excess shows up there. For N ≥ 4 the central weight is ≥ 2 and F50's weight-1 count remains 2N for all topologies tested. **K_3 N=3 is not a special algebraic phenomenon; it is the small-N manifestation of the universal "central-weight excess in high-symmetry topologies" pattern.**

**Why the excess is palindromic:** the F1 Π-conjugation palindrome `d_real(w) = d_real(N − w)` (proven for the full Liouvillian spectrum) forces the per-weight excess to be palindromic too. The conjugation pairing weight-w ↔ weight-(N−w) commutes with H's action on operators (via the standard Π = Z⊗N involution), so any centralizer dimension at weight w has a partner at weight N−w.

**Connection to the WEIGHT2_KERNEL.md observations:** the +23 at K_4 weight-2 and similar topology-dependent counts at weight-2 across N=4..6 (their original table) are the SAME phenomenon as today's K_3 N=3 weight-1 finding, different value of N, same "central weight excess + F1 palindrome". WEIGHT2_KERNEL had documented the topology dependence at weight-2 ≥ 4 weeks ago; we now understand the K_3 N=3 case as the same pattern with N=3's central weight happening to land at w=1.

**What remains genuinely open:**
- A closed-form formula for the excess as a function of (graph G, N, weight w). The values 2, 3, 9, 10, 23 (N=4 w=2 across topologies) and 8, 40 (N=5 K_n w=2) don't fit an obvious combinatorial family.
- The micro-structure of the central-weight extras: what specific S_N-irrep + spin alignment produces them, and why K_N has the largest count.
- Higher N: does the pattern hold at N ≥ 6? (Untested today but no reason to expect deviation.)

### F94/F96/F97 cross-check (do today's other formulas help?)

Asked: can F94's sym3 closed form (4/3 Q²K³ for |0+0+⟩ pair (0,2)) or F96's universal subdominant slopes (-16/9, -8/3) detect K_4 vs Ring topology dependence?

Test ([`simulations/f94_topology_visibility_probe.py`](../../simulations/f94_topology_visibility_probe.py)): compute F94's sym3 matrix element for all 4 outcomes on |0+0+⟩ pair (0,2) at chain, ring, K_4.

| Outcome | Chain | Ring | K_4 |
|---------|-------|------|-----|
| sym3 \|00⟩ | +5 | +8 | **+8** (= Ring) |
| sym3 \|01⟩ | -4 | -4 | **-4** (= Ring) |
| sym3 \|10⟩ | -1 | -4 | **-4** (= Ring) |
| sym3 \|11⟩ | 0 | 0 | **0** (= Ring) |

K_4 and Ring give bit-identical Dyson matrix elements across all 4 outcomes. F96's slopes are also identical (chain, ring, K_4 all give -16/9 and -8/3 for the |01⟩ and |11⟩ subdominant slopes). F94/F96 are **blind to the K_4 vs Ring topology distinction** at this canonical lens.

**Why the blindness:** the two extra K_4 bonds vs Ring are (0,2) and (1,3), the "diagonal" bonds. For |0+0+⟩ initial state and pair (0,2) measurement:
- Bond (0,2) connects the two kept-pair sites (both prepared as |0⟩, Z eigenstates), so [Z_0 Z_2, |0,·,0,·⟩] = 0 and the (0,2) bond's contributions vanish.
- Bond (1,3) connects the two traced-out sites (both prepared as |+⟩, X eigenstates), so its contributions trace out to zero on the pair (0,2) observable.

Both extra K_4 bonds fall in the **symmetric blind spots** of this specific (state, pair) lens.

For asymmetric initial states, F94's sym3 DOES detect K_4 vs Ring:

| Lens | Ring sym3 | K_4 sym3 | K_4 − Ring |
|------|-----------|----------|------------|
| \|++00⟩ pair (0,2) | (-0.5, -0.5, 1.5, -0.5) | (-0.5, -3.5, 4.5, -0.5) | (0, **-3**, **+3**, 0) |
| \|10+0⟩ pair (0,1) | (-2, -2, 5, -1) | (-3, -2, 6, -1) | (**-1**, 0, **+1**, 0) |

The K_4 − Ring difference is an **antisymmetric integer shift between outcomes differing by one bit flip**. Magnitudes 1-3, much smaller than F50's +23 central-weight excess at K_4 N=4 weight-2.

**Conclusion:** F94/F96 and F50 are **complementary lenses** on the same Heisenberg + Z-deph topology dependence:

| Lens | What it tracks | K_4 vs Ring sensitivity |
|------|----------------|------------------------|
| F94 sym3 (\|0+0+⟩, pair (0,2)) | Specific observable on specific state | Blind (symmetric blind spots) |
| F94 sym3 (\|++00⟩ etc.) | Same observable on asymmetric state | Visible: integer shifts (1-3) |
| F50 weight-1 ker | Full operator algebra, state-independent | Invisible at w=1 for N ≥ 4 (= 2N) |
| F50 weight-w ker (w = central) | Full operator algebra, state-independent | Visible: large excess (+23 at K_4 N=4 w=2) |

No direct closed-form transfer from F94/F96 to F50's central-weight excess. But the structural lesson is consistent: K_4's extra bonds produce **bond-commutator cancellation patterns** that show up as small shifts in observable-specific lenses (F94 with asymmetric states) and as large excesses in operator-algebra-level counts (F50 central weight).

The F94/F96/F97 framework is the right toolkit for tracking specific Heisenberg + Z-deph observables; F50 is the right toolkit for tracking the operator algebra's universal degeneracy structure. Both are needed to fully characterize the Heisenberg + Z-deph topology family.

### Spin-isotypic decomposition: partial closed-form

Decomposing the per-weight ker `dim(ker[H_G, ·] | weight-w)` by spin sector splits the count into **single-block** (operators supported entirely in one H-eigenspace) and **multi-block** (block-diagonal operators that span multiple spin sectors with cross-spin coordination):

    dim(ker[H_K_N, ·] | weight-w) = Σ_S single_block(S, w) + multi_block(w)

For K_N, the spin sectors are indexed by total spin S ∈ {N/2, N/2 − 1, ...}, with block dimension `dim_S = m_S · (2S+1)` where `m_S` is the SU(2) multiplicity (= number of irreducible copies). Empirical decomposition ([`simulations/f50_spin_isotypic_decomposition.py`](../../simulations/f50_spin_isotypic_decomposition.py)):

| Setup | (m_S, 2S+1, dim_S) | (w=0, w=1, ..., w=N) single-block pattern |
|-------|--------------------|-------------------------------------------|
| K_3 S=3/2 max | (1, 4, 4) | (2, 4, 4, 2), sum 12 = 4N |
| K_3 S=1/2 | (2, 2, 4) | (0, 2, 2, 0), central only |
| K_4 S=2 max | (1, 5, 5) | (2, 4, 4, 4, 2), sum 16 = 4N |
| K_4 S=1 | (3, 3, 9) | (0, 0, 26, 0, 0), central only |
| K_4 S=0 | (2, 1, 2) | (0, 0, 1, 0, 0), central only |
| K_5 S=5/2 max | (1, 6, 6) | (2, 4, 4, 4, 4, 2), sum 20 = 4N |
| K_5 S=3/2 | (4, 4, 16) | (0, 0, 22, 22, 0, 0), central pair |
| K_5 S=1/2 | (5, 2, 10) | (0, 0, 8, 8, 0, 0), central pair |
| K_6 S=3 max | (1, 7, 7) | (2, 4, 4, 4, 4, 4, 2), sum 24 = 4N |
| K_6 S=2 | (5, 5, 25) | (0, 0, 38, **0**, 38, 0, 0), even-w only (parity!) |
| K_6 S=1 | (9, 3, 27) | (0, 0, 30, 124, 30, 0, 0), central triple |
| K_6 S=0 | (5, 1, 5) | (0, 0, 0, 0, 0, 0, 0), **vanishes** |

**Two universal structural facts:**

1. **Max-spin block (S = N/2, m_S = 1, dim = N+1) gives the universal palindromic pattern (2, 4, 4, ..., 4, 2)** with `(N − 1)` interior 4s and edge 2s, total sum `4N` for all N ≥ 3. This is the SWAP-invariant `T_c^{(a)}` contribution from the original F50 proof, those 2N operators live entirely in the max-spin (fully symmetric) eigenspace.

2. **Sub-max spin blocks (S < N/2) contribute pure-weight operators ONLY at central weights**, with zero contribution at edge weights. The "central window" width and parity depend on (m_S, 2S+1) in a non-trivial way: at K_5 sub-max contributes at w ∈ {2, 3}; at K_6 S=2 contributes at w ∈ {2, 4} (skipping the true center w=3 due to a parity selection rule); at K_6 S=0 the contribution vanishes entirely.

**Consequence: central-weight excess decomposes cleanly.**

    central-weight-excess(K_N) = Σ_{S < N/2} single_block(S, central w) + multi-block contributions

The max-spin pattern is **N-uniform and palindromic in w**, so it contributes equally to central and neighbor weights, max-spin alone does NOT create central-weight excess. The excess comes entirely from sub-max spin sectors concentrating their pure-weight mass at central weights.

For K_3 N=3 weight-1: single-block excess vs chain = 6 − 4 = 2 (from S=1/2 block contributing 2 at w=1). Multi-block contribution is **identical for K_3 and chain** (both = 2 at w=1). The K_3 N=3 "anomaly" is **entirely** a single-block phenomenon: the K_3 spin-1/2 block has 2 pure-weight-1 operators absent from chain's denser eigenspace structure.

For K_4 N=4 weight-2: single-block contribution = 4 (chain) → 31 (K_4) = +27 single-block excess. Multi-block contribution = 9 (chain) → 5 (K_4) = −4 multi-block diff. Net excess = +23 ✓.

**Closed-form for max-spin contribution (Tier 1 derived):**

For the max-spin block (S = N/2, m_S = 1, dim N+1 = symmetric subspace), the pure-weight-w sym-supported operator space has a complete closed-form, bit-exact for N = 2..5 across all w:

    single_block(S = N/2, w) = { 2  if w = 0 or w = N
                                  4  if 1 ≤ w ≤ N-1 }

    Σ_w single_block(S = N/2, w) = 4N
    Multi-weight in max-spin block = (N+1)² − 4N = (N−1)²

The explicit operators are the **Dicke-endpoint ladder rungs** between |D_0⟩ = |0⟩^N and |D_N⟩ = |1⟩^N:

| w | Basis of pure-weight-w max-spin ops | Dim |
|---|-------------------------------------|-----|
| 0 | `|D_0⟩⟨D_0|, |D_N⟩⟨D_N|` (diagonal endpoint projectors) | 2 |
| 1..N-1 | `|D_0⟩⟨D_w| ± h.c., |D_{N-w}⟩⟨D_N| ± h.c.` (two endpoint-anchored rungs, real + imag) | 4 |
| N | `|D_0⟩⟨D_N| + h.c., i(|D_0⟩⟨D_N| − h.c.)` (full-ladder jump) | 2 |

Explicit identification of the w=0 endpoint projectors (derived from `Π_i (I ± Z_i)`):
- `|D_0⟩⟨D_0| = (1/2^N) · Σ_k e_k(Z_1, ..., Z_N) = (1/2^N) · Π_i (I + Z_i)`
- `|D_N⟩⟨D_N| = (1/2^N) · Σ_k (-1)^k e_k(Z_1, ..., Z_N) = (1/2^N) · Π_i (I − Z_i)`

where e_k is the k-th elementary symmetric polynomial in Z-operators (verified bit-exact for N = 2..5).

**Multi-weight identification:** the (N−1)² multi-weight ops in M(N+1) correspond to the (N−1)² **middle-Dicke transitions** |D_k⟩⟨D_l| for k, l ∈ {1, ..., N-1}. These intrinsically mix multiple Pauli weights when expanded (whereas endpoint-anchored rungs are pure-weight). Count check: 4N pure-weight + (N−1)² multi-weight = N² + 2N + 1 = (N+1)² ✓.

**Structural reason:** the extreme Dicke states |D_0⟩ = |0⟩^N and |D_N⟩ = |1⟩^N are SINGLE computational basis states (not superpositions), so operators involving them have clean Pauli weight expansions. Middle Dicke states |D_k⟩ for 1 ≤ k ≤ N-1 are superpositions of C(N, k) computational states, and off-diagonal ops between middle states generate multi-weight Pauli expansions.

**Consequence for F50:** the max-spin contribution to ker(K_N, w) is weight-uniform (always 2 or 4 per weight, never different at central vs edge weights). Therefore the **central-weight excess of K_N over chain is entirely a sub-max-spin phenomenon**, confirming the empirical observation that the K_N anomaly comes from the lower-spin sectors concentrating their mass at central weights.

**Forward link.** The same X⊗N-eigenbasis / Dicke-endpoint mechanism powers the closed-form proof of the 3/8 K-intermediate anchor for the Dicke superposition (|D_n⟩ + |D_{n+1}⟩)/√2 in F86b, see [`PROOF_F86B_UNIVERSAL_SHAPE.md`](PROOF_F86B_UNIVERSAL_SHAPE.md) §Statement 2 closed-form proof block. Both are instances of the orthogonal-symmetry-decomposition pattern: under the involution X⊗N (or equivalently the popcount-mirror structure), Pauli operators split into ±1 eigenspaces whose pure-weight content has clean algebraic descriptions.

**Partial closed-form for sub-max via antisymmetric Pauli orbits (Tier 1 for small cases):**

For sub-max spin sectors at K_N, the **sign-rep (antisymmetric) contribution** is identified as totally antisymmetric Pauli tensors over distinct-letter multisets:

    sign-rep contribution at K_N (S, w) = # distinct-letter Pauli multisets of weight w on N sites

A "distinct-letter multiset" at weight w on N sites means w letters from {X, Y} and N-w letters from {Z, I}, with ALL N letters distinct (no repeats). Since the letter alphabet has size 4 = |{X, Y, Z, I}|, distinct multisets exist only when w ≤ 2 AND N-w ≤ 2:

    distinct_multisets(N, w) = C(2, w) · C(2, N-w)   [zero if w > 2 or N-w > 2]

| N | w=0 | w=1 | w=2 | w=3 | w=4 |
|---|-----|-----|-----|-----|-----|
| 2 | 1   | 4   | 1   | 0   | 0   |
| 3 | 0   | 2   | 2   | 0   | 0   |
| 4 | 0   | 0   | 1   | 0   | 0   |
| ≥5 | 0  | 0   | 0   | 0   | 0   |

**Verified bit-exact** (the antisymmetric op is supported entirely in the corresponding sub-max block):

- **K_3 S=1/2 pure-weight pattern (0, 2, 2, 0)** = exactly the distinct-multiset count ✓. Explicit ops: antisym(X, Z, I), antisym(Y, Z, I) at w=1; antisym(X, Y, Z), antisym(X, Y, I) at w=2.
- **K_4 S=0 pure-weight pattern (0, 0, 1, 0, 0)** = exactly the distinct-multiset count ✓. The 1 op = `ε_{ijkl} · σ_i σ_j σ_k σ_l` summed over permutations of (I, X, Y, Z) across the 4 sites = totally antisymmetric Levi-Civita Pauli tensor.
- **K_6 S=0 pure-weight pattern (0, 0, 0, 0, 0, 0, 0)** = exactly zero from distinct-multiset count (impossible to have 4+ distinct letters in 6 positions) ✓. Predicts the empirically observed vanishing!

The general pattern is consistent with **S_N representation theory**: the sub-max block algebra `M(V_S^{(m_S)}) = V_SN ⊗ V_SN^*` decomposes under S_N into irreps, and the sign-rep [1^N] slot is filled by antisymmetric Pauli orbits which exist only with distinct letters.

**What the antisym rule DOES NOT capture (closed-form remains incomplete):**

The pattern fails for sub-max sectors at N ≥ 4 with non-zero contributions beyond the sign-rep:

- K_4 S=1 pure-weight-2 = 26 (sign-rep gives 0; the 26 comes from [4], [2,2], [3,1], [2,1,1] irrep contributions in the block algebra)
- K_5 S=1/2 pure-weight = 8 each at w=2, 3 (sign-rep gives 0 for N=5)
- K_5 S=3/2 = 22, K_6 S=1 = 30/124/30, K_6 S=2 = 38/0/38

These come from **non-sign S_N-irrep components** that are non-zero for non-distinct-letter multisets. The general formula:

    single_block(K_N, S, w) = Σ_{ρ ∈ S_N-irreps in V_SN ⊗ V_SN^*}
                                (block-mult of ρ) · (# weight-w Pauli orbits in ρ)

For each S_N-irrep ρ:
  - Block multiplicity of ρ = `dim Hom_{S_N}(V_SN, V_SN ⊗ ρ)` (Schur-style)
  - Weight-w Pauli orbits in ρ = sum over multisets of `dim Hom_{S_N}(orbit_multiset, ρ)`

This is the **Schur-Weyl decomposition on operator space**, evaluating Frobenius reciprocity for each (ρ, multiset) pair. Tractable per case via character tables, but a single uniform closed-form for arbitrary (N, S, w) requires summing many irrep contributions.

**Structural insight:** the K_N central-weight excess over chain is driven by:
1. Antisym (sign-rep) ops: # distinct-letter multisets at weight w (clean closed-form, dominates at small N)
2. Other S_N-irrep contributions: more complex, dominate at N ≥ 5 where distinct-letter count vanishes

For K_4 N=4 w=2: excess = 23 = 1 (sign-rep antisym from S=0) + 22 (other irreps from S=1, S=0, S=2). The "antisym" mechanism is a clean but minority piece of the total excess at K_4.

**Status:** 
- Max-spin closed-form: Tier 1 derived (Dicke endpoint ladder rungs, complete identification).
- Sub-max antisym contribution: Tier 1 derived for sign-rep [1^N] (= distinct-letter multiset count), bit-exact verified for K_3 S=1/2, K_4 S=0, K_6 S=0 = 0.
- Sub-max non-sign-rep contributions: Tier 2 empirical, requires per-case S_N-character calculations on weight-w Pauli orbit decompositions.

A complete closed-form formula `f(N, S, w)` exists in principle via Schur-Weyl + character theory, but writing it as a single closed expression rather than a sum-over-irreps remains open.

### Open questions (refined post-resolution; resolved Q1/Q2/Q4 retired 2026-05-17 evening)

> **Q1 (matrix-commutator vs left-action gap), RESOLVED 2026-05-17 morning.** The gap
> is fully characterized in the appendix above (§ "Where Step 5 of the original
> proof breaks down" and § "Matrix-commutator framework: the right structural
> angle") and in `docs/ANALYTICAL_FORMULAS.md` F50 section. Briefly: `[H, A] = 0`
> reduces to `Σ_b [SWAP_b, A] = 0` (matrix-commutator sum), not to
> `Σ_b (SWAP_b(A) − A) = 0` (conjugation-action sum); the centralizer of H decomposes
> as `⊕_λ M(d_λ)` on each H-eigenspace, and the weight-w intersection is a graph-and-N
> specific algebraic question. Explicit falsification of the naive class-sum conjecture
> at K_4 N=4: `simulations/f50_irrep_decomposition_probe.py`.

> **Q2 (why K_3 N=3 specifically), RESOLVED 2026-05-17 evening.** K_3 N=3 is NOT
> algebraically special. It is the small-N manifestation of the universal
> "central-weight excess in high-symmetry topologies" pattern (this proof §
> "Resolution of the open question (full topology sweep)"). The
> structural cause: sub-max-spin sectors concentrate pure-weight mass at central
> weights only, while max-spin baseline is weight-uniform. N=3's central weight = 1
> coincides with F50's tracked weight, so the excess shows there; N ≥ 4 has central
> weight ≥ 2 and the excess shifts there too (K_4 w=2 excess 23, K_5 w=2 excess 40,
> etc.). Simulation: `f50_spin_isotypic_decomposition.py`, `f50_max_spin_closed_form.py`.

1. **Sub-max-spin single-block closed-form** (Q3, refined). For max-spin `S = N/2`,
   `single_block(S, w) = 2 if w ∈ {0, N}, else 4` is the universal Dicke-endpoint
   ladder-rungs pattern (Tier 1 derived; this proof § "Max-spin closed-form via Dicke
   endpoint ladder rungs"). For sub-max spin `S < N/2`, the corresponding `f(N, S, w)`
   in closed form is open. Empirical patterns:
   - K_3 S=1/2: (0, 2, 2, 0), the antisym piece, now closed-form via Q3' below.
   - K_4 S=0: (0, 0, 1, 0, 0); K_4 S=1: (0, 0, 26, 0, 0); K_6 S=1: (0, 0, 30, 124, 30, 0, 0).
   - K_6 S=0 = (0, 0, 0, 0, 0, 0, 0), predicted to vanish from antisym; confirmed.

2. **Q3' Sub-max non-sign-rep contributions** (partial). The sign-rep [1^N]
   contribution to `single_block(S, w)` is closed-form via distinct-letter Pauli
   multisets: `C(2, w) · C(2, N−w)` (this proof § "Partial closed-form for sub-max
   via antisymmetric Pauli orbits", 2026-05-17 night, Tier 1 derived for small cases).
   For K_4 N=4 w=2 excess = 23 = 1 (sign-rep antisym) + 22 (other irreps from S=1, S=0,
   S=2). The Schur-Weyl per-S_N-irrep decomposition of the 22 piece is open and would
   complete the closed-form via Frobenius reciprocity character calculations.

> **Q4 (does any other high-Aut graph beyond K_3 N=3 anomalize), RESOLVED 2026-05-17
> evening.** The n_XY=1 commutant sweep via efficient rank-of-`[H, ·]`-on-n_XY=1
> subspace (`simulations/f50_weight1_commutant_efficient.py`, bypasses the
> 4^N × 4^N dense Liouvillian via the observation that L A = -2γ A forces
> n_XY(A) = 1, and `[H, ·]` restricted to the n_XY=1 subspace is a 2N · 2^(N−1)
> × O(N²) · 2^(N−1) sparse matrix). Tested high-Aut graphs all show δ = 0:
>
> | graph | N | Aut size | dim ker | δ |
> |-------|---|----------|---------|---|
> | K_{3,3} (bipartite, edge-transitive) | 6 | 72 | 12 | 0 |
> | Q_3 hypercube (edge-transitive) | 8 | 48 | 16 | 0 |
> | Möbius-Kantor M_8 (Cayley Z_8 with {±1, ±3}) | 8 | 96 | 16 | 0 |
> | Petersen (vertex+edge transitive) | 10 | 120 | 20 | 0 |
>
> Consistent with the Q2 resolution: K_3 N=3 is uniquely anomalous because central
> weight = 1 coincides with the F50-tracked weight ONLY at N=3. For N ≥ 4 the
> central weight is ≥ 2 (where the excess does live, e.g. K_4 w=2 excess = 23),
> so weight-1 has no excess regardless of graph automorphism. K_3 N=3 remains the
> unique weight-1 anomaly in this sweep.

4. **N=3 algebraic-forcing meta-pattern (Q5).** Several N=3 specialties are recorded
   in the repo: F33's exact-rational decay rates (chain N=3 only; D10's
   `cos(πk/3) = ±1/2`), F69's irreducible sextic over ℚ (GHZ_3 + W_3 slice),
   `hypotheses/THE_BOOT_SCRIPT.md` OQ-294 (V-Effect 14 broken combos at N=3),
   `hypotheses/NESTED_MIRROR_STRUCTURE.md`'s 12-class refraction at N=3,
   `experiments/PRIMORDIAL_QUBIT_ALGEBRA.md`'s `w_M = 3/4` (N=3-chain only). K_3
   weight-1 extras (now resolved as the small-N face of central-weight excess) are
   another facet of "N=3 is algebraically distinguished". A synthesis tying these
   together is open.

### Status update

The F50 formula `d_real(Re = −2γ) = 2N` should be read with the K_3 N=3 caveat:
- **Tier 1 lower bound `≥ 2N`**: rigorously proven via SWAP-invariant construction (Steps 1-4 of the original proof).
- **Tier 2 upper bound `≤ 2N`**: empirically verified for chain N=2..7 and most other connected graphs at N ≥ 4, but **violated at N=3 K_3 by 2** (one S_3 standard 2-dim irrep contribution). The proof's Step 5 derivation has a matrix-commutator vs conjugation-action gap that explains the missed K_3 case.

