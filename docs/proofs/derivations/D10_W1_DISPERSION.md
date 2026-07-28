# D10: Weight-1 Dispersion Relation

**Derives:** ω_k = 4J(1 − cos(πk/N)), k = 1, ..., N−1, on the (0,1) coherence block
**From:** Heisenberg Hamiltonian structure + Z-dephasing diagonal in Pauli basis
**Convention:** the Pauli normalisation H = J·Σ_⟨i,j⟩ (X_iX_j + Y_iY_j + Z_iZ_j).
The spin convention H = J·Σ S_i·S_j used by the ring and star Im-max proofs is
this one divided by 4, so its frequencies are a quarter of these.
**Status:** PROVEN

---

## Statement

For the N-qubit Heisenberg XXX chain with uniform Z-dephasing (rate γ
per qubit), the Liouvillian restricted to the (0,1) coherence block,
spanned by the |0⟩⟨j| between the all-up state and the single-excitation
states, carries N modes: one at zero frequency and exactly N−1 oscillating
ones, with frequencies

    ω_k = 4J(1 − cos(πk/N)),    k = 1, ..., N−1

and uniform decay rate d = 2γ. The eigenvectors are Neumann standing waves
with amplitudes proportional to cos(πk(j − ½)/N) at site j.

This block sits inside the XY-weight-1 sector but is not all of it; Step 6
gives the scope.

## Definitions

**XY-weight-1 sector:** Pauli strings of the form σ₁ ⊗ ... ⊗ σ_N where
exactly one factor is X or Y, and all others are I or Z. There are
2N · 2^(N−1) such strings: choose which site carries X or Y (N choices),
choose X or Y (2 choices), choose I or Z for each remaining site
(2^(N−1) choices).

**Liouvillian in Pauli basis:** L = L_H + L_D where L_H = −i[H, ·] and
L_D = Σ_k γ_k D_{Z_k}. The dissipator is diagonal: L_D(σ) = −2γ · n_XY(σ) · σ,
where n_XY(σ) is the XY-weight. For w=1 strings, L_D(σ) = −2γ · σ.

## Proof

### Step 1: The dissipator is a scalar on the w=1 sector

For any Pauli string σ with XY-weight 1:

    L_D(σ) = −2γ · σ

This follows from the Z-dephasing dissipator eigenvalue formula
(proven in the Mirror Symmetry Proof, Step 1): D_{Z_k} acting on σ
contributes −2γ_k if site k has X or Y, and 0 if site k has I or Z.
For uniform γ and exactly one XY site, the total is −2γ.

Therefore L_D restricted to the w=1 sector is −2γ · I. All w=1 modes
decay at the same rate 2γ. The frequencies are determined entirely
by L_H.

### Step 2: what is invariant is the popcount grading, not the XY weight

The XY weight is the right label for the dissipator, by Step 1, but it is
NOT conserved by L_H, so it cannot by itself carve out an invariant block.
One counterexample settles it. Take σ = X₀Z₁ at N=3, which has XY weight 1,
and the bond term X₁X₂:

    [X₁X₂, X₀Z₁] = X₀ ⊗ [X₁, Z₁] ⊗ X₂ = −2i · X₀Y₁X₂,

a string of XY weight 3. The leak is not marginal: the largest matrix
element of L_H out of the w=1 sector is 2.0 at N = 3, 4 and 5. What Step 3
below does hold for is the bare string X_j with identities elsewhere, and
that case is worth writing out, but the Z-dressed strings of the same
sector do not follow it.

What IS conserved is the total excitation number: H = J Σ (XX + YY + ZZ)
commutes with Σ_l Z_l, so every coherence block between the popcount-p and
popcount-q sectors is invariant under L_H, exactly (verified: the largest
element of H between different popcounts is 0). The dephasing is diagonal
on that grading too. The block this derivation uses is (p, q) = (0, 1).

For orientation, here is the bare-string commutator algebra that gives the
hopping amplitude. Consider a basis string with X at site j and I elsewhere
(the Y case is analogous), for the Heisenberg Hamiltonian
H = J Σ_{⟨i,j⟩} (X_i X_j + Y_i Y_j + Z_i Z_j).

The key commutators of Pauli matrices are:

    [X, Y] = 2iZ,  [Y, Z] = 2iX,  [Z, X] = 2iY

For the bond (j, j+1) with X at site j and I at site j+1:

    [X_j X_{j+1}, X_j ⊗ I_{j+1}] = [X_j, X_j] ⊗ X_{j+1} = 0
    [Y_j Y_{j+1}, X_j ⊗ I_{j+1}] = [Y_j, X_j] ⊗ Y_{j+1} = −2iZ_j ⊗ Y_{j+1}
    [Z_j Z_{j+1}, X_j ⊗ I_{j+1}] = [Z_j, X_j] ⊗ Z_{j+1} = 2iY_j ⊗ Z_{j+1}

Similarly, for the bond (j−1, j):

    [X_{j-1} X_j, I_{j-1} ⊗ X_j] = X_{j-1} ⊗ [X_j, X_j] = 0
    [Y_{j-1} Y_j, I_{j-1} ⊗ X_j] = Y_{j-1} ⊗ [Y_j, X_j] = −2i Y_{j-1} ⊗ Z_j
    [Z_{j-1} Z_j, I_{j-1} ⊗ X_j] = Z_{j-1} ⊗ [Z_j, X_j] = 2i Z_{j-1} ⊗ Y_j

So [H, ·] on the bare string X_j moves the excitation to the neighbouring
sites j±1 and dresses the vacated site with a Z: the outputs Z_j ⊗ Y_{j+1}
and Y_{j-1} ⊗ Z_j still have XY weight 1. This is the hopping that Step 3
recovers on the coherence block, with amplitude 2J per bond.

It is exactly the Z-dressed strings, whose weight the next commutator does
raise, that make the XY-weight sector the wrong invariant object. The
popcount block is the right one, and it is where the derivation continues.

### Step 3: On the (0,1) block the generator is the graph Laplacian

Restrict to the coherence block spanned by |0⟩⟨j|, where |0⟩ is the
all-up state and |j⟩ carries the single excitation at site j. This is the
repository's own object: the joint-popcount block (p, q) = (0, 1) of
[`compute/MirrorWorld/Block.cs`](../../../compute/MirrorWorld/Block.cs), of
size C(N,0)·C(N,1) = N. It is closed under L, and it is the block the
dispersion speaks about; Step 6 says why the rest of the w=1 sector is not.

Since |0⟩ sits in the popcount-0 sector and |j⟩ in popcount-1, L_H acts
on this block through the one-magnon Hamiltonian alone. In the site basis,
the XXX chain gives

    H|j⟩ = 2J (|j−1⟩ + |j+1⟩) + J [(N−1) − 2·deg(j)] |j⟩,

the hopping 2J from the XX + YY terms of each incident bond, and the
diagonal from the ZZ terms: a bond not touching j contributes +J, a bond
touching j contributes −J, and site j has deg(j) incident bonds. The
all-up state has H|0⟩ = J(N−1)|0⟩, so the oscillation frequency of
|0⟩⟨j| is the gap below the ferromagnet:

    E_ferro·Id − H₁  =  2J · [ deg(j)·δ_{jj'} − A_{jj'} ]  =  2J · 𝓛,

where A is the adjacency matrix of the chain and 𝓛 = D − A its graph
Laplacian. The generator on the block is therefore

    L|₍₀,₁₎  =  −2iJ · 𝓛  −  2γ · Id,

the −2γ being Step 1: |0⟩ and |j⟩ differ in exactly one bit.

**The boundary condition is Neumann, not Dirichlet.** The end sites are
not places where the amplitude is forced to vanish; they are sites with
one neighbour instead of two, so the Laplacian diagonal reads
(1, 2, ..., 2, 1), not 2 throughout. That is what puts N rather than N+1
in the denominator below.

### Step 4: Exactly N−1 of the N modes oscillate

The path graph on N vertices has Laplacian spectrum

    λ_k = 2 − 2cos(πk/N),    k = 0, 1, ..., N−1,

so the block carries N modes, all decaying at the same rate 2γ, with
frequencies ω_k = 2J·λ_k = 4J(1 − cos(πk/N)).

**Why exactly N−1 oscillate.** A graph Laplacian is positive semidefinite
and satisfies 𝓛·𝟙 = 0. For a CONNECTED graph its kernel is exactly the
constants, hence one-dimensional; the chain is connected, so λ_0 = 0 is
simple and every other λ_k is strictly positive. Exactly one mode sits at
zero frequency and the remaining N−1 oscillate:

    ω_k = 4J(1 − cos(πk/N)),    k = 1, ..., N−1.    ∎

The kernel mode is not an accident of the algebra. Its eigenvector is the
uniform one, Σ_j |0⟩⟨j|, which is the total lowering operator S⁻ applied
to the ferromagnet; SU(2) makes it degenerate with |0⟩ itself, so it
cannot oscillate. The count N−1 is the connectivity of the chain speaking.

### Step 5: Eigenvectors are Neumann standing waves

The eigenvectors of the path Laplacian are

    v_k(j) = cos( πk·(j − ½) / N ),    j = 1, ..., N,   k = 0, ..., N−1,

the Neumann standing waves of Step 3, with antinodes at the two open ends
and k nodes in the interior. Mode 0 is the constant, the kernel of Step 4.

The Dirichlet family sin(πkj/N) is NOT the eigenbasis here: it belongs to
the boundary condition that forces the amplitude to vanish beyond the ends,
which an open chain does not impose.

### Step 6: Scope, and the rest of the XY-weight-1 sector

The Pauli w=1 sector of the Definitions is much larger than this block: it
is 2N·2^(N−1)-dimensional and holds the coherences between popcount p and
p+1 for EVERY p, not only p = 0. Those higher blocks carry frequencies of
their own, which are not given by ω_k. At N=5 the full w=1 sector shows 21
distinct oscillation frequencies while the (0,1) block shows 5, four of
them nonzero. This derivation is the (0,1) statement; it does not describe
the XY-weight-1 sector as a whole.

## Corollaries

**F7** (Q-factor spectrum): Q_k = ω_k / d = 4J(1−cos(πk/N)) / (2γ)
= 2J/γ · (1 − cos(πk/N)).

**F41** (Palindromic time): t_Π = 2π/ω_min = 2π/(4J(1−cos(π/N)))
= π/(4J sin²(π/(2N))), since 1−cos(π/N) = 2 sin²(π/(2N)) and so
ω_min = 8J sin²(π/(2N)).

**D01** (Bandwidth): BW = ω_max − ω_min = 4J(cos(π/N) − cos(π(N−1)/N))
= 8J cos(π/N).

**D07** (Q distribution): follows from the density of cosine-spaced
eigenvalues, yielding an arcsine distribution.

## Numerical Verification

Verified to machine precision (relative error < 10⁻¹²) for N = 2
through N = 7 by comparison with full Liouvillian eigendecomposition.
CV (coefficient of variation) of the residuals: 0.000 for all tested N.

See: [`simulations/analytical_spectrum_verify.py`](../../../simulations/analytical_spectrum_verify.py)
