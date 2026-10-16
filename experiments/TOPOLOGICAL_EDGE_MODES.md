# Topological Edge Modes: The Localization Is Geometric

<!-- Keywords: topological edge modes Liouvillian, SSH model palindromic
spectrum, chiral block structure winding number, SSH-type diagnostic
open quantum, standing wave mode localization, sacrifice zone mechanism
geometric not topological, biorthogonal Berry phase Lindbladian,
R=CPsi2 topological edge modes -->

**Status:** Computationally verified (Tier 2). Negative result: not topological.
**Date:** April 1, 2026
**Script:** [topological_edge_modes.py](../simulations/topological_edge_modes.py)
**Data:** [topological_edge_modes.txt](../simulations/results/topological_edge_modes.txt)
**Depends on:**
- [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994, the question)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (Pi gives a sectorwise P relation; full SRP class open)
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (Poisson-like spacing statistics)
- [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) ((0,1) coherence block dispersion)

---

## What this document is about

The cavity mode localization result (r = 0.994) looks like it could be
"topological protection," meaning protected by a mathematical invariant that
cannot change smoothly, the way a coffee mug cannot become a sphere without
tearing a hole. This document runs five checks. Three tests (the SSH
structure, the rank of the chiral blocks, a sweep for a phase boundary)
and the exposed mechanism point the same way; a Berry phase cannot decide
on this kind of path, and a mode count is no invariant either way.
The answer is no: the localization is geometric, not topological. Standing wave
patterns on a chain naturally place some modes at the edges and others at
the center; the noise profile simply selects which subset survives. No
topological invariant is involved. This is a negative result.

---

## Abstract

The cavity mode localization result (r = 0.994, profile [0.52, 0.63,
0.70, 0.63, 0.52] invariant under noise profiles) looks topological.
Three tests and the mechanism show it is not; a Berry phase cannot decide:

1. **SSH analogy fails.** The (0,1) coherence block is tight-binding with
   uniform hopping (no alternating structure). No SSH topology.

2. **Chiral blocks have full rank (N=5).** The off-diagonal blocks of
   L_c in the Pi-eigenbasis have no zero singular values: no zero
   modes. (A finite square block always has index 0, so this rank test
   reads zero modes, not a winding.) (N=4 shows 4+3 zero SVs,
   the even-N centre modes of Phase 2, absent at N = 3, 5.)

3. **A Berry phase cannot decide here.** The parameter path (uniform
   to sacrifice noise) is an open segment of an interval, on which a
   Berry phase is defined only by a closing rule and carries no
   topological meaning. Summed step by step, the phase grows with the
   number of steps; it measures the eigensolver, not the system.

4. **No sharp phase boundary.** The protection factor is continuous as
   the sacrifice moves from edge to center (its largest step shrinks when
   the grid is refined) and the slowest rate never closes. No
   topological phase transition.

5. **The mechanism is exposed.** Under edge sacrifice: center-heavy
   modes survive. Under center sacrifice: modes dark at the center survive.
   The same Hamiltonian modes exist in both cases; the noise selects
   which ones dominate. The sacrifice zone exploits geometry, not
   topology.

---

## Phase 1: SSH Analogy

The SSH model (Su-Schrieffer-Heeger, 1979) requires alternating
coupling t1-t2. The Heisenberg chain has uniform coupling J. The
dispersion on the (0,1) coherence block is:

    omega_k = 4J(1 - cos(pi*k/N)),    k = 1, ..., N-1

This is tight-binding with hopping 2J, quantized on an open chain.
The eigenstates are the NEUMANN standing waves
psi_k(j) proportional to cos(pi*k*(j - 1/2)/N), the eigenbasis of the open
chain's Laplacian (D10 Step 5). The Dirichlet family sin(pi*k*j/(N+1)) is a
different object's, the XY single-excitation Hamiltonian's (F2b); it is not
the eigenbasis here, and the denominator is the tell.

For the highest-frequency mode (k = N-1): amplitude peaks at the
chain center and is smallest at the edges. For the lowest (k = 1):
amplitude is distributed toward the edges.

The non-uniform gamma profile (sacrifice zone) does NOT create
alternating hopping. It adds site-dependent DAMPING to the pre-existing
standing wave modes. The spatial profile of each mode is determined
by H alone; gamma only determines which modes survive longest.

**No SSH topology.** Uniform hopping gives winding number zero at any
gamma profile.

---

## Phase 2: Chiral Block Structure

The conjugation operator Pi has eigenvalues {+1, -1, +i, -i}
(fourth roots of unity, since Pi^4 = I). The anti-commutation
{Pi, L_c} = 0 forces L_c to be block-off-diagonal between paired
eigenspaces:

    L_c: V_{+1} <-> V_{-1}    (off-diagonal blocks A1, B1)
    L_c: V_{+i} <-> V_{-i}    (off-diagonal blocks A2, B2)

The on-diagonal blocks are zero to machine precision (verified: residual
< 1e-15), confirming the chiral decomposition.

### Singular value analysis

If the off-diagonal block A has zero singular values, the system has
topological zero modes (analogous to SSH edge modes). If A has full
rank, the topology is trivial.

| N | Pair | Block size | Min SV | Zero SVs | Verdict |
|---|------|-----------|--------|----------|---------|
| 4 | +/-1 | 64x64 | 3.8e-17 | 4 | Zero modes present |
| 4 | +/-i | 64x64 | 3.7e-16 | 3 | Zero modes present |
| 5 | +/-1 | 256x256 | 2.4e-03 | 0 | **Full rank (trivial)** |
| 5 | +/-i | 256x256 | 2.4e-03 | 0 | **Full rank (trivial)** |

**N=4 has 7 zero singular values (4+3).** A zero singular value is a mode
of L exactly at the spectral centre −Σγ, which is XY-weight N/2. At even N
that weight holds, among others, the orbit sums of letter multisets with N/2
letters X or Y, each of which commutes with the Heisenberg H and decays at
exactly the centre rate (at N = 4 they give 9 of the 14; the rest are
commutant elements of the corner rows and one weight-mixing mode); at odd N the weight N/2 is not an integer, so no pure-weight mode can sit
there, and none does at N = 3, 5. So the zero modes are an even-N law, not a finite-size artifact:
dim ker(L + Σγ) is 0, 14, 0, 16 at N = 3, 4, 5, 6 (the centre counts of
[Degeneracy Palindrome](DEGENERACY_PALINDROME.md)), and the table lists the
A blocks only, 7 of the 14 at N = 4. Testing N = 5 alone could not see them.

**N=5 is the system where the r = 0.994 localization was measured.**
Full rank, no zero modes. A winding number (an integer invariant that counts how many times the system's parameters "wrap around" a cycle; nonzero means topologically protected edge states) needs a periodic parameter, which a finite open chain does not have; what the rank shows is that there are no zero modes for such an invariant to protect.

---

## Phase 3: Berry Phase

Parametrize the gamma profile continuously from uniform (theta = 0) to
sacrifice-edge (theta = 1):

    gamma_k(theta) = (1-theta) * gamma_mean + theta * sacrifice_k

Follow the slowest oscillating mode. Coherences split into blocks by how
many excitations the ket and the bra carry (their popcounts), and H and
Z dephasing never move a coherence from one block to another. The mode
lives in the (1,0) block, one excitation in the ket and none in the bra,
which at N = 5 is 5×5. Inside it the mode is simple, with a gap of 2.0
along the whole path. In the full space it has an exact twin, its image
under flipping every qubit (X^N), which lives in the (4,5) block, so a
full-space eigensolver returns an arbitrary mix of the two. At θ = 0 all
five block modes share the rate 2γ_mean = 0.216, and only continuity
from θ = 1 says which one is followed.

A Berry phase is the geometric phase picked up when the parameters go
around a closed loop; quantized values (0 or π) can signal topological
protection, as the Zak phase does when the parameter runs around a
periodic Brillouin zone. Here θ runs along an interval, and the path is
an open segment. An open path has a gauge-invariant phase only once a
rule closes it, and its value then depends on that rule; a closed loop
on an interval only retraces itself. Because H is real, the block is
also complex symmetric (equal to its own transpose): the left
eigenvector is the transposed right one and the biorthogonal connection
vanishes, so a loop in several γ parameters would read 0 or π, π only
when it encircles an exceptional point. This path encircles nothing.
**The test cannot decide topology on this path, in either direction.**

Reading it anyway shows the trap. Two readings are computed: the
per-step sum of the phases between neighbouring left and right
eigenvectors (both are needed, since the Liouvillian is not Hermitian),

    φ = −Im Σ_n log ⟨ψ_L(θ_{n+1}) | ψ_R(θ_n)⟩,

with each eigenvector in the eigensolver's own phase convention; and the
phase of the product with the pair normalized to ⟨ψ_L|ψ_R⟩ = 1 at every
step and closed back to the start (the Pancharatnam closing rule),

    Π_n ⟨ψ_L(θ_n) | ψ_R(θ_{n+1})⟩ · ⟨ψ_L(1) | ψ_R(0)⟩,

in which every vector's arbitrary phase cancels:

| steps | per-step sum | closed normalized phase |
|-------|--------------|-------------------------|
| 60 | −0.772 | −4.06·10⁻⁴ |
| 600 | −7.83 | −4.00·10⁻⁴ |
| 6000 | −78.5 | −4.00·10⁻⁴ |

The per-step sum grows in proportion to the number of steps. The
eigensolver normalizes left and right vectors separately, so each step
adds the phase between the two at that point of the path (from 0 at
θ = 0 to 0.025 rad at θ = 1, 0.013 on average), a property of the
solver's normalization and not of the path. The closed product converges, to
−4.0·10⁻⁴ rad: gauge-invariant under the closing rule, and neither
quantized nor a topological invariant. The verdict rests on the other
tests.

---

## Phase 5: Robustness Sweep

Sweep the gamma profile from edge sacrifice (theta = +1) through uniform
(theta = 0) to center sacrifice (theta = -1):

| theta | Profile type | Min rate | Protection factor |
|-------|-------------|----------|-------------------|
| +1.00 | Edge sacrifice | 0.056 | 19.3 |
| +0.50 | Edge-biased | 0.137 | 7.9 |
| 0.00 | Uniform | 0.216 | 5.0 |
| -0.50 | Center-biased | 0.118 | 9.2 |
| -1.00 | Center sacrifice | 0.020 | 54.0 |

The protection factor is the largest rate over the slowest, and it is
**continuous**. Its steepest stretch is the tail toward center
sacrifice, where the slowest rate falls toward 2γ_low = 0.02; there the
step between neighbouring grid points is largest (17.8 near θ = −0.98 on
41 points), but doubling the grid to 81 points shrinks it by a factor
1.67. A discontinuity would keep its size; a continuous curve halves it
once the grid resolves the tail. The slowest rate never reaches zero. No sharp jump, no topological phase boundary.

### The mechanism exposed

The localization profile of the SLOWEST modes changes with the noise:

| Noise profile | Slowest mode profile | Rate |
|--------------|---------------------|------|
| Edge sacrifice | [0.86, 1.05, **1.17**, 1.05, 0.87] (center-heavy) | 0.056 |
| Uniform | [1.00, 1.00, 1.00, 1.00, 1.00] (no selection: all 16 tie) | 0.216 |
| Center sacrifice | [**1.04**, **1.04**, 0.83, **1.04**, **1.04**] (center-dark) | 0.020 |

Each row averages every oscillating mode at the slowest rate, each
eigenspace read in a basis-free way (a single vector of a degenerate
eigenvalue would depend on the eigensolver's basis). Under uniform noise
nothing is selected: all sixteen oscillating modes decay at 2γ, and they
carry the four standing-wave profiles of the frequencies ω = 7.24, 5.24,
2.76 and 0.76, which average to flat. The ω = 7.24 mode among them
carries [0.87, 1.05, 1.17, 1.05, 0.87], the profile edge sacrifice
selects. Center sacrifice selects the modes at the two frequencies ω = 5.24 and 0.76,
both dark at the center: the center qubit carries their least weight
(0.83, against 0.95 to 1.13 elsewhere).

Under edge sacrifice: center-heavy modes survive (they avoid the noisy
edge). Under center sacrifice: center-dark modes survive (they avoid the
noisy center). The SAME set of Hamiltonian modes exists in all three
cases. The noise profile selects which subset dominates.

This is the reading [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)
rests on. Each mode keeps nearly the same qubit-weight pattern under
sacrifice noise and uniform noise (the ω = 7.24 mode's weights move by at
most 0.003 between the two), because the pattern comes from the chain's
geometry; which mode is SLOWEST depends on the noise profile.

---

## Phase 4: Mode Counting

| N | Center-localized / total oscillating | Fraction |
|---|-------------------------------------|----------|
| 3 | 20 / 40 | 50.0% |
| 4 | 178 / 202 | 88.1% |
| 5 | 356 / 904 | 39.4% |

A mode is center-localized when a middle site carries the highest
weight: the single center qubit at odd N, the middle pair at even N
(the reflection makes the two middle weights exactly equal in every
mode, so at N = 4 no single site could be the center). Equality is read
to rounding; the nearest genuine near-tie is 8.6·10⁻⁵ (at N = 5), and a
coarser threshold counts more modes there (372 within 10⁻³, 404 within
10⁻²). The label belongs to a single vector, and at N = 3 and 5 every
oscillating eigenvalue of the full spectrum has an exact twin (its X^N
image in the mirrored popcount block), so each mode is read in its own
popcount block, the sector-pure basis. There every eigenvalue is simple
at N = 3 and 5; at N = 4 the in-block pairs kept their label under 20
random changes of basis each.

The fraction is NOT determined by a topological invariant (which would
be an integer, stable across N). It moves with what counts as the
middle and with where the threshold is drawn, which no invariant could.

---

## Connection to Previous Results

### RMT (finite-N pooled spacings)

The Poisson-like statistics do not prove integrability or topology.
Invariant sectors can suppress repulsion in pooled spacings; the topology
verdict here rests on its own block/rank and edge-mode diagnostics, not on RMT.

### PT-Symmetry Analysis (chiral, order 4)

Pi provides the P-type relation used in the block decomposition (Phase 2).
Its global order-four structure is resolved through Pi² = U_X; within a
fixed parity sector, a phase-normalized Pi is involutive. Both tested
subsystems have full rank at N=5 (trivial), but this does not assign a
global or fully reduced SRP class.

### Sacrifice Zone Formula

The sacrifice zone formula (gamma_edge = N*gamma_base - (N-1)*epsilon)
is NOT a topological phase boundary condition. It is the noise profile
that maximizes the GEOMETRIC advantage: concentrate noise where
edge-heavy modes live, protect where center-heavy modes live. The
360x improvement at N=5 comes from the standing wave structure of the
chain, not from a topological invariant.

---

## What This Establishes

1. **The localization is geometric.** Standing wave patterns on a 1D
   chain create modes with varying spatial profiles (center-heavy to
   edge-heavy). This is a property of the Hamiltonian, independent of
   noise.

2. **The sacrifice zone exploits geometry.** By placing noise at the
   chain edge, center-heavy modes (which naturally avoid the edge)
   survive. The optimization is geometric matching of noise to mode
   structure.

3. **No topological invariant.** No SSH structure, no zero modes in the
   chiral blocks at N = 5, no sharp phase boundary. The three tests that
   can decide agree; the Berry phase cannot decide on an open
   one-parameter path.

4. **The N=4 zero SVs are the even-N centre modes.** They sit exactly at
   −Σγ, which no pure-weight mode can reach at odd N, and they return at N = 6 (16 modes).
   They come from operators commuting with H at XY-weight N/2, not
   from a topological invariant.

5. **The profile reverses with noise.** Edge sacrifice protects center
   modes. Center sacrifice protects center-dark modes. The Hamiltonian modes
   are fixed; the noise selects the survivors. Confirmed independently
   in proton water chains (5.1x at N=5) and DNA G-C base pairs (3.8x)
   ([Proton Water Chain](PROTON_WATER_CHAIN.md),
   [DNA Base Pairing](DNA_BASE_PAIRING.md)).

## What This Does NOT Rule Out

- Topological effects at larger N (only tested N=3, 4, 5)
- Topological effects in ring or star topology (only tested chain)
- Topological effects under non-Z dephasing (only tested Z-dephasing)
- A more exotic topological invariant beyond the tested SSH-type winding number

---

## Reproducibility

| Component | Location |
|-----------|----------|
| Script | [`simulations/topological_edge_modes.py`](../simulations/topological_edge_modes.py) |
| Output | [`simulations/results/topological_edge_modes.txt`](../simulations/results/topological_edge_modes.txt) |

---

## References

- Su, W.P., Schrieffer, J.R., Heeger, A.J. (1979). "Solitons in
  Polyacetylene." PRL 42, 1698.
- Lieu, S., McGinley, M., Cooper, N.R. (2020). "Tenfold Way for
  Quadratic Lindbladians." PRL 124, 040401.
- Kawabata, K. et al. (2019). "Symmetry and Topology in Non-Hermitian
  Physics." PRX 9, 041015.
- Bianco, R., Resta, R. (2011). "Mapping topological order in
  coordinate space." PRB 84, 241106(R).
- Song, F., Yao, S., Wang, Z. (2019). "Non-Hermitian Skin Effect and
  Chiral Damping." PRL 123, 170401.

---

*The data screamed "topological." The math says "geometric." Three tests
and one mechanism, one answer: the sacrifice zone exploits the standing wave structure of
the chain, not a topological invariant. The mode profiles are determined
by the Hamiltonian. The noise selects the survivors. No topology needed.*
