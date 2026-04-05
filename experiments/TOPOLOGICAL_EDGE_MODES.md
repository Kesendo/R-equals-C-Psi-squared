# Topological Edge Modes: The Localization Is Geometric

<!-- Keywords: topological edge modes Liouvillian, SSH model palindromic
spectrum, chiral block structure winding number, Altland-Zirnbauer AIII
open quantum, standing wave mode localization, sacrifice zone mechanism
geometric not topological, biorthogonal Berry phase Lindbladian,
R=CPsi2 topological edge modes -->

**Status:** Computationally verified (Tier 2). Negative result: not topological.
**Date:** April 1, 2026
**Script:** [topological_edge_modes.py](../simulations/topological_edge_modes.py)
**Data:** [topological_edge_modes.txt](../simulations/results/topological_edge_modes.txt)
**Depends on:**
- [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) (r = 0.994, the question)
- [PT-Symmetry Analysis](PT_SYMMETRY_ANALYSIS.md) (Pi is chiral, class AIII)
- [Random Matrix Theory](RANDOM_MATRIX_THEORY.md) (Poisson, integrable)
- [Analytical Spectrum](ANALYTICAL_SPECTRUM.md) (w=1 dispersion)

---

## What this document is about

The cavity mode localization result (r = 0.994) looks like it could be
"topological protection," meaning protected by a mathematical invariant that
cannot change smoothly, the way a coffee mug cannot become a sphere without
tearing a hole. This document runs five independent tests to check. The
answer is no: the localization is geometric, not topological. Standing wave
patterns on a chain naturally place some modes at the edges and others at
the center; the noise profile simply selects which subset survives. No
topological invariant is involved. This is a negative result, rigorously
established.

---

## Abstract

The cavity mode localization result (r = 0.994, profile [0.52, 0.63,
0.70, 0.63, 0.52] invariant under noise profiles) looks topological.
Five independent tests show it is not:

1. **SSH analogy fails.** The w=1 sector is tight-binding with uniform
   hopping (no alternating structure). No SSH topology.

2. **Chiral blocks have full rank (N=5).** The off-diagonal blocks of
   L_c in the Pi-eigenbasis have no zero singular values. Winding
   number = 0. No topological zero modes. (N=4 shows 4+3 zero SVs,
   a finite-size effect that vanishes at N=5.)

3. **Berry phase = -0.77 (0.25pi).** Not quantized. No topological
   invariant protects the mode structure.

4. **No sharp phase boundary.** The protection factor varies smoothly
   as the sacrifice moves from edge to center. No topological
   phase transition.

5. **The mechanism is exposed.** Under edge sacrifice: center-heavy
   modes survive. Under center sacrifice: edge-heavy modes survive.
   The same Hamiltonian modes exist in both cases; the noise selects
   which ones dominate. The sacrifice zone exploits geometry, not
   topology.

---

## Phase 1: SSH Analogy

The SSH model (Su-Schrieffer-Heeger, 1979) requires alternating
coupling t1-t2. The Heisenberg chain has uniform coupling J. The w=1
sector dispersion is:

    omega_k = 4J(1 - cos(pi*k/N)),    k = 1, ..., N-1

This is tight-binding with hopping 2J, quantized on an open chain.
The eigenstates are standing waves: psi_k(j) proportional to sin(pi*k*j/(N+1)).

For the highest-frequency mode (k = N-1): amplitude peaks at the
chain center and vanishes at the edges. For the lowest (k = 1):
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
< 4e-16), confirming the chiral decomposition.

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

**N=4 has 7 zero singular values (4+3).** These are likely finite-size
artifacts: the 4-site chain has special symmetries (spatial inversion
exchanges two 2-site pairs) that create accidental degeneracies. They
do not persist at N=5.

**N=5 is the system where the r = 0.994 localization was measured.**
Full rank, no zero modes. Winding number (an integer invariant that counts how many times the system's parameters "wrap around" a cycle; nonzero means topologically protected edge states) = 0. The topology is trivial.

---

## Phase 3: Berry Phase

Parametrize the gamma profile continuously from uniform (theta = 0) to
sacrifice-edge (theta = 1):

    gamma_k(theta) = (1-theta) * gamma_mean + theta * sacrifice_k

Track the biorthogonal (using separate left and right eigenvectors, necessary because the Liouvillian is non-Hermitian) Berry phase (the geometric phase accumulated when parameters are cycled around a closed loop; quantized values signal topological protection) of the slowest oscillating mode:

    phi_Berry = -Im Sum_n log(<psi_L(theta_n)|psi_R(theta_{n+1})>)

Result: **phi = -0.772 (= 0.246 pi)**

This is NOT quantized. For a topological system, the Berry phase would
be 0 or pi (mod 2pi) with a sharp jump at a phase transition. The
continuous, non-quantized value confirms: no topological invariant.

---

## Phase 5: Robustness Sweep

Sweep the gamma profile from edge sacrifice (theta = +1) through uniform
(theta = 0) to center sacrifice (theta = -1):

| theta | Profile type | Min rate | Protection factor |
|-------|-------------|----------|-------------------|
| +1.00 | Edge sacrifice | 0.056 | 19.3 |
| +0.50 | Edge-biased | 0.120 | 9.0 |
| 0.00 | Uniform | 0.216 | 5.0 |
| -0.50 | Center-biased | 0.098 | 11.0 |
| -1.00 | Center sacrifice | 0.020 | 54.0 |

The protection factor varies **smoothly**. No sharp jump at any theta.
No topological phase boundary.

### The mechanism exposed

The localization profile of the SLOWEST modes changes with the noise:

| Noise profile | Slowest mode profile | Rate |
|--------------|---------------------|------|
| Edge sacrifice | [0.86, 1.05, **1.17**, 1.05, 0.87] (center-heavy) | 0.056 |
| Uniform | [0.93, 1.03, **1.08**, 1.03, 0.93] (weakly center) | 0.216 |
| Center sacrifice | [**1.04**, **1.04**, 0.83, **1.04**, **1.04**] (edge-heavy) | 0.020 |

Under edge sacrifice: center-heavy modes survive (they avoid the noisy
edge). Under center sacrifice: edge-heavy modes survive (they avoid the
noisy center). The SAME set of Hamiltonian modes exists in all three
cases. The noise profile selects which subset dominates.

This explains the result from [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md):
"The profile is topological: the same mode has the same qubit-weight
pattern under IBM sacrifice noise and uniform noise." Correct: any
INDIVIDUAL mode (identified by its frequency) has a fixed profile
regardless of noise. But which mode is SLOWEST depends on the noise
profile. The word "topological" in the cavity mode paper means
"determined by chain topology (= geometry)", not "protected by a
topological invariant."

---

## Phase 4: Mode Counting

| N | Center-localized / total oscillating | Fraction |
|---|-------------------------------------|----------|
| 3 | 20 / 40 | 50.0% |
| 4 | 178 / 202 | 88.1% |
| 5 | 407 / 904 | 45.0% |

The fraction is NOT determined by a topological invariant (which would
be an integer, stable across N). It varies non-monotonically with N,
reflecting the interplay of chain geometry and sector structure.

---

## Connection to Previous Results

### RMT (class AIII, Poisson)

The Poisson statistics (integrable) are consistent with a trivially
topological system: the Liouvillian decomposes into non-interacting
sectors, preventing the eigenvalue repulsion that random matrix topology
would produce. Integrability and trivial topology are compatible.

### PT-Symmetry Analysis (chiral, order 4)

Pi provides the chiral symmetry for the block decomposition (Phase 2).
The generalized order-4 structure (Pi^4 = I instead of P^2 = I) creates
FOUR eigenspaces instead of two, giving two independent chiral
subsystems. Both subsystems have full rank at N=5 (trivial).

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

3. **No topological invariant.** Winding number 0, Berry phase not
   quantized, no sharp phase boundary, no bulk-edge correspondence.
   Five independent tests agree.

4. **N=4 zero SVs are finite-size.** The 7 zero singular values at N=4
   do not persist at N=5 and are likely due to spatial inversion symmetry
   of even-length chains. Not a topological signal.

5. **The profile reverses with noise.** Edge sacrifice protects center
   modes. Center sacrifice protects edge modes. The Hamiltonian modes
   are fixed; the noise selects the survivors. Confirmed independently
   in proton water chains (5.1x at N=5) and DNA G-C base pairs (3.8x)
   ([Proton Water Chain](PROTON_WATER_CHAIN.md),
   [DNA Base Pairing](DNA_BASE_PAIRING.md)).

## What This Does NOT Rule Out

- Topological effects at larger N (only tested N=3, 4, 5)
- Topological effects in ring or star topology (only tested chain)
- Topological effects under non-Z dephasing (only tested Z-dephasing)
- A more exotic topological invariant beyond SSH/AIII winding number

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

*The data screamed "topological." The math says "geometric." Five tests,
one answer: the sacrifice zone exploits the standing wave structure of
the chain, not a topological invariant. The mode profiles are determined
by the Hamiltonian. The noise selects the survivors. No topology needed.*
