# Analytical Spectrum of the Heisenberg Chain's (0,1) Coherence Block

**Status:** Verified N=2-6 against a full Liouvillian eigendecomposition, N=2-10 by block closure
**Date:** March 31, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [analytical_spectrum_verify.py](../simulations/analytical_spectrum_verify.py)
**Results:** [analytical_spectrum_verify.txt](../simulations/results/analytical_spectrum_verify.txt)
**Depends on:** [Thermal Breaking](THERMAL_BREAKING.md),
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)

---

## What this document is about

One block of the Liouvillian has a closed-form frequency spectrum, and this
document is that spectrum and what follows from it: the Q-factors, the
bandwidth, the mode density, the large-N scaling. It replaces full matrix
diagonalization, which scales as 4^{3N}, with a single cosine evaluation per
mode. Which block, and why it and not the larger sector it sits in, is the
next section.

---

## The Result

The (0,1) coherence block of the Liouvillian, spanned by the |0⟩⟨j| between
the ferromagnet and the single excitations, has an exact dispersion relation (the formula relating
oscillation frequency to mode index) for the Heisenberg chain:

    ω_k = 4J · (1 - cos(πk/N)),    k = 1, ..., N-1

This gives N-1 distinct frequencies for an N-qubit chain, N modes with one of
them at zero. The block sits INSIDE the Pauli w=1 sector (strings where exactly
one qubit carries X or Y, the rest I or Z) and is far smaller than it, 5 against
160 at N=5; that sector is not L-invariant and has no spectrum of its own, so it
is the block and not the sector that this formula is about (D10 Step 6). Under
uniform Z-dephasing every mode of the block decays at the same rate 2γ, so the
Q-factor of each mode is:

    Q_k = ω_k / (2γ) = 2J/γ · (1 - cos(πk/N))

---

## Verification

Tested against numerical Liouvillian eigenvalues. Six candidate
dispersion relations were compared. Only one matches:

| N | k | Predicted | Numerical | Error |
|:--|:--|:----------|:----------|:------|
| 2 | 1 | 4.000000 | 4.000000 | 0 |
| 3 | 1 | 2.000000 | 2.000000 | 0 |
| 3 | 2 | 6.000000 | 6.000000 | 0 |
| 4 | 1 | 1.171573 | 1.171573 | 0 |
| 4 | 2 | 4.000000 | 4.000000 | 0 |
| 4 | 3 | 6.828427 | 6.828427 | 0 |
| 5 | 1 | 0.763932 | 0.763932 | 0 |
| 5 | 2 | 2.763932 | 2.763932 | 0 |
| 5 | 3 | 5.236068 | 5.236068 | 0 |
| 5 | 4 | 7.236068 | 7.236068 | 0 |
| 6 | 1 | 0.535898 | 0.535898 | 0 |
| 6 | 2 | 2.000000 | 2.000000 | 0 |
| 6 | 3 | 4.000000 | 4.000000 | 0 |
| 6 | 4 | 6.000000 | 6.000000 | 0 |
| 6 | 5 | 7.464102 | 7.464102 | 0 |

15/15 frequencies match, at the 10⁻⁶ resolution the comparison rounds to.
The machine-precision figure is the block-closure route's, in
[D10](../docs/proofs/derivations/D10_W1_DISPERSION.md): entry-wise 0.00e+00
through N=10. Five other candidate
formulas (including 4J·sin(πk/(2N)), 4J(1-cos(πk/(N+1))), and
2J(1-cos(πk/N))) all fail.

---

## What This Replaces

Previously, computing these frequencies required:
1. Build Hamiltonian (2^N × 2^N matrix)
2. Build Liouvillian (4^N × 4^N matrix)
3. Diagonalize (O(4^{3N}) operations)
4. Filter by decay rate, |Re λ| = 2γ

For N=7: 16384² matrix, 92 minutes eigendecomposition.
For N=8: 65536² matrix, 10.6 hours.

Now: evaluate cos(πk/N) for k=1..N-1. Instant for any N.

---

## Derived Results

### A. The V-Effect gain for every mode

The [V-Effect gain](THERMAL_BREAKING.md) V(N) = 1+cos(π/N) was the
ratio of maximum frequencies. The full spectrum gives a gain for each
mode index k (comparing the k-th frequency at chain length N to that
at N=2):

    ω_k(N) / ω_1(2) = (1 - cos(πk/N)) / (1 - cos(π/2)) = 1 - cos(πk/N)

Since ω_1(2) = 4J (the only N=2 frequency), the gain for mode k at
chain length N relative to the single N=2 mode is simply 1-cos(πk/N).
The maximum (at k=N-1) gives V(N) = 1+cos(π/N) as before.

### B. Q-factor spectrum for any N

Under uniform Z-dephasing with rate γ per qubit:

    Q_k = 2J/γ · (1 - cos(πk/N))

| Property | Formula |
|:---------|:--------|
| Maximum Q | Q_{N-1} = 2J/γ · (1+cos(π/N)) |
| Minimum Q | Q_1 = 2J/γ · (1-cos(π/N)) → 2Jπ²/(γN²) for large N |
| Mean Q | Q_mean = 2J/γ (exactly, from Σcos = 0) |
| Q spread | Q_{max}/Q_{min} = (1+cos(π/N))/(1-cos(π/N)) = cot²(π/(2N)) |

The Q spread grows as ~N²/π² for large N. Longer chains have a wider
range of mode lifetimes.

### C. Frequency bandwidth

    Bandwidth = ω_{N-1} - ω_1 = 4J · (cos(π/N) - cos(π(N-1)/N)) = 8J · cos(π/N)

For large N: bandwidth → 8J. The block's band saturates at 8J, covering
frequencies from ~0 to ~8J.

### D. Mode density at large N

For N → ∞, the frequencies fill the interval [0, 8J] with density:

    ρ(ω) = N/(π · √(8Jω - ω²))

This is the density of states for a cosine band (van Hove singularities,
the divergences in mode density that occur at band edges where the
group velocity vanishes). Confirmed: the N=6 frequencies visually cluster near
the edges (0.54 and 7.46) relative to the interior (2, 4, 6).

### E. Scaling N → ∞

| Quantity | N=5 | N=10 | N=100 | N→∞ |
|:---------|:----|:-----|:------|:----|
| Number of oscillating modes | 4 | 9 | 99 | N-1 |
| ω_max = 4J(1−cos(π(N−1)/N)) | 7.236 | 7.804 | 7.998 | 8J |
| ω_min = 4J(1−cos(π/N)) | 0.764 | 0.196 | 0.00197 | 0 |
| V(N) = Q_max/Q_mean = 1+cos(π/N) | 1.809 | 1.951 | 2.000 | 2 |
| Q spread = cot²(π/2N) | 9.47 | 39.86 | 4052 | ∞ |

---

## What This Does NOT Cover

The (0,1) block carries N-1 nonzero frequencies. The full Liouvillian has
up to 4^N eigenvalues distributed across sectors w=0, w=1, ..., w=N.
This formula covers only that one block, which sits inside w=1.

- **w=0 sector:** All {I,Z} operators. Stationary (rate=0, freq=0).
  Not oscillating.
- **Coherences at Hamming distance 2, 3, ...:** higher decay rates (4γ, 6γ, ...).
  Their frequencies are determined by multi-magnon excitations, which
  are NOT simple cosine bands for the Heisenberg model (magnon-magnon
  interactions create complex spectra).
- **Mode localization:** The spatial profile of the block's modes (center vs
  edge weight) is NOT given by this formula. It requires the
  eigenvectors, not just the eigenvalues. The
  [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) result
  (r=0.994) was computed from eigenvectors.
- **Non-uniform dephasing:** under a sacrifice-zone profile the block's
  generator is −2iJ·𝓛 − 2·diag(γ), still exactly closed, but the frequency
  formula does NOT survive. diag(γ) and the Laplacian no longer commute, so
  the two cannot be diagonalized together, the eigenvalues are not
  −2γ_j − 2iJ·μ_m, and the operator is no longer normal. At N=5, J=1, on the IBM Torino
  sacrifice profile γ = [2.336, 0.099, 0.050, 0.072, 0.051] of
  [Concentrator Geometry](CONCENTRATOR_GEOMETRY.md), the |Im λ| read
  [0.238, 1.666, 2.068, 4.884, 7.144] against ω_k =
  [0, 0.764, 2.764, 5.236, 7.236]: the zero mode acquires a frequency and
  every other one moves. What the modes gain is separated decay rates,
  which is what makes a slowest one exist at all; see the
  `site_resolved_vacuum_block` open arc.

---

## Connection to Known Physics

The formula ω_k = 4J(1-cos(πk/N)) is the dispersion relation of
the nearest-neighbor tight-binding model (the simplest lattice model
where a particle hops between adjacent sites, with hopping amplitude 2J) with
the quantization k_n = πn/N. This quantization corresponds to a
chain with specific boundary conditions.

For the Heisenberg XXX model, the single-magnon sector is equivalent
to a tight-binding hopping problem. The (0,1) coherence block inherits
this structure because the ferromagnet is an eigenstate, so the block
sees the one-magnon Hamiltonian alone: a |0⟩⟨j| coherence oscillates at
E_magnon − E_vac, and ω_k IS that energy, with no doubling.

So the 4J is not a commutator effect, and it is worth saying because it
would be an easy one to assume. The Liouvillian acts on |0⟩⟨j| from
both sides, but one of those sides is the vacuum and contributes only
the constant E_vac; the frequency is a single energy difference, not
twice one. The 4 is two independent 2s from the Hamiltonian itself:
XX + YY each contribute J to the hop, so the hopping amplitude is 2J,
and the open-path Laplacian's eigenvalues are 2(1 − cos(πk/N)). Their
product is 4J·(1 − cos(πk/N)). The ZZ term is what supplies the degree
diagonal that makes it a Laplacian rather than an adjacency matrix.
The equivalent entry-wise statement, that L restricted to the block is
exactly −2iJ·𝓛 − 2γ·Id, is gated at N=2 through 10 by
[`simulations/d10_block_closure_verify.py`](../simulations/d10_block_closure_verify.py).

**Proven.** The dispersion relation is derived analytically in
[D10](../docs/proofs/derivations/D10_W1_DISPERSION.md) via reduction
of the (0,1) coherence block to a nearest-neighbour tight-binding
hopping problem on N sites. The proof shows that L restricted to that
block is −2iJ·𝓛 − 2γ·Id with 𝓛 the graph Laplacian, giving the cosine
eigenvalues for an open chain. It is NOT the w=1 Pauli sector that
reduces this way: D10 Step 2 gives an explicit counterexample showing
L_H does not close on that sector. Verified for N=2−6 against a full
Liouvillian eigendecomposition, and through N=10 by block closure.

---

## Tier Assessment

- Dispersion relation ω_k = 4J(1-cos(πk/N)): **Tier 1** (proven
  analytically in D10; verified N=2-6 against a full Liouvillian
  eigendecomposition, N=2-10 by block closure)
- Q-factor spectrum formulas: **Tier 1** (algebraic consequences of
  the dispersion relation + the proven 2γ decay rate on the block, at
  uniform γ)
- Mode density and large-N scaling: **Tier 2** (standard band theory
  applied to the verified dispersion relation)
