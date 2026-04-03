# Analytical w=1 Spectrum of the Heisenberg Chain

**Status:** Verified to machine precision (N=2-6, zero error)
**Date:** March 31, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [analytical_spectrum_verify.py](../simulations/analytical_spectrum_verify.py)
**Results:** [analytical_spectrum_verify.txt](../simulations/results/analytical_spectrum_verify.txt)
**Depends on:** [Thermal Breaking](THERMAL_BREAKING.md),
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)

---

## What this document is about

The w=1 sector of the Liouvillian (modes where exactly one qubit carries
an X or Y Pauli operator) has an exact closed-form frequency formula:
ω_k = 4J·(1 − cos(πk/N)). This replaces full matrix diagonalization
(which scales as 4^{3N}) with a single cosine evaluation per mode.
Verified to machine precision for N=2 through 6 (15 frequencies, zero error).

---

## The Result

The w=1 sector of the Liouvillian (modes where exactly one qubit
carries X or Y, all others I or Z) has an exact dispersion relation (the formula relating oscillation
frequency to mode index) for the Heisenberg chain:

    ω_k = 4J · (1 - cos(πk/N)),    k = 1, ..., N-1

This gives N-1 distinct frequencies for an N-qubit chain. Under uniform
Z-dephasing, all w=1 modes decay at the same rate 2γ, so the Q-factor
of each mode is:

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

15/15 frequencies match at machine precision. Five other candidate
formulas (including 4J·sin(πk/(2N)), 4J(1-cos(πk/(N+1))), and
2J(1-cos(πk/N))) all fail.

---

## What This Replaces

Previously, computing w=1 frequencies required:
1. Build Hamiltonian (2^N × 2^N matrix)
2. Build Liouvillian (4^N × 4^N matrix)
3. Diagonalize (O(4^{3N}) operations)
4. Filter for w=1 modes

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

For large N: bandwidth → 8J. The w=1 band saturates at 8J, covering
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
| Number of w=1 modes | 4 | 9 | 99 | N-1 |
| ω_max | 7.24 | 7.90 | 8.00 | 8J |
| ω_min | 0.76 | 0.10 | 0.001 | 0 |
| V(N) = Q_max ratio | 1.81 | 1.95 | 2.00 | 2 |
| Q spread | 9.47 | 76.4 | ~N² | ∞ |

---

## What This Does NOT Cover

The w=1 sector contains N-1 frequencies. The full Liouvillian has
up to 4^N eigenvalues distributed across sectors w=0, w=1, ..., w=N.
This formula covers only w=1.

- **w=0 sector:** All {I,Z} operators. Stationary (rate=0, freq=0).
  Not oscillating.
- **w=2, w=3, ... sectors:** Higher decay rates (4γ, 6γ, ...).
  Their frequencies are determined by multi-magnon excitations, which
  are NOT simple cosine bands for the Heisenberg model (magnon-magnon
  interactions create complex spectra).
- **Mode localization:** The spatial profile of w=1 modes (center vs
  edge weight) is NOT given by this formula. It requires the
  eigenvectors, not just the eigenvalues. The
  [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) result
  (r=0.994) was computed from eigenvectors.
- **Non-uniform dephasing:** Under sacrifice-zone profiles, w=1 modes
  acquire different decay rates. The frequency formula still holds
  (frequencies are Hamiltonian-determined), but the Q-factors change.

---

## Connection to Known Physics

The formula ω_k = 4J(1-cos(πk/N)) is the dispersion relation of
the nearest-neighbor tight-binding model (the simplest lattice model
where a particle hops between adjacent sites, with hopping amplitude 2J) with
the quantization k_n = πn/N. This quantization corresponds to a
chain with specific boundary conditions.

For the Heisenberg XXX model, the single-magnon sector is equivalent
to a tight-binding hopping problem. The w=1 Liouvillian sector
inherits this structure because [H, σ_+^(j)] creates single-magnon
transitions. The factor 4J (instead of the standard 2J for a
tight-binding chain) comes from the Liouvillian being a commutator
(double action of H).

A formal proof that ω_k = 4J(1-cos(πk/N)) holds for ALL N remains
open (see [Thermal Breaking, Open Question 7](THERMAL_BREAKING.md)).
The numerical evidence (zero error for N=2-6, 15 frequencies) is
strong but not a proof. A third independent validation comes from the
spectral form factor: the SFF modulation peak matches ω_min =
4J(1-cos(π/N)) to <1% for N=2-4 and N=6, confirming the dispersion
relation in the time domain
([Spectral Form Factor](SPECTRAL_FORM_FACTOR.md)).

---

## Tier Assessment

- Dispersion relation ω_k = 4J(1-cos(πk/N)): **Tier 1-2** (verified
  to machine precision N=2-6, motivated by tight-binding analogy,
  formal proof for all N: open)
- Q-factor spectrum formulas: **Tier 1** (algebraic consequences of
  the dispersion relation + the proven 2γ decay rate for w=1)
- Mode density and large-N scaling: **Tier 2** (standard band theory
  applied to the verified dispersion relation)
