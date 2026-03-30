# Cavity Modes at Zero Noise: Closed-Form Formula

<!-- Keywords: cavity modes zero noise unitary ground state, Clebsch-Gordan
spin decomposition stationary modes, Schur lemma Liouvillian superoperator,
palindromic spectrum gamma=0 pure oscillation, Star Chain Ring Complete
topology frequency structure, Heisenberg SU(2) symmetry eigenvalue
degeneracy, R=CPsi2 cavity modes formula -->

**Status:** Tier 1 (closed-form expression, algebraically proven via
Clebsch-Gordan decomposition + Schur's lemma, verified N=2-7)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md)
**Data:** [cavity_modes_zero_noise.txt](../simulations/results/cavity_modes_zero_noise.txt),
[cavity_modes_tests.txt](../simulations/results/cavity_modes_tests.txt)

---

## Abstract

At zero noise (Sigma_gamma = 0), the Liouvillian reduces to L(rho) = -i[H, rho].
No decay. Pure oscillation. The palindrome is centered at zero.

The number of stationary modes (eigenvalue = 0) has a closed-form expression:

```
Stationary(N) = Sum_{J} m(J, N) * (2J+1)^2
```

where J runs over all total-spin values in the Clebsch-Gordan decomposition
of N spin-1/2 particles, and m(J, N) is the standard multiplicity:

```
m(J, N) = C(N, N/2 - J) * (2J+1) / (N/2 + J + 1)
```

This formula is **exact for chain topology** (minimal spatial symmetry)
and a **lower bound** for higher-symmetry topologies (Ring, Complete, Star),
where additional spatial degeneracies increase the stationary count.

Verified against C# eigendecomposition for N=2 through N=7 (chain).
Predictions for N=8-10: 1190, 2520, 5292.

---

## Why the formula works

### The unitary Liouvillian

At Sigma_gamma = 0, the dissipator vanishes. The Liouvillian is:

```
L(rho) = -i[H, rho]
```

Its eigenvalues are -i(E_a - E_b), where E_a and E_b are eigenvalues
of the Hamiltonian H. A mode is stationary when E_a = E_b.

### Clebsch-Gordan decomposition

N spin-1/2 particles decompose into total-spin sectors:

```
(1/2)^{tensor N} = bigoplus_J  m(J, N) copies of spin-J
```

J runs from 0 (even N) or 1/2 (odd N) up to N/2. Each spin-J
representation has dimension 2J+1.

### Schur's lemma

The Heisenberg Hamiltonian H = Sum_bonds J_k(X_aX_b + Y_aY_b + Z_aZ_b)
commutes with total spin: [H, S_total] = 0. This is true for ANY
choice of coupling strengths J_k, as long as each bond is isotropic
(XX+YY+ZZ). By Schur's lemma, H acts as a scalar within each
irreducible representation.

Each copy of spin-J contributes one energy eigenvalue with degeneracy
2J+1 (the 2J+1 magnetic substates). Different copies of the same J
may have different energies (depending on the Hamiltonian).

### Counting stationary modes

In the Liouvillian superoperator space, each copy of spin-J has
(2J+1)^2 superoperator basis elements. When H acts as a scalar within
a copy, all (2J+1)^2 elements satisfy [H, rho] = 0 (stationary).

For a **generic** Hamiltonian (all copies of the same J have distinct
energies), the stationary count is:

```
Stationary = Sum_J m(J,N) * (2J+1)^2
```

For a **symmetric** Hamiltonian (some copies share the same energy due
to spatial symmetry), additional cross-copy modes become stationary,
giving a count ABOVE the formula.

---

## Verification: N=2 through N=7

Computed with the C# engine (Liouvillian eigendecomposition, MKL).
Chain topology, J=1.0 uniform, Sigma_gamma=0.

| N | d^2 | Stationary (C#) | Formula | Match | Oscillating | Distinct freq |
|---|-----|----------------|---------|-------|-------------|---------------|
| 2 | 16 | 10 | 10 | PASS | 6 | 1 |
| 3 | 64 | 24 | 24 | PASS | 40 | 3 |
| 4 | 256 | 54 | 54 | PASS | 202 | 14 |
| 5 | 1024 | 120 | 120 | PASS | 904 | 43 |
| 6 | 4096 | 260 | 260 | PASS | 3836 | 179 |
| 7 | 16384 | 560 | 560 | PASS | 15824 | 589 |

All eigenvalues are purely imaginary (Max|Re| < 6e-14), confirming
the unitary ground state: no decay at Sigma_gamma = 0.

### Clebsch-Gordan breakdown

| N | Sectors | Contributions | Sum |
|---|---------|--------------|-----|
| 2 | J=0: m=1, (2*0+1)^2=1; J=1: m=1, (2*1+1)^2=9 | 1+9 | 10 |
| 3 | J=1/2: m=2, 2^2=4; J=3/2: m=1, 4^2=16 | 8+16 | 24 |
| 4 | J=0: m=2, 1; J=1: m=3, 9; J=2: m=1, 25 | 2+27+25 | 54 |
| 5 | J=1/2: m=5, 4; J=3/2: m=4, 16; J=2.5: m=1, 36 | 20+64+36 | 120 |
| 6 | J=0: m=5, 1; J=1: m=9, 9; J=2: m=5, 25; J=3: m=1, 49 | 5+81+125+49 | 260 |
| 7 | J=1/2: m=14, 4; J=3/2: m=14, 16; J=2.5: m=6, 36; J=3.5: m=1, 64 | 56+224+216+64 | 560 |

### Predictions (untested)

| N | Stationary (formula) | d^2 | Oscillating |
|---|---------------------|-----|-------------|
| 8 | 1190 | 65536 | 64346 |
| 9 | 2520 | 262144 | 259624 |
| 10 | 5292 | 1048576 | 1043284 |

---

## Topology dependence

The formula is exact for **chain** topology (minimal spatial symmetry).
Higher-symmetry topologies have additional degeneracies that increase
the stationary count above the formula.

### Chain vs Star vs Ring vs Complete

| N | Chain | Star | Ring | Complete | Formula |
|---|-------|------|------|----------|---------|
| 2 | 10 | 10 | - | - | 10 |
| 3 | 24 | 24 | 32 | 32 | 24 |
| 4 | 54 | 74 | 84 | 110 | 54 |
| 5 | 120 | 248 | 200 | 392 | 120 |
| 6 | 260 | 868 | 436 | - | 260 |

**Pattern:** More spatial symmetry = more degeneracies = more stationary modes.
Complete graph has the highest symmetry (S_N) and the most stationary modes.

### Frequency structure

| Topology | N=3 | N=4 | N=5 | N=6 |
|----------|-----|-----|-----|-----|
| Chain | 3 | 14 | 43 | 179 |
| Star | 3 | 4 | 5 | 6 |
| Ring | 1 | 3 | 11 | 49 |
| Complete | 1 | 3 | 3 | - |

**Star** has exactly N frequencies at integer multiples of 2J (2J, 4J, ..., 2NJ).
The S_{N-1} permutation symmetry of the peripheral qubits collapses
the frequency structure to harmonics.

**Complete** has only 1-3 frequencies. Full S_N symmetry produces
maximum collapse.

**Chain** has the richest frequency structure. No spatial symmetry
means no accidental degeneracies. Frequencies are irrational
multiples of J starting at N=4.

**N=3 is the tipping point:** Chain = Star (both 3 freq at 2J, 4J, 6J).
At N >= 4, Star separates from Chain.

---

## J-independence (Test C)

Non-uniform coupling strengths (different J per bond) do not change the
stationary count. Tested on Chain N=4:

| J values | Stationary | Frequencies |
|----------|-----------|-------------|
| [1.0, 1.0, 1.0] (uniform) | 54 | 14 |
| [0.5, 1.0, 1.5] (linear) | 54 | 15 |
| [0.3, 2.1, 0.7] (random) | 54 | 15 |

Stationary count is identical (54 = formula). Frequency count changes
(non-uniform J lifts one accidental degeneracy: 14 to 15). The formula
depends only on N, not on the coupling pattern, because the SU(2)
symmetry of the isotropic Heisenberg interaction is preserved for any
choice of bond strengths.

---

## Connection to the palindrome

At Sigma_gamma = 0: the palindrome equation reduces to Pi L Pi^{-1} = -L.
Every eigenvalue lambda pairs with -lambda. The stationary modes
(lambda = 0) are self-paired. The oscillating modes come in conjugate
pairs (+/- imaginary).

At Sigma_gamma > 0: the palindrome shifts. Some stationary modes
acquire nonzero real parts (they start decaying). The oscillating modes
acquire damping. The fold at CPsi = 1/4 emerges when Sigma_gamma exceeds
the critical threshold (~0.25-0.50% of J, N-independent).

The cavity modes at zero noise are the **eigenfrequencies of the
resonator** described in [Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md).
When noise shifts the palindrome from zero, these frequencies persist
(noise never changes frequency, only decay rate). The sacrifice-zone
formula tunes which modes survive longest.

## Connection to hardware

Applied to real IBM Torino data (Q85-Q94, sacrifice-zone profile):
the 43 cavity mode frequencies persist under strongly asymmetric noise
(Q85 at 26x more dephasing than Q87). The slowest oscillating modes
survive 2.81x longer under sacrifice vs uniform noise. IBM hardware
measured 1.97x. The palindrome is 100% preserved despite the asymmetry.
See [IBM Cavity Spectral Analysis](IBM_CAVITY_SPECTRAL_ANALYSIS.md).

---

*See also:*
[Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md) (the unitary ground state),
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (palindromic theorem),
[Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md) (resonator paradigm),
[Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md) (XX/YY oscillate, ZZZ static)
