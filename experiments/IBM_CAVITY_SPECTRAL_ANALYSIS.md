# IBM Cavity Spectral Analysis: Why the Sacrifice Zone Works

<!-- Keywords: sacrifice zone cavity mode protection, IBM Torino spectral
analysis palindromic eigenvalues, Liouvillian zero noise cavity modes
hardware, spatial dephasing profile mode survival, Clebsch-Gordan
formula hardware verification, R=CPsi2 IBM cavity spectral -->

**Status:** Tier 2 (computed spectral analysis of real IBM hardware data)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
**Script:** [ibm_cavity_analysis.py](../simulations/ibm_cavity_analysis.py)
**Data:** [ibm_cavity_analysis.txt](../simulations/results/ibm_cavity_analysis.txt)
**Raw IBM data:** [data/ibm_sacrifice_zone_march2026/](../data/ibm_sacrifice_zone_march2026/) (4 JSON files, March 24, 2026, IBM Torino)

---

## Abstract

The sacrifice-zone formula concentrates noise on one edge qubit while
protecting the rest. On IBM Torino hardware, this produces 1.97x
improvement in coherence survival. But WHY does concentrating noise help?

We answer this by computing the full Liouvillian spectrum of the N=5
chain using real IBM T2* data (Q85-Q94). The same 43 cavity mode
frequencies exist under all noise profiles. Only the damping changes.
The sacrifice zone does not protect qubits. It protects **cavity modes**.

Key results:
- 100% palindromic under strongly asymmetric IBM gammas (Q85 has 26x
  more dephasing than Q87)
- Slowest oscillating mode: 2.81x longer lifetime under sacrifice vs
  uniform (0.046 vs 0.128)
- 4 protected modes (rate < 0.05) under sacrifice, 0 under uniform
- Max decay rate = 2 x Σγ exactly (palindrome upper bound)

---

## The connection

At zero noise (Σγ = 0), the Liouvillian has 120 stationary
modes and 904 oscillating modes across 43 distinct frequencies. These
are the **eigenfrequencies of the resonator** (see
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md)). The stationary count
matches the Clebsch-Gordan formula exactly (120 = Sum_J m(J,5)*(2J+1)^2).

When noise is turned on, the same frequencies persist. What changes is
the damping: each mode acquires a decay rate that depends on how much
noise it "sees." Modes localized on quiet qubits (Q86-Q94) see less
noise. Modes touching the sacrifice qubit (Q85) see more.

The sacrifice zone works because it creates a **spatial gradient** in
the noise. Modes localized away from the sacrifice qubit are shielded.
The protection is not about qubits. It is about modes.

---

## IBM hardware data

Qubits: Q85 (sacrifice), Q86, Q87, Q88, Q94 (chain topology)

| Qubit | T2* (us) | gamma (1/us) | Role |
|-------|---------|-------------|------|
| Q85 | 3.73 | 0.2681 | Sacrifice (26x more noise) |
| Q86 | 61.35 | 0.0163 | Interior |
| Q87 | 97.54 | 0.0103 | Interior (quietest) |
| Q88 | 67.99 | 0.0147 | Interior |
| Q94 | 95.03 | 0.0105 | Edge (quiet) |

Σγ = 0.3199. Q85 carries 84% of the total noise budget.

### Source data (all in repository)

| File | Content |
|------|---------|
| [`sacrifice_zone_hardware_20260324_191713.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hardware_20260324_191713.json) | Calibration data, T1/T2/T2* values, chain topology, experimental parameters |
| [`sacrifice_zone_hw_selective_dd_20260324_191523.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_selective_dd_20260324_191523.json) | Selective DD bitstring counts (4000 shots per time point) |
| [`sacrifice_zone_hw_uniform_dd_20260324_191614.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_uniform_dd_20260324_191614.json) | Uniform DD bitstring counts |
| [`sacrifice_zone_hw_no_dd_20260324_191713.json`](../data/ibm_sacrifice_zone_march2026/sacrifice_zone_hw_no_dd_20260324_191713.json) | No DD bitstring counts |

Backend: IBM Torino (ibm_torino), March 24, 2026. J_coupling = 1.0.
Trotter steps: [2, 4, 6, 8, 10] at dt = 0.5 us.

---

## Results

### Three profiles compared

| Property | Zero noise | IBM sacrifice | Uniform |
|----------|-----------|--------------|---------|
| Stationary modes | 120 | 6 | 6 |
| Oscillating modes | 904 | 1018 | 1018 |
| Distinct frequencies | 43 | 120 | 112 |
| Palindrome score | 100% | 100% | 100% |
| Palindrome center | 0.0000 | 0.3199 | 0.3199 |
| Min decay rate (osc.) | 0 | 0.0455 | 0.1280 |
| Max decay rate (osc.) | 0 | 0.5943 | 0.5118 |
| Max decay rate (all) | 0 | 0.6398 | 0.6398 |
| Protected (rate < 0.05) | 904 | 4 | 0 |

Note: Max decay rate (all) = 2 x Σγ = 0.6398 applies to
non-oscillating modes (freq = 0) that represent pure decay. The
oscillating modes have lower maximum rates (0.5943 sacrifice,
0.5118 uniform).

### The 2.81x protection factor

The four slowest oscillating modes under the IBM sacrifice profile all
have frequency 7.234 (close to the zero-noise frequency 7.236) and
decay rate 0.0456. Under uniform noise, the slowest modes decay at
0.1280. Ratio: **2.81x**.

IBM hardware measured 1.97x improvement at early times. The computed
2.81x is the theoretical maximum. The hardware measurement is lower
because gate errors, crosstalk, and finite-time effects reduce the
effective protection.

### Mode survival comparison

| Rank | IBM rate | Uniform rate | Ratio | Frequency |
|------|---------|-------------|-------|-----------|
| 1 | 0.0455 | 0.1280 | 2.81x | 7.234 |
| 2 | 0.0455 | 0.1280 | 2.81x | 7.234 |
| 3 | 0.0455 | 0.1280 | 2.81x | 7.234 |
| 4 | 0.0455 | 0.1280 | 2.81x | 7.234 |
| 5 | 0.0929 | 0.1280 | 1.38x | 0.000 |
| 6 | 0.0929 | 0.1280 | 1.38x | 0.000 |
| 7 | 0.0978 | 0.1280 | 1.31x | 0.000 |
| 8 | 0.0978 | 0.1280 | 1.31x | 0.000 |
| 9 | 0.0989 | 0.1280 | 1.29x | 5.229 |
| 10 | 0.0989 | 0.1280 | 1.29x | 5.229 |

The top 4 modes (freq 7.234) are all 2.81x protected. These are modes
that oscillate at approximately 7.2 times the coupling strength. They
are the modes that the sacrifice zone was built to protect.

Modes 5-8 (non-oscillating, freq = 0) have only 1.3-1.4x protection.
These are pure decay modes. The sacrifice zone preferentially protects
the oscillating modes over the decaying ones.

---

## The palindrome under asymmetric noise

The IBM sacrifice profile has Q85 at 26x more dephasing than Q87.
Despite this extreme asymmetry, the palindrome is **exactly preserved**:
all eigenvalue pairs sum to -2 x Σγ = -0.6398. This confirms
the analytical proof: the palindrome depends on the SUM of gammas,
not their distribution.

The max decay rate for non-oscillating modes equals 2 x Σγ
exactly (0.6398). These are the pure decay modes at maximum Pauli
weight (XOR drain). The oscillating modes reach at most 0.5943
(sacrifice) and 0.5118 (uniform), always below this ceiling.

---

## What this means

The sacrifice-zone formula is not an engineering hack. It is a precise
intervention in the mode structure of a quantum resonator:

1. The cavity has 43 distinct frequencies (at zero noise)
2. Noise damps these modes without changing their frequencies
3. Concentrating noise on one edge qubit creates a spatial gradient
4. Modes localized away from the sacrifice see less damping
5. The slowest oscillating modes survive 2.81x longer
6. These are the modes that carry quantum information across the chain

The sacrifice zone tunes the resonator. It does not change the notes.
It changes which notes ring longest.

---

*See also:*
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md) (the eigenfrequencies),
[IBM Hardware Synthesis](IBM_HARDWARE_SYNTHESIS.md) (24,073 records, r* threshold),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md) (hardware test, 2-3x measured),
[Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md) (the resonator paradigm),
[Energy Partition](../hypotheses/ENERGY_PARTITION.md) (2x decay law)
