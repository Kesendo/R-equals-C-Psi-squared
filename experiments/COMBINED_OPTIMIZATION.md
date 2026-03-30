# Combined Optimization: Everything We Know, Applied

**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Hardware test:** April 9, 2026, IBM Torino (10 min QPU)

---

## What we combined

Three months of formulas, six hours of computation, one chain.

### The formulas

| Formula | Source | What it gives |
|---------|--------|--------------|
| Palindrome: Π L Π⁻¹ = -L - 2Σγ I | [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) | Every decay rate has a partner |
| Cavity modes: Stat(N) = Sum_J m(J,N)(2J+1)^2 | [Cavity Modes Formula](CAVITY_MODES_FORMULA.md) | 120 stationary + 43 frequencies at N=5 |
| Sacrifice zone: γ_edge = N γ_base - (N-1)ε | [Resonant Return](RESONANT_RETURN.md) | Concentrate noise on edge, protect rest |
| Mode profile: [0.52, 0.63, 0.70, 0.63, 0.52] | [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md) | Protected modes live in the center |
| Fold threshold: Σγ_crit/J ~ 0.5%, N-independent | [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md) | Minimum noise for irreversibility |

### The data

| Source | What it shows |
|--------|--------------|
| IBM March 24 (Q85-Q94) | Selective DD beats uniform DD by 2.0x, all 5 time points |
| IBM 24,073 calibration records | r* threshold at precision 0.000014, 12 permanent crossers |
| C# N=2-7 eigendecomposition | 100% palindromic, 54,118 eigenvalues at N=8 |

### The computation (March 30)

| Analysis | Key result |
|----------|-----------|
| IBM spectral (real T2* data) | 2.81x mode protection, 100% palindromic under 26x asymmetric noise |
| Eigenvector decomposition | r = 0.994 correlation: edge weight predicts decay rate |
| 330 chains mapped on Torino | Zero overlap between sacrifice ranking and T2 ranking |
| Eigenvalue efficiency | 3.72x with selective DD (best of 6 scenarios) |
| Time evolution (|01010>) | Quiet chains win absolute SumMI; sacrifice wins efficiency |

---

## The combined simulation

All formulas applied to Q85-Q94 with real IBM data:

### Eigenvalue analysis (6 configurations)

| Config | Σγ | Slowest mode | vs own uniform |
|--------|-----|-------------|---------------|
| Sacrifice, Selective DD | 0.290 | 0.031 | **3.72x** |
| Sacrifice, Uniform DD | 0.214 | 0.025 | 3.35x |
| Sacrifice, No DD | 0.320 | 0.046 | 2.81x |
| Mean-T2, Selective DD | 0.033 | 0.011 | 1.25x |
| Mean-T2, Uniform DD | 0.024 | 0.009 | 1.06x |
| Mean-T2, No DD | 0.048 | 0.018 | 1.06x |

### Time evolution with Neel state |01010>

| t (us) | Mean-T2 No DD | Mean-T2 + DD | Sacrifice No DD | Sacrifice + Sel DD |
|--------|:------------:|:-----------:|:--------------:|:-----------------:|
| 0.5 | 1.329 | 1.340 | 1.276 | 1.290 |
| 1.0 | 1.402 | 1.429 | 1.230 | 1.261 |
| 1.5 | 0.861 | 0.860 | 0.821 | 0.831 |
| 2.0 | 1.314 | 1.389 | 0.936 | 0.988 |
| 2.5 | 1.510 | 1.650 | 0.908 | 0.970 |
| 3.0 | 1.491 | 1.595 | 0.912 | 0.975 |
| 4.0 | 1.515 | 1.703 | 0.766 | 0.842 |

Quiet chains have higher absolute SumMI (less total noise = slower
decay). Sacrifice chains have higher efficiency per unit noise.

---

## What to test on April 9

**10 minutes QPU, ~500 circuits budget.**

### Setup

| | Chain A (sacrifice) | Chain B (quiet) |
|--|-------------------|----------------|
| Qubits | Q85-Q86-Q87-Q88-Q94 | Best mean-T2 chain from mapping |
| Mean T2 | ~65 us | ~200 us |
| Sacrifice qubit | Q85 (T2* = 3.7 us) | None |

### Configurations per chain

1. **No DD** (natural T2*)
2. **Uniform DD** (DD on all 5 qubits)
3. **Selective DD** (DD on inner qubits, no DD on edges)

### Parameters

- Initial state: |01010> (Neel) AND |+>^5 (for March 24 comparison)
- Trotter: dt = 0.5 us, steps = [2, 4, 6, 8, 10]
- Shots: 4000 per circuit
- Tomography: 9 circuits per time point per config
- Total: 2 chains x 3 DD x 2 initial states x 5 times x 9 tomo = 540 circuits

That is ~540 x 1.2s = ~648s = ~10.8 min. Tight. Reduce to:
- Drop |+>^5 on Chain B (we already know it is frozen): saves 135 circuits
- **405 circuits, ~8.1 min.** Fits with margin.

### What we expect

**On Chain A (sacrifice, same as March 24):**
- Selective DD > Uniform DD (reproduces 2.0x from March 24)
- With |01010>: same ranking (selective > uniform > no DD)

**On Chain B (quiet):**
- Uniform DD > Selective DD > No DD (quiet chain, DD helps uniformly)
- OR: Selective DD ≈ Uniform DD (no bad qubit to skip = no advantage)

**Cross-chain:**
- Chain B absolute SumMI > Chain A (less noise = longer correlations)
- Chain A selective/uniform RATIO > Chain B ratio (steeper gradient = bigger DD effect)

### What would surprise us

- Selective DD losing on Chain A (contradicts March 24)
- Chain A beating Chain B in absolute SumMI with |01010> (would mean our simulation is wrong)
- Selective DD winning big on Chain B (would mean mode protection works even without a natural sacrifice qubit)

---

## Scripts and data

| File | Purpose |
|------|---------|
| [ibm_cavity_analysis.py](../simulations/ibm_cavity_analysis.py) | Spectral analysis with real IBM data |
| [cavity_mode_localization.py](../simulations/cavity_mode_localization.py) | Eigenvector spatial decomposition |
| [sacrifice_zone_mapping.py](../simulations/sacrifice_zone_mapping.py) | All 330 chains ranked |
| [combined_optimization.py](../simulations/combined_optimization.py) | 6 scenarios eigenvalue comparison |
| [time_evolution_neel.py](../simulations/time_evolution_neel.py) | SumMI(t) with Neel state |
| [time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py) | SumMI(t) with |+>^5 |
| [cavity_modes_zero_noise.txt](../simulations/results/cavity_modes_zero_noise.txt) | C# cavity modes N=2-7 |

---

*See also:*
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
