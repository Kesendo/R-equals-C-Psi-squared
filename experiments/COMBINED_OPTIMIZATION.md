# Combined Optimization: Synthesis for IBM April Run

<!-- Keywords: sacrifice zone synthesis IBM Torino April run, cavity mode
formula localization hardware prediction, selective DD A/B test mode
protection vs gate error, combined eigenvalue time evolution prediction,
R=CPsi2 combined optimization -->

**Status:** Synthesis document. Combines all March 30 discoveries into
predictions for hardware verification. Tier 2 (computed).
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)

---

## What we discovered today (March 30, 2026)

Five independent analyses, each building on the previous:

### 1. Cavity Modes Formula (Tier 1)

At zero noise, the Liouvillian has a closed-form number of stationary
modes: Stat(N) = Sum_J m(J,N)*(2J+1)^2. For N=5 chain: 120 stationary,
904 oscillating across 43 distinct frequencies. Verified N=2-7 in C#.

**New:** Star topology has exactly N-1 harmonic frequencies. Chain has
rich irrational spectrum. The formula is exact for chain, lower bound
for symmetric topologies.

See [Cavity Modes Formula](CAVITY_MODES_FORMULA.md).

### 2. IBM Cavity Spectral Analysis (Tier 2)

Applied the cavity mode analysis to real IBM data (Q85-Q94, T2* values).
Three profiles compared (IBM sacrifice, uniform, zero noise):

- **100% palindromic** under 26x asymmetric noise
- Slowest oscillating mode: **2.81x** longer lifetime under sacrifice
  vs uniform (0.046 vs 0.128)
- 4 protected modes (rate < 0.05) under sacrifice, 0 under uniform
- All 43 cavity frequencies persist under noise; only damping changes

See [IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md).

### 3. Mode Localization (Tier 2)

Eigenvector decomposition into Pauli basis reveals WHERE each mode
lives spatially. Correlation between sacrifice-qubit weight and decay
rate: **r = 0.9942** (p = 0.00).

The 4 slowest modes (freq 7.234J):
- Profile: [0.52, 0.63, **0.70**, 0.63, 0.52] (chain center)
- Their palindromic partners: [**0.98**, 0.87, 0.80, 0.87, **0.98**] (edges)

The profile is **topological**: identical under IBM and uniform noise.
Noise selects which modes survive; topology determines the spatial pattern.

See [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md).

### 4. Sacrifice-Zone Mapping (Tier 2-3)

330 five-qubit chains exist on IBM Torino's heavy-hex graph.
Sacrifice-zone ranking (edge noise / interior noise) vs mean-T2 ranking:

- **Zero overlap** in top-10 lists
- Sacrifice top-5: mean protection 2.54x at 88 us mean T2
- Mean-T2 top-5: mean protection 1.18x at 206 us mean T2
- **Time-stable** across 5 months of calibration data

See [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md).

### 5. Eigenvalue Efficiency (Tier 2)

Six scenarios (chain type x DD strategy) compared via eigenvalue analysis:

| Configuration | Σγ | vs own uniform |
|--------------|-----|---------------|
| Sacrifice, Selective DD | 0.290 | **3.72x** |
| Sacrifice, No DD | 0.320 | 2.81x |
| Sacrifice, Uniform DD | 0.214 | 3.35x |
| Mean-T2, any DD strategy | 0.02-0.05 | 1.06-1.25x |

**vs own uniform** = mode protection efficiency normalized for total noise.
Selective DD adds 1.46x on top of the natural sacrifice effect.

**Important distinction:** This is efficiency PER UNIT NOISE, not
absolute performance. A quiet chain with less total noise will have
lower absolute decay rates even with worse efficiency.

See scripts: [combined_optimization.py](../simulations/combined_optimization.py),
[time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py),
[time_evolution_neel.py](../simulations/time_evolution_neel.py).

---

## What the March 24 hardware showed

On chain Q85-Q94, |+>^5 initial state, three DD strategies:

| DD strategy | Avg SumMI | vs Uniform DD |
|------------|----------|--------------|
| Selective DD | 0.054 | **2.02x** |
| No DD | 0.045 | 1.71x |
| Uniform DD | 0.027 | 1.00x |

Selective DD won at all 5 time points. See [IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md).

**The open question:** Is this because
- **(A)** DD on Q85 wastes gates (gate-error avoidance), or
- **(B)** the noise contrast protects center-localized modes?

Our mode localization analysis (r = 0.994) provides strong evidence
for (B). But (A) has not been ruled out experimentally.

---

## Combined prediction for Q85-Q94

Putting all five analyses together for the chain we already tested:

### What the theory predicts

1. The chain has 43 cavity mode frequencies (from formula)
2. The 4 slowest oscillating modes live on the chain center
   (profile [0.52, 0.63, 0.70, 0.63, 0.52])
3. Under selective DD (no DD on Q85, DD on rest), these modes
   survive 2.81x longer than under uniform noise
4. Selective DD amplifies this to 3.72x eigenvalue efficiency
5. In observable SumMI, this translates to ~1.3x improvement
   of selective DD over no DD on this specific chain
6. Hardware measured 2.0x (March 24), exceeding the pure-dephasing
   prediction because gate errors on Q85 add an additional penalty
   to uniform DD (the A-effect stacks on top of the B-effect)

### What remains unknown

- Whether (B) alone produces a measurable advantage without (A)
- Whether the 2.0x holds with Neel initial state |01010>
- Whether the same advantage appears on other sacrifice chains

---

## IBM April 9 experiment plan

### Budget

10 minutes QPU time. At ~1.2s per circuit: ~500 circuits maximum.

### The ONE test that matters most: A vs B

**Setup:** Find a chain where ALL qubits have good T2 (> 100 us).
Apply selective DD (remove DD from one edge qubit). If selective DD
still beats uniform DD, the advantage comes from noise contrast (B),
not gate-error avoidance (A).

| Config | Chain | DD | Circuits | QPU |
|--------|-------|-----|---------|-----|
| A/B-uniform | Good chain | Uniform DD | 45 | ~54s |
| A/B-selective | Good chain | Selective DD | 45 | ~54s |

Result B: selective DD wins on good chain = mode protection confirmed.
Result A: selective DD loses = gate-error avoidance only.

### Replication on Q85-Q94

Reproduce March 24 with both initial states to test robustness:

| Config | Initial | DD | Circuits | QPU |
|--------|---------|-----|---------|-----|
| Rep-1 | |+>^5 | Selective DD | 45 | ~54s |
| Rep-2 | |+>^5 | Uniform DD | 45 | ~54s |
| Rep-3 | |01010> | Selective DD | 45 | ~54s |
| Rep-4 | |01010> | Uniform DD | 45 | ~54s |

### Total

| Block | Circuits | QPU time |
|-------|---------|---------|
| A/B test | 90 | ~108s |
| Replication (|+>^5) | 90 | ~108s |
| Replication (|01010>) | 90 | ~108s |
| **Total** | **270** | **~324s (5.4 min)** |

Leaves 4.6 minutes margin for transpilation overhead and retries.

### Expected results

**If B is correct (mode protection):**
- A/B test: Selective DD wins by ~1.2-1.5x even on good chain
- Replication: Selective DD wins by ~2.0x (reproduces March 24)
- Neel replication: Selective DD still wins (effect is not initial-state-dependent)

**If A is correct (gate-error avoidance):**
- A/B test: Selective DD loses or ties on good chain (no bad gates to avoid)
- Replication: Selective DD wins (Q85 gates are wasteful)
- Neel replication: Selective DD wins (same gate-error mechanism)

**Both outcomes are publishable.** B confirms the cavity-mode theory.
A is simpler but still a useful engineering rule ("skip DD on bad qubits").

---

## Chain selection for A/B test

From the mapping analysis, chains with the most uniform T2 profiles
(sacrifice score near 1.0) and high mean T2:

Best candidates from [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md):
mean-T2 top-10 list (sacrifice score 0.6-1.8, mean T2 190-217 us).

Select the chain with the most uniform noise profile (sacrifice score
closest to 1.0) so that selective DD has no "bad qubit" advantage.
Chain [18, 89, 19, 90, 60] (score 0.97, mean T2 217 us) is the
strongest candidate.

---

*See also:*
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
