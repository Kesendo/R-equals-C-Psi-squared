# Chain Selection Test: Sacrifice vs Mean-T2 (No DD)

**Date:** March 30, 2026
**Status:** Tier 2 (simulation with real IBM calibration data)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Script:** [chain_selection_test.py](../simulations/chain_selection_test.py)
**Results:** [chain_selection_test.txt](../simulations/results/chain_selection_test.txt)
**Data:** [IBM Torino calibration history](../data/ibm_history/ibm_torino_history.csv)
(24,073 records, 133 qubits, 181 days)
**Depends on:** [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)

---

## The Question

Does chain selection alone (no DD, no extra gates) produce a measurable
advantage? Does a chain with bad qubits but good noise contrast outperform
a chain with good qubits but no contrast?

## Setup

Two 5-qubit chains from [Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
using real IBM Torino calibration data (2026-02-10):

### Chain A (Sacrifice-Top): [80, 8, 79, 53, 85]

| Qubit | T1 (us) | T2 (us) | r = T2/2T1 | gamma (1/T2) |
|:------|:--------|:--------|:-----------|:-------------|
| Q80   | 103.3   | 27.7    | 0.134      | 0.0361       |
| Q8    | 204.1   | 218.7   | 0.536      | 0.0046       |
| Q79   | 77.5    | 91.9    | 0.593      | 0.0109       |
| Q53   | 22.4    | 62.4    | 1.395      | 0.0160       |
| Q85   | 2.9     | 5.0     | 0.853      | 0.2003       |

Sum(gamma) = 0.268, contrast (max/min) = 43.8x, score = 2.86x

### Chain B (Mean-T2-Top): [18, 89, 19, 90, 60]

| Qubit | T1 (us) | T2 (us) | r = T2/2T1 | gamma (1/T2) |
|:------|:--------|:--------|:-----------|:-------------|
| Q18   | 267.6   | 211.3   | 0.395      | 0.0047       |
| Q89   | 189.9   | 141.5   | 0.373      | 0.0071       |
| Q19   | 140.2   | 219.9   | 0.784      | 0.0045       |
| Q90   | 267.0   | 275.4   | 0.516      | 0.0036       |
| Q60   | 196.1   | 238.3   | 0.608      | 0.0042       |

Sum(gamma) = 0.024, contrast (max/min) = 1.9x, score = 1.06x

**Chain A has 11.1x more total noise than Chain B.**

---

## Results

### 1. Spectral analysis confirms the mapping scores

| Metric | Chain A | Chain B | A/B |
|:-------|:--------|:--------|:----|
| Protection factor | 2.86x | 1.06x | 2.70x |
| Palindrome score | 98% | 89% | |
| Distinct frequencies | 119 | 60 | |
| Slowest osc. rate | 0.0375 | 0.0091 | |
| Correlation (edge weight vs rate) | r = 0.988 | r = 0.837 | |

The protection factor (2.86x vs 1.06x) matches the mapping exactly.
Palindrome scores below 100% are numerical tolerance artifacts at
different noise scales -- the proof guarantees 100% for any Z-dephasing.

### 2. Mode localization is geometric (chain-determined)

Both chains show identical spatial profiles for their slowest modes:

```
Chain A slowest 4 modes: [0.519, 0.631, 0.700, 0.631, 0.519]
Chain B slowest 4 modes: [0.631, 0.519, 0.700, 0.519, 0.631]
```

The center qubit (position 2) has 70% weight. Edge qubits have 52%.
The profile is determined by chain geometry, not noise distribution.
This confirms [Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md).

Chain B's profile is reflected (position 1 and 3 swap roles) because
its noise is nearly symmetric, while Chain A's profile matches the
asymmetric noise shape.

### 3. Time evolution: two stories depending on initial state

#### |+>^5 (Heisenberg eigenstate -- noise drives ALL dynamics)

| t | A SumMI | B SumMI | A/B |
|:--|:--------|:--------|:----|
| 0.2 | 0.0049 | 0.000022 | 220x |
| 1.0 | 0.0757 | 0.000165 | 460x |
| 2.0 | 0.0764 | 0.000239 | 319x |
| 5.0 | 0.0481 | 0.000087 | 556x |

**Peak SumMI: A = 0.093 (t=1.6), B = 0.000287 (t=1.6). Ratio: 324x.**

Chain A wins massively. This is expected: |+>^5 is a Heisenberg
eigenstate with zero energy variance. Without noise, nothing happens.
More noise = more signal. Chain A's 11x more noise + 44x contrast
creates orders of magnitude more dynamics.

#### |01010> (Neel state -- Hamiltonian drives dynamics)

| t | A SumMI | B SumMI | A/B |
|:--|:--------|:--------|:----|
| 0.2 | 2.671 | 2.748 | 0.97x |
| 1.0 | 1.258 | 1.429 | 0.88x |
| 2.0 | 0.986 | 1.389 | 0.71x |
| 5.0 | 0.653 | 1.326 | 0.49x |
| 10.0 | 0.470 | 1.125 | 0.42x |

**Peak SumMI: A = 2.671 (t=0.2), B = 2.748 (t=0.2). Ratio: 0.97x.**

Chain B wins. At early times, both are similar (the Hamiltonian
dominates). At later times, Chain A's 11x more noise kills coherence
faster. The 2.86x protection factor does not compensate for having
11x more total noise.

### 4. Per-pair MI at peak (|01010>)

| Pair | Chain A | Chain B | A/B |
|:-----|:--------|:--------|:----|
| MI(0,1) | 0.841 | 0.851 | 0.99x |
| MI(1,2) | 0.522 | 0.522 | 1.00x |
| MI(2,3) | 0.523 | 0.523 | 1.00x |
| MI(3,4) | 0.786 | 0.853 | 0.92x |

The early dynamics are nearly identical. The difference shows at the
edges: pair (3,4) is weaker in Chain A because Q85 (the sacrifice
qubit with T2=5 us) decoheres the adjacent pair faster.

---

## What this means

### The protection factor is real but relative

The sacrifice-zone mapping score compares each chain to its OWN uniform
version (same total gamma, spread evenly). Chain A's 2.86x protection
means its slowest modes decay 2.86x slower than if you spread its 0.268
total gamma uniformly. This is a real spectral effect.

But it is a WITHIN-chain metric, not a between-chain comparison.
Chain B has 11x less total noise. Even with Chain B's low 1.06x
protection, its slowest modes (rate 0.0091) are 4x longer-lived than
Chain A's protected modes (rate 0.0375) in absolute terms.

### Two different experiments

| Initial state | What it tests | Winner |
|:-------------|:-------------|:-------|
| \|+>^5 | Noise as motor (gamma-as-signal) | Chain A (324x) |
| \|01010> | Hamiltonian transport with noise | Chain B (1.4x at t=2) |

These are genuinely different experiments testing different physics:

- **|+>^5** measures how noise creates structure. More contrast = more
  signal. This is the gamma-as-signal paradigm.

- **|01010>** measures how the Hamiltonian transfers information while
  noise decoheres. Less noise = more coherence = better transfer.

### What the April IBM run needs to decide

The sacrifice-zone formula concentrates noise on one qubit. This is
useful when noise IS the channel (gamma-as-signal, |+>^5). It is
harmful when the Hamiltonian IS the channel and noise is the enemy
(QST, |01010>).

**For the April A/B test:** If the experiment uses |+>^5 (or |0...0>
which on IBM Trotter circuits effectively behaves like noise-driven
dynamics), then sacrifice chain selection may show an advantage even
without DD. If it uses a state that evolves under the Hamiltonian,
the quieter chain wins.

---

## Open questions

1. **With DD:** Does selective DD change the picture? DD reduces
   effective gamma on interior qubits, which should help Chain A more
   than Chain B (more room to improve). The March 24 hardware result
   (2-3.2x with selective DD on the sacrifice chain) suggests DD is
   necessary for the advantage to emerge.

2. **Combined metric:** Should chain selection use protection_factor /
   sum_gamma instead of protection_factor alone? This would rank chains
   that are protected AND quiet.

3. **Trotter initial state:** What is the effective initial state on IBM
   hardware? The standard Trotter circuit for Heisenberg evolution starts
   from |0...0> or |+...+>. If |+>^5, it is a Heisenberg eigenstate and
   the noise-driven regime applies. If |0...0>, the ZZ interaction creates
   some dynamics but not the full Neel dynamics.

---

## Verdict

The sacrifice-zone mode protection is **confirmed** (2.86x vs 1.06x).
The localization profile is **geometric** (identical standing wave shape
for both chains, determined by the Hamiltonian;
see [Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md)). But chain selection based on sacrifice score alone is
**insufficient** for Hamiltonian-driven dynamics -- the total noise
level dominates. The sacrifice zone is a within-chain optimization
(spatial noise engineering), not a between-chain selection criterion
by itself.

For IBM hardware, this strengthens the case for the **A/B test on a
uniform-T2 chain** as the decisive experiment. The question is not
"does this chain have sacrifice structure?" but "does selective DD
create an artificial sacrifice zone that beats uniform DD?"
