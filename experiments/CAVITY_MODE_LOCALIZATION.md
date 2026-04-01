# Cavity Mode Localization: Where the Protected Modes Live

<!-- Keywords: Liouvillian eigenvector spatial localization, Pauli basis
decomposition qubit weight profile, sacrifice zone mode protection
mechanism, palindromic partner edge center topology, chain interior
modes survive dephasing, R=CPsi2 cavity mode localization -->

**Status:** Tier 2 (computed eigenvector analysis, three noise profiles, N=5)
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md)
**Script:** [cavity_mode_localization.py](../simulations/cavity_mode_localization.py)
**Data:** [cavity_mode_localization.txt](../simulations/results/cavity_mode_localization.txt)

---

## Abstract

The sacrifice zone protects cavity modes, not qubits
([IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md)). But which
modes are protected, and why?

By decomposing Liouvillian eigenvectors into the Pauli basis, we
compute the spatial profile of each mode: how much of the mode "lives"
on each qubit. The correlation between sacrifice-qubit weight and decay
rate is **r = 0.9942** (p = 0.00). Modes concentrated on the chain
center survive longest. Modes concentrated on the edges die fastest.

The spatial profile is **geometric**: the same mode has the same
qubit-weight pattern under IBM sacrifice noise and uniform noise. The
geometry of the chain (standing wave patterns) determines WHERE modes
live. The noise profile determines WHICH modes survive. The sacrifice
zone exploits a pre-existing geometric structure. (Not topologically
protected in the Altland-Zirnbauer sense: winding number = 0, Berry
phase not quantized. See
[Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md).)

---

## Method

For each Liouvillian eigenvector v (1024 components for N=5), we:

1. Reshape v into a 32x32 operator V
2. Decompose V into the Pauli basis: c_j = Tr(P_j^dag V) / d
3. Compute per-qubit weight: weight[k] = Sum |c_j|^2 over all Pauli
   strings where qubit k is active (not identity)
4. Normalize to a probability distribution over qubits

---

## Results

### The 7.234J modes: center-localized

The 4 slowest oscillating modes (rate = 0.0456, freq = 7.234J) have
the profile:

| Q0 (sacrifice) | Q1 | Q2 (center) | Q3 | Q4 |
|----------------|------|------------|------|------|
| 0.519 | 0.631 | **0.700** | 0.631 | 0.519 |

Weight peaks at the chain center (Q2 = 0.700) and is minimal at the
edges (Q0 = Q4 = 0.519). These modes are **center-localized**.

Their palindromic partners (rate = 0.594, fastest oscillating) have
the complementary profile:

| Q0 | Q1 | Q2 | Q3 | Q4 |
|------|------|------|------|------|
| **0.981** | 0.869 | 0.800 | 0.869 | **0.981** |

Weight peaks at the edges. The slow mode lives in the center. The
fast mode lives on the edges. **Palindromic partners are spatially
complementary.**

### Correlation: Q0 weight vs decay rate

| Profile | Pearson r | p-value |
|---------|----------|---------|
| IBM sacrifice | 0.9942 | 0.00 |
| Uniform | 0.8588 | 2e-297 |
| Zero noise | 0.0169 | 0.61 |

Under IBM sacrifice noise: near-perfect correlation. Modes with more
weight on Q0 (the sacrifice qubit) decay faster. Under uniform noise:
still strong (the topology effect persists but all qubits contribute
equally to damping). At zero noise: no correlation (no damping, no
selection).

### Quartile comparison

**IBM sacrifice:**
- Slowest quartile (254 modes): [0.666, 0.711, 0.716, 0.712, 0.675]
- Fastest quartile (254 modes): [0.839, 0.795, 0.786, 0.795, 0.830]

Slow modes avoid the edges. Fast modes concentrate on the edges. The
difference between slowest and fastest is 0.17 on Q0 (25% relative).

**Zero noise:**
- Slowest quartile: [0.750, 0.750, 0.746, 0.750, 0.750]
- Fastest quartile: [0.748, 0.749, 0.750, 0.749, 0.748]

No spatial preference at all. The profiles are flat.

### The profile is geometric (not topologically protected)

The 4 slowest 7.234J modes have profile [0.519, 0.631, 0.700, 0.631, 0.519]
under **both** IBM sacrifice and uniform noise. The profile does not
depend on the noise distribution. It depends on the chain geometry.

The chain creates modes with different spatial profiles through
standing wave patterns (sin(πkj/N) eigenstates of the tight-binding
w=1 sector). Some modes concentrate on the center, others on the
edges. This is a property of the Hamiltonian, not the noise. When
asymmetric noise is applied, it selectively damps the edge-heavy
modes (because the sacrifice qubit is on the edge). The center-heavy
modes survive. Crucially, which modes are "slowest" REVERSES when
noise moves from edge to center: edge sacrifice selects center-heavy
modes, center sacrifice selects edge-heavy modes
([Topological Edge Modes](TOPOLOGICAL_EDGE_MODES.md)).

---

## The mechanism

1. The Heisenberg chain creates modes with varying spatial profiles
   (center-heavy to edge-heavy)
2. At zero noise, all modes have zero decay rate (no selection)
3. When noise is turned on, each mode's decay rate depends on how
   much of it touches noisy qubits
4. Under uniform noise: all modes see the same noise (no differential
   protection, r = 0.86 from topology alone)
5. Under sacrifice noise: edge modes see much more noise than center
   modes (r = 0.994, strong differential protection)
6. The slowest modes happen to be center-localized at the 7.234J
   frequency (a topological property of the chain)

The sacrifice zone does not create the spatial structure. It
**exploits** it. The chain topology provides the modes. The sacrifice
zone provides the selection pressure. Together: 2.81x protection
for the modes that matter most.

---

*See also:*
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md) (the 2.81x result),
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md) (the eigenfrequencies),
[Resonant Return](RESONANT_RETURN.md) (the sacrifice-zone formula),
[Energy Partition](../hypotheses/ENERGY_PARTITION.md) (2x decay law)
