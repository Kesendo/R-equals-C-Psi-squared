# One-Site Marginal-Similarity Landscape C_ij(t) under Bond Perturbation

**Status:** Tier 1 for raw `dC_ij/dJ` landscape data (three-N scan,
reproducible). `C_ij` is a similarity of one-site marginals, not a two-site
correlation or information-flow observable; propagation and modal attribution
are not established by this scan.
**Date:** 2026-04-20 (evening)
**Authors:** Tom, Claude Opus 4.7 (1M)
**Relates to:** [the orthogonality-selection family](ORTHOGONALITY_SELECTION_FAMILY.md) (Step 2 of §6.2 plan), [standing wave theory](../docs/STANDING_WAVE_THEORY.md), [the relay protocol](RELAY_PROTOCOL.md)

---

## Observable

```
C_ij(t) = Tr(rho_i(t) · rho_j(t))              one-site marginal similarity
dC_ij/dJ = (C_ij^(B+) - C_ij^(B-)) / (2 dJ)    response to bond perturbation
```

For `i=j` this reduces to per-site purity. For `i≠j` it compares two reduced
one-site density matrices. It is not the expectation of a joint observable and
does not measure correlation or information shared between the sites.

Initial state: PTF bonding `(|vac⟩ + |psi_1⟩)/sqrt(2)`. Perturbation: bond 0, delta_J = ±0.01. Scanned N = 4, 5, 6.

---

## Results

Peak |dC_ij/dJ| and peak time, summarised over all site pairs at distance `d = |i-j|`:

```
                N=4              N=5              N=6
d   peak     time    peak     time    peak     time
0   0.137    2.80    0.116    4.40    0.102    5.20
1   0.088    2.80    0.099    4.00    0.079    5.60
2   0.081    2.80    0.068    1.60    0.068    2.40
3   0.040    4.40    0.080    4.00    0.069    2.80
4                    0.051    1.60    0.069    5.60
5                                     0.046    2.00
```

Every distance has a non-zero derivative at `t = 0.4`, the first sampled time.
Because the initial probe is delocalized and `C_ij` is a marginal-similarity
quantity, this observation is not a test for a Lieb-Robinson front.

---

## Three structural observations

### 1. The delocalized probe does not test propagation

PTF's bonding state `(|vac⟩ + |psi_1⟩)/sqrt(2)` contains `|psi_1⟩` as a
**delocalized sine mode** spanning all sites. A bond perturbation can therefore
change one-site marginals across the chain without launching a local excitation.

The scan therefore reports sensitivity of marginal similarity to the bond. A
front test would need a localized preparation and a propagation observable; no
front speed follows from these data.

### 2. Peak times cluster in two groups

Clearest at N=5: peaks at `t ≈ 1.6` (d=2, d=4) and at `t ≈ 4.0` (d=0, d=1, d=3). A similar split at N=6: `t ≈ 2.0-2.8` vs `t ≈ 5.2-5.6`. At N=4 the trend is less pronounced (time resolution dt=0.4 may be too coarse).

The peak-time clusters are roughly consistent with pair differences among the
single-excitation frequencies Im(λ) = ±E_k = ±2J·cos(π·k/(N+1)):

- At N=5: E_1 = sqrt(3), E_2 = 1, E_3 = 0, E_4 = -1, E_5 = -sqrt(3)
- Pair-difference frequencies E_k - E_m give the oscillation period of C_{ij}(t). The "fast cluster" corresponds to |E_1 - E_5| = 2·sqrt(3) ≈ 3.46, period ≈ 1.82 (close to 1.6). The "slow cluster" corresponds to |E_2 - E_4| = 2, period ≈ 3.14 (close to 4.0). Rough match.

This is a rough numerical comparison, not a modal attribution. An explicit
mode decomposition was not performed.

### 3. Spatial-reflection asymmetry decays in the measured scan

The bond-0 perturbation is spatially asymmetric (endpoint perturbation), so the immediate dC_ij(t) does NOT satisfy R-symmetry dC_ij = dC_{N-1-i, N-1-j}. But the residual **decays in time**:

```
N=5, bond=0:
t = 20:  max residual = 2.2e-2
t = 40:  max residual = 1.4e-3
t = 60:  max residual = 1.3e-4
t = 80:  max residual = 5.7e-6
```

The sampled reductions are about 10–22× per 20 time units. Four late-time
samples do not determine an exponential rate; its dependence on `N` and
`γ_0` remains open.

**Interpretation:** F71 gives a mirror-symmetric `c_1` bond profile for
reflection-symmetric probes. In this endpoint-perturbation scan, the measured
reflection asymmetry decays as the trajectory approaches its symmetric steady
state.

What is visible here is the decay of reflection asymmetry in `dC_ij(t)`. The
simultaneous nonzero derivatives do not test a propagating front, and neither
these data nor the Pi pairing identify forward/backward spatial modes or a
standing wave; that reading would require a localized preparation and a
propagation or current observable.

---

## Relation to the Meta-Theorem

The Meta-Theorem from [the orthogonality-selection family](ORTHOGONALITY_SELECTION_FAMILY.md) says: any measurement M projects onto a detector basis; conservation laws + basis diagonality produce blind channels.

What this scan records:

- **Measurements have a time axis.** The response contains multiple temporal
  scales roughly consistent with eigenvalue-difference frequencies; an explicit
  mode decomposition has not been performed.
- **The static and dynamical observations are distinct.** F71 supplies the
  reflection-symmetric profile for symmetric probes. This scan additionally
  observes late-time decay of one endpoint perturbation's reflection asymmetry,
  without establishing a universal attractor, invariant manifold, or decay rate.

---

## Status

**Raw derivative landscape reproduced at N = 4, 5, 6 (walltime 167s combined).**

**Open sub-questions:**
- Explicit mode-decomposition of the peak-time clustering (would sharpen observation 2).
- Dependence of the reflection-residual decay on N and `γ_0`.
- A separate propagation test with a localized preparation and a suitable
  local observable.

**Files:**
- `simulations/eq018_info_flow.py`
- `simulations/results/eq018_info_flow/info_flow_N{4,5,6}.json`
- `simulations/results/eq018_info_flow/run.log`

**Next step:** Step 3 from the plan - test at the Liouvillian mode level whether
paired modes (`alpha_fast + alpha_slow = 2 Sigma gamma`) carry any
oppositely-directed current in the mode-resolved `dC_ij` landscape. The spectral
identity itself supplies no conserved-flux statement.

---

*On this delocalized probe, the response starts globally and its measured
reflection asymmetry decays.*
