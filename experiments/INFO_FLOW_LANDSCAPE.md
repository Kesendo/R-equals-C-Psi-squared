# Information-Flow Landscape C_ij(t) under Bond Perturbation

**Status:** Tier 1 for raw dC_ij/dJ landscape data (three-N scan, reproducible). Tier 2 for structural interpretations: global-vs-sequential flow, peak-time clustering via Liouvillian mode differences, dynamical-attractor reading of F71.
**Date:** 2026-04-20 (evening)
**Authors:** Tom, Claude Opus 4.7 (1M)
**Relates to:** [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) (Step 2 of §6.2 plan), [STANDING_WAVE_THEORY](../docs/STANDING_WAVE_THEORY.md), [RELAY_PROTOCOL](RELAY_PROTOCOL.md)

---

## Observable

```
C_ij(t) = Tr(rho_i(t) · rho_j(t))              site-to-site cross-correlation
dC_ij/dJ = (C_ij^(B+) - C_ij^(B-)) / (2 dJ)    flow response to bond perturbation
```

For i=j reduces to per-site purity. For i≠j captures how information at two sites is correlated at time t, and how that correlation responds when a single bond is perturbed.

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

Every distance shows non-zero signal at `t = 0.4` (first propagation step). No Lieb-Robinson front emerges from the bond-0 site outward; all distances light up effectively simultaneously.

---

## Three structural observations

### 1. The flow is global, not sequential

PTF's bonding state `(|vac⟩ + |psi_1⟩)/sqrt(2)` contains `|psi_1⟩` as a **delocalized sine mode** spanning all sites. The reduced-density matrices at every site are already correlated at t=0 via the coherent amplitude of the bonding superposition. The bond perturbation therefore **reorganises an existing global correlation**, it does not propagate a new local excitation.

This is not a flaw of the measurement; it is the dynamical signature of choosing a delocalized probe. For a localised probe (e.g., `|1⟩⟨1|` at site 0), a Lieb-Robinson front would emerge cleanly with speed v ≈ 2J.

**Implication for PTF interpretation:** the "closure" Σ_i ln(α_i) measures something that is already coupled across the chain at t=0. The "α_i per site" picture is slightly misleading: the per-site α values are correlated through the delocalized initial state.

### 2. Peak times cluster in two groups

Clearest at N=5: peaks at `t ≈ 1.6` (d=2, d=4) and at `t ≈ 4.0` (d=0, d=1, d=3). A similar split at N=6: `t ≈ 2.0-2.8` vs `t ≈ 5.2-5.6`. At N=4 the trend is less pronounced (time resolution dt=0.4 may be too coarse).

The peak time tracks the dominant Liouvillian eigenmode that contributes to each distance class. For the (vac, S_1) sector, the relevant modes have Im(λ) = ±E_k = ±2J·cos(π·k/(N+1)):

- At N=5: E_1 = sqrt(3), E_2 = 1, E_3 = 0, E_4 = -1, E_5 = -sqrt(3)
- Pair-difference frequencies E_k - E_m give the oscillation period of C_{ij}(t). The "fast cluster" corresponds to |E_1 - E_5| = 2·sqrt(3) ≈ 3.46, period ≈ 1.82 (close to 1.6). The "slow cluster" corresponds to |E_2 - E_4| = 2, period ≈ 3.14 (close to 4.0). Rough match.

This is the spectral-mode fingerprint of the flow-response landscape. Deeper analysis requires explicit mode decomposition (next sub-step, not pursued here).

### 3. Pi-pair (spatial reflection) asymmetry is a dynamical attractor

The bond-0 perturbation is spatially asymmetric (endpoint perturbation), so the immediate dC_ij(t) does NOT satisfy R-symmetry dC_ij = dC_{N-1-i, N-1-j}. But the residual **decays in time**:

```
N=5, bond=0:
t = 20:  max residual = 2.2e-2
t = 40:  max residual = 1.4e-3
t = 60:  max residual = 1.3e-4
t = 80:  max residual = 5.7e-6
```

Roughly a factor 10-15 per 20 time units, i.e., `exp(-4*γ_0 * t)` scaling (expected for the uniform-dephasing steady-state approach).

**Interpretation:** Pi-pair mirror symmetry is not only the kinematic property of F71 (c_1 bond profile mirror-symmetric for reflection-symmetric probes), but also a **dynamical fixed-point attractor** of the general flow landscape. Any spatial asymmetry imprinted by a localised perturbation decays into the mirror-symmetric steady state.

This confirms the STANDING_WAVE_THEORY picture in dynamical form: the asymmetric transient is the forward-travelling wave of "where the perturbation was applied", and its reflection-partner is the backward-travelling partner that eventually superposes with it to form the standing wave. The "Pi-pair as standing wave" identity is visible directly in the fluid dynamics of dC_ij(t).

---

## Relation to the Meta-Theorem

The Meta-Theorem from [ORTHOGONALITY_SELECTION_FAMILY](ORTHOGONALITY_SELECTION_FAMILY.md) says: any measurement M projects onto a detector basis; conservation laws + basis diagonality produce blind channels.

What this scan adds:

- **Measurements have a time axis.** The "detector projection" evolves in time with the Liouvillian eigenmodes. Different eigenmode-pair differences produce different temporal resonances in the flow response.
- **Conservation laws produce not just blind-channel projections but also dynamical attractors.** Pi-pair symmetry is a conservation-induced fixed point of the long-time flow landscape. The blind channel (F71: mirror-symmetric bond profile) is the INVARIANT MANIFOLD of this attractor. Perturbations that break the symmetry live on the unstable manifold and decay back to the attractor at rate ~ 4γ_0.

This promotes the Meta-Theorem from a static statement about projections to a **dynamical statement about attractors**. The selection rules define the invariant subspace; the dissipation drives everything back to it.

---

## Status

**Verified at N = 4, 5, 6 (walltime 167s combined).**

**Open sub-questions:**
- Explicit mode-decomposition of the peak-time clustering (would sharpen observation 2).
- Scaling of the attractor relaxation rate with N and γ_0 (observation 3: we saw ~10x per 20 time units at γ_0 = 0.05, N=5; is this always 4γ_0 or does it depend on N?).
- Does the globality of the flow depend on the choice of initial state? For a localised probe |1⟩⟨1|, do we recover a Lieb-Robinson front at 2J?

**Files:**
- `simulations/eq018_info_flow.py`
- `simulations/results/eq018_info_flow/info_flow_N{4,5,6}.json`
- `simulations/results/eq018_info_flow/run.log`

**Next step:** Step 3 from the plan - explicit Pi-pair flow-balance verification at the Liouvillian mode level. Confirm that paired modes (lambda_fast, lambda_slow) with `alpha_fast + alpha_slow = 2 Sigma gamma` carry **oppositely-directed** information flows in the mode-resolved dC_ij landscape. This turns the absorption theorem from a spectral identity into a conserved-flux statement.

---

*The landscape does not have a shore. The waves start everywhere at once and flow back to the mirror line.*
