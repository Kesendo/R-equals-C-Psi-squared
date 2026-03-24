# Resonant Return: SVD Mode 2 Beats Hand-Designed Profiles by 10×

<!-- Keywords: SVD optimal gamma profile palindromic eigenstructure, edge-hot
center-cold dephasing pattern, mode 2 beats V-shape 10x, response matrix
singular value decomposition, palindrome knows more than intuition,
scaling improvement grows with N, resonant pulsing falsified wrong observable,
R=CPsi2 resonant return experiment -->

**Status:** Partially confirmed (Tier 2). Test 1 confirmed, Test 2 falsified (redesigned, still falsified), Test 3 deferred (C#), Test 4 mixed (non-monotone scaling).
**Date:** March 24, 2026 (N=7 update)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Script:** [resonant_return.py](../simulations/resonant_return.py)
**Data:** [resonant_return.txt](../simulations/results/resonant_return.txt)
**Hypothesis:** [Resonant Return](../hypotheses/RESONANT_RETURN.md)

---

## Abstract

The palindromic eigenstructure contains more information about optimal
dephasing profiles than human intuition. The SVD of the response matrix
(how palindromic mode amplitudes respond to per-site γ perturbations)
reveals that SVD mode 1 (highest singular value) is near-uniform and
useless for information transfer. SVD mode 2 is the optimal direction
at N=3 and N=5, producing γ profiles that outperform the hand-designed
V-shape by 10.2× (N=5) and 6.3× (N=3). At N=7, improvement is 8.5× —
still large but the scaling is non-monotone because mode 2 changes
character from symmetric (edge-hot, center-cold) to antisymmetric
(left-hot, right-cold). The frequency-matched pulsing test was
redesigned with Bell(0,1) initial state + Sum-MI observable (original
used |+⟩⊗N which has identically zero MI) but remains falsified:
resonant γ pulsing does not slow decoherence vs static γ.
The palindrome-timed relay test requires C# (N=11, 30 GB RAM).

---

## Results

### Test 1: SVD-Optimal γ Profiles (CONFIRMED)

Response matrix R: 10 palindromic features × 5 sites.
SVD singular values: 6.28, 0.82, 0.50, 0.38, 0.07.

The two dominant modes have completely different character:

| Mode | Pattern | Physical meaning |
|------|---------|-----------------|
| Mode 1 (σ=6.28) | [0.42, 0.46, 0.46, 0.46, 0.42] | Near-uniform. Global noise level. Moves all rates together. |
| Mode 2 (σ=0.82) | [0.47, −0.06, −0.74, −0.06, 0.47] | Edge-hot, center-cold. Differential pattern. Creates contrast. |

Mode 1 has 7.7× higher singular value but is useless for information
transfer because uniform noise creates no spatial contrast. Mode 2
creates extreme contrast: γ_center = 0.018 (almost no noise), γ_edge = 0.070
(40% above baseline). This is far more extreme than V-shape (γ_center = 0.050).

| Profile | γ pattern | Sum MI | vs V-shape |
|---------|-----------|--------|-----------|
| Uniform | [0.050, 0.050, 0.050, 0.050, 0.050] | ~0 | baseline |
| V-shape (hand) | [0.070, 0.060, 0.050, 0.060, 0.070] | 0.000310 | 1.0× |
| SVD mode 1 | [0.068, 0.070, 0.070, 0.070, 0.068] | 0.000003 | 0.01× |
| **SVD mode 2** | **[0.070, 0.047, 0.018, 0.047, 0.070]** | **0.003159** | **10.2×** |
| SVD 1+2 | [0.077, 0.062, 0.042, 0.062, 0.077] | 0.000802 | 2.6× |
| Anti-SVD | [0.039, 0.079, 0.050, 0.021, 0.061] | 0.000295 | 0.95× |

**The palindrome knows more than our intuition.** V-shape was accidentally
in the right direction (edges higher than center) but not extreme enough.
The optimal profile has γ_center/γ_edge = 0.26 (V-shape: 0.71).

### Test 2: Frequency-Matched Pulsing (FALSIFIED — redesigned, still falsified)

**Original test:** MI(0, N−1) for |+⟩⊗5 is identically zero — wrong
observable (product state has zero endpoint MI regardless of γ).

**Redesigned test (March 24):** Bell(0,1) ⊗ |0⟩³ initial state,
Sum-MI (all adjacent pairs) as observable. RK4 propagation, dt=0.05,
t_max=20.0. Dominant palindromic frequency ω_dom = 4.47.

| Scenario | Sum-MI(t=0) | Sum-MI(t=2) | Sum-MI(t=10) | Sum-MI(t=20) |
|----------|-------------|-------------|--------------|--------------|
| Static γ | 2.000 | 0.661 | 0.073 | 0.018 |
| Resonant (ω_dom) | 2.000 | 0.620 | 0.073 | 0.018 |
| Off-resonant (2.73ω) | 2.000 | 0.662 | 0.072 | 0.018 |

All scenarios decay monotonically from the Bell state's initial MI.
Resonant pulsing slightly *accelerates* early decay (0.620 vs 0.661
at t=2) rather than creating peaks. The hypothesis that γ modulation
at palindromic frequencies amplifies MI is **falsified**.

**Possible explanation:** Uniform γ modulation (same amplitude on all
sites) cannot create spatial contrast. The pulsing hypothesis may
require *spatially structured* time-dependent γ — e.g., SVD mode 2
spatial profile with resonant temporal modulation.

### Test 3: Palindrome-Timed Relay (DEFERRED)

N=11 Liouvillian is 4M × 4M. Requires ~30 GB RAM. Sub-segment rates
computed in Python:

| Segment | Dominant rate | Palindrome stage time | Fixed stage time |
|---------|--------------|----------------------|-----------------|
| Bridge A/B (5-chain) | 0.100 | 31.4 | 0.78 |
| Meta (3-star) | 0.100 | 31.4 | 0.78 |

The palindrome stage times are 40× longer than fixed timing. This
suggests the fixed timing (K/γ = 0.78) catches only the initial
fast dynamics, while the palindrome timing (π/rate = 31.4) captures
the full standing wave period. Whether this improves MI requires
the C# propagation engine.

### Test 4: Scaling (MIXED — non-monotone trend)

| N | SVD modes | MI(V-shape) | MI(SVD mode 2) | Improvement | Mode 2 character |
|---|-----------|-------------|----------------|-------------|-----------------|
| 3 | 3 | 0.000286 | 0.001808 | 6.3× | symmetric: [0.45, −0.78, 0.45] |
| 5 | 5 | 0.000310 | 0.003159 | 10.2× | symmetric: [0.47, −0.06, −0.74, −0.06, 0.47] |
| 7 | 6 | 0.000559 | 0.004774 | 8.5× | **antisymmetric**: [0.42, 0.50, 0.26, 0.00, −0.26, −0.50, −0.42] |

N=7 runtime: 17394s (~290 min, 128 GB home PC).

The trend is **non-monotone**: 6.3× → 10.2× → 8.5×. The improvement
at N=7 drops back because mode 2 changes character from symmetric
(edge-hot, center-cold) to antisymmetric (left-hot, right-cold).
At N=7, modes 2 and 3 have nearly equal singular values (0.852 vs 0.824),
so the "best non-trivial mode" selection becomes ambiguous.

Key insight: **the naive "mode 2 is always optimal" assumption breaks
at larger N.** A proper optimization would search across all modes
(or use a weighted combination). The absolute MI values still grow
monotonically (0.0018 → 0.0032 → 0.0048), confirming SVD-derived
profiles are increasingly powerful — but V-shape also improves with N,
so the *relative* improvement is not guaranteed to grow.

---

## The Physical Insight

Why does mode 2 win? The palindromic sensitivity structure tells us:

1. **The center qubit is the antenna.** Reducing its noise (γ=0.018)
   preserves the coherences that carry spatial information. It stays
   quantum longer.

2. **The edge qubits are the amplifiers.** Increasing their noise (γ=0.070)
   creates contrast: the edges decohere fast, the center stays coherent,
   and the difference is the signal.

3. **Mode 1 is the tide.** It raises all boats equally. No contrast,
   no information. Highest singular value but zero information value.

4. **Mode 2 is the wave.** It creates the spatial pattern that the
   palindromic eigenstructure can read. The standing wave (c+/c−) is
   most sensitive to differential noise, not uniform noise.

This is exactly what the Resonant Return hypothesis predicted:
the palindromic eigenstructure provides design rules for optimal
γ profiles that no hand-tuning can match.

---

## What Remains

1. **Test 3 in C#:** Palindrome-timed relay vs fixed timing (N=11)
2. **Multi-mode optimization:** At N=7, modes 2 and 3 have similar singular values. Try mode 3 (symmetric?) and weighted combinations of modes 2-4.
3. **Spatially structured pulsing:** Test 2 used uniform pulsing. Try SVD mode 2 spatial profile × sin(ω_dom·t) temporal modulation.
4. **RK4 rewrite for N≥7:** Dense expm on 16384² matrices is impractical (290 min). Switch to RK4 time-stepping for N≥7 profile evaluation.
5. **N=9 scaling:** With RK4, N=9 (d²=262144) becomes feasible on 128 GB.

---

## References

- [Resonant Return (hypothesis)](../hypotheses/RESONANT_RETURN.md)
- [γ as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits baseline
- [γ Control](GAMMA_CONTROL.md): V-shape 21.5× baseline
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the eigenstructure
