# Resonant Return: SVD Mode 2 Beats Hand-Designed Profiles by 10×

<!-- Keywords: SVD optimal gamma profile palindromic eigenstructure, edge-hot
center-cold dephasing pattern, mode 2 beats V-shape 10x, response matrix
singular value decomposition, palindrome knows more than intuition,
scaling improvement grows with N, resonant pulsing falsified wrong observable,
R=CPsi2 resonant return experiment -->

**Status:** Partially confirmed (Tier 2). Test 1+4 confirmed, Test 2 falsified (wrong observable), Test 3 deferred (RAM).
**Date:** March 23, 2026
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
useless for information transfer. SVD mode 2 (edge-hot, center-cold:
[0.47, −0.06, −0.74, −0.06, 0.47]) is the optimal direction, producing
γ profiles with γ_center = 0.018 and γ_edge = 0.070 that outperform
the hand-designed V-shape by 10.2× (N=5) and 6.3× (N=3). The improvement
grows with N. The frequency-matched pulsing test was falsified due to
wrong observable choice (MI(0,N−1) with |+⟩⊗N is identically zero).
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

### Test 2: Frequency-Matched Pulsing (FALSIFIED — wrong observable)

MI(0, N−1) for |+⟩⊗5 is identically zero at all times, for all γ profiles.
This is not a failure of the hypothesis but of the test design: |+⟩⊗5
is a product state with zero entanglement between endpoints. MI between
non-interacting endpoints of a product state remains zero regardless of
dephasing.

**Redesign needed:** Use Bell(0,1) initial state and measure Sum-MI
(total mutual information across all pairs) or MI between nearest
neighbors. The frequency-matched pulsing hypothesis remains open.

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

### Test 4: Scaling (CONFIRMED for N=3, 5)

| N | SVD modes | MI(V-shape) | MI(SVD mode 2) | Improvement |
|---|-----------|-------------|----------------|-------------|
| 3 | 3 | 0.000286 | 0.001808 | 6.3× |
| 5 | 5 | 0.000310 | 0.003159 | 10.2× |
| 7 | — | — | — | (needs 128 GB) |

Improvement grows with N: 6.3× → 10.2×. Trend consistent with
prediction that larger systems have more modes to exploit and
V-shape becomes increasingly suboptimal.

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

1. **Test 2 redesign:** Bell initial state + Sum-MI metric
2. **Test 3 in C#:** Palindrome-timed relay vs fixed timing (N=11)
3. **N=7 scaling:** Run at home with 128 GB RAM
4. **Higher modes:** What do SVD modes 3-5 contribute?
5. **Combined optimization:** Mode 2 + frequency pulsing simultaneously

---

## References

- [Resonant Return (hypothesis)](../hypotheses/RESONANT_RETURN.md)
- [γ as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits baseline
- [γ Control](GAMMA_CONTROL.md): V-shape 21.5× baseline
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the eigenstructure
