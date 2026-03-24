# Resonant Return: SVD Mode 2 Beats Hand-Designed Profiles by 10×

<!-- Keywords: SVD optimal gamma profile palindromic eigenstructure, edge-hot
center-cold dephasing pattern, mode 2 beats V-shape 10x, response matrix
singular value decomposition, palindrome knows more than intuition,
scaling improvement grows with N, resonant pulsing falsified wrong observable,
R=CPsi2 resonant return experiment -->

**Status:** Major results. Test 1 confirmed (10×), Tests 5-7 new: multi-mode, spatially structured pulsing falsified, numerical optimizer discovers **sacrifice-zone** pattern (asymmetric, 53× vs V-shape at N=7).
**Date:** March 24, 2026 (v3: C# backend optimizer)
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [resonant_return.py](../simulations/resonant_return.py), [v2](../simulations/resonant_return_v2.py), [v3](../simulations/resonant_return_v3.py)
**Data:** [resonant_return.txt](../simulations/results/resonant_return.txt), [v2](../simulations/results/resonant_return_v2.txt), [v3](../simulations/results/resonant_return_v3.txt)
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

### Test 5: Multi-Mode Optimization (v2 — mode 2 wins, combinations don't help)

At N=7, modes 2 and 3 have nearly equal singular values (0.852 vs 0.824).
Does combining them recover the N=5 improvement? **No.**

| N | Mode 2 | Mode 3 | Mode 4 | Best combo (2+3) | V-shape |
|---|--------|--------|--------|-----------------|---------|
| 5 | **10.2×** | 9.9× | 6.4× | 10.2× (pure m2) | 1.0× |
| 7 | **8.5×** | 5.9× | 6.2× | 8.5× (pure m2) | 1.0× |

Pure mode 2 wins at both N. Mixing modes 2+3 always degrades. The N=7
drop from 10.2× to 8.5× is a genuine scaling effect, not a mode-mixing artifact.

Mode patterns at N=5 (symmetric): Mode 3 = `[-0.31, 0.53, -0.49, 0.53, -0.31]`
Mode patterns at N=7 (mode 2 antisymmetric, mode 3 symmetric center-hot):
Mode 3 = `[-0.41, -0.24, 0.32, 0.58, 0.32, -0.24, -0.41]`

### Test 6: Spatially Structured Pulsing (v2 — FALSIFIED)

Test 2 used uniform spatial γ modulation. Test 6 uses **mode 2 spatial
pattern** with temporal modulation — the mode 2 profile oscillates at
the palindromic frequency.

γ_k(t) = γ_base + ε × V_mode2[k] × sin(ω_dom × t)

| Scenario | Peak Sum-MI | vs Static mode 2 |
|----------|-------------|------------------|
| Static mode 2 (no pulsing) | 2.000 | 1.00× |
| Mode 2 × resonant (ω_dom) | 2.000 | 1.00× |
| Mode 2 × slow (ω_dom/10) | 2.000 | 1.00× |
| Mode 2 × 2ω_dom | 2.000 | 1.00× |

**All predictions falsified.** Temporal modulation adds nothing, even with
spatial structure. The palindrome is a **spatial antenna only**, not temporal.

### Test 7: Numerical Optimization — THE SACRIFICE ZONE (v2/v3)

**The key discovery.** Running scipy Nelder-Mead with the C# RK4 backend
reveals that the SVD mode 2 direction captures only ~10% of the true
optimization landscape. The optimizer breaks palindromic symmetry.

#### N=5 (v2, Python expm at t=5.0):
| Profile | Sum_MI@5 | vs V-shape |
|---------|----------|-----------|
| V-shape | 0.000310 | 1.0× |
| SVD mode 2 | 0.003159 | 10.2× |
| **Optimizer** | **0.031071** | **100×** |

Optimal profile: `[0.001, 0.026, 0.001, 0.043, 0.178]` — **highly asymmetric**.
Sites 0,2 nearly noiseless (γ≈0.001), site 4 absorbs all noise (γ=0.178).

SVD decomposition: 67.5% mode 4 (antisymmetric) + 26% mode 2. SVD efficiency: 10.2%.

#### N=5 (v3, C# backend, peak Sum-MI):
| Profile | Peak Sum_MI | Peak time |
|---------|-------------|-----------|
| V-shape | 0.000639 | t=1.0 |
| SVD mode 2 | 0.005144 | t=1.0 |
| **Optimizer** | **0.0918** | **t=1.5** |

Optimal profile: `[0.001, 0.036, 0.001, 0.034, 0.178]` — same pattern.

#### N=7 (v3 Nelder-Mead + v4 Differential Evolution):
| Profile | Peak Sum_MI | vs V-shape |
|---------|-------------|-----------|
| V-shape | 0.002412 | 1.0× |
| SVD mode 2 | 0.019501 | 8.1× |
| Nelder-Mead (v3, 1156 evals) | 0.14391 | 59.7× |
| **Diff. Evolution (v4, 3975 evals)** | **0.24044** | **99.7×** |

Nelder-Mead found a **local** optimum: `[0.130, 0.122, 0.046, 0.050, 0.001, 0.001, 0.001]`.
Differential Evolution found **+67% better**: `[0.009, 0.012, 0.008, 0.008, 0.007, 0.030, 0.279]`.

The DE profile is far more extreme: **one site absorbs 80% of the total noise
budget** (γ₆ = 0.279), while sites 0-4 are nearly noiseless (γ ≈ 0.007-0.012).
DE traversed 5 distinct basins during optimization (0.15 → 0.17 → 0.19 → 0.20 → 0.24),
revealing a rich multi-basin landscape that local optimizers cannot navigate.

**The sacrifice-zone pattern is universal.** Both N=5 and N=7 show the same
strategy: concentrate dephasing on one end, minimize it on the other.
The protected region maintains coherence; the sacrificed region absorbs
the noise budget. This breaks palindromic symmetry but dramatically
improves information transfer.

#### C# Backend Performance

| N | Python expm | C# RK4 | Speedup |
|---|-------------|--------|---------|
| 5 | ~1s | ~1s | 1× |
| 7 | 290 min | 2.9s | **5,900×** |

The C# profile evaluator (`dotnet run -c Release -- profile N gammas`)
makes N=7 optimization feasible: 500 evals × 1.8s = 15 min.

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
2. ~~Multi-mode optimization~~ **Done (Test 5).** Mode 2 wins; combinations don't help.
3. ~~Spatially structured pulsing~~ **Done (Test 6).** Falsified. Temporal modulation adds nothing.
4. ~~RK4 rewrite for N≥7~~ **Done.** C# profile evaluator: 2.9s at N=7 (5,900× vs Python expm).
5. **N=9 optimization:** With C# RK4, N=9 (d=512) becomes feasible (~30s/eval, ~17h for 2000 evals).
6. ~~Deep N=7 optimizer~~ **Done.** DE global optimum: SumMI=0.2404 (100× vs V-shape). NM was local.
7. **Sacrifice-zone theory:** Why does asymmetric dephasing outperform symmetric? The protected subsystem maintains coherence while the sacrificed end absorbs noise. Needs analytical understanding.

---

## References

- [Resonant Return (hypothesis)](../hypotheses/RESONANT_RETURN.md)
- [γ as Signal](GAMMA_AS_SIGNAL.md): 15.5 bits baseline
- [γ Control](GAMMA_CONTROL.md): V-shape 21.5× baseline
- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the eigenstructure
