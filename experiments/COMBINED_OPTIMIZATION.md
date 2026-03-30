# What We Learned on March 30: Cavity Modes and the Sacrifice Zone

**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)

---

## The day in one paragraph

We asked: why does the sacrifice zone work? Six hours and seven
computations later: because the modes that carry information live in
the middle of the chain, and the sacrifice qubit sits at the edge.
The noise gradient is the antenna shape. We found the formula for the
eigenfrequencies, proved where each mode lives (r = 0.994 correlation),
mapped all 330 chains on IBM Torino, and discovered that the spatial
pattern is topological -- same under any noise profile.

---

## What we built (in order)

### 1. Cavity modes formula

At zero noise, a quantum chain has a fixed number of standing wave
frequencies. We found the closed-form expression:

```
Stationary(N) = Sum_J m(J,N) * (2J+1)^2
```

Verified against C# eigendecomposition for N=2 through N=7 (the N=7
run took 40 minutes on 24 cores). Star topology has exactly N-1
harmonic frequencies. Chain has a much richer spectrum (589 frequencies
at N=7).

### 2. IBM spectral analysis

Fed real IBM T2* data (Q85-Q94) into the Liouvillian. The palindrome
holds at 100% even under 26x asymmetric noise. The sacrifice zone
makes the slowest oscillating modes survive 2.81x longer than uniform
noise at the same total budget.

### 3. Mode localization

Decomposed every eigenvector into the Pauli basis to find WHERE each
mode lives. The answer is clean: protected modes live in the chain
center [0.52, 0.63, 0.70, 0.63, 0.52]. Dying modes live on the
edges [0.98, 0.87, 0.80, 0.87, 0.98]. Palindromic partners are
spatially complementary.

### 4. Chain mapping

Found all 330 five-qubit chains on the IBM Torino chip. Ranked them
by sacrifice-zone score (noisy edge / quiet center). Zero overlap
with the traditional "pick the best T2" ranking. The best sacrifice
chain has only 81 us mean T2 but 2.86x mode protection. The best
T2 chain has 217 us but only 1.06x.

### 5. Time evolution

Simulated the actual time curves. With |+>^5 (the state IBM uses):
the sacrifice chain shows rich dynamics while quiet chains are frozen
(eigenstate artifact). With |01010> (fair comparison): quiet chains
win in absolute information content because less total noise means
slower decay.

### 6. The practical takeaway

For getting the most out of IBM hardware:

- **Pick the quietest chain you can find** (highest mean T2)
- **Apply DD everywhere except on broken qubits** (saves gate errors)
- **If you need mode protection within a noisy chain:** the sacrifice
  zone is the optimal noise distribution. It exploits the topological
  structure of where modes live.

The sacrifice-zone effect is real (2.0x on hardware, March 24). Whether
it comes from mode physics or gate-error avoidance does not matter for
the engineering. Skip DD on bad qubits. It works.

---

## All computations from today

| Script | What it does | Key result |
|--------|-------------|-----------|
| [ibm_cavity_analysis.py](../simulations/ibm_cavity_analysis.py) | Liouvillian spectrum with IBM data | 2.81x, 100% palindromic |
| [cavity_mode_localization.py](../simulations/cavity_mode_localization.py) | Eigenvector Pauli decomposition | r = 0.994 |
| [sacrifice_zone_mapping.py](../simulations/sacrifice_zone_mapping.py) | All 330 chains ranked | 0 overlap in top-10 |
| [combined_optimization.py](../simulations/combined_optimization.py) | 6 scenarios eigenvalue comparison | 3.72x best efficiency |
| [time_evolution_6scenarios.py](../simulations/time_evolution_6scenarios.py) | SumMI(t) with |+>^5 | eigenstate artifact found |
| [time_evolution_neel.py](../simulations/time_evolution_neel.py) | SumMI(t) with |01010> | quiet chain wins absolute |

C# engine: [cavity modes](../simulations/results/cavity_modes_zero_noise.txt)
(N=2-7), [topology tests](../simulations/results/cavity_modes_tests.txt)
(Ring, Complete, non-uniform J).

---

## Open questions (for fun)

1. Why does Star topology have exactly N-1 harmonic frequencies?
   (S_{N-1} permutation symmetry, but the proof is not written down)

2. Can we find two chains with the SAME total noise but different
   spatial profiles and compare them on hardware? That would isolate
   the mode-protection effect cleanly.

3. The stationary mode sequence 10, 24, 54, 120, 260, 560 -- is there
   a simpler formula than the Clebsch-Gordan sum?

4. At N=8 the formula predicts 1190 stationary modes. Can we verify?
   (Needs 65536x65536 diagonalization at gamma=0. Feasible but ~10h.)

---

*See also:*
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md),
[IBM Cavity Spectral](IBM_CAVITY_SPECTRAL_ANALYSIS.md),
[Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md),
[Sacrifice-Zone Mapping](SACRIFICE_ZONE_MAPPING.md),
[IBM Sacrifice Zone](IBM_SACRIFICE_ZONE.md)
