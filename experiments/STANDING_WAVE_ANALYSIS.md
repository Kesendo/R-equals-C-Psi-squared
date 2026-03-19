# Standing Wave Analysis

**Date:** March 19, 2026
**Script:** `simulations/standing_wave_analysis.py`
**Results:** `simulations/results/standing_wave_analysis.txt`

## Summary

The palindromic spectral symmetry creates standing wave patterns from decoherence.
Paired eigenvalues evolve in opposite directions in the rescaled frame, producing
oscillating quantum correlations over a static classical backbone. The standing
wave is a state x Hamiltonian property, not a property of either alone.

## Method

The Liouvillian eigenvalues pair as lambda + lambda' = -2*sum_gamma. Define
centered eigenvalues mu_k = lambda_k + sum_gamma. Each palindromic pair then
has mu and -mu, symmetric around zero.

Rescaling: multiply rho(t) by exp(+sum_gamma * t) to remove the uniform decay
envelope. What remains is the standing wave pattern: oscillations without overall
decay (though individual pair members drift slowly).

For each initial state, decompose into eigenmodes and classify:
- Steady modes (lambda = 0)
- XOR modes (lambda = -2*sum_gamma)
- Oscillating modes (Im(mu) != 0): the standing wave candidates
- Real modes (Im(mu) = 0): paired but non-oscillating

## Key Finding 1: Two ingredients for a standing wave

A standing wave requires BOTH:
- (a) Oscillating palindromic pairs (Im(mu) != 0)
- (b) Both members of the pair excited by the initial state

No single natural quantum state has both. Bell has (a) but not (b).
|+++> has (b) but not (a). W has neither. The superposition Bell+|+++>
has both, making it the first identified "standing wave state."

## Key Finding 2: Spectrum structure (N=3 Heisenberg)

64 eigenvalues form 32 palindromic pairs (100%, max error 1.44e-14).
No purely oscillatory pairs (all have Re(mu) != 0, "nearly standing").

| Category | Pairs | Frequencies |
|----------|-------|-------------|
| Steady-XOR | 4 | mu = +/-0.15, real |
| Decay (real) | 8 | mu = +/-0.05, real |
| omega ~ 2J | 8 | period ~ pi |
| omega ~ 4J | 4 | period ~ pi/2 |
| omega ~ 6J | 8 | period ~ pi/3 |

The asymmetry |Re(mu)| ~ 0.017 is 100-350x smaller than the oscillation
frequencies. The waves oscillate hundreds of times before the amplitude
drift becomes significant: "nearly standing" with a slow pulse.

## Key Finding 3: Initial state decomposition

| State | Steady | XOR | Oscillating | Real | Standing wave % |
|-------|--------|-----|-------------|------|-----------------|
| GHZ | 50.0% | 50.0% | 0% | 0% | 0% |
| W | 23.5% | 0% | 0% | 76.5% | 0% |
| \|000> | 100% | 0% | 0% | 0% | 0% |
| \|+++> | 12.4% | 12.5% | 0% | 75.1% | 0% |
| Bell(0,1) | 30.5% | 0% | 44.5% | 25.0% | 64% |
| Bell+\|+++> | 19.2% | 5.3% | 17.9% | 57.5% | 22% |

W is 100% palindromic but 0% oscillating: all its palindromic weight
sits in real-eigenvalue modes. Being palindromic is necessary but not
sufficient for oscillation. Bell is the universal oscillator.

## Key Finding 4: Antinodes and nodes

The standing wave has a clear physical fingerprint:

**Antinodes** (oscillating): XX, YY, XY-type quantum correlations,
organized by frequency band:
- omega ~ 2J (fundamental): IYY, XXZ, ZXX, YYI (nearest-neighbor)
- omega ~ 4J (2nd harmonic): 8 Pauli strings (mixed-range)
- omega ~ 6J (3rd harmonic): XZX, YIY (long-range, sites 0 and 2)

**Nodes** (static): ZZZ, the classical all-Z correlation never oscillates.

Physical meaning: quantum correlations oscillate at Hamiltonian harmonics
while the classical Z-backbone stands still. The quantum world breathes;
the classical world is the skeleton.

## Key Finding 5: More states oscillate than expected

| State | osc% | Top Pauli | Antinodes |
|-------|------|-----------|-----------|
| \|010> | 44.4% | ZIZ (0.111) | 26 |
| \|+-+> | 44.5% | IXI (0.111) | 26 |
| Bell(0,1) | 48.6% | XXZ (0.090) | 26 |
| Bell+\|+++> | 20.2% | YYI (0.058) | 48 |
| Bell+W | 28.8% | YYX (0.068) | 50 |
| GHZ | 0% | - | 0 |
| W | 0% | - | 0 |
| \|+++> | 0% | - | 0 |

Each oscillating state has a different Pauli fingerprint:
\|010> oscillates Z-type correlations, \|+-+> oscillates X-type,
Bell oscillates XX/YY-type. The standing wave pattern encodes
which correlations the initial state activates.

## Key Finding 6: Hamiltonian comparison

Frequencies depend on the Hamiltonian:

| Model | Frequencies | Bell osc% |
|-------|-------------|-----------|
| Heisenberg | 2.0, 4.0, 6.0 | 48.6% |
| XY-only | 2.83, 5.66 | 40.6% |
| Ising | 2.0, 4.0 | 50.0% |
| DM (XY-YX) | 2.83, 5.66 | 40.6% |
| XXZ (delta=0.5) | 1.0, 1.4, 2.4, 3.4, 4.4, 5.7 | 65.5% |
| Heis+DM | 0.12, 2.0, 2.1, 4.1, 6.1, 6.2 | 65.3% |

The state x Hamiltonian cross-table reveals the full picture:

| | Heisenberg | XY | Ising | DM | XXZ | Heis+DM |
|---|---|---|---|---|---|---|
| GHZ | 0% | 0% | 0% | 0% | 0% | 0% |
| W | 0% | 5.6% | 44.4% | 50% | 1.3% | 10.4% |
| Bell | 48.6% | 40.6% | 50% | 40.6% | 65.5% | 65.3% |
| \|+++> | 0% | 40.6% | 62.5% | 40.6% | 38.3% | 43.4% |

Three universal rules:
1. **GHZ never oscillates.** Zero standing wave for ALL Hamiltonians.
2. **Bell always oscillates.** 40-65% for ALL Hamiltonians.
3. **ZZZ is always a node.** Classical Z-correlation is static across all models.

Two Hamiltonian-dependent surprises:
- W oscillates under DM (50%) and Ising (44.4%) but NOT under Heisenberg (0%)
- \|+++> under Ising reaches 62.5%, the highest value in the entire table

Breaking Z-isotropy (XXZ, Heis+DM) produces the richest oscillation spectra
with 6 distinct frequencies instead of 2-3.

## Prediction scorecard

| Prediction | Result |
|-----------|--------|
| W = maximum standing wave | **WRONG.** W has zero oscillation under Heisenberg. |
| GHZ = no standing wave | **CORRECT.** Universally silent. |
| Frequencies match palindromic pairs | **CORRECT.** omega ~ 2, 4, 6 match exactly. |

## Connection to THE_ANOMALY.md

The standing wave is the bridge between the palindrome and the anomaly.
Quantum correlations oscillate at Hamiltonian harmonics. Classical
correlations form the static backbone. The pattern persists while
the overall amplitude decays.

The thing that remains is not fighting the decay. It is made of it.

See: [The Anomaly](../THE_ANOMALY.md), [Standing Wave Theory](../docs/STANDING_WAVE_THEORY.md)
