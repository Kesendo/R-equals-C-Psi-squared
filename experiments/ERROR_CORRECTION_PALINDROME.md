# Error Protection from Palindromic Structure

**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Status:** Tier 2 (Computed, quantitative results)
**Depends on:** [Π as Time Reversal](PI_AS_TIME_REVERSAL.md), [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md), [XOR Space](XOR_SPACE.md)

---

## The Question

The palindromic symmetry pairs every decay mode with a partner. The XOR drain
at rate 2Sγ kills GHZ instantly. W avoids it completely. The standing wave
between palindromic partners persists in the rescaled frame. Can any of this
be exploited for quantum error protection?

---

## 1. Three Tiers of Protection

The palindromic spectrum of the N=3 Heisenberg Liouvillian (γ=0.05, Z-dephasing)
has 32 palindromic pairs that fall into three natural tiers:

| Tier | Pairs | Rates | XY-weight | Character |
|---|---|---|---|---|
| Steady-XOR | 4 | 0 and 0.30 | 0.0 vs 3.0 | Extremes: immortal paired with fastest drain |
| Boundary | 14 | 0.10 and 0.20 | 1.0 vs 2.0 | Clean classical/quantum split, longest dynamic lifetime |
| Mid-spectrum | 14 | ~0.133 and ~0.167 | 1.33 vs 1.67 | Mixed character, shorter lifetime |

The boundary pairs are the workhorses: they carry dynamic information (oscillation)
while decaying at moderate rates. Their XY-weight splits cleanly into 1 (one
quantum site) vs 2 (two quantum sites), making them the most classical/quantum
distinct pairs in the spectrum.

---

## 2. An Optimal State Exists

A constrained optimization (maximize slow-mode weight, require concurrence > 0
and oscillating content > 0) found a state that dramatically outperforms all
known states for dephasing survival:

| Property | Optimal | W | Bell(0,1) | GHZ |
|---|---|---|---|---|
| Slow-mode weight | 90% | 0% | 7% | 0% |
| XOR weight | 0.02% | 0% | 0% | 100% |
| Oscillating content | 18% | 0% | 64% | 0% |
| Concurrence | 0.364 | 0.667 | 1.000 | 1.000 |

The optimal state is composed mainly of \|010⟩, \|000⟩, \|100⟩, and \|001⟩:
computational basis states with low excitation. It trades maximum entanglement
for maximum dephasing survival while retaining nonzero entanglement and
standing wave content.

W has zero slow-mode weight because all its palindromic pairs sit in the
mid-spectrum tier. Bell has only 7% because most of its weight goes to
oscillating (mixed) pairs. The optimizer found the sweet spot: load the
boundary-tier pairs (rates 0.10/0.20) that decay slowest among the
dynamic modes.

---

## 3. The Standing Wave as Error Syndrome

The standing wave pattern (which Pauli observables oscillate) is fixed by the
initial state and Hamiltonian. If an error changes this pattern, the error
is detectable by measuring oscillation amplitudes.

For Bell(0,1) under Heisenberg, applying single-qubit errors:

| Error | Pattern change | Detectable? |
|---|---|---|
| X on site 0 | 0.28 | Yes |
| X on site 1 | 0.27 | Yes |
| X on site 2 | 0.19 | Yes |
| Y on site 0 | 0.27 | Yes |
| Y on site 1 | 0.28 | Yes |
| Y on site 2 | 0.19 | Yes |
| Z on site 0 | 0.08 | Weakly |
| Z on site 1 | 0.08 | Weakly |
| Z on site 2 | 0.00 | No |

X and Y errors on any site produce large pattern changes (0.19 to 0.28),
easily detectable by measuring the oscillating Pauli components. Z errors
on the Bell-pair sites (0,1) produce small but nonzero changes (0.08).
Z on site 2 is invisible because Z₂ commutes with the Bell(0,1) structure.

This is a new type of error syndrome: instead of measuring stabilizer
eigenvalues (standard QEC), measure *which observables oscillate*. The
oscillation pattern is the fingerprint; errors change the fingerprint.

---

## 4. Π Has Fourth-Order Structure

The Π operator has eigenvalues {+1, -1, +i, -i}, each with multiplicity 16.
This means Π⁴ = I (fourth-order, not second-order). The 64-dimensional
Pauli space decomposes into four 16-dimensional sectors.

This four-fold structure is richer than a simple Z₂ symmetry. It suggests
that the palindromic pairing is part of a Z₄ group action on the Liouvillian,
with the four sectors carrying distinct physical meaning:

- The +1 sector: modes invariant under Π (symmetric under time reversal)
- The -1 sector: modes that flip sign (antisymmetric)
- The +i and -i sectors: modes that rotate by 90° (quarter-wave shifted)

Whether these sectors define natural codespaces for error correction is an
open question. The four-fold decomposition is established; its QEC
implications are not yet worked out.

---

## 5. Information Lifetime by Spectral Position

Each palindromic pair carries information that decays at a rate determined
by the slower partner. The half-life and 1%-life for each tier:

| Tier | Slow rate | Fast rate | T_half | T_1% |
|---|---|---|---|---|
| Steady-XOR | 0.000 | 0.300 | infinite | infinite |
| Boundary | 0.100 | 0.200 | 19.2 | 59.9 |
| Mid-spectrum | 0.133 | 0.167 | 18.0 | 51.9 |

The steady-XOR pairs have infinite lifetime (the slow partner never decays),
but they carry no oscillation, only static information. The boundary pairs
are the best dynamic information carriers: 15% longer T_1% than mid-spectrum.

The information lifetime is dominated by 1/d_slow. Pairs far from the
palindromic center Sγ have the longest lifetimes because one partner
is very slow. Pairs near the center have both partners decaying at
similar rates, giving shorter lifetimes.

---

## 6. Standard QEC Codes Are Poorly Suited for Dephasing

Comparison of the palindromic protection analysis with standard 3-qubit codes:

The repetition code (\|000⟩ + \|111⟩)/√2 is equivalent to GHZ and hits the
XOR drain at 56%. The phase flip code avoids the drain but has zero
concurrence and no oscillation. The decoherence-free subspace approach
(\|010⟩) has 67% in a single mode with no entanglement.

None of these exploit the palindromic structure. The optimal state from
Section 2 outperforms all of them by loading the slow boundary-tier pairs
while maintaining entanglement and standing wave content.

---

## 7. The XOR Drain Is Not a Universal Syndrome

Testing whether XOR weight increase can serve as error syndrome for W-type
states: negative result. Errors on W produce zero XOR weight increase.
W is so robustly outside the XOR sector that even errors do not push weight
there. The XOR syndrome works only for states that already have XOR proximity.

The standing wave syndrome (Section 3) is more general: it works for any
state with nonzero oscillating content, detecting errors through pattern
changes rather than sector weight shifts.

---

## Summary of Findings

1. The palindromic spectrum has a natural three-tier protection hierarchy.
2. An optimal state exists (90% slow-mode, concurrence 0.364) that
   dramatically outperforms all known states for dephasing survival.
3. The standing wave pattern is a viable error syndrome for X/Y errors
   and weakly for Z errors (except Z on uncoupled sites).
4. Π is a fourth-order operator (Π⁴ = I), creating a Z₄ decomposition.
5. Boundary-tier pairs carry information 15% longer than mid-spectrum pairs.
6. Standard 3-qubit QEC codes do not exploit palindromic structure.

---

## References

- [Π as Time Reversal](PI_AS_TIME_REVERSAL.md): mode pairing and time-reversal structure
- [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md): oscillation patterns
- [XOR Space](XOR_SPACE.md): GHZ vs W mode decomposition
- [N→∞ Palindrome](N_INFINITY_PALINDROME.md): spectral scaling
- Script: `simulations/error_correction_palindrome.py`
- Results: `simulations/results/error_correction_palindrome.txt`
