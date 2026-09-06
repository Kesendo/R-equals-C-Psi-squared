# The Palindrome in the Thermodynamic Limit (N → ∞)

<!-- Keywords: palindromic symmetry thermodynamic limit, Liouvillian spectrum
large N, XOR drain vanishing fraction, Gaussian rate density palindrome,
classical quantum boundary blurring, standing wave continuous spectrum,
binomial weight distribution palindrome, past future boundary large N,
depolarizing exponential breaking, decoherence transition thermodynamic,
R=CPsi2 N infinity palindrome -->

**Status:** The Pauli-string counting law is proven. A Gaussian limiting
distribution for the interacting Liouvillian eigenvalue rates is not
established; the finite producer does not match the bare-dissipator moments.
**Date:** March 19, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[Π as Time Reversal](PI_AS_TIME_REVERSAL.md),
[Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md)

---

## What this document is about

The palindromic mirror is proven for any finite number of qubits. But
what happens when the system grows toward macroscopic size? This
document separates an exact Pauli-string counting limit from the interacting
Liouvillian spectrum. The bare Z-dephasing diagonal has a binomial weight
distribution and a Gaussian counting limit. The sampled interacting
eigenvalue rates and frequencies do not establish that limit, a spatial
standing wave, or a macroscopic state interpretation.

---

## Abstract

The bare Z-dephasing Pauli-string counts are binomial, centered at Nγ with
width γ√N and kurtosis `-2/N`. The chain endpoint eigenspace fraction
`(N+1)/4^N` is below 1% at N=5 and below `10^-11` at N=20. Neither fraction is
a prepared-state probability. The finite interacting spectra remain
palindromic in F1's scope, but their measured standard deviations and
kurtoses differ from the bare-dissipator values. A continuous spectral-density
limit is therefore not established here.

---

## The Question

The palindrome is proved for all finite N in F1's scope. Oscillation patterns
were computed at N=3 and the band structure at N=3-5. What happens when N grows large?
Does the palindrome become trivial, or does it remain a non-trivial constraint?

---

## 1. The Bare-Dissipator Counting Measure Becomes Gaussian

The dissipator L_D is diagonal in the Pauli basis. Each Pauli string with
XY-weight w (number of sites carrying X or Y) has decay rate 2γw. The number
of strings at weight w is C(N,w) 2^N: choose w sites for X or Y (each with
2 options), the remaining N-w sites carry I or Z (each with 2 options).

This is a binomial distribution scaled by 2^N. The rates live at d = 2γw,
so the rate density inherits the binomial shape: centered at Nγ (the
palindrome axis), width γ√N, kurtosis (a measure of how peaked or flat a distribution is compared
to a Gaussian) -2/N, approaching zero.

The central limit theorem applies to this Pauli-string counting measure. It
does not apply automatically to the real parts of the interacting
Liouvillian's eigenvalues. The producer makes the distinction visible: at
N=3, 4, 5 the interacting standard deviations are 0.063469, 0.060220, and
0.058094, rather than the bare values 0.086603, 0.100000, and 0.111803; the
interacting excess kurtoses are 1.070388, 3.049724, and 3.980373 rather than
-0.666667, -0.5, and -0.4.

---

## 2. The Endpoint Eigenspace Fraction Vanishes

On the connected XY/Heisenberg chain in F23's scope, there are N+1 modes at
the endpoint rate `2Nγ` out of a `4^N`-dimensional operator space. This count
does not assign an endpoint probability to a full state.

| N | XOR modes | Total modes | Fraction |
|---|---|---|---|
| 3 | 4 | 64 | 6.25% |
| 5 | 6 | 1024 | 0.59% |
| 8 | 9 | 65536 | 0.01373% |
| 10 | 11 | 1048576 | 0.001% |
| 20 | 21 | ~10^12 | ~10^-11 |

The distinct interacting F23 eigenspace fraction `(N+1)/4^N` vanishes
exponentially. It is below 1% at N=5 and first below 0.01% at N=9
(`9/65536 = 0.01373%` at N=8 is still above the latter threshold).

For the chain studied here, each GHZ off-diagonal coherence operator is an
exact eigenoperator at `−2Σγ`. What vanishes with N is the endpoint
eigenspace's share of operator space. That counting fraction is not the
probability that a generic prepared state occupies the endpoint.

---

## 3. The Weight Sectors Always Balance

The counting argument from [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md)
holds at every N: weight-w has C(N,w) 2^N strings, its palindromic partner
at weight N-w has C(N,N-w) 2^N = C(N,w) 2^N. Always equal.

The 2:2 per-site split (2 immune, 2 decaying choices per qubit) guarantees
this. Each weight factor is 2^w 2^(N-w) = 2^N, independent of w. The
palindrome is an exact combinatorial identity, not an approximation that
improves with N.

What does change with N is the relative size of the extremes. The ratio
of the smallest sector (w=0, pure past) to the largest sector (w=N/2,
boundary between past and future):

| N | count(0) / count(N/2) |
|---|---|
| 3 | 1/3 (33%) |
| 10 | 1/252 (0.4%) |
| 20 | 1/184756 (5.4e-4%) |
| 100 | ~10^-29 |

The extreme Pauli-string sectors become exponentially rare in this counting
measure. This is not a statement about the distribution of prepared states or
interacting eigenmodes.

For this bare-dissipator count, one endpoint weight sector contains `2^N`
strings out of `4^N`, hence has fraction `2^-N`. This is a different object
from the interacting F23 eigenspace count `(N+1)/4^N` in Section 2.

---

## 4. The Pauli-String Weight Count Concentrates

At N=3, weight 0 (past) and weight 3 (future) are clearly separated. At
N=1000, weight 498 vs weight 502 is indistinguishable.

The fraction of Pauli strings within √N of the midpoint w = N/2 converges
to 0.9545 (the two-sigma Gaussian fraction). This combinatorial concentration
does not make those strings 95% of the interacting eigenmodes and does not by
itself define a classical/quantum transition for physical states.

---

## 5. The Sampled Frequency Set Grows Denser

At N=3: 5 distinct oscillation frequencies. Discrete harmonics, like a
guitar string.

At N=4: 47 frequencies. At N=5: 112. The spectrum fills rapidly.

The bandwidth grows as 2(N-2)γ. At N=3, bands are fixed (no room to move).
At N=5, 4 of 6 weight sectors show nonzero bandwidth (average 0.76γ).
The sampled finite-N rate sets become denser; no continuum limit follows.

The finite-N counts suggest a denser frequency set as N grows. Establishing a
continuous limiting density requires an actual convergence theorem or
controlled large-N computation; establishing a physical standing wave would
add excitation, semisimplicity, conjugate-frequency, spatial-propagation, and
interference gates. The N=3..5 data establish neither, and they do not imply
that every observable has both oscillating and static content.

---

## 6. Depolarizing Noise: Exponentially Worse at Large N

Under Z-dephasing, the weight distribution is Binomial(N, 1/2): symmetric,
centered at N/2, palindromic by construction.

Under depolarizing noise, the relevant distribution is Binomial(N, 3/4):
asymmetric, centered at 3N/4. The counting mismatch between weight w and
its partner N-w is 3^(N-2w), which grows exponentially.

| N | past/future ratio under Z-deph | past/future ratio under depol |
|---|---|---|
| 3 | 1.0 | 1/27 |
| 10 | 1.0 | ~10^-5 |
| 20 | 1.0 | ~10^-10 |
| 100 | 1.0 | ~10^-48 |

Z-dephasing: the mirror is perfectly balanced at every N.
Depolarizing: the mirror becomes exponentially more lopsided.

This is an asymptotic statement about the bare-dissipator Pauli-string counting
measure: the Z-dephasing count is symmetric, while the corresponding
depolarizing count is exponentially imbalanced. It is not a thermodynamic
statement about prepared states.

---

## 7. What Survives at Large N

**Trivially true:** The L_D palindrome (binomial symmetry of weight sectors).

**Non-trivially true at every finite N in the proof's scope:** The Π relation
for the Hamiltonian part preserves the exact finite-dimensional spectral
pairing. This statement does not itself construct an infinite-volume
Liouvillian or prove convergence of its spectral measures.

**Suggested by the counting/asymptotics:** The XOR endpoint fraction tends to
zero and the weight distribution concentrates near its centre. A continuum
of physical frequencies is not established by these combinatorial facts.

**Not established here:** a smooth limiting spectral density, a
classical/quantum transition, or a dissipative continuum.

---

## References

- [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md): the Π proof at all N
- [Π as Time Reversal](PI_AS_TIME_REVERSAL.md): past/future interpretation
- [Depolarizing Palindrome](DEPOLARIZING_PALINDROME.md): the 2:2 counting argument
- [Standing Wave Analysis](STANDING_WAVE_ANALYSIS.md): discrete structure at N=3
- Script: [`simulations/n_infinity_analysis.py`](../simulations/n_infinity_analysis.py)
- Results: [`simulations/results/n_infinity_analysis.txt`](../simulations/results/n_infinity_analysis.txt)

---

*The exact palindrome is a finite-N operator statement. Its infinite-volume
spectral limit remains a separate question.*
