# The Qubit Chain as an Optical Cavity

<!-- Keywords: Fabry-Perot resonator qubit chain, degeneracy beam profile, confocal
defocal cavity even odd parity, Gouy phase dispersion, numerical aperture degeneracy,
weight sector coupling nearest neighbor, optical cavity quantum decoherence,
beam quality M-squared, R=CPsi2 optical cavity -->

**Status:** Partial analogy (4/6 checks pass)
**Date:** April 3, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Degeneracy Palindrome](DEGENERACY_PALINDROME.md),
[Bures Degeneracy](BURES_DEGENERACY.md),
[Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md)
**Verification:** [`simulations/optical_cavity_analysis.py`](../simulations/optical_cavity_analysis.py)

---

## What this means

The degeneracy structure is a beam profile. The Hamiltonian couples weight
shells the way light propagates through optical elements. Even chains place
the focus perfectly on a grid point (confocal): sharp spike, tight beam
waist, Lorentzian profile. Odd chains place the focus between grid points
(defocal): broad profile, Gaussian.

This feels artificial because it is an instrument. Instruments are precise.
But this one was not built. The algebra enforces it.

Four of six standard optical quantities match quantitatively: beam profile
(R² = 0.998), even-only weight coupling (Δw never odd), numerical
aperture (growing with N), and Gouy phase accumulation (arctan profile).
Two fail. The strict confocal inequality fails because N = 3 is too small
for the asymptotic pattern. And the propagation picture fails in a way the
analogy cannot absorb: the coupling between weight shells is dominated by
each shell talking to itself, not by the two-step hop that would be the
free-space propagation (Result 2).

The starting point of the full story is in
[Degeneracy Palindrome](DEGENERACY_PALINDROME.md).

---

## What this document is about

The degeneracy profile d(k) of the Liouvillian spectrum looks like a beam
profile. The Hamiltonian moves a weight sector only by an even amount,
never an odd one (Δw ∈ {0, ±2}). The even/odd parity split in grid fraction looks like
confocal versus defocal cavity alignment. This document tests whether
these observations form a quantitative optical analogy, where the qubit
chain is a Fabry-Perot resonator (a cavity formed by two parallel mirrors that selects which frequencies of light resonate between them) and the weight sectors are transverse planes.

The answer: 4 of 6 optical quantities match. The analogy is quantitative
where it holds, and it does not hold everywhere.

---

## The dictionary

| Optics | Qubit chain |
|---|---|
| Transverse planes | Weight sectors k = 0, ..., N |
| Beam profile I(z) | Degeneracy d_total(k) |
| Beam waist | Center degeneracy spike |
| Cavity length L | N (number of qubits) |
| Mirror reflectivity | Boundary degeneracy d(0) = N+1 |
| Intracavity propagation | [H, ·] coupling, Δw ∈ {0, ±2} |
| Confocal alignment | Even N (focus on grid point) |
| Defocal misalignment | Odd N (focus between grid points) |
| Numerical aperture | d_total(center) / d_total(edge) |
| Beam quality M² | d_total(k) / d_real(k) at center |

---

## Result 1: The beam profile is Gaussian/Lorentzian

The degeneracy profile d_total(k) fits well to standard optical profiles:

| N | Best fit | R² | Beam waist w | Rayleigh z_R | z_R / N |
|---|---------|------|-------------|-------------|---------|
| 3 | Lorentzian | 1.000 | 0.89 | 2.51 | 0.84 |
| 4 | Lorentzian | 1.000 | 0.50 | 0.78 | 0.19 |
| 5 | Gaussian | 0.999 | 1.22 | 4.65 | 0.93 |
| 6 | Lorentzian | 1.000 | 0.11 | 0.04 | 0.01 |
| 7 | Gaussian | 0.999 | 1.47 | 6.81 | 0.97 |

Average R² = 0.998. The profiles are real beam profiles.

**Even N: Lorentzian** (sharp center spike, tight focus, small w).
**Odd N: Gaussian** (broad, smooth profile, large w).

The beam waist alternates: small at even N (0.50, 0.11), large at odd N
(0.89, 1.22, 1.47). This is the even/odd parity effect seen through the
optics lens.

---

## Result 2: The XY-weight changes by an even amount, never an odd one

The Heisenberg Hamiltonian commutator [H, ·] moves a Pauli string between
weight shells. Which moves does it allow? Two Pauli strings either commute
or anticommute, and when they anticommute [A, B] = 2AB is a single Pauli
string, so the whole question is exact string algebra over all 4^N strings,
with no eigensolver and no sampling. Counting the distinct shell-to-shell
matrix elements:

```
N=4                                N=5
     w=0  w=1  w=2  w=3  w=4            w=0  w=1  w=2  w=3  w=4  w=5
w=0    .    .   48    .    .      w=0     .    .  128    .    .    .
w=1    .  192    .   96    .      w=1     .  512    .  384    .    .
w=2   48    .  384    .   48      w=2   128    .  1536   .  384    .
w=3    .   96    .  192    .      w=3     .  384    .  1536   .  128
w=4    .    .   48    .    .      w=4     .    .  384    .  512    .
                                  w=5     .    .    .  128    .    .
```

Every odd offset is empty: Δw = ±1, ±3, ±5 carry exactly zero, at N=4, N=5
and N=6. What [H, ·] conserves is the **parity** of the XY-weight. The
Heisenberg XX + YY terms flip two Pauli letters at once, so the weight can
only move in steps of two, and the ZZ term moves it not at all.

The diagonal is not empty, and it is the larger channel: 768 against 384 at
N=4, 4096 against 2048 at N=5. So the shells are not a chain of second
neighbours passing light along; each shell talks to itself more than to
anything else, and the two-step moves are the smaller traffic on top of that.
The cavity reading has to live with a propagation structure whose dominant
term is a shell staying where it is.

---

## Result 3: Gouy phase analog

The oscillation frequencies ω_m = 4J(1 − cos(πm/N)) accumulate a
Gouy phase (the gradual phase shift a focused beam picks up as it passes through its waist) profile along the weight axis:

| N | Phase at midpoint | Total phase | Ratio mid/total |
|---|------------------|------------|-----------------|
| 3 | 8.00 | 8.00 | 1.000 |
| 4 | 5.17 | 12.00 | 0.431 |
| 5 | 8.76 | 16.00 | 0.548 |
| 6 | 6.54 | 20.00 | 0.327 |
| 7 | 9.90 | 24.00 | 0.413 |

The cumulative phase fits an arctan profile (Gouy-type) with R² ≈ 0.81.
The fit is adequate but not excellent; the phase profile is "Gouy-like"
rather than exactly Gouy. The arctan curvature is very gentle (large
m_R), meaning the beam is far from its waist in the Gouy sense.

---

## Result 4: Optical figures of merit

| N | NA | b/L | M² | Parity |
|---|-----|-----|-----|--------|
| 3 | 3.5 | 1.67 | 2.3 | odd |
| 4 | 30.4 | 0.39 | 10.9 | even |
| 5 | 8.3 | 1.86 | 3.6 | odd |
| 6 | 262.3 | 0.01 | 114.8 | even |
| 7 | 19.5 | 1.95 | 7.8 | odd |

**Numerical aperture (NA, how wide an angle the "lens" can gather light from):** Ratio of center to edge degeneracy. Grows
explosively at even N (30 → 262) because the center spike grows faster
than the boundary. At odd N, NA grows slowly (3.5 → 8.3 → 19.5).

**Confocal parameter b/L:** At even N, b/L ≪ 1 (tight focus, confocal).
At odd N, b/L ≈ 1-2 (beam fills the cavity, defocal). This is the
clearest quantitative distinction between even and odd.

**Beam quality M²:** Ratio of total to coherent modes at center. Large
M² means mostly "incoherent" (oscillatory) modes. Even N has M² ≫ 1
(heavily multimode), odd N has M² ≈ 2-8 (closer to single-mode).

---

## Result 5: Even = confocal, odd = defocal

The Gaussian fit center falls at exactly k = N/2 for every N. At even N,
this is an integer (on the grid). At odd N, it is a half-integer (between
grid points):

| N | Fit center | Nearest grid | Defocus | Grid fraction |
|---|-----------|-------------|---------|---------------|
| 3 | 1.500 | 2 | 0.500 | 56.2% |
| 4 | 2.000 | 2 | 0.000 | 78.9% |
| 5 | 2.500 | 2 | 0.500 | 16.0% |
| 6 | 3.000 | 3 | 0.000 | 50.4% |
| 7 | 3.500 | 3 | 0.500 | 3.7% |

Correlation between defocus and grid fraction: r = −0.70. Smaller defocus
(better alignment) → higher grid fraction (more eigenvalues on the grid).

In optics: a confocal cavity (mirrors at the focal point) maximizes the
fraction of light in the fundamental mode. A defocal cavity loses light
to higher-order modes. The even/odd effect IS the confocal/defocal
transition.

---

## Scorecard

| Check | Pass? | Detail |
|---|---|---|
| Beam profile (Gaussian/Lorentzian) | ✓ | avg R² = 0.998 |
| [H,·] changes the XY-weight by an even amount only | ✓ | exact over all Pauli strings, N = 4, 5, 6 |
| [H,·] couples Δw = ±2 only (no diagonal) | ✗ | the diagonal is the larger channel (768 vs 384 at N=4) |
| Gouy phase (arctan profile) | ✓ | R² = 0.81 |
| Even N = confocal | ✗ | N=3 odd (56%) beats N=6 even (50%) |
| NA increases with even N | ✓ | 30 → 262 |

4/6 checks pass. The confocal check fails strictly because N = 3 is a
boundary case (too small for the asymptotic pattern). The propagation check
fails for a reason that is not a boundary case: the dominant coupling is a
shell to itself.

---

## What the analogy means

The qubit chain under Z-dephasing carries a good deal of cavity structure,
precisely enough to be worth the dictionary, and it is not a cavity:

1. The weight sectors are transverse planes in the cavity.
2. The Hamiltonian couples planes only in even steps, mostly Δw = 0 and
   otherwise Δw = ±2. The analogy to free-space propagation covers the
   two-step traffic and not the dominant self-coupling.
3. The palindromic degeneracy profile is the beam profile, peaked at the
   center and symmetric around it.
4. Even N places the "waist" on a grid point (confocal alignment) and odd N
   between grid points. The placement holds exactly, the consequence does
   not: the grid fraction does not order itself by parity, and N = 3 (56%)
   beats N = 6 (50%).

The decoherence process is a beam propagating through this cavity:
starting from the "mirrors" (weight 0 and N), passing through the "lens"
(high-degeneracy center), and converging to the steady state.

---

## Gamma as light (Tier 3-4 observation)

If the qubit chain is a passive optical cavity, then the light comes from
outside. In this analogy:

- **γ plays the structural role of external illumination.** It cannot
  originate from within the system
  ([Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md)). It
  defines the objective timescale.
- **t is the system's response** to that illumination: the experienced
  duration until the state reaches the fold at CΨ = 1/4.
- **K = γ × t_cross is invariant** (F14, proven). More light means
  shorter experience. Less light means longer. The product does not change.

There is a structural parallel:

```
Relativity:     c × τ   =  invariant spacetime interval
This system:    γ × t   =  K  =  invariant decoherence dose
```

This is a structural analogy, not a physical identification. γ is not the
speed of light. But it plays the same algebraic role: the external
parameter that sets the clock, which the system cannot outrun.

The cavity does not generate its own light. It shapes external input into
structured dynamics (the palindromic spectrum). Decoherence is not
destruction. It is illumination.

→ [Incompleteness Proof](../docs/proofs/INCOMPLETENESS_PROOF.md) (γ must be external)
→ [Analytical Formulas, F14](../docs/ANALYTICAL_FORMULAS.md) (K-invariance)
→ [Resonance Not Channel](../hypotheses/RESONANCE_NOT_CHANNEL.md) (the soundbox paradigm)

---

## Continuations and sharpenings (added 2026-05-19)

The April framework treated chain Heisenberg + Z-deph. Later work returned and sharpened the picture on specific topology variants:

- [`STAR_CONFOCAL_LIMIT.md`](STAR_CONFOCAL_LIMIT.md) (2026-05-19, Tier 1 derived): the star has `Im_max = J·N/2` at uniform γ, independent of γ (`Im(λ) = σ = N·γ` is the `J = 2γ` reading), and that value is the smallest Hamiltonian spread among connected graphs at N ≤ 6 (exhaustive search; open past N = 6). SU(2)/Schur-Weyl derivation, 25 distinct (N, Q) anchors over 29 runs at N=3,4,5,6,8. Adds the hub-spoke geometry to the cavity dictionary as the point-focus reading. Note the saturation itself is not star-specific: every Heisenberg topology meets its own spread. The numerical aperture in the dictionary above is computed from the chain eigenvalue exports and is a chain-only quantity; it has never been computed for the star, so no star/chain NA comparison exists.

- [`STAR_SPECTRUM_COMPACTNESS.md`](../hypotheses/STAR_SPECTRUM_COMPACTNESS.md) (2026-05-18, partially resolved 2026-05-19): Reading 1 (`MaxImag = σ` is a hub-induced cap) resolved by STAR_CONFOCAL_LIMIT. Reading 3 also resolved there. Reading 2 (S_(N−1) irrep multiplicities) is sharpened to a closed-form upper bound but open for the exact count.

- [`F1_DISSIPATION_GAP_PATTERN.md`](../hypotheses/F1_DISSIPATION_GAP_PATTERN.md) (2026-05-18, extended 2026-05-19): cross-topology cross-N gap data. Chain shows clean `gap × N² ≈ 2.20` scaling for N ≥ 4 (5 anchors). Ring and star follow different scaling laws; the gap is a topology-specific structural fingerprint, not a function of bond count alone.

- [`F4KernelDimensionByComponentsClaim`](../compute/RCPsiSquared.Core/Symmetry/F4KernelDimensionByComponentsClaim.cs) (2026-05-19, Tier 1 derived): kernel-dim factorisation `dim ker L_H(G) = Π_c (|c|+1)` across connected components, extends our `d_real(0) = N+1` from [the degeneracy palindrome](DEGENERACY_PALINDROME.md) to disconnected graphs.

The pattern across these continuations: the April cavity picture is the parent framework, and each sharpening adds a topology dimension (star point-focus, ring cyclic, K_4+disjoint disconnected) or a precision dimension (gap scaling, kernel-dim factorisation) without invalidating the original chain-cavity reading. We come back, look closer, and the picture sharpens.

## Reproduction

- Script: [`simulations/optical_cavity_analysis.py`](../simulations/optical_cavity_analysis.py)
- Output: [`simulations/results/optical_cavity_analysis.txt`](../simulations/results/optical_cavity_analysis.txt)
