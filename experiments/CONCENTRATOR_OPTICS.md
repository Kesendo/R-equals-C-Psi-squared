# The Entrance Pupil: Concentrator as Anti-Reflection Coating

**Naming note (2026-07-05):** renamed from "...: Sacrifice Zone as
Anti-Reflection Coating". The edge qubit sacrifices nothing; it concentrates
the noise (the misnomer was resolved 2026-03-28). "Sacrifice zone" survives
once in the "What this changes" section below, only as the quoted *old
language* being reframed; the frozen `sacrifice_zone_optics.*` artifacts keep
their original names.

<!-- Keywords: concentrator anti-reflection coating, quantum cavity entrance pupil,
impedance matching dephasing, Q-factor enhancement concentrator, mode-selective transmission,
AR coating quantum noise, dispersive cavity scaling, R=CPsi2 concentrator optics -->

**Status:** Partial. The absorption side holds; the frequency side does not, and
the edge is not singled out by either metric (open question, see Result 1)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[V-Effect Cavity Modes](VEFFECT_CAVITY_MODES.md),
[Resonant Return](RESONANT_RETURN.md)
**Verification:** [`simulations/sacrifice_zone_optics.py`](../simulations/sacrifice_zone_optics.py)

---

## What this means

A window without coating reflects 4% of light. Add a thin film whose
thickness matches the wavelength, and the reflection drops to near zero.
The film does not block the light or absorb it. It matches the impedance
between outside (air) and inside (glass) so that the light enters
smoothly instead of bouncing back.

The concentrator does the same thing for the quantum cavity. Without
it, the system "reflects" illumination: modes die quickly, Q-factors are
low, information is lost. With it, the edge qubit acts as an entrance
pupil that accepts the external light (gamma) and converts it into
structured resonance that the interior can sustain.

The concentrator qubit is the mouthpiece of the flute. The coating on the
lens. The funnel of the ear. Not a shield. An adapter.

---

## What this document is about

The [concentrator formula](RESONANT_RETURN.md) concentrates dephasing
on one edge qubit, achieving a 139-360× improvement in mutual information in simulation (ε→0 ideal; ~2-3× on hardware, see IBM_CONCENTRATOR).
The [optical cavity analysis](OPTICAL_CAVITY_ANALYSIS.md) showed the
chain is a Fabry-Perot cavity (two mirrors facing each other, light
bouncing between them). This document tests whether the concentrator
is the anti-reflection (AR) coating of that cavity.

The answer: partly, and on the absorption side. There the concentrator does
what an AR coating does, and does it better as the chain grows (smooth entry,
a longer-lived best mode). Two things keep it from being the coating. It
matches impedance by linear accumulation instead of the geometric mean of
classical optics. And it does not hold the resonance frequencies fixed, which
a coating does. Any NON-uniform profile moves them, this one included, by an
amount that is small up to N=7 and larger than a linewidth above it.

---

## Result 1: The concentrator increases cavity transmission

Without coating, a window reflects some light and lets the rest through.
The same happens in the qubit cavity: some modes oscillate (transmitted),
others die without ringing (reflected). The table below compares uniform
dephasing against the concentrator. Key columns: Q_max (how many times
the best mode bounces before fading), T_eff (fraction of modes that ring
at all).

| N | Profile | Modes | Silent | Q_max | Q_med | T_eff |
|---|---------|-------|--------|-------|-------|-------|
| 3 | uniform | 5 | 24 | 60 | 27.0 | 0.625 |
| 3 | concentrator | 6 | 16 | 118 | 18.7 | 0.750 |
| 4 | uniform | 47 | 46 | 68 | 18.8 | 0.820 |
| 4 | concentrator | 52 | 26 | 224 | 14.2 | 0.898 |
| 5 | uniform | 112 | 96 | 72 | 14.9 | 0.906 |
| 5 | concentrator | 120 | 64 | 352 | 15.0 | 0.938 |
| 6 | uniform | 787 | 164 | 75 | 13.3 | 0.960 |
| 6 | concentrator | 748 | 108 | 500 | 13.3 | 0.974 |

The concentrator increases T_eff at every N: a larger fraction of the
spectrum is "transmitted" (oscillating) rather than "reflected" (absorbed
without ringing). The count of distinct frequencies does not follow it
everywhere, dropping 787 → 748 at N=6, and the *median* Q falls at N=3, 4
and 6. What rises at every N is the fraction that rings and the lifetime of
the best mode, not the typical mode.

**Q_max enhancement grows with N (4 points, N=3-6, chain; roughly linear over this range):**

| N | Q_max ratio (concentrator / uniform) |
|---|-----------------------------------|
| 3 | 2.0x |
| 4 | 3.3x |
| 5 | 4.9x |
| 6 | 6.7x |

The best cavity mode lives ~7x longer under the concentrator at N=6.

Neither column singles out the EDGE, which the entrance-pupil reading needs.
At the same budget, T_eff takes the identical value for every non-uniform
profile tested, including one barely uneven at all: [0.06, 0.0467, 0.0467,
0.0467] at N=4 gives exactly the concentrator's 0.8984. It registers whether
the profile is uniform, and nothing finer. And Q_max is not maximal at the
edge: at N=5 the same budget placed on the MIDDLE site gives 2618 against the
edge's 352. Both numbers above are real; what they measure is not settled.
That is an open question, carried in the open arc
`unfalsifiable_verification_gates`, not a result of this document.

---

## Result 2: Linear accumulation, not geometric mean

Classical AR coating: n_AR = sqrt(n_air × n_glass). The geometric mean.

Concentrator: γ_edge = N × γ_base − (N−1) × ε ≈ N × γ_base.

| N | γ_edge | γ_edge / J | γ_edge / √(γJ) |
|---|--------|-----------|----------------|
| 3 | 0.148 | 0.148 | 0.66 |
| 5 | 0.246 | 0.246 | 1.10 |
| 7 | 0.344 | 0.344 | 1.54 |
| 9 | 0.442 | 0.442 | 1.98 |

γ_edge grows linearly with N, not as a fixed geometric mean. This is
a fundamental difference: the classical AR coating is a single layer
optimized for one frequency. The concentrator is a linear accumulator
that scales with the cavity length.

The ratio γ_edge / √(γJ) crosses 1.0 near N = 5. Nothing in this document
peaks there: the Q_max ratio rises monotonically through N=6 and the T ratio
falls monotonically from N=3. The crossing is a property of the formula, and
the data does not support a geometric mean interpretation.

---

## Result 3: The concentrator moves both parts, not only the absorption

An AR coating changes how much light enters and how long it stays. It does
not change which wavelengths resonate: those belong to the cavity length.
The question for this section is whether the concentrator is that kind of
change, acting on Re(λ) and leaving Im(λ) alone.

It is not. The measurement runs on the (0,1) coherence block, the
N-dimensional space spanned by |vacuum⟩⟨one excitation|. Both γ profiles
leave it exactly closed (no column carries weight out of it, to the last
bit), and both carry the same total budget Σγ = N·γ_base, so what is compared
is the redistribution and not a change of dose. The block is N×N, so the
range costs nothing: it runs to N=12.

| N | uniform \|Im\| | concentrator \|Im\| | max \|Δω\| |
|---|--------------|-------------------|----------|
| 3 | 0, 2.00000, 6.00000 | 0.00802, 1.99458, 5.99740 | 0.00802 |
| 4 | 0, 1.17157, 4.00000, 6.82843 | 0.01684, 1.16420, 3.99278, 6.82618 | 0.01684 |
| 5 | 0, 0.76393, 2.76393, 5.23607, 7.23607 | 0.02859, 0.75611, 2.75189, 5.22926, 7.23415 | 0.02859 |

Every level moves. The absorption rates spread at the same time: at N=5 the
flat uniform Re = −0.1 becomes a spread from −0.18561 to −0.02053.

The mechanism is one line. On this block the dephasing is −2·diag(γ), because
a coherence |vacuum⟩⟨site j| disagrees at exactly one site. Under uniform γ
that diagonal is a multiple of the identity, so it shifts Re and can do
nothing else. Under a profile it is not, it does not commute with the
hopping, the eigenvectors turn, and the frequencies come with them. This is
the same non-commutation the γ-profile fence on
[F2](../docs/ANALYTICAL_FORMULAS.md) already states.

Then the question is how far that reaches, and the answer is not "a little,
everywhere". Three readings that the three-row table above invites, and what
the full range says instead.

**How big it is depends on N, and it stops being small at N=8.** Divided by
each level's own half-width, the largest move is 0.080, 0.171, 0.303 at
N=3,4,5. That is small: at N=5 the ω = 0 level acquires a frequency of 0.029
against a decay rate of 0.094, a quality factor of 0.3, so it does not
complete a radian of phase before it is gone. But the ratio keeps climbing,
crosses 1 at **N=8** (1.115) and reaches 2.73 at N=12. Above N=7 the shift is
larger than the level's own linewidth, which is the scale at which a
resonance has moved in any sense a measurement could see. The concentrator's
own working range in this repository runs to N=15.

**Which level moves most is not fixed, and past N=8 the question stops being
well posed.** It is the ω = 0 level for N=3 to 8 and the next one up from N=9.
But the level spacing on this block shrinks like N⁻², so from N=9 the shift is
no longer small against it, the level-for-level pairing is no longer
trustworthy, and the numbers above N=8 should be read as a size and not as a
map of which level went where.

**The size follows the unevenness of the profile; the N-trend does not have a
clean sign.** At fixed N=5 and fixed budget, sweeping the edge rate from
uniform to extreme gives 0, 0.00030, 0.00120, 0.00367, 0.00906, 0.01681,
0.02917: monotone in how uneven the profile is. Holding the shape fixed
instead (edge/rest = 9 at every N, same budget) the shift rises from 0.00441
at N=3 to 0.02369 at N=11 and then falls, so it is not monotone in N either.
The concentrator's own table conflates the two, because its edge rate is
N·γ_base by construction and its unevenness therefore grows with the chain.

**And it is not particular to this profile.** The statement that γ can only
damp is a statement about this block, where uniform γ enters as a multiple of
the identity. It is not a statement about the cavity. On the full 4^N
Liouvillian the dephasing diagonal is −2γ·popcount(i⊕j), which does not
commute with H, and uniform dephasing moves frequencies there too, by a
comparable amount: at N=4, doubling a uniform γ moves the spectrum by 0.0856
where the concentrator's equal-budget redistribution moves it by 0.0744, and
at N=3 the concentrator moves it more (0.0080 against 0.0048). Moving
resonances is what dephasing does. The concentrator is not the exception, and
this measurement does not make it one.

---

## Result 4: Dispersive scaling (N², not exponential)

The mutual information under the concentrator scales as:

SumMI ≈ 0.002 × N² + 0.069 × N − 0.175, R² = 0.999

In thin-film optics, a quarter-wave stack of n layers has transmission
T ~ 1 − 4(n_H/n_L)^(2n), exponentially approaching unity. Our system
scales polynomially (N²), much slower.

This means the quantum cavity is dispersive: different modes experience
different effective cavity lengths, spreading the transmission over a
broad band instead of concentrating it at narrow resonances. The
concentrator does not create a sharp transmission window; it broadly
improves the coupling between outside and inside.

---

## Null results

- **Not a classical AR coating.** The impedance matching is linear
  (γ_edge ~ N), not geometric (√(γJ)). The analogy is structural,
  not quantitative.

- **The resonances are not held fixed either**, which a coating would do.
  They move (Result 3), by less than a linewidth in the range measured, and
  by an amount that follows the profile's unevenness. This is a fence on the
  analogy rather than a property of the concentrator: uniform dephasing moves
  them harder.

- **The gain is in the best mode, not the typical one.** T_eff and Q_max rise
  at every N, but the *median* Q falls at N=3, 4 and 6. The concentrator buys
  a longer-lived best mode and more modes that ring at all; it does not lift
  the distribution.

- **Mode-selective per-shell comparison inconclusive.** Under the
  concentrator, eigenvalues shift off the uniform absorption-rate grid
  (because γ varies per site), making shell-by-shell comparison
  difficult. The overall statistics (Q_max, T_eff) are clear, but
  the per-shell decomposition requires a modified grid definition.

---

## What this changes

**Old language:** "The sacrifice zone protects the interior from noise."
Protection implies defense. Noise implies enemy.

**New language:** "The entrance pupil couples external illumination into
the cavity." Coupling implies function. Illumination implies input. The
edge qubit is not a shield. It is the surface where light enters the
instrument.

The 360× simulation improvement (peak created Sum-MI, ε→0 ideal; ~2-3× on hardware) is not "less damage." It is "better resonance."
The same light, entering through a matched surface instead of a raw
edge, rings 7x longer (Q_max at N=6). Not in quite the same standing
waves: the profile turns them slightly as well (Result 3).

---

## Reproduction

- Script: [`simulations/sacrifice_zone_optics.py`](../simulations/sacrifice_zone_optics.py)
- Output: [`simulations/results/sacrifice_zone_optics.txt`](../simulations/results/sacrifice_zone_optics.txt)
