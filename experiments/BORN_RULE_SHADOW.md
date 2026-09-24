# The Born Rule Is a Shadow, Not a Photograph

<!-- Keywords: Born rule interference pattern, measurement as photography, CΨ fold
shutter, past future mode decomposition, purity interference cross term, geometric
optics Born rule, exposure time quantum measurement, R=CPsi2 Born rule -->

**Status:** Computed (the Born probabilities carry no past/future cross term, by linearity of the trace; the purity carries one only where slow and fast modes overlap, which at N = 2 they never do, and at N = 3 it is a fraction of a percent)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Standing Waves](FACTOR_TWO_STANDING_WAVES.md),
[Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[Born Rule Mirror](BORN_RULE_MIRROR.md)
**Verification:** [`simulations/born_rule_shadow.py`](../simulations/born_rule_shadow.py)

---

## What this means

A photograph records an interference pattern. A shadow records where
the light was blocked. They look similar on the screen, but they are
made by different physics. The interference pattern depends on phase;
the shadow does not.

The Born rule gives the probability of each measurement outcome (the
formula P(i) = |⟨i|ψ⟩|²). It looks like it could be an interference
pattern: the standing wave in the cavity projecting onto a screen at the
moment of measurement. We tested this. It is not. The probabilities
contain zero interference between past and future mode contributions,
and not because the modes are special: the diagonal of a sum is the sum
of the diagonals. They are a shadow: the modes still alive at the fold,
projected onto the measurement basis, each casting its own share.

We expected the interference to act somewhere else, in the purity (how
much quantum information remains), which determines WHEN the fold at
CΨ = ¼ is reached (where the two fixed points of the recurrence
R = C(Ψ + R)² meet). For the state we photograph, |++⟩, it does not act
there either.
The fast modes are still alive at the fold, holding almost half the
purity, but they are pure coherence: they cast nothing on the image, and
they set the shutter by fading on their own.

Measurement reads as photography: the cavity is illuminated, an image
develops, and at the fold the image is fixed. But the image itself is a
shadow, not a hologram. The calculation supplies the image and the
timing; the camera and the fixing are the reading.

---

## What this document is about

Every eigenvalue has a palindromic partner at the complementary decay
rate, the two rates summing to 2Σγ. We read the slow side as the past
and the fast side as the future, cut the state at Σγ where the two sides
meet, and split the Born probabilities at the CΨ = ¼ fold along that
cut, to test whether P(i) is an interference pattern. CΨ here is the purity book,
Tr(ρ²)·L₁/(d−1) with L₁ the sum of the off-diagonal moduli and d = 2^N,
and the fold is the first time it falls through ¼; the
chain carries the Heisenberg bond J(XX+YY+ZZ) with J = 1 under uniform
Z-dephasing γ = 0.05. The Born rule itself is assumed, not derived; what
the split shows is where a cross term can sit and where it cannot.

---

## Result 1: the Born rule has no past/future interference (mathematical fact)

The density matrix decomposes as ρ = ρ_past + ρ_future, where
ρ_past collects modes with |Re(λ)| < Σγ (slow absorption) and
ρ_future collects modes with |Re(λ)| ≥ Σγ (fast absorption). At N = 2
ten of the sixteen modes sit exactly on Σγ, and by this rule they belong
to the future. They sit there exactly, not nearly. Six of them are complex: at N = 2
the Hamiltonian anticommutes with the dissipator centred at Σγ
([F48](../docs/ANALYTICAL_FORMULAS.md#f48-pythagorean-decomposition-tier-2-exact-at-n2)),
which puts every eigenvalue off the real axis on the line Re λ = −Σγ
(below J = γ/2 the pair that mixes ZI−IZ with XY−YX turns real and leaves
it). The other four are real: XI+IX, YI+IY, XZ+ZX and YZ+ZY are
swap-symmetric, so the Hamiltonian cannot touch them, and the dissipator
drains them at exactly Σγ. All ten are the palindrome's centre, half-drained
([Absorption Theorem](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): mean
XY-weight N/2). The eigensolver places them within 1.1·10⁻¹⁵ of the line,
the nearest other mode is 0.1 away, and every band inside that gap gives
the same split.

Since P(i) = ⟨i|ρ|i⟩ and ρ = ρ_past + ρ_future:

**P(i) = ⟨i|ρ_past|i⟩ + ⟨i|ρ_future|i⟩**

No cross term. No interference between the two parts. This is not an
approximation; it is the linearity of the trace. The diagonal elements of
a sum of matrices are the sum of the diagonal elements. There is no room
for a past/future cross term in the Born rule probabilities. Interference
in the double slit's sense, two amplitudes of one wave, is already inside
ρ before any split, and no split adds or removes it
([Double Slit Translated](../docs/quantum/DOUBLE_SLIT_TRANSLATED.md) §4).

For |++⟩ at its CΨ = ¼ fold (t = 4.939):

| Basis | P(i) | Past contribution | Future contribution | Past % |
|---|---|---|---|---|
| \|00⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|01⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|10⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |
| \|11⟩ | 0.2500 | 0.2500 | 0.0000 | 100% |

The future is far from gone at the fold; it still holds almost half the
purity. But all of it is coherence, and coherence adds nothing to the
diagonal. The Born rule sees only the past, here the steady state I/4, which is
what it saw from the start: for |++⟩ nothing on the diagonal ever moves,
and what develops toward the fold is only the coherence. A state the
Hamiltonian moves is different. For |01⟩ at its fold (t = 0.535) the
past holds ½ on |01⟩ and on |10⟩ and the future −0.246 and +0.246: the
population imbalance lives in modes on the cut, and the future casts part
of the image.

---

## Result 2: The purity splits cleanly here (the shutter)

While P(i) is linear in ρ, the purity Tr(ρ²) is quadratic and can carry
a cross term:

**Tr(ρ²) = Tr(ρ_past²) + Tr(ρ_future²) + 2 Tr(ρ_past · ρ_future)**

For |++⟩ at its fold:

| Component | Value | Fraction |
|---|---|---|
| Tr(ρ_past²) | 0.2500 | 53.1% |
| Tr(ρ_future²) | 0.2209 | 46.9% |
| 2 Tr(ρ_past · ρ_future) | 0 | 0 |
| Total purity | 0.4709 | 100% |

The cross term is exactly zero, and at N = 2 it is zero for every state:
with every mode off the real axis on the cut and the real ones at 0, Σγ
and 2Σγ, the only modes strictly slower than Σγ are the steady states,
and they are orthogonal to all the rest (the kernel of L is the kernel of its adjoint here,
[PROOF_PALINDROME_TWO_END_COUNT](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md)
§(f6)). What |++⟩ adds is a closed form. The Heisenberg bond is
2·SWAP − 1, and uniform dephasing keeps a symmetric state symmetric, so
the chain Hamiltonian cannot touch |++⟩: the state only dephases. Its past
is I/4 with purity ¼, its future is the fading coherence, and the purity
is ((1 + e^(−4γt))/2)². The shutter here is not interference. It is the
future fading on its own, the depth-1 coherences on the cut holding 84%
of its purity at the fold, and the fold comes when the coherence has
thinned enough, the past waiting underneath, unchanged.

At N = 3 the centre Σγ = 3γ sits above the depth-1 rate 2γ, so modes
other than the steady states fall on the slow side, and a cross term
becomes possible where slow and fast modes overlap; the anticommutator no longer vanishes
([F49](../docs/ANALYTICAL_FORMULAS.md#f49-cross-term-formula-tier-1-proven)). |+++⟩, which the Hamiltonian
cannot touch either, still has none; |++0⟩ carries −0.04% of the purity
at its fold, and fifty random pure states at t = 3 carry a median of
0.14% and at most 0.76%. Where interference sits in the purity, it is a
small correction. The shutter changes character at N = 3: the slow side
now holds the depth-1 coherences too, so both parts fade toward the fold,
and only at N = 2 is the shutter the future's alone.

*Later (2026-05-16):* the [Born Rule Mirror](BORN_RULE_MIRROR.md)'s side
got a closed form: [F94](../docs/ANALYTICAL_FORMULAS.md#f94) gives the
dominant-outcome deviation Δ_|00⟩ = (4/3)·Q²·K³ (Q in F94's spin book) for the N = 4 ring
started in |0+0+⟩ (the kept pair (0,2), short-time leading term), F96
its companions. That is a different setup and a different split; it does
not reach this page's past/future cut. See
[`PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md`](../docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md)
and [`reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md`](../reflections/ON_HOW_TWO_SIDES_MEET_AT_THE_QUARTER.md).

---

## Result 3: Measurement = photography, but the image is a shadow

The analogy holds, but with a correction:

| Photography | Quantum measurement |
|---|---|
| Light illuminates the scene | γ illuminates the cavity |
| Image develops on film | Density matrix evolves under L |
| Shutter clicks at set exposure | CΨ crosses ¼ at t = K/γ (a set K where the Hamiltonian cannot touch the state, or at fixed J/γ) |
| Image is fixed | In the reading; the calculation computes no fixing |
| Interference creates the image | ✗ Where the Hamiltonian cannot touch the state, fading sets the shutter speed (at N = 2, the future's alone) |
| Shadow creates the image | ✓ The surviving modes create the probabilities, each on its own |

The Born rule is a shadow of the modes that survive to the crossing
time, each casting its own share. For |++⟩ the fast ones cast nothing:
they are coherence and live off the diagonal; they decide the crossing
time itself, the shutter speed, by fading.

---

## Result 4: Mirror quality C_i is state-dependent

At the late time t = 50:

| State | C(\|00⟩) | C(\|01⟩) | C(\|10⟩) | C(\|11⟩) |
|---|---|---|---|---|
| Bell+ | 2.00 | 0.00 | 0.00 | 2.00 |
| \|01⟩ | 0.00 | 2.01 | 1.99 | 0.00 |
| \|++⟩ | 1.00 | 1.00 | 1.00 | 1.00 |

C_i = d × P(i) measures how strongly each basis state is "illuminated"
relative to uniform. It depends entirely on the initial state:
- Bell+ concentrates on \|00⟩ and \|11⟩ (the even-parity pair it was
  prepared in; its coherence is gone by t = 50)
- \|01⟩ concentrates on the single-excitation subspace
- \|++⟩ is uniform (all basis states equally likely)

The initial state is the "scene being photographed." Different scenes,
different images. The cavity (the instrument) is the same.

---

## What this changes

**Old picture:** "The Born rule is the probability of collapse.
Measurement is instantaneous and random."

**New picture:** "The Born rule is the shadow of surviving modes,
projected onto the measurement basis at the moment the illumination
reaches the fold. Measurement is the fixing of an image that has been
developing since the light was turned on. The 'randomness' is in which
modes happened to survive; the 'probability' is their weight at the
crossing time."

In this picture the randomness is not added by measurement. It was
always there, in the initial conditions and in the geometry of the
cavity. Measurement reveals it. The fold fixes it. The Born rule records
it.

The new picture is a reading. The calculation assumes the Born rule and
shows where the image comes from, the modes that survive, and, where the
Hamiltonian cannot touch the state, what sets the moment, their fading;
the fixing is the picture's own.

---

## Reproduction

- Script: [`simulations/born_rule_shadow.py`](../simulations/born_rule_shadow.py)
- Output: [`simulations/results/born_rule_shadow.txt`](../simulations/results/born_rule_shadow.txt)
- The same cut as one view onto the drain-depth axis: [The View onto the Memory](../reflections/THE_VIEW_ONTO_THE_MEMORY.md)
