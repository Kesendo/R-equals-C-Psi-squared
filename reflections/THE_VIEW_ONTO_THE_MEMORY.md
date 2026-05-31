# The View Onto the Memory

**Date:** 2026-05-28
**Status:** Reflection / synthesis. Ties existing Tier-1 results into one reading. No new claim, one recognition.
**Authors:** Thomas Wicht, Claude (Opus 4.7)

> The carbon "painter" figure is not a picture of a relaxation time. It is the view onto the memory: how a system stores what it has been, and how it carries that forward.

![The N=4 carbon ring's eight slowest Liouvillian modes, each with its decay rate and its "popcount-mean" storage reading. The slow modes (small rate) are the memory; the fast modes are flushed first.](../visualizations/memory.jpg)

---

## Start here, if the words are new

A molecule, a spin, any small quantum system sits in a noisy world that constantly nudges it. Some of what the system holds about its own past is washed out almost at once. Some of it lingers for a long time. The part that lingers is its **memory**: the shape it carries forward from the past into the future.

This project found the exact rule for which part is kept and which is lost. The picture above reads that rule off one real example, a four-site carbon ring. Each row is a **mode**, one pattern the system can be in. The number beside it, `Re(λ)`, is how fast that pattern fades: small means long-lived, large means gone quickly. The slowest patterns are the memory. The fastest are forgotten first.

Everything below is one idea seen from several sides: **how a system keeps what it has been**. You do not need the math to take the idea. If you want to go deeper, every claim links to where it is proven.

---

## One axis, six names

The trunk is the **Absorption Theorem** ([`PROOF_ABSORPTION_THEOREM.md`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md)):

    |Re(λ_k)| = 2γ · ⟨popcount(i XOR j)⟩    for every mode k

`popcount(i XOR j)` is the **drain depth** of a piece of the state: the number of sites at which two basis configurations differ, which is how many places the noise can grip. Depth zero (the diagonal) means rate zero, never fading, stored forever. High depth means flushed fastest, at up to `2Σγ`. Every other name in this project is a reading of that one axis:

| Lens | slow end (low depth) | fast end (high depth) | Where it lives |
|------|----------------------|-----------------------|----------------|
| Drain depth | popcount 0, never fades | popcount N, fastest | Absorption Theorem |
| Spectral mode | small \|Re(λ)\|, long-lived | large \|Re(λ)\|, lost | F3, D6 (gap = 2γ) |
| Time | {I,Z}, decided, classical, **past** | {X,Y}, undecided, quantum, **future** | [`PI_AS_TIME_REVERSAL.md`](../experiments/PI_AS_TIME_REVERSAL.md) |
| Born shadow | ρ_past (\|Re(λ)\| < Σγ), **97%** | ρ_future (≥ Σγ), ~3% with interference | [`BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md) |
| Memory | static part (kernel of L) | dynamical part that fades | F88b, [`MemoryAxisRho.cs`](../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs) |
| Information | what survives the palindrome | the XOR drain, stores nothing | [`XOR_SPACE.md`](../experiments/XOR_SPACE.md) |

These are not analogies. They are the **same split of the state**, made by the same operator L, read through six lenses. The conjugation Π makes the pairing exact (`Π·L·Π⁻¹ = −L − 2Σγ·I`): every slow mode has a fast partner whose rates sum to `2Σγ`. The reflection [`ON_TWO_TIMES.md`](ON_TWO_TIMES.md) already says it in words: "memory is the shape that survives long enough to be re-recognized; the mode with the smaller |Re(λ)| determines how far back the envelope remembers; the XOR drain stores nothing at all."

---

## The picture is this axis, mode by mode

The figure sorts the carbon ring's slow modes by exactly this depth. Reading the state-side popcount distribution of each mode (verified in [`_carbon_painter_xor_and_depth.py`](../simulations/_carbon_painter_xor_and_depth.py)):

- the steady state is **100% depth-0**, pure stored past;
- the slowest mortal mode (rate 0.172) is **91.8% depth-0**, almost all memory, with a thin tail that will fade;
- the fastest modes (rate 2.127) are **89.9% depth-1**, future, flushed first.

And the figure's "popcount-mean identity" is the Absorption Theorem read state-side: the mean drain depth of each mode equals its decay rate, bit-exact across all eight modes (0.1723, 0.2190, 0.5969, 0.9014, 2.0670, 2.1273, twice). The figure is the storage map.

The two painter towers (Y content versus non-Y content) are then just **which part of the past the molecule holds in which transverse axis**. The reported T2 anisotropy, `T2(x)/T2(y) = 0.2190/0.1723 = 1.27`, is the statement that x-memory and y-memory fade through different rungs of the same drain. It is exact as a ratio of two mode rates (so exact as a ratio of two mean depths, by the Absorption Theorem), but it is not a simple closed-form fraction: the once-conjectured 4/3, 8/7, 14/13, 20/19 sequence does not hold.

---

## The bridge into the classical world

This is why substrate inheritance is more than inheriting decay rates. A substrate mapped onto the Liouvillian (a two-state unit per site, a site-local noise bath, a coupling graph) inherits the **whole past, future, and memory split**, because the Absorption Theorem depends only on the noise, not on the Hamiltonian. It holds for any Hermitian Hamiltonian, real or complex. Carbon, water, a neural Jacobian: each, once mapped, carries the same split.

And the split is the bridge. The **classical world is the slow, low-depth, {I,Z}, stored part**: the decided past that survives long enough to be re-recognized. The **quantum world is the fast, high-depth, {X,Y}, undecided part**, drained at `2Σγ`. "The bridge into the classical world" is not a separate construction. It is reading the slow rim of this one axis. A molecule's classical, measurable, remembered identity is the depth-0 floor of its own Liouvillian. Everything above it is the quantum tide the bath is already pulling out.

---

## Honest seam: the 97/3 number

The literal **97 / 3** lives in the Born-rule branch, and it surfaces there in two windows: as `Tr(ρ_past²) = 97.1%` for the |++⟩ pair at the CΨ = 1/4 crossing, the rest being future plus interference ([`BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md)), and as the 97%-Hamiltonian / 3%-decoherence-correction split of the |0+0+⟩ Heisenberg pair's Born probabilities ([`BORN_RULE_MIRROR.md`](../experiments/BORN_RULE_MIRROR.md)); both closed as **F94**, `Δ = (4/3)·Q²·K³`. The carbon painter system shows the structurally identical split, a slow majority over a fast minority (about 98/2 by mode count, the generic Heisenberg tally in [`XOR_SPACE.md`](../experiments/XOR_SPACE.md); the figure's own slowest mortal mode is 91.8% depth-0 pure) but the exact 97/3 numeral is the sibling instance, not this system. The painter figure's own "3%" is unrelated: it is the finite-time fitting gap on the T2 ratio. Same split, two windows onto it. The recognition is right, and the number has a specific home.

---

## What this figure also held (seen 2026-05-29)

We drew this figure as one molecule's storage map and filed it. A more general view was sitting in it the whole time.

Each mode here already comes with a **storage portfolio** read in percentages: how much of its content sits in the Y axis versus the non-Y axes (the two painter towers), how much at depth 0 versus depth 1. And the state-side identity prints the rate as the mean of that storage, `|Re(λ_k)| = 2γ·⟨popcount(i XOR j)⟩`. We read the depth as one number against one noise rate γ, because the carbon bath is uniform: every site dephases at the same γ, so the vector of rates lay flat and the percentages looked like a single count.

The deeper law lets the noise be a **vector**, one rate per channel, and reads any mode's rate as the γ-weighted version of that same percentage portfolio:

    −Re(λ_k) = 2 · Σ_x γ_x · ⟨Δ_x⟩_k,   ⟨Δ_x⟩_k ∈ [0, 1]

`⟨Δ_x⟩` is the percent of its difference the mode stores in channel x (on a basis coherence, whether bra and ket disagree there); `popcount` with a single γ is the all-channels-equal case. The carrier is a vector; the rate is its pairing with the mode's portfolio.

Where the clocks differ, the portfolio becomes the whole story. A ³¹P donor in silicon carries channels spanning decades, electron spin fast, nuclear spin very slow and protected, charge sensitive, valley medium. There a coherence lives as long as it stores its difference only in the slow channels: the nuclear spin's long memory is that, left alone, its coherence is almost entirely nuclear difference, and hyperfine coupling that mixes even a few percent of the fast electron channel into the mode costs decades. Clock transitions and sweet spots are the same instruction in the portfolio's own words, hold zero percent in the fast channels. (`simulations/_sip_carrier_channels.py`; the per-mode law verified in `simulations/_absorption_gamma_vector.py`; the carrier-vector reading now in [`AbsorptionTheoremClaim`](../compute/RCPsiSquared.Core/Symmetry/AbsorptionTheoremClaim.cs).)

The honest part is that the figure handed us the perspective and we took it for a relaxation time. The per-mode percentages were printed here, the state-side mean was printed here, even the leftover sat beside it as a percent. The percentages were the lens the whole time. We learned to read them as the general view, the carrier as a vector and each mode as a portfolio, only after looking at a substrate whose channels refused to weigh the same. The structure was always on the screen; the seeing is what was new.

---

## Seen again 2026-05-31: the axis has a rail, and a currency

We came back to this axis from the other end, the post-EP flow
([THE_FLOW_BETWEEN_TWO_SINGULARITIES](../experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md)), and it
handed back two things the figure did not yet name: a **rail** running along the axis, and a
**currency** that prices every rung. None of the pieces is new; the project found each before. What
was hard to gather is that they are one ladder.

**The rail: parity.** The drain depth popcount(i⊕j) is, for any basis coherence |i⟩⟨j|, exactly the
light content n_XY (the sites that differ are the sites carrying {X,Y}); one axis, two names. And the
axis carries a parity, depth mod 2, which is the change in particle number:

- **even rungs** (depth 0, 2, 4 …) are NUMBER-CONSERVING: the diagonal memory and the |i⟩⟨j|
  coherences a fixed-number state holds. The dark, intra-sector rail the single-excitation flow rides.
- **odd rungs** (depth 1, 3 …) are NUMBER-CHANGING: the |n⟩⟨n±1| coherences coupling adjacent
  particle-number sectors. The bright rail an optical dipole grips, the birth channel.

Because the rate is the depth (the Absorption Theorem, |Re λ| = 2γ⟨depth⟩), the rail is a ladder of
rates. depth-0 (even) is the kernel, the steady state, the memory this whole figure is about. depth-1
(odd) is the slowest MORTAL mode, rate 2γ, the longest-lived number-changing coherence, the birth
channel. depth-2 (even) is the flow's own coherence, rate 4γ, faster than the birth channel. So the
single-excitation flow, riding the even rungs, relaxes through depth-2 and never touches the depth-1
birth channel: its weight on odd parity is identically 0, and the C# GameObject measures its
between-block content as machine zero
([PostEpFlowField](../compute/RCPsiSquared.Diagnostics/Foundation/PostEpFlowField.cs)). The birth
channel is the slowest mortal mode of the whole Liouvillian, but a latent door: number-conserving
dynamics never opens it.

**The currency: the bilinear p(1−p).** Every rung is priced by one shape. The per-site light is
2·p(1−p) with p the single-site occupation, peaking at ½ on the Bloch equator. The maximal
between-block coherence is C_block = p(1−p) with p the weight between two adjacent number sectors,
peaking at ¼ for the Dicke superposition (|D_n⟩+|D_{n+1}⟩)/√2
([BlockCoherenceContent](../compute/RCPsiSquared.Core/F86/BlockCoherenceContent.cs), Theorem 2). Same
bilinear x(1−x); the recurring ½ and ¼ that thread this project, carbon's half-filling, the polarity
triple, the coherence ceiling, are its doubled apex and its apex, with p wearing two hats: an
occupation, or a sector weight.

**The crown switches rails.** Which rail holds the longest-lived mode depends on the coupling. The
even (flow) modes' rate rises with Q = J/γ (stronger coupling forgets the occupation faster); the odd
birth channel stays near 2γ. They cross: below the crossing the longest memory is the even occupation,
above it the odd birth channel. For the uniform chain the crossing sits near the exceptional point,
but that coincidence is uniform-specific, a non-uniform γ-profile breaks it (the tempting "EP rate =
2γ" and "⟨n_XY⟩ = 1 at the EP" are uniform-only, refuted at peaked-V; recorded honestly rather than
kept).

**The carbon mirror.** bright/dark is the rail seen optically: the bright band edge is the odd,
number-changing rung a photon couples to; the dark 2Ag triplet pair is the even, number-conserving
rung, the singlet-fission object, the same bound pair we isolate
([SINGLET_FISSION_AND_THE_TWO_CLOCKS](../docs/carbon/SINGLET_FISSION_AND_THE_TWO_CLOCKS.md)). What the
chemist calls dark is what this axis sees clearly.

So the full lock: depth that is light that is rate, split even and odd into flow and birth, priced
everywhere by p(1−p), mirrored slow to fast by Π (the palindrome, rates summing to 2Σγ), with the
memory at its foot. The six lenses above all read the one axis; the rail is the perpendicular split,
and the bilinear is what sits on each rung. Probes:
[`flow_depth_parity.py`](../simulations/flow_depth_parity.py),
[`light_content.py`](../simulations/light_content.py),
[`bound_pair_light.py`](../simulations/bound_pair_light.py), and this session's parity and
birth-channel sweeps.

---

## What this holds

Nothing here is a new theorem. It is the statement that the painter alternation, the popcount / XOR storage, the Born shadow's 97/3, F88b's static-versus-memory decomposition, Π's past and future, and the classical bridge are **one axis**, the drain-depth axis the Absorption Theorem quantizes, seen from six sides. The figure that looked like a relaxation measurement is the clearest single picture of it: the view onto the memory, mode by mode, of how a state keeps what it has been.

---

## Go deeper

- [`PROOF_ABSORPTION_THEOREM.md`](../docs/proofs/PROOF_ABSORPTION_THEOREM.md): the trunk, `Re(λ) = −2γ⟨popcount⟩`, holds for any Hermitian H.
- [`PI_AS_TIME_REVERSAL.md`](../experiments/PI_AS_TIME_REVERSAL.md): {I,Z} is past, {X,Y} is future, Π maps one to the other.
- [`BORN_RULE_SHADOW.md`](../experiments/BORN_RULE_SHADOW.md), [`BORN_RULE_MIRROR.md`](../experiments/BORN_RULE_MIRROR.md): ρ_past and ρ_future, the 97/3 split, closed form F94.
- [`PROOF_F86B_UNIVERSAL_SHAPE.md`](../docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md) (F88b) and [`MemoryAxisRho.cs`](../compute/RCPsiSquared.Diagnostics/Foundation/MemoryAxisRho.cs): the static-versus-memory decomposition of the state.
- [`ON_TWO_TIMES.md`](ON_TWO_TIMES.md): memory as a standing wave, and the two times (the noise time that flows, the felt time with a horizon set by the slowest mode).
- [`XOR_SPACE.md`](../experiments/XOR_SPACE.md), [`GLOSSARY.md`](../docs/GLOSSARY.md): where information lives in the palindrome, and the XOR fraction as how fast a state is drained.
- [`PAINTER_ALTERNATION_NMR_BRIDGE.md`](../docs/carbon/PAINTER_ALTERNATION_NMR_BRIDGE.md) and [`_carbon_painter_xor_and_depth.py`](../simulations/_carbon_painter_xor_and_depth.py): the carbon figure and its state-side storage read.
