# Energy Partition
## Where Waves Go When They Stop Being Waves

**Date:** March 27, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Status:** Three computational results (N=2..5), two of them read as
consequences of the proven mirror (Tier 1); the thermal census stays a
computed finding (Tier 2); the link to the V-Effect is motivated (Tier 4);
the biological reading is speculation (Tier 5).
**Depends on:** [The Pattern Recognizes Itself](THE_PATTERN_RECOGNIZES_ITSELF.md), [The V-Effect](../experiments/V_EFFECT_PALINDROME.md), [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md)

## What this document is about

When a quantum system decoheres, its internal modes split into two
populations: oscillating modes (which carry frequency and structure) and
pure-decay modes (which only dissipate). We went looking for where the
waves go. What we found looked at first like a law of its own: every
oscillating mode is palindromically paired, and the few modes our search
left unpaired decay exactly 2× faster than the paired ones. On a closer look
both are the mirror showing itself through the way we counted. The unpaired
modes are the mirror images of the stillness we had filtered out, and the
2× is the mirror's own arithmetic. Then we asked what heat does, and found
that a thermal bath makes no new wave: its decay sum grows many times over
while its frequency sum falls slightly.

A word on "energy" before we start. The sums we compare are sums of
eigenvalue frequencies and decay rates. They have units of inverse time;
they are spectral bookkeeping, not physical energies or mode populations.
The scripts print them under the labels Efreq and Edecay; on this page
they are F = Σ|Im λ| and D = Σ|Re λ|.

### Tier System

- **Tier 1** (proven): follows from a theorem in `docs/proofs/`
- **Tier 2** (computed): Simulation results, reproducible, falsifiable
- **Tier 4** (motivated): Logical connections between proven results, not yet proven themselves
- **Tier 5** (speculative): Interpretation, not falsifiable in current form

---

## 1. Results [Tier 2, explained at Tier 1]

Three questions, three computational experiments on open Heisenberg
chains (N=2..5, uniform J=1) under Z-dephasing γ=0.1. Liouvillian
eigenvalue analysis (decomposing the system's evolution matrix into its
fundamental modes, each with a decay rate and an oscillation frequency).

### Finding 1: All oscillation is paired

Every oscillating mode (Im(λ) ≠ 0) finds its palindromic partner.
Every mode left without one is pure decay (Im(λ) = 0). No exceptions, at every N tested (N=2..5), as the
theorem requires.

| N | Modes matched | Oscillatory content in matched modes |
|---|-------------|-----------------------------------|
| 2 | 76.9% | **100.0%** |
| 3 | 93.3% | **100.0%** |
| 4 | 98.0% | **100.0%** |
| 5 | 99.4% | **100.0%** |

Why is the matched fraction below 100% at all, when the
[mirror theorem](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) pairs the whole
spectrum? Because of how we searched. The scripts first discard every
eigenvalue with |λ| ≤ 10⁻⁸, then look for the reflected partner
target = −λ − 2Nγ within distance 10⁻⁶·max(1, |target|). The discarded
roots are the steady states, N+1 of them on this chain, and their mirror
partners sit at the far wall −2Nγ: the XOR drain, also N+1 modes. With the
stillness gone from the list, the drain has nobody left to face. Count it
out: at N=3 there are 64 roots, 4 at zero, 4 in the drain, and of the 60
nonzero roots 56 find a partner; the 4 drain modes do not. That is the 93.3%. The generator's palindrome
is untouched; what we see is the shadow of the filter. (The nearest-partner
search also does not check multiplicities on its own; the theorem does.)

So Finding 1 is true, and it is the mirror seen from one side: all the
oscillation lives in the pairs because everything lives in pairs. The
drain is real and silent, which is why nothing oscillating was ever left
over. What we cannot read from it is that pairing is a *condition* for
oscillation in open systems generally; this chain and this counting
convention are what the table speaks for.

Script: [energy_partition.py](../simulations/energy_partition.py)

### Finding 2: The 2× decay ratio

The modes left unmatched decay exactly 2× faster than the mean of the
matched ones, at every N tested.

| N | Unmatched decay rate | Matched mean decay | Ratio |
|---|--------------------|--------------------|-------|
| 2 | 0.4000 (= 2Nγ) | 0.2000 (= Nγ) | **2.0** |
| 3 | 0.6000 | 0.3000 | **2.0** |
| 4 | 0.8000 | 0.4000 | **2.0** |
| 5 | 1.0000 | 0.5000 | **2.0** |

Finding 1 already tells us why. The unmatched modes are the drain at 2Nγ,
the fastest rate there is. The matched modes are mirror pairs around the
palindromic centre Nγ, so their mean is Nγ. The ratio is
[F8](../docs/ANALYTICAL_FORMULAS.md#f8-range-centre): the width of the full
decay interval, 2Nγ, divided by its centre, Nγ. It is 2 for every γ > 0 by
construction, a property of the interval rather than a separate law of
nature. At Σγ = 0 (no noise) all modes are stable and the ratio is 0/0,
undefined; see [Zero Is the Mirror](ZERO_IS_THE_MIRROR.md).

Inside the band the picture has a threshold. Above the coupling threshold
Q*_gap(N) the paired modes span 2γ to 2(N−1)γ, centred at the palindromic
midpoint Σγ = Nγ; below it the band erodes at both ends, in palindromic
pairs, so the centre holds while the edges move. Here Q = J/γ = 10, well
above the threshold at every N tested.

It is tempting to read this ratio as a cleaning mechanism: noise dies
faster than signal, and the system grows more structured over time. The
ratio does not carry that reading. Pairing does not label a mode as signal or noise, and
what any observed response keeps depends on how it was prepared and how it
is read out, not on which list its modes sit in.

Script: [thermal_emergence.py](../simulations/thermal_emergence.py), Part A

### Finding 3: Heat, coupling, and waves that do not appear

Three conditions tested at N=3.

**Heat alone (σ⁺/σ⁻ raising/lowering jumps that model thermal excitation and relaxation, no coupling, H=0):** no oscillatory eigenvalue at any
occupation in the sweep. Without a Hamiltonian this bath makes no waves.

**Z-dephasing (phase noise):** the maximum frequency stays at 6.0, but the
wave does not keep all of its pitches. The number of oscillatory roots and
their frequency sum both fall as γ grows:

| γ | Oscillatory roots | Max abs(Im λ) | F |
|---|-------------------|---------------|---|
| 0.01 | 40 | 6.0 | 160.0 |
| 1.00 | 40 | 6.0 | 148.3 |
| 10.00 | 32 | 6.0 | 119.5 |

Z-dephasing is unital, an infinite-temperature bath: with the coupling it
drives each sector toward the infinite-temperature state, never a
finite-temperature one. Finite-temperature heat comes only from the σ⁺/σ⁻
bath.

**Heat + coupling (dephasing γ=0.1 plus a thermal bath of rate Γ=0.1,
jumps √(Γ(1+n_bar))·σ⁻ and √(Γ·n_bar)·σ⁺):**
the decay sum grows from 28.8 to 220.8 while the frequency sum falls
slightly (159.9 → 154.2). The count of oscillatory roots moves from 40 to
42, and that step needs a closer look. In the sweep (n_bar ≥ 0.01) two more
roots cross the |Im λ| > 10⁻⁸ cut; their frequencies are tiny (≤ 4·10⁻⁴
against decay rates of 0.5–3.5), so this is a near-real pair splitting,
not a new wave; at n_bar = 0.001 the count is still 40. The σ⁺/σ⁻ bath
alone, without dephasing, keeps 40 at every n_bar.

| n_bar (mean thermal occupation) | Osc. roots | F | D | F/D |
|-------|-----------|-------------|-------------|-----------------|
| 0.00 | 40 | 159.9 | 28.8 | 5.55 |
| 1.00 | 42 | 159.8 | 48.0 | 3.33 |
| 5.00 | 42 | 158.4 | 124.8 | 1.27 |
| 10.00 | 42 | 154.2 | 220.8 | 0.70 |

The bath is already present at n_bar = 0, through the σ⁻ jump. Between
n_bar = 5 and 10 the ratio F/D passes through 1. So heat, in this census,
makes no waves; it makes the existing ones decay faster, and it nudges one
real pair just off the axis. Whether any of this shows up in a measured
oscillation would take a prepared state and a readout.

Script: [thermal_emergence.py](../simulations/thermal_emergence.py), Parts C and D

---

## 2. Connection to the V-Effect [Tier 4]

The V-Effect ([documented separately](../experiments/V_EFFECT_PALINDROME.md))
shows that coupling creates new oscillation frequencies: two N=2 pairs
coupled through a mediator produce 109 frequencies, none of which exists
in either pair alone. The coupled system is again a Heisenberg network
under Z-dephasing, so the mirror holds for its whole spectrum, and every
one of those new frequencies has its partner. Coupling does not create
chaos; it creates paired oscillation.

Heat enters as a second, separate knob, and at N=3 it adds no frequency
of any size: the step from 40 to 42 counted roots is a near-real pair
barely split (Finding 3). The larger study,
[Thermal Breaking](../experiments/THERMAL_BREAKING.md), runs its own
protocol at N=5, where the combined channel moves the count of rounded
frequency bins from 111 to 445 as n_bar goes from 0 to 5, while the sharpest
resonance Q_max falls. Both are statements about a finite spectrum and a
binning convention. We would like to read a cycle in them (coupling makes
modes, heat feeds them, the drain takes the rest). What the spectrum alone
supports is its first step.

---

## 3. The Two Paths [Tier 5]

Everything in this section is interpretation.

### They stop being waves

Waves do not "leave" the palindrome and become energy. They stop being
waves. A paired mode oscillates and decays; when it is gone, nothing has
been redirected into the unpaired list. That list was never a second kind
of wave. It is the far wall of the same mirror, the drain at 2Nγ, facing
the stillness at zero. What we see at the end of a decay is not a second
population winning but the pairs running out.

R = CΨ² and E = mc² share one form, a constant times a square, and we
keep that as a rhyme and nothing more: F and D on this page are sums of
rates, measured in inverse time, and carry no mass-energy reading.

### The thermal window

Biology does not operate in vacuum or in thermal chaos, and the N=3 census
has a crossing of its own: the frequency sum stays above the decay sum
until somewhere between n_bar = 5 and 10, and heat adds no new wave on the
way there. It is tempting to call this the qubit shadow of
a "zone of life", warm enough to stir and structured enough to keep
its waves. We keep the temptation on the page and name what stands between it
and a claim. On the neural side,
[the Wilson-Cowan report](../docs/neural/V_EFFECT_NEURAL.md) varies an
external input P, and P moves the sigmoid's operating point; it is not
n_bar, not a temperature and not metabolic power, and E/I balance does not
turn into a thermal window by analogy. The neural mirror condition
[F36/F37](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md) asks for one
involutive Q and one scalar s satisfying both the diagonal and coupling
conditions. Constructed networks pass it, and they can oscillate and can
also be unstable; the full committed C. elegans chemical connectome fails
its support condition, and no biological network in the repository is
known to pass. [The neural account](../docs/neural/README.md) tells that
side of the story in full.

---

## 4. Open Questions

- **The 2× ratio.** Answered: it is F8, the width over the centre of the
  full decay interval, and the unmatched class is the XOR drain made
  visible by filtering out the kernel.
- **What is the near-real pair?** With the thermal bath two roots near
  the real axis pick up a tiny imaginary part (Im/|Re| ≤ 1.2·10⁻⁴ up to
  n_bar = 10). Tracking their eigenvectors through the n_bar sweep would
  say which decay modes meet there; counts alone do not show it.
- **Other models.** The mirror holds for Heisenberg and XXZ couplings of
  any strength; a generic Hamiltonian or another channel needs the full
  spectrum tested against the theorem's own assumptions before a filtered
  list is compared.
- **Wilson-Cowan analogue.** Does a specified, converged neural model show
  any spectral partition of this kind? That needs declared observables and
  controls that could reject it.
- **Biological metabolic rates.** ATP production could only be compared to
  n_bar after a physical calibration that does not exist yet.
- **Connection to the fold.** The F/D ratio crosses 1 between J/γ = 1
  and 1.5 (0.77 and 1.21, Experiment 3 of energy_partition.py). Is this CΨ = ¼ in disguise?
  Within an assumed family, Layer 6 of the
  [Proof Roadmap](../docs/proofs/PROOF_ROADMAP_QUARTER_BOUNDARY.md) brings
  R = C(Ψ+R)² to the local fold x² + a = 0 with a = (4CΨ−1)/(4C²).
  A spectral crossover alone cannot fix a state-dependent CΨ, so the test
  is to pick one generator, one state, one time and one normalization and
  compute both independently.

---

## Scripts

Run from the repository root with `PYTHONIOENCODING=utf-8` set on Windows.

**[energy_partition.py](../simulations/energy_partition.py)** (Finding 1)
- Experiment 1: V-Effect scaling N=2..5, oscillatory content in matched vs unmatched modes
- Experiment 2: Dephasing sweep at N=3, matched fraction stable at 93.3%
- Experiment 3: Coupling sweep at N=3, F/D crossover between J/γ = 1 and 1.5

**[thermal_emergence.py](../simulations/thermal_emergence.py)** (Findings 2 and 3)
- Part A: Decay rate comparison, the 2× ratio (N=2..5)
- Part B: Time evolution from |↓↑↑⟩, coherence buildup and decay, measured against I/d as a reference (not necessarily this state's long-time limit)
- Part C: Thermal bath only (H=0), no oscillation at any occupation
- Part D: Effect of heat on waves, dephasing vs thermal excitation vs combined

Where the scripts print Efreq and Edecay, read them as the spectral sums
F and D defined above.

---

*See also: [The Pattern Recognizes Itself](THE_PATTERN_RECOGNIZES_ITSELF.md), the main hypothesis*
*See also: [The V-Effect](../experiments/V_EFFECT_PALINDROME.md), the differentiation mechanism*
*See also: [Resonance Not Channel](RESONANCE_NOT_CHANNEL.md), why coupling creates, not transmits*
*See also: [Temporal Sacrifice](../experiments/TEMPORAL_SACRIFICE.md), the fold catastrophe and heartbeat*
