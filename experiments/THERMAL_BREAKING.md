# Thermal Breaking: A Warmer Bath Trades Q-Factor for Frequency Diversity

<!-- Keywords: thermal amplitude damping sigma minus sigma plus finite
occupation spectral census F137 palindrome centre frequency bins V(N)
1+cos(pi/N) golden ratio F6 Q gain sacrifice profile -->

**Status:** V(N) = 1 + cos(π/N) is exact on the (0,1) coherence block under
uniform Z-dephasing and proven for all N through D10; that the block also holds
the Liouvillian's best Q is measured at N = 2 to 6 and open beyond. The thermal
palindrome is exact at N = 2 and 3 (rational characteristic polynomial) and at
the eigensolver floor to N = 5. The self-consistent occupation is excluded by
an identity for β > 0. The frequency counts are a four-decimal binning protocol
on a finite grid; n̄ is an external channel parameter throughout.
**Date:** March 30, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Producers:** [v_effect_gamma_sweep.py](../simulations/v_effect_gamma_sweep.py),
[v_effect_thermal.py](../simulations/v_effect_thermal.py),
[thermal_palindrome_centre.py](../simulations/thermal_palindrome_centre.py),
[self_heating_fixpoint.py](../simulations/self_heating_fixpoint.py)
**Outputs:** [v_effect_gamma_sweep.txt](../simulations/results/v_effect_gamma_sweep.txt),
[v_effect_thermal.txt](../simulations/results/v_effect_thermal.txt),
[thermal_palindrome_centre.txt](../simulations/results/thermal_palindrome_centre.txt),
[self_heating_fixpoint.txt](../simulations/results/self_heating_fixpoint.txt)
**Depends on:** [V-Effect](V_EFFECT_PALINDROME.md),
[D10, the (0,1) block dispersion](../docs/proofs/derivations/D10_W1_DISPERSION.md),
[F137](../docs/ANALYTICAL_FORMULAS.md#f137)

---

## What this document is about

A quantum resonator has two qualities one might want: a high Q-factor (how
many times a mode oscillates before it dies) and frequency diversity (how many
different frequencies the spectrum carries). This document follows three knobs
on a Heisenberg chain and finds that they pull these two qualities apart.
Coupling sets how the best Q grows with the chain, through a closed form,
V(N) = 1 + cos(π/N), in which the dephasing rate cancels. Z-dephasing splits
degenerate frequencies and leaves the palindromic pairing exact. A bath with a finite thermal occupation
n̄ multiplies the frequency count and lowers Q, and the sacrifice profile's
Q advantage fades as n̄ grows. The palindrome survives the bath at a centre that
moves with the total rate, and breaks only where amplitude damping shares its
axis with Z-dephasing.

Q-factor here is Q = |Im λ| / |Re λ| of a Liouvillian eigenvalue, and
Q_max is its largest value over the spectrum. A "frequency" is one distinct
value of |Im λ| after rounding to four decimals, over the modes with nonzero
rate and frequency. The counts are properties of that protocol and of the
finite grid, not resolution-free mode counts.

## The channel that was run

The chain is an open isotropic Heisenberg chain, H = J·Σ (XX + YY + ZZ) with
J = 1. The bath acts on each qubit through

```text
L_down = √(γ_amp·(n̄ + 1)) · σ⁻
L_up   = √(γ_amp·n̄)       · σ⁺
σ⁻ = |0⟩⟨1|,  σ⁺ = |1⟩⟨0|
```

in the repository's convention, |0⟩ the ground state. Spontaneous emission
remains at n̄ = 0, where absorption vanishes. A one-qubit endpoint gate checks
the direction directly: at n̄ = 0, |0⟩⟨0| is fixed and |1⟩⟨1| loses excited
population at γ_amp; at n̄ > 0, population rises from |0⟩⟨0| at γ_amp·n̄.
Swapping the two operators fails the same gate with residual 0.2, while the
implemented direction has residual 0.0 in the spectral producer and 2.78e-17 in
the channel audit.

The swap leaves every eigenvalue unchanged. The Heisenberg coupling and the
Z-dephasing are both invariant under conjugation by X^⊗N, which exchanges σ⁻
and σ⁺ on every site, so the swapped generator is the exact similarity
(X^⊗N ⊗ X^⊗N)·L·(X^⊗N ⊗ X^⊗N) of the correct one, with the entrywise
difference 0.0 at N = 2 and 3 (n̄ = 0.7, with and without Z-dephasing). Every Q
and every frequency count on this page is therefore the same under either
convention. The F1 distances are a different kind of number: the canonical
scorer is greedy, so one spectrum has an orbit of readings, and the value
printed depends on the order in which the eigensolver returns the eigenvalues,
which the permutation changes. The amplitude-only rows stay at the floor and
the co-axial rows stay broken under either label; their digits are not
spectral. Among the page's readings, only the sign of the steady-state ⟨Z⟩
tells the two conventions apart; a reading that fixes a basis state (an initial
or ancilla state) would as well.

The occupation n̄ is supplied from outside. No producer here contains a law by
which decay raises n̄, so none computes self-heating.

---

## The Q gain V(N) and why γ cancels

The Q-factor gain Q_max(N=5) / Q_max(N=2) was measured across 19 uniform
dephasing rates, γ/J = 0.001 to 5.0. All 19 rows read the same ratio:

| γ/J | Q_max (N=2) | Q_max (N=5) | ratio |
|:----|:------------|:------------|:------|
| 0.001 | 2000.0 | 3618.0 | 1.81 |
| 0.01 | 200.0 | 361.8 | 1.81 |
| 0.1 | 20.0 | 36.2 | 1.81 |
| 1.0 | 2.0 | 3.6 | 1.81 |
| 5.0 | 0.4 | 0.7 | 1.81 |

The absolute Q scales as J/γ; the ratio does not move. It is
V(5) = (5 + √5)/4 = 1.80902…, exactly, for the following reason.

The best Q sits on the chain's **(0,1) coherence block**: the N-dimensional
span of the |0⟩⟨j| between the ferromagnetic vacuum |0…0⟩ and the single
excitations |j⟩, the object F2's dispersion describes
([D10](../docs/proofs/derivations/D10_W1_DISPERSION.md)). Every |0⟩⟨j|
disagrees with the vacuum at exactly one site, so uniform Z-dephasing acts on
the whole block as the scalar −2γ (D10 Step 1): every mode of the block decays
at 2γ, whatever its frequency. Therefore

    Q_max = ω_max / (2γ)

and the factor 2γ, the same at every N, cancels from the ratio:

    V(N) = Q_max(N) / Q_max(2) = ω_max(N) / ω_max(2)

The gain is a frequency ratio. It is independent of the dephasing RATE; it is
not independent of the dephasing PROFILE, as "Uniform γ only" below shows.

The same number is also a ratio inside one chain. The block's mean Q over its
oscillating modes is 2J/γ, because Σ_{k=1}^{N−1} cos(πk/N) = 0, and that is
exactly Q_max(2). So V(N) = Q_max/Q_mean at fixed N as well, which is the form
[F6](../docs/ANALYTICAL_FORMULAS.md#f6) records (D2 in the same registry): the
edge of the Q spectrum against its mean.

The block is not the XY-weight-1 Pauli sector, which is far larger (160
dimensions against 5 at N = 5), is not L-invariant, and has no spectrum of its
own; D10 Step 6 carries that scope.

### The exact formula

The block's ω_max for the Heisenberg chain, computed at each N:

| N | ω_max | exact form | V(N) = ω_max/ω_max(2) |
|:--|:------|:-----------|:-----------------|
| 2 | 4.0000 | 4J | 1.000 |
| 3 | 6.0000 | 4J + 4J·cos(π/3) = 6J | 1.500 |
| 4 | 6.8284 | 4J + 4J·cos(π/4) = (4 + 2√2)J | 1.707 |
| 5 | 7.2361 | 4J + 4J·cos(π/5) = (5 + √5)J | **1.809** |
| 6 | 7.4641 | 4J + 4J·cos(π/6) = (4 + 2√3)J | 1.866 |

The eigensolver on the N × N block reproduces the closed form to within
5.3e-15 at every row. The pattern is

    ω_max(N) = 4J·(1 + cos(π/N)) = 8J·cos²(π/(2N))
    V(N)     = 1 + cos(π/N)       = 2·cos²(π/(2N))

and it is not a pattern but a theorem: ω_max is the k = N−1 member of F2's
dispersion ω_k = 4J·(1 − cos(πk/N)), derived from the block's tight-binding
reduction in D10 for every N. The block's generator is 2J times the chain's
graph Laplacian, and the Laplacian is the one-magnon Hamiltonian measured from
the ferromagnetic vacuum, so V(N) is a ratio of single-magnon energies.

For N = 5 the golden ratio appears: cos(π/5) = φ/2 with φ = (1 + √5)/2, so

    V(5) = 1 + cos(36°) = 1 + φ/2 = (5 + √5)/4 ≈ 1.80902

For N → ∞, cos(π/N) → 1 and V → 2: the gain saturates at exactly 2 for an
infinite chain.

The value 1.809 is not a topological invariant. It is not an integer, it
depends on N and on the coupling (Heisenberg here), and it is a smooth function
of 1/N approaching 2. It is a geometric constant of the Heisenberg chain's
one-magnon band. Its angles π/(2N) are read as wave angles in
[Off-Niven as wave breaking](../docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md), and
the block's whole spectrum, not only its edge, is written out in
[Analytical Spectrum](ANALYTICAL_SPECTRUM.md).

### Uniform γ only

The cancellation holds for uniform dephasing and for that alone. Under a
non-uniform profile the dephasing stops being a scalar on the block, diag(γ)
stops commuting with the block's Laplacian, and both halves of the argument go
at once. The decay rates spread
([Cavity Mode Localization](CAVITY_MODE_LOCALIZATION.md)), and the frequencies
themselves move, level by level, growing with the profile's unevenness
([Concentrator Optics](CONCENTRATOR_OPTICS.md), Result 3). At N = 5 the block's
ω_max is 7.236016 at the edge profile [0.05, 0.01, 0.01, 0.01, 0.01],
7.234910 at [0.2, 0.01, …] and 7.228595 at [0.5, 0.01, …], against the uniform
7.236068. The best-Q mode is that same shifted mode in each case. So the 1.81
ratio does not survive by retreating to the extremal mode: ω_max moves too, and
V(N) = 1 + cos(π/N) is a uniform-γ statement.

Zero-temperature amplitude damping (n̄ = 0) does not break the ratio either.
σ⁻ annihilates the vacuum, so on each |0⟩⟨j| it only adds the decay γ_amp/2 of
the one excited site, again a scalar on the block: with γ_amp = 0.1 and no
Z-dephasing, Q_max(5)/Q_max(2) = 144.7/80.0 = 1.81, and with Z-dephasing 0.1
plus γ_amp = 0.05 it is 32.2/17.8 = 1.81.

---

## Dephasing lifts degeneracies

At zero noise the N = 5 chain has 43 distinct frequencies, the distinct
nonzero |E_i − E_j| of H in four-decimal bins (the
[Cavity Modes Formula](CAVITY_MODES_FORMULA.md) count). As uniform
Z-dephasing rises:

| γ/J | N = 5 frequencies |
|:----|:------------------|
| 0.001 | 50 |
| 0.01 | 78 |
| 0.1 | 111 |
| 0.15 | **112** |
| 0.3 | **112** |
| 1.0 | 109 |
| 5.0 | 103 |

Dephasing splits modes that shared a frequency at γ = 0. The count is largest
around γ/J ≈ 0.15 to 0.3, at 112, up from 43; above that it falls again as the
modes broaden. The bin width is a choice, and a coarser bin lowers every count;
the rows belong to this protocol. N = 2 stays at one or two frequencies
throughout, so the new frequencies are the coupled chain's.

The palindromic pairing stays exact at every γ, since Z-dephasing is the
channel F1 covers.

---

## A warmer bath trades Q for frequencies

Adding the finite-occupation channel on top of Z-dephasing (γ_z = 0.1,
γ_amp = 0.05):

| channels | n̄ | Q_max (N=2) | Q_max (N=5) | ratio | frequencies, N=5 |
|---|---:|---:|---:|---:|---:|
| Z-dephasing γ_z = 0.1 | 0 | 20.0 | 36.2 | 1.81 | 111 |
| Z-dephasing 0.1 + amplitude 0.05 | 0 | 17.8 | 32.2 | 1.81 | 111 |
| same combined channel | 0.5 | 14.5 | 21.0 | 1.44 | 403 |
| same combined channel | 2.0 | 9.4 | 12.5 | 1.33 | 423 |
| same combined channel | 5.0 | 5.5 | 7.1 | 1.29 | 445 |

The ratio leaves 1.81 at the first nonzero occupation: σ⁺ lifts population out
of the vacuum, so the (0,1) block is no longer closed under the generator and
its scalar decay no longer sets Q_max. The frequency count roughly quadruples,
from 111 to 445, while Q_max at N = 5 falls from 32.2 to 7.1.
N = 2 stays at two frequencies at every n̄. These rows record how this finite
spectrum and this binning respond to the channel parameter; they do not show
energy being converted into frequencies, and they do not define a temperature.
The registry keeps this section as [F20](../docs/ANALYTICAL_FORMULAS.md#f20).
The same channel at a named bath occupation runs in
[DNA Base Pairing](DNA_BASE_PAIRING.md), and
[Energy Partition](../hypotheses/ENERGY_PARTITION.md) reads this census beside its
own N = 3 protocol.

## The sacrifice profile's advantage fades as n̄ grows

The [sacrifice profile](RESONANT_RETURN.md) concentrates the dephasing on one
edge qubit (γ_z = 0.5 there, 0.01 on the interior four), compared here with the
uniform profile of equal total, [0.108]·5, both with γ_amp = 0.05 at N = 5:

| n̄ | Q_max sacrifice | Q_max uniform | Q ratio | frequencies sacrifice | frequencies uniform | excess |
|:--|:----|:----|:----|:----|:----|:----|
| 0.00 | 89.3 | 30.0 | **2.97** | 120 | 111 | 8.1% |
| 0.05 | 77.6 | 28.6 | 2.71 | 374 | 313 | 19.5% |
| 0.20 | 55.7 | 25.0 | 2.23 | 439 | 377 | 16.4% |
| 0.50 | 35.8 | 20.0 | 1.79 | 456 | 406 | 12.3% |
| 1.00 | 22.6 | 16.1 | 1.40 | 469 | 420 | 11.7% |
| 5.00 | 7.5 | 7.0 | 1.08 | 479 | 443 | 8.1% |
| 10.0 | 4.2 | 4.1 | **1.02** | 488 | 462 | 5.6% |

At n̄ = 0 the sacrifice profile holds three times the uniform Q_max; by
n̄ = 10 the advantage is gone, the bath's rates having overtaken the profile's
spatial structure. The frequency count tells a different story: the sacrifice
profile carries more frequencies than the uniform one in all eleven sampled
rows (the table shows seven; all eleven are in the results file), by 5.6% to 19.5%, largest at small n̄ (19.5% at n̄ = 0.05) and smallest
at the hottest row. The profile's asymmetry lifts degeneracies the bath alone
does not.

This is a within-chain comparison, one profile against another on the same
chain at the same total γ. Between chains the total noise level has to be
accounted for as well ([Chain Selection Test](CHAIN_SELECTION_TEST.md)). It is
a sweep of two specified generators and identifies n̄ with no physical
temperature.

---

## The palindrome under the bath

F1's palindromic pairing is proven for Z-dephasing
([Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)). The bath
brings σ⁻ and σ⁺, which lie outside that proof, and the pairing survives them
at a shifted centre:

| channel | palindrome centre | status |
|:--------|:------------------|:-------|
| Z-dephasing | −Σγ | F1, proven |
| amplitude damping, n̄ = 0 | −Σγ/2 | [F137](../docs/ANALYTICAL_FORMULAS.md#f137) |
| thermal bath, σ⁻ and σ⁺ | −Σ(γ↓ + γ↑)/2 | F137 extended |
| Z-dephasing **and** amplitude damping on the same axis | none exists | breaks |

Dephasing pays the full shift, amplitude damping half of it, and the thermal
bath half of the **total** per-site rate γ↓ + γ↑ = γ_amp·(2n̄ + 1).

The centre is not estimated; it is read off the trace. If a multiset is closed
under λ ↦ 2c − λ, every pair sums to 2c, so Σλ = n·c and c = mean(λ) exactly.
There is one candidate centre, and the commutator part of L is traceless, so it
never depends on H. That turns "is it palindromic?" into a check: a large
distance at that centre means no centre works. The
[Label Map](../docs/quantum/THE_LABEL_MAP.md) records this as the answer to
what "the centre" means for a palindrome.

Where the rates are rational the characteristic polynomial is exact, and the
palindrome is the identity p(2c − x) ≡ p(x). It holds with the Heisenberg H at
N = 2 and N = 3, no eigensolver anywhere, and the same test returns "neither"
for N = 2 with the thermal bath beside Z-dephasing, so that break is proven
there as well. At N = 2 to 5, with independently drawn per-site rates, the
canonical F1 distance at the trace centre sits at the eigensolver floor
([`thermal_palindrome_centre.py`](../simulations/thermal_palindrome_centre.py)).
The scorer is the one [Concentrator Mapping](CONCENTRATOR_MAPPING.md) made
canonical for every producer, read at the trace centre rather than at the
Z-dephasing centre.
In the spectral census the pure amplitude channel reads 51 to 113 in units of
machine epsilon times the spectral radius at N = 5 across all eleven
occupations, the floor there, while co-axial Z-dephasing plus amplitude reads
1.7e14 to 2.4e15: not paired.

The break is sharper than "two channels are worse than one". At N = 4 the
thermal bath beside transverse dephasing (X or Y) keeps the pairing, and only
co-axial Z breaks it (4.4e14 against 27 and 26); MIRROR_SYMMETRY_PROOF carries
the composition rule. The three additions carry the same total rate, so all
three share one centre to every digit, and only the distance tells them apart:
the centre must be computed **and** the distance checked; neither decides
alone.

Where this leaves temperature: the palindrome's existence does not depend on n̄,
only the position of its centre does. There is no critical temperature to look
for. And the hardware operating point is not the hot one. A flown two-leg
protocol on ibm_kingston, repeated on ibm_marrakesh
(`fw.Confirmations.lookup('f84_heating_leg_attribution_kingston_july2026')`),
found the bath cold: flat up-legs putting the thermal populations at 0.23% to
0.83% on the three Kingston qubits and at most 0.5% on Marrakesh, with γ↑
pinned near zero. The regime that matters on real machines is T1 beside
co-axial Z, the one combination that fails.

Still open: a proof of the H ≠ 0 pairing at general N.

---

## What the second producer establishes

For the selected local rates, one qubit has stationary excited-state
population

```text
p_excited = n̄ / (2n̄ + 1)
```

The tensor product of that local state is stationary for the tested uniform
Heisenberg chains, because it commutes with the Hamiltonian and each local
channel fixes it. The audit measures ‖L ρ_target‖ from 0.0 to 5.56e-17 at
N = 3 and N = 5, including one non-uniform Z profile. A wrong one-site
population is the negative control and gives residual 0.141.

This local-channel target is not, in general, the Gibbs state of the
interacting Heisenberg Hamiltonian, and comparing their energies does not close
a feedback equation.

## The self-consistent occupation does not exist

A self-heating fixed point would be an occupation n̄ at which the channel's
steady energy equals the thermal energy that same n̄ names, so the loop closes
on itself. The loop is a model's closing rule, not part of the generator, which
takes n̄ as an input; the registry keeps what follows as
[F21](../docs/ANALYTICAL_FORMULAS.md#f21). The fixed point does not exist, and this needs no sweep: both sides are known in
closed form on this branch.

The steady state is the product of the local targets diag(1 − p, p) with
p = n̄/(2n̄ + 1), the object the stationarity gate above already pins, so its
energy is

    Tr(H ρ_steady) = (N − 1)·⟨Z⟩² = (N − 1) / (2n̄ + 1)² ≥ 0

which the producer reproduces in every one of its 27 printed rows, three
configurations at nine occupations each, with the worst residual 7.55e-15.
Meanwhile tr H = 0 exactly, so Tr(H ρ_β) = 0 at β = 0, and since
d/dβ Tr(H ρ_β) = −Var_β(H) < 0 for this non-constant H, any β > 0 gives
Tr(H ρ_Gibbs) < 0 strictly. The gap is a non-negative number minus a negative
one and cannot vanish, at any N, any rate, and any map from n̄ to a temperature
with β > 0.

That is stronger than a sweep, and it is also what a sweep could never have
shown: a scan reporting "no sign change over five decades" is reporting an
identity it cannot see. The two energies belong to different objects, one a
local emission/absorption target and one a Gibbs state of the interacting
chain, and no temperature makes a state of one into a state of the other.

---

## Open question

Does any other joint-popcount block of the Liouvillian reach a higher Q than
the (0,1) coherence block? The V(N) argument reads Q_max off that block alone,
and D10 proves the block's ω_max for every N; the step from "best on the block"
to "best in the Liouvillian" is measured, not proven. At γ = 0.05 the full
Liouvillian's Q_max equals the block's at N = 2 to 6, 40, 60, 68.284271,
72.360680 and 74.641016, carried by the block's own eigenvalue −2γ ± i·ω_max
(the block is closed under L, so that eigenvalue is exactly one of L's); the
19-row γ sweep agrees at N = 2, 3 and 5 at every rate to its printed precision.
Outside the block and its X^⊗N image (N−1, N), which carries the same
spectrum and so ties it exactly, the nearest competitor stands at 0.9997, 0.750,
0.693, 0.809 and 0.871 of the block's Q_max at N = 2 to 6, a margin that narrows from N = 4
on, so the measurement cannot stand in for the general-N proof that is
missing. It is the one step the V(N) reading needs to speak for the whole
spectrum. (The thermal and sacrifice-profile Q_max
values above run under amplitude damping or a γ profile, where neither the
block's closed form nor this question applies as stated.)

---

## Reproduction

```bash
python simulations/v_effect_gamma_sweep.py
python simulations/v_effect_thermal.py
python simulations/thermal_palindrome_centre.py
python simulations/self_heating_fixpoint.py
```

The V(N) section reads the uniform-γ rows of the first; the thermal and
sacrifice tables, Part 2 and Part 4 of the second.
