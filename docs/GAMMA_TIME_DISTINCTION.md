# The Three Levels of Time: Dephasing Removes the Return

<!-- Keywords: gamma dephasing three levels of time, parameter time vs
observable change vs direction, recurrence removed by dephasing, Bell+
eigenstate stationary, |01> recurrence period pi/(2J), tau=gamma*t
homogeneity, collapse at fixed Q not at fixed J, J provides content gamma
provides direction reading, R=CPsi2 gamma time -->

**Status:** Computed trajectory comparison at N = 2 (Tier 2); the reading of γ as the direction of time is a reading (Tier 5)
**Date:** March 22, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Scripts:** [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py) ([output](../simulations/results/gamma_is_time_proof.txt)), [two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py), [gamma_unit_scaling_gate.py](../simulations/gamma_unit_scaling_gate.py), [disprove_gamma_is_time.py](../simulations/disprove_gamma_is_time.py)

---

## What this document is about

Think about what "time" means in everyday life. Not the number on a
clock face, but the feeling that things move forward: that milk spills
and does not unspill, that you age and do not un-age, that decisions,
once made, stay made. Physicists call this the "arrow of time," the
difference between past and future.

We asked what, in a small quantum system, looks like that arrow, and
the word "time" turned out to hide three different things. The
mathematical parameter t (just a number in an equation) is there in
every model, noisy or not. Oscillation (things swinging back and
forth, like a pendulum) looks like change but comes back to where it
started. And dephasing noise (γ) does something neither of the other
two does: it takes the return away.

The experiment is simple and a little dramatic. Take two qubits with
no noise: depending on how you prepare them, they either sit perfectly
still or swing in circles, returning exactly to where they began, again
and again. Turn on a small amount of noise, and the returns stop; the
state settles and does not come back. That much is measured. That this
is what gives time its direction is how we read it, and the page keeps
the two apart.

---

## Abstract

Time has three levels here: (1) the formal parameter t in d/dt, present
with or without γ; (2) observable change, which without γ is either
stillness or exact recurrence; and (3) direction, before and after,
which is where we read γ in. On the N = 2 Heisenberg chain at γ = 0,
Bell+ does not move at all over t = 0 to 50, because it is an
eigenstate of the Hamiltonian, and |01⟩ returns to itself with period
π/(2J), its CΨ crossing ¼ 63 times downward and 64 times upward. At
γ = 0.05 the returns are gone, purity falls from 1 to ½, and the state
comes to rest. The generator is homogeneous, L(J, γ)·t = τ·L(J/γ, 1)
with τ = γt, so every trajectory is a function of τ **and** of
Q = J/γ: holding Q, the curves collapse onto τ exactly; holding J, they
do not, because that sweep moves Q. A trajectory here needs two
numbers, and γ supplies one of them. We read this as γ giving the
direction and J giving the content.

---

## The Three Levels

The word "time" hides three different things. This table separates
them.

| Level | Without γ | With γ |
|-------|-----------|--------|
| t as the symbol in d/dt | Yes (the parameter is always there) | Yes |
| t as observable change | Stillness or exact recurrence | Damping; the returns stop |
| t as direction (before and after) | No direction in the tested cases: nothing accumulates | The reading: γ supplies the direction |

γ has no bearing on Level 1: Hamiltonian evolution runs in t at γ = 0
just as it does at γ > 0. The difference between the columns lives at
Level 2, where it is measured, and our reading carries it up to
Level 3. Even there γ is not the whole story: γ sets the scale on which
the returns stop, J sets what happens along the way, and Part 3 below
shows that a trajectory needs both.

---

## The Evidence

The following runs show what "time" looks like with and without noise
on two qubits coupled by a Heisenberg bond (J = 1) under local Z
dephasing. Bell+ is a maximally entangled state (two qubits perfectly
correlated). |01⟩ is a simple product state (one qubit up, one qubit
down). CΨ is our composite diagnostic of how quantum the state is (see
[What We Found](WHAT_WE_FOUND.md)).

### Bell+ at γ = 0: the parameter runs, and nothing moves

The parameter t goes from 0 to 50. Every observable is constant:
CΨ = 0.3333, concurrence = 1.0, trace distance (the standard measure of
how distinguishable two quantum states are) = 0.0000 at every time
step. From inside, t = 0 and t = 50 are the same state. The reason is
the preparation, not the absence of noise: Bell+ is an eigenstate of
the Heisenberg Hamiltonian, the commutator [H, ρ_Bell+] vanishes in
every entry, and a stationary state reports on itself rather than on
whether anything could move. The next run moves plenty.

### |01⟩ at γ = 0: the parameter runs in circles

|01⟩ is not an eigenstate, and it swings. The Hamiltonian's spectrum is
{−3J, J}, a single gap of 4J, so the state returns to itself with
period π/(2J), exactly, forever. Over t = 0 to 50 CΨ crosses ¼ 63 times
downward and 64 times upward. (On the script's 0.1 grid the first
sampled return falls at t = 11.0, which is seven periods; the returns
before it fall between grid points, and the grid, not the period,
decides which one it samples first.) There is
change, but no direction: no observable accumulates, and the ¼ line is
crossed freely in both directions. The γ = 0 Liouvillian spectrum is
purely imaginary (0 ten times and ±4i three times each) and pairs 16 of 16 about its own centre, zero.

### |01⟩ at γ = 0.05: the returns stop

Switch on γ = 0.05 and the trace distance never comes back below 0.01.
Purity falls from 1.0 to 0.5000. The rate of change ‖dρ/dt‖ at t = 50
drops from 2.81 at γ = 0 to 0.019: the state has all but stopped. The
spectrum acquires real parts, −0.1 and −0.2, the decay rates. CΨ still
crosses ¼ twice in each direction before it settles, because the
Hamiltonian's swing rides on the decay; Bell+ under the same noise
crosses ¼ once, downward, and stays below for the rest of the window.
Either way, what the dephasing removes is the recurrence. Things happen
that do not unhappen. This is where we read a direction in: past and
future become distinguishable by how far the decay has gone.

---

## The Tests Reread

This page began as an attempt to disprove "γ is time"
([disprove_gamma_is_time.py](../simulations/disprove_gamma_is_time.py)).
Its numbers hold, and read with the distinction below they say more
than "no". One distinction runs through all of them: a frequency
(swinging back and forth) is not a direction (counting forward and
never coming back).

**Test 1 (γ = 0 evolution):** The parameter t exists without noise, and
Hamiltonian evolution uses it. Bell+ at γ = 0 shows no observable
change, because it is stationary; |01⟩ changes and recurs. Nothing
accumulates.

**Test 2 (same τ, different t):** At the same τ = γt, two systems carry
different Hamiltonian phases J·t. That phase is periodic, a rotation
and not an accumulation. It is the second number of Part 3 below, seen
from the other side.

**Test 3 (the Hamiltonian's frequency):** The peak frequency of the
σ_x ⊗ σ_y signal is 0.5994 at every γ from 0 to 0.10; only its
amplitude falls, from 1.0000 to 0.9267. The Hamiltonian has its own
frequency, independent of the noise, and a frequency comes back.

**Test 4 (several γ at once):** With a γ profile [0.01, 0.10, 0.01] on
three qubits, all sites share one t. That is because t is a coordinate
of the whole equation, not a property of a site. What the local γ sets
is the local decay scale.

**Test 6 (Π reversal):** Π maps the Liouvillian onto its mirror image
exactly (‖ΠLΠ⁻¹ + L + 2Σγ·I‖ = 0), yet applying Π at t = 5 and evolving
on leaves the state at t = 10 a trace distance 0.388 from where it
began. Π reverses the spectrum's structure about its centre, not the
decay. The irreversible part is not undone.

---

## The Computation in Three Parts

Script: [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py).
Data: [gamma_is_time_proof.txt](../simulations/results/gamma_is_time_proof.txt).
Part 1 asks what switching γ on changes. Part 2 lays the γ-only and the
J-only corner side by side. Part 3 tests whether τ = γt carries the
whole trajectory, and finds exactly when it does.

### Part 1: What switching γ on changes

| Property | \|01⟩ γ=0 | \|01⟩ γ=0.05 | Bell+ γ=0 | Bell+ γ=0.05 |
|----------|---------|------------|-----------|-------------|
| S(ρ_A) non-decreasing after t = 1 | No | No | Yes (trivially: S = 1 constant) | Yes (trivially: S = 1 constant) |
| trace distance returns (< 0.01) | Yes (recurrence) | No | Yes (trivially: always 0) | No |
| CΨ crossings of ¼, down/up | 63/64 | 2/2 | 0/0 | 1/0 |
| ‖dρ/dt‖ at t = 50 | 2.81 | 0.019 | 0.000 | 0.000006 |

Without γ: recurrence (|01⟩) or stillness (Bell+). With γ: the state
converges and stays. The entropy of |01⟩ at γ = 0.05 still does not
rise monotonically, because the Hamiltonian's oscillation modulates
the decay; pure monotonicity belongs to J = 0, and Bell+ passes only
trivially, its one-site entropy pinned at 1.

### Part 2: The γ-only corner and the J-only corner

| Configuration | S non-decreasing | D does not return | CΨ crosses one way only |
|---------------|------------------|-------------------|-------------------------|
| J=0.1, γ=0.05 | No | Yes | Yes |
| J=1.0, γ=0.05 | No | Yes | No |
| J=10, γ=0.05 | No | Yes | No |
| **J=0, γ=0.05** | **Yes** | **Yes** | **Yes** |
| J=1.0, γ=0 | No | No | No |

**J = 0, γ > 0:** pure decay, no Hamiltonian. Entropy rises
monotonically, the distance never returns, CΨ crosses once, downward.

**J > 0, γ = 0:** oscillation and recurrence. No diagnostic says yes.

On this state and this grid the diagnostics separate the pure-decay
corner from the unitary one, and every row with γ > 0 loses its return.
They are diagnostics of damping. That damping is the direction we
experience as time is the reading we bring to them, not something a
table of yes and no could establish.

### Part 3: The τ = γt collapse, and what it actually measures

The generator is homogeneous: L(J, γ)·t = τ·L(J/γ, 1) with τ = γt. So
every observable is a function of τ **and** of Q = J/γ, and a collapse
test has to be run twice or it says nothing. Spread is the largest gap
between two curves at equal τ, range is how far the observable itself
travels, and collapse is called only under 5% of the range, so an
observable cannot pass by standing still.

Holding J = 1 and sweeping γ from 0.01 to 0.20, which sweeps Q from 100
to 5:

| Observable | spread | range | spread/range | collapses? |
|---|---|---|---|---|
| S(ρ_A) | 0.789752 | 1.000000 | 0.7898 | no |
| Tr(ρ²) | 0.057319 | 0.491263 | 0.1167 | no |
| CΨ | 0.259139 | 0.307348 | 0.8431 | no |
| Concurrence | 0.861100 | 0.959735 | 0.8972 | no |

Holding Q = 20 instead and sweeping the same γ values, with J moving
along:

| Observable | spread | range | spread/range | collapses? |
|---|---|---|---|---|
| S(ρ_A) | 0.000000 | 0.999986 | 0.0000 | yes |
| Tr(ρ²) | 0.000000 | 0.490770 | 0.0000 | yes |
| CΨ | 0.000000 | 0.307348 | 0.0000 | yes |
| Concurrence | 0.000000 | 0.959735 | 0.0000 | yes |

Both arms come out as the identity requires, and the producer raises
if either one does not. At fixed Q the curves land on each other to
machine precision: τ = γt **is** this generator's own time. At fixed J
they do not, and the reason is not that τ is the wrong time variable.
Holding J while sweeping γ moves Q, so those five runs are five
different systems compared at matched τ, and what the failure measures
is the second knob. The generator's scaling symmetry is **joint**,
L(λJ, λγ) = λ·L(J, γ) ([Q Scale Three Bands](../experiments/Q_SCALE_THREE_BANDS.md),
Tier 1), and the left table alone, read as "irreversible observables do
not scale with τ", mistakes a change of system for a failure of the
time variable.

One observable is left out of these tables on purpose. arg(ρ₀₁) looks
like a phase and is not: ρ₀₁ links popcount 0 to popcount 1, and both
the Hamiltonian and the Z-dephasing conserve popcount, so for this
preparation the element is identically zero at every γ and every t.
Its neighbour arg(ρ₁₂), the |01⟩⟨10| coherence this state does
populate, is no better as a row: it takes three values, so its spread
equals its range by construction and it could only ever print a
failure. Both are observables whose verdict is fixed before the run.

### Why Bell+ collapses at fixed J anyway: the state

[F14](ANALYTICAL_FORMULAS.md) reports K = γ·t_cross constant across a γ
sweep at fixed J ([Crossing Taxonomy](../experiments/CROSSING_TAXONOMY.md)),
and the left table above reports the trajectory not collapsing in τ.
Both are right, and what separates them is the initial state. F14 runs
**Bell+**, the eigenstate: Z-dephasing keeps it in its sector, J never
enters its trajectory, and its observables really are functions of τ
alone, with K = ln(4/3)/8 in the concurrence book. Part 3 runs **|01⟩**,
which is not an eigenstate. The mechanism is measured rather than
argued: the commutator [H, ρ_Bell+] is zero in every entry, against a
largest magnitude of 2.0 for |01⟩, and at fixed J the purity curves at
matched τ collapse to 6.7·10⁻¹⁶ for Bell+ against 4.6·10⁻² for |01⟩
([`gamma_unit_scaling_gate.py`](../simulations/gamma_unit_scaling_gate.py),
part F). Where the state gives Q nothing to act on, moving Q moves
nothing.

---

## What the Two Knobs Say Together

Read together, the three parts say something narrow and firm: a
trajectory here needs two numbers, γ and J, and γ supplies one of them.
γ sets the scale on which the returns stop; J sets what swings while
they do. A quantity that sets the time variable only alongside a second
knob is not what "γ is time" would claim, so that sentence does not
survive its own tables. What survives is a division of labour, and our
reading of it: **γ gives the direction, J gives the content, and
neither alone makes a trajectory.**

This is also why the repository keeps the word "clock" for a specific
object: the two-handed clock of [F95](ANALYTICAL_FORMULAS.md#f95),
one hand set by γ and the other by J. This page meets the same
two-handedness from the trajectory side, as the pair (τ, Q); that the
two are one object is an adjacency we note here, not a result this page
derives.

## Where γ comes from

The [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md) certifies
that a generator with a nonzero dissipative part is open. It does not
say what is outside: a source of decay written as a Lindblad dissipator
already assumes an environment, and one written unitarily is closed and
carries no decay, so inside the formalism the origin question cannot be
settled. That γ comes from an outside is a reading; whether any part of it could come
from within stays open. Nothing on this page depends on the answer.

---

## Precise Language

| Statement | Status |
|-----------|--------|
| The formal parameter t exists without γ | Correct (Level 1): Hamiltonian evolution runs in t at γ = 0 |
| γ sets a dissipative timescale in the generator | Established by the model definition |
| At γ = 0, the two tested states are stationary (Bell+, an eigenstate) or exactly recurrent (\|01⟩, period π/(2J)) | Measured at N = 2 |
| With γ > 0, the tested trajectories stop returning and come to rest | Measured at N = 2 |
| The Hamiltonian provides a frequency | Correct (0.5994, independent of γ); a frequency comes back |
| τ = γt alone determines the trajectory at fixed J | False; fails on all four observables of Part 3 |
| τ = γt determines the trajectory at fixed Q = J/γ | True; collapses exactly |
| γ gives time its direction, J its content | Our reading (Level 3); the computation establishes the decay scale |
| γ is the condition for experienced time | A reading, not established; nothing here measures experience |
| A nonzero dissipator identifies its microscopic source | Not established ([Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md)) |
| The ¼ line is a one-way door | Not in general: \|01⟩ at γ = 0.05 crosses 2/2; Bell+ crosses once, downward |

---

## Reproduction

```bash
python simulations/gamma_is_time_proof.py
python simulations/two_qubits_no_noise.py
python simulations/gamma_unit_scaling_gate.py
```

Every number above belongs to its stated system (N = 2, Heisenberg,
local Z dephasing), preparation, time window, and observable set.

---

## References

- [Incompleteness Proof](proofs/INCOMPLETENESS_PROOF.md): the system is certified open; the origin of γ cannot be settled inside the formalism
- [The Bridge Was Always Open](THE_BRIDGE_WAS_ALWAYS_OPEN.md): the γ-as-time reading and how far it reaches
- [Q Belongs to No Substance](Q_BELONGS_TO_NO_SUBSTANCE.md): the same scaling identity, L(J, γ) = γ·L₁(Q)
- [two_qubits_no_noise.py](../simulations/two_qubits_no_noise.py): what time looks like without γ
- [disprove_gamma_is_time.py](../simulations/disprove_gamma_is_time.py): the original tests, reread above
- [gamma_is_time_proof.py](../simulations/gamma_is_time_proof.py): the three parts (what γ changes, the two corners, the two-armed τ collapse)
- [gamma_unit_scaling_gate.py](../simulations/gamma_unit_scaling_gate.py): the joint scaling and the Bell+ / |01⟩ control
- [gamma_is_time_proof.txt](../simulations/results/gamma_is_time_proof.txt): raw results
