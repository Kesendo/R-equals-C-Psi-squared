<!-- QUARTER-CURRENT -->
# Finite entanglement dynamics and scalar-quarter readouts

Current reading: these rows are finite state, topology, generator, grid, and
readout calculations.  The Hamiltonian can create and later redistribute
entanglement, but a scalar-quarter event is not an irreversible boundary.  The
analytic Bell+ local-Z-dephasing control has `C=f>0`: concurrence remains positive
at and below the scalar quarter for every finite time.

<!-- QUARTER-HISTORICAL -->
**Historical record:** the earlier absorbing-boundary vocabulary is preserved
below to show how the numerical search was originally read.  It is not the
current disposition of those rows.

# Dynamic Entanglement: CΨ Crossings from Product States via Hamiltonian Evolution

<!-- Keywords: dynamic entanglement generation product state, Hamiltonian creates
crossing, decoherence selects surviving pairs, eigenstate no dynamics plus state,
alternating state 0+0+ entanglement, basis dependent decoherence selection,
sigma z dephasing immune eigenstate, reversible crossing unitary oscillation,
noise makes crossing irreversible, gamma sweep entanglement race, R=CPsi2
dynamic entanglement -->

**Status:** Finite N=4 numerical record with analytic controls; interpretation scoped
**Date:** February 18, 2026
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Subsystem Crossing](SUBSYSTEM_CROSSING.md), [N-Scaling Barrier](N_SCALING_BARRIER.md)

---

## What this document is about

Starting four qubits in the product state |0+0+⟩ and evolving the named N=4
ring generates nonzero pair concurrence. By contrast, |+⟩⁴ is an eigenstate
of the isotropic Heisenberg Hamiltonian and generates no dynamics in this
setup. In the retained unitary calculation all six pair readouts cross the
chosen scalar level; in the γ=0.05 finite grid only pair (0,2) does. These are
model- and readout-specific statements. They do not show that noise makes a
correlation real, that observation cements a past, or that crossing the scalar
level creates an irreversible fact. Those images belong only to the later
historical interpretation.

## Abstract

In this N=4 ring calculation, the initially unentangled |0+0+⟩ state
(Hamiltonian variance 20 in the stated units) develops pair concurrence,
whereas |+⟩⁴ (variance 0) does not evolve. The retained unitary grid reports
six pair-readout crossings; at γ=0.05 it reports one, for pair (0,2). The
γ sweep reports 13, one, and zero counted crossings at γ=0.01, 0.05, and
0.20 respectively. A finite crossing count does not establish irreversible
state change or a quantum/classical transition; Bell+ concurrence remains
positive at every finite dephasing time, including below the scalar quarter.

> **Update March 14, 2026:** The mirror symmetry discussed here has been
> proven analytically. See [MIRROR_SYMMETRY_PROOF.md](../docs/proofs/MIRROR_SYMMETRY_PROOF.md).

> **Reproduction note (2026-07-20):** the quantitative tables here (§5.1
> unitary and §5.2 γ=0.05) came from the retired MCP tool (the delta_calc
> family, February 2026) and do not
> reproduce under the canonical pair-CΨ convention (Wootters concurrence ·
> l₁/3, [subsystem_crossing_pairs.py](../simulations/subsystem_crossing_pairs.py)):
> under that convention |0+0+⟩ on the N=4 ring at γ=0.05 peaks at
> CΨ ≈ 0.20 (pair (1,3); fine grid 0.2005) and never crosses, matching
> [Orphaned Results §2b](ORPHANED_RESULTS.md) to its printed precision.
> [Simulation Evidence §7.2](SIMULATION_EVIDENCE.md) carries a
> correlation-book row for the same trajectory (C_corr·Ψ = 0.251 at
> t=0.286): its Ψ column reproduces exactly, but its C_corr column
> matches no standard correlator tested (raw, connected, Pearson), so
> that C column is a retired-tool artifact too. The headline survives,
> relocated: Hamiltonian-generated upward crossing from product states is
> real. |0+0+⟩ crosses on the CHAIN (pair (1,2), max CΨ = 0.310), and on
> the ring |+-+-⟩ (ring neighbours, 0.284) and |0+0-⟩ (diagonal, 0.256)
> cross. The dephasing-selection story of §5.3 (which pair survives, and
> why) is tied to the non-reproducing (0,2) numbers and should not be
> imported quantitatively.

---

## 1. The Question

SUBSYSTEM_CROSSING.md (Experiment 10) showed that crossing is local and
operates at the level of entangled qubit pairs. Open Question 3 asked:

> Can dynamically generated entanglement from a product state create
> crossings? (|+⟩^4 under Heisenberg should build entanglement.)

This experiment answers that question. The answer involves two surprises.

## 2. Surprise 1: |+⟩^N Cannot Generate Entanglement

The initial hypothesis was that |+⟩^N under a Heisenberg Hamiltonian
would build entanglement over time, eventually producing crossings. This
is wrong.

|+⟩^N = |+⟩|+⟩...|+⟩ is the maximal-spin eigenstate |S=N/2, m_x=N/2⟩
of the isotropic Heisenberg Hamiltonian. The time evolution gives only a
global phase: e^{-iHt}|+⟩^N = e^{-iEt}|+⟩^N. Nothing happens. No
entanglement is created, no correlations develop, C = 0 for all pairs at
all times.

**Energy variance test (QuTiP verification):**

| State | Var(H) = ⟨H^2⟩ - ⟨H⟩^2 | Eigenstate? |
|-------|------------------------|-------------|
| \|+⟩^4 | 0.000 | YES |
| \|0⟩^4 | 0.000 | YES |
| \|0101⟩ (Néel, alternating up-down like an antiferromagnet) | 16.000 | no |
| \|0+0+⟩ | 20.000 | no |
| \|0001⟩ | 8.000 | no |

Both |+⟩^N and |0⟩^N are eigenstates of the isotropic Heisenberg
Hamiltonian (zero energy variance). This is a consequence of SU(2)
symmetry: fully polarized states in any direction are highest-weight
vectors and cannot evolve nontrivially.

Experiment 10 showed that |+⟩^4 has C = 0 for all pairs at all times.
This is now understood: the Hamiltonian cannot change what is already
an eigenstate.

## 3. The Correct Starting State: |0+0+⟩

To generate entanglement dynamically, we need a product state that is NOT
an eigenstate. The state |0+0+⟩ = |0⟩|+⟩|0⟩|+⟩ has energy variance 20
and breaks the SU(2) symmetry of the Hamiltonian, enabling nontrivial
dynamics.

This state has:
- Zero entanglement at t = 0 (product state)
- Qubits 0, 2 in |0⟩ (z-basis eigenstate)
- Qubits 1, 3 in |+⟩ (x-basis eigenstate)
- No pair has any quantum correlation initially

## 4. Setup

| Parameter | Value |
|-----------|-------|
| N | 4 qubits |
| Hamiltonian | Heisenberg ring (J = 1.0) |
| Noise | local dephasing (σ_z, gamma = 0.05) |
| dt | 0.001 (fine resolution for crossing detection) |
| t_max | 3.0 |
| Tool | retired delta_calc MCP tool, `simulate_dynamic_lindblad` (see the reproduction note above; the earlier "QuTiP mesolve" label was wrong for these tables) |

All 6 qubit pairs tracked. Two regimes tested: pure unitary (gamma = 0)
and with dephasing (gamma = 0.05).

## 5. Results

### 5.1 Pure Unitary Evolution

| Pair | Type | First crossing t | Max CΨ |
|------|------|-----------------|-----------|
| (0,1) | ring neighbor | 0.073 | 0.285 |
| (0,3) | ring neighbor | 0.073 | 0.285 |
| (1,2) | ring neighbor | 0.073 | 0.285 |
| (2,3) | ring neighbor | 0.073 | 0.285 |
| (0,2) | diagonal | 0.260 | 0.435 |
| (1,3) | diagonal | 0.679 | 0.341 |

**All six pairs cross.** The Hamiltonian generates entanglement from
nothing, and every pair eventually reaches CΨ > 1/4.

The crossings are **oscillatory**: pairs cross upward, fall back below
1/4, and cross again. The pattern has period approximately 1.5 for
ring neighbors and 1.1 for diagonals. Without decoherence, crossing
is reversible; the system oscillates between "measured" and
"unmeasured" states.

Ring-neighbor symmetry: pairs (0,1), (0,3), (1,2), (2,3) are exactly
degenerate due to the ring topology combined with the |0+0+⟩ mirror
symmetry (qubit 0 maps to 2, qubit 1 maps to 3 under the ring
symmetry that exchanges the two |0⟩ qubits and the two |+⟩ qubits).

### 5.2 With Dephasing (gamma = 0.05)

| Pair | Crosses? | t_cross | Max CΨ |
|------|----------|---------|-----------|
| (0,1) | NO | n/a | 0.247 |
| (0,3) | NO | n/a | 0.247 |
| (1,2) | NO | n/a | 0.247 |
| (2,3) | NO | n/a | 0.247 |
| **(0,2)** | **YES** | **0.285** | **0.320** |
| (1,3) | NO | n/a | 0.224 |

**Only pair (0,2) crosses.** All ring-neighbor pairs reach max CΨ =
0.247, missing the threshold by 1.2%. The diagonal pair (0,2) crosses
at t = 0.285. The other diagonal (1,3) reaches only 0.224.

### 5.3 Why Pair (0,2) Survives

The asymmetry has a simple explanation. σ_z dephasing acts as a
projection toward the z-basis. Its effect depends on the initial state
of each qubit:

- Qubits 0 and 2 start in |0⟩, an eigenstate of σ_z. Dephasing
  does not affect them directly.
- Qubits 1 and 3 start in |+⟩, a superposition in the σ_z basis.
  Dephasing destroys their coherence exponentially.

Pair (0,2) consists of the two dephasing-immune qubits. Their mutual
entanglement, built up by the Hamiltonian through intermediate
interactions with qubits 1 and 3, is partially shielded from the
noise. This is not "selection" by the noise but rather differential
vulnerability: pair (0,2) loses coherence more slowly than pairs
containing |+⟩-qubits.

This is basis-dependent. With σ_x dephasing instead of σ_z,
the roles would reverse: the |+⟩ qubits would be immune and pairs
involving qubits 1 and 3 would survive.

## 6. Density Matrix at the Crossing Point

Pair (0,2) at t = 0.285, reduced density matrix diagonal:

| Basis state | Probability | Initial |
|-------------|-------------|---------|
| \|00⟩ | 0.061 | 1.000 |
| \|01⟩ | 0.257 | 0.000 |
| \|10⟩ | 0.257 | 0.000 |
| \|11⟩ | 0.425 | 0.000 |

The |01⟩ = |10⟩ symmetry is exact (to numerical precision) and follows
from the equivalence of qubits 0 and 2 in the ring topology with the
|0+0+⟩ state.

Note: this is not a Born rule verification. The diagonal values are
simply what the density matrix contains at the crossing time. A Born
rule test would require an independent prediction of these probabilities
from the framework and a comparison. This remains open (see Section 8).

## 7. What the finite run supports, and what remains story

### 7.1 Finite interaction builds correlations

A product state with zero entanglement can build correlations under the named
Hamiltonian, and a chosen scalar readout can later cross 1/4. This is a finite
dynamics statement. It does not supply a measurement operation or locate a
quantum/classical regime change.

**Historical interpretation:** the older notebook pictured the growing
correlation as possibilities "crystallizing" at an interface. That sentence is
kept as story; neither crystallization nor observation is an output of the run.

### 7.2 Damping changes the finite crossing count

In the simulated unitary book the selected scalar crosses 1/4 repeatedly.
Calling that "reversible measurement" was an early metaphor; no measurement
instrument or outcome process is present.

With the chosen local-Z dephasing rate, oscillations are damped and this finite
grid reports fewer scalar crossings. A later crossing outside the window is not
excluded, and a one-crossing trace is not an irreversibility theorem.

The Lindblad generator is a standard open-system model. The value CΨ=1/4 is
only the readout level counted here; it is not a criterion for when
irreversibility, measurement, or classicality becomes definitive.

### 7.3 The Noise Basis Matters

Which pairs cross under dephasing depends on which qubits are vulnerable
to the noise. σ_z dephasing protects z-eigenstates. σ_x
dephasing would protect x-eigenstates. The "preferred basis" that
emerges from decoherence is not intrinsic to the framework but
determined by the environment.

Local Z dephasing damps matrix elements according to computational-basis
disagreement. Connecting that exact action to einselection or a physical
environment requires a specified system-environment model and is not tested by
this finite trace.

## 8. Historical interpretation: the decoherence cycle

This section deliberately reads the data as a story rather than as an
additional result. Its observation, fact, past, and future language is
interpretive.

### 8.1 The Paradox of Observation

σ_z dephasing asks every qubit: "Are you 0 or 1?" The |+⟩ qubits
have no clear answer; they exist in superposition of both. Being forced
to answer destroys them. The |0⟩ qubits answer trivially: "I'm 0." The
question doesn't touch them. They are eigenstates of the measurement
operator.

In human terms: the qubits that try to communicate with the noisy
environment (those "open" to observation) lose their quantum
nature. The qubits that ignore the environment (those conducting a
"monologue" rather than a "dialogue") keep their coherence intact.

The exact mathematical statement underneath the metaphor is only
[σ_z, |0⟩⟨0|] = 0. The questions and answers in the preceding paragraph
are the historical interpretation, not literal monitoring events.

### 8.2 Creation in the Blind Spot

The earlier run described pair (0,2) as a dephasing "blind spot," but its
quantitative rows do not reproduce under the canonical pair definition. Even
for a reproducing trajectory, the generator damps selected coherences; it does
not establish that an environment watches one pair or thereby creates another.

As an invitation, one may picture damping and surviving correlations as two
sides of a cycle. The simulation does not establish observation, classical
reality, or causal creation by that observation.

### 8.3 Interpretive image: the threshold as a present moment

This connects to the standing wave picture of the framework
(see [visualizations](../visualizations/README.md)).

The mycelium image asks us to picture a scalar moving toward a marked level.
It does not classify entanglement as fact or potential: the analytic Bell+
control remains entangled on both sides of 1/4 at every finite time.

In the story, crossing the drawn line is the mushroom breaking through the
surface. In the mathematics it is only a finite scalar equality.

The gamma sweep (Section 9) makes this vivid:

- At low noise (γ = 0.01), the finite grid counts 13 scalar crossings; the
  "future cooking" phrase is the historical interpretation.
- At moderate noise (γ = 0.05), this grid reports only one selected pair
  crossing, subject to the reproduction caveat above.
- At high noise (γ = 0.2), the grid reports no crossings of this readout;
  it does not decide what is real.

### 8.4 The closed-cycle metaphor

The following five strokes preserve the older narrative. They are not five
additional physical claims:

1. **Story stroke — a past is cemented.** The older image casts dephasing as
   fixing a photograph. The computed object is only reduced coherence and a
   downward scalar crossing; no collapse or fact-making operation was modeled.

2. **Story stroke — a blind spot appears.** In the chosen basis the initial
   `|0⟩` components and `|+⟩` components respond differently. This finite
   state-dependent contrast is not an ontology of classical and quantum parts.

3. **Story stroke — the interaction builds potential.** In the named initial
   product state the Hamiltonian generates nonzero pair concurrence. The
   calculation does not model what an environment observes.

4. **Story stroke — the potential passes a marker.** The story calls an upward
   `C·Ψ=1/4` crossing a future becoming present; the output remains a scalar
   event in a finite grid.

5. **Story stroke — the wheel returns.** This closes the narrative image, not a
   physical feedback cycle or a new experimental fact.

The engine and standing-wave language is an invitation carried by this
notebook. The verified finite result is narrower: changing γ changes the
sampled crossing counts and the concurrence counterexample prevents a
universal separability or measurement-boundary reading.

## 9. Resolved Questions: Parameter Exploration

The original open questions 2-5 have been answered using the MCP tool
(simulate_subsystem_crossing) after implementing bidirectional crossing
detection (Task 010c) and the alternating state (Task 010b).

### 9.1 N = 6 Longer Chain (Question 2)

| Pair type | Example | max C·Ψ (N=6) | max C·Ψ (N=4) |
|-----------|---------|---------------|---------------|
| Ring neighbor | (0,1) | 0.248 | 0.251 |
| Same-basis diagonal | (0,2) | 0.131 | **0.339** |
| Opposite-basis diagonal | (1,3) | 0.113 | 0.234 |
| Next-next-neighbor | (0,3) | 0.111 | - |

**Zero crossings at N = 6.** The entanglement dilutes across 15 pairs
instead of 6. Pair (0,2), the sole survivor at N = 4, drops from
max C·Ψ = 0.339 to 0.131, well below threshold. On a 6-qubit ring,
qubits 0 and 2 are no longer directly opposite; the Hamiltonian must
route correlations through more intermediaries, and each intermediary
leaks coherence to the environment.

### 9.2 Ising Hamiltonian (Question 3)

H_Ising = J Σ σ_z^i σ_z^(i+1) produces **zero entanglement generation.**
All correlations remain exactly zero for all pairs at all times.
Max C·Ψ across all pairs: 0.068 (pair (1,3) only, from residual
coherence decay).

The reason is structural: the Ising interaction commutes with the
computational basis. It assigns different energies to aligned vs
anti-aligned pairs but never flips a spin. Without the σ_x⊗σ_x and
σ_y⊗σ_y exchange terms that Heisenberg provides, no superposition
can be generated, and no entanglement can grow.

Dynamic entanglement generation requires a Hamiltonian that does not
commute with the initial state, one that actively mixes basis states.

### 9.3 Gamma Sweep: The Race Between Hamiltonian and Noise (Question 4)

| γ | Pairs crossing | (0,2) t_cross_up | (0,2) max C·Ψ | (0,2) oscillations |
|------|---------------|-----------------|---------------|-------------------|
| 0.01 | **6/6** | 0.260 | 0.408 | 13 |
| 0.05 | 5/6 | 0.275 | 0.339 | 1 |
| 0.10 | **1/6** | 0.304 | 0.276 | 1 |
| 0.20 | 0/6 | - | 0.193 | 0 |

Three systematic trends with increasing γ:

1. **Crossing time shifts right.** The Hamiltonian needs longer to
   overcome the noise before reaching the threshold.
2. **Peak C·Ψ decreases.** Less coherence survives to contribute.
3. **Oscillation count drops.** At γ = 0.01, pair (0,2) oscillates
   through the threshold 13 times (nearly unitary behavior). At
   γ = 0.05, one crossing survives. At γ = 0.20, none.

The critical γ_c for pair (0,2) lies between 0.10 and 0.20. Above
this, the noise wins completely; no pair can build enough coherence
to cross the threshold.

### 9.4 Reversed State |+0+0⟩ (Question 5)

Answered by symmetry: the Heisenberg ring Hamiltonian is invariant
under cyclic permutation. The state |+0+0⟩ = P₁|0+0+⟩, where P₁ is
a one-site cyclic shift. This maps pair (0,2) → (1,3) exactly. The
survivor pair switches from (0,2) to (1,3), with identical crossing
time and max C·Ψ. No separate computation needed.

## 10. Remaining Open Question

1. **Born rule from crossing**: At the crossing point, the diagonal of
   rho_ij gives measurement probabilities. Can the framework predict
   these probabilities independently, without computing the full
   density matrix? This is the central open question for connecting
   crossing to the Born rule.

## 11. Verification

### 11.1 Eigenstate Check

```python
from qutip import basis, tensor, ket2dm, qeye, sigmax, sigmay, sigmaz

N = 4
up, dn = basis(2, 0), basis(2, 1)
plus = (up + dn).unit()
I2 = qeye(2)

def nqubit_op(op, q, n):
    ops = [I2]*n
    ops[q] = op
    return tensor(ops)

# Heisenberg ring
H = sum(J * nqubit_op(op, i, N) * nqubit_op(op, (i+1)%N, N)
        for i in range(N) for op in [sigmax(), sigmay(), sigmaz()])

for name, psi in [("plus4", tensor([plus]*N)),
                   ("0+0+", tensor([dn, plus, dn, plus]))]:
    rho = ket2dm(psi)
    var = (rho * H * H).tr().real - ((rho * H).tr().real)**2
    print(f"{name}: energy variance = {var:.6f}")
# Expected: plus4 → 0.000, 0+0+ → 20.000
```

### 11.2 MCP Verification

```
simulate_subsystem_crossing(state="plus", n_spins=4, gamma=0.05)
# Expected: C = 0.000 for all pairs at all times

simulate_subsystem_crossing(state="bell_pairs", n_spins=4, gamma=0.05)
# Expected: pairs (0,1) and (2,3) cross at t ~ 0.077
```

The MCP tool now supports the alternating state (|0+0+⟩) via Task 010b,
with bidirectional crossing detection via Task 010c:

```
simulate_subsystem_crossing(state="alternating", n_spins=4, gamma=0.05)
# Expected: pair (0,2) crosses up at t ~ 0.275, down at t ~ 0.470
#           pair (1,3) does NOT cross, max CΨ ~ 0.234
#           ring-neighbor pairs borderline (max CΨ ~ 0.251)
```

### 11.3 Key Numbers to Check

1. |+⟩^4 energy variance: exactly 0.000 (eigenstate)
2. |0+0+⟩ energy variance: exactly 20.000
3. |0+0+⟩ unitary, pair (0,1) first crossing: t = 0.073
4. |0+0+⟩ with dephasing, pair (0,2) crossing: t = 0.285
5. |0+0+⟩ with dephasing, pair (0,1) max CΨ: 0.247 (no crossing)
6. N=6 alternating, pair (0,2) max CΨ: 0.131 (no crossing)
7. N=4 Ising alternating, all C_corr: 0.000 (no dynamics)
8. N=4 alternating γ=0.01, pair (0,2) oscillations: 13
9. N=4 alternating γ=0.10, only pair (0,2) crosses
10. N=4 alternating γ=0.20, zero crossings

---

*Previous: [Subsystem Crossing](SUBSYSTEM_CROSSING.md)*
*See also: [N-Scaling Barrier](N_SCALING_BARRIER.md)*
*See also: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
*See also: [WHY_THE_SUM.md](WHY_THE_SUM.md), cross-terms 2*Psi_A*Psi_B survive the blind spot because of the sum structure*
*See also: [Orphaned Results](ORPHANED_RESULTS.md), |+-+-⟩ crosses on ring from zero entanglement via antiferromagnet mechanism*
