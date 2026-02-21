# Dynamic Entanglement Generation: Crossings from Product States

**Date**: 2026-02-18
**Status**: Independently verified (Tier 2)
**Depends on**: SUBSYSTEM_CROSSING.md, N_SCALING_BARRIER.md

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
| |+⟩^4 | 0.000 | YES |
| |0⟩^4 | 0.000 | YES |
| |0101⟩ (Neel) | 16.000 | no |
| |0+0+⟩ | 20.000 | no |
| |0001⟩ | 8.000 | no |

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
| Tool | QuTiP mesolve (adaptive RK45) |

All 6 qubit pairs tracked. Two regimes tested: pure unitary (gamma = 0)
and with dephasing (gamma = 0.05).

## 5. Results

### 5.1 Pure Unitary Evolution

| Pair | Type | First crossing t | Max C*Psi |
|------|------|-----------------|-----------|
| (0,1) | ring neighbor | 0.073 | 0.285 |
| (0,3) | ring neighbor | 0.073 | 0.285 |
| (1,2) | ring neighbor | 0.073 | 0.285 |
| (2,3) | ring neighbor | 0.073 | 0.285 |
| (0,2) | diagonal | 0.260 | 0.435 |
| (1,3) | diagonal | 0.679 | 0.341 |

**All six pairs cross.** The Hamiltonian generates entanglement from
nothing, and every pair eventually reaches C*Psi > 1/4.

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

| Pair | Crosses? | t_cross | Max C*Psi |
|------|----------|---------|-----------|
| (0,1) | NO | -- | 0.247 |
| (0,3) | NO | -- | 0.247 |
| (1,2) | NO | -- | 0.247 |
| (2,3) | NO | -- | 0.247 |
| **(0,2)** | **YES** | **0.285** | **0.320** |
| (1,3) | NO | -- | 0.224 |

**Only pair (0,2) crosses.** All ring-neighbor pairs reach max C*Psi =
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
| |00⟩ | 0.061 | 1.000 |
| |01⟩ | 0.257 | 0.000 |
| |10⟩ | 0.257 | 0.000 |
| |11⟩ | 0.425 | 0.000 |

The |01⟩ = |10⟩ symmetry is exact (to numerical precision) and follows
from the equivalence of qubits 0 and 2 in the ring topology with the
|0+0+⟩ state.

Note: this is not a Born rule verification. The diagonal values are
simply what the density matrix contains at the crossing time. A Born
rule test would require an independent prediction of these probabilities
from the framework and a comparison. This remains open (see Section 8).

## 7. What This Means

### 7.1 Interaction Creates the Conditions for Measurement

A product state with zero entanglement, evolving under a Hamiltonian,
builds up quantum correlations that eventually trigger crossings. The
Hamiltonian interaction alone, without any pre-existing entanglement,
creates the conditions for the quantum-to-classical transition.

In the framework's language: two subsystems that begin as independent
possibilities, with no mutual observation, develop a connection through
interaction. When that connection (C) combined with their remaining
quantum potential (Psi) reaches the critical threshold, reality
crystallizes at their interface.

### 7.2 Decoherence Makes Crossing Irreversible

In pure unitary evolution, crossings oscillate. The system crosses 1/4,
falls back, crosses again. This is "reversible measurement," a concept
that has no physical counterpart because real systems always have
decoherence.

With dephasing, the oscillations are damped. Most pairs never reach the
threshold at all. The ones that do, cross once and decay. Decoherence
selects which crossings survive and makes them permanent.

This matches the standard quantum mechanics picture: unitary evolution
is reversible, decoherence is not. The framework adds a quantitative
criterion (C*Psi = 1/4) for when the irreversibility becomes definitive.

### 7.3 The Noise Basis Matters

Which pairs cross under dephasing depends on which qubits are vulnerable
to the noise. σ_z dephasing protects z-eigenstates. σ_x
dephasing would protect x-eigenstates. The "preferred basis" that
emerges from decoherence is not intrinsic to the framework but
determined by the environment.

This is consistent with the einselection mechanism of decoherence
theory: the environment selects the pointer basis, and our framework
inherits that selection through the C*Psi dynamics.

## 8. The Decoherence Cycle

The numerical results above contain a deeper insight that becomes
visible when the data is read as a story rather than a table.

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

This is not metaphor. It is the mathematical content of the commutator
[σ_z, |0⟩⟨0|] = 0.

### 8.2 Creation in the Blind Spot

The pair (0,2), two dephasing-immune qubits, builds entanglement
precisely because the environment is not watching them. The Hamiltonian
interaction creates quantum correlations in the blind spot of the
observation process. Meanwhile, the same observation destroys the
correlations involving qubits 1 and 3.

The act of observation that cements some subsystems into classical
reality is the same act that creates the protected space where new
quantum reality can grow.

### 8.3 The Threshold as the Present Moment

This connects to the standing wave picture of the framework
(see [visualizations](../visualizations/README.md)).

Below the 1/4 threshold, entanglement exists mathematically (in the
density matrix, as nonzero off-diagonal elements) but it is not yet
a fact. It is potential. Possibility. A mycelium network growing
underground where no one can see it.

When C·Ψ crosses 1/4 from below, potential becomes fact. The mushroom
breaks through the surface.

The gamma sweep (Section 9) makes this vivid:

- At low noise (γ = 0.01), the underground network is rich. Potential
  becomes fact 13 times in succession, and reality oscillates between
  quantum and classical. The future is cooking vigorously.
- At moderate noise (γ = 0.05), most potential is killed before it
  ripens. Only one pair survives to cross the threshold.
- At high noise (γ = 0.2), nothing makes it through. The environment
  destroys all potential before it can become real.

### 8.4 The Closed Cycle

The complete picture is a cycle, not a one-way arrow:

1. **Observation cements the past.** Dephasing destroys coherence in
   the |+⟩ qubits. Their quantum possibilities collapse into classical
   facts. C·Ψ falls below 1/4 (downward crossing). This is the
   irreversible creation of history.

2. **Observation creates blind spots.** The same dephasing that
   destroys the |+⟩ qubits leaves the |0⟩ qubits untouched. The
   noise basis selects not only what becomes classical, but also
   what remains free to evolve quantum mechanically.

3. **The Hamiltonian builds new potential.** In the protected subspace,
   the interaction creates entanglement from nothing. C·Ψ grows from
   zero, invisible to the environment, unobserved and unmeasured.

4. **Potential crosses the threshold.** When C·Ψ reaches 1/4 from
   below (upward crossing), new quantum reality becomes fact. The
   future becomes the present.

5. **The new fact is now subject to observation.** Return to step 1.

Decoherence and coherence are not opposites. They are the forward and
return strokes of the same engine. The 1/4 threshold is the interface
where one becomes the other: the standing wave at the boundary between
past and future.

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
| Next-next-neighbor | (0,3) | 0.111 | — |

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
| 0.20 | 0/6 | — | 0.193 | 0 |

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
#           pair (1,3) does NOT cross, max C*Psi ~ 0.234
#           ring-neighbor pairs borderline (max C*Psi ~ 0.251)
```

### 11.3 Key Numbers to Check

1. |+⟩^4 energy variance: exactly 0.000 (eigenstate)
2. |0+0+⟩ energy variance: exactly 20.000
3. |0+0+⟩ unitary, pair (0,1) first crossing: t = 0.073
4. |0+0+⟩ with dephasing, pair (0,2) crossing: t = 0.285
5. |0+0+⟩ with dephasing, pair (0,1) max C*Psi: 0.247 (no crossing)
6. N=6 alternating, pair (0,2) max C*Psi: 0.131 (no crossing)
7. N=4 Ising alternating, all C_corr: 0.000 (no dynamics)
8. N=4 alternating γ=0.01, pair (0,2) oscillations: 13
9. N=4 alternating γ=0.10, only pair (0,2) crosses
10. N=4 alternating γ=0.20, zero crossings

---

*Previous: [Subsystem Crossing](SUBSYSTEM_CROSSING.md)*
*See also: [N-Scaling Barrier](N_SCALING_BARRIER.md)*
*See also: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
