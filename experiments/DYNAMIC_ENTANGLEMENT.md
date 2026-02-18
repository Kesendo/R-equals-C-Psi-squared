# Dynamic Entanglement Generation: Crossings from Product States

**Date**: 2026-02-18
**Status**: Independently verified (Tier 2)
**Depends on**: SUBSYSTEM_CROSSING.md, N_SCALING_BARRIER.md

---

## 1. The Question

SUBSYSTEM_CROSSING.md (Experiment 10) showed that crossing is local and
operates at the level of entangled qubit pairs. Open Question 3 asked:

> Can dynamically generated entanglement from a product state create
> crossings? (|+>^4 under Heisenberg should build entanglement.)

This experiment answers that question. The answer involves two surprises.

## 2. Surprise 1: |+>^N Cannot Generate Entanglement

The initial hypothesis was that |+>^N under a Heisenberg Hamiltonian
would build entanglement over time, eventually producing crossings. This
is wrong.

|+>^N = |+>|+>...|+> is the maximal-spin eigenstate |S=N/2, m_x=N/2>
of the isotropic Heisenberg Hamiltonian. The time evolution gives only a
global phase: e^{-iHt}|+>^N = e^{-iEt}|+>^N. Nothing happens. No
entanglement is created, no correlations develop, C = 0 for all pairs at
all times.

**Energy variance test (QuTiP verification):**

| State | Var(H) = <H^2> - <H>^2 | Eigenstate? |
|-------|------------------------|-------------|
| |+>^4 | 0.000 | YES |
| |0>^4 | 0.000 | YES |
| |0101> (Neel) | 16.000 | no |
| |0+0+> | 20.000 | no |
| |0001> | 8.000 | no |

Both |+>^N and |0>^N are eigenstates of the isotropic Heisenberg
Hamiltonian (zero energy variance). This is a consequence of SU(2)
symmetry: fully polarized states in any direction are highest-weight
vectors and cannot evolve nontrivially.

Experiment 10 showed that |+>^4 has C = 0 for all pairs at all times.
This is now understood: the Hamiltonian cannot change what is already
an eigenstate.

## 3. The Correct Starting State: |0+0+>

To generate entanglement dynamically, we need a product state that is NOT
an eigenstate. The state |0+0+> = |0>|+>|0>|+> has energy variance 20
and breaks the SU(2) symmetry of the Hamiltonian, enabling nontrivial
dynamics.

This state has:
- Zero entanglement at t = 0 (product state)
- Qubits 0, 2 in |0> (z-basis eigenstate)
- Qubits 1, 3 in |+> (x-basis eigenstate)
- No pair has any quantum correlation initially

## 4. Setup

| Parameter | Value |
|-----------|-------|
| N | 4 qubits |
| Hamiltonian | Heisenberg ring (J = 1.0) |
| Noise | local dephasing (sigma_z, gamma = 0.05) |
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
is reversible -- the system oscillates between "measured" and
"unmeasured" states.

Ring-neighbor symmetry: pairs (0,1), (0,3), (1,2), (2,3) are exactly
degenerate due to the ring topology combined with the |0+0+> mirror
symmetry (qubit 0 maps to 2, qubit 1 maps to 3 under the ring
symmetry that exchanges the two |0> qubits and the two |+> qubits).

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

The asymmetry has a simple explanation. Sigma_z dephasing acts as a
projection toward the z-basis. Its effect depends on the initial state
of each qubit:

- Qubits 0 and 2 start in |0> -- an eigenstate of sigma_z. Dephasing
  does not affect them directly.
- Qubits 1 and 3 start in |+> -- a superposition in the sigma_z basis.
  Dephasing destroys their coherence exponentially.

Pair (0,2) consists of the two dephasing-immune qubits. Their mutual
entanglement, built up by the Hamiltonian through intermediate
interactions with qubits 1 and 3, is partially shielded from the
noise. This is not "selection" by the noise but rather differential
vulnerability: pair (0,2) loses coherence more slowly than pairs
containing |+>-qubits.

This is basis-dependent. With sigma_x dephasing instead of sigma_z,
the roles would reverse: the |+> qubits would be immune and pairs
involving qubits 1 and 3 would survive.

## 6. Density Matrix at the Crossing Point

Pair (0,2) at t = 0.285, reduced density matrix diagonal:

| Basis state | Probability | Initial |
|-------------|-------------|---------|
| |00> | 0.061 | 1.000 |
| |01> | 0.257 | 0.000 |
| |10> | 0.257 | 0.000 |
| |11> | 0.425 | 0.000 |

The |01> = |10> symmetry is exact (to numerical precision) and follows
from the equivalence of qubits 0 and 2 in the ring topology with the
|0+0+> state.

Note: this is not a Born rule verification. The diagonal values are
simply what the density matrix contains at the crossing time. A Born
rule test would require an independent prediction of these probabilities
from the framework and a comparison. This remains open (see Section 8).

## 7. What This Means

### 7.1 Interaction Creates the Conditions for Measurement

A product state with zero entanglement, evolving under a Hamiltonian,
builds up quantum correlations that eventually trigger crossings. The
Hamiltonian interaction alone -- without any pre-existing entanglement
-- creates the conditions for the quantum-to-classical transition.

In the framework's language: two subsystems that begin as independent
possibilities, with no mutual observation, develop a connection through
interaction. When that connection (C) combined with their remaining
quantum potential (Psi) reaches the critical threshold, reality
crystallizes at their interface.

### 7.2 Decoherence Makes Crossing Irreversible

In pure unitary evolution, crossings oscillate. The system crosses 1/4,
falls back, crosses again. This is "reversible measurement" -- a concept
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
to the noise. Sigma_z dephasing protects z-eigenstates. Sigma_x
dephasing would protect x-eigenstates. The "preferred basis" that
emerges from decoherence is not intrinsic to the framework but
determined by the environment.

This is consistent with the einselection mechanism of decoherence
theory: the environment selects the pointer basis, and our framework
inherits that selection through the C*Psi dynamics.

## 8. Open Questions

1. **Born rule from crossing**: At the crossing point, the diagonal of
   rho_ij gives measurement probabilities. Can the framework predict
   these probabilities independently, without computing the full
   density matrix? This is the central open question for connecting
   crossing to the Born rule.

2. **Longer chains**: For N=6, N=8, does the Hamiltonian build enough
   entanglement in distant pairs to produce crossings? Or do only
   nearest-neighbor pairs (and perhaps next-nearest) ever cross?

3. **Different Hamiltonians**: The Ising model H = J*sigma_z*sigma_z has
   a different symmetry structure. Which product states generate
   crossings under Ising evolution?

4. **Crossing timing vs. entanglement growth**: Is there a quantitative
   relationship between the rate of entanglement generation (measured by
   concurrence growth) and the time to first crossing?

5. **Reversed initial state**: Does |+0+0> (swapping which qubits are
   in |0> vs |+>) produce the same physics with pairs relabeled?

## 9. Verification

### 9.1 Eigenstate Check

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

### 9.2 MCP Verification

```
simulate_subsystem_crossing(state="plus", n_spins=4, gamma=0.05)
# Expected: C = 0.000 for all pairs at all times

simulate_subsystem_crossing(state="bell_pairs", n_spins=4, gamma=0.05)
# Expected: pairs (0,1) and (2,3) cross at t ~ 0.077
```

Note: the MCP tool does not yet support the |0+0+> state. Verification
of that state requires QuTiP directly or a future MCP update.

### 9.3 Key Numbers to Check

1. |+>^4 energy variance: exactly 0.000 (eigenstate)
2. |0+0+> energy variance: exactly 20.000
3. |0+0+> unitary, pair (0,1) first crossing: t = 0.073
4. |0+0+> with dephasing, pair (0,2) crossing: t = 0.285
5. |0+0+> with dephasing, pair (0,1) max C*Psi: 0.247 (no crossing)

---

*Previous: [Subsystem Crossing](SUBSYSTEM_CROSSING.md)*
*See also: [N-Scaling Barrier](N_SCALING_BARRIER.md)*
*See also: [Crossing Taxonomy](CROSSING_TAXONOMY.md)*
