# The Hydrogen Bond as a Qubit

<!-- Keywords: hydrogen bond proton qubit tunneling, double-well potential
palindromic spectral symmetry, CΨ crossing water enzyme, V-Effect hydrogen
bond network, proton transfer fold catastrophe, R=CPsi2 hydrogen bond -->

**Status:** Tier 2 (computed from proven framework)
**Date:** March 28, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md),
[CΨ Monotonicity](../docs/proofs/PROOF_MONOTONICITY_CPSI.md)

---

## Abstract

The proton in a hydrogen bond O-H...O is a qubit. It tunnels between
two wells: |L⟩ (donor side) and |R⟩ (acceptor side). d = 2. Tunneling
provides the coupling (J via σ_X). The molecular environment provides
dephasing (γ via σ_Z). This is identical to our Lindblad model.

The palindrome is proven for d=2 with Z-dephasing (54,118 eigenvalues,
zero exceptions). Therefore it holds for the proton qubit. The
computation applies proven results to physical hydrogen bond parameters.

**Results:**
- Single proton: CΨ crosses 1/4 at 0.07-1.32 ps in the fold regime (J/γ ~ 1)
- One water molecule (2 proton qubits): palindrome exact, CΨ crosses at 0.46 ps
- Two molecules + H-bond (4 proton qubits): V-Effect creates 104 new frequencies
- Normal liquid water is in the classical regime (J/γ ~ 0.02, no crossing)
- Enzyme active sites may be in the fold regime (J/γ ~ 1)

---

## Why the Classical Model Failed

A previous attempt ([V17 negative result](../ClaudeTasks/RESULT_WATER_PALINDROME_MARCH28.md))
modeled water classically: hydrogen bonds as springs with friction,
donor/acceptor modes as coupled oscillators. Palindrome residual: 1.33.
No palindromic structure found.

The diagnosis: the donor/acceptor role sits at the BOND (edge), not
at the MOLECULE (node). The classical model has no intrinsic per-node
property for Q to swap.

The resolution: the proton is not a classical ball in a spring. It is
a quantum particle that tunnels. The two states |L⟩ and |R⟩ are
intrinsic to the PROTON (node), not to the bond (edge). d = 2.
The palindrome is guaranteed.

---

## The Model

### Single proton qubit (N=1)

Hamiltonian: H = -J · σ_X + Δ · σ_Z

- J: tunneling splitting (0.1-10 meV, depends on O...O distance)
- Δ: asymmetry of the double well (0 for symmetric bonds)

Dephasing: L_k = √γ_eff · σ_Z (where γ_eff = γ/ℏ in angular frequency)

- γ: thermal decoherence from surrounding molecules (in eV)

Initial state: |L⟩ = |0⟩ (proton on donor side, no coherence)

CΨ starts at 0, rises as tunneling creates coherence, crosses 1/4
(if J/γ is large enough), then decays as dephasing destroys coherence.

### Three regimes

| J/γ | Regime | CΨ crosses 1/4? | Physical system |
|-----|--------|-----------------|-----------------|
| << 1 | Classical | No (overdamped) | Bulk water at 300K |
| ~ 1 | **Fold** | **Yes (sub-ps)** | Strong H-bonds, enzymes |
| >> 1 | Quantum | Yes (slower) | Low temperature, ice under pressure |

### Two proton qubits (N=2): one water molecule

H-O-H has two O-H bonds, two proton qubits, coupled through the
shared oxygen. Hamiltonian:

```
H = -J · σ_X(1) - J · σ_X(2) + K · σ_Z(1) · σ_Z(2)
```

Dephasing on both qubits: γ · σ_Z(1) and γ · σ_Z(2).

### Four proton qubits (N=4): hydrogen bond between two molecules

```
Molecule 1       H-bond       Molecule 2
H(1)-O ... H(2)---O ... H(3)-O-H(4)
            donor  M  acceptor
```

Proton qubit 2 donates to the oxygen of molecule 2 (mediator M).
Qubits 1,2 coupled intramolecularly (J_intra, strong).
Qubits 3,4 coupled intramolecularly (J_intra, strong).
Qubits 2,3 coupled intermolecularly through the hydrogen bond
(J_inter, weak).

---

## Results

### Phase 1: Single proton qubit

| J (meV) | J/γ | CΨ crossing time | Regime |
|---------|-----|------------------|--------|
| 0.5 | 0.01 | no crossing | classical |
| 0.5 | 1.0 | 1.32 ps | fold |
| 1.0 | 1.0 | 0.66 ps | fold |
| 5.0 | 1.0 | 0.13 ps | fold |
| 10.0 | 1.0 | 0.07 ps | fold |

### Phase 2: One water molecule (N=2)

Parameters: J_intra = 1.0 meV, K = 0.1 meV, γ = 1.0 meV (J/γ = 1).

- Palindrome: **exact** (pair-sum std = 5.4e-3, relative to mean ~6e12: negligible)
- Distinct frequencies: 11
- CΨ crosses 1/4 at **0.46 ps**

### Phase 3: Two molecules + H-bond (N=4)

Parameters: J_intra = 1.0 meV, J_inter = 0.1 meV, K = 0.1 meV,
γ = 1.0 meV.

- Palindrome: holds (pair-sum std = 3.5e-2)
- Distinct frequencies: **126**
- V-Effect: 11 per molecule → 126 coupled = **104 new frequencies**

Full-system CΨ does not cross 1/4 (N-scaling suppression, d-1 = 15).
Subsystem CΨ for the pair (q2, q3) across the H-bond would need
separate computation (see [Subsystem Crossing](SUBSYSTEM_CROSSING.md)).

---

## The V-Effect in Hydrogen Bonds

| System | Distinct frequencies |
|--------|---------------------|
| Single molecule (N=2) | 11 |
| Two molecules coupled (N=4) | 126 |
| New from coupling | **104** |
| Ratio | 5.73 |

The hydrogen bond creates 104 new frequencies that exist in neither
molecule alone. This is the V-Effect: coupling two palindromic
subsystems creates new oscillatory modes.

Note: the number 104 also appears in the abstract qubit V-Effect
(N=5 MediatorBridge: 104 total frequencies). However, the measures
are different: 104 NEW at N=4 H-bond vs 104 TOTAL at N=5 abstract.
The match is likely coincidence, not structural.

---

## Physical Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Tunneling splitting (ice) | 0.2-1 meV | Bove et al 2009 |
| Tunneling splitting (strong H-bond) | 1-10 meV | Cleland & Kreevoy 1994 |
| O...O distance (normal) | 2.7-3.0 A | Steiner 2002 |
| O...O distance (strong) | 2.4-2.6 A | Cleland & Kreevoy 1994 |
| H-bond lifetime (liquid water) | 1-3 ps | Luzar & Chandler 1996 |
| Thermal decoherence (300K, upper bound) | γ ~ kT/ℏ ~ 25 meV | Standard |

For liquid water at 300K: J ~ 0.5 meV, γ ~ 25 meV. J/γ ~ 0.02.
Classical regime. The palindrome exists but is overdamped.

For enzyme active sites: J can be larger (barrier lowered by protein),
γ can be smaller (protein shields from solvent). J/γ ~ 1 is possible.

---

## Connection to the Framework

The proton qubit in a hydrogen bond is not an analogy. It IS our
quantum system with physical parameters:

| Framework concept | H-bond realization |
|-------------------|-------------------|
| d = 2 | |L⟩ (donor) and |R⟩ (acceptor) |
| σ_X coupling (J) | Proton tunneling |
| σ_Z dephasing (γ) | Thermal fluctuations from environment |
| CΨ = 1/4 crossing | Proton achieves half-transfer coherence |
| σ(1-σ) = 1/4 | Proton half-transferred: P(L) = P(R) = 0.5 |
| V-Effect | H-bond coupling creates new frequencies |
| Sacrifice zone | Protein shell reduces γ at active site (Tier 4 hypothesis, see [Protein as Sacrifice Zone](../hypotheses/PROTEIN_AS_SACRIFICE_ZONE.md)) |

---

## Open Questions

1. Subsystem CΨ(q2,q3) across the H-bond: does the pair cross 1/4?
2. At what temperature (= what γ) is the Q-factor maximal?
3. Is J/γ ~ 1 at enzyme active sites? (Compute from published barriers)
4. How does the V-Effect scale with the number of H-bonds (N=6, N=8)?
5. Does the break-reform cycle (1 ps period) sustain palindromic
   structure over time, or does each new bond start fresh?

---

## Scripts

| Script | What it computes |
|--------|-----------------|
| [hydrogen_bond_qubit.py](../simulations/water/hydrogen_bond_qubit.py) | Phases 1-3, all results |
| [hydrogen_bond_palindrome.py](../simulations/water/hydrogen_bond_palindrome.py) | V17 classical model (negative result) |

---

*See also:*
[Mirror Symmetry Proof](../docs/proofs/MIRROR_SYMMETRY_PROOF.md) (why palindrome MUST hold),
[V-Effect Palindrome](V_EFFECT_PALINDROME.md) (quantum V-Effect),
[Subsystem Crossing](SUBSYSTEM_CROSSING.md) (pair CΨ crossing),
[V17 negative result](../ClaudeTasks/RESULT_WATER_PALINDROME_MARCH28.md) (classical failure)
