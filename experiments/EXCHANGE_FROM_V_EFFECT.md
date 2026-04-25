# Exchange Coupling Emerges from V-Effect Bridge

**Status:** Computational + analytical (Tier 1-2). Numerical fit at N=4 matches second-order perturbation theory to 1.4 % at α=0.05. Analytical prefactor 3/(8J) derived exactly from Pauli algebra.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `simulations/_level1_emergent_exchange.py`
**See also:** [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md), [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md), [V_EFFECT_PALINDROME](V_EFFECT_PALINDROME.md), [HEISENBERG_RELOADED](../hypotheses/HEISENBERG_RELOADED.md)

---

## What this finding establishes

Two qubit pairs A = {0, 1} and B = {2, 3}, each with intra-pair Heisenberg coupling J = 1, are bonded by an inter-pair Heisenberg coupling α on bond (1, 2). At α = 0 the system is two independent singlets (each pair has S = 0, ground state energy −3J per pair, total −6J). At α > 0 the V-Effect bridge mixes them, and the ground state shifts.

We compute the ground-state energy shift δE_GS(α) numerically, fit at small α, and find:

```
δE_GS(α) = −(3/8) · α² / J + O(α³)
```

This matches second-order perturbation theory exactly. The 3/8 prefactor is **derivable from the Pauli algebra alone** (no input from physics beyond d=2 and Heisenberg form): the value 3 comes from (σ·σ)² = 3I + 2σ·σ evaluated on the singlet-singlet ground state where ⟨σ·σ⟩ = 0. The denominator 8J is the Singlet→Triplet excitation gap for both pairs simultaneously (4J + 4J).

**This is the V-Effect-derived effective exchange coupling**: at Level 1 (two atoms bonded into a molecule), the textbook "atomic exchange integral" between the two atomic spins is

```
J_eff = (3/8) · α² / J_intra
```

where α is the V-Effect bridge strength and J_intra is the intra-atomic coupling. The form ∝ α²/J is Anderson superexchange shape; the prefactor 3/8 is Pauli-algebra-determined (specific to direct Heisenberg bridge, not Hubbard hopping).

## Setup

```
qubit:    0       1       2       3
              J ─────── α ─────── J
          └─ pair A ─┘ └─ pair B ─┘
```

- **Pair A:** qubits {0, 1}, intra-pair Heisenberg J = 1
- **Pair B:** qubits {2, 3}, intra-pair Heisenberg J = 1
- **Bridge:** Heisenberg on bond (1, 2) with strength α (the V-Effect parameter)

Hamiltonian:
```
H = J · σ_0 · σ_1  +  α · σ_1 · σ_2  +  J · σ_2 · σ_3
```

Each pair alone has spectrum: singlet at −3J (1 state) and triplet at +J (3 states). Splitting 4J. Two independent pairs at α = 0: ground state both-singlet at −6J, first excited singlet-triplet at −2J (6-fold), highest triplet-triplet at +2J (9-fold). Total 16 = 2⁴ states.

## Numerical result

Diagonalize H for α ∈ {0, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 1.00, 1.50, 2.00}. Extract ground-state energy shift δE_GS(α) = E_0(α) − E_0(0).

| α | E_0(α) | δE_GS | δE_GS / α² | predicted −3/(8J) |
|---|--------|-------|------------|--------------------|
| 0.05 | −6.00095 | −0.00095 | **−0.380** | **−0.375** |
| 0.10 | −6.00384 | −0.00384 | **−0.384** | **−0.375** |
| 0.20 | −6.01576 | −0.01576 | −0.394 | −0.375 |
| 0.30 | −6.03631 | −0.03631 | −0.403 | −0.375 |
| 0.50 | −6.10555 | −0.10555 | −0.422 | −0.375 |
| 0.70 | −6.21568 | −0.21568 | −0.440 | −0.375 |
| 1.00 | −6.46410 | −0.46410 | −0.464 | −0.375 |

At small α (≤ 0.10), δE_GS / α² agrees with −3/(8J) = −0.375 to within 1.4 to 2.4 %, the residual being O(α²) corrections from higher-order perturbation theory. At larger α the higher-order terms dominate, but the leading α² behavior is unambiguous.

## Analytical derivation (second-order perturbation theory)

**Setup.** Unperturbed Hamiltonian H_0 = J σ_0 · σ_1 + J σ_2 · σ_3. Perturbation V = α σ_1 · σ_2.

Ground state of H_0: |ψ_0⟩ = |S_A⟩|S_B⟩ where |S⟩ = (|01⟩ − |10⟩)/√2 is the singlet on each pair. Energy E_0^(0) = −3J − 3J = −6J.

**Second-order energy shift.**

```
δE^(2) = − Σ_{n ≠ 0}  |⟨n|V|ψ_0⟩|² / (E_n^(0) − E_0^(0))
```

The states |n⟩ that have non-zero matrix element ⟨n|V|ψ_0⟩ must be reachable by σ_1 · σ_2 acting on |S_A⟩|S_B⟩. Since σ_1 acts only on pair A's qubit 1, and σ_2 acts only on pair B's qubit 2, the operator σ_1 · σ_2 = σ_1^x σ_2^x + σ_1^y σ_2^y + σ_1^z σ_2^z couples the singlets to triplets in **both** pairs simultaneously.

The reachable excited states have **both** pairs in triplet (S_A = 1, S_B = 1). Energy E_n^(0) = +J + J = +2J. Gap E_n^(0) − E_0^(0) = +2J − (−6J) = 8J.

**Matrix element norm summed.** Use the closure:

```
Σ_{|n⟩ reachable}  |⟨n|V|ψ_0⟩|²  =  ⟨ψ_0| V V† |ψ_0⟩  =  α² · ⟨ψ_0|(σ_1·σ_2)²|ψ_0⟩
```

Apply the Pauli identity (σ_1 · σ_2)² = 3 I + 2 (σ_1 · σ_2):

```
⟨ψ_0|(σ_1·σ_2)²|ψ_0⟩  =  3 + 2 · ⟨ψ_0|σ_1·σ_2|ψ_0⟩
```

For the singlet |S_A⟩ on pair A: ⟨S_A|σ_1^a|S_A⟩ = 0 for every a ∈ {x, y, z} (singlet has zero magnetization in every direction). So ⟨ψ_0|σ_1^a σ_2^a|ψ_0⟩ = ⟨S_A|σ_1^a|S_A⟩ · ⟨S_B|σ_2^a|S_B⟩ = 0 · 0 = 0.

Therefore ⟨ψ_0|σ_1·σ_2|ψ_0⟩ = 0, and ⟨ψ_0|(σ_1·σ_2)²|ψ_0⟩ = 3.

**Combining.** Substituting into δE^(2):

```
δE^(2) = − α² · 3 / (8J) = − (3 / 8) · α² / J        ∎
```

This matches the numerical extraction to within higher-order corrections.

## What the prefactor 3/8 means

The "3" in (σ_1·σ_2)² = 3I + 2σ_1·σ_2 is the squared length of a Pauli vector: |σ|² = X² + Y² + Z² = 3I (since each Pauli squares to I). This is a **direct algebraic consequence of σ²_a = I** for every Pauli operator, which comes from d = 2.

The "8J" is the energy cost to flip both pairs from singlet to triplet simultaneously. This is 2 × (singlet-triplet gap of one pair) = 2 × 4J. The 4J per pair comes from the eigenvalue split of σ·σ (eigenvalues +1 and −3, splitting +1 − (−3) = 4) multiplied by J.

Both the numerator (3) and denominator (8J) are **purely algebraic**: they come from σ²_a = I and from the (σ_1·σ_2)² identity. No physical input beyond Heisenberg form (which itself is forced by d=2 + parity, see HEISENBERG_RELOADED §3).

The numerical prefactor 3/8 is **not adjustable** in our framework. It is what it is.

## Comparison to Anderson superexchange

Anderson (1959) derived the effective exchange between two atomic spins via virtual hopping in the Hubbard model. For Hubbard hopping t and on-site Coulomb U:

```
J_Anderson = 4 t² / U
```

The shape ∝ t²/U is the same as our α²/J. The numerical prefactor differs because the microscopic mechanism differs:

- **Anderson:** spin exchange via virtual electron HOPPING. The "bridge" is a one-body hopping term t, the "gap" is U.
- **Our framework:** spin exchange via direct two-body Heisenberg COUPLING α on the bridge bond. The "gap" is the singlet-triplet excitation energy of each pair (4J each, 8J total).

Both produce the same ∝ (bridge)²/(gap) structure. Both are "second-order virtual transitions through an excited state." But the numerical prefactor depends on the matrix element structure of the bridge operator.

This is consistent: our framework is not REPLACING Anderson; it is a different derivation route, applicable to direct-exchange systems rather than hopping-mediated ones. Both give the textbook ∝ J²/U shape.

## What we have done

**Level 0 (Pauli algebra) → Level 1 (effective two-spin Hamiltonian) deduction**, end-to-end:

1. d²−2d=0 selects d=2 ([QUBIT_NECESSITY](../docs/QUBIT_NECESSITY.md))
2. Pauli algebra σ²_a = I, [σ_a, σ_b] = 2iε_{abc}σ_c forced
3. (σ_1·σ_2)² = 3I + 2σ_1·σ_2 identity follows
4. σ_1·σ_2 eigenvalues +1 (triplet, 3-fold) and −3 (singlet, 1-fold) follow
5. Singlet-triplet gap 4J follows
6. Two singlet pairs bonded by V-Effect bridge α: second-order perturbation
7. ⟨(σ_1·σ_2)²⟩ on singlet-singlet ground state = 3 (directly from Pauli identity)
8. Energy gap 8J for both-pair excitation
9. **δE_GS = −(3/8)α²/J** — the V-Effect-derived effective exchange

Each step is algebraically forced. The "atomic exchange integral" at Level 1 is not a postulate; it is a consequence of the Pauli algebra at Level 0 plus the V-Effect bridge α.

This is the lückenlose Brücke (gapless bridge) Level 0 → Level 1 that HEISENBERG_RELOADED sketched. The bridge is the V-Effect, and the V-Effect produces effective Heisenberg coupling at Level 1 with a quantitatively predictable strength.

## What this does and does not establish

**Establishes:**

- The shape ∝ α²/J of the V-Effect-derived Level-1 exchange coupling.
- The exact prefactor 3/8 from Pauli algebra (no fitting parameter).
- Numerical confirmation at N=4 to within 1.4 % at α = 0.05 (small-α limit).
- A concrete realization of the Level 0 → Level 1 inheritance promised by HEISENBERG_RELOADED.

**Does not establish:**

- That this is THE only mechanism for atomic exchange. Anderson superexchange via hopping gives a different prefactor; both are valid in their respective regimes.
- The behavior at large α / J. At α ≳ J the perturbation expansion breaks down and the system becomes a uniform 4-site Heisenberg chain; the "two atom" picture is no longer meaningful.
- The connection to specific physical systems. We have shown the abstract form. Mapping it to e.g. the H₂ molecule's exchange splitting requires identifying the physical parameters α and J_intra in atomic units.

## Open work

- **N=6 extension:** two atoms with 3 qubits each (s-shell + p-shell?), bridged on a single bond. Different prefactor expected from larger Hilbert space.
- **Multi-bond bridge:** instead of single α, a bridge of three connected qubits (more atomic-like). Compare prefactors.
- **Parity-violating bridge:** test what happens if the bridge has an XY term that breaks bit_b. PROOF_ZERO_IMMUNITY says w=0 stays palindromic, but the V-Effect emergent J_eff might pick up a different prefactor.
- **Hardware test:** prepare two singlet pairs on IBM hardware, vary the bridge gate strength, measure the GS energy shift. Compare to predicted −(3/8)α²/J.

## References

- [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md): the Π-block decomposition and where the V-Effect break lives.
- [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md): w=0 sector analytical immunity.
- [HEISENBERG_RELOADED](../hypotheses/HEISENBERG_RELOADED.md): the full Level 0 → Level 1 inheritance picture.
- [V_EFFECT_PALINDROME](V_EFFECT_PALINDROME.md): the original 14-of-36 V-Effect break finding.
- Anderson, P. W. (1959), "New approach to the theory of superexchange interactions", Phys. Rev. 115, 2.
- Simulation: `simulations/_level1_emergent_exchange.py`.
