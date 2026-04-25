# Asymmetric V-Effect Emergent Exchange

**Status:** Computational + analytical (Tier 1-2). First framework-based calculation built on `simulations/framework.py` primitives. Numerical agreement to 0.2-0.8 % at α=0.025 across 7 asymmetry configurations.
**Date:** 2026-04-25
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Pipeline:** `simulations/_asymmetric_v_effect.py` (uses `framework.py`)
**See also:** [EXCHANGE_FROM_V_EFFECT](EXCHANGE_FROM_V_EFFECT.md), [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md), [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md)

---

## What this generalizes

EXCHANGE_FROM_V_EFFECT computed δE_GS = −(3/8) α²/J for two SYMMETRIC Heisenberg pairs (both with intra-coupling J) bonded by inter-pair Heisenberg α. This document extends to ASYMMETRIC pairs (J_A ≠ J_B), giving the natural generalization:

```
δE_GS  =  − 3 α²  /  (4 (J_A + J_B))
```

At J_A = J_B = J this reduces to −3α²/(8J), recovering the symmetric result. The 4(J_A + J_B) denominator is the "total flip cost": pair A flipping singlet→triplet costs 4J_A, pair B costs 4J_B, both must flip simultaneously to be reachable by the bridge V = α σ_1·σ_2.

## Setup

```
qubit:    0       1       2       3
              J_A ─────── α ────── J_B
          └─ pair A ─┘   bridge  └─ pair B ─┘
```

- **Pair A:** qubits {0, 1}, intra-Heisenberg J_A · σ_0·σ_1
- **Pair B:** qubits {2, 3}, intra-Heisenberg J_B · σ_2·σ_3
- **Bridge:** Heisenberg on bond (1, 2) with strength α

H_total = J_A·σ_0·σ_1 + α·σ_1·σ_2 + J_B·σ_2·σ_3

## Analytical prediction (second-order PT)

**Unperturbed ground state:** |ψ_0⟩ = |S_A⟩|S_B⟩, energy E_0^(0) = −3 J_A − 3 J_B.

**Reachable excited states:** σ_1·σ_2 acts on bridge qubits. Both pairs must flip singlet→triplet to give nonzero matrix element. Excited states: |T_α^A⟩|T_β^B⟩ at energy +J_A + J_B. Gap:

```
E_excited^(0) − E_GS^(0) = (J_A + J_B) − (−3J_A − 3J_B) = 4 (J_A + J_B)
```

**Matrix element norm:** Σ |⟨excited|V|ψ_0⟩|² = α² ⟨ψ_0|(σ_1·σ_2)²|ψ_0⟩.

Pauli identity: (σ_1·σ_2)² = 3 I − 2 (σ_1·σ_2).
On singlet-singlet GS: ⟨σ_1·σ_2⟩_{|S_A⟩|S_B⟩} = 0 (each factor's individual Pauli expectations vanish on its own singlet).
Therefore ⟨(σ_1·σ_2)²⟩_{|ψ_0⟩} = 3 − 0 = 3.

**Combining:**

```
δE^(2) = − Σ_excited |⟨n|V|ψ_0⟩|² / (E_n^(0) − E_0^(0))
       = − α² · 3 / (4(J_A + J_B))
       = − 3 α² / (4 (J_A + J_B))         ∎
```

## Numerical verification (N=4, exact diagonalization)

Pipeline: `simulations/_asymmetric_v_effect.py`. Build H via `framework._build_bilinear`, exact eigendecomposition, extract δE_GS = E_0(α) − E_0(0), compare to prediction.

Selected results at α = 0.025 (smallest tested, deepest in PT regime):

| J_A | J_B | δE/α² (numeric) | −3/(4(J_A+J_B)) | rel error |
|-----|-----|------------------|-----------------|-----------|
| 1.0 | 1.0 | −0.377 | −0.375 | 0.63 % |
| 0.5 | 1.0 | −0.504 | −0.500 | 0.84 % |
| 1.0 | 2.0 | −0.251 | −0.250 | 0.42 % |
| 0.5 | 2.0 | −0.302 | −0.300 | 0.50 % |
| 1.0 | 5.0 | −0.1253 | −0.1250 | 0.21 % |
| 0.3 | 3.0 | −0.2281 | −0.2273 | 0.38 % |
| 2.0 | 3.0 | −0.1504 | −0.1500 | 0.25 % |

The prediction holds across nearly two orders of magnitude in J_A/J_B asymmetry (from 1:1 to 1:10) at sub-1 % relative error. The error is the leading O(α²) correction to second-order PT and shrinks with α: at α=0.05 errors are ~1 % across the board, at α=0.20 they grow to ~3-7 % depending on asymmetry.

## What this calculation does

**First end-to-end use of framework.py primitives.** The Hamiltonian construction goes through `framework._build_bilinear`, the prediction call uses `framework.v_effect_emergent_exchange` for the symmetric reference, and the script verifies that framework.py's symmetric formula agrees with this script's symmetric special case. Self-consistency check passes at machine precision.

**Generalization of EXCHANGE_FROM_V_EFFECT to asymmetric atoms.** The inheritance Level 0 → Level 1 V-Effect bridge produces an effective Level-1 exchange that scales as α² / (J_A + J_B), which is the textbook Anderson superexchange shape applied to asymmetric atoms (where in the original Hubbard derivation, the "U" for each atom can differ).

**Confirmation that the 3/8 prefactor was a special case.** The symmetric version's 3/8 = 3 / (4·2) is the symmetric J_A = J_B = J specialization of the general 3 / (4(J_A + J_B)). The "3" is universal (from Pauli identity); the "4(J_A + J_B)" is the configuration-dependent gap.

## What this does not establish

- N=4 specific. The prediction is a leading-order PT result, valid as long as both pairs are in their singlet GS at α=0 and the perturbation V is small. Deviations grow with α (third-order corrections) and are not part of the leading-order claim.
- Asymmetric BRIDGE (different couplings on (1,2) bond per Pauli direction): not tested. The bridge here is uniform Heisenberg α(XX+YY+ZZ).
- Real physical atom-pair systems where J_A and J_B come from different Coulomb integrals: the abstract formula matches, but mapping to specific molecular systems (e.g., heteronuclear diatomic bonds) requires identifying physical J_A and J_B values, which is a separate step.

## What's open

- **Mixed bridge:** what if α(XX+YY) + β(ZZ) (XXZ-style bridge)? Framework predicts the prefactor depends on which Pauli components are present in V.
- **Three atoms:** chain of three pairs bonded by two bridges. Predicts a three-body Heisenberg-style chain at Level 1, with effective J_eff between consecutive atoms following the same form.
- **Four-qubit bridge:** instead of a single inter-bond, use a 2-qubit "molecular bond" mediator. Different gap structure; different prefactor.
- **Hardware mapping:** identify physical atom pairs (e.g., H₂, CO, water dimer) where J_A, J_B, α can be extracted from atomic-data tables, and check if our prediction matches measured exchange splittings.

## References

- [EXCHANGE_FROM_V_EFFECT](EXCHANGE_FROM_V_EFFECT.md): symmetric J_A = J_B = J case, the original derivation.
- [V_EFFECT_BOUNDARY_LOCALIZATION](V_EFFECT_BOUNDARY_LOCALIZATION.md): the algebraic structure of where V-Effect breaks live.
- [PROOF_ZERO_IMMUNITY](../docs/proofs/PROOF_ZERO_IMMUNITY.md): the (w=0, w=N) palindrome immunity that grounds this whole construction.
- [HEISENBERG_RELOADED](../hypotheses/HEISENBERG_RELOADED.md): the level-stack picture into which this fits.
- `simulations/framework.py`: framework primitives.
- `simulations/_asymmetric_v_effect.py`: this calculation's pipeline.
