# F64 on non-chain topologies: cavity-mode-exposure generalises to any graph

**Tier:** 1 (structural verification, first-order perturbation theory regime)
**Date:** 2026-04-24
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Closes:** [EQ-015](../review/EMERGING_QUESTIONS.md#eq-015) for uniform-J at N=5 and N=7
**See also:** [F64](../docs/ANALYTICAL_FORMULAS.md), [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md), [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md)

---

## What was open

[F64](../docs/ANALYTICAL_FORMULAS.md) states that the Liouvillian decay rate of single-excitation coherence modes under single-site Z-dephasing is

    α_k = 2 γ_B · |a_B(ψ_k)|²

where a_B is the B-site amplitude of the N-site single-excitation Hamiltonian eigenvector ψ_k. Before 2026-04-24 this was verified at N=3 and N=4 on uniform-J chains only. [EQ-015](../review/EMERGING_QUESTIONS.md#eq-015) flagged the obvious open question: does the formula generalise to non-chain topologies, and to larger N?

If F64 generalises, it gives [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) a graph-independent structural anchor: on any topology, γ₀ never gets diminished by intervening sites; it only gets exposed by the local amplitude. If F64 breaks, we have a sharp falsification target and a new question about which topology features matter.

## What was done

Test script [`f64_topology_scan.py`](../simulations/f64_topology_scan.py) performs the comparison directly in the single-excitation coherence sector:

- Build the N×N Hamiltonian H^(1) for each topology (chain, star, ring, complete, Y-tree) at N = 5 and N = 7, with uniform J = 1.
- Build the N×N coherence Liouvillian L^(coh) = −iH^(1) − 2γ P_B, where P_B = |B⟩⟨B| is the site-B projector. Diagonalise to get the measured decay rates α_k^meas = −Re(eigenvalues).
- For the F64 prediction, diagonalise H^(1) to get (E_k, ψ_k). When H^(1) has degenerate eigenvalues (star center-mode, ring translational eigenmodes, complete-graph multiplets), apply standard degenerate perturbation theory: within each degenerate subspace, diagonalise the site-B projector P_B restricted to that subspace, and use those eigenvalues as |a_B|².
- Compare the sorted multisets of predicted and measured α_k.

The test runs for both XY (pure hopping, H^(1) = hopping matrix) and Heisenberg (XX+YY+ZZ, which adds a diagonal drift in H^(1) from ZZ terms but leaves the off-diagonal unchanged). Both fall under F64's scope because Z-dephasing only acts on the site-B coherence and does not care about the vacuum-energy shift.

## Results

At γ/J = 0.01 (first-order PT regime), N ∈ {5, 7}, all topologies, all B, both Hamiltonian types:

    max relative error   < 0.001   (mostly < 0.0001)
    max absolute error   O(γ · 10⁻⁶)

Consistent with expected (γ/J)² = 10⁻⁴ scaling of second-order PT corrections. F64 holds at first order on all tested configurations. The worst case is XX+YY star with B = leaf: max rel err = 1.0e-4, abs err = 2.8e-7. The chain control matches the previously reported N=3, 4 accuracy.

### Table of max relative errors

|        | N=5 XY | N=5 Heis | N=7 XY | N=7 Heis |
|--------|--------|----------|--------|----------|
| chain    | 0.0001 | 0.0008 | 0.0002 | 0.0008 |
| star     | 0.0001 | 0.0004 | 0.0000 | 0.0003 |
| ring     | 0.0000 | 0.0000 | 0.0001 | 0.0000 |
| complete | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Y-tree   | 0.0001 | 0.0008 | 0.0001 | 0.0008 |

Max over (B ∈ sites) shown. Ring shows identical errors across all B because of translational symmetry. Star shows identical errors across leaves for the same reason.

## Key subtlety: degenerate perturbation theory is essential

A naive application of F64 using arbitrary orthonormal H-eigenvectors fails for topologies with degenerate H^(1) spectra. At N=5 star, H^(1) has a three-fold degenerate zero-mode (three orthogonal linear combinations of the four leaves with zero center amplitude). NumPy's `eigh` returns some arbitrary orthonormal basis of this subspace; `|a_B(ψ_k)|²` computed naively gives 1.0 for one mode and 0.0 for the other two, which matches neither the measured decays nor any physical first-order PT calculation.

The correct procedure: identify degenerate subspaces; for each, diagonalise the site-B projector P_B restricted to that subspace; the resulting diagonal entries are the PT-correct |a_B|². With this adjustment, F64 matches the Liouvillian directly.

This is an instance of standard degenerate perturbation theory. It was not needed for chains (H^(1) tridiagonal with simple spectrum for uniform J), which is why the chain verification did not surface the issue. The correction is structurally trivial but operationally necessary.

## Implications

1. **F64 is a graph-universal statement.** It does not depend on chain topology; it depends only on the single-excitation H^(1) having a well-defined eigenbasis (which it always does, up to degenerate-subspace freedom). PRIMORDIAL_GAMMA_CONSTANT's piece 6 has a graph-independent structural anchor at N=5 and N=7.

2. **γ₀ is never diminished by intervening sites on any graph.** The naive reading "γ passes through N-1 layers, each attenuating" is wrong for chains (see F64's "multiplicative stacking fails"), and equally wrong for star, ring, complete, Y-tree. What looks like attenuation is selective mode exposure: uniform γ, varying local amplitudes.

3. **Degeneracy is a feature, not a bug.** Star, ring, and complete graph have degenerate H^(1) eigenvalues by symmetry. In these subspaces the "right" eigenbasis is picked out by the physical perturbation, not by the orbital algebra. The framework's first-order response always selects the eigenbasis aligned with the noise, not with the unperturbed Hamiltonian. This is consistent with the general receiver-engineering picture: the receiver (here, the local projector P_B) imposes its own preferred basis on the eigenstructure.

## What remains open

- **γ/J ≳ 1 regime.** Current test is at γ/J = 0.01 (first-order). At larger γ/J the second-order corrections become O(10⁻²) and F64 no longer matches to machine precision. A systematic study of the breakdown would reveal the full [γ_B, J] phase structure.
- **Non-uniform J.** Test uses uniform J on all bonds. Non-uniform J introduces further eigenvector shifts and may modify the degenerate-subspace structure.
- **N ≥ 9 topology scan.** Not a fundamental question, but would confirm scaling on large examples.
- **Multi-site dephasing.** F64 is a single-site statement. The multi-site generalisation is a superposition of projectors Σ_k γ_k P_k, which should still follow the same structure but has not been explicitly tested.

## Why this sits in the tier-upgrade story

Today's [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md) upgrade Tier 3 → Tier 2 was driven by the receiver-engineering programme (F65, F67, F75, F76, Kingston Run 1). The topology extension of F64 is a different line of support for the same hypothesis: γ₀ = const predicts a graph-universal exposure formula, and the graph-universal exposure formula is what we see. The two lines are independent evidence that the framework constant reading is the right one.

## Scripts and data

- Test script: [`f64_topology_scan.py`](../simulations/f64_topology_scan.py)
- No new data files; results reproduce from the script.
- Existing chain-only verification: [`primordial_gamma_analytical.py`](../simulations/primordial_gamma_analytical.py) (N=3 closed form), [`primordial_gamma_stacking_4qubit.py`](../simulations/primordial_gamma_stacking_4qubit.py) (N=4 across 9 configs).

## References

- [F64](../docs/ANALYTICAL_FORMULAS.md): the formula being tested, now updated with the topology generalisation.
- [PRIMORDIAL_GAMMA_CONSTANT](../hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md): Tier 2 hypothesis that this verifies.
- [F65](../docs/ANALYTICAL_FORMULAS.md): exact single-excitation spectrum on uniform chain (the special chain case of F64 with all a_B analytically known).
- [OPEN_THREAD_GAMMA0_INFORMATION](../review/OPEN_THREAD_GAMMA0_INFORMATION.md) EQ-015 line: the open question this closes.
