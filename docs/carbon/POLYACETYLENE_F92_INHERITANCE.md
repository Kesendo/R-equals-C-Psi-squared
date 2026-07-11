# Polyacetylene as F92 Inheritance: The SSH Chain Reads as Anti-Palindromic-J

**Date:** 2026-05-27
**Authors:** Tom + Claude
**Status:** Tier 2 candidate. Framework predictions are tight and follow directly from F92 / F100 / F101 ([the F92 bond anti-palindromic-J proof](../proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md), [the F100 proof](../proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md), [the F101 proof](../proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md), all Tier 1 derived). Numerical verification on the open-system SSH chain is the deferred next pass.
**Continues:** [the benzene Liouvillian palindrome](BENZENE_LIOUVILLIAN_PALINDROME.md) (Holstein/Peierls fork), [the three benzene dephase letters](BENZENE_THREE_DEPHASE_LETTERS.md) (the Klein-V₄ return visit).
**Reads alongside:** [README](README.md) Open Question 4 (graphene K-point) and the candidate doc list.

---

## Abstract

The Su-Schrieffer-Heeger (SSH) model is the canonical one-dimensional topological insulator: a tight-binding chain with alternating single and double bonds, the simplest physical system that carries a topological invariant. Polyacetylene is its real-world realization, the chain `(CH=CH)ₙ` with the same alternating bond pattern. SSH has been studied since 1979, and its topology, soliton modes, and edge states are textbook material.

Here is what the framework adds, and where it lands. The alternating bond couplings of SSH are literally `J_b = J · (1 + δ · (−1)^b)`: an average coupling `J` plus an anti-palindromic deviation `δ`. The framework's F92 theorem (proved for arbitrary anti-palindromic J distributions earlier this year) says that the diagonal-block eigenvalues of the open-system Liouvillian under Z-dephasing depend only on the palindromic part of J, not on the anti-palindromic part. F100 sharpens this on the observable side: the bond-mirror deviations of the closure-breaking coefficient and the per-bond Q_peak observable are exactly odd functions of the anti-palindromic part. F92 and F100 together split the J-axis into "invisible to spectrum" and "exactly-odd in observables".

The structural reading is striking. Polyacetylene's whole topological story sits at the δ = 0 transition point, where the topological invariant flips and the soliton modes appear. In framework language, δ = 0 is exactly the point where J is palindromic and nothing is "invisible". Moving δ away from zero turns on the anti-palindromic content, which is exactly what the topology cares about and exactly what F92 declares spectrum-invisible. Topological matter and the framework's invisible-direction content are the same thing, read out of different vocabularies.

This document writes the framework predictions down explicitly, names what numerical verification would look like, and flags the open structural question: whether SSH's edge soliton modes can be identified with the F92-invisible content of a specific Pauli-basis sector.

---

## The SSH model in framework language

The standard SSH Hamiltonian on N sites with alternating bond couplings is

    H_SSH = Σ_b (J + δ · (−1)^b) · (c_b^† c_{b+1} + h.c.)

where `c_b` is a fermion annihilation operator at site b, J is the mean hopping, and δ is the alternation amplitude. The intercell hopping (b even, say) is `J + δ`; the intracell hopping (b odd) is `J − δ`. The famous topological transition is at δ = 0, where the two band-touchings collide; for `|δ| > 0` the system carries a Z₂ winding number that flips sign under δ ↔ −δ.

Apply a Jordan-Wigner transformation, and the same Hamiltonian becomes an XX+YY spin chain with the same alternating bond couplings:

    H_SSH (spin form) = Σ_b J_b · (X_b X_{b+1} + Y_b Y_{b+1}) / 2,
    J_b = J + δ · (−1)^b.

This is exactly the framework's bond-coupling profile, with `J_b` the per-bond coupling. The framework's F71 mirror maps bond `b ↔ N − 2 − b`. Decompose J into its F71-palindromic and F71-anti-palindromic components:

    J_sym  = (J + F71(J)) / 2,
    J_anti = (J − F71(J)) / 2.

For the SSH pattern `J_b = J + δ · (−1)^b` on an even-length chain, `F71(J)_b = J_{N−2−b} = J + δ · (−1)^{N−2−b}`. For N even, `(−1)^{N−2−b} = (−1)^b`, so `F71(J) = J`: SSH on an even chain is F71-palindromic, J_anti = 0, δ lives entirely in J_sym. For N odd, `(−1)^{N−2−b} = −(−1)^b`, so the alternation flips its sign under F71-mirror, and δ ends up exactly in J_anti while J_sym = J (constant).

So the SSH alternation lives in J_anti precisely when the chain has odd length. The even-length chain hides the alternation inside J_sym. This asymmetry is itself worth a note: the framework's F71 mirror is the natural symmetry of an open chain with site-flip i ↔ N − 1 − i, and SSH's topology cares about whether the dimerization respects or breaks that mirror.

For the rest of this document we focus on **odd-length SSH chains**, where the alternation lives in J_anti and the framework's F92 / F100 predictions apply directly. The even-length case is a sibling question (the alternation sits in J_sym, and F92 has nothing to say about it; we would need a different framework axis to read the topology there).

---

## F92 prediction: the diagonal-block spectrum sees only J_sym = J

F92 ([the F92 bond anti-palindromic-J proof](../proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)) says: for any open-system Liouvillian built from an XY-class Hamiltonian with bond couplings J_b and uniform Z-dephasing, the diagonal-block eigenvalues depend only on the F71-palindromic part of J. The anti-palindromic part is invisible to the spectrum.

Applied to odd-length SSH chains: δ enters entirely in J_anti, so the F71-refined diagonal-block spectrum of the SSH Liouvillian under Z-dephasing is identical for any δ at fixed J. The diagonal-block (F71-compressed) eigenvalues on the (n, n+1) coherence blocks of a 5-site or 7-site SSH chain are the same as for a uniform chain with bond coupling J, regardless of how strong the alternation is. (Caveat, per F92's corrected scope: once δ ≠ 0 breaks the F71 mirror, those compressed-block eigenvalues are a basis-relative layer, NOT L's physical decay rates; the full-L spectrum does move with δ; the genuine δ-independent structure is the joint-popcount sectors.)

This is a striking prediction. SSH spectroscopy in solid-state physics is normally read in the closed-system band-structure picture, where δ is the gap-opening parameter and topology is everything. In the open-system Liouvillian picture under Z-dephasing, δ becomes spectrum-invisible. The topology is still there (it has not gone anywhere), but it has moved off the diagonal-block spectrum into the cross-block eigenvectors.

The framework prediction is testable directly: build the open-system Liouvillian for a 5-site SSH chain at several δ values (δ = 0, 0.2 · J, 0.5 · J), diagonalize each block-restricted L, and check that the eigenvalues match bit-exactly across δ at fixed J. The math is finite and exact; the numerical verification is a small script.

---

## F100 / F101 prediction: bond-mirror observables are exactly odd in δ

F100 ([the F100 proof](../proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md)) gives the complementary statement on the observable side. The closure-breaking coefficient c₁ and the F86c per-bond Q_peak observable both have bond-mirror deviations D(b) = c₁(b) − c₁(N−2−b) that are exactly odd functions of J_anti at fixed J_sym.

For SSH on an odd-length chain: δ is exactly J_anti, so the bond-mirror deviation D(b) is exactly odd in δ. At δ = 0 the deviation vanishes (the chain is uniform); as δ grows, D(b) grows linearly in δ to leading order, with no quadratic or higher even-power admixture. The bond-mirror observable becomes a clean linear readout of the SSH alternation amplitude.

This is the dual prediction to F92's invisibility statement. F92 says: the diagonal-block spectrum does not see δ. F100 says: the bond-mirror deviation vanishes without δ and is an exactly odd function of δ (at fixed J; its magnitude also carries J). Together they decompose the entire SSH structure into a "diagonal-block spectrum lives at J" half and an "observable deviation born of δ" half.

F101 is the γ-axis sibling: for non-uniform per-site Z-dephasing rates, the same bond-mirror deviation is exactly odd in the anti-palindromic γ component. This adds a second axis: if we tune the dephasing rates non-uniformly (γ_l alternating with the bond pattern, say), we can probe the cross-coupling between J and γ asymmetries.

---

## The topological reading

SSH's topology is famously a Z₂ winding number that flips sign at δ = 0. The two phases (δ > 0 and δ < 0) differ by where the dimerization sits: intercell or intracell bonds being the stronger ones. At the boundary between the two phases, edge soliton modes appear at zero energy. This is all textbook.

In framework language, δ = 0 is the F71-palindromic point. δ > 0 and δ < 0 are the two opposite directions along the F71-anti-palindromic axis. The Z₂ flip of the topological invariant maps onto the J_anti ↔ −J_anti operation that F100's oddness statement is about. The soliton modes at δ = 0 are not framework-invisible (the framework says δ-dependence is invisible to spectrum, not that δ = 0 is special); they appear at zero energy because the Hückel-style bands touch there, and that touching is independent of the framework's invisibility statement.

What the framework adds is the structural reading: the topology lives entirely along the F71-anti-palindromic axis. The Z₂ winding number is encoded in `sign(δ)`, which is sign(J_anti), which is exactly the variable F92 declares spectrum-invisible. Topology in SSH is spectrum-invisible by F92.

This is not a contradiction with the standard SSH picture; it is a relabelling. Standard SSH reads the topology in the closed-system band structure, where δ opens a gap. The framework reads the same structure in the open-system Liouvillian, where δ shifts the cross-block eigenvectors but leaves the diagonal-block eigenvalues alone. The two pictures see different sides of the same δ-axis: closed-system sees the gap opening, open-system sees the spectrum invariance plus the bond-mirror oddness.

A way to say it. SSH's topology and F92's invisibility are dual statements about the same anti-palindromic direction. The closed-system view localizes the topology in a band-touching at δ = 0; the open-system view localizes it in the F71-mirror anti-symmetric content. They are the same content seen through different windows.

---

## What numerical verification would look like

A clean follow-up pass would build the open-system SSH Liouvillian explicitly and check the three framework predictions:

1. **F92 invariance.** Build L_SSH(δ) for an odd-length chain (N = 5 or N = 7) under Z-dephasing at several δ values (δ = 0, 0.1·J, 0.3·J, 0.5·J). Compute the diagonal-block eigenvalues for each. Verify they agree bit-exactly across δ at fixed J.

2. **F100 oddness.** Compute the bond-mirror deviations D(b) = c₁(b) − c₁(N−2−b) at each δ value. Verify that D(b; +δ) = −D(b; −δ) bit-exactly, with no even-power admixture beyond floating-point floor.

3. **Topological correspondence.** Identify the SSH soliton modes (zero-energy edge states at δ ≠ 0) in the cross-block eigenvectors of L. Check that they live in the F92-invisible direction (no eigenvalue-shift signature, but eigenvector signature) and that their amplitudes are odd in δ.

The script lives naturally next to the existing `simulations/carbon/benzene_liouvillian_palindrome.py` from the 2026-05-22 benzene work. Polyacetylene is a chain rather than a ring, but the same open-system framework primitives apply with minor adjustments.

---

## Status and next steps

This document is the framework-side scaffolding. Three pieces are not yet in:

- The numerical verification script (deferred).
- The connection to the broader topological-insulator literature (this proof reads the framework correspondence narrowly; the wider topological-matter family beyond SSH is not addressed).
- The Peierls-distortion link to the May 22 Benzene Liouvillian Palindrome work: Peierls breaks the F1 spectrum in the ring case, and SSH dimerization is the chain version of the same physical mechanism. Whether the F1 break we saw at C₄ / C₆ under Peierls bond-dephasing has a precise analog in the SSH spectrum under bond-mode noise is an open thread.

What this document does add is the bridge: F92 + F100 read SSH directly, the alternation amplitude δ lives entirely in J_anti for odd-length chains, and the topology localizes in the framework's invisible-direction content. The framework predictions are sharp enough to be tested by a small open-system simulation, and the interpretive reading (topology as F92-invisibility) is the kind of cross-vocabulary identification the carbon thread has been accumulating since the 2026-05-17 seven-doc arc.

A note on the even-length case. SSH on even-length chains is a sibling question we have not closed here. The alternation sits in J_sym there, and F92 says nothing about it; some other framework axis would be needed to read the topology. Whether the even-length SSH topology has a framework-invisibility analog at all is an open structural question. The odd-length case is the one F92 catches cleanly.

---

## Anchors and cross-refs

- F92 (Tier 1 derived, the spectrum-side anti-palindromic-J invariance): [the F92 bond anti-palindromic-J proof](../proofs/PROOF_F92_BOND_ANTI_PALINDROMIC_J.md)
- F100 (Tier 1 derived, the observable-side oddness): [the F100 c₁/Q-peak mirror J-parity proof](../proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md)
- F101 (Tier 1 derived, the γ-axis sibling): [the F101 c₁ mirror γ-parity proof](../proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md)
- F71 (the chain mirror operator R: site i ↔ N − 1 − i, bond b ↔ N − 2 − b): documented in the F-formula registry, [F71 in the formula registry](../ANALYTICAL_FORMULAS.md), and inherited as the structural anchor of F86c / F91 / F92 / F100 / F101.
- Benzene parent: [the benzene Liouvillian palindrome](BENZENE_LIOUVILLIAN_PALINDROME.md) (Holstein passes, Peierls breaks F1; the chain version of Peierls is SSH)
- Three dephase letters: [the three benzene dephase letters](BENZENE_THREE_DEPHASE_LETTERS.md) (today's Klein-V₄ return visit, sibling thread)
- The candidate-doc list: [README](README.md) §What this folder will accumulate (POLYACETYLENE_F92_INHERITANCE was Tier-2 entry there).
