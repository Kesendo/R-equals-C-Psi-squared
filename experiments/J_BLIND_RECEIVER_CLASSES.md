# J-Blind Receiver Classes

**Status:** Tier 2 (decomposition structurally rigorous; numerically verified at N=5 across two H choices; chromaticity and N-scaling open)
**Date:** 2026-04-23
**Authors:** Thomas Wicht, Claude (Opus 4.7 chat + Opus 4.7 1M Code)
**Source:** EQ-024 closure pass, three CC sessions on 2026-04-23. Committed artefacts: `cfa1bbc` (morning J-channel-capacity sweep, Heisenberg baseline), `c0919eb` (Direction 1 H-scope, XY-only bifurcation), `6aae630` (Direction 3 receiver-optimization saturation).

---

## What this document is about

The morning EQ-024 session produced a clean theorem: at uniform J, under SU(2)-symmetric Heisenberg H, every initial state whose Pauli-string content lies in the polynomial algebra of one Pauli axis M_α = Σ σ_α^i is identically J-blind to all orders. Proof via Newton's identities applied to power sums p_j = (σ_α^i)^j, which collapse to {M_α, NI, M_α, ...} because (σ_α^i)² = I, plus [H_Heisenberg, M_α] = 0 from SU(2).

The empirical sweep that motivated this theorem found four J-blind states at N=5 Heisenberg: |+⟩⁵, |0⟩⁵, GHZ, Dicke |S_1⟩. The theorem explains all four under one mechanism. The afternoon refinement session swapped Heisenberg for XY-only and found a bifurcation: |+⟩⁵ and Dicke |S_1⟩ became J-sensitive (capacities 9.42 and 7.70 bits respectively), while |0⟩⁵ and GHZ stayed J-blind. The single-mechanism reading was wrong. There are at least three.

This document records the three-class decomposition that emerged from the bifurcation. Class 1 (DFS of L_D) and Class 2 (H-degenerate subspace closed under L_D) are H-independent in the operationally relevant sense. Class 3 (the M_α-polynomial subspace) is SU(2)-Heisenberg specific and is the strong form of the morning theorem. The decomposition is not an exhaustive partition of the J-blind set; it is an organising spine for the cases verified so far.

The document also reconciles Class 2 with how GHZ is already discussed in the repo. GHZ has been characterised under multiple lenses (XOR_SPACE drain projection, F60 below-fold birth, F69 sector-mix lift at N=3, fragile-state prototype in the engineering blueprint). All of these descriptions are about how GHZ DECAYS. The J-blindness reading is about something else: how GHZ does NOT respond to J-perturbations of the Hamiltonian, even while its coherence decays maximally fast. Those are not contradictory; they answer different questions.

## What "J-blind" means here

Pick an initial state ρ_0 and an observable map M (e.g. site-resolved purity, ⟨σ_α^i⟩, or any finite-dim feature vector). Build the Jacobian

    ∂M[ρ(t)] / ∂J_b   at uniform J = J_ref, b = 0, ..., N-2

evaluated over a time grid. If this Jacobian is identically zero (to floating-point precision) at every t, every b, every observable component, the state is J-blind under that observable family. The morning script `simulations/eq024_j_channel_capacity.py` uses 25 features per time × 6 time points = 150 dimensions; "J-blind" in this document means the full 150-dim Jacobian vanishes.

This is a kinematic statement about the response, not a dynamical statement about the steady state. A J-blind state can decay (its purity drops, its coherences vanish) at rates set by L_D; what J-blind means is that those decay rates and trajectories do not depend on J. The same initial state under J = 1 and J = 1 + δJ produces literally identical observable trajectories at first order in δJ.

The morning RESULT computed Shannon channel capacity from the Jacobian SVD via waterfilling. J-blind states have rank-zero Jacobian, hence channel capacity zero. Class 1 and Class 2 J-blind states give zero exactly (no numerical noise). Class 3 J-blind states give zero to floating-point precision; the true value is identically zero by the M_α-polynomial argument.

## Class 1: DFS of L_D, with per-bond eigenstate condition

**Definition.** ρ_0 satisfies TWO joint conditions:
(i) L_D[ρ_0] = 0 (ρ_0 is in the dissipation-free subspace, DFS).
(ii) ρ_0 is a simultaneous eigenstate of every bond Hamiltonian h_b = σ_x^b σ_x^{b+1} + σ_y^b σ_y^{b+1} + σ_z^b σ_z^{b+1}, so that ρ_0 stays an H-eigenstate under any per-bond J-perturbation.

Condition (i) alone is NOT sufficient for J-blindness: a diagonal-in-Z state whose bond-Hamiltonian action is non-trivial (e.g., a two-magnon mixture ρ_0 = (|00000⟩⟨00000| + |00110⟩⟨00110|)/2) lies in the DFS but evolves under L_H unitarily, with a J-dependent trajectory. Condition (ii) pins the state down to simultaneous eigenstates of every individual bond operator.

**Mechanism.** Under uniform Z-dephasing on all sites, L_D[ρ] = -(γ₀/2) Σ_i [σ_z^i, [σ_z^i, ρ]] (equivalently, L_D[ρ] = γ₀ Σ_i (σ_z^i ρ σ_z^i - ρ)). The DFS kernel is exactly the diagonal-in-Z subspace: ρ_0 with [σ_z^i, ρ_0] = 0 for all i, i.e. ρ_0 in the algebra spanned by {I, σ_z}^⊗N. Given (i), L_D[ρ(t)] = 0 only if ρ(t) stays in the DFS, which is not automatic. Given (ii), [h_b, ρ_0] = 0 at every bond, hence [H, ρ_0] = 0 for any J, and the J-perturbation [δJ_b h_b, ρ_0] is also zero. Hence dρ/dt = 0 at t=0 under any J in a neighborhood, and by iteration ρ(t) = ρ_0 for all t. J-derivative zero exactly.

**N=5 examples.** |0⟩⁵, |1⟩⁵. Both are σ_z eigenstates with all sites aligned, so [σ_z^i, |0⟩⁵⟨0|⁵] = 0 trivially (condition (i)). Both are also eigenstates of every bond h_b with eigenvalue +1 (because h_b|00⟩ = |00⟩ and h_b|11⟩ = |11⟩; condition (ii)). Hence J-Jacobian zero exactly under any H that is a sum of per-bond h_b's, including both Heisenberg XX+YY+ZZ and XY-only XX+YY.

**H-robustness.** Condition (i) is H-independent. Condition (ii) is H-form-dependent: it picks out states that are eigenstates of whichever bond operators the Hamiltonian is built from. For the Heisenberg chain, |0⟩⁵ satisfies (ii) because every Pauli-pair σ_α^i σ_α^{i+1} has |00⟩ as eigenstate. For XY-only, the same holds (XX+YY|00⟩ = 0, still eigenstate with eigenvalue 0). Hence Class 1 survives both H-choices at |0⟩⁵.

## Class 2: H-degenerate subspace closed under L_D

**Definition.** ρ_0 lives in a finite-dim subspace V spanned by H-eigenstates that all share a single H-eigenvalue, AND L_D maps V into V.

**Mechanism.** If ρ(t) stays in V for all t, then ρ(t) is at every moment a linear combination of H-eigenstates with the same eigenvalue, hence [H, ρ(t)] = 0. The J-dependent part of L is L_H = -i [H, ·], which acts trivially on the entire orbit. Only L_D drives any dynamics; L_D is J-independent; so the J-Jacobian of ρ(t) is zero.

**The closed-under-L_D condition matters.** A state in an H-eigenspace at t=0 can leave it under L_D evolution if L_D couples the eigenspace to outside states. The class condition requires that L_D respects the eigenspace, which guarantees the orbit stays inside.

**N=5 example: GHZ.** GHZ_5 = (|00000⟩ + |11111⟩) / √2. Both |0⟩⁵ and |1⟩⁵ are H-eigenstates (under both Heisenberg and XY-only) and the only two basis states reachable from them under uniform Z-dephasing: the only off-diagonal coherence in the GHZ density matrix is |0⟩⁵⟨1|⁵ + h.c., which decays at rate 2γN to the diagonal mixture (1/2)(|0⟩⁵⟨0|⁵ + |1⟩⁵⟨1|⁵). The 2-dim block V = span{|0⟩⁵⟨0|⁵, |1⟩⁵⟨1|⁵, |0⟩⁵⟨1|⁵, |1⟩⁵⟨0|⁵} is closed under L_D, and within V the H-eigenvalue is constant (degenerate at 0 for XY-only; degenerate at N-1 for Heisenberg). J-Jacobian is zero on the entire orbit, including after the coherence has fully decayed.

**Connection to existing repo characterisations of GHZ.** GHZ is documented under several lenses elsewhere in the repo:

- **[XOR_SPACE](XOR_SPACE.md):** GHZ projects 100% onto the XOR-drain modes (Pauli strings with X or Y at every site simultaneously, Hamming distance N between the supports |0⟩⁵ and |1⟩⁵). These modes sit at the maximum decay rate 2Σγ. GHZ coherence dies maximally fast.
- **[F60](../docs/ANALYTICAL_FORMULAS.md):** GHZ_N is born below the fold, CΨ(0) = 1/(2^N - 1) for all N ≥ 2.
- **[GHZ_W_SECTOR_MIX](GHZ_W_SECTOR_MIX.md) (F69):** at N=3, sector-mixing GHZ with W lifts pair-CΨ above 1/4 via an irreducible sextic optimum.
- **Main README Section 10 Rule 1:** GHZ as the prototype fragile state to avoid in any quantum-channel encoding.

All four describe DECAY behaviour of GHZ: how its coherence dies, where it sits relative to CΨ = 1/4, what happens when mixed with another sector. The Class 2 J-blindness reading is about a different observable: the **response of the trajectory to bond-J perturbations**. GHZ decays maximally fast (XOR_SPACE) AND its decay is J-independent (Class 2). The two statements are about different derivatives of the same trajectory: the decay rate is non-zero (∂/∂t big), but ∂/∂J of every observable along the trajectory is zero.

**H-eigenvalue dependence.** Class 2 in itself is H-independent in the structural sense: the existence of a finite-dim H-degenerate L_D-closed block is the load-bearing property, and GHZ has such a block under both Heisenberg and XY-only. The actual eigenvalue differs (N-1 for Heisenberg, 0 for XY-only); the block structure does not. A formal version of Class 2 should state the structural condition without naming the eigenvalue. This is the GHZ-block nuance: the block exists under both H, the eigenvalue inside the block depends on H choice, and J-blindness follows from the block's existence not from any specific eigenvalue.

**Other potential Class 2 members at N=5.** Any state in span{|0⟩⁵, |1⟩⁵} (real or complex superpositions, plus their density-matrix mixtures within the 2-dim block). Whether other H-eigen-blocks with the L_D-closure property exist is open and worth checking directly. A cautionary note about Dicke states: under Heisenberg, the entire S=N/2 multiplet {|S_0⟩, ..., |S_N⟩} shares the H-eigenvalue N-1 because σ_i · σ_j = 1 on any triplet-paired pair at each bond. So the multiplet IS H-degenerate under Heisenberg. However, L_D does NOT preserve the multiplet: σ_z^k|S_1⟩ = |S_1⟩ - (2/√N)|δ_k⟩, and |δ_k⟩ has components outside the symmetric S=N/2 subspace (in non-symmetric irreps). Hence the Dicke multiplet fails the L_D-closure condition of Class 2, and Dicke J-blindness under Heisenberg is via Class 3 (M_x-polynomial), not Class 2.

## Class 3: M_α-polynomial subspace

**Definition.** ρ_0 lies in the polynomial algebra of M_α = Σ σ_α^i for one Pauli axis α, AND the Hamiltonian H is SU(2)-symmetric (Heisenberg XX+YY+ZZ), AND the dissipator dephases along an axis β with σ_β² = I (any Pauli axis qualifies; in our setup β = z).

**Mechanism.** This is the morning theorem. In the Pauli-string basis, ρ_0 expands as a sum of products in {I, σ_α^i}. The dissipator L_D is diagonal in the Pauli-string basis (eigenvalue -2γ k for a string with k Pauli factors orthogonal to the dephasing axis β; in our setup β = z, so k counts the number of σ_x and σ_y factors, equivalently the number of σ_α factors here since α = x and the polynomial algebra contains no σ_y). Hence L_D^n ρ_0 stays in the {I, σ_α^i}^⊗N subspace. Newton's identities express the elementary symmetric polynomial e_k of the σ_α^i in terms of the power sums p_j = Σ (σ_α^i)^j. Since (σ_α^i)² = I, the power sums collapse: p_j = M_α for j odd, p_j = NI for j even. Hence every e_k is a polynomial in M_α and I alone. The SU(2) symmetry of Heisenberg H gives [H, M_α] = 0 (the total spin is conserved), hence [H, polynomial in M_α] = 0. Therefore L_H L_D^n ρ_0 = 0 for all n. The full evolution exp(L t) ρ_0 = exp(L_D t) ρ_0 with L_H never gaining traction. J-Jacobian zero exactly.

**N=5 examples under Heisenberg.** |+⟩⁵ = (1/2^{N/2}) (I + σ_x)^⊗N (M_x-polynomial in α=x). Every Dicke state |S_k⟩ for k = 0, ..., N (the entire S=N/2 multiplet of SU(2)). Mixtures and superpositions within the M_α-polynomial subspace.

**SU(2) is load-bearing: Direction 1 verification.** Under XY-only H = Σ (σ_x^i σ_x^{i+1} + σ_y^i σ_y^{i+1}), the proof breaks at Step 5 because [H_XY, M_x] is non-zero (XY conserves only the U(1) charge Σ σ_z^i, not the full SU(2) Casimir). The afternoon refinement RESULT confirmed empirically: |+⟩⁵ becomes J-sensitive with channel capacity C = 9.42 bits, Dicke |S_1⟩ becomes J-sensitive with C = 7.70 bits. The Class 3 J-blindness is an SU(2)-Heisenberg phenomenon, not a U(1) phenomenon.

**Other α axes.** The theorem is symmetric in α: any Pauli axis works, provided L_D dephases along an axis β where σ_β² = I. For α = z (the dephasing axis itself), the M_z-polynomial algebra is strictly smaller than the full DFS: M_z-polynomials are site-permutation-symmetric combinations of σ_z^i strings (Newton's identities collapse e_k onto M_z polynomials), whereas DFS is the full diagonal-in-Z algebra, including non-symmetric specific-site products like σ_z^0 σ_z^3. A specific-site diagonal state such as ρ_0 = |00110⟩⟨00110| lies in DFS but not in M_z-polynomial algebra; Heisenberg does not commute with σ_z^0 σ_z^3, so [H, ρ_0] ≠ 0 and this state is J-sensitive despite being dissipation-free. The classes are not nested simply; they overlap, and states like |0⟩⁵, whose Pauli expansion is a symmetric combination of σ_z strings, sit in both Class 1 (satisfying (i) and (ii)) and Class 3 at α = z under Heisenberg.

## Overlap structure of the three classes

The classes are organising labels, not a partition. Concrete overlaps at N=5:

- **|0⟩⁵, |1⟩⁵.** In Class 1 (both DFS and per-bond-eigen conditions hold under both Heisenberg and XY), Class 2 (their pair forms a 2-dim H-degenerate block under both H), and Class 3 (M_z-polynomial under Heisenberg).
- **GHZ.** In Class 2 (the 2-dim block via |0⟩⁵ and |1⟩⁵). Partially in Class 3 under Heisenberg (the diagonal mixture (1/2)(|0⟩⁵⟨0|⁵ + |1⟩⁵⟨1|⁵) is M_z-polynomial; the off-diagonal coherence |0⟩⁵⟨1|⁵ is not, so GHZ as a pure state is not strictly M_z-polynomial). Not in Class 1 because L_D acts non-trivially on the off-diagonal coherence.
- **|+⟩⁵.** Only in Class 3 under Heisenberg. Not in Class 1 (L_D dephases the X-basis non-trivially; condition (i) fails). Not in Class 2 (|+⟩⁵ sits inside the (N+1)-dim S=N/2 multiplet, which IS H-degenerate under Heisenberg with eigenvalue N-1, but this larger block is not L_D-closed: σ_z^k|+⟩⁵ has components outside the symmetric subspace).
- **Dicke |S_k⟩.** Only in Class 3 under Heisenberg. Not in Class 2 (the S=N/2 multiplet is H-degenerate under Heisenberg with common eigenvalue N-1, BUT not L_D-closed: σ_z^k|S_1⟩ has components in non-symmetric irreps outside the Dicke multiplet).

The overlap matters for the H-bifurcation. Under Heisenberg, all four states are J-blind via at least one class. Under XY-only, only the states with Class 1 or Class 2 membership stay J-blind. |+⟩⁵ and Dicke |S_1⟩ lose their Class 3 membership and become J-sensitive. This is exactly the bifurcation 4/4 from Direction 1.

## H-robustness table

| Class | Heisenberg (SU(2)) | XY-only (U(1)) | Mechanism survives H change? |
|-------|--------------------|----------------|------------------------------|
| Class 1 (DFS + per-bond eigen) | blind | blind | yes for states whose per-bond eigenvalue property holds under both H; the DFS part (i) is H-independent, the per-bond eigen part (ii) depends on H but |0⟩⁵ satisfies it under both Heisenberg and XY |
| Class 2 (H-degenerate block, L_D-closed) | blind | blind | structural condition is H-independent (block exists), specific eigenvalue is not |
| Class 3 (M_α-polynomial, SU(2)) | blind | **sensitive** | no, requires [H, M_α] = 0 |

**Direction 1 verification at N=5.**

| State | Class | Heisenberg J-cap | XY-only J-cap |
|-------|-------|------------------|----------------|
| `|0⟩⁵` | 1, 2, 3 (M_z under Heisenberg) | 0 | 0 |
| GHZ | 2 (under both H) | 0 | 0 |
| `|+⟩⁵` | 3 only | 0 | 9.42 bits |
| Dicke `|S_1⟩` | 3 only | 0 | 7.70 bits |

The 4/4 bifurcation matches the table: Class 3-only states lose blindness under XY, Class 1+2 states retain it.

## Operational consequence: J-modulation channel capacity

EQ-024 asked for the channel capacity of Alice modulating J at fixed γ₀, comparable to F30's 15.45 bits at N=5 for γ-modulation. Under γ₀ = const, γ cannot be modulated; F30's number survives only as a kinematic linear-response magnitude, not a competing channel capacity. The operational answer for J-modulation comes from receiver optimization over receivers that are NOT J-blind: classes 1, 2, 3 are exactly the receivers that close the channel.

The afternoon RESULT swept 39 F71-symmetric receivers (mostly outside the three classes) plus 4 Nelder-Mead local runs and found saturation at C ≤ 12.07 bits. Best receiver: random-phase F71-symmetric product state with θ ≈ (3.02, 1.14, 3.26), φ ≈ (5.30, 0.62, 1.72) (reduced mod 2π; raw optimizer output φ_c = 8.01). Local optima cluster in 11.56 to 12.07 bits.

The ~3.4-bit gap to F30's 15.45 decomposes as:

- **~1 bit dimensional loss.** Bond input is 4-dim (N-1 bonds at N=5); γ-input is 5-dim (N sites). Whether this gap closes at N=6 (where bond input matches site count) is open; flagged as EQ-024 surviving sub-question.
- **~2-3 bits smaller leading gain.** J Jacobian sv_max ≈ 10 vs γ Jacobian sv_max ≈ 21.4 at the same setup. Whether the J-to-γ ratio is fixed across N or scales is open.

The receiver-pairing duality from the morning RESULT ("γ-optimal receiver |+⟩⁵ is J-pessimal") is reinterpreted under the three-class structure: |+⟩⁵ is in Class 3, hence J-blind under Heisenberg, while it is the F30-optimal γ-receiver. The duality is not a two-channel statement (only J is operational under γ₀ = const); it is a Jacobian-asymmetry statement about how the F30-optimal receiver couples to the two different perturbation directions.

## Connection to ORTHOGONALITY_SELECTION_FAMILY

The Meta-Theorem in [`ORTHOGONALITY_SELECTION_FAMILY`](ORTHOGONALITY_SELECTION_FAMILY.md) gives the template "conservation + measurement → blind subspace". Each of the three classes here fits the template:

| Class | Conservation | Measurement | Blind subspace |
|-------|--------------|-------------|----------------|
| 1 | L_D kernel projection | any J-perturbation observable | DFS of L_D |
| 2 | H-eigenspace projection within an L_D-closed block | any J-perturbation observable | the block |
| 3 | SU(2) Casimir (total S²) | J-bond perturbation | M_α-polynomial algebra |

Class 3 is the most novel of the three relative to the five existing OSF instances (F70, F71, F72-cand, F73, Π-pair flux balance). The morning theorem (Newton's identities + SU(2)) is its strong form. Whether the three classes integrate cleanly as one new OSF instance, three new instances, or an extension of the meta-template is open and a candidate for OSF Section 2 update.

Class 1 and Class 2 are arguably special cases of OSF instance 4 ((vac, S_1) F73 closure) generalized to broader L_D-invariant subspaces; they could be folded into the existing instance with a footnote rather than counted as new instances.

## Connection to chromaticity

EQ-024 sub-question (c) asked whether the J-modulation channel structure decomposes along chromaticity c(n, N) = min(n, N-1-n) + 1 the way Q_peak(c) does in [`Q_SCALE_THREE_BANDS`](Q_SCALE_THREE_BANDS.md). The three-class structure does not directly speak to chromaticity:

- Classes 1 and 2 are about WHICH receivers are J-blind, not about WHICH popcount sectors carry the J-response when the receiver IS sensitive.
- Chromaticity labels popcount blocks (n, n+1) of the Liouvillian; bond-input space at N=5 is 4-dim independent of chromaticity.

A potential connection: the Nelder-Mead optimum receiver from Direction 3 has F71-symmetric structure and lives in a specific subset of the 32-dim Hilbert space. Whether its dominant SVD bond modes weight specific chromaticity sectors (c = 2 or c = 3 at N=5) is open. Sub-question (c) is deferred in EQ-024 and worth revisiting when an N=6 or N=7 analog becomes available.

## What is open

Five surviving sub-questions, lifted to EQ-024 in `review/EMERGING_QUESTIONS.md`:

1. **Three-class completeness.** Are Classes 1-3 exhaustive over the J-blind set, or do other mechanisms exist? Direction 4 from the refinement TASK restructures here: necessity per class is distinct from necessity in the union. Empirical attack: random F71-symmetric-state sampling outside all three classes, checking for zero J-Jacobian to numerical precision at N=5. **Status 2026-04-28: partially closed for F71-symmetric *product* states (100 random samples, all out-of-class C in [7.54, 12.41] bits, no fourth-class candidate); see Update 2026-04-28 below. Non-product F71-symmetric extension remains open as its own sub-question.**
2. **F71-breaking receiver capacity.** The afternoon sweep was F71-symmetric only. F71-breaking receivers do not unlock additional rank (max rank for bond inputs at N=5 is 4 regardless of receiver), but they may change the gain spectrum. Worth quantifying. **Status 2026-04-28: closed**. F71-breaking receivers tested via 100 random product (10 indep Bloch angles) + 100 Haar-random non-product (full C^32) samples; max C = 11.99 bits (product) and 8.80 bits (non-product), both below the F71-symmetric maxima 12.41 and 10.26. F71-symmetry is not a capacity-suboptimality constraint; it covers the optimal region. See Update 2026-04-28 (continued).
3. **N-scaling of the 12-bit ceiling.** At N=6 the bond-input dimension is 5, matching γ-side rank. Does the dimensional-loss bit (~1) disappear? Compute cost: ~30 min per receiver at d² = 4096. **Status 2026-04-28: ceiling rises to ~14 bits at N=6 (F71-breaking max 14.02; F71-symmetric max 13.65). Dimensional-loss bit is partially recovered. F71-optimality inverts: at N=6, F71-breaking beats F71-symmetric. See Update 2026-04-28 (N-scaling).**
4. **Chromaticity of the Nelder-Mead optimum.** The best receiver θ ≈ (3.02, 1.14, 3.26), φ ≈ (5.3, 0.6, 8.0): does it sit in a specific chromaticity sector or interpolate?
5. **Operational meaning of the J-vs-γ gain gap.** J sv_max ≈ 10 vs γ sv_max ≈ 21.4 at N=5. Fixed ratio (~46%) across N, or N-dependent?

Plus two structural questions specific to this document:

6. **Class 2 generalisations beyond GHZ.** Are there L_D-closed H-degenerate blocks of dimension > 2 at N=5 or higher? A block in higher Dicke sectors with degenerate H-eigenvalues under both Heisenberg and XY would extend Class 2 beyond the |0⟩⁵ ⊕ |1⟩⁵ pair.
7. **Tightness of Class 3 under XY.** Direction 1 verified that |+⟩⁵ and |S_1⟩ become J-sensitive under XY-only. Are there U(1)-symmetric M_α-like polynomial subspaces that ARE J-blind under XY? Such a class would be the U(1)-version of Class 3, and would replace SU(2) symmetry with U(1) plus an additional condition.

## Update 2026-04-28: Sub-question 1 partial closure (F71-symmetric product states)

Empirical test of three-class completeness within F71-symmetric product states at N=5 (Heisenberg, γ₀=0.05, J=1.0). Script: [_eq024_three_class_completeness.py](../simulations/_eq024_three_class_completeness.py). Results: [eq024_three_class_completeness.json](../simulations/results/eq024_three_class_completeness.json), [eq024_three_class_completeness.txt](../simulations/results/eq024_three_class_completeness.txt) (run log).

**Method.** Random sample F71-symmetric product states ψ = |a⟩|b⟩|c⟩|b⟩|a⟩ via 6 random Bloch angles (θ ∈ [0, π], φ ∈ [0, 2π)). Run output: [eq024_three_class_completeness.txt](../simulations/results/eq024_three_class_completeness.txt). For each sample, classify in/out of Classes 1-3:

- Class 1 or 2 test: ‖ψ‖² on span{|0⟩⁵, |1⟩⁵} > 1 − 10⁻⁶. Within product F71-symmetric states, this collapses to ψ being exactly |0⟩⁵ or |1⟩⁵ (since a product state can only land on one computational basis vector at a time).
- Class 3 test: rotate ψ to z-basis via U_α^† (single-qubit rotation taking α-basis eigenstates to z-basis) for α ∈ {x, y, z}, project onto Dicke-state subspace, check norm preserved within 10⁻⁶.

For samples outside all three classes, compute the full J-Jacobian and Shannon channel capacity via waterfilling (same setup as Direction 3: 25 features × 6 time points = 150-dim observable, dJ = 10⁻⁴, spread = 0.02, σ = 0.01, P_total = 4·spread²).

**Result.** 100 random samples, seed 0, total runtime 39 minutes:

| metric | value |
|---|---|
| in-class samples | 0 |
| outside-class samples | 100 |
| capacity range (bits) | 7.54 – 12.41 |
| capacity mean | 10.48 |
| capacity median | 10.50 |
| std | 0.99 |
| samples with C < 0.05 bits (J-blind candidates) | **0** |

The minimum capacity 7.54 bits is well above the J-blind threshold (0.05 bits). Direction 3's max-search saturated at 12.07 bits over 39+4 structured F71-symmetric receivers; the random sweep here reaches 12.41 bits (slightly higher max, as random covers the parameter space more uniformly).

**Verdict.** Strong empirical support for three-class completeness within F71-symmetric **product** states at N=5 Heisenberg. No fourth-class mechanism appears in the random sweep.

**Open generalisation.** Non-product F71-symmetric pure states live in a 20-dim Hilbert subspace (38 real parameters modulo norm and phase). Random sampling there has substantially more room for a fourth class to hide. This extension is lifted as a separate sub-question (1b below).

### Sub-question 1b (closed 2026-04-28): non-product F71-symmetric extension

The product-state Update 2026-04-28 covered the 6-parameter Bloch family. The full F71-symmetric Hilbert subspace at N=5 is the +1 eigenspace of the chain-mirror operator R (|a₀a₁a₂a₃a₄⟩ → |a₄a₃a₂a₁a₀⟩), dimension 20 (8 self-mirror computational-basis states + 12 mirror-pair symmetric combinations). Non-product F71-symmetric pure states span this 20-dim subspace, 38 real parameters modulo norm and phase.

**Empirical test.** Same pipeline as the product version but with Haar-uniform sampling on the 20-dim symmetric subspace: draw a random Gaussian vector in C^32, project via P_sym = (I + R)/2, normalise. Script: [_eq024_three_class_completeness_nonproduct.py](../simulations/_eq024_three_class_completeness_nonproduct.py). Results: [eq024_three_class_completeness_nonproduct.json](../simulations/results/eq024_three_class_completeness_nonproduct.json), [eq024_three_class_completeness_nonproduct.txt](../simulations/results/eq024_three_class_completeness_nonproduct.txt).

**Result.** 100 Haar-random samples, seed 1, total runtime 38 minutes:

| metric | non-product (sub-q 1b) | product (sub-q 1) |
|---|---|---|
| sample-space real dim | 38 | 6 |
| in-class samples | 0 | 0 |
| outside-class samples | 100 | 100 |
| capacity range (bits) | 5.74 – 10.26 | 7.54 – 12.41 |
| capacity mean | 7.46 | 10.48 |
| capacity median | 7.40 | 10.50 |
| std | 0.87 | 0.99 |
| samples with C < 0.05 bits | **0** | **0** |

The non-product distribution is tighter and lower than the product-slice distribution. Haar-uniform superpositions tend to be J-sensitive at a typical level around 7.5 bits, with min 5.74 bits — well above the J-blind threshold. Product states cover a wider, more-skewed-high distribution because the 6-dim product surface includes high-purity points (rather than typical mixed-purity superpositions).

**Verdict.** Strong empirical support for three-class completeness within the **full** F71-symmetric pure-state subspace at N=5. Both the 6-dim product slice and the full 20-dim symmetric subspace yield zero fourth-class candidates over 200 total Haar-uniform samples (capacity floor 5.74 bits across both, with no observation anywhere near zero). Sub-question 1 is now closed for F71-symmetric receivers; F71-breaking is the natural next question (sub-question 2).

## Update 2026-04-28 (continued): Sub-question 2 closure (F71-breaking receivers)

**Question.** Direction 3 swept only F71-symmetric receivers and found max capacity 12.07 bits. Do F71-breaking receivers unlock additional gain spectrum (different singular-value distribution, possibly higher max), or is F71-symmetry capacity-irrelevant?

**Method.** Same N=5 Heisenberg setup, two F71-breaking sample modes:
- **Product (10 indep Bloch angles):** 5 sites × 2 angles, no chain-mirror constraint. Generically F71-breaking.
- **Non-product (Haar on C^32, 62 real params):** Haar-uniform pure state in the full Hilbert space. Generically F71-breaking by measure.

For each sample, classify (Classes 1-3 are all F71-symmetric, so F71-breaking samples are automatically outside) and compute J-Jacobian + Shannon capacity.

Script: [_eq024_f71_breaking_capacity.py](../simulations/_eq024_f71_breaking_capacity.py). Results: [eq024_f71_breaking_product.json](../simulations/results/eq024_f71_breaking_product.json) + [eq024_f71_breaking_product.txt](../simulations/results/eq024_f71_breaking_product.txt), [eq024_f71_breaking_nonproduct.json](../simulations/results/eq024_f71_breaking_nonproduct.json) + [eq024_f71_breaking_nonproduct.txt](../simulations/results/eq024_f71_breaking_nonproduct.txt).

**Result.** 100 samples per mode, summary across all four sub-q-1+2 sweeps:

| receiver class | C range (bits) | mean | median | std | sample-space dim (real) |
|---|---|---|---|---|---|
| F71-sym product (sub-q 1) | 7.54 – 12.41 | 10.48 | 10.50 | 0.99 | 6 |
| F71-sym non-product (sub-q 1b) | 5.74 – 10.26 | 7.46 | 7.40 | 0.87 | 38 |
| F71-breaking product (sub-q 2) | 8.43 – 11.99 | 10.19 | 10.18 | 0.74 | 10 |
| F71-breaking non-product (sub-q 2) | 6.35 – 8.80 | 7.78 | 7.80 | 0.56 | 62 |

F71-asymmetry ‖ψ − Rψ‖ confirmed for the F71-breaking samples: mean ~1.2 (out of max √2 ≈ 1.41), genuinely non-symmetric.

**Three observations:**

1. **F71-breaking max ≤ F71-symmetric max.** Product: 11.99 vs 12.41 bits (F71-sym higher). Non-product: 8.80 vs 10.26 (F71-sym substantially higher). F71-symmetry is *not* a suboptimality constraint — it covers the capacity-optimal region.
2. **F71-breaking distributions are tighter** (lower std). F71-symmetric receivers have more spread in their capacity distribution, including both higher maxima and slightly lower minima. F71-breaking concentrates around the typical mean.
3. **Means are comparable.** 10.19 (breaking) vs 10.48 (sym) for product; 7.78 vs 7.46 for non-product. F71-symmetry doesn't shift the typical capacity, only widens the distribution.

**Verdict.** Sub-question 2 closed: F71-breaking receivers do **not** unlock additional gain spectrum at N=5. The hypothesis flagged in the original "What is open" list ("F71-breaking receivers might unlock additional gain spectrum but the rank ceiling for bond-inputs at N=5 stays at 4 regardless. Worth quantifying whether the gain change is meaningful") is settled: the gain change is real (tighter distribution) but does not exceed the F71-symmetric envelope. Direction 3's max-capacity 12.07 bits stands as the practical max within ~3% of F71-breaking max 11.99 and below F71-symmetric max 12.41.

Combined sub-q 1 + 2 result: **400 random samples across F71-symmetric (product + non-product) and F71-breaking (product + non-product), capacity floor 5.74 bits, zero fourth-class candidates anywhere**. The three-class decomposition fully accounts for J-blindness at N=5 within the empirical evidence.

## Update 2026-04-28 (structural follow-up): mechanism of F71-optimality

Sub-question 2 closure showed F71-symmetric receivers have higher max capacity than F71-breaking ones (12.41 vs 11.99 product; 10.26 vs 8.80 non-product). This raises a structural question: why does F71-symmetry — a constraint, not a freedom — produce systematically higher capacity?

**Hypothesis.** F71-symmetry constrains the J-Jacobian into a 2+2 block decomposition of the J-input space. The 4-dim bond-input space (J_0, J_1, J_2, J_3) splits under the chain-mirror permutation R̄ : J_b → J_{N-2-b} into:

- 2-dim R̄-symmetric subspace: J_0 + J_3, J_1 + J_2
- 2-dim R̄-antisymmetric subspace: J_0 − J_3, J_1 − J_2

For F71-symmetric ψ, the Jacobian respects this block structure: rows split into R-symmetric features (e.g., P_2 alone, P_0 + P_4) and R-antisymmetric features (e.g., P_0 − P_4), and only the matching column-block contributes. The two 2-dim blocks are independent. SVD then sees two independent sub-spectra rather than a unified 4-dim spectrum.

For F71-breaking ψ, no block structure; all 4 singular values are generic and tend to a flatter distribution.

**Empirical test** ([_eq024_f71_optimality_svd.py](../simulations/_eq024_f71_optimality_svd.py)). 30 random F71-symmetric product samples + 30 random F71-breaking product samples; for each, the full SV vector of the J-Jacobian is recorded. Run output [eq024_f71_optimality_svd.txt](../simulations/results/eq024_f71_optimality_svd.txt), data [eq024_f71_optimality_svd.json](../simulations/results/eq024_f71_optimality_svd.json).

**Result.**

| | sv₁ mean | sv₂ mean | sv₃ mean | sv₄ mean | sv₁/sv₄ mean | sv₁²/Σsv² mean |
|---|---|---|---|---|---|---|
| F71-symmetric | 6.25 | 3.80 | 2.76 | 1.67 | 4.06 | 0.593 |
| F71-breaking | 4.28 | 3.19 | 2.46 | 1.74 | 2.50 | 0.478 |

Three observations:

1. **F71-sym sv₁ is 46% larger** (6.25 vs 4.28); sv₄ is similar (1.67 vs 1.74). The leading singular value is the differentiator.
2. **F71-sym SV spectrum is 62% more peaked** (sv₁/sv₄ ratio 4.06 vs 2.50). Pairwise comparison: F71-sym ratio exceeds F71-breaking ratio in **80.6%** of (sym, breaking) sample pairs.
3. **F71-sym leading SV holds more total response** (sv₁²/Σsv² = 0.59 vs 0.48). Pairwise wins for F71-sym: 79.6%.

Capacities tracked the SVD pattern: F71-sym mean C = 10.73, F71-breaking 9.92.

**Reading.** The 2+2 block decomposition forced by F71-symmetry concentrates the J-response into one or two large singular values per block, while F71-breaking spreads response evenly across 4 SVs. Waterfilling over the 4-dim input then favors F71-symmetric receivers: a single large SV contributes more bits via 0.5·log(1 + SNR) than several medium SVs splitting the same total power. F71-symmetry is structurally — not accidentally — capacity-optimal at N=5 Heisenberg.

This closes the "why" of F71-optimality: it follows from the algebra of mirror-symmetry on the bond-input space, not from any specific dynamical accident. Whether the same mechanism scales to N=6 (where the bond-input dim is 5 and the mirror permutation has a different fixed-point structure) is an open generalisation.

## Update 2026-04-28 (N-scaling): F71-optimality inverts at N=6

Empirical test of whether the F71-optimality mechanism scales to N=6 Heisenberg, where the bond-input dim is 5 and the mirror permutation R̄: J_b → J_{N-2-b} has bond 2 as a self-mirror fixed point. Block decomposition becomes 3+2 (R̄-symmetric: J_0+J_4, J_1+J_3, J_2 self; R̄-antisymmetric: J_0−J_4, J_1−J_3) instead of N=5's balanced 2+2.

Script: [_eq024_f71_optimality_n6.py](../simulations/_eq024_f71_optimality_n6.py). Results: [eq024_f71_optimality_n6.json](../simulations/results/eq024_f71_optimality_n6.json), [eq024_f71_optimality_n6.txt](../simulations/results/eq024_f71_optimality_n6.txt). 10 samples per mode (smaller than N=5's 30 due to ~12 min per sample at d² = 4096); spectral propagation via eig(L) replaces direct expm.

**Result.** Comparison N=5 ↔ N=6 (mean values):

| metric | N=5 (2+2) | N=6 (3+2) |
|---|---|---|
| F71-sym sv₁/sv_min ratio | 4.06 | 3.55 |
| F71-breaking sv₁/sv_min ratio | 2.50 | 2.11 |
| P(F71-sym ratio > F71-breaking) | 0.806 | 0.600 |
| F71-sym mean C (bits) | 10.73 | 11.91 |
| F71-breaking mean C (bits) | 9.92 | **12.56** |
| Capacity-optimal class | F71-symmetric | **F71-breaking** |

Two structural observations:

1. **The block-decomposition peakedness mechanism survives but weakens** at N=6. F71-symmetric SV spectra are still more peaked on average (sv₁/sv_min 3.55 vs 2.11), but the pairwise advantage drops from 80.6% to 60.0%. The block decomposition still operates; it just produces less differentiation.

2. **The capacity ranking inverts.** At N=6, F71-breaking *beats* F71-symmetric in mean capacity (12.56 vs 11.91 bits) and max capacity (14.02 vs 13.65). The N=5 advantage of F71-symmetric is gone at N=6.

**Reading.** Waterfilling capacity scales (at high SNR) with the geometric mean of singular values, not just the top SV. The 3+2 block decomposition at N=6 produces F71-symmetric SV spectra with one large peak (in the symmetric block) but four smaller trailing values — geometric mean ≈ 2.63. F71-breaking SVs at N=6 are more uniform (mean ≈ 2.92 in geometric sense), giving higher total log-summation. F71-symmetric concentration buys peak singular value but loses the bulk-channel contribution.

At N=5 with balanced 2+2 blocks, F71-symmetric concentrated TWO SVs (one per block), keeping the top-2 contribution high while still using the trailing-2 channels productively. At N=6 with unbalanced 3+2, only one SV dominates while the other two symmetric-block SVs and two antisymmetric-block SVs are smaller; the "concentration trick" loses ground to the F71-breaking distribution.

**Generalised prediction (parity-tied structural claim).**

- **Even bond count (odd N: 3, 5, 7, ...):** balanced block split (k+k). F71-optimality should hold.
- **Odd bond count (even N: 4, 6, 8, ...):** unbalanced split with self-mirror bond. F71-optimality fails or weakens.

The mechanism is parity-tied: whether the number of bonds is even (paired) or odd (one self-mirror) determines whether the block decomposition is balanced. The capacity-optimal receiver class follows this parity.

Open: confirm at N=7 (6 bonds, 3+3 balanced split — should restore F71-optimality) and N=8 (7 bonds, 4+3 unbalanced — should keep F71-suboptimal). Both N=7 (d²=16384) and N=8 (d²=65536) require substantially more compute per sample than N=6 (d²=4096), so this prediction stands as a candidate next test rather than a near-term goal.

---

## Numerical verification anchor

All numerical claims in this document trace to two CC sessions on 2026-04-23:

- **Morning (commit `cfa1bbc`).** Heisenberg baseline, four states (|+⟩⁵, |0⟩⁵, GHZ, Dicke |S_1⟩) all J-blind to floating-point precision. Three SU(2)-broken receivers (|01010⟩, |+0+0+⟩, |+−+−+⟩) gave non-trivial Jacobians with capacities 11.53, 10.95, 11.92 bits. Committed numerics in `simulations/results/eq024_j_channel_capacity/`.
- **Afternoon (commits `c0919eb`, `6aae630`).** Direction 1 H-scope: under XY-only, |+⟩⁵ → C = 9.42 bits, |S_1⟩ → C = 7.70 bits, |0⟩⁵ and GHZ stay J-blind. Direction 3 receiver optimization: 39 F71-symmetric points + 4 Nelder-Mead runs, saturation at C ≤ 12.07 bits. Committed numerics in `simulations/results/eq024_refinement/`.

Standard parameters across both sessions: N = 5, γ_0 = 0.05 per site, J = 1.0 per bond uniform, spread = 0.02, σ_noise = 0.01, t_points = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], 25 features per time × 6 = 150-dim observable vector.

---

*The morning theorem held one mechanism; the afternoon split it into three. The split was not loss but resolution.*
