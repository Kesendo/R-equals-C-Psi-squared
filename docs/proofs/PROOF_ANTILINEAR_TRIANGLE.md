# PROOF: the antilinear triangle (θ, conj, †: one Klein four-group, five proofs, one engine)

**Status:** Tier 1 derived (one-line algebraic identities, universal N; verified machine-exact at N = 1..3, the closure block exact at N = 2; §7 exact at N = 1..5 in `AntilinearTriangleTests`). The five corollary legs are existing, independently derived repository proofs; this document does not re-prove them, it exhibits the one engine they share.
**Date:** 2026-06-11
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [The mirror group factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md): the mirror group D₄ = ⟨R, D⟩, the transpose mirror D, and the cube of characters this proof extends by one antilinear unit.
- [F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) Lemmas A + B: the dagger leg.
- [F113](PROOF_F113_COEFFICIENT_DERIVATION.md) Lemma C and [F87 windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2: the Hermitian-conjugacy leg.
- [F87 windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4: the reversal kill, the word-level leg.
- [the chiral mirror rate law](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) §5/§8 and [the F86b obstruction](PROOF_F86B_OBSTRUCTION.md): the K₁/K_b antiunitary mirror T = Σ₁∘conj, the conjugation leg.
- F114 and F118 in [Analytical Formulas](../ANALYTICAL_FORMULAS.md): the transpose sign law and the cube of characters.

## Abstract

Operator space carries three involutions that every working physicist knows separately: the **transpose** θ(A) = Aᵀ, the **entrywise conjugation** conj(A) = Ā, and the **adjoint** †(A) = A†. Each is the product of the other two, † = θ∘conj = conj∘θ, so together with the identity they form a Klein four-group. What organizes them is a pair of ±1 gradings: the **linearity character** ℓ (is the map ℂ-linear or antilinear?) and the **multiplicativity character** m (does it preserve or reverse products?). The three maps occupy the three nontrivial cells: θ = (linear, reversing), conj = (antilinear, preserving), † = (antilinear, reversing). On a Pauli string σ all of this collapses to one number: θ(σ) = conj(σ) = (−1)^{n_Y(σ)}·σ and †(σ) = σ, because Y is the only Pauli that is both antisymmetric and imaginary.

The engine this proof banks is the **transport law**: for any operator H, Hermitian or not, the commutator superoperator L_H = −i[H, ·] transports through each map as

  **μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)},**

the sign being the *product of the two characters*. One sign comes from the −i (paid only by antilinear maps), one from the commutator's order (paid only by reversing maps); θ and conj each pay exactly one of the two, † pays both and so pays nothing. Five results in this repository, derived independently over five weeks, are five faces of this one line: the [F114](../ANALYTICAL_FORMULAS.md) transpose sign law (the θ face), the girth-ladder reversal kill (θ at word length), [F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md)'s dagger isometry and anti-Hermiticity (the † face at the Hilbert-Schmidt pairing), [F113](PROOF_F113_COEFFICIENT_DERIVATION.md)'s Hermitian conjugacy (the fixed-point collapse: on Hermitian operators θ = conj), and the K_b mode mirror T = Σ₁∘conj (the conj face, dressed). In the Pauli basis the triangle docks onto the [mirror group](PROOF_PI_FACTORS_AS_R_TIMES_D.md): θ *is* the D₄ reflection D, † is the antilinear unit 𝒦 (pure coefficient conjugation), and the full closure ⟨R, D, 𝒦⟩ is the **antilinear double** D₄ × Z₂ of order 16, eight of whose members are antiunitary, and two of the triangle's vertices, θ and conj, are the reflections of the mirror group's rotation dial.

One level up (§7), the transpose of the generator's matrix is θ again, and it carries the Lindbladian to the Lindbladian of the transposed system exactly when the generator is unital. Read through the trace pairing this is reciprocity: for a real H with Pauli-string jumps, preparing a real A and reading a real B gives what preparing B and reading A gives, whatever the couplings and the dephasing profile, and with a flux the swap lands on the system with the flux reversed. What trades places are the preparation and readout operators; moving an experiment to the other end of a chain is covered where the two coincide.

The computational anchor is [`simulations/antilinear_triangle.py`](../../simulations/antilinear_triangle.py), self-validating; §7 is pinned exactly in `AntilinearTriangleTests` (`compute/MirrorWorld.Tests`).

## §1 The triangle

Work on the operator space of N qubits. Define three maps:

  θ(A) = Aᵀ,   conj(A) = Ā,   †(A) = A† = (Aᵀ)‾.

All three are involutions, and any two compose to the third; {id, θ, conj, †} ≅ V₄. Two characters grade this group:

- **ℓ (linearity):** θ is ℂ-linear (ℓ = +1); conj and † are antilinear (ℓ = −1).
- **m (multiplicativity):** conj is an automorphism, conj(AB) = conj(A)conj(B) (m = +1); θ and † are antiautomorphisms, θ(AB) = θ(B)θ(A) (m = −1).

The pair (ℓ, m) identifies each element uniquely: id = (+,+), θ = (+,−), conj = (−,+), † = (−,−). So V₄ ≅ Z₂(ℓ) × Z₂(m), and the product character ℓ·m, which §2 shows is the transport sign, takes the value −1 exactly on θ and conj.

On a Pauli string σ with n_Y(σ) letters Y: σᵀ = (−1)^{n_Y}σ and σ̄ = (−1)^{n_Y}σ (Y is the only antisymmetric Pauli and the only imaginary one, and it is both), so σ† = σ. In the Pauli basis, where A = Σ c_σ σ:

  θ(A) = Σ c_σ(−1)^{n_Y}σ,  conj(A) = Σ c̄_σ(−1)^{n_Y}σ,  †(A) = Σ c̄_σ σ.

The first line says **θ is the mirror D** = diag((−1)^{n_Y}) of the [mirror group](PROOF_PI_FACTORS_AS_R_TIMES_D.md), as a ℂ-linear map on coherence space; the third says **† is the antilinear unit 𝒦**, pure coefficient conjugation; the second says conj = D∘𝒦. The triangle is the D-axis of the mirror group, thickened by antilinearity.

## §2 The transport law (the engine)

**Theorem.** For every H and every μ ∈ {id, θ, conj, †},

  μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)},   L_H = −i[H, ·].

*Proof.* Write μ(L_H(μ(ρ))) = μ(−i(H·μ(ρ) − μ(ρ)·H)). The scalar −i passes through μ as ℓ(μ)·(−i) (antilinear maps conjugate it). The products transport as μ(H·μ(ρ)) = μ(H)·ρ when m = +1 and = ρ·μ(H) when m = −1 (using μ² = id), so the bracket becomes [μ(H), ρ] when m = +1 and −[μ(H), ρ] when m = −1. Collecting, μ∘L_H∘μ = ℓ(μ)m(μ)·(−i)[μ(H), ρ]. ∎

The four instances:

  id: +L_H.   θ: **−L_{Hᵀ}**.   conj: **−L_{H̄}**.   †: **+L_{H†}** (the two signs cancel).

The same bookkeeping runs for any scalar-times-bracket superoperator; nothing used Hermiticity, dimension, or the Pauli structure. The dephasing dissipator, diagonal in the string basis with real rates, is fixed by all three transports, which is why every leg below extends from L_H to the full dephased Lindbladian without further argument.

## §3 The five legs

**§3.1 The θ face: F114, and the reversal kill at word length.** Specialize μ = θ to a Pauli string H = σ: θ∘L_σ∘θ = −L_{σᵀ} = (−1)^{n_Y(σ)+1}·L_σ, which is [F114](../ANALYTICAL_FORMULAS.md) verbatim (typed as `CommutatorDConjugationSign`). The same antiautomorphism applied *inside a trace* gives the word face: for Pauli factors A_1, …, A_j,

  Tr(A_1⋯A_j) = Tr((A_1⋯A_j)ᵀ) = Tr(θ(A_j)⋯θ(A_1)) = (−1)^{Σᵢ n_Y(A_i)} · Tr(A_j⋯A_1),

the [reversal kill](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) of the girth ladder: word reversal = transpose × (−1)^{n_Y}, a sign-reversing involution exactly at odd total y-parity. F114 is θ at word length 2 (a commutator); the kill is θ at word length j.

**§3.2 The † face: F112 Lemmas A and B.** At the Hilbert-Schmidt pairing ⟨A, B⟩ = Tr(A†B), the adjoint of a commutator superoperator is (L_H)\* = −L_{H†}: a one-line trace computation, or equivalently the transport law read through the pairing. For Hermitian H this is the anti-Hermiticity (L_H)\* = −L_H of [F112 Lemma B](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md). Lemma A is the same vertex one level up: the dagger map on superoperators, S ↦ S\*, is an *antilinear isometry*, so it conjugates conjugation-eigenvalues, (Π S Π⁻¹)\* = Π S\* Π⁻¹ sends the Π-eigenvalue λ to λ̄ and swaps the ±i eigenspaces norm-preservingly. Antilinearity is the load-bearing property in both lemmas, and it is the ℓ = −1 of the † vertex. The same skew/self-adjoint split under this pairing (L_H skew, the dephasing dissipator self-adjoint) is the structural floor under the [Absorption Theorem](PROOF_ABSORPTION_THEOREM.md)'s Rayleigh reading: imaginary parts from H alone, real parts from the dissipator alone.

**§3.3 The fixed-point collapse: F113 Lemma C and the F117 Hermitian conjugacy.** Because † = θ∘conj, the fixed-point set of any one vertex is exactly the locus where the other two vertices *agree*:

  H = H†  ⟺  Hᵀ = H̄,   A = Aᵀ  ⟺  Ā = A†,   A = Ā  ⟺  Aᵀ = A†.

The first of these is the engine of [F113 Lemma C](PROOF_F113_COEFFICIENT_DERIVATION.md) and of the [windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)'s Hermitian-conjugacy step (F117's leg factorization): for Hermitian H the ket leg, built from Hᵀ, is the entrywise conjugate of the bra leg, built from H, at the same indices, T̃^{(l)}_{β} = conj(T^{(l)}_{β}) (l and β are leg tuples, written without an arrow). Nothing about moments is special here; the moment identity is the collapse evaluated inside a trace.

**§3.4 The conj face: antiunitary intertwiners, K₁ and K_b.** A generator whose matrix is *real in the Pauli basis* commutes with the conj-transport, and an antilinear symmetry pairs eigendata as (λ, v) ↦ (λ̄, 𝒯v): the spectrum closes under conjugation with intertwined eigenvectors. That is the undressed form. The [K_b mode mirror](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) T = Σ₁∘conj is the dressed form: a signed diagonal Σ₁ (the K₁ shadow on the coherence block) composed with conjugation, chosen so that T commutes with the block generator *and* the perturbation direction; antilinearity is what lets the pairing land at λ̄ and hence at the mirrored mode index, K_b(b; mode k) = K_b(b; mode N+1−k) pointwise. The [K₁ trajectory identity](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) (typed as `ChiralMirrorTrajectoryClaim`) is the same vertex on the propagation side. In transport-law terms the conj face reads conj∘L_H∘conj = −L_{H̄}: for real H this is −L_H, the chiral flip these mirrors exploit, with the dissipator invariant.

## §4 The antilinear double, and the reflections of the rotation dial

In the Pauli basis the triangle's generators dock onto the mirror group: D (= θ) is already a D₄ reflection, and 𝒦 (= †) is new. The closure is

  **⟨R, D, 𝒦⟩ ≅ D₄ × Z₂,  order 16,**

verified exactly at N = 2 by closure enumeration: eight members are linear (the mirror group), eight are antilinear, and the triangle {id, D, 𝒦, D𝒦} sits inside as the D-axis Klein subgroup. The direct-product structure deserves care, because the first guess is wrong: 𝒦 does *not* commute with R. R carries ±i phases (Y·X = −iZ per site), so 𝒦R𝒦 = conj(R) = D·R·D = 𝓕R: conjugation by 𝒦 swaps the two edge mirrors, ρ·F ↔ F·ρ, which is exactly the side swap an antiautomorphism must perform. The element that *is* central is **conj = D𝒦**, and physically it could not be otherwise: entrywise conjugation commutes with multiplication by the real matrix F = X^⊗N on either side, while the adjoint trades the sides. So the Z₂ factor is generated by conj, the double's centre is {I, 𝓕} × {I, conj}, and the group is the **antilinear double** of the mirror group: every mirror acquires an antiunitary twin.

The bridge this builds leads to the rotation circle of [the mirror group factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5's addendum, and a circle wants reflections. The triangle supplies them with a clean division of labor, the **dial trio**:

  θ∘Ad_{R_z(φ)}∘θ = Ad_{R_z(−φ)},  conj∘Ad_{R_z(φ)}∘conj = Ad_{R_z(−φ)},  †∘Ad_U∘† = Ad_U for every unitary U.

The two ℓ·m = −1 vertices invert the R_z dial (they are the O(2) reflections of that circle; a real rotation, R_y, they fix, since θ∘Ad_U∘θ = conj∘Ad_U∘conj = Ad_Ū); the † vertex commutes with every rotation dial (the antiunitary that respects every rotation). What is settled here is which algebra the circle's reflections come from.

## §5 What is ours and what is the home

The three maps are textbook, and antiunitary symmetry is Wigner's. What this document banks is the viewpoint, assembled independently and recognized afterwards: the two-character grading with the **product character as the transport sign**, the identification θ = D docking the triangle onto our mirror group, the antilinear double of D₄, and above all the recognition that five repository proofs, written weeks apart for different theorems (a sign law, a moment kill, a balance, a coefficient formula, a mode mirror), were each holding one face of the same four-element group. None of the five proofs changes; they gain a shared parent. The typed claim is `AntilinearTriangleClaim` (Core/Symmetry, parents `MirrorGroupD4Claim`, `CommutatorDConjugationSign`, `LindbladBitBPiBalance`; the `ChiralMirrorTrajectoryClaim` edge is carried here in prose, across the layer boundary).

The verification anchor is [`simulations/antilinear_triangle.py`](../../simulations/antilinear_triangle.py): the V₄ table and gradings, the transport law for all four vertices on non-Hermitian H, the five legs each re-derived from the engine and cross-checked against its home formulation, the order-16 closure, and the dial trio. Every block raises on failure; the process exits 0 only if the whole ledger holds.

## §6 The qudit generalization: the Weyl-Heisenberg lattice (added 2026-06-11, same day)

The triangle above is a qubit statement: on a Pauli string θ(σ) = conj(σ) = (−1)^{n_Y}σ, †(σ) = σ. That sign is the d = 2 shadow of a lattice action the Weyl-Heisenberg algebra carries at every local dimension d.

Replace the Pauli basis by the **Weyl-Heisenberg operators** P_{a,b} = X^a Z^b (a, b ∈ Z_d), where X is the shift |x⟩ ↦ |x+1⟩ and Z = diag(ω^x) the clock operator, with ω = e^{2πi/d}, and ZX = ωXZ. They span the d² operators per site and generalize {I, X, Y, Z} = {P_{0,0}, P_{1,0}, P_{1,1}, P_{0,1}} (d = 2, ω = −1, Y = XZ up to phase). The three involutions act on the labels with a **symplectic phase**:

  θ(P_{a,b}) = ω^{−ab} P_{−a, b},  conj(P_{a,b}) = P_{a, −b},  †(P_{a,b}) = ω^{ab} P_{−a, −b},

verified machine-exact at d = 2, 3, 4, 5, with † = θ∘conj on the labels (the Klein four-group). The transport law μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is **basis-free and holds verbatim at every d** (the −i and the commutator order know nothing of the dimension; verified machine-zero d = 2..5 on non-Hermitian H).

Two things become clear at once. First, **the qubit (−1)^{n_Y} is the degeneration of the symplectic phase**: at d = 2, ω = −1 and the label flip a ↦ −a is trivial (every label is its own inverse), so θ(P_{a,b}) = (−1)^{ab}P_{a,b}, and (−1)^{ab} = −1 only on (1,1) = Y. The textbook "only Y is both antisymmetric and imaginary" is ω^{ab} evaluated at the one label where a and b are both nonzero, at the one ω where the flip collapses. Second, **for d > 2 the involutions genuinely move the labels** (a ↦ −a is nontrivial for a ≠ 0, d/2): the triangle is no longer a sign on each operator but a reflection of the Z_d × Z_d lattice dressed by the symplectic phase, and the Klein four-group {id, θ, conj, †} acts as the sign-flip group {(±a, ±b)}.

This dovetails with the [qudit mirror group](PROOF_QUDIT_PARTIAL_PALINDROME.md) (F121's Z_d ≀ Z₂): the mirror generator Π_d = ρᵀ·Shift^{⊗N} is the lattice **translation** a ↦ a + 1 composed with θ, while θ and conj are its **reflections**, so ⟨Π_d, D⟩ acts as the symmetry group of the Weyl-Heisenberg lattice, with the qubit D₄ = ⟨R, D⟩ its d = 2 cell.

**Open (the arc this opens, named for re-entry):** the **antilinear double** ⟨Π_d, D, 𝒦⟩ over general d (F119's D₄ × Z₂ is the d = 2 case), and whether the Weyl-Heisenberg clock operator Z_d (not the F95 clock), as d → ∞, becomes the continuous rotation circle of [the mirror group factorization](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5's addendum: the discrete phase-space lattice Z_d × Z_d limiting to the continuous torus, the mirror group's two axes (local dimension and rotation) becoming one.

The verification anchor is [`simulations/qudit_mirror_group_family.py`](../../simulations/qudit_mirror_group_family.py) (self-validating: the symplectic-phase label action d = 2..5, the d = 2 collapse to (−1)^{n_Y}, the universal transport law, and the d > 2 label permutation).

## §7 One level up: the generator's transpose, and the two ends of a channel

The triangle acts on operators. Superoperators carry the same three maps one level up, and the θ vertex there answers, in algebra, a question this repository has carried in words since its beginning, about the two ends of a channel: what is prepared, and what is read out. Does it matter which end does which?

**The sweep behind this section.** The F-registry in [Analytical Formulas](../ANALYTICAL_FORMULAS.md) holds no reciprocity statement in this sense; "detailed balance" occurs there only in thermal and T1 settings (F84, F113, F120, F137), and [KMS and detailed balance](../KMS_DETAILED_BALANCE.md) names the detailed-balance conditions that carry an antiunitary time reversal without working one out. The word reciprocity occurs in `docs/proofs/` only as Frobenius reciprocity and in `experiments/` only in the photographic sense, exposure against brightness ([K-dosimetry](../../experiments/K_DOSIMETRY.md), [Beer-Lambert](../../experiments/BEER_LAMBERT_BREAKDOWN.md), [the cusp's slowing](../../experiments/CRITICAL_SLOWING_AT_THE_CUSP.md)). The complex symmetry this section reads is used widely as a tool (listed at the end). Every transfer in `experiments/` with an asymmetric layout ran in one direction ([the star topology](../../experiments/STAR_TOPOLOGY_OBSERVERS.md), [the QST bridge](../../experiments/QST_BRIDGE.md), [the relay protocol](../../experiments/RELAY_PROTOCOL.md)), apart from a reversed pair in the scaling curve that reading 3 takes up; [the XOR space](../../experiments/XOR_SPACE.md) leaves open whether A→B and B→A excite different partners, which would make the channel directional, [what's inside the windows](../../experiments/WHATS_INSIDE_THE_WINDOWS.md) lists the reverse direction as untested, and [the bridge closure](../../experiments/BRIDGE_CLOSURE.md) is the null result for a channel at J = 0. In `hypotheses/` the question is asked in words: [the quantum transistor](../../hypotheses/MEDIATOR_AS_QUANTUM_TRANSISTOR.md) argues the reversed direction from the A↔B permutation, [the other side](../../hypotheses/THE_OTHER_SIDE.md) §17 and [resonance, not channel](../../hypotheses/RESONANCE_NOT_CHANNEL.md) carry the search for a two-way bridge, and [the bidirectional bridge](../historical/THE_BIDIRECTIONAL_BRIDGE.md) is its historical form; the OpenArcs entry `the_gate_that_does_not_gate` records it as the question the March work named its starting point. `fw.Confirmations` holds no flight with sender and receiver swapped, the glossary nothing on reciprocity or on the generator's transpose, `docs/CAUGHT_ERRORS.md` the complex symmetry only in the fragile-bridge entry named at the end, and the typed layer no claim stating the transpose law or the reciprocity.

**The map.** Write the matrix of a superoperator S in the matrix-unit basis E_ij = |i⟩⟨j|, with the entry S(E_ij)[k,l] at row (k,l) and column (i,j), and let Sᵀ be the superoperator whose matrix is the transpose. Sᵀ is linear and reverses composition, (ST)ᵀ = TᵀSᵀ: the (linear, reversing) cell, θ one level up. It factors through the other two vertices, Sᵀ = conj∘S\*∘conj with S\* the Hilbert-Schmidt adjoint of §3.2 (its matrix is the conjugate transpose, and conj conjugates every entry back, the E_ij being real), so it is the † face of §3.2 composed with the conj face of §3.4.

**Theorem (the transpose law).** Let L = L_H + Σ_k r_k D_{J_k}, with D_J(ρ) = JρJ† − ½{J†J, ρ}, any H, any jumps J_k and any real rates r_k. Then

  **Lᵀ = L_{Hᵀ} + Σ_k r_k D_{J_kᵀ}**

exactly when L is unital, L(I) = Σ_k r_k [J_k, J_k†] = 0: the transpose of the generator is then the generator of the transposed system. Normal jumps (J†J = JJ†) are unital term by term.

*Proof.* Compare entries. L_H(E_ij)[k,l] = −i(H_ki δ_jl − δ_ki H_jl), and L_{Hᵀ}(E_kl)[i,j] = −i((Hᵀ)_ik δ_lj − δ_ik (Hᵀ)_lj) is the same number. The sandwich gives (J E_ij J†)[k,l] = J_ki·J̄_lj and (Jᵀ E_kl (Jᵀ)†)[i,j] = (Jᵀ)_ik·conj((Jᵀ)_jl) = J_ki·J̄_lj, for every jump with no condition. What remains are the anticommutators, and the transpose of X ↦ {M, X} is X ↦ {Mᵀ, X} by the same entry count: {K, ·} on one side with K = Σ_k r_k (J_k†J_k)ᵀ = Σ_k r_k J_kᵀJ̄_k, and {K′, ·} on the other with K′ = Σ_k r_k (J_kᵀ)†J_kᵀ = Σ_k r_k J̄_kJ_kᵀ. A map X ↦ {M, X} vanishes only for M = 0, so the two agree exactly when K = K′, whose complex conjugate reads Σ_k r_k J_k†J_k = Σ_k r_k J_kJ_k†. ∎

Unitality is also what keeps the transposed generator trace-preserving: the pairing identity below, at A = I, gives Tr(Lᵀ(Bᵀ)) = Tr(B·L(I)).

Every Pauli string is normal (it is unitary) and satisfies Pᵀ = ±P, the sign squaring away in D_P. So for Pauli-string jumps (dephasing along X, Y or Z at any site profile, depolarizing) **Lᵀ = L_{Hᵀ} + D**: the dissipator is its own transpose, and only H can tell the two ends apart. In general, once L is unital, the dissipator is its own transpose when every jump is symmetric up to a phase, or when the jumps come in transposed pairs at equal rates; such a pair (J, Jᵀ) is unital exactly when [J, J†] is real. σ⁻ and σ⁺ are such a pair: alone, σ⁻ is not unital and breaks the law; with σ⁺ at the same rate, the infinite-temperature bath, it holds; at unequal rates, a finite temperature, it breaks. Symmetric alone is not enough (a symmetric jump that is not normal, such as [[1, i], [i, −1]], breaks the law), and neither is real (a real, normal, non-symmetric jump such as a cyclic shift keeps the law but not Lᵀ = L). The transpose belongs to a basis, here the computational one; dephasing along a tilted axis such as (X + Y)/√2 is not its own transpose in this basis, and a common rotation about Z turns it into X dephasing.

**Corollary (reciprocity).** For all operators A and B and every t,

  Tr(B · e^{tL}(A)) = Tr(Aᵀ · e^{tLᵀ}(Bᵀ)),

and for symmetric A and B (Aᵀ = A, Bᵀ = B) this reads

  **Tr(B · e^{tL}(A)) = Tr(A · e^{tLᵀ}(B)).**

*Proof.* With M the matrix of S, Tr(B·S(A)) = Σ B_lk·M[(k,l),(i,j)]·A_ij, and Tr(Aᵀ·Sᵀ(Bᵀ)) = Σ (Aᵀ)_ji·M[(k,l),(i,j)]·(Bᵀ)_kl = Σ A_ij·M[(k,l),(i,j)]·B_lk, term by term the same sum. Apply this to S = e^{tL}, whose transpose is e^{tLᵀ}. ∎

Read physically, an experiment is a pair: A is what is prepared, B is what is read out. Its swap prepares B and reads A; for a positive readout operator that is not a state, that means preparing B/Tr B, and the swapped signal carries the factor (for a readout that is not positive the identity stands without a swapped experiment). A Hermitian operator is symmetric exactly when it is real, so A and B may be populations, the Bell state Φ⁺ = (|00⟩ + |11⟩)/√2, |+⟩ states, any real density matrix or real observable. A coherence such as |0⟩ + i|1⟩ is not symmetric, but Aᵀ = Ā is again a state, and the general identity pairs it with |0⟩ − i|1⟩.

**Three readings.**

1. *Real H: preparation and readout trade places.* A Hermitian H with Hᵀ = H is real, so with Pauli-string jumps Lᵀ = L: the generator is complex symmetric, and for real A and B, preparing A and reading B gives the signal of preparing B and reading A. Nothing about the layout enters: the couplings may be as uneven as one likes and every site may carry its own dephasing rate. The Heisenberg and XXZ chains of this repository, under Z dephasing, are this case. The complex symmetry also follows inside this document: §2 and §3.2 give θ∘L∘θ = −L_H + D = L\* for real H, so Lᵀ = conj∘L\*∘conj = †∘L∘†, and a Hermiticity-preserving L commutes with †.

2. *A flux: the swap lands on the reversed flux.* A Hermitian H that is not real (a hopping phase around a loop that no gauge removes) has Hᵀ = H̄. The corollary then pairs an experiment with its swap on the conjugate system: preparing A and reading B at flux +φ gives what preparing B and reading A gives at flux −φ. One fixed system can now tell the two ends apart; the swapped experiment still exists, on the mirrored system.

3. *Operators, and when positions follow.* The swap exchanges the preparation operator with the readout operator. Moving an experiment to the other end of a chain while keeping its habits (prepare over the rest of the chain in |0⟩, read with the rest ignored) is a different act, and the corollary covers it when the readout, restricted to what the dynamics can reach from the preparation, is the mirrored preparation. On the chains of this repository that is easy to read, because a number-conserving H under Z dephasing keeps every (ket, bra) excitation block apart.

   *One excitation, end to end.* Inside the one-excitation block, "end site b excited, the rest ignored" is exactly the projector on that excitation, so sending one excitation from one end site to the other is its own swap and gives the same signal both ways, whatever the couplings and rates. The same holds for any single-qubit state sent end to end and read at the far end site, since its coherence travels in the (0, 1) block, whose own generator is complex symmetric. Read at an interior site instead, the positional mirror is no longer the swap, and it differs.

   *A Bell pair.* Preparing Φ⁺ on the near pair over the rest in |0⟩ and reading Φ⁺ on the far pair with the rest ignored mixes the 0 and 2 excitation blocks. Two of the far pair's elements are transfers the corollary fixes, the |11⟩ population (|11 0…0⟩ carried to |0…0 11⟩) and the |00⟩⟨11| coherence; they agree from both ends. The |00⟩ population is different: in the two-excitation block it collects every state with the far pair empty, starting with the preparation's own state, so it holds a return probability, and the corollary does not tie that to the return probability of the mirrored experiment. The |01⟩ and |10⟩ populations and the |01⟩⟨10| coherence are not fixed either. So the Bell readout is not its own swap, and where it agrees from both ends a further symmetry does the work: at N = 4 the spin flip X^⊗N, which carries |1100⟩ to |0011⟩ and which an uneven longitudinal field breaks; on the XX chain (no ZZ term) the free fermions of Jordan-Wigner, since Z dephasing is quadratic in them and the one-body correlations ⟨c_i†c_j⟩ then evolve on their own under the one-excitation generator, the corollary fixing the far pair's total occupation and with it the |00⟩ population and the Φ⁺ fidelity, while |01⟩ and |10⟩ still differ one by one. On uneven Heisenberg and XXZ chains with a site between the pairs, and on an XX star (whose Jordan-Wigner strings break the closure), the Bell readout differs between the two ends.

   The N = 11 pull chain of [the scaling curve](../../experiments/SCALING_CURVE.md) is a Bell pair on a Heisenberg chain. Its rows "2:1 all (cascading pull)" and "2:1 reversed (push)" (in [`pull_principle.txt`](../../simulations/results/pull_principle.txt) "2:1 all pull" and "2:1 reversed", with six digits) are one chain run from its two ends, as the doc's own closing note says, and the sampled peak mutual informations differ, 0.120625 against 0.102039. Nonlinearity is not needed for such a difference: the far pair's Φ⁺ fidelity, a linear readout, can already differ through the populations the corollary does not fix. For the quantum transistor's "the coupling ratio J_AM/J_MB controls directionality" this cuts both ways: for one excitation or any single-qubit state sent from end site to end site over the static chain, the corollary rules out a direction bias in any readout of the far end site; for a Bell-pair payload or a pair receiver it does not decide. The same answer meets the XOR space's open question in the sense of the operator swap. A nonlinear readout (mutual information, CΨ) is outside the corollary altogether.

**Where the repository already stood on this.** The complex symmetry of reading 1 has worked for months as a rigidity tool, always derived from "H real symmetric":
- [the weight-one degeneracy](PROOF_WEIGHT1_DEGENERACY.md) Lemma 1, the transpose as a Krein form, τB†τ = B, which carries F50's Theorems A to D (typed in `ExceptionalCouplingSetClaim`);
- [the frozen divisor](PROOF_R90_FROZEN_DIVISOR.md) Lemma 6 (the cofactor determinant as a Gram determinant) and its isotropic eigenvector heading a Jordan chain, and the transpose parity that grades [the frozen band](PROOF_FROZEN_BAND_SO4.md)'s perturbation orders (Proposition 5.2);
- the complex-symmetric Kᵣ of [the missing phase](PROOF_MISSING_PHASE_SLOW_READOUT.md) (`MissingPhaseReadoutAlgebra`, `MissingPhaseSlowReadoutWitness`);
- the complex-symmetric reductions of [Route B at N = 4](PROOF_ROUTE_B_N4_SELF_FOLD.md) and [N = 6](PROOF_ROUTE_B_N6_UNFOLDING.md), and the measured max|L − Lᵀ| = 0.0 of [the N = 5 real-q diabolic](PROOF_N5_REAL_Q_DIABOLIC.md) (`RouteBN5RealQSemisimpleClaim`);
- the pencil M(z) = C + z·A₀ of [the exceptional couplings](../../experiments/THE_EXCEPTIONAL_COUPLINGS.md);
- the self-orthogonal coalescing eigenvector rᵀr → 0 of [the β-exotic genericity study](../../experiments/F89_BETA_EXOTIC_GENERICITY.md), the vᵀv gauge of the eigenvector holonomy (`EigenvectorHolonomy` in Core's numerics, `SeedHolonomyClaim`, `SeedHolonomyWitness`, and the Core.Tests pin `WeightCoherenceBlock_IsComplexSymmetric_TheVTvGaugePrecondition`), and the "complex-symmetric plus pseudo-Hermitian" rigidity of the F89 open arc;
- the pseudo-Hermiticity of iL in [the PT-symmetry analysis](../../experiments/PT_SYMMETRY_ANALYSIS.md) and in the Krein signatures of [the fragile bridge](../../hypotheses/FRAGILE_BRIDGE.md)'s exceptional point (certified in the 2026-09-25 entry of [the caught errors](../CAUGHT_ERRORS.md));
- the real-H step of [the pairing condition](../../experiments/THE_PAIRING_CONDITION.md), where for a real symmetric H the transpose upgrades the metric relation to a similarity, the same precondition under which F1 writes Π_Z(ρ) = ρᵀ·X^⊗N;
- single-particle, the bipartite conjugation S·conj(M)·S⁻¹ = M of [the edge block](PROOF_EDGE_BLOCK_DEFECTIVE_UNDER_PROFILE.md) at Δ = 0, whose Hermitian-hop form S·conj(M)·S⁻¹ = Mᵀ the open arc `site_resolved_vacuum_block` records with the reason that a path is a tree and its hop phases carry no flux: reading 2 seen from the one-excitation side.

What this section adds is the reading on the flow: the same Lᵀ = L says that preparation and readout can trade places. The Hilbert-Schmidt adjoint L\* = −L_H + L_D, written for Hamiltonian plus dephasing, already stands in [KMS and detailed balance](../KMS_DETAILED_BALANCE.md), which compares it with Π; the reciprocity is not drawn there, and whether θ∘L∘θ = L\* fills one of the Sá-Ribeiro-Prosen symmetry slots that document leaves open is not decided here.

Checked for adjacency beyond this list: F29's receiver/sender asymmetry on the star (about 1.7 at low partner rates, inverting at 0.17292; F29's table and the star topology's γ_A against γ_B scan) compares placements of γ read through CΨ, a positional and nonlinear comparison, so it neither tests nor contradicts the corollary. The textbook neighbour is closer than a resemblance. Take complex conjugation K in the computational basis as the time reversal: a real H is K-invariant, and for Hermitian A and B the corollary reads Tr(B·e^{tL}(A)) = Tr(Ā·e^{tL}(B̄)). That is the time-reversal form of quantum detailed balance at the infinite-temperature state I/d, the relation Onsager-Casimir reciprocity is drawn from, with the flux reversal as Casimir's field reversal; it is one of the antiunitary detailed-balance conditions the KMS document names, worked out at β = 0. The time reversal here is K, not the spin time reversal iσ_yK. Lorentz reciprocity in a lossy medium, broken by a magnetic bias as εᵀ(B) = ε(−B), is the same shape in electromagnetism.

The verification anchor for this section is `AntilinearTriangleTests` in `compute/MirrorWorld.Tests`, on integer inputs; every identity is compared exactly (== 0.0), every break is an exactly computed residual above 0.1:
- **The transpose law:** at N = 1, 2, 3 on a non-Hermitian H, with Z on every site, the all-X string and a string carrying Y and, from N = 2, a phased cycle (normal, neither symmetric nor unitary) as jumps. On the cycle, carrying the jumps back unchanged, conjugated or daggered instead of transposed breaks the law.
- **The condition:** an added σ⁻ breaks the law; σ⁻ with σ⁺ at equal rates keeps it, at unequal rates breaks it. The Lindbladian used preserves the trace and equals the world's dephasing dissipator at H = 0.
- **Complex symmetry:** Lᵀ = L for a dense real symmetric H at N = 4, and the flux ring breaks it.
- **Reciprocity, order by order:** Tr(B·Lᵖ(A)) = Tr(A·Lᵖ(B)) for p = 0..4 on a real H, broken by a Hermitian, non-real payload. On a dense real H (not number-conserving) the operator swap agrees while the positional mirror differs.
- **When positions follow:** on uneven Heisenberg and XXZ chains at N = 5 with unequal rates, through order t¹⁰: one excitation end to end and a |+⟩ agree from both ends (the excitation's signal first appears at t⁸), an interior readout differs, and the Bell pair differs (first at t⁶). The Bell pair agrees on the XX chain, with or without a longitudinal field, and at N = 4 under the spin flip; the field breaks the N = 4 agreement, and on the XX star the Bell pair differs while one excitation leaf to leaf still agrees.
- **The flux reversal:** on a three-site ring carrying a quarter flux, Hᵀ is the reversed ring, and the swap is exact against the reversed ring and broken against the same ring.
