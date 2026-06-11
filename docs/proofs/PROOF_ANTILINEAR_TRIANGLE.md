# PROOF: the antilinear triangle (θ, conj, †: one Klein four-group, five proofs, one engine)

**Status:** Tier 1 derived (one-line algebraic identities, universal N; verified machine-exact at N = 1..3, the closure block exact at N = 2). The five corollary legs are existing, independently derived repository proofs; this document does not re-prove them, it exhibits the one engine they share.
**Date:** 2026-06-11
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [PROOF_PI_FACTORS_AS_R_TIMES_D](PROOF_PI_FACTORS_AS_R_TIMES_D.md): the mirror group D₄ = ⟨R, D⟩, the transpose mirror D, and the cube of characters this proof extends by one antilinear unit.
- [PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md) Lemmas A + B: the dagger leg.
- [PROOF_F113_COEFFICIENT_DERIVATION](PROOF_F113_COEFFICIENT_DERIVATION.md) Lemma C and [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2: the Hermitian-conjugacy leg.
- [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §4: the reversal kill, the word-level leg.
- [PROOF_PTF_CHIRAL_MIRROR_RATE_LAW](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) §5/§8 and [PROOF_F86B_OBSTRUCTION](PROOF_F86B_OBSTRUCTION.md): the K₁/K_b antiunitary mirror T = Σ₁∘conj, the conjugation leg.
- F114 and F118 in [ANALYTICAL_FORMULAS](../ANALYTICAL_FORMULAS.md): the transpose sign law and the cube of characters.

## Abstract

Operator space carries three involutions that every working physicist knows separately: the **transpose** θ(A) = Aᵀ, the **entrywise conjugation** conj(A) = Ā, and the **adjoint** †(A) = A†. Each is the product of the other two, † = θ∘conj = conj∘θ, so together with the identity they form a Klein four-group. What organizes them is a pair of ±1 gradings: the **linearity character** ℓ (is the map ℂ-linear or antilinear?) and the **multiplicativity character** m (does it preserve or reverse products?). The three maps occupy the three nontrivial cells: θ = (linear, reversing), conj = (antilinear, preserving), † = (antilinear, reversing). On a Pauli string σ all of this collapses to one number: θ(σ) = conj(σ) = (−1)^{n_Y(σ)}·σ and †(σ) = σ, because Y is the only Pauli that is both antisymmetric and imaginary.

The engine this proof banks is the **transport law**: for any operator H, Hermitian or not, the commutator superoperator L_H = −i[H, ·] transports through each map as

  **μ ∘ L_H ∘ μ = ℓ(μ)·m(μ) · L_{μ(H)},**

the sign being the *product of the two characters*. One sign comes from the −i (paid only by antilinear maps), one from the commutator's order (paid only by reversing maps); θ and conj each pay exactly one of the two, † pays both and so pays nothing. Five results in this repository, derived independently over five weeks, are five faces of this one line: the [F114](../ANALYTICAL_FORMULAS.md) transpose sign law (the θ face), the girth-ladder reversal kill (θ at word length), [F112](PROOF_F112_LINDBLAD_BIT_B_PI_BALANCE.md)'s dagger isometry and anti-Hermiticity (the † face at the Hilbert-Schmidt pairing), [F113](PROOF_F113_COEFFICIENT_DERIVATION.md)'s Hermitian conjugacy (the fixed-point collapse: on Hermitian operators θ = conj), and the K_b mode mirror T = Σ₁∘conj (the conj face, dressed). In the Pauli basis the triangle docks onto the [mirror group](PROOF_PI_FACTORS_AS_R_TIMES_D.md): θ *is* the D₄ reflection D, † is the antilinear unit 𝒦 (pure coefficient conjugation), and the full closure ⟨R, D, 𝒦⟩ is the **antilinear double** D₄ × Z₂ of order 16, eight of whose members are antiunitary, the supply the circle-thickened completion was waiting for.

The computational anchor is [`simulations/antilinear_triangle.py`](../../simulations/antilinear_triangle.py), self-validating.

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

The first of these is the engine of [F113 Lemma C](PROOF_F113_COEFFICIENT_DERIVATION.md) and of the [windowed converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md)'s Hermitian-conjugacy step (F117's leg factorization): for Hermitian H the ket leg, built from Hᵀ, is the entrywise conjugate of the bra leg, built from H, at the same indices, T̃^{(l⃗)}_{β⃗} = conj(T^{(l⃗)}_{β⃗}). Nothing about moments is special here; the moment identity is the collapse evaluated inside a trace.

**§3.4 The conj face: antiunitary intertwiners, K₁ and K_b.** A generator whose matrix is *real in the Pauli basis* commutes with the conj-transport, and an antilinear symmetry pairs eigendata as (λ, v) ↦ (λ̄, 𝒯v): the spectrum closes under conjugation with intertwined eigenvectors. That is the undressed form. The [K_b mode mirror](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) T = Σ₁∘conj is the dressed form: a signed diagonal Σ₁ (the K₁ shadow on the coherence block) composed with conjugation, chosen so that T commutes with the block generator *and* the perturbation direction; antilinearity is what lets the pairing land at λ̄ and hence at the mirrored mode index, K_b(b; mode k) = K_b(b; mode N+1−k) pointwise. The [K₁ trajectory identity](PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md) (typed as `ChiralMirrorTrajectoryClaim`) is the same vertex on the propagation side. In transport-law terms the conj face reads conj∘L_H∘conj = −L_{H̄}: for real H this is −L_H, the chiral flip these mirrors exploit, with the dissipator invariant.

## §4 The antilinear double, and where the circle completion lives

In the Pauli basis the triangle's generators dock onto the mirror group: D (= θ) is already a D₄ reflection, and 𝒦 (= †) is new. The closure is

  **⟨R, D, 𝒦⟩ ≅ D₄ × Z₂,  order 16,**

verified exactly at N = 2 by closure enumeration: eight members are linear (the mirror group), eight are antilinear, and the triangle {id, D, 𝒦, D𝒦} sits inside as the D-axis Klein subgroup. The direct-product structure deserves care, because the first guess is wrong: 𝒦 does *not* commute with R. R carries ±i phases (Y·X = −iZ per site), so 𝒦R𝒦 = conj(R) = D·R·D = 𝓕R: conjugation by 𝒦 swaps the two edge mirrors, ρ·F ↔ F·ρ, which is exactly the side swap an antiautomorphism must perform. The element that *is* central is **conj = D𝒦**, and physically it could not be otherwise: entrywise conjugation commutes with multiplication by the real matrix F = X^⊗N on either side, while the adjoint trades the sides. So the Z₂ factor is generated by conj, the double's centre is {I, 𝓕} × {I, conj}, and the group is the **antilinear double** of the mirror group: every mirror acquires an antiunitary twin. (The agent-built verifier caught this; the original scout had checked only the order and the antilinear count, not centrality.)

The bridge this builds, named and left open in [PROOF_PI_FACTORS_AS_R_TIMES_D](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5: the discrete completion ⟨r, d, h⟩ must thicken each letter-stabilizing transposition into a circle of rotations, and a circle wants reflections. The triangle supplies them with a clean division of labor, the **dial trio**:

  θ∘Ad_{R_z(φ)}∘θ = Ad_{R_z(−φ)},  conj∘Ad_{R_z(φ)}∘conj = Ad_{R_z(−φ)},  †∘Ad_U∘† = Ad_U for every unitary U.

The two ℓ·m = −1 vertices invert every dial (they are the O(2) reflections of each thickened circle); the † vertex commutes with all of them (the antiunitary that respects every rotation). The S₃ ⋉ D₄ completion itself stays open; what is settled here is which algebra its reflections come from.

## §5 What is ours and what is the home

The three maps are textbook, and antiunitary symmetry is Wigner's. What this document banks is the viewpoint, assembled independently and recognized afterwards: the two-character grading with the **product character as the transport sign**, the identification θ = D docking the triangle onto our mirror group, the antilinear double of D₄, and above all the recognition that five repository proofs, written weeks apart for different theorems (a sign law, a moment kill, a balance, a coefficient formula, a mode mirror), were each holding one face of the same four-element group. None of the five proofs changes; they gain a shared parent. The typed claim is `AntilinearTriangleClaim` (Core/Symmetry, parents `MirrorGroupD4Claim`, `CommutatorDConjugationSign`, `LindbladBitBPiBalance`; the `ChiralMirrorTrajectoryClaim` edge is carried here in prose, across the layer boundary).

The verification anchor is [`simulations/antilinear_triangle.py`](../../simulations/antilinear_triangle.py): the V₄ table and gradings, the transport law for all four vertices on non-Hermitian H, the five legs each re-derived from the engine and cross-checked against its home formulation, the order-16 closure, and the dial trio. Every block raises on failure; the process exits 0 only if the whole ledger holds.

## §6 The qudit generalization: the Weyl-Heisenberg lattice (added 2026-06-11, same day)

The triangle above is a qubit statement: on a Pauli string θ(σ) = conj(σ) = (−1)^{n_Y}σ, †(σ) = σ. That sign is the d = 2 shadow of a lattice action the Weyl-Heisenberg algebra carries at every local dimension d.

Replace the Pauli basis by the **Weyl-Heisenberg operators** P_{a,b} = X^a Z^b (a, b ∈ Z_d), where X is the clock shift |x⟩ ↦ |x+1⟩, Z = diag(ω^x) with ω = e^{2πi/d}, and ZX = ωXZ. They span the d² operators per site and generalize {I, X, Y, Z} = {P_{0,0}, P_{1,0}, P_{1,1}, P_{0,1}} (d = 2, ω = −1, Y = XZ up to phase). The three involutions act on the labels with a **symplectic phase**:

  θ(P_{a,b}) = ω^{−ab} P_{−a, b},  conj(P_{a,b}) = P_{a, −b},  †(P_{a,b}) = ω^{ab} P_{−a, −b},

verified machine-exact at d = 2, 3, 4, 5, with † = θ∘conj on the labels (the Klein four-group). The transport law μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is **basis-free and holds verbatim at every d** (the −i and the commutator order know nothing of the dimension; verified machine-zero d = 2..5 on non-Hermitian H).

Two things become clear at once. First, **the qubit (−1)^{n_Y} is the degeneration of the symplectic phase**: at d = 2, ω = −1 and the label flip a ↦ −a is trivial (every label is its own inverse), so θ(P_{a,b}) = (−1)^{ab}P_{a,b}, and (−1)^{ab} = −1 only on (1,1) = Y. The textbook "only Y is both antisymmetric and imaginary" is ω^{ab} evaluated at the one label where a and b are both nonzero, at the one ω where the flip collapses. Second, **for d > 2 the involutions genuinely move the labels** (a ↦ −a is nontrivial for a ≠ 0, d/2): the triangle is no longer a sign on each operator but a reflection of the Z_d × Z_d lattice dressed by the symplectic phase, and the Klein four-group {id, θ, conj, †} acts as the sign-flip group {(±a, ±b)}.

This dovetails with the [qudit mirror group](PROOF_QUDIT_PARTIAL_PALINDROME.md) (F121's Z_d ≀ Z₂): the mirror generator Π_d = ρᵀ·Shift^{⊗N} is the lattice **translation** a ↦ a + 1 composed with θ, while θ and conj are its **reflections**, so ⟨Π_d, D⟩ acts as the symmetry group of the Weyl-Heisenberg lattice, with the qubit D₄ = ⟨R, D⟩ its d = 2 cell.

**Open (the arc this opens, named for re-entry):** the **antilinear double** ⟨Π_d, D, 𝒦⟩ over general d (F119's D₄ × Z₂ is the d = 2 case), and whether the clock Z_d, as d → ∞, becomes the continuous rotation circle the [S₃ ⋉ D₄ completion](PROOF_PI_FACTORS_AS_R_TIMES_D.md) §5 must thicken: the discrete phase-space lattice Z_d × Z_d limiting to the continuous torus, the mirror group's two axes (local dimension and rotation) becoming one.

The verification anchor is [`simulations/qudit_mirror_group_family.py`](../../simulations/qudit_mirror_group_family.py) (self-validating: the symplectic-phase label action d = 2..5, the d = 2 collapse to (−1)^{n_Y}, the universal transport law, and the d > 2 label permutation).
