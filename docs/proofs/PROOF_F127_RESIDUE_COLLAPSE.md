# PROOF: F127 by residue collapse, the sheet lattice, the core identity, and the mirror anchor

**Status:** derived chain over ℚ(i), machine-verified end to end (six exact gate scripts, total runtime about ten minutes; the one non-looped link is the textbook simple-pole change-of-uniformizer rule, §4). The overall F127 grade stays Tier1Candidate per the registry: the code-trust caveat (bespoke single implementations) applies to these gates as it does to the wall, and the two proofs discharge it for each other only partially, their trust surfaces being disjoint. The grid+CRT wall (527 deterministic items, 17 primes) remains in force as the independent certificate of the same theorem; this proof supersedes it as the primary derivation and explains it. Verifier scripts are being promoted with this document; the C# witness extension (`inspect --root crosstriple`) is owed.
**Date:** 2026-07-14
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- The F127 registry entry in [Analytical Formulas](../ANALYTICAL_FORMULAS.md) (~§F127): the theorem statement and the grid+CRT certificate this proof replaces as primary.
- [experiments/F89_SEED_EXISTENCE_REDUCTION.md](../../experiments/F89_SEED_EXISTENCE_REDUCTION.md), sections "The cross-triple orthogonality" and "The variety identity, proved over ℚ(i): the grid+CRT wall": where 𝔉 comes from (the twinning of resonant levels) and how the wall was built.
- [`simulations/cross_triple_orthogonality.py`](../../simulations/cross_triple_orthogonality.py): the committed definition of `cross_form` (𝔉) and the ℚ-proof of the mirror specialisation, this proof's anchor (§6).
- [`simulations/residue_assembly_close.py`](../../simulations/residue_assembly_close.py): the committed w₃-eliminated frame; its structure results (simple poles, Laurent window) are load-bearing inputs (§5).
- [`simulations/halfangle_residue_proof.py`](../../simulations/halfangle_residue_proof.py): the committed pre-elimination residue results; §4's transport reconciles the two frames.
- [reflections/ON_LEAVING_THE_CIRCLE.md](../../reflections/ON_LEAVING_THE_CIRCLE.md): the fourth and fifth visits, which predicted this document ("inside there is something simple and fragile, a one-line reason of the same family as: the length cancels").

## Names imported (one line each, so this document parses cold)

**𝔉** = `cross_form(a; b)`, the committed six-angle cotangent function (never plain F; F₀/F₁ are its w₃-split coefficients below). **V** = the variety {cos a₁+cos a₂+cos a₃ = 0} × {cos b₁+cos b₂+cos b₃ = 0}; **Ca**, **Cb** name the two constraints separately. **F127** = the theorem 𝔉 ≡ 0 on V. **z_k = e^{ia_k}, w_k = e^{ib_k}** are the monomial coordinates; in them Ca reads **Qz** = z₃² + S_z·z₃ + 1 = 0 with S_z = z₁+1/z₁+z₂+1/z₂, and Cb reads **Qw** = w₃² + S_w·w₃ + 1 = 0 with S_w = w₁+1/w₁+w₂+1/w₂. **_pieces convention**: for slot i with ascending complement {j < l}, psum = angle_j+angle_l, pdif = angle_l−angle_j, α_i = sin(angle_l)−sin(angle_j), β_i = −(sin(angle_l)+sin(angle_j)). An **atom** = one of the 72 elementary two-pole terms of §1, indexed (i, j, ξ, υ, e). A **sheet** = a hypersurface {L·(a,b) ≡ 0 mod 2π} with L ∈ {±1}⁶ canonical (first coefficient +1); an **event** = one atom-pole incidence on a sheet. **T** = the nine-term core function of §3. **K** = ℚ(i)(z₁,z₂,w₂)[z₃]/(Qz), the committed frame's coefficient field; **F₀, F₁** = the coefficients of 𝔉 = F₀(w₁) + F₁(w₁)·w₃ over K (the w₃-split; unique because K(w₁)[w₃]/(Qw) is free of rank 2 on {1, w₃}, Qw being monic of degree 2). A **norm** = a w₃-free denominator factor N(c) = c² + c·S_w + 1 produced by rationalizing a w₃-carrying cotangent factor (M−1), M = c·w₃^s with c a monomial and s = ±1. **The wall** = the grid+CRT certificate (527 items, 17 primes). **LEM** = the half-angle lemma cot(μ/2)[cot((μ+ρ)/2)+cot((μ−ρ)/2)] = 2(1+cos μ)/(cos ρ−cos μ); **tanLEM** its tangent twin. **Res** = the z₃-resultant of §3, not a residue; residues are always written Res_{w₁=r} or "residue". **The mirror specialisation** = b_j = π − a_{4−j}, i.e. (w₁,w₂,w₃) = (−1/z₃, −1/z₂, −1/z₁). **S1..S6** = the six exact gate statements of the exactify script (§8). **STEP_FIELD** = the committed gate in [`residue_assembly_close.py`](../../simulations/residue_assembly_close.py) proving Qz irreducible over ℚ(i)(z₁,z₂), so K is a field.

## Preface

The wall took ninety hours of grids. The note left at the door before it was finished said: if the folding pattern of this arc holds one more time, the ninety hours are not the proof, they are the pile, and inside there is something simple and fragile. The fifth visit sharpened the guess: the residues collapsing the moment both cosine sums vanish.

That is exactly what was inside. The 288 pole events of 𝔉 sort themselves onto 32 sheets, nine to a sheet, and on every sheet the nine residues sum to one and the same nine-term expression: a cotangent bilinear form in cofactor coefficients. That expression vanishes precisely when both cosine sums do, and the reason is one divisibility: the resultant of the two constraint quadratics divides the paired combination of the reduced numerator. The mirror quadratic of the old single-triple proof, doubled. Everything after that is bookkeeping a symmetry does for free: the function is odd under every single-angle flip, so its Laurent window collapses, and the one evaluation that pins the constant is the mirror case, which was proved over ℚ months before anyone knew it would become the anchor of the general theorem.

## Abstract

**Theorem (F127).** 𝔉 ≡ 0 on V.

**Proof chain (each algebraic step exact over ℚ(i), modulo one textbook analytic link, §4):** (§1) 𝔉 is a signed sum of 72 elementary atoms, each a difference of two simple poles (gate S1, per atom). (§2) The atoms' 288 pole events lie on 32 sheets, nine events each, and on every sheet the residue equals −(ΠL)/4 · T(L∘x) per event (gate S2), so one core function T covers all sheets through sign flips. (§3) The core identity: T = 0 on Ca ∧ Cb ∧ {Σa+Σb ≡ 0 mod 2π}, proved by a 540-term numerator and one exact divisibility, Res_{z₃}(Qz, P) | E; the degenerate sublocus where the divisibility is silent is closed separately. (§4) The transport lemma: in the committed w₃-eliminated frame the residue of F₀ and of F₁ at every norm root is a unit multiple of a sheet residue, hence zero; the three sparse norms outside the sheet lattice vanish by the b₁↔b₃ oddness (gate S6). (§5) With the committed structure results (all poles simple, Laurent window ⊆ {−1,0,1}, pinned again with w₂ distinguished) and the exact oddness bijections (gates S3, S5, S6), F₀ and F₁ collapse to g·(w₁−1/w₁)(w₂−1/w₂) with g free of both. (§6) The anchor: pulling back along the mirror morphism (gate S4 plus the committed ℚ-proof of the mirror case) and along its b₃-flip kills g on both Qw-branches: g = 0, so F₀ = F₁ = 0, so 𝔉 ≡ 0 on V. ∎

The mechanism in one line: **the residues collapse exactly when both cosine sums vanish, as divisibility by the resultant of the two constraint quadratics.**

## §1 The atoms

Write cot(μ/2)·Xh(μ; ξ, υ) for the committed building block of 𝔉 (Xh is the four-cotangent kernel; `cross_form` sums 36 of these with signs and coefficients). Two partial-fraction steps, both instances of LEM, turn each block into atoms:

  𝔉 = ¼ Σ_{i,j,ξ,υ,e} (−1)^{i+j} (−e) · c_ξ c_υ · (1+cos φ) · [ 1/(cos φ − cos(a_i+b_j)) − 1/(cos φ − cos(a_i−b_j)) ],

with φ = ξ + eυ; the index i runs over the a-slots, j over the b-slots, ξ over {psum_i, pdif_i} of the a-triple with c_ξ = α_i (psum) or β_i (pdif), υ the same on the b-side, e = ±1. That is 3·3·2·2·2 = 72 atoms, each carrying two simple poles, 144 poles, 288 pole events (each pole factor has two sheets, §2).

**Gate S1** proves the decomposition exactly, atom by atom (72 identities over ℚ(i) in monomial coordinates). The proof route matters: LEM folds the cotangent block into a common cos-pole form first; comparing raw cotangent and atom forms directly explodes (documented dead end, §8).

## §2 The sheet lattice

A pole of an atom is the condition cos φ = cos(a_i ± b_j), i.e. φ ∓ (a_i ± b_j) ≡ 0 mod 2π. Because φ = ξ + eυ mixes the two non-i a-angles and the two non-j b-angles with coefficients ±1, and a_i, b_j enter with coefficient ∓1, every pole condition is a linear form L·(a₁,a₂,a₃,b₁,b₂,b₃) ≡ 0 with **all six coefficients in {±1}**. Up to global sign that is at most 2⁶/2 = 32 sheets, and all 32 occur, with exactly **nine events each** (one per (i,j) block; machine-asserted in gate S2 and reproduced by an independent enumeration in review).

Near a sheet with canonical coordinate ℓ, cos φ − cos(a_i ± b_j) = −ε·sin(φ)·ℓ + O(ℓ²) for both sheet signs (ε = the canonicalization sign), so each event has residue −(atom prefactor)·c_s·ε/sin φ in ℓ, with c_s = ±1 the pole's sign in the bracket.

**Gate S2** proves, exactly and per event (288 identities): on the sheet, the event residue equals −(Π_k L_k)/4 times the corresponding (i,j) term of T(L∘x), where (L∘x)_k = L_k·x_k. Since cosine is even, Ca and Cb are invariant under every sign flip L, and the sheet condition says precisely that the flipped point satisfies Σa′+Σb′ ≡ 0. **One identity therefore covers all 32 sheets.**

## §3 The core identity

Define, in the _pieces convention,

  **T(a; b) = Σ_{i,j=1..3} (−1)^{i+j} α_i(a) α_j(b) · cot((a_i+b_j)/2).**

**Core identity.** T = 0 on Ca ∧ Cb ∧ {a₁+a₂+a₃+b₁+b₂+b₃ ≡ 0 mod 2π}.

The third constraint is load-bearing: on Ca ∧ Cb alone T is generically of order one (falsified numerically; the seed-pinned gate in `f127_core_identity.py` prints max |T| ≈ 35). The sheet relation is what turns the nine cotangents into the sheet's residue sum.

**Proof.** In monomial coordinates the sheet relation is z₁z₂z₃w₁w₂w₃ = 1, so w₃ = 1/(z₁z₂z₃w₁w₂) eliminates monomially, with no radicals and no coupling cost. Assemble T's numerator N in product form (nine terms, each a product of small factors; 540 terms expanded, degree 9 in z₃). Reduce N mod Qz: N ≡ r₀ + r₁·z₃. The second constraint transports to P = A·z₃² + B·z₃ + C with A = m₀², B = S_w·m₀, C = 1, m₀ = z₁z₂w₁w₂. On a common root of Qz and P, (B − A·S_z)·z₃ + (C − A) = 0, so define E := r₀(B − A·S_z) − r₁(C − A). Then:

1. **Res_{z₃}(Qz, P) divides E exactly** (E has 3128 terms, Res has 51; the division over ℚ(i) leaves remainder zero). Since {Res = 0} contains the projection of the whole variety, and E additionally vanishes on all 32 complex branches of {Res = 0} (checked numerically before the exact division), the divisibility proves N = 0 at every common root with B − A·S_z ≠ 0: no component caveat, no torus-only step.
2. **The degenerate sublocus.** E is silent where B = A·S_z; on the variety that forces A = C, i.e. m₀² = 1, and then S_z = m₀·S_w, where P ≡ Qz (the two quadratics coincide). Two exact facts close it: (i) gcd(Res, S_w − m₀·S_z) = 1, so no branch of {Res = 0} lies inside the degenerate locus, and Zariski closure extends N = 0 from the dense good locus to every branch; (ii) independently, on the sublocus itself every irreducible factor of S_z − εS_w (ε = m₀ = ±1) divides both numerators r₀ and r₁ (exact ideal membership), so r₀ = r₁ = 0 there outright.
3. **One reading caution.** On the ε = +1 sublocus the cotangent form of T is genuinely singular (z₃w₃ = 1 makes the (3,3) denominator vanish); "T = 0" is throughout the statement r₀ = r₁ = 0 in K[z₃]/(Qz), the pole-free witness, not a pointwise cotangent evaluation. ∎

This is the fragile thing. The mirror sub-case's proof mechanism (one quadratic divides one numerator) appears here doubled: two constraint quadratics, one resultant, one divisibility.

## §4 The transport lemma

Steps 1-3 kill the residues in the angle chart. The Laurent argument of §5 lives in the committed w₃-eliminated frame, so the vanishing must be transported there. The committed inventory has 43 denominator factors in w₁: 6 w₁-free (no w₁-pole), 6 linear (the outer poles a_i + b₁, whose grouped residues are the empty polynomial, committed), and 31 quadratic norms N(c) = c² + c·S_w + 1.

**Incidence (exact).** Each norm is norm(c) for a unique monomial c (machine-matched against the committed factor list). The 32 sheets collapse **2:1 onto 16 monomials** (M and its w₃-conjugate share c); 28 norms have c carrying w₁; the remaining **3 sparse norms norm(z_i)** come from the a_i + b₃ outer poles, which live outside the ±1 sheet lattice; the 12 root-sharing quad-quad pairs of the committed component plan are exactly the pairs {N(c), N(1/c)}, which differ by the unit c². There is no 32-vs-31 mismatch. The sheets collapse 2:1 onto 16 conjugation classes of monomials; the norm count is separately 31 = 28 dense + 3 sparse; the two bookkeepings meet factor by factor in the transport script's incidence table.

**Squarefreeness (exact, all 31).** Every norm has cleared w₁-degree exactly 2 and w₁-discriminant a nonzero element of K (reduced mod Qz and checked exactly). The tangency locus is {c² = 1} = {w₃ = ±1}, a proper sublocus; off it every pole of every norm is simple, sharpening the committed within-term coprimality.

**Transport (the linchpin).** Fix a norm N(c) with branch sign s and a root r. Two exact, loop-verified identities:

- (i) (M−1)(M̄−1) ≡ N mod Qw, for all 62 distinct M. At a simple root, differentiating along w₁ gives N′(r) = (M−1)′·(c²−1) with c²−1 a unit by squarefreeness: the change-of-uniformizer factor between the eliminated frame and the branch.
- (ii) For every incident term (384 incidences): A_t + B_t·c^s ≡ 0 mod (Qz, N), because the elimination multiplied the conjugate factor (M̄−1) into the numerator and M̄ = 1 on the branch.

Writing RF₀, RF₁ for the root-grouped residues of F₀, F₁ at r and G for the pulled-back sheet residue, (ii) gives RF₁ = −c^{−s}·RF₀ and (i) gives RF₀ + c^{−s}·RF₁ = G, hence

  **RF₀ · (1 − c^{−2s}) = G.**

The factor is a unit (c² ≠ 1 on the simple locus), and G = 0 by §2 + §3. So **RF₀ = RF₁ = 0 separately** on the simple locus; since RF₀, RF₁ are elements of K vanishing on a dense open subset, they vanish identically, at every root of every dense norm. The three sparse norms carry the a_i + b₃ outer poles; the b₁↔b₃ transposition is an exact oddness of 𝔉 (gate S6), so their grouped residues equal minus the a_i + b₁ outer residues, which are identically zero (committed). The single non-looped link in this section is the simple-pole product rule itself, standard calculus on an algebraic curve; every algebraic input to it is loop-proved. F₀ and F₁ therefore have **no finite w₁-pole over ℚ(i)**.

## §5 The window and the symmetries

The committed structure results now apply with their vanishing hypotheses supplied: F₀ and F₁, pole-free in w₁ with Laurent window ⊆ {−1,0,1} (committed, exact), are Laurent polynomials c₋₁/w₁ + c₀ + c₁·w₁ over K.

**Oddness (gates S3, S5, S6).** 𝔉 is odd under each single-angle flip and under the b-transpositions; each oddness is realized by an exact sign-reversing bijection of the 72 atoms (S3: b₁ → −b₁, υ-swap, 24 fixed atoms; b₁ ↔ b₂, j-swap with e-flip, 12 fixed; S5: b₃ → −b₃, υ-swap with e-flip, 24 fixed; S6: b₁ ↔ b₃). The flip b₁ → −b₁ is w₁ → 1/w₁; it preserves Qz and Qw (S_w is inversion-invariant), fixes w₃ (it does not touch b₃), and fixes every element of K. By uniqueness of the {1, w₃} coordinates, F₀ and F₁ are each odd: c₀ = 0 and c₋₁ = −c₁, so

  F_which = c·(w₁ − 1/w₁), c ∈ K.

**Relabel to w₂.** The same structure holds with w₂ distinguished; this is not argued by symmetry alone but literally rerun (the swap σ: w₁↔w₂ applied to the built term set, then the unmodified committed step functions: identical inventory 43 = 6+6+31, all poles simple, window ⊆ {−1,0,1}). With the b₁↔b₂ oddness, F_which = d·(w₂ − 1/w₂) likewise. Comparing coefficients in the common field:

  **F_which = g_which · (w₁ − 1/w₁)(w₂ − 1/w₂), g_which ∈ ℚ(i)(z₁,z₂)[z₃]/(Qz).**

## §6 The anchor

The mirror substitution φ*: (w₁, w₂, w₃) ↦ (−1/z₃, −1/z₂, −1/z₁) is a morphism into the variety defined over ℚ(i) (it satisfies Cb on Ca, since cos(π−a) = −cos a). Gate S4 proves exactly that φ*(𝔉) = ½·mirror_form(a) as rational functions, and mirror_form ≡ 0 on Ca is the committed ℚ-proof. So

  φ*(F₀) + φ*(F₁)·(−1/z₁) = 0 in the domain ℚ(i)(z₁,z₂)[z₃]/(Qz).

The factor φ*((w₁−1/w₁)(w₂−1/w₂)) = (z₃−1/z₃)(z₂−1/z₂) (each factor directly: w−1/w = z−1/z under w = −1/z) is a nonzero element of the domain: z₃² = 1 contradicts Qz's irreducibility (committed STEP_FIELD), and z₂ is a free variable. Dividing: g₀ − g₁/z₁ = 0. The second Qw-branch at the same (w₁, w₂) is w₃′ = −z₁ (the Qw-roots multiply to 1); it is reached from the mirror point by b₃ → −b₃, and gate S5's oddness makes 𝔉 vanish there too: g₀ − g₁·z₁ = 0. Subtracting, g₁(z₁ − 1/z₁) = 0, and z₁ is free, so **g₁ = 0, then g₀ = 0**. Hence F₀ = F₁ = 0, i.e. **𝔉 ≡ 0 on V**. ∎

The circle closes on itself: the mirror sub-case, proved over ℚ back when it was a lone specialisation gate, turns out to be the single evaluation the general theorem needs.

## §7 What the wall keeps

The wall (527 items, 17 primes, CRT with the realness guard) stays committed and untouched: it is an independent, fully deterministic certificate of the same statement, computed before this chain existed, and it gates this chain in return (the committed GF(p) residue certificate confirms §4's conclusion end to end at 764 root-point pairs). The honest trust surface of this proof: exact sympy/small-field arithmetic throughout the gates; one textbook analytic link (§4's product rule); the committed structure results it builds on, each itself exact. The wall's trust surface is disjoint (modular arithmetic, height bounds), which is exactly what makes the pair strong.

## §8 Verification index

All gates run from the repo root; outputs land next to the scripts. Runtimes on the dev machine.

| script | proves | runtime |
|---|---|---|
| `simulations/f127_core_identity.py` | §3: naked-T falsification, N assembly (540 terms), Res \| E exact, branch check | ~4 min |
| `simulations/f127_core_locus_patch.py` | §3.2: forcing chain, gcd(Res, S_w−m₀S_z) = 1, sublocus r₀ = r₁ = 0 | 19 s |
| `simulations/f127_exact_gates.py` | S1 (72 atoms), S2 (288 events), S3/S5/S6 (oddness bijections), S4 (mirror gate); LEM/tanLEM | 40 s |
| `simulations/f127_transport_lemma.py` | §4: incidence table, 31 discriminants, norm identity (62), step-2 divisibility (384) | 126 s |
| `simulations/f127_relabel_pin.py` | §5: the literal w₂-distinguished rerun (inventory, simple, window) | 8 s |
| `simulations/f127_sheet_lattice.py` | §2 numerics: the 32×9 lattice survey, residue bookkeeping gates | ~3 min |

Documented dead ends, kept so nobody retries them: sympy cancel/together on the unassembled six-variable sum (hangs; the product-form assembly and LEM folding are the cure); the naked core identity without the sheet relation (false); comparing cotangent and atom forms without the common cos-pole folding (explodes).
