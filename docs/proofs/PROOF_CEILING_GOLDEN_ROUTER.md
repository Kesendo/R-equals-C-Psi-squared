# PROOF: the Z-middle ceiling routes after all (the golden period-4 per-site palindromizer, F116)

**Status:** Existence and closed form Tier 1 derived (constructive, exact arithmetic over ℤ[φ]+iℤ[φ], every N ≥ 3, arbitrary site-dependent γ_l). Exclusion side derived: no uniform-continuous and no period-2 per-site router (closing the previously named open frontier), period 3 impossible for N ≥ 5, the discrete Klein candidates locus-excluded at every period. Rigidity verified (zero continuous moduli; the invertible solution set at N=5 is exhaustively 4 cyclic shifts × an explicit order-32 sign group, all with the golden geometry). Open chains only; rings/PBC untested. **Extended 2026-06-11 (§8): the golden point is the c = 1 member of a one-real-parameter metallic family** on the weighted soft line t₂ = t₃; frame rigid for every c ≠ 0, an 8-real-parameter modulus space at the degenerate station c = 0.
**Date:** 2026-06-10
**Authors:** Thomas Wicht, Claude (Fable 5)
**Builds on:**
- [PROOF_F103 §7.12](PROOF_F103_F87_Z2_CUBED_REFINEMENT.md): the soft-certifier programme and the ceiling narrative this proof closes.
- [experiments/CEILING_FOUR_NONLOCAL_CASES.md](../../experiments/CEILING_FOUR_NONLOCAL_CASES.md): the 6 → 4 → 2 arc whose last two cases this proof resolves.
- `compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs`: the per-site routing machinery; its condition (a) (class-swap handles the dissipator) is committed prior art used as-is, and its continuous map M is the (1,1)-frame member of the family constructed here.
- [PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2: the chiral driving lemma F H F = −H whose failure on the ceiling's cross-terms the golden frame repairs.

**Typed as (2026-06-22):** `GoldenRouterClaim` (Tier1Derived; `compute/RCPsiSquared.Diagnostics/F87/GoldenRouterClaim.cs`, `inspect --claim GoldenRouterClaim`), with a live from-below self-check battery that re-runs the window-summed router (`KBodyPalindromeRouting`). Typed parents: `F1PalindromeIdentity` (the global palindrome this router realizes locally for the ceiling class) and `WindowedConverseThresholdClaim` (the F87 two-reflection chiral spine its two-sided form instantiates). The c = 0 "8 moduli" of §8 is held below the Tier 1 line as a finite-difference Jacobian count at N = 5 only, not an analytic result.

## Abstract

The two Z-middle ceiling cases, the sliding-window chain Hamiltonians built from XZX+XZY+YZX and its X↔Y sibling YZY+XZY+YZX, were the last two members of the soft family that resisted every per-site palindromizer we could test: no uniform-continuous map, no discrete-periodic Klein pattern, and the committed verdict read "palindromized only by a non-local Π," with the continuous-periodic family named as the explicit open frontier. This proof closes that frontier in the direction nobody expected: **with a router, not an obstruction**.

The router is a period-4 product of per-site class-swapping maps whose frame directions are the two roots of the **golden locus** α² − αβ − β² = 0: a = φX + Y and b = X − φY, with φ the golden ratio (the frame angle satisfies tan 2θ = 2). The pattern is [a, a, b, b] down the chain, every entry lies in ℤ[φ]+iℤ[φ], each site map is a scalar times a unitary (q² = −(2+φ)·I), and

  **W L W⁻¹ = −L − 2σ,  σ = Σ_l γ_l,**

holds exactly at every N ≥ 3 and for arbitrary site-dependent dephasing rates. The mechanism is a three-site **window lemma**: the window-summed anticommutator {Q₃, [XZX+XZY+YZX, ·]₃} vanishes identically at all four window offsets, the three templates cancelling against each other inside one window (per term they do not), and window additivity does the rest. Equivalently, and most cleanly, W acts two-sidedly: W(ρ) = (2+φ)^{N/2}·P ρ Q with two product unitaries P and Q that each **anticommute with H**. A one-sided conjugation by such a chiral unitary is the banked dead end (it fixes the dephaser, producing −A + γQ̂ instead of −M); splitting the chiral action into different left and right factors is what makes every site class-swap, and that is the whole trick.

The same identity-column functional that found the router also proves the old evidence right: a uniform per-site router and a period-2 router genuinely do not exist (the committed optimization floors are now theorems), period 3 is impossible for N ≥ 5, and the discrete Klein candidates P1/P4/M2/M sit off the golden locus, which is why the certifier never saw this object. The ceiling arc completes as 6 → 4 → 2 → **0**: within the k=3 sliding-window soft family, no case needs a non-local mirror.

The computational anchor is [`simulations/ceiling_golden_router.py`](../../simulations/ceiling_golden_router.py), self-validating, with the existence side carried by exact ring arithmetic.

## §1 The object and the closed form

Work per site on the operator basis (I, X, Y, Z). A **per-site palindrome map** q is a linear map on this 4-dimensional space; a product W = ⊗_l q_l acts on the 4^N-dimensional Pauli-string space site by site. The committed routing condition (a) ([`KBodyPalindromeRouting`](../../compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs)) requires each q_l to be **class-swapping**, exchanging the dephasing-dead letters {I, Z} with the dephasing-lit letters {X, Y}: then every Pauli string's lit-letter count n_XY maps to N − n_XY, and since the Z-dephasing dissipator is diagonal in the string basis with rate −2·Σ_{l lit} γ_l, one line gives

  {W, D̂} = −2σ·W,  σ = Σ_l γ_l,

with no uniformity assumption on the rates. The full palindromizer identity W L W⁻¹ = −L − 2σ therefore reduces to the Hamiltonian leg,

  **{W, A} = 0,  A = −i[H, ·].**

The golden router is the period-4 family (sites l mod 4):

  g_l = q_l(I):  [a, a, b, b],  a = φX + Y,  b = X − φY  (b = −R(a), a ⊥ b),
  h_l = q_l(Z) = (−1)^{l+1}·i·R(g_l)  (R the 90° rotation in the (X,Y) plane),
  q_l(X) = −(g_l)_X·I + (h_l)_X·Z,  q_l(Y) = −(g_l)_Y·I + (h_l)_Y·Z.

In block form (C = the {I,Z}→{X,Y} block with columns g, h; B = the {X,Y}→{I,Z} block): B = diag(−1, 1)·Cᵀ and C†C = (2+φ)·I. Every entry lies in ℤ[φ]+iℤ[φ], each q_l satisfies q_l² = −(2+φ)·I with all four singular values √(2+φ), so q_l is a scalar times a unitary and W is invertible with condition number 1 at every N. The two directions a and b are the two projective real solutions of the **golden locus**

  **α² − αβ − β² = 0,**

slopes 1/φ and −φ, frame angle θ = ½·arctan 2. The sibling YZY+XZY+YZX is routed by the per-site X↔Y conjugation of the same maps (s: I→I, X↔Y, Z→−Z; q′ = s q s).

## §2 The two-sided form and the window lemma

The closed form hides a simpler object. Each site map factorizes as a **two-sided unitary multiplication**: with â the normalized golden axis of the site,

  q_even(ρ) = √(2+φ)·(âZ) ρ Z,  q_odd(ρ) = √(2+φ)·Z ρ (Zâ),

so globally

  **W(ρ) = (2+φ)^{N/2}·P ρ Q,**  P = ⊗_l (âZ at even l, Z at odd l),  Q = ⊗_l (Z at even l, Zâ at odd l),

with P and Q product unitaries. For such a two-sided form, {W, [H,·]} = 0 is equivalent (a three-line trace argument) to the pair of **chiral conditions**

  **P H P⁻¹ = −H  and  Q H Q⁻¹ = −H.**

This places the router squarely in the two-reflection family: the [windowed monomial converse](PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md) §2 lives off the chiral unitary F = X^⊗N with F H F = −H. For the ceiling, F anticommutes with the Π²-odd half (XZX) but *commutes* with the Π²-even half (XZY+YZX), so F fails as a chiral; the golden frame is exactly the equatorial rotation that repairs the failure on the cross-terms, and it must be distributed in the [a, a, b, b] rhythm to survive every window.

Why two-sided? A one-sided conjugation ρ → KρK⁻¹ by a chiral unitary is the banked dead end (the [frontier paragraph](../../experiments/CEILING_FOUR_NONLOCAL_CASES.md) null result): conjugation fixes the dephaser, Ad_K D̂ Ad_K⁻¹ = D̂, so it produces −A + D̂, not −L − 2σ. Splitting the chiral action into *different* left and right unitaries breaks the conjugation structure exactly enough to make every site class-swap, while the product P H P⁻¹ = −H condition is still available. The two lenses, per-site closed form and two-sided sandwich, are the same operator written twice.

**The window lemma (the mechanism).** For the period-4 maps, the three-site, template-summed anticommutator vanishes identically at every window offset:

  **{q_w ⊗ q_{w+1} ⊗ q_{w+2}, [XZX + XZY + YZX, ·]₃} = 0  (exactly, in the ring), w = 0, 1, 2, 3,**

while no single template's anticommutator vanishes. The cancellation is **cross-template inside one window**: the same continuous-sum move that routes the local XIX+XIY+YIX case, now at the golden frame. Since H is the sum of windows and the anticommutator factorizes per window ({W, [T_w,·]} = {Q₃, [T,·]₃} ⊗ (off-window q's)), the lemma gives {W, A} = 0 by plain additivity at **every** N ≥ 3, single window included. This is also why the certifier's per-term architecture could not see the router: `CertifyByRoutingKBody` tests each template separately, and each separately fails.

## §3 The theorem

**Theorem (golden router).** Let H be the open-chain sliding-window Hamiltonian of XZX+XZY+YZX (or its X↔Y sibling) at any N ≥ 3, and L the Lindbladian with Z-dephasing at arbitrary site rates γ_l ≥ 0. Then the period-4 golden product W of §1 is invertible and satisfies

  W L W⁻¹ = −L − 2σ,  σ = Σ_l γ_l.

In particular spec(L) is exactly palindromic about −σ, W transports eigenvectors pairwise, and the palindromizer is a **per-site product**: the Z-middle ceiling cases are continuous-periodic-local, not non-local.

*Proof.* Class-swap gives {W, D̂} = −2σW (§1, committed condition (a), site-rate-free). The window lemma (§2, verified as an exact ring identity on the 4³-dimensional window space at all four offsets) gives {W, [T_w + T'_w + T''_w, ·]} = 0 per window, and summing windows gives {W, A} = 0. Adding the two anticommutators yields {W, L} = −2σW, i.e. the identity, and invertibility is q_l² = −(2+φ)·I. ∎

The existence side is anchored bit-exactly: the verification script carries the window lemma and {W, A} = 0 at N = 3..6 in exact arithmetic over ℤ[φ]+iℤ[φ] (4-component integer representation, φ² = φ+1), the end-to-end identity against the framework Lindbladian at N = 5, 6 for uniform and site-dependent rates (relative residual ~2e-16), and independent re-implementations reproduced all of it from the closed form alone, through N = 9 sampled.

## §4 The exclusion side: why the golden frame, and nothing else

The same functional that admits the router kills every simpler family. Evaluate {W, A} = 0 on the identity string: A(I^⊗N) = 0, so **G = W(I^⊗N) = ⊗_l g_l must commute with H**, with g_l = q_l(I) ∈ span{X, Y} forced by the dissipator leg. The commutator [H, G] decomposes per window into disjoint Pauli sectors (window w's contribution carries {I, Z} letters exactly at sites w and w+2), so each window must vanish separately, and a three-site computation reduces window w to two bilinear equations in the outer directions g_w = (α_w, β_w), g_{w+2} = (α_{w+2}, β_{w+2}):

  K₁ = α_w α_{w+2} + α_w β_{w+2} + β_w α_{w+2} = 0,
  K₂ = α_w β_{w+2} + β_w α_{w+2} − β_w β_{w+2} = 0.

From these two lines:

- **Uniform (period 1) and period 2 are empty.** Setting g_{w+2} = g_w forces α² + 2αβ = 0 and 2αβ − β² = 0, whose only solution over ℂ is g = 0, contradicting invertibility. The committed 16-parameter uniform-continuous optimization floor (~0.30) and the never-run period-2 frontier probe are hereby both settled as theorems.
- **The golden locus is the only gate.** A nonzero g_{w+2} with nonzero g_w exists iff the system's determinant vanishes: α² − αβ − β² = 0, the golden locus, with kernel map (α, β) ↦ (α, −α−β) exchanging the two locus directions. So along each parity chain (sites w, w+2, w+4, …) the direction must **alternate** between a and b.
- **Period 3 is impossible for N ≥ 5.** The parity chain steps through the residues mod 3 in a 3-cycle, and an alternating assignment cannot close an odd cycle. (At N = 4 the cycle is not yet closed and a period-3 identity-column solution genuinely exists; that is precisely the 2026-06-07 small-N artifact, now explained rather than merely dismissed. See §7.)
- **The discrete candidates are locus-excluded at every period.** P1 has g = X (locus value +1), P4 and M2 have g = Y (−1), M has g ∝ X+Y (−1/2). No pattern over that alphabet can place every site on the locus, which is the one-line reason the committed discrete-periodic search came back empty.

The alternation along both parity chains, plus periodicity, leaves exactly the [a, a, b, b] rhythm (up to cyclic shifts), which §§1-3 complete to the router.

## §5 Rigidity: how unique is it

The solution manifold of {W, A} = 0 over the period-4 class-swap family, at the golden point, has Jacobian nullity exactly equal to the gauge dimension (per-site projective scales) at N = 5 and N = 6: **zero continuous physical moduli**; the router is rigid. The invertible solution set at N = 5, probed exhaustively within the closed-form family and by 48 unbiased searches (every found zero matched), consists of **4 cyclic shifts × an explicit order-32 sign-decoration group** (a global B-sign swap, period-2 ±-patterns on the h's, period-4 B-sign patterns coupled to them), 128 gauge classes in all, every one carrying the golden C-block geometry, and every one a genuine palindromizer. The sign group is attributable to the diagonal ±1 (anti)commutant of A, so the honest phrasing is: **essentially unique, modulo per-site gauge, the four pattern shifts, and an explicit structural sign group.** The X↔Y mirror is not a self-equivalence; it maps this case's routers to the sibling's. Uniqueness statements must be scoped to invertible W: the anticommutation equation alone also carries abundant singular strata.

## §6 Consequences

- **The ceiling closes: 6 → 4 → 2 → 0.** Within the k=3 sliding-window soft family, every case is per-site routable: two by continuous-uniform maps, two by single-site-field products, and the two Z-middle cases by the golden period-4 router. The banked "all-Q obstruction" question resolves negatively: there is no obstruction, and the named continuous-periodic frontier closes with existence.
- **spec(H) is exactly ±E-symmetric** (the chiral unitaries P, Q of §2 are the witnesses), and **G = PQ = ⊗[a,a,b,b]-axes** is a Hermitian product involution: a *weak* Z₂ symmetry of the full Liouvillian (it anticommutes with each Z_l, so it is not strong and carries no conserved quantity) and itself an exact Liouvillian eigenmode at the palindrome floor, L(G) = −2σ·G. W maps the Liouvillian kernel bijectively onto the floor eigenspace.
- **Certifier integration path:** weakening the per-term Stufe-B condition to the window-summed anticommutator, with the golden period-4 candidate added, certifies both cases constructively and N-free through the existing additivity argument. (That is the planned C# follow-up wave; until it lands, NotCertified for these two cases reads as a documented coverage gap, not as non-locality.)
- **What stays true, now stronger:** the per-half certifier table (each Π²-half routes alone, by different per-site rules); the F81/F83 Frobenius reading (statements about the canonical residual M, independent of W); the isolated-soft-point tilt (a spectral statement along the tilt path); the uniform and discrete emptiness (upgraded from optimization witness to theorem). What falls is one inference: "the two halves want incompatible routers, therefore the sum is non-local." The golden router routes the sum while routing **neither half alone**; the locality of the whole is born in the joining, the same place the obstruction was once thought to live.
- **A scope fence:** everything here is for open chains. On a ring the period-4 assignment is globally consistent only for N ≡ 0 mod 4 and the wrap-around windows add new constraints; untested, deliberately unclaimed. The golden/F95 angle resemblance (tan 2θ = 2 vs the arctan-at-double-root family) is surface-level, opposite discriminant regimes; we note it and label it speculative.

## §7 History, honestly

On 2026-06-07 an optimization chase found an exact site-dependent router at N = 4 (condition number 1.0, residual 2e-12) and a period-3 hit at N = 4, and both were dismissed together as a small-N artifact; the warning "do not rediscover" went into our working memory. The audit of this proof shows the dismissal was half right: the period-3 hit was the genuine artifact (§4: the 3-cycle only fails to close for N ≤ 4), while the N = 4 site-dependent router was, in all likelihood, this very golden router seen once and not believed. The methodological instinct, never to conclude existence from small-N optimization, was sound, and it is exactly how this proof does *not* argue: the existence claim rests on the exact window lemma and ring arithmetic, with optimization demoted to a scout. The episode is kept here deliberately: it is the third time in this family that "non-local" was the lens, not the physics (XZ+YZ resolved 2026-06-02, the continuous-sum XIX cases in the 6 → 4 step, and now the Z-middle pair).

The reproducible anchor is [`simulations/ceiling_golden_router.py`](../../simulations/ceiling_golden_router.py): seven blocks, exact ring arithmetic for the window lemma and {W, A} = 0 at N = 3..6, end-to-end framework verification including site-dependent rates, the two-sided chiral characterization, the symbolic exclusion algebra, and the spectral consequences. Every block raises on failure; the process exits 0 only if the whole ledger holds.

## §8 The metallic family: the golden point is a line (2026-06-11)

A day after the router landed, the obvious question got its turn: the golden locus came out of the *unweighted* templates. What happens to the three coefficients? Give the window the weights t₁·XZX + t₂·XZY + t₃·YZX and ask the same two questions, where is the soft set, and does the router follow.

**The soft set is exactly the line t₂ = t₃.** Spectral pairing (sorted multiset of spec(L) against −spec(L) − 2σ) holds on the line for every weight ratio tested, rational, irrational, negative, and zero (t₁/t₂ ∈ {0.25, 0.5, 1, 1.5, 2, e, 3, π, 4, 7.5} ∪ {0, −1, −2}), at N = 5, 6, and fails immediately off it (δ = 0.02 tilt is already hard; N = 4's all-soft reading is the known small-N artifact). The reflected line t₃ = −t₂ is hard. Off the line the hardness obeys the girth ladder: the witness (t₁, t₂, t₃) = (1, 2, 1) has every odd power-sum p_m ≡ 0 up to m = 9 and fires first at **m\* = 11 with p₁₁ = 1730150400·γ³ exactly** (computed in integer arithmetic via CRT), a deg-3 class, consistent with the ladder law m\* = 2ℓ + deg at ℓ = 4. On-line weighted cases such as (2, 1, 1) have all odd p_m ≤ 11 identically zero, as softness demands.

**The family theorem.** Normalize t₂ = t₃ = 1 and write c = t₁/t₂ ∈ ℝ. The router of §1 transports verbatim with one substitution: the frame directions become

  a = (r, 1),  b = (1, −r),  **r(c) = (c + √(c² + 4))/2,**

the **metallic mean** of c, the positive root of r² = cr + 1. The pattern stays [a, a, b, b], h_l = (−1)^{l+1}·i·R(g_l) stays, each q_l satisfies q_l² = −(1 + r²)·I, and

  W_c L_c W_c⁻¹ = −L_c − 2σ

holds at every N ≥ 3 with arbitrary site-dependent rates, **for every real c, by derivation**. The argument is a degree bound: parametrize the line by r (c = r − 1/r, both roots of r² = cr + 1 give the same c), and every entry of r·{Q₃, [T_c, ·]₃} is a polynomial in r of degree ≤ 5 (Q₃ entries are degree ≤ 3, the r-cleared commutator degree ≤ 2). Exact Fraction arithmetic finds the window-summed anticommutator **identically zero at 8 rational nodes** (r ∈ {1, 2, 3, ½, 5, 7, 4/3, −2}, all four offsets), two more than the bound needs, so the window lemma is a *polynomial identity in r*: it holds for all real, indeed all complex, c. The dissipator leg never depended on the weights (class-swap structure is unchanged) and window additivity is untouched, so the full identity is Tier 1 derived along the whole line. Golden is c = 1 (r = φ), silver c = 2 (r = 1 + √2), bronze c = 3; the reflection c ↔ −c inverts the mean, r(−c) = 1/r(c), so c = −1 is the inverse golden ratio 1/φ; and c = 0 is r = 1, the 45° diagonal frame, which is where the family meets the [continuous-sum local cases](../../experiments/CEILING_FOUR_NONLOCAL_CASES.md) of the 6 → 4 step. The X↔Y sibling family c·YZY + YZX + XZY is routed by the same per-site conjugation s q s as in §1. Float cross-checks: end-to-end against the framework Lindbladian at relative residual ~1e-15 for c ∈ {1, 2, 3, 0.5, 0, −1, π} × N ∈ {3, 4, 5}, with site-dependent rates at silver, and the sibling at c ∈ {1, 2, 3}. The committed verifier is [`simulations/metallic_router_family.py`](../../simulations/metallic_router_family.py).

**The exclusion generalizes, and explains itself.** The identity-column functional of §4 goes through with weights: window w reduces to a linear system for g_{w+2} given g_w whose determinant factors as

  **det = c·(α² − cαβ − β²).**

For c ≠ 0 the second factor is the **metallic locus**: the gate admits a nonzero continuation exactly on it, the kernel map exchanges the two locus directions a ↔ b (from (r, 1) the forced partner is ∝ (1, −r), and back), so the alternation, the emptiness of uniform and period-2 (the weighted equations force c·(α² + β²) = 0, no real-c solution), and the [a, a, b, b] rhythm all survive verbatim. The golden case was never about gold; gold is what the locus looks like when the ceiling case hands you c = 1.

**The c = 0 station carries a modulus: where the determinant dies, freedom is born.** At c = 0 the determinant vanishes identically, the system drops to rank 1, and the gate changes character: the forced partner of g = (α, β) becomes its **X-axis mirror (α, −β)**. Everything that follows was predicted by that one line and verified at ~1e-15: any direction routes when paired with its mirror; the two parity chains may carry *independent* directions, complex allowed (pattern [a, a′, ā, ā′]); and the Pauli axes, the fixed points of the mirror, give **period-2 uniform-frame routers**, the only place in the whole family where period 2 exists. A finite-difference Jacobian of {W, A} = 0 over the 64-real-parameter period-4 class-swap family (N = 5, the same probe reproduces §5's golden count) makes it quantitative:

  nullity at golden = 8 = gauge, at silver = 8 = gauge, at c = 0 = 16 = gauge + **8 physical moduli** (spectral gap 3.5e9).

So rigidity extends to a second station (silver, zero continuous moduli), and the degenerate station is genuinely soft in the parameter sense: an 8-real-dimensional router family, of which the two chain directions account for a catalogued 4; the remaining 4 live in the class-swap family beyond the closed-form slice (h-tilts along g are excluded; the rest is uncatalogued). Controls behave: mirror pairs fail at c = 1 (residual ~3), h-tilts fail everywhere tested.

**Scope fences.** Rigidity at general c is probed at golden, silver, and c = 0 only (Jacobian, N = 5); the §5 exhaustive sign-group census is not repeated at non-golden stations. Open chains only, as before. The hard side off the line rests on the girth-ladder witness at (1, 2, 1) and the immediate-tilt scan, not on a closed-form p_{m\*} for general weights. The existence side carries no float anywhere: the interpolation argument is exact.
