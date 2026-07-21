# PROOF F139: the seam identity, the F134 wall as a Chebyshev divisor

*2026-07-21. Consumes [PROOF_F133_W_SYMPLECTIC_CLOSED_FORM](PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md) (the letter system and the read-off) and [PROOF_F134_TWO_ROW_REFLECTION_LAW](PROOF_F134_TWO_ROW_REFLECTION_LAW.md) (the law this identity derives, and the column polynomials P_k it recovers as quotients). Grade: a derivation chain of six hand-provable lemmas ending in a finite exact-ℤ polynomial division, machine-checked with a trust-disjoint C# recomputation; the F127-wall epistemic class, now for a mechanism rather than a bare law.*

Before the machinery: what this feels like. F134 ended with one open thing, the identity performing the reflection. A season of hunting had proved what the mechanism is not; this document banks what it is. The whole two-row table, read through the right chain of collapses, is a single polynomial in a cosine variable, and the affine wall at eleven is not a symmetry of anything: it is a divisor. Dividing by the vanishing polynomial of the eleventh-roots cosine lattice returns the committed column polynomials as quotients (up to the sign (−1)^k) and a small remainder, and the reflection law is nothing but the statement that the remainder stays below the window. The breaks at a third row of two, unexplained since the day they were measured, are the remainder overflowing. The mystery mapped in F134 §5 closes here.

## §1 Definitions

All objects derive from F133's letter system. There, X = Δ̂·∏₃₁ is a Laurent polynomial in six half-angle units t₁, t′ = (t₂, …, t₆), and the committed coefficients n_λ are read off by the signed sum n_raw(λ) = Σ_{ε∈{±1}⁶} sgn(ε)·[X]_{ε∘2(λ+ρ₆)} = 2 n_λ, ρ₆ = (6,5,4,3,2,1). A two-row λ is (λ₁, k, 0, 0, 0, 0), written n_(λ₁,k); μ₁ = λ₁ + 6.

**Units.** Work in x = t′² (five variables). Write ρ₅ = (5,4,3,2,1), δ = (2,1,0,−1,−2) = ρ₅ − 3·𝟙, and a_δ(x) = Σ_{σ∈S₅} sgn(σ)·x^{σδ} (the A₄ alternant at δ; equal to the half-angle Vandermonde Δ̂₅ = ∏_{u<v}(x_u^{1/2}x_v^{−1/2} − x_v^{1/2}x_u^{−1/2}), integral because each variable meets four factors).

**The letter system.** Λ₃₆ := {e₁, …, e₅} ∪ ({±1}⁵ ∖ {−𝟙}), the x-exponent vectors of F133's 36 t₁-carrying factor coefficients (F133 [R3] gives X = Δ̂·∏₃₁; the 36 = 5 + 31 factorization in t₁ is re-derived here in Lemma 1). Split it as

- **P** := {e₁, …, e₅, 𝟙}, the 6-letter *window* (the five vector letters and the top cube vertex), and
- **Core** := {±1}⁵ ∖ {±𝟙}, the 30 remaining cube letters: negation-closed, total sum 0, falling into 15 antipodal pairs {M, −M} with *cosines* C_M := x^M + x^{−M}.

**The read-off functional.** For a Laurent polynomial G in x and a one-row target ν = (k, 0, 0, 0, 0) with μ = ν + ρ₅,

  R⁻[G] := Σ_{ε∈{±1}⁵} Σ_{σ∈S₅} sgn(ε)·sgn(σ)·[G]_{ε∘μ + 𝟙 − σδ}.

(This is the coefficient of x^{ε∘μ}, summed with signs, of a_δ·x^{−𝟙}·G.) Two notational warnings: the bare 5-vector μ = ν + ρ₅ here is not the six-variable μ = λ + ρ₆ above, whose first coordinate μ₁ = λ₁ + 6 survives as the slice label; and ψ, φ below depend on k through R⁻'s target ν, a subscript we suppress (Φ_k carries it for both).

**The seam objects.**

  ψ_{i,t} := R⁻[e_i(P)·e_t({C_M})]  (i = 0..6, t = 0..15; elementary symmetric functions of the window letters and of the 15 cosines),
  φ_i(y) := Σ_t ψ_{i,t}·y^{15−t},
  **Φ_k(y) := (y² − 1)·φ₀(y) + y·φ₁(y) + φ₂(y)**  (a polynomial of degree ≤ 15, Lemma 5).

**Chebyshev conventions.** S_m(y) is the Chebyshev polynomial of the second kind in the 2cos normalization: S₀ = 1, S₁ = y, S_m = y·S_{m−1} − S_{m−2}, so S_m(2cosθ) = sin((m+1)θ)/sinθ and z^{m+1} − z^{−(m+1)} = (z − z⁻¹)·S_m(z + z⁻¹). In particular

  **S₁₀(y) = y¹⁰ − 9y⁸ + 28y⁶ − 35y⁴ + 15y² − 1 = ∏_{r=1}^{10} (y − 2cos(rπ/11))**,

the vanishing polynomial of the Niven cosine lattice at 11 (over ℚ it splits into the two quintic minimal polynomials of 2cos(π/11) and 2cos(2π/11)), and (z − z⁻¹)·S₁₀(z + z⁻¹) = z¹¹ − z⁻¹¹: **the F134 wall factor** (t²² − t⁻²² in F133's units, z = t²).

**The committed column polynomials** (F134 §2c, written in the S-basis via C_a = t^a + t^{−a}, z = t²): P₀ = S₁ − S₅, P₁ = S₄ − S₂ − 2S₀, P₂ = 3S₁, P₃ = −(S₂ + 2S₀), P₄ = S₁, P₅ = 0.

## §2 The statement

**Theorem (F139, the seam identity).** For k = 0..5, polynomial division of Φ_k by S₁₀ gives

  **Φ_k = S₁₀ · (−1)^k P_k + R_k,  with deg R_k ≤ 4 + k, and R₀ = 0**

(so Φ₀ = S₁₀·(S₁ − S₅) exactly). Moreover the a-priori chain reproduces the committed table: writing Φ_k = Σ_m b_m S_m,

  **n_(λ₁,k) = (−1)^{λ₁} · b_{λ₁+5}**  on every parity-live strip cell (gate G4, 27/27 including zeros).

Here n_(λ₁,k) is F134's two-row coefficient n_(j,k) with j = λ₁, the *strip* is the dominant two-row index set {0 ≤ k ≤ min(λ₁, 5), λ₁ ≤ 10}, and *parity-live* means λ₁ + k even (odd |λ| vanishes identically in F133); the strip has 27 such cells, and the gate checks all of them, zeros included.

**Corollary (F134).** n_(λ₁,k) = n_(10−λ₁,k) on the two-row slice. *Proof from the theorem:* S₁₀·S_j = S_{10+j} + S_{10+j−2} + … + S_{10−j} (Chebyshev product rule), so any multiple S₁₀·Q with deg_S Q ≤ 5 has S-coefficients symmetric about m = 10 throughout m ∈ [5, 15]; deg Φ_k ≤ 15 forces deg_S Q_k ≤ 5. The remainder feeds only m ≤ 4 + k, strictly below the k-window's low edge m = 5 + k (μ₁ = 6 + k). Hence b_{10+u} = b_{10−u} for every m = 10 ± u in the window [5+k, 15−k], which is the law read through n = (−1)^{λ₁}b_{λ₁+5}. ∎

The wall 11 is therefore not a hidden symmetry of the coefficient array (F134's centroid theorem forbids that reading); it is a **divisor**: the window palindromy is inherited from S₁₀, the beyond-wall shadow of F134 §2c is the remainder, and Θ_k ↔ S₁₀·(−1)^k P_k under (z − z⁻¹)·S₁₀ = z¹¹ − z⁻¹¹. The F134 Θ-decomposition, banked there as "equivalent to the law, not a proof", is now the division algorithm.

## §3 The derivation chain

Each step is a small lemma; together they build Φ_k from F133's letters with no input from the coefficient table. The committed gate ([simulations/f139_seam_identity.py](../../simulations/f139_seam_identity.py)) checks each machine-checkable step; the table enters only as the final gate.

**Lemma 1 (the slice).** [t₁^{2μ₁}] X = (−1)^d · Δ̂₅ · t′^{−2𝟙} · e_d(Λ₃₆-doubled), d = 18 − μ₁, where the letters appear doubled in t′-units (2e_v and 2M), i.e. exactly Λ₃₆ in x = t′² units with the shift x^{−𝟙}. *Proof:* X = Δ̂₅(t′)·∏_{i=1}^{36}(a_i t₁ − a_i^{−1} t₁^{−1}) (F133's factorization, collected in t₁; the 36 letters a_i are the five difference factors and the 31 sheets). Choosing the +t₁ branch from p factors gives t₁-degree 2p − 36 = 2μ₁ and coefficient (−1)^{36−p}·(∏ a_i^{−1})·e_p(a²); with ∏ a_i^{−1} = t′^{2𝟙} and e_p(a²) = t′^{−4𝟙}·e_{36−p}(a^{−2}), 36 − p = d. Gate G1 samples this against the committed F133 read-off at μ₁ ∈ {6, 9, 11, 13, 16}. ∎

**Lemma 2 (one slice suffices).** n_(λ₁,ν) equals the 5-variable read-off of the positive slice alone: the ε₁ = −1 half of F133's read-off is the x-inversion image of the ε₁ = +1 half, and inversion carries sign (−1)⁵ through the 5-variable signed sum, so the two halves contribute equally; the factor 2 cancels against n_raw = 2n. Hence, with c_d := [z^d] R⁻[∏_{ℓ∈Λ₃₆}(1 + z·x^ℓ)],

  n_(λ₁,ν) = (−1)^d · c_d,  d = 12 − λ₁. ∎

**Lemma 3 (the pair factorization).** Per Core pair, (1 + z x^M)(1 + z x^{−M}) = 1 + z·C_M + z² = z·(y + C_M) with y := z + z⁻¹. Hence E_Core(z) = z^{15}·∏_{15 pairs}(y + C_M): the Core's whole z-dependence lives in y. Expanding ∏(y + C_M) = Σ_t e_t({C_M})·y^{15−t} gives the cosine matrix: c-contributions ψ_{i,t} with the pure kinematic factor y^{15−t}. ∎

**Lemma 4 (the window skew).** ∏_{p∈P} x^p = x^{2𝟙}, so e_{6−i}(P)(x) = x^{2𝟙}·e_i(P)(x⁻¹); and R⁻ is inversion-odd (a_δ is inversion-invariant since −δ is the even reversal of δ; the coefficient extraction at ε∘μ summed over ε with sgn(ε) picks up (−1)⁵ under x ↦ x⁻¹; the x^{±𝟙} twists exchange). Hence for any inversion-invariant H (in particular H = ∏(y + C_M), Core being negation-closed):

  R⁻[e_{6−i}(P)·H] = −R⁻[e_i(P)·H],  so ψ_{6−i,t} = −ψ_{i,t} and ψ_{3,t} = 0.

Gate G2 checks the full grid. ∎

**Lemma 5 (the Chebyshev fold).** Assembling Lemmas 2 to 4: c(z) = z^{15} Σ_{i=0}^{6} z^i φ_i(y) = z^{15}·Σ_{i=0}^{2} (z^i − z^{6−i})·φ_i(y) = −z^{18}·[(z³ − z⁻³)φ₀ + (z² − z⁻²)φ₁ + (z − z⁻¹)φ₂] = **−z^{18}·(z − z⁻¹)·Φ(y)** by z^j − z^{−j} = (z − z⁻¹)S_{j−1}(y) and S₂ = y² − 1, S₁ = y, S₀ = 1. Extracting [z^d] via z^{m+1} − z^{−(m+1)} = (z − z⁻¹)S_m(y) gives, for the S-coefficients b_m of Φ, c_d = −f_{d−18} with f_e = b_{e−1} (e ≥ 1), f_e = −b_{−e−1} (e ≤ −1); at e = −μ₁ this is n_(λ₁,k) = (−1)^{λ₁}·b_{λ₁+5}. The degree bound deg Φ_k ≤ 15 is gate G3 (the top ψ entries cancel in the (y²−1)-assembly). ∎

**Lemma 6 (the division).** Finite exact-ℤ polynomial division of each Φ_k by S₁₀, gate G5: the quotients equal (−1)^k P_k coefficient for coefficient, the remainders obey deg R_k ≤ 4 + k, and R₀ = 0. ∎

The chain's shape, in the arc's older vocabulary: Lemma 3 is where the affine translation by 22 becomes visible (it acts as the y-grading of the Core), Lemma 4 is the second exact mirror, and Lemma 6's remainder-degree bound is the one fact that is verified rather than derived (the division itself is mechanical; Lemma 5's degree bound is likewise gate-checked but has a derivable cancellation reason); §5 records that residue.

## §4 The domain fence

The same chain runs at any ν (replace μ = ν + ρ₅). Measured through the division (gate G6):

- **l = 1** (ν = (k₁, 1), where F134's reflection holds on the three-row slice): remainder degree ≤ 4, below every window; the l = 1 hold is remainder smallness.
- **l = 2** (ν = (2,2), (3,2), (4,2), the slices carrying F134's 8 breaks): remainder degree 8 to 9, **overflowing the window**: the breaks are exactly the remainder reaching into the palindromy range. The break atlas of F134 §4 is the support of the remainder tower.

The two-row law's own out-of-window continuations behave the same way: for k = 2 and k = 4 the b-sequence fails palindromy only at non-partition edge positions (λ₁ < k), which no table cell reads.

## §5 What remains open

The identity reduces F134's mechanism question to one line: **why is deg R_k ≤ 4 + k?** The division verifies it; a structural derivation (a reason the shadow is small, without performing the division) is the remaining beauty. Two doors, both consistent with everything measured: the values Φ_k(2cos(rπ/11)) as root-of-unity specializations z = e^{irπ/11} of the letter product (the factorization y_r + C_M = ζ^{−r}(ζ^r + x^M)(ζ^r + x^{−M}) turns Φ(y_r) into an evaluation of ∏_{ℓ}(1 + ζ^r x^ℓ), a Weyl-denominator-type object), and the θ₁/Appell-Lerch reading of the remainder tower (F134 §3's (−,+) characteristic; the truncation-with-shadow shape is the partial-theta signature). Neither is needed for the theorem.

## §6 Verification

- Python gate: [simulations/f139_seam_identity.py](../../simulations/f139_seam_identity.py) (G1 slice identity vs the committed F133 read-off, G2 skew, G3 degree lemma, G4 the a-priori table 27/27, G5 the division with Q_k = (−1)^k P_k and the remainder bounds, G6 the fence, G7 corruption control; ~2-3 min, exact ℤ throughout; imports only the committed F133 gate).
- C# second implementation: `WSymplecticClosedForm.AnalyzeSeamIdentity` rebuilds ψ, Φ, and the division on an independent code path (own packing, own Chebyshev algebra, embedded expected quotients) and cross-checks the a-priori table against the embedded F134 two-row table; surfaced by the `CrossTripleOrthogonality` witness family (see the F139 registry entry for the test count).
- The coefficient table: [simulations/results/f133_w_closed_form/chiC_coeffs.txt](../../simulations/results/f133_w_closed_form/chiC_coeffs.txt) (committed with F133; used here only as the gate).
