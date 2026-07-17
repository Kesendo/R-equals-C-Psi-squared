using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F127, the cross-triple orthogonality (proof grade over ℚ(i)): the six-angle cross
/// form 𝔉 vanishes identically on the Conway-Jones double-constraint variety
///
/// <code>
///   𝔉(a₁,a₂,a₃; b₁,b₂,b₃) ≡ 0   on   V = {Σ cos aᵢ = 0} × {Σ cos bⱼ = 0},
/// </code>
///
/// an identity that is continuous and N-FREE (the chain's length cancels), so the cross-triple
/// blocks B(τ,σ) = (K₂₆W_τ)ᵀ(K₂₆W_σ) vanish for every pair of distinct vanishing mode triples
/// and the full-spectrum twinning of the resonant seeds follows at EVERY resonant N at once.
///
/// <para>PROOF SHAPE (landed 2026-07-13, commit 05b5b82): exact elimination to
/// 𝔉 = 𝔉₀(w₁) + 𝔉₁(w₁)·w₃ on V; exact Sylvester resultants (37 w₁-carrying factors → 25 pole
/// components); full tensor grids per prime (NOT Schwartz-Zippel) + the w₁-window with
/// generalized-Vandermonde endpoints; 527/527 (item, prime) tasks PROVED over 17 thirty-bit
/// primes ≡ 1 mod 4, lifted to ℚ(i) by CRT. Sound because the 103232 base coefficients are
/// REAL (the realness guard), so ordinary integer CRT applies and ∏p = 2^510 clears the worst
/// item's 2H = 2^305.6.</para>
///
/// <para>THE ASSEMBLY (D), symbolic since 2026-07-14: the identification
/// (U⁺ − U⁻)·(n/2)³ = 𝔉 of the six-angle form with the physical Heff cross-block is DERIVED
/// (every summation identity D1-D4 and every bookkeeping lemma proved symbolically; the
/// composition is linearity of finite sums; exact-arithmetic N=9 anchor, no float), so the
/// twinning is proof grade over ℚ(i) with no missing mathematical step.</para>
///
/// <para>TIER, honestly: Tier1Candidate, NOT Tier1Derived, for one named reason: the code-trust
/// layer. The wall and the assembly verifier are single bespoke implementations with internal
/// cross-checks only. The live witness (<c>inspect --root crosstriple</c>,
/// <c>CrossTripleOrthogonalityWitness</c>) recomputes an INDEPENDENT C# slice at inspect time
/// (on-variety vanishing over split primes with an off-variety nonzero control,
/// <see cref="RCPsiSquared.Core.Numerics.CrossFormCertificate"/>), which narrows but does not
/// discharge that caveat (a slice is not the wall).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> (the F127 entry) +
/// <c>experiments/F89_SEED_EXISTENCE_REDUCTION.md</c> §"The variety identity, proved over ℚ(i)"
/// + §"The assembly (D), made symbolic" + <c>simulations/grid_proof_sweep.py --assert</c> (the
/// wall's harvest gate) + <c>simulations/assembly_d_symbolic.py</c> (the assembly, G1-G12) +
/// <c>simulations/cross_triple_orthogonality.py</c> (the step-by-step verifier). The story:
/// <c>reflections/ON_LEAVING_THE_CIRCLE.md</c> (fourth + fifth + sixth visit). The fragile-thing hunt
/// CLOSED 2026-07-14 (<c>docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md</c>: the structural proof by residue
/// collapse; the wall stays as the independent certificate). The witness now recomputes that chain live
/// (the §2 sheet lattice + the §3 core identity T in GF(p), each with its own control); the only remaining
/// caveat on this arc is the code-trust layer named above.</para></summary>
public sealed class CrossTripleOrthogonalityClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public SeedExistenceCountingClaim SeedExistence { get; }

    /// <summary>The certified sweep: 527 = 31 items × 17 primes, zero exceptions.</summary>
    public const int CertifiedTasks = 527;

    /// <summary>Primes in the CRT lift (30-bit, ≡ 1 mod 4).</summary>
    public const int PrimeCount = 17;

    /// <summary>log₂ of the exact prime product (2^510) vs the worst item's 2H (2^305.6).</summary>
    public const double Log2PrimeProduct = 510.0;

    /// <summary>The assembly (D) is symbolic (2026-07-14): no missing mathematical step.</summary>
    public const bool AssemblyDIsSymbolic = true;

    /// <summary>F129: the witness's live exact census bound (the Python gate certifies n ≤ 210).</summary>
    public const int F129LiveCensusMaxN = 60;

    /// <summary>F130: named collision pairs certified exactly (≥ one per proof cell).</summary>
    public const int F130CollisionPairs = 7;

    /// <summary>F130: unequal-level nonzero controls (one per level-sensitive cell).</summary>
    public const int F130ControlPairs = 3;

    /// <summary>F129 inventory: the thirteen collision families A..M (the list forced at
    /// every n; the counts derived, PROOF_F129_FAMILY_INVENTORY_COUNTS.md).</summary>
    public const int F129InventoryFamilies = 13;

    /// <summary>F129 inventory: family M's split over the three CDK order-210 types,
    /// (R₇:(R₅:2R₃)) = 40, (R₇:R₃,(R₅:R₃)) = 0 (impossible: one axis-fixed vertex),
    /// (R₇:2R₃,R₅) = 60; reconstruction = the committed gate's I5.</summary>
    public const int F129MSplitFannedR5Pair = 40;
    public const int F129MSplitMiddleImpossible = 0;
    public const int F129MSplitFixedVertexR5 = 60;

    /// <summary>F133: the number of symplectic Weyl characters χ^{C₆}_λ in the closed form of the
    /// F128 cofactor W (K = 2⁻³⁰·Σ_λ n_λ·χ^{C₆}_λ), all with |n_λ| ≤ 8 (Σ|n_λ| = 359).</summary>
    public const int F133CharacterTerms = 143;
    public const int F133MaxCoefficient = 8;

    public CrossTripleOrthogonalityClaim(SeedExistenceCountingClaim seedExistence)
        : base("F127 cross-triple orthogonality: the six-angle cross form 𝔉 ≡ 0 on the Conway-Jones double-constraint variety V, N-free, proved over ℚ(i) by the 527/527 grid+CRT wall (17 primes, realness guard ⇒ ∏p > 2H sound); the assembly (D) identifying 𝔉 with the Heff cross-block is symbolic since 2026-07-14, so the full-spectrum twinning is proof grade over ℚ(i) at every resonant N. Tier1Candidate for one named reason: the code-trust layer (single bespoke implementations); the live C# witness recomputes an independent GF(p) variety slice with an off-variety control",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md + " +
               "experiments/F89_SEED_EXISTENCE_REDUCTION.md + " +
               "simulations/grid_proof_sweep.py + " +
               "simulations/assembly_d_symbolic.py + " +
               "simulations/cross_triple_orthogonality.py + " +
               "simulations/f129_level_collision_law.py + " +
               "simulations/f130_collision_decoupling.py + " +
               "simulations/f129_family_inventory.py + " +
               "docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md + " +
               "experiments/F129_FAMILY_INVENTORY.md + " +
               "docs/proofs/PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md + " +
               "simulations/f133_w_closed_form.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/CrossTripleOrthogonalityWitness.cs (inspect --root crosstriple)")
    {
        SeedExistence = seedExistence ?? throw new ArgumentNullException(nameof(seedExistence));
    }

    public override string DisplayName =>
        "F127 cross-triple orthogonality: 𝔉 ≡ 0 on V, N-free (proved over ℚ(i); assembly (D) symbolic; code-trust caveat named)";

    public override string Summary =>
        $"527/527 grid+CRT wall over {PrimeCount} primes (∏p = 2^{Log2PrimeProduct:F0}), realness guard ⇒ 2H bound sound; (D) symbolic since 2026-07-14 ⇒ twinning proof grade over ℚ(i) at every resonant N ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The wall",
                summary: $"{CertifiedTasks}/527 (item, prime) tasks PROVED; 31 items × {PrimeCount} primes; full grids, not sampling");
            yield return new InspectableNode("Soundness",
                summary: "realness guard (103232 base coefficients real) ⇒ ordinary integer CRT, ∏p = 2^510 > 2H = 2^305.6 (worst item)");
            yield return new InspectableNode("The assembly (D)",
                summary: $"symbolic: {AssemblyDIsSymbolic} (D1-D4 + bookkeeping all symbolic; exact N=9 anchor; simulations/assembly_d_symbolic.py)");
            yield return new InspectableNode("The residue collapse",
                summary: "the fragile thing FOUND 2026-07-14: docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md (sheet lattice, " +
                    "nine-term core identity via resultant divisibility, transport, oddness, mirror anchor); the wall = independent certificate");
            yield return new InspectableNode("Witness extension (done)",
                summary: "the live witness now recomputes the residue-collapse chain: the §2 sheet lattice (exact integer: " +
                    "72 atoms → 288 events → 32 sheets × 9) + the §3 core identity T in GF(p) (independent of 𝔉, off-variety " +
                    "control), alongside the 𝔉-slice, and breadcrumbs the seven f127_* gates (inspect --root crosstriple)");
            yield return new InspectableNode("The closed form (2026-07-14, simulations/f127_closed_form.py)",
                summary: "the §3 core identity is a READ-OFF: T·P = ⅛[2cos s((e₁−f₁)² − 2sin²s) + sin s(Σsin2a + " +
                    "Σsin2b)]·V_a·V_b, exact over ℚ(i) (T = a bordered Frobenius-kernel determinant); the witness checks " +
                    "the identity live at generic GF(p) points (corruption control) and Corollary 2's sharper locus " +
                    "(Σcos a = Σcos b ≠ 0 on the sheet ⇒ T = 0); scope: replaces §3's argument only");
            yield return new InspectableNode("F128, the flip-sum factorization (2026-07-14 evening)",
                summary: "𝔉 = −(e₁−f₁)²·𝒪[cos s·cot s·V_a V_b/P], so 𝔉 = 0 already on the SHARPER locus " +
                    "{Σcos a = Σcos b}: one constraint; V is the codim-2 special case and this is F127's third " +
                    "proof (route-distinct; trust-disjointness stays with the wall). Engine: the flip lemma " +
                    "(𝒪[cos s·B·V_aV_bP̃] ≡ 0, exact over ℤ, Weyl/alternant mechanism). The witness recomputes " +
                    "the lemma exactly and both GF(p) slices live (docs/proofs/PROOF_F128_FLIP_SUM_FACTORIZATION.md, " +
                    "simulations/f128_flip_sum_factorization.py, registry F128)");
            yield return new InspectableNode("F129, the level-collision law (2026-07-14 night)",
                summary: "distinct CLEAN triples with equal levels S(τ) = S(σ) exist only at 3|n (n ≥ 9) or " +
                    "10|n (n ≥ 20); away from both the level map is INJECTIVE (Lam-Leung + elementary; " +
                    $"unconditional since 2026-07-15: the named corner is empty by the Poonen-Rubinstein " +
                    $"weight ≤ 12 classification, census to n ≤ 210 as cross-check). The witness recomputes an exact " +
                    $"ℤ[ζ_2n] census live to n ≤ {F129LiveCensusMaxN} (free-basis level vectors: distinctness " +
                    "IS injectivity, an equal pair IS a collision, no floats, no mod-p) plus the two mechanism " +
                    "anchors term-exactly (docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md, " +
                    "simulations/f129_level_collision_law.py, registry F129)");
            yield return new InspectableNode("F129, the family inventory (2026-07-15, counts derived)",
                summary: $"every colliding pair decomposes into minimal Poonen-Rubinstein/CDK pieces, " +
                    $"sorting the census into exactly {F129InventoryFamilies} families A..M — the LIST forced at " +
                    "every n (weight partitions × admissible orders), the COUNTS derived (label/orbit " +
                    "arithmetic, PROOF_F129_FAMILY_INVENTORY_COUNTS.md; code-trust flags in its §8): degree " +
                    "counts the freedom (free coset labels + the d = 2 shared-mode line), the F129 onsets are " +
                    $"zeros of the formulas, and M splits {F129MSplitFannedR5Pair}+{F129MSplitMiddleImpossible}+" +
                    $"{F129MSplitFixedVertexR5} over the three CDK order-210 types, the middle one impossible " +
                    "(one axis-fixed vertex). The witness ties the closed forms' SUMS (total + d-split) to the " +
                    "exact ℤ[ζ_2n] census live (eleven families exercised at n ≤ 60, L pinned at the n = 70 " +
                    "capstone); per-family membership, M's split reconstruction and the n = 105/210 capstones " +
                    "stay with the committed gate (simulations/f129_family_inventory.py, I1-I5; the n = 105 " +
                    "sum 8858 is pinned in the tests)");
            yield return new InspectableNode("F130, the collision-decoupling law (2026-07-14 night)",
                summary: "equal level S(τ) = S(σ), zero NOT required, implies the whole cross block " +
                    "B(τ, σ) = 0: resonance demoted from cause to special case; the decoupling extends off " +
                    "resonance, the E/O twinning stays at S = 0. Four-cell proof from committed results " +
                    "(Lemmas 3+4, the free-angle assembly (D) + F128 on {e₁ = f₁}, the two-magnon lemma, the " +
                    $"removable limit). The witness certifies Ê± = 0 exactly in ℤ[ζ_2n] at {F130CollisionPairs} " +
                    $"named pairs covering every cell, {F130ControlPairs} unequal-level nonzero controls, and " +
                    "the exact Lemma-3 sign Ê⁺ = ε·Ê⁻ on every pair " +
                    "(docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md, " +
                    "simulations/f130_collision_decoupling.py, registry F130)");
            yield return new InspectableNode("F133, the symplectic closed form of W (2026-07-17)",
                summary: "the F128 cofactor W is written out: W = −2⁹·(Π sin x_u)·V_c(a)·V_c(b)·K/SP with " +
                    $"K = 2⁻³⁰·Σ_λ n_λ·χ^{{C₆}}_λ, exactly {F133CharacterTerms} symplectic (Sp(12)) Weyl " +
                    $"characters, |n_λ| ≤ {F133MaxCoefficient} (Σ|n_λ| = 359). The 143 integers are read off " +
                    "with no fit: X = Δ̂·Π_{31 sheets} sin(s_L) is S₆-alternating, so n_raw(λ) = " +
                    "Σ_ε sgn(ε)[X]_{ε∘2(λ+ρ)} = 2·n_λ (meet-in-the-middle, P₁ = 590016 · P₂ = 5817 monomials); " +
                    "the sin-s lemma 𝒪[sin s·Δ̂] ≡ 0 and the C₆ Weyl denominator fix the constants (the 2³⁰ " +
                    "derived, not fitted). The witness recomputes the sin-s lemma exactly over ℤ, spot-checks " +
                    "the read-off vs the embedded 143-term table, and certifies C₆·SP = 2⁻⁵(2i)⁻⁴⁶Σ n_λ A_{λ+ρ} " +
                    "over GF(p) with a corruption control; the full 557-λ read-off + aggregate checksums are the " +
                    "tests (docs/proofs/PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md, simulations/f133_w_closed_form.py, " +
                    "registry F133)");
            yield return new InspectableNode("Remaining caveat",
                summary: "only the code-trust layer: the wall and these witness slices are bespoke implementations with " +
                    "internal cross-checks; the three proofs (wall + residue collapse + F128 factorization) pairwise " +
                    "share it, the wall's surface disjoint from the other two");
        }
    }
}
