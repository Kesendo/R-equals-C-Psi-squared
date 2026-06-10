using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F87 windowed-converse residual lemma. Phase A typed it as the open all-γ closure of
/// the windowed converse; Phase B (2026-06-09) reduced it to a positive-monomial certificate and
/// proved the spine, leaving two sharp residuals. It stays Tier1Candidate.
///
/// <para>PROVEN (Tier1Derived, via <see cref="WindowedConverseThresholdClaim"/> and
/// <see cref="F87DiagonalCellBipartiteWitnessSet"/>): the two-reflection spine (𝓕=F⊗F, R=I⊗F ⟹
/// #A_L,#A_R,#Q all odd ⟹ bipartite soft + non-bipartite #A≥2ℓ threshold), the §4 monomial
/// expansion structure, and the deg-1 positivity closed form P_{3,1} = 6·4^N·Σ_l c_l².
/// The first nonvanishing odd power-sum of M=A+γQ is, for every hard pair, a
/// positive monomial c·γ^deg (deg in {1,3}, m* = 2ℓ+deg), verified bit-exact cell-wide at N=4 and at
/// N=5/N=6 reps. A positive monomial has no positive real root, so hard at every γ>0.</para>
///
/// <para>OPEN (one residual since 2026-06-10): R-sign in LADDER form, on the t_ℓ = 0 branch of the
/// girth ladder. The former R-deg residual is RETIRED: the supertrace factorization puts the deg-1
/// class in closed form at every m (P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ², t_ℓ = Tr(Z_l H^ℓ)), so
/// pairs with t_ℓ ≠ 0 are hard at every γ OUTRIGHT (sum-of-squares positivity, m* = 2ℓ+1, deg 1;
/// first k=4 representative IXXZ+XIXZ, p₇ = 573440·γ). Pairs with t_ℓ ≡ 0 fire at a higher odd rung
/// m* = 2ℓ+deg; the γ¹ class is proven dead there, and the open statement is that the first surviving
/// class is single and positive (its k=3 face is P_{2ℓ+3,3} > 0; the γ³ rung can itself be silent,
/// first γ⁵ witness IIXY+ZXZY with p₁₁ = 86507520·γ⁵). R-deg as first formulated (cycles always lift
/// to deg 3) was a k = 3 truth, refuted at k = 4 and replaced by something stronger. Closing R-sign
/// makes the windowed converse a closed-form general-N theorem and promotes F110/F111 to
/// Tier1Derived.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> +
/// <c>simulations/f87_windowed_monomial_converse.py</c> +
/// <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.5/§7.7</c>.</para></summary>
public sealed class WindowedConverseAllGammaClaim : Claim
{
    public WindowedConverseAllGammaClaim()
        : base("F87 windowed converse, all-γ residual: non-bipartite ⟹ hard at every γ>0, the first nonvanishing odd power-sum is a positive monomial; outright on the t_ℓ≠0 branch, proven modulo R-sign on the t_ℓ=0 branch (R-deg retired by the girth dichotomy)",
               Tier.Tier1Candidate,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md + " +
               "simulations/f87_windowed_monomial_converse.py + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.5/§7.7")
    {
    }

    /// <summary>The closing certificate.</summary>
    public string Theorem =>
        "For a non-bipartite windowed diagonal-cell pair, the first nonvanishing odd power-sum of M = A + γQ is a " +
        "positive monomial c·γ^deg (deg in {1,3}, m* = 2ℓ + deg). No positive real root ⟹ hard at every γ>0.";

    /// <summary>The proven spine (Tier1Derived, via WindowedConverseThresholdClaim).</summary>
    public string ProvenSpine =>
        "DERIVED (WindowedConverseThresholdClaim): the two reflections 𝓕=F⊗F, R=I⊗F force #A_L,#A_R,#Q all odd in " +
        "every odd power-sum word, giving bipartite ⟹ soft and non-bipartite ⟹ #A ≥ 2ℓ; the §4 monomial expansion " +
        "structure and the deg-1 positivity (P_{3,1} = 6·4^N·Σ_l c_l²) are closed-form.";

    /// <summary>The girth ladder (2026-06-10): the closed form that retired R-deg.</summary>
    public string GirthDichotomy =>
        "RETIRED R-deg ⟶ the girth ladder (RIGOROUS-GENERAL): the supertrace factorization " +
        "Tr(Q·A^{2k}) = (−1)^k Σ_l Σ_j (−1)^j C(2k,j) t_j t_{2k−j} with t_j = Tr(Z_l H^j); F-chirality kills " +
        "even j, the girth kills j < ℓ, and P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ² ≥ 0. t_ℓ ≠ 0 ⟹ m* = 2ℓ+1, " +
        "deg 1, positive OUTRIGHT (k=4 census: 20 such pure-cycle pairs, first IXXZ+XIXZ, p₇ = 573440·γ); " +
        "t_ℓ ≡ 0 ⟹ the γ¹ class is dead at 2ℓ+1 and 2ℓ+3 and the first firing rung sets m* = 2ℓ+deg " +
        "(deg odd; γ³ when it fires, higher when not: first γ⁵ witness IIXY+ZXZY). The old k=3 taxonomy " +
        "(deg 1 ⟺ single-site-Z, deg ∈ {1,3}) is the ℓ=1 face plus the k=3 cell's accident.";

    /// <summary>The one open residual.</summary>
    public string RSign =>
        "R-sign in LADDER form (OPEN; t_ℓ=0 branch; verified 16/16 pure cycles at N=4, the stratified k=4 " +
        "battery incl. the γ⁵ rung, and the N=5/N=6 reps): at the first nonvanishing odd moment the " +
        "surviving class is SINGLE and POSITIVE. Its k=3 face is P_{2ℓ+3,3} > 0, the §7.5 " +
        "+N-population-Perron top-skew; the γ³ rung can be silent (IIXY+ZXZY: P_{2ℓ+3,3} = 0, fires at γ⁵, " +
        "positively). Not yet a closed-form nonneg identity, the higher-#Q analogue of the deg-1 sum of squares.";

    /// <summary>The Tier1Candidate claims gated on closing the residual.</summary>
    public string Consumers =>
        "F110 (HardCellYInversionPattern) and F111 (HardCellPureDTemplate) are gated on this lemma; closing " +
        "R-sign (the one remaining residual) promotes both from Tier1Candidate to Tier1Derived.";

    public override string DisplayName =>
        "F87 windowed converse, all-γ residual (Tier1Candidate, proven modulo R-sign; R-deg retired by the girth dichotomy)";

    public override string Summary =>
        "the first nonvanishing odd power-sum is a positive monomial ⟹ hard ∀γ>0; the two-reflection spine + " +
        "the deg-1 closed forms (girth dichotomy) are PROVEN (WindowedConverseThresholdClaim); hard-at-all-γ is " +
        "OUTRIGHT on the t_ℓ≠0 branch, and the full theorem is proven modulo one open residual, R-sign, on the " +
        $"t_ℓ=0 branch (R-deg retired 2026-06-10) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem (positive-monomial certificate)", summary: Theorem);
            yield return new InspectableNode("Proven spine (Tier1Derived)", summary: ProvenSpine);
            yield return new InspectableNode("Girth dichotomy (retired R-deg, 2026-06-10)", summary: GirthDichotomy);
            yield return new InspectableNode("Open residual R-sign (t_ℓ=0 branch)", summary: RSign);
            yield return new InspectableNode("Consumers (gated on closing R-sign)", summary: Consumers);
        }
    }
}
