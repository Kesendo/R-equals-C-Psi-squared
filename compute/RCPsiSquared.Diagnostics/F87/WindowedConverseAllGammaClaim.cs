using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F87 windowed-converse all-γ theorem, CLOSED. Phase A typed it as the open all-γ
/// closure of the windowed converse; Phase B (2026-06-09) reduced it to a positive-monomial
/// certificate with two residuals (R-deg + R-sign); 2026-06-10 retired R-deg in the morning wave
/// (the girth dichotomy) and resolved R-sign the same day (the Pascal-Gram positivity theorem).
/// Tier1Derived, no residual.
///
/// <para>THE THEOREM: for a non-bipartite windowed diagonal-cell pair, every γ-coefficient of the
/// first nonvanishing odd power-sum of M = A + γQ is non-negative; each surviving #Q = d class is
/// the equal-leg-total Pascal-Gram sum of squares P_{m*,d} = (m*/d)·Σ_{l⃗}Σ_{k⃗}|U^{(l⃗)}_{k⃗}|²
/// (U the binomial transform of the d-leg moments T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁}···Z_{l_d}H^{α_d})),
/// and every other class vanishes exactly (parity mod 4, girth, or the cascade). At least one class
/// is positive, so p_{m*}(γ) > 0 for every γ > 0: hard at one γ is hard at ALL γ.</para>
///
/// <para>The proof chain: cyclic decomposition → leg factorization (A_L/A_R binomial + supertrace
/// split) → Hermitian conjugacy (ket leg = conj bra leg, the F113-Lemma-C transpose-trick sibling)
/// → leg parity (F-chirality: odd totals only) + leg girth (totals ≥ ℓ) → Vandermonde assembly
/// (C(α+β,β) = Σ_k C(α,k)C(β,k), prefactor (−i)^u(+i)^u = +1) → slice inversion (U at |k⃗| = u IS
/// T) → cascade induction (p_m ≡ 0 below m* kills all lower-total moments). Selection rule
/// corollary: classes fire only for d ≡ m*−2 (mod 4), d ≤ m*−2ℓ, so monomiality is DERIVED for
/// deg ≤ 3 (every k=3 case); from deg = 5 two classes may coexist and positivity carries alone
/// (the γ⁵ witness IIXY+ZXZY is single because t₅ = 0 too, but nothing needs that).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md</c> §5 +
/// <c>simulations/f87_pascal_gram_positivity.py</c> (all five proof steps verified exactly:
/// d = 1/3/5 Pascal-Gram == exact CRT coefficient on all five branch representatives) +
/// <c>simulations/f87_girth_dichotomy.py</c> + <c>simulations/f87_windowed_monomial_converse.py</c>.
/// </para></summary>
public sealed class WindowedConverseAllGammaClaim : Claim
{
    public WindowedConverseAllGammaClaim()
        : base("F87 windowed converse, all-γ theorem: non-bipartite ⟹ hard at every γ>0; every coefficient of the first nonvanishing odd power-sum is a Pascal-Gram sum of squares or exactly zero (R-deg retired + R-sign resolved 2026-06-10, no residual)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F87_WINDOWED_MONOMIAL_CONVERSE.md + " +
               "simulations/f87_pascal_gram_positivity.py + " +
               "simulations/f87_windowed_monomial_converse.py + " +
               "docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md §7.5/§7.7")
    {
    }

    /// <summary>The closing certificate.</summary>
    public string Theorem =>
        "For a non-bipartite windowed diagonal-cell pair, every γ-coefficient of the first nonvanishing odd " +
        "power-sum of M = A + γQ is ≥ 0 (each surviving class an equal-total Pascal-Gram sum of squares), at " +
        "least one is > 0. Non-negative coefficients have no positive real root ⟹ hard at every γ>0.";

    /// <summary>The proven spine (Tier1Derived, via WindowedConverseThresholdClaim).</summary>
    public string ProvenSpine =>
        "DERIVED (WindowedConverseThresholdClaim): the two reflections 𝓕=F⊗F, R=I⊗F force #A_L,#A_R,#Q all odd in " +
        "every odd power-sum word, giving bipartite ⟹ soft and non-bipartite ⟹ #A ≥ 2ℓ; the §4 monomial expansion " +
        "structure and the deg-1 positivity (P_{3,1} = 6·4^N·Σ_l c_l²) are closed-form.";

    /// <summary>The girth ladder (2026-06-10 morning): the closed form that retired R-deg.</summary>
    public string GirthDichotomy =>
        "RETIRED R-deg ⟶ the girth ladder (RIGOROUS-GENERAL): the supertrace factorization " +
        "Tr(Q·A^{2k}) = (−1)^k Σ_l Σ_j (−1)^j C(2k,j) t_j t_{2k−j} with t_j = Tr(Z_l H^j); F-chirality kills " +
        "even j, the girth kills j < ℓ, and P_{2ℓ+1,1} = (2ℓ+1)·C(2ℓ,ℓ)·Σ_l t_ℓ² ≥ 0. t_ℓ ≠ 0 ⟹ m* = 2ℓ+1, " +
        "deg 1, positive OUTRIGHT (k=4 census: 20 such pure-cycle pairs, first IXXZ+XIXZ, p₇ = 573440·γ); " +
        "t_ℓ ≡ 0 ⟹ the γ¹ class is dead at 2ℓ+1 and 2ℓ+3 and the first firing rung sets m* = 2ℓ+deg " +
        "(deg odd; γ³ when it fires, higher when not: first γ⁵ witness IIXY+ZXZY). The old k=3 taxonomy " +
        "(deg 1 ⟺ single-site-Z, deg ∈ {1,3}) is the ℓ=1 face plus the k=3 cell's accident.";

    /// <summary>The Pascal-Gram positivity theorem (2026-06-10, same day): R-sign resolved.</summary>
    public string PascalGramPositivity =>
        "RESOLVED R-sign ⟶ Pascal-Gram positivity (RIGOROUS-GENERAL): every coefficient P_{m,d} factorizes " +
        "through d-leg moments T^{(l⃗)}_{α⃗} = Tr(Z_{l₁}H^{α₁}···Z_{l_d}H^{α_d}); at m* the cascade (p_m ≡ 0 " +
        "below m* + slice inversion of the unitriangular binomial transform) kills every unequal-total block, " +
        "and the equal-total block is Σ|U|² with prefactor (−i)^u(+i)^u = +1. Verified exactly at d=1 " +
        "(IXXZ+XIXZ 573440), d=3 (K3 2064384, flux 589824, multi-Z 61440), d=5 (IIXY+ZXZY 86507520), plus the " +
        "cascade's forced zeros and the selection rule (d ≡ m*−2 mod 4, d ≤ m*−2ℓ: singleness DERIVED for " +
        "deg ≤ 3). The old k=3 face P_{2ℓ+3,3} > 0 (the §7.5 Perron top-skew) is the d=3 instance.";

    /// <summary>Claims that were gated on this lemma's residual.</summary>
    public string Consumers =>
        "F110 (HardCellYInversionPattern) and F111 (HardCellPureDTemplate) were gated on R-sign; with the " +
        "residual resolved the gate is open and both promote from Tier1Candidate to Tier1Derived.";

    public override string DisplayName =>
        "F87 windowed converse, all-γ theorem (Tier1Derived, no residual; R-deg retired + R-sign resolved 2026-06-10)";

    public override string Summary =>
        "every coefficient of the first nonvanishing odd power-sum is a Pascal-Gram sum of squares or exactly " +
        "zero ⟹ p_{m*}(γ) > 0 ∀γ>0 ⟹ hard ∀γ>0; the two-reflection spine + the girth dichotomy + the " +
        "Pascal-Gram positivity theorem close the windowed converse with NO residual (R-deg retired and R-sign " +
        $"resolved 2026-06-10) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem (Pascal-Gram positivity certificate)", summary: Theorem);
            yield return new InspectableNode("Proven spine (Tier1Derived)", summary: ProvenSpine);
            yield return new InspectableNode("Girth dichotomy (retired R-deg, 2026-06-10)", summary: GirthDichotomy);
            yield return new InspectableNode("Pascal-Gram positivity (resolved R-sign, 2026-06-10)", summary: PascalGramPositivity);
            yield return new InspectableNode("Consumers (gate now open: F110/F111 promote)", summary: Consumers);
        }
    }
}
