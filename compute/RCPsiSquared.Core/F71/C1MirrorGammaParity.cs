using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>F101 (Tier 1 derived, 2026-05-21): the F71 c₁ bond-mirror deviation
///
/// <code>
///     D(b)  :=  c₁(b) − c₁(N−2−b)
/// </code>
///
/// is an EXACTLY ODD function of the F71-anti-palindromic component of the per-site
/// Z-dephasing profile γ. F71 (<see cref="C1MirrorIdentity"/>) was proven for R-symmetric
/// γ; this is its graceful-breakdown extension to non-uniform γ, and closes the F71
/// "non-uniform γ_i" open question. Observable-side twin of F91 (the γ spectral-invariance
/// claim); the J-side counterpart is F100 (<see cref="C1QPeakMirrorJParity"/>).
///
/// <para>Decompose γ (a per-site profile of length N) at the site-mirror l ↔ N−1−l,
/// γ = γ_sym + γ_anti with γ_sym = (γ + F71(γ))/2 and γ_anti = (γ − F71(γ))/2. Then:</para>
///
/// <list type="bullet">
///   <item>D ≡ 0 for every palindromic γ (γ_anti = 0), however non-uniform γ_sym is. F71
///         never required uniform γ; it requires palindromic (R-symmetric) γ.</item>
///   <item>D(γ_sym, −γ_anti) = −D(γ_sym, γ_anti): exactly odd in γ_anti, all orders. The
///         Taylor series of D in γ_anti has odd powers only, so the breakdown is graceful
///         (leading-order linear in the per-site asymmetry), never a hard jump.</item>
///   <item>The leading coefficient κ_γ generically DEPENDS on γ_sym (Tier 2 empirical):
///         the parity argument fixes the oddness, not the coefficient. κ_γ admits no
///         closed form: it is ∂c₁/∂γ and c₁ is itself a non-closed numerical fit. See the
///         κ-obstruction section of <c>docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md</c>.</item>
/// </list>
///
/// <para>Mechanism (cf. <c>docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md</c>): the F71
/// chain-mirror R is R-equivariant on the whole c₁ pipeline — the dephasing dissipator
/// D[ρ] = Σ_l γ_l (Z_l ρ Z_l − ρ) is linear in each γ_l and R relabels site l ↔ N−1−l —
/// so c₁(b; γ) = c₁(N−2−b; F71(γ)); hence D(b; F71(γ)) = −D(b; γ). Substituting
/// F71(γ) = γ_sym − γ_anti yields the oddness. This is the observable-side twin of F91:
/// F91 says the F71-refined diagonal-block spectrum depends only on γ_sym; F101 says the
/// c₁ bond-mirror deviation lives in γ_anti.</para>
///
/// <para>Scope: covers c₁ (closure-breaking coefficient, vac+SE probe states). The F86c
/// Q_peak observable is NOT covered: its Q-axis Q = J/γ₀ is defined against a scalar γ₀,
/// which a non-uniform γ profile does not provide; the γ_avg-anchored Q_peak route is a
/// separable extension noted in the proof. Numerically witnessed with residuals ≤ 5e-9 at
/// N=3,4,5 (<c>simulations/f71_nonuniform_gamma_verification.py</c>).</para>
/// </summary>
public sealed class C1MirrorGammaParity : Claim
{
    public C1MirrorGammaParity()
        : base("F101 c₁ mirror-deviation parity: D(b) = c₁(b) − c₁(N−2−b) is odd in γ_anti, zero for palindromic γ",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F101 + docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md")
    { }

    /// <summary>F71-palindromic component γ_sym = (γ + F71(γ))/2, where F71 maps site
    /// l ↔ N−1−l. Palindromic by construction (γ_sym[l] = γ_sym[N−1−l]). The c₁
    /// bond-mirror deviation D is blind to this component.</summary>
    public static double[] PalindromicComponent(IReadOnlyList<double> siteGamma)
    {
        if (siteGamma is null) throw new ArgumentNullException(nameof(siteGamma));
        int n = siteGamma.Count;
        var gSym = new double[n];
        for (int l = 0; l < n; l++)
            gSym[l] = 0.5 * (siteGamma[l] + siteGamma[n - 1 - l]);
        return gSym;
    }

    /// <summary>F71-anti-palindromic component γ_anti = (γ − F71(γ))/2. Anti-palindromic
    /// by construction (γ_anti[l] = −γ_anti[N−1−l]); the deviation D is exactly odd in it,
    /// and zero iff γ_anti is zero.</summary>
    public static double[] AntiPalindromicComponent(IReadOnlyList<double> siteGamma)
    {
        if (siteGamma is null) throw new ArgumentNullException(nameof(siteGamma));
        int n = siteGamma.Count;
        var gAnti = new double[n];
        for (int l = 0; l < n; l++)
            gAnti[l] = 0.5 * (siteGamma[l] - siteGamma[n - 1 - l]);
        return gAnti;
    }

    /// <summary>Norm ‖γ_anti‖ = sqrt(Σ_l γ_anti[l]²) of the anti-palindromic component:
    /// the deviation of γ from being F71-palindromic. Zero iff γ_l = γ_{N−1−l} for all l.
    /// The c₁ bond-mirror deviation D vanishes iff this is zero.</summary>
    public static double PalindromicDeviation(IReadOnlyList<double> siteGamma)
    {
        double sum = 0.0;
        foreach (double gAnti in AntiPalindromicComponent(siteGamma))
            sum += gAnti * gAnti;
        return Math.Sqrt(sum);
    }

    /// <summary>True iff γ is F71-palindromic (γ_l = γ_{N−1−l} ∀l) within tolerance: the
    /// regime in which the c₁ bond-mirror identity D(b) = 0 holds exactly.</summary>
    public static bool IsPalindromic(IReadOnlyList<double> siteGamma, double tolerance = 1e-10)
        => PalindromicDeviation(siteGamma) < tolerance;

    public override string DisplayName =>
        "F101: D(b) = c₁(b) − c₁(N−2−b) is odd in γ_anti (graceful breakdown of F71)";

    public override string Summary =>
        "F71 c₁ bond-mirror deviation is exactly odd in the anti-palindromic γ-component; "
        + $"zero for palindromic γ, leading-order linear; numerically verified N=3,4,5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parity mechanism",
                summary: "R-equivariance gives c₁(b;γ) = c₁(N−2−b;F71(γ)); hence D(b;F71(γ)) = −D(b;γ), exactly odd in γ_anti (all orders)");
            yield return new InspectableNode("palindromic survival",
                summary: "D ≡ 0 for every palindromic γ (γ_anti = 0), however non-uniform γ_sym; F71 needs palindromic γ, not uniform γ");
            yield return new InspectableNode("graceful breakdown",
                summary: "Taylor of D in γ_anti has odd powers only → leading-order linear in the per-site asymmetry; no hard violation");
            yield return new InspectableNode("γ_sym-dependence (Tier 2 empirical)",
                summary: "the leading coefficient κ_γ depends on γ_sym; parity fixes the oddness, not κ_γ; κ_γ has no closed form (∂c₁/∂γ, c₁ itself a non-closed numerical fit; see PROOF_F101)");
            yield return new InspectableNode("scope: c₁ only",
                summary: "covers c₁ (closure-breaking coefficient); F86c Q_peak is not covered (its Q-axis needs a scalar γ₀); the γ_avg-anchored Q_peak route is a separable extension");
            yield return new InspectableNode("twin of F91",
                summary: "observable-side counterpart of F91 (diagonal-block spectrum depends only on γ_sym); J-side counterpart is F100");
            yield return new InspectableNode("verified",
                summary: "c₁ via the α-rescaling pipeline, ψ_1+vac / ψ_2+vac, N=3,4,5: oddness residual ≤ 5e-9, palindromic survival ≤ 5e-9 (simulations/f71_nonuniform_gamma_verification.py)");
        }
    }
}
