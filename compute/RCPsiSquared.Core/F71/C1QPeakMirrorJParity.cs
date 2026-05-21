using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>F100 (Tier 1 derived, 2026-05-20): the F71 c₁/Q_peak bond-mirror deviation
///
/// <code>
///     D(b)  :=  c₁(b) − c₁(N−2−b)
/// </code>
///
/// is an EXACTLY ODD function of the F71-anti-palindromic component of the bond-coupling
/// profile J. F71 (<see cref="C1MirrorIdentity"/>) was proven for uniform J; this is its
/// graceful-breakdown extension to non-uniform J, and closes the F71 "non-uniform J_b"
/// open question.
///
/// <para>Decompose J into its F71-palindromic and F71-anti-palindromic parts,
/// J = J_sym + J_anti with J_sym = (J + F71(J))/2 and J_anti = (J − F71(J))/2 (F71 maps
/// bond b ↔ N−2−b). Then:</para>
///
/// <list type="bullet">
///   <item>D ≡ 0 for every palindromic J (J_anti = 0), however non-uniform J_sym is. F71
///         never required uniform J; it requires palindromic J. Uniform is merely the
///         simplest palindromic profile.</item>
///   <item>D(J_sym, −J_anti) = −D(J_sym, J_anti): exactly odd in J_anti, all orders. The
///         Taylor series of D in J_anti has odd powers only, so the breakdown is graceful
///         (leading-order linear in B_b = J_b − J_{N−2−b}), never a hard jump.</item>
///   <item>The leading coefficient κ_b generically DEPENDS on J_sym (Tier 2 empirical):
///         the parity argument fixes the oddness, not the coefficient. κ_b admits no
///         closed form: it is ∂c₁/∂J and c₁ is itself a non-closed numerical fit. See the
///         κ-obstruction section of <c>docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md</c>.</item>
/// </list>
///
/// <para>Mechanism (cf. <c>docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md</c>): the F71
/// chain-mirror R is R-equivariant on the whole c₁ pipeline, so c₁(b; J) = c₁(N−2−b; F71(J));
/// hence D(b; F71(J)) = −D(b; J). Substituting F71(J) = J_sym − J_anti (F92's J_sym/J_anti
/// split, H linear in J) yields the oddness. This is the observable-side twin of F92: F92
/// says the F71-refined diagonal-block spectrum depends only on J_sym; F100 says the
/// c₁/Q_peak bond-mirror deviation lives in J_anti.</para>
///
/// <para>Scope: covers c₁ (closure-breaking coefficient, vac+SE) and Q_peak (F86c per-bond
/// observable; the identical F71-conjugation argument). Both numerically witnessed with
/// residuals ≤ 1e-9: c₁ at N=3,4,5 (<c>simulations/_f71_nonuniform_j_verification.py</c>),
/// Q_peak at N=4,5,6 (<c>simulations/_f100_qpeak_nonuniform_j_verification.py</c>).</para>
/// </summary>
public sealed class C1QPeakMirrorJParity : Claim
{
    public C1QPeakMirrorJParity()
        : base("F100 c₁/Q_peak mirror-deviation parity: D(b) = c₁(b) − c₁(N−2−b) is odd in J_anti, zero for palindromic J",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F100 + docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md")
    { }

    /// <summary>F71-palindromic component J_sym = (J + F71(J))/2, where F71 maps bond
    /// b ↔ N−2−b. Palindromic by construction (J_sym[b] = J_sym[N−2−b]). The c₁/Q_peak
    /// bond-mirror deviation D is blind to this component.</summary>
    public static double[] PalindromicComponent(IReadOnlyList<double> bondJ)
    {
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        int n = bondJ.Count;
        var jSym = new double[n];
        for (int b = 0; b < n; b++)
            jSym[b] = 0.5 * (bondJ[b] + bondJ[n - 1 - b]);
        return jSym;
    }

    /// <summary>F71-anti-palindromic component J_anti = (J − F71(J))/2. Anti-palindromic
    /// by construction (J_anti[b] = −J_anti[N−2−b]); the deviation D is exactly odd in it,
    /// and zero iff J_anti is zero.</summary>
    public static double[] AntiPalindromicComponent(IReadOnlyList<double> bondJ)
    {
        if (bondJ is null) throw new ArgumentNullException(nameof(bondJ));
        int n = bondJ.Count;
        var jAnti = new double[n];
        for (int b = 0; b < n; b++)
            jAnti[b] = 0.5 * (bondJ[b] - bondJ[n - 1 - b]);
        return jAnti;
    }

    /// <summary>Norm ‖J_anti‖ = sqrt(Σ_b J_anti[b]²) of the anti-palindromic component:
    /// the deviation of J from being F71-palindromic. Zero iff J_b = J_{N−2−b} for all b.
    /// The c₁/Q_peak bond-mirror deviation D vanishes iff this is zero; it is the palindromic-axis
    /// twin of F92's AntiPalindromicDeviation.</summary>
    public static double PalindromicDeviation(IReadOnlyList<double> bondJ)
    {
        double sum = 0.0;
        foreach (double jAnti in AntiPalindromicComponent(bondJ))
            sum += jAnti * jAnti;
        return Math.Sqrt(sum);
    }

    /// <summary>True iff J is F71-palindromic (J_b = J_{N−2−b} ∀b) within tolerance: the
    /// regime in which the c₁/Q_peak bond-mirror identity D(b) = 0 holds exactly.</summary>
    public static bool IsPalindromic(IReadOnlyList<double> bondJ, double tolerance = 1e-10)
        => PalindromicDeviation(bondJ) < tolerance;

    public override string DisplayName =>
        "F100: D(b) = c₁(b) − c₁(N−2−b) is odd in J_anti (graceful breakdown of F71)";

    public override string Summary =>
        "F71 c₁/Q_peak bond-mirror deviation is exactly odd in the anti-palindromic J-component; "
        + $"zero for palindromic J, leading-order linear; c₁ numerically verified N=3,4,5 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("parity mechanism",
                summary: "R-equivariance gives c₁(b;J) = c₁(N−2−b;F71(J)); hence D(b;F71(J)) = −D(b;J), exactly odd in J_anti (all orders)");
            yield return new InspectableNode("palindromic survival",
                summary: "D ≡ 0 for every palindromic J (J_anti = 0), however non-uniform J_sym; F71 needs palindromic J, not uniform J");
            yield return new InspectableNode("graceful breakdown",
                summary: "Taylor of D in J_anti has odd powers only → leading-order linear in B_b = J_b − J_{N−2−b}; no hard violation");
            yield return new InspectableNode("J_sym-dependence (Tier 2 empirical)",
                summary: "the leading coefficient κ_b depends on J_sym; parity fixes the oddness, not κ; 62–143% relative κ spread across J_sym profiles, N=3,4,5; κ has no closed form (∂c₁/∂J, c₁ itself a non-closed numerical fit; see PROOF_F100)");
            yield return new InspectableNode("scope: c₁ and Q_peak",
                summary: "covers c₁ (closure-breaking coefficient) and Q_peak (F86c per-bond observable); the identical F71-conjugation argument");
            yield return new InspectableNode("twin of F92",
                summary: "observable-side counterpart of F92 (diagonal-block spectrum depends only on J_sym); F92 and F100 are the two faces of the J_sym/J_anti split");
            yield return new InspectableNode("verified",
                summary: "c₁ via the α-rescaling pipeline, ψ_1+vac / ψ_2+vac, N=3,4,5: oddness residual ≤ 1.0e-9, palindromic survival ≤ 4.0e-10 (simulations/_f71_nonuniform_j_verification.py)");
        }
    }
}
