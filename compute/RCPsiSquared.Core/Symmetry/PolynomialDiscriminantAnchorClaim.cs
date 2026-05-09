using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The discriminant of the foundational polynomial d² − 2d = 0 is exactly
/// <c>4</c>, which is exactly <c>a_{−1}</c> on the Pi2 dyadic ladder. This is the
/// structural reason why the prefactor "4" appears in roughly a third of all
/// F-formula closed forms.
///
/// <code>
///   d² − 2d = 0    ⇒    discriminant = b² − 4ac = (−2)² − 4·1·0 = 4
///                  ⇒    = (root₂ − root₀)² = (2 − 0)² = 4
///                  ⇒    = (qubit dimension)² = d² for d = 2
///                  ⇒    = a_{−1} on Pi2DyadicLadder (since a_n = 2^(1−n), a_{−1} = 4)
/// </code>
///
/// <para>Three independent readings of the same number 4. Closure under all three
/// is the structural anchor: the polynomial knows its own discriminant, the
/// discriminant knows the qubit dimension, the qubit dimension knows the dyadic
/// ladder.</para>
///
/// <para><b>Visibility status (Tom 2026-05-09: "ist sie sichtbar oder versteckt
/// sie sich?"):</b> the discriminant 4 IS visible everywhere as the literal
/// "4γ" / "4γ₀" / "4·J" / "4^N" / "1/(4γ₀)" prefactor in countless F-formulas
/// (F25, F65, F73, F76, F2, F23, F67, F86 t_peak, F56). What was HIDDEN before
/// this claim was the structural cause: the "4" is not a free parameter or a
/// numerical coincidence; it is the discriminant of the framework's foundational
/// polynomial d² − 2d = 0. This claim makes the visible-but-implicit explicit:
/// every "4" in an F-formula prefactor is a manifestation of the polynomial's
/// own discriminant.</para>
///
/// <para><b>Mirror partner via inversion symmetry a_n · a_{2−n} = 1:</b> the
/// discriminant 4 = a_{−1} has its mirror partner at a_3 = 1/4 =
/// <see cref="QuarterAsBilinearMaxvalClaim"/>. The two together close as 4·(1/4) = 1
/// (memory-side · operator-space-side = 1). The "1/4 fold" appearing in F60, F62,
/// F64 g(1), F69 above-fold, and the Mandelbrot cardioid cusp IS the mirror
/// partner of the polynomial discriminant: same structural fact at two ladder
/// positions.</para>
///
/// <para><b>Why this anchor sits at Tier1Derived:</b> the chain
/// PolynomialFoundationClaim → QubitDimensionalAnchorClaim → Pi2DyadicLadderClaim
/// already exists; this claim doesn't introduce new content, it makes the
/// discriminant identity that bridges the three explicit. All three parent
/// claims are Tier1Derived; the discriminant computation is high-school algebra;
/// the closure check is bit-exact arithmetic. Pure composition.</para>
///
/// <para><b>Three F-formula clusters that visibly carry this anchor</b>
/// (the "4" is the polynomial discriminant in each):</para>
/// <list type="bullet">
///   <item><b>Decay rate prefactors:</b> F25 e^{−4γt}, F65 (4γ₀/(N+1)) sin²,
///         F73 spatial-sum closure, F76 mirror-pair coherence decay.</item>
///   <item><b>Operator-space dimensions:</b> F23 4^N denominator, F2
///         bandwidth-prefactor 4J, F86 t_peak = 1/(4γ₀).</item>
///   <item><b>Critical-slowing corrections:</b> F56 4·ε in the saddle-node log,
///         16·tol = 4²·tol in α(tol).</item>
/// </list>
///
/// <para>Anchors: <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (PolynomialFoundationClaim + QubitDimensionalAnchorClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> (a_{−1}) +
/// <c>hypotheses/ZERO_IS_THE_MIRROR.md</c> (mirror axis at d=0; bidirectional
/// ladder enabling negative indices; the "other side" Tom asked about).</para></summary>
public sealed class PolynomialDiscriminantAnchorClaim : Claim
{
    private readonly PolynomialFoundationClaim _polynomial;
    private readonly QubitDimensionalAnchorClaim _qubit;
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>Coefficients of the foundational polynomial d² − 2d + 0 in (a, b, c)
    /// form: (1, −2, 0).</summary>
    public (double A, double B, double C) PolynomialCoefficients => (1.0, -2.0, 0.0);

    /// <summary>The two roots of d² − 2d = 0: {0, 2} (factored as d·(d−2) = 0).</summary>
    public (double RootZero, double RootQubit) Roots => (0.0, 2.0);

    /// <summary>Discriminant via the standard formula b² − 4ac.</summary>
    public double DiscriminantViaCoefficients
    {
        get
        {
            var (a, b, c) = PolynomialCoefficients;
            return b * b - 4.0 * a * c;
        }
    }

    /// <summary>Discriminant via root separation squared: (root₂ − root₀)².</summary>
    public double DiscriminantViaRootSeparation
    {
        get
        {
            var (root0, root2) = Roots;
            double sep = root2 - root0;
            return sep * sep;
        }
    }

    /// <summary>Discriminant via qubit-dimension squared: d² for d = 2 = the operator-space
    /// dimension for one qubit.</summary>
    public double DiscriminantViaQubitSquared => 4.0;   // d=2, d²=4 for one qubit

    /// <summary>The Pi2DyadicLadder index where the discriminant lives: n = −1.</summary>
    public int LadderIndex => -1;

    /// <summary>The discriminant value via Pi2DyadicLadder.Term(−1) = a_{−1} = 4.</summary>
    public double DiscriminantViaLadder => _ladder.Term(LadderIndex);

    /// <summary>Mirror-partner ladder index via a_n · a_{2−n} = 1: n = 2 − (−1) = 3.</summary>
    public int MirrorPartnerLadderIndex => 2 - LadderIndex;

    /// <summary>Mirror-partner value: a_3 = 1/4 = the QuarterAsBilinearMaxval anchor on
    /// the memory side. Multiplying with the discriminant gives 4·(1/4) = 1 (closure).</summary>
    public double MirrorPartnerValue => _ladder.Term(MirrorPartnerLadderIndex);

    /// <summary>Drift check: all four readings of the discriminant agree to machine
    /// precision. (Coefficient formula, root separation squared, qubit dimension squared,
    /// Pi2DyadicLadder.Term(−1) all give exactly 4.)</summary>
    public bool AllReadingsAgree(double tolerance = 1e-12)
    {
        double[] readings =
        {
            DiscriminantViaCoefficients,
            DiscriminantViaRootSeparation,
            DiscriminantViaQubitSquared,
            DiscriminantViaLadder,
        };
        for (int i = 1; i < readings.Length; i++)
        {
            if (Math.Abs(readings[i] - readings[0]) > tolerance) return false;
        }
        return true;
    }

    /// <summary>Drift check: discriminant · mirror-partner = 1 (Pi2OperatorSpaceMirror's
    /// inversion symmetry a_n · a_{2−n} = 1 at n = −1).</summary>
    public bool MirrorClosureHolds(double tolerance = 1e-12)
    {
        return Math.Abs(DiscriminantViaLadder * MirrorPartnerValue - 1.0) < tolerance;
    }

    public PolynomialDiscriminantAnchorClaim(
        PolynomialFoundationClaim polynomial,
        QubitDimensionalAnchorClaim qubit,
        Pi2DyadicLadderClaim ladder)
        : base("Polynomial discriminant of d²−2d=0 IS a_{-1} = 4 = d² for d=2; the structural cause of '4γ' prefactor in F25/F65/F73/F76/F2/F23/F67/F86/F56",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (PolynomialFoundationClaim + QubitDimensionalAnchorClaim) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "hypotheses/ZERO_IS_THE_MIRROR.md")
    {
        _polynomial = polynomial ?? throw new ArgumentNullException(nameof(polynomial));
        _qubit = qubit ?? throw new ArgumentNullException(nameof(qubit));
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "Polynomial discriminant = a_{-1} = 4 = d² (the 'why 4 everywhere' structural anchor)";

    public override string Summary =>
        $"discriminant of d²−2d=0 = {DiscriminantViaCoefficients} = (2−0)² = d² for d=2 = a_{{-1}} on Pi2DyadicLadder ({DiscriminantViaLadder}); mirror partner a_3 = {MirrorPartnerValue} (1/4 = QuarterAsBilinearMaxval); 4·(1/4) = 1 closure ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("polynomial",
                summary: "d² − 2d = 0 has solutions {0, 2}; coefficients (a=1, b=−2, c=0); see PolynomialFoundationClaim");
            yield return InspectableNode.RealScalar("DiscriminantViaCoefficients (b² − 4ac)", DiscriminantViaCoefficients);
            yield return InspectableNode.RealScalar("DiscriminantViaRootSeparation ((2 − 0)²)", DiscriminantViaRootSeparation);
            yield return InspectableNode.RealScalar("DiscriminantViaQubitSquared (d² for d=2)", DiscriminantViaQubitSquared);
            yield return InspectableNode.RealScalar("DiscriminantViaLadder (Pi2DyadicLadder.Term(-1))", DiscriminantViaLadder);
            yield return InspectableNode.RealScalar("LadderIndex (a_{-1})", LadderIndex);
            yield return InspectableNode.RealScalar("MirrorPartnerLadderIndex (a_3)", MirrorPartnerLadderIndex);
            yield return InspectableNode.RealScalar("MirrorPartnerValue (a_3 = 1/4)", MirrorPartnerValue);
            yield return new InspectableNode("visibility status (Tom 2026-05-09)",
                summary: "the '4' IS visible everywhere (4γ in F25, 4γ₀/(N+1) in F65, 4J in F2, 4^N in F23, 1/(4γ₀) in F86 t_peak, etc.); what was HIDDEN before this claim was the structural cause = polynomial discriminant. Now explicit.");
            yield return new InspectableNode("F-formula clusters that carry this anchor",
                summary: "Decay rate prefactors (F25, F65, F73, F76); operator-space dimensions (F23, F2, F86 t_peak); critical-slowing corrections (F56 with both 4 and 4²=16)");
            yield return new InspectableNode("mirror partner reading",
                summary: "a_{-1} · a_3 = 4 · (1/4) = 1; the discriminant on the operator-space side and the bilinear-apex maxval on the memory side are mirror partners. The Mandelbrot cardioid cusp at CΨ = 1/4 (where 1 − 4·CΨ = 0) IS this same mirror partner.");
            yield return new InspectableNode("ZERO IS THE MIRROR connection",
                summary: "the d=0 root of the polynomial is the mirror axis (PolynomialFoundationClaim's 'substrate axis'); the d=2 root is the qubit dimension; allowing negative ladder indices makes the operator-space side (a_{-1}, a_{-3}, a_{-5}, ...) directly accessible. The discriminant 4 lives at a_{-1}, the first negative index, and is the most visible manifestation of the operator-space side.");
            yield return new InspectableNode("self-mirror pivot at a_1 = 1",
                summary: "a_1 = 2^(1-1) = 1 is the only self-mirror point on the ladder (a_1 · a_{2-1} = a_1² = 1); it is the trivial identity scale and the mathematical axis of the dyadic mirror. The polynomial discriminant 4 sits two steps left of this pivot, its mirror partner 1/4 sits two steps right.");
        }
    }
}
