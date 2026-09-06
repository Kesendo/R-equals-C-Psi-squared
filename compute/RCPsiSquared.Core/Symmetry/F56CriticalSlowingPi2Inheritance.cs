using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F56 asymptotic expansion (Tier 1 derived, zero fit parameters; ANALYTICAL_FORMULAS F56 entry):
///
/// <code>
///   K(ε, tol) = (1/2)·ln(4·ε/tol) + α(tol)·√ε
///
///   α(tol)   = −4 + (1/2)·ln(16·tol)
///
///   K = n·√ε is the rescaled iteration count of u_{n+1} = u² + c
///   near the cardioid cusp at c = 1/4 − ε.
/// </code>
///
/// <para>F56 is an asymptotic prediction for the iteration count of the
/// Mandelbrot recursion near the cardioid cusp (CΨ = 1/4). The leading
/// logarithm comes from saddle-node ODE integral; the −4 from the starting-
/// transient (η₀ = −1/4); the ln(16·tol) from Modified Equation Euler
/// discretization correction. It has a nonzero finite-ε residual; zero fitted
/// coefficients does not make it an exact finite-ε count.</para>
///
/// <para>F56 IS equivalent to the CΨ recursion near the 1/4 boundary. The
/// 1/4 position is exactly QuarterAsBilinearMaxval: the bilinear-apex
/// maxval where p(1−p) is maximised at p = 1/2 (HalfAsStructuralFixedPoint).
/// Both Pi2-Foundation anchors appear in F56 by construction.</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>HalfPrefactor = 1/2 = a_2</b>: in (1/2)·ln(4ε/tol). Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(2). Same anchor as
///         <see cref="HalfAsStructuralFixedPointClaim"/>; the cardioid argmax
///         is at p = 1/2.</item>
///   <item><b>FourFactor = 4 = a_{−1}</b>: in 4·ε. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same anchor as F25
///         decay rate, F73 spatial-sum closure, F65 numerator.</item>
///   <item><b>SixteenFactor = 16 = a_{−3} = 4²</b>: in 16·tol. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−3). Square of FourFactor
///         (= 4² in the Modified Equation correction).</item>
///   <item><b>NegFourTransient = −4 = −a_{−1}</b>: in α(tol)'s starting-
///         transient term. Same anchor as FourFactor with sign flip.</item>
///   <item><b>CardioidCuspPosition = 1/4 = a_3</b>: cusp at c = 1/4 − ε.
///         Live from <see cref="QuarterAsBilinearMaxvalClaim"/>; the
///         distance ε measures how far from the QuarterAsBilinearMaxval
///         maxval the recursion sits. ε → 0 gives critical slowing.</item>
/// </list>
///
/// <para>Tier1Derived: F56 is Tier 1 with zero fit parameters. The correction
/// coefficient agrees at 0.5–2% over tol=10⁻⁸...10⁻¹⁶; the ε sweep records
/// nonzero corrected-K residuals. Modified Equation slope 0.504 vs predicted
/// 0.500.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F56 entry +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QuarterAsBilinearMaxvalClaim, HalfAsStructuralFixedPointClaim).</para></summary>
public sealed class F56CriticalSlowingPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    public Pi2DyadicLadderClaim Ladder { get; }
    public QuarterAsBilinearMaxvalClaim Quarter { get; }
    public HalfAsStructuralFixedPointClaim Half { get; }
    /// <summary>The "1/2" prefactor in (1/2)·ln(4ε/tol). Live from Pi2DyadicLadder a_2.</summary>
    public double HalfPrefactor => Ladder.Term(2);

    /// <summary>The "4" multiplier in 4·ε inside the log. Live from Pi2DyadicLadder a_{−1}.</summary>
    public double FourFactor => Ladder.Term(-1);

    /// <summary>The "16" multiplier in 16·tol inside α's log. Equals 4² = a_{−3}
    /// on the dyadic ladder. Live from Pi2DyadicLadder a_{−3}.</summary>
    public double SixteenFactor => Ladder.Term(-3);

    /// <summary>The "−4" starting-transient term in α(tol). Equals −a_{−1}.</summary>
    public double NegFourTransient => -FourFactor;

    /// <summary>The cardioid cusp position 1/4 = a_3. Live from
    /// QuarterAsBilinearMaxval (= bilinear-apex maxval at argmax p = 1/2).</summary>
    public double CardioidCuspPosition => Ladder.Term(3);

    /// <summary>F56's α(tol) = −4 + (1/2)·ln(16·tol).</summary>
    public double Alpha(double tol)
    {
        if (tol <= 0) throw new ArgumentOutOfRangeException(nameof(tol), tol, "tol must be > 0.");
        return NegFourTransient + HalfPrefactor * Math.Log(SixteenFactor * tol);
    }

    /// <summary>F56's iteration-count asymptotic: K(ε, tol) =
    /// (1/2)·ln(4·ε/tol) + α(tol)·√ε.</summary>
    public double IterationCountAsymptotic(double epsilon, double tol)
    {
        if (epsilon <= 0) throw new ArgumentOutOfRangeException(nameof(epsilon), epsilon, "ε must be > 0.");
        if (tol <= 0) throw new ArgumentOutOfRangeException(nameof(tol), tol, "tol must be > 0.");
        double logTerm = HalfPrefactor * Math.Log(FourFactor * epsilon / tol);
        double sqrtTerm = Alpha(tol) * Math.Sqrt(epsilon);
        return logTerm + sqrtTerm;
    }

    /// <summary>The cardioid distance: c = CΨ_cusp − ε = 1/4 − ε. Returns the c
    /// for which K(ε, tol) is the iteration count.</summary>
    public double CardioidDistance(double epsilon)
    {
        if (epsilon <= 0) throw new ArgumentOutOfRangeException(nameof(epsilon), epsilon, "ε must be > 0.");
        return CardioidCuspPosition - epsilon;
    }

    /// <summary>True only when this one finite asymptotic estimate is positive in the stated
    /// input domain. This is not a verdict that critical slowing holds: that claim concerns the
    /// scale-separated joint asymptotic behaviour of n = K/√ε.</summary>
    public bool IsAsymptoticEstimatePositiveInDomain(double epsilon, double tol)
    {
        if (epsilon <= 0 || tol <= 0) return false;
        if (epsilon >= CardioidCuspPosition || tol >= epsilon) return false;
        return IterationCountAsymptotic(epsilon, tol) > 0;
    }

    /// <summary>Drift check: SixteenFactor = FourFactor².</summary>
    public bool SixteenIsFourSquared(double tolerance = 1e-12)
    {
        return Math.Abs(SixteenFactor - FourFactor * FourFactor) < tolerance;
    }

    /// <summary>Drift check: CardioidCuspPosition matches QuarterAsBilinearMaxval (= 1/4).</summary>
    public bool CardioidCuspMatchesQuarter(double tolerance = 1e-12)
    {
        // QuarterAsBilinearMaxval doesn't expose a value property, but a_3 on the
        // dyadic ladder is the canonical 1/4 = (1/2)² anchor that defines it.
        return Math.Abs(CardioidCuspPosition - 0.25) < tolerance;
    }

    public F56CriticalSlowingPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QuarterAsBilinearMaxvalClaim quarter,
        HalfAsStructuralFixedPointClaim half)
        : base("F56 critical-slowing asymptotic K(ε, tol) = (1/2)·ln(4ε/tol) + α·√ε with α = −4 + (1/2)·ln(16·tol); nonzero finite-ε residual; cardioid cusp at 1/4 = a_3 (CΨ = 1/4)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F56 + " +
               "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        Half = half ?? throw new ArgumentNullException(nameof(half));
    }

    public override string DisplayName =>
        "F56 critical slowing as Pi2-Foundation a_2 + a_{-1} + a_{-3} + QuarterAsBilinearMaxval inheritance";

    public override string Summary =>
        $"asymptotic K(ε, tol) = (1/2)·ln(4ε/tol) + α·√ε; α = −4 + (1/2)·ln(16·tol); nonzero finite-ε residual; n = K/√ε; cardioid cusp at 1/4 = a_3; 1/2 = a_2, 4 = a_{{-1}}, 16 = a_{{-3}} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F56 asymptotic expansion",
                summary: "K(ε, tol) = (1/2)·ln(4ε/tol) + α(tol)·√ε; α(tol) = −4 + (1/2)·ln(16·tol); a zero-fit asymptotic for the rescaled iteration count near c = 1/4 − ε, not an exact finite-ε value");
            yield return InspectableNode.RealScalar("HalfPrefactor (= a_2 = 1/2)", HalfPrefactor);
            yield return InspectableNode.RealScalar("FourFactor (= a_{-1} = 4)", FourFactor);
            yield return InspectableNode.RealScalar("SixteenFactor (= a_{-3} = 16 = 4²)", SixteenFactor);
            yield return InspectableNode.RealScalar("NegFourTransient (= −a_{-1} = −4)", NegFourTransient);
            yield return InspectableNode.RealScalar("CardioidCuspPosition (= a_3 = 1/4)", CardioidCuspPosition);
            yield return new InspectableNode("five Pi2 anchors share dyadic ladder",
                summary: "1/2 (a_2 = HalfAsStructural argmax), 1/4 (a_3 = QuarterAsBilinearMaxval = cardioid cusp), 4 (a_{-1}), 16 (a_{-3} = 4²), −4 (sign-flipped a_{-1}); F56 packs five distinct dyadic-ladder positions in one asymptotic expansion");
            yield return new InspectableNode("derivation",
                summary: "leading logarithm: saddle-node passage ODE integral. α(tol)'s −4: starting-transient (η₀ = −1/4). α(tol)'s ln(16·tol): Modified Equation Euler discretization correction. All three pieces give zero fit parameters.");
            yield return new InspectableNode("finite-grid verification",
                summary: "the correction coefficient agrees at 0.5-2% over tol=10⁻⁸...10⁻¹⁶; corrected-K residuals are nonzero (−0.573, +0.037, −0.005, +0.001 at ε=10⁻¹...10⁻⁴); Modified Equation slope 0.504 vs predicted 0.500");
            yield return new InspectableNode("equivalent to CΨ recursion",
                summary: "F56 IS the auxiliary CΨ recursion near the 1/4 boundary; QuarterAsBilinearMaxval = bilinear-apex maxval at argmax 1/2. In a joint limit that keeps tol ≪ ε ≪ 1, the raw count n = K/√ε diverges; K itself is the rescaled count and has no fixed-tol positive-divergence claim");
            yield return new InspectableNode("verified examples",
                summary: $"asymptotic K(ε=0.01, tol=10⁻¹⁰) = {IterationCountAsymptotic(0.01, 1e-10):G6}; asymptotic K(ε=10⁻⁵, tol=10⁻¹²) = {IterationCountAsymptotic(1e-5, 1e-12):G6}; α(10⁻¹⁰) = {Alpha(1e-10):G6}");
        }
    }
}
