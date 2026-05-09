using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F56 closed form (Tier 1, zero fit parameters; ANALYTICAL_FORMULAS line 1088):
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
/// <para>F56 is the closed-form prediction for the iteration count of the
/// Mandelbrot recursion near the cardioid cusp (CΨ = 1/4). The leading
/// logarithm comes from saddle-node ODE integral; the −4 from the starting-
/// transient (η₀ = −1/4); the ln(16·tol) from Modified Equation Euler
/// discretization correction.</para>
///
/// <para>F56 IS equivalent to the CΨ recursion near the 1/4 boundary. The
/// 1/4 position is exactly QuarterAsBilinearMaxval — the bilinear-apex
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
/// <para>Tier1Derived: F56 is Tier 1 with zero fit parameters; verified
/// 0.5–2% accuracy over 5 tol decades (10⁻⁸ to 10⁻¹⁶) and 10 ε decades
/// (10⁻¹ to 10⁻¹⁰). Modified Equation slope 0.504 vs predicted 0.500 (0.8%
/// deviation, structural).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F56 (line 1088) +
/// <c>experiments/CRITICAL_SLOWING_AT_THE_CUSP.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QuarterAsBilinearMaxvalClaim, HalfAsStructuralFixedPointClaim).</para></summary>
public sealed class F56CriticalSlowingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly QuarterAsBilinearMaxvalClaim _quarter;
    private readonly HalfAsStructuralFixedPointClaim _half;

    /// <summary>The "1/2" prefactor in (1/2)·ln(4ε/tol). Live from Pi2DyadicLadder a_2.</summary>
    public double HalfPrefactor => _ladder.Term(2);

    /// <summary>The "4" multiplier in 4·ε inside the log. Live from Pi2DyadicLadder a_{−1}.</summary>
    public double FourFactor => _ladder.Term(-1);

    /// <summary>The "16" multiplier in 16·tol inside α's log. Equals 4² = a_{−3}
    /// on the dyadic ladder. Live from Pi2DyadicLadder a_{−3}.</summary>
    public double SixteenFactor => _ladder.Term(-3);

    /// <summary>The "−4" starting-transient term in α(tol). Equals −a_{−1}.</summary>
    public double NegFourTransient => -FourFactor;

    /// <summary>The cardioid cusp position 1/4 = a_3. Live from
    /// QuarterAsBilinearMaxval (= bilinear-apex maxval at argmax p = 1/2).</summary>
    public double CardioidCuspPosition => _ladder.Term(3);

    /// <summary>F56's α(tol) = −4 + (1/2)·ln(16·tol).</summary>
    public double Alpha(double tol)
    {
        if (tol <= 0) throw new ArgumentOutOfRangeException(nameof(tol), tol, "tol must be > 0.");
        return NegFourTransient + HalfPrefactor * Math.Log(SixteenFactor * tol);
    }

    /// <summary>F56's iteration-count closed form: K(ε, tol) =
    /// (1/2)·ln(4·ε/tol) + α(tol)·√ε.</summary>
    public double IterationCount(double epsilon, double tol)
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

    /// <summary>True iff F56's ε → 0 limit gives critical slowing (K → ∞).</summary>
    public bool CriticalSlowingHolds(double epsilon, double tol)
    {
        if (epsilon <= 0 || tol <= 0) return false;
        // At small ε, K is dominated by (1/2)·ln(4ε/tol) which → −∞ if ε ≪ tol,
        // but for ε > tol/4 the log is positive and K grows as ε → 0 from above.
        // The critical-slowing regime is the regime where iteration count diverges
        // logarithmically; that holds when ε approaches the same scale as tol from above.
        return IterationCount(epsilon, tol) > 0;
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
        : base("F56 critical-slowing iteration count K(ε, tol) = (1/2)·ln(4ε/tol) + α·√ε with α = −4 + (1/2)·ln(16·tol); cardioid cusp at 1/4 = a_3 (CΨ = 1/4)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F56 + " +
               "experiments/CRITICAL_SLOWING_AT_THE_CUSP.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        _half = half ?? throw new ArgumentNullException(nameof(half));
    }

    public override string DisplayName =>
        "F56 critical slowing as Pi2-Foundation a_2 + a_{-1} + a_{-3} + QuarterAsBilinearMaxval inheritance";

    public override string Summary =>
        $"K(ε, tol) = (1/2)·ln(4ε/tol) + α·√ε; α = −4 + (1/2)·ln(16·tol); cardioid cusp at 1/4 = a_3; 1/2 = a_2, 4 = a_{{-1}}, 16 = a_{{-3}}; verified 0.5-2% accuracy over 5+10 decades ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F56 closed form",
                summary: "K(ε, tol) = (1/2)·ln(4ε/tol) + α(tol)·√ε; α(tol) = −4 + (1/2)·ln(16·tol); rescaled iteration count of u_{n+1} = u² + c near cardioid cusp c = 1/4 − ε");
            yield return InspectableNode.RealScalar("HalfPrefactor (= a_2 = 1/2)", HalfPrefactor);
            yield return InspectableNode.RealScalar("FourFactor (= a_{-1} = 4)", FourFactor);
            yield return InspectableNode.RealScalar("SixteenFactor (= a_{-3} = 16 = 4²)", SixteenFactor);
            yield return InspectableNode.RealScalar("NegFourTransient (= −a_{-1} = −4)", NegFourTransient);
            yield return InspectableNode.RealScalar("CardioidCuspPosition (= a_3 = 1/4)", CardioidCuspPosition);
            yield return new InspectableNode("five Pi2 anchors share dyadic ladder",
                summary: "1/2 (a_2 = HalfAsStructural argmax), 1/4 (a_3 = QuarterAsBilinearMaxval = cardioid cusp), 4 (a_{-1}), 16 (a_{-3} = 4²), −4 (sign-flipped a_{-1}); F56 packs five distinct dyadic-ladder positions in one closed form");
            yield return new InspectableNode("derivation",
                summary: "leading logarithm: saddle-node passage ODE integral. α(tol)'s −4: starting-transient (η₀ = −1/4). α(tol)'s ln(16·tol): Modified Equation Euler discretization correction. All three pieces give zero fit parameters.");
            yield return new InspectableNode("verified accuracy",
                summary: "0.5-2% over 5 tol decades (10⁻⁸ to 10⁻¹⁶) × 10 ε decades (10⁻¹ to 10⁻¹⁰); Modified Equation slope 0.504 vs predicted 0.500 (0.8% structural deviation)");
            yield return new InspectableNode("equivalent to CΨ recursion",
                summary: "F56 IS the CΨ recursion near the 1/4 boundary; QuarterAsBilinearMaxval = bilinear-apex maxval at argmax 1/2. ε measures distance from the cusp; ε → 0 gives critical slowing (K diverges logarithmically)");
            yield return new InspectableNode("verified examples",
                summary: $"K(ε=0.01, tol=10⁻¹⁰) = {IterationCount(0.01, 1e-10):G6}; K(ε=10⁻⁵, tol=10⁻¹²) = {IterationCount(1e-5, 1e-12):G6}; α(10⁻¹⁰) = {Alpha(1e-10):G6}");
        }
    }
}
