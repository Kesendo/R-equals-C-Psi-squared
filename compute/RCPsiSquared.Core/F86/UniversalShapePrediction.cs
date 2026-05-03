using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 Statement 2 (Tier 1 candidate, 2026-05-02): per bond class, the K-curve
/// HWHM_left/Q_peak ratio is universal — Interior ≈ 0.756, Endpoint ≈ 0.770 across the
/// tested (c, N, γ₀) grid.
///
/// <para>Each prediction carries:</para>
/// <list type="bullet">
///   <item>The bond <see cref="BondClass"/> it applies to.</item>
///   <item>The expected ratio (the named structural fingerprint).</item>
///   <item>A tolerance (the empirical scatter across witnesses).</item>
///   <item>A list of <see cref="UniversalShapeWitness"/> data points that establish the
///         Tier-1-candidate status (numerical multi-(c, N) stability).</item>
/// </list>
///
/// <para>Methodological caveat (per the PTF retraction lesson): the ratio is a numerical
/// witness, not a derived constant. Promotion to Tier 1 requires deriving f_class(x) from
/// the 2-level EP analytics. Candidate closed forms (3/4, 2/√7, …) must come from the
/// algebra, not from least-squares fits.</para>
/// </summary>
public sealed class UniversalShapePrediction : F86Claim
{
    public BondClass BondClass { get; }
    public double ExpectedHwhmOverQPeak { get; }
    public double Tolerance { get; }
    public IReadOnlyList<UniversalShapeWitness> Witnesses { get; }

    public UniversalShapePrediction(BondClass bondClass, double expectedRatio, double tolerance,
        IReadOnlyList<UniversalShapeWitness> witnesses)
        : base($"universal HWHM/Q for {bondClass}",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F86 + docs/proofs/PROOF_F86_QPEAK.md Statement 2")
    {
        BondClass = bondClass;
        ExpectedHwhmOverQPeak = expectedRatio;
        Tolerance = tolerance;
        Witnesses = witnesses;
    }

    /// <summary>Compare a measured peak against this prediction. <c>Within</c> is true iff
    /// |measured − expected| ≤ tolerance.</summary>
    public PredictionMatch CompareTo(PeakResult measured)
    {
        if (measured.HwhmLeftOverQPeak is not { } actual)
            return new PredictionMatch(this, null, false, "no HWHM available");
        double delta = actual - ExpectedHwhmOverQPeak;
        bool within = Math.Abs(delta) <= Tolerance;
        return new PredictionMatch(this, actual, within,
            $"actual = {actual:F4}, expected = {ExpectedHwhmOverQPeak:F4}, Δ = {Formatting.SignedDelta(delta)}");
    }

    public override string DisplayName => $"{BondClass} HWHM_left/Q_peak ≈ {ExpectedHwhmOverQPeak:F4}";

    public override string Summary =>
        $"= {ExpectedHwhmOverQPeak:F4} ± {Tolerance:F4} ({Tier.Label()}, {Witnesses.Count} witness(es))";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return InspectableNode.RealScalar("expected HWHM_left/Q_peak", ExpectedHwhmOverQPeak, "F4");
            yield return InspectableNode.RealScalar("tolerance", Tolerance, "F4");
            yield return InspectableNode.Group("witnesses",
                Witnesses.Cast<IInspectable>().ToArray());
        }
    }

    public override InspectablePayload Payload =>
        new InspectablePayload.Real("HWHM_left/Q_peak (predicted)", ExpectedHwhmOverQPeak, "F4");
}

/// <summary>Result of comparing a measured <see cref="PeakResult"/> to a
/// <see cref="UniversalShapePrediction"/>. Use as the witness side of the OOP knowledge —
/// "here is what F86 says vs. here is what we measured today".</summary>
public sealed record PredictionMatch(
    UniversalShapePrediction Prediction,
    double? Actual,
    bool Within,
    string Description) : IInspectable
{
    public string DisplayName => $"vs measured ({(Within ? "within tolerance" : "OUTSIDE tolerance")})";
    public string Summary => Description;
    public IEnumerable<IInspectable> Children { get; } = Array.Empty<IInspectable>();
    public InspectablePayload Payload => Actual is { } a
        ? new InspectablePayload.Real("measured HWHM_left/Q_peak", a, "F4")
        : InspectablePayload.Empty;
}
