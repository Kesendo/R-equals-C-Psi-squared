using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Core.F86;

/// <summary>F86 per-block Q_peak (Q_SCALE convention, relative-J derivative ΔJ = 0.05·J):
/// Tier 2 empirical chromaticity-specific N-invariant constants.
///
/// <para>Q_peak(c=3) = 1.6, Q_peak(c=4) = 1.8, Q_peak(c=5) = 1.8 — saturates at 1.8 for c ≥ 4.
/// c=2 is finite-size-sensitive (1.4..1.6 across N=4..9) and not a clean framework constant.</para>
///
/// <para>Verified across: c=3 N∈{5,6,7,8,9}, c=4 N∈{7,8,9}, c=5 at N=9 (commit 4612468 in
/// Q_SCALE_THREE_BANDS). Used by the γ₀-extraction protocol: <c>γ₀ ≈ J*/Q_peak(c)</c>.</para>
/// </summary>
public sealed class PerBlockQPeakClaim : F86Claim
{
    public int Chromaticity { get; }
    public double QPeakValue { get; }
    public string VerifiedRange { get; }
    public bool Saturated { get; }
    public string? Caveat { get; }

    public PerBlockQPeakClaim(int chromaticity, double qPeak, string verifiedRange,
        bool saturated, string? caveat = null)
        : base($"per-block Q_peak (c={chromaticity})",
               caveat is null ? Tier.Tier2Empirical : Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F86 + experiments/Q_SCALE_THREE_BANDS.md")
    {
        Chromaticity = chromaticity;
        QPeakValue = qPeak;
        VerifiedRange = verifiedRange;
        Saturated = saturated;
        Caveat = caveat;
    }

    /// <summary>γ₀ extraction: given measured peak J*, what is γ₀ predicted by this block's Q_peak?</summary>
    public double ExtractGammaZero(double measuredJStar) => measuredJStar / QPeakValue;

    public override string DisplayName => $"Q_peak(c={Chromaticity}) = {QPeakValue:F2}";
    public override string Summary =>
        $"= {QPeakValue:F2} ({Tier.Label()}, verified at {VerifiedRange}" +
        (Caveat is null ? ")" : $", caveat: {Caveat})");

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("chromaticity", Chromaticity);
            yield return InspectableNode.RealScalar("Q_peak (per-block, Q_SCALE)", QPeakValue, "F2");
            yield return new InspectableNode("verified range", summary: VerifiedRange);
            yield return new InspectableNode("saturated?", summary: Saturated ? "yes" : "no");
            if (Caveat is not null)
                yield return new InspectableNode("caveat", summary: Caveat);
        }
    }

    public override InspectablePayload Payload =>
        new InspectablePayload.Real($"Q_peak(c={Chromaticity})", QPeakValue, "F2");

    /// <summary>The four standard per-block Q_peak claims from Q_SCALE_THREE_BANDS.</summary>
    public static IReadOnlyList<PerBlockQPeakClaim> Standard { get; } = new[]
    {
        new PerBlockQPeakClaim(2, 1.5, "N=4..9 (wobbles 1.4..1.6)", saturated: false,
            caveat: "finite-size sensitive, not a clean framework constant"),
        new PerBlockQPeakClaim(3, 1.6, "N=5..9", saturated: true),
        new PerBlockQPeakClaim(4, 1.8, "N=7..9", saturated: true),
        new PerBlockQPeakClaim(5, 1.8, "N=9 (commit 4612468)", saturated: true),
    };
}
