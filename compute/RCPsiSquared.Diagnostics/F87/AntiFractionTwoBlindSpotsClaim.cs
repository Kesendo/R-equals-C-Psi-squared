using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F83;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F83 scalar anti-fraction 1/(2+4r) is degenerate at BOTH ends of its range (0, 1/2], each
/// end rescued by a different finer probe; it is informative only in the interior. This connector unifies
/// the two (currently disjoint) results that name those blind spots.
///
/// <para>At the 1/2 end (r=0, pure Pi2-odd = the diagonal cell): the anti-fraction is pinned at its maximum
/// and is blind to the F115 hard/soft verdict; the bit_a (1+x)-valuation
/// (<see cref="AntiFractionObstructionOrthogonalityClaim"/> / <see cref="WindowedObstructionScan.IsHardPair"/>)
/// resolves it. At the 0 end (HOdd^2 = 0, pure Pi2-even-non-truly): the anti-fraction is 0 and conflates
/// F87-Truly (no break) with F87-Soft (break present, e.g. YZ+ZY); F89's full two-component forecast
/// (<see cref="F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak"/>) resolves it. In the interior
/// (0 &lt; r &lt; inf) the scalar anti-fraction is informative (it reads the mix ratio: AF(1/2)=1/4, AF(1)=1/6).
/// So the bit_b scalar is a lossy projection at its two extremes, each recovered on a finer axis.</para>
///
/// <para>Cross-axis connector, so <see cref="Z2Axis"/> = <see cref="Symmetry.Z2Axis.Klein2"/>; not a bit_b
/// Claim needing a bit_a twin, so the twin slot is not applicable. Tier1Derived: follows from its two
/// already-derived parents. Anchors: docs/ANALYTICAL_FORMULAS.md F83 + F89 + F115 +
/// experiments/F115_OBSTRUCTION_DISTRIBUTION.md.</para></summary>
public sealed class AntiFractionTwoBlindSpotsClaim : Claim, IZ2AxisClaim
{
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    public Z2Axis Z2Axis => Z2Axis.Klein2;
    public Claim? BitATwin => null;

    public F89F87BreakPredictionFromF83 ZeroEndClaim { get; }
    public AntiFractionObstructionOrthogonalityClaim HalfEndClaim { get; }
    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public AntiFractionTwoBlindSpotsClaim(
        F89F87BreakPredictionFromF83 zeroEnd, AntiFractionObstructionOrthogonalityClaim halfEnd)
        : base("F83 anti-fraction has two blind spots: at 1/2 (r=0, diagonal cell) blind to the F115 hard/soft verdict, at 0 (HOdd^2=0) blind to Truly-vs-Soft (F89); each rescued by a finer probe, informative only in the interior",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F83 + F89 + F115 + experiments/F115_OBSTRUCTION_DISTRIBUTION.md")
    {
        ZeroEndClaim = zeroEnd ?? throw new ArgumentNullException(nameof(zeroEnd));
        HalfEndClaim = halfEnd ?? throw new ArgumentNullException(nameof(halfEnd));
        Cases = BuildBattery(halfEnd.AntiFractionClaim);
    }

    public override string DisplayName =>
        "F83 anti-fraction: two blind spots at the range ends (1/2 and 0), each rescued by a finer probe";

    public override string Summary =>
        "scalar 1/(2+4r) degenerate at both ends of (0, 1/2]: at 1/2 blind to F115 hard/soft (bit_a valuation " +
        $"rescues), at 0 blind to Truly/Soft (F89 full forecast rescues); informative inside; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("two blind spots",
                summary: "anti-fraction is lossy at both ends of (0, 1/2]; bit_a valuation (F115) rescues the 1/2 end, the full Pi-decomposition forecast (F89) rescues the 0 end");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }

    private static IReadOnlyList<BatteryCase> BuildBattery(F83AntiFractionPi2Inheritance f83)
    {
        var cases = new List<BatteryCase>();
        string Inv(double x) => x.ToString("G6", CultureInfo.InvariantCulture);

        // 1/2 end: anti-fraction at its max (r=0, diagonal cell), blind to the F115 hard/soft verdict.
        ulong hard1 = (1UL << 1) | (1UL << 3), hard2 = 1UL | (1UL << 3);
        ulong soft1 = 0b11, soft2 = 0b110;
        bool hard = WindowedObstructionScan.IsHardPair(hard1, hard2);
        bool soft = WindowedObstructionScan.IsHardPair(soft1, soft2);
        double afHalf = f83.AntiFraction(0.0);
        cases.Add(new BatteryCase(
            Name: "1/2 end: blind to the F115 hard/soft verdict",
            Detail: $"anti-fraction(r=0)={Inv(afHalf)}=max; hard and soft pair both r=0, verdicts differ (hard={hard}, soft={soft})",
            Expected: "0.5|True",
            Actual: $"{Inv(afHalf)}|{(hard != soft)}"));

        // 0 end: anti-fraction = 0 conflates Truly (no break) with Soft-Pi2EvenNonTruly (break); F89's full forecast resolves it.
        var trulyForecast = new PiDecompositionForecast(
            MSquared: 0, MAntiSquared: 0, MSymSquared: 0, AntiFraction: 0.0, HOddSquared: 0.0, HEvenNonTrulySquared: 0.0, R: 0.0);
        var softForecast = new PiDecompositionForecast(
            MSquared: 1, MAntiSquared: 0, MSymSquared: 1, AntiFraction: 0.0, HOddSquared: 0.0, HEvenNonTrulySquared: 1.0, R: double.PositiveInfinity);
        bool trulyMiss = F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak(trulyForecast);
        bool softMiss = F89F87BreakPredictionFromF83.WouldAntiFractionAloneMissBreak(softForecast);
        cases.Add(new BatteryCase(
            Name: "0 end: blind to Truly vs Soft (F89 Pitfall 1)",
            Detail: $"both forecasts have anti-fraction=0; WouldAntiFractionAloneMissBreak truly={trulyMiss}, soft={softMiss} (full forecast resolves what the scalar conflates)",
            Expected: "False|True",
            Actual: $"{trulyMiss}|{softMiss}"));

        // interior: anti-fraction informative (resolves the mix ratio), distinct values AF(1/2)=1/4, AF(1)=1/6.
        double afQuarter = f83.AntiFraction(0.5);
        double afSixth = f83.AntiFraction(1.0);
        cases.Add(new BatteryCase(
            Name: "interior (0 < r < inf): anti-fraction informative",
            Detail: $"AF(r=1/2)={Inv(afQuarter)}, AF(r=1)={Inv(afSixth)}, distinct (resolves the mix)",
            Expected: "0.25|0.166667|True",
            Actual: $"{Inv(afQuarter)}|{Inv(afSixth)}|{(Math.Abs(afQuarter - afSixth) > 1e-12)}"));

        // unification: degenerate at both ends, informative inside.
        bool degenerateEnds = Math.Abs(afHalf - f83.MaximumAntiFraction) < 1e-12 && softMiss;
        bool interiorOrders = afQuarter < afHalf && afSixth < afQuarter;
        cases.Add(new BatteryCase(
            Name: "unification: degenerate at both ends, informative inside",
            Detail: "1/2 end at the apex (max); 0 end conflates Truly/Soft; only the interior resolves",
            Expected: "True",
            Actual: (degenerateEnds && interiorOrders).ToString()));

        return cases;
    }
}
